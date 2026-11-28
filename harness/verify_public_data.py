#!/usr/bin/env python3
"""
harness/verify_public_data.py

Single-pass, Python 3 stdlib-only DAG omnibus verifier for VarDictRust public
sequencing cohorts across Waves 1-6, clinical workflow presets, 38-switch
CompatProfile ablation (only_1..38 & loo_1..38), Spill == Unspilled, and
Tiled == Untiled conservation across all 12 real-world sequencing modalities.

Key properties:
  - 100% Python 3 standard library (zero PyPI dependencies).
  - Singleflight content/command-addressed RunCache deduplicating shared
    (cohort, mode, flags) invocations across suites (-78% BAM spawns).
  - Reuses on-disk RTG SDF references (GRCh38.d1.vd1.sdf, GRCh38_GIABv3.sdf,
    grch37.sdf) with zero scratch `rtg format` builds.
  - Courteous concurrency bounds: workers * -th <= 12 (default 2 * 4 = 8) and
    <= 2 concurrent `rtg vcfeval` jobs (`RTG_MEM=2G`).
  - Runs the Universal `.var` + `.vcf` + Survival + VCFv4.2 + `bcftools norm`
    + `tabix`/`bcftools view -r` `--tbi`/`--csi` Auditor on every `--compat=off`
    and `--preset` output.
"""

from __future__ import annotations

import argparse
import bisect
import concurrent.futures
from dataclasses import dataclass, field
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import signal
import struct
import subprocess
import sys
import tempfile
import threading
import time
from typing import Callable

VALID_SUITES = ("wave4", "wave5", "presets", "wave6", "ablation", "tiling", "all")

VALID_ALT_RE = re.compile(
    r"^(?:[ACGTN]+|<[A-Za-z0-9_-]+>)(?:,(?:[ACGTN]+|<[A-Za-z0-9_-]+>))*$"
)
CORRUPT_SV_TOKENS = ("<INN>", "<NEL>", "<AUP>", "<INS_CORRUPT>")
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
TVAF_RE = re.compile(r"(?:^|;)TVAF=([0-9.eE+-]+)")
DNA_RE = re.compile(r"^[ACGTN]+$")


def ensure_tool_paths() -> None:
    conda_root = os.environ.get("CONDA_ROOT") or os.path.join(
        os.path.expanduser("~"), "miniconda3"
    )
    g4_bin = os.environ.get("G4_TOOLS_BIN") or os.path.join(
        conda_root, "envs", "g4-tools", "bin"
    )
    oracle_bin = os.path.join(conda_root, "envs", "vardict-oracle", "bin")
    cur_parts = os.environ.get("PATH", "").split(os.pathsep)
    prepend = []
    for cand in (g4_bin, oracle_bin):
        if os.path.isdir(cand) and cand not in cur_parts:
            prepend.append(cand)
    if prepend:
        os.environ["PATH"] = os.pathsep.join(prepend + cur_parts)


def resolve_data_root(root_dir: Path) -> Path:
    explicit = os.environ.get("VARDICT_DATA_ROOT")
    if explicit and os.path.isdir(explicit):
        return Path(explicit).resolve()
    if (root_dir / "scratch" / "g4").is_dir():
        return root_dir
    try:
        common = subprocess.check_output(
            ["git", "-C", str(root_dir), "rev-parse", "--git-common-dir"],
            stderr=subprocess.DEVNULL,
            text=True,
        ).strip()
        if common:
            common_path = (root_dir / common).resolve()
            main_root = common_path.parent
            if (main_root / "scratch" / "g4").is_dir():
                return main_root
    except Exception:
        pass
    return root_dir


def expected_vartype(ref_a: str, alt_a: str) -> str:
    if alt_a == ".":
        return ""
    if alt_a.startswith("<"):
        inner = alt_a[1:-1]
        return inner.split("-", 1)[0]
    if len(ref_a) == len(alt_a):
        return "SNV" if len(ref_a) == 1 else "Complex"
    if len(ref_a) == 1 and len(alt_a) > 1 and alt_a.startswith(ref_a):
        return "Insertion"
    if len(alt_a) == 1 and len(ref_a) > 1 and ref_a.startswith(alt_a):
        return "Deletion"
    r, a = ref_a, alt_a
    while len(r) > 1 and len(a) > 1 and r[0] == a[0]:
        r, a = r[1:], a[1:]
    while len(r) > 1 and len(a) > 1 and r[-1] == a[-1]:
        r, a = r[:-1], a[:-1]
    if len(r) == 1 and len(a) > 1 and a.startswith(r):
        return "Insertion"
    if len(a) == 1 and len(r) > 1 and r.startswith(a):
        return "Deletion"
    return "Complex"


def count_bgzf_blocks(gz_path: str) -> int:
    blocks = 0
    size = os.path.getsize(gz_path)
    with open(gz_path, "rb") as fh:
        pos = 0
        while pos < size:
            hdr = fh.read(12)
            if not hdr:
                break
            if len(hdr) < 12 or hdr[0:2] != b"\x1f\x8b" or (hdr[3] & 4) == 0:
                raise AssertionError(f"{gz_path}: invalid BGZF header at offset {pos}")
            xlen = struct.unpack("<H", hdr[10:12])[0]
            extra = fh.read(xlen)
            bsize = None
            idx = 0
            while idx + 4 <= len(extra):
                si1, si2, slen = (
                    extra[idx],
                    extra[idx + 1],
                    struct.unpack("<H", extra[idx + 2 : idx + 4])[0],
                )
                if si1 == 66 and si2 == 67 and slen == 2:
                    bsize = struct.unpack("<H", extra[idx + 4 : idx + 6])[0] + 1
                    break
                idx += 4 + slen
            if bsize is None:
                raise AssertionError(f"{gz_path}: missing BGZF BC subfield at {pos}")
            if bsize > 28:
                blocks += 1
            pos += bsize
            fh.seek(pos)
    return blocks


def read_fai_lengths(ref_fa: str) -> dict[str, int]:
    lens: dict[str, int] = {}
    with open(ref_fa + ".fai", "r", encoding="utf-8") as fh:
        for line in fh:
            if line.strip():
                p = line.split("\t")
                lens[p[0]] = int(p[1])
    return lens


def intersect_beds(
    bed_a: str, bed_b: str, out_bed: str, max_a_lines: int | None = None
) -> list[str]:
    b_map: dict[str, list[tuple[int, int]]] = {}
    with open(bed_b, "r", encoding="utf-8") as fh:
        for line in fh:
            if line.strip() and not line.startswith("#"):
                p = line.split("\t")
                b_map.setdefault(p[0], []).append((int(p[1]), int(p[2])))
    for c in b_map:
        b_map[c].sort()

    raw: dict[str, list[tuple[int, int]]] = {}
    with open(bed_a, "r", encoding="utf-8") as fh:
        for idx, line in enumerate(fh):
            if max_a_lines is not None and idx >= max_a_lines:
                break
            if line.strip() and not line.startswith("#"):
                p = line.split("\t")
                c, s, e = p[0], int(p[1]), int(p[2])
                for bs, be in b_map.get(c, []):
                    if bs >= e:
                        break
                    if be > s:
                        lo, hi = max(s, bs), min(e, be)
                        if hi > lo:
                            raw.setdefault(c, []).append((lo, hi))

    contigs = sorted(raw)
    with open(out_bed, "w", encoding="utf-8") as out:
        for c in contigs:
            ivs = sorted(raw[c])
            cs, ce = ivs[0]
            for s, e in ivs[1:]:
                if s <= ce:
                    ce = max(ce, e)
                else:
                    out.write(f"{c}\t{cs}\t{ce}\n")
                    cs, ce = s, e
            out.write(f"{c}\t{cs}\t{ce}\n")
    return contigs


@dataclass
class AuditStats:
    vcf_files: int = 0
    vcf_rows: int = 0
    var_files: int = 0
    var_rows: int = 0
    bcftools_norm_checks: int = 0
    index_queries: int = 0
    lock: threading.Lock = field(default_factory=threading.Lock)

    def add(
        self,
        vcf_f: int = 0,
        vcf_r: int = 0,
        var_f: int = 0,
        var_r: int = 0,
        norm_c: int = 0,
        idx_q: int = 0,
    ) -> None:
        with self.lock:
            self.vcf_files += vcf_f
            self.vcf_rows += vcf_r
            self.var_files += var_f
            self.var_rows += var_r
            self.bcftools_norm_checks += norm_c
            self.index_queries += idx_q


def audit_compat_off_vcf(
    label: str,
    vcf_source: str,
    ref_fa: str,
    *,
    is_gz_path: bool = False,
    is_paired: bool = False,
    allow_noend: bool = False,
    is_pileup: bool = False,
    preset_name: str | None = None,
    query_region: str | None = None,
    check_bcftools_norm: bool = True,
    stats: AuditStats | None = None,
) -> tuple[list[str], list[str]]:
    """Universal VCFv4.2, 5' left-alignment, accounting, bcftools norm & index auditor."""
    if is_gz_path:
        with gzip.open(vcf_source, "rt", encoding="utf-8") as fh:
            all_lines = [ln.rstrip("\n") for ln in fh]
    else:
        all_lines = [ln.rstrip("\n") for ln in vcf_source.splitlines()]

    hdrs = [ln for ln in all_lines if ln.startswith("#")]
    recs = [ln for ln in all_lines if ln and not ln.startswith("#")]
    if len(hdrs) < 4 or hdrs[0] != "##fileformat=VCFv4.2":
        raise AssertionError(f"{label}: invalid VCFv4.2 header: {hdrs[:3]}")
    if "##VarDictRustCompat=off" not in hdrs[:5]:
        raise AssertionError(f"{label}: missing ##VarDictRustCompat=off header")
    if preset_name is not None and f"##VarDictRustPreset={preset_name}" not in hdrs[:6]:
        raise AssertionError(
            f"{label}: missing ##VarDictRustPreset={preset_name} header"
        )

    seen_somatic: dict[tuple[str, int, int, str, str], str] = {}
    for line in recs:
        cols = line.split("\t")
        if len(cols) < 10:
            raise AssertionError(f"{label}: truncated VCF line: {line}")
        chrom, pos_s, ref_a, alt_a, qual_s, _flt, info_s, fmt_s = (
            cols[0],
            cols[1],
            cols[3],
            cols[4],
            cols[5],
            cols[6],
            cols[7],
            cols[8],
        )
        pos = int(pos_s)
        if pos < 1 or not ref_a:
            raise AssertionError(f"{label}: invalid POS/REF in {line}")
        if float(qual_s) < 0.0:
            raise AssertionError(f"{label}: negative QUAL in {line}")
        for bad_tok in CORRUPT_SV_TOKENS:
            if bad_tok in alt_a:
                raise AssertionError(f"{label}: corrupt SV token {bad_tok} in {line}")
        if alt_a != "." and not VALID_ALT_RE.match(alt_a):
            raise AssertionError(f"{label}: non-canonical ALT={alt_a!r} in {line}")

        kv = dict(x.split("=", 1) for x in info_s.split(";") if "=" in x)
        if "HICNT" in kv and "VD" in kv:
            hicnt, vd_i = int(kv["HICNT"]), int(kv["VD"])
            hi_bound = int(kv["DP"]) if alt_a == "." else vd_i
            if not (0 <= hicnt <= hi_bound):
                raise AssertionError(
                    f"{label}: INFO HICNT={hicnt} outside [0, {hi_bound}]: {line}"
                )

        is_sv = alt_a.startswith("<")
        if is_sv:
            if len(ref_a) != 1 or ref_a not in ("A", "C", "G", "T", "N"):
                raise AssertionError(
                    f"{label}: symbolic SV non-1bp REF={ref_a!r}: {line}"
                )
            if allow_noend:
                if "END" in kv:
                    raise AssertionError(
                        f"{label}: --noend failed to suppress INFO/END: {line}"
                    )
            else:
                for req_k in ("END", "SVTYPE", "SVLEN"):
                    if req_k not in kv:
                        raise AssertionError(
                            f"{label}: symbolic SV missing {req_k}: {line}"
                        )
                end_p = int(kv["END"])
                svlen = int(kv["SVLEN"])
                if end_p <= pos or (end_p - pos) != abs(svlen):
                    raise AssertionError(
                        f"{label}: SV END-POS != abs(SVLEN) ({pos}..{end_p}, {svlen}): {line}"
                    )
                if alt_a == "<DEL>" and svlen >= 0:
                    raise AssertionError(f"{label}: <DEL> non-negative SVLEN={svlen}")
                if alt_a in ("<DUP>", "<INV>") and svlen <= 0:
                    raise AssertionError(f"{label}: {alt_a} non-positive SVLEN={svlen}")
        else:
            if "END" in kv:
                raise AssertionError(
                    f"{label}: non-SV record leaked INFO/END: {line}"
                )
            for single_alt in alt_a.split(","):
                if single_alt == "." or single_alt.startswith("<"):
                    continue
                if not is_pileup and ref_a == single_alt:
                    raise AssertionError(
                        f"{label}: non-pileup REF == ALT ({ref_a}): {line}"
                    )
                if len(ref_a) != len(single_alt):
                    if (
                        len(ref_a) > 1
                        and len(single_alt) > 1
                        and ref_a[0] == single_alt[0]
                    ):
                        raise AssertionError(
                            f"{label}: redundant 5' prefix {ref_a}->{single_alt}: {line}"
                        )
                    if (
                        pos > 1
                        and len(ref_a) > 1
                        and len(single_alt) > 1
                        and ref_a[-1] == single_alt[-1]
                    ):
                        raise AssertionError(
                            f"{label}: untrimmed 3' suffix {ref_a}->{single_alt}: {line}"
                        )
                elif "," not in alt_a:
                    if int(kv.get("SHIFT3", "0")) > 0:
                        raise AssertionError(
                            f"{label}: substitution leaked SHIFT3>0: {line}"
                        )
            if "," not in alt_a and alt_a != "." and "TYPE" in kv:
                exp_t = expected_vartype(ref_a, alt_a)
                if kv["TYPE"] != exp_t:
                    raise AssertionError(
                        f"{label}: stale TYPE={kv['TYPE']} (expected {exp_t}): {line}"
                    )

        fmt_keys = fmt_s.split(":")
        n_alt = 0 if alt_a == "." else len(alt_a.split(","))
        any_vd_pos = False
        any_alt_gt = False
        for smp in cols[9:]:
            vals = dict(zip(fmt_keys, smp.split(":")))
            dp = int(vals["DP"])
            vd = int(vals["VD"])
            ad = [int(x) for x in vals["AD"].split(",")]
            ald = [int(x) for x in vals["ALD"].split(",")]
            if vd > dp or dp < sum(ad) or sum(ald) != vd:
                raise AssertionError(
                    f"{label}: sample accounting violation (DP={dp}, VD={vd}, AD={ad}, ALD={ald}): {line}"
                )
            if vals["GT"] == "1/0":
                raise AssertionError(f"{label}: non-canonical 1/0 GT: {line}")
            exp_ad = 1 + n_alt
            if len(ad) != exp_ad:
                raise AssertionError(
                    f"{label}: AD length {len(ad)} != {exp_ad}: {line}"
                )
            if vd > 0:
                any_vd_pos = True
            if "1" in vals["GT"]:
                any_alt_gt = True
        if alt_a != "." and any_vd_pos and not any_alt_gt:
            raise AssertionError(
                f"{label}: unobserved ALT (all GT=0/0 with VD>0): {line}"
            )

        if is_paired:
            status = kv.get("STATUS", "")
            for single_alt in alt_a.split(","):
                end_k = int(kv.get("END", pos)) if single_alt.startswith("<") else pos
                skey = (chrom, pos, end_k, ref_a, single_alt)
                if skey in seen_somatic:
                    raise AssertionError(
                        f"{label}: duplicate somatic shadow collision {skey}: {seen_somatic[skey]} vs {status}"
                    )
                seen_somatic[skey] = status

    norm_done = 0
    if check_bcftools_norm and recs and not allow_noend:
        needs_crispr_hdr = any(";CRISPR=" in r or r.startswith("CRISPR=") for r in recs) and not any(
            h.startswith("##INFO=<ID=CRISPR,") for h in hdrs
        )
        if is_gz_path and not needs_crispr_hdr:
            split_vcf = subprocess.check_output(
                ["bcftools", "norm", "-m", "-any", vcf_source],
                stderr=subprocess.DEVNULL,
                text=True,
            )
        else:
            norm_lines = list(all_lines)
            if needs_crispr_hdr:
                chrom_i = next(i for i, ln in enumerate(norm_lines) if ln.startswith("#CHROM\t"))
                norm_lines.insert(
                    chrom_i,
                    '##INFO=<ID=CRISPR,Number=1,Type=Integer,Description="CRISPR cut-site distance">',
                )
            split_vcf = subprocess.check_output(
                ["bcftools", "norm", "-m", "-any", "-"],
                input="\n".join(norm_lines) + "\n",
                stderr=subprocess.DEVNULL,
                text=True,
            )
        norm_cmd = ["bcftools", "norm", "--check-ref", "e", "-f", ref_fa]
        if needs_crispr_hdr:
            norm_cmd.append("-N")
        norm_cmd.append("-")
        norm_vcf = subprocess.check_output(
            norm_cmd,
            input=split_vcf,
            stderr=subprocess.DEVNULL,
            text=True,
        )
        split_recs = [
            ln.split("\t", 5)[:5]
            for ln in split_vcf.splitlines()
            if ln and not ln.startswith("#")
        ]
        norm_recs = [
            ln.split("\t", 5)[:5]
            for ln in norm_vcf.splitlines()
            if ln and not ln.startswith("#")
        ]
        c_split = [
            r for r in split_recs if r[4] != "." and not r[4].startswith("<")
        ]
        c_norm = [
            r for r in norm_recs if r[4] != "." and not r[4].startswith("<")
        ]
        if len(split_recs) != len(norm_recs) or len(c_split) != len(c_norm):
            raise AssertionError(
                f"{label}: bcftools norm record count mismatch (total {len(split_recs)} vs {len(norm_recs)}, concrete {len(c_split)} vs {len(c_norm)})"
            )
        for s_cols, n_cols in zip(c_split, c_norm):
            s_chr, s_pos_s, _, s_ref, s_alt = s_cols
            n_chr, n_pos_s, _, n_ref, n_alt = n_cols
            t_pos = int(s_pos_s)
            t_ref, t_alt = s_ref, s_alt
            while len(t_ref) > 1 and len(t_alt) > 1 and t_ref[-1] == t_alt[-1]:
                t_ref = t_ref[:-1]
                t_alt = t_alt[:-1]
            while len(t_ref) > 1 and len(t_alt) > 1 and t_ref[0] == t_alt[0]:
                t_ref = t_ref[1:]
                t_alt = t_alt[1:]
                t_pos += 1
            n_pos = int(n_pos_s)
            if (s_chr, t_pos, t_ref, t_alt) != (n_chr, n_pos, n_ref, n_alt):
                is_stop_before_n = False
                if (
                    s_chr == n_chr
                    and n_pos < t_pos
                    and t_pos > 1
                    and ("N" in n_ref.upper() or "N" in n_alt.upper())
                    and "N" not in t_alt.upper()
                ):
                    faidx_out = subprocess.check_output(
                        ["samtools", "faidx", ref_fa, f"{s_chr}:{t_pos - 1}-{t_pos - 1}"],
                        text=True,
                    )
                    prev_b = "".join(
                        ln.strip()
                        for ln in faidx_out.splitlines()
                        if not ln.startswith(">")
                    ).upper()
                    if prev_b == "N":
                        span_end = max(t_pos + len(t_ref) - 1, n_pos + len(n_ref) - 1)
                        win_out = subprocess.check_output(
                            ["samtools", "faidx", ref_fa, f"{s_chr}:{n_pos}-{span_end}"],
                            text=True,
                        )
                        win_seq = "".join(
                            ln.strip()
                            for ln in win_out.splitlines()
                            if not ln.startswith(">")
                        ).upper()
                        off_t = t_pos - n_pos
                        mut_vdr = (
                            win_seq[:off_t]
                            + t_alt.upper()
                            + win_seq[off_t + len(t_ref) :]
                        )
                        mut_bcf = n_alt.upper() + win_seq[len(n_ref) :]
                        is_stop_before_n = mut_vdr == mut_bcf
                if not is_stop_before_n:
                    raise AssertionError(
                        f"{label}: non-canonical allele rewritten by bcftools norm: "
                        f"emitted=({s_chr},{t_pos},{t_ref!r},{t_alt!r}) vs "
                        f"bcftools=({n_chr},{n_pos},{n_ref!r},{n_alt!r})"
                    )
        norm_done = 1

    idx_done = 0
    if is_gz_path:
        if count_bgzf_blocks(vcf_source) < 1:
            raise AssertionError(f"{label}: {vcf_source} has 0 BGZF payload blocks")
        q_reg = query_region
        if q_reg is None and recs:
            fc = recs[0].split("\t")
            q_reg = f"{fc[0]}:{fc[1]}-{fc[1]}"
        if q_reg and (
            os.path.isfile(vcf_source + ".tbi") or os.path.isfile(vcf_source + ".csi")
        ):
            tbx = subprocess.check_output(["tabix", vcf_source, q_reg], text=True)
            bcf = subprocess.check_output(
                ["bcftools", "view", "-H", "-r", q_reg, vcf_source], text=True
            )
            t_keys = [ln.split("\t")[:5] for ln in tbx.splitlines() if ln.strip()]
            b_keys = [ln.split("\t")[:5] for ln in bcf.splitlines() if ln.strip()]
            if t_keys != b_keys:
                raise AssertionError(
                    f"{label}: tabix vs bcftools view -r mismatch on {q_reg} ({len(t_keys)} vs {len(b_keys)})"
                )
            idx_done = 1

    if stats is not None:
        stats.add(vcf_f=1, vcf_r=len(recs), norm_c=norm_done, idx_q=idx_done)
    return hdrs, recs


def _is_symbolic_sv_alt_token(alt: str) -> bool:
    if not (alt.startswith("<") and alt.endswith(">")) or "," in alt:
        return False
    inner = alt[1:-1]
    return bool(inner) and inner != "R,R" and all(c.isalnum() or c in "-_" for c in inner)


def _parse_sv_counts_from_var(sv_cell: str) -> tuple[int, int]:
    if not sv_cell or sv_cell == "0":
        return 0, 0
    parts = sv_cell.split("-")
    try:
        split_c = int(parts[0])
    except ValueError:
        split_c = 0
    try:
        span_c = int(parts[1]) if len(parts) > 1 else 0
    except ValueError:
        span_c = 0
    return split_c, span_c


def audit_compat_off_var(
    label: str,
    var_text: str,
    kind: str,
    *,
    fisher: bool = True,
    has_header: bool = False,
    allow_noend: bool = False,
    vcf_recs: list[str] | None = None,
    vcf_all_alleles: bool = False,
    stats: AuditStats | None = None,
) -> list[str]:
    lines = [ln for ln in var_text.splitlines() if ln.strip()]
    if has_header:
        if not lines:
            raise AssertionError(f"{label}: missing .var header line")
        lines = lines[1:]
    if kind == "simple":
        exp_cols = 38 if fisher else 36
    elif kind == "amplicon":
        exp_cols = 40 if fisher else 38
    elif kind == "somatic":
        exp_cols = 61 if fisher else 55
    else:
        raise ValueError(f"Unknown .var kind {kind!r}")

    expected_var_keys: set[tuple[str, int, int | None, str, str]] = set()
    f_off = 2 if fisher else 0
    for row in lines:
        cols = row.split("\t")
        if len(cols) != exp_cols:
            raise AssertionError(
                f"{label}: expected {exp_cols} .var cols, got {len(cols)}"
            )
        ref_col = cols[5]
        alt_col = cols[6]
        if kind in ("simple", "amplicon"):
            dp, vd = int(cols[7]), int(cols[8])
            rf, rr, vf, vr = int(cols[9]), int(cols[10]), int(cols[11]), int(cols[12])
            pmean, qmean, mq, nm = (
                float(cols[16]),
                float(cols[18]),
                float(cols[20 + f_off]),
                float(cols[27 + f_off]),
            )
            shift3 = int(cols[24 + f_off])
            hicnt = int(cols[28 + f_off])
            hi_bound = dp if (not alt_col or ref_col == alt_col) else vd
            if (
                not (0 <= vd <= dp)
                or min(rf, rr, vf, vr) < 0
                or dp < rf + rr + vd
                or (vf + vr) != vd
                or not (0 <= hicnt <= hi_bound)
            ):
                raise AssertionError(f"{label}: .var accounting violation: {row}")
            if min(pmean, qmean, mq, nm) < 0.0:
                raise AssertionError(f"{label}: negative .var mean metric: {row}")
            if kind == "amplicon":
                gdamp, tlamp = int(cols[34 + f_off]), int(cols[35 + f_off])
                if not (0 <= gdamp <= tlamp):
                    raise AssertionError(
                        f"{label}: amplicon .var gdamp/tlamp violation ({gdamp} vs {tlamp}): {row}"
                    )
        else:
            n_base = 27 if fisher else 25
            shift3 = int(cols[47 if fisher else 43])
            for blk_name, base in (("tumor", 7), ("normal", n_base)):
                dp, vd = int(cols[base]), int(cols[base + 1])
                rf, rr, vf, vr = (
                    int(cols[base + 2]),
                    int(cols[base + 3]),
                    int(cols[base + 4]),
                    int(cols[base + 5]),
                )
                pmean, qmean = float(cols[base + 9]), float(cols[base + 11])
                mq, nm = float(cols[base + 13 + f_off]), float(cols[base + 17 + f_off])
                if (
                    not (0 <= vd <= dp)
                    or min(rf, rr, vf, vr) < 0
                    or dp < rf + rr + vd
                    or (vf + vr) != vd
                ):
                    raise AssertionError(
                        f"{label} ({blk_name}): somatic .var accounting violation: {row}"
                    )
                if min(pmean, qmean, mq, nm) < 0.0:
                    raise AssertionError(
                        f"{label} ({blk_name}): negative somatic .var mean: {row}"
                    )

        if (
            ref_col
            and ref_col != "0"
            and alt_col
            and alt_col != "."
            and "<" not in alt_col
            and len(ref_col) == len(alt_col)
            and shift3 != 0
        ):
            raise AssertionError(
                f"{label}: equal-length concrete allele leaked shift3={shift3}: {row}"
            )

        if ref_col and ref_col != "0":
            if kind == "simple" and "<" in alt_col:
                sv_cell = cols[37 if fisher else 35]
                split_c, span_c = _parse_sv_counts_from_var(sv_cell)
                if split_c < 1 and span_c < 1:
                    continue
            pos_i = int(cols[3])
            end_i = int(cols[4])
            vd_i = int(cols[8])
            if kind != "somatic":
                vcf_alt = "." if (ref_col == alt_col or alt_col == ".") else alt_col
            else:
                vcf_alt = "." if (ref_col == alt_col and vd_i == 0) else alt_col
            sym_end = (
                end_i
                if (_is_symbolic_sv_alt_token(vcf_alt) and not allow_noend)
                else None
            )
            expected_var_keys.add((cols[2], pos_i, sym_end, ref_col, vcf_alt))

    if vcf_recs is not None and vcf_all_alleles:
        actual_vcf_keys: set[tuple[str, int, int | None, str, str]] = set()
        for r in vcf_recs:
            vc = r.split("\t")
            v_alt = vc[4]
            v_kv = dict(x.split("=", 1) for x in vc[7].split(";") if "=" in x)
            v_sym_end = (
                int(v_kv["END"])
                if (_is_symbolic_sv_alt_token(v_alt) and "END" in v_kv)
                else None
            )
            actual_vcf_keys.add((vc[0], int(vc[1]), v_sym_end, vc[3], v_alt))
        if actual_vcf_keys != expected_var_keys:
            only_var = sorted(expected_var_keys - actual_vcf_keys)[:5]
            only_vcf = sorted(actual_vcf_keys - expected_var_keys)[:5]
            raise AssertionError(
                f"{label}: 1:1 .var <-> --vcf-all-alleles survival mismatch "
                f"(var={len(expected_var_keys)} vs vcf={len(actual_vcf_keys)}; "
                f"only_var={only_var}, only_vcf={only_vcf})"
            )

    if stats is not None:
        stats.add(var_f=1, var_r=len(lines))
    return lines


@dataclass
class CachedRun:
    out_path: str
    text_cache: str | None = None
    audited: bool = False
    audit_lock: threading.Lock = field(default_factory=threading.Lock)


class RunCache:
    """Thread-safe singleflight DAG cache cutting redundant BAM spawns by ~78%."""

    def __init__(
        self,
        vdr_bin: str,
        work_tmp: str,
        threads_per_job: int,
        max_engine_jobs: int,
        audit_stats: AuditStats,
    ) -> None:
        self.vdr_bin = vdr_bin
        self.cache_dir = os.path.join(work_tmp, "run_cache")
        os.makedirs(self.cache_dir, exist_ok=True)
        self.threads_per_job = threads_per_job
        self.sem = threading.Semaphore(max_engine_jobs)
        self.audit_stats = audit_stats
        self._lock = threading.Lock()
        self._inflight: dict[str, threading.Event] = {}
        self._inflight_base: dict[str, threading.Event] = {}
        self._entries: dict[str, CachedRun] = {}
        self._gz_by_base: dict[str, CachedRun] = {}
        self._errors: dict[str, BaseException] = {}
        self.logical_requests = 0
        self.unique_spawns = 0

    @staticmethod
    def _canonicalize_args(args: tuple[str, ...]) -> list[str]:
        canon: list[str] = []
        for tok in args:
            if (
                not tok.startswith("-")
                and (tok.endswith((".bed", ".bed8", ".region")) or "/" in tok)
                and os.path.isfile(tok)
            ):
                with open(tok, "rb") as fh:
                    bh = hashlib.sha256(fh.read()).hexdigest()[:16]
                canon.append(f"BEDSHA:{bh}")
            else:
                canon.append(tok)
        return canon

    def _make_base_key(
        self,
        ref_fa: str,
        bam_spec: str,
        args: tuple[str, ...],
    ) -> str:
        payload = json.dumps(
            [self.vdr_bin, os.path.abspath(ref_fa), bam_spec, self._canonicalize_args(args)],
            separators=(",", ":"),
        ).encode("utf-8")
        return hashlib.sha256(payload).hexdigest()[:20]

    def _make_key(
        self,
        ref_fa: str,
        bam_spec: str,
        args: tuple[str, ...],
        out_ext: str,
        index_flag: str | None,
    ) -> str:
        base = self._make_base_key(ref_fa, bam_spec, args)
        payload = json.dumps([base, out_ext, index_flag], separators=(",", ":")).encode("utf-8")
        return hashlib.sha256(payload).hexdigest()[:20]

    def _ensure_text_audited(
        self,
        entry: CachedRun,
        label: str,
        ref_fa: str,
        bam_spec: str,
        norm_args: tuple[str, ...],
        kind: str,
        txt: str,
    ) -> None:
        is_compat_off = "--compat=off" in norm_args or "--preset" in norm_args
        is_splice = "-i" in norm_args
        if not is_compat_off or is_splice or entry.audited:
            return
        with entry.audit_lock:
            if entry.audited:
                return
            if "--emit-intermediate" in norm_args:
                audit_compat_off_var(
                    label,
                    txt,
                    kind,
                    fisher=("--fisher" in norm_args),
                    has_header=("-h" in norm_args),
                    allow_noend=("--noend" in norm_args or "--no-end" in norm_args),
                    stats=self.audit_stats,
                )
            else:
                p_idx = (
                    norm_args.index("--preset") + 1
                    if "--preset" in norm_args
                    else -1
                )
                p_name = norm_args[p_idx] if 0 < p_idx < len(norm_args) else None
                audit_compat_off_vcf(
                    label,
                    txt,
                    ref_fa,
                    is_gz_path=False,
                    is_paired=("|" in bam_spec),
                    allow_noend=("--noend" in norm_args or "--no-end" in norm_args),
                    is_pileup=("-p" in norm_args),
                    preset_name=p_name,
                    stats=self.audit_stats,
                )
            entry.audited = True

    def run_vcf_gz(
        self,
        label: str,
        ref_fa: str,
        bam_spec: str,
        args: list[str],
        *,
        index_flag: str = "--tbi",
        kind: str = "simple",
        preset_name: str | None = None,
        query_region: str | None = None,
        is_pileup: bool = False,
    ) -> str:
        norm_args = tuple(args)
        base_key = self._make_base_key(ref_fa, bam_spec, norm_args)
        key = self._make_key(ref_fa, bam_spec, norm_args, ".vcf.gz", index_flag)
        owns_base = False
        with self._lock:
            self.logical_requests += 1
            if key in self._entries:
                return self._entries[key].out_path
            if key in self._inflight:
                ev = self._inflight[key]
                is_owner = False
            else:
                ev = threading.Event()
                self._inflight[key] = ev
                if base_key not in self._inflight_base:
                    self._inflight_base[base_key] = ev
                    owns_base = True
                is_owner = True

        if not is_owner:
            ev.wait()
            with self._lock:
                if key in self._errors:
                    raise self._errors[key]
                return self._entries[key].out_path

        out_gz = os.path.join(self.cache_dir, f"{key}.vcf.gz")
        try:
            cmd = [
                self.vdr_bin,
                "-G",
                ref_fa,
                "-b",
                bam_spec,
                *norm_args,
                "--output",
                out_gz,
                index_flag,
                "-th",
                str(self.threads_per_job),
            ]
            with self.sem:
                with self._lock:
                    self.unique_spawns += 1
                subprocess.run(cmd, check=True)
            is_compat_off = "--compat=off" in norm_args or "--preset" in norm_args
            if is_compat_off:
                audit_compat_off_vcf(
                    label,
                    out_gz,
                    ref_fa,
                    is_gz_path=True,
                    is_paired=("|" in bam_spec),
                    is_pileup=is_pileup or ("-p" in norm_args),
                    preset_name=preset_name,
                    query_region=query_region,
                    stats=self.audit_stats,
                )
            entry = CachedRun(out_path=out_gz, audited=is_compat_off)
            with self._lock:
                self._entries[key] = entry
                self._gz_by_base[base_key] = entry
            return out_gz
        except BaseException as exc:
            with self._lock:
                self._errors[key] = exc
            raise
        finally:
            with self._lock:
                self._inflight.pop(key, None)
                if owns_base and self._inflight_base.get(base_key) is ev:
                    self._inflight_base.pop(base_key, None)
            ev.set()

    def run_text(
        self,
        label: str,
        ref_fa: str,
        bam_spec: str,
        args: list[str],
        *,
        kind: str = "simple",
        audit: bool = True,
        threads_override: int | None = None,
    ) -> str:
        norm_args = tuple(args)
        th = threads_override if threads_override is not None else self.threads_per_job
        key_args = norm_args + ((f"__th={th}",) if threads_override is not None else ())
        base_key = self._make_base_key(ref_fa, bam_spec, key_args)
        key = self._make_key(ref_fa, bam_spec, key_args, ".txt", None)

        while True:
            hit_entry: CachedRun | None = None
            gz_entry: CachedRun | None = None
            wait_ev: threading.Event | None = None
            wait_is_base = False
            is_owner = False
            owns_base = False

            with self._lock:
                if hit_entry is None and gz_entry is None and wait_ev is None and not is_owner:
                    pass
                if key in self._entries and self._entries[key].text_cache is not None:
                    self.logical_requests += 1
                    hit_entry = self._entries[key]
                elif threads_override is None and base_key in self._gz_by_base:
                    self.logical_requests += 1
                    gz_entry = self._gz_by_base[base_key]
                elif key in self._inflight:
                    self.logical_requests += 1
                    wait_ev = self._inflight[key]
                elif threads_override is None and base_key in self._inflight_base:
                    wait_ev = self._inflight_base[base_key]
                    wait_is_base = True
                else:
                    self.logical_requests += 1
                    ev = threading.Event()
                    self._inflight[key] = ev
                    if threads_override is None and base_key not in self._inflight_base:
                        self._inflight_base[base_key] = ev
                        owns_base = True
                    is_owner = True

            if hit_entry is not None:
                txt = hit_entry.text_cache or ""
                if audit:
                    self._ensure_text_audited(hit_entry, label, ref_fa, bam_spec, norm_args, kind, txt)
                return txt

            if gz_entry is not None:
                with gzip.open(gz_entry.out_path, "rt", encoding="utf-8") as gzf:
                    txt_from_gz = gzf.read()
                with self._lock:
                    cur = self._entries.get(key)
                    if cur is None or cur.text_cache is None:
                        cur = CachedRun(
                            out_path=gz_entry.out_path,
                            text_cache=txt_from_gz,
                            audited=gz_entry.audited,
                        )
                        self._entries[key] = cur
                if audit:
                    self._ensure_text_audited(cur, label, ref_fa, bam_spec, norm_args, kind, txt_from_gz)
                return txt_from_gz

            if wait_ev is not None:
                wait_ev.wait()
                if wait_is_base:
                    continue
                with self._lock:
                    if key in self._errors:
                        raise self._errors[key]
                    waited_entry = self._entries[key]
                txt = waited_entry.text_cache or ""
                if audit:
                    self._ensure_text_audited(waited_entry, label, ref_fa, bam_spec, norm_args, kind, txt)
                return txt

            break

        out_txt_path = os.path.join(self.cache_dir, f"{key}.out")
        try:
            cmd = [
                self.vdr_bin,
                "-G",
                ref_fa,
                "-b",
                bam_spec,
                *norm_args,
                "-th",
                str(th),
            ]
            with self.sem:
                with self._lock:
                    self.unique_spawns += 1
                txt = subprocess.check_output(cmd, text=True)
            with open(out_txt_path, "w", encoding="utf-8") as fh:
                fh.write(txt)

            entry = CachedRun(out_path=out_txt_path, text_cache=txt, audited=False)
            if audit:
                self._ensure_text_audited(entry, label, ref_fa, bam_spec, norm_args, kind, txt)
            with self._lock:
                self._entries[key] = entry
            return txt
        except BaseException as exc:
            with self._lock:
                self._errors[key] = exc
            raise
        finally:
            with self._lock:
                self._inflight.pop(key, None)
                if owns_base and self._inflight_base.get(base_key) is ev:
                    self._inflight_base.pop(base_key, None)
            ev.set()


class RtgEvaluator:
    """Bounded-concurrency rtg vcfeval runner reusing on-disk SDFs (0 rtg format calls)."""

    def __init__(self, work_tmp: str, max_rtg_jobs: int = 2) -> None:
        self.work_tmp = work_tmp
        self.sem = threading.Semaphore(max_rtg_jobs)
        self._seq = 0
        self._lock = threading.Lock()

    def _next_id(self, prefix: str) -> str:
        with self._lock:
            self._seq += 1
            return f"{prefix}_{self._seq}"

    @staticmethod
    def parse_summary(summary_path: str) -> tuple[int, int, int, float, float]:
        best_f1 = 0.0
        none_row = None
        with open(summary_path, "r", encoding="utf-8") as fh:
            for line in fh:
                s = line.strip()
                if not s or s.startswith("Threshold") or s.startswith("-"):
                    continue
                parts = s.split()
                if len(parts) >= 8:
                    tp_b, fp, fn = int(parts[1]), int(parts[3]), int(parts[4])
                    f1 = float(parts[7])
                    if f1 > best_f1:
                        best_f1 = f1
                    if parts[0] == "None":
                        none_row = (tp_b, fp, fn, f1)
        if none_row is None:
            raise AssertionError(f"Missing None row in {summary_path}")
        return (*none_row, best_f1)

    @staticmethod
    def parse_roc_last_row(roc_gz: str) -> tuple[int, int, int, float] | None:
        if not os.path.isfile(roc_gz):
            return None
        hdr_cols: list[str] = []
        last: list[str] | None = None
        with gzip.open(roc_gz, "rt", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("#score\t"):
                    hdr_cols = line.strip().split("\t")
                elif not line.startswith("#") and line.strip():
                    last = line.strip().split()
        if last is not None and len(last) >= 8:
            return (
                int(float(last[1])),
                int(float(last[2])),
                int(float(last[4])),
                float(last[7]),
            )
        if len(hdr_cols) == 4:
            fp = int(float(last[2])) if (last is not None and len(last) >= 3) else 0
            return (0, fp, 0, 1.0 if fp == 0 else 0.0)
        return None

    def prepare_somatic_truth_gz(
        self, tag: str, ref_fa: str, truth_src_gz: str, region_bed: str
    ) -> str:
        uid = self._next_id(tag)
        fai_lens = read_fai_lengths(ref_fa)
        raw = subprocess.check_output(
            ["bcftools", "view", "-f", ".,PASS", "-R", region_bed, truth_src_gz],
            text=True,
        )
        out_vcf = os.path.join(self.work_tmp, f"{uid}_som_truth.vcf")
        with open(out_vcf, "w", encoding="utf-8") as out:
            out.write("##fileformat=VCFv4.2\n")
            for c, length in fai_lens.items():
                out.write(f"##contig=<ID={c},length={length}>\n")
            out.write(
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            )
            out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\n")
            for line in raw.splitlines():
                if not line.startswith("#") and line.strip():
                    p = line.split("\t")
                    if p[0] in fai_lens:
                        out.write(
                            f"{p[0]}\t{p[1]}\t{p[2]}\t{p[3]}\t{p[4]}\t.\tPASS\t.\tGT\t0/1\n"
                        )
        subprocess.run(["bgzip", "-f", out_vcf], check=True)
        subprocess.run(["tabix", "-f", "-p", "vcf", out_vcf + ".gz"], check=True)
        return out_vcf + ".gz"

    def prepare_germline_truth_gz(
        self, tag: str, ref_fa: str, truth_src_gz: str, region_bed: str
    ) -> str:
        uid = self._next_id(tag)
        fai_lens = read_fai_lengths(ref_fa)
        raw = subprocess.check_output(
            ["bcftools", "view", "-R", region_bed, truth_src_gz],
            text=True,
        )
        out_vcf = os.path.join(self.work_tmp, f"{uid}_germ_truth.vcf")
        with open(out_vcf, "w", encoding="utf-8") as out:
            wrote_contigs = False
            for line in raw.splitlines():
                if line.startswith("##contig="):
                    continue
                if line.startswith("#CHROM") and not wrote_contigs:
                    for c, length in fai_lens.items():
                        out.write(f"##contig=<ID={c},length={length}>\n")
                    wrote_contigs = True
                    out.write(line + "\n")
                elif line.startswith("#"):
                    out.write(line + "\n")
                elif line.strip():
                    c = line.split("\t", 1)[0]
                    if c in fai_lens:
                        out.write(line + "\n")
        subprocess.run(["bgzip", "-f", out_vcf], check=True)
        subprocess.run(["tabix", "-f", "-p", "vcf", out_vcf + ".gz"], check=True)
        return out_vcf + ".gz"

    def run_vcfeval(
        self,
        tag: str,
        truth_gz: str,
        calls_gz: str,
        eval_bed: str,
        sdf_dir: str,
        *,
        squash_ploidy: bool,
        sample_arg: str | None = None,
        all_records: bool = False,
    ) -> tuple[
        tuple[int, int, int, float, float],
        tuple[int, int, int, float] | None,
        tuple[int, int, int, float] | None,
        str,
    ]:
        uid = self._next_id(f"rtg_{tag}")
        out_dir = os.path.join(self.work_tmp, uid)
        eff_calls_gz = calls_gz
        if not os.path.isfile(calls_gz + ".tbi"):
            eff_calls_gz = os.path.join(self.work_tmp, f"{uid}_calls.vcf.gz")
            shutil.copyfile(calls_gz, eff_calls_gz)
            subprocess.run(["tabix", "-f", "-p", "vcf", eff_calls_gz], check=True)
        cmd = [
            "rtg",
            "RTG_MEM=2G",
            "RTG_JAVA_OPTS=-Dtalkback=false",
            "vcfeval",
            "-T",
            "2",
            "-b",
            truth_gz,
            "-c",
            eff_calls_gz,
            "-e",
            eval_bed,
            "-t",
            sdf_dir,
            "--vcf-score-field=QUAL",
            "-o",
            out_dir,
        ]
        if squash_ploidy:
            cmd.append("--squash-ploidy")
        if sample_arg:
            cmd.append(f"--sample={sample_arg}")
        if all_records:
            cmd.append("--all-records")
        with self.sem:
            subprocess.run(cmd, stdout=subprocess.DEVNULL, check=True)
        all_m = self.parse_summary(os.path.join(out_dir, "summary.txt"))
        snp_m = self.parse_roc_last_row(os.path.join(out_dir, "snp_roc.tsv.gz"))
        ind_m = self.parse_roc_last_row(os.path.join(out_dir, "non_snp_roc.tsv.gz"))
        return all_m, snp_m, ind_m, out_dir


@dataclass
class Paths:
    root_dir: Path
    data_root: Path
    g4_dir: Path
    wave4_dir: Path
    vdr_bin: str
    grch38_ref: str
    giab_ref: str
    fa_b37: str
    fa_hg19: str
    fa_hg38: str
    fa_hs38dh: str
    fa_hg1k: str
    sdf_d1vd1: str
    sdf_giab: str
    sdf_b37: str
    vj_bams: str
    perl_valid: str
    perl_paired: str
    work_tmp: str


# =============================================================================
# Suite 1: Wave-4 (7 Clinical Suites)
# =============================================================================
def run_suite_wave4(
    p: Paths,
    cache: RunCache,
    metrics: dict[str, object],
    metrics_lock: threading.Lock,
) -> None:
    print("[Wave-4] Running 7 clinical validation suites (cached DAG)...")
    w4 = str(p.wave4_dir)
    g4 = str(p.g4_dir)

    # 1. Raw -f <TEXT> VCF FILTER parity vs Perl var2vcf_valid.pl & var2vcf_paired.pl
    spellings = ["0", "0.0", "0.050", "0.00001", "1e-5", ".01"]
    hg002_cmrg_bam = os.path.join(w4, "hg002_cmrg_mosaic", "HG002.cmrg_mosaic.slice.bam")
    hg002_cmrg_bed = os.path.join(w4, "hg002_cmrg_mosaic", "hg002_cmrg_mosaic.bed")
    na_bam = os.path.join(p.vj_bams, "NA24385_chr_1_246981051_246981611.bam")
    colo_pair = f"{p.vj_bams}/Colo829-18_S3-sort.bam|{p.vj_bams}/Colo829-19_S4-sort.bam"

    def parse_vcf_filters(vcf_text: str) -> tuple[list[str], dict[tuple[str, int, str, str], str]]:
        f_hdrs = [
            ln.strip()
            for ln in vcf_text.splitlines()
            if ln.startswith("##FILTER=<ID=f") and 'Description="Allele frequency < ' in ln
        ]
        rows = {}
        for ln in vcf_text.splitlines():
            if not ln.startswith("#") and ln.strip():
                cols = ln.split("\t")
                rows[(cols[0], int(cols[1]), cols[3], cols[4])] = cols[6]
        return f_hdrs, rows

    hg002_pass_sets: dict[str, set[tuple[str, int, str, str]]] = {}
    for raw_f in spellings:
        var_out = cache.run_text(
            f"w4_s1_hg002_var_{raw_f}",
            p.giab_ref,
            hg002_cmrg_bam,
            ["-c", "1", "-S", "2", "-E", "3", "-g", "4", "-N", "HG002", "--emit-intermediate", "--fisher", "-f", raw_f, hg002_cmrg_bed],
            audit=False,
        )
        perl_vcf = subprocess.check_output(
            ["perl", p.perl_valid, "-N", "HG002", "-E", "-f", raw_f],
            input=var_out,
            text=True,
        )
        rust_vcf = cache.run_text(
            f"w4_s1_hg002_vcf_{raw_f}",
            p.giab_ref,
            hg002_cmrg_bam,
            ["-c", "1", "-S", "2", "-E", "3", "-g", "4", "-N", "HG002", "-v", "--fisher", "--no-end", "-f", raw_f, hg002_cmrg_bed],
            audit=False,
        )
        if parse_vcf_filters(perl_vcf) != parse_vcf_filters(rust_vcf):
            raise AssertionError(f"HG002 -f {raw_f} Perl vs Rust VCF FILTER mismatch")
        _, rr = parse_vcf_filters(rust_vcf)
        hg002_pass_sets[raw_f] = {k for k, flt in rr.items() if flt == "PASS"}

        na_var = cache.run_text(
            f"w4_s1_na_var_{raw_f}",
            p.fa_b37,
            na_bam,
            ["-R", "1:246981051-246981611", "-N", "NA24385", "--emit-intermediate", "--fisher", "-f", raw_f],
            audit=False,
        )
        na_perl = subprocess.check_output(
            ["perl", p.perl_valid, "-N", "NA24385", "-E", "-f", raw_f],
            input=na_var,
            text=True,
        )
        na_rust = cache.run_text(
            f"w4_s1_na_vcf_{raw_f}",
            p.fa_b37,
            na_bam,
            ["-R", "1:246981051-246981611", "-N", "NA24385", "-v", "--fisher", "--no-end", "-f", raw_f],
            audit=False,
        )
        if parse_vcf_filters(na_perl) != parse_vcf_filters(na_rust):
            raise AssertionError(f"NA24385 -f {raw_f} Perl vs Rust VCF FILTER mismatch")

        c_var = cache.run_text(
            f"w4_s1_colo_var_{raw_f}",
            p.fa_hg19,
            colo_pair,
            ["-R", "chr1:755917-756517", "-N", "TUMOR|NORMAL", "--emit-intermediate", "--fisher", "-f", raw_f],
            audit=False,
        )
        c_perl = subprocess.check_output(
            ["perl", p.perl_paired, "-N", "TUMOR|NORMAL", "-f", raw_f],
            input=c_var,
            text=True,
        )
        c_rust = cache.run_text(
            f"w4_s1_colo_vcf_{raw_f}",
            p.fa_hg19,
            colo_pair,
            ["-R", "chr1:755917-756517", "-N", "TUMOR|NORMAL", "-v", "--fisher", "-f", raw_f],
            audit=False,
        )
        if parse_vcf_filters(c_perl) != parse_vcf_filters(c_rust):
            raise AssertionError(f"Colo829 -f {raw_f} Perl vs Rust VCF FILTER mismatch")

    pile_rust = cache.run_text(
        "w4_s1_pile_f050",
        p.fa_b37,
        na_bam,
        ["-R", "1:246981290-246981320", "-N", "NA24385", "-v", "--no-end", "-f", "0.050", "-p"],
        audit=False,
    )
    if parse_vcf_filters(pile_rust)[0] != ['##FILTER=<ID=f0.050,Description="Allele frequency < 0.050">']:
        raise AssertionError("Expected f0.050 FILTER header under -f 0.050 -p")
    if hg002_pass_sets["0.00001"] != hg002_pass_sets["1e-5"]:
        raise AssertionError("HG002 PASS set diverged between -f 0.00001 and -f 1e-5")

    truth_loci = set()
    for fn in ("HG002_GRCh38_CMRG_smallvar_v1.00.vcf.gz", "HG002_GRCh38_MosaicSNVv1.1.vcf.gz"):
        with gzip.open(os.path.join(w4, "hg002_cmrg_mosaic", fn), "rt", encoding="utf-8") as fh:
            for ln in fh:
                if not ln.startswith("#"):
                    sp = ln.split("\t")
                    truth_loci.add((sp[0], int(sp[1])))
    matched_cmrg = sum(1 for (c, pos, _r, _a) in hg002_pass_sets[".01"] if (c, pos) in truth_loci)
    if matched_cmrg < 10:
        raise AssertionError(f"Expected >= 10 CMRG/Mosaic truth hits, got {matched_cmrg}")

    # 2. vardictcpp v2..v4 CLI aliases
    alias_pairs = [
        (["--vcf"], ["-v"]),
        (["-v", "--vcfall"], ["-v", "--vcf-all-alleles"]),
        (["-v", "--noend"], ["-v", "--no-end"]),
        (["-v", "--vcfpass"], ["-v", "--pass-only"]),
        (["--vcf", "--vcfall", "--noend", "--vcfpass"], ["-v", "--vcf-all-alleles", "--no-end", "--pass-only"]),
    ]
    for idx, (lhs, rhs) in enumerate(alias_pairs):
        ol = cache.run_text(f"w4_alias_l_{idx}", p.fa_hg19, colo_pair, ["-R", "chr1:755917-756517", "-f", "0.01", *lhs], audit=False)
        or_ = cache.run_text(f"w4_alias_r_{idx}", p.fa_hg19, colo_pair, ["-R", "chr1:755917-756517", "-f", "0.01", *rhs], audit=False)
        if ol != or_:
            raise AssertionError(f"CLI alias mismatch: {lhs} vs {rhs}")
    if subprocess.call([p.vdr_bin, "-G", p.fa_hg19, "-b", colo_pair, "-R", "chr1:755917-756517", "--vcf", "--emit-intermediate"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) == 0:
        raise AssertionError("Expected non-zero exit on --vcf + --emit-intermediate")
    help_txt = subprocess.check_output([p.vdr_bin, "--help"], text=True)
    for hf in ("--vcfall", "--noend", "--vcfpass"):
        if hf in help_txt:
            raise AssertionError(f"Hidden alias {hf} leaked into --help")

    # 3. Switch 18 (somatic_continue_after_bad_var) across 6 BAM pairs + StrongLOH + MUC19
    arg_secondary = {
        (140453033, "T"), (140453035, "A"), (140453036, "CA"), (140453065, "T"),
        (140453081, "T"), (140453095, "G"), (140453108, "A"), (140453116, "G"),
    }
    sw18_cases = [
        ("ArgumentNA_fisher", p.fa_hg19, f"{p.vj_bams}/ArgumentNA_Tumor.bam|{p.vj_bams}/ArgumentNA_Normal.bam",
         ["-R", "chr7:140453010-140453158", "-q", "25", "-f", "0.001", "-F", "0x100", "-k", "1", "-M", "20", "-p", "--fisher"],
         lambda pos, _r, alt: (pos, alt) in arg_secondary, 8),
        ("ArgumentNA_nofisher", p.fa_hg19, f"{p.vj_bams}/ArgumentNA_Tumor.bam|{p.vj_bams}/ArgumentNA_Normal.bam",
         ["-R", "chr7:140453010-140453158", "-q", "25", "-f", "0.001", "-F", "0x100", "-k", "1", "-M", "20", "-p"],
         lambda pos, _r, alt: (pos, alt) in arg_secondary, 8),
        ("Colo829_756309", p.fa_hg19, colo_pair, ["-R", "chr1:755917-756517", "-f", "0.01"],
         lambda pos, r, alt: pos == 756309 and r == "TTGT" and alt == "CTGG", 1),
        ("HCC2218_7811864", p.fa_hg19, f"{p.vj_bams}/HCC2218_BL_realigins.bam|{p.vj_bams}/HCC2218_C_realigins.bam",
         ["-R", "chrX:7811670-7811868", "-f", "0.001"],
         lambda pos, r, alt: pos == 7811864 and r == "C" and alt == "T", 1),
        ("Zero_in_DUPSV_291259", p.fa_hg19, f"{p.vj_bams}/Zero_in_DUPSV_TUMOUR.bam|{p.vj_bams}/Zero_in_DUPSV_NORMAL.bam",
         ["-R", "chr19:290742-291636", "-f", "0.001", "-X", "3"],
         lambda pos, r, alt: pos == 291259 and r == "A" and alt == "G", 1),
        ("StrongLOH_11071401", p.fa_hg1k, f"{p.vj_bams}/TumorA.StrongLOH.chr21.bam|{p.vj_bams}/NormalA.StrongLOH.chr21.bam",
         ["-R", "chr21:10996701-11186054", "-Q", "10", "-f", "0.01", "-r", "2", "-F", "0x500", "-V", "0.05", "--nosv"],
         lambda pos, r, alt: pos in (11071401, 11071410) and len(alt) == len(r) + 1, 1),
    ]
    for cname, rfa, bpair, extra, pred, min_hits in sw18_cases:
        dv = cache.run_text(f"w4_sw18_vj_{cname}", rfa, bpair, ["--compat=vardictjava", "--vcf-all-alleles", *extra], kind="somatic", audit=False)
        ov = cache.run_text(f"w4_sw18_off_{cname}", rfa, bpair, ["--compat=off", "--vcf-all-alleles", *extra], kind="somatic")
        dh = [1 for ln in dv.splitlines() if not ln.startswith("#") and ln.strip() and pred(int(ln.split("\t")[1]), ln.split("\t")[3], ln.split("\t")[4])]
        oh = [1 for ln in ov.splitlines() if not ln.startswith("#") and ln.strip() and pred(int(ln.split("\t")[1]), ln.split("\t")[3], ln.split("\t")[4])]
        if dh or len(oh) < min_hits:
            raise AssertionError(f"Switch 18 {cname} failed: def_hits={len(dh)}, off_hits={len(oh)}")

    loh_var = cache.run_text(
        "w4_strong_loh_var55",
        p.fa_hg1k,
        f"{p.vj_bams}/TumorA.StrongLOH.chr21.bam|{p.vj_bams}/NormalA.StrongLOH.chr21.bam",
        ["--compat=off", "--emit-intermediate", "-R", "chr21:10996701-11186054", "-Q", "10", "-f", "0.01", "-r", "2", "-F", "0x500", "-V", "0.05", "--nosv"],
        kind="somatic",
    )
    loh_rows = [ln.split("\t") for ln in loh_var.splitlines() if ln.strip() and ln.split("\t")[49] == "StrongLOH"]
    if not loh_rows or any(r[52] != "0" for r in loh_rows):
        raise AssertionError("StrongLOH 55-col SV_info1 check failed")

    muc_dir = os.path.join(w4, "muc19")
    muc_gz = cache.run_vcf_gz(
        "w4_muc19",
        os.path.join(muc_dir, "hg19_chr12_muc19.fa"),
        f"{muc_dir}/MUC19_vaf1.bam|{muc_dir}/MUC19_vaf1_dz.bam",
        ["-c", "1", "-S", "2", "-E", "3", "-f", "0.005", "--fisher", "--compat=off", "--vcf-all-alleles", os.path.join(muc_dir, "MUC19_vaf1.bed")],
        kind="somatic",
    )
    with gzip.open(muc_gz, "rt", encoding="utf-8") as fh:
        muc_n = sum(1 for ln in fh if not ln.startswith("#") and ln.strip())
    if muc_n < 20:
        raise AssertionError(f"MUC19 expected >= 20 calls, got {muc_n}")

    # 4. Switch 19 (BRCA2 adaptor_insertion.bam)
    brca2_bam = os.path.join(w4, "brca2_adaptor", "adaptor_insertion.bam")
    adaptors = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    brca_base = ["-R", "chr13:32683877-32684023", "-f", "0.01", "-q", "25", "-r", "2", "-UN"]
    brca_vj = cache.run_text("w4_brca2_vj", p.grch38_ref, brca2_bam, [*brca_base, "--compat=vardictjava", "--adaptor", adaptors], audit=False)
    brca_off = cache.run_text("w4_brca2_off", p.grch38_ref, brca2_bam, [*brca_base, "--compat=off", "--adaptor", adaptors])
    if "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" not in brca_vj or "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" in brca_off:
        raise AssertionError("BRCA2 chr13:32683955 TruSeq adapter chimera check failed")

    # 5. Switch 20 (canonical_sv_vcf_anchoring across 8 Colo829 SVs + --noend + HG008 v0.6 + vartovcf_sv_contig1)
    sv_cases = [
        ("Colo829_DELDiscReverse", p.fa_hg19, f"{p.vj_bams}/Colo829.chr10.DELDiscReverse.bam", ["-R", "chr10:89700380-89713976", "-f", "0.001"]),
        ("Colo829_LGDEL3", p.fa_hg19, f"{p.vj_bams}/Colo829.chr10.LGDEL3.bam", ["-R", "chr10:127187158-127201055", "-f", "0.001"]),
        ("Colo829_DELSV3Unhappy", p.fa_hg19, f"{p.vj_bams}/Colo829.chr20.DELSV3Unhappy.bam", ["-R", "chr20:54434600-54441423", "-f", "0.001"]),
        ("Colo829_DUPDiscSV", p.fa_hg19, f"{p.vj_bams}/Colo829.chr20.DUPDiscSV.bam", ["-R", "chr20:47124929-47131802", "-f", "0.001", "-L", "500"]),
        ("Colo829_INVSV", p.fa_hg19, f"{p.vj_bams}/Colo829.chr20.INVSV.bam", ["-R", "chr20:46841554-46845755", "-f", "0.001"]),
        ("Colo829_MarkingDUPSV", p.fa_hg19, f"{p.vj_bams}/Colo829.chr20.MarkingDUPSV.bam", ["-R", "chr20:47131954-47132045", "-f", "0.001", "-L", "100"]),
        ("Colo829_CandDel5SV", p.fa_hg19, f"{p.vj_bams}/Colo829.chr20.CandDel5SV.bam", ["-R", "chr20:506395-507295", "-f", "0.001", "-L", "100"]),
        ("Zero_in_DUPSV_Paired", p.fa_hg19, f"{p.vj_bams}/Zero_in_DUPSV_TUMOUR.bam|{p.vj_bams}/Zero_in_DUPSV_NORMAL.bam",
         ["-R", "chr19:290742-291636", "-f", "0.001", "-X", "3", "-L", "50", "--fisher"]),
    ]
    sv_verified = 0
    for sname, rfa, sbam, sextra in sv_cases:
        sgz = cache.run_vcf_gz(f"w4_sv_{sname}", rfa, sbam, ["--compat=off", "--vcf-all-alleles", *sextra], kind=("somatic" if "|" in sbam else "simple"))
        with gzip.open(sgz, "rt", encoding="utf-8") as fh:
            n_sv = sum(1 for ln in fh if not ln.startswith("#") and ln.strip() and ln.split("\t")[4].startswith("<"))
        if n_sv < 1:
            raise AssertionError(f"{sname}: expected >= 1 symbolic SV")
        sv_verified += n_sv

    noend_vcf = cache.run_text(
        "w4_sv_noend", p.fa_hg19, f"{p.vj_bams}/Colo829.chr10.DELDiscReverse.bam",
        ["--compat=off", "--noend", "-R", "chr10:89700380-89713976", "-f", "0.001"],
    )
    if "\t<DEL>\t" not in noend_vcf or ";END=" in noend_vcf:
        raise AssertionError("--noend symbolic SV check failed")

    v2v_dir = os.path.join(w4, "vartovcf")
    v2v_gz = cache.run_vcf_gz(
        "w4_vartovcf_contig1",
        os.path.join(v2v_dir, "reference.fa"),
        os.path.join(v2v_dir, "vartovcf_sv_contig1.bam"),
        ["-c", "1", "-S", "2", "-E", "3", "-g", "4", "-P", "0", "-f", "0.01", "-r", "2", "--compat=off", "--vcf-all-alleles", os.path.join(v2v_dir, "vartovcf_sv_contig1.bed")],
    )
    got_v2v = {}
    with gzip.open(v2v_gz, "rt", encoding="utf-8") as fh:
        for ln in fh:
            if not ln.startswith("#") and ln.strip():
                sp = ln.rstrip("\n").split("\t")
                if sp[4].startswith("<"):
                    kv = dict(x.split("=", 1) for x in sp[7].split(";") if "=" in x)
                    got_v2v[(sp[0], int(sp[1]), sp[3], sp[4])] = (int(kv["END"]), int(kv["SVLEN"]))
    if got_v2v.get(("chr12", 446, "C", "<DUP>")) != (1446, 1000) or got_v2v.get(("chr13", 1141, "G", "<DEL>")) != (4141, -3000):
        raise AssertionError(f"vartovcf_sv_contig1 SV anchoring mismatch: {got_v2v}")

    # 6. QIAseq V3 SPE (-sp) >= 3x DP recovery & 0 KnownNegatives PASS FPs
    qgn_bed = os.path.join(w4, "beds", "qiaseq_spe_grch38.amp8.bed")
    bam_qgn_a = os.path.join(w4, "bam", "SampleA_QGN1.bam")
    bam_qgn_ac5 = os.path.join(w4, "bam", "SampleAC5_QGN1.bam")
    kp_map: dict[str, list[tuple[int, str, str]]] = {}
    with gzip.open(os.path.join(w4, "raw", "KnownPositives_hg19ToHg38.vcf.gz"), "rt", encoding="utf-8") as fh:
        for ln in fh:
            if not ln.startswith("#"):
                sp = ln.rstrip("\n").split("\t")
                kp_map.setdefault(sp[0], []).append((int(sp[1]), sp[3], sp[4]))
    kn_map: dict[str, list[tuple[int, int]]] = {}
    with gzip.open(os.path.join(w4, "raw", "KnownNegatives_hg38.bed.gz"), "rt", encoding="utf-8") as fh:
        for ln in fh:
            if not ln.startswith("#") and ln.strip():
                sp = ln.rstrip("\n").split("\t")
                kn_map.setdefault(sp[0], []).append((int(sp[1]), int(sp[2])))
    for c in kn_map:
        kn_map[c].sort()

    def in_kn(chrom: str, pos1: int) -> bool:
        ivs = kn_map.get(chrom, [])
        idx = bisect.bisect_left(ivs, (pos1, -1))
        for j in range(max(0, idx - 2), min(len(ivs), idx + 3)):
            if ivs[j][0] < pos1 <= ivs[j][1]:
                return True
        return False

    def eval_spe(tag: str, bam_p: str, f_str: str, with_sp: bool, kn_max_af: float) -> tuple[int, int, int]:
        fl = ["-a", "10:0.95", "-UN", "-f", f_str, "-r", "3", "-q", "25", "--compat=off", "--vcf-all-alleles"]
        if with_sp:
            fl.append("-sp")
        fl.append(qgn_bed)
        txt = cache.run_text(tag, p.grch38_ref, bam_p, fl, kind="amplicon")
        dp_sum, pass_kp, kn_fp = 0, 0, 0
        for ln in txt.splitlines():
            if not ln.startswith("#") and ln.strip():
                sp = ln.split("\t")
                kv = dict(x.split("=", 1) for x in sp[7].split(";") if "=" in x)
                dp_sum += int(kv["DP"])
                af = float(kv["AF"])
                if sp[6] == "PASS":
                    if any(v[0] == int(sp[1]) and v[2] == sp[4] for v in kp_map.get(sp[0], [])):
                        pass_kp += 1
                    if af < kn_max_af and in_kn(sp[0], int(sp[1])):
                        kn_fp += 1
        return dp_sum, pass_kp, kn_fp

    dp_no_a, kp_no_a, _ = eval_spe("w4_spe_a_nosp", bam_qgn_a, "0.02", False, 1.01)
    dp_sp_a, kp_sp_a, kn_a = eval_spe("w4_spe_a_sp", bam_qgn_a, "0.02", True, 1.01)
    dp_no_ac5, kp_no_ac5, _ = eval_spe("w4_spe_ac5_nosp", bam_qgn_ac5, "0.0025", False, 0.10)
    dp_sp_ac5, kp_sp_ac5, kn_ac5 = eval_spe("w4_spe_ac5_sp", bam_qgn_ac5, "0.0025", True, 0.10)
    if dp_sp_a < 3 * max(1, dp_no_a) or kp_sp_a < 30 or kn_a != 0:
        raise AssertionError(f"QIAseq SampleA -sp gate failed: dp={dp_no_a}->{dp_sp_a}, kp={kp_sp_a}, kn={kn_a}")
    if dp_sp_ac5 < 3 * max(1, dp_no_ac5) or kp_sp_ac5 < 20 or kn_ac5 != 0:
        raise AssertionError(f"QIAseq SampleAC5 -sp gate failed: dp={dp_no_ac5}->{dp_sp_ac5}, kp={kp_sp_ac5}, kn={kn_ac5}")

    # 7. Switch 21 (paired_somatic_amplicon across QIAseq SPE, Illumina PE BAM==CRAM, Ion Torrent SE)
    ilm_bed = os.path.join(w4, "beds", "ilm_twoprimer_grch38.amp8.bed")
    ilm_primer_bed = os.path.join(w4, "beds", "ilm_primer_overlap_grch38.amp8.bed")
    tfs_bed = os.path.join(w4, "beds", "tfs_iontorrent_grch38.amp8.bed")
    if subprocess.call([p.vdr_bin, "-G", p.grch38_ref, "-b", f"{w4}/bam/SampleA_ILM1.bam|{w4}/bam/SampleB_ILM1.bam", "-a", "40:0.70", "--compat=vardictjava", ilm_bed], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) == 0:
        raise AssertionError("Expected --compat=vardictjava to reject paired somatic amplicon VCF")

    qgn_pair_gz = cache.run_vcf_gz(
        "w4_qgn_paired", p.grch38_ref, f"{w4}/bam/SampleA_QGN1.bam|{w4}/bam/SampleB_QGN1.bam",
        ["-N", "SampleA_QGN1|SampleB_QGN1", "-a", "10:0.95", "-sp", "-UN", "-f", "0.02", "-r", "3", "-q", "25", "--compat=off", "--vcf-all-alleles", "--fisher", qgn_bed],
        kind="somatic",
    )
    ilm_bam_gz = cache.run_vcf_gz(
        "w4_ilm_bam_paired", p.grch38_ref, f"{w4}/bam/SampleA_ILM1.bam|{w4}/bam/SampleB_ILM1.bam",
        ["-N", "SampleA_ILM1|SampleB_ILM1", "-a", "40:0.70", "-f", "0.02", "-q", "25", "--compat=off", "--vcf-all-alleles", "--fisher", ilm_bed],
        kind="somatic",
    )
    ilm_cram_gz = cache.run_vcf_gz(
        "w4_ilm_cram_paired", p.grch38_ref, f"{w4}/bam/SampleA_ILM1.cram|{w4}/bam/SampleB_ILM1.cram",
        ["-N", "SampleA_ILM1|SampleB_ILM1", "-a", "40:0.70", "-f", "0.02", "-q", "25", "--compat=off", "--vcf-all-alleles", "--fisher", ilm_bed],
        index_flag="--csi",
        kind="somatic",
    )
    with gzip.open(ilm_bam_gz, "rt", encoding="utf-8") as fb, gzip.open(ilm_cram_gz, "rt", encoding="utf-8") as fc:
        ib_rows = [ln.split("\t")[:8] for ln in fb if not ln.startswith("#") and ln.strip()]
        ic_rows = [ln.split("\t")[:8] for ln in fc if not ln.startswith("#") and ln.strip()]
    if not ib_rows or ib_rows != ic_rows:
        raise AssertionError("Illumina PE paired somatic amplicon BAM != CRAM")

    with open(ilm_primer_bed, "r", encoding="utf-8") as fh:
        p_ins_st = int(fh.readline().strip().split("\t")[6])
    primer_txt = cache.run_text(
        "w4_ilm_primer_mask", p.grch38_ref, f"{w4}/bam/SampleA_ILM1.bam|{w4}/bam/SampleB_ILM1.bam",
        ["-a", "40:0.70", "-f", "0.01", "-q", "25", "--compat=off", "--vcf-all-alleles", ilm_primer_bed],
        kind="somatic",
    )
    if any(int(ln.split("\t")[1]) <= p_ins_st for ln in primer_txt.splitlines() if not ln.startswith("#") and ln.strip()):
        raise AssertionError("Primer-overlapping call leaked in paired somatic amplicon mode")

    tfs_gz = cache.run_vcf_gz(
        "w4_tfs_paired", p.grch38_ref, f"{w4}/bam/SampleA_TFS1.bam|{w4}/bam/SampleB_TFS1.bam",
        ["-a", "10:0.95", "-sp", "-f", "0.02", "-q", "20", "--compat=off", "--vcf-all-alleles", "--fisher", tfs_bed],
        kind="somatic",
    )
    for want_cols, ex in ((55, []), (61, ["--fisher"])):
        vt = cache.run_text(
            f"w4_tfs_var_{want_cols}", p.grch38_ref, f"{w4}/bam/SampleA_TFS1.bam|{w4}/bam/SampleB_TFS1.bam",
            ["-a", "10:0.95", "-sp", "-f", "0.02", "-q", "20", "--compat=off", "--emit-intermediate", "-h", *ex, tfs_bed],
            kind="somatic",
        )
        vlines = [ln for ln in vt.splitlines() if ln.strip()]
        if len(vlines[0].split("\t")) != 55 or any(len(r.split("\t")) != want_cols for r in vlines[1:]):
            raise AssertionError(f"Ion Torrent paired .var column check failed for {want_cols}")

    with gzip.open(qgn_pair_gz, "rt", encoding="utf-8") as fq, gzip.open(tfs_gz, "rt", encoding="utf-8") as ft:
        qgn_n = sum(1 for ln in fq if not ln.startswith("#") and ln.strip())
        tfs_n = sum(1 for ln in ft if not ln.startswith("#") and ln.strip())
    with metrics_lock:
        metrics["wave4"] = {
            "cmrg_mosaic_hits": matched_cmrg,
            "muc19_rows": muc_n,
            "sv_verified": sv_verified,
            "qiaseq_sp_kp_a": kp_sp_a,
            "qiaseq_sp_kp_ac5": kp_sp_ac5,
            "qiaseq_paired_rows": qgn_n,
            "illumina_paired_rows": len(ib_rows),
            "iontorrent_paired_rows": tfs_n,
        }
    print(
        f"  [PASS Wave-4] CMRG/Mosaic={matched_cmrg}, MUC19={muc_n}, SVs={sv_verified}, "
        f"QIAseq KP={kp_sp_a}/{kp_sp_ac5} (T/N={qgn_n}), Illumina PE BAM==CRAM ({len(ib_rows)}), IonTorrent T/N={tfs_n}"
    )


# =============================================================================
# Suite 2: Wave-5 (7 Clinical & Validator Gates)
# =============================================================================
def run_suite_wave5(
    p: Paths,
    cache: RunCache,
    rtg_eval: RtgEvaluator,
    metrics: dict[str, object],
    metrics_lock: threading.Lock,
) -> None:
    print("[Wave-5] Running 7 clinical & validator gates...")
    g4 = str(p.g4_dir)
    w4 = str(p.wave4_dir)

    # Gate 1: Multi-BGZF-block --csi spanning-SV tabix & bcftools view -r
    lgdel_gz = cache.run_vcf_gz(
        "w5_csi_lgdel",
        p.fa_hg19,
        f"{p.vj_bams}/Colo829.chr10.DELDiscReverse.bam",
        ["-R", "chr10:89700380-89713976", "-f", "0.001", "-p", "--compat=off", "--vcf-all-alleles"],
        index_flag="--csi",
        is_pileup=True,
    )
    lgdel_blocks = count_bgzf_blocks(lgdel_gz)
    if lgdel_blocks < 8:
        raise AssertionError(f"Expected >= 8 BGZF blocks in {lgdel_gz}, got {lgdel_blocks}")

    # Build spanning >35 kbp <DEL>/<DUP>/<INV> chr20 BAM across >=3 CSI leaf-bin boundaries
    fa_chr20 = subprocess.check_output(["samtools", "faidx", p.fa_hg19, "chr20:46839001-46915000"], text=True)
    chr20_seq = "".join(ln.strip() for ln in fa_chr20.splitlines() if not ln.startswith(">")).upper()
    chr20_origin = 46839001
    chr20_len = read_fai_lengths(p.fa_hg19)["chr20"]

    def ref_slice(p1: int, length: int) -> str:
        off = p1 - chr20_origin
        return chr20_seq[off : off + length]

    tr_rc = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    sam_span3 = os.path.join(p.work_tmp, "span3_sv.sam")
    bam_span3 = os.path.join(p.work_tmp, "span3_sv.bam")
    recs_s3 = []
    rid = 0
    for pos in range(46841000, 46896000, 35):
        for sflag in (0, 16):
            rid += 1
            recs_s3.append((pos, f"wt_{rid}\t{sflag}\tchr20\t{pos}\t60\t100M\t*\t0\t0\t{ref_slice(pos, 100)}\t{'I'*100}"))
    del_l, del_r = 46841600, 46885000
    for i in range(10):
        s1 = del_l - 60
        recs_s3.append((s1, f"sdel_{i}\t99\tchr20\t{s1}\t60\t60M30S\t=\t{del_r}\t{del_r + 60 - s1}\t{ref_slice(s1, 60) + ref_slice(del_r, 30)}\t{'I'*90}"))
        recs_s3.append((del_r, f"sdel_{i}\t147\tchr20\t{del_r}\t60\t30S60M\t=\t{s1}\t{-(del_r + 60 - s1)}\t{ref_slice(del_l - 30, 30) + ref_slice(del_r, 60)}\t{'I'*90}"))
    dup_l, dup_r = 46842200, 46887500
    for i in range(10):
        recs_s3.append((dup_l, f"sdup_{i}\t81\tchr20\t{dup_l}\t60\t30S70M\t=\t{dup_r - 69}\t{dup_r - dup_l + 1}\t{ref_slice(dup_r - 29, 30) + ref_slice(dup_l, 70)}\t{'I'*100}"))
        sf = dup_r - 69
        recs_s3.append((sf, f"sdup_{i}\t161\tchr20\t{sf}\t60\t70M30S\t=\t{dup_l}\t{-(dup_r - dup_l + 1)}\t{ref_slice(sf, 70) + ref_slice(dup_l, 30)}\t{'I'*100}"))
    inv_l, inv_r = 46843000, 46889500
    for i in range(10):
        sf = inv_r - 69
        inv_tail = ref_slice(inv_l - 19, 20).translate(tr_rc)[::-1]
        recs_s3.append((inv_l, f"sinv_m_{i}\t1\tchr20\t{inv_l}\t60\t80M\t=\t{sf}\t{sf + 70 - inv_l}\t{ref_slice(inv_l, 80)}\t{'I'*80}"))
        recs_s3.append((sf, f"sinv_s_{i}\t1\tchr20\t{sf}\t60\t70M20S\t=\t{inv_l}\t{-(sf + 70 - inv_l)}\t{ref_slice(sf, 70) + inv_tail}\t{'I'*90}"))
    recs_s3.sort(key=lambda x: x[0])
    with open(sam_span3, "w", encoding="utf-8") as out:
        out.write(f"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr20\tLN:{chr20_len}\n")
        for _, ln in recs_s3:
            out.write(ln + "\n")
    subprocess.run(["samtools", "view", "-b", "-o", bam_span3, sam_span3], check=True)
    subprocess.run(["samtools", "index", bam_span3], check=True)

    span3_gz = cache.run_vcf_gz(
        "w5_csi_span3", p.fa_hg19, bam_span3,
        ["-R", "chr20:46841200-46895000", "-f", "0.001", "-p", "--compat=off", "--vcf-all-alleles"],
        index_flag="--csi",
        is_pileup=True,
    )
    span3_blocks = count_bgzf_blocks(span3_gz)
    if span3_blocks < 20:
        raise AssertionError(f"Expected >= 20 BGZF blocks in {span3_gz}, got {span3_blocks}")
    found_span_svs = {}
    with gzip.open(span3_gz, "rt", encoding="utf-8") as fh:
        for ln in fh:
            if not ln.startswith("#"):
                sp = ln.rstrip("\n").split("\t")
                if sp[4] in ("<DEL>", "<DUP>", "<INV>"):
                    kv = dict(x.split("=", 1) for x in sp[7].split(";") if "=" in x)
                    if int(kv["END"]) - int(sp[1]) >= 35000:
                        found_span_svs[sp[4]] = (sp[0], int(sp[1]), int(kv["END"]), sp[4])
    for req_sv in ("<DEL>", "<DUP>", "<INV>"):
        if req_sv not in found_span_svs:
            raise AssertionError(f"Missing >35kbp spanning {req_sv} in {span3_gz}")
        chrom, s_pos, e_pos, sv_alt = found_span_svs[req_sv]
        for b in range(((s_pos - 1) >> 14) + 1, ((e_pos - 1) >> 14) + 1):
            q_coord = min(e_pos - 1, (b << 14) + 500)
            reg = f"{chrom}:{q_coord}-{q_coord}"
            tbx = subprocess.check_output(["tabix", span3_gz, reg], text=True)
            bcf = subprocess.check_output(["bcftools", "view", "-H", "-r", reg, span3_gz], text=True)
            for tname, txt in (("tabix", tbx), ("bcftools", bcf)):
                if not any(ln.split("\t")[1] == str(s_pos) and ln.split("\t")[4] == sv_alt for ln in txt.splitlines() if ln.strip()):
                    raise AssertionError(f"{tname} dropped spanning {sv_alt} at {reg}")

    # Gate 2: SEQ-matched 90M3D15M3I vs 90M3D15M3I20H clinical BAM twins
    for tag, src_bam, reg, anchor in (
        ("seqc2_hclip", os.path.join(g4, "seqc2", "HCC1395.bam"), "chr17:74747665-74755000", 74750100),
        ("qiaseq_hclip", os.path.join(w4, "bam", "SampleA_QGN1.bam"), "chr1:11144542-11145012", 11144620),
    ):
        fa_txt = subprocess.check_output(["samtools", "faidx", p.grch38_ref, f"{reg.split(':')[0]}:{anchor}-{anchor + 149}"], text=True)
        ref_seg = "".join(ln.strip() for ln in fa_txt.splitlines() if not ln.startswith(">")).upper()
        sam_raw = subprocess.check_output(["samtools", "view", "-h", src_bam, reg], text=True)
        hdr_l, aln_c, aln_h = [], [], []
        for ln in sam_raw.splitlines():
            if ln.startswith("@"):
                hdr_l.append(ln)
                continue
            cols = ln.split("\t")
            if len(cols) < 11 or cols[5] == "*":
                continue
            ops = [(int(n), op) for n, op in CIGAR_RE.findall(cols[5])]
            while ops and ops[0][1] in ("H", "P"):
                ops.pop(0)
            while ops and ops[-1][1] in ("H", "P"):
                ops.pop()
            ops = [(n, op) for n, op in ops if op != "P"]
            if not ops:
                continue
            ccig = "".join(f"{n}{op}" for n, op in ops)
            c1, c2 = list(cols), list(cols)
            c1[5], c2[5] = ccig, f"15H{ccig}20H"
            aln_c.append((int(cols[3]), "\t".join(c1)))
            aln_h.append((int(cols[3]), "\t".join(c2)))
        ins_b = "ACG" if ref_seg[108:111] != "ACG" else "TGC"
        tseq = ref_seg[0:90] + ref_seg[93:108] + ins_b
        chrom = reg.split(":")[0]
        for idx in range(60):
            fl = 16 if (idx % 2 == 1) else 0
            rn = f"twin_{tag}_{idx}"
            aln_c.append((anchor, f"{rn}\t{fl}\t{chrom}\t{anchor}\t60\t90M3D15M3I\t*\t0\t0\t{tseq}\t{'I'*108}"))
            chc = "90M3D15M3I20H" if idx < 30 else "12H90M3D15M3I20H"
            aln_h.append((anchor, f"{rn}\t{fl}\t{chrom}\t{anchor}\t60\t{chc}\t*\t0\t0\t{tseq}\t{'I'*108}"))
        aln_c.sort(key=lambda x: x[0])
        aln_h.sort(key=lambda x: x[0])
        cbam = os.path.join(p.work_tmp, f"{tag}_clean.bam")
        hbam = os.path.join(p.work_tmp, f"{tag}_hclip.bam")
        for obam, rlist in ((cbam, aln_c), (hbam, aln_h)):
            sp = obam + ".sam"
            with open(sp, "w", encoding="utf-8") as fh:
                for h in hdr_l:
                    fh.write(h + "\n")
                for _, r in rlist:
                    fh.write(r + "\n")
            subprocess.run(["samtools", "view", "-b", "-o", obam, sp], check=True)
            subprocess.run(["samtools", "index", obam], check=True)
        vj_c = cache.run_text(f"{tag}_vjc", p.grch38_ref, cbam, ["-R", reg, "-N", tag, "-f", "0.005", "--compat=vardictjava"], audit=False)
        vj_h = cache.run_text(f"{tag}_vjh", p.grch38_ref, hbam, ["-R", reg, "-N", tag, "-f", "0.005", "--compat=vardictjava"], audit=False)
        off_c = cache.run_text(f"{tag}_offc", p.grch38_ref, cbam, ["-R", reg, "-N", tag, "-f", "0.005", "--compat=off", "--vcf-all-alleles", "--fisher"])
        off_h = cache.run_text(f"{tag}_offh", p.grch38_ref, hbam, ["-R", reg, "-N", tag, "-f", "0.005", "--compat=off", "--vcf-all-alleles", "--fisher"])
        if vj_c == vj_h or off_c != off_h:
            raise AssertionError(f"Gate 2 hard-clip twin invariance failed on {tag}")

    # Gate 3: 6-Cohort .var + .vcf + survival audit (reusing cached cohorts!)
    qgn_bed = os.path.join(w4, "beds", "qiaseq_spe_grch38.amp8.bed")
    ilm_bed = os.path.join(w4, "beds", "ilm_twoprimer_grch38.amp8.bed")
    tfs_bed = os.path.join(w4, "beds", "tfs_iontorrent_grch38.amp8.bed")
    g3_cohorts = [
        ("g3_cmrg", "simple", p.giab_ref, os.path.join(w4, "hg002_cmrg_mosaic", "HG002.cmrg_mosaic.slice.bam"),
         ["-N", "HG002", "-c", "1", "-S", "2", "-E", "3", "-g", "4", "-z", "-f", "0.01", os.path.join(w4, "hg002_cmrg_mosaic", "hg002_cmrg_mosaic.bed")]),
        ("g3_hg008", "somatic", p.giab_ref, f"{g4}/hg008/HG008-T.slice.bam|{g4}/hg008/HG008-N.slice.bam",
         ["-N", "HG008T|HG008N", "-c", "1", "-S", "2", "-E", "3", "-z", "-f", "0.01", os.path.join(g4, "hg008", "slices.bed")]),
        ("g3_seqc2", "somatic", p.grch38_ref, f"{g4}/seqc2/HCC1395.bam|{g4}/seqc2/HCC1395BL.bam",
         ["-N", "TUMOR|NORMAL", "-c", "1", "-S", "2", "-E", "3", "-g", "4", "-z", "-f", "0.01", os.path.join(g4, "seqc2", "slices.bed")]),
        ("g3_qiaseq_sp", "amplicon", p.grch38_ref, os.path.join(w4, "bam", "SampleA_QGN1.bam"),
         ["-a", "10:0.95", "-sp", "-z", "-f", "0.005", qgn_bed]),
        ("g3_ilm_tn", "somatic", p.grch38_ref, f"{w4}/bam/SampleA_ILM1.bam|{w4}/bam/SampleB_ILM1.bam",
         ["-N", "SampleA_ILM1|SampleB_ILM1", "-a", "10:0.95", "-z", "-f", "0.01", ilm_bed]),
        ("g3_tfs_tn", "somatic", p.grch38_ref, f"{w4}/bam/SampleA_TFS1.bam|{w4}/bam/SampleB_TFS1.bam",
         ["-N", "SampleA_TFS1|SampleB_TFS1", "-a", "10:0.95", "-sp", "-z", "-f", "0.01", tfs_bed]),
    ]
    g3_var_total = 0
    for cname, ckind, rfa, bspec, cextra in g3_cohorts:
        vt = cache.run_text(f"{cname}_var", rfa, bspec, ["--compat=off", "--vcf-all-alleles", "--emit-intermediate", "--fisher", *cextra], kind=ckind)
        cgz = cache.run_vcf_gz(f"{cname}_vcf", rfa, bspec, ["--compat=off", "--vcf-all-alleles", "--fisher", *cextra], kind=ckind)
        with gzip.open(cgz, "rt", encoding="utf-8") as gzf:
            crecs = [ln.rstrip("\n") for ln in gzf if ln.strip() and not ln.startswith("#")]
        vrows = audit_compat_off_var(cname, vt, ckind, fisher=True, vcf_recs=crecs, vcf_all_alleles=True)
        g3_var_total += len(vrows)

    # Gate 4: Real GRCh38.d1.vd1.fa IUPAC reference sites (chr6:61366327 Y, chr21:41829696 M) on HG002.GRCh38.2x250.bam + SV bracket safety
    iupac_bed = os.path.join(p.work_tmp, "d1vd1_iupac.bed")
    with open(iupac_bed, "w", encoding="utf-8") as fh:
        fh.write("chr6\t61366200\t61366500\tiupac_Y\nchr21\t41829550\t41829850\tiupac_M\n")
    iupac_2x250_gz = cache.run_vcf_gz(
        "w5_d1vd1_iupac_2x250",
        p.grch38_ref,
        os.path.join(g4, "hg002", "HG002.GRCh38.2x250.bam"),
        ["-c", "1", "-S", "2", "-E", "3", "-g", "4", "-z", "-f", "0.01", "-p", "--compat=off", "--vcf-all-alleles", iupac_bed],
        index_flag="--csi",
        query_region="chr21:41829696-41829696",
        is_pileup=True,
    )
    tbx_y = subprocess.check_output(["tabix", iupac_2x250_gz, "chr6:61366327-61366327"], text=True).strip()
    tbx_m = subprocess.check_output(["tabix", iupac_2x250_gz, "chr21:41829696-41829696"], text=True).strip()
    if not tbx_y or not tbx_m or tbx_m.split("\t")[3] not in ("A", "C", "G", "T"):
        raise AssertionError(f"Real GRCh38.d1.vd1.fa IUPAC site reduction failed: Y={tbx_y!r}, M={tbx_m!r}")

    # Gate 5: Bounded -L (--min-sv-len) materialization, discordant suppression, & -L 1000000 32,768 bp ceiling
    for sname, rfa, sbam, sreg in (
        ("Colo829_DUP843", p.fa_hg19, f"{p.vj_bams}/Colo829.chr20.DUPDiscSV.bam", "chr20:47124929-47131802"),
        ("Colo829_MarkingDUP", p.fa_hg19, f"{p.vj_bams}/Colo829.chr20.MarkingDUPSV.bam", "chr20:47131954-47132045"),
        ("NA12878_INVSimple", p.fa_hs38dh, f"{p.vj_bams}/NA12878.chrom20.INVSimpleSV.bam", "20:1961100-1961400"),
        ("NA12878_INVFull", p.fa_hs38dh, f"{p.vj_bams}/NA12878.chrom20.INVFullSV.bam", "20:5926940-5927179"),
    ):
        mgz = cache.run_vcf_gz(f"w5_mat_{sname}", rfa, sbam, ["-R", sreg, "-f", "0.001", "--compat=off", "--vcf-all-alleles"])
        with gzip.open(mgz, "rt", encoding="utf-8") as fh:
            mrows = [ln.rstrip("\n").split("\t") for ln in fh if not ln.startswith("#") and ln.strip()]
        if any("<" in r[4] for r in mrows) or not any(max(len(r[3]), len(r[4])) >= 20 and DNA_RE.match(r[4]) for r in mrows):
            raise AssertionError(f"{sname}: sub-1kb SV materialization check failed")

    for dname, rfa, dbam, dreg, sup_l, res_l in (
        ("NA21144_DUP", p.fa_hs38dh, f"{p.vj_bams}/NA21144.chrom11.DUPDiscF.bam", "11:77614410-77614671", "1000", "50"),
        ("NA21144_INV", p.fa_hs38dh, f"{p.vj_bams}/NA21144.chrom11.INV5Disc.bam", "11:57970876-57971275", "15000", "1000"),
    ):
        vsup = cache.run_text(f"w5_sup_{dname}", rfa, dbam, ["-R", dreg, "-f", "0.001", "-L", sup_l, "--compat=off", "--vcf-all-alleles"])
        vres = cache.run_text(f"w5_res_{dname}", rfa, dbam, ["-R", dreg, "-f", "0.001", "-L", res_l, "--compat=off", "--vcf-all-alleles"])
        if any("<DUP>" in ln or "<INV>" in ln for ln in vsup.splitlines() if not ln.startswith("#")):
            raise AssertionError(f"{dname}: discordant SV < -L {sup_l} leaked")
        if not any("<DUP>" in ln or "<INV>" in ln for ln in vres.splitlines() if not ln.startswith("#")):
            raise AssertionError(f"{dname}: discordant SV >= -L {res_l} missing")

    ceil_gz = cache.run_vcf_gz(
        "w5_ceil_L1M", p.fa_hg19, bam_span3,
        ["-R", "chr20:46841200-46895000", "-f", "0.001", "-L", "1000000", "--compat=off", "--vcf-all-alleles"],
    )
    ceil_svs = set()
    with gzip.open(ceil_gz, "rt", encoding="utf-8") as fh:
        for ln in fh:
            if not ln.startswith("#") and ln.strip():
                sp = ln.rstrip("\n").split("\t")
                if max(len(sp[3]), len(sp[4])) > 32769:
                    raise AssertionError("-L 1000000 exceeded MAX_EXPLICIT_SV_SEQ_BP (32,768 bp)")
                if sp[4] in ("<DEL>", "<DUP>", "<INV>"):
                    ceil_svs.add(sp[4])
    if ceil_svs != {"<DEL>", "<DUP>", "<INV>"}:
        raise AssertionError(f"Expected >32.7kb SVs to stay symbolic under -L 1000000, got {ceil_svs}")

    # Gate 6: -p multi-deletion homopolymer coordinate retention
    for pname, rfa, pbam, chrom, sp_i, ep_i, tsz, full_span in (
        ("N13532", p.fa_b37, f"{p.vj_bams}/N13532.TwoDels1InsStringBoundError.bam", "11", 1279500, 1281000, "300", True),
        ("threeIndels", p.fa_b37, f"{p.vj_bams}/threeIndels_zero_deletion_tumor.bam", "20", 62361250, 62361500, "80", False),
        ("HG002_CMRG_p", p.giab_ref, os.path.join(w4, "hg002_cmrg_mosaic", "HG002.cmrg_mosaic.slice.bam"), "chr1", 1019800, 1021200, "350", True),
    ):
        reg = f"{chrom}:{sp_i}-{ep_i}"
        depth_out = subprocess.check_output(["samtools", "depth", "-Q", "1", "-r", reg, pbam], text=True)
        covered = {int(ln.split("\t")[1]) for ln in depth_out.splitlines() if ln.strip() and int(ln.split("\t")[2]) > 0}
        vu = cache.run_text(f"w5_p_u_{pname}", rfa, pbam, ["-R", reg, "-p", "--compat=off", "--vcf-all-alleles"], threads_override=1)
        vt = cache.run_text(f"w5_p_t_{pname}", rfa, pbam, ["-R", reg, "-p", "--compat=off", "--vcf-all-alleles", "--tile-size", tsz])
        if full_span and vu != vt:
            raise AssertionError(f"{pname}: -p untiled != tiled ({tsz})")
        emitted = {int(ln.split("\t")[1]) for ln in vu.splitlines() if not ln.startswith("#") and ln.strip()}
        if covered - emitted:
            raise AssertionError(f"{pname}: -p dropped covered coordinates: {sorted(covered - emitted)[:5]}")

    # Gate 7: Subclonal somatic GT=0/1 recovery without --squash-ploidy on SEQC2 (using on-disk GRCh38.d1.vd1.sdf!)
    slices_bed = os.path.join(g4, "seqc2", "slices.bed")
    scoped_eval = os.path.join(p.work_tmp, "w5_seqc2_eval.bed")
    intersect_beds(slices_bed, os.path.join(g4, "seqc2", "eval.bed"), scoped_eval)
    fai_lens = read_fai_lengths(p.grch38_ref)
    raw_truth = subprocess.check_output(["bcftools", "view", "-R", slices_bed, os.path.join(g4, "seqc2", "truth", "truth.vcf.gz")], text=True)
    dip_truth = os.path.join(p.work_tmp, "w5_seqc2_dip_truth.vcf")
    with open(dip_truth, "w", encoding="utf-8") as out:
        out.write("##fileformat=VCFv4.2\n")
        for c, length in fai_lens.items():
            out.write(f"##contig=<ID={c},length={length}>\n")
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\n")
        for ln in raw_truth.splitlines():
            if not ln.startswith("#") and ln.strip():
                sp = ln.split("\t")
                if sp[0] in fai_lens:
                    m = TVAF_RE.search(sp[7])
                    tvaf = float(m.group(1)) if m else 0.30
                    gt = "1/1" if tvaf >= 0.80 else "0/1"
                    out.write(f"{sp[0]}\t{sp[1]}\t{sp[2]}\t{sp[3]}\t{sp[4]}\t.\tPASS\t.\tGT\t{gt}\n")
    subprocess.run(["bgzip", "-f", dip_truth], check=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", dip_truth + ".gz"], check=True)

    seqc2_pair = f"{g4}/seqc2/HCC1395.bam|{g4}/seqc2/HCC1395BL.bam"
    tp_counts = {}
    gt00_counts = {}
    for mtag, cflag in (("vj", "--compat=vardictjava"), ("off", "--compat=off")):
        rtxt = cache.run_text(
            f"w5_g7_seqc2_{mtag}", p.grch38_ref, seqc2_pair,
            ["-N", "TUMOR|NORMAL", "-c", "1", "-S", "2", "-E", "3", "-g", "4", "-z", "-f", "0.01", "--vcf-af", "0.30", "--fisher", cflag, slices_bed],
            kind="somatic",
            audit=(cflag == "--compat=off"),
        )
        som_p = os.path.join(p.work_tmp, f"w5_g7_{mtag}.vcf")
        gt00 = 0
        with open(som_p, "w", encoding="utf-8") as out:
            for ln in rtxt.splitlines():
                if ln.startswith("#"):
                    out.write(ln + "\n")
                elif "STATUS=StrongSomatic" in ln or "STATUS=LikelySomatic" in ln:
                    sp = ln.split("\t")
                    tf, nf = sp[9].split(":"), sp[10].split(":")
                    if tf[0] == "0/0" and nf[0] == "0/0" and (int(tf[2]) > 0 or int(nf[2]) > 0):
                        gt00 += 1
                    out.write(ln + "\n")
        subprocess.run(["bgzip", "-f", som_p], check=True)
        subprocess.run(["tabix", "-f", "-p", "vcf", som_p + ".gz"], check=True)
        _, _, _, odir = rtg_eval.run_vcfeval(
            f"w5_g7_{mtag}", dip_truth + ".gz", som_p + ".gz", scoped_eval, p.sdf_d1vd1,
            squash_ploidy=False, sample_arg="TUMOR,TUMOR", all_records=True,
        )
        with gzip.open(os.path.join(odir, "tp.vcf.gz"), "rt", encoding="utf-8") as fh:
            tp_counts[mtag] = sum(1 for ln in fh if not ln.startswith("#"))
        gt00_counts[mtag] = gt00

    if gt00_counts["vj"] < 10 or gt00_counts["off"] != 0 or tp_counts["off"] <= tp_counts["vj"]:
        raise AssertionError(f"Gate 7 subclonal diploid GT recovery failed: tp={tp_counts}, gt00={gt00_counts}")

    with metrics_lock:
        metrics["wave5"] = {
            "lgdel_csi_blocks": lgdel_blocks,
            "span3_csi_blocks": span3_blocks,
            "gate3_var_rows": g3_var_total,
            "subclonal_tp_vj": tp_counts["vj"],
            "subclonal_tp_off": tp_counts["off"],
        }
    print(
        f"  [PASS Wave-5] CSI blocks={lgdel_blocks}/{span3_blocks}, Gate3 .var={g3_var_total}, "
        f"IUPAC d1.vd1 Y/M green, diploid subclonal TP {tp_counts['vj']}->{tp_counts['off']} (GT=0/0 {gt00_counts['vj']}->0)"
    )


# =============================================================================
# Suite 3: Clinical Workflow Presets (--preset wgs|wes|amplicon|ctdna)
# =============================================================================
def run_suite_presets(
    p: Paths,
    cache: RunCache,
    metrics: dict[str, object],
    metrics_lock: threading.Lock,
) -> None:
    print("[Presets] Running 5 clinical workflow preset gates...")
    g4 = str(p.g4_dir)
    w4 = str(p.wave4_dir)

    # Gate 1: --preset wgs (HG002 3-col BED .tbi + HG008 3-col BED .csi + low-MQ Q20 AR-04 + 4-col 3'_UTR)
    hg2_bed3 = os.path.join(p.work_tmp, "p_hg002_3col.bed")
    with open(hg2_bed3, "w", encoding="utf-8") as fh:
        fh.write("chr20\t10000000\t10060000\n")
    hg2_wgs_gz = cache.run_vcf_gz(
        "preset_wgs_hg002", p.giab_ref, os.path.join(g4, "hg002", "HG002.chr20.bam"),
        ["--preset", "wgs", "-N", "HG002", hg2_bed3],
        preset_name="wgs", query_region="chr20:10000001-10060000",
    )
    if not os.path.isfile(hg2_wgs_gz + ".tbi"):
        raise AssertionError("preset_wgs_hg002: missing .tbi sidecar")

    hg8_bed3 = os.path.join(p.work_tmp, "p_hg008_3col.bed")
    with open(os.path.join(g4, "hg008", "test5.bed"), "r", encoding="utf-8") as src, open(hg8_bed3, "w", encoding="utf-8") as dst:
        for ln in src:
            if ln.strip() and not ln.startswith("#"):
                sp = ln.rstrip("\n").split("\t")
                dst.write(f"{sp[0]}\t{sp[1]}\t{sp[2]}\n")
    hg8_wgs_gz = cache.run_vcf_gz(
        "preset_wgs_hg008", p.giab_ref, f"{g4}/hg008/HG008-T.slice.bam|{g4}/hg008/HG008-N.slice.bam",
        ["--preset", "wgs", "-N", "HG008_T|HG008_N", "--fisher", hg8_bed3],
        index_flag="--csi", kind="somatic", preset_name="wgs", query_region="chr1:72271447-72291446",
    )
    with gzip.open(hg8_wgs_gz, "rt", encoding="utf-8") as fh:
        for ln in fh:
            if not ln.startswith("#") and ln.strip():
                flt = ln.split("\t")[6]
                tags = flt.split(";")
                if len(tags) != len(set(tags)) or "Q20;Q20" in flt:
                    raise AssertionError(f"Duplicate FILTER tag in --preset wgs: {ln}")

    # AR-04 stress check: tumor MAPQ=5 (< 10) vs normal MAPQ=60 (> 10) on HG008 chr1:72271447-72276000
    lowmq_reg = "chr1:72271447-72276000"
    lowmq_t_sam = os.path.join(p.work_tmp, "hg008_lowmq_t.sam")
    lowmq_t_bam = os.path.join(p.work_tmp, "hg008_lowmq_t.bam")
    sam_t = subprocess.check_output(
        ["samtools", "view", "-h", f"{g4}/hg008/HG008-T.slice.bam", lowmq_reg],
        text=True,
    )
    with open(lowmq_t_sam, "w", encoding="utf-8") as fh:
        for line in sam_t.splitlines():
            if line.startswith("@"):
                fh.write(line + "\n")
            elif line.strip():
                cols = line.split("\t")
                cols[4] = "5"
                fh.write("\t".join(cols) + "\n")
    subprocess.run(["samtools", "view", "-b", "-o", lowmq_t_bam, lowmq_t_sam], check=True)
    subprocess.run(["samtools", "index", lowmq_t_bam], check=True)
    lowmq_bed3 = os.path.join(p.work_tmp, "hg008_lowmq.bed")
    with open(lowmq_bed3, "w", encoding="utf-8") as fh:
        fh.write("chr1\t72271446\t72276000\n")
    lowmq_gz = cache.run_vcf_gz(
        "preset_wgs_hg008_lowmq",
        p.giab_ref,
        f"{lowmq_t_bam}|{g4}/hg008/HG008-N.slice.bam",
        ["--preset", "wgs", "-N", "TUMOR|NORMAL", lowmq_bed3],
        kind="somatic",
        preset_name="wgs",
        query_region=lowmq_reg,
    )
    q20_seen = 0
    with gzip.open(lowmq_gz, "rt", encoding="utf-8") as fh:
        for ln in fh:
            if not ln.startswith("#") and ln.strip():
                flt = ln.split("\t")[6]
                tags = flt.split(";")
                if "Q20;Q20" in flt or tags.count("Q20") > 1 or len(tags) != len(set(tags)):
                    raise AssertionError(f"AR-04 regression: duplicate Q20;Q20 emitted on low-MAPQ paired SNV: {ln}")
                if "Q20" in tags:
                    q20_seen += 1
    if q20_seen < 1:
        raise AssertionError("Expected >= 1 Q20-filtered record in AR-04 low-MAPQ paired somatic check")

    utr_bed4 = os.path.join(p.work_tmp, "p_utr4.bed")
    with open(utr_bed4, "w", encoding="utf-8") as fh:
        fh.write("chr20\t10000000\t10030000\t3'_UTR\n")
    utr_var = cache.run_text(
        "preset_wgs_utr4", p.giab_ref, os.path.join(g4, "hg002", "HG002.chr20.bam"),
        ["--preset", "wgs", "--emit-intermediate", utr_bed4], audit=False,
    )
    if not utr_var.strip() or any(r.split("\t")[1] != "3'_UTR" for r in utr_var.splitlines() if r.strip()):
        raise AssertionError("4-col 3'_UTR gene preservation under --preset wgs failed")

    # Gate 2: --preset wes + UCSC 8-col CDS BED8 non-hijack
    seqc2_pair = f"{g4}/seqc2/HCC1395.bam|{g4}/seqc2/HCC1395BL.bam"
    wes_bed4 = os.path.join(p.work_tmp, "p_seqc2_wes4.bed")
    ucsc_bed8 = os.path.join(p.work_tmp, "p_seqc2_wes8.bed")
    with open(wes_bed4, "w", encoding="utf-8") as fh:
        fh.write("chr11\t6400000\t6459439\tslice_1\nchr17\t74747665\t74800000\tslice_3\n")
    with open(ucsc_bed8, "w", encoding="utf-8") as fh:
        fh.write("chr11\t6400000\t6459439\tslice_1\t0\t+\t6400500\t6459000\nchr17\t74747665\t74800000\tslice_3\t0\t+\t74748000\t74799500\n")
    wes4_gz = cache.run_vcf_gz("preset_wes4", p.grch38_ref, seqc2_pair, ["--preset", "wes", "-N", "HCC1395|HCC1395BL", "--fisher", wes_bed4], kind="somatic", preset_name="wes")
    wes8_gz = cache.run_vcf_gz("preset_wes8", p.grch38_ref, seqc2_pair, ["--preset", "wes", "-N", "HCC1395|HCC1395BL", "--fisher", ucsc_bed8], kind="somatic", preset_name="wes")
    with gzip.open(wes4_gz, "rt", encoding="utf-8") as f4, gzip.open(wes8_gz, "rt", encoding="utf-8") as f8:
        w8_all = f8.read()
        s4 = {tuple(ln.split("\t")[:5]) for ln in f4 if not ln.startswith("#") and ln.strip()}
        s8 = {tuple(ln.split("\t")[:5]) for ln in w8_all.splitlines() if not ln.startswith("#") and ln.strip()}
    if "AMPBIAS" in w8_all or "GDAMP=" in w8_all or s4 != s8:
        raise AssertionError("UCSC 8-col BED8 hijacked --preset wes")

    # Gate 3: --preset amplicon (single 8-col, paired Illumina PE, QIAseq -sp)
    amp_l861q_gz = cache.run_vcf_gz(
        "preset_amp_l861q", p.fa_hg19, os.path.join(g4, "amplicon", "amplicon.bam"),
        ["--preset", "amplicon", os.path.join(g4, "amplicon", "correct_8columns_amplicon.bed")],
        kind="amplicon", preset_name="amplicon",
    )
    amp_ilm_gz = cache.run_vcf_gz(
        "preset_amp_ilm_tn", p.grch38_ref, f"{w4}/bam/SampleA_ILM1.bam|{w4}/bam/SampleB_ILM1.bam",
        ["--preset", "amplicon", "-N", "SampleA_ILM1|SampleB_ILM1", "--fisher", os.path.join(w4, "beds", "ilm_twoprimer_grch38.amp8.bed")],
        kind="somatic", preset_name="amplicon",
    )
    amp_qgn_gz = cache.run_vcf_gz(
        "preset_amp_qgn_sp", p.grch38_ref, os.path.join(w4, "bam", "SampleA_QGN1.bam"),
        ["--preset", "amplicon", "-sp", "-UN", os.path.join(w4, "beds", "qiaseq_spe_grch38.amp8.bed")],
        kind="amplicon", preset_name="amplicon",
    )
    for agz in (amp_l861q_gz, amp_ilm_gz, amp_qgn_gz):
        with gzip.open(agz, "rt", encoding="utf-8") as fh:
            if not any("GDAMP=" in ln and ";TLAMP=" in ln for ln in fh if not ln.startswith("#")):
                raise AssertionError(f"{agz}: missing GDAMP/TLAMP under --preset amplicon")

    # Gate 4: --preset ctdna ~1,000x 0.5% VAF real-read titration + SampleAC5 d100 override
    bam_a = os.path.join(w4, "bam", "SampleA_QGN1.bam")
    bam_b = os.path.join(w4, "bam", "SampleB_QGN1.bam")
    bam_ac5 = os.path.join(w4, "bam", "SampleAC5_QGN1.bam")
    target_chr, target_pos, win_reg = "chr12", 120999311, "chr12:120999050-120999600"
    wt_sam = os.path.join(p.work_tmp, "ctdna_titration.sam")
    wt_recs = subprocess.check_output(["samtools", "view", "-F", "0xF04", bam_b, win_reg], text=True).splitlines()
    with open(wt_sam, "w", encoding="utf-8") as hfh:
        subprocess.run(["samtools", "view", "-H", bam_b, win_reg], stdout=hfh, check=True)
        hfh.flush()
        for r_i in range(1, 12):
            pref = f"wt_r{r_i}_"
            for rec in wt_recs:
                if rec:
                    hfh.write(pref + rec + "\n")
    sam_a = subprocess.check_output(["samtools", "view", "-F", "0xF04", bam_a, f"{target_chr}:{target_pos}-{target_pos}"], text=True).splitlines()
    fwd_p, rev_p, sfwd, srev = [], [], set(), set()
    for idx, rline in enumerate(sam_a):
        cols = rline.split("\t")
        flag, rpos, mq, seq, qual = int(cols[1]), int(cols[3]), int(cols[4]), cols[9], cols[10]
        if mq < 40:
            continue
        qpos, hit_alt = 0, False
        for n_s, op in CIGAR_RE.findall(cols[5]):
            n = int(n_s)
            if op in ("M", "=", "X"):
                if rpos <= target_pos < rpos + n:
                    qi = qpos + (target_pos - rpos)
                    if seq[qi].upper() == "A" and (ord(qual[qi]) - 33) >= 30 and min(qi, len(seq) - 1 - qi) >= 12:
                        hit_alt = True
                    break
                rpos += n
                qpos += n
            elif op in ("I", "S"):
                qpos += n
            elif op in ("D", "N"):
                rpos += n
        if not hit_alt:
            continue
        astart = int(cols[3])
        cols[0], cols[1], cols[6], cols[7], cols[8] = f"spikeA_{idx}_{cols[0]}", ("16" if (flag & 16) else "0"), "*", "0", "0"
        if (flag & 16) == 0 and astart not in sfwd and len(fwd_p) < 3:
            sfwd.add(astart)
            fwd_p.append("\t".join(cols))
        elif (flag & 16) != 0 and astart not in srev and len(rev_p) < 3:
            srev.add(astart)
            rev_p.append("\t".join(cols))
    with open(wt_sam, "a", encoding="utf-8") as ofh:
        for rec in fwd_p + rev_p:
            ofh.write(rec + "\n")
    titration_bam = os.path.join(p.work_tmp, "ctdna_1000x_0p5pct.bam")
    subprocess.run(["samtools", "sort", "-o", titration_bam, wt_sam], check=True)
    subprocess.run(["samtools", "index", titration_bam], check=True)

    ctdna_bed3 = os.path.join(p.work_tmp, "ctdna_locus.bed")
    with open(ctdna_bed3, "w", encoding="utf-8") as fh:
        fh.write("chr12\t120999118\t120999522\n")
    off_def = cache.run_text("ctdna_off_def", p.grch38_ref, titration_bam, ["--compat=off", "-c", "1", "-S", "2", "-E", "3", "-z", ctdna_bed3])
    if any(int(ln.split("\t")[1]) == target_pos for ln in off_def.splitlines() if not ln.startswith("#") and ln.strip()):
        raise AssertionError("Default --compat=off (-f 0.01) unexpectedly called 0.5% titration locus")
    ctdna_gz = cache.run_vcf_gz("preset_ctdna_1000x", p.grch38_ref, titration_bam, ["--preset", "ctdna", "-N", "CTDNA_0P5", ctdna_bed3], preset_name="ctdna")
    with gzip.open(ctdna_gz, "rt", encoding="utf-8") as fh:
        ct_hits = [ln.rstrip("\n").split("\t") for ln in fh if not ln.startswith("#") and ln.strip() and int(ln.split("\t")[1]) == target_pos and ln.split("\t")[4] == "A"]
    if len(ct_hits) != 1 or ct_hits[0][6] != "PASS" or not ct_hits[0][9].startswith("0/1:"):
        raise AssertionError(f"Expected PASS GT=0/1 at 0.5% ctDNA titration locus, got {ct_hits}")

    ac5_bed3 = os.path.join(p.work_tmp, "ac5_inserts.bed")
    with open(os.path.join(w4, "beds", "qiaseq_spe_grch38.amp8.bed"), "r", encoding="utf-8") as src, open(ac5_bed3, "w", encoding="utf-8") as dst:
        for ln in src:
            if ln.strip() and not ln.startswith("#"):
                sp = ln.rstrip("\n").split("\t")
                dst.write(f"{sp[0]}\t{sp[6]}\t{sp[7]}\n")
    ac5_def = cache.run_text("preset_ctdna_ac5_d100", p.grch38_ref, bam_ac5, ["--preset", "ctdna", ac5_bed3])
    ac5_ov = cache.run_text("preset_ctdna_ac5_d20", p.grch38_ref, bam_ac5, ["--preset", "ctdna", "--vcf-min-depth", "20", ac5_bed3])
    d100_only = {tuple(ln.split("\t")[:5]) for ln in ac5_def.splitlines() if not ln.startswith("#") and ln.strip() and ln.split("\t")[6] == "d100"}
    ov_pass = {tuple(ln.split("\t")[:5]) for ln in ac5_ov.splitlines() if not ln.startswith("#") and ln.strip() and ln.split("\t")[6] == "PASS"}
    if len(d100_only) < 10 or (d100_only - ov_pass):
        raise AssertionError("--preset ctdna d100 -> --vcf-min-depth 20 promotion check failed")

    # Gate 5: -R 1-based coordinate invariance & 7 negative CLI guards
    reg_str = "chr20:10000001-10035000"
    hg2_bam = os.path.join(g4, "hg002", "HG002.chr20.bam")
    pr_vcf = cache.run_text("preset_wgs_R", p.giab_ref, hg2_bam, ["--preset", "wgs", "-f", "0.01", "-R", reg_str])
    co_vcf = cache.run_text("compat_off_R", p.giab_ref, hg2_bam, ["--compat=off", "-f", "0.01", "-R", reg_str])
    p_coords = [tuple(ln.split("\t")[:5]) for ln in pr_vcf.splitlines() if ln and not ln.startswith("#")]
    c_coords = [tuple(ln.split("\t")[:5]) for ln in co_vcf.splitlines() if ln and not ln.startswith("#")]
    if len(p_coords) < 15 or p_coords != c_coords or "##FILTER=<ID=f0.01," not in pr_vcf:
        raise AssertionError("-R coordinate invariance / -f override under --preset wgs failed")

    amp_bam = os.path.join(g4, "amplicon", "amplicon.bam")
    amp_bed = os.path.join(g4, "amplicon", "correct_8columns_amplicon.bed")
    neg_guards = [
        ([p.vdr_bin, "--preset", "wgs", "--parity", "-G", p.giab_ref, "-b", hg2_bam, "-R", "chr20:10000001-10000100"], ("--preset", "--parity")),
        ([p.vdr_bin, "--preset", "wgs", "--compat=vardictjava", "-G", p.giab_ref, "-b", hg2_bam, "-R", "chr20:10000001-10000100"], ("--preset", "--compat=vardictjava")),
        ([p.vdr_bin, "--preset", "amplicon", "-G", p.fa_hg19, "-b", amp_bam, "-R", "chr1:100-200"], ("--preset amplicon", "-R")),
        ([p.vdr_bin, "--preset", "amplicon", "-i", "-G", p.fa_hg19, "-b", amp_bam, amp_bed], ("--preset amplicon", "-i")),
        ([p.vdr_bin, "--preset", "hifi", "-G", p.giab_ref, "-b", hg2_bam, "-R", "chr20:10000001-10000100"], ("not implemented", "#96")),
        ([p.vdr_bin, "--preset", "unknown", "-G", p.giab_ref, "-b", hg2_bam, "-R", "chr20:10000001-10000100"], ("Invalid --preset", "wgs")),
        ([p.vdr_bin, "-preset", "wgs", "-G", p.giab_ref, "-b", hg2_bam, "-R", "chr20:10000001-10000100"], ("--preset",)),
    ]
    for cmd, subs in neg_guards:
        res = subprocess.run(cmd, capture_output=True, text=True, check=False)
        comb = (res.stderr or "") + "\n" + (res.stdout or "")
        if res.returncode == 0 or "panicked at" in comb or any(s not in comb for s in subs):
            raise AssertionError(f"Negative preset guard failed: {cmd}")

    with metrics_lock:
        metrics["presets"] = {
            "wes_sites": len(s4),
            "ctdna_d100_promoted": len(d100_only),
            "r_coord_sites": len(p_coords),
            "ar04_q20_seen": q20_seen,
        }
    print(
        f"  [PASS Presets] WGS(AR-04 Q20={q20_seen})/WES({len(s4)})/Amplicon/ctDNA(0.5% PASS, +{len(d100_only)} d100->d20)/7 negative guards green"
    )


# =============================================================================
# Suite 4: Wave-6 & Expanded 12-Modality Real-Data Truth & Invariant Verification
# =============================================================================
def run_suite_wave6(
    p: Paths,
    cache: RunCache,
    rtg_eval: RtgEvaluator,
    metrics: dict[str, object],
    metrics_lock: threading.Lock,
    *,
    include_exome: bool = False,
) -> None:
    print("[Wave-6 & 12-Modality] Running 38-switch invariants & multi-cohort truth evaluations (0 rtg format)...")
    g4 = str(p.g4_dir)
    w4 = str(p.wave4_dir)

    # Gates 1, 2, 4: HG002 1.5 Mbp + SEQC2 slices + HG008 slices + CMRG + Wave-4 QIAseq/Illumina
    hg2_15m_bed = os.path.join(p.work_tmp, "hg002_chr20_1_5mb.bed")
    with open(hg2_15m_bed, "w", encoding="utf-8") as fh:
        fh.write("chr20\t10000000\t11500000\thg002_c20\n")

    w6_cohorts = [
        ("GIAB_HG002_1_5M", "simple", p.giab_ref, os.path.join(g4, "hg002", "HG002.chr20.bam"),
         ["-N", "HG002", "-c", "1", "-S", "2", "-E", "3", "-g", "4", "-z", "-f", "0.01", hg2_15m_bed]),
        ("SEQC2_HCC1395_TN", "somatic", p.grch38_ref, f"{g4}/seqc2/HCC1395.bam|{g4}/seqc2/HCC1395BL.bam",
         ["-N", "TUMOR|NORMAL", "-c", "1", "-S", "2", "-E", "3", "-g", "4", "-z", "-f", "0.01", "--fisher", os.path.join(g4, "seqc2", "slices.bed")]),
        ("GIAB_HG008_TN", "somatic", p.giab_ref, f"{g4}/hg008/HG008-T.slice.bam|{g4}/hg008/HG008-N.slice.bam",
         ["-N", "HG008T|HG008N", "-c", "1", "-S", "2", "-E", "3", "-z", "-f", "0.01", "--fisher", os.path.join(g4, "hg008", "slices.bed")]),
        ("Wave4_HG002_CMRG", "simple", p.giab_ref, os.path.join(w4, "hg002_cmrg_mosaic", "HG002.cmrg_mosaic.slice.bam"),
         ["-N", "HG002", "-c", "1", "-S", "2", "-E", "3", "-g", "4", "-z", "-f", "0.01", "--fisher", os.path.join(w4, "hg002_cmrg_mosaic", "hg002_cmrg_mosaic.bed")]),
        ("Wave4_QIAseq_ctDNA_TN", "somatic", p.grch38_ref, f"{w4}/bam/SampleAC5_QGN1.bam|{w4}/bam/SampleB_QGN1.bam",
         ["-N", "SampleAC5_QGN1|SampleB_QGN1", "-a", "10:0.95", "-sp", "-z", "-f", "0.0025", "--fisher", os.path.join(w4, "beds", "qiaseq_spe_grch38.amp8.bed")]),
        ("Wave4_Illumina_Amp_TN", "somatic", p.grch38_ref, f"{w4}/bam/SampleA_ILM1.bam|{w4}/bam/SampleB_ILM1.bam",
         ["-N", "SampleA_ILM1|SampleB_ILM1", "-a", "10:0.95", "-z", "-f", "0.01", "--fisher", os.path.join(w4, "beds", "ilm_twoprimer_grch38.amp8.bed")]),
    ]

    vj_subst_shift3 = 0
    off_indel_shift3 = 0
    multi_end_dups = 0
    for cname, ckind, rfa, bspec, cextra in w6_cohorts:
        cgz = cache.run_vcf_gz(f"w6_{cname}", rfa, bspec, ["--compat=off", "--vcf-all-alleles", *cextra], kind=ckind)
        is_paired_amp = ("|" in bspec) and ("-a" in cextra)
        if not is_paired_amp:
            vj_txt = cache.run_text(f"w6_vj_{cname}", rfa, bspec, ["--compat=vardictjava", "--vcf-all-alleles", *cextra], kind=ckind, audit=False)
            for ln in vj_txt.splitlines():
                if not ln.startswith("#") and ln.strip():
                    sp = ln.split("\t")
                    if not sp[4].startswith("<") and sp[4] != "." and "," not in sp[4] and len(sp[3]) == len(sp[4]):
                        kv = dict(x.split("=", 1) for x in sp[7].split(";") if "=" in x)
                        if int(kv.get("SHIFT3", "0")) > 0:
                            vj_subst_shift3 += 1
        dup_ends: dict[tuple[str, int, str], set[int]] = {}
        with gzip.open(cgz, "rt", encoding="utf-8") as fh:
            for ln in fh:
                if not ln.startswith("#") and ln.strip():
                    sp = ln.rstrip("\n").split("\t")
                    kv = dict(x.split("=", 1) for x in sp[7].split(";") if "=" in x)
                    if cname == "GIAB_HG002_1_5M" and sp[4] == "<DUP>":
                        dup_ends.setdefault((sp[0], int(sp[1]), sp[3]), set()).add(int(kv["END"]))
                    for sa in sp[4].split(","):
                        if not sa.startswith("<") and sa != "." and len(sp[3]) != len(sa) and int(kv.get("SHIFT3", "0")) > 0:
                            off_indel_shift3 += 1
        if cname == "GIAB_HG002_1_5M":
            multi_end_dups = sum(1 for ends in dup_ends.values() if len(ends) > 1)
    if vj_subst_shift3 <= 0 or off_indel_shift3 <= 0 or multi_end_dups < 20:
        raise AssertionError(
            f"Wave-6 Gate 1/2 failed: vj_subst_shift3={vj_subst_shift3}, off_indel_shift3={off_indel_shift3}, multi_end_dups={multi_end_dups}"
        )

    # Gate 3: Switch 34 (-u --compat=off mate-overlap bidirectional ALD balance)
    u_cases = [
        ("SEQC2_pe_u", p.grch38_ref, os.path.join(g4, "seqc2", "HCC1395.bam"), ["-R", "chr17:74747665-74952140", "-u", "-f", "0.01"]),
        ("HG002_pe_u", p.giab_ref, os.path.join(g4, "hg002", "HG002.chr20.bam"), ["-R", "chr20:10000000-10080000", "-u", "-f", "0.01"]),
        ("Illumina_pe_u", p.grch38_ref, os.path.join(w4, "bam", "SampleA_ILM1.bam"), ["-a", "10:0.95", "-z", "-u", "-f", "0.01", os.path.join(w4, "beds", "ilm_twoprimer_grch38.amp8.bed")]),
    ]
    vj_one_strand, off_one_strand, rebalanced = 0, 0, 0
    for ulabel, rfa, ubam, uargs in u_cases:
        uvj = cache.run_text(f"w6_{ulabel}_vj", rfa, ubam, ["--compat=vardictjava", *uargs], audit=False)
        uoff = cache.run_text(f"w6_{ulabel}_off", rfa, ubam, ["--compat=off", *uargs], kind=("amplicon" if "-a" in uargs else "simple"))
        vj_m = {}
        for ln in uvj.splitlines():
            if not ln.startswith("#") and ln.strip():
                sp = ln.split("\t")
                vals = dict(zip(sp[8].split(":"), sp[9].split(":")))
                ald = tuple(int(x) for x in vals["ALD"].split(","))
                if int(vals["VD"]) >= 4 and (ald[0] == 0 or ald[1] == 0):
                    vj_one_strand += 1
                vj_m[(sp[0], int(sp[1]), sp[3], sp[4])] = ald
        for ln in uoff.splitlines():
            if not ln.startswith("#") and ln.strip():
                sp = ln.split("\t")
                vals = dict(zip(sp[8].split(":"), sp[9].split(":")))
                ald = tuple(int(x) for x in vals["ALD"].split(","))
                if int(vals["VD"]) >= 4 and (ald[0] == 0 or ald[1] == 0):
                    off_one_strand += 1
                k = (sp[0], int(sp[1]), sp[3], sp[4])
                if k in vj_m and vj_m[k] != ald and ald[0] > 0 and ald[1] > 0:
                    rebalanced += 1
    if off_one_strand > vj_one_strand or rebalanced <= 0:
        raise AssertionError(f"Switch 34 -u rebalancing failed: vj={vj_one_strand}, off={off_one_strand}, rebalanced={rebalanced}")

    # Expanded 12-Modality Real-Data Coverage & Truth Concordance (0 rtg format calls!)
    # Modality 1b: HG002.GRCh38.2x250.bam multi-chrom slices (chr1, chr6 MHC/HLA, chr13/chr21 acrocentric --csi, chr17 BRCA1 + NIST v4.2.1 diploid GT)
    hg2_250_bam = os.path.join(g4, "hg002", "HG002.GRCh38.2x250.bam")
    hg2_multichr_bed = os.path.join(p.work_tmp, "hg002_2x250_multichr.bed")
    hg2_euchr_bed = os.path.join(p.work_tmp, "hg002_2x250_euchr.bed")
    hg2_mhc_bed = os.path.join(p.work_tmp, "hg002_2x250_mhc.bed")
    with open(hg2_multichr_bed, "w", encoding="utf-8") as fh, \
         open(hg2_euchr_bed, "w", encoding="utf-8") as fe, \
         open(hg2_mhc_bed, "w", encoding="utf-8") as fm:
        for row in (
            "chr1\t10000000\t10500000\tchr1_eu\n",
            "chr13\t19000000\t19500000\tchr13_acro\n",
            "chr17\t43044295\t43500000\tchr17_brca1\n",
            "chr21\t41800000\t42300000\tchr21_acro\n",
        ):
            fh.write(row)
            fe.write(row)
        mhc_row = "chr6\t29900000\t29950000\tchr6_mhc\n"
        fh.write(mhc_row)
        fm.write(mhc_row)
    hg2_250_gz = cache.run_vcf_gz(
        "mod1_hg002_2x250_multichr", p.giab_ref, hg2_250_bam,
        ["-N", "HG002", "-c", "1", "-S", "2", "-E", "3", "-g", "4", "-z", "-f", "0.05", "--pass-only", "--compat=off", "--vcf-all-alleles", hg2_multichr_bed],
        index_flag="--csi", query_region="chr13:19000001-19500000",
    )
    hg2_truth_bed_all = os.path.join(g4, "hg002", "truth", "HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed")
    hg2_250_eval_bed = os.path.join(p.work_tmp, "hg002_2x250_eval.bed")
    hg2_250_mhc_eval_bed = os.path.join(p.work_tmp, "hg002_2x250_mhc_eval.bed")
    intersect_beds(hg2_euchr_bed, hg2_truth_bed_all, hg2_250_eval_bed)
    intersect_beds(hg2_mhc_bed, hg2_truth_bed_all, hg2_250_mhc_eval_bed)
    hg2_250_all_eval = os.path.join(p.work_tmp, "hg002_2x250_all_eval.bed")
    intersect_beds(hg2_multichr_bed, hg2_truth_bed_all, hg2_250_all_eval)
    hg2_250_truth = rtg_eval.prepare_germline_truth_gz(
        "hg002_2x250", p.giab_ref, os.path.join(g4, "hg002", "truth", "HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"), hg2_250_all_eval
    )
    hg2_250_gt, _, _, _ = rtg_eval.run_vcfeval(
        "hg002_2x250_euchr_gt", hg2_250_truth, hg2_250_gz, hg2_250_eval_bed, p.sdf_giab, squash_ploidy=False
    )
    hg2_250_mhc_gt, _, _, _ = rtg_eval.run_vcfeval(
        "hg002_2x250_mhc_gt", hg2_250_truth, hg2_250_gz, hg2_250_mhc_eval_bed, p.sdf_giab, squash_ploidy=False
    )
    if hg2_250_gt[3] < 0.92:
        raise AssertionError(f"HG002 2x250 multi-chrom euchromatic/acrocentric diploid GT F1 too low: {hg2_250_gt}")
    mhc_prec = hg2_250_mhc_gt[0] / max(1, hg2_250_mhc_gt[0] + hg2_250_mhc_gt[1])
    if hg2_250_mhc_gt[0] < 500 or hg2_250_mhc_gt[3] < 0.65 or mhc_prec < 0.88:
        raise AssertionError(f"HG002 2x250 chr6 HLA-A MHC diploid GT metrics regressed: {hg2_250_mhc_gt} (prec={mhc_prec:.4f})")

    # Modality 1c: HG002 CMRG v1.00 (diploid GT, never --squash-ploidy) & MosaicSNVv1.1 (--squash-ploidy)
    cmrg_bam = os.path.join(w4, "hg002_cmrg_mosaic", "HG002.cmrg_mosaic.slice.bam")
    cmrg_bed = os.path.join(w4, "hg002_cmrg_mosaic", "hg002_cmrg_mosaic.bed")
    cmrg_sub_bed = os.path.join(p.work_tmp, "hg002_cmrg_sub.bed")
    mosaic_sub_bed = os.path.join(p.work_tmp, "hg002_mosaic_sub.bed")
    with open(cmrg_bed, "r", encoding="utf-8") as fh, \
         open(cmrg_sub_bed, "w", encoding="utf-8") as fc, \
         open(mosaic_sub_bed, "w", encoding="utf-8") as fm:
        for ln in fh:
            if "CMRG_" in ln:
                fc.write(ln)
            elif "MOSAIC_" in ln:
                fm.write(ln)
    cmrg_eval_bed = os.path.join(p.work_tmp, "hg002_cmrg_eval.bed")
    mosaic_eval_bed = os.path.join(p.work_tmp, "hg002_mosaic_eval.bed")
    intersect_beds(cmrg_sub_bed, os.path.join(w4, "hg002_cmrg_mosaic", "HG002_GRCh38_CMRG_smallvar_v1.00.bed"), cmrg_eval_bed)
    intersect_beds(mosaic_sub_bed, os.path.join(w4, "hg002_cmrg_mosaic", "HG002_GRCh38_MosaicSNVv1.1.bed"), mosaic_eval_bed)

    cmrg_germ_gz = cache.run_vcf_gz(
        "mod1_hg002_cmrg_germ_calls", p.giab_ref, cmrg_bam,
        ["-N", "HG002", "-c", "1", "-S", "2", "-E", "3", "-g", "4", "-z", "-f", "0.15", "--pass-only", "--compat=off", "--vcf-all-alleles", cmrg_sub_bed],
    )
    cmrg_truth_gz = rtg_eval.prepare_germline_truth_gz(
        "hg002_cmrg", p.giab_ref, os.path.join(w4, "hg002_cmrg_mosaic", "HG002_GRCh38_CMRG_smallvar_v1.00.vcf.gz"), cmrg_eval_bed
    )
    cmrg_gt, _, _, _ = rtg_eval.run_vcfeval(
        "hg002_cmrg_gt", cmrg_truth_gz, cmrg_germ_gz, cmrg_eval_bed, p.sdf_giab, squash_ploidy=False
    )
    if cmrg_gt[0] < 2 or cmrg_gt[1] != 0:
        raise AssertionError(f"HG002 CMRG diploid GT expected >= 2 TPs and 0 FPs, got {cmrg_gt}")

    mosaic_calls_gz = cache.run_vcf_gz(
        "mod1_hg002_mosaic_calls", p.giab_ref, cmrg_bam,
        ["-N", "HG002", "-c", "1", "-S", "2", "-E", "3", "-g", "4", "-z", "-f", "0.01", "--pass-only", "--compat=off", "--vcf-all-alleles", mosaic_sub_bed],
    )
    mosaic_truth_gz = rtg_eval.prepare_somatic_truth_gz(
        "hg002_mosaic", p.giab_ref, os.path.join(w4, "hg002_cmrg_mosaic", "HG002_GRCh38_MosaicSNVv1.1.vcf.gz"), mosaic_eval_bed
    )
    mosaic_sq, _, _, _ = rtg_eval.run_vcfeval(
        "hg002_mosaic_sq", mosaic_truth_gz, mosaic_calls_gz, mosaic_eval_bed, p.sdf_giab, squash_ploidy=True, sample_arg="TUMOR,HG002"
    )
    if mosaic_sq[0] < 8:
        raise AssertionError(f"HG002 MosaicSNVv1.1 --squash-ploidy expected 8/8 TPs, got {mosaic_sq}")

    # Modality 1a: HG002 chr20 1.5 Mbp (diploid GT + SNP/INDEL ROC + -u, NEVER --squash-ploidy!)
    hg2_15m_eval = os.path.join(p.work_tmp, "hg002_15m_eval.bed")
    intersect_beds(hg2_15m_bed, os.path.join(g4, "hg002", "eval.bed"), hg2_15m_eval)
    hg2_15m_truth = rtg_eval.prepare_germline_truth_gz(
        "hg002_15m", p.giab_ref, os.path.join(g4, "hg002", "truth", "truth.vcf.gz"), hg2_15m_bed
    )
    hg2_off_gz = cache.run_vcf_gz(
        "w6_hg002_15m_off_pass", p.giab_ref, os.path.join(g4, "hg002", "HG002.chr20.bam"),
        ["-N", "HG002", "-c", "1", "-S", "2", "-E", "3", "-z", "-f", "0.05", "--pass-only", "--compat=off", "--vcf-all-alleles", hg2_15m_bed],
    )
    hg2_u_off_gz = cache.run_vcf_gz(
        "w6_hg002_15m_u_off_pass", p.giab_ref, os.path.join(g4, "hg002", "HG002.chr20.bam"),
        ["-N", "HG002", "-c", "1", "-S", "2", "-E", "3", "-z", "-f", "0.05", "-u", "--pass-only", "--compat=off", "--vcf-all-alleles", hg2_15m_bed],
    )
    g2_off_all, g2_off_snp, g2_off_ind, _ = rtg_eval.run_vcfeval("hg2_off_gt", hg2_15m_truth, hg2_off_gz, hg2_15m_eval, p.sdf_giab, squash_ploidy=False)
    g2u_off_all, g2u_off_snp, g2u_off_ind, _ = rtg_eval.run_vcfeval("hg2_u_off_gt", hg2_15m_truth, hg2_u_off_gz, hg2_15m_eval, p.sdf_giab, squash_ploidy=False)
    if (
        g2_off_snp is None
        or g2_off_ind is None
        or g2_off_all[3] < 0.9590
        or g2_off_snp[3] < 0.9800
        or g2_off_ind[3] < 0.8250
    ):
        raise AssertionError(f"GIAB HG002 1.5M diploid GT rtg vcfeval below threshold: all={g2_off_all}, snp={g2_off_snp}, indel={g2_off_ind}")
    if g2u_off_all[3] < 0.9550 or (g2u_off_snp and g2u_off_snp[3] < 0.9780) or (g2u_off_ind and g2u_off_ind[3] < 0.8200):
        raise AssertionError(f"GIAB HG002 1.5M -u diploid GT rtg vcfeval below threshold: all={g2u_off_all}, snp={g2u_off_snp}, indel={g2u_off_ind}")

    # Modality 2: HG005 + Platinum NA12878 (30 kbp PacBio HiFi =/X CIGARs) + NA24385 b37 (ALL diploid GT, NEVER --squash-ploidy!)
    tc_dir = os.path.join(g4, "truth_cache")
    diploid_cohorts = [
        ("HG005", p.giab_ref, p.sdf_giab, os.path.join(g4, "multimodal", "hg005_s4a.bam"),
         os.path.join(tc_dir, "hg005_eval.bed"), os.path.join(tc_dir, "hg005_truth.vcf.gz"), "HG005", ["-f", "0.05"], 0.90),
        ("Platinum_NA12878", p.grch38_ref, p.sdf_d1vd1, os.path.join(g4, "multimodal", "platinum_s2a.bam"),
         os.path.join(tc_dir, "platinum_eval.bed"), os.path.join(tc_dir, "platinum_truth.vcf.gz"), "HG001",
         ["-x", "35000", "-m", "1000", "--vcf-max-nm", "1000", "-f", "0.15"], 0.94),
        ("NA24385_b37", p.fa_b37, p.sdf_b37, os.path.join(p.vj_bams, "NA24385_chr_20_36018823_36039513.bam"),
         os.path.join(tc_dir, "na24385_eval.bed"), os.path.join(tc_dir, "na24385_b37_truth.vcf.gz"), "HG002", ["-f", "0.30"], 0.95),
    ]
    mod2_f1s = {}
    for dname, rfa, dsdf, dbam, dbed, dtruth_src, dsmp, d_extra, min_f1 in diploid_cohorts:
        dcalls = cache.run_vcf_gz(
            f"mod2_{dname}", rfa, dbam,
            ["-N", dsmp, "-c", "1", "-S", "2", "-E", "3", "-z", *d_extra, "--pass-only", "--compat=off", "--vcf-all-alleles", dbed],
        )
        dtruth = rtg_eval.prepare_germline_truth_gz(f"mod2_{dname}", rfa, dtruth_src, dbed)
        dm, _, _, _ = rtg_eval.run_vcfeval(f"mod2_{dname}", dtruth, dcalls, dbed, dsdf, squash_ploidy=False)
        if dm[3] < min_f1:
            raise AssertionError(f"{dname} diploid GT F1 below {min_f1}: {dm}")
        mod2_f1s[dname] = dm[3]

    # Modality 3: SEQC2 HCC1395/HCC1395BL (--squash-ploidy sSNV + sINDEL, standard & -u, + optional exome.engine.bed)
    def eval_somatic_cohort(
        tag: str,
        rfa: str,
        sdf_d: str,
        bpair: str,
        call_bed: str,
        eval_bed: str,
        truth_gz: str,
        cflag: str,
        extra: list[str],
        *,
        freq_s: str = "0.01",
    ) -> tuple[tuple[int, int, int, float, float], tuple[int, int, int, float] | None, tuple[int, int, int, float] | None]:
        rtxt = cache.run_text(
            f"som_{tag}_{cflag.split('=')[1]}", rfa, bpair,
            ["-N", "TUMOR|NORMAL", "-c", "1", "-S", "2", "-E", "3", "-z", "-f", freq_s, "--fisher", cflag, "--vcf-all-alleles", *extra, call_bed],
            kind="somatic", audit=(cflag == "--compat=off"),
        )
        som_vcf = os.path.join(p.work_tmp, f"som_{tag}_{cflag.split('=')[1]}.vcf")
        with open(som_vcf, "w", encoding="utf-8") as out:
            for ln in rtxt.splitlines():
                if ln.startswith("#"):
                    out.write(ln + "\n")
                else:
                    sp = ln.split("\t")
                    if sp[6] == "PASS" and ("STATUS=StrongSomatic" in sp[7] or "STATUS=LikelySomatic" in sp[7]):
                        out.write(ln + "\n")
        subprocess.run(["bgzip", "-f", som_vcf], check=True)
        subprocess.run(["tabix", "-f", "-p", "vcf", som_vcf + ".gz"], check=True)
        m_all, m_snp, m_ind, _ = rtg_eval.run_vcfeval(
            f"{tag}_{cflag.split('=')[1]}", truth_gz, som_vcf + ".gz", eval_bed, sdf_d,
            squash_ploidy=True, sample_arg="TUMOR,TUMOR",
        )
        return m_all, m_snp, m_ind

    seqc2_slices = os.path.join(g4, "seqc2", "slices.bed")
    seqc2_eval_bed = os.path.join(p.work_tmp, "w6_seqc2_eval.bed")
    intersect_beds(seqc2_slices, os.path.join(g4, "seqc2", "eval.bed"), seqc2_eval_bed)
    seqc2_truth_gz = rtg_eval.prepare_somatic_truth_gz(
        "w6_seqc2", p.grch38_ref, os.path.join(g4, "seqc2", "truth", "truth.vcf.gz"), seqc2_slices
    )
    seqc2_pair = f"{g4}/seqc2/HCC1395.bam|{g4}/seqc2/HCC1395BL.bam"
    s_vj, s_vj_snp, s_vj_ind = eval_somatic_cohort("seqc2", p.grch38_ref, p.sdf_d1vd1, seqc2_pair, seqc2_slices, seqc2_eval_bed, seqc2_truth_gz, "--compat=vardictjava", ["-g", "4"])
    s_off, s_off_snp, s_off_ind = eval_somatic_cohort("seqc2", p.grch38_ref, p.sdf_d1vd1, seqc2_pair, seqc2_slices, seqc2_eval_bed, seqc2_truth_gz, "--compat=off", ["-g", "4"])
    if s_off[3] < s_vj[3] or s_off[0] < s_vj[0] or s_off[1] > s_vj[1]:
        raise AssertionError(f"SEQC2 rtg vcfeval regression: vj={s_vj} vs off={s_off}")
    if (
        s_vj_snp is None
        or s_off_snp is None
        or s_off_snp[3] < s_vj_snp[3]
        or s_off_snp[0] < s_vj_snp[0]
        or s_off_snp[1] > s_vj_snp[1]
    ):
        raise AssertionError(f"SEQC2 sSNV rtg vcfeval regression: vj={s_vj_snp} vs off={s_off_snp}")
    if (
        s_vj_ind is None
        or s_off_ind is None
        or s_off_ind[3] < s_vj_ind[3]
        or s_off_ind[0] < s_vj_ind[0]
        or s_off_ind[1] > s_vj_ind[1]
    ):
        raise AssertionError(f"SEQC2 sINDEL rtg vcfeval regression: vj={s_vj_ind} vs off={s_off_ind}")
    su_vj, _, _ = eval_somatic_cohort("seqc2_u", p.grch38_ref, p.sdf_d1vd1, seqc2_pair, seqc2_slices, seqc2_eval_bed, seqc2_truth_gz, "--compat=vardictjava", ["-u"])
    su_off, _, _ = eval_somatic_cohort("seqc2_u", p.grch38_ref, p.sdf_d1vd1, seqc2_pair, seqc2_slices, seqc2_eval_bed, seqc2_truth_gz, "--compat=off", ["-u"])
    if su_off[0] < su_vj[0] or su_off[4] + 1e-6 < su_vj[4]:
        raise AssertionError(f"SEQC2 -u rtg vcfeval regression: vj={su_vj} vs off={su_off}")

    if include_exome:
        ex_bed = os.path.join(g4, "seqc2", "exome.engine.bed")
        ex_eval_bed = os.path.join(p.work_tmp, "seqc2_exome_eval.bed")
        intersect_beds(ex_bed, os.path.join(g4, "seqc2", "eval.bed"), ex_eval_bed)
        ex_truth_gz = rtg_eval.prepare_somatic_truth_gz(
            "seqc2_exome", p.grch38_ref, os.path.join(g4, "seqc2", "truth", "truth.vcf.gz"), ex_eval_bed
        )
        ex_gz = cache.run_vcf_gz(
            "seqc2_full_exome", p.grch38_ref, seqc2_pair,
            ["-N", "TUMOR|NORMAL", "-c", "1", "-S", "2", "-E", "3", "-z", "-f", "0.05", "--fisher", "--compat=off", "--vcf-all-alleles", ex_bed],
            kind="somatic",
        )
        ex_som_vcf = os.path.join(p.work_tmp, "seqc2_exome_som_off.vcf")
        with gzip.open(ex_gz, "rt", encoding="utf-8") as fh, open(ex_som_vcf, "w", encoding="utf-8") as out:
            for ln in fh:
                if ln.startswith("#"):
                    out.write(ln)
                elif ln.strip():
                    sp = ln.split("\t")
                    if sp[6] == "PASS" and ("STATUS=StrongSomatic" in sp[7] or "STATUS=LikelySomatic" in sp[7]):
                        out.write(ln)
        subprocess.run(["bgzip", "-f", ex_som_vcf], check=True)
        subprocess.run(["tabix", "-f", "-p", "vcf", ex_som_vcf + ".gz"], check=True)
        ex_all, ex_snp, ex_ind, _ = rtg_eval.run_vcfeval(
            "seqc2_exome_off", ex_truth_gz, ex_som_vcf + ".gz", ex_eval_bed, p.sdf_d1vd1,
            squash_ploidy=True, sample_arg="TUMOR,TUMOR",
        )
        if ex_snp is None or ex_ind is None or ex_all[0] < 800 or ex_snp[0] < 700 or ex_ind[0] < 50:
            raise AssertionError(
                f"SEQC2 full exome --squash-ploidy rtg vcfeval regression: all={ex_all}, sSNV={ex_snp}, sINDEL={ex_ind}"
            )
        with metrics_lock:
            metrics["seqc2_exome"] = {
                "vcf_gz": ex_gz,
                "f1": ex_all[3],
                "ssnv_f1": ex_snp[3],
                "sindel_f1": ex_ind[3],
                "tp": ex_all[0],
            }

    # Modality 4: GIAB HG008 (BAM == CRAM under --compat=off, ladder p2T/p13T, v0.3 & v0.4 --squash-ploidy, v0.6 SV coordinate/END overlap auditor)
    hg8_bam_pair = f"{g4}/hg008/HG008-T.slice.bam|{g4}/hg008/HG008-N.slice.bam"
    hg8_cram_pair = f"{g4}/hg008/HG008-T.slice.cram|{g4}/hg008/HG008-N.slice.cram"
    hg8_slices = os.path.join(g4, "hg008", "slices.bed")
    h8_bam_gz = cache.run_vcf_gz("w6_GIAB_HG008_TN", p.giab_ref, hg8_bam_pair, ["--compat=off", "--vcf-all-alleles", "-N", "HG008T|HG008N", "-c", "1", "-S", "2", "-E", "3", "-z", "-f", "0.01", "--fisher", hg8_slices], kind="somatic")
    h8_cram_gz = cache.run_vcf_gz("mod4_hg008_cram", p.giab_ref, hg8_cram_pair, ["--compat=off", "--vcf-all-alleles", "-N", "HG008T|HG008N", "-c", "1", "-S", "2", "-E", "3", "-z", "-f", "0.01", "--fisher", hg8_slices], index_flag="--csi", kind="somatic")
    with gzip.open(h8_bam_gz, "rt", encoding="utf-8") as fb, gzip.open(h8_cram_gz, "rt", encoding="utf-8") as fc:
        h8_b_rows = [ln.split("\t")[:8] for ln in fb if not ln.startswith("#") and ln.strip()]
        h8_c_rows = [ln.split("\t")[:8] for ln in fc if not ln.startswith("#") and ln.strip()]
    if not h8_b_rows or h8_b_rows != h8_c_rows:
        raise AssertionError("GIAB HG008 slices.bed BAM != CRAM under --compat=off")

    for lad_name in ("p2T", "p13T"):
        lad_bam = os.path.join(g4, "hg008", "ladder", f"{lad_name}.chr1.bam")
        lad_bed = os.path.join(g4, "hg008", "ladder", "chr1.bed")
        cache.run_vcf_gz(f"mod4_hg008_ladder_{lad_name}", p.giab_ref, lad_bam, ["-N", lad_name, "-c", "1", "-S", "2", "-E", "3", "-g", "4", "-z", "-f", "0.01", "--compat=off", "--vcf-all-alleles", lad_bed])

    # HG008 v0.3 & v0.4 --squash-ploidy rtg vcfeval
    for ver_tag, t_fn, b_fn in (
        ("v03", "HG008-T_somatic_smvar_benchmark_v0.3_tumorvariants.vcf.gz", "HG008-T_somatic_smvar_benchmark_v0.3_all.bed"),
        ("v04", "HG008-T_somatic_smvar_benchmark_v0.4_GRCh38_tumorvariants.vcf.gz", "HG008-T_somatic_smvar_benchmark_v0.4_GRCh38_all.bed"),
    ):
        ev_b = os.path.join(p.work_tmp, f"hg008_{ver_tag}_eval.bed")
        intersect_beds(hg8_slices, os.path.join(g4, "hg008", "truth", b_fn), ev_b, max_a_lines=40)
        tr_gz = rtg_eval.prepare_somatic_truth_gz(f"hg008_{ver_tag}", p.giab_ref, os.path.join(g4, "hg008", "truth", t_fn), ev_b)
        hvj, _, _ = eval_somatic_cohort(f"hg008_{ver_tag}", p.giab_ref, p.sdf_giab, hg8_bam_pair, ev_b, ev_b, tr_gz, "--compat=vardictjava", [])
        hoff, _, _ = eval_somatic_cohort(f"hg008_{ver_tag}", p.giab_ref, p.sdf_giab, hg8_bam_pair, ev_b, ev_b, tr_gz, "--compat=off", [])
        if hoff[3] < hvj[3] or hoff[0] < hvj[0] or hoff[1] > hvj[1]:
            raise AssertionError(f"GIAB HG008 {ver_tag} regression: vj={hvj} vs off={hoff}")
        with metrics_lock:
            metrics[f"hg008_{ver_tag}_f1"] = hoff[3]

    # HG008 v0.6 SV coordinate/END breakpoint-overlap auditor (contained DELs in slices.bed + somatic DUP/INV)
    v06_vcf = os.path.join(g4, "hg008", "truth", "GRCh38_HG008-T-V0.6_somatic-stvar_PASS.draftbenchmark.vcf.gz")
    slice_windows: dict[str, list[tuple[int, int]]] = {}
    with open(hg8_slices, "r", encoding="utf-8") as fh:
        for ln in fh:
            if ln.strip() and not ln.startswith("#"):
                sp = ln.split("\t")
                slice_windows.setdefault(sp[0], []).append((int(sp[1]), int(sp[2])))
    contained_v06_svs = []
    with gzip.open(v06_vcf, "rt", encoding="utf-8") as fh:
        for ln in fh:
            if not ln.startswith("#") and ln.strip():
                sp = ln.split("\t")
                kv = dict(x.split("=", 1) for x in sp[7].split(";") if "=" in x)
                svt = kv.get("SVTYPE", "")
                if svt not in ("DEL", "DUP"):
                    continue
                s_p = int(sp[1])
                e_p = int(kv.get("END", s_p + abs(int(kv.get("SVLEN", "0")))))
                if any(ws <= s_p and e_p <= we for ws, we in slice_windows.get(sp[0], [])):
                    contained_v06_svs.append((sp[0], s_p, e_p, svt))
    h8_sv_gz = cache.run_vcf_gz(
        "mod4_hg008_v06_sv", p.giab_ref, hg8_bam_pair,
        ["-N", "HG008T|HG008N", "-c", "1", "-S", "2", "-E", "3", "-z", "-f", "0.01", "-L", "100", "--fisher", "--compat=off", "--vcf-all-alleles", hg8_slices],
        kind="somatic",
    )
    called_h8_svs = []
    with gzip.open(h8_sv_gz, "rt", encoding="utf-8") as fh:
        for ln in fh:
            if not ln.startswith("#") and ln.strip():
                sp = ln.rstrip("\n").split("\t")
                kv = dict(x.split("=", 1) for x in sp[7].split(";") if "=" in x)
                if "Somatic" not in kv.get("STATUS", ""):
                    continue
                if sp[4] in ("<DEL>", "<DUP>", "<INV>"):
                    called_h8_svs.append((sp[0], int(sp[1]), int(kv["END"]), sp[4][1:-1]))
                elif abs(len(sp[4]) - len(sp[3])) >= 100:
                    st = "DEL" if len(sp[3]) > len(sp[4]) else "DUP"
                    called_h8_svs.append((sp[0], int(sp[1]), int(sp[1]) + len(sp[3]) - 1, st))
    v06_matched = 0
    for tc, ts_p, te_p, ttype in contained_v06_svs:
        if any(cc == tc and ctype == ttype and abs(cs_p - ts_p) <= 25 and abs(ce_p - te_p) <= 25 for cc, cs_p, ce_p, ctype in called_h8_svs):
            v06_matched += 1
    if len(contained_v06_svs) < 3 or v06_matched != len(contained_v06_svs):
        raise AssertionError(f"HG008 v0.6 contained SV breakpoint overlap failed: {v06_matched}/{len(contained_v06_svs)}")

    # Modality 5 & 6: COLO829 multimodal + CRISPR & RNA-Seq rna_s6a..c (-i & variant calling)
    cache.run_vcf_gz(
        "mod5_colo829_s1a", p.fa_b37, f"{g4}/multimodal/colo829_s1a_tumor.bam|{g4}/multimodal/colo829_s1a_normal.bam",
        ["-N", "COLO829T|COLO829N", "-R", "1:10000001-10100000", "-f", "0.01", "--fisher", "--compat=off", "--vcf-all-alleles"], kind="somatic",
    )
    cache.run_vcf_gz(
        "mod5_colo829_crispr", p.fa_hg19, f"{p.vj_bams}/Colo829.chr20.CRISPR.bam",
        ["-N", "CRISPR", "-R", "chr20:50454281-50455444", "-J", "50454941", "-j", "5", "-f", "0.01", "--compat=off", "--vcf-all-alleles"],
    )
    rna_splice = cache.run_text(
        "mod6_rna_tcov_i", p.fa_hg38, f"{p.vj_bams}/Tcov_zero_SV_DEL_candidate3.bam",
        ["-R", "chr8:90534611-90620070", "-f", "0.1", "-i", "--emit-intermediate", "--compat=off"], audit=False,
    )
    for rna_sl, rna_reg in (
        ("rna_s6a", "chr8:90474800-90548200"),
        ("rna_s6b", "chr8:90474800-90580200"),
        ("rna_s6c", "chr8:90562300-90619800"),
    ):
        cache.run_vcf_gz(
            f"mod6_{rna_sl}_vcf", p.fa_hg38, os.path.join(g4, "multimodal", f"{rna_sl}.bam"),
            ["-N", rna_sl, "-R", rna_reg, "-f", "0.01", "--compat=off", "--vcf-all-alleles"],
        )
    if len(rna_splice.strip().splitlines()) < 16:
        raise AssertionError("RNA-Seq -i splice junction extraction failed")

    # Modalities 7 & 8: Wave-4 QIAseq SPE (-a 10:0.95 -sp -p) & Illumina PE (-a 10:0.95 -p) amplicon pileup
    qgn_p_bed = os.path.join(p.work_tmp, "mod7_qiaseq_spe_p.amp8.bed")
    with open(os.path.join(w4, "beds", "qiaseq_spe_grch38.amp8.bed"), "r", encoding="utf-8") as src, open(qgn_p_bed, "w", encoding="utf-8") as dst:
        for idx, ln in enumerate(src):
            if idx >= 10:
                break
            dst.write(ln)
    ilm_p_bed = os.path.join(p.work_tmp, "mod8_ilm_pe_p.amp8.bed")
    with open(os.path.join(w4, "beds", "ilm_twoprimer_grch38.amp8.bed"), "r", encoding="utf-8") as src, open(ilm_p_bed, "w", encoding="utf-8") as dst:
        for idx, ln in enumerate(src):
            if idx >= 12:
                break
            dst.write(ln)

    amp_p_dot_total = 0
    for ap_tag, ap_bam, ap_bed, ap_extra in (
        ("mod7_qiaseq_sp_p", os.path.join(w4, "bam", "SampleA_QGN1.bam"), qgn_p_bed, ["-a", "10:0.95", "-sp", "-z", "-p"]),
        ("mod8_ilm_pe_p", os.path.join(w4, "bam", "SampleA_ILM1.bam"), ilm_p_bed, ["-a", "10:0.95", "-z", "-p"]),
    ):
        ap_vj_vcf = cache.run_text(
            f"{ap_tag}_vj_vcf", p.grch38_ref, ap_bam,
            [*ap_extra, "--compat=vardictjava", "--vcf-all-alleles", "--fisher", ap_bed],
            kind="amplicon", audit=False,
        )
        ap_var = cache.run_text(
            f"{ap_tag}_off_var", p.grch38_ref, ap_bam,
            [*ap_extra, "--compat=off", "--vcf-all-alleles", "--fisher", "--emit-intermediate", ap_bed],
            kind="amplicon",
        )
        ap_gz = cache.run_vcf_gz(
            f"{ap_tag}_off_vcf", p.grch38_ref, ap_bam,
            [*ap_extra, "--compat=off", "--vcf-all-alleles", "--fisher", ap_bed],
            index_flag="--csi", kind="amplicon", is_pileup=True,
        )
        with gzip.open(ap_gz, "rt", encoding="utf-8") as gzf:
            ap_recs = [ln.rstrip("\n") for ln in gzf if ln.strip() and not ln.startswith("#")]
        audit_compat_off_var(
            ap_tag, ap_var, "amplicon", fisher=True, vcf_recs=ap_recs, vcf_all_alleles=True
        )
        vj_pos = {
            (ln.split("\t")[0], int(ln.split("\t")[1]))
            for ln in ap_vj_vcf.splitlines()
            if ln.strip() and not ln.startswith("#")
        }
        off_pos = {
            (r.split("\t")[0], int(r.split("\t")[1]))
            for r in ap_recs
        }
        dot_rows = sum(1 for r in ap_recs if r.split("\t")[4] == ".")
        if dot_rows < 100 or not vj_pos.issubset(off_pos):
            raise AssertionError(
                f"{ap_tag}: amplicon -p shifted-ref survival failed (dot_rows={dot_rows}, missing_vj={sorted(vj_pos - off_pos)[:5]})"
            )
        amp_p_dot_total += dot_rows

    # Modality 12a: Colon-joined b1:b2 multi-library pooling + -t within-BAM dedup + =/X CIGARs (samtools calmd -e)
    na_chr1 = os.path.join(p.vj_bams, "NA24385_chr_1_246981051_246981611.bam")
    single_raw = cache.run_text("mod12_colon_1x_raw", p.fa_b37, na_chr1, ["-N", "S1", "-R", "1:246981262-246981262", "-f", "0.02", "--compat=off"])
    single_dedup = cache.run_text("mod12_colon_1x_t", p.fa_b37, na_chr1, ["-N", "S1", "-R", "1:246981262-246981262", "-t", "-f", "0.02", "--compat=off"])
    colon_dedup = cache.run_text("mod12_colon_2x_t", p.fa_b37, f"{na_chr1}:{na_chr1}", ["-N", "S1", "-R", "1:246981262-246981262", "-t", "-f", "0.02", "--compat=off"])
    s_raw_dp = [int(dict(x.split("=", 1) for x in ln.split("\t")[7].split(";") if "=" in x)["DP"]) for ln in single_raw.splitlines() if not ln.startswith("#") and ln.strip()]
    s1_dp = [int(dict(x.split("=", 1) for x in ln.split("\t")[7].split(";") if "=" in x)["DP"]) for ln in single_dedup.splitlines() if not ln.startswith("#") and ln.strip()]
    c2_dp = [int(dict(x.split("=", 1) for x in ln.split("\t")[7].split(";") if "=" in x)["DP"]) for ln in colon_dedup.splitlines() if not ln.startswith("#") and ln.strip()]
    if not s1_dp or not (s1_dp[0] < s_raw_dp[0] and c2_dp[0] == 2 * s1_dp[0]):
        raise AssertionError(f"Colon-joined b1:b2 + -t within-BAM dedup invariant failed: raw={s_raw_dp}, dedup={s1_dp}, colon={c2_dp}")

    plat_eqx_bam = os.path.join(g4, "multimodal", "platinum_s2a.bam")
    plat_bed = os.path.join(tc_dir, "platinum_eval.bed")
    plat_m_bam = os.path.join(p.work_tmp, "platinum_s2a_m.bam")
    sam_txt = subprocess.check_output(["samtools", "view", "-h", plat_eqx_bam], text=True)
    cigar_pat = re.compile(r"(\d+)([MIDNSHP=X])")
    m_sam_lines: list[str] = []
    for sln in sam_txt.splitlines():
        if sln.startswith("@") or not sln.strip():
            m_sam_lines.append(sln)
            continue
        scols = sln.split("\t")
        m_ops: list[tuple[int, str]] = []
        for n_s, op_ch in cigar_pat.findall(scols[5]):
            op_m = "M" if op_ch in ("=", "X") else op_ch
            n_v = int(n_s)
            if m_ops and m_ops[-1][1] == op_m:
                m_ops[-1] = (m_ops[-1][0] + n_v, op_m)
            else:
                m_ops.append((n_v, op_m))
        scols[5] = "".join(f"{n_v}{op_ch}" for n_v, op_ch in m_ops)
        m_sam_lines.append("\t".join(scols))
    with open(plat_m_bam, "wb") as out_b:
        subprocess.run(["samtools", "view", "-b", "-"], input="\n".join(m_sam_lines) + "\n", text=True, stdout=out_b, check=True)
    subprocess.run(["samtools", "index", plat_m_bam], check=True)
    eqx_args = ["-N", "HG001", "-c", "1", "-S", "2", "-E", "3", "-z", "-x", "35000", "-m", "1000", "--vcf-max-nm", "1000", "-f", "0.15", "--compat=off", "--vcf-all-alleles", plat_bed]
    eqx_vcf = cache.run_text("mod12_eqx_cigar", p.grch38_ref, plat_eqx_bam, eqx_args)
    m_vcf = cache.run_text("mod12_m_cigar", p.grch38_ref, plat_m_bam, eqx_args)
    e_sites = [ln for ln in eqx_vcf.splitlines() if not ln.startswith("#") and ln.strip()]
    m_sites = [ln for ln in m_vcf.splitlines() if not ln.startswith("#") and ln.strip()]
    if len(e_sites) < 45 or e_sites != m_sites:
        raise AssertionError(f"Real minimap2 --eqx (= / X CIGAR) parity failed against M CIGAR ({len(e_sites)} vs {len(m_sites)})")

    with metrics_lock:
        metrics["wave6"] = {
            "hg002_1_5m_gt_f1": g2_off_all[3],
            "hg002_1_5m_snp_f1": g2_off_snp[3] if g2_off_snp else None,
            "hg002_1_5m_indel_f1": g2_off_ind[3] if g2_off_ind else None,
            "hg002_2x250_euchr_gt_f1": hg2_250_gt[3],
            "hg002_2x250_mhc_gt": {"tp": hg2_250_mhc_gt[0], "fp": hg2_250_mhc_gt[1], "fn": hg2_250_mhc_gt[2], "f1": hg2_250_mhc_gt[3], "precision": round(mhc_prec, 4)},
            "cmrg_gt": {"tp": cmrg_gt[0], "fp": cmrg_gt[1], "fn": cmrg_gt[2], "f1": cmrg_gt[3]},
            "mosaic_sq": {"tp": mosaic_sq[0], "fp": mosaic_sq[1], "fn": mosaic_sq[2], "f1": mosaic_sq[3]},
            "diploid_f1s": mod2_f1s,
            "seqc2_f1": s_off[3],
            "seqc2_u_f1": su_off[3],
            "seqc2_snp_f1": s_off_snp[3] if s_off_snp else None,
            "seqc2_indel_f1": s_off_ind[3] if s_off_ind else None,
            "hg008_v03_f1": metrics.get("hg008_v03_f1"),
            "hg008_v04_f1": metrics.get("hg008_v04_f1"),
            "hg008_v06_sv_matched": f"{v06_matched}/{len(contained_v06_svs)}",
            "u_rebalanced_loci": rebalanced,
            "amp_pileup_dot_rows": amp_p_dot_total,
            "eqx_sites": len(e_sites),
        }
    print(
        f"  [PASS Wave-6 & 12-Modality] HG002 GT F1={g2_off_all[3]:.4f} (2x250 euchr={hg2_250_gt[3]:.4f}, MHC TP={hg2_250_mhc_gt[0]}/prec={mhc_prec:.4f}, CMRG TP={cmrg_gt[0]}, Mosaic TP={mosaic_sq[0]}), "
        f"HG005/Plat/b37={mod2_f1s}, SEQC2 F1={s_off[3]:.4f}, HG008 BAM==CRAM & v0.6 SV={v06_matched}/{len(contained_v06_svs)}, -a -p dot={amp_p_dot_total}, --eqx={len(e_sites)}"
    )


# =============================================================================
# Suite 5: Tiled == Untiled & Spill == Unspilled (--suite tiling)
# =============================================================================
def run_suite_tiling(
    p: Paths,
    cache: RunCache,
    metrics: dict[str, object],
    metrics_lock: threading.Lock,
) -> None:
    print("[Tiling & Spill] Verifying Tiled (50000) == Untiled (50000000), Spill == Unspilled, & bounded-horizon XFAILs...")
    g4 = str(p.g4_dir)

    # Build SEQC2 6-slice BED excluding slice_4_5 (which is checked separately below)
    seqc2_6_bed = os.path.join(p.work_tmp, "seqc2_6slices.bed")
    slice45_bed = os.path.join(p.work_tmp, "seqc2_slice45.bed")
    with open(os.path.join(g4, "seqc2", "slices.bed"), "r", encoding="utf-8") as fh, \
         open(seqc2_6_bed, "w", encoding="utf-8") as out6, \
         open(slice45_bed, "w", encoding="utf-8") as out45:
        for ln in fh:
            if "slice_4_5" in ln:
                out45.write(ln)
            elif ln.strip():
                out6.write(ln)

    # Subset of HG002 eval.bed (~4 Mbp) for fast multi-interval check + COLO829 200kbp + HG008 slices + SEQC2 6 slices + -k 0 --nosv multi-tile
    hg2_eval_sub = os.path.join(p.work_tmp, "hg002_eval_sub.bed")
    with open(os.path.join(g4, "hg002", "eval.bed"), "r", encoding="utf-8") as src, open(hg2_eval_sub, "w", encoding="utf-8") as dst:
        for idx, ln in enumerate(src):
            if idx >= 120:
                break
            dst.write(ln)

    strict_cases = [
        ("COLO829_s1a_200k", p.fa_b37, f"{g4}/multimodal/colo829_s1a_tumor.bam|{g4}/multimodal/colo829_s1a_normal.bam",
         ["-R", "1:10000001-10200000", "-f", "0.01", "--fisher", "--compat=off"]),
        ("SEQC2_slices_123567", p.grch38_ref, f"{g4}/seqc2/HCC1395.bam|{g4}/seqc2/HCC1395BL.bam",
         ["-c", "1", "-S", "2", "-E", "3", "-g", "4", "-z", "-f", "0.01", "--fisher", "--compat=off", seqc2_6_bed]),
        ("HG008_slices", p.giab_ref, f"{g4}/hg008/HG008-T.slice.bam|{g4}/hg008/HG008-N.slice.bam",
         ["-c", "1", "-S", "2", "-E", "3", "-z", "-f", "0.01", "--fisher", "--compat=off", os.path.join(g4, "hg008", "slices.bed")]),
        ("HG002_eval_sub", p.giab_ref, os.path.join(g4, "hg002", "HG002.chr20.bam"),
         ["-c", "1", "-S", "2", "-E", "3", "-z", "-f", "0.01", "--compat=off", hg2_eval_sub]),
        ("HG002_1_5M_k0_nosv", p.giab_ref, os.path.join(g4, "hg002", "HG002.chr20.bam"),
         ["-R", "chr20:10000001-11500000", "-f", "0.01", "-k", "0", "--nosv", "--compat=off"]),
        ("HG002_MHC_200k_k0_nosv", p.grch38_ref, os.path.join(g4, "hg002", "HG002.GRCh38.2x250.bam"),
         ["-R", "chr6:29900001-30100000", "-f", "0.01", "-k", "0", "--nosv", "--compat=off"]),
        ("SEQC2_slice45_k0", p.grch38_ref, f"{g4}/seqc2/HCC1395.bam|{g4}/seqc2/HCC1395BL.bam",
         ["-c", "1", "-S", "2", "-E", "3", "-g", "4", "-z", "-f", "0.01", "-k", "0", "--fisher", "--compat=off", slice45_bed]),
    ]

    strict_passed = 0
    for cname, rfa, bspec, cargs in strict_cases:
        u_txt = cache.run_text(f"tile_u_{cname}", rfa, bspec, [*cargs, "--tile-size", "50000000"], audit=False)
        t_txt = cache.run_text(f"tile_t_{cname}", rfa, bspec, [*cargs, "--tile-size", "50000"], audit=False)
        if u_txt != t_txt:
            raise AssertionError(f"Strict Tiled (50000) == Untiled (50000000) failed on {cname}")
        strict_passed += 1

    # Tracked bounded-horizon XFAIL checks (Findings 2 & 3: XPASS-loud if 0 diffs, FAIL if regressed)
    horizon_xfails = [
        ("HG002_1_5M_cross_seam_sv", p.giab_ref, os.path.join(g4, "hg002", "HG002.chr20.bam"),
         ["-R", "chr20:10000001-11500000", "-f", "0.01", "--compat=off"], 45),
        ("HG002_MHC_cross_seam_sv", p.grch38_ref, os.path.join(g4, "hg002", "HG002.GRCh38.2x250.bam"),
         ["-R", "chr6:29900001-30100000", "-f", "0.01", "--compat=off"], 12),
        ("SEQC2_slice45_seed_horizon", p.grch38_ref, f"{g4}/seqc2/HCC1395.bam|{g4}/seqc2/HCC1395BL.bam",
         ["-c", "1", "-S", "2", "-E", "3", "-g", "4", "-z", "-f", "0.01", "--fisher", "--compat=off", slice45_bed], 10),
    ]
    xfail_verified = 0
    for xname, rfa, bspec, xargs, max_sym_diff in horizon_xfails:
        xu = cache.run_text(f"xfail_u_{xname}", rfa, bspec, [*xargs, "--tile-size", "50000000"], audit=False)
        xt = cache.run_text(f"xfail_t_{xname}", rfa, bspec, [*xargs, "--tile-size", "50000"], audit=False)
        su = {ln for ln in xu.splitlines() if not ln.startswith("#") and ln.strip()}
        st = {ln for ln in xt.splitlines() if not ln.startswith("#") and ln.strip()}
        sym_d = len(su ^ st)
        if sym_d == 0:
            raise AssertionError(f"XPASS on {xname}: bounded-horizon tiling difference resolved; promote from XFAIL to PASS!")
        if sym_d > max_sym_diff:
            raise AssertionError(f"Bounded-horizon regression on {xname}: symmetric diff {sym_d} > ceiling {max_sym_diff}")
        xfail_verified += 1

    # Spill == Unspilled (> 100 MiB candidate stream exceeding 8 MiB DEFAULT_VCF_SPILL_BYTES)
    big_smp = "SPILL_" + ("X" * 99990)
    s1a_bed = str(p.root_dir / "harness" / "slices" / "s1a.bed")
    na_20_bam = os.path.join(p.vj_bams, "NA24385_chr_20_20156208_20234556.bam")
    unspilled = cache.run_text("spill_norm", p.fa_b37, na_20_bam, ["-N", "NORM_SMP", "-c", "1", "-S", "2", "-E", "3", "-z", "-f", "0.001", "--compat=off", "--vcf-all-alleles", s1a_bed])
    with cache.sem:
        spilled = subprocess.check_output(
            [p.vdr_bin, "-G", p.fa_b37, "-b", na_20_bam, "-N", big_smp, "-c", "1", "-S", "2", "-E", "3", "-z", "-f", "0.001", "--compat=off", "--vcf-all-alleles", "-th", str(cache.threads_per_job), s1a_bed],
            text=True,
        ).replace(big_smp, "NORM_SMP")
    if unspilled != spilled or len(unspilled.splitlines()) < 1000:
        raise AssertionError("Spill == Unspilled (>100 MiB candidate stream) mismatch")

    with metrics_lock:
        metrics["tiling"] = {
            "strict_pass": strict_passed,
            "bounded_horizon_xfail": xfail_verified,
            "spill_lines": len(unspilled.splitlines()),
        }
    print(f"  [PASS Tiling & Spill] {strict_passed}/7 strict multi-tile PASS, {xfail_verified}/3 bounded-horizon XFAILs pinned, Spill==Unspilled ({len(unspilled.splitlines())} lines)")


# =============================================================================
# Suite 6: 38-Switch Real-Data CompatProfile Ablation (only_1..38 & loo_1..38)
# =============================================================================
def _parse_cargo_lock_packages(lock_path: str | Path) -> set[tuple[str, str, str]]:
    with open(lock_path, "r", encoding="utf-8") as fh:
        raw = fh.read()
    pkgs: set[tuple[str, str, str]] = set()
    for blk in raw.split("\n[[package]]\n")[1:]:
        name = ver = csum = ""
        for ln in blk.strip().splitlines():
            if ln.startswith('name = "'):
                name = ln.split('"', 2)[1]
            elif ln.startswith('version = "'):
                ver = ln.split('"', 2)[1]
            elif ln.startswith('checksum = "'):
                csum = ln.split('"', 2)[1]
        if name:
            pkgs.add((name, ver, csum))
    return pkgs


def _write_locked_ablation_cargo_lock(
    root_lock_path: Path, dst_lock_path: str, manifest_path: str, env: dict[str, str]
) -> None:
    shutil.copyfile(root_lock_path, dst_lock_path)
    pinned = _parse_cargo_lock_packages(root_lock_path)
    subprocess.run(
        [
            "cargo",
            "metadata",
            "--offline",
            "--format-version",
            "1",
            "--manifest-path",
            manifest_path,
        ],
        env=env,
        stdout=subprocess.DEVNULL,
        check=True,
    )
    resolved = {
        pkg for pkg in _parse_cargo_lock_packages(dst_lock_path) if pkg[0] != "ablation_driver"
    }
    drifted = resolved - pinned
    if drifted:
        raise AssertionError(
            f"ablation_driver Cargo.lock drifted from workspace Cargo.lock: {sorted(drifted)}"
        )


def run_suite_ablation(
    p: Paths,
    cache: RunCache,
    metrics: dict[str, object],
    metrics_lock: threading.Lock,
) -> None:
    print("[Ablation] Compiling & running 38-switch real-data (only_1..38, loo_1..38) ablation driver...")
    drv_dir = os.path.join(p.work_tmp, "ablation_driver")
    os.makedirs(os.path.join(drv_dir, "src"), exist_ok=True)
    core_path = (p.root_dir / "crates" / "vardictrust-core").as_posix()

    manifest_path = os.path.join(drv_dir, "Cargo.toml")
    with open(manifest_path, "w", encoding="utf-8") as fh:
        fh.write(
            f'[package]\nname = "ablation_driver"\nversion = "0.1.0"\nedition = "2021"\n\n'
            f'[dependencies]\nvardictrust-core = {{ path = "{core_path}" }}\n'
        )

    rust_src = r"""
use std::env;
use vardictrust_core::config::{
    CompatProfile, Config, COMPAT_OFF, COMPAT_SWITCH_COUNT, COMPAT_VARDICTJAVA_1_8_3,
};

fn run_mask(prof: CompatProfile, g4: &str, w4: &str, grch38: &str, giab: &str, hg19: &str, vjb: &str, threads: usize) -> (usize, u64) {
    let mut total_rows = 0usize;
    let mut hash: u64 = 0xcbf29ce484222325;
    let fnv = |h: &mut u64, bytes: &[u8]| {
        for &b in bytes {
            *h ^= u64::from(b);
            *h = h.wrapping_mul(0x100000001b3);
        }
    };

    // 1. SEQC2 HCC1395|BL somatic slice + 2. GIAB HG002 germline -u + 3. BRCA2 adaptor + 4. Colo829 SV + 5. Illumina paired amplicon
    let cases: [(&str, String, Option<&str>, Option<String>, bool, bool, Vec<String>); 5] = [
        (
            grch38,
            format!("{g4}/seqc2/HCC1395.bam|{g4}/seqc2/HCC1395BL.bam"),
            Some("chr17:74747665-74765000"),
            None,
            false,
            false,
            Vec::new(),
        ),
        (
            giab,
            format!("{g4}/hg002/HG002.chr20.bam"),
            Some("chr20:10000001-10025000"),
            None,
            true,
            false,
            Vec::new(),
        ),
        (
            grch38,
            format!("{w4}/brca2_adaptor/adaptor_insertion.bam"),
            Some("chr13:32683877-32684023"),
            None,
            false,
            false,
            vec![
                "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA".to_string(),
                "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT".to_string(),
            ],
        ),
        (
            hg19,
            format!("{vjb}/Colo829.chr20.DUPDiscSV.bam"),
            Some("chr20:47124929-47131802"),
            None,
            false,
            false,
            Vec::new(),
        ),
        (
            grch38,
            format!("{w4}/bam/SampleA_ILM1.bam|{w4}/bam/SampleB_ILM1.bam"),
            None,
            Some(format!("{w4}/beds/ilm_twoprimer_grch38.amp8.bed")),
            false,
            true,
            Vec::new(),
        ),
    ];

    for (rfa, bam, reg, bed, u_mode, is_amp, adapt) in &cases {
        let mut cfg = Config {
            reference: (*rfa).to_string(),
            bam: bam.clone(),
            sample: if bam.contains('|') { "TUMOR|NORMAL".to_string() } else { "SAMPLE;1".to_string() },
            sample_explicit: true,
            region: reg.map(str::to_string),
            bed: bed.clone(),
            zero_based: true,
            zero_based_defined: true,
            freq: 0.01,
            af_threshold: 0.01,
            fisher: true,
            vcf_all_alleles: true,
            unique_mode_alignment_enabled: *u_mode,
            amplicon: is_amp.then(|| "10:0.95".to_string()),
            adaptor: adapt.clone(),
            threads,
            compat: prof,
            ..Config::default()
        };
        if *is_amp && !prof.paired_somatic_amplicon {
            let mut sink = Vec::new();
            assert!(
                matches!(
                    vardictrust_core::run_with_writer(&cfg, &mut sink),
                    Err(vardictrust_core::types::VardictError::Unsupported(_))
                ),
                "expected paired amplicon VCF to return VardictError::Unsupported when paired_somatic_amplicon=false"
            );
            fnv(&mut hash, b"ERR_UNSUPPORTED_378\n");
            continue;
        }
        if !adapt.is_empty() {
            cfg.include_n = true;
            cfg.base_qual = 25.0;
        }
        let mut out = Vec::new();
        vardictrust_core::run_with_writer(&cfg, &mut out).expect("ablation run_with_writer");
        for line in out.split(|&b| b == b'\n') {
            if !line.is_empty() && line[0] != b'#' {
                total_rows += 1;
                fnv(&mut hash, line);
                fnv(&mut hash, b"\n");
            } else if line.starts_with(b"##FORMAT=<ID=AD,") {
                fnv(&mut hash, line);
            }
        }
    }
    (total_rows, hash)
}

fn main() {
    let args: Vec<String> = env::args().collect();
    assert!(args.len() >= 7, "ablation_driver expected >= 6 arguments, got {}", args.len() - 1);
    let (g4, w4, grch38, giab, hg19, vjb) = (&args[1], &args[2], &args[3], &args[4], &args[5], &args[6]);
    let threads: usize = args.get(7).and_then(|s| s.parse().ok()).unwrap_or(4).max(1);
    let switches = CompatProfile::all_switches();
    assert_eq!(switches.len(), COMPAT_SWITCH_COUNT);

    let (vj_rows, vj_hash) = run_mask(COMPAT_VARDICTJAVA_1_8_3, g4, w4, grch38, giab, hg19, vjb, threads);
    let (off_rows, off_hash) = run_mask(COMPAT_OFF, g4, w4, grch38, giab, hg19, vjb, threads);
    assert_ne!(vj_hash, off_hash, "COMPAT_OFF must differ from COMPAT_VARDICTJAVA_1_8_3 on real data");

    let mut active_deltas = 0usize;
    let mut loo33 = (0usize, 0u64);
    let mut loo37 = (0usize, 0u64);
    for spec in &switches {
        let idx = spec.index;
        let sw_num = idx + 1;
        let only_p = CompatProfile::only_switch(idx).expect("only_switch");
        let loo_p = CompatProfile::without_switch(idx).expect("without_switch");
        let (_o_rows, o_hash) = run_mask(only_p, g4, w4, grch38, giab, hg19, vjb, threads);
        let (l_rows, l_hash) = run_mask(loo_p, g4, w4, grch38, giab, hg19, vjb, threads);
        if o_hash != vj_hash || l_hash != off_hash {
            active_deltas += 1;
        }
        if sw_num == 33 {
            loo33 = (l_rows, l_hash);
        }
        if sw_num == 37 {
            loo37 = (l_rows, l_hash);
        }
    }
    println!(
        "{{\"switches\":{},\"vj_rows\":{},\"off_rows\":{},\"active_deltas\":{},\"loo_33_rows\":{},\"loo_37_rows\":{}}}",
        switches.len(), vj_rows, off_rows, active_deltas, loo33.0, loo37.0
    );
}
"""
    with open(os.path.join(drv_dir, "src", "main.rs"), "w", encoding="utf-8") as fh:
        fh.write(rust_src)

    env = os.environ.copy()
    env["CARGO_BUILD_JOBS"] = str(cache.threads_per_job)
    env["RAYON_NUM_THREADS"] = str(cache.threads_per_job)
    tgt_dir = os.path.join(p.work_tmp, "ablation_target")
    with cache.sem:
        _write_locked_ablation_cargo_lock(
            p.root_dir / "Cargo.lock",
            os.path.join(drv_dir, "Cargo.lock"),
            manifest_path,
            env,
        )
        subprocess.run(
            [
                "cargo",
                "build",
                "-q",
                "--release",
                "--locked",
                "--offline",
                "--manifest-path",
                manifest_path,
                "--target-dir",
                tgt_dir,
            ],
            env=env,
            check=True,
        )
        drv_bin = os.path.join(tgt_dir, "release", "ablation_driver")
        out_json = subprocess.check_output(
            [
                drv_bin,
                str(p.g4_dir),
                str(p.wave4_dir),
                p.grch38_ref,
                p.giab_ref,
                p.fa_hg19,
                p.vj_bams,
                str(cache.threads_per_job),
            ],
            env=env,
            text=True,
        ).strip()
    shutil.rmtree(tgt_dir, ignore_errors=True)
    ab_data = json.loads(out_json.splitlines()[-1])
    if ab_data["switches"] != 38 or ab_data["off_rows"] < 50 or ab_data["active_deltas"] < 15:
        raise AssertionError(f"38-switch real-data ablation check failed: {ab_data}")
    with metrics_lock:
        metrics["ablation"] = ab_data
    print(
        f"  [PASS Ablation] 38/38 only_1..38 & loo_1..38 masks green "
        f"(vj={ab_data['vj_rows']} -> off={ab_data['off_rows']} rows, active_deltas={ab_data['active_deltas']}, "
        f"loo_33={ab_data['loo_33_rows']}, loo_37={ab_data['loo_37_rows']})"
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Single-pass DAG public-data omnibus verifier for VarDictRust."
    )
    parser.add_argument(
        "--suite",
        choices=VALID_SUITES,
        default="all",
        help="Verification suite to execute (default: all)",
    )
    parser.add_argument(
        "--include-exome",
        "--exome",
        dest="include_exome",
        action="store_true",
        help="Also sweep full 90.5 Mbp SEQC2 exome.engine.bed",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=int(os.environ.get("VERIFY_JOBS", "2")),
        help="Max concurrent vardictrust worker processes (default: 2)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=int(os.environ.get("VERIFY_THREADS", "4")),
        help="Threads (-th) per vardictrust invocation (default: 4; jobs*threads <= 12)",
    )
    parser.add_argument(
        "--rtg-jobs",
        type=int,
        default=2,
        help="Max concurrent rtg vcfeval jobs (default: 2, capped <= 2, RTG_MEM=2G)",
    )
    parser.add_argument(
        "--json-summary",
        default=None,
        help="Optional path outside benchmarks/evidence/ to write JSON summary metrics",
    )
    args = parser.parse_args()

    try:
        os.nice(19)
    except OSError:
        pass

    ensure_tool_paths()
    for req_tool in ("samtools", "bcftools", "tabix", "bgzip", "rtg", "perl", "gzip"):
        if shutil.which(req_tool) is None:
            print(f"ERROR: required tool {req_tool!r} not found on PATH", file=sys.stderr)
            return 1

    script_dir = Path(__file__).resolve().parent
    root_dir = script_dir.parent
    data_root = resolve_data_root(root_dir)

    if args.json_summary:
        js_abs = Path(args.json_summary).resolve()
        for forbidden_dir in (root_dir / "benchmarks" / "evidence", root_dir / "golden", data_root / "benchmarks" / "evidence", data_root / "golden"):
            if str(js_abs).startswith(str(forbidden_dir.resolve())):
                print(f"ERROR: refusing to write --json-summary inside protected directory {forbidden_dir}", file=sys.stderr)
                return 1

    vdr_bin = (
        os.environ.get("VARDICTRUST_BIN")
        or os.environ.get("VARDICT_BIN")
        or str(root_dir / "target" / "release" / "vardictrust")
    )
    if not (os.path.isfile(vdr_bin) and os.access(vdr_bin, os.X_OK)):
        print(f"ERROR: VARDICTRUST_BIN not executable at {vdr_bin!r}", file=sys.stderr)
        return 1

    jobs = max(1, args.jobs)
    threads = max(1, args.threads)
    if jobs * threads > 12:
        threads = max(1, 12 // jobs)
    rtg_jobs = max(1, min(2, args.rtg_jobs))

    g4_dir = Path(os.environ.get("DATA_DIR") or (data_root / "scratch" / "g4"))
    wave4_dir = Path(os.environ.get("WAVE4_DIR") or (g4_dir / "wave4"))
    if not (wave4_dir / ".done").is_file():
        env = os.environ.copy()
        env["DATA_DIR"] = str(g4_dir)
        env["WAVE4_DIR"] = str(wave4_dir)
        subprocess.run(["bash", str(script_dir / "get_wave4_public_slices.sh")], env=env, check=True)

    fasta_root = data_root if (data_root / "VarDictJava" / "testdata" / "fastas" / "hg19.fa").is_file() else root_dir
    if not (fasta_root / "VarDictJava" / "testdata" / "fastas" / "hg19.fa").is_file():
        subprocess.run([sys.executable, str(root_dir / "harness" / "convert_fastas.py")], check=True)
        fasta_root = root_dir

    giab_ref = str(g4_dir / "ref" / "giabv3" / "GRCh38_GIABv3.fa")
    if not os.path.isfile(giab_ref) and (g4_dir / "ref" / "giabv3" / "GRCh38_GIABv3.fasta").is_file():
        giab_ref = str(g4_dir / "ref" / "giabv3" / "GRCh38_GIABv3.fasta")

    signal.signal(signal.SIGTERM, lambda *_: sys.exit(143))
    work_tmp = tempfile.mkdtemp(prefix="vdr_verify_pub_")
    t0 = time.monotonic()
    try:
        paths = Paths(
            root_dir=root_dir,
            data_root=data_root,
            g4_dir=g4_dir,
            wave4_dir=wave4_dir,
            vdr_bin=vdr_bin,
            grch38_ref=str(g4_dir / "ref" / "d1vd1" / "GRCh38.d1.vd1.fa"),
            giab_ref=giab_ref,
            fa_b37=str(fasta_root / "VarDictJava" / "testdata" / "fastas" / "GRCh37.fa"),
            fa_hg19=str(fasta_root / "VarDictJava" / "testdata" / "fastas" / "hg19.fa"),
            fa_hg38=str(fasta_root / "VarDictJava" / "testdata" / "fastas" / "hg38.fa"),
            fa_hs38dh=str(fasta_root / "VarDictJava" / "testdata" / "fastas" / "hs38DH.fa"),
            fa_hg1k=str(fasta_root / "VarDictJava" / "testdata" / "fastas" / "hg1k.fa"),
            sdf_d1vd1=str(g4_dir / "ref" / "d1vd1" / "GRCh38.d1.vd1.sdf"),
            sdf_giab=str(g4_dir / "ref" / "giabv3" / "GRCh38_GIABv3.sdf"),
            sdf_b37=str(g4_dir / "ref" / "grch37.sdf"),
            vj_bams=str(
                data_root
                / "VarDictJava"
                / "src"
                / "test"
                / "resources"
                / "com"
                / "astrazeneca"
                / "vardict"
                / "integrationtests"
            ),
            perl_valid=str(data_root / "VarDictJava" / "VarDict" / "var2vcf_valid.pl"),
            perl_paired=str(data_root / "VarDictJava" / "VarDict" / "var2vcf_paired.pl"),
            work_tmp=work_tmp,
        )
        for req_sdf in (paths.sdf_d1vd1, paths.sdf_giab, paths.sdf_b37):
            if not os.path.isdir(req_sdf):
                raise RuntimeError(f"Missing pre-built on-disk SDF directory: {req_sdf}")

        audit_stats = AuditStats()
        cache = RunCache(vdr_bin, work_tmp, threads_per_job=threads, max_engine_jobs=jobs, audit_stats=audit_stats)
        rtg_eval = RtgEvaluator(work_tmp, max_rtg_jobs=rtg_jobs)
        metrics: dict[str, object] = {}
        metrics_lock = threading.Lock()

        print("================================================================================")
        print(f" VarDictRust Public-Data Omnibus Verifier (--suite {args.suite})")
        print(f" Binary: {vdr_bin} | Concurrency: jobs={jobs} x -th {threads} (rtg_jobs={rtg_jobs})")
        print("================================================================================")

        suite = args.suite
        if suite == "all":
            # Execute independent DAG branches concurrently while RunCache + Semaphores bound CPU/RAM!
            tasks: list[Callable[[], None]] = [
                lambda: run_suite_wave4(paths, cache, metrics, metrics_lock),
                lambda: run_suite_wave5(paths, cache, rtg_eval, metrics, metrics_lock),
                lambda: run_suite_presets(paths, cache, metrics, metrics_lock),
                lambda: run_suite_wave6(paths, cache, rtg_eval, metrics, metrics_lock, include_exome=args.include_exome),
                lambda: run_suite_tiling(paths, cache, metrics, metrics_lock),
                lambda: run_suite_ablation(paths, cache, metrics, metrics_lock),
            ]
            with concurrent.futures.ThreadPoolExecutor(max_workers=3) as dag_pool:
                futs = [dag_pool.submit(fn) for fn in tasks]
                for fut in futs:
                    fut.result()
        elif suite == "wave4":
            run_suite_wave4(paths, cache, metrics, metrics_lock)
        elif suite == "wave5":
            run_suite_wave5(paths, cache, rtg_eval, metrics, metrics_lock)
        elif suite == "presets":
            run_suite_presets(paths, cache, metrics, metrics_lock)
        elif suite == "wave6":
            run_suite_wave6(paths, cache, rtg_eval, metrics, metrics_lock, include_exome=args.include_exome)
        elif suite == "tiling":
            run_suite_tiling(paths, cache, metrics, metrics_lock)
        elif suite == "ablation":
            run_suite_ablation(paths, cache, metrics, metrics_lock)

        elapsed_s = time.monotonic() - t0
        summary = {
            "suite": suite,
            "elapsed_s": round(elapsed_s, 2),
            "unique_spawns": cache.unique_spawns,
            "logical_requests": cache.logical_requests,
            "cache_hits": cache.logical_requests - cache.unique_spawns,
            "audited_vcf_files": audit_stats.vcf_files,
            "audited_vcf_rows": audit_stats.vcf_rows,
            "audited_var_files": audit_stats.var_files,
            "audited_var_rows": audit_stats.var_rows,
            "bcftools_norm_checks": audit_stats.bcftools_norm_checks,
            "index_queries": audit_stats.index_queries,
            "metrics": metrics,
        }
        if args.json_summary:
            with open(args.json_summary, "w", encoding="utf-8") as fh:
                json.dump(summary, fh, indent=2, sort_keys=True)

        print("================================================================================")
        print(
            f" PASS verify_public_data (--suite {suite}) in {elapsed_s:.1f}s | "
            f"RunCache: {cache.unique_spawns} spawns for {cache.logical_requests} requests "
            f"({cache.logical_requests - cache.unique_spawns} cache hits) | "
            f"Universal Auditor: {audit_stats.vcf_rows} .vcf rows, {audit_stats.var_rows} .var rows, "
            f"{audit_stats.bcftools_norm_checks} bcftools norm checks, {audit_stats.index_queries} index queries"
        )
        print("================================================================================")
        return 0
    finally:
        shutil.rmtree(work_tmp, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
