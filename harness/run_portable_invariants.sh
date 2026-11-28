#!/usr/bin/env bash
# =============================================================================
# harness/run_portable_invariants.sh
# Fast (<12s), zero-JVM, zero-scratch/g4 portable CI invariant gate executing
# strictly against git-tracked VarDictJava test resources & harness/slices/*.bed.
#
# Enforces 4 end-to-end invariants on every PR and differential.yml run:
#   1. --compat=off --vcf-all-alleles --output out.vcf.gz --csi +
#      --preset {wgs,wes,amplicon,ctdna} + Universal .var/.vcf schema &
#      depth/ALD accounting + bcftools norm --check-ref e + tabix/bcftools view -r
#   2. Tiled (--tile-size 200/300) == Untiled (--tile-size 50000000) across
#      Simple, Paired Somatic, Amplicon (-a & -sp), Pileup (-p), Splice (-i),
#      Spanning-D seams, and SV (<DEL>, <INV>) under --compat=vardictjava|off
#   3. Live portable BAM == CRAM (samtools view -C -T) across single & paired
#      modes under both --compat=vardictjava and --compat=off
#   4. --vcf-all-alleles vs Perl var2vcf_valid.pl -A / var2vcf_paired.pl -A
#      under --compat=vardictjava (piping vardictrust --emit-intermediate)
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(CDPATH='' cd -- "$(dirname -- "$0")" && pwd)"
ROOT_DIR="$(CDPATH='' cd -- "$SCRIPT_DIR/.." && pwd)"
. "$SCRIPT_DIR/lib/portable.sh"

ensure_g4_tools_path
if [ -n "${CONDA_ROOT:-}" ] && [ -d "$CONDA_ROOT/envs/vardict-oracle/bin" ]; then
  case ":$PATH:" in
    *":$CONDA_ROOT/envs/vardict-oracle/bin:"*) ;;
    *) export PATH="$CONDA_ROOT/envs/vardict-oracle/bin:$PATH" ;;
  esac
elif [ -n "${HOME:-}" ] && [ -d "$HOME/miniconda3/envs/vardict-oracle/bin" ]; then
  case ":$PATH:" in
    *":$HOME/miniconda3/envs/vardict-oracle/bin:"*) ;;
    *) export PATH="$HOME/miniconda3/envs/vardict-oracle/bin:$PATH" ;;
  esac
fi

for _req_tool in bcftools tabix samtools perl; do
  if ! command -v "$_req_tool" >/dev/null 2>&1; then
    echo "ERROR: required tool '$_req_tool' not found on PATH" >&2
    exit 1
  fi
done
unset _req_tool

VARDICTRUST_BIN="${VARDICTRUST_BIN:-${VARDICT_BIN:-$ROOT_DIR/target/release/vardictrust}}"
if [ ! -x "$VARDICTRUST_BIN" ]; then
  cargo build -q --release -p vardictrust --manifest-path "$ROOT_DIR/Cargo.toml"
fi

ASSET_ROOT="$ROOT_DIR"
DATA_ROOT="$(resolve_data_root "$ROOT_DIR")"
if [ ! -f "$ASSET_ROOT/VarDictJava/testdata/fastas/hg19.fa" ] \
  || [ ! -d "$ASSET_ROOT/VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests" ] \
  || [ ! -f "$ASSET_ROOT/VarDictJava/VarDict/var2vcf_valid.pl" ]; then
  if [ -f "$DATA_ROOT/VarDictJava/testdata/fastas/hg19.fa" ] && [ -f "$DATA_ROOT/VarDictJava/testdata/fastas/hs38DH.fa" ]; then
    ASSET_ROOT="$DATA_ROOT"
  else
    python3 "$ROOT_DIR/harness/convert_fastas.py"
  fi
fi

WORK_TMP="$(mktemp -d "${TMPDIR:-/tmp}/vdr_portable_inv.XXXXXX")"
trap 'rm -rf "$WORK_TMP"' EXIT INT TERM

T0_MS="$(now_ms)"

python3 - "$VARDICTRUST_BIN" "$ROOT_DIR" "$ASSET_ROOT" "$WORK_TMP" <<'PY'
import concurrent.futures
import gzip
import os
import re
import shutil
import struct
import subprocess
import sys

vdr, root_dir, asset_root, tmp = sys.argv[1:5]
sys.path.insert(0, os.path.join(root_dir, "harness"))
import compare_vcf  # noqa: E402

for req_tool in ("bcftools", "tabix", "samtools", "perl"):
    if shutil.which(req_tool) is None:
        raise SystemExit(f"ERROR: required tool {req_tool!r} not found on PATH")

fa_b37 = os.path.join(asset_root, "VarDictJava", "testdata", "fastas", "GRCh37.fa")
fa_hg19 = os.path.join(asset_root, "VarDictJava", "testdata", "fastas", "hg19.fa")
fa_hg38 = os.path.join(asset_root, "VarDictJava", "testdata", "fastas", "hg38.fa")
fa_hs38dh = os.path.join(asset_root, "VarDictJava", "testdata", "fastas", "hs38DH.fa")

vj_bams_rel = ("VarDictJava", "src", "test", "resources", "com", "astrazeneca", "vardict", "integrationtests")
vj_bams = os.path.join(root_dir, *vj_bams_rel)
if not os.path.isdir(vj_bams):
    vj_bams = os.path.join(asset_root, *vj_bams_rel)

slices_dir = os.path.join(root_dir, "harness", "slices")

perl_valid = os.path.join(root_dir, "VarDictJava", "VarDict", "var2vcf_valid.pl")
if not os.path.isfile(perl_valid):
    perl_valid = os.path.join(asset_root, "VarDictJava", "VarDict", "var2vcf_valid.pl")

perl_paired = os.path.join(root_dir, "VarDictJava", "VarDict", "var2vcf_paired.pl")
if not os.path.isfile(perl_paired):
    perl_paired = os.path.join(asset_root, "VarDictJava", "VarDict", "var2vcf_paired.pl")

if not os.path.isdir(vj_bams):
    raise SystemExit(f"ERROR: missing portable BAM directory {vj_bams}")
for req in (fa_b37, fa_hg19, fa_hg38, fa_hs38dh, perl_valid, perl_paired):
    if not os.path.isfile(req):
        raise SystemExit(f"ERROR: missing portable asset {req}")

valid_alt_re = re.compile(r"^(?:[ACGTN]+|<[A-Za-z0-9_-]+>)(?:,(?:[ACGTN]+|<[A-Za-z0-9_-]+>))*$")


def verify_bgzf_header(gz_path: str) -> int:
    blocks = 0
    size = os.path.getsize(gz_path)
    with open(gz_path, "rb") as fh:
        pos = 0
        while pos < size:
            hdr = fh.read(12)
            if not hdr:
                break
            if len(hdr) < 12 or hdr[0:2] != b"\x1f\x8b" or (hdr[3] & 4) == 0:
                raise AssertionError(f"{gz_path}: invalid BGZF header at {pos}")
            xlen = struct.unpack("<H", hdr[10:12])[0]
            extra = fh.read(xlen)
            bsize = None
            idx = 0
            while idx + 4 <= len(extra):
                si1, si2, slen = extra[idx], extra[idx + 1], struct.unpack("<H", extra[idx + 2 : idx + 4])[0]
                if si1 == 66 and si2 == 67 and slen == 2:
                    bsize = struct.unpack("<H", extra[idx + 4 : idx + 6])[0] + 1
                    break
                idx += 4 + slen
            if bsize is None:
                raise AssertionError(f"{gz_path}: missing BGZF BC subfield")
            if bsize > 28:
                blocks += 1
            pos += bsize
            fh.seek(pos)
    return blocks


def audit_universal(label: str, vcf_text: str, ref_fa: str, kind: str, vcf_gz_path: str = "", query_reg: str = "", var_text: str = "") -> int:
    hdrs = []
    recs = []
    for line in vcf_text.splitlines():
        if line.startswith("#"):
            hdrs.append(line)
        elif line.strip():
            recs.append(line)
    if len(hdrs) < 4 or hdrs[0] != "##fileformat=VCFv4.2" or "##VarDictRustCompat=off" not in hdrs[:5]:
        raise AssertionError(f"{label}: missing VCFv4.2 / ##VarDictRustCompat=off header")
    if not recs:
        raise AssertionError(f"{label}: expected >= 1 VCF record")

    is_paired = kind == "somatic"
    seen_somatic = set()
    for line in recs:
        cols = line.split("\t")
        chrom, pos_s, ref_a, alt_a, qual_s, _flt, info_s, fmt_s = (
            cols[0], cols[1], cols[3], cols[4], cols[5], cols[6], cols[7], cols[8]
        )
        pos = int(pos_s)
        if pos < 1 or float(qual_s) < 0.0:
            raise AssertionError(f"{label}: invalid POS/QUAL in {line}")
        if alt_a != "." and not valid_alt_re.match(alt_a):
            raise AssertionError(f"{label}: invalid ALT {alt_a!r} in {line}")
        kv = dict(x.split("=", 1) for x in info_s.split(";") if "=" in x)
        is_sv = alt_a.startswith("<")
        if is_sv:
            end_p, svlen = int(kv["END"]), int(kv["SVLEN"])
            if len(ref_a) != 1 or end_p <= pos or (end_p - pos) != abs(svlen):
                raise AssertionError(f"{label}: malformed symbolic SV anchoring in {line}")
        else:
            if "END" in kv:
                raise AssertionError(f"{label}: non-SV leaked INFO/END in {line}")
            for single_alt in alt_a.split(","):
                if single_alt == "." or single_alt.startswith("<"):
                    continue
                if len(ref_a) != len(single_alt):
                    if len(ref_a) > 1 and len(single_alt) > 1 and ref_a[0] == single_alt[0]:
                        raise AssertionError(f"{label}: redundant 5' prefix in {line}")
                    if pos > 1 and len(ref_a) > 1 and len(single_alt) > 1 and ref_a[-1] == single_alt[-1]:
                        raise AssertionError(f"{label}: untrimmed 3' suffix in {line}")
                elif "," not in alt_a and int(kv.get("SHIFT3", "0")) > 0:
                    raise AssertionError(f"{label}: substitution leaked SHIFT3>0 in {line}")

        fmt_keys = fmt_s.split(":")
        n_alt = 0 if alt_a == "." else len(alt_a.split(","))
        any_vd, any_gt1 = False, False
        for smp in cols[9:]:
            vals = dict(zip(fmt_keys, smp.split(":")))
            dp, vd = int(vals["DP"]), int(vals["VD"])
            ad = [int(x) for x in vals["AD"].split(",")]
            ald = [int(x) for x in vals["ALD"].split(",")]
            if vd > dp or dp < sum(ad) or sum(ald) != vd:
                raise AssertionError(f"{label}: depth/AD/ALD accounting violation in {line}")
            if vals["GT"] == "1/0":
                raise AssertionError(f"{label}: non-canonical 1/0 GT in {line}")
            if len(ad) != (1 + n_alt):
                raise AssertionError(f"{label}: invalid AD arity in {line}")
            if vd > 0:
                any_vd = True
            if "1" in vals["GT"]:
                any_gt1 = True
        if alt_a != "." and any_vd and not any_gt1:
            raise AssertionError(f"{label}: unobserved ALT with GT=0/0 across all samples in {line}")
        if is_paired:
            for single_alt in alt_a.split(","):
                end_k = int(kv.get("END", pos)) if single_alt.startswith("<") else pos
                sk = (chrom, pos, end_k, ref_a, single_alt)
                if sk in seen_somatic:
                    raise AssertionError(f"{label}: duplicate paired somatic shadow record {sk}")
                seen_somatic.add(sk)

    if var_text:
        v_rows = [r for r in var_text.splitlines() if r.strip()]
        exp_cols = 61 if is_paired else (40 if kind == "amplicon" else 38)
        for row in v_rows:
            c = row.split("\t")
            if len(c) != exp_cols:
                raise AssertionError(f"{label}: expected {exp_cols} .var cols, got {len(c)}")
            offsets = (7, 27) if is_paired else (7,)
            for base in offsets:
                dp, vd = int(c[base]), int(c[base + 1])
                rf, rr, vf, vr = int(c[base + 2]), int(c[base + 3]), int(c[base + 4]), int(c[base + 5])
                if dp < rf + rr + vd or (vf + vr) != vd:
                    raise AssertionError(f"{label}: .var accounting violation in {row}")

    if vcf_gz_path:
        if verify_bgzf_header(vcf_gz_path) < 1:
            raise AssertionError(f"{label}: empty BGZF stream {vcf_gz_path}")
        split_out = subprocess.check_output(
            ["bcftools", "norm", "-m", "-any", vcf_gz_path],
            stderr=subprocess.DEVNULL,
            text=True,
            timeout=120,
        )
        norm_out = subprocess.check_output(
            ["bcftools", "norm", "--check-ref", "e", "-f", ref_fa, "-"],
            input=split_out,
            stderr=subprocess.DEVNULL,
            text=True,
            timeout=120,
        )
        split_recs = [l.split("\t", 5)[:5] for l in split_out.splitlines() if l and not l.startswith("#")]
        norm_recs = [l.split("\t", 5)[:5] for l in norm_out.splitlines() if l and not l.startswith("#")]
        c_split = [r for r in split_recs if r[4] != "." and not r[4].startswith("<")]
        c_norm = [r for r in norm_recs if r[4] != "." and not r[4].startswith("<")]
        if len(split_recs) != len(norm_recs) or len(c_split) != len(c_norm):
            raise AssertionError(
                f"{label}: bcftools norm record count mismatch ({len(split_recs)}->{len(norm_recs)}, concrete {len(c_split)}->{len(c_norm)})"
            )
        for s_cols, n_cols in zip(c_split, c_norm):
            s_chr, s_pos_s, _, s_ref, s_alt = s_cols
            n_chr, n_pos_s, _, n_ref, n_alt = n_cols
            t_pos = int(s_pos_s)
            t_ref, t_alt = s_ref, s_alt
            while len(t_ref) > 1 and len(t_alt) > 1 and t_ref[-1] == t_alt[-1]:
                t_ref, t_alt = t_ref[:-1], t_alt[:-1]
            while len(t_ref) > 1 and len(t_alt) > 1 and t_ref[0] == t_alt[0]:
                t_ref, t_alt = t_ref[1:], t_alt[1:]
                t_pos += 1
            if (s_chr, t_pos, t_ref, t_alt) != (n_chr, int(n_pos_s), n_ref, n_alt):
                raise AssertionError(
                    f"{label}: non-canonical allele rewritten by bcftools norm: {(s_chr, t_pos, t_ref, t_alt)} != {(n_chr, int(n_pos_s), n_ref, n_alt)}"
                )
        if query_reg:
            bcf_q = subprocess.check_output(
                ["bcftools", "view", "-H", "-r", query_reg, vcf_gz_path],
                text=True,
                timeout=120,
            )
            if not bcf_q.strip():
                raise AssertionError(f"{label}: empty bcftools view -r {query_reg} on {vcf_gz_path}")
            tbx_q = subprocess.check_output(
                ["tabix", vcf_gz_path, query_reg],
                text=True,
                timeout=120,
            )
            if [l.split("\t")[:5] for l in bcf_q.splitlines() if l.strip()] != [l.split("\t")[:5] for l in tbx_q.splitlines() if l.strip()]:
                raise AssertionError(f"{label}: tabix vs bcftools view -r mismatch on {query_reg}")
    return len(recs)


def gate1_compat_off_and_presets() -> str:
    na_bam = os.path.join(vj_bams, "NA24385_chr_1_246981051_246981611.bam")
    hcc_pair = f"{os.path.join(vj_bams, 'HCC2218_BL_realigins.bam')}|{os.path.join(vj_bams, 'HCC2218_C_realigins.bam')}"
    sv_bam = os.path.join(vj_bams, "Colo829.chr10.DELDiscReverse.bam")

    # 1a) --compat=off --vcf-all-alleles --output .vcf.gz --csi + .var across Simple, Somatic, SV
    total = 0
    for tag, kind, fa, bam, reg, extra in [
        ("g1_simple", "simple", fa_b37, na_bam, "1:246981051-246981611", ["-f", "0.01"]),
        ("g1_somatic", "somatic", fa_hg19, hcc_pair, "chrX:7811670-7811868", ["-f", "0.001"]),
        ("g1_sv", "simple", fa_hg19, sv_bam, "chr10:89700380-89713976", ["-f", "0.001"]),
    ]:
        gz = os.path.join(tmp, f"{tag}.vcf.gz")
        subprocess.run(
            [vdr, "-G", fa, "-b", bam, "-R", reg, "--compat=off", "--vcf-all-alleles", "--fisher", "--output", gz, "--csi", "-th", "2", *extra],
            check=True,
            timeout=120,
        )
        if not os.path.isfile(gz + ".csi"):
            raise AssertionError(f"{tag}: missing .csi sidecar")
        var_txt = subprocess.check_output(
            [vdr, "-G", fa, "-b", bam, "-R", reg, "--compat=off", "--vcf-all-alleles", "--fisher", "--emit-intermediate", "-th", "2", *extra],
            text=True,
            timeout=120,
        )
        with gzip.open(gz, "rt", encoding="utf-8") as fh:
            vcf_txt = fh.read()
        total += audit_universal(tag, vcf_txt, fa, kind, gz, reg, var_txt)

    # 1b) --preset {wgs, wes, amplicon, ctdna} with --output .vcf.gz --tbi / --csi
    bed_wgs = os.path.join(tmp, "p_wgs.bed")
    with open(bed_wgs, "w", encoding="utf-8") as fh:
        fh.write("1\t246981050\t246981611\n")
    bed_wes4 = os.path.join(tmp, "p_wes.bed")
    with open(bed_wes4, "w", encoding="utf-8") as fh:
        fh.write("chrX\t7811669\t7811868\t3'_UTR\n")
    bed_wes8 = os.path.join(tmp, "p_wes8.bed")
    with open(bed_wes8, "w", encoding="utf-8") as fh:
        fh.write("chrX\t7811669\t7811868\tGENE1\t0\t+\t7811680\t7811860\n")

    preset_checks = [
        ("wgs", "simple", fa_b37, na_bam, bed_wgs, "1:246981051-246981611", "--tbi", []),
        ("wes", "somatic", fa_hg19, hcc_pair, bed_wes4, "chrX:7811670-7811868", "--csi", ["--fisher"]),
        ("amplicon", "amplicon", fa_hg19, os.path.join(vj_bams, "Colo829-18_S3-sort.bam"), os.path.join(slices_dir, "s6.bed"), "chr1:933867-934466", "--csi", []),
        ("ctdna", "simple", fa_b37, na_bam, bed_wgs, "1:246981051-246981611", "--tbi", []),
    ]
    for pname, pkind, fa, bam, bed_p, qreg, idx_flag, pextra in preset_checks:
        pgz = os.path.join(tmp, f"preset_{pname}.vcf.gz")
        idx_ext = ".tbi" if idx_flag == "--tbi" else ".csi"
        subprocess.run(
            [vdr, "--preset", pname, *pextra, "-G", fa, "-b", bam, "--output", pgz, idx_flag, "-th", "2", bed_p],
            check=True,
            timeout=120,
        )
        if not os.path.isfile(pgz + idx_ext):
            raise AssertionError(f"preset_{pname}: missing {idx_ext} sidecar")
        with gzip.open(pgz, "rt", encoding="utf-8") as fh:
            ptxt = fh.read()
        if f"##VarDictRustPreset={pname}\n" not in ptxt:
            raise AssertionError(f"preset {pname} missing ##VarDictRustPreset={pname}")
        total += audit_universal(f"preset_{pname}", ptxt, fa, pkind, pgz, qreg)

    # Verify UCSC 8-col BED8 non-hijack under --preset wes
    wes8_txt = subprocess.check_output(
        [vdr, "--preset", "wes", "--fisher", "-G", fa_hg19, "-b", hcc_pair, "-th", "2", bed_wes8],
        text=True,
        timeout=120,
    )
    if "AMPBIAS" in wes8_txt or "GDAMP=" in wes8_txt:
        raise AssertionError("UCSC 8-col BED8 hijacked --preset wes into PCR amplicon schema")
    return f"Gate 1 PASS ({total} audited --compat=off & --preset records, .csi/.tbi + bcftools norm green)"


def gate2_tiled_vs_untiled() -> str:
    na_bam = os.path.join(vj_bams, "NA24385_chr_1_246981051_246981611.bam")

    cases = [
        ("simple", fa_b37, na_bam, ["-R", "1:246981051-246981611", "-f", "0.01"], "200", 10),
        ("somatic_s4", fa_hg19, f"{vj_bams}/HCC2218_BL_realigins.bam|{vj_bams}/HCC2218_C_realigins.bam",
         ["-c", "1", "-S", "2", "-E", "3", "-g", "4", "-z", "-f", "0.001", os.path.join(slices_dir, "s4.bed")], "200", 14),
        ("somatic_s5", fa_hg19, f"{vj_bams}/Zero_in_DUPSV_TUMOUR.bam|{vj_bams}/Zero_in_DUPSV_NORMAL.bam",
         ["-c", "1", "-S", "2", "-E", "3", "-g", "4", "-z", "-f", "0.001", "-X", "3", os.path.join(slices_dir, "s5.bed")], "200", 60),
        ("amplicon_a", fa_hg19, f"{vj_bams}/Colo829-18_S3-sort.bam",
         ["-a", "10:0.95", "-z", os.path.join(slices_dir, "s6.bed")], "200", 4),
        ("amplicon_sp", fa_b37, f"{vj_bams}/MRE11a.bam",
         ["-a", "10:0.95", "-sp", "-p", "-z", os.path.join(slices_dir, "s7.bed")], "200", 90),
        ("pileup_200", fa_b37, na_bam, ["-R", "1:246981051-246981611", "-p"], "200", 500),
        ("pileup_300", fa_b37, f"{vj_bams}/N13532.TwoDels1InsStringBoundError.bam",
         ["-R", "11:1279500-1281000", "-p"], "300", 900),
        ("splice_i", fa_hg38, f"{vj_bams}/Tcov_zero_SV_DEL_candidate3.bam",
         ["-R", "chr8:90534611-90620070", "-f", "0.1", "-i", "--emit-intermediate"], "200", 16),
        ("sv_del", fa_hg19, f"{vj_bams}/Colo829.chr10.DELDiscReverse.bam",
         ["-R", "chr10:89712340-89713976", "-f", "0.001"], "200", 1),
        ("sv_inv", fa_hg19, f"{vj_bams}/Colo829.chr20.INVSV.bam",
         ["-R", "chr20:46842755-46845755", "-f", "0.001"], "200", 10),
    ]
    checked = 0
    for name, fa, bam, args, ts, min_rows in cases:
        for compat in ("--compat=vardictjava", "--compat=off"):
            u = subprocess.check_output([vdr, "-G", fa, "-b", bam, compat, "--tile-size", "50000000", "-th", "1", *args], text=True, timeout=120)
            t = subprocess.check_output([vdr, "-G", fa, "-b", bam, compat, "--tile-size", ts, "-th", "2", *args], text=True, timeout=120)
            rows = [l for l in u.splitlines() if l.strip() and not l.startswith("#")]
            if len(rows) < min_rows:
                raise AssertionError(f"Gate 2 {name} ({compat}): expected >= {min_rows} rows, got {len(rows)}")
            if u != t:
                raise AssertionError(f"Gate 2 {name} ({compat}): Tiled (--tile-size {ts}) != Untiled")
            checked += 1

    # Spanning-D & far-breakpoint SV split-tile seam parity (Colo829 chr1:756320 35D & NA12878.INVSimpleSV 20:1961170)
    seam_cases = [
        ("colo829_35d_seam", fa_hg19, f"{vj_bams}/Colo829-18_S3-sort.bam|{vj_bams}/Colo829-19_S4-sort.bam",
         ["-R", "chr1:755917-756517", "-f", "0.01"]),
        ("na12878_inv_seam", fa_hs38dh, f"{vj_bams}/NA12878.chrom20.INVSimpleSV.bam",
         ["-R", "20:1961050-1961650", "-f", "0.01"]),
    ]
    for sname, sfa, sbam, sargs in seam_cases:
        for compat in ("--compat=vardictjava", "--compat=off"):
            su = subprocess.check_output([vdr, "-G", sfa, "-b", sbam, compat, "--tile-size", "50000000", "-th", "1", *sargs], text=True, timeout=120)
            st = subprocess.check_output([vdr, "-G", sfa, "-b", sbam, compat, "--tile-size", "200", "-th", "2", *sargs], text=True, timeout=120)
            if su != st:
                raise AssertionError(f"Gate 2 {sname} ({compat}): spanning-deletion / SV tile-seam Tiled != Untiled")
            checked += 1
    return f"Gate 2 PASS ({checked}/24 modal Tiled==Untiled checks byte-identical)"


def make_scoped_bam_and_cram(src_bam: str, ref_fa: str, chrom: str, prefix: str) -> tuple[str, str]:
    fai_len = None
    with open(ref_fa + ".fai", "r", encoding="utf-8") as fh:
        for line in fh:
            p = line.split("\t")
            if p[0] == chrom:
                fai_len = int(p[1])
                break
    if fai_len is None:
        raise AssertionError(f"{chrom} missing from {ref_fa}.fai")
    sam_raw = subprocess.check_output(["samtools", "view", "-h", src_bam, chrom], text=True, timeout=120)
    sam_clean = os.path.join(tmp, f"{prefix}.sam")
    with open(sam_clean, "w", encoding="utf-8") as out:
        out.write(f"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:{chrom}\tLN:{fai_len}\n")
        for ln in sam_raw.splitlines():
            if ln.startswith("@RG") or ln.startswith("@PG"):
                out.write(ln + "\n")
            elif not ln.startswith("@") and ln.strip():
                cols = ln.split("\t")
                if cols[6] not in ("*", "=", chrom):
                    cols[6], cols[7], cols[8] = "*", "0", "0"
                out.write("\t".join(cols) + "\n")
    out_bam = os.path.join(tmp, f"{prefix}.bam")
    out_cram = os.path.join(tmp, f"{prefix}.cram")
    subprocess.run(["samtools", "view", "-b", "-o", out_bam, sam_clean], check=True, timeout=120)
    subprocess.run(["samtools", "index", out_bam], check=True, timeout=120)
    subprocess.run(["samtools", "view", "-C", "-T", ref_fa, "-o", out_cram, out_bam], check=True, timeout=120)
    subprocess.run(["samtools", "index", out_cram], check=True, timeout=120)
    return out_bam, out_cram


def gate3_portable_bam_vs_cram() -> str:
    na_bam_src = os.path.join(vj_bams, "NA24385_chr_1_246981051_246981611.bam")
    na_bam, na_cram = make_scoped_bam_and_cram(na_bam_src, fa_b37, "1", "na24385_cram")

    hcc_t_src = os.path.join(vj_bams, "HCC2218_BL_realigins.bam")
    hcc_n_src = os.path.join(vj_bams, "HCC2218_C_realigins.bam")
    hcc_t_bam, hcc_t_cram = make_scoped_bam_and_cram(hcc_t_src, fa_hg19, "chrX", "hcc2218_t")
    hcc_n_bam, hcc_n_cram = make_scoped_bam_and_cram(hcc_n_src, fa_hg19, "chrX", "hcc2218_n")

    for compat in ("--compat=vardictjava", "--compat=off"):
        b_out = subprocess.check_output(
            [vdr, "-G", fa_b37, "-b", na_bam, "-N", "NA24385", "-R", "1:246981051-246981611", "-f", "0.01", compat, "-th", "2"],
            text=True,
            timeout=120,
        )
        c_out = subprocess.check_output(
            [vdr, "-G", fa_b37, "-b", na_cram, "-N", "NA24385", "-R", "1:246981051-246981611", "-f", "0.01", compat, "-th", "2"],
            text=True,
            timeout=120,
        )
        b_rows = [l for l in b_out.splitlines() if not l.startswith("#")]
        c_rows = [l for l in c_out.splitlines() if not l.startswith("#")]
        if len(b_rows) < 10 or b_rows != c_rows:
            raise AssertionError(f"Gate 3 single-sample BAM != CRAM under {compat}")

        pb_out = subprocess.check_output(
            [vdr, "-G", fa_hg19, "-b", f"{hcc_t_bam}|{hcc_n_bam}", "-N", "TUMOR|NORMAL",
             "-c", "1", "-S", "2", "-E", "3", "-g", "4", "-z", "-f", "0.001", compat, "-th", "2",
             os.path.join(slices_dir, "s4.bed")],
            text=True,
            timeout=120,
        )
        pc_out = subprocess.check_output(
            [vdr, "-G", fa_hg19, "-b", f"{hcc_t_cram}|{hcc_n_cram}", "-N", "TUMOR|NORMAL",
             "-c", "1", "-S", "2", "-E", "3", "-g", "4", "-z", "-f", "0.001", compat, "-th", "2",
             os.path.join(slices_dir, "s4.bed")],
            text=True,
            timeout=120,
        )
        pb_rows = [l for l in pb_out.splitlines() if not l.startswith("#")]
        pc_rows = [l for l in pc_out.splitlines() if not l.startswith("#")]
        if len(pb_rows) < 14 or pb_rows != pc_rows:
            raise AssertionError(f"Gate 3 paired somatic BAM != CRAM under {compat}")
    return "Gate 3 PASS (4/4 single & paired live BAM==CRAM checks byte-identical across vardictjava & off)"


def gate4_perl_all_alleles_parity() -> str:
    na_bam = os.path.join(vj_bams, "NA24385_chr_1_246981051_246981611.bam")
    colo_pair = f"{vj_bams}/Colo829-18_S3-sort.bam|{vj_bams}/Colo829-19_S4-sort.bam"
    hcc_pair = f"{vj_bams}/HCC2218_BL_realigins.bam|{vj_bams}/HCC2218_C_realigins.bam"

    perl_cases = [
        (
            "perl_valid_simple",
            [vdr, "-G", fa_b37, "-b", na_bam, "-R", "1:246981051-246981611", "-N", "NA24385", "-f", "0.01", "--compat=vardictjava", "--fisher"],
            ["perl", perl_valid, "-N", "NA24385", "-A", "-E", "-f", "0.01"],
            ["--vcf-all-alleles", "--no-end"],
        ),
        (
            "perl_valid_amplicon",
            [vdr, "-G", fa_hg19, "-b", f"{vj_bams}/Colo829-18_S3-sort.bam", "-N", "AMP_SMP", "-a", "10:0.95", "-z", "-f", "0.01",
             "--compat=vardictjava", "--fisher", os.path.join(slices_dir, "s6.bed")],
            ["perl", perl_valid, "-N", "AMP_SMP", "-A", "-E", "-f", "0.01"],
            ["--vcf-all-alleles", "--no-end"],
        ),
        (
            "perl_paired_colo829",
            [vdr, "-G", fa_hg19, "-b", colo_pair, "-N", "TUMOR|NORMAL", "-c", "1", "-S", "2", "-E", "3", "-g", "4", "-z",
             "-f", "0.01", "--compat=vardictjava", "--fisher", os.path.join(slices_dir, "s3.bed")],
            ["perl", perl_paired, "-N", "TUMOR|NORMAL", "-A", "-f", "0.01"],
            ["--vcf-all-alleles"],
        ),
        (
            "perl_paired_hcc2218",
            [vdr, "-G", fa_hg19, "-b", hcc_pair, "-N", "TUMOR|NORMAL", "-c", "1", "-S", "2", "-E", "3", "-g", "4", "-z",
             "-f", "0.001", "--compat=vardictjava", "--fisher", os.path.join(slices_dir, "s4.bed")],
            ["perl", perl_paired, "-N", "TUMOR|NORMAL", "-A", "-f", "0.001"],
            ["--vcf-all-alleles"],
        ),
    ]
    total_compared = 0
    for cname, base_cmd, perl_cmd, rust_vcf_flags in perl_cases:
        var_out = subprocess.check_output([*base_cmd, "--emit-intermediate"], text=True, timeout=120)
        perl_vcf = subprocess.check_output(perl_cmd, input=var_out, text=True, timeout=120)
        rust_vcf = subprocess.check_output([*base_cmd, *rust_vcf_flags], text=True, timeout=120)
        p_path = os.path.join(tmp, f"{cname}_perl.vcf")
        r_path = os.path.join(tmp, f"{cname}_rust.vcf")
        with open(p_path, "w", encoding="utf-8") as fp:
            fp.write(perl_vcf)
        with open(r_path, "w", encoding="utf-8") as fr:
            fr.write(rust_vcf)
        _, lj = compare_vcf.load_vcf(p_path)
        _, lr = compare_vcf.load_vcf(r_path)
        gj, gr = compare_vcf.group(lj), compare_vcf.group(lr)
        errs = []
        if set(gj) != set(gr) or any(len(gj[k]) != len(gr[k]) for k in gj) or not lj:
            raise AssertionError(f"Gate 4 {cname} key-set mismatch: perl={len(lj)} vs rust={len(lr)}")
        for k in sorted(gj):
            for jr, rr in zip(gj[k], gr[k]):
                compare_vcf.compare_pair(k, jr, rr, 1e-6, errs, [])
        if errs:
            raise AssertionError(f"Gate 4 {cname} --vcf-all-alleles vs Perl -A mismatch: {errs[:5]}")
        total_compared += len(lj)
    return f"Gate 4 PASS (4/4 single, amplicon & paired --vcf-all-alleles == var2vcf_*.pl -A, {total_compared} alleles)"


max_workers = min(4, max(1, (os.cpu_count() or 2) // 2))
with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as pool:
    futs = [
        pool.submit(gate1_compat_off_and_presets),
        pool.submit(gate2_tiled_vs_untiled),
        pool.submit(gate3_portable_bam_vs_cram),
        pool.submit(gate4_perl_all_alleles_parity),
    ]
    for fut in futs:
        print("  " + fut.result())
PY

ELAPSED_MS=$(( $(now_ms) - T0_MS ))
echo "PASS run_portable_invariants: all 4 portable invariant gates green (${ELAPSED_MS}ms)"
