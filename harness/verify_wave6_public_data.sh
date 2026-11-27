#!/usr/bin/env bash
# =============================================================================
# harness/verify_wave6_public_data.sh
# Local $0.00 public truth-set and invariant verification suite for Wave-6
# CompatProfile accuracy switches 33..38 (bringing --compat=off to 38 switches)
# and zero-diff contiguous-slice realignment fast-paths.
#
# Cross-Ecosystem Provenance & Citation Registry (owner/repo#N + DOIs):
# - Switch 33 (sclip_majority_count_first):
#   VariationUtils.java:90, palakpsheth/VarDictRust#160 (Item 4).
# - Switch 34 (mate_overlap_both_strands):
#   Clinical-Genomics/BALSAMIC#1396, AstraZeneca-NGS/VarDictJava#239,
#   palakpsheth/VarDictRust#160 (Item 8).
# - Switch 35 (canonical_vcf_complex_left_align):
#   AstraZeneca-NGS/VarDictJava#264, AstraZeneca-NGS/VarDict#134,
#   mskcc/vcf2maf#279, samtools/bcftools#1379, palakpsheth/VarDictRust#160 (Item 5);
#   Tan et al., Bioinformatics 2015 (DOI: 10.1093/bioinformatics/btv112).
# - Switch 36 (somatic_reconcile_split_alleles):
#   clintval/vartovcf#191, clintval/vartovcf#195, palakpsheth/VarDictRust#160 (Item 11).
# - Switch 37 (adj_snv_full_sclip_match):
#   clintval/vartovcf#181, clintval/vartovcf#186, AstraZeneca-NGS/VarDict#174,
#   palakpsheth/VarDictRust#160 (Item 10A).
# - Switch 38 (indel_only_shift3):
#   AstraZeneca-NGS/VarDictJava#263, clintval/vartovcf#117,
#   palakpsheth/VarDictRust#160 (Item 12B);
#   Cleary et al., bioRxiv 2015 (DOI: 10.1101/023754, rtg vcfeval);
#   Fang et al., Nat. Biotechnol. 2021 (DOI: 10.1038/s41587-021-00993-6, SEQC2).
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
source "$SCRIPT_DIR/lib/portable.sh"

ensure_g4_tools_path
for _conda_env in "${CONDA_ROOT:-${HOME:-}/miniconda3}/envs/g4-tools/bin" \
                  "${CONDA_ROOT:-${HOME:-}/miniconda3}/envs/vardict-oracle/bin"; do
    if [ -d "$_conda_env" ]; then
        case ":$PATH:" in
            *":$_conda_env:"*) ;;
            *) export PATH="$_conda_env:$PATH" ;;
        esac
    fi
done
unset _conda_env

for _req_tool in samtools bcftools tabix bgzip rtg python3 gzip; do
    if ! command -v "$_req_tool" >/dev/null 2>&1; then
        echo "ERROR: required tool '$_req_tool' not found on PATH ($PATH)" >&2
        exit 1
    fi
done
unset _req_tool

VARDICTRUST_BIN="${VARDICTRUST_BIN:-${VARDICT_BIN:-$ROOT_DIR/target/release/vardictrust}}"
if [ ! -x "$VARDICTRUST_BIN" ]; then
    echo "ERROR: VARDICTRUST_BIN is not executable at '$VARDICTRUST_BIN' (run 'cargo build --release' first)" >&2
    exit 1
fi

DATA_ROOT="$(resolve_data_root "$ROOT_DIR")"
G4_DIR="${DATA_DIR:-$DATA_ROOT/scratch/g4}"
WAVE4_DIR="${WAVE4_DIR:-$G4_DIR/wave4}"

GRCH38_REF="$G4_DIR/ref/d1vd1/GRCh38.d1.vd1.fa"
GIAB_REF="$G4_DIR/ref/giabv3/GRCh38_GIABv3.fa"
if [ ! -f "$GIAB_REF" ] && [ -f "$G4_DIR/ref/giabv3/GRCh38_GIABv3.fasta" ]; then
    GIAB_REF="$G4_DIR/ref/giabv3/GRCh38_GIABv3.fasta"
fi

for _req_file in "$GRCH38_REF" "${GRCH38_REF}.fai" "$GIAB_REF" "${GIAB_REF}.fai" \
                 "$G4_DIR/hg002/HG002.chr20.bam" \
                 "$G4_DIR/hg002/eval.bed" "$G4_DIR/hg002/truth/truth.vcf.gz" \
                 "$G4_DIR/seqc2/HCC1395.bam" "$G4_DIR/seqc2/HCC1395BL.bam" \
                 "$G4_DIR/seqc2/slices.bed" "$G4_DIR/seqc2/eval.bed" \
                 "$G4_DIR/seqc2/truth/truth.vcf.gz"; do
    if [ ! -s "$_req_file" ]; then
        echo "ERROR: missing Wave-6 public dataset prerequisite: '$_req_file'" >&2
        exit 1
    fi
done
unset _req_file

WORK_TMP="$(mktemp -d "${TMPDIR:-/tmp}/vdr_wave6_verify.XXXXXX")"
trap 'rm -rf "$WORK_TMP"' EXIT

T0_MS="$(now_ms)"
echo "================================================================================"
echo " Wave-6 5-Gate Public Truth-Set & Invariant Suite (38 CompatProfile Switches)"
echo " Binary: $VARDICTRUST_BIN"
echo " Data:   $G4_DIR"
echo "================================================================================"

# =============================================================================
# Gates 1..4: Execute & Audit Cached Public Cohorts (HG002, SEQC2, HG008, Wave-4)
#   1. bcftools norm --check-ref e == 0 mismatches + 0 unanchored/untrimmed alleles
#   2. 0 SNV / equal-length MNV records carry SHIFT3 > 0 under --compat=off
#      (vs > 0 in --compat=vardictjava)
#   3. -u --compat=off preserves bidirectional ALD balance without VD > DP or
#      ALD[0] + ALD[1] != VD
#   4. Paired somatic --compat=off has 0 duplicate (CHROM, POS, REF, ALT)
#      StrongSomatic + StrongLOH shadow collisions
# =============================================================================
echo "[1/5] Running --compat=vardictjava vs --compat=off across HG002, SEQC2, HG008, & Wave-4..."
nice -n 19 python3 - "$VARDICTRUST_BIN" "$GRCH38_REF" "$GIAB_REF" "$G4_DIR" "$WAVE4_DIR" "$WORK_TMP" <<'PY'
import gzip
import os
import re
import subprocess
import sys

vdr, grch38_ref, giab_ref, g4_dir, w4_dir, tmp = sys.argv[1:7]
try:
    os.nice(19)
except OSError:
    pass

valid_alt_re = re.compile(r"^(?:[ACGTN]+|<[A-Za-z0-9_-]+>)(?:,(?:[ACGTN]+|<[A-Za-z0-9_-]+>))*$")

hg002_bed = os.path.join(tmp, "hg002_chr20_1_5mb.bed")
with open(hg002_bed, "w", encoding="utf-8") as fh:
    fh.write("chr20\t10000000\t11500000\thg002_c20\n")

seqc2_slices = os.path.join(g4_dir, "seqc2", "slices.bed")
cohorts = [
    (
        "GIAB_HG002_chr20",
        False,
        giab_ref,
        os.path.join(g4_dir, "hg002", "HG002.chr20.bam"),
        ["-c", "1", "-S", "2", "-E", "3", "-g", "4", "-z", "-f", "0.01", "-th", "4", hg002_bed],
    ),
    (
        "SEQC2_HCC1395_TN",
        True,
        grch38_ref,
        f"{g4_dir}/seqc2/HCC1395.bam|{g4_dir}/seqc2/HCC1395BL.bam",
        ["-c", "1", "-S", "2", "-E", "3", "-g", "4", "-z", "-f", "0.01", "--fisher", "-th", "4", seqc2_slices],
    ),
]

hg008_t = os.path.join(g4_dir, "hg008", "HG008-T.slice.bam")
hg008_n = os.path.join(g4_dir, "hg008", "HG008-N.slice.bam")
hg008_bed = os.path.join(g4_dir, "hg008", "slices.bed")
if os.path.isfile(hg008_t) and os.path.isfile(hg008_n) and os.path.isfile(hg008_bed):
    cohorts.append(
        (
            "GIAB_HG008_TN",
            True,
            giab_ref,
            f"{hg008_t}|{hg008_n}",
            ["-c", "1", "-S", "2", "-E", "3", "-z", "-f", "0.01", "--fisher", "-th", "4", hg008_bed],
        )
    )

if os.path.isdir(w4_dir):
    cmrg_bam = os.path.join(w4_dir, "hg002_cmrg_mosaic", "HG002.cmrg_mosaic.slice.bam")
    cmrg_bed = os.path.join(w4_dir, "hg002_cmrg_mosaic", "hg002_cmrg_mosaic.bed")
    if os.path.isfile(cmrg_bam) and os.path.isfile(cmrg_bed):
        cohorts.append(
            (
                "Wave4_HG002_CMRG",
                False,
                giab_ref,
                cmrg_bam,
                ["-c", "1", "-S", "2", "-E", "3", "-g", "4", "-z", "-f", "0.01", "-th", "4", cmrg_bed],
            )
        )
    qgn_ac5 = os.path.join(w4_dir, "bam", "SampleAC5_QGN1.bam")
    qgn_b = os.path.join(w4_dir, "bam", "SampleB_QGN1.bam")
    qgn_bed = os.path.join(w4_dir, "beds", "qiaseq_spe_grch38.amp8.bed")
    if os.path.isfile(qgn_ac5) and os.path.isfile(qgn_b) and os.path.isfile(qgn_bed):
        cohorts.append(
            (
                "Wave4_QIAseq_ctDNA_TN",
                True,
                grch38_ref,
                f"{qgn_ac5}|{qgn_b}",
                ["-a", "10:0.95", "-sp", "-z", "-f", "0.0025", "--fisher", "-th", "4", qgn_bed],
            )
        )
    ilm_a = os.path.join(w4_dir, "bam", "SampleA_ILM1.bam")
    ilm_b = os.path.join(w4_dir, "bam", "SampleB_ILM1.bam")
    ilm_bed = os.path.join(w4_dir, "beds", "ilm_twoprimer_grch38.amp8.bed")
    if os.path.isfile(ilm_a) and os.path.isfile(ilm_b) and os.path.isfile(ilm_bed):
        cohorts.append(
            (
                "Wave4_Illumina_Amp_TN",
                True,
                grch38_ref,
                f"{ilm_a}|{ilm_b}",
                ["-a", "10:0.95", "-z", "-f", "0.01", "--fisher", "-th", "4", ilm_bed],
            )
        )

total_off_rows = 0
vj_subst_shift3_pos = 0
off_subst_shift3_pos = 0
off_indel_shift3_pos = 0
somatic_cohorts_checked = 0

for name, is_paired, ref_fa, bam_arg, extra_args in cohorts:
    off_gz = os.path.join(tmp, f"{name}_off.vcf.gz")
    subprocess.run(
        [vdr, "-G", ref_fa, "-b", bam_arg, "--compat=off", "--vcf-all-alleles", "--output", off_gz, "--tbi", *extra_args],
        check=True,
    )
    # Invariant 1a: compare concrete split+unpadded (CHROM, POS, REF, ALT) alleles before vs after
    # `bcftools norm --check-ref e -f ref_fa` so right-shifted indels (e.g. non-canonical `AA>A`)
    # that `bcftools norm` rewrites instead of rejecting fail immediately (allowing ONLY the
    # documented VCFv4.2 stop-before-'N' exception where `bcftools norm -f` rolls onto an 'N' base).
    split_vcf = subprocess.check_output(
        ["bcftools", "norm", "-m", "-any", off_gz],
        stderr=subprocess.DEVNULL,
        text=True,
    )
    norm_vcf = subprocess.check_output(
        ["bcftools", "norm", "--check-ref", "e", "-f", ref_fa, "-"],
        input=split_vcf,
        stderr=subprocess.DEVNULL,
        text=True,
    )
    all_split_recs = [
        l.split("\t", 5)[:5]
        for l in split_vcf.splitlines()
        if l and not l.startswith("#")
    ]
    all_norm_recs = [
        l.split("\t", 5)[:5]
        for l in norm_vcf.splitlines()
        if l and not l.startswith("#")
    ]
    if len(all_split_recs) != len(all_norm_recs):
        raise AssertionError(
            f"{name}: bcftools norm record count mismatch: split={len(all_split_recs)} vs norm={len(all_norm_recs)}"
        )
    # Filter symbolic `<...>` SVs before pairing concrete alleles (`bcftools norm -f` shifts symbolic `<DUP>`
    # across homologous endpoints and re-sorts across adjacent SNVs, while Switch 35 targets concrete alleles).
    split_recs = [r for r in all_split_recs if r[4] != "." and not r[4].startswith("<")]
    norm_recs = [r for r in all_norm_recs if r[4] != "." and not r[4].startswith("<")]
    if len(split_recs) != len(norm_recs):
        raise AssertionError(
            f"{name}: concrete bcftools norm count mismatch: split={len(split_recs)} vs norm={len(norm_recs)}"
        )
    for s_cols, n_cols in zip(split_recs, norm_recs):
        s_chr, s_pos_s, _, s_ref, s_alt = s_cols
        n_chr, n_pos_s, _, n_ref, n_alt = n_cols
        t_pos = int(s_pos_s)
        t_ref, t_alt = s_ref, s_alt
        # Trim shared multiallelic 3' suffix then 5' prefix down to a 1-bp anchor at existing `t_pos`
        # (NO reference left-rolling!).
        while len(t_ref) > 1 and len(t_alt) > 1 and t_ref[-1] == t_alt[-1]:
            t_ref = t_ref[:-1]
            t_alt = t_alt[:-1]
        while len(t_ref) > 1 and len(t_alt) > 1 and t_ref[0] == t_alt[0]:
            t_ref = t_ref[1:]
            t_alt = t_alt[1:]
            t_pos += 1
        n_pos = int(n_pos_s)
        if (s_chr, t_pos, t_ref, t_alt) != (n_chr, n_pos, n_ref, n_alt):
            # Allow ONLY VCFv4.2 stop-before-'N' where VarDictRust halted immediately 3' of 'N'
            # while `bcftools norm -f` rolled onto an 'N' base in `n_ref`/`n_alt`.
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
                prev_base = "".join(
                    ln.strip() for ln in faidx_out.splitlines() if not ln.startswith(">")
                ).upper()
                if prev_base == "N":
                    span_end = max(t_pos + len(t_ref) - 1, n_pos + len(n_ref) - 1)
                    win_out = subprocess.check_output(
                        ["samtools", "faidx", ref_fa, f"{s_chr}:{n_pos}-{span_end}"],
                        text=True,
                    )
                    win_seq = "".join(
                        ln.strip() for ln in win_out.splitlines() if not ln.startswith(">")
                    ).upper()
                    off_t = t_pos - n_pos
                    mut_vdr = win_seq[:off_t] + t_alt.upper() + win_seq[off_t + len(t_ref) :]
                    mut_bcf = n_alt.upper() + win_seq[len(n_ref) :]
                    is_stop_before_n = mut_vdr == mut_bcf
            if not is_stop_before_n:
                raise AssertionError(
                    f"{name}: non-canonical allele rewritten by bcftools norm -f: "
                    f"emitted=({s_chr}, {t_pos}, {t_ref!r}, {t_alt!r}) vs "
                    f"bcftools=({n_chr}, {n_pos}, {n_ref!r}, {n_alt!r})"
                )

    # Run --compat=vardictjava on non-paired-amplicon cohorts (paired amplicon is unsupported in v1.8.3 VCF mode)
    is_paired_amp = is_paired and "-a" in extra_args
    if not is_paired_amp:
        vj_txt = subprocess.check_output(
            [vdr, "-G", ref_fa, "-b", bam_arg, "--compat=vardictjava", "--vcf-all-alleles", "-v", *extra_args],
            text=True,
        )
        for line in vj_txt.splitlines():
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split("\t")
            ref_a, alt_a = cols[3], cols[4]
            if not alt_a.startswith("<") and alt_a != "." and "," not in alt_a and len(ref_a) == len(alt_a):
                kv = dict(x.split("=", 1) for x in cols[7].split(";") if "=" in x)
                if int(kv.get("SHIFT3", "0")) > 0:
                    vj_subst_shift3_pos += 1

    seen_keys = {}
    hg002_dup_ends = {}
    with gzip.open(off_gz, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            total_off_rows += 1
            cols = line.rstrip("\n").split("\t")
            chrom, pos_s, ref_a, alt_a, info_s, fmt_s = cols[0], cols[1], cols[3], cols[4], cols[7], cols[8]
            pos = int(pos_s)
            if pos < 1 or not ref_a:
                raise AssertionError(f"{name}: invalid POS/REF in {line}")
            if alt_a != "." and not valid_alt_re.match(alt_a):
                raise AssertionError(f"{name}: invalid/unanchored ALT allele {alt_a!r} in {line}")

            kv = dict(x.split("=", 1) for x in info_s.split(";") if "=" in x)
            if name == "GIAB_HG002_chr20" and alt_a == "<DUP>":
                hg002_dup_ends.setdefault((chrom, pos, ref_a), set()).add(int(kv["END"]))
            for single_alt in alt_a.split(","):
                if single_alt == "." or single_alt.startswith("<"):
                    continue
                if len(ref_a) != len(single_alt):
                    # Invariant 1b: parsimonious 5' anchoring (no redundant 5' prefix & no untrimmed 3' suffix when POS > 1)
                    if len(ref_a) > 1 and len(single_alt) > 1 and ref_a[0] == single_alt[0]:
                        raise AssertionError(f"{name}: redundant 5' prefix in {ref_a}->{single_alt}: {line}")
                    if pos > 1 and len(ref_a) > 1 and len(single_alt) > 1 and ref_a[-1] == single_alt[-1]:
                        raise AssertionError(f"{name}: untrimmed 3' suffix in {ref_a}->{single_alt}: {line}")
                    if int(kv.get("SHIFT3", "0")) > 0:
                        off_indel_shift3_pos += 1
                elif "," not in alt_a:
                    # Invariant 2: 0 SNV or equal-length MNV records carry SHIFT3 > 0 under --compat=off
                    if int(kv.get("SHIFT3", "0")) > 0:
                        off_subst_shift3_pos += 1
                        raise AssertionError(f"{name}: equal-length substitution leaked SHIFT3>0: {line}")

            # Accounting check across all sample columns
            fmt_keys = fmt_s.split(":")
            for smp in cols[9:]:
                vals = dict(zip(fmt_keys, smp.split(":")))
                dp = int(vals["DP"])
                vd = int(vals["VD"])
                ad = [int(x) for x in vals["AD"].split(",")]
                ald = [int(x) for x in vals["ALD"].split(",")]
                if vd > dp or dp < sum(ad) or sum(ald) != vd:
                    raise AssertionError(f"{name}: depth/ALD accounting violation in {line}")

            # Invariant 4: 0 duplicate (CHROM, POS, END, REF, ALT) shadow collisions in paired somatic --compat=off
            if is_paired:
                status = kv.get("STATUS", "")
                for single_alt in alt_a.split(","):
                    end_key = int(kv.get("END", pos)) if single_alt.startswith("<") else pos
                    key = (chrom, pos, end_key, ref_a, single_alt)
                    if key in seen_keys:
                        raise AssertionError(
                            f"{name}: duplicate somatic (CHROM, POS, END, REF, ALT) collision {key}: "
                            f"prev={seen_keys[key]} vs cur={status}"
                        )
                    seen_keys[key] = status
    if name == "GIAB_HG002_chr20":
        multi_end_dups = sum(1 for ends in hg002_dup_ends.values() if len(ends) > 1)
        if multi_end_dups < 20:
            raise AssertionError(
                f"Expected >= 20 distinct-END symbolic <DUP> loci preserved at shared POS on GIAB HG002 1.5 Mbp, got {multi_end_dups}"
            )
    if is_paired:
        somatic_cohorts_checked += 1

if vj_subst_shift3_pos <= 0:
    raise AssertionError("Expected --compat=vardictjava to exhibit > 0 SNVs/MNVs with spurious SHIFT3 > 0")
if off_subst_shift3_pos != 0:
    raise AssertionError(f"Expected 0 SNVs/MNVs with SHIFT3 > 0 under --compat=off, got {off_subst_shift3_pos}")
if off_indel_shift3_pos <= 0:
    raise AssertionError("Expected true indels under --compat=off to retain SHIFT3 > 0")

print(
    f"  PASS Gate 1 (bcftools norm --check-ref e & 5' left-alignment): {len(cohorts)} cohorts, "
    f"{total_off_rows} --compat=off records (0 ref mismatches, 0 unanchored indels)"
)
print(
    f"  PASS Gate 2 (Switch 38 indel-only SHIFT3): SNV/MNV SHIFT3>0 {vj_subst_shift3_pos} (--compat=vardictjava) "
    f"-> {off_subst_shift3_pos} (--compat=off); {off_indel_shift3_pos} indels retained SHIFT3>0"
)
print(
    f"  PASS Gate 4 (Switch 36 somatic split reconciliation): {somatic_cohorts_checked} paired somatic cohorts "
    f"(0 duplicate CHROM/POS/REF/ALT or StrongSomatic+StrongLOH shadow collisions)"
)
PY

# =============================================================================
# Gate 3: Paired-End Short-Insert `-u --compat=off` Bidirectional ALD Balance
#         (Switch 34: Clinical-Genomics/BALSAMIC#1396, VarDictJava#239)
# =============================================================================
echo "[3/5] Verifying -u --compat=off mate-overlap bidirectional ALD balance & DP/VD conservation..."
nice -n 19 python3 - "$VARDICTRUST_BIN" "$GRCH38_REF" "$GIAB_REF" "$G4_DIR" "$WAVE4_DIR" "$WORK_TMP" <<'PY'
import os
import subprocess
import sys

vdr, grch38_ref, giab_ref, g4_dir, w4_dir, tmp = sys.argv[1:7]
try:
    os.nice(19)
except OSError:
    pass

u_cases = [
    (
        "SEQC2_HCC1395_pe_u",
        grch38_ref,
        os.path.join(g4_dir, "seqc2", "HCC1395.bam"),
        ["-R", "chr17:74747665-74952140", "-u", "-f", "0.01", "-th", "4"],
    ),
    (
        "GIAB_HG002_pe_u",
        giab_ref,
        os.path.join(g4_dir, "hg002", "HG002.chr20.bam"),
        ["-R", "chr20:10000000-10080000", "-u", "-f", "0.01", "-th", "4"],
    ),
]

ilm_a = os.path.join(w4_dir, "bam", "SampleA_ILM1.bam")
ilm_bed = os.path.join(w4_dir, "beds", "ilm_twoprimer_grch38.amp8.bed")
if os.path.isfile(ilm_a) and os.path.isfile(ilm_bed):
    u_cases.append(
        (
            "Wave4_Illumina_PE_u",
            grch38_ref,
            ilm_a,
            ["-a", "10:0.95", "-z", "-u", "-f", "0.01", "-th", "4", ilm_bed],
        )
    )

checked_u_rows = 0
vj_one_strand_ge4 = 0
off_one_strand_ge4 = 0
rebalanced_loci = 0

for label, ref_fa, bam_p, args in u_cases:
    vj_vcf = subprocess.check_output([vdr, "-G", ref_fa, "-b", bam_p, "--compat=vardictjava", "-v", *args], text=True)
    off_vcf = subprocess.check_output([vdr, "-G", ref_fa, "-b", bam_p, "--compat=off", "-v", *args], text=True)

    vj_map = {}
    for line in vj_vcf.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        c = line.split("\t")
        vals = dict(zip(c[8].split(":"), c[9].split(":")))
        ald = [int(x) for x in vals["ALD"].split(",")]
        vd = int(vals["VD"])
        if vd >= 4 and (ald[0] == 0 or ald[1] == 0):
            vj_one_strand_ge4 += 1
        vj_map[(c[0], int(c[1]), c[3], c[4])] = tuple(ald)

    for line in off_vcf.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        checked_u_rows += 1
        c = line.split("\t")
        vals = dict(zip(c[8].split(":"), c[9].split(":")))
        dp = int(vals["DP"])
        vd = int(vals["VD"])
        ad = [int(x) for x in vals["AD"].split(",")]
        ald = [int(x) for x in vals["ALD"].split(",")]
        if vd > dp or dp < sum(ad) or (ald[0] + ald[1]) != vd:
            raise AssertionError(f"{label}: -u --compat=off accounting violation (DP={dp}, VD={vd}, AD={ad}, ALD={ald}): {line}")
        if vd >= 4 and (ald[0] == 0 or ald[1] == 0):
            off_one_strand_ge4 += 1
        key = (c[0], int(c[1]), c[3], c[4])
        if key in vj_map:
            vj_ald = vj_map[key]
            if vj_ald != tuple(ald) and ald[0] > 0 and ald[1] > 0:
                rebalanced_loci += 1

if off_one_strand_ge4 > vj_one_strand_ge4:
    raise AssertionError(
        f"Expected -u --compat=off to reduce or equal single-strand VD>=4 starvation: vj={vj_one_strand_ge4}, off={off_one_strand_ge4}"
    )
if rebalanced_loci <= 0:
    raise AssertionError("Expected >= 1 overlapping paired-end locus rebalanced to bidirectional ALD under -u --compat=off")

print(
    f"  PASS Gate 3 (Switch 34 -u mate-overlap strand balance): {checked_u_rows} rows audited "
    f"(0 VD>DP or ALD!=VD; {rebalanced_loci} loci rebalanced; single-strand VD>=4 {vj_one_strand_ge4} -> {off_one_strand_ge4})"
)
PY

# =============================================================================
# Gate 5: `rtg vcfeval` Ground-Truth Concordance Non-Inferiority / Improvement
#         across SEQC2 HCC1395/BL, GIAB HG008 v0.4 (if present), & GIAB HG002
# =============================================================================
echo "[5/5] Running rtg vcfeval truth concordance on SEQC2, GIAB HG008 v0.4, & GIAB HG002..."
nice -n 19 python3 - "$VARDICTRUST_BIN" "$GRCH38_REF" "$GIAB_REF" "$G4_DIR" "$WORK_TMP" <<'PY'
import os
import subprocess
import sys

vdr, grch38_ref, giab_ref, g4_dir, tmp = sys.argv[1:6]
try:
    os.nice(19)
except OSError:
    pass

def intersect_beds(bed_a, bed_b, out_bed, max_a_lines=None):
    b_map = {}
    with open(bed_b, "r", encoding="utf-8") as fh:
        for line in fh:
            if line.strip() and not line.startswith("#"):
                p = line.split("\t")
                b_map.setdefault(p[0], []).append((int(p[1]), int(p[2])))
    for c in b_map:
        b_map[c].sort()

    raw = {}
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

def build_scoped_sdf(ref_fa, contigs, prefix):
    scoped_fa = os.path.join(tmp, f"{prefix}_scoped.fa")
    with open(scoped_fa, "w", encoding="utf-8") as fh:
        subprocess.run(["samtools", "faidx", ref_fa, *contigs], stdout=fh, check=True)
    subprocess.run(["samtools", "faidx", scoped_fa], check=True)
    sdf_dir = os.path.join(tmp, f"{prefix}_scoped.sdf")
    subprocess.run(
        ["rtg", "RTG_MEM=4G", "RTG_JAVA_OPTS=-Dtalkback=false", "format", "-o", sdf_dir, scoped_fa],
        stdout=subprocess.DEVNULL,
        check=True,
    )
    return scoped_fa, sdf_dir

def parse_roc_last_row(roc_gz):
    import gzip as _gz
    if not os.path.isfile(roc_gz):
        return None
    last = None
    with _gz.open(roc_gz, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            last = line.strip().split()
    if last is None or len(last) < 8:
        return None
    return (int(float(last[1])), int(float(last[2])), int(float(last[4])), float(last[7]))

def parse_vcfeval_summary(summary_path):
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

def parse_vcfeval_dir(out_dir):
    all_row = parse_vcfeval_summary(os.path.join(out_dir, "summary.txt"))
    snp_row = parse_roc_last_row(os.path.join(out_dir, "snp_roc.tsv.gz"))
    indel_row = parse_roc_last_row(os.path.join(out_dir, "non_snp_roc.tsv.gz"))
    return all_row, snp_row, indel_row

# 5a) SEQC2 HCC1395|HCC1395BL somatic truth concordance (standard & paired-end -u)
seqc2_slices = os.path.join(g4_dir, "seqc2", "slices.bed")
seqc2_eval_src = os.path.join(g4_dir, "seqc2", "eval.bed")
seqc2_eval_bed = os.path.join(tmp, "seqc2_eval_scoped.bed")
seqc2_contigs = intersect_beds(seqc2_slices, seqc2_eval_src, seqc2_eval_bed)
seqc2_fa, seqc2_sdf = build_scoped_sdf(grch38_ref, seqc2_contigs, "seqc2")

seqc2_truth_raw = subprocess.check_output(
    ["bcftools", "view", "-R", seqc2_slices, os.path.join(g4_dir, "seqc2", "truth", "truth.vcf.gz")],
    text=True,
)
with open(seqc2_fa + ".fai", encoding="utf-8") as fai_fh:
    seqc2_lens = {l.split("\t")[0]: int(l.split("\t")[1]) for l in fai_fh}
seqc2_truth_vcf = os.path.join(tmp, "seqc2_truth.vcf")
with open(seqc2_truth_vcf, "w", encoding="utf-8") as out:
    out.write("##fileformat=VCFv4.2\n")
    for c in seqc2_contigs:
        out.write(f"##contig=<ID={c},length={seqc2_lens[c]}>\n")
    out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\n")
    for line in seqc2_truth_raw.splitlines():
        if not line.startswith("#") and line.strip():
            p = line.split("\t")
            if p[0] in seqc2_lens:
                out.write(f"{p[0]}\t{p[1]}\t{p[2]}\t{p[3]}\t{p[4]}\t.\tPASS\t.\tGT\t0/1\n")
subprocess.run(["bgzip", "-f", seqc2_truth_vcf], check=True)
subprocess.run(["tabix", "-f", "-p", "vcf", seqc2_truth_vcf + ".gz"], check=True)

def eval_somatic_rtg(tag, ref_fa, sdf_dir, bam_pair, bed_path, eval_bed, truth_gz, compat_flag, extra_flags=None):
    extra_flags = extra_flags or []
    raw_vcf = os.path.join(tmp, f"{tag}_{compat_flag.split('=')[1]}_raw.vcf")
    som_vcf = os.path.join(tmp, f"{tag}_{compat_flag.split('=')[1]}_som.vcf")
    subprocess.run(
        [vdr, "-G", ref_fa, "-b", bam_pair, "-N", "TUMOR|NORMAL",
         "-c", "1", "-S", "2", "-E", "3", "-z", "-f", "0.01", "--fisher",
         compat_flag, "--vcf-all-alleles", "--output", raw_vcf, "-th", "4", *extra_flags, bed_path],
        check=True,
    )
    with open(raw_vcf, "r", encoding="utf-8") as inp, open(som_vcf, "w", encoding="utf-8") as out:
        for line in inp:
            if line.startswith("#"):
                out.write(line)
            else:
                p = line.split("\t")
                if p[6] == "PASS" and ("STATUS=StrongSomatic" in p[7] or "STATUS=LikelySomatic" in p[7]):
                    out.write(line)
    subprocess.run(["bgzip", "-f", som_vcf], check=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", som_vcf + ".gz"], check=True)
    out_dir = os.path.join(tmp, f"rtg_{tag}_{compat_flag.split('=')[1]}")
    subprocess.run(
        ["rtg", "RTG_MEM=4G", "RTG_JAVA_OPTS=-Dtalkback=false", "vcfeval",
         "-T", "4", "-b", truth_gz, "-c", som_vcf + ".gz", "-e", eval_bed, "-t", sdf_dir,
         "--squash-ploidy", "--vcf-score-field=QUAL", "--sample=ALT,TUMOR", "-o", out_dir],
        stdout=subprocess.DEVNULL,
        check=True,
    )
    return parse_vcfeval_summary(os.path.join(out_dir, "summary.txt"))

seqc2_pair = f"{g4_dir}/seqc2/HCC1395.bam|{g4_dir}/seqc2/HCC1395BL.bam"
s_vj = eval_somatic_rtg("seqc2", seqc2_fa, seqc2_sdf, seqc2_pair, seqc2_slices, seqc2_eval_bed, seqc2_truth_vcf + ".gz", "--compat=vardictjava")
s_off = eval_somatic_rtg("seqc2", seqc2_fa, seqc2_sdf, seqc2_pair, seqc2_slices, seqc2_eval_bed, seqc2_truth_vcf + ".gz", "--compat=off")
if s_off[3] < s_vj[3] or s_off[0] < s_vj[0] or s_off[1] > s_vj[1]:
    raise AssertionError(f"SEQC2 rtg vcfeval regression: vj={s_vj} vs off={s_off}")

su_vj = eval_somatic_rtg("seqc2_u", seqc2_fa, seqc2_sdf, seqc2_pair, seqc2_slices, seqc2_eval_bed, seqc2_truth_vcf + ".gz", "--compat=vardictjava", ["-u"])
su_off = eval_somatic_rtg("seqc2_u", seqc2_fa, seqc2_sdf, seqc2_pair, seqc2_slices, seqc2_eval_bed, seqc2_truth_vcf + ".gz", "--compat=off", ["-u"])
# Note: under `-u`, v1.1.0 Switch 18 (`somatic_continue_after_bad_var`) evaluates secondary alleles at multi-allelic
# loci (+1 candidate vs v1.8.3 early-return); require zero TP loss and non-inferior best-threshold F1!
if su_off[0] < su_vj[0] or su_off[4] + 1e-6 < su_vj[4]:
    raise AssertionError(f"SEQC2 -u rtg vcfeval TP/ROC regression: vj={su_vj} vs off={su_off}")

# 5b) GIAB HG002 chr20:10000000-11500000 (1.5 Mbp) germline truth concordance (SNP + INDEL stratified, diploid GT + squash-ploidy, standard + -u)
hg2_call_bed = os.path.join(tmp, "hg002_1_5mb.bed")
with open(hg2_call_bed, "w", encoding="utf-8") as fh:
    fh.write("chr20\t10000000\t11500000\n")
hg2_eval_bed = os.path.join(tmp, "hg002_1_5mb_eval.bed")
intersect_beds(hg2_call_bed, os.path.join(g4_dir, "hg002", "eval.bed"), hg2_eval_bed)
hg2_fa, hg2_sdf = build_scoped_sdf(giab_ref, ["chr20"], "hg002_c20")
hg2_truth_src = os.path.join(g4_dir, "hg002", "truth", "truth.vcf.gz")
hg2_truth_scoped = os.path.join(tmp, "hg002_c20_truth.vcf.gz")
subprocess.run(
    ["bcftools", "view", "-R", hg2_call_bed, "-Oz", "-o", hg2_truth_scoped, hg2_truth_src],
    check=True,
)
subprocess.run(["tabix", "-f", "-p", "vcf", hg2_truth_scoped], check=True)

def eval_hg002_rtg(tag, compat_flag, extra_flags=None, check_diploid=False):
    extra_flags = extra_flags or []
    out_gz = os.path.join(tmp, f"{tag}_{compat_flag.split('=')[1]}.vcf.gz")
    subprocess.run(
        [vdr, "-G", hg2_fa, "-b", os.path.join(g4_dir, "hg002", "HG002.chr20.bam"),
         "-N", "HG002", "-c", "1", "-S", "2", "-E", "3", "-z", "-f", "0.05",
         "--pass-only", compat_flag, "--vcf-all-alleles", "--output", out_gz, "--tbi", "-th", "4", *extra_flags, hg2_call_bed],
        check=True,
    )
    sq_dir = os.path.join(tmp, f"rtg_{tag}_{compat_flag.split('=')[1]}_sq")
    subprocess.run(
        ["rtg", "RTG_MEM=4G", "RTG_JAVA_OPTS=-Dtalkback=false", "vcfeval",
         "-T", "4", "-b", hg2_truth_scoped, "-c", out_gz, "-e", hg2_eval_bed, "-t", hg2_sdf,
         "--squash-ploidy", "--vcf-score-field=QUAL", "-o", sq_dir],
        stdout=subprocess.DEVNULL,
        check=True,
    )
    sq_res = parse_vcfeval_dir(sq_dir)
    gt_res = None
    if check_diploid:
        gt_dir = os.path.join(tmp, f"rtg_{tag}_{compat_flag.split('=')[1]}_gt")
        subprocess.run(
            ["rtg", "RTG_MEM=4G", "RTG_JAVA_OPTS=-Dtalkback=false", "vcfeval",
             "-T", "4", "-b", hg2_truth_scoped, "-c", out_gz, "-e", hg2_eval_bed, "-t", hg2_sdf,
             "--vcf-score-field=QUAL", "-o", gt_dir],
            stdout=subprocess.DEVNULL,
            check=True,
        )
        gt_res = parse_vcfeval_dir(gt_dir)
    return sq_res, gt_res

(g2_vj, _), (g2_off, g2_off_gt) = (
    eval_hg002_rtg("hg002", "--compat=vardictjava"),
    eval_hg002_rtg("hg002", "--compat=off", check_diploid=True),
)
for stratum_name, vj_m, off_m in [("ALL", g2_vj[0], g2_off[0]), ("SNP", g2_vj[1], g2_off[1]), ("INDEL", g2_vj[2], g2_off[2])]:
    if vj_m is None or off_m is None:
        raise AssertionError(f"Missing GIAB HG002 {stratum_name} ROC stratum")
    if off_m[3] + 1e-6 < vj_m[3] or off_m[0] < vj_m[0] or off_m[1] > vj_m[1]:
        raise AssertionError(f"GIAB HG002 1.5M {stratum_name} rtg vcfeval regression: vj={vj_m} vs off={off_m}")
if (
    g2_off_gt is None
    or g2_off_gt[0] is None
    or g2_off_gt[1] is None
    or g2_off_gt[2] is None
    or g2_off_gt[0][3] < 0.9590
    or g2_off_gt[1][3] < 0.9800
    or g2_off_gt[2][3] < 0.8250
):
    raise AssertionError(f"GIAB HG002 1.5M diploid GT rtg vcfeval below threshold: {g2_off_gt}")

(g2u_vj, _), (g2u_off, _) = (
    eval_hg002_rtg("hg002_u", "--compat=vardictjava", ["-u"]),
    eval_hg002_rtg("hg002_u", "--compat=off", ["-u"]),
)
for stratum_name, vj_m, off_m in [("ALL", g2u_vj[0], g2u_off[0]), ("SNP", g2u_vj[1], g2u_off[1]), ("INDEL", g2u_vj[2], g2u_off[2])]:
    if vj_m is None or off_m is None:
        raise AssertionError(f"Missing GIAB HG002 -u {stratum_name} ROC stratum")
    if off_m[3] + 1e-6 < vj_m[3] or off_m[0] < vj_m[0] or off_m[1] > vj_m[1]:
        raise AssertionError(f"GIAB HG002 1.5M -u {stratum_name} rtg vcfeval regression: vj={vj_m} vs off={off_m}")

# 5c) GIAB HG008 somatic v0.4 benchmark (standard + paired-end -u, if present locally)
hg8_truth_v04 = os.path.join(g4_dir, "hg008", "truth", "HG008-T_somatic_smvar_benchmark_v0.4_GRCh38_tumorvariants.vcf.gz")
hg8_bed_v04 = os.path.join(g4_dir, "hg008", "truth", "HG008-T_somatic_smvar_benchmark_v0.4_GRCh38_all.bed")
hg8_slices = os.path.join(g4_dir, "hg008", "slices.bed")
hg8_msg = "skipped (not cached)"
if os.path.isfile(hg8_truth_v04) and os.path.isfile(hg8_bed_v04) and os.path.isfile(hg8_slices):
    hg8_eval_bed = os.path.join(tmp, "hg008_v04_eval.bed")
    hg8_contigs = intersect_beds(hg8_slices, hg8_bed_v04, hg8_eval_bed, max_a_lines=40)
    hg8_fa, hg8_sdf = build_scoped_sdf(giab_ref, hg8_contigs, "hg008_v04")
    with open(hg8_fa + ".fai", encoding="utf-8") as fai_fh:
        hg8_lens = {l.split("\t")[0]: int(l.split("\t")[1]) for l in fai_fh}
    hg8_raw = subprocess.check_output(["bcftools", "view", "-f", ".,PASS", "-R", hg8_eval_bed, hg8_truth_v04], text=True)
    hg8_truth_vcf = os.path.join(tmp, "hg008_v04_scoped.vcf")
    with open(hg8_truth_vcf, "w", encoding="utf-8") as out:
        out.write("##fileformat=VCFv4.2\n")
        for c in hg8_contigs:
            out.write(f"##contig=<ID={c},length={hg8_lens[c]}>\n")
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\n")
        for line in hg8_raw.splitlines():
            if not line.startswith("#") and line.strip():
                p = line.split("\t")
                if p[0] in hg8_lens:
                    out.write(f"{p[0]}\t{p[1]}\t{p[2]}\t{p[3]}\t{p[4]}\t.\tPASS\t.\tGT\t0/1\n")
    subprocess.run(["bgzip", "-f", hg8_truth_vcf], check=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", hg8_truth_vcf + ".gz"], check=True)

    hg8_pair = f"{g4_dir}/hg008/HG008-T.slice.bam|{g4_dir}/hg008/HG008-N.slice.bam"
    h8_vj = eval_somatic_rtg("hg008v04", hg8_fa, hg8_sdf, hg8_pair, hg8_eval_bed, hg8_eval_bed, hg8_truth_vcf + ".gz", "--compat=vardictjava")
    h8_off = eval_somatic_rtg("hg008v04", hg8_fa, hg8_sdf, hg8_pair, hg8_eval_bed, hg8_eval_bed, hg8_truth_vcf + ".gz", "--compat=off")
    if h8_off[3] < h8_vj[3] or h8_off[0] < h8_vj[0] or h8_off[1] > h8_vj[1]:
        raise AssertionError(f"GIAB HG008 v0.4 rtg vcfeval regression: vj={h8_vj} vs off={h8_off}")
    h8u_vj = eval_somatic_rtg("hg008v04_u", hg8_fa, hg8_sdf, hg8_pair, hg8_eval_bed, hg8_eval_bed, hg8_truth_vcf + ".gz", "--compat=vardictjava", ["-u"])
    h8u_off = eval_somatic_rtg("hg008v04_u", hg8_fa, hg8_sdf, hg8_pair, hg8_eval_bed, hg8_eval_bed, hg8_truth_vcf + ".gz", "--compat=off", ["-u"])
    if h8u_off[3] < h8u_vj[3] or h8u_off[0] < h8u_vj[0] or h8u_off[1] > h8u_vj[1]:
        raise AssertionError(f"GIAB HG008 v0.4 -u rtg vcfeval regression: vj={h8u_vj} vs off={h8u_off}")
    hg8_msg = (
        f"HG008 v0.4 F1 {h8_vj[3]:.4f}->{h8_off[3]:.4f} (FP {h8_vj[1]}->{h8_off[1]}), "
        f"-u F1 {h8u_vj[3]:.4f}->{h8u_off[3]:.4f} (FP {h8u_vj[1]}->{h8u_off[1]})"
    )

print(
    f"  PASS Gate 5 (rtg vcfeval standard & -u, SNP/INDEL & diploid GT): "
    f"SEQC2 F1 {s_vj[3]:.4f}->{s_off[3]:.4f} (-u TP={su_off[0]}), "
    f"HG002 1.5M F1 {g2_vj[0][3]:.4f}->{g2_off[0][3]:.4f} (diploid GT F1={g2_off_gt[0][3]:.4f}, Indel F1={g2_off[2][3]:.4f}, -u F1 {g2u_vj[0][3]:.4f}->{g2u_off[0][3]:.4f}), "
    f"{hg8_msg}"
)
PY

ELAPSED_S=$(( ( $(now_ms) - T0_MS ) / 1000 ))
echo "================================================================================"
echo " ALL 5 WAVE-6 PUBLIC TRUTH-SET & INVARIANT GATES PASSED (${ELAPSED_S}s)"
echo "================================================================================"
