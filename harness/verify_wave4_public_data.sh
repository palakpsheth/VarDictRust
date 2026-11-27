#!/usr/bin/env bash
# =============================================================================
# harness/verify_wave4_public_data.sh
# Comprehensive 6-cohort real-world clinical validation & stress-test suite for
# Wave-4 features, CompatProfile switches 18..21, single-primer amplicon (-sp),
# raw -f <TEXT> VCF FILTER parity, and vardictcpp v2..v4 CLI aliases.
#
# Cross-Ecosystem Provenance & Citation Registry (owner/repo#N + DOIs):
# - Item 1 (Raw -f <TEXT> VCF FILTER ID & Header Parity):
#   AstraZeneca-NGS/VarDict@009e017 (var2vcf_valid.pl, var2vcf_paired.pl);
#   MHH-Bioinformatics-Hematology/vardictcpp#4 (2f6d811c23);
#   Lai et al., Nucleic Acids Res. 2016 (DOI: 10.1093/nar/gkw227).
# - Item 2 (Drop-In vardictcpp v2..v4 CLI Aliases --vcfall/--noend/--vcfpass/--vcf):
#   MHH-Bioinformatics-Hematology/vardictcpp#4 (a689e43877),
#   MHH-Bioinformatics-Hematology/vardictcpp#1 (4f7e242e0a),
#   bioconda/bioconda-recipes#68008, bioconda/bioconda-recipes#71981,
#   bioconda/bioconda-recipes#71875, nf-core/modules#13830,
#   AstraZeneca-NGS/VarDictJava#65, AstraZeneca-NGS/VarDict#155,
#   palakpsheth/VarDictRust#94, palakpsheth/VarDictRust#92,
#   palakpsheth/VarDictRust#156.
# - Item 4 (CompatProfile switch 18: somatic_continue_after_bad_var):
#   joachimwolff/VarDictJava#4, clintval/vartovcf#139, clintval/vartovcf#149,
#   AstraZeneca-NGS/VarDictJava#376; Fang et al., Nat. Biotechnol. 2021
#   (DOI: 10.1038/s41587-021-00993-6, SEQC2 HCC1395); NIST GIAB HG008 v0.4/v0.6.
# - Item 5 (CompatProfile switch 19: adaptor_all_consensus_dirs):
#   AstraZeneca-NGS/VarDictJava#240, AstraZeneca-NGS/VarDictJava#241,
#   AstraZeneca-NGS/VarDict#119,
#   PolinaBevad/VarDictJava@92c62cf6c3520b877c394a4c9fd5fe0d371cd972.
# - Item 6 (CompatProfile switch 20: canonical_sv_vcf_anchoring):
#   AstraZeneca-NGS/VarDictJava#398, AstraZeneca-NGS/VarDictJava#225,
#   AstraZeneca-NGS/VarDictJava#229, AstraZeneca-NGS/VarDict#111,
#   AstraZeneca-NGS/VarDictJava#382, clintval/vartovcf#192 (baaaddb415),
#   clintval/vartovcf#107, clintval/vartovcf#108,
#   clintval/vartovcf@3b7bd742926c62988ec2fd8420b5f237f43db6ff,
#   samtools/hts-specs VCFv4.2 S1.4.5 & S5.4.
# - Item 7 (--amplicon-single-primer / -sp):
#   AstraZeneca-NGS/VarDictJava#301, AstraZeneca-NGS/VarDictJava#303,
#   AstraZeneca-NGS/VarDictJava#326, AstraZeneca-NGS/VarDict#191;
#   Gong et al., Sci. Data 2022 (DOI: 10.1038/s41597-022-01359-6, ENA SRP295113,
#   Figshare DOI: 10.6084/m9.figshare.19128005.v2 &
#   DOI: 10.6084/m9.figshare.13511829.v4).
# - Item 8 (CompatProfile switch 21: paired_somatic_amplicon):
#   AstraZeneca-NGS/VarDictJava#378, AstraZeneca-NGS/VarDictJava#155,
#   AstraZeneca-NGS/VarDict#109, AstraZeneca-NGS/VarDictJava#402,
#   MHH-Bioinformatics-Hematology/vardictcpp#1,
#   MHH-Bioinformatics-Hematology/vardictcpp#3,
#   palakpsheth/VarDictRust#95, palakpsheth/VarDictRust#94.
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

for _req_tool in samtools bcftools tabix perl python3 gzip; do
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

# Ensure all 6 clinical cohorts are staged and SHA-256 verified
DATA_DIR="$G4_DIR" WAVE4_DIR="$WAVE4_DIR" bash "$SCRIPT_DIR/get_wave4_public_slices.sh"

FA_HG19="$DATA_ROOT/VarDictJava/testdata/fastas/hg19.fa"
FA_B37="$DATA_ROOT/VarDictJava/testdata/fastas/GRCh37.fa"
FA_HG1K="$DATA_ROOT/VarDictJava/testdata/fastas/hg1k.fa"
if [ ! -f "$FA_HG19" ] || [ ! -f "$FA_B37" ] || [ ! -f "$FA_HG1K" ]; then
    python3 "$ROOT_DIR/harness/convert_fastas.py"
fi

GRCH38_REF="$G4_DIR/ref/d1vd1/GRCh38.d1.vd1.fa"
GIAB_REF="$G4_DIR/ref/giabv3/GRCh38_GIABv3.fa"
if [ ! -f "$GIAB_REF" ] && [ -f "$G4_DIR/ref/giabv3/GRCh38_GIABv3.fasta" ]; then
    GIAB_REF="$G4_DIR/ref/giabv3/GRCh38_GIABv3.fasta"
fi

PERL_VALID="$DATA_ROOT/VarDictJava/VarDict/var2vcf_valid.pl"
PERL_PAIRED="$DATA_ROOT/VarDictJava/VarDict/var2vcf_paired.pl"
VJ_BAMS="$DATA_ROOT/VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests"

if ! { [ -s "$GIAB_REF" ] && [ -s "${GIAB_REF}.fai" ] && [ -f "$PERL_VALID" ] && [ -f "$PERL_PAIRED" ]; }; then
    echo "ERROR: missing Wave-4 verification prerequisite(s): require non-empty GIAB_REF ('$GIAB_REF' + .fai) and Perl scripts ('$PERL_VALID', '$PERL_PAIRED')" >&2
    exit 1
fi

WORK_TMP="$(mktemp -d "${TMPDIR:-/tmp}/vdr_wave4_verify.XXXXXX")"
trap 'rm -rf "$WORK_TMP"' EXIT

T0_MS="$(now_ms)"
echo "================================================================================"
echo " Wave-4 6-Cohort Clinical Stress-Test & Verification Suite"
echo " Binary: $VARDICTRUST_BIN"
echo " Data:   $WAVE4_DIR"
echo "================================================================================"

# =============================================================================
# Suite 1: Raw -f <TEXT> VCF FILTER ID & Header Parity vs Live Perl var2vcf_*.pl
#          (AstraZeneca-NGS/VarDict@009e017, MHH-Bioinformatics-Hematology/vardictcpp#4)
# =============================================================================
echo "[1/7] Verifying raw -f <TEXT> VCF FILTER parity vs live Perl var2vcf_valid.pl & var2vcf_paired.pl..."
python3 - "$VARDICTRUST_BIN" "$PERL_VALID" "$PERL_PAIRED" "$GIAB_REF" "$FA_B37" "$FA_HG19" "$VJ_BAMS" "$WAVE4_DIR" "$WORK_TMP" <<'PY'
import gzip, os, subprocess, sys

vdr, perl_valid, perl_paired, giab_ref, fa_b37, fa_hg19, vj_bams, w4, tmp = sys.argv[1:10]
spellings = ["0", "0.0", "0.050", "0.00001", "1e-5", ".01"]

def parse_vcf_filters(vcf_text):
    f_headers = [
        l.strip() for l in vcf_text.splitlines()
        if l.startswith("##FILTER=<ID=f") and 'Description="Allele frequency < ' in l
    ]
    rows = {}
    for l in vcf_text.splitlines():
        if l.startswith("#") or not l.strip():
            continue
        p = l.split("\t")
        rows[(p[0], int(p[1]), p[3], p[4])] = p[6]
    return f_headers, rows

# 1a) Single-sample HG002 CMRG + MosaicSNVv1.1 & NA24385
hg002_bam = os.path.join(w4, "hg002_cmrg_mosaic", "HG002.cmrg_mosaic.slice.bam")
hg002_bed = os.path.join(w4, "hg002_cmrg_mosaic", "hg002_cmrg_mosaic.bed")
na_bam = os.path.join(vj_bams, "NA24385_chr_1_246981051_246981611.bam")

hg002_pass_sets = {}
for raw_f in spellings:
    # HG002 single-sample vs var2vcf_valid.pl (--emit-intermediate --fisher produces the 38-col table var2vcf_valid.pl reads)
    var_out = subprocess.check_output(
        [vdr, "-G", giab_ref, "-b", hg002_bam, "-c", "1", "-S", "2", "-E", "3", "-g", "4", "-N", "HG002", "--emit-intermediate", "--fisher", "-f", raw_f, hg002_bed],
        text=True,
    )
    perl_vcf = subprocess.check_output(
        ["perl", perl_valid, "-N", "HG002", "-E", "-f", raw_f],
        input=var_out,
        text=True,
    )
    rust_vcf = subprocess.check_output(
        [vdr, "-G", giab_ref, "-b", hg002_bam, "-c", "1", "-S", "2", "-E", "3", "-g", "4", "-N", "HG002", "-v", "--fisher", "--no-end", "-f", raw_f, hg002_bed],
        text=True,
    )
    ph, pr = parse_vcf_filters(perl_vcf)
    rh, rr = parse_vcf_filters(rust_vcf)
    if ph != rh:
        raise AssertionError(f"HG002 -f {raw_f} header mismatch: perl={ph} vs rust={rh}")
    if pr != rr:
        raise AssertionError(f"HG002 -f {raw_f} FILTER row mismatch: perl_n={len(pr)} rust_n={len(rr)}")
    hg002_pass_sets[raw_f] = {k for k, flt in rr.items() if flt == "PASS"}

    # NA24385 single-sample vs var2vcf_valid.pl
    na_var = subprocess.check_output(
        [vdr, "-G", fa_b37, "-b", na_bam, "-R", "1:246981051-246981611", "-N", "NA24385", "--emit-intermediate", "--fisher", "-f", raw_f],
        text=True,
    )
    na_perl = subprocess.check_output(["perl", perl_valid, "-N", "NA24385", "-E", "-f", raw_f], input=na_var, text=True)
    na_rust = subprocess.check_output(
        [vdr, "-G", fa_b37, "-b", na_bam, "-R", "1:246981051-246981611", "-N", "NA24385", "-v", "--fisher", "--no-end", "-f", raw_f],
        text=True,
    )
    if parse_vcf_filters(na_perl) != parse_vcf_filters(na_rust):
        raise AssertionError(f"NA24385 -f {raw_f} Perl vs Rust VCF FILTER mismatch")

    # Colo829 paired somatic vs var2vcf_paired.pl (--emit-intermediate --fisher produces the 61-col table var2vcf_paired.pl reads)
    colo_pair = f"{os.path.join(vj_bams, 'Colo829-18_S3-sort.bam')}|{os.path.join(vj_bams, 'Colo829-19_S4-sort.bam')}"
    c_var = subprocess.check_output(
        [vdr, "-G", fa_hg19, "-b", colo_pair, "-R", "chr1:755917-756517", "-N", "TUMOR|NORMAL", "--emit-intermediate", "--fisher", "-f", raw_f],
        text=True,
    )
    c_perl = subprocess.check_output(["perl", perl_paired, "-N", "TUMOR|NORMAL", "-f", raw_f], input=c_var, text=True)
    c_rust = subprocess.check_output(
        [vdr, "-G", fa_hg19, "-b", colo_pair, "-R", "chr1:755917-756517", "-N", "TUMOR|NORMAL", "-v", "--fisher", "-f", raw_f],
        text=True,
    )
    if parse_vcf_filters(c_perl) != parse_vcf_filters(c_rust):
        raise AssertionError(f"Colo829 paired -f {raw_f} Perl vs Rust VCF FILTER mismatch")

# Verify -f 0.050 followed by -p (pileup) preserves pre-pileup f0.050 label
pile_rust = subprocess.check_output(
    [vdr, "-G", fa_b37, "-b", na_bam, "-R", "1:246981290-246981320", "-N", "NA24385", "-v", "--no-end", "-f", "0.050", "-p"],
    text=True,
)
ph_pile, _ = parse_vcf_filters(pile_rust)
if ph_pile != ['##FILTER=<ID=f0.050,Description="Allele frequency < 0.050">']:
    raise AssertionError(f"Expected f0.050 header under -f 0.050 -p, got {ph_pile}")

# Verify 0 precision/recall shift across equivalent numeric spellings on HG002 CMRG + MosaicSNVv1.1
if hg002_pass_sets["0.00001"] != hg002_pass_sets["1e-5"]:
    raise AssertionError("HG002 PASS set diverged between -f 0.00001 and -f 1e-5")
rust_01 = subprocess.check_output(
    [vdr, "-G", giab_ref, "-b", hg002_bam, "-c", "1", "-S", "2", "-E", "3", "-g", "4", "-N", "HG002", "-v", "--fisher", "--no-end", "-f", "0.01", hg002_bed],
    text=True,
)
_, rr_01 = parse_vcf_filters(rust_01)
if {k for k, flt in rr_01.items() if flt == "PASS"} != hg002_pass_sets[".01"]:
    raise AssertionError("HG002 PASS set diverged between -f 0.01 and -f .01")

# Verify HG002 CMRG v1.00 + MosaicSNVv1.1 truth overlap is non-empty
truth_loci = set()
for fn in ("HG002_GRCh38_CMRG_smallvar_v1.00.vcf.gz", "HG002_GRCh38_MosaicSNVv1.1.vcf.gz"):
    with gzip.open(os.path.join(w4, "hg002_cmrg_mosaic", fn), "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.split("\t")
            truth_loci.add((p[0], int(p[1])))
matched_truth = sum(1 for (c, pos, _r, _a) in hg002_pass_sets[".01"] if (c, pos) in truth_loci)
if matched_truth < 10:
    raise AssertionError(f"Expected >= 10 HG002 CMRG/MosaicSNVv1.1 truth calls in slice, got {matched_truth}")
print(f"  PASS Suite 1: 6 -f spellings x 3 cohorts matched live Perl; HG002 CMRG/Mosaic truth hits={matched_truth}")
PY

# =============================================================================
# Suite 2: Drop-In vardictcpp v2..v4 CLI Aliases (--vcfall/--noend/--vcfpass/--vcf)
#          (MHH-Bioinformatics-Hematology/vardictcpp#4, palakpsheth/VarDictRust#94)
# =============================================================================
echo "[2/7] Verifying vardictcpp v2..v4 CLI alias byte-equivalence (--vcfall, --noend, --vcfpass, --vcf)..."
python3 - "$VARDICTRUST_BIN" "$GIAB_REF" "$FA_HG19" "$VJ_BAMS" "$WAVE4_DIR" <<'PY'
import os, subprocess, sys

vdr, giab_ref, fa_hg19, vj_bams, w4 = sys.argv[1:6]
hg002_bam = os.path.join(w4, "hg002_cmrg_mosaic", "HG002.cmrg_mosaic.slice.bam")
hg002_bed = os.path.join(w4, "hg002_cmrg_mosaic", "hg002_cmrg_mosaic.bed")
colo_pair = f"{os.path.join(vj_bams, 'Colo829-18_S3-sort.bam')}|{os.path.join(vj_bams, 'Colo829-19_S4-sort.bam')}"

# Note: NEVER pass bare -A/-E/-S as boolean VCF flags (-E is BED end column, -S is BED start column!)
alias_pairs = [
    (["--vcf"], ["-v"]),
    (["-v", "--vcfall"], ["-v", "--vcf-all-alleles"]),
    (["-v", "--noend"], ["-v", "--no-end"]),
    (["-v", "--vcfpass"], ["-v", "--pass-only"]),
    (["--vcf", "--vcfall", "--noend", "--vcfpass"], ["-v", "--vcf-all-alleles", "--no-end", "--pass-only"]),
]

for lhs, rhs in alias_pairs:
    out_l = subprocess.check_output(
        [vdr, "-G", giab_ref, "-b", hg002_bam, "-c", "1", "-S", "2", "-E", "3", "-g", "4", "-f", "0.01", *lhs, hg002_bed]
    )
    out_r = subprocess.check_output(
        [vdr, "-G", giab_ref, "-b", hg002_bam, "-c", "1", "-S", "2", "-E", "3", "-g", "4", "-f", "0.01", *rhs, hg002_bed]
    )
    if out_l != out_r:
        raise AssertionError(f"Single-sample CLI alias byte mismatch: {lhs} vs {rhs}")

    p_l = subprocess.check_output([vdr, "-G", fa_hg19, "-b", colo_pair, "-R", "chr1:755917-756517", "-f", "0.01", *lhs])
    p_r = subprocess.check_output([vdr, "-G", fa_hg19, "-b", colo_pair, "-R", "chr1:755917-756517", "-f", "0.01", *rhs])
    if p_l != p_r:
        raise AssertionError(f"Paired CLI alias byte mismatch: {lhs} vs {rhs}")

# Verify --vcf + --emit-intermediate conflict errors cleanly
bad_rc = subprocess.call(
    [vdr, "-G", fa_hg19, "-b", colo_pair, "-R", "chr1:755917-756517", "--vcf", "--emit-intermediate"],
    stdout=subprocess.DEVNULL,
    stderr=subprocess.DEVNULL,
)
if bad_rc == 0:
    raise AssertionError("Expected non-zero exit when combining --vcf with --emit-intermediate")

help_txt = subprocess.check_output([vdr, "--help"], text=True)
for hidden_flag in ("--vcfall", "--noend", "--vcfpass"):
    if hidden_flag in help_txt:
        raise AssertionError(f"Hidden alias {hidden_flag} leaked into --help output")
if "--amplicon-single-primer" not in help_txt:
    raise AssertionError("--amplicon-single-primer (-sp) missing from --help output")
print("  PASS Suite 2: all vardictcpp v2..v4 CLI aliases byte-identical across single & paired cohorts")
PY

# =============================================================================
# Suite 3: CompatProfile Switch 18 (somatic_continue_after_bad_var)
#          (joachimwolff/VarDictJava#4, clintval/vartovcf#139 & clintval/vartovcf#149,
#           AstraZeneca-NGS/VarDictJava#376, SEQC2 HCC1395, GIAB HG008 v0.4)
# =============================================================================
echo "[3/7] Verifying CompatProfile switch 18 (somatic_continue_after_bad_var) on 6 BAM pairs + MUC19 + SEQC2/HG008..."
python3 - "$VARDICTRUST_BIN" "$FA_HG19" "$FA_HG1K" "$GRCH38_REF" "$GIAB_REF" "$VJ_BAMS" "$G4_DIR" "$WAVE4_DIR" "$WORK_TMP" <<'PY'
import gzip, os, subprocess, sys

vdr, fa_hg19, fa_hg1k, grch38_ref, giab_ref, vj_bams, g4_dir, w4, tmp = sys.argv[1:10]

# 3a) All 6 joachimwolff/VarDictJava#4 tumor/normal multi-allelic pairs
# Note: in ArgumentNA_Tumor|Normal, 140453065 C>A (Q=7.3 < -q 25) and 140453033 G>A (SN=1.0 < 1.5)
# are the failing primary candidates that prematurely aborted VarDictJava v1.8.3's somatic loop;
# somatic_continue_after_bad_var recovers 16 passing secondary alleles across the slice.
arg_secondary = {
    (140453033, "T"),
    (140453035, "A"),
    (140453036, "CA"),
    (140453065, "T"),
    (140453081, "T"),
    (140453095, "G"),
    (140453108, "A"),
    (140453116, "G"),
}
cases = [
    (
        "ArgumentNA_fisher",
        fa_hg19,
        f"{vj_bams}/ArgumentNA_Tumor.bam|{vj_bams}/ArgumentNA_Normal.bam",
        ["-R", "chr7:140453010-140453158", "-q", "25", "-f", "0.001", "-F", "0x100", "-k", "1", "-M", "20", "-p", "--fisher"],
        lambda pos, ref, alt: (pos, alt) in arg_secondary,
        8,
    ),
    (
        "ArgumentNA_nofisher",
        fa_hg19,
        f"{vj_bams}/ArgumentNA_Tumor.bam|{vj_bams}/ArgumentNA_Normal.bam",
        ["-R", "chr7:140453010-140453158", "-q", "25", "-f", "0.001", "-F", "0x100", "-k", "1", "-M", "20", "-p"],
        lambda pos, ref, alt: (pos, alt) in arg_secondary,
        8,
    ),
    (
        "Colo829_chr1_756309",
        fa_hg19,
        f"{vj_bams}/Colo829-18_S3-sort.bam|{vj_bams}/Colo829-19_S4-sort.bam",
        ["-R", "chr1:755917-756517", "-f", "0.01"],
        lambda pos, ref, alt: pos == 756309 and ref == "TTGT" and alt == "CTGG",
        1,
    ),
    (
        "HCC2218_chrX_7811864",
        fa_hg19,
        f"{vj_bams}/HCC2218_BL_realigins.bam|{vj_bams}/HCC2218_C_realigins.bam",
        ["-R", "chrX:7811670-7811868", "-f", "0.001"],
        lambda pos, ref, alt: pos == 7811864 and ref == "C" and alt == "T",
        1,
    ),
    (
        "Zero_in_DUPSV_chr19_291259",
        fa_hg19,
        f"{vj_bams}/Zero_in_DUPSV_TUMOUR.bam|{vj_bams}/Zero_in_DUPSV_NORMAL.bam",
        ["-R", "chr19:290742-291636", "-f", "0.001", "-X", "3"],
        lambda pos, ref, alt: pos == 291259 and ref == "A" and alt == "G",
        1,
    ),
    (
        "StrongLOH_chr21_11071401",
        fa_hg1k,
        f"{vj_bams}/TumorA.StrongLOH.chr21.bam|{vj_bams}/NormalA.StrongLOH.chr21.bam",
        ["-R", "chr21:10996701-11186054", "-Q", "10", "-f", "0.01", "-r", "2", "-F", "0x500", "-V", "0.05", "--nosv"],
        lambda pos, ref, alt: pos in (11071401, 11071410) and len(alt) == len(ref) + 1,
        1,
    ),
]

for name, ref_fa, bam_pair, extra, pred, min_off_hits in cases:
    def_vcf = subprocess.check_output([vdr, "-G", ref_fa, "-b", bam_pair, "--compat=vardictjava", "--vcfall", "-v", *extra], text=True)
    off_vcf = subprocess.check_output([vdr, "-G", ref_fa, "-b", bam_pair, "--compat=off", "--vcfall", "-v", *extra], text=True)
    def_hits = [
        (int(l.split("\t")[1]), l.split("\t")[3], l.split("\t")[4])
        for l in def_vcf.splitlines() if not l.startswith("#") and pred(int(l.split("\t")[1]), l.split("\t")[3], l.split("\t")[4])
    ]
    off_hits = [
        (int(l.split("\t")[1]), l.split("\t")[3], l.split("\t")[4])
        for l in off_vcf.splitlines() if not l.startswith("#") and pred(int(l.split("\t")[1]), l.split("\t")[3], l.split("\t")[4])
    ]
    if def_hits:
        raise AssertionError(f"{name}: expected default --compat=vardictjava to drop secondary allele, got {def_hits}")
    if len(off_hits) < min_off_hits:
        raise AssertionError(f"{name}: expected --compat=off to recover >= {min_off_hits} secondary alleles, got {off_hits}")

# Also check 55-col StrongLOH .var rows keep unpolluted tumor SV_info1 ("0" at col index 52, VarLabel at col index 49)
loh_var = subprocess.check_output(
    [vdr, "-G", fa_hg1k, "-b", f"{vj_bams}/TumorA.StrongLOH.chr21.bam|{vj_bams}/NormalA.StrongLOH.chr21.bam",
     "--compat=off", "--emit-intermediate", "-R", "chr21:10996701-11186054", "-Q", "10", "-f", "0.01", "-r", "2", "-F", "0x500", "-V", "0.05", "--nosv"],
    text=True,
)
strong_loh_checked = 0
for row in loh_var.splitlines():
    cols = row.split("\t")
    if len(cols) == 55 and cols[49] == "StrongLOH":
        strong_loh_checked += 1
        if cols[52] != "0":
            raise AssertionError(f"StrongLOH row polluted tumor SV_info1 column (cols[52]): {cols[52]!r}")
if strong_loh_checked < 1:
    raise AssertionError("Expected >= 1 StrongLOH row in 55-column .var output")

# 3b) Clinical MUC19 paralog somatic pair (AstraZeneca-NGS/VarDictJava#376)
muc_dir = os.path.join(w4, "muc19")
muc_vcf = os.path.join(tmp, "muc19_off.vcf.gz")
subprocess.run(
    [vdr, "-G", os.path.join(muc_dir, "hg19_chr12_muc19.fa"),
     "-b", f"{muc_dir}/MUC19_vaf1.bam|{muc_dir}/MUC19_vaf1_dz.bam",
     "-c", "1", "-S", "2", "-E", "3", "-f", "0.005", "--fisher", "--compat=off", "--vcfall",
     "--output", muc_vcf, "--tbi", os.path.join(muc_dir, "MUC19_vaf1.bed")],
    check=True,
)
subprocess.run(["bcftools", "norm", "--check-ref", "e", "-f", os.path.join(muc_dir, "hg19_chr12_muc19.fa"), muc_vcf],
               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
with gzip.open(muc_vcf, "rt", encoding="utf-8") as fh:
    muc_rows = [l.rstrip("\n").split("\t") for l in fh if not l.startswith("#")]
if len(muc_rows) < 20:
    raise AssertionError(f"MUC19 expected >= 20 variants under --compat=off -f 0.005, got {len(muc_rows)}")
for r in muc_rows:
    fmt = r[8].split(":")
    for smp in (r[9], r[10]):
        vals = dict(zip(fmt, smp.split(":")))
        dp, vd, af = int(vals["DP"]), int(vals["VD"]), float(vals["AF"])
        if vd > dp or af > 1.000001:
            raise AssertionError(f"MUC19 VD > DP / AF > 1 invariant violated at {r[0]}:{r[1]} ({vals})")

# 3c) Quantitative TP_21 >= TP_17 & FP_21 <= FP_17 on SEQC2 HCC1395|BL and GIAB HG008-T|N
def intersect_beds(bed_a, bed_b, out_bed, max_a_lines=None):
    b_map = {}
    with open(bed_b, "r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split("\t")
            b_map.setdefault(p[0], []).append((int(p[1]), int(p[2])))
    for c in b_map:
        b_map[c].sort()
    out_ivs = []
    with open(bed_a, "r", encoding="utf-8") as fh:
        for idx, line in enumerate(fh):
            if max_a_lines is not None and idx >= max_a_lines:
                break
            if line.startswith("#") or not line.strip():
                continue
            p = line.split("\t")
            c, s, e = p[0], int(p[1]), int(p[2])
            for bs, be in b_map.get(c, []):
                if bs >= e:
                    break
                if be > s:
                    out_ivs.append((c, max(s, bs), min(e, be)))
    with open(out_bed, "w", encoding="utf-8") as out:
        for c, s, e in out_ivs:
            if e > s:
                out.write(f"{c}\t{s}\t{e}\n")

def eval_somatic_tp_fp(ref_fa, bam_pair, eval_bed, truth_vcf):
    raw_vcf = os.path.join(tmp, "eval_calls.vcf")
    somatic_vcf = os.path.join(tmp, "eval_somatic.vcf")
    with open(raw_vcf, "w", encoding="utf-8") as out:
        subprocess.run(
            [vdr, "-G", ref_fa, "-b", bam_pair, "-c", "1", "-S", "2", "-E", "3",
             "-f", "0.05", "--fisher", "--compat=off", "-v", eval_bed],
            stdout=out,
            check=True,
        )
    with open(raw_vcf, "r", encoding="utf-8") as fh, open(somatic_vcf, "w", encoding="utf-8") as out:
        for line in fh:
            if line.startswith("#"):
                out.write(line)
                continue
            p = line.split("\t")
            if p[6] == "PASS" and ("STATUS=StrongSomatic" in p[7] or "STATUS=LikelySomatic" in p[7]):
                out.write(line)
    view_calls = subprocess.run(
        ["bcftools", "view", "-T", eval_bed, somatic_vcf],
        capture_output=True,
        text=True,
        check=True,
    )
    norm_calls = subprocess.run(
        ["bcftools", "norm", "-m", "-any", "-f", ref_fa],
        input=view_calls.stdout,
        capture_output=True,
        text=True,
        check=True,
    )
    view_truth = subprocess.run(
        ["bcftools", "view", "-f", ".,PASS", "-T", eval_bed, truth_vcf],
        capture_output=True,
        text=True,
        check=True,
    )
    norm_truth = subprocess.run(
        ["bcftools", "norm", "-m", "-any", "-f", ref_fa],
        input=view_truth.stdout,
        capture_output=True,
        text=True,
        check=True,
    )
    calls = {
        tuple(l.split("\t")[i] for i in (0, 1, 3, 4))
        for l in norm_calls.stdout.splitlines() if not l.startswith("#") and l.strip()
    }
    truth = {
        tuple(l.split("\t")[i] for i in (0, 1, 3, 4))
        for l in norm_truth.stdout.splitlines() if not l.startswith("#") and l.strip()
    }
    if not truth:
        raise AssertionError(f"Loaded 0 truth variants from {truth_vcf} over {eval_bed}")
    return len(calls & truth), len(calls - truth)

seqc2_dir = os.path.join(g4_dir, "seqc2")
seqc2_t = os.path.join(seqc2_dir, "HCC1395.bam")
seqc2_n = os.path.join(seqc2_dir, "HCC1395BL.bam")
if not os.path.isfile(seqc2_t) and os.path.isfile(os.path.join(seqc2_dir, "HCC1395.slice.bam")):
    seqc2_t = os.path.join(seqc2_dir, "HCC1395.slice.bam")
    seqc2_n = os.path.join(seqc2_dir, "HCC1395BL.slice.bam")
if os.path.isfile(seqc2_t):
    seqc2_bed = os.path.join(tmp, "seqc2_hc.bed")
    intersect_beds(os.path.join(seqc2_dir, "slices.bed"), os.path.join(seqc2_dir, "truth", "HC.bed"), seqc2_bed)
    tp_seqc2, fp_seqc2 = eval_somatic_tp_fp(
        grch38_ref,
        f"{seqc2_t}|{seqc2_n}",
        seqc2_bed,
        os.path.join(seqc2_dir, "truth", "truth.vcf.gz"),
    )
    # Pinned 17-switch baseline on SEQC2 HCC1395|HCC1395BL: TP_17 = 43, FP_17 = 134
    if tp_seqc2 < 43 or fp_seqc2 > 134:
        raise AssertionError(f"SEQC2 HCC1395 regression vs 17-switch baseline (43/134): TP={tp_seqc2}, FP={fp_seqc2}")
else:
    raise AssertionError(f"Missing SEQC2 HCC1395 slice at {seqc2_dir}")

hg008_dir = os.path.join(g4_dir, "hg008")
if os.path.isfile(os.path.join(hg008_dir, "HG008-T.slice.bam")):
    hg008_bed = os.path.join(tmp, "hg008_hc40.bed")
    intersect_beds(
        os.path.join(hg008_dir, "slices.bed"),
        os.path.join(hg008_dir, "truth", "HG008-T_somatic_smvar_benchmark_v0.4_GRCh38_all.bed"),
        hg008_bed,
        max_a_lines=40,
    )
    tp_hg008, fp_hg008 = eval_somatic_tp_fp(
        giab_ref,
        f"{hg008_dir}/HG008-T.slice.bam|{hg008_dir}/HG008-N.slice.bam",
        hg008_bed,
        os.path.join(hg008_dir, "truth", "HG008-T_somatic_smvar_benchmark_v0.4_GRCh38_tumorvariants.vcf.gz"),
    )
    # Pinned 17-switch baseline on GIAB HG008-T|N (slices[1..40] x v0.4): TP_17 = 97, FP_17 = 27
    if tp_hg008 < 97 or fp_hg008 > 27:
        raise AssertionError(f"GIAB HG008 regression vs 17-switch baseline (97/27): TP={tp_hg008}, FP={fp_hg008}")
else:
    raise AssertionError(f"Missing GIAB HG008 slice at {hg008_dir}")

print(
    f"  PASS Suite 3: 6/6 multi-allelic pairs recovered, MUC19={len(muc_rows)} rows, "
    f"SEQC2 TP={tp_seqc2}/FP={fp_seqc2} (baseline 43/134), HG008 TP={tp_hg008}/FP={fp_hg008} (baseline 97/27)"
)
PY

# =============================================================================
# Suite 4: CompatProfile Switch 19 (adaptor_all_consensus_dirs)
#          (AstraZeneca-NGS/VarDictJava#240 & AstraZeneca-NGS/VarDictJava#241,
#           AstraZeneca-NGS/VarDict#119, PolinaBevad/VarDictJava@92c62cf6)
# =============================================================================
echo "[4/7] Verifying CompatProfile switch 19 (adaptor_all_consensus_dirs) on patient BRCA2 BAM..."
python3 - "$VARDICTRUST_BIN" "$GRCH38_REF" "$WAVE4_DIR" <<'PY'
import os, subprocess, sys

vdr, grch38_ref, w4 = sys.argv[1:4]
bam = os.path.join(w4, "brca2_adaptor", "adaptor_insertion.bam")
adaptors = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
base_cmd = [vdr, "-G", grch38_ref, "-b", bam, "-R", "chr13:32683877-32684023", "-f", "0.01", "-q", "25", "-r", "2", "-UN", "-v"]

def_vcf = subprocess.check_output([*base_cmd, "--compat=vardictjava", "--adaptor", adaptors], text=True)
off_vcf = subprocess.check_output([*base_cmd, "--compat=off", "--adaptor", adaptors], text=True)
no_ad_vcf = subprocess.check_output([*base_cmd, "--compat=vardictjava"], text=True)

def find_adapter_calls(vcf_text):
    hits = []
    for l in vcf_text.splitlines():
        if l.startswith("#") or not l.strip():
            continue
        p = l.split("\t")
        if int(p[1]) == 32683955 and ("AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" in p[4] or len(p[4]) > 100):
            hits.append((p[0], int(p[1]), p[3], p[4]))
    return hits

def_hits = find_adapter_calls(def_vcf)
off_hits = find_adapter_calls(off_vcf)
if len(def_hits) == 0:
    raise AssertionError("Expected --compat=vardictjava to emit BRCA2 chr13:32683955 TruSeq adapter insertion")
if len(off_hits) != 0:
    raise AssertionError(f"Expected --compat=off to suppress BRCA2 adapter insertion to 0 calls, got {off_hits}")
if not find_adapter_calls(no_ad_vcf):
    raise AssertionError("Expected omitting --adaptor to preserve raw chimeric call")
print(f"  PASS Suite 4: BRCA2 chr13:32683955 adapter chimera emitted in default ({len(def_hits[0][3]) - 1} bp) and suppressed to 0 under --compat=off")
PY

# =============================================================================
# Suite 5: CompatProfile Switch 20 (canonical_sv_vcf_anchoring)
#          (AstraZeneca-NGS/VarDictJava#398, AstraZeneca-NGS/VarDictJava#225,
#           AstraZeneca-NGS/VarDictJava#229, AstraZeneca-NGS/VarDict#111,
#           clintval/vartovcf#107, clintval/vartovcf#192, GIAB HG008 v0.6 SVs)
# =============================================================================
echo "[5/7] Verifying CompatProfile switch 20 (canonical_sv_vcf_anchoring) across 7+ SV BAMs, HG008 v0.6, & clintval/vartovcf..."
python3 - "$VARDICTRUST_BIN" "$FA_HG19" "$GIAB_REF" "$VJ_BAMS" "$G4_DIR" "$WAVE4_DIR" "$WORK_TMP" <<'PY'
import gzip, os, subprocess, sys

vdr, fa_hg19, giab_ref, vj_bams, g4_dir, w4, tmp = sys.argv[1:8]

def verify_sv_vcf_gz(label, ref_fa, vcf_gz, min_svs=1):
    subprocess.run(
        ["bcftools", "norm", "--check-ref", "e", "-f", ref_fa, vcf_gz],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=True,
    )
    sv_count = 0
    with gzip.open(vcf_gz, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            chrom, pos, ref, alt, info = p[0], int(p[1]), p[3], p[4], p[7]
            if pos < 1:
                raise AssertionError(f"{label}: invalid POS < 1 at {chrom}:{pos}")
            if not alt.startswith("<"):
                continue
            sv_count += 1
            if len(ref) != 1 or ref not in ("A", "C", "G", "T"):
                raise AssertionError(f"{label}: symbolic SV {chrom}:{pos} {alt} has non-1bp canonical REF={ref!r}")
            kv = dict(item.split("=", 1) for item in info.split(";") if "=" in item)
            for req_key in ("END", "SVTYPE", "SVLEN"):
                if req_key not in kv:
                    raise AssertionError(f"{label}: symbolic SV {chrom}:{pos} {alt} missing {req_key} in INFO={info}")
            end_pos = int(kv["END"])
            svlen = int(kv["SVLEN"])
            if end_pos <= pos or (end_pos - pos) != abs(svlen):
                raise AssertionError(
                    f"{label}: {chrom}:{pos} {alt} violated END - POS == abs(SVLEN) (POS={pos}, END={end_pos}, SVLEN={svlen})"
                )
            if alt == "<DEL>" and svlen >= 0:
                raise AssertionError(f"{label}: <DEL> at {chrom}:{pos} must have negative SVLEN, got {svlen}")
            if alt in ("<DUP>", "<INV>") and svlen <= 0:
                raise AssertionError(f"{label}: {alt} at {chrom}:{pos} must have positive SVLEN, got {svlen}")
            tbx = subprocess.check_output(["tabix", vcf_gz, f"{chrom}:{pos}-{end_pos}"], text=True)
            if not any(int(tl.split("\t")[1]) == pos and tl.split("\t")[4] == alt for tl in tbx.splitlines() if tl.strip()):
                raise AssertionError(f"{label}: tabix query {chrom}:{pos}-{end_pos} missed {alt}")
    if sv_count < min_svs:
        raise AssertionError(f"{label}: expected >= {min_svs} symbolic SVs, got {sv_count}")
    return sv_count

# 5a) All 8 Colo829 / Zero_in_DUPSV SV BAMs (covering <DEL>, complex &/# <DEL>, <INV>, and all 3 <DUP> origins + paired SV)
sv_cases = [
    ("Colo829_DELDiscReverse", fa_hg19, f"{vj_bams}/Colo829.chr10.DELDiscReverse.bam", ["-R", "chr10:89700380-89713976", "-f", "0.001"]),
    ("Colo829_LGDEL3", fa_hg19, f"{vj_bams}/Colo829.chr10.LGDEL3.bam", ["-R", "chr10:127187158-127201055", "-f", "0.001"]),
    ("Colo829_DELSV3Unhappy", fa_hg19, f"{vj_bams}/Colo829.chr20.DELSV3Unhappy.bam", ["-R", "chr20:54434600-54441423", "-f", "0.001"]),
    ("Colo829_DUPDiscSV", fa_hg19, f"{vj_bams}/Colo829.chr20.DUPDiscSV.bam", ["-R", "chr20:47124929-47131802", "-f", "0.001", "-L", "500"]),
    ("Colo829_INVSV", fa_hg19, f"{vj_bams}/Colo829.chr20.INVSV.bam", ["-R", "chr20:46841554-46845755", "-f", "0.001"]),
    ("Colo829_MarkingDUPSV", fa_hg19, f"{vj_bams}/Colo829.chr20.MarkingDUPSV.bam", ["-R", "chr20:47131954-47132045", "-f", "0.001", "-L", "100"]),
    ("Colo829_CandDel5SV", fa_hg19, f"{vj_bams}/Colo829.chr20.CandDel5SV.bam", ["-R", "chr20:506395-507295", "-f", "0.001", "-L", "100"]),
    ("Zero_in_DUPSV_Paired", fa_hg19, f"{vj_bams}/Zero_in_DUPSV_TUMOUR.bam|{vj_bams}/Zero_in_DUPSV_NORMAL.bam",
     ["-R", "chr19:290742-291636", "-f", "0.001", "-X", "3", "-L", "50", "--fisher"]),
]

total_svs = 0
for name, ref_fa, bam_arg, extra in sv_cases:
    out_gz = os.path.join(tmp, f"{name}.vcf.gz")
    subprocess.run([vdr, "-G", ref_fa, "-b", bam_arg, "--compat=off", "--vcfall", "--output", out_gz, "--tbi", *extra], check=True)
    total_svs += verify_sv_vcf_gz(name, ref_fa, out_gz, min_svs=1)

# Verify --noend suppresses END on symbolic SVs (asserting >= 1 <DEL> row is actually present)
noend_vcf = subprocess.check_output(
    [vdr, "-G", fa_hg19, "-b", f"{vj_bams}/Colo829.chr10.DELDiscReverse.bam", "--compat=off", "-v", "--noend", "-R", "chr10:89700380-89713976", "-f", "0.001"],
    text=True,
)
noend_rows = [l for l in noend_vcf.splitlines() if not l.startswith("#") and l.strip()]
if not any(l.split("\t")[4] == "<DEL>" for l in noend_rows):
    raise AssertionError("Expected >= 1 <DEL> row in --noend SV check")
if any(";END=" in l or l.split("\t")[7].startswith("END=") for l in noend_rows):
    raise AssertionError("Expected --noend to suppress END= on symbolic SV records")

# 5b) GIAB HG008-T|N pancreatic cancer v0.6 somatic SV loci (both single-sample & paired somatic VCFs)
hg008_dir = os.path.join(g4_dir, "hg008")
hg008_sv_bed = os.path.join(tmp, "hg008_v06_sv.bed")
with open(hg008_sv_bed, "w", encoding="utf-8") as fh:
    fh.write("chr3\t101702500\t101707500\nchr10\t123542500\t123544500\nchr13\t39535800\t39539000\n")

for mode_name, bam_arg in [
    ("HG008_v06_single", f"{hg008_dir}/HG008-T.slice.bam"),
    ("HG008_v06_paired", f"{hg008_dir}/HG008-T.slice.bam|{hg008_dir}/HG008-N.slice.bam"),
]:
    out_gz = os.path.join(tmp, f"{mode_name}.vcf.gz")
    subprocess.run(
        [vdr, "-G", giab_ref, "-b", bam_arg, "-c", "1", "-S", "2", "-E", "3",
         "-f", "0.01", "-L", "100", "--fisher", "--compat=off", "--vcfall", "--output", out_gz, "--tbi", hg008_sv_bed],
        check=True,
    )
    total_svs += verify_sv_vcf_gz(mode_name, giab_ref, out_gz, min_svs=2)

# 5c) clintval/vartovcf@3b7bd742 reference suite & contig-start POS == 1 check
v2v_dir = os.path.join(w4, "vartovcf")
v2v_ref = os.path.join(v2v_dir, "reference.fa")
for fn in ("calls.vcf", "calls.contig-ends.vcf"):
    subprocess.run(
        ["bcftools", "norm", "--check-ref", "e", "-f", v2v_ref, os.path.join(v2v_dir, fn)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=True,
    )
v2v_gz = os.path.join(tmp, "vartovcf_sv_contig1.vcf.gz")
subprocess.run(
    [vdr, "-G", v2v_ref, "-b", os.path.join(v2v_dir, "vartovcf_sv_contig1.bam"),
     "-c", "1", "-S", "2", "-E", "3", "-g", "4", "-P", "0", "-f", "0.01", "-r", "2",
     "--compat=off", "--vcfall", "--output", v2v_gz, "--tbi", os.path.join(v2v_dir, "vartovcf_sv_contig1.bed")],
    check=True,
)
verify_sv_vcf_gz("vartovcf_sv_contig1", v2v_ref, v2v_gz, min_svs=2)
with gzip.open(v2v_gz, "rt", encoding="utf-8") as fh:
    got_sv = {}
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        p = line.rstrip("\n").split("\t")
        if p[4].startswith("<"):
            kv = dict(item.split("=", 1) for item in p[7].split(";") if "=" in item)
            got_sv[(p[0], int(p[1]), p[3], p[4])] = (int(kv["END"]), int(kv["SVLEN"]))
if got_sv.get(("chr12", 446, "C", "<DUP>")) != (1446, 1000):
    raise AssertionError(f"clintval/vartovcf chr12:446 <DUP> mismatch: got {got_sv}")
if got_sv.get(("chr13", 1141, "G", "<DEL>")) != (4141, -3000):
    raise AssertionError(f"clintval/vartovcf chr13:1141 <DEL> mismatch: got {got_sv}")
print(f"  PASS Suite 5: {total_svs} clinical/fixture symbolic SVs passed bcftools norm --check-ref e + tabix + clintval/vartovcf exact match")
PY

# =============================================================================
# Suite 6: --amplicon-single-primer (-sp) on QIAseq V3 SPE Tumor & 0.5% ctDNA
#          (AstraZeneca-NGS/VarDictJava#301, AstraZeneca-NGS/VarDictJava#303,
#           AstraZeneca-NGS/VarDictJava#326, AstraZeneca-NGS/VarDict#191;
#           Gong et al. Sci. Data 2022, DOI: 10.1038/s41597-022-01359-6)
# =============================================================================
echo "[6/7] Verifying --amplicon-single-primer (-sp) on QIAseq V3 SPE solid tumor & 0.5% ctDNA..."
python3 - "$VARDICTRUST_BIN" "$GRCH38_REF" "$WAVE4_DIR" "$WORK_TMP" <<'PY'
import bisect, gzip, os, subprocess, sys

vdr, grch38_ref, w4, tmp = sys.argv[1:5]
qgn_bed = os.path.join(w4, "beds", "qiaseq_spe_grch38.amp8.bed")
bam_a = os.path.join(w4, "bam", "SampleA_QGN1.bam")
bam_ac5 = os.path.join(w4, "bam", "SampleAC5_QGN1.bam")

# Verify -sp without -a on a 3-col non-amplicon BED errors cleanly (CmdParser.java:239)
plain_bed = os.path.join(tmp, "plain3.bed")
with open(plain_bed, "w", encoding="utf-8") as fh:
    fh.write("chr1\t11150069\t11150407\n")
bad_sp_rc = subprocess.call(
    [vdr, "-G", grch38_ref, "-b", bam_a, "-c", "1", "-S", "2", "-E", "3", "-sp", "-v", plain_bed],
    stdout=subprocess.DEVNULL,
    stderr=subprocess.DEVNULL,
)
if bad_sp_rc == 0:
    raise AssertionError("Expected -sp without -a on a 3-col BED to exit non-zero")

kp = {}
with gzip.open(os.path.join(w4, "raw", "KnownPositives_hg19ToHg38.vcf.gz"), "rt", encoding="utf-8") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        kp.setdefault(p[0], []).append((int(p[1]), p[3], p[4]))

kn = {}
with gzip.open(os.path.join(w4, "raw", "KnownNegatives_hg38.bed.gz"), "rt", encoding="utf-8") as fh:
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        p = line.rstrip("\n").split("\t")
        kn.setdefault(p[0], []).append((int(p[1]), int(p[2])))
for c in kn:
    kn[c].sort()

def in_kn(c, pos1):
    ivs = kn.get(c, [])
    idx = bisect.bisect_left(ivs, (pos1, -1))
    for j in range(max(0, idx - 2), min(len(ivs), idx + 3)):
        if ivs[j][0] < pos1 <= ivs[j][1]:
            return True
    return False

def analyze_spe_run(bam_path, f_str, with_sp, kn_max_af=1.01):
    cmd = [vdr, "-G", grch38_ref, "-b", bam_path, "-a", "10:0.95", "-UN", "-f", f_str, "-r", "3", "-q", "25", "--compat=off", "--vcfall", "-v"]
    if with_sp:
        cmd.append("-sp")
    cmd.append(qgn_bed)
    vcf = subprocess.check_output(cmd, text=True)
    dps, pass_kp, low_vaf_kp, kn_fps = [], 0, 0, []
    for l in vcf.splitlines():
        if l.startswith("#") or not l.strip():
            continue
        p = l.split("\t")
        c, pos, ra, va, flt, info = p[0], int(p[1]), p[3], p[4], p[6], p[7]
        kv = dict(item.split("=", 1) for item in info.split(";") if "=" in item)
        dp = int(kv["DP"])
        af = float(kv["AF"])
        dps.append(dp)
        if flt == "PASS":
            if any(v[0] == pos and v[2] == va for v in kp.get(c, [])):
                pass_kp += 1
                if af < 0.10:
                    low_vaf_kp += 1
            if af < kn_max_af and in_kn(c, pos):
                kn_fps.append((c, pos, ra, va, dp, af))
    return vcf, sum(dps), (sum(dps) / len(dps) if dps else 0.0), pass_kp, low_vaf_kp, kn_fps

_, sum_no_a, mean_no_a, kp_no_a, _, _ = analyze_spe_run(bam_a, "0.02", False, kn_max_af=1.01)
_, sum_sp_a, mean_sp_a, kp_sp_a, _, kn_fp_a = analyze_spe_run(bam_a, "0.02", True, kn_max_af=1.01)

dp_ratio_a = (sum_sp_a / max(1, sum_no_a)) if sum_no_a > 0 else float(sum_sp_a)
if dp_ratio_a < 3.0 and (mean_sp_a / max(1.0, mean_no_a)) < 3.0:
    raise AssertionError(
        f"SampleA_QGN1 -sp failed >= 3x DP recovery gate: no_sp=(sum={sum_no_a}, mean={mean_no_a:.1f}), "
        f"sp=(sum={sum_sp_a}, mean={mean_sp_a:.1f})"
    )
if kp_sp_a < 30 or kp_sp_a <= kp_no_a:
    raise AssertionError(f"SampleA_QGN1 -sp KnownPositives recovery too low: sp={kp_sp_a} vs no_sp={kp_no_a}")
if kn_fp_a:
    raise AssertionError(f"SampleA_QGN1 -sp produced KnownNegatives PASS FPs: {kn_fp_a[:5]}")

# SampleAC5_QGN1 is 95% Sample C gDNA + 5% Sample A ctDNA dilution: enforce 0 KnownNegatives PASS FPs across all depths in the < 0.10 ctDNA AF band
_, sum_no_ac5, mean_no_ac5, kp_no_ac5, _, _ = analyze_spe_run(bam_ac5, "0.0025", False, kn_max_af=0.10)
vcf_ac5, sum_sp_ac5, mean_sp_ac5, kp_sp_ac5, low_ac5, kn_fp_ac5 = analyze_spe_run(bam_ac5, "0.0025", True, kn_max_af=0.10)
if '##FILTER=<ID=f0.0025,Description="Allele frequency < 0.0025">' not in vcf_ac5:
    raise AssertionError("SampleAC5_QGN1 missing raw f0.0025 FILTER header")
if (sum_sp_ac5 / max(1, sum_no_ac5)) < 3.0 and (mean_sp_ac5 / max(1.0, mean_no_ac5)) < 3.0:
    raise AssertionError(f"SampleAC5_QGN1 -sp failed >= 3x DP recovery gate: {sum_no_ac5} vs {sum_sp_ac5}")
if kp_sp_ac5 < 20 or kp_sp_ac5 <= kp_no_ac5:
    raise AssertionError(f"SampleAC5_QGN1 -sp KnownPositives recovery too low: sp={kp_sp_ac5} vs no_sp={kp_no_ac5}")
if kn_fp_ac5:
    raise AssertionError(f"SampleAC5_QGN1 -sp produced KnownNegatives PASS FPs: {kn_fp_ac5[:5]}")

print(
    f"  PASS Suite 6: QIAseq SPE -sp DP gain A={sum_sp_a / max(1, sum_no_a):.1f}x (mean {mean_no_a:.1f}->{mean_sp_a:.1f}, KP {kp_no_a}->{kp_sp_a}, KN_FP=0), "
    f"ctDNA AC5={sum_sp_ac5 / max(1, sum_no_ac5):.1f}x (KP {kp_no_ac5}->{kp_sp_ac5}, lowVAF={low_ac5}, KN_FP=0)"
)
PY

# =============================================================================
# Suite 7: CompatProfile Switch 21 (paired_somatic_amplicon) Across QIAseq SPE,
#          Illumina PE (BAM + CRAM + TBI + CSI + Primer Mask), Ion Torrent SE,
#          and Colo829 (AstraZeneca-NGS/VarDictJava#378, palakpsheth/VarDictRust#95)
# =============================================================================
echo "[7/7] Verifying CompatProfile switch 21 (paired_somatic_amplicon) across QIAseq SPE, Illumina PE BAM/CRAM, Ion Torrent SE, & Colo829..."
python3 - "$VARDICTRUST_BIN" "$GRCH38_REF" "$FA_HG19" "$VJ_BAMS" "$WAVE4_DIR" "$WORK_TMP" <<'PY'
import gzip, os, subprocess, sys

vdr, grch38_ref, fa_hg19, vj_bams, w4, tmp = sys.argv[1:7]
qgn_bed = os.path.join(w4, "beds", "qiaseq_spe_grch38.amp8.bed")
ilm_bed = os.path.join(w4, "beds", "ilm_twoprimer_grch38.amp8.bed")
ilm_primer_bed = os.path.join(w4, "beds", "ilm_primer_overlap_grch38.amp8.bed")
tfs_bed = os.path.join(w4, "beds", "tfs_iontorrent_grch38.amp8.bed")

# 7a) Default --compat=vardictjava must still reject paired + amplicon VCF (AstraZeneca-NGS/VarDictJava#378)
def_rc = subprocess.call(
    [vdr, "-G", grch38_ref, "-b", f"{w4}/bam/SampleA_ILM1.bam|{w4}/bam/SampleB_ILM1.bam",
     "-a", "40:0.70", "--compat=vardictjava", "-v", ilm_bed],
    stdout=subprocess.DEVNULL,
    stderr=subprocess.DEVNULL,
)
if def_rc == 0:
    raise AssertionError("Expected default --compat=vardictjava to reject paired somatic amplicon VCF")

kp_loci = set()
with gzip.open(os.path.join(w4, "raw", "KnownPositives_hg19ToHg38.vcf.gz"), "rt", encoding="utf-8") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        kp_loci.add((p[0], int(p[1]), p[4]))

def verify_paired_amp_vcf(label, ref_fa, vcf_gz):
    subprocess.run(
        ["bcftools", "norm", "--check-ref", "e", "-f", ref_fa, vcf_gz],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=True,
    )
    headers = []
    rows = []
    with gzip.open(vcf_gz, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                headers.append(line.strip())
            elif line.strip():
                rows.append(line.rstrip("\n").split("\t"))
    hdr_blob = "\n".join(headers)
    for req_hdr in ("##INFO=<ID=GDAMP,", "##INFO=<ID=TLAMP,", "##INFO=<ID=NCAMP,", "##INFO=<ID=AMPFLAG,", "##FILTER=<ID=AMPBIAS,"):
        if req_hdr not in hdr_blob:
            raise AssertionError(f"{label}: missing required paired-amplicon header {req_hdr}")

    somatic_kp, germline_n, tlamp2_n, ampbias_n = 0, 0, 0, 0
    for r in rows:
        chrom, pos, _id, _ref, alt, _q, flt, info, fmt_s, t_s, n_s = r[:11]
        pos_i = int(pos)
        if "pSTD" in flt.split(";"):
            raise AssertionError(f"{label}: pSTD filter must be suppressed in paired somatic amplicon mode ({chrom}:{pos})")
        kv = dict(item.split("=", 1) for item in info.split(";") if "=" in item)
        for req_k in ("STATUS", "GDAMP", "TLAMP", "NCAMP", "AMPFLAG"):
            if req_k not in kv:
                raise AssertionError(f"{label}: missing {req_k} in INFO at {chrom}:{pos}: {info}")
        gdamp, tlamp, ncamp, ampflag = int(kv["GDAMP"]), int(kv["TLAMP"]), int(kv["NCAMP"]), int(kv["AMPFLAG"])
        if tlamp >= 2:
            tlamp2_n += 1
        # Ports var2vcf_valid.pl:185: (($gamp < $tamp - $ncamp) || $ampflag), evaluated outside normal_clean
        is_amp_discordant = (gdamp < (tlamp - ncamp)) or (ampflag != 0)
        flt_tags = flt.split(";")
        if is_amp_discordant and (flt == "PASS" or "AMPBIAS" not in flt_tags):
            raise AssertionError(
                f"{label}: AMPBIAS discordant call {chrom}:{pos} (GDAMP={gdamp}, TLAMP={tlamp}, NCAMP={ncamp}, AMPFLAG={ampflag}) missing AMPBIAS: FILTER={flt}"
            )
        if "AMPBIAS" in flt_tags:
            ampbias_n += 1
        if kv["STATUS"] == "StrongSomatic" and (chrom, pos_i, alt) in kp_loci:
            somatic_kp += 1
        if kv["STATUS"] == "Germline":
            germline_n += 1
        fmt = fmt_s.split(":")
        for smp in (t_s, n_s):
            vals = dict(zip(fmt, smp.split(":")))
            dp, vd, af = int(vals["DP"]), int(vals["VD"]), float(vals["AF"])
            if vd > dp or af > 1.000001:
                raise AssertionError(f"{label}: VD > DP / AF > 1 at {chrom}:{pos} ({vals})")
    if rows:
        first_c, first_p = rows[0][0], int(rows[0][1])
        tbx = subprocess.check_output(["tabix", vcf_gz, f"{first_c}:{first_p}-{first_p}"], text=True)
        if not tbx.strip():
            raise AssertionError(f"{label}: tabix query failed at {first_c}:{first_p}")
    return rows, somatic_kp, germline_n, tlamp2_n, ampbias_n

def verify_var_cols(label, ref_fa, bam_pair, bed_path, amp_args):
    for want_cols, extra in [(55, []), (61, ["--fisher"])]:
        var_out = subprocess.check_output(
            [vdr, "-G", ref_fa, "-b", bam_pair, *amp_args, "--compat=off", "--emit-intermediate", "-h", *extra, bed_path],
            text=True,
        )
        lines = [l for l in var_out.splitlines() if l.strip()]
        if len(lines) < 2:
            raise AssertionError(f"{label}: expected header + >= 1 .var row, got {len(lines)} lines")
        hdr_cols = len(lines[0].split("\t"))
        if hdr_cols != 55:
            raise AssertionError(f"{label}: -h header line expected 55 columns (SomaticMode.java:134-150), got {hdr_cols}")
        for idx, l in enumerate(lines[1:], start=2):
            ncols = len(l.split("\t"))
            if ncols != want_cols:
                raise AssertionError(f"{label}: data line {idx} expected {want_cols} columns, got {ncols}")

# 7b) QIAseq V3 SPE Paired Somatic (SampleA_QGN1 | SampleB_QGN1)
qgn_gz = os.path.join(tmp, "qgn_paired.vcf.gz")
subprocess.run(
    [vdr, "-G", grch38_ref, "-b", f"{w4}/bam/SampleA_QGN1.bam|{w4}/bam/SampleB_QGN1.bam",
     "-N", "SampleA_QGN1|SampleB_QGN1",
     "-a", "10:0.95", "-sp", "-UN", "-f", "0.02", "-r", "3", "-q", "25",
     "--compat=off", "--vcfall", "--fisher", "--output", qgn_gz, "--tbi", qgn_bed],
    check=True,
)
qgn_rows, qgn_som_kp, qgn_germ, _, _ = verify_paired_amp_vcf("QIAseq_SPE_Paired", grch38_ref, qgn_gz)
if qgn_som_kp < 5 or qgn_germ < 5:
    raise AssertionError(f"QIAseq_SPE_Paired expected StrongSomatic KP >= 5 & Germline >= 5, got som={qgn_som_kp}, germ={qgn_germ}")

# 7c) Illumina Paired-End Two-Primer (BAM vs CRAM + .tbi + .csi + 55/61 .var + primer-overlap suppression)
ilm_bam_gz = os.path.join(tmp, "ilm_bam.vcf.gz")
ilm_cram_gz = os.path.join(tmp, "ilm_cram.vcf.gz")
subprocess.run(
    [vdr, "-G", grch38_ref, "-b", f"{w4}/bam/SampleA_ILM1.bam|{w4}/bam/SampleB_ILM1.bam",
     "-N", "SampleA_ILM1|SampleB_ILM1",
     "-a", "40:0.70", "-f", "0.02", "-q", "25", "--compat=off", "--vcfall", "--fisher",
     "--output", ilm_bam_gz, "--tbi", ilm_bed],
    check=True,
)
subprocess.run(
    [vdr, "-G", grch38_ref, "-b", f"{w4}/bam/SampleA_ILM1.cram|{w4}/bam/SampleB_ILM1.cram",
     "-N", "SampleA_ILM1|SampleB_ILM1",
     "-a", "40:0.70", "-f", "0.02", "-q", "25", "--compat=off", "--vcfall", "--fisher",
     "--output", ilm_cram_gz, "--csi", ilm_bed],
    check=True,
)
ilm_rows, ilm_som_kp, ilm_germ, ilm_tlamp2, ilm_ampbias = verify_paired_amp_vcf("Illumina_PE_BAM", grch38_ref, ilm_bam_gz)
cram_rows, _, _, _, _ = verify_paired_amp_vcf("Illumina_PE_CRAM", grch38_ref, ilm_cram_gz)
if [(r[0], r[1], r[3], r[4], r[6], r[7]) for r in ilm_rows] != [(r[0], r[1], r[3], r[4], r[6], r[7]) for r in cram_rows]:
    raise AssertionError("Illumina paired amplicon BAM vs CRAM VCF output diverged")
if ilm_som_kp < 3 or ilm_germ < 3 or ilm_tlamp2 < 5:
    raise AssertionError(
        f"Illumina_PE expected StrongSomatic KP >= 3, Germline >= 3, TLAMP>=2 >= 5; "
        f"got som={ilm_som_kp}, germ={ilm_germ}, tlamp2={ilm_tlamp2}"
    )
verify_var_cols(
    "Illumina_PE_VAR",
    grch38_ref,
    f"{w4}/bam/SampleA_ILM1.bam|{w4}/bam/SampleB_ILM1.bam",
    ilm_bed,
    ["-a", "40:0.70", "-f", "0.02", "-q", "25"],
)

# Verify primer-overlap suppression on ilm_primer_overlap_grch38.amp8.bed
with open(ilm_primer_bed, "r", encoding="utf-8") as fh:
    p_cols = fh.readline().strip().split("\t")
p_chr, p_st, p_ins_st = p_cols[0], int(p_cols[1]), int(p_cols[6])
primer_vcf = subprocess.check_output(
    [vdr, "-G", grch38_ref, "-b", f"{w4}/bam/SampleA_ILM1.bam|{w4}/bam/SampleB_ILM1.bam",
     "-a", "40:0.70", "-f", "0.01", "-q", "25", "--compat=off", "--vcfall", "-v", ilm_primer_bed],
    text=True,
)
leaked_primer = [
    l for l in primer_vcf.splitlines()
    if not l.startswith("#") and l.strip() and int(l.split("\t")[1]) <= p_ins_st
]
if leaked_primer:
    raise AssertionError(f"Primer-overlapping variant leaked into paired somatic amplicon VCF: {leaked_primer[0]}")

# 7d) Ion Torrent S5 XL Single-End Homopolymer Amplicon (SampleA_TFS1 | SampleB_TFS1)
tfs_gz = os.path.join(tmp, "tfs_paired.vcf.gz")
subprocess.run(
    [vdr, "-G", grch38_ref, "-b", f"{w4}/bam/SampleA_TFS1.bam|{w4}/bam/SampleB_TFS1.bam",
     "-a", "10:0.95", "-sp", "-f", "0.02", "-q", "20", "--compat=off", "--vcfall", "--fisher",
     "--output", tfs_gz, "--tbi", tfs_bed],
    check=True,
)
tfs_rows, _, _, _, _ = verify_paired_amp_vcf("IonTorrent_SE_Paired", grch38_ref, tfs_gz)
verify_var_cols(
    "IonTorrent_SE_VAR",
    grch38_ref,
    f"{w4}/bam/SampleA_TFS1.bam|{w4}/bam/SampleB_TFS1.bam",
    tfs_bed,
    ["-a", "10:0.95", "-sp", "-f", "0.02", "-q", "20"],
)

# 7e) Colo829-18 | Colo829-19 Overlapping 8-Col Amplicon Pair
colo_amp_bed = os.path.join(tmp, "colo829_amp8.bed")
with open(colo_amp_bed, "w", encoding="utf-8") as fh:
    fh.write("chr1\t756240\t756380\tCOLO_AMP_1\t0\t+\t756260\t756360\n")
    fh.write("chr1\t756255\t756395\tCOLO_AMP_2\t0\t+\t756275\t756375\n")
colo_gz = os.path.join(tmp, "colo_amp.vcf.gz")
subprocess.run(
    [vdr, "-G", fa_hg19, "-b", f"{vj_bams}/Colo829-18_S3-sort.bam|{vj_bams}/Colo829-19_S4-sort.bam",
     "-a", "40:0.60", "-f", "0.01", "--compat=off", "--vcfall", "--fisher",
     "--output", colo_gz, "--tbi", colo_amp_bed],
    check=True,
)
colo_rows, _, _, _, _ = verify_paired_amp_vcf("Colo829_Paired_Amplicon", fa_hg19, colo_gz)

print(
    f"  PASS Suite 7: QIAseq_SPE={len(qgn_rows)} rows, Illumina_PE={len(ilm_rows)} rows "
    f"(BAM==CRAM, TLAMP>=2={ilm_tlamp2}, AMPBIAS={ilm_ampbias}, primer-mask=0 leaks), "
    f"IonTorrent_SE={len(tfs_rows)} rows, Colo829={len(colo_rows)} rows"
)
PY

T1_MS="$(now_ms)"
ELAPSED_S=$(( (T1_MS - T0_MS) / 1000 ))
echo "================================================================================"
echo " PASS verify_wave4_public_data: all 7 Wave-4 clinical suites green (${ELAPSED_S}s)"
echo "================================================================================"
