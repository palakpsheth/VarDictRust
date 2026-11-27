#!/usr/bin/env bash
# =============================================================================
# harness/verify_presets_public_data.sh
# Comprehensive 5-gate real-world clinical validation suite for v1.1.0
# --preset wgs|wes|amplicon|ctdna (Issue #94, AR-01..AR-06, AR-10..AR-12).
#
# Clinical Cohorts & Provenance Registry:
# - Gate 1 (--preset wgs): NIST GIAB HG002 germline WGS (chr20) + GIAB HG008
#   PDAC paired T/N somatic WGS (--tbi & --csi, 3-col 0-based BED + 4-col
#   3'_UTR BED annotation AR-01 + paired Q20 single-tag check AR-04).
# - Gate 2 (--preset wes): FDA SEQC2 HCC1395|HCC1395BL paired T/N somatic WES
#   (Fang et al., Nat. Biotechnol. 2021, DOI: 10.1038/s41587-021-00993-6) +
#   UCSC 8-column CDS BED8 non-hijack verification (AR-02).
# - Gate 3 (--preset amplicon): Targeted 8-col PCR amplicon panel, paired
#   somatic Illumina PE two-primer amplicon, & QIAseq V3 single-primer (-sp)
#   solid tumor panel (Gong et al., Sci. Data 2022, DOI: 10.1038/s41597-022-01359-6).
# - Gate 4 (--preset ctdna): Deterministic ~1,000x real-read 0.5% VAF QIAseq
#   ctDNA titration (SampleB_QGN1 wild-type background + SampleA_QGN1 KnownPositive
#   chr12:120999311 G>A spike-in, AR-06) + SampleAC5_QGN1 d100 soft-filtering &
#   explicit --vcf-min-depth 20 override precedence.
# - Gate 5 (-R 1-Based Coordinate Invariance & Negative CLI Guards):
#   --preset wgs -f 0.01 -R vs --compat=off -f 0.01 -R coordinate identity
#   (AR-03, AR-10) + loud non-zero rejection of invalid/conflicting combinations
#   (AR-11, AR-12, Issue #96).
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

for _req_tool in samtools bcftools tabix python3 gzip; do
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

if [ ! -f "$WAVE4_DIR/.done" ]; then
    DATA_DIR="$G4_DIR" WAVE4_DIR="$WAVE4_DIR" bash "$SCRIPT_DIR/get_wave4_public_slices.sh"
fi

FA_HG19="$DATA_ROOT/VarDictJava/testdata/fastas/hg19.fa"
if [ ! -f "$FA_HG19" ]; then
    python3 "$ROOT_DIR/harness/convert_fastas.py"
fi

GRCH38_REF="$G4_DIR/ref/d1vd1/GRCh38.d1.vd1.fa"
GIAB_REF="$G4_DIR/ref/giabv3/GRCh38_GIABv3.fa"
if [ ! -f "$GIAB_REF" ] && [ -f "$G4_DIR/ref/giabv3/GRCh38_GIABv3.fasta" ]; then
    GIAB_REF="$G4_DIR/ref/giabv3/GRCh38_GIABv3.fasta"
fi

SEQ_T="$G4_DIR/seqc2/HCC1395.bam"
SEQ_N="$G4_DIR/seqc2/HCC1395BL.bam"
if [ ! -f "$SEQ_T" ] && [ -f "$G4_DIR/seqc2/HCC1395.slice.bam" ]; then
    SEQ_T="$G4_DIR/seqc2/HCC1395.slice.bam"
    SEQ_N="$G4_DIR/seqc2/HCC1395BL.slice.bam"
fi

for _req_file in "$FA_HG19" "${FA_HG19}.fai" \
                 "$GRCH38_REF" "${GRCH38_REF}.fai" \
                 "$GIAB_REF" "${GIAB_REF}.fai" \
                 "$G4_DIR/hg002/HG002.chr20.bam" \
                 "$G4_DIR/hg008/HG008-T.slice.bam" "$G4_DIR/hg008/HG008-N.slice.bam" \
                 "$SEQ_T" "$SEQ_N" \
                 "$G4_DIR/amplicon/amplicon.bam" "$G4_DIR/amplicon/correct_8columns_amplicon.bed" \
                 "$WAVE4_DIR/bam/SampleA_QGN1.bam" "$WAVE4_DIR/bam/SampleB_QGN1.bam" \
                 "$WAVE4_DIR/bam/SampleAC5_QGN1.bam" "$WAVE4_DIR/bam/SampleA_ILM1.bam" \
                 "$WAVE4_DIR/bam/SampleB_ILM1.bam" "$WAVE4_DIR/beds/qiaseq_spe_grch38.amp8.bed" \
                 "$WAVE4_DIR/beds/ilm_twoprimer_grch38.amp8.bed"; do
    if [ ! -s "$_req_file" ]; then
        echo "ERROR: missing clinical preset verification prerequisite: '$_req_file'" >&2
        exit 1
    fi
done
unset _req_file

_TMP_BASE="${TMPDIR:-/tmp}"
WORK_TMP="$(mktemp -d "${_TMP_BASE%/}/vdr_presets_verify.XXXXXX")"
unset _TMP_BASE
trap 'rm -rf "$WORK_TMP"' EXIT

T0_MS="$(now_ms)"
echo "================================================================================"
echo " Clinical Workflow Presets (--preset wgs|wes|amplicon|ctdna) 5-Gate Suite"
echo " Binary: $VARDICTRUST_BIN"
echo " Data:   $G4_DIR"
echo "================================================================================"

# =============================================================================
# Gate 1: --preset wgs on GIAB HG002 Germline + GIAB HG008 Paired T/N Somatic
#         (3-col 0-based BED, --tbi & --csi, 3'_UTR 4-col BED AR-01, Q20 AR-04)
# =============================================================================
echo "[1/5] Verifying --preset wgs on GIAB HG002 germline & GIAB HG008 paired somatic..."
python3 - "$VARDICTRUST_BIN" "$GIAB_REF" "$G4_DIR" "$WORK_TMP" <<'PY'
import gzip, os, subprocess, sys

vdr, giab_ref, g4_dir, tmp = sys.argv[1:5]

def check_preset_vcf_gz(vcf_gz, ref_fa, preset_name, expected_filters, query_region):
    subprocess.run(
        ["bcftools", "norm", "--check-ref", "e", "-f", ref_fa, vcf_gz],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=True,
    )
    tbx_out = subprocess.check_output(["tabix", vcf_gz, query_region], text=True)
    with gzip.open(vcf_gz, "rt", encoding="utf-8") as fh:
        all_lines = [l.rstrip("\n") for l in fh]
    hdr_lines = [l for l in all_lines if l.startswith("#")]
    rec_lines = [l for l in all_lines if l and not l.startswith("#")]
    if len(hdr_lines) < 5:
        raise AssertionError(f"{vcf_gz}: truncated VCF header")
    if hdr_lines[2] != "##VarDictRustCompat=off":
        raise AssertionError(f"{vcf_gz}: expected line 3 '##VarDictRustCompat=off', got {hdr_lines[2]!r}")
    if hdr_lines[3] != f"##VarDictRustPreset={preset_name}":
        raise AssertionError(f"{vcf_gz}: expected line 4 '##VarDictRustPreset={preset_name}', got {hdr_lines[3]!r}")
    hdr_blob = "\n".join(hdr_lines)
    for flt_id in expected_filters:
        if f"##FILTER=<ID={flt_id}," not in hdr_blob:
            raise AssertionError(f"{vcf_gz}: missing expected FILTER header ID={flt_id}")
    return hdr_lines, rec_lines, tbx_out

# 1a) GIAB HG002 germline on 3-column 0-based BED (NO -c/-S/-E/-z passed!) + --tbi
hg002_bed3 = os.path.join(tmp, "hg002_wgs_3col.bed")
with open(hg002_bed3, "w", encoding="utf-8") as fh:
    fh.write("chr20\t10000000\t10060000\n")

hg002_gz = os.path.join(tmp, "hg002_wgs.vcf.gz")
subprocess.run(
    [vdr, "--preset", "wgs", "-G", giab_ref, "-b", os.path.join(g4_dir, "hg002", "HG002.chr20.bam"),
     "-N", "HG002", "--output", hg002_gz, "--tbi", "-th", "4", hg002_bed3],
    check=True,
)
_, hg002_recs, hg002_tbx = check_preset_vcf_gz(
    hg002_gz, giab_ref, "wgs", ("f0.05", "d10", "v3", "Q20"), "chr20:10000001-10060000"
)
if len(hg002_recs) < 20 or not hg002_tbx.strip():
    raise AssertionError(f"HG002 --preset wgs expected >= 20 records and non-empty tabix query, got {len(hg002_recs)}")

# 1b) GIAB HG008 paired T/N somatic on 3-column 0-based BED (NO -c/-S/-E/-z passed!) + --csi
hg008_bed3 = os.path.join(tmp, "hg008_wgs_3col.bed")
with open(os.path.join(g4_dir, "hg008", "test5.bed"), "r", encoding="utf-8") as src, \
     open(hg008_bed3, "w", encoding="utf-8") as dst:
    for line in src:
        if line.strip() and not line.startswith("#"):
            p = line.rstrip("\n").split("\t")
            dst.write(f"{p[0]}\t{p[1]}\t{p[2]}\n")

hg008_pair = f"{g4_dir}/hg008/HG008-T.slice.bam|{g4_dir}/hg008/HG008-N.slice.bam"
hg008_gz = os.path.join(tmp, "hg008_wgs.vcf.gz")
subprocess.run(
    [vdr, "--preset", "wgs", "-G", giab_ref, "-b", hg008_pair, "-N", "HG008_T|HG008_N",
     "--fisher", "--output", hg008_gz, "--csi", "-th", "4", hg008_bed3],
    check=True,
)
_, hg008_recs, hg008_tbx = check_preset_vcf_gz(
    hg008_gz, giab_ref, "wgs", ("f0.05", "d10", "v3", "Q20"), "chr1:72271447-72291446"
)
if len(hg008_recs) < 10 or not hg008_tbx.strip():
    raise AssertionError(f"HG008 --preset wgs expected >= 10 records and non-empty CSI tabix query, got {len(hg008_recs)}")

for r in hg008_recs:
    flt = r.split("\t")[6]
    flt_tags = flt.split(";")
    if len(flt_tags) != len(set(flt_tags)) or "Q20;Q20" in flt:
        raise AssertionError(f"Duplicate FILTER tag (AR-04) in HG008 paired somatic record: {r}")

# Also stress-test AR-04 directly on a real-read HG008 slice where tumor MAPQ is lowered to 5 (< 10)
# while normal MAPQ stays 60 (> 10), triggering both mapq_low (< 20) and mapq_snv_low (mapq1 < 10 && eff_mapq > 10 && SNV)!
lowmq_reg = "chr1:72271447-72276000"
lowmq_t_sam = os.path.join(tmp, "hg008_lowmq_t.sam")
lowmq_t_bam = os.path.join(tmp, "hg008_lowmq_t.bam")
sam_t = subprocess.check_output(["samtools", "view", "-h", f"{g4_dir}/hg008/HG008-T.slice.bam", lowmq_reg], text=True)
with open(lowmq_t_sam, "w", encoding="utf-8") as fh:
    for line in sam_t.splitlines():
        if line.startswith("@"):
            fh.write(line + "\n")
        else:
            cols = line.split("\t")
            cols[4] = "5"
            fh.write("\t".join(cols) + "\n")
subprocess.run(["samtools", "view", "-b", "-o", lowmq_t_bam, lowmq_t_sam], check=True)
subprocess.run(["samtools", "index", lowmq_t_bam], check=True)

lowmq_bed3 = os.path.join(tmp, "hg008_lowmq.bed")
with open(lowmq_bed3, "w", encoding="utf-8") as fh:
    fh.write("chr1\t72271446\t72276000\n")
lowmq_gz = os.path.join(tmp, "hg008_lowmq_wgs.vcf.gz")
subprocess.run(
    [vdr, "--preset", "wgs", "-G", giab_ref, "-b", f"{lowmq_t_bam}|{g4_dir}/hg008/HG008-N.slice.bam",
     "-N", "TUMOR|NORMAL", "--output", lowmq_gz, "--tbi", lowmq_bed3],
    check=True,
)
_, lowmq_recs, _ = check_preset_vcf_gz(
    lowmq_gz, giab_ref, "wgs", ("f0.05", "d10", "v3", "Q20"), lowmq_reg
)
q20_seen = 0
for r in lowmq_recs:
    flt = r.split("\t")[6]
    tags = flt.split(";")
    if "Q20;Q20" in flt or tags.count("Q20") > 1:
        raise AssertionError(f"AR-04 regression: duplicate Q20;Q20 emitted on low-MAPQ paired SNV: {r}")
    if "Q20" in tags:
        q20_seen += 1
if q20_seen < 1:
    raise AssertionError("Expected >= 1 Q20-filtered record in low-MAPQ paired somatic stress check")

# 1c) 4-column BED with 3'_UTR gene annotation (AR-01: column_for_chromosome=Some(0) prevents thick_start=3 parse error)
utr_bed4 = os.path.join(tmp, "hg002_utr_4col.bed")
with open(utr_bed4, "w", encoding="utf-8") as fh:
    fh.write("chr20\t10000000\t10030000\t3'_UTR\n")
utr_var = subprocess.check_output(
    [vdr, "--preset", "wgs", "-G", giab_ref, "-b", os.path.join(g4_dir, "hg002", "HG002.chr20.bam"),
     "--emit-intermediate", utr_bed4],
    text=True,
)
utr_rows = [l.split("\t") for l in utr_var.splitlines() if l.strip()]
if not utr_rows or any(r[1] != "3'_UTR" for r in utr_rows):
    raise AssertionError(f"4-col 3'_UTR BED under --preset wgs failed to preserve gene column: rows={len(utr_rows)}")

print(
    f"[PASS] Gate 1/5 (--preset wgs): HG002={len(hg002_recs)} calls (.tbi), HG008={len(hg008_recs)} calls (.csi, "
    f"low-MQ Q20={q20_seen}, 0 duplicate Q20;Q20), 4-col 3'_UTR={len(utr_rows)} rows"
)
PY

# =============================================================================
# Gate 2: --preset wes on FDA SEQC2 HCC1395 T/N Somatic + UCSC BED8 Non-Hijack
#         (AR-02: 8-col UCSC CDS BED8 must NOT trigger PCR amplicon mode!)
# =============================================================================
echo "[2/5] Verifying --preset wes on SEQC2 HCC1395|HCC1395BL & UCSC 8-col CDS BED8 non-hijack..."
python3 - "$VARDICTRUST_BIN" "$GRCH38_REF" "$SEQ_T" "$SEQ_N" "$WORK_TMP" <<'PY'
import gzip, os, subprocess, sys

vdr, grch38_ref, seq_t, seq_n, tmp = sys.argv[1:6]
bam_pair = f"{seq_t}|{seq_n}"

wes_bed4 = os.path.join(tmp, "seqc2_wes_4col.bed")
with open(wes_bed4, "w", encoding="utf-8") as fh:
    fh.write("chr11\t6400000\t6459439\tslice_1\nchr17\t74747665\t74800000\tslice_3\n")

wes_gz = os.path.join(tmp, "seqc2_wes_4col.vcf.gz")
subprocess.run(
    [vdr, "--preset", "wes", "-G", grch38_ref, "-b", bam_pair, "-N", "HCC1395|HCC1395BL",
     "--fisher", "--output", wes_gz, "--tbi", "-th", "4", wes_bed4],
    check=True,
)
subprocess.run(
    ["bcftools", "norm", "--check-ref", "e", "-f", grch38_ref, wes_gz],
    stdout=subprocess.DEVNULL,
    stderr=subprocess.DEVNULL,
    check=True,
)

with gzip.open(wes_gz, "rt", encoding="utf-8") as fh:
    wes_lines = [l.rstrip("\n") for l in fh]
wes_hdr = [l for l in wes_lines if l.startswith("#")]
wes_recs = [l for l in wes_lines if l and not l.startswith("#")]
if wes_hdr[2] != "##VarDictRustCompat=off" or wes_hdr[3] != "##VarDictRustPreset=wes":
    raise AssertionError(f"Unexpected --preset wes header prefix: {wes_hdr[:5]}")
hdr_blob = "\n".join(wes_hdr)
for req_flt in ("f0.02", "d20", "v4", "Q20"):
    if f"##FILTER=<ID={req_flt}," not in hdr_blob:
        raise AssertionError(f"--preset wes missing FILTER header {req_flt}")
if len(wes_recs) < 15:
    raise AssertionError(f"Expected >= 15 SEQC2 WES somatic records, got {len(wes_recs)}")

# UCSC 8-column CDS BED8 (thickStart >= start && thickEnd <= end, which in default mode auto-triggers
# 8-col PCR amplicon mode and drops sheared capture reads!): under --preset wes WITHOUT -a (AR-02),
# it must stay in somatic capture mode with 0 AMPBIAS/GDAMP and call identical somatic loci!
ucsc_bed8 = os.path.join(tmp, "seqc2_ucsc_cds.bed8")
with open(ucsc_bed8, "w", encoding="utf-8") as fh:
    fh.write("chr11\t6400000\t6459439\tslice_1\t0\t+\t6400500\t6459000\n")
    fh.write("chr17\t74747665\t74800000\tslice_3\t0\t+\t74748000\t74799500\n")

bed8_gz = os.path.join(tmp, "seqc2_wes_bed8.vcf.gz")
subprocess.run(
    [vdr, "--preset", "wes", "-G", grch38_ref, "-b", bam_pair, "-N", "HCC1395|HCC1395BL",
     "--fisher", "--output", bed8_gz, "--tbi", "-th", "4", ucsc_bed8],
    check=True,
)
subprocess.run(
    ["bcftools", "norm", "--check-ref", "e", "-f", grch38_ref, bed8_gz],
    stdout=subprocess.DEVNULL,
    stderr=subprocess.DEVNULL,
    check=True,
)
with gzip.open(bed8_gz, "rt", encoding="utf-8") as fh:
    bed8_all = [l.rstrip("\n") for l in fh]
bed8_hdr_blob = "\n".join(l for l in bed8_all if l.startswith("#"))
bed8_recs = [l for l in bed8_all if l and not l.startswith("#")]

if "##FILTER=<ID=AMPBIAS," in bed8_hdr_blob or "##INFO=<ID=GDAMP," in bed8_hdr_blob:
    raise AssertionError("AR-02 regression: UCSC 8-col BED8 hijacked --preset wes into PCR Amplicon VCF schema!")
for r in bed8_recs:
    cols = r.split("\t")
    if "AMPBIAS" in cols[6] or "GDAMP=" in cols[7] or "TLAMP=" in cols[7]:
        raise AssertionError(f"AR-02 regression: amplicon tag leaked into --preset wes BED8 record: {r}")

sites_4col = {(r.split("\t")[0], int(r.split("\t")[1]), r.split("\t")[3], r.split("\t")[4]) for r in wes_recs}
sites_bed8 = {(r.split("\t")[0], int(r.split("\t")[1]), r.split("\t")[3], r.split("\t")[4]) for r in bed8_recs}
if sites_4col != sites_bed8:
    raise AssertionError(
        f"AR-02 regression: UCSC BED8 dropped/altered capture calls vs 4-col BED ({len(sites_bed8)} vs {len(sites_4col)})"
    )

print(
    f"[PASS] Gate 2/5 (--preset wes): SEQC2 HCC1395={len(wes_recs)} somatic records; "
    f"UCSC 8-col BED8 non-hijack matched {len(sites_bed8)}/{len(sites_4col)} calls (0 AMPBIAS/GDAMP)"
)
PY

# =============================================================================
# Gate 3: --preset amplicon on 8-Col Panel, Paired T/N Amplicon, & QIAseq -sp
#         (NO -a and NO -z required!)
# =============================================================================
echo "[3/5] Verifying --preset amplicon on single-sample 8-col, paired Illumina PE, & QIAseq V3 -sp..."
python3 - "$VARDICTRUST_BIN" "$FA_HG19" "$GRCH38_REF" "$G4_DIR" "$WAVE4_DIR" "$WORK_TMP" <<'PY'
import gzip, os, subprocess, sys

vdr, fa_hg19, grch38_ref, g4_dir, w4_dir, tmp = sys.argv[1:7]

def audit_amplicon_preset(label, ref_fa, bam_arg, bed_path, extra_args):
    out_gz = os.path.join(tmp, f"{label}.vcf.gz")
    subprocess.run(
        [vdr, "--preset", "amplicon", *extra_args, "-G", ref_fa, "-b", bam_arg,
         "--output", out_gz, "--tbi", "-th", "4", bed_path],
        check=True,
    )
    subprocess.run(
        ["bcftools", "norm", "--check-ref", "e", "-f", ref_fa, out_gz],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=True,
    )
    with gzip.open(out_gz, "rt", encoding="utf-8") as fh:
        lines = [l.rstrip("\n") for l in fh]
    hdrs = [l for l in lines if l.startswith("#")]
    recs = [l for l in lines if l and not l.startswith("#")]
    if hdrs[2] != "##VarDictRustCompat=off" or hdrs[3] != "##VarDictRustPreset=amplicon":
        raise AssertionError(f"{label}: unexpected header lines 3..4: {hdrs[:5]}")
    hdr_blob = "\n".join(hdrs)
    for req in (
        "##FILTER=<ID=f0.01,",
        "##FILTER=<ID=d50,",
        "##FILTER=<ID=v5,",
        "##FILTER=<ID=Q20,",
        "##FILTER=<ID=AMPBIAS,",
        "##INFO=<ID=GDAMP,",
        "##INFO=<ID=TLAMP,",
    ):
        if req not in hdr_blob:
            raise AssertionError(f"{label}: missing required amplicon preset header {req}")
    if not recs:
        raise AssertionError(f"{label}: expected >= 1 emitted VCF record")
    for r in recs:
        info = r.split("\t")[7]
        if "GDAMP=" not in info or ";TLAMP=" not in info:
            raise AssertionError(f"{label}: record missing GDAMP;TLAMP in INFO: {r}")
    return recs

# 3a) Single-sample 8-column amplicon ($G4_DIR/amplicon/amplicon.bam + correct_8columns_amplicon.bed, NO -a, NO -z)
amp_single_recs = audit_amplicon_preset(
    "amp_single_l861q",
    fa_hg19,
    os.path.join(g4_dir, "amplicon", "amplicon.bam"),
    os.path.join(g4_dir, "amplicon", "correct_8columns_amplicon.bed"),
    [],
)

# 3b) Paired tumor-normal Illumina PE two-primer amplicon (-b "T|N", NO -a, NO -z)
amp_paired_recs = audit_amplicon_preset(
    "amp_paired_ilm",
    grch38_ref,
    f"{w4_dir}/bam/SampleA_ILM1.bam|{w4_dir}/bam/SampleB_ILM1.bam",
    os.path.join(w4_dir, "beds", "ilm_twoprimer_grch38.amp8.bed"),
    ["-N", "SampleA_ILM1|SampleB_ILM1", "--fisher"],
)
if len(amp_paired_recs) < 10:
    raise AssertionError(f"Expected >= 10 paired Illumina PE amplicon records, got {len(amp_paired_recs)}")

# 3c) QIAseq V3 single-primer extension (--preset amplicon -sp, NO -a, NO -z)
qgn_recs = audit_amplicon_preset(
    "amp_qiaseq_sp",
    grch38_ref,
    os.path.join(w4_dir, "bam", "SampleA_QGN1.bam"),
    os.path.join(w4_dir, "beds", "qiaseq_spe_grch38.amp8.bed"),
    ["-sp", "-UN"],
)
kp_loci = set()
with gzip.open(os.path.join(w4_dir, "raw", "KnownPositives_hg19ToHg38.vcf.gz"), "rt", encoding="utf-8") as fh:
    for line in fh:
        if not line.startswith("#"):
            p = line.rstrip("\n").split("\t")
            kp_loci.add((p[0], int(p[1]), p[4]))
qgn_pass_kp = sum(
    1 for r in qgn_recs
    if r.split("\t")[6] == "PASS"
    and (r.split("\t")[0], int(r.split("\t")[1]), r.split("\t")[4]) in kp_loci
)
if qgn_pass_kp < 15:
    raise AssertionError(f"QIAseq --preset amplicon -sp expected >= 15 PASS KnownPositives, got {qgn_pass_kp}")

print(
    f"[PASS] Gate 3/5 (--preset amplicon): L861Q={len(amp_single_recs)} calls, "
    f"Illumina_PE_T/N={len(amp_paired_recs)} calls, QIAseq_SPE(-sp)={len(qgn_recs)} calls (PASS KP={qgn_pass_kp})"
)
PY

# =============================================================================
# Gate 4: --preset ctdna on ~1,000x Real-Read 0.5% VAF Titration + d100 Override
#         (AR-06: SampleB_QGN1 WT + SampleA_QGN1 KnownPositive chr12:120999311 G>A)
# =============================================================================
echo "[4/5] Verifying --preset ctdna on ~1,000x real-read 0.5% VAF titration & d100 override precedence..."
python3 - "$VARDICTRUST_BIN" "$GRCH38_REF" "$WAVE4_DIR" "$WORK_TMP" <<'PY'
import gzip, os, re, subprocess, sys

vdr, grch38_ref, w4_dir, tmp = sys.argv[1:5]
bam_a = os.path.join(w4_dir, "bam", "SampleA_QGN1.bam")
bam_b = os.path.join(w4_dir, "bam", "SampleB_QGN1.bam")
bam_ac5 = os.path.join(w4_dir, "bam", "SampleAC5_QGN1.bam")

# Target KnownPositive SNV: chr12:120999311 G>A (inside QGN_KP_77 insert chr12:120999118-120999522)
target_chr = "chr12"
target_pos = 120999311
win_reg = "chr12:120999050-120999600"

# Replicate wild-type SampleB_QGN1.bam alignments 11x with unique read names in a single pass
wt_sam = os.path.join(tmp, "ctdna_titration.sam")
wt_records = subprocess.check_output(
    ["samtools", "view", "-F", "0xF04", bam_b, win_reg], text=True
).splitlines()
with open(wt_sam, "w", encoding="utf-8") as hdr_fh:
    subprocess.run(["samtools", "view", "-H", bam_b, win_reg], stdout=hdr_fh, check=True)
    hdr_fh.flush()
    for rep_idx in range(1, 12):
        prefix = f"wt_r{rep_idx}_"
        for rec in wt_records:
            if rec:
                hdr_fh.write(prefix + rec + "\n")

# Deterministically select 6 real variant-supporting (G>A) reads (3 forward + 3 reverse with distinct
# alignment positions so PSTD=1, PMEAN>=12, BIAS=2:2, QUAL>=30, MQ>=40) from SampleA_QGN1.bam!
cig_re = re.compile(r"(\d+)([MIDNSHP=X])")
sam_a_lines = subprocess.check_output(["samtools", "view", "-F", "0xF04", bam_a, f"{target_chr}:{target_pos}-{target_pos}"], text=True).splitlines()
fwd_picks, rev_picks = [], []
seen_fwd_pos, seen_rev_pos = set(), set()

for idx, raw_line in enumerate(sam_a_lines):
    cols = raw_line.split("\t")
    flag, rpos, mapq, seq, qual = int(cols[1]), int(cols[3]), int(cols[4]), cols[9], cols[10]
    if mapq < 40:
        continue
    qpos = 0
    hit_alt = False
    for n_s, op in cig_re.findall(cols[5]):
        n = int(n_s)
        if op in ("M", "=", "X"):
            if rpos <= target_pos < rpos + n:
                qi = qpos + (target_pos - rpos)
                dist_end = min(qi, len(seq) - 1 - qi)
                if seq[qi].upper() == "A" and (ord(qual[qi]) - 33) >= 30 and dist_end >= 12:
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
    aln_start = int(cols[3])
    # Clear paired flags on the 6 spiked single-molecule reads so each has an independent template name
    cols[0] = f"spikeA_{idx}_{cols[0]}"
    cols[1] = "16" if (flag & 16) else "0"
    cols[6] = "*"
    cols[7] = "0"
    cols[8] = "0"
    if (flag & 16) == 0 and aln_start not in seen_fwd_pos and len(fwd_picks) < 3:
        seen_fwd_pos.add(aln_start)
        fwd_picks.append("\t".join(cols))
    elif (flag & 16) != 0 and aln_start not in seen_rev_pos and len(rev_picks) < 3:
        seen_rev_pos.add(aln_start)
        rev_picks.append("\t".join(cols))

if len(fwd_picks) != 3 or len(rev_picks) != 3:
    raise AssertionError(f"Expected 3 fwd + 3 rev real SampleA A-reads at {target_chr}:{target_pos}, got {len(fwd_picks)}+{len(rev_picks)}")

with open(wt_sam, "a", encoding="utf-8") as out_fh:
    for rec in fwd_picks + rev_picks:
        out_fh.write(rec + "\n")

titration_bam = os.path.join(tmp, "ctdna_1000x_0p5pct.bam")
subprocess.run(["samtools", "sort", "-o", titration_bam, wt_sam], check=True)
subprocess.run(["samtools", "index", titration_bam], check=True)

ctdna_bed3 = os.path.join(tmp, "ctdna_locus_3col.bed")
with open(ctdna_bed3, "w", encoding="utf-8") as fh:
    fh.write("chr12\t120999118\t120999522\n")

# 4a) Default --compat=off (default -f 0.01) MUST emit 0 calls at the 0.5% VAF locus
off_default_vcf = subprocess.check_output(
    [vdr, "--compat=off", "-G", grch38_ref, "-b", titration_bam, "-c", "1", "-S", "2", "-E", "3", "-z", ctdna_bed3],
    text=True,
)
off_hits = [
    l for l in off_default_vcf.splitlines()
    if not l.startswith("#") and l.strip() and int(l.split("\t")[1]) == target_pos
]
if off_hits:
    raise AssertionError(f"Expected default --compat=off (-f 0.01) to emit 0 calls at 0.5% locus, got: {off_hits}")

# 4b) --preset ctdna (NO -c/-S/-E/-z and NO -f!) MUST call the 0.5% allele as FILTER=PASS with GT=0/1
ctdna_gz = os.path.join(tmp, "ctdna_1000x.vcf.gz")
subprocess.run(
    [vdr, "--preset", "ctdna", "-G", grch38_ref, "-b", titration_bam, "-N", "CTDNA_0P5",
     "--output", ctdna_gz, "--tbi", ctdna_bed3],
    check=True,
)
subprocess.run(
    ["bcftools", "norm", "--check-ref", "e", "-f", grch38_ref, ctdna_gz],
    stdout=subprocess.DEVNULL,
    stderr=subprocess.DEVNULL,
    check=True,
)
with gzip.open(ctdna_gz, "rt", encoding="utf-8") as fh:
    ct_lines = [l.rstrip("\n") for l in fh]
ct_hdrs = [l for l in ct_lines if l.startswith("#")]
ct_recs = [l for l in ct_lines if l and not l.startswith("#")]
if ct_hdrs[2] != "##VarDictRustCompat=off" or ct_hdrs[3] != "##VarDictRustPreset=ctdna":
    raise AssertionError(f"Unexpected --preset ctdna header lines: {ct_hdrs[:5]}")
ct_hdr_blob = "\n".join(ct_hdrs)
for req_flt in ("f0.002", "d100", "v4", "q25", "Q20"):
    if f"##FILTER=<ID={req_flt}," not in ct_hdr_blob:
        raise AssertionError(f"--preset ctdna missing FILTER header {req_flt}")

locus_hits = [r.split("\t") for r in ct_recs if int(r.split("\t")[1]) == target_pos and r.split("\t")[4] == "A"]
if len(locus_hits) != 1:
    raise AssertionError(f"Expected exact chr12:{target_pos} G>A call under --preset ctdna, got: {ct_recs}")
hit = locus_hits[0]
info_kv = dict(x.split("=", 1) for x in hit[7].split(";") if "=" in x)
fmt_kv = dict(zip(hit[8].split(":"), hit[9].split(":")))
dp_val, vd_val, af_val = int(info_kv["DP"]), int(info_kv["VD"]), float(info_kv["AF"])
if dp_val < 1000 or vd_val != 6 or not (0.003 <= af_val <= 0.008):
    raise AssertionError(f"Unexpected titration metrics at chr12:{target_pos}: DP={dp_val}, VD={vd_val}, AF={af_val}")
if hit[6] != "PASS" or fmt_kv.get("GT") != "0/1":
    raise AssertionError(f"Expected FILTER=PASS and GT=0/1 on 0.5% ctDNA call, got FILTER={hit[6]!r}, FORMAT={fmt_kv}")

# 4c) Raw SampleAC5_QGN1.bam (max depth < 100 across QIAseq inserts): verify d100 soft-filtering
#     (accounting for var2vcf_valid.pl:183 hicnt * hiaf >= 0.5 germline rescue) and explicit
#     --vcf-min-depth 20 override promotion!
ac5_bed3 = os.path.join(tmp, "ac5_inserts_3col.bed")
with open(os.path.join(w4_dir, "beds", "qiaseq_spe_grch38.amp8.bed"), "r", encoding="utf-8") as src, \
     open(ac5_bed3, "w", encoding="utf-8") as dst:
    for line in src:
        if line.strip() and not line.startswith("#"):
            p = line.rstrip("\n").split("\t")
            dst.write(f"{p[0]}\t{p[6]}\t{p[7]}\n")

ac5_default = subprocess.check_output(
    [vdr, "--preset", "ctdna", "-G", grch38_ref, "-b", bam_ac5, "-th", "4", ac5_bed3],
    text=True,
)
ac5_rows = [l.split("\t") for l in ac5_default.splitlines() if l and not l.startswith("#")]
ac5_pass = [r for r in ac5_rows if r[6] == "PASS"]
ac5_d100 = [r for r in ac5_rows if "d100" in r[6].split(";")]
d100_only_sites = {(r[0], r[1], r[3], r[4]) for r in ac5_rows if r[6] == "d100"}
low_af_pass = []
for r in ac5_pass:
    kv = dict(x.split("=", 1) for x in r[7].split(";") if "=" in x)
    if int(kv["HICNT"]) * float(kv["HIAF"]) < 0.5:
        low_af_pass.append((r[0], r[1], r[3], r[4], kv.get("DP"), kv.get("AF")))

if len(ac5_d100) < 20 or len(d100_only_sites) < 10 or low_af_pass:
    raise AssertionError(
        f"Expected raw SampleAC5_QGN1 under --preset ctdna to have >=20 d100 rows, >=10 d100-only sites, "
        f"and 0 low-AF (hicnt*hiaf<0.5) PASS rows; got d100={len(ac5_d100)}, "
        f"d100_only={len(d100_only_sites)}, low_af_pass={low_af_pass[:5]}"
    )

ac5_override = subprocess.check_output(
    [vdr, "--preset", "ctdna", "--vcf-min-depth", "20", "-G", grch38_ref, "-b", bam_ac5, "-th", "4", ac5_bed3],
    text=True,
)
if "##FILTER=<ID=d20," not in ac5_override or "##FILTER=<ID=d100," in ac5_override:
    raise AssertionError("Expected --vcf-min-depth 20 to replace d100 header with d20")
ac5_ov_rows = [
    l.split("\t") for l in ac5_override.splitlines()
    if l and not l.startswith("#")
]
if any("d100" in r[6].split(";") for r in ac5_ov_rows):
    raise AssertionError("Unexpected residual d100 FILTER tag under --vcf-min-depth 20")
ac5_ov_pass = [r for r in ac5_ov_rows if r[6] == "PASS"]
ac5_ov_pass_sites = {(r[0], r[1], r[3], r[4]) for r in ac5_ov_pass}
missing_promoted = d100_only_sites - ac5_ov_pass_sites
if missing_promoted or (len(ac5_ov_pass) - len(ac5_pass)) < 10:
    raise AssertionError(
        f"Expected all {len(d100_only_sites)} d100-only sites (>=10) to promote to PASS under --vcf-min-depth 20; "
        f"missing={sorted(missing_promoted)[:5]}, pass_delta={len(ac5_ov_pass) - len(ac5_pass)}"
    )

print(
    f"[PASS] Gate 4/5 (--preset ctdna): chr12:{target_pos} G>A DP={dp_val} VD={vd_val} AF={af_val:.4f} "
    f"(--compat=off=0 calls vs ctdna=PASS GT=0/1); SampleAC5 d100={len(ac5_d100)} "
    f"(d100_only={len(d100_only_sites)}) -> d20 PASS={len(ac5_ov_pass)} (+{len(ac5_ov_pass) - len(ac5_pass)})"
)
PY

# =============================================================================
# Gate 5: -R 1-Based Coordinate Invariance & Negative CLI Guards
#         (AR-03, AR-10, AR-11, AR-12, Issue #96)
# =============================================================================
echo "[5/5] Verifying -R 1-based coordinate invariance (-f override) & 7 negative CLI guards..."
python3 - "$VARDICTRUST_BIN" "$GIAB_REF" "$FA_HG19" "$G4_DIR" "$WORK_TMP" <<'PY'
import os, subprocess, sys

vdr, giab_ref, fa_hg19, g4_dir, tmp = sys.argv[1:6]
hg002_bam = os.path.join(g4_dir, "hg002", "HG002.chr20.bam")
amp_bam = os.path.join(g4_dir, "amplicon", "amplicon.bam")
amp_bed = os.path.join(g4_dir, "amplicon", "correct_8columns_amplicon.bed")

# 5a) --preset wgs -f 0.01 -R vs --compat=off -f 0.01 -R:
#     Must produce identical 1-based (CHROM, POS, REF, ALT) coordinates (0 1-bp shift, AR-03)
#     and emit ##FILTER=<ID=f0.01,...> instead of f0.05 (-f cleared preset vcf_af=0.05, AR-10)!
reg_str = "chr20:10000001-10035000"
preset_r_vcf = subprocess.check_output(
    [vdr, "--preset", "wgs", "-f", "0.01", "-G", giab_ref, "-b", hg002_bam, "-R", reg_str],
    text=True,
)
compat_r_vcf = subprocess.check_output(
    [vdr, "--compat=off", "-f", "0.01", "-G", giab_ref, "-b", hg002_bam, "-R", reg_str],
    text=True,
)

if '##FILTER=<ID=f0.01,Description="Allele frequency < 0.01">' not in preset_r_vcf or "##FILTER=<ID=f0.05," in preset_r_vcf:
    raise AssertionError("AR-10 regression: explicit -f 0.01 did not override --preset wgs FILTER=f0.05")

preset_coords = [
    tuple(l.split("\t")[:5])
    for l in preset_r_vcf.splitlines() if l and not l.startswith("#")
]
compat_coords = [
    tuple(l.split("\t")[:5])
    for l in compat_r_vcf.splitlines() if l and not l.startswith("#")
]
if len(preset_coords) < 15 or preset_coords != compat_coords:
    raise AssertionError(
        f"AR-03 regression: --preset wgs -R coordinates diverged from --compat=off -R "
        f"({len(preset_coords)} vs {len(compat_coords)})"
    )

# 5b) All 7 negative guards must exit non-zero cleanly with actionable error messages (0 panics)
negative_cases = [
    (
        "--preset wgs --parity",
        [vdr, "--preset", "wgs", "--parity", "-G", giab_ref, "-b", hg002_bam, "-R", "chr20:10000001-10000100"],
        ("--preset", "--parity"),
    ),
    (
        "--preset wgs --compat=vardictjava",
        [vdr, "--preset", "wgs", "--compat=vardictjava", "-G", giab_ref, "-b", hg002_bam, "-R", "chr20:10000001-10000100"],
        ("--preset", "--compat=vardictjava"),
    ),
    (
        "--preset amplicon -R chr1:100-200",
        [vdr, "--preset", "amplicon", "-G", fa_hg19, "-b", amp_bam, "-R", "chr1:100-200"],
        ("--preset amplicon", "-R"),
    ),
    (
        "--preset amplicon -i",
        [vdr, "--preset", "amplicon", "-i", "-G", fa_hg19, "-b", amp_bam, amp_bed],
        ("--preset amplicon", "-i"),
    ),
    (
        "--preset hifi",
        [vdr, "--preset", "hifi", "-G", giab_ref, "-b", hg002_bam, "-R", "chr20:10000001-10000100"],
        ("not implemented", "#96"),
    ),
    (
        "--preset unknown",
        [vdr, "--preset", "unknown", "-G", giab_ref, "-b", hg002_bam, "-R", "chr20:10000001-10000100"],
        ("Invalid --preset", "wgs"),
    ),
    (
        "-preset wgs",
        [vdr, "-preset", "wgs", "-G", giab_ref, "-b", hg002_bam, "-R", "chr20:10000001-10000100"],
        ("--preset",),
    ),
]

for label, cmd, expected_substrings in negative_cases:
    res = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if res.returncode == 0:
        raise AssertionError(f"Negative guard '{label}' unexpectedly exited 0")
    combined = (res.stderr or "") + "\n" + (res.stdout or "")
    if "panicked at" in combined or "RUST_BACKTRACE" in combined:
        raise AssertionError(f"Negative guard '{label}' panicked:\n{combined}")
    for sub in expected_substrings:
        if sub not in combined:
            raise AssertionError(f"Negative guard '{label}' missing expected substring {sub!r} in:\n{combined}")

print(
    f"[PASS] Gate 5/5 (-R Coordinate Invariance & Negative Guards): "
    f"{len(preset_coords)}/{len(compat_coords)} 1-based -R sites identical (FILTER=f0.01), 7/7 negative guards rejected cleanly"
)
PY

ELAPSED_S=$(( ( $(now_ms) - T0_MS ) / 1000 ))
echo "================================================================================"
echo " ALL 5 CLINICAL PRESET PUBLIC-DATA GATES PASSED (${ELAPSED_S}s)"
echo "================================================================================"
