#!/usr/bin/env bash
# =============================================================================
# harness/get_wave4_public_slices.sh
# Idempotent, SHA-256-pinned ingestion and alignment of 6 open-access clinical
# cancer, somatic SV, and liquid-biopsy benchmarking cohorts for Wave-4
# real-world stress testing (zero binary fixtures staged in git).
#
# Clinical Cohorts & Upstream Provenance Registry:
# 1. FDA SEQC2 Oncopanel Sequencing Working Group (Gong et al., Sci. Data 2022,
#    DOI: 10.1038/s41597-022-01359-2; ENA SRP295113; Figshare 19128005.v2 &
#    13511829.v4):
#    - Cohort A: QIAseq V3 Single-Primer Extension (SPE) solid tumor Sample A
#      (SRR13224707), matched normal Sample B (SRR13224695), and 0.5-1% ctDNA
#      liquid-biopsy dilution Sample AC5 (SRR13224675). Exercises
#      AstraZeneca-NGS/VarDictJava#301, AstraZeneca-NGS/VarDictJava#303,
#      AstraZeneca-NGS/VarDictJava#326, AstraZeneca-NGS/VarDict#191 (-sp), and
#      AstraZeneca-NGS/VarDictJava#378 / palakpsheth/VarDictRust#95.
#    - Cohort B: Illumina paired-end two-primer oncology panel Sample A
#      (SRR13289946) & Sample B (SRR13289937), staged as both BAM and CRAM.
#      Exercises AstraZeneca-NGS/VarDictJava#155, AstraZeneca-NGS/VarDict#109,
#      MHH-Bioinformatics-Hematology/vardictcpp#1,
#      MHH-Bioinformatics-Hematology/vardictcpp#3, and
#      palakpsheth/VarDictRust#95.
#    - Cohort C: Ion Torrent S5 XL single-end homopolymer amplicon panel
#      Sample A (SRR13169349) & Sample B (SRR13169319). Exercises
#      AstraZeneca-NGS/VarDictJava#402 + reconcile_vd_dp + -sp + paired somatic
#      amplicon calling.
# 2. Clinical BRCA2 Patient Adapter-Chimera Tumor BAM
#    (PolinaBevad/VarDictJava@92c62cf6c3520b877c394a4c9fd5fe0d371cd972;
#    AstraZeneca-NGS/VarDictJava#240, AstraZeneca-NGS/VarDictJava#241,
#    AstraZeneca-NGS/VarDict#119) + Clinical MUC19 Somatic Paralog Pair
#    (AstraZeneca-NGS/VarDictJava#376).
# 3. clintval/vartovcf@3b7bd742926c62988ec2fd8420b5f237f43db6ff SV & contig-end
#    reference suite (clintval/vartovcf#107, clintval/vartovcf#108,
#    clintval/vartovcf#139, clintval/vartovcf#149, clintval/vartovcf#192,
#    AstraZeneca-NGS/VarDictJava#398, AstraZeneca-NGS/VarDictJava#225,
#    AstraZeneca-NGS/VarDictJava#229, AstraZeneca-NGS/VarDict#111,
#    AstraZeneca-NGS/VarDictJava#382).
# 4. NIST GIAB HG002 Clinical Medical Relevant Genes (CMRG v1.00) &
#    MosaicSNVv1.1 truth + 2x250 bp GRCh38 slice (AstraZeneca-NGS/VarDict@009e017,
#    MHH-Bioinformatics-Hematology/vardictcpp#4, palakpsheth/VarDictRust#94,
#    AstraZeneca-NGS/VarDictJava#65, AstraZeneca-NGS/VarDict#155), alongside
#    on-disk GIAB HG008-T|N (Fang et al. / NIST v0.4 & v0.6) and SEQC2
#    HCC1395|HCC1395BL (DOI: 10.1038/s41587-021-00993-6; joachimwolff/VarDictJava#4).
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

for _req_tool in bwa samtools bcftools tabix perl python3 curl gzip; do
    if ! command -v "$_req_tool" >/dev/null 2>&1; then
        echo "ERROR: required tool '$_req_tool' not found on PATH ($PATH)" >&2
        exit 1
    fi
done
unset _req_tool

DATA_ROOT="$(resolve_data_root "$ROOT_DIR")"
G4_DIR="${DATA_DIR:-$DATA_ROOT/scratch/g4}"
OUT_DIR="${WAVE4_DIR:-$G4_DIR/wave4}"
RAW_DIR="$OUT_DIR/raw"
ENA_DIR="$RAW_DIR/ena"
FASTQ_DIR="$OUT_DIR/fastq"
BAM_DIR="$OUT_DIR/bam"
BEDS_DIR="$OUT_DIR/beds"
BRCA2_DIR="$OUT_DIR/brca2_adaptor"
MUC19_DIR="$OUT_DIR/muc19"
VARTOVCF_DIR="$OUT_DIR/vartovcf"
BEDZIP_DIR="$OUT_DIR/bed_zip"
HG002_DIR="$OUT_DIR/hg002_cmrg_mosaic"

GRCH38_REF="$G4_DIR/ref/d1vd1/GRCh38.d1.vd1.fa"
GIAB_REF="$G4_DIR/ref/giabv3/GRCh38_GIABv3.fa"
if [ ! -f "$GIAB_REF" ] && [ -f "$G4_DIR/ref/giabv3/GRCh38_GIABv3.fasta" ]; then
    GIAB_REF="$G4_DIR/ref/giabv3/GRCh38_GIABv3.fasta"
fi

WAVE4_SCHEMA_VERSION="v1.1-wave4-6cohorts-r1"
DONE_STAMP="$OUT_DIR/.done"
MANIFEST_FILE="$OUT_DIR/.wave4_manifest.sha256"
T0_MS="$(now_ms)"

# -----------------------------------------------------------------------------
# Fast-path idempotent cache check (<2s when already staged & SHA-256 verified)
# -----------------------------------------------------------------------------
if [ "${FORCE_RESTAGE:-0}" != "1" ] && [ -f "$DONE_STAMP" ] && [ -s "$MANIFEST_FILE" ] \
   && [ "$(cat "$DONE_STAMP" 2>/dev/null || true)" = "$WAVE4_SCHEMA_VERSION" ]; then
    if python3 - "$OUT_DIR" "$MANIFEST_FILE" "$GRCH38_REF" <<'PY'
import hashlib, os, subprocess, sys

out_dir, manifest_path, grch38_ref = sys.argv[1], sys.argv[2], sys.argv[3]
if not os.path.isfile(grch38_ref) or not os.path.isfile(grch38_ref + ".fai"):
    sys.exit(1)

with open(manifest_path, "r", encoding="utf-8") as fh:
    for line in fh:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        want_sha, rel_path = line.split("\t", 1)
        full_path = os.path.join(out_dir, rel_path)
        if not os.path.isfile(full_path):
            sys.exit(1)
        h = hashlib.sha256()
        with open(full_path, "rb") as rf:
            for chunk in iter(lambda: rf.read(1 << 20), b""):
                h.update(chunk)
        if h.hexdigest() != want_sha:
            sys.exit(1)

required_nonempty = [
    "brca2_adaptor/adaptor_insertion.bam",
    "brca2_adaptor/adaptor_insertion.bam.bai",
    "brca2_adaptor/brca2_adaptor.bed",
    "muc19/MUC19_vaf1.bam",
    "muc19/MUC19_vaf1.bam.bai",
    "muc19/MUC19_vaf1_dz.bam",
    "muc19/MUC19_vaf1_dz.bam.bai",
    "muc19/MUC19_vaf1.bed",
    "muc19/hg19_chr12_muc19.fa",
    "muc19/hg19_chr12_muc19.fa.fai",
    "vartovcf/calls.var",
    "vartovcf/calls.vcf",
    "vartovcf/calls.contig-ends.var",
    "vartovcf/calls.contig-ends.vcf",
    "vartovcf/reference.fa",
    "vartovcf/reference.fa.fai",
    "vartovcf/vartovcf_sv_contig1.bam",
    "vartovcf/vartovcf_sv_contig1.bam.bai",
    "vartovcf/vartovcf_sv_contig1.bed",
    "bam/SampleA_QGN1.bam",
    "bam/SampleA_QGN1.bam.bai",
    "bam/SampleB_QGN1.bam",
    "bam/SampleB_QGN1.bam.bai",
    "bam/SampleAC5_QGN1.bam",
    "bam/SampleAC5_QGN1.bam.bai",
    "bam/SampleA_ILM1.bam",
    "bam/SampleA_ILM1.bam.bai",
    "bam/SampleA_ILM1.cram",
    "bam/SampleA_ILM1.cram.crai",
    "bam/SampleB_ILM1.bam",
    "bam/SampleB_ILM1.bam.bai",
    "bam/SampleB_ILM1.cram",
    "bam/SampleB_ILM1.cram.crai",
    "bam/SampleA_TFS1.bam",
    "bam/SampleA_TFS1.bam.bai",
    "bam/SampleB_TFS1.bam",
    "bam/SampleB_TFS1.bam.bai",
    "hg002_cmrg_mosaic/HG002.cmrg_mosaic.slice.bam",
    "hg002_cmrg_mosaic/HG002.cmrg_mosaic.slice.bam.bai",
    "hg002_cmrg_mosaic/hg002_cmrg_mosaic.bed",
    "beds/qiaseq_spe_grch38.amp8.bed",
    "beds/ilm_twoprimer_grch38.amp8.bed",
    "beds/ilm_primer_overlap_grch38.amp8.bed",
    "beds/tfs_iontorrent_grch38.amp8.bed",
]
for rel in required_nonempty:
    p = os.path.join(out_dir, rel)
    if not os.path.isfile(p) or os.path.getsize(p) == 0:
        sys.exit(1)

bams_to_check = [
    os.path.join(out_dir, r) for r in required_nonempty if r.endswith((".bam", ".cram"))
]
rc = subprocess.call(["samtools", "quickcheck", *bams_to_check])
sys.exit(0 if rc == 0 else 1)
PY
    then
        T1_MS="$(now_ms)"
        ELAPSED_MS=$(( T1_MS - T0_MS ))
        echo "[cached] Wave-4 6-cohort clinical datasets verified at $OUT_DIR (${ELAPSED_MS} ms)"
        exit 0
    fi
fi

mkdir -p "$RAW_DIR" "$ENA_DIR" "$FASTQ_DIR" "$BAM_DIR" "$BEDS_DIR" \
         "$BRCA2_DIR" "$MUC19_DIR" "$VARTOVCF_DIR" "$BEDZIP_DIR" "$HG002_DIR"

echo "================================================================================"
echo " Staging Wave-4 6-Cohort Clinical Datasets"
echo " Destination: $OUT_DIR"
echo "================================================================================"

# -----------------------------------------------------------------------------
# Reference Verification (GRCh38.d1.vd1.fa + BWA index + chr13 LN:114364328)
# -----------------------------------------------------------------------------
echo "[0/7] Verifying on-disk GRCh38.d1.vd1.fa reference & BWA index..."
for _ext in "" ".fai" ".amb" ".ann" ".bwt" ".pac" ".sa"; do
    if [ ! -s "${GRCH38_REF}${_ext}" ]; then
        echo "ERROR: missing required reference file ${GRCH38_REF}${_ext}" >&2
        exit 1
    fi
done
unset _ext

CHR13_LEN="$(awk '$1 == "chr13" {print $2; exit}' "${GRCH38_REF}.fai")"
if [ "$CHR13_LEN" != "114364328" ]; then
    echo "ERROR: expected chr13 length 114364328 in ${GRCH38_REF}.fai, got '${CHR13_LEN}'" >&2
    exit 1
fi

# -----------------------------------------------------------------------------
# Asset 1: Clinical BRCA2 Patient Adapter-Chimera BAM (AstraZeneca-NGS/VarDictJava#240)
# -----------------------------------------------------------------------------
echo "[1/7] Staging BRCA2 patient adapter-chimera BAM (PolinaBevad/VarDictJava@92c62cf6)..."
ADAPTOR_BASE_URL="https://raw.githubusercontent.com/PolinaBevad/VarDictJava/92c62cf6c3520b877c394a4c9fd5fe0d371cd972/src/test/resources/com/astrazeneca/vardict/integrationtests"
ADAPTOR_BAM_SHA="86f4dcaddacc40bc18ce4acd5a854a6c65e71b5460ab037f52d4c0e6d5165f8d"
ADAPTOR_BAI_SHA="855ed867215c799568b01548819383c6f1e66a3a203a14883e9824d63e604afc"

stage_verified_sha256 "$ADAPTOR_BASE_URL/adaptor_insertion.bam" "$ADAPTOR_BAM_SHA" "$RAW_DIR/adaptor_insertion.bam"
stage_verified_sha256 "$ADAPTOR_BASE_URL/adaptor_insertion.bam.bai" "$ADAPTOR_BAI_SHA" "$RAW_DIR/adaptor_insertion.bam.bai"

cp -f "$RAW_DIR/adaptor_insertion.bam" "$BRCA2_DIR/adaptor_insertion.bam"
cp -f "$RAW_DIR/adaptor_insertion.bam.bai" "$BRCA2_DIR/adaptor_insertion.bam.bai"
BAM_CHR13_LEN="$(samtools view -H "$BRCA2_DIR/adaptor_insertion.bam" | awk -F'\t' '$1=="@SQ" && !found {sn=""; ln=""; for(i=2;i<=NF;i++){if($i~/^SN:/) sn=substr($i,4); if($i~/^LN:/) ln=substr($i,4)} if(sn=="chr13") {print ln; found=1}}')"
if [ "$BAM_CHR13_LEN" != "114364328" ]; then
    echo "ERROR: adaptor_insertion.bam chr13 LN ($BAM_CHR13_LEN) != GRCh38.d1.vd1.fa chr13 LN (114364328)" >&2
    exit 1
fi
printf 'chr13\t32683876\t32684023\tBRCA2_ADAPTOR\n' > "$BRCA2_DIR/brca2_adaptor.bed"

# -----------------------------------------------------------------------------
# Asset 2: Clinical MUC19 Somatic Paralog Fixture (AstraZeneca-NGS/VarDictJava#376)
# -----------------------------------------------------------------------------
echo "[2/7] Staging MUC19 somatic paralog fixture (AstraZeneca-NGS/VarDictJava#376)..."
MUC19_ZIP_URL="https://github.com/AstraZeneca-NGS/VarDictJava/files/10389155/MUC19_vaf1.zip"
MUC19_ZIP_SHA="1874bae4622ea30c0e7065e39899ef5ec850e17688762be4896e302740c2edc4"
MUC19_REF_WIN_URL="https://rest.ensembl.org/sequence/region/human/12:40870001..40880000:1?coord_system_version=GRCh37;content-type=text/plain"
MUC19_REF_WIN_SHA="13048d713c664c8a10b39d2e60ca9854fe89822ae1a0bf768db6d4e59adc8857"

stage_verified_sha256 "$MUC19_ZIP_URL" "$MUC19_ZIP_SHA" "$RAW_DIR/MUC19_vaf1.zip"
stage_verified_sha256 "$MUC19_REF_WIN_URL" "$MUC19_REF_WIN_SHA" "$RAW_DIR/grch37_chr12_40870001_40880000.txt"

python3 - "$RAW_DIR/MUC19_vaf1.zip" "$RAW_DIR/grch37_chr12_40870001_40880000.txt" "$MUC19_DIR" <<'PY'
import os, sys, zipfile

zip_path, win_path, out_dir = sys.argv[1], sys.argv[2], sys.argv[3]
with zipfile.ZipFile(zip_path, "r") as zf:
    for member in zf.infolist():
        if member.is_dir():
            continue
        base = os.path.basename(member.filename)
        if not base or base.startswith("."):
            continue
        dest = os.path.join(out_dir, base)
        if not os.path.isfile(dest) or os.path.getsize(dest) != member.file_size:
            with zf.open(member, "r") as src, open(dest, "wb") as dst:
                dst.write(src.read())

fa_path = os.path.join(out_dir, "hg19_chr12_muc19.fa")
total_len = 133851895
if not os.path.isfile(fa_path) or os.path.getsize(fa_path) < total_len:
    with open(win_path, "r", encoding="utf-8") as fh:
        win_seq = "".join(fh.read().split()).upper()
    if len(win_seq) != 10000:
        raise RuntimeError(f"Expected 10000 bp GRCh37 chr12 window, got {len(win_seq)}")
    tmp_fa = fa_path + ".tmp"
    # 1-based 40870001..40880000 == 0-based [40870000, 40880000), exact multiples of 1000!
    line_n = ("N" * 1000 + "\n").encode("ascii")
    with open(tmp_fa, "wb") as out:
        out.write(b">chr12\n")
        for _ in range(40870):
            out.write(line_n)
        for off in range(0, 10000, 1000):
            out.write((win_seq[off : off + 1000] + "\n").encode("ascii"))
        rem = total_len - 40880000
        full_lines, tail = divmod(rem, 1000)
        for _ in range(full_lines):
            out.write(line_n)
        if tail:
            out.write(("N" * tail + "\n").encode("ascii"))
    os.replace(tmp_fa, fa_path)
PY
if [ ! -s "$MUC19_DIR/hg19_chr12_muc19.fa.fai" ]; then
    samtools faidx "$MUC19_DIR/hg19_chr12_muc19.fa"
fi

# -----------------------------------------------------------------------------
# Asset 3: clintval/vartovcf@3b7bd742 SV & Contig-End Reference Suite
# -----------------------------------------------------------------------------
echo "[3/7] Staging clintval/vartovcf@3b7bd742 SV & contig-end fixtures..."
V2V_BASE_URL="https://raw.githubusercontent.com/clintval/vartovcf/3b7bd742926c62988ec2fd8420b5f237f43db6ff/tests"
V2V_CALLS_VAR_SHA="97559956a0df08a1cc9af1841e51692b3e8e443b78f0d30610c5b41d04e12113"
V2V_CALLS_VCF_SHA="1e12f5f7ab5f46bd1014cff3c103954d364a2009b5869dc0e79d3e0d037f0888"
V2V_REF_FA_SHA="9454726c6559c3994864b9d3860bd4bf22f384d8dcc952134172cad352153d09"
V2V_CENDS_VAR_SHA="09a31b8002212ac5af3cffa820abb47c09a479a92570e6b3dbfca0b95bb8076a"
V2V_CENDS_VCF_SHA="c848f3a9d12a710a80abd3ef20463b963f551f2bab9fd0d05ff661f5f35614be"

stage_verified_sha256 "$V2V_BASE_URL/calls.var" "$V2V_CALLS_VAR_SHA" "$RAW_DIR/vartovcf_calls.var"
stage_verified_sha256 "$V2V_BASE_URL/calls.vcf" "$V2V_CALLS_VCF_SHA" "$RAW_DIR/vartovcf_calls.vcf"
stage_verified_sha256 "$V2V_BASE_URL/reference.fa" "$V2V_REF_FA_SHA" "$RAW_DIR/vartovcf_reference.fa"
stage_verified_sha256 "$V2V_BASE_URL/calls.contig-ends.var" "$V2V_CENDS_VAR_SHA" "$RAW_DIR/vartovcf_calls.contig-ends.var"
stage_verified_sha256 "$V2V_BASE_URL/calls.contig-ends.vcf" "$V2V_CENDS_VCF_SHA" "$RAW_DIR/vartovcf_calls.contig-ends.vcf"

cp -f "$RAW_DIR/vartovcf_calls.var" "$VARTOVCF_DIR/calls.var"
cp -f "$RAW_DIR/vartovcf_calls.vcf" "$VARTOVCF_DIR/calls.vcf"
cp -f "$RAW_DIR/vartovcf_reference.fa" "$VARTOVCF_DIR/reference.fa"
cp -f "$RAW_DIR/vartovcf_calls.contig-ends.var" "$VARTOVCF_DIR/calls.contig-ends.var"
cp -f "$RAW_DIR/vartovcf_calls.contig-ends.vcf" "$VARTOVCF_DIR/calls.contig-ends.vcf"
samtools faidx "$VARTOVCF_DIR/reference.fa"

# Build deterministic BAM reproducing clintval/vartovcf@3b7bd742 SV breakpoints
# (chr12:446 C <DUP> END=1446 SVLEN=1000, chr13:1141 G <DEL> END=4141 SVLEN=-3000,
# and chr1:1 contig-start pileup locus) directly on vartovcf/reference.fa.
if [ ! -s "$VARTOVCF_DIR/vartovcf_sv_contig1.bam" ] || [ ! -s "$VARTOVCF_DIR/vartovcf_sv_contig1.bam.bai" ]; then
    python3 - "$VARTOVCF_DIR/reference.fa" "$VARTOVCF_DIR/vartovcf_sv_contig1.sam" "$VARTOVCF_DIR/vartovcf_sv_contig1.bed" <<'PY'
import sys

ref_fa, out_sam, out_bed = sys.argv[1], sys.argv[2], sys.argv[3]
seqs = {}
cur = None
with open(ref_fa, "r", encoding="utf-8") as fh:
    for line in fh:
        line = line.strip()
        if line.startswith(">"):
            cur = line[1:].split()[0]
            seqs[cur] = []
        elif cur:
            seqs[cur].append(line.upper())
seqs = {k: "".join(v) for k, v in seqs.items()}

sam = [
    "@HD\tVN:1.6\tSO:coordinate",
    f"@SQ\tSN:chr1\tLN:{len(seqs['chr1'])}",
    f"@SQ\tSN:chr12\tLN:{len(seqs['chr12'])}",
    f"@SQ\tSN:chr13\tLN:{len(seqs['chr13'])}",
]
rid = 0
def add_read(chrom, pos1, cigar, seq, flag=0):
    global rid
    rid += 1
    sam.append(f"v2v_{rid}\t{flag}\t{chrom}\t{pos1}\t60\t{cigar}\t*\t0\t0\t{seq}\t{'I' * len(seq)}\tNM:i:0")

# Contig-start chr1:1 coverage + chr1:1..1000 split-read locus
for _ in range(12):
    add_read("chr1", 1, "100M", seqs["chr1"][0:100])
for i in range(6):
    off = i % 3
    m1, sc = 60 - off, 40 + off
    p1 = 1001 - m1
    add_read("chr1", p1, f"{m1}M{sc}S", seqs["chr1"][p1 - 1 : 1000] + seqs["chr1"][0:sc])
    add_read("chr1", 1, f"{sc}S{m1}M", seqs["chr1"][1000 - sc : 1000] + seqs["chr1"][0:m1])

# chr12:447..1446 <DUP> -> anchored chr12:446 C <DUP> END=1446 SVLEN=1000
for i in range(16):
    add_read("chr12", 400 + i * 3, "100M", seqs["chr12"][399 + i * 3 : 499 + i * 3])
    add_read("chr12", 1390 + i * 3, "100M", seqs["chr12"][1389 + i * 3 : 1489 + i * 3])
for i in range(6):
    off = i % 3
    m1, sc = 60 - off, 40 + off
    p1 = 1447 - m1
    add_read("chr12", p1, f"{m1}M{sc}S", seqs["chr12"][p1 - 1 : 1446] + seqs["chr12"][446 : 446 + sc])
    add_read("chr12", 447, f"{sc}S{m1}M", seqs["chr12"][1446 - sc : 1446] + seqs["chr12"][446 : 446 + m1])

# chr13:1141..4141 <DEL> -> anchored chr13:1141 G <DEL> END=4141 SVLEN=-3000
for i in range(16):
    add_read("chr13", 1090 + i * 3, "100M", seqs["chr13"][1089 + i * 3 : 1189 + i * 3])
    add_read("chr13", 4100 + i * 3, "100M", seqs["chr13"][4099 + i * 3 : 4199 + i * 3])
for i in range(6):
    off = i % 3
    m1, sc = 60 - off, 40 + off
    p1 = 1142 - m1
    add_read("chr13", p1, f"{m1}M{sc}S", seqs["chr13"][p1 - 1 : 1141] + seqs["chr13"][4141 : 4141 + sc])
    add_read("chr13", 4142, f"{sc}S{m1}M", seqs["chr13"][1141 - sc : 1141] + seqs["chr13"][4141 : 4141 + m1])

with open(out_sam, "w", encoding="utf-8") as fh:
    fh.write("\n".join(sam) + "\n")
with open(out_bed, "w", encoding="utf-8") as fh:
    fh.write("chr1\t0\t1200\tBRINP3\nchr12\t167\t1600\tPTPN11\nchr13\t1000\t4300\tFLT3\n")
PY
    samtools view -bS "$VARTOVCF_DIR/vartovcf_sv_contig1.sam" \
        | samtools sort -o "$VARTOVCF_DIR/vartovcf_sv_contig1.bam"
    rm -f "$VARTOVCF_DIR/vartovcf_sv_contig1.sam"
    samtools index "$VARTOVCF_DIR/vartovcf_sv_contig1.bam"
    samtools quickcheck "$VARTOVCF_DIR/vartovcf_sv_contig1.bam"
fi

# -----------------------------------------------------------------------------
# Asset 4: Official FDA SEQC2 Oncopanel GRCh38 Truth & Panel BEDs
#          (Gong et al., Sci. Data 2022, DOI: 10.1038/s41597-022-01359-2)
# -----------------------------------------------------------------------------
echo "[4/7] Staging FDA SEQC2 Oncopanel GRCh38 truth & panel BEDs (Figshare 19128005.v2 & 13511829.v4)..."
FIGSHARE_BED_ZIP_SHA="ff41bd6deb6cd2fc88e2f8f34c22a265d4f91cc3ae8fa31185fe27818e4dc082"
FIGSHARE_KP_VCF_SHA="66201ebd74452b7c0bd60750e49ff371c272c619736804f2d927682192a01bd6"
FIGSHARE_CTR_BED_SHA="bc1fb53cf7952312c464b44907daeda9fd84b91243ec08cc3cfd50f2f2374c31"
FIGSHARE_KN_BED_SHA="8891eb16d40756c254a3ae960ffab23577a6c906f8c3e134b0285e0240d4a102"

stage_verified_sha256 "https://ndownloader.figshare.com/files/33984191" "$FIGSHARE_BED_ZIP_SHA" "$RAW_DIR/BED.zip"
stage_verified_sha256 "https://ndownloader.figshare.com/files/41951184" "$FIGSHARE_KP_VCF_SHA" "$RAW_DIR/KnownPositives_hg19ToHg38.vcf.gz"
stage_verified_sha256 "https://ndownloader.figshare.com/files/25927461" "$FIGSHARE_CTR_BED_SHA" "$RAW_DIR/CTR_hg38.bed.gz"
stage_verified_sha256 "https://ndownloader.figshare.com/files/41951178" "$FIGSHARE_KN_BED_SHA" "$RAW_DIR/KnownNegatives_hg38.bed.gz"

python3 - "$RAW_DIR/BED.zip" "$BEDZIP_DIR" <<'PY'
import os, sys, zipfile

zip_path, out_dir = sys.argv[1], sys.argv[2]
with zipfile.ZipFile(zip_path, "r") as zf:
    for member in zf.infolist():
        if member.is_dir():
            continue
        rel = os.path.normpath(member.filename).lstrip("/\\")
        if ".." in rel.split(os.sep):
            raise RuntimeError(f"Unsafe zip member {member.filename}")
        dest = os.path.join(out_dir, rel)
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        if not os.path.isfile(dest) or os.path.getsize(dest) != member.file_size:
            with zf.open(member, "r") as src, open(dest, "wb") as dst:
                dst.write(src.read())
PY

# -----------------------------------------------------------------------------
# Assets 5-7: FDA SEQC2 ENA SRP295113 Real Clinical Sequencing Runs
#             (Bounded HTTPS Byte-Range Slices, SHA-256 Verified on Raw Bytes)
# -----------------------------------------------------------------------------
echo "[5/7] Staging ENA SRP295113 clinical FASTQ range slices (12 streams, SHA-256 verified)..."

ENA_SPECS=(
    # Cohort A: QIAseq V3 SPE Solid Tumor, Matched Normal, & 0.5% ctDNA Liquid Biopsy
    "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR132/007/SRR13224707/SRR13224707_1.fastq.gz|8388607|d99a6bf10e49bafe9b0b79d410ac21bd4a7ee4775d6bfbcedacf6d41bfd178d3|SRR13224707_1.range_0_8388607.fastq.gz"
    "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR132/007/SRR13224707/SRR13224707_2.fastq.gz|10485759|ba3b2937b99ef122de16ee6b452d36629086757a4dad6ec5966733600be72135|SRR13224707_2.range_0_10485759.fastq.gz"
    "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR132/095/SRR13224695/SRR13224695_1.fastq.gz|8388607|e31d932f49844af5a0da3984060d544fc22663d4331497ec425a53f8e71540fd|SRR13224695_1.range_0_8388607.fastq.gz"
    "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR132/095/SRR13224695/SRR13224695_2.fastq.gz|10485759|a7b2b561888246dfe4b00e7dffb105e9899ab3acfa340d43a4f288bd3c27590d|SRR13224695_2.range_0_10485759.fastq.gz"
    "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR132/075/SRR13224675/SRR13224675_1.fastq.gz|8388607|02e5187d566f8f9c7741540d40c589137c7e810820bdf2e2cbf34e30a2a0f698|SRR13224675_1.range_0_8388607.fastq.gz"
    "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR132/075/SRR13224675/SRR13224675_2.fastq.gz|10485759|3fa014156fd03b3c9f468f198355092c32afcb937fed0e6fd40cea1077887d7a|SRR13224675_2.range_0_10485759.fastq.gz"
    # Cohort B: Illumina Paired-End Oncology Panel (SampleA_ILM1 | SampleB_ILM1)
    "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR132/046/SRR13289946/SRR13289946_1.fastq.gz|8388607|09b60319dee644ff096e35ea1ad2e21f3ca37989398e97d13fc9defc7ba777e7|SRR13289946_1.range_0_8388607.fastq.gz"
    "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR132/046/SRR13289946/SRR13289946_2.fastq.gz|8388607|e414e4eba16554c504299630be827837a4ee6970768d5cc5a4a1e54cb7b4010b|SRR13289946_2.range_0_8388607.fastq.gz"
    "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR132/037/SRR13289937/SRR13289937_1.fastq.gz|8388607|96c3d024f1d389912d723b897d72487d13c31f49af3a14c28b118a6f0c013b51|SRR13289937_1.range_0_8388607.fastq.gz"
    "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR132/037/SRR13289937/SRR13289937_2.fastq.gz|8388607|50c9fb063065e893055a191de3a9d60f941a02c8316e5903c0e55a4a1324069f|SRR13289937_2.range_0_8388607.fastq.gz"
    # Cohort C: Ion Torrent S5 XL Single-End Homopolymer Amplicon (SampleA_TFS1 | SampleB_TFS1)
    "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR131/049/SRR13169349/SRR13169349.fastq.gz|8388607|0deed25a7fcd75e37415fd6efcba0c5f79c8567fbab852d69ab4d608dd2d5f74|SRR13169349.range_0_8388607.fastq.gz"
    "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR131/019/SRR13169319/SRR13169319.fastq.gz|8388607|5d90c9c9d9c0b621b2677308eb4e0b3381cdca998b99c1f4596c0c049045dd92|SRR13169319.range_0_8388607.fastq.gz"
)

fetch_ena_range_one() {
    local url="$1"
    local end_byte="$2"
    local want_sha="$3"
    local dest="$4"
    local want_len=$(( end_byte + 1 ))
    if [ -f "$dest" ]; then
        local cur_sha
        cur_sha="$(portable_sha256 "$dest")"
        if [ "$cur_sha" = "$want_sha" ]; then
            return 0
        fi
    fi
    local part="${dest}.part"
    rm -f "$part"
    echo "  [downloading] $(basename "$dest") (bytes 0-${end_byte})..."
    if ! curl -fsSL --retry 5 --retry-delay 5 -r "0-${end_byte}" -o "$part" "$url"; then
        rm -f "$part"
        echo "ERROR: ENA range download failed for $url" >&2
        return 1
    fi
    local got_len got_sha
    got_len="$(portable_filesize "$part")"
    got_sha="$(portable_sha256 "$part")"
    if [ "$got_len" != "$want_len" ] || [ "$got_sha" != "$want_sha" ]; then
        rm -f "$part"
        echo "ERROR: ENA range verification failed for $(basename "$dest") (len=$got_len want=$want_len, sha=$got_sha want=$want_sha)" >&2
        return 1
    fi
    mv -f "$part" "$dest"
}

ENA_PIDS=()
for _spec in "${ENA_SPECS[@]}"; do
    IFS='|' read -r _url _end _sha _fname <<< "$_spec"
    if [ -f "$ENA_DIR/$_fname" ] && [ "$(portable_sha256 "$ENA_DIR/$_fname")" = "$_sha" ]; then
        echo "  [verified] $_fname"
    else
        fetch_ena_range_one "$_url" "$_end" "$_sha" "$ENA_DIR/$_fname" &
        ENA_PIDS+=("$!")
    fi
done
unset _spec _url _end _sha _fname

if [ "${#ENA_PIDS[@]}" -gt 0 ]; then
    for _pid in "${ENA_PIDS[@]}"; do
        wait "$_pid"
    done
fi
unset ENA_PIDS _pid

# Decompress bounded gzip streams via Python zlib.decompressobj(16 + zlib.MAX_WBITS),
# verify read-ID sync, trim QIAseq V3 5' 23 bp UMI+ATTGGAGTCCT & 3' read-through
# adapters, and align to GRCh38.d1.vd1.fa with bwa mem.
align_clinical_sample() {
    local mode="$1"
    local sample="$2"
    local r1_gz="$3"
    local r2_gz="$4"
    local out_bam="$BAM_DIR/${sample}.bam"
    local fq1="$FASTQ_DIR/${sample}_R1.fastq"
    local fq2="$FASTQ_DIR/${sample}_R2.fastq"

    if [ -s "$out_bam" ] && [ -s "${out_bam}.bai" ] && samtools quickcheck "$out_bam" 2>/dev/null; then
        echo "  [cached] ${sample}.bam"
    else
        echo "  [aligning] ${sample} (${mode}) -> ${sample}.bam..."
        python3 - "$mode" "$sample" "$r1_gz" "$r2_gz" "$fq1" "$fq2" <<'PY'
import sys, zlib

mode, sample, r1_gz, r2_gz, fq1_out, fq2_out = sys.argv[1:7]

def read_partial_fastq_gz(path):
    dec = zlib.decompressobj(16 + zlib.MAX_WBITS)
    with open(path, "rb") as fh:
        raw = fh.read()
    text = dec.decompress(raw).decode("ascii", errors="replace")
    lines = text.splitlines()
    # Drop trailing partial line/record at byte-range boundary
    if not text.endswith("\n") and lines:
        lines = lines[:-1]
    n_recs = len(lines) // 4
    recs = []
    for i in range(n_recs):
        h, s, plus, q = lines[4 * i : 4 * i + 4]
        if not h.startswith("@") or not plus.startswith("+") or len(s) != len(q):
            break
        qname = h[1:].split()[0].split("/")[0]
        recs.append((qname, h, s, q))
    return recs

def trim_3p_motif(seq, qual, motif):
    idx = seq.find(motif)
    if idx != -1:
        return seq[:idx], qual[:idx]
    for k in range(len(motif) - 1, 7, -1):
        if seq.endswith(motif[:k]):
            return seq[:-k], qual[:-k]
    return seq, qual

if mode == "iontorrent_se":
    r1 = read_partial_fastq_gz(r1_gz)
    kept = 0
    with open(fq1_out, "w", encoding="ascii") as out1:
        for qname, _h, s, q in r1:
            if len(s) >= 35:
                out1.write(f"@{qname}\n{s}\n+\n{q}\n")
                kept += 1
    if kept < 10000:
        raise RuntimeError(f"{sample}: too few single-end reads ({kept})")
    print(f"    {sample}: {kept} SE reads extracted")
else:
    r1 = read_partial_fastq_gz(r1_gz)
    r2 = read_partial_fastq_gz(r2_gz)
    k = min(len(r1), len(r2))
    if k < 20000:
        raise RuntimeError(f"{sample}: too few synced read pairs ({k})")
    motif = "ATTGGAGTCCT"
    exact_hits = 0
    hd_le2_hits = 0
    kept = 0
    with open(fq1_out, "w", encoding="ascii") as out1, open(fq2_out, "w", encoding="ascii") as out2:
        for i in range(k):
            qn1, _h1, s1, q1 = r1[i]
            qn2, _h2, s2, q2 = r2[i]
            if qn1 != qn2:
                raise RuntimeError(f"{sample}: R1/R2 QNAME mismatch at record {i}: {qn1} != {qn2}")
            if mode == "qiaseq_spe":
                if len(s2) < 58:
                    continue
                seg = s2[12:23]
                hd = sum(1 for a, b in zip(seg, motif) if a != b)
                if hd == 0:
                    exact_hits += 1
                if hd <= 2:
                    hd_le2_hits += 1
                s2, q2 = s2[23:], q2[23:]
                s1, q1 = trim_3p_motif(s1, q1, "AGGACTCCAAT")
                s2, q2 = trim_3p_motif(s2, q2, "CAAAACGCAATACTGTACATT")
            if len(s1) >= 35 and len(s2) >= 35:
                out1.write(f"@{qn1}\n{s1}\n+\n{q1}\n")
                out2.write(f"@{qn2}\n{s2}\n+\n{q2}\n")
                kept += 1
    if mode == "qiaseq_spe" and (exact_hits < 500 or hd_le2_hits < 15000):
        raise RuntimeError(
            f"{sample}: QIAseq V3 5' R2 ATTGGAGTCCT motif check failed "
            f"(exact={exact_hits}, hd<=2={hd_le2_hits})"
        )
    print(f"    {sample}: {kept}/{k} PE pairs synced (motif exact={exact_hits}, hd<=2={hd_le2_hits})")
PY
        local pl="ILLUMINA"
        if [ "$mode" = "iontorrent_se" ]; then
            pl="IONTORRENT"
            bwa mem -v 1 -t 4 -R "@RG\\tID:${sample}\\tSM:${sample}\\tPL:${pl}" "$GRCH38_REF" "$fq1" \
                | samtools sort -@ 4 -o "${out_bam}.tmp"
        else
            bwa mem -v 1 -t 4 -R "@RG\\tID:${sample}\\tSM:${sample}\\tPL:${pl}" "$GRCH38_REF" "$fq1" "$fq2" \
                | samtools sort -@ 4 -o "${out_bam}.tmp"
        fi
        mv -f "${out_bam}.tmp" "$out_bam"
        samtools index -@ 4 "$out_bam"
        samtools quickcheck "$out_bam"
    fi

    if [ "$mode" = "illumina_pe" ]; then
        local out_cram="$BAM_DIR/${sample}.cram"
        if [ -s "$out_cram" ] && [ -s "${out_cram}.crai" ] && samtools quickcheck "$out_cram" 2>/dev/null; then
            echo "  [cached] ${sample}.cram"
        else
            echo "  [encoding] ${sample}.cram..."
            samtools view -@ 4 -C -T "$GRCH38_REF" -o "${out_cram}.tmp" "$out_bam"
            mv -f "${out_cram}.tmp" "$out_cram"
            samtools index -@ 4 "$out_cram"
            samtools quickcheck "$out_cram"
        fi
    fi
}

align_clinical_sample "qiaseq_spe" "SampleA_QGN1" \
    "$ENA_DIR/SRR13224707_1.range_0_8388607.fastq.gz" \
    "$ENA_DIR/SRR13224707_2.range_0_10485759.fastq.gz"
align_clinical_sample "qiaseq_spe" "SampleB_QGN1" \
    "$ENA_DIR/SRR13224695_1.range_0_8388607.fastq.gz" \
    "$ENA_DIR/SRR13224695_2.range_0_10485759.fastq.gz"
align_clinical_sample "qiaseq_spe" "SampleAC5_QGN1" \
    "$ENA_DIR/SRR13224675_1.range_0_8388607.fastq.gz" \
    "$ENA_DIR/SRR13224675_2.range_0_10485759.fastq.gz"

align_clinical_sample "illumina_pe" "SampleA_ILM1" \
    "$ENA_DIR/SRR13289946_1.range_0_8388607.fastq.gz" \
    "$ENA_DIR/SRR13289946_2.range_0_8388607.fastq.gz"
align_clinical_sample "illumina_pe" "SampleB_ILM1" \
    "$ENA_DIR/SRR13289937_1.range_0_8388607.fastq.gz" \
    "$ENA_DIR/SRR13289937_2.range_0_8388607.fastq.gz"

align_clinical_sample "iontorrent_se" "SampleA_TFS1" \
    "$ENA_DIR/SRR13169349.range_0_8388607.fastq.gz" ""
align_clinical_sample "iontorrent_se" "SampleB_TFS1" \
    "$ENA_DIR/SRR13169319.range_0_8388607.fastq.gz" ""

# -----------------------------------------------------------------------------
# Deterministic 8-Column Amplicon BED Generation Across All 3 Clinical Panels
# -----------------------------------------------------------------------------
echo "[6/7] Generating deterministic 8-column GRCh38 amplicon BEDs..."
python3 - "$RAW_DIR" "$BAM_DIR" "$BEDS_DIR" <<'PY'
import bisect, gzip, os, re, subprocess, sys
from collections import Counter

raw_dir, bam_dir, beds_dir = sys.argv[1:4]

kp = {}
with gzip.open(os.path.join(raw_dir, "KnownPositives_hg19ToHg38.vcf.gz"), "rt", encoding="utf-8") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        kp.setdefault(p[0], []).append((int(p[1]), p[3], p[4]))

ctr = {}
with gzip.open(os.path.join(raw_dir, "CTR_hg38.bed.gz"), "rt", encoding="utf-8") as fh:
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        p = line.rstrip("\n").split("\t")
        ctr.setdefault(p[0], []).append((int(p[1]), int(p[2])))
for c in ctr:
    ctr[c].sort()

kn = {}
with gzip.open(os.path.join(raw_dir, "KnownNegatives_hg38.bed.gz"), "rt", encoding="utf-8") as fh:
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        p = line.rstrip("\n").split("\t")
        kn.setdefault(p[0], []).append((int(p[1]), int(p[2])))
for c in kn:
    kn[c].sort()

def overlaps_sorted(bed_map, chrom, start0, end0):
    ivs = bed_map.get(chrom, [])
    idx = bisect.bisect_left(ivs, (start0, -1))
    for j in range(max(0, idx - 2), min(len(ivs), idx + 5)):
        if max(start0, ivs[j][0]) < min(end0, ivs[j][1]):
            return True
    return False

cigar_re = re.compile(r"(\d+)([MIDNSHP=X])")
def aligned_len(cigar):
    return sum(int(n) for n, op in cigar_re.findall(cigar) if op in "MDN=X")

def collect_spe_clusters(bam_path, r1_only=True):
    clusters = Counter()
    spans = {}
    cmd = ["samtools", "view", "-F", "2308"]
    if r1_only:
        cmd.extend(["-f", "64"])
    cmd.append(bam_path)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    assert proc.stdout is not None
    for line in proc.stdout:
        f = line.split("\t")
        fl, c, pos0, tlen = int(f[1]), f[2], int(f[3]) - 1, int(f[8])
        if not c.startswith("chr"):
            continue
        end0 = pos0 + aligned_len(f[5])
        if (fl & 16) == 0:
            key = (c, "+", pos0)
            pe = max(end0, pos0 + tlen) if 0 < tlen < 600 else end0
            clusters[key] += 1
            s0, e0 = spans.get(key, (pos0, pe))
            spans[key] = (min(s0, pos0), max(e0, pe))
        else:
            key = (c, "-", end0)
            ps = min(pos0, end0 + tlen) if -600 < tlen < 0 else pos0
            clusters[key] += 1
            s0, e0 = spans.get(key, (ps, end0))
            spans[key] = (min(s0, ps), max(e0, end0))
    if proc.wait() != 0:
        raise RuntimeError(f"samtools view failed on {bam_path}")
    return clusters, spans

def chrom_sort_key(c):
    body = c[3:] if c.startswith("chr") else c
    if body.isdigit():
        return (0, int(body), "")
    return (1, 99, body)

# 1) QIAseq V3 SPE 8-col BED:
#    + strand R1 GSP: region = [gsp_5p, max_r2_end], insert = [gsp_5p + 22, max_r2_end]
#    - strand R1 GSP: region = [min_r2_start, gsp_3p], insert = [min_r2_start, gsp_3p - 22]
ca, sa = collect_spe_clusters(os.path.join(bam_dir, "SampleA_QGN1.bam"), r1_only=True)
cb, sb = collect_spe_clusters(os.path.join(bam_dir, "SampleB_QGN1.bam"), r1_only=True)
cac5, sac5 = collect_spe_clusters(os.path.join(bam_dir, "SampleAC5_QGN1.bam"), r1_only=True)

selected_qgn = {}
def try_add_qgn(key, tag_prefix):
    if key in selected_qgn:
        return False
    c, strand, _ = key
    s_list = [ m[key][0] for m in (sa, sb, sac5) if key in m ]
    e_list = [ m[key][1] for m in (sa, sb, sac5) if key in m ]
    if not s_list:
        return False
    st, en = min(s_list), max(e_list)
    ins_s, ins_e = (st + 22, en) if strand == "+" else (st, en - 22)
    if ins_e - ins_s < 45 or en - st > 480:
        return False
    if not overlaps_sorted(ctr, c, ins_s, ins_e) and not any(ins_s < v[0] <= ins_e for v in kp.get(c, [])):
        return False
    selected_qgn[key] = (c, st, en, f"{tag_prefix}_{len(selected_qgn) + 1}", strand, ins_s, ins_e)
    return True

# Top SampleA KnownPositives SPE clusters
for key, cnt in sorted(ca.items(), key=lambda kv: (-kv[1], chrom_sort_key(kv[0][0]), kv[0][2], kv[0][1])):
    if cnt < 18 or len(selected_qgn) >= 130:
        break
    c, strand, _ = key
    st, en = sa[key]
    ins_s, ins_e = (st + 22, en) if strand == "+" else (st, en - 22)
    if any(ins_s < v[0] <= ins_e for v in kp.get(c, [])):
        try_add_qgn(key, "QGN_KP")

# SampleAC5 ctDNA KnownPositives SPE clusters
for key, cnt in sorted(cac5.items(), key=lambda kv: (-kv[1], chrom_sort_key(kv[0][0]), kv[0][2], kv[0][1])):
    if cnt < 16 or len(selected_qgn) >= 165:
        break
    c, strand, _ = key
    st, en = sac5[key]
    ins_s, ins_e = (st + 22, en) if strand == "+" else (st, en - 22)
    if any(ins_s < v[0] <= ins_e for v in kp.get(c, [])):
        try_add_qgn(key, "QGN_AC5")

# High-depth KnownNegatives-only control SPE clusters
for key, cnt in sorted(ca.items(), key=lambda kv: (-kv[1], chrom_sort_key(kv[0][0]), kv[0][2], kv[0][1])):
    if cnt < 35 or len(selected_qgn) >= 190:
        break
    if cb[key] < 25 or cac5[key] < 22:
        continue
    c, strand, _ = key
    st, en = sa[key]
    ins_s, ins_e = (st + 22, en) if strand == "+" else (st, en - 22)
    if any(st <= v[0] <= en for v in kp.get(c, [])):
        continue
    if overlaps_sorted(kn, c, ins_s, ins_e):
        try_add_qgn(key, "QGN_NEG")

qgn_rows = sorted(selected_qgn.values(), key=lambda r: (chrom_sort_key(r[0]), r[1], r[2], r[3]))
with open(os.path.join(beds_dir, "qiaseq_spe_grch38.amp8.bed"), "w", encoding="utf-8") as out:
    for c, st, en, name, strand, ins_s, ins_e in qgn_rows:
        out.write(f"{c}\t{st}\t{en}\t{name}\t0\t{strand}\t{ins_s}\t{ins_e}\n")

# 2) Illumina Two-Primer Overlapping 8-Col Amplicon BED + Primer-Overlap Boundary BED
def count_kp_depths(bam_path, kp_positions_by_chrom):
    counts = Counter()
    proc = subprocess.Popen(["samtools", "view", "-F", "2308", bam_path], stdout=subprocess.PIPE, text=True)
    assert proc.stdout is not None
    for line in proc.stdout:
        f = line.split("\t")
        c = f[2]
        pts = kp_positions_by_chrom.get(c)
        if not pts:
            continue
        pos0 = int(f[3]) - 1
        end0 = pos0 + aligned_len(f[5])
        lo = bisect.bisect_right(pts, pos0)
        hi = bisect.bisect_right(pts, end0)
        for j in range(lo, hi):
            counts[(c, pts[j])] += 1
    if proc.wait() != 0:
        raise RuntimeError(f"samtools view failed on {bam_path}")
    return counts

kp_pos_sorted = {c: sorted({v[0] for v in vs}) for c, vs in kp.items()}
ilm_dp_a = count_kp_depths(os.path.join(bam_dir, "SampleA_ILM1.bam"), kp_pos_sorted)
ilm_dp_b = count_kp_depths(os.path.join(bam_dir, "SampleB_ILM1.bam"), kp_pos_sorted)
ilm_tiles = [
    (ca_dp + ilm_dp_b[(c, pos1)], c, pos1)
    for (c, pos1), ca_dp in ilm_dp_a.items()
    if ca_dp >= 35 and ilm_dp_b[(c, pos1)] >= 25
]

ilm_tiles.sort(key=lambda t: (-t[0], chrom_sort_key(t[1]), t[2]))
chosen_ilm = []
seen_loci = set()
for _dp, c, pos1 in ilm_tiles:
    if any(ec == c and abs(ep - pos1) < 200 for ec, ep in seen_loci):
        continue
    seen_loci.add((c, pos1))
    # Overlapping two-primer amplicon pair with 22 bp synthetic 5'/3' primers
    s1, e1 = pos1 - 65, pos1 + 65
    s2, e2 = pos1 - 53, pos1 + 77
    chosen_ilm.append((c, s1, e1, f"ILM_{c}_{pos1}_T1", "+", s1 + 22, e1 - 22))
    chosen_ilm.append((c, s2, e2, f"ILM_{c}_{pos1}_T2", "+", s2 + 22, e2 - 22))
    if len(chosen_ilm) >= 80:
        break

chosen_ilm.sort(key=lambda r: (chrom_sort_key(r[0]), r[1], r[2], r[3]))
with open(os.path.join(beds_dir, "ilm_twoprimer_grch38.amp8.bed"), "w", encoding="utf-8") as out:
    for c, st, en, name, strand, ins_s, ins_e in chosen_ilm:
        out.write(f"{c}\t{st}\t{en}\t{name}\t0\t{strand}\t{ins_s}\t{ins_e}\n")

# Dedicated primer-overlap BED placing high-depth KnownPositive hotspot inside the 22 bp 5' primer
# (EGFR chr7:55146655 sits inside [55146645, 55146667) primer margin -> must be suppressed in T & N!)
primer_pos = 55146655
primer_chr = "chr7"
if ("chr7", 55146655) not in seen_loci and seen_loci:
    primer_chr, primer_pos = sorted(seen_loci, key=lambda x: (chrom_sort_key(x[0]), x[1]))[0]
p_st = primer_pos - 10
p_en = primer_pos + 120
with open(os.path.join(beds_dir, "ilm_primer_overlap_grch38.amp8.bed"), "w", encoding="utf-8") as out:
    out.write(f"{primer_chr}\t{p_st}\t{p_en}\tILM_PRIMER_MASK_{primer_chr}_{primer_pos}\t0\t+\t{p_st + 22}\t{p_en - 22}\n")

# 3) Ion Torrent S5 XL Single-End 8-Col Amplicon BED
ctfs_a, stfs_a = collect_spe_clusters(os.path.join(bam_dir, "SampleA_TFS1.bam"), r1_only=False)
ctfs_b, stfs_b = collect_spe_clusters(os.path.join(bam_dir, "SampleB_TFS1.bam"), r1_only=False)
tfs_rows = []
for key, cnt in sorted(ctfs_a.items(), key=lambda kv: (-kv[1], chrom_sort_key(kv[0][0]), kv[0][2], kv[0][1])):
    if cnt < 20 or len(tfs_rows) >= 80:
        break
    if ctfs_b[key] < 15:
        continue
    c, strand, pos5 = key
    st = min(stfs_a[key][0], stfs_b[key][0])
    en = max(stfs_a[key][1], stfs_b[key][1])
    ins_s, ins_e = (st + 22, en) if strand == "+" else (st, en - 22)
    if ins_e - ins_s < 40 or en - st > 350:
        continue
    tfs_rows.append((c, st, en, f"TFS_{strand}_{pos5}", strand, ins_s, ins_e))

tfs_rows.sort(key=lambda r: (chrom_sort_key(r[0]), r[1], r[2], r[3]))
with open(os.path.join(beds_dir, "tfs_iontorrent_grch38.amp8.bed"), "w", encoding="utf-8") as out:
    for c, st, en, name, strand, ins_s, ins_e in tfs_rows:
        out.write(f"{c}\t{st}\t{en}\t{name}\t0\t{strand}\t{ins_s}\t{ins_e}\n")

print(f"  Generated BEDs: QGN={len(qgn_rows)} rows, ILM={len(chosen_ilm)} rows (+1 primer-mask), TFS={len(tfs_rows)} rows")
PY

# -----------------------------------------------------------------------------
# Cohort 6: GIAB HG002 CMRG v1.00 & MosaicSNVv1.1 Truth + 2x250 bp GRCh38 Slice
# -----------------------------------------------------------------------------
echo "[7/7] Staging GIAB HG002 CMRG v1.00 + MosaicSNVv1.1 truth & BAM slice..."
CMRG_BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/SmallVariant"
MOSAIC_BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/mosaic_v1.1/GRCh38/SNV"

CMRG_VCF_SHA="80f18efa283d953cfc615d30dd01b1ff80debcc46f52042b942da783acd48da2"
CMRG_TBI_SHA="d3ddfd35884bf9d8dcde4fdba0ae10edbb539b87207df5d4e5db2f7aec925231"
CMRG_BED_SHA="dc77803bb474802fd3213d2f3649076dc14675d3b3c0ae3284036be09d2813f7"
MOSAIC_VCF_SHA="c30fabbb23e25f3eb0bf08c9b39edc3bb731db612b83f8874225d0c1c53a48d0"
MOSAIC_TBI_SHA="10c895bfe25cddec1818cfb9936ae93c2a2cc2146749a81aa9285e45a639b2e9"
MOSAIC_BED_SHA="d80eae760c5581f22205cb4845db99f1f5501662e02be8a03ad3a2aa3873502f"

stage_verified_sha256 "$CMRG_BASE/HG002_GRCh38_CMRG_smallvar_v1.00.vcf.gz" "$CMRG_VCF_SHA" "$HG002_DIR/HG002_GRCh38_CMRG_smallvar_v1.00.vcf.gz"
stage_verified_sha256 "$CMRG_BASE/HG002_GRCh38_CMRG_smallvar_v1.00.vcf.gz.tbi" "$CMRG_TBI_SHA" "$HG002_DIR/HG002_GRCh38_CMRG_smallvar_v1.00.vcf.gz.tbi"
stage_verified_sha256 "$CMRG_BASE/HG002_GRCh38_CMRG_smallvar_v1.00.bed" "$CMRG_BED_SHA" "$HG002_DIR/HG002_GRCh38_CMRG_smallvar_v1.00.bed"
stage_verified_sha256 "$MOSAIC_BASE/HG002_GRCh38_MosaicSNVv1.1.vcf.gz" "$MOSAIC_VCF_SHA" "$HG002_DIR/HG002_GRCh38_MosaicSNVv1.1.vcf.gz"
stage_verified_sha256 "$MOSAIC_BASE/HG002_GRCh38_MosaicSNVv1.1.vcf.gz.tbi" "$MOSAIC_TBI_SHA" "$HG002_DIR/HG002_GRCh38_MosaicSNVv1.1.vcf.gz.tbi"
stage_verified_sha256 "$MOSAIC_BASE/HG002_GRCh38_MosaicSNVv1.1.bed" "$MOSAIC_BED_SHA" "$HG002_DIR/HG002_GRCh38_MosaicSNVv1.1.bed"

python3 - "$HG002_DIR" <<'PY'
import gzip, os, sys

d = sys.argv[1]
intervals = []
with open(os.path.join(d, "HG002_GRCh38_CMRG_smallvar_v1.00.bed"), "r", encoding="utf-8") as fh:
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        p = line.rstrip("\n").split("\t")
        s, e = int(p[1]), int(p[2])
        if p[0] == "chr1" and 500 <= (e - s) <= 5000:
            intervals.append((p[0], s, e, f"CMRG_{len(intervals) + 1}"))
        if len(intervals) >= 8:
            break

mosaic_n = 0
with gzip.open(os.path.join(d, "HG002_GRCh38_MosaicSNVv1.1.vcf.gz"), "rt", encoding="utf-8") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if p[0] in ("chr1", "chr2"):
            pos = int(p[1])
            intervals.append((p[0], pos - 200, pos + 200, f"MOSAIC_{mosaic_n + 1}"))
            mosaic_n += 1
        if mosaic_n >= 8:
            break

with open(os.path.join(d, "hg002_cmrg_mosaic.bed"), "w", encoding="utf-8") as out:
    for c, s, e, name in intervals:
        out.write(f"{c}\t{s}\t{e}\t{name}\n")
PY

HG002_SLICE_BAM="$HG002_DIR/HG002.cmrg_mosaic.slice.bam"
LOCAL_HG002_FULL="$G4_DIR/hg002/HG002.GRCh38.2x250.bam"
if [ -s "$HG002_SLICE_BAM" ] && [ -s "${HG002_SLICE_BAM}.bai" ] && samtools quickcheck "$HG002_SLICE_BAM" 2>/dev/null; then
    echo "  [cached] HG002.cmrg_mosaic.slice.bam"
elif [ -s "$LOCAL_HG002_FULL" ] && [ -s "${LOCAL_HG002_FULL}.bai" ]; then
    echo "  [slicing] HG002.cmrg_mosaic.slice.bam from local $LOCAL_HG002_FULL..."
    samtools view -@ 4 -b -M -L "$HG002_DIR/hg002_cmrg_mosaic.bed" "$LOCAL_HG002_FULL" \
        | samtools sort -@ 4 -o "${HG002_SLICE_BAM}.tmp"
    mv -f "${HG002_SLICE_BAM}.tmp" "$HG002_SLICE_BAM"
    samtools index -@ 4 "$HG002_SLICE_BAM"
    samtools quickcheck "$HG002_SLICE_BAM"
else
    HG002_S3_BAM="https://giab.s3.amazonaws.com/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams/HG002.GRCh38.2x250.bam"
    HG002_S3_BAI="https://giab.s3.amazonaws.com/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams/HG002.GRCh38.2x250.bam.bai"
    HG002_S3_BAI_SHA="723f69de14b8adbbf1a60f23f9af007310011b51b13ee01d18a2a8daf7282f36"
    stream_verified_bam_slice --bed "$HG002_DIR/hg002_cmrg_mosaic.bed" --threads 4 \
        "HG002_cmrg_mosaic" "$HG002_S3_BAM" "$HG002_S3_BAI" "$HG002_S3_BAI_SHA" \
        "$HG002_DIR/.HG002.GRCh38.2x250.bam.bai" "" "$HG002_SLICE_BAM"
fi

# -----------------------------------------------------------------------------
# Write Pinned SHA-256 Manifest & Completion Stamp
# -----------------------------------------------------------------------------
cat > "${MANIFEST_FILE}.tmp" <<EOF
# Wave-4 Clinical Dataset Pinned Raw Asset SHA-256 Manifest (${WAVE4_SCHEMA_VERSION})
${ADAPTOR_BAM_SHA}	raw/adaptor_insertion.bam
${ADAPTOR_BAI_SHA}	raw/adaptor_insertion.bam.bai
${MUC19_ZIP_SHA}	raw/MUC19_vaf1.zip
${MUC19_REF_WIN_SHA}	raw/grch37_chr12_40870001_40880000.txt
${V2V_CALLS_VAR_SHA}	raw/vartovcf_calls.var
${V2V_CALLS_VCF_SHA}	raw/vartovcf_calls.vcf
${V2V_REF_FA_SHA}	raw/vartovcf_reference.fa
${V2V_CENDS_VAR_SHA}	raw/vartovcf_calls.contig-ends.var
${V2V_CENDS_VCF_SHA}	raw/vartovcf_calls.contig-ends.vcf
${FIGSHARE_BED_ZIP_SHA}	raw/BED.zip
${FIGSHARE_KP_VCF_SHA}	raw/KnownPositives_hg19ToHg38.vcf.gz
${FIGSHARE_CTR_BED_SHA}	raw/CTR_hg38.bed.gz
${FIGSHARE_KN_BED_SHA}	raw/KnownNegatives_hg38.bed.gz
d99a6bf10e49bafe9b0b79d410ac21bd4a7ee4775d6bfbcedacf6d41bfd178d3	raw/ena/SRR13224707_1.range_0_8388607.fastq.gz
ba3b2937b99ef122de16ee6b452d36629086757a4dad6ec5966733600be72135	raw/ena/SRR13224707_2.range_0_10485759.fastq.gz
e31d932f49844af5a0da3984060d544fc22663d4331497ec425a53f8e71540fd	raw/ena/SRR13224695_1.range_0_8388607.fastq.gz
a7b2b561888246dfe4b00e7dffb105e9899ab3acfa340d43a4f288bd3c27590d	raw/ena/SRR13224695_2.range_0_10485759.fastq.gz
02e5187d566f8f9c7741540d40c589137c7e810820bdf2e2cbf34e30a2a0f698	raw/ena/SRR13224675_1.range_0_8388607.fastq.gz
3fa014156fd03b3c9f468f198355092c32afcb937fed0e6fd40cea1077887d7a	raw/ena/SRR13224675_2.range_0_10485759.fastq.gz
09b60319dee644ff096e35ea1ad2e21f3ca37989398e97d13fc9defc7ba777e7	raw/ena/SRR13289946_1.range_0_8388607.fastq.gz
e414e4eba16554c504299630be827837a4ee6970768d5cc5a4a1e54cb7b4010b	raw/ena/SRR13289946_2.range_0_8388607.fastq.gz
96c3d024f1d389912d723b897d72487d13c31f49af3a14c28b118a6f0c013b51	raw/ena/SRR13289937_1.range_0_8388607.fastq.gz
50c9fb063065e893055a191de3a9d60f941a02c8316e5903c0e55a4a1324069f	raw/ena/SRR13289937_2.range_0_8388607.fastq.gz
0deed25a7fcd75e37415fd6efcba0c5f79c8567fbab852d69ab4d608dd2d5f74	raw/ena/SRR13169349.range_0_8388607.fastq.gz
5d90c9c9d9c0b621b2677308eb4e0b3381cdca998b99c1f4596c0c049045dd92	raw/ena/SRR13169319.range_0_8388607.fastq.gz
${CMRG_VCF_SHA}	hg002_cmrg_mosaic/HG002_GRCh38_CMRG_smallvar_v1.00.vcf.gz
${CMRG_TBI_SHA}	hg002_cmrg_mosaic/HG002_GRCh38_CMRG_smallvar_v1.00.vcf.gz.tbi
${CMRG_BED_SHA}	hg002_cmrg_mosaic/HG002_GRCh38_CMRG_smallvar_v1.00.bed
${MOSAIC_VCF_SHA}	hg002_cmrg_mosaic/HG002_GRCh38_MosaicSNVv1.1.vcf.gz
${MOSAIC_TBI_SHA}	hg002_cmrg_mosaic/HG002_GRCh38_MosaicSNVv1.1.vcf.gz.tbi
${MOSAIC_BED_SHA}	hg002_cmrg_mosaic/HG002_GRCh38_MosaicSNVv1.1.bed
EOF
mv -f "${MANIFEST_FILE}.tmp" "$MANIFEST_FILE"
printf '%s\n' "$WAVE4_SCHEMA_VERSION" > "$DONE_STAMP"

T1_MS="$(now_ms)"
ELAPSED_MS=$(( T1_MS - T0_MS ))
echo "================================================================================"
echo " Wave-4 6-Cohort Clinical Datasets Ready ($OUT_DIR, ${ELAPSED_MS} ms)"
echo "================================================================================"
