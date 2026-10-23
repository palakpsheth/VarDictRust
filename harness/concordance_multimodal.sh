#!/usr/bin/env bash
# =============================================================================
# harness/concordance_multimodal.sh
# Automated Tier 4 Expanded Biological Ground Truth Concordance Suite
# Evaluates dual-engine parity and RTG vcfeval concordance across expanded
# consortium benchmarks: GIAB HG005 (GRCh38) and GIAB NA24385 (GRCh37/b37).
#
# Adheres to Host Policy A3:
# - Minimal network bandwidth: streams small intervals to persistent truth_cache/
# - Maximum 8 threads (-th 8, -T 8)
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
CACHE_DIR="$ROOT_DIR/scratch/g4/truth_cache"
REF_DIR="$ROOT_DIR/scratch/g4/ref"
OUT_BASE="$ROOT_DIR/scratch/g4/concordance_multimodal"
THREADS=8

if ! command -v rtg &>/dev/null || ! command -v bcftools &>/dev/null; then
  export PATH="/home/palak/miniconda3/envs/g4-tools/bin:$PATH"
fi

mkdir -p "$CACHE_DIR" "$OUT_BASE"

# Ensure GRCh37 SDF exists
FA_B37="$ROOT_DIR/VarDictJava/testdata/fastas/GRCh37.fa"
SDF_B37="$REF_DIR/grch37.sdf"
if [[ ! -d "$SDF_B37" ]]; then
  echo "Formatting GRCh37 SDF..."
  rtg format -o "$SDF_B37" "$FA_B37"
fi

FA_GIAB="$ROOT_DIR/scratch/g4/ref/giabv3/GRCh38_GIABv3.fasta"
SDF_GIAB="$REF_DIR/giabv3/GRCh38_GIABv3.sdf"

TMPDIR=$(mktemp -d /tmp/vardict_multimodal_concordance_XXXXXX)
trap 'rm -rf "$TMPDIR"' EXIT

echo "================================================================================"
echo " VarDictRust v1.0.0 Expanded Tier 4 Truth Concordance Suite"
echo " Evaluating HG005 (GRCh38) and NA24385 (GRCh37/b37)"
echo " Cache: $CACHE_DIR | Temp: $TMPDIR"
echo "================================================================================"

# Helper: Ensure truth slice exists in cache (isolated execution in CACHE_DIR to avoid root .tbi litter)
ensure_truth_slice() {
  local id="$1"
  local url="$2"
  local interval="$3"
  local vcf_out="$CACHE_DIR/${id}_truth.vcf.gz"

  if [[ -f "$vcf_out" && -f "$vcf_out.tbi" ]] && [[ $(zgrep -v '^#' "$vcf_out" | head -n 1 | wc -l) -gt 0 ]]; then
    echo "  [cached truth] $id ($interval)"
    return 0
  fi

  echo -n "  [fetching truth] $id ($interval)... "
  (
    cd "$CACHE_DIR"
    tabix -H "$url" > "$TMPDIR/${id}_truth.vcf"
    tabix "$url" "$interval" >> "$TMPDIR/${id}_truth.vcf"
    rm -f ./*.tbi 2>/dev/null || true
  )
  bgzip -f "$TMPDIR/${id}_truth.vcf"
  mv "$TMPDIR/${id}_truth.vcf.gz" "$vcf_out"
  tabix -f -p vcf "$vcf_out"
  local count=$(zgrep -v '^#' "$vcf_out" | wc -l)
  echo "done ($count truth variants)"
}

# -----------------------------------------------------------------------------
# 1. NIST GIAB HG005 (Han Chinese Diversity WGS, GRCh38)
# -----------------------------------------------------------------------------
echo "[1/2] Benchmarking GIAB HG005 Han Chinese WGS vs NIST v4.2.1..."
URL_HG005="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/ChineseTrio/HG005_NA24631_son/NISTv4.2.1/GRCh38/HG005_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
ensure_truth_slice "hg005" "$URL_HG005" "chr3:40000000-40010000"

printf "chr3\t40000000\t40010000\n" > "$CACHE_DIR/hg005_eval.bed"
BAM_HG005="$ROOT_DIR/scratch/g4/multimodal/hg005_s4a.bam"

echo "  Calling variants (Oracle & Rust)..."
bash "$ROOT_DIR/harness/golden_store.sh" oracle conc-hg005 "$FA_GIAB" "$BAM_HG005" "$CACHE_DIR/hg005_eval.bed" "$TMPDIR/hg005_java.vcf" "HG005" -z -th "$THREADS" -f 0.05 2>/dev/null
"$ROOT_DIR/target/release/vardictrust" -G "$FA_GIAB" -N "HG005" -b "$BAM_HG005" -c 1 -S 2 -E 3 -g 4 -z -th "$THREADS" -f 0.05 "$CACHE_DIR/hg005_eval.bed" > "$TMPDIR/hg005_rust.vcf" 2>/dev/null

# NOTE: no conda-pin probe here by design (presence via run_cpp.sh exit-3 gate; version keys via staple_cpp_extra.sh/concordance_germline.sh); keep pin format in sync with run_extended_coverage.sh + run_multimodal_slices.sh (CPP_PIN_STR).
# Third arm: vardictcpp (advisory only — same flags as the gate arms; a missing
# pin or engine failure degrades to a loud SKIP and the gate below never sees it).
CPP_HG005="$TMPDIR/hg005_cpp.vcf"
if [ "${CPP_SKIP:-0}" = "1" ] || [ "${CPP_SKIP:-0}" = "true" ]; then
  echo "  CPP: SKIP: disabled via CPP_SKIP=1"; CPP_HG005=""
elif bash "$ROOT_DIR/harness/golden_store.sh" cpp conc-hg005 "$FA_GIAB" "$BAM_HG005" "$CACHE_DIR/hg005_eval.bed" "$CPP_HG005" "HG005" -z -th "$THREADS" -f 0.05 2>"$TMPDIR/hg005_cpp.err"; then
  :
else
  cpp_rc=$?
  if [ "$cpp_rc" -eq 3 ]; then reason=$(head -1 "$TMPDIR/hg005_cpp.err"); echo "  CPP: SKIP: ${reason#SKIP: }";
  else
    first_err=$(grep -vE '^(ERROR conda\.cli|cpp_chunk=|[[:space:]]*$)' "$TMPDIR/hg005_cpp.err" 2>/dev/null | head -n 1 | tr -d '\r' | cut -c 1-140 || true)
    echo "  CPP: SKIP: cpp engine failure (exit $cpp_rc${first_err:+: $first_err})"
  fi
  CPP_HG005=""
fi

echo "  Verifying live oracle parity..."
python3 "$ROOT_DIR/harness/compare_vcf.py" --strict-bytes "$TMPDIR/hg005_java.vcf" "$TMPDIR/hg005_rust.vcf" --eps 1e-6 >/dev/null

echo "  Running rtg vcfeval concordance..."
bash "$ROOT_DIR/harness/concordance_germline.sh" "hg005-truth" \
  "$TMPDIR/hg005_java.vcf" "$TMPDIR/hg005_rust.vcf" \
  "$CACHE_DIR/hg005_truth.vcf.gz" "$CACHE_DIR/hg005_eval.bed" \
  "$SDF_GIAB" "$OUT_BASE/hg005" "HG005" "$CPP_HG005"

# -----------------------------------------------------------------------------
# 2. NIST GIAB NA24385 b37 (Ashkenazi Son Deep WGS, GRCh37/b37)
# -----------------------------------------------------------------------------
echo "[2/2] Benchmarking GIAB NA24385 b37 WGS vs NIST v4.2.1..."
URL_NA24385="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh37/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"
ensure_truth_slice "na24385_b37" "$URL_NA24385" "20:36018823-36039513"

printf "20\t36018823\t36039513\n" > "$CACHE_DIR/na24385_eval.bed"
BAM_NA24385="$ROOT_DIR/VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/NA24385_chr_20_36018823_36039513.bam"

echo "  Calling variants (Oracle & Rust, -N HG002)..."
bash "$ROOT_DIR/harness/golden_store.sh" oracle conc-na24385 "$FA_B37" "$BAM_NA24385" "$CACHE_DIR/na24385_eval.bed" "$TMPDIR/na24385_java.vcf" "HG002" -z -th "$THREADS" -f 0.001 2>/dev/null
"$ROOT_DIR/target/release/vardictrust" -G "$FA_B37" -N "HG002" -b "$BAM_NA24385" -c 1 -S 2 -E 3 -g 4 -z -th "$THREADS" -f 0.001 "$CACHE_DIR/na24385_eval.bed" > "$TMPDIR/na24385_rust.vcf" 2>/dev/null

# Third arm: vardictcpp (advisory only — same flags as the gate arms).
CPP_NA24385="$TMPDIR/na24385_cpp.vcf"
if [ "${CPP_SKIP:-0}" = "1" ] || [ "${CPP_SKIP:-0}" = "true" ]; then
  echo "  CPP: SKIP: disabled via CPP_SKIP=1"; CPP_NA24385=""
elif bash "$ROOT_DIR/harness/golden_store.sh" cpp conc-na24385 "$FA_B37" "$BAM_NA24385" "$CACHE_DIR/na24385_eval.bed" "$CPP_NA24385" "HG002" -z -th "$THREADS" -f 0.001 2>"$TMPDIR/na24385_cpp.err"; then
  :
else
  cpp_rc=$?
  if [ "$cpp_rc" -eq 3 ]; then reason=$(head -1 "$TMPDIR/na24385_cpp.err"); echo "  CPP: SKIP: ${reason#SKIP: }";
  else
    first_err=$(grep -vE '^(ERROR conda\.cli|cpp_chunk=|[[:space:]]*$)' "$TMPDIR/na24385_cpp.err" 2>/dev/null | head -n 1 | tr -d '\r' | cut -c 1-140 || true)
    echo "  CPP: SKIP: cpp engine failure (exit $cpp_rc${first_err:+: $first_err})"
  fi
  CPP_NA24385=""
fi

echo "  Verifying live oracle parity..."
python3 "$ROOT_DIR/harness/compare_vcf.py" --strict-bytes "$TMPDIR/na24385_java.vcf" "$TMPDIR/na24385_rust.vcf" --eps 1e-6 >/dev/null

echo "  Running rtg vcfeval concordance..."
bash "$ROOT_DIR/harness/concordance_germline.sh" "na24385-b37-truth" \
  "$TMPDIR/na24385_java.vcf" "$TMPDIR/na24385_rust.vcf" \
  "$CACHE_DIR/na24385_b37_truth.vcf.gz" "$CACHE_DIR/na24385_eval.bed" \
  "$SDF_B37" "$OUT_BASE/na24385" "HG002" "$CPP_NA24385"

echo "================================================================================"
echo " All Expanded Multi-Modal Truth Benchmarks Verified: EQUAL (0 site diffs)"
echo " Evidence recorded in benchmarks/evidence/"
echo "================================================================================"
