#!/usr/bin/env bash
# harness/recertify.sh — Full release parity recertification & evidence updater
# Governed by rewrites.bio Manifesto §0, §0b, §5 and AGENTS.md §15.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$ROOT_DIR"
COMMON_GIT_DIR="$(git -C "$ROOT_DIR" rev-parse --path-format=absolute --git-common-dir 2>/dev/null || true)"
if [[ -n "$COMMON_GIT_DIR" && -d "$(dirname "$COMMON_GIT_DIR")/VarDictJava" ]]; then
    DATA_ROOT="$(cd "$(dirname "$COMMON_GIT_DIR")" && pwd)"
else
    DATA_ROOT="$ROOT_DIR"
fi

if ! command -v rtg &>/dev/null || ! command -v bcftools &>/dev/null; then
    export PATH="/home/palak/miniconda3/envs/g4-tools/bin:$PATH"
fi

ts() { date '+[%H:%M:%S]'; }
RECERT_T0=$SECONDS
echo "================================================================================"
echo " VarDictRust Live Dual-Engine Release Parity Recertification"
echo "================================================================================"

echo "$(ts) [Recertify 1/5] Building release binary..."
cargo build --release -p vardictrust
export VARDICT_BIN="$ROOT_DIR/target/release/vardictrust"

echo "$(ts) [Recertify 2/5] Materializing reference FASTAs..."
python3 "$ROOT_DIR/harness/convert_fastas.py"

FP=$(bash "$ROOT_DIR/harness/engine_fingerprint.sh")
PERF_FP=$(bash "$ROOT_DIR/harness/engine_fingerprint.sh" --perf)
echo "$(ts) Current Engine Fingerprint: $FP (perf_fingerprint: $PERF_FP)"
export GOLDEN_RUN_ID="recertify-$(date +%s)"

needs_stamp() {
    local name="$1"
    python3 - "$ROOT_DIR/benchmarks/evidence/${name}.json" "$FP" <<'PY'
import json, sys
p, fp = sys.argv[1], sys.argv[2]
try:
    d = json.load(open(p))
except Exception:
    sys.exit(0)
if d.get("fingerprint") == fp and d.get("result") in ("MATCH", "EQUAL", "PASS"):
    sys.exit(1)
sys.exit(0)
PY
}

# Verify scale.json is valid under either FP or PERF_FP; if stale, run scale_bench.sh under media_stack minimal
if python3 - "$ROOT_DIR/benchmarks/evidence/scale.json" "$FP" "$PERF_FP" <<'PY'
import json, sys
p, fp, perf = sys.argv[1], sys.argv[2], sys.argv[3]
try:
    d = json.load(open(p))
except Exception:
    sys.exit(1)
x = d.get("extra") or {}
sys.exit(0 if (d.get("fingerprint") == fp or x.get("perf_fingerprint") == perf) and d.get("result") in ("MATCH", "EQUAL", "PASS") else 1)
PY
then
    echo "$(ts) [REUSE] scale.json: certified under perf_fingerprint $PERF_FP"
else
    echo "$(ts) [RUN] scale.json: perf_fingerprint changed ($PERF_FP); running 24-cell scale benchmark under media_stack minimal..."
    trap 'python3 "$ROOT_DIR/harness/media_stack_lifecycle.py" up || true' EXIT
    python3 "$ROOT_DIR/harness/media_stack_lifecycle.py" minimal
    EST_MIN=55 bash "$ROOT_DIR/harness/courtesy.sh" bash "$ROOT_DIR/harness/scale_bench.sh" \
      scratch/g4/ref/d1vd1/GRCh38.d1.vd1.fa \
      scratch/chr1/HG002.chr1.bam \
      scratch/chr1/chr1.bed \
      HG002 scratch/g4/scale
    python3 "$ROOT_DIR/harness/record_evidence.py" scale PASS scratch/g4/scale \
      -- env EST_MIN=55 bash harness/courtesy.sh bash harness/scale_bench.sh scratch/g4/ref/d1vd1/GRCh38.d1.vd1.fa scratch/chr1/HG002.chr1.bam scratch/chr1/chr1.bed HG002 scratch/g4/scale
    python3 "$ROOT_DIR/harness/media_stack_lifecycle.py" up
    trap - EXIT
fi

echo "$(ts) [Recertify 3/5] Running live dual-engine regression & cohort certification suites..."
if needs_stamp "amplicon-tier2"; then
    step_t0=$SECONDS
    echo "$(ts) [RUN] amplicon-tier2: Certifying Targeted Amplicon Tier-2..."
    KEEP_DIR=scratch/g4/amplicon/cmp MIN_VARIANTS=2 bash harness/compare.sh \
      scratch/g4/amplicon/ref.fa \
      scratch/g4/amplicon/amplicon.bam \
      scratch/g4/amplicon/correct_8columns_amplicon.bed \
      SAMPLE -a 10:0.95 -th 8
    python3 harness/record_evidence.py amplicon-tier2 MATCH scratch/g4/amplicon/cmp \
      -- env KEEP_DIR=scratch/g4/amplicon/cmp MIN_VARIANTS=2 bash harness/compare.sh scratch/g4/amplicon/ref.fa scratch/g4/amplicon/amplicon.bam scratch/g4/amplicon/correct_8columns_amplicon.bed SAMPLE -a 10:0.95 -th 8
    echo "$(ts) [DONE] amplicon-tier2 (${SECONDS-step_t0}s)"
else
    echo "$(ts) [REUSE] amplicon-tier2: already certified for $FP"
fi

if needs_stamp "hg008-tier2" || needs_stamp "hg008-truth"; then
    step_t0=$SECONDS
    echo "$(ts) [RUN] hg008-tier2 & hg008-truth: Certifying GIAB HG008 Somatic..."
    KEEP_DIR=scratch/g4/hg008/cmp MIN_VARIANTS=18426 bash harness/compare.sh \
      scratch/g4/ref/giabv3/GRCh38_GIABv3.fa \
      "scratch/g4/hg008/HG008-T.slice.bam|scratch/g4/hg008/HG008-N.slice.bam" \
      scratch/g4/hg008/engine.bed \
      "HG008-T|HG008-N" -th 4
    python3 harness/record_evidence.py hg008-tier2 MATCH scratch/g4/hg008/cmp \
      -- env KEEP_DIR=scratch/g4/hg008/cmp MIN_VARIANTS=18426 bash harness/compare.sh scratch/g4/ref/giabv3/GRCh38_GIABv3.fa "scratch/g4/hg008/HG008-T.slice.bam|scratch/g4/hg008/HG008-N.slice.bam" scratch/g4/hg008/engine.bed "HG008-T|HG008-N" -th 4
    bash harness/concordance_pair.sh hg008-truth \
      scratch/g4/hg008/cmp/java.vcf \
      scratch/g4/hg008/cmp/rust.vcf \
      scratch/g4/hg008/truth/HG008-T_somatic_smvar_benchmark_v0.3_tumorvariants.vcf.gz \
      scratch/g4/hg008/eval.bed \
      scratch/g4/ref/giabv3/GRCh38_GIABv3.sdf \
      scratch/g4/hg008/truth \
      HG008-T \
      scratch/g4/hg008/cmp/cpp.vcf
    echo "$(ts) [DONE] hg008-tier2 & hg008-truth ($((SECONDS - step_t0))s)"
else
    echo "$(ts) [REUSE] hg008-tier2 & hg008-truth: already certified for $FP"
fi

if needs_stamp "seqc2-slices-tier2"; then
    step_t0=$SECONDS
    echo "$(ts) [RUN] seqc2-slices-tier2: Certifying SEQC2 Somatic Slices..."
    KEEP_DIR=scratch/g4/seqc2/cmp-slices MIN_VARIANTS=4151 bash harness/compare.sh \
      scratch/g4/ref/d1vd1/GRCh38.d1.vd1.fa \
      "scratch/g4/seqc2/HCC1395.bam|scratch/g4/seqc2/HCC1395BL.bam" \
      scratch/g4/seqc2/slices.engine.bed \
      "HCC1395|HCC1395BL" -th 8
    python3 harness/record_evidence.py seqc2-slices-tier2 MATCH scratch/g4/seqc2/cmp-slices \
      -- env KEEP_DIR=scratch/g4/seqc2/cmp-slices MIN_VARIANTS=4151 bash harness/compare.sh scratch/g4/ref/d1vd1/GRCh38.d1.vd1.fa "scratch/g4/seqc2/HCC1395.bam|scratch/g4/seqc2/HCC1395BL.bam" scratch/g4/seqc2/slices.engine.bed "HCC1395|HCC1395BL" -th 8
    echo "$(ts) [DONE] seqc2-slices-tier2 ($((SECONDS - step_t0))s)"
else
    echo "$(ts) [REUSE] seqc2-slices-tier2: already certified for $FP"
fi

if needs_stamp "seqc2-exome-tier2" || needs_stamp "seqc2-truth"; then
    step_t0=$SECONDS
    echo "$(ts) [RUN] seqc2-exome-tier2 & seqc2-truth: Certifying SEQC2 Somatic Exome..."
    KEEP_DIR=scratch/g4/seqc2/cmp-exome MIN_VARIANTS=623444 bash harness/compare.sh \
      scratch/g4/ref/d1vd1/GRCh38.d1.vd1.fa \
      "scratch/g4/seqc2/HCC1395.bam|scratch/g4/seqc2/HCC1395BL.bam" \
      scratch/g4/seqc2/exome.engine.bed \
      "HCC1395|HCC1395BL" -th 4
    python3 harness/record_evidence.py seqc2-exome-tier2 MATCH scratch/g4/seqc2/cmp-exome \
      -- env KEEP_DIR=scratch/g4/seqc2/cmp-exome MIN_VARIANTS=623444 bash harness/compare.sh scratch/g4/ref/d1vd1/GRCh38.d1.vd1.fa "scratch/g4/seqc2/HCC1395.bam|scratch/g4/seqc2/HCC1395BL.bam" scratch/g4/seqc2/exome.engine.bed "HCC1395|HCC1395BL" -th 4
    bash harness/concordance_pair.sh seqc2-truth \
      scratch/g4/seqc2/cmp-exome/java.vcf \
      scratch/g4/seqc2/cmp-exome/rust.vcf \
      scratch/g4/seqc2/truth/truth.vcf.gz \
      scratch/g4/seqc2/eval.bed \
      scratch/g4/ref/d1vd1/GRCh38.d1.vd1.sdf \
      scratch/g4/seqc2/truth \
      HCC1395 \
      scratch/g4/seqc2/cmp-exome/cpp.vcf
    echo "$(ts) [DONE] seqc2-exome-tier2 & seqc2-truth ($((SECONDS - step_t0))s)"
else
    echo "$(ts) [REUSE] seqc2-exome-tier2 & seqc2-truth: already certified for $FP"
fi

if needs_stamp "hg002-tier2" || needs_stamp "hg002-truth"; then
    step_t0=$SECONDS
    echo "$(ts) [RUN] hg002-tier2 & hg002-truth: Certifying GIAB HG002..."
    KEEP_DIR=scratch/g4/hg002/cmp MIN_VARIANTS=227562 bash harness/compare.sh \
      scratch/g4/ref/giabv3/GRCh38_GIABv3.fasta \
      scratch/g4/hg002/HG002.chr20.bam \
      scratch/g4/hg002/eval.bed \
      HG002 -f 0.01 -th 8
    python3 harness/record_evidence.py hg002-tier2 MATCH scratch/g4/hg002/cmp \
      -- env KEEP_DIR=scratch/g4/hg002/cmp MIN_VARIANTS=227562 bash harness/compare.sh scratch/g4/ref/giabv3/GRCh38_GIABv3.fasta scratch/g4/hg002/HG002.chr20.bam scratch/g4/hg002/eval.bed HG002 -f 0.01 -th 8
    bash harness/concordance_germline.sh hg002-truth \
      scratch/g4/hg002/cmp/java.vcf \
      scratch/g4/hg002/cmp/rust.vcf \
      scratch/g4/hg002/truth/truth.vcf.gz \
      scratch/g4/hg002/eval.bed \
      scratch/g4/ref/giabv3/GRCh38_GIABv3.sdf \
      scratch/g4/hg002/truth \
      HG002 \
      scratch/g4/hg002/cmp/cpp.vcf
    echo "$(ts) [DONE] hg002-tier2 & hg002-truth ($((SECONDS - step_t0))s)"
else
    echo "$(ts) [REUSE] hg002-tier2 & hg002-truth: already certified for $FP"
fi

if needs_stamp "hg008-cram-twin"; then
    step_t0=$SECONDS
    echo "$(ts) [RUN] hg008-cram-twin: Certifying HG008 CRAM Twin..."
    KEEP_DIR=scratch/g4/hg008/cram CMP_TIMEOUT=12h bash harness/cram_twin.sh
    python3 harness/record_evidence.py hg008-cram-twin MATCH scratch/g4/hg008/cram \
      -- env KEEP_DIR=scratch/g4/hg008/cram CMP_TIMEOUT=12h bash harness/cram_twin.sh
    echo "$(ts) [DONE] hg008-cram-twin ($((SECONDS - step_t0))s)"
else
    echo "$(ts) [REUSE] hg008-cram-twin: already certified for $FP"
fi

if needs_stamp "extended-coverage-tier2"; then
    echo "$(ts) [RUN] extended-coverage-tier2: Running Extended Coverage Suite (17 modalities)..."
    bash "$ROOT_DIR/harness/run_extended_coverage.sh"
else
    echo "$(ts) [REUSE] extended-coverage-tier2: already certified for $FP"
fi

if [ ! -d "$DATA_ROOT/scratch/g4/multimodal" ]; then
    echo "ERROR: scratch/g4/multimodal is required for full recertification." >&2
    exit 1
fi
if needs_stamp "multimodal-slices-tier2"; then
    echo "$(ts) [RUN] multimodal-slices-tier2: Running Multi-Modal Multi-Slice Suite (18 biological slices)..."
    bash "$ROOT_DIR/harness/run_multimodal_slices.sh"
else
    echo "$(ts) [REUSE] multimodal-slices-tier2: already certified for $FP"
fi

if needs_stamp "hg005-truth" || needs_stamp "na24385-b37-truth"; then
    step_t0=$SECONDS
    echo "$(ts) [RUN] hg005-truth & na24385-b37-truth: Certifying Multi-Modal Truth Concordance..."
    bash "$ROOT_DIR/harness/concordance_multimodal.sh"
    echo "$(ts) [DONE] hg005-truth & na24385-b37-truth ($((SECONDS - step_t0))s)"
else
    echo "$(ts) [REUSE] hg005-truth & na24385-b37-truth: already certified for $FP"
fi

if needs_stamp "perf-budget" || needs_stamp "hg008-tier2" || needs_stamp "seqc2-exome-tier2" || needs_stamp "hg002-tier2" || needs_stamp "amplicon-tier2"; then
    echo "$(ts) [RUN] perf-budget: Certifying 4-Target Performance Budget..."
    mkdir -p scratch/g4/perf-budget
    python3 harness/perf_budget.py hg008-tier2 seqc2-exome-tier2 hg002-tier2 amplicon-tier2 > scratch/g4/perf-budget/extra.json
    python3 harness/record_evidence.py perf-budget PASS scratch/g4/perf-budget \
      -- python3 harness/perf_budget.py hg008-tier2 seqc2-exome-tier2 hg002-tier2 amplicon-tier2
else
    echo "$(ts) [REUSE] perf-budget: already certified for $FP"
fi

if needs_stamp "musl-ci"; then
    echo "$(ts) [RUN] musl-ci: Certifying Static musl Binary..."
    cargo test --target x86_64-unknown-linux-musl
    python3 harness/record_evidence.py musl-ci PASS - -- cargo test --target x86_64-unknown-linux-musl
else
    echo "$(ts) [REUSE] musl-ci: already certified for $FP"
fi

echo "$(ts) [RUN] Micro vardictcpp third-arm cell (advisory; gate stays java==rust) ---"
KEEP_DIR="$ROOT_DIR/scratch/recertify-micro" MIN_VARIANTS=1 \
  bash "$ROOT_DIR/harness/compare.sh" fixtures/micro/micro.ref.fa \
  fixtures/micro/micro.bam fixtures/micro/micro.bed

echo "$(ts) --- Golden store: pinned-baseline reuse in this run ---"
bash "$ROOT_DIR/harness/golden_store.sh" stats --run "$GOLDEN_RUN_ID" ${GOLDEN_MIN_HIT:+--min-hit "$GOLDEN_MIN_HIT"}
# Unset GOLDEN_RUN_ID before acceptance.sh so item 3 stamps each suite's own 33-cell ledger
unset GOLDEN_RUN_ID

echo "$(ts) [Recertify 4/5] Synchronizing documentation fingerprint & SHA-256 manifest tables..."
python3 "$ROOT_DIR/harness/sync_evidence_docs.py" --repo-root "$ROOT_DIR"

echo "$(ts) [Recertify 5/5] Running Master Acceptance Gate (8/8)..."
bash "$ROOT_DIR/harness/acceptance.sh"

echo "================================================================================"
echo " $(ts) RELEASE RECERTIFICATION COMPLETE ($((SECONDS - RECERT_T0))s): All gates PASS under fingerprint $FP"
echo "================================================================================"
