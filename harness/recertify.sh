#!/usr/bin/env bash
# harness/recertify.sh — Full pre-merge parity recertification & evidence updater
# Governed by rewrites.bio Manifesto §0, §0b, §5 and AGENTS.md §15.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

echo "================================================================================"
echo " VarDictRust Pre-Merge Live Dual-Engine Parity Recertification"
echo "================================================================================"

echo "[1/4] Building release binary..."
cargo build --release -p vardictrust

echo "[2/4] Materializing reference FASTAs..."
python3 "$ROOT_DIR/harness/convert_fastas.py"

FP=$(bash "$ROOT_DIR/harness/engine_fingerprint.sh")
echo "Current Engine Fingerprint: $FP"
export GOLDEN_RUN_ID="recertify-$(date +%s)"

echo "[3/4] Running live dual-engine regression suites..."
echo "--- Running Extended Coverage Suite (17 modalities)..."
bash "$ROOT_DIR/harness/run_extended_coverage.sh"

if [ ! -d "$ROOT_DIR/scratch/g4/multimodal" ]; then
    echo "ERROR: scratch/g4/multimodal is required for full pre-merge recertification." >&2
    exit 1
fi
echo "--- Running Multi-Modal Multi-Slice Suite (18 biological slices)..."
bash "$ROOT_DIR/harness/run_multimodal_slices.sh"

# Fingerprints asserted (not forced) by acceptance.sh ev(); runners stamp at
# generation (per maintainer ruling: overwriting a runner-stamped fingerprint
# is Class-0).
echo "--- Micro vardictcpp third-arm cell (advisory; gate stays java==rust) ---"
KEEP_DIR="$ROOT_DIR/scratch/recertify-micro" MIN_VARIANTS=1 \
  bash "$ROOT_DIR/harness/compare.sh" "$ROOT_DIR/fixtures/micro/micro.ref.fa" \
  "$ROOT_DIR/fixtures/micro/micro.bam" "$ROOT_DIR/fixtures/micro/micro.bed" micro

echo "--- Golden store: pinned-baseline reuse in this run ---"
bash "$ROOT_DIR/harness/golden_store.sh" stats --run "$GOLDEN_RUN_ID" ${GOLDEN_MIN_HIT:+--min-hit "$GOLDEN_MIN_HIT"}
# Unset GOLDEN_RUN_ID before acceptance.sh so item 3 stamps each suite's own 33-cell ledger
unset GOLDEN_RUN_ID

echo "[4/4] Running Master Acceptance Gate (8/8)..."
bash "$ROOT_DIR/harness/acceptance.sh"

echo "================================================================================"
echo " PRE-MERGE RECERTIFICATION COMPLETE: All gates PASS under fingerprint $FP"
echo "================================================================================"
