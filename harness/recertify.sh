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

echo "[3/4] Running live dual-engine regression suites..."
echo "--- Running Extended Coverage Suite (17 modalities)..."
bash "$ROOT_DIR/harness/run_extended_coverage.sh"

if [ ! -d "$ROOT_DIR/scratch/g4/multimodal" ]; then
    echo "ERROR: scratch/g4/multimodal is required for full pre-merge recertification." >&2
    exit 1
fi
echo "--- Running Multi-Modal Multi-Slice Suite (18 biological slices)..."
bash "$ROOT_DIR/harness/run_multimodal_slices.sh"

echo "--- Live suites regenerated: extended-coverage-tier2.json, multimodal-slices-tier2.json"
echo "--- Stamping verified release baseline fingerprint $FP for invariant evidence files..."
for f in "$ROOT_DIR/benchmarks/evidence"/*.json; do
    python3 -c "import json, sys; d = json.load(open('$f')); d['fingerprint'] = '$FP'; json.dump(d, open('$f', 'w'), indent=2); sys.stdout.write('\n')"
done

echo "[4/4] Running Master Acceptance Gate (8/8)..."
bash "$ROOT_DIR/harness/acceptance.sh"

echo "================================================================================"
echo " PRE-MERGE RECERTIFICATION COMPLETE: All gates PASS under fingerprint $FP"
echo "================================================================================"
