#!/usr/bin/env bash
# =============================================================================
# harness/test_bioconda_recipe.sh
# Local & containerized builder/tester for the vardictrust Bioconda recipe.
#
# Adheres strictly to AGENTS.md Rule 0/0b:
#   - Runs real verification commands.
#   - Never outputs a fake MATCH or mocks a passing build.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
RECIPE_DIR="$ROOT_DIR/packaging/bioconda/recipes/vardictrust"

echo "================================================================================"
echo " VarDictRust Bioconda Recipe Verification"
echo " Recipe Path: $RECIPE_DIR"
echo "================================================================================"

# 1. Run in-tree Python validator (static schema, required sections, licenses)
python3 "$SCRIPT_DIR/tests/test_bioconda_recipe.py"
echo "[PASS] In-tree recipe structure, metadata, and license checks passed."

# 2. Check for Docker for real containerized build
if [ "${1:-}" = "--docker" ] || [ "${DOCKER_BUILD:-0}" = "1" ]; then
  if command -v docker &>/dev/null && docker info &>/dev/null; then
    echo "[INFO] Running containerized bioconda build test:"
    docker run --rm -v "$ROOT_DIR:/workspace" -w /workspace quay.io/bioconda/bioconda-utils-build-env:latest \
      bioconda-utils build --docker --packages vardictrust
    echo "[PASS] Containerized bioconda build succeeded."
  else
    echo "ERROR: --docker requested but docker daemon is not available." >&2
    exit 1
  fi
elif command -v conda-build &>/dev/null; then
  echo "[INFO] conda-build detected. Building recipe locally:"
  conda-build "$RECIPE_DIR" --no-anaconda-upload
  echo "[PASS] Local conda-build succeeded."
else
  echo "[INFO] Container build skipped (run with --docker to execute full container build)."
  echo "[PASS] Static schema, licensing, and file contracts verified."
fi

echo "================================================================================"
echo " Recipe Verification Check Complete"
echo "================================================================================"
