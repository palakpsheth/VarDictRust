#!/usr/bin/env bash
# =============================================================================
# harness/verify_wave5_public_data.sh
# Compatibility shim delegating Wave-5 public clinical & validator verification
# (--suite wave5) to the single-pass omnibus runner harness/verify_public_data.py.
# Provenance: clintval/vartovcf#161; Cleary et al., bioRxiv 2015 (DOI: 10.1101/023754);
# Fang et al., Nat. Biotechnol. 2021 (DOI: 10.1038/s41587-021-00993-6).
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(CDPATH='' cd -- "$(dirname -- "$0")" && pwd)"
exec python3 "$SCRIPT_DIR/verify_public_data.py" --suite wave5 "$@"
