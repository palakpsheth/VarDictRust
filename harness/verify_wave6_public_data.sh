#!/usr/bin/env bash
# =============================================================================
# harness/verify_wave6_public_data.sh
# Compatibility shim delegating Wave-6 38-switch public truth-set & invariant
# verification (--suite wave6) to the single-pass omnibus harness/verify_public_data.py.
# Provenance: Tan et al., Bioinformatics 2015 (DOI: 10.1093/bioinformatics/btv112);
# Fang et al., Nat. Biotechnol. 2021 (DOI: 10.1038/s41587-021-00993-6).
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(CDPATH='' cd -- "$(dirname -- "$0")" && pwd)"
exec python3 "$SCRIPT_DIR/verify_public_data.py" --suite wave6 "$@"
