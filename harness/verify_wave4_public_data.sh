#!/usr/bin/env bash
# =============================================================================
# harness/verify_wave4_public_data.sh
# Compatibility shim delegating Wave-4 public clinical verification (--suite wave4)
# to the single-pass omnibus runner harness/verify_public_data.py.
# Provenance: Gong et al., Sci. Data 2022 (DOI: 10.1038/s41597-022-01359-6);
# bioconda/bioconda-recipes#71875; Fang et al., Nat. Biotechnol. 2021
# (DOI: 10.1038/s41587-021-00993-6).
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(CDPATH='' cd -- "$(dirname -- "$0")" && pwd)"
exec python3 "$SCRIPT_DIR/verify_public_data.py" --suite wave4 "$@"
