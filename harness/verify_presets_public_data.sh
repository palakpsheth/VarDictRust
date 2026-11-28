#!/usr/bin/env bash
# =============================================================================
# harness/verify_presets_public_data.sh
# Compatibility shim delegating clinical workflow preset verification
# (--suite presets: --preset wgs|wes|amplicon|ctdna) to harness/verify_public_data.py.
# Provenance: Gong et al., Sci. Data 2022 (DOI: 10.1038/s41597-022-01359-6);
# Fang et al., Nat. Biotechnol. 2021 (DOI: 10.1038/s41587-021-00993-6).
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(CDPATH='' cd -- "$(dirname -- "$0")" && pwd)"
exec python3 "$SCRIPT_DIR/verify_public_data.py" --suite presets "$@"
