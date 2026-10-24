#!/usr/bin/env bash
# Thin shim so suites call `bash harness/golden_store.sh <verb> ...`; all logic is in golden_store.py (stdlib only).
exec python3 "$(dirname "${BASH_SOURCE[0]}")/golden_store.py" "$@"
