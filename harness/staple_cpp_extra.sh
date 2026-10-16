#!/usr/bin/env bash
# Usage: bash harness/staple_cpp_extra.sh <keepdir> <cpp_line>
# Merge-only persistence of third-arm evidence keys for record_evidence.py
# (ev_cpp presence gate): success verdicts -> extra.cpp_version, any SKIP ->
# extra.cpp_skip_reason. The two keys are mutually exclusive (stale one dropped
# so KEEP_DIR reuse never reports both). Never touches any gate exit code.
set -euo pipefail
OUT=$1; cpp_line=$2
case "$cpp_line" in
  "CPP: SKIP: "*)
    cpp_key="cpp_skip_reason"; cpp_val="${cpp_line#CPP: SKIP: }" ;;
  *)
    cpp_key="cpp_version"
    cpp_val=$(conda list -n vardict-oracle vardictcpp 2>/dev/null | awk '/^vardictcpp[[:space:]]/{print "vardictcpp " $2 " " $3 " bioconda"}')
    [ -n "$cpp_val" ] || cpp_val=$(conda run -n vardict-oracle vardictcpp --version 2>/dev/null | head -n 1) ;;
esac
if [ -n "${cpp_val:-}" ]; then
  # CPP-5 chunk pin: record the effective --chunk beside the version key
  # (protocol: default frozen in prereg, recorded in extra). run_cpp.sh logs
  # `cpp_chunk=<val>` to its stderr ($OUT/cpp.err under KEEP_DIR); fall back
  # to the same ${CPP_CHUNK:-50000} default the runner used.
  cpp_chunk=""
  if [ -f "$OUT/cpp.err" ]; then
    cpp_chunk=$(grep -m1 -o '^cpp_chunk=[0-9]*' "$OUT/cpp.err" | cut -d= -f2 || true)
  fi
  [ -n "$cpp_chunk" ] || cpp_chunk="${CPP_CHUNK:-50000}"
  python3 - "$OUT/extra.json" "$cpp_key" "$cpp_val" "$cpp_chunk" <<'PYEOF' || true
import json, sys
p, k, v, chunk = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
try:
    with open(p) as f: d = json.load(f)
    if not isinstance(d, dict): d = {}
except Exception: d = {}
if k and v:
    d[k] = v
    d.pop("cpp_version" if k == "cpp_skip_reason" else "cpp_skip_reason", None)
if chunk:
    d["cpp_chunk"] = chunk
with open(p, "w") as f: json.dump(d, f, indent=2, sort_keys=True)
PYEOF
fi
