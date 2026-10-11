#!/usr/bin/env bash
# Content hash of everything that can change engine output or the Tier-2 gate.
set -euo pipefail
git diff --quiet HEAD -- crates harness Cargo.lock \
  || { echo "engine_fingerprint: uncommitted changes under crates/ harness/ Cargo.lock" >&2; exit 1; }
{ git rev-parse HEAD:crates/vardict-core/src HEAD:crates/vardictrust/src \
    HEAD:harness/compare.sh HEAD:harness/run_oracle.sh HEAD:harness/compare_vcf.py
  cargo metadata --locked --format-version 1 | python3 -c 'import json,sys
m=json.load(sys.stdin); print("\n".join(sorted(p["name"]+"@"+p["version"] for p in m["packages"] if p.get("source"))))'
} | sha256sum | cut -c1-16
