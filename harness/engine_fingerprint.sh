#!/usr/bin/env bash
# Content hash of everything that can change engine output or the Tier-2 gate.
set -euo pipefail
GATE_FILES=(
  Cargo.toml
  Cargo.lock
  crates/vardictrust-core/Cargo.toml
  crates/vardictrust/Cargo.toml
  crates/vardictrust-core/src
  crates/vardictrust/src
  harness/compare.sh
  harness/run_oracle.sh
  harness/run_cpp.sh
  benchmarks/env/vardict-oracle.yml
  harness/setup_oracle.sh
  harness/compare_vcf.py
  harness/run_extended_coverage.sh
  harness/run_multimodal_slices.sh
  harness/slices.tsv
  fixtures/micro
  VarDictJava
)

git diff --quiet HEAD -- "${GATE_FILES[@]}" \
  || { echo "engine_fingerprint: uncommitted changes under engine sources, manifests, comparators, runners, fixtures, or oracle" >&2; exit 1; }

{
  for f in "${GATE_FILES[@]}"; do
    git rev-parse "HEAD:$f"
  done
  cargo metadata --locked --format-version 1 | python3 -c 'import json,sys
m=json.load(sys.stdin); print("\n".join(sorted(p["name"]+"@"+p["version"] for p in m["packages"] if p.get("source"))))'
} | sha256sum | cut -c1-16
