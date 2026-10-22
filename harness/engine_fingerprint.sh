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
  harness/cram_twin.sh
  harness/staple_cpp_extra.sh
  harness/golden_store.sh
  harness/golden_store.py
  benchmarks/env/vardict-oracle.yml
  harness/setup_oracle.sh
  harness/compare_vcf.py
  harness/run_extended_coverage.sh
  harness/run_multimodal_slices.sh
  harness/slices.tsv
  fixtures/micro
  VarDictJava
)

# --perf: only what can change speed or RSS (engine sources, manifests, lockfile) — the scale-bench carry key.
if [ "${1:-}" = "--perf" ]; then
  GATE_FILES=(Cargo.toml Cargo.lock crates/vardictrust-core/Cargo.toml crates/vardictrust/Cargo.toml crates/vardictrust-core/src crates/vardictrust/src)
fi

git diff --quiet HEAD -- "${GATE_FILES[@]}" \
  || { echo "engine_fingerprint: uncommitted changes under engine sources, manifests, comparators, runners, fixtures, or oracle" >&2; exit 1; }

normalize_gate_blob() {
  case "$1" in
    Cargo.toml)
      python3 -c 'import re,sys; s=sys.stdin.read(); s,n=re.subn(r"(\[workspace\.package\]\nversion = \")[0-9]+\.[0-9]+\.[0-9]+(?:-[^\"]+)?(?:\+[^\"]+)?(\")", r"\g<1>1.0.0\2", s, count=1); assert n==1; sys.stdout.write(s)'
      ;;
    crates/vardictrust/Cargo.toml)
      python3 -c 'import re,sys; s=sys.stdin.read(); s,n=re.subn(r"(vardictrust-core = \{ version = \")[0-9]+\.[0-9]+\.[0-9]+(?:-[^\"]+)?(?:\+[^\"]+)?(\", path = \"\.\./vardictrust-core\" \})", r"\g<1>1.0.0\2", s, count=1); assert n==1; sys.stdout.write(s)'
      ;;
    Cargo.lock)
      python3 -c 'import re,sys; s=sys.stdin.read(); s,n=re.subn(r"(\[\[package\]\]\nname = \"vardictrust(?:-core)?\"\nversion = \")[0-9]+\.[0-9]+\.[0-9]+(?:-[^\"]+)?(?:\+[^\"]+)?(\"\ndependencies = \[)", r"\g<1>1.0.0\2", s, count=2); assert n==2; sys.stdout.write(s)'
      ;;
    *)
      cat
      ;;
  esac
}

blob_oid() {
  case "$1" in
    Cargo.toml|crates/vardictrust/Cargo.toml|Cargo.lock)
      git show "HEAD:$1" | normalize_gate_blob "$1" | git hash-object -t blob --stdin
      ;;
    *)
      git rev-parse "HEAD:$1"
      ;;
  esac
}

{
  for f in "${GATE_FILES[@]}"; do
    blob_oid "$f"
  done
  cargo metadata --locked --format-version 1 | python3 -c 'import json,sys
m=json.load(sys.stdin); print("\n".join(sorted(p["name"]+"@"+p["version"] for p in m["packages"] if p.get("source"))))'
} | sha256sum | cut -c1-16

