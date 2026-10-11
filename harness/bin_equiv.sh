#!/usr/bin/env bash
# harness/bin_equiv.sh
# Validates that the static musl binary (target/x86_64-unknown-linux-musl/release/vardictrust)
# produces byte-identical VCF output to the standard glibc binary across representative benchmarks.
set -euo pipefail

GNU_BIN="target/release/vardictrust"
MUSL_BIN="target/x86_64-unknown-linux-musl/release/vardictrust"

echo "Building GNU binary..."
cargo build --release -p vardictrust

echo "Building MUSL binary..."
cargo build --release --target x86_64-unknown-linux-musl -p vardictrust

echo "=== 1. Verifying MUSL Static Linkage ==="
file_info=$(file "$MUSL_BIN")
echo "$file_info"
if echo "$file_info" | grep -E -q "static.*linked" || ldd "$MUSL_BIN" 2>&1 | grep -q "statically linked"; then
  echo "✓ Confirmed: x86_64-unknown-linux-musl binary is statically linked (no glibc dependency)."
else
  echo "ERROR: $MUSL_BIN is not statically linked!" >&2
  exit 1
fi

TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TMP_DIR"' EXIT

echo "=== 2. Micro Fixture Byte-for-Byte Comparison ==="
"$GNU_BIN" -G fixtures/micro/micro.ref.fa -b fixtures/micro/micro.bam -c 1 -S 2 -E 3 -g 4 -z fixtures/micro/micro.bed > "$TMP_DIR/gnu_micro.vcf"
"$MUSL_BIN" -G fixtures/micro/micro.ref.fa -b fixtures/micro/micro.bam -c 1 -S 2 -E 3 -g 4 -z fixtures/micro/micro.bed > "$TMP_DIR/musl_micro.vcf"

diff -u "$TMP_DIR/gnu_micro.vcf" "$TMP_DIR/musl_micro.vcf" || {
  echo "ERROR: Micro fixture VCF differs between GNU and MUSL binaries!" >&2
  exit 2
}
echo "✓ Micro fixture: 100% byte-identical."

echo "=== 3. Representative Slices Equivalence Sweep ==="
SLICES=("s1a" "tiny" "s3" "s6" "s9")
for s in "${SLICES[@]}"; do
  row=$(awk -F'\t' -v id="$s" '$1==id {print $0}' harness/slices.tsv)
  if [ -z "$row" ]; then
    echo "Skipping unknown slice $s"
    continue
  fi
  IFS=$'\t' read -r id ref bam bed sample min_variants extra <<< "$row"
  echo "Testing slice: $id..."
  
  "$GNU_BIN" -G "$ref" -b "$bam" -N "$sample" -c 1 -S 2 -E 3 -g 4 -z ${extra:-} "$bed" > "$TMP_DIR/gnu_${id}.vcf"
  "$MUSL_BIN" -G "$ref" -b "$bam" -N "$sample" -c 1 -S 2 -E 3 -g 4 -z ${extra:-} "$bed" > "$TMP_DIR/musl_${id}.vcf"
  
  diff -u "$TMP_DIR/gnu_${id}.vcf" "$TMP_DIR/musl_${id}.vcf" || {
    echo "ERROR: Slice $id VCF differs between GNU and MUSL binaries!" >&2
    exit 3
  }
  echo "✓ Slice $id: 100% byte-identical."
done

echo "================================================================="
echo "=== MUSL STATIC BINARY EQUIVALENCE: 100% PASS (BYTE-IDENTICAL) ==="
echo "================================================================="
