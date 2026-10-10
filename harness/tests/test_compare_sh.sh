#!/usr/bin/env bash
set -uo pipefail
cargo build -q --release -p vardictrust
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT; fail=0
M="fixtures/micro/micro.ref.fa fixtures/micro/micro.bam fixtures/micro/micro.bed"
mkdir -p "$T/bin"; cat > "$T/bin/cargo" <<EOF
#!/usr/bin/env bash
while [ "\$1" != "--" ]; do shift; done; shift
"\$PWD/target/release/vardictrust" "\$@" | head -n -1
EOF
chmod +x "$T/bin/cargo"
out=$(PATH="$T/bin:$PATH" KEEP_DIR="$T/k1" bash harness/compare.sh $M 2>&1); rc=$?
[ $rc -ne 0 ] || { echo "FAIL1 rc=0"; fail=1; }
grep -qx MISMATCH <<<"$out" || { echo "FAIL1 no MISMATCH line"; fail=1; }
grep -q '^[<>]' <<<"$out" || { echo "FAIL1 no diff"; fail=1; }
[ -s "$T/k1/java.vcf" ] && [ -f "$T/k1/rust.vcf" ] || { echo "FAIL1 files not kept"; fail=1; }
[ "$(grep -vc '^#' "$T/k1/java.vcf")" -ge 1 ] || { echo "FAIL1 fixture vacuous"; fail=1; }
KEEP_DIR="$T/k2" MIN_VARIANTS=1 bash harness/compare.sh $M > "$T/o2" 2>&1 || { echo "FAIL2 rc"; fail=1; }
[ "$(tail -1 "$T/o2")" = MATCH ] && [ -f "$T/k2/java.time" ] || { echo "FAIL2 output"; fail=1; }
MIN_VARIANTS=100000 bash harness/compare.sh $M >/dev/null 2>&1 && { echo "FAIL3 floor ignored"; fail=1; }
mkdir -p "$T/tmp"; PATH="$T/bin:$PATH" TMPDIR="$T/tmp" bash harness/compare.sh $M >/dev/null 2>&1
[ -z "$(ls -A "$T/tmp")" ] || { echo "FAIL4 temp dir leaked"; fail=1; }
[ $fail -eq 0 ] && echo "PASS test_compare_sh"; exit $fail
