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
grep -qx MATCH "$T/o2" && [ -f "$T/k2/java.time" ] || { echo "FAIL2 output"; fail=1; }
MIN_VARIANTS=100000 bash harness/compare.sh $M >/dev/null 2>&1 && { echo "FAIL3 floor ignored"; fail=1; }
mkdir -p "$T/tmp"; PATH="$T/bin:$PATH" TMPDIR="$T/tmp" bash harness/compare.sh $M >/dev/null 2>&1
[ -z "$(ls -A "$T/tmp")" ] || { echo "FAIL4 temp dir leaked"; fail=1; }
# CPP5: cpp missing -> SKIP, gate verdict intact (stub conda fails vardictcpp only).
REAL_CONDA=$(command -v conda)
mkdir -p "$T/nocppbin"; cat > "$T/nocppbin/conda" <<EOF
#!/usr/bin/env bash
for a in "\$@"; do [ "\$a" = vardictcpp ] && { echo "stub conda: no vardictcpp" >&2; exit 1; }; done
exec "$REAL_CONDA" "\$@"
EOF
chmod +x "$T/nocppbin/conda"
out5=$(PATH="$T/nocppbin:$PATH" KEEP_DIR="$T/k5" bash harness/compare.sh $M 2>&1); rc5=$?
[ $rc5 -eq 0 ] || { echo "FAIL5 gate rc=$rc5"; fail=1; }
grep -qx MATCH <<<"$out5" || { echo "FAIL5 no gate MATCH"; fail=1; }
grep -q '^CPP: SKIP:' <<<"$out5" || { echo "FAIL5 no CPP SKIP"; fail=1; }
grep -q '^CPP: SKIP:' "$T/k5/compare_cpp.txt" || { echo "FAIL5 no compare_cpp SKIP"; fail=1; }
# CPP6: empty BED + MIN_VARIANTS=1 -> CPP MISMATCH (never vacuous), header-only fallback.
: > "$T/empty.bed"
E="fixtures/micro/micro.ref.fa fixtures/micro/micro.bam $T/empty.bed"
out6=$(MIN_VARIANTS=1 KEEP_DIR="$T/k6" bash harness/compare.sh $E 2>&1); rc6=$?
[ $rc6 -ne 0 ] || { echo "FAIL6 rc=0"; fail=1; }
grep -q '^CPP: MISMATCH' <<<"$out6" || { echo "FAIL6 no CPP MISMATCH"; fail=1; }
[ -s "$T/k6/cpp.vcf" ] || { echo "FAIL6 cpp.vcf not header-only"; fail=1; }
[ "$(grep -vc '^#' "$T/k6/cpp.vcf")" -eq 0 ] || { echo "FAIL6 cpp.vcf has records"; fail=1; }
grep -q '^#CHROM' "$T/k6/cpp.vcf" || { echo "FAIL6 cpp.vcf no #CHROM"; fail=1; }
# CPP7: micro third arm records a verdict; cpp agrees on calls but its native
# -vcf serialization differs (no ##reference/##contig headers, extra END INFO
# key), so verbatim compare_vcf.py reports MISMATCH — recorded, never gating.
out7=$(MIN_VARIANTS=1 KEEP_DIR="$T/k7" bash harness/compare.sh $M 2>&1); rc7=$?
[ $rc7 -eq 0 ] || { echo "FAIL7 gate rc=$rc7"; fail=1; }
grep -qx MATCH <<<"$out7" || { echo "FAIL7 no gate MATCH"; fail=1; }
grep -q '^CPP: MISMATCH' <<<"$out7" || { echo "FAIL7 no CPP MISMATCH"; fail=1; }
grep -q 'END' "$T/k7/compare_cpp.txt" || { echo "FAIL7 no END diff in compare_cpp"; fail=1; }
[ -s "$T/k7/cpp.vcf" ] && [ -f "$T/k7/cpp.time" ] || { echo "FAIL7 cpp artifacts"; fail=1; }
# CPP8: CPP_SKIP=1 escapes the third arm; gate verdict unchanged.
out8=$(CPP_SKIP=1 KEEP_DIR="$T/k8" bash harness/compare.sh $M 2>&1); rc8=$?
[ $rc8 -eq 0 ] || { echo "FAIL8 gate rc=$rc8"; fail=1; }
grep -q '^CPP: SKIP:' <<<"$out8" || { echo "FAIL8 no CPP SKIP"; fail=1; }
[ ! -e "$T/k8/cpp.vcf" ] || { echo "FAIL8 cpp ran despite SKIP"; fail=1; }
# CPP9: first-comparator MISMATCH + second MATCH must still report CPP: MISMATCH
# (both-rc regression for a7854d5; single-rc code saw only rc2 and misreported MATCH).
# Stub python3 synthesizes java!=cpp (rc1=1) and cpp==rust (rc2=0); the gate call
# (java+rust, no cpp) delegates to the real python3 so the gate stays honest.
REAL_PY3=$(command -v python3)
mkdir -p "$T/pybin"; cat > "$T/pybin/python3" <<EOF
#!/usr/bin/env bash
for a in "\$@"; do case "\$a" in *java.vcf) hj=1;; *cpp.vcf) hc=1;; *rust.vcf) hr=1;; esac; done
if [ "\${hj:-0}" -eq 1 ] && [ "\${hc:-0}" -eq 1 ]; then echo "MISMATCH: synthetic java!=cpp (CPP9)"; exit 1; fi
if [ "\${hc:-0}" -eq 1 ] && [ "\${hr:-0}" -eq 1 ]; then echo "MATCH (synthetic cpp==rust CPP9)"; exit 0; fi
exec "$REAL_PY3" "\$@"
EOF
chmod +x "$T/pybin/python3"
out9=$(PATH="$T/pybin:$PATH" MIN_VARIANTS=1 KEEP_DIR="$T/k9" bash harness/compare.sh $M 2>&1); rc9=$?
[ $rc9 -eq 0 ] || { echo "FAIL9 gate rc=$rc9"; fail=1; }
grep -qx MATCH <<<"$out9" || { echo "FAIL9 no gate MATCH"; fail=1; }
grep -q '^CPP: MISMATCH' <<<"$out9" || { echo "FAIL9 no CPP MISMATCH (rc1=1 rc2=0)"; fail=1; }
grep -q 'synthetic java!=cpp' "$T/k9/compare_cpp.txt" || { echo "FAIL9 no rc1 evidence"; fail=1; }
grep -q 'synthetic cpp==rust' "$T/k9/compare_cpp.txt" || { echo "FAIL9 no rc2 evidence"; fail=1; }
[ $fail -eq 0 ] && echo "PASS test_compare_sh"; exit $fail
