#!/usr/bin/env bash
# harness/tests/test_golden_wiring.sh — no harness script may invoke an oracle without going through
# harness/golden_store.sh (issue #31); new suites are covered automatically (every harness/*.sh is scanned).
# Exempt per line: comments, `--version` probes, golden_store.sh calls, and `# golden:bypass <reason>`.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
# Scripts that legitimately name the oracles: the wrappers, the fingerprint file list, and trace_compare.sh.
ALLOW=(run_oracle.sh run_cpp.sh engine_fingerprint.sh trace_compare.sh golden_store.sh)
# Suites whose wiring task has not landed yet — each wiring task deletes its own entry; empty when the branch is done.
PENDING=()
in_list() { local x=$1; shift; local y; for y in "$@"; do [ "$x" = "$y" ] && return 0; done; return 1; }
scan() {
  grep -nE '(^|[^A-Za-z_./-])bash[[:space:]]+"?[^ ]*run_(oracle|cpp)\.sh|^[[:space:]]*"?[^ ]*harness/run_(oracle|cpp)\.sh|conda run[^|;&]*vardict(-java|cpp)([[:space:]]|$)|\$VARDICTCPP_BIN' "$1" \
    | grep -vE '^[0-9]+:[[:space:]]*#' | grep -vE 'golden_store\.sh|golden:bypass|--version' || true
}
fail=0
for f in harness/*.sh; do
  b=$(basename "$f")
  in_list "$b" "${ALLOW[@]}" && continue
  in_list "$b" "${PENDING[@]}" && continue
  bad=$(scan "$f")
  if [ -n "$bad" ]; then echo "FAIL: un-wired oracle invocation in $f:"; echo "$bad"; fail=1; fi
done
grep -q 'GOLDEN_LABEL="slice-' harness/run_slice.sh || { echo "FAIL: run_slice.sh lost its GOLDEN_LABEL"; fail=1; }
grep -q 'golden_store.sh cpp cram-twin-bam' harness/cram_twin.sh || { echo "FAIL: cram_twin.sh cpp arm not wired"; fail=1; }
grep -q 'GOLDEN_STORE_DIR' harness/run_self_tests.sh || { echo "FAIL: run_self_tests.sh does not isolate the store"; fail=1; }
grep -q '"golden": \$GOLDEN_JSON' harness/run_extended_coverage.sh || { echo "FAIL: extended evidence lost its golden block"; fail=1; }
grep -q "'golden': json.loads" harness/run_multimodal_slices.sh || { echo "FAIL: multimodal evidence lost its golden block"; fail=1; }
for g in harness/golden_store.sh harness/golden_store.py; do
  grep -qE "^  ${g}\$" harness/engine_fingerprint.sh || { echo "FAIL: $g is not in engine_fingerprint.sh GATE_FILES"; fail=1; }
done
grep -q 'golden_store.sh verify' harness/acceptance.sh || { echo "FAIL: acceptance.sh item 3 lost the store verify"; fail=1; }
grep -q 'GOLDEN_RUN_ID' harness/recertify.sh || { echo "FAIL: recertify.sh lost its hit-rate print"; fail=1; }
pf=$(bash harness/engine_fingerprint.sh --perf 2>/dev/null) || pf=dirty
if [ "$pf" != dirty ]; then
  { [[ "$pf" =~ ^[0-9a-f]{16}$ ]] && [ "$pf" != "$(bash harness/engine_fingerprint.sh 2>/dev/null)" ]; } \
    || { echo "FAIL: --perf fingerprint malformed or equal to the full fingerprint"; fail=1; }
fi
[ -z "$(git ls-files docs/release/waivers.md)" ] || { echo "FAIL: docs/release/waivers.md still exists"; fail=1; }
! grep -q 'waived' harness/acceptance.sh || { echo "FAIL: acceptance.sh still has a waiver path"; fail=1; }
stale=$(git grep -nIl 'release/waivers' -- . ':!VarDictJava' ':!harness/tests/test_golden_wiring.sh' || true)
[ -z "$stale" ] || { echo "FAIL: references to the deleted waiver registry remain: $stale"; fail=1; }
# Task 10b: no bare 'cpp engine exit $VAR' or 'engine failure (exit $VAR)' without stderr reason interpolation
bare_err=$(grep -nE 'cpp engine exit \$[A-Za-z0-9_]+\)|engine failure \(exit \$[A-Za-z0-9_]+\)' harness/*.sh || true)
[ -z "$bare_err" ] || { echo "FAIL: bare cpp exit reason without stderr detail: $bare_err"; fail=1; }
grep -q 'golden_store.sh" cpp ext-m2' harness/run_extended_coverage.sh || { echo "FAIL: ext-m2 cpp arm missing"; fail=1; }
grep -q 'golden_store.sh" cpp ext-m3' harness/run_extended_coverage.sh || { echo "FAIL: ext-m3 cpp arm missing"; fail=1; }
grep -q 'golden_store.sh" cpp "mm-\$slice_id"' harness/run_multimodal_slices.sh || { echo "FAIL: multimodal vcf cpp arm missing"; fail=1; }
grep -q 'cpp_f1_normalized' harness/concordance_germline.sh || { echo "FAIL: concordance_germline missing cpp_f1_normalized"; fail=1; }
grep -q 'cpp_f1_normalized' harness/concordance_pair.sh || { echo "FAIL: concordance_pair missing cpp_f1_normalized"; fail=1; }

[ $fail -eq 0 ] && echo "PASS test_golden_wiring"; exit $fail
