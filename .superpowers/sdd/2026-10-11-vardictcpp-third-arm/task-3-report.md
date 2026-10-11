# Task 3 report — `compare.sh` third arm (gate unchanged)

Commit: `e279407` (`T3A-03: compare.sh third-arm (gate unchanged)`, branch `feature/vardictcpp-reproducibility`).
Files: `harness/compare.sh` (+40), `harness/tests/test_compare_sh.sh` (+37/−1). No other files touched.

## What was built

`compare.sh` runs the gate exactly as before (oracle → rust → `compare_vcf.py`, `MATCH`/`MISMATCH`, `exit $rc`),
then runs an advisory vardictcpp arm that can never change the verdict or exit code:

- `CPP_SKIP=1` (or `true`) → `CPP: SKIP: disabled via CPP_SKIP=1`, no engine run.
- `conda` absent → `CPP: SKIP: conda not on PATH`.
- Otherwise `LC_ALL=C /usr/bin/time -v -o cpp.time timeout $TO bash harness/run_cpp.sh <ref> <bam> <bed> <cpp.vcf> <sample> -z "$@"`
  (`-z` added exactly like the oracle call so cpp runs zero-based; `-th`/paired `*|*` flow through `$@` untouched).
- `run_cpp.sh` exit 3 → `CPP: SKIP: <run_cpp.sh reason>`; any other failure → `ENGINE-FAIL: cpp` on stdout
  plus `CPP: SKIP: engine failure (exit N)`. Both keep the java/rust exit code.
- Success → `compare_vcf.py [$STRICT_ARG] --min-variants $MIN_VARIANTS` (verbatim flags, default eps 1e-6)
  run twice (java vs cpp, cpp vs rust) into `compare_cpp.txt`, then `CPP: MATCH` (both rc 0) else `CPP: MISMATCH`
  appended to both stdout and `compare_cpp.txt`.
- Carried Task-2 fix (ledger Important item): 0-byte `cpp.vcf` (empty BED / variant-free) gets a header-only
  fallback (`grep '^#' java.vcf`, minimal fileformat/`#CHROM` if that is ever empty) so `compare_vcf.py` sees a
  valid empty VCF. Verified live: empty BED → oracle/rust emit 65-line header-only VCFs (4901 B each).
- Per-`KEEP_DIR` artifacts: `cpp.vcf`, `cpp.time`, `cpp.err`, `compare_cpp.txt`. Without `KEEP_DIR` all vanish
  with the existing trap (FAIL4 temp-leak check still passes).
- Existing `FAIL2` assertion updated (`tail -1 == MATCH` → `grep -qx MATCH`) since stdout now ends with the
  `CPP:` line; gate intent unchanged.

## Tests (fail-first done)

Added CPP5–CPP8 to `harness/tests/test_compare_sh.sh`. Ran the file against unmodified `compare.sh` first:
all four new cases failed as required (no `CPP:` lines, no `cpp.vcf`/`compare_cpp.txt`); after the fix the
whole file prints `PASS test_compare_sh` (run twice: pre- and post-commit, second run is the committed state).

- CPP5 cpp-missing (stub `conda` failing only `vardictcpp` args, delegating the rest to real conda):
  gate rc 0 + `MATCH` intact, `CPP: SKIP:` on stdout and in `compare_cpp.txt`.
- CPP6 empty BED + `MIN_VARIANTS=1`: rc ≠ 0, `CPP: MISMATCH` (floor, never vacuous), `cpp.vcf` non-empty,
  zero records, `#CHROM` present (fallback path).
- CPP7 micro + `MIN_VARIANTS=1`: gate rc 0 + `MATCH`; `CPP: MISMATCH` recorded; `compare_cpp.txt` names `END`;
  `cpp.vcf`/`cpp.time` kept.
- CPP8 `CPP_SKIP=1`: gate rc 0, `CPP: SKIP:`, no `cpp.vcf` created.

## Micro gate (brief Step 4)

`MIN_VARIANTS=1 KEEP_DIR=<tmp> bash harness/compare.sh fixtures/micro/micro.ref.fa fixtures/micro/micro.bam
fixtures/micro/micro.bed SAMPLE -th 8` → stdout `MATCH (2 variants)` / `MATCH` / `CPP: MISMATCH`, rc 0.
`grep -E '^(MATCH|MISMATCH|CPP:)' $KEEP_DIR/compare*.txt` → `MATCH (2 variants)`, two `MISMATCH:` blocks,
`CPP: MISMATCH`.

## CONCERN (brief expectation corrected — NOT implemented as written)

The brief expected micro → `CPP: MATCH`. Live evidence proves that is unachievable with the mandated
verbatim `compare_vcf.py` (eps 1e-6, no normalization), so I recorded the honest `CPP: MISMATCH` instead of
manufacturing a MATCH (§0b). Characterization (record-by-record diff of `java.vcf` vs cpp output on micro):

- Variant sets identical (2/2, same CHROM/POS/REF/ALT); FILTER, QUAL, FORMAT, all sample fields, and every
  INFO value identical. cpp agrees on all calls.
- Exactly two serialization diffs, both structural to vardictcpp's native `-vcf` emitter:
  1. Headers: java 65 lines vs cpp 62 — cpp omits `##reference=` + 2× `##contig=` (the perl `var2vcf_valid.pl`
     provenance lines; note AGENTS.md §11/#320 context — cpp emits no contig lengths at all).
  2. Records: cpp INFO carries `END=101` / `END=153`; java/rust INFO key-set has no `END`
     (`compare_vcf.py` requires exact INFO key-set equality).
- No normalization was applied: `cpp.vcf` is byte-verbatim engine output and `compare_cpp.txt` shows both
  diffs, so Tasks 4–7 and any future concordance work consume truthful signal. If the plan wants a
  call-level `CPP: MATCH` later, that needs an explicit normalized-metric decision — out of scope for this task.

## Self-review

- Gate preamble byte-identical in behavior (oracle/rust ENGINE-FAIL exit-2 paths untouched; callers
  `acceptance.sh`/`pr_check.sh`/`rgate.sh`/`run_slice.sh` key off exit code, safe).
- `set -euo pipefail` audited: every fallible cpp-arm command is inside an `elif` condition, an `||` guard,
  or an `if` test. `bash -n` clean on both files; pre-commit (gitleaks) passed.
- No subagents used. No data/target/envs committed (`.opencode/` left untracked).

## For the ledger

- Task-2 carried Important item (0-byte cpp OUT) is FIXED by the header-only fallback in `compare.sh`.
- Downstream note for Task 7: `compare_cpp.txt` always exists when the arm runs or skips (SKIP lines included);
  `cpp.vcf`/`cpp.time` exist only when the engine was attempted. `extra.cpp_*` can parse the trailing
  `^CPP: ` line from stdout or `compare_cpp.txt` — same string in both.

## Fix round 1/5 — capture both cpp comparator exit codes (review Important item)

Fix commit: `a7854d5` (`T3A-03-fix1: capture both cpp comparator exit codes separately`).
Report commit: (this commit; see `git log --oneline -1`).

Finding (verbatim): the cpp arm block `{ python java-vs-cpp; echo; python cpp-vs-rust; }
>compare_cpp.txt || cpp_rc=$?` captured only the SECOND comparator's exit, so a
first-comparator MISMATCH with a second-comparator MATCH (reachable when the gate MISMATCHes,
e.g. cpp==rust≠java) wrongly reported `CPP: MATCH`.

Fix (`harness/compare.sh` only; gate preamble and everything else identical): run both
comparators as separate commands writing to the same `compare_cpp.txt` (`>` then `>>`, separator
`echo` in between), recording each exit separately (`|| rc1=$?` / `|| rc2=$?`, both
`set -e`-safe), and report `CPP: MATCH` only if `rc1 == 0 && rc2 == 0`, else `CPP: MISMATCH`.
Both comparators still always run; output ordering in `compare_cpp.txt` unchanged
(java-vs-cpp output, separator, cpp-vs-rust output, trailing `CPP:` line).

Tests: `bash -n harness/compare.sh` clean; full `bash harness/tests/test_compare_sh.sh`
prints `PASS test_compare_sh` (covers CPP5–CPP8 incl. CPP7 MISMATCH, gate rc intact).
No subagents used. Controller ruling noted: CPP7 asserting `CPP: MISMATCH` on micro is
correct as written — untouched by this fix.
