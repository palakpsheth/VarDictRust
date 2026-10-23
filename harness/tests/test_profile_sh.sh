#!/usr/bin/env bash
# harness/tests/test_profile_sh.sh: harness/profile.sh plan/guards/execution, with stub profiler binaries (no build, no
# real data run). Run from the repo root (run_self_tests.sh contract: <10 s, no persistent artifacts).
set -uo pipefail
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT; fail=0
PS=harness/profile.sh
ck() { if ! eval "$2"; then echo "FAIL test_profile_sh: $1" >&2; fail=1; fi; }

# stub engine: honours the VARDICT_PROFILE_* env contract of crates/vardictrust/src/profile.rs
cat > "$T/engine" <<'ENG'
#!/usr/bin/env bash
echo "$@" >> "$STUB_LOG"
if [ -n "${VARDICT_PROFILE_CPU:-}" ]; then
  printf 'main;run;vardictrust_core::stub::hot 5\n' > "$VARDICT_PROFILE_CPU/cpu-collapsed.txt"
  printf '{"ticks_total":6,"samples_kept":5}\n' > "$VARDICT_PROFILE_CPU/cpu-meta.json"
fi
if [ -n "${VARDICT_PROFILE_DHAT:-}" ]; then
  printf '{"pps":[{"tb":10,"tbk":1,"mb":10,"fs":[1,2]}],"ftbl":["[root]","a::alloc","vardictrust_core::stub::hot"]}\n' > "$VARDICT_PROFILE_DHAT/dhat-heap.json"
fi
printf '##fileformat=VCFv4.2\n'
ENG
chmod +x "$T/engine"
export STUB_LOG="$T/engine.log"

# ---- dry run: plan only, nothing written, no build -------------------------------------------------
out=$(PROFILE_OUT="$T/prof" PROFILE_BIN="$T/engine" bash $PS baseline --only tiny --dry-run 2>&1); rc=$?
ck "dry-run exits 0" "[ $rc -eq 0 ]"
ck "dry-run names the cell" "grep -q '^PLAN tiny' <<<\"\$out\""
ck "dry-run shows courtesy wrapper" "grep -q 'harness/courtesy.sh' <<<\"\$out\""
ck "dry-run shows all three profiler envs" "grep -q 'VARDICT_PROFILE_CPU=' <<<\"\$out\" && grep -q 'VARDICT_PROFILE_HEAP=' <<<\"\$out\" && grep -q 'VARDICT_PROFILE_DHAT=' <<<\"\$out\""
ck "dry-run uses the stub binary and time -v" "grep -q \"$T/engine\" <<<\"\$out\" && grep -q '/usr/bin/time -v' <<<\"\$out\""
ck "dry-run caps threads at 8" "grep -q -- '-th 8' <<<\"\$out\""
ck "dry-run writes nothing" "[ ! -e '$T/prof' ] && [ ! -e '$T/engine.log' ]"
# dhat is slice-only: the WGS and HG002 plans have no dhat modality
out=$(PROFILE_OUT="$T/prof" PROFILE_BIN="$T/engine" bash $PS optimized --only wgs-scale --dry-run 2>&1)
ck "wgs-scale plan has no dhat" "! grep -q 'VARDICT_PROFILE_DHAT' <<<\"\$out\" && grep -q -- '-R chr1:1-10000000' <<<\"\$out\" && grep -q 'systemd-run' <<<\"\$out\""
out=$(PROFILE_OUT="$T/prof" PROFILE_BIN="$T/engine" bash $PS optimized --only rna-i --dry-run 2>&1)
ck "rna-i plan is the extended m1 command" "grep -q -- '-R chr8:90534611-90620070 -f 0.1 -i' <<<\"\$out\""
out=$(PROFILE_OUT="$T/prof" PROFILE_BIN="$T/engine" bash $PS optimized --only s3 --dry-run 2>&1)
ck "s3 plan passes tumor|normal as one -b" "grep -q 'Colo829-18_S3-sort.bam.|.*Colo829-19_S4-sort.bam' <<<\"\$out\""

# default build trees must be git-ignored (a `git add -A` after the first real run must never stage them)
out=$(env -u PROFILE_BIN -u PROFILE_BIN_CPU -u PROFILE_BIN_HEAP -u PROFILE_BIN_DHAT -u CARGO_TARGET_DIR -u PROFILE_TARGET_BASE \
  PROFILE_OUT="$T/prof" bash $PS baseline --only tiny --dry-run 2>&1)
n=0
while read -r b; do
  n=$((n + 1)); ck "default profile binary dir is git-ignored: $b" "git check-ignore -q '$b'"
done < <(grep -o '[^ ]*/profiling/vardictrust' <<<"$out" | sort -u)
ck "three default profile binaries resolved" "[ $n -eq 3 ]"
ck "legacy /target-*/ dirs are git-ignored too" "git check-ignore -q target-cpu/profiling/vardictrust"

# ---- argument / safety errors --------------------------------------------------------------------------
out=$(bash $PS baseline --only nope 2>&1); rc=$?
ck "unknown cell errors clearly" "[ $rc -ne 0 ] && grep -q \"unknown cell 'nope'\" <<<\"\$out\""
bash $PS >/dev/null 2>&1; ck "missing tag errors" "[ $? -ne 0 ]"
bash $PS bogus >/dev/null 2>&1; ck "bad tag errors" "[ $? -ne 0 ]"
mkdir -p "$T/alias" && ln -s "$PWD/benchmarks/evidence" "$T/alias/ev"
for o in benchmarks/evidence "$PWD/benchmarks/evidence/x" "$T/alias/ev"; do
  out=$(PROFILE_OUT="$o" PROFILE_BIN="$T/engine" bash $PS baseline --only tiny --dry-run 2>&1); rc=$?
  ck "refuses evidence dir $o" "[ $rc -ne 0 ] && grep -q 'benchmarks/evidence' <<<\"\$out\""
done
PROFILE_THREADS=64 PROFILE_BIN="$T/engine" bash $PS baseline --only tiny --dry-run >/dev/null 2>&1
ck "thread cap >8 refused" "[ $? -ne 0 ]"
out=$(PROFILE_OUT="$T/prof" PROFILE_BIN_CPU="$T/nonexistent" PROFILE_BIN_HEAP="$T/engine" PROFILE_BIN_DHAT="$T/engine" \
  bash $PS baseline --only tiny 2>&1); rc=$?
ck "missing binary refused outside dry-run" "[ $rc -ne 0 ] && grep -q 'missing/not executable' <<<\"\$out\" && [ ! -e '$T/prof' ]"

# ---- skip semantics: a mirror checkout with no data ------------------------------------------------------
M="$T/mirror"; mkdir -p "$M/harness" "$M/benchmarks/evidence"
cp harness/profile.sh harness/courtesy.sh harness/profile_catalog.py harness/slices.tsv "$M/harness/"
out=$(cd "$M" && MIN_FREE_GB=0 MIN_AVAIL_GB=0 COURTESY_NOW=0600 PROFILE_BIN="$T/engine" bash harness/profile.sh baseline --only hg002 2>&1); rc=$?
ck "missing input is an error when named via --only" "[ $rc -ne 0 ] && grep -q 'missing input' <<<\"\$out\""
out=$(cd "$M" && MIN_FREE_GB=0 MIN_AVAIL_GB=0 COURTESY_NOW=0600 PROFILE_BIN="$T/engine" bash harness/profile.sh baseline 2>&1); rc=$?
ck "missing inputs SKIP (not fail) on a full run" "[ $rc -eq 0 ] && grep -q '^SKIP hg002: missing input' <<<\"\$out\" && grep -q '^SKIP wgs-scale' <<<\"\$out\""

# ---- real execution through courtesy.sh with the stub engine (tiny slice inputs live in the repo) -----------
if [ -f VarDictJava/testdata/fastas/GRCh37.fa ] && [ -f harness/slices/tiny.bed ]; then
  mkdir -p "$T/bin"
  cat > "$T/bin/systemd-run" <<'SD'
#!/usr/bin/env bash
echo "$@" >> "$SD_LOG"
while [ $# -gt 0 ]; do case "$1" in --user|--scope) shift ;; -p) shift 2 ;; *) break ;; esac; done
exec "$@"
SD
  chmod +x "$T/bin/systemd-run"; export SD_LOG="$T/sd.log"
  CE="env MIN_FREE_GB=0 MIN_AVAIL_GB=0 COURTESY_NOW=0600"
  out=$($CE PATH="$T/bin:$PATH" PROFILE_OUT="$T/prof" PROFILE_BIN="$T/engine" bash $PS baseline --only tiny 2>&1); rc=$?
  d="$T/prof/baseline/tiny"
  ck "stub run exits 0" "[ $rc -eq 0 ]"
  ck "three runs through the stub with -th" "[ \$(wc -l < '$T/engine.log') -eq 3 ] && grep -q -- '-th 8' '$T/engine.log'"
  ck "time -v recorded for every run" "[ \$(grep -c 'Exit status: 0' '$d/time.txt') -eq 3 ]"
  ck "cmd.txt records the commands" "grep -q 'modality: dhat' '$d/cmd.txt' && grep -q courtesy.sh '$d/cmd.txt'"
  ck "dhat run is inside the MemoryMax guard" "grep -q 'MemoryMax=14G' '$T/sd.log' && grep -q 'MemorySwapMax=0' '$T/sd.log'"
  ck "slice cpu/heap runs are not guarded" "[ \$(wc -l < '$T/sd.log') -eq 2 ]"   # probe + the dhat run
  ck "profile outputs land in the cell dir" "[ -s '$d/cpu-collapsed.txt' ] && [ -s '$d/dhat-heap.json' ]"
  ck "catalog.md built from them" "grep -q 'vardictrust_core::stub::hot' '$T/prof/baseline/catalog.md'"
  # systemd-run --user unavailable: still runs, says so
  printf '#!/usr/bin/env bash\nexit 1\n' > "$T/bin/systemd-run"
  out=$($CE PATH="$T/bin:$PATH" PROFILE_OUT="$T/prof2" PROFILE_BIN="$T/engine" bash $PS optimized --only tiny 2>&1); rc=$?
  ck "no systemd-run: warns and proceeds" "[ $rc -eq 0 ] && grep -q 'WITHOUT the MemoryMax' <<<\"\$out\" && [ -s '$T/prof2/optimized/tiny/dhat-heap.json' ]"
  # courtesy refusal is final, never bypassed
  out=$(env MIN_FREE_GB=0 MIN_AVAIL_GB=0 COURTESY_NOW=0400 PROFILE_OUT="$T/prof3" PROFILE_BIN="$T/engine" bash $PS baseline --only tiny 2>&1); rc=$?
  ck "courtesy maintenance-window refusal aborts" "[ $rc -ne 0 ] && grep -q 'maintenance' <<<\"\$out\" && [ ! -e '$T/prof3' ]"
  # a failing engine fails the run
  printf '#!/usr/bin/env bash\nexit 7\n' > "$T/bad"; chmod +x "$T/bad"
  out=$($CE PATH="$T/bin:$PATH" PROFILE_OUT="$T/prof4" PROFILE_BIN="$T/bad" bash $PS baseline --only tiny 2>&1); rc=$?
  ck "engine failure propagates" "[ $rc -ne 0 ] && grep -q 'FAILED tiny' <<<\"\$out\""
else
  echo "test_profile_sh: tiny slice inputs absent; real-execution section skipped" >&2
fi

[ "$fail" -eq 0 ] && echo "test_profile_sh: OK"
exit "$fail"
