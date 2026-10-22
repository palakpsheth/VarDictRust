#!/usr/bin/env bash
# harness/tests/test_scale_bench.sh — scale_bench.sh measurement protocol (issue #31). Stub engine, no data.
set -uo pipefail
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT; fail=0
SB=harness/scale_bench.sh
printf '>chr1\nACGT\n' > "$T/ref.fa"; printf 'chr1\t120000000\t6\t60\t61\n' > "$T/ref.fa.fai"
: > "$T/x.bam"; printf 'chr1\t1\t100\tg\n' > "$T/x.bed"
mkdir -p "$T/stub"; export STUB_DIR="$T/stub"
cat > "$T/engine" <<'ENG'
#!/usr/bin/env bash
th=1; for ((i=1;i<=$#;i++)); do if [ "${!i}" = "-th" ]; then j=$((i+1)); th=${!j}; fi; done
f="$STUB_DIR/n_$th"; n=$(cat "$f" 2>/dev/null || echo 0); echo $((n+1)) > "$f"
if [ "$th" = 1 ]; then case $n in 0) sleep 0.2;; 1) sleep 0.8;; *) sleep 0.5;; esac; else sleep 0.05; fi
if [ -n "${STUB_FAIL_TH:-}" ] && [ "$th" = "$STUB_FAIL_TH" ]; then echo boom >&2; exit 9; fi
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
ENG
chmod +x "$T/engine"
CW="env MIN_FREE_GB=0 MIN_AVAIL_GB=0 COURTESY_NOW=0600 EST_MIN=1"
# run [args...]: through the real courtesy.sh (nice 19 + COURTESY_ACTIVE=1), 1 cpu pinned unless the caller overrides
run() { $CW bash harness/courtesy.sh env VARDICT_BIN="$T/engine" bash $SB "$@" "$T/ref.fa" "$T/x.bam" "$T/x.bed" S "$T/out" 2>"$T/err"; }
bare() { env VARDICT_BIN="$T/engine" bash $SB "$@" "$T/ref.fa" "$T/x.bam" "$T/x.bed" S "$T/out" 2>"$T/err"; }
reset() { rm -rf "$T/out" "$T/stub"/*; }
jget() { python3 -c "import json,sys; d=json.load(open('$T/out/extra.json')); print(eval(sys.argv[1]))" "$1"; }

# 1. courtesy gate: refused without the wrapper; the escape hatch runs but never certifies; a spoofed marker is not enough
reset; if ( unset COURTESY_ACTIVE; SCALE_CPUS=0 bare --no-cpp ); then echo "FAIL1 ran without courtesy"; fail=1; fi
grep -q 'courtesy.sh' "$T/err" || { echo "FAIL1 no courtesy hint"; fail=1; }
reset; ( unset COURTESY_ACTIVE; SCALE_ALLOW_UNCOURTEOUS=1 SCALE_RUNS=1 SCALE_CPUS=0 bare --no-cpp ) || { echo "FAIL1b escape hatch rc"; cat "$T/err"; fail=1; }
[ "$(jget "d['courtesy']")" = False ] && [ "$(jget "d['certifying']")" = False ] || { echo "FAIL1b courtesy/certifying"; fail=1; }
reset; COURTESY_ACTIVE=1 SCALE_ALLOW_UNCOURTEOUS=1 SCALE_RUNS=1 SCALE_CPUS=0 bare --no-cpp || { echo "FAIL1c rc"; cat "$T/err"; fail=1; }
[ "$(jget "d['courtesy']")" = False ] || { echo "FAIL1c a spoofed COURTESY_ACTIVE without nice must not count"; fail=1; }

# 2. median of 3 runs (th1 sleeps 0.2/0.8/0.5s -> median 0.5), run logs, host facts, --no-cpp
reset; SCALE_CPUS=0 run --no-cpp || { echo "FAIL2 rc"; cat "$T/err"; fail=1; }
w=$(jget "d['wall_s']['1']"); python3 -c "import sys; w=float('$w'); sys.exit(0 if 0.4 <= w <= 0.75 else 1)" || { echo "FAIL2 median wall=$w"; fail=1; }
[ "$(jget "d['runs']")" = 3 ] && [ "$(jget "d['courtesy']")" = True ] || { echo "FAIL2 runs/courtesy"; fail=1; }
for k in governor loadavg_start loadavg_end cpus perf_fingerprint; do jget "d['$k']" >/dev/null 2>&1 || { echo "FAIL2 missing $k"; fail=1; }; done
[ -f "$T/out/th1.run3.time" ] && [ ! -e "$T/out/cpp_th1.run1.time" ] || { echo "FAIL2 run logs / cpp skipped"; fail=1; }
grep -q 'CPP: SKIP: disabled via CPP_SKIP=1' "$T/err" || { echo "FAIL2 no CPP SKIP"; fail=1; }

# 3. certification rules. One CPU can never certify (the frozen bar 2.4 is defined at P=4)
[ "$(jget "d['P']")" = 1 ] && [ "$(jget "d['pinned']")" = True ] && [ "$(jget "d['certifying']")" = False ] || { echo "FAIL3 P=1 must not certify"; fail=1; }
[ "$(jget "d['bar']")" = 0.6 ] || { echo "FAIL3 bar"; fail=1; }
jget "d['maxrss_kb']['120000000']" >/dev/null && [ "$(jget "len(d['maxrss_kb'])")" = 4 ] || { echo "FAIL3 span from .fai"; fail=1; }
if [ "$(nproc --all)" -ge 8 ]; then   # default pin 0-7: certifying must equal the documented rule, whatever this host's topology
  reset; SCALE_RUNS=3 run --no-cpp || { echo "FAIL3b rc"; cat "$T/err"; fail=1; }
  [ "$(jget "d['cpus']")" = "0-7" ] && [ "$(jget "d['pinned']")" = True ] || { echo "FAIL3b default pin 0-7"; fail=1; }
  [ "$(jget "d['certifying'] == (d['P'] >= 4 and d['perf_fingerprint'] != 'dirty' and d['courtesy'] and d['runs'] >= 3)")" = True ] || { echo "FAIL3b certifying rule"; fail=1; }
fi
if [ "$(nproc --all)" -ge 2 ]; then
  reset; SCALE_CPUS= SCALE_RUNS=1 run --no-cpp || { echo "FAIL3c rc"; cat "$T/err"; fail=1; }
  [ "$(jget "d['pinned']")" = False ] && [ "$(jget "d['certifying']")" = False ] || { echo "FAIL3c unpinned must not certify"; fail=1; }
fi
reset; SCALE_RUNS=2 SCALE_CPUS=0 run --no-cpp; [ "$(jget "d['certifying']")" = False ] || { echo "FAIL3d <3 runs must not certify"; fail=1; }
reset; SCALE_CPUS=0-6:2 run --no-cpp && { echo "FAIL3e stride accepted"; fail=1; }
grep -q 'plain list/range' "$T/err" || { echo "FAIL3e no message"; fail=1; }

# 4. named cell failure keeps the partial time log, exits non-zero, and leaves no stale verdict
reset; SCALE_CPUS=0 SCALE_RUNS=1 run --no-cpp; [ -f "$T/out/extra.json" ] || { echo "FAIL4 setup"; fail=1; }
STUB_FAIL_TH=4 SCALE_CPUS=0 SCALE_RUNS=1 run --no-cpp && { echo "FAIL4 rc 0"; fail=1; }
grep -q 'CELL FAILED: th4 run 1 (exit 9)' "$T/err" || { echo "FAIL4 message: $(cat "$T/err")"; fail=1; }
[ -f "$T/out/th4.run1.time" ] || { echo "FAIL4 partial log not kept"; fail=1; }
[ ! -e "$T/out/extra.json" ] || { echo "FAIL4 stale extra.json survived a failed run"; fail=1; }

# 5. default outdir is under scratch/ (textual: running it would write into the repo)
grep -q 'OUT=${5:-scratch/g4/scale}' $SB || { echo "FAIL5 default OUT"; fail=1; }

# 6. both cpp arms route through run_cpp.sh (stub conda answers `run -n vardict-oracle vardictcpp ...`)
mkdir -p "$T/bin"
cat > "$T/bin/conda" <<'CND'
#!/usr/bin/env bash
shift 3
[ "${1:-}" = vardictcpp ] || exit 0
shift
if [ "${1:-}" = "--version" ]; then echo "vardictcpp stub"; exit 0; fi
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
CND
chmod +x "$T/bin/conda"
reset; PATH="$T/bin:$PATH" SCALE_RUNS=1 SCALE_CPUS=0 run || { echo "FAIL6 rc"; cat "$T/err"; fail=1; }
[ "$(jget "d['cpp_wall_s']['1'] is not None")" = True ] || { echo "FAIL6 cpp BED arm not recorded"; fail=1; }
[ "$(jget "d['cpp_maxrss_kb']['10000000'] is not None")" = True ] || { echo "FAIL6 cpp region arm not recorded"; fail=1; }
[ "$(cat "$T/out/R10000000.bed")" = "$(printf 'chr1\t0\t10000000\tR10000000')" ] || { echo "FAIL6 region BED"; fail=1; }

# 7. courtesy.sh marks the run so scale_bench can tell it was wrapped
[ "$($CW bash harness/courtesy.sh bash -c 'echo $COURTESY_ACTIVE')" = 1 ] || { echo "FAIL7 courtesy.sh does not export COURTESY_ACTIVE"; fail=1; }
[ $fail -eq 0 ] && echo "PASS test_scale_bench"; exit $fail
