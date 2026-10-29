#!/usr/bin/env bash
# Red-green test for evidence recording and fingerprint staleness check (T06)
set -euo pipefail

fail=0
EVID_DIR="benchmarks/evidence"
mkdir -p "$EVID_DIR"

cleanup() {
  rm -f "$EVID_DIR/t-test-stale.json" "$EVID_DIR/t-test-badres.json" "$EVID_DIR/t-test.json" \
        "$EVID_DIR/t-test-cpp-both.json" "$EVID_DIR/t-test-cpp-empty.json" "$EVID_DIR/t-test-cpp-ok.json" \
        "$EVID_DIR/t-test-reuse.json" "$EVID_DIR/t-test-scale.json"
  [ -n "${T_KEEP:-}" ] && rm -rf "$T_KEEP" || true
}
trap cleanup EXIT

# 0. Test ev_cpp() from acceptance.sh (mutual exclusivity + non-empty string)
eval "$(sed -n '/^ev_cpp()/,/^}/p' harness/acceptance.sh)"

cat > "$EVID_DIR/t-test-cpp-both.json" <<'EOF'
{
  "name": "t-test-cpp-both",
  "result": "MATCH",
  "extra": {
    "cpp_version": "vardictcpp 3 hc52dbad_0 bioconda",
    "cpp_skip_reason": "disabled via CPP_SKIP=1"
  }
}
EOF
if ev_cpp t-test-cpp-both >/dev/null 2>&1; then
  echo "FAIL: ev_cpp accepted record with BOTH cpp_version and cpp_skip_reason" >&2
  fail=1
fi

cat > "$EVID_DIR/t-test-cpp-empty.json" <<'EOF'
{
  "name": "t-test-cpp-empty",
  "result": "MATCH",
  "extra": {
    "cpp_version": "   "
  }
}
EOF
if ev_cpp t-test-cpp-empty >/dev/null 2>&1; then
  echo "FAIL: ev_cpp accepted record with empty/whitespace cpp_version" >&2
  fail=1
fi

cat > "$EVID_DIR/t-test-cpp-ok.json" <<'EOF'
{
  "name": "t-test-cpp-ok",
  "result": "MATCH",
  "extra": {
    "cpp_version": "vardictcpp 3 hc52dbad_0 bioconda"
  }
}
EOF
if ! ev_cpp t-test-cpp-ok >/dev/null 2>&1; then
  echo "FAIL: ev_cpp rejected valid cpp_version record" >&2
  fail=1
fi

# 1. Obtain current real fingerprint
FP=$(bash harness/engine_fingerprint.sh)
if [[ ! "$FP" =~ ^[0-9a-f]{16}$ ]]; then
  echo "FAIL: invalid fingerprint format '$FP'" >&2
  exit 1
fi

# Define the exact ev function from acceptance.sh
ev() {
  python3 - "$EVID_DIR/$1.json" "$FP" <<'PY'
import json, sys
p, fp = sys.argv[1], sys.argv[2]
try: d = json.load(open(p))
except Exception as e: sys.exit(f"{p}: {e}")
if d.get("fingerprint") != fp: sys.exit(f"{p}: stale fingerprint {d.get('fingerprint')} != {fp}")
if d.get("result") not in ("MATCH", "EQUAL", "PASS"): sys.exit(f"{p}: result={d.get('result')}")
PY
}

# 2. Verify stale/mismatched fingerprint 'deadbeef' is rejected
cat > "$EVID_DIR/t-test-stale.json" <<EOF
{
  "name": "t-test-stale",
  "result": "MATCH",
  "fingerprint": "deadbeef"
}
EOF

if err=$(ev t-test-stale 2>&1); then
  echo "FAIL: ev accepted stale fingerprint deadbeef" >&2
  fail=1
else
  if [[ "$err" != *"stale fingerprint deadbeef !="* ]]; then
    echo "FAIL: expected stale fingerprint error message, got: $err" >&2
    fail=1
  fi
fi

# Also verify invalid result is rejected
cat > "$EVID_DIR/t-test-badres.json" <<EOF
{
  "name": "t-test-badres",
  "result": "MISMATCH",
  "fingerprint": "$FP"
}
EOF

if err=$(ev t-test-badres 2>&1); then
  echo "FAIL: ev accepted MISMATCH result" >&2
  fail=1
else
  if [[ "$err" != *"result=MISMATCH"* ]]; then
    echo "FAIL: expected result=MISMATCH error message, got: $err" >&2
    fail=1
  fi
fi

# 3. Verify record_evidence.py writes an evidence JSON containing the real fingerprint
python3 harness/record_evidence.py t-test MATCH - -- true >/dev/null

if [ ! -f "$EVID_DIR/t-test.json" ]; then
  echo "FAIL: $EVID_DIR/t-test.json was not created" >&2
  fail=1
else
  recorded_fp=$(python3 -c 'import json; print(json.load(open("benchmarks/evidence/t-test.json"))["fingerprint"])')
  if [ "$recorded_fp" != "$FP" ]; then
    echo "FAIL: recorded fingerprint $recorded_fp != expected $FP" >&2
    fail=1
  fi
  # Verify ev accepts the newly recorded evidence
  if ! ev t-test; then
    echo "FAIL: ev rejected newly recorded evidence t-test" >&2
    fail=1
  fi
fi

# 4. Verify record_evidence.py ignores stale compare_cpp.txt and cpp.time when extra.json has a pre-execution skip
T_KEEP=$(mktemp -d)
echo "MATCH (5 variants)" > "$T_KEEP/compare.txt"
echo "CPP: MATCH (java==cpp, cpp==rust)" > "$T_KEEP/compare_cpp.txt"
printf '\tElapsed (wall clock) time (h:mm:ss or m:ss): 0:01.00\n\tMaximum resident set size (kbytes): 1000\n\tExit status: 0\n' > "$T_KEEP/cpp.time"
echo '{"cpp_skip_reason": "disabled via CPP_SKIP=1"}' > "$T_KEEP/extra.json"
python3 harness/record_evidence.py t-test-reuse MATCH "$T_KEEP" -- true >/dev/null
if ! python3 -c '
import json
d = json.load(open("benchmarks/evidence/t-test-reuse.json"))
assert d["cpp_time"] is None, d["cpp_time"]
assert d["extra"]["cpp_head"] == ["CPP: SKIP: disabled via CPP_SKIP=1"], d["extra"]
'; then
  echo "FAIL: record_evidence.py leaked stale compare_cpp.txt or cpp.time on KEEP_DIR reuse" >&2
  fail=1
fi

# 5. Verify engine_fingerprint.sh hashes raw HEAD git blob OIDs directly (no normalize_gate_blob workaround)
if grep -q 'normalize_gate_blob' harness/engine_fingerprint.sh; then
  echo "FAIL: harness/engine_fingerprint.sh still contains obsolete normalize_gate_blob()" >&2
  fail=1
fi
fp_default=$(bash harness/engine_fingerprint.sh)
fp_perf=$(bash harness/engine_fingerprint.sh --perf)
if ! [[ "$fp_default" =~ ^[0-9a-f]{16}$ && "$fp_perf" =~ ^[0-9a-f]{16}$ ]]; then
  echo "FAIL: engine_fingerprint.sh produced invalid digest(s): default='$fp_default' perf='$fp_perf'" >&2
  fail=1
fi

# 6. ev_scale (acceptance.sh item 5): the scale bench is carried only while the perf-relevant code is unchanged
#    (extra.perf_fingerprint); there is no waiver. P>=4 / bar>=2.4 is the pre-registered bar.
eval "$(sed -n '/^ev_scale()/,/^}/p' harness/acceptance.sh)"
PERF_FP=$(bash harness/engine_fingerprint.sh --perf 2>/dev/null || echo 1298c1d1091a0fc8)
scale_case() {  # <ok|fail> <result> <fingerprint> <extra-json> <message>
  local want=$1 res=$2 fpv=$3 extra=$4 msg=$5 rc=0
  printf '{"name":"t-test-scale","result":"%s","fingerprint":"%s","extra":%s}\n' "$res" "$fpv" "$extra" > "$EVID_DIR/t-test-scale.json"
  ev_scale t-test-scale >/dev/null 2>&1 || rc=$?
  if [ "$want" = ok ] && [ "$rc" -ne 0 ]; then echo "FAIL: ev_scale rejected: $msg" >&2; fail=1; fi
  if [ "$want" = fail ] && [ "$rc" -eq 0 ]; then echo "FAIL: ev_scale accepted: $msg" >&2; fail=1; fi
}
GOOD='"pass":true,"certifying":true,"P":4,"bar":2.4'
scale_case ok   PASS     "$FP"  "{$GOOD,\"perf_fingerprint\":\"other\"}"        "matching engine fingerprint"
scale_case ok   PASS     stale  "{$GOOD,\"perf_fingerprint\":\"$PERF_FP\"}"     "stale engine fingerprint but unchanged perf paths (the carry rule)"
scale_case fail PASS     stale  "{$GOOD,\"perf_fingerprint\":\"changed\"}"      "perf-relevant code changed since the run"
scale_case fail PASS     stale  "{$GOOD}"                                       "evidence without a perf_fingerprint"
scale_case fail PASS     "$FP"  '{"pass":true,"certifying":true,"P":1,"bar":0.6}'   "a one-cpu bar of 0.6 is not the pre-registered bar"
scale_case fail PASS     "$FP"  '{"pass":false,"certifying":true,"P":4,"bar":2.4}'  "extra.pass is false"
scale_case fail PASS     "$FP"  '{"pass":true,"certifying":false,"P":4,"bar":2.4}'  "non-certifying run"
scale_case fail MISMATCH "$FP"  "{$GOOD}"                                       "non-PASS result"

if [ $fail -eq 0 ]; then
  echo "PASS test_evidence"
  exit 0
else
  echo "FAIL test_evidence"
  exit 1
fi

