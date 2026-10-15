#!/usr/bin/env bash
# Red-green test for evidence recording and fingerprint staleness check (T06)
set -euo pipefail

fail=0
EVID_DIR="benchmarks/evidence"
mkdir -p "$EVID_DIR"

cleanup() {
  rm -f "$EVID_DIR/t-test-stale.json" "$EVID_DIR/t-test-badres.json" "$EVID_DIR/t-test.json" \
        "$EVID_DIR/t-test-cpp-both.json" "$EVID_DIR/t-test-cpp-empty.json" "$EVID_DIR/t-test-cpp-ok.json" \
        "$EVID_DIR/t-test-reuse.json"
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

# 5. Verify normalize_gate_blob() is invariant to workspace SemVer bumps and sensitive to dependency/profile mutations
eval "$(sed -n '/^normalize_gate_blob()/,/^}/p' harness/engine_fingerprint.sh)"
if ! declare -F normalize_gate_blob >/dev/null 2>&1; then
  echo "FAIL: normalize_gate_blob() not defined in harness/engine_fingerprint.sh" >&2
  fail=1
else
  for mf in Cargo.toml crates/vardictrust/Cargo.toml Cargo.lock; do
    base_oid=$(git show "HEAD:$mf" | normalize_gate_blob "$mf" | git hash-object -t blob --stdin)
    v100_raw_oid=$(git show "HEAD:$mf" | python3 -c 'import re,sys; s=sys.stdin.read(); s=re.sub(r"(\[workspace\.package\]\r?\nversion\s*=\s*\")[^\"]+(\")", r"\g<1>1.0.0\2", s, count=1); s=re.sub(r"(vardictrust-core\s*=\s*\{[^\n}]*?\bversion\s*=\s*\")[^\"]+(\")", r"\g<1>1.0.0\2", s, count=1); s=re.sub(r"(\[\[package\]\]\r?\nname = \"vardictrust(?:-core)?\"\r?\nversion = \")[^\"]+(\")", r"\g<1>1.0.0\2", s, count=2); sys.stdout.write(s)' | git hash-object -t blob --stdin)
    if [ "$base_oid" != "$v100_raw_oid" ]; then
      echo "FAIL: normalize_gate_blob($mf) did not match canonical 1.0.0 OID ($base_oid != $v100_raw_oid)" >&2
      fail=1
    fi
    bumped_oid=$(git show "HEAD:$mf" | python3 -c 'import re,sys; s=sys.stdin.read(); s=re.sub(r"(\[workspace\.package\]\r?\nversion\s*=\s*\")[^\"]+(\")", r"\g<1>9.9.9-rc.1\2", s, count=1); s=re.sub(r"(vardictrust-core\s*=\s*\{[^\n}]*?\bversion\s*=\s*\")[^\"]+(\")", r"\g<1>9.9.9-rc.1\2", s, count=1); s=re.sub(r"(\[\[package\]\]\r?\nname = \"vardictrust(?:-core)?\"\r?\nversion = \")[^\"]+(\")", r"\g<1>9.9.9-rc.1\2", s, count=2); sys.stdout.write(s)' | normalize_gate_blob "$mf" | git hash-object -t blob --stdin)
    if [ "$base_oid" != "$bumped_oid" ]; then
      echo "FAIL: normalize_gate_blob($mf) changed OID on workspace SemVer bump ($base_oid != $bumped_oid)" >&2
      fail=1
    fi
  done
  dep_oid=$(git show "HEAD:Cargo.toml" | sed 's/noodles = { version = "0.117"/noodles = { version = "0.118"/' | normalize_gate_blob Cargo.toml | git hash-object -t blob --stdin)
  base_cargo_oid=$(git show "HEAD:Cargo.toml" | normalize_gate_blob Cargo.toml | git hash-object -t blob --stdin)
  if [ "$dep_oid" = "$base_cargo_oid" ]; then
    echo "FAIL: normalize_gate_blob(Cargo.toml) did not change OID when dependency version changed" >&2
    fail=1
  fi
  prof_oid=$(git show "HEAD:Cargo.toml" | sed 's/opt-level = 3/opt-level = 2/' | normalize_gate_blob Cargo.toml | git hash-object -t blob --stdin)
  if [ "$prof_oid" = "$base_cargo_oid" ]; then
    echo "FAIL: normalize_gate_blob(Cargo.toml) did not change OID when release profile changed" >&2
    fail=1
  fi
  lock_dep_oid=$(git show "HEAD:Cargo.lock" | sed -E 's/"1\.0\.[2-9][0-9]*"/"1.0.999"/' | normalize_gate_blob Cargo.lock | git hash-object -t blob --stdin)
  base_lock_oid=$(git show "HEAD:Cargo.lock" | normalize_gate_blob Cargo.lock | git hash-object -t blob --stdin)
  if [ "$lock_dep_oid" = "$base_lock_oid" ]; then
    echo "FAIL: normalize_gate_blob(Cargo.lock) did not change OID when external 1.0.x crate version changed" >&2
    fail=1
  fi
fi

if [ $fail -eq 0 ]; then
  echo "PASS test_evidence"
  exit 0
else
  echo "FAIL test_evidence"
  exit 1
fi

