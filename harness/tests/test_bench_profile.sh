#!/usr/bin/env bash
# harness/tests/test_bench_profile.sh
# Verifies that .cargo/config.toml configures [profile.bench] with debuginfo and unstripped
# symbols without polluting release profiles or engine_fingerprint.sh.
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT_DIR"

CFG=".cargo/config.toml"
if [ ! -f "$CFG" ]; then
  echo "FAIL test_bench_profile: missing $CFG" >&2
  exit 1
fi

# 1. Assert .cargo/config.toml only defines [profile.bench*] sections and no rustflags/release overrides.
python3 - "$CFG" <<'PY'
import re, sys

path = sys.argv[1]
with open(path, "r", encoding="utf-8") as fh:
    lines = [line.split("#", 1)[0].strip() for line in fh]

sections = []
kv = {}
cur = None
for line in lines:
    if not line:
        continue
    m = re.match(r"^\[([^\]]+)\]$", line)
    if m:
        cur = m.group(1).strip()
        sections.append(cur)
        kv.setdefault(cur, {})
        continue
    if "=" in line:
        if cur is None:
            sys.exit(f"FAIL: top-level key outside section in {path}: {line}")
        k, v = [x.strip().strip('"').strip("'") for x in line.split("=", 1)]
        kv[cur][k] = v

if not sections:
    sys.exit(f"FAIL: no sections found in {path}")
for sec in sections:
    if not (sec == "profile.bench" or sec.startswith("profile.bench.")):
        sys.exit(f"FAIL: forbidden section [{sec}] in {path}; only [profile.bench*] allowed")

bench = kv.get("profile.bench", {})
if bench.get("inherits") != "release":
    sys.exit(f"FAIL: [profile.bench] inherits={bench.get('inherits')!r}, expected 'release'")
if bench.get("strip") != "false":
    sys.exit(f"FAIL: [profile.bench] strip={bench.get('strip')!r}, expected 'false'")
if bench.get("debug") not in ("1", "2", "true", "line-tables-only"):
    sys.exit(f"FAIL: [profile.bench] debug={bench.get('debug')!r}, expected debuginfo enabled")
PY

# 2. Verify cargo bench --no-run produces an unstripped bench binary with debuginfo != 0
TMP_JSON=$(mktemp)
trap 'rm -f "$TMP_JSON"' EXIT
cargo bench --no-run --message-format=json -p vardictrust-core > "$TMP_JSON"
BENCH_EXE=$(python3 - "$TMP_JSON" <<'PY'
import json, sys

exe = None
with open(sys.argv[1], "r", encoding="utf-8") as fh:
    for line in fh:
        if not line.startswith("{"):
            continue
        msg = json.loads(line)
        if msg.get("reason") != "compiler-artifact":
            continue
        target = msg.get("target") or {}
        if target.get("name") == "engine" and "bench" in (target.get("kind") or []):
            prof = msg.get("profile") or {}
            dbg = prof.get("debuginfo")
            if not dbg:
                sys.exit(f"FAIL: engine bench profile.debuginfo={dbg!r}, expected non-zero")
            exe = msg.get("executable")
if not exe:
    sys.exit("FAIL: could not locate engine bench executable in cargo bench --message-format=json")
print(exe)
PY
)

if [ ! -x "$BENCH_EXE" ]; then
  echo "FAIL test_bench_profile: bench executable not found at $BENCH_EXE" >&2
  exit 1
fi

SYM_COUNT=$(nm "$BENCH_EXE" | grep -c "vardictrust_core" || true)
if [ "$SYM_COUNT" -le 0 ]; then
  echo "FAIL test_bench_profile: expected vardictrust_core symbols in $BENCH_EXE, found 0" >&2
  exit 1
fi

echo "PASS test_bench_profile"
