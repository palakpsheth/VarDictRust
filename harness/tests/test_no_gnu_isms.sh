#!/usr/bin/env bash
# harness/tests/test_no_gnu_isms.sh
# Static scan of harness/*.sh that fails on raw GNU/Linux-only constructs
# (/usr/bin/time -v, date +%s%3N, date +%s%N, /proc/, mktemp --suffix, stat -c,
# df -BG, ionice, md5sum, sha256sum) and hardcoded /home/ paths outside harness/lib/portable.sh
# (excluding compare.sh and run_oracle.sh which are owned by Issue B / Task 3).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$ROOT_DIR"

fail=0

for f in harness/*.sh; do
  base="$(basename "$f")"
  # compare.sh and run_oracle.sh belong to Task 3 (Issue B)
  if [[ "$base" == "compare.sh" || "$base" == "run_oracle.sh" ]]; then
    continue
  fi
  # Strip comment-only lines before scanning
  non_comment="$(grep -n -v '^[[:space:]]*#' "$f" || true)"
  for pat in \
    '/usr/bin/time -v' \
    'date +%s%3N' \
    'date +%s%N' \
    '/proc/' \
    'mktemp --suffix' \
    'stat -c' \
    'df -BG' \
    'ionice' \
    'md5sum' \
    'sha256sum' \
    '/home/'
  do
    if matches="$(printf '%s\n' "$non_comment" | grep -F "$pat" || true)" && [[ -n "$matches" ]]; then
      echo "FAIL test_no_gnu_isms: $f contains forbidden construct '$pat':" >&2
      printf '%s\n' "$matches" | sed 's/^/  /' >&2
      fail=1
    fi
  done
done

# Ensure zero hardcoded /home/palak references anywhere in harness/ or .githooks/
palak_hits="$(grep -rn '/home/palak' harness/ .githooks/ 2>/dev/null | grep -v 'harness/tests/test_no_gnu_isms.sh' || true)"
if [[ -n "$palak_hits" ]]; then
  echo "FAIL test_no_gnu_isms: hardcoded /home/palak found:" >&2
  printf '%s\n' "$palak_hits" | sed 's/^/  /' >&2
  fail=1
fi

if [ "$fail" -eq 0 ]; then
  echo "PASS test_no_gnu_isms"
  exit 0
else
  exit 1
fi
