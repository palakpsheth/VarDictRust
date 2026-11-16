#!/usr/bin/env bash
# harness/tests/test_no_gnu_isms.sh
# Static scan of harness/ and packaging/ shell scripts that fails on raw GNU/Linux-only
# or bash-4-only constructs (/usr/bin/time -v, date +%s%3N, date +%s%N, /proc/,
# mktemp --suffix, stat -c, df -BG, ionice, md5sum, sha256sum, sed -i, realpath -m,
# mapfile, zcat, raw timeout, ${var^^}, unguarded free -h) and hardcoded /home/ paths
# outside harness/lib/portable.sh.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$ROOT_DIR"

fail=0

SCAN_FILES=()
while IFS= read -r f; do
  [ -n "$f" ] && SCAN_FILES+=("$f")
done < <(find harness packaging -type f -name "*.sh" ! -path "harness/tests/*" ! -path "harness/lib/portable.sh" | LC_ALL=C sort)

for f in "${SCAN_FILES[@]}"; do
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
    'sed -i' \
    'realpath -m' \
    'mapfile' \
    'zcat' \
    '/home/'
  do
    if matches="$(printf '%s\n' "$non_comment" | grep -F "$pat" || true)" && [[ -n "$matches" ]]; then
      echo "FAIL test_no_gnu_isms: $f contains forbidden construct '$pat':" >&2
      printf '%s\n' "$matches" | sed 's/^/  /' >&2
      fail=1
    fi
  done

  if matches="$(printf '%s\n' "$non_comment" | grep -E '\$\{[A-Za-z0-9_]+\^\^\}' || true)" && [[ -n "$matches" ]]; then
    echo "FAIL test_no_gnu_isms: $f contains bash-4 uppercase expansion '\${var^^}':" >&2
    printf '%s\n' "$matches" | sed 's/^/  /' >&2
    fail=1
  fi

  if matches="$(printf '%s\n' "$non_comment" | grep -E '([0-9]+:|;|&&|\|\||\{|\(|[[:space:]])timeout[[:space:]]+(\$|"|\\\$|\\"|'\''|[0-9])' || true)" && [[ -n "$matches" ]]; then
    echo "FAIL test_no_gnu_isms: $f contains raw 'timeout' invocation (use TIMEOUT_CMD / portable.sh --timeout):" >&2
    printf '%s\n' "$matches" | sed 's/^/  /' >&2
    fail=1
  fi

  if matches="$(printf '%s\n' "$non_comment" | grep -F 'free -h' | grep -v 'command -v free' | grep -v 'echo ' || true)" && [[ -n "$matches" ]]; then
    echo "FAIL test_no_gnu_isms: $f contains unguarded 'free -h' (guard with command -v free):" >&2
    printf '%s\n' "$matches" | sed 's/^/  /' >&2
    fail=1
  fi
done

# Ensure zero hardcoded /home/palak references anywhere in harness/, packaging/, or .githooks/
palak_hits="$(grep -rn '/home/palak' harness/ packaging/ .githooks/ 2>/dev/null | grep -v 'harness/tests/test_no_gnu_isms.sh' || true)"
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
