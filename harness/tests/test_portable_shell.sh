#!/usr/bin/env bash
# harness/tests/test_portable_shell.sh
# Tests every helper in harness/lib/portable.sh under a minimal `env -i PATH=/usr/bin:/bin` environment.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
PORTABLE_SH="$ROOT_DIR/harness/lib/portable.sh"

[ -f "$PORTABLE_SH" ] || { echo "FAIL test_portable_shell: $PORTABLE_SH does not exist" >&2; exit 1; }
bash -n "$PORTABLE_SH" || { echo "FAIL test_portable_shell: bash -n failed on $PORTABLE_SH" >&2; exit 1; }

TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

printf 'hello world\n' > "$TMP_DIR/sample.txt"

env -i PATH="/usr/bin:/bin" HOME="$TMP_DIR" TMPDIR="$TMP_DIR" bash -c '
  set -euo pipefail
  source "$1"
  tmp="$2"

  # 1. now_ms
  t1="$(now_ms)"
  [[ "$t1" =~ ^[0-9]{10,}$ ]] || { echo "FAIL: now_ms non-numeric: $t1" >&2; exit 1; }
  t2="$(now_ms)"
  [[ "$t2" =~ ^[0-9]{10,}$ ]] && [ "$t2" -ge "$t1" ] || { echo "FAIL: now_ms monotonic check failed: $t1 -> $t2" >&2; exit 1; }

  # 2. mem_avail_gb
  mem="$(mem_avail_gb)"
  [[ "$mem" =~ ^[0-9]+$ ]] || { echo "FAIL: mem_avail_gb non-numeric: $mem" >&2; exit 1; }

  # 3. disk_free_gb & disk_total_gb
  dfree="$(disk_free_gb "$tmp")"
  [[ "$dfree" =~ ^[0-9]+$ ]] && [ "$dfree" -ge 1 ] || { echo "FAIL: disk_free_gb invalid: $dfree" >&2; exit 1; }
  dtotal="$(disk_total_gb "$tmp")"
  [[ "$dtotal" =~ ^[0-9]+$ ]] && [ "$dtotal" -ge "$dfree" ] || { echo "FAIL: disk_total_gb invalid: $dtotal < $dfree" >&2; exit 1; }

  # 4. portable_md5 / md5
  got_md5="$(portable_md5 "$tmp/sample.txt")"
  [ "$got_md5" = "6f5902ac237024bdd0c176cb93063dc4" ] || { echo "FAIL: portable_md5 got $got_md5" >&2; exit 1; }
  [ "$(md5 "$tmp/sample.txt")" = "$got_md5" ] || { echo "FAIL: md5 alias mismatch" >&2; exit 1; }

  # 5. portable_sha256 / sha256 (file and stdin)
  want_sha="a948904f2f0f479b8f8197694b30184b0d2ed1c1cd2a1ec0fb85d299a192a447"
  got_sha="$(portable_sha256 "$tmp/sample.txt")"
  [ "$got_sha" = "$want_sha" ] || { echo "FAIL: portable_sha256(file) got $got_sha" >&2; exit 1; }
  got_sha_stdin="$(portable_sha256 < "$tmp/sample.txt")"
  [ "$got_sha_stdin" = "$want_sha" ] || { echo "FAIL: portable_sha256(stdin) got $got_sha_stdin" >&2; exit 1; }
  [ "$(sha256 "$tmp/sample.txt")" = "$want_sha" ] || { echo "FAIL: sha256 alias mismatch" >&2; exit 1; }

  # 6. portable_filesize / filesize
  fsize="$(portable_filesize "$tmp/sample.txt")"
  [ "$fsize" = "12" ] || { echo "FAIL: portable_filesize expected 12, got $fsize" >&2; exit 1; }
  [ "$(filesize "$tmp/sample.txt")" = "12" ] || { echo "FAIL: filesize alias mismatch" >&2; exit 1; }

  # 7. nice_idle
  out_nice="$(nice_idle echo ok)"
  [ "$out_nice" = "ok" ] || { echo "FAIL: nice_idle output: $out_nice" >&2; exit 1; }

  # 8. portable_time_v and TIME_V_CMD
  portable_time_v -o "$tmp/t1.time" true
  grep -q "Elapsed (wall clock) time" "$tmp/t1.time" || { echo "FAIL: t1.time missing wall clock" >&2; exit 1; }
  grep -q "Maximum resident set size" "$tmp/t1.time" || { echo "FAIL: t1.time missing max rss" >&2; exit 1; }
  grep -q "Exit status: 0" "$tmp/t1.time" || { echo "FAIL: t1.time missing exit status" >&2; exit 1; }

  "${TIME_V_CMD[@]}" -a -o "$tmp/t1.time" true
  [ "$(grep -c "Exit status: 0" "$tmp/t1.time")" -eq 2 ] || { echo "FAIL: TIME_V_CMD -a append did not record second run" >&2; exit 1; }

  "$1" --time-v -o "$tmp/t2.time" true
  grep -q "Exit status: 0" "$tmp/t2.time" || { echo "FAIL: portable.sh --time-v entrypoint failed" >&2; exit 1; }

  # 9. CONDA_PREFIX named environment root resolution
  conda_sub_g4="$(CONDA_PREFIX="/opt/conda/envs/vardict-oracle" G4_TOOLS_BIN="" bash -c "source \"\$1\"; printf \"%s\" \"\$G4_TOOLS_BIN\"" _ "$1")"
  [ "$conda_sub_g4" = "/opt/conda/envs/g4-tools/bin" ] || { echo "FAIL: named CONDA_PREFIX resolved to $conda_sub_g4" >&2; exit 1; }
  conda_explicit_g4="$(CONDA_PREFIX="/opt/conda/envs/vardict-oracle" G4_TOOLS_BIN="/custom/g4/bin" bash -c "source \"\$1\"; printf \"%s\" \"\$G4_TOOLS_BIN\"" _ "$1")"
  [ "$conda_explicit_g4" = "/custom/g4/bin" ] || { echo "FAIL: explicit G4_TOOLS_BIN resolved to $conda_explicit_g4" >&2; exit 1; }

  # 10. Checkout path with spaces + fallback TIME_V_CMD + stderr preservation
  spaced_dir="$tmp/path with spaces"
  mkdir -p "$spaced_dir"
  cp "$1" "$spaced_dir/portable.sh"
  chmod +x "$spaced_dir/portable.sh"
  stderr_out="$("$spaced_dir/portable.sh" --time-v -o "$tmp/t_space.time" sh -c "echo child_stderr_preserved >&2" 2>&1)"
  [ "$stderr_out" = "child_stderr_preserved" ] || { echo "FAIL: child stderr lost or polluted: $stderr_out" >&2; exit 1; }
  grep -q "Exit status: 0" "$tmp/t_space.time" || { echo "FAIL: spaced portable.sh --time-v did not write t_space.time" >&2; exit 1; }
' _ "$PORTABLE_SH" "$TMP_DIR"

echo "PASS test_portable_shell"
