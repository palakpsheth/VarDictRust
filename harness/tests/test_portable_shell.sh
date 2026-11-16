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

  # 11. Deterministic Python timing backend: stderr preservation, piped stdin, failing exit status, Darwin RSS conversion
  py_stdin_out="$(printf "piped_payload" | PORTABLE_TIME_BACKEND=python portable_time_v -o "$tmp/t_py.time" sh -c "cat; echo py_stderr_msg >&2" 2>"$tmp/t_py.stderr")"
  [ "$py_stdin_out" = "piped_payload" ] || { echo "FAIL: python backend lost piped stdin: $py_stdin_out" >&2; exit 1; }
  [ "$(cat "$tmp/t_py.stderr")" = "py_stderr_msg" ] || { echo "FAIL: python backend stderr mismatch" >&2; exit 1; }
  grep -q "Exit status: 0" "$tmp/t_py.time" || { echo "FAIL: python backend missing Exit status: 0" >&2; exit 1; }
  py_rc=0
  PORTABLE_TIME_BACKEND=python portable_time_v -o "$tmp/t_py_fail.time" sh -c "exit 17" || py_rc=$?
  [ "$py_rc" -eq 17 ] || { echo "FAIL: python backend did not propagate exit 17 (got $py_rc)" >&2; exit 1; }
  grep -q "Exit status: 17" "$tmp/t_py_fail.time" || { echo "FAIL: python backend did not record Exit status: 17" >&2; exit 1; }

  # 12. Deterministic BSD timing backend (mock time -l): stderr preservation, piped stdin, failing exit status, Darwin bytes->kbytes RSS
  mock_bsd_time="$tmp/mock_bsd_time"
  cat > "$mock_bsd_time" <<'\''EOF'\''
#!/bin/sh
if [ "${1:-}" = "-l" ]; then
  shift
  "$@"
  rc=$?
  printf "        0.01 real         0.00 user         0.00 sys\n  104857600  maximum resident set size\n" >&2
  exit $rc
fi
exit 1
EOF
  chmod +x "$mock_bsd_time"
  bsd_stdin_out="$(printf "bsd_payload" | PORTABLE_TIME_BACKEND=bsd PORTABLE_TIME_BIN="$mock_bsd_time" PORTABLE_UNAME_S=Darwin portable_time_v -o "$tmp/t_bsd.time" sh -c "cat; echo bsd_stderr_msg >&2" 2>"$tmp/t_bsd.stderr")"
  [ "$bsd_stdin_out" = "bsd_payload" ] || { echo "FAIL: bsd backend lost piped stdin: $bsd_stdin_out" >&2; exit 1; }
  [ "$(cat "$tmp/t_bsd.stderr")" = "bsd_stderr_msg" ] || { echo "FAIL: bsd backend stderr mismatch: $(cat "$tmp/t_bsd.stderr")" >&2; exit 1; }
  grep -q "Maximum resident set size (kbytes): 102400" "$tmp/t_bsd.time" || { echo "FAIL: bsd backend Darwin RSS conversion failed: $(cat "$tmp/t_bsd.time")" >&2; exit 1; }
  bsd_rc=0
  PORTABLE_TIME_BACKEND=bsd PORTABLE_TIME_BIN="$mock_bsd_time" portable_time_v -o "$tmp/t_bsd_fail.time" sh -c "exit 23" || bsd_rc=$?
  [ "$bsd_rc" -eq 23 ] || { echo "FAIL: bsd backend did not propagate exit 23 (got $bsd_rc)" >&2; exit 1; }
  grep -q "Exit status: 23" "$tmp/t_bsd_fail.time" || { echo "FAIL: bsd backend did not record Exit status: 23" >&2; exit 1; }

  # 13. exec_nice_idle replaces shell process
  exec_out="$(bash -c "source \"\$1\"; exec_nice_idle echo exec_ok" _ "$1")"
  [ "$exec_out" = "exec_ok" ] || { echo "FAIL: exec_nice_idle output: $exec_out" >&2; exit 1; }

  # 14. portable_timeout / --timeout CLI dispatch + Python fallback + process-group cleanup
  to_rc=0
  "$1" --timeout 0.2s sh -c '\''sleep 30 & echo $! > "$1"; wait'\'' _ "$tmp/to_child.pid" || to_rc=$?
  [ "$to_rc" -eq 124 ] || { echo "FAIL: --timeout 0.2s expected exit 124, got $to_rc" >&2; exit 1; }
  sleep 0.05
  if [ -s "$tmp/to_child.pid" ] && kill -0 "$(cat "$tmp/to_child.pid")" 2>/dev/null; then
    kill -9 "$(cat "$tmp/to_child.pid")" 2>/dev/null || true
    echo "FAIL: --timeout left orphan child process $(cat "$tmp/to_child.pid")" >&2
    exit 1
  fi
  "$1" --timeout 1h true || { echo "FAIL: --timeout 1h true failed" >&2; exit 1; }
  portable_time_v -o "$tmp/t_to.time" "${TIMEOUT_CMD[@]}" 30s true || { echo "FAIL: portable_time_v + TIMEOUT_CMD failed" >&2; exit 1; }
  grep -q "Exit status: 0" "$tmp/t_to.time" || { echo "FAIL: t_to.time missing Exit status: 0" >&2; exit 1; }

  no_to_bin="$tmp/no_timeout_bin"
  mkdir -p "$no_to_bin"
  for util in sh bash python3 sleep true cat kill rm mkdir dirname basename pwd uname tr awk cut wc date nice; do
    if u_path="$(type -P "$util" 2>/dev/null)" && [ -n "$u_path" ]; then
      ln -sf "$u_path" "$no_to_bin/$util"
    fi
  done
  py_to_rc=0
  PATH="$no_to_bin" "$1" --timeout 0.2s sh -c '\''sleep 30 & echo $! > "$1"; wait'\'' _ "$tmp/py_to_child.pid" || py_to_rc=$?
  [ "$py_to_rc" -eq 124 ] || { echo "FAIL: Python fallback --timeout 0.2s expected 124, got $py_to_rc" >&2; exit 1; }
  sleep 0.05
  if [ -s "$tmp/py_to_child.pid" ] && kill -0 "$(cat "$tmp/py_to_child.pid")" 2>/dev/null; then
    kill -9 "$(cat "$tmp/py_to_child.pid")" 2>/dev/null || true
    echo "FAIL: Python fallback --timeout left orphan child $(cat "$tmp/py_to_child.pid")" >&2
    exit 1
  fi
  PATH="$no_to_bin" "$1" --timeout 1h true || { echo "FAIL: Python fallback --timeout 1h true failed" >&2; exit 1; }
  py_exit_rc=0
  PATH="$no_to_bin" "$1" --timeout 10 sh -c "exit 19" || py_exit_rc=$?
  [ "$py_exit_rc" -eq 19 ] || { echo "FAIL: Python fallback did not preserve exit 19 (got $py_exit_rc)" >&2; exit 1; }

  # 15. find_gitleaks under minimal PATH=/usr/bin:/bin
  fake_home="$tmp/fake_home"
  mkdir -p "$fake_home/miniforge3/envs/vardict-oracle/bin" "$fake_home/.local/bin"
  printf "#!/bin/sh\nexit 0\n" > "$fake_home/miniforge3/envs/vardict-oracle/bin/gitleaks"
  chmod +x "$fake_home/miniforge3/envs/vardict-oracle/bin/gitleaks"
  got_gl1="$(HOME="$fake_home" CONDA_PREFIX="" G4_TOOLS_BIN="$tmp/nope" find_gitleaks)"
  [ "$got_gl1" = "$fake_home/miniforge3/envs/vardict-oracle/bin/gitleaks" ] || { echo "FAIL: find_gitleaks miniforge3 fallback got $got_gl1" >&2; exit 1; }
  rm -f "$fake_home/miniforge3/envs/vardict-oracle/bin/gitleaks"
  printf "#!/bin/sh\nexit 0\n" > "$fake_home/.local/bin/gitleaks"
  chmod +x "$fake_home/.local/bin/gitleaks"
  got_gl2="$(HOME="$fake_home" CONDA_PREFIX="" G4_TOOLS_BIN="$tmp/nope" find_gitleaks)"
  [ "$got_gl2" = "$fake_home/.local/bin/gitleaks" ] || { echo "FAIL: find_gitleaks .local/bin fallback got $got_gl2" >&2; exit 1; }

  # 16. resolve_data_root (explicit env, direct scratch/g4, and linked git worktree)
  got_dr_env="$(DATA_ROOT="/explicit/data/root" resolve_data_root "$tmp")"
  [ "$got_dr_env" = "/explicit/data/root" ] || { echo "FAIL: resolve_data_root explicit env got $got_dr_env" >&2; exit 1; }
  if command -v git >/dev/null 2>&1; then
    mkdir -p "$tmp/main_repo/scratch/g4"
    (
      cd "$tmp/main_repo"
      git init -q -b main
      git config user.email "t@example.com"
      git config user.name "T"
      echo ok > README.md
      git add README.md
      git commit -q -m init
      git worktree add -q -b wt "$tmp/wt_repo"
    )
    mkdir -p "$tmp/wt_repo/scratch/acceptance"
    got_dr_wt="$(DATA_ROOT="" resolve_data_root "$tmp/wt_repo")"
    want_dr_main="$(cd "$tmp/main_repo" && pwd)"
    [ "$got_dr_wt" = "$want_dr_main" ] || { echo "FAIL: resolve_data_root worktree got $got_dr_wt, want $want_dr_main" >&2; exit 1; }
  fi
' _ "$PORTABLE_SH" "$TMP_DIR"

echo "PASS test_portable_shell"
