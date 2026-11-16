#!/usr/bin/env bash
# harness/tests/test_git_hooks.sh — Self-test for .githooks/pre-commit, .githooks/pre-push,
# harness/setup_hooks.sh, and the Git-2.34-compatible 3-way mergeability guard.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$ROOT_DIR"

fail=0
ck() {
  local desc="$1"
  shift
  if ! "$@"; then
    echo "FAIL test_git_hooks: $desc" >&2
    fail=1
  fi
}

# 1. Executable bits & syntax checks on hooks and setup script
ck "setup_hooks.sh exists and passes bash -n" bash -n harness/setup_hooks.sh
ck "pre-commit passes bash -n" bash -n .githooks/pre-commit
ck "pre-push passes bash -n" bash -n .githooks/pre-push
ck "pr_check.sh passes bash -n" bash -n harness/pr_check.sh
ck "pre-commit is executable" test -x .githooks/pre-commit
ck "pre-push is executable" test -x .githooks/pre-push

# 2. setup_hooks.sh contract: chmods both hooks and mentions both hooks + worktrees
ck "setup_hooks.sh chmods both pre-commit and pre-push" grep -qE 'chmod \+x.*\.githooks/pre-commit.*\.githooks/pre-push' harness/setup_hooks.sh
setup_out=$(bash harness/setup_hooks.sh 2>&1)
ck "setup_hooks.sh output mentions pre-commit and pre-push" grep -q 'pre-commit' <<<"$setup_out"
ck "setup_hooks.sh output mentions pre-push" grep -q 'pre-push' <<<"$setup_out"
ck "setup_hooks.sh output mentions worktree" grep -qi 'worktree' <<<"$setup_out"

# 3. Timestamp helpers & courteous thread caps in pre-commit, pre-push, and pr_check.sh
for f in .githooks/pre-commit .githooks/pre-push harness/pr_check.sh; do
  ck "$f defines ts() helper" grep -q "ts() { date '+\\[%H:%M:%S\\]'; }" "$f"
  ck "$f exports RUST_TEST_THREADS" grep -q 'RUST_TEST_THREADS="\${RUST_TEST_THREADS:-4}"' "$f"
  ck "$f exports RAYON_NUM_THREADS" grep -q 'RAYON_NUM_THREADS="\${RAYON_NUM_THREADS:-4}"' "$f"
done
ck "pre-push exports RUST_THREADS=1" grep -q 'RUST_THREADS="\${RUST_THREADS:-1}"' .githooks/pre-push
ck "pr_check.sh exports RUST_THREADS=1" grep -q 'RUST_THREADS="\${RUST_THREADS:-1}"' harness/pr_check.sh

# 4. pre-commit fast staged checks & ACMRD filter
ck "pre-commit resolves ROOT_DIR via git rev-parse --show-toplevel" grep -q 'git rev-parse --show-toplevel' .githooks/pre-commit
ck "pre-commit uses --diff-filter=ACMRD for Rust/Cargo detection" grep -q -- '--diff-filter=ACMRD' .githooks/pre-commit
ck "pre-commit syntax-checks staged shell scripts with bash -n" grep -q 'bash -n' .githooks/pre-commit
ck "pre-commit compiles staged Python files with py_compile" grep -q 'py_compile' .githooks/pre-commit
ck "pre-commit runs test_stdlib_purity.py on staged Python changes" grep -q 'test_stdlib_purity\.py' .githooks/pre-commit
ck "pre-commit runs check_release_consistency.py on workflow/release surface changes" grep -q 'check_release_consistency\.py' .githooks/pre-commit
ck "pre-commit runs TestMarkdownStaticHygiene on staged Markdown/docs changes" grep -q 'TestMarkdownStaticHygiene' .githooks/pre-commit
ck "pre-commit blocks partially staged documentation files" grep -q 'UNSTAGED_MD=' .githooks/pre-commit
ck "pre-commit matches zensical.toml in staged doc regex" grep -Fq 'zensical\.toml' .githooks/pre-commit
ck "pre-commit matches .readthedocs.yaml in staged doc regex" grep -Fq '\.readthedocs\.yaml' .githooks/pre-commit
ck "pre-commit matches ^docs/ in staged doc regex" grep -Fq '^docs/' .githooks/pre-commit
ck "pre-push matches zensical.toml in DOC_CHANGES" grep -Fq 'zensical\.toml' .githooks/pre-push
ck "pre-push matches .readthedocs.yaml in DOC_CHANGES" grep -Fq '\.readthedocs\.yaml' .githooks/pre-push
ck "pre-push matches ^docs/ in DOC_CHANGES" grep -Fq '^docs/' .githooks/pre-push
ck "pre-push runs prepare_docs_site.py on doc changes" grep -q 'prepare_docs_site\.py' .githooks/pre-push
ck "pre-push exports pushed commit via git archive for doc build" grep -q 'git archive "\$local_oid"' .githooks/pre-push
ck "pre-push uses ZENSICAL_CMD bash array" grep -q 'ZENSICAL_CMD=()' .githooks/pre-push
ck "pre-push runs zensical build --clean --strict via array" grep -Fq '"${ZENSICAL_CMD[@]}" build --clean --strict' .githooks/pre-push
ck "pre-push fails fast when zensical is missing on doc changes" grep -q "Install Zensical before pushing documentation updates" .githooks/pre-push

# 5. pre-push runs engine_behavior and robustness in dev profile; acceptance.sh Item 4 runs them in --release
ck "pre-push runs engine_behavior (dev)" grep -q 'cargo test -p vardictrust-core --test engine_behavior' .githooks/pre-push
ck "pre-push runs robustness (dev)" grep -q 'cargo test -p vardictrust --test robustness' .githooks/pre-push
ck "acceptance.sh runs engine_behavior (--release)" grep -q 'cargo test --release -p vardictrust-core --test engine_behavior' harness/acceptance.sh
ck "acceptance.sh runs robustness (--release)" grep -q 'cargo test --release -p vardictrust --test robustness' harness/acceptance.sh

# 6. pre-push expanded path detection & diff base fallback
for pat in '\.github/' '\.githooks/' '\.agents/' 'fixtures/' 'golden/' 'benchmarks/evidence/' 'deny\.toml' '\.cargo/'; do
  ck "pre-push CODE_CHANGES matches $pat" grep -Fq "$pat" .githooks/pre-push
done
ck "pre-push has main...local_oid fallback" grep -q 'main\.\.\."\$local_oid"' .githooks/pre-push

# 6b. Supply-chain checks (cargo deny & cargo audit) in pre-commit, pr_check.sh, and acceptance.sh
for f in .githooks/pre-commit harness/pr_check.sh harness/acceptance.sh; do
  ck "$f invokes cargo deny --all-features check advisories licenses bans sources" grep -q -- '--all-features check advisories licenses bans sources' "$f"
  ck "$f invokes cargo audit" grep -q 'cargo audit' "$f"
done

# 7. Git-2.34-compatible 3-way merge-tree conflict guard in pre-push and pr_check.sh
for f in .githooks/pre-push harness/pr_check.sh; do
  ck "$f uses 3-way git merge-tree conflict check" grep -q 'git merge-tree "\$base" origin/main HEAD' "$f"
  ck "$f checks for ^+<<<<<<< \\.our$" grep -q '\^+<<<<<<< \\\.our\$' "$f"
  ck "$f instructs git fetch origin main && git rebase origin/main" grep -q 'git fetch origin main && git rebase origin/main' "$f"
done

# 8. Functional test of the 3-way git merge-tree conflict guard in an isolated temp repo
TMP_REPO=$(mktemp -d)
trap 'rm -rf "$TMP_REPO"' EXIT
(
  cd "$TMP_REPO"
  git init -q -b main origin.git --bare
  git clone -q origin.git work 2>/dev/null
  cd work
  git config user.email "test@example.com"
  git config user.name "Test"
  echo "line1" > file.txt
  git add file.txt
  git commit -q -m "initial"
  git push -q origin main

  # Create feature branch modifying file.txt
  git checkout -q -b feat/conflict
  echo "feature-change" > file.txt
  git commit -q -a -m "feature change"

  # Advance origin/main with a conflicting change
  git checkout -q main
  echo "main-change" > file.txt
  git commit -q -a -m "main change"
  git push -q origin main

  # Switch back to feat/conflict and verify 3-way git merge-tree detects conflict
  git checkout -q feat/conflict
  base=$(git merge-base origin/main HEAD)
  if ! git merge-tree "$base" origin/main HEAD | grep -q '^+<<<<<<< \.our$'; then
    echo "FAIL: git merge-tree did not detect expected conflict" >&2
    exit 1
  fi

  # Verify pre-push fails fast on conflict with origin/main before running cargo
  local_oid=$(git rev-parse HEAD)
  zero_oid="0000000000000000000000000000000000000000"
  set +e
  push_out=$(printf "refs/heads/feat/conflict %s refs/heads/feat/conflict %s\n" "$local_oid" "$zero_oid" | bash "$ROOT_DIR/.githooks/pre-push" 2>&1)
  push_rc=$?
  set -e
  if [ "$push_rc" -eq 0 ]; then
    echo "FAIL: pre-push succeeded despite merge conflict with origin/main" >&2
    exit 1
  fi
  if ! grep -q 'git fetch origin main && git rebase origin/main' <<<"$push_out"; then
    echo "FAIL: pre-push conflict output missing rebase instruction: $push_out" >&2
    exit 1
  fi
) || fail=1

# 9. Functional test of pre-push ZENSICAL_CMD=() with spaces in checkout path + missing-binary fast-fail
(
  cd "$TMP_REPO"
  git init -q -b main origin_docs.git --bare
  git clone -q origin_docs.git "checkout with spaces" 2>/dev/null
  cd "checkout with spaces"
  git config user.email "test@example.com"
  git config user.name "Test"

  mkdir -p .githooks harness
  cp "$ROOT_DIR/.githooks/pre-push" .githooks/pre-push
  chmod +x .githooks/pre-push

  cat > harness/prepare_docs_site.py <<'PYEOF'
import pathlib
p = pathlib.Path(__file__).resolve().parent.parent / ".docs_staged"
p.mkdir(parents=True, exist_ok=True)
(p / "README.md").write_text("# Test\n", encoding="utf-8")
PYEOF
  cat > zensical.toml <<'TOMLEOF'
[project]
site_name = "TestDocs"
docs_dir = ".docs_staged"
site_dir = "site"
TOMLEOF
  echo "# Initial" > README.md

  git add .githooks/pre-push harness/prepare_docs_site.py zensical.toml README.md
  git commit -q -m "initial docs setup"
  git push -q origin main

  git checkout -q -b docs/spaced-path-check
  echo "# Updated" >> README.md
  git commit -q -a -m "docs: update readme"

  local_oid=$(git rev-parse HEAD)
  zero_oid="0000000000000000000000000000000000000000"
  mkdir -p "$TMP_REPO/empty_home"

  # 9a. Missing zensical binary must fail fast with clear error
  set +e
  missing_out=$(printf "refs/heads/docs/spaced-path-check %s refs/heads/docs/spaced-path-check %s\n" "$local_oid" "$zero_oid" | env -i PATH="/usr/bin:/bin" HOME="$TMP_REPO/empty_home" bash "$TMP_REPO/checkout with spaces/.githooks/pre-push" 2>&1)
  missing_rc=$?
  set -e
  if [ "$missing_rc" -eq 0 ]; then
    echo "FAIL: pre-push succeeded when zensical binary was missing" >&2
    exit 1
  fi
  if ! grep -q "Install Zensical before pushing documentation updates" <<<"$missing_out"; then
    echo "FAIL: pre-push missing-zensical output lacked actionable install hint: $missing_out" >&2
    exit 1
  fi

  # 9b. .venv_docs/bin/zensical inside a checkout path with spaces must execute cleanly
  mkdir -p ".venv_docs/bin"
  MARKER_FILE="$TMP_REPO/zensical_invoked.log"
  cat > ".venv_docs/bin/zensical" <<SH_EOF
#!/usr/bin/env bash
printf '%s\n' "\$*" >> "$MARKER_FILE"
exit 0
SH_EOF
  chmod +x ".venv_docs/bin/zensical"

  set +e
  spaced_out=$(printf "refs/heads/docs/spaced-path-check %s refs/heads/docs/spaced-path-check %s\n" "$local_oid" "$zero_oid" | env -i PATH="/usr/bin:/bin" HOME="$TMP_REPO/empty_home" bash "$TMP_REPO/checkout with spaces/.githooks/pre-push" 2>&1)
  spaced_rc=$?
  set -e
  if [ "$spaced_rc" -ne 0 ]; then
    echo "FAIL: pre-push failed in spaced checkout path: $spaced_out" >&2
    exit 1
  fi
  if [ ! -f "$MARKER_FILE" ] || ! grep -q '^build --clean --strict$' "$MARKER_FILE"; then
    echo "FAIL: .venv_docs/bin/zensical was not invoked with 'build --clean --strict' in spaced path" >&2
    exit 1
  fi
) || fail=1

if [ "$fail" -eq 0 ]; then
  echo "PASS test_git_hooks"
  exit 0
else
  exit 1
fi
