#!/usr/bin/env bash
# Configure Git to use repository-managed hooks in .githooks/ (main checkout and git worktrees)
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

chmod +x "$ROOT_DIR/.githooks/pre-commit" "$ROOT_DIR/.githooks/pre-push"
git config core.hooksPath .githooks

echo "Success: Git pre-commit and pre-push hooks configured via 'git config core.hooksPath .githooks' (active across main checkout and git worktrees)."
echo "To test manually, run: .githooks/pre-commit or .githooks/pre-push"
echo "To revert to default hooks, run: git config --unset core.hooksPath"
