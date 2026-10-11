#!/usr/bin/env bash
# Configure Git to use repository-managed hooks in .githooks/
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

chmod +x .githooks/pre-commit
git config core.hooksPath .githooks

echo "Success: Git pre-commit hook configured via 'git config core.hooksPath .githooks'."
echo "To test manually, run: .githooks/pre-commit"
echo "To revert to default hooks, run: git config --unset core.hooksPath"
