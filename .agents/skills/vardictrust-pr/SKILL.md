---
name: vardictrust-pr
description: Use when opening, updating, reviewing, recertifying, or merging pull requests in VarDictRust, filling .github/pull_request_template.md, or replying inline to automated and human code review threads.
---

# VarDictRust Pull Request Lifecycle & Governance (`vardictrust-pr`)

## Overview

Every change to `VarDictRust` lands via a feature branch off `origin/main` and a pull request completing `.github/pull_request_template.md` (`AGENTS.md` §§0, 0b, 14, 15, 18). Direct pushes to `main` are blocked by `.githooks/pre-push`.

**Core principle:** Evidence before completion (`AGENTS.md` §0) and Class-0 integrity (`AGENTS.md` §0b). Never hand-edit `benchmarks/evidence/*.json`, never relax `--eps` or `--cols`, and never mark engine checkboxes `N/A` when engine, `Cargo.*`, `harness/`, or workflow files change.

---

## Quick Reference: Code-Change Filter & Required Gates

`.githooks/pre-push` inspects `git diff --name-only` against `^(crates/|Cargo\.toml$|Cargo\.lock$|harness/|\.github/workflows/)`:

| PR Scope | `.githooks/pre-push` Classification | `engine_fingerprint.sh` (`FP` / `PERF_FP`) | Required Verification & PR Template Policy |
| :--- | :--- | :--- | :--- |
| Docs/Wiki/Skills only (`docs/`, `wiki/`, `README.md`, `.agents/skills/`) | Docs/infra-only (`CODE_CHANGES` empty) | Unchanged | `cargo fmt`, `cargo clippy`, `bash harness/run_self_tests.sh`, `mkdocs build --strict`. Engine checkboxes in `.github/pull_request_template.md` may be marked `N/A` with a 1-line justification. |
| `Cargo.toml` / `Cargo.lock` workspace version bump or `harness/` test edits | `CODE_CHANGES` non-empty | Unchanged (`normalize_gate_blob` strips workspace SemVer) | Full `bash harness/pr_check.sh` (6/6) + `bash harness/acceptance.sh` (8/8). **All engine checkboxes in `.github/pull_request_template.md` MUST be checked `[x]`** (never `N/A`). Revert timestamp-only dirt in `extended-coverage-tier2.json` and `multimodal-slices-tier2.json` after `acceptance.sh`. |
| `crates/vardictrust-core/`, `crates/vardictrust/`, external deps, or gated `harness/` runners | `CODE_CHANGES` non-empty | **Changed** (`FP` and/or `PERF_FP` updates) | Run `bash harness/recertify.sh` (or `bash harness/courtesy.sh bash harness/recertify.sh`) to regenerate live evidence in `benchmarks/evidence/*.json`, update SHA-256 table in `docs/validation/omnibus-empirical-validation.md` §13, then run `pr_check.sh` + `acceptance.sh`. |

---

## End-to-End Pull Request Workflow

### 1. Branch & Hook Initialization
```bash
bash harness/setup_hooks.sh
git fetch origin
git checkout -b <type>/<slug> origin/main
```
- Allowed `<type>` prefixes: `feat`, `fix`, `perf`, `refactor`, `test`, `docs`, `chore`, `release`.
- **Commit hygiene (`AGENTS.md` §18):** Atomic commits (`git commit -m "<type>: <short desc>"`). **Never** include `Co-Authored-By` trailers.
- **Validation statement (`AGENTS.md` §15.5):** Engine commits/PRs must state:
  `validated vs VarDictJava 1.8.3, dataset <name>, cmd <…>`

### 2. Pre-PR Verification Sequence
Commit your changes first (`engine_fingerprint.sh` requires a clean working tree on gated paths), then run:
```bash
# 1. Format, strict Clippy lints, and 277-test workspace suite (zero warnings)
cargo fmt --all -- --check
cargo clippy --all-targets --all-features -- -D warnings
RUST_TEST_THREADS=4 cargo test --workspace

# 2. Check engine and performance fingerprints
bash harness/engine_fingerprint.sh
bash harness/engine_fingerprint.sh --perf

# 3. If FP or PERF_FP changed, recertify live evidence honestly (NEVER hand-edit JSONs)
# bash harness/courtesy.sh bash harness/recertify.sh

# 4. Run all 24 harness self-tests and strict MkDocs build check
bash harness/run_self_tests.sh
python3 harness/prepare_docs_site.py && mkdocs build --strict && rm -rf .docs_staged site

# 5. Run 6-stage pre-PR compound gate, 8-suite acceptance gate, and secret scan
bash harness/pr_check.sh
bash harness/acceptance.sh
# If FP did not change, restore timestamp-only updates from live Tier-2 runners:
git checkout -- benchmarks/evidence/extended-coverage-tier2.json benchmarks/evidence/multimodal-slices-tier2.json
gitleaks detect --source . -v
```

### 3. Push & Open Pull Request (`gh pr create`)
Complete every section of `.github/pull_request_template.md`:
```bash
git push -u origin <type>/<slug>
gh pr create --repo palakpsheth/VarDictRust \
  --base main \
  --head <type>/<slug> \
  --title "<type>(<scope>): <concise summary>" \
  --body-file /tmp/pr_body.md
```

### 4. Inline Code Review Resolution & CI Polling
When automated reviewers (`gemini-code-assist`, Copilot, CodeRabbit) or humans post inline review comments:
1. Fetch all review threads and inline comments:
   ```bash
   gh api repos/palakpsheth/VarDictRust/pulls/<PR>/comments \
     --jq '.[] | {id: .id, path: .path, line: .line, user: .user.login, body: .body}'
   ```
2. Verify each suggestion technically before implementing. Push a surgical fix commit, re-run verification gates, and reply **inline** to each comment ID with the commit SHA and verification receipt:
   ```bash
   gh api -X POST repos/palakpsheth/VarDictRust/pulls/<PR>/comments/<COMMENT_ID>/replies \
     -f body="Fixed in \`<commit-sha>\` — <concise technical explanation and verification receipt>."
   ```
3. Confirm all CI checks are green before squash-merging:
   ```bash
   gh pr checks <PR>
   gh pr merge <PR> --squash --delete-branch
   ```

---

## Rationalization Counter-Table

| Excuse / Rationalization | Reality / Mandatory Action |
| :--- | :--- |
| "Only a comment or variable name changed in `crates/`; I'll manually update `engine_fingerprint` in `benchmarks/evidence/*.json`." | **Class-0 Integrity Violation (`AGENTS.md` §0b).** Never hand-edit `benchmarks/evidence/*.json`. Run `bash harness/recertify.sh`. |
| "I only bumped `Cargo.toml` and `harness/tests/*`, so I'll mark the engine checkboxes `N/A` in the PR body." | `.githooks/pre-push` classifies `Cargo.*` and `harness/` as `CODE_CHANGES`. Run `pr_check.sh` + `acceptance.sh` and check `[x]` every box. |
| "I'll post one top-level PR comment instead of replying to each inline review comment." | Reply directly to each inline review thread (`/pulls/<PR>/comments/<id>/replies`) so reviewers and auditors can trace every resolution. |

## Red Flags — STOP Immediately If You Catch Yourself

- Editing any file under `benchmarks/evidence/*.json` or `golden/` by hand.
- Adding `unwrap()`, `expect()`, `panic!()`, `unreachable!()`, or `unsafe` in `crates/vardictrust-core/src/**` outside `#[cfg(test)]`.
- Including `Co-Authored-By` trailers in commit messages or squash-merge descriptions.
- Pushing binary genomic files (`.bam`, `.cram`, `.fa`, `.vcf` outside `golden/`), `target/`, `.docs_staged/`, or `site/`.
