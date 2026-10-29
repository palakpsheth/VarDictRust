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

`.githooks/pre-push` inspects `git diff --name-only` against `\.rs$|Cargo\.(toml|lock)$|^crates/|^harness/|^\.github/|^\.githooks/|^\.agents/|^fixtures/|^golden/|^benchmarks/evidence/`:

| PR Scope | `.githooks/pre-push` Classification | `engine_fingerprint.sh` (`FP` / `PERF_FP`) | Required Verification & PR Template Policy |
| :--- | :--- | :--- | :--- |
| Docs/Wiki only (`docs/`, `wiki/`, `README.md`) | Docs/infra-only (`CODE_CHANGES` empty) | Unchanged | `cargo fmt`, `cargo clippy`, `bash harness/run_self_tests.sh`, `mkdocs build --strict`. Engine checkboxes in `.github/pull_request_template.md` may be marked `N/A` with a 1-line justification. |
| Feature, Bugfix, Harness, Workflow, or Skill PR (`crates/`, `Cargo.*`, `harness/`, `.github/`, `.githooks/`, `.agents/`) | `CODE_CHANGES` non-empty | Unchanged or Changed | **Tier A Per-PR Gate (`<25s`):** Run `bash harness/pr_check.sh` (6/6: workspace tests, 92-fixture corpus, `micro` + `s1a..s9`, `run_extended_coverage.sh --no-record`, `run_multimodal_slices.sh --no-record`, `perf_budget.py`, `gitleaks`). Do **not** rewrite `benchmarks/evidence/*.json` on normal PRs; mark the Release Recertification checkbox `N/A — deferred to release cut per two-tier testing policy`. |
| Release Branch (`release/vX.Y.Z`) or Milestone Certification PR | `CODE_CHANGES` non-empty | Certified at `HEAD` | **Tier B Release Gate:** Run `bash harness/recertify.sh` (which refreshes all 15 non-scale `benchmarks/evidence/*.json` files, runs `python3 harness/sync_evidence_docs.py`, and executes `bash harness/acceptance.sh` 8/8) + `bash harness/pr_check.sh` (6/6). (`scale_bench.sh` runs only if `PERF_FP` changed since `scale.json`.) |

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
Commit your changes first, then run:
```bash
# 1. Run all 25 harness self-tests and strict MkDocs build check
bash harness/run_self_tests.sh
python3 harness/prepare_docs_site.py && mkdocs build --strict && rm -rf .docs_staged site

# 2. Tier A (Every PR): Run 6-stage pre-PR compound gate (<25s, --no-record on Tier-2 suites)
bash harness/pr_check.sh

# 3. Tier B (Release branch release/vX.Y.Z or milestone certification ONLY):
# bash harness/recertify.sh
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

### 4. Mandatory Multi-Skill Review, Copilot Review Wait & Thread Resolution (`AGENTS.md` §15.8)
Every PR in `VarDictRust` **must** complete both the multi-skill review gate and the Copilot review thread-resolution gate before merging:
1. **Run Multi-Skill Review (`/ponytail-review`, `/rust-router`, and `/pragmatic-rust-guidelines`):**
   Audit `git diff origin/main...HEAD` through all three skills (`/ponytail-review`, `/rust-router`, and `/pragmatic-rust-guidelines`) and fix all findings using Red/Green TDD before merging.
2. **Wait for `copilot-pull-request-reviewer[bot]` Review Completion:**
   **Never** merge a PR solely because CI checks are green while `copilot-pull-request-reviewer[bot]` review is still pending. Poll until `copilot-pull-request-reviewer[bot]` posts its review on the PR:
   ```bash
   gh api --paginate repos/palakpsheth/VarDictRust/pulls/<PR>/reviews \
     --jq '.[] | {id: .id, user: .user.login, state: .state, commit_id: .commit_id}'
   gh api --paginate repos/palakpsheth/VarDictRust/pulls/<PR>/comments \
     --jq '.[] | {id: .id, path: .path, line: .line, user: .user.login, body: .body}'
   ```
3. **Refute or Fix, Reply Inline, and Resolve Every Review Thread (`resolveReviewThread`):**
   Verify each review comment technically. Fix valid findings (or refute invalid ones with code/oracle evidence), push and re-verify gates, reply **inline** to each comment ID, and resolve every review thread via GraphQL before merging:
   ```bash
   # Reply inline to each review comment:
   gh api -X POST repos/palakpsheth/VarDictRust/pulls/<PR>/comments/<COMMENT_ID>/replies \
     -f body="Fixed in \`<commit-sha>\` — <concise technical explanation and verification receipt>."

   # Query all review thread node IDs with GraphQL cursor pagination and resolve each thread:
   gh api graphql --paginate -f query='
     query($pr: Int!, $endCursor: String) {
       repository(owner: "palakpsheth", name: "VarDictRust") {
         pullRequest(number: $pr) {
           reviewThreads(first: 100, after: $endCursor) {
             pageInfo { hasNextPage endCursor }
             nodes { id isResolved comments(first: 1) { nodes { databaseId path } } }
           }
         }
       }
     }' -F pr=<PR>

   gh api graphql -f query='
     mutation($threadId: ID!) {
       resolveReviewThread(input: {threadId: $threadId}) {
         thread { id isResolved }
       }
     }' -f threadId="<THREAD_NODE_ID>"
   ```
4. **Confirm All CI Checks Green & Zero Unresolved Threads Before Squash-Merging:**
   ```bash
   gh pr checks <PR>
   gh pr merge <PR> --squash --delete-branch
   ```

---

## Rationalization Counter-Table

| Excuse / Rationalization | Reality / Mandatory Action |
| :--- | :--- |
| "Only a comment or variable name changed in `crates/`; I'll manually update `engine_fingerprint` in `benchmarks/evidence/*.json`." | **Class-0 Integrity Violation (`AGENTS.md` §0b).** Never hand-edit `benchmarks/evidence/*.json`. Routine feature/bugfix PRs run `bash harness/pr_check.sh` (`<25s`, `--no-record`) without touching `benchmarks/evidence/*.json`; release certification PRs (`release/vX.Y.Z`) run `bash harness/recertify.sh`. |
| "I only bumped `Cargo.toml` and `harness/tests/*`, so I'll mark the engine checkboxes `N/A` in the PR body." | `.githooks/pre-push` classifies `Cargo.*` and `harness/` as `CODE_CHANGES`. Run `bash harness/pr_check.sh` (6/6) and check `[x]` the Tier A engine checkboxes (marking the Tier B Release Recertification checkbox `N/A — deferred to release cut per two-tier testing policy` on routine PRs, or running `recertify.sh` + `acceptance.sh` on `release/vX.Y.Z` PRs). |
| "I'll post one top-level PR comment instead of replying to each inline review comment." | Reply directly to each inline review thread (`/pulls/<PR>/comments/<id>/replies`) and resolve every review thread via `resolveReviewThread` so reviewers and auditors can trace every resolution. |
| "All 11 CI checks are green, so I can merge right away without waiting for Copilot or running the review skills." | **Forbidden (`AGENTS.md` §15.8).** Every PR must run `/ponytail-review`, `/rust-router`, and `/pragmatic-rust-guidelines`, wait for `copilot-pull-request-reviewer[bot]` to complete its review, refute/fix + reply inline to every comment, and resolve all review threads (`resolveReviewThread`) before squash-merging. |

## Red Flags — STOP Immediately If You Catch Yourself

- Editing any file under `benchmarks/evidence/*.json` or `golden/` by hand.
- Adding `unwrap()`, `expect()`, `panic!()`, `unreachable!()`, or `unsafe` in `crates/vardictrust-core/src/**` outside `#[cfg(test)]`.
- Including `Co-Authored-By` trailers in commit messages or squash-merge descriptions.
- Pushing binary genomic files (`.bam`, `.cram`, `.fa`, `.vcf` outside `golden/`), `target/`, `.docs_staged/`, or `site/`.
