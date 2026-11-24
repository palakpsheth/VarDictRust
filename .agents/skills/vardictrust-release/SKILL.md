---
name: vardictrust-release
description: Use when preparing, verifying, tagging, or publishing a VarDictRust release across GitHub Releases (cargo-dist SLSA Level 3), Crates.io, Bioconda, Nextflow nf-core/modules, and ReadTheDocs.
---

# VarDictRust Multi-Channel Release Engineering (`vardictrust-release`)

## Overview

VarDictRust ships across **four binary/package channels** plus **ReadTheDocs**:
1. **GitHub Releases (`cargo-dist` + SLSA Level 3 provenance):** Prebuilt archives for `x86_64-unknown-linux-gnu`, `x86_64-unknown-linux-musl`, `aarch64-unknown-linux-musl`, `aarch64-apple-darwin`, and `x86_64-apple-darwin`, shell installer (`vardictrust-installer.sh`), and `sha256.sum`.
2. **Crates.io:** Sequential publish of `vardictrust-core` followed by `vardictrust`.
3. **Bioconda (`bioconda/bioconda-recipes`):** Conda package built from the tagged GitHub source tarball (`packaging/bioconda/`).
4. **Nextflow (`nf-core/modules`):** `vardictrust` module (`packaging/nf-core/modules/vardictrust/`) with multi-arch Wave Singularity/Docker containers built from Bioconda (`bioconda::vardictrust` + `bioconda::htslib`).
5. **ReadTheDocs (`vardictrust.readthedocs.io`):** Strictly built Zensical site (`zensical.toml` + `harness/prepare_docs_site.py` + `harness/trigger_rtd.sh`).

**Core principle:** Zero version drift across all 20 version-bearing files (`references/version-surface-map.md`) and zero fake SHA-256s or evidence edits (`AGENTS.md` §0b).

---

## Quick Reference: Release Phases & Dependency Order

| Phase | Branch / Context | Validator Flag | Key Deliverables & Upstream Dependency |
| :--- | :--- | :--- | :--- |
| **Phase 1 (`pre-tag`)** | `release/vX.Y.Z` PR $\to$ `main` | `--phase pre-tag --require-head-fingerprint` | Bump 18 pre-tag version surfaces (including `.github/SECURITY.md` `X.Y.x` when `major.minor` advances), author `docs/release/vX.Y.Z-pr.md`, run `bash harness/recertify.sh` (refreshing all 16 `benchmarks/evidence/*.json` files + `sync_evidence_docs.py` + `acceptance.sh` 8/8) + `pr_check.sh` (6/6), squash-merge to `main`. |
| **Phase 2 (`post-tag`)** | Annotated tag `vX.Y.Z` on `main` | `--phase post-tag` | Push `vX.Y.Z` tag $\to$ `.github/workflows/release.yml` publishes GitHub Release + SLSA provenance + Crates.io (`vardictrust-core`, `vardictrust`). |
| **Phase 3A (`post-bioconda`)** | Post-GitHub tarball `vX.Y.Z.tar.gz` | `--phase post-bioconda` | Run `bash packaging/bioconda/prepare_bioconda_pr.sh X.Y.Z` + `bash harness/test_bioconda_recipe.sh`, sync `meta.yaml` to `VarDictRust` `main` (`chore/bioconda-vX.Y.Z`), open PR on `bioconda/bioconda-recipes`. |
| **Phase 3B (`post-nfcore`)** | Post-Anaconda `linux-64` + `linux-aarch64` upload | `--phase post-nfcore` | Once `api.anaconda.org/package/bioconda/vardictrust` lists both `linux-64` and `linux-aarch64` tarballs for `X.Y.Z`, run `bash packaging/nf-core/prepare_nfcore_pr.sh X.Y.Z`, sync `packaging/nf-core/`, and push to `nf-core/modules`. |
| **Phase 3C (`docs` & portals)** | Post-merge on `main` | `--phase post-nfcore` | Trigger `bash harness/trigger_rtd.sh build vardictrust latest` and verify ReadTheDocs project metadata + ELIXIR `bio.tools` (`biotools:vardictrust`). |

---

## Step-by-Step Execution Runbook

### Phase 1: Pre-Tag Release Branch & PR (`--phase pre-tag`)
1. Create branch `git checkout -b release/vX.Y.Z origin/main`.
2. Update all 18 pre-tag files listed in [`references/version-surface-map.md`](references/version-surface-map.md) (updating `.github/SECURITY.md` Supported Versions to `X.Y.x` `:white_check_mark:` whenever `major.minor` changes) and run `cargo update -w`.
3. Create `docs/release/vX.Y.Z-pr.md` and add a superseded note at the top of the prior `docs/release/v*-pr.md`.
4. Commit the version surfaces (`git commit -m "release: prepare VarDictRust vX.Y.Z"`, zero `Co-Authored-By` trailers), then run `bash harness/recertify.sh` to refresh all 15 non-scale `benchmarks/evidence/*.json` files to the release `FP`, auto-sync `README.md` / `omnibus-empirical-validation.md` / `DIVERGENCES.md` / `test-matrix-and-ci-audit.md` / `benchmarks/README.md` / `docs/ecosystem-survey.md` / `docs/release/vX.Y.Z-pr.md` / `version-surface-map.md` via `python3 harness/sync_evidence_docs.py`, re-compile `paper/` manuscript artifacts via `python3 paper/build_manuscript.py`, and verify the full gate suite:
   ```bash
   # 1. Recertify all 15 non-scale evidence JSONs + sync documentation SHA-256/FP tables + rebuild paper/ + run acceptance.sh (8/8)
   bash harness/recertify.sh
   python3 paper/build_manuscript.py
   git add benchmarks/evidence/ README.md docs/ DIVERGENCES.md benchmarks/README.md .agents/skills/vardictrust-release/references/version-surface-map.md paper/
   git commit --amend --no-edit

   # 2. Run bundled consistency checker across all 18 pre-tag surfaces + 16 evidence JSONs + live HEAD fingerprint
   python3 .agents/skills/vardictrust-release/scripts/check_release_consistency.py \
     --version X.Y.Z --phase pre-tag --require-head-fingerprint

   # 3. Verify self-tests, strict Zensical build, Cargo packaging, and 6-stage PR gate
   bash harness/run_self_tests.sh
   python3 harness/prepare_docs_site.py && zensical build --clean --strict && rm -rf .docs_staged site .cache
   cargo publish -p vardictrust-core --dry-run
   cargo package --list -p vardictrust
   bash harness/pr_check.sh
   bash harness/acceptance.sh
   gitleaks detect --source . -v
   ```
5. Push `release/vX.Y.Z`, open PR completing `.github/pull_request_template.md` with **all engine checkboxes checked `[x]`**, resolve any inline review comments, and squash-merge to `main`.

### Phase 2: Tag & Automated CI Release (`--phase post-tag`)
```bash
git checkout main && git pull --ff-only origin main
git tag -a vX.Y.Z -m "VarDictRust vX.Y.Z"
git push origin vX.Y.Z

# Monitor .github/workflows/release.yml (plan -> build-local-artifacts -> build-global-artifacts -> host -> attest-release -> publish-crates)
# Note: v* tags are immutable under tag-protection-v*; if release.yml needs re-dispatch from main after a CI/runner fix (with only CI/test/skill diffs vs tag; all 5 downstream jobs check out refs/tags/vX.Y.Z), run:
#   gh workflow run release.yml --repo palakpsheth/VarDictRust --ref main -f tag=vX.Y.Z
gh run list --workflow=release.yml --limit 1
gh run watch <run-id> --exit-status
python3 .agents/skills/vardictrust-release/scripts/check_release_consistency.py --version X.Y.Z --phase post-tag
```

### Phase 3A: Bioconda Recipe Update (`--phase post-bioconda`)
```bash
bash packaging/bioconda/prepare_bioconda_pr.sh X.Y.Z
bash harness/test_bioconda_recipe.sh
python3 .agents/skills/vardictrust-release/scripts/check_release_consistency.py --version X.Y.Z --phase post-bioconda

# Sync updated meta.yaml to VarDictRust main via PR, then push branch to bioconda-recipes fork
git checkout -b chore/bioconda-vX.Y.Z origin/main
git add packaging/bioconda/recipes/vardictrust/meta.yaml
git commit -m "chore(bioconda): update in-tree meta.yaml to vX.Y.Z sha256"
git push -u origin chore/bioconda-vX.Y.Z

# Once all 4 Bioconda CI checks pass (Lint, Linux Tests, OSX-64 Tests, CircleCI build and test (ARM)):
# 1. Issue the official BiocondaBot command to apply the "please review & merge" label:
gh pr comment <bioconda-pr> --repo bioconda/bioconda-recipes --body "@BiocondaBot please add label"
# 2a. If PR was opened by BiocondaBot (autobump), approve as recipe maintainer:
gh pr review <bioconda-pr> --repo bioconda/bioconda-recipes --approve
# 2b. If PR was opened by palakpsheth from a fork (GitHub blocks self-APPROVE), post maintainer sign-off and ping @bioconda/core (which BiocondaBot reposts):
gh pr review <bioconda-pr> --repo bioconda/bioconda-recipes --comment --body "Maintainer sign-off (@palakpsheth): all 4 multi-arch CI checks are green. Ready for merge."
gh pr comment <bioconda-pr> --repo bioconda/bioconda-recipes --body "All CI checks (Lint, Linux x86_64, Linux ARM64, macOS x86_64) are green. Ready for review and merge from @bioconda/core. Thank you!"
```

### Phase 3B & 3C: Nextflow `nf-core/modules` & ReadTheDocs (`--phase post-nfcore`)
```bash
# Run once bioconda::vardictrust=X.Y.Z has both linux-64 and linux-aarch64 tarballs on anaconda.org
NFCORE_MODULES_DIR=/path/to/nf-core-modules bash packaging/nf-core/prepare_nfcore_pr.sh X.Y.Z
python3 .agents/skills/vardictrust-release/scripts/check_release_consistency.py --version X.Y.Z --phase post-nfcore

# Once the upstream nf-core/modules PR passes all CI checks (conda, docker, singularity, nf-core lint, pre-commit)
# and copilot-pull-request-reviewer[bot] review:
# Etiquette rule (per @famosab on nf-core/modules#13830): DO NOT mass-ping @nf-core/maintainers
# or request "team_reviewers":["maintainers"] on routine nf-core/modules version-bump PRs!
# 1. Apply the "Ready for Review" and "update module" labels:
gh api repos/nf-core/modules/issues/<nfcore-pr>/labels -X POST --input - <<< '{"labels":["Ready for Review","update module"]}'
# 2. Request review in the nf-core Slack #request-review channel (or request a specific prior module
#    maintainer such as @nvnieuwk or @SPPearce individually only when appropriate, without @nf-core/maintainers):
#    e.g. post the PR URL in #request-review on nf-core Slack.
# 3. Once 1 approving review lands, enqueue into the nf-core/modules merge queue:
PR_NODE_ID="$(gh api repos/nf-core/modules/pulls/<nfcore-pr> --jq '.node_id')"
gh api graphql -f query='mutation($id: ID!) { enqueuePullRequest(input: {pullRequestId: $id}) { mergeQueueEntry { id state position } } }' -f id="$PR_NODE_ID"

# Trigger ReadTheDocs latest build and sync external registry metadata:
bash harness/trigger_rtd.sh build vardictrust latest
# Verify ReadTheDocs admin project metadata (description & project_url) at https://app.readthedocs.org/projects/vardictrust/
# and update ELIXIR bio.tools (https://bio.tools/vardictrust) version X.Y.Z, SPDX license, and dual-mode description.
```

---

## Rationalization Counter-Table

| Excuse / Rationalization | Reality / Mandatory Action |
| :--- | :--- |
| "I'll bump `Cargo.toml` and `CHANGELOG.md` only, like `v1.0.2` did." | `v1.0.2` left 15 files at `1.0.1` (`crates/vardictrust/Cargo.toml`, `CITATION.cff`, `README.md`, etc.). Always run `bash harness/recertify.sh` and `check_release_consistency.py --version X.Y.Z --phase pre-tag --require-head-fingerprint`. |
| "I'll compute `meta.yaml`'s SHA-256 from a local `git archive` before pushing the `vX.Y.Z` tag." | GitHub's generated source tarball can differ from local `git archive` metadata. Only run `prepare_bioconda_pr.sh X.Y.Z` in Phase 3A after `vX.Y.Z` is tagged on GitHub. |
| "I'll run `nf-core modules containers create` as soon as the Bioconda PR merges." | Seqera Wave requires both `linux-64` and `linux-aarch64` packages to be indexed on `anaconda.org`. `prepare_nfcore_pr.sh` verifies both architectures before invoking Wave. |
| "I'll tag `@nf-core/maintainers` or pass `"team_reviewers":["maintainers"]` on the `nf-core/modules` bump PR." | Per `@famosab` (`nf-core/modules#13830`), never mass-ping `@nf-core/maintainers` on routine module version bumps. Apply the `"Ready for Review"` and `"update module"` labels and post in `nf-core` Slack `#request-review` (or request a prior module maintainer individually without `@nf-core/maintainers`). |

## Red Flags — STOP Immediately If You Catch Yourself

- Pushing `vX.Y.Z` tag before the `release/vX.Y.Z` PR has been squash-merged to `main`.
- Running `cargo publish -p vardictrust --dry-run` before `vardictrust-core = X.Y.Z` is published on `crates.io` (use `cargo publish -p vardictrust-core --dry-run && cargo package --list -p vardictrust` in Phase 1).
- Editing `version`, `sha256`, or `build.number` in `packaging/bioconda/recipes/vardictrust/meta.yaml` or container/snapshot pins in `packaging/nf-core/modules/vardictrust/tests/main.nf.test.snap` by hand at release time instead of via `prepare_bioconda_pr.sh` and `prepare_nfcore_pr.sh`.
