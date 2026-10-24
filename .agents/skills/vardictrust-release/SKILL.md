---
name: vardictrust-release
description: Use when preparing, verifying, tagging, or publishing a VarDictRust release across GitHub Releases (cargo-dist SLSA Level 3), Crates.io, Bioconda, Nextflow nf-core/modules, and ReadTheDocs.
---

# VarDictRust Multi-Channel Release Engineering (`vardictrust-release`)

## Overview

VarDictRust ships across **four binary/package channels** plus **ReadTheDocs**:
1. **GitHub Releases (`cargo-dist` + SLSA Level 3 provenance):** Prebuilt archives for `x86_64-unknown-linux-gnu`, `x86_64-unknown-linux-musl`, `aarch64-unknown-linux-musl`, and `aarch64-apple-darwin`, shell installer (`vardictrust-installer.sh`), and `sha256.sum`.
2. **Crates.io:** Sequential publish of `vardictrust-core` followed by `vardictrust`.
3. **Bioconda (`bioconda/bioconda-recipes`):** Conda package built from the tagged GitHub source tarball (`packaging/bioconda/`).
4. **Nextflow (`nf-core/modules`):** `vardictrust` module (`packaging/nf-core/modules/vardictrust/`) with multi-arch Wave Singularity/Docker containers built from Bioconda (`bioconda::vardictrust` + `bioconda::htslib`).
5. **ReadTheDocs (`vardictrust.readthedocs.io`):** Strictly built MkDocs site (`harness/prepare_docs_site.py` + `harness/trigger_rtd.sh`).

**Core principle:** Zero version drift across all 19 version-bearing files (`references/version-surface-map.md`) and zero fake SHA-256s or evidence edits (`AGENTS.md` §0b).

---

## Quick Reference: Release Phases & Dependency Order

| Phase | Branch / Context | Validator Flag | Key Deliverables & Upstream Dependency |
| :--- | :--- | :--- | :--- |
| **Phase 1 (`pre-tag`)** | `release/vX.Y.Z` PR $\to$ `main` | `--phase pre-tag` | Bump 17 pre-tag version surfaces, author `docs/release/vX.Y.Z-pr.md`, verify `FP` & `PERF_FP`, run `pr_check.sh` (6/6) + `acceptance.sh` (8/8), squash-merge to `main`. |
| **Phase 2 (`post-tag`)** | Annotated tag `vX.Y.Z` on `main` | `--phase post-tag` | Push `vX.Y.Z` tag $\to$ `.github/workflows/release.yml` publishes GitHub Release + SLSA provenance + Crates.io (`vardictrust-core`, `vardictrust`). |
| **Phase 3A (`post-bioconda`)** | Post-GitHub tarball `vX.Y.Z.tar.gz` | `--phase post-bioconda` | Run `bash packaging/bioconda/prepare_bioconda_pr.sh X.Y.Z` + `bash harness/test_bioconda_recipe.sh`, sync `meta.yaml` to `VarDictRust` `main` (`chore/bioconda-vX.Y.Z`), open PR on `bioconda/bioconda-recipes`. |
| **Phase 3B (`post-nfcore`)** | Post-Anaconda `linux-64` + `linux-aarch64` upload | `--phase post-nfcore` | Once `api.anaconda.org/package/bioconda/vardictrust` lists both `linux-64` and `linux-aarch64` tarballs for `X.Y.Z`, run `bash packaging/nf-core/prepare_nfcore_pr.sh X.Y.Z`, sync `packaging/nf-core/`, and push to `nf-core/modules`. |
| **Phase 3C (`docs`)** | Post-merge on `main` | `--phase post-nfcore` | Trigger `bash harness/trigger_rtd.sh build vardictrust latest`. |

---

## Step-by-Step Execution Runbook

### Phase 1: Pre-Tag Release Branch & PR (`--phase pre-tag`)
1. Create branch `git checkout -b release/vX.Y.Z origin/main`.
2. Update all 17 pre-tag files listed in [`references/version-surface-map.md`](references/version-surface-map.md) and run `cargo update -w`.
3. Create `docs/release/vX.Y.Z-pr.md` and add a superseded note at the top of the prior `docs/release/v*-pr.md`.
4. Commit changes (`git commit -m "release: prepare VarDictRust vX.Y.Z"`, zero `Co-Authored-By` trailers) and run the full gate suite:
   ```bash
   cargo fmt --all -- --check
   cargo clippy --all-targets --all-features -- -D warnings
   RUST_TEST_THREADS=4 cargo test --workspace
   cargo build --release
   ./target/release/vardictrust --version

   # Run bundled consistency checker across all 17 pre-tag surfaces + 16 evidence JSONs
   python3 .agents/skills/vardictrust-release/scripts/check_release_consistency.py --version X.Y.Z --phase pre-tag

   bash harness/engine_fingerprint.sh
   bash harness/engine_fingerprint.sh --perf
   bash harness/run_self_tests.sh
   python3 harness/prepare_docs_site.py && mkdocs build --strict && rm -rf .docs_staged site
   cargo publish -p vardictrust-core --dry-run
   cargo package --list -p vardictrust
   bash harness/pr_check.sh
   bash harness/acceptance.sh
   git checkout -- benchmarks/evidence/extended-coverage-tier2.json benchmarks/evidence/multimodal-slices-tier2.json
   gitleaks detect --source . -v
   ```
5. Push `release/vX.Y.Z`, open PR completing `.github/pull_request_template.md` with **all engine checkboxes checked `[x]`**, resolve any inline review comments, and squash-merge to `main`.

### Phase 2: Tag & Automated CI Release (`--phase post-tag`)
```bash
git checkout main && git pull --ff-only origin main
git tag -a vX.Y.Z -m "VarDictRust vX.Y.Z"
git push origin vX.Y.Z

# Monitor .github/workflows/release.yml (plan -> build-local-artifacts -> build-global-artifacts -> host -> attest-release -> publish-crates)
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
bash harness/trigger_rtd.sh build vardictrust latest
```

---

## Rationalization Counter-Table

| Excuse / Rationalization | Reality / Mandatory Action |
| :--- | :--- |
| "I'll bump `Cargo.toml` and `CHANGELOG.md` only, like `v1.0.2` did." | `v1.0.2` left 15 files at `1.0.1` (`crates/vardictrust/Cargo.toml`, `CITATION.cff`, `README.md`, etc.). Always run `check_release_consistency.py --version X.Y.Z --phase pre-tag`. |
| "I'll compute `meta.yaml`'s SHA-256 from a local `git archive` before pushing the `vX.Y.Z` tag." | GitHub's generated source tarball can differ from local `git archive` metadata. Only run `prepare_bioconda_pr.sh X.Y.Z` in Phase 3A after `vX.Y.Z` is tagged on GitHub. |
| "I'll run `nf-core modules containers create` as soon as the Bioconda PR merges." | Seqera Wave requires both `linux-64` and `linux-aarch64` packages to be indexed on `anaconda.org`. `prepare_nfcore_pr.sh` verifies both architectures before invoking Wave. |

## Red Flags — STOP Immediately If You Catch Yourself

- Pushing `vX.Y.Z` tag before the `release/vX.Y.Z` PR has been squash-merged to `main`.
- Running `cargo publish -p vardictrust --dry-run` before `vardictrust-core = X.Y.Z` is published on `crates.io` (use `cargo publish -p vardictrust-core --dry-run && cargo package --list -p vardictrust` in Phase 1).
- Editing `packaging/bioconda/recipes/vardictrust/meta.yaml` or `packaging/nf-core/modules/vardictrust/tests/main.nf.test.snap` by hand instead of via `prepare_bioconda_pr.sh` and `prepare_nfcore_pr.sh`.
