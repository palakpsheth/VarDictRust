# VarDictRust Version Surface Map (`19` Tracked Surfaces)

Reference companion for [`.agents/skills/vardictrust-release/SKILL.md`](../SKILL.md) and [`scripts/check_release_consistency.py`](../scripts/check_release_consistency.py).

---

## 1. Pre-Tag Surfaces (`--phase pre-tag`, 17 Files Updated in Phase 1 PR)

`normalize_gate_blob()` in `harness/engine_fingerprint.sh` normalizes workspace SemVer in surfaces #1–#3 to `"1.0.0"` before hashing, so bumping `X.Y.Z` preserves `FP` and `PERF_FP`. When feature or bugfix PRs merged since the prior release have advanced `FP` or `PERF_FP`, running `bash harness/recertify.sh` in Phase 1 refreshes the 16 `benchmarks/evidence/*.json` files and invokes `python3 harness/sync_evidence_docs.py` to auto-sync the `Certified Release Engine Fingerprint`, `perf_fingerprint`, and 16-row §13 SHA-256 table across surfaces #10, #12, #14, #15, and #16.

| # | File Path | Exact Version / Metadata Pattern | Purpose |
| :--- | :--- | :--- | :--- |
| 1 | `Cargo.toml` | `[workspace.package]` $\to$ `version = "X.Y.Z"` | Workspace package version inherited by both crates |
| 2 | `crates/vardictrust/Cargo.toml` | `vardictrust-core = { path = "../vardictrust-core", version = "X.Y.Z" }` | Pin required by `cargo publish -p vardictrust` on Crates.io |
| 3 | `Cargo.lock` | `name = "vardictrust"` & `name = "vardictrust-core"` $\to$ `version = "X.Y.Z"` | Updated via `cargo update -w` |
| 4 | `CITATION.cff` | `version: X.Y.Z` and `date-released: "YYYY-MM-DD"` | Zenodo / scholarly citation metadata |
| 5 | `CHANGELOG.md` | `## [Unreleased]`, `## [X.Y.Z] - YYYY-MM-DD`, `[Unreleased]: .../vX.Y.Z...HEAD`, `[X.Y.Z]: .../vPrev...vX.Y.Z` | Keep-a-Changelog release history & diff links |
| 6 | `.github/ISSUE_TEMPLATE/bug_report.yml` | `placeholder: "vardictrust X.Y.Z (commit ...)"` | Issue template version example |
| 7 | `packaging/bioconda/prepare_bioconda_pr.sh` | `VERSION="${1:-X.Y.Z}"` | Default version fallback for Phase 3A Bioconda script |
| 8 | `packaging/nf-core/prepare_nfcore_pr.sh` | `VERSION="${1:-X.Y.Z}"` | Default version fallback for Phase 3B Nextflow script |
| 9 | `harness/run_self_tests.sh` | `# VarDictRust vX.Y.Z — Harness Self-Test Runner` | Header version banner |
| 10 | `docs/release/vX.Y.Z-pr.md` | `# Release vX.Y.Z:` | Release certification dossier (`FP`, `PERF_FP`, test counts, 4-channel checklist) |
| 11 | `README.md` | `--tag vX.Y.Z`, `version = {X.Y.Z}`, and `416` test breakdown (`252 + 110 + 25 + 2 + 2 + 17 + 8`) | Main repository landing page |
| 12 | `DIVERGENCES.md` | `VarDictRust_vX.Y.Z`, `416 PASS` (`252 PASS` core lib), `3.48x speedup, 1.22x flat RSS`, `d925892846918c2a` | Empirical divergence & certification summary |
| 13 | `CONTRIBUTING.md` | `416` workspace tests (`§5.4`), `aarch64-unknown-linux-musl` (`§9.1`), and `.agents/skills/` runbook links | Contributor & release guide |
| 14 | `benchmarks/README.md` | `Certified in \`vX.Y.Z\`` and link to `docs/release/vX.Y.Z-pr.md` | Benchmark scoreboard header |
| 15 | `docs/testing/test-matrix-and-ci-audit.md` | `vX.Y.Z`, `416` workspace tests (`252` core lib + `110` `corpus_parity` + `25` `engine_behavior` + `2` `amplicon_stream_mem` + `2` `cram_stack` + `17` bin + `8` CLI integration), `d925892846918c2a` | Testing matrix & CI audit specification |
| 16 | `docs/validation/omnibus-empirical-validation.md` | `Release \`vX.Y.Z\`` header, `d925892846918c2a`, and §13 16-file SHA-256 manifest | Omnibus empirical validation report |
| 17 | `mkdocs.yml`, `wiki/Home.md`, & `wiki/_Sidebar.md` | `Release vX.Y.Z Certification: docs/release/vX.Y.Z-pr.md`, `vX.Y.Z` links, and `416` workspace tests | ReadTheDocs navigation & GitHub Wiki |

---

## 2. Post-Tag & Downstream Surfaces (Updated in Phases 3A & 3B)

These two surfaces depend on live upstream artifacts and must **not** be hand-edited in Phase 1 (`AGENTS.md` §0b):

| # | File Path | Phase Updated | Exact Pattern | Upstream Dependency |
| :--- | :--- | :--- | :--- | :--- |
| 18 | `packaging/bioconda/recipes/vardictrust/meta.yaml` | **Phase 3A (`--phase post-bioconda`)** | `{% set version = "X.Y.Z" %}` + live `sha256:` + `build: number: 0` | Requires `https://github.com/palakpsheth/VarDictRust/archive/refs/tags/vX.Y.Z.tar.gz` created in Phase 2 |
| 19 | `packaging/nf-core/modules/vardictrust/{environment.yml,main.nf,.conda-lock/*.txt,tests/main.nf.test.snap}` | **Phase 3B (`--phase post-nfcore`)** | `bioconda::vardictrust=X.Y.Z`, Wave container URLs, `.conda-lock/` lockfiles, and `["VARDICTRUST", "vardictrust", "X.Y.Z"]` in `main.nf.test.snap` | Requires `bioconda::vardictrust=X.Y.Z` (`linux-64` + `linux-aarch64`) live on `https://api.anaconda.org/package/bioconda/vardictrust` |
