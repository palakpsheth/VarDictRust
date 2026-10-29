## Summary of Changes
<!-- Provide a clear, concise description of changes, motivation, and references to tracking issues. -->

## Type of Change
- [ ] Bug fix (non-breaking change restoring upstream parity or robustness)
- [ ] Performance optimization (scaling or memory footprint improvement)
- [ ] Test harness / benchmarking expansion
- [ ] Documentation update
- [ ] Infrastructure / CI workflow update

---

## Parity Validation Receipt ([AGENTS.md §15 Item 5](../AGENTS.md#15-definition-of-done-every-change))
<!-- PRs modifying calling logic or statistical filters MUST provide empirical validation evidence: -->
> **Validation Statement:** `validated vs VarDictJava 1.8.3, dataset <name>, cmd <...>`
```
<!-- Paste comparator output or live runner evidence here (e.g. compare.sh MATCH receipt) -->
```

---

## Mandatory PR Regression Suite Checklist
All checks must pass before PR approval. Run locally via `bash harness/pr_check.sh`:
Docs/infra-only PRs whose diff touches none of `*.rs`, `Cargo.toml`/`Cargo.lock`, `crates/`, `harness/` may mark engine boxes below N/A with a one-line justification instead.
- [ ] **Unit & Workspace Tests ([Level 1](../docs/testing/test-matrix-and-ci-audit.md#level-1)):** `cargo test --workspace` passes with zero failures.
- [ ] **Integration Corpus ([Level 2 / Tier-1](../docs/testing/test-matrix-and-ci-audit.md#level-2)):** `harness/run_integration.sh` reports `PASS: 91  XFAIL: 1  FAIL: 0  SKIP: 0  XPASS: 0` (UB-4 upstream fixture defect; 100% live oracle parity).
- [ ] **Crash Robustness Suite ([Level 2](../docs/testing/test-matrix-and-ci-audit.md#level-2)):** `cargo test -p vardictrust --test robustness` passes (Hostile CIGAR [PR #402](https://github.com/AstraZeneca-NGS/VarDictJava/pull/402), missing BAI [#187](https://github.com/AstraZeneca-NGS/VarDictJava/issues/187)).
- [ ] **Differential Parity ([Level 3 / Tier-2](../docs/testing/test-matrix-and-ci-audit.md#tier-2)):** `harness/compare.sh` and `pr_check.sh` Stage 4 (`micro`, `s1a..s9`, `run_extended_coverage.sh --no-record`, `run_multimodal_slices.sh --no-record`) executed with 0 diffs ($\varepsilon \le 10^{-6}$).
- [ ] **Release Recertification ([AGENTS.md §15 Item 4](../AGENTS.md#15-definition-of-done-every-change)):** If this is a release branch (`release/vX.Y.Z`) or milestone certification PR, run `bash harness/recertify.sh` and `bash harness/acceptance.sh` (`N/A — deferred to release cut per two-tier testing policy` on normal feature/fix PRs).

---

## Mandatory Performance & Memory Verification
Every PR modifying variant calling, pileup, or realignment logic must verify performance:
- [ ] **Performance Budget:** `python3 harness/perf_budget.py hg008-tier2 seqc2-exome-tier2 hg002-tier2 amplicon-tier2` passes.
- [ ] **Runtime Delta:** PR runtime $\le$ Baseline runtime (within 5% tolerance).
- [ ] **Peak RSS Delta:** PR peak RSS $\le$ Baseline peak RSS (no memory regressions).
- [ ] **Bounded Memory:** Verified that memory usage remains $O(1)$ across wide intervals via region tiling.

---

## Security, Hygiene & Cluster Etiquette
- [ ] `cargo fmt --all -- --check` clean.
- [ ] `cargo clippy --all-targets --all-features -- -D warnings` clean (zero warnings).
- [ ] `gitleaks detect --source .` clean (0 leaks).
- [ ] Cluster resource etiquette respected (concurrency capped at $\le 8$ threads or `RUST_TEST_THREADS=4`).
- [ ] Zero `Co-Authored-By` commit trailers per commit policy.
- [ ] Engine fingerprint verified (`bash harness/engine_fingerprint.sh`).
- [ ] Multi-skill review completed and findings fixed (`/ponytail-review`, `/rust-router`, and `/pragmatic-rust-guidelines` per [AGENTS.md §15.8](../AGENTS.md#15-definition-of-done-every-change)).
- [ ] Waited for `copilot-pull-request-reviewer[bot]` review to finish; every review comment refuted or fixed with an inline reply and every review thread resolved (`resolveReviewThread`).
- [ ] Security disclosures follow [SECURITY.md](SECURITY.md) (contact: `palakpsheth@gmail.com`).
