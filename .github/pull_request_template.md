## Summary of Changes
<!-- Provide a clear, concise description of changes, motivation, and references to tracking issues. -->

## Type of Change
- [ ] Bug fix (non-breaking change restoring upstream parity or robustness)
- [ ] Performance optimization (scaling or memory footprint improvement)
- [ ] Test harness / benchmarking expansion
- [ ] Documentation update
- [ ] Infrastructure / CI workflow update

---

## Mandatory PR Regression Suite Checklist
All checks must pass before PR approval. Run locally via `bash harness/pr_check.sh`:
- [ ] **Unit & Workspace Tests:** `cargo test --workspace` passes (all 208 workspace tests).
- [ ] **Tier-1 Integration Corpus:** `harness/run_integration.sh` reports `PASS: 91  FAIL: 1  SKIP: 0` (sole expected failure is UB-4 fixture reference defect; 100% live oracle parity).
- [ ] **Robustness Crash Suite:** `cargo test -p vardictrust --test robustness` passes (Issue #341 family, PR #402, missing BAI #187).
- [ ] **Tier-2 Differential Parity:** `harness/compare.sh` executed on affected modalities with 0 diffs ($\varepsilon = 10^{-6}$).

---

## Mandatory Performance & Memory Verification
Every PR modifying variant calling, pileup, or realignment logic must verify performance:
- [ ] **Performance Budget:** `python3 harness/perf_budget.py hg008-tier2 seqc2-exome-tier2 hg002-tier2 amplicon-tier2` passes.
- [ ] **Runtime Delta:** PR runtime $\le$ Baseline runtime (within 5% tolerance).
- [ ] **Peak RSS Delta:** PR peak RSS $\le$ Baseline peak RSS (no memory regressions).
- [ ] **Bounded Memory:** Verified that memory usage remains $O(1)$ across wide intervals.

---

## Security, Hygiene & Compliance
- [ ] `cargo fmt --all -- --check` clean.
- [ ] `cargo clippy --all-targets --all-features -- -D warnings` clean (zero warnings).
- [ ] `gitleaks detect --source .` clean (0 leaks).
- [ ] Zero `Co-Authored-By` commit trailers per commit policy.
- [ ] Engine fingerprint verified (`bash harness/engine_fingerprint.sh`).
- [ ] Security disclosures follow [SECURITY.md](SECURITY.md) (contact: `palakpsheth@gmail.com`).
