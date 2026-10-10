# Contributing to VarDictRust

Thank you for contributing to **VarDictRust** (`vardictrust`)! VarDictRust is an idiomatic Rust rewrite of [VarDictJava](https://github.com/AstraZeneca-NGS/VarDictJava) v1.8.3, built under the [rewrites.bio](https://rewrites.bio) manifesto.

Our core commitment is to deliver a faster, more robust variant caller that produces output **identical to VarDictJava v1.8.3 on valid inputs ($\varepsilon = 10^{-6}$, 0 variant diffs)**.

---

## 1. Core Engineering Principles

### The Prime Directive (rewrites.bio §5)
Output must match **VarDictJava v1.8.3** on valid inputs:
- **Exact Parity:** The variant set (CHROM, POS, REF, ALT), every FILTER decision, INFO/FORMAT tags, genotypes, variant TYPE, and somatic STATUS must match upstream VarDictJava.
- **Float Tolerance:** P-values, odds ratios, and allele frequencies must match within $\varepsilon = 10^{-6}$, adhering strictly to IEEE 754 half-even banker's rounding (`round_half_even`).
- **Intermediate Column Counts:** Intermediate `.var` tables must strictly adhere to upstream column counts:
  - **Simple Mode:** 36 columns (default) / 38 columns (`--fisher`).
  - **Amplicon Mode:** 38 columns (default) / 40 columns (`--fisher`).
  - **Somatic Mode:** 55 columns (default) / 61 columns (`--fisher`).
- **Bug Emulation Policy (v1.0):** If v1.8.3 has an output quirk, **replicate it** (see [`docs/upstream-bugs.md`](docs/upstream-bugs.md)). Corrections are deferred to v1.1.0+ behind opt-in `--fix <id>` flags per [`ROADMAP.md`](ROADMAP.md). Where v1.8.3 crashes (e.g. `modifyCigar` #341), VarDictRust produces the correct output or fails fast without panicking.

### Integrity — Evidence Before Assertions (§0 & §0b)
- Never claim a task or PR is done without verifiable evidence.
- Never edit or skip failing tests or fixtures in `golden/` or `VarDictJava/testdata/` to turn a mismatch into a match. If a test fails against live Java, your code or understanding is wrong: fix the port, not the test.

### Coordinate Indexing (§6)
- Noodles BAM/SAM positions are **0-based, half-open** (`Pos0`).
- VCF output and VarDict-internal coordinates are **1-based inclusive** (`Pos1`).
- Convert **only** via `.to_one_based()` / `.to_zero_based()`. Coordinate errors are silent and catastrophic.

---

## 2. Standard Pull Request Workflow

Every pull request must pass the automated PR regression and performance gate prior to review.

### Step 1: Run the Automated PR Check
Before submitting a PR, run the local PR check script:
```bash
bash harness/pr_check.sh
```
This script automatically executes:
1. **Hygiene & Unit Tests:** `cargo fmt --check`, `cargo clippy -D warnings`, and `cargo test --workspace` (all 202 tests).
2. **Crash Robustness Suite:** `cargo test -p vardictrust --test robustness` (validates PR #402, missing BAI #187, and issue #341 family).
3. **Tier-1 Upstream Integration Corpus:** `harness/run_integration.sh` (validates all 92 fixtures in `VarDictJava/testdata/integrationtestcases/`, asserting the expected `91 PASS / 1 FAIL` signature where the sole static diff is UB-4).
4. **Tier-2 Micro Differential Parity:** `harness/compare.sh` on `fixtures/micro/` verifying bitwise equivalence against live Java.
5. **Performance Budget:** `harness/perf_budget.py` asserting that runtime and peak RSS memory meet pre-registered performance budgets.
6. **Secret Scanning:** `gitleaks detect --source . -v` verifying 0 credential leaks.

### Step 2: Fill Out the PR Template
Fill out all checkboxes in [`.github/pull_request_template.md`](.github/pull_request_template.md), providing runtime and memory benchmarks for any changes touching variant calling, pileup, or realignment logic.

---

## 3. Commit Guidelines & Host Policies

- **Atomic Commits:** Structure changes into logical atomic commits with task prefixes (e.g., `Txxx: <short description>`).
- **No Co-Authored-By Trailers:** In accordance with repository Host Policy A10, **do not** include `Co-Authored-By` commit trailers. AI disclosures are maintained transparently in `README.md` and `docs/release/v1.0.0-pr.md`.
- **Courtesy Resource Limits:** For long-running or multi-threaded jobs on shared systems, constrain executions to $\le 8$ threads and wrap in `EST_MIN=<min> bash harness/courtesy.sh`.
