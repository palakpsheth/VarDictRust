# VarDictRust Documentation Portal

Welcome to the **VarDictRust** documentation portal.

VarDictRust is an idiomatic Rust rewrite of **VarDictJava v1.8.3** (a somatic and germline variant caller for NGS), consolidating its fragmented Java→R→Perl pipeline into a single, statically-compiled, memory-safe binary that scales to deep whole-genome sequencing (WGS). Built under the **[rewrites.bio](https://rewrites.bio)** manifesto.

---

## Quick Navigation

### 1. Getting Started
- **[README](https://github.com/palakpsheth/VarDictRust#readme)** — Executive overview, performance scoreboards, quickstart commands, and installation.
- **[Installation Guide](https://github.com/palakpsheth/VarDictRust#installation)** — Pre-built static musl binaries, cargo install, and verification.
- **[Contributing Guidelines](https://github.com/palakpsheth/VarDictRust/blob/main/CONTRIBUTING.md)** — Development prerequisites, code style, and test harnesses.

### 2. Empirical Validation & Parity
- **[Omnibus Empirical Validation](https://github.com/palakpsheth/VarDictRust/blob/main/docs/validation/omnibus-empirical-validation.md)** — Comprehensive validation suite across GIAB HG002, SEQC2 HCC1395, and 92 integration fixtures.
- **[NIST v4.2.1 HG002 Benchmark](https://github.com/palakpsheth/VarDictRust/blob/main/docs/validation/nist-v4.2.1-hg002-benchmark.md)** — Precision, recall, and F1 concordance benchmarks evaluated with RTG `vcfeval`.
- **[Extended Modality Validation](https://github.com/palakpsheth/VarDictRust/blob/main/docs/validation/extended-coverage-validation.md)** — RNA-seq, ultra-deep amplicon, WGS, and edge-case parity evidence.
- **[Multi-Modal Multi-Slice Validation](https://github.com/palakpsheth/VarDictRust/blob/main/docs/validation/multimodal-multislice-validation.md)** — Cross-dataset concordance on HG005, NA24385 (b37), and NA12878.

### 3. Architecture & Methodology
- **[Divergences & Bug Tracking](https://github.com/palakpsheth/VarDictRust/blob/main/DIVERGENCES.md)** — Categorized audit of VarDictJava quirks, edge cases, and deterministic HashMap emulation.
- **[Upstream Bugs](https://github.com/palakpsheth/VarDictRust/blob/main/docs/upstream-bugs.md)** — Consolidated catalog of upstream VarDictJava v1.8.3 defects preserved for parity.
- **[Ecosystem Survey](https://github.com/palakpsheth/VarDictRust/blob/main/docs/ecosystem-survey.md)** — Comparative survey of modern NGS variant callers.

### 4. Testing & Quality Assurance
- **[Comprehensive Test Matrix & CI Audit](https://github.com/palakpsheth/VarDictRust/blob/main/docs/testing/test-matrix-and-ci-audit.md)** — 7-tier test architecture, 208 cargo unit/integration tests, and CI workflows.
- **[Release Waivers & Parity Audit](https://github.com/palakpsheth/VarDictRust/blob/main/docs/release/waivers.md)** — Formal justification of intentional behavioral divergences.
- **[v1.0.0 Release PR Specification](https://github.com/palakpsheth/VarDictRust/blob/main/docs/release/v1.0.0-pr.md)** — Release manifest, cryptographic fingerprints, and certification.

### 5. Project Roadmap
- **[ROADMAP](https://github.com/palakpsheth/VarDictRust/blob/main/ROADMAP.md)** — Governance for post-1.0 opt-in `--fix` flags, UMI deduplication, and long-read support.
