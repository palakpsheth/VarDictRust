# Research: VarDictRust Rewrite

This document details the architectural decisions, rationales, and alternatives considered for the Rust port of VarDict.

## 1. Crate Layout and Dependencies

- **Decision**: Structuring the project as a Cargo workspace containing two crates:
  - `crates/vardict-core`: Engine library implementing all BAM parsing, pileup logic, realignment, statistics, and formatting. Enforces `#![forbid(unsafe_code)]` and guarantees no panics.
  - `crates/vardict-cli`: The CLI frontend wrapping the core engine and presenting the user-facing interface via `clap`.
  - Use `noodles` for pure-Rust alignment and indexing I/O.
- **Rationale**: Isolates CLI concern from core genomics logic. Forbidding unsafe code and eliminating panics in the core library prevents runtime crashes on malformed data. `noodles` compiles statically without linking to external C libraries, easing cross-compilation.
- **Alternatives Considered**:
  - *HTSlib Rust Bindings (`rust-htslib`)*: Rejected. While high-performance, it binds directly to the C `libhts` library. This introduces complex link dependencies, making static compilation targeting `musl` difficult, and violates the pure-Rust/safe-memory design.
  - *Single Monolithic Crate*: Rejected. Monoliths make it difficult to expose the core variant caller API for integration into other Rust pipelines or services.

## 2. Odds Ratio Conditional MLE Calculation

- **Decision**: Direct translation of `FisherExact.java` and `UnirootZeroIn.java` from the oracle reference codebase to a custom statistics module (`crates/vardict-core/src/stats/fisher.rs`).
- **Rationale**: Replicating R's two-sided Fisher's Exact test and Brent's root-finding conditional Maximum Likelihood Estimation (MLE) of the odds ratio is necessary for exact float equivalence within $\epsilon = 10^{-5}$ tolerance.
- **Alternatives Considered**:
  - *Standard `fishers_exact` crate*: Considered, but it only computes p-values, not the conditional MLE odds ratio. We use `fishers_exact` as a unit-test sanity cross-check, but implement the odds ratio Brent solver manually.

## 3. Somatic STATUS Classification

- **Decision**: Port the Java-based postprocessing code (`SomaticPostProcessModule.java` and related modules) to Rust.
- **Rationale**: Replaces the legacy Perl/R post-filtering scripts (`var2vcf_somatic.pl`, etc.) by embedding their exact heuristics directly into the core engine. Replicating the threshold evaluation ensures somatic classifications (e.g. `StrongSomatic`, `LOH`, `Germline`, `SampleSpecific`) are computed efficiently and identically to the oracle.
- **Alternatives Considered**:
  - *Maintain Perl/R scripts wrapper*: Rejected. Shelling out to Perl or R violates the single statically-compiled binary requirement, requires users to configure runtime interpreters, and degrades performance.

## 4. Region Tiling Strategy

- **Decision**: Partition genome regions into tiles with a halo overlap region (equivalent to maximum read length) to permit local realignments without edge artifacts.
- **Rationale**: Restricts active memory usage to the current tile's pileup size. Ensures that memory remains bounded even at extreme coverage depths (e.g., >100,000x) and prevents OOM crashes.
- **Alternatives Considered**:
  - *Whole-Chromosome Loading*: Rejected. In deep WGS datasets, loading entire chromosome read blocks into memory quickly leads to OOM crashes.
