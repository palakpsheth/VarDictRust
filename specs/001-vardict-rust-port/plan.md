# Implementation Plan: VarDictRust Port

**Branch**: `init-scaffold` | **Date**: 2026-06-19 | **Spec**: [specs/001-vardict-rust-port/spec.md](spec.md)

**Input**: Feature specification from `specs/001-vardict-rust-port/spec.md`

## Summary

The VarDictRust Port project is a rewrite of VarDictJava v1.8.3 in idiomatic Rust as a single, statically-compiled binary. The primary goal is to achieve exact output parity (variant set, filters, genotypes, somatic STATUS, intermediate columns) on valid inputs, while eliminating runtime crashes, reducing memory usage via region tiling, and improving performance through Rayon thread work-stealing.

## Technical Context

**Language/Version**: Rust 1.75+

**Primary Dependencies**: `noodles (0.75+)`, `bio`, `block-aligner`, `statrs`, `fishers_exact`, `rayon`, `clap`, `indexmap`, `thiserror`, `anyhow`

**Storage**: Local file I/O (BAM, CRAM, FASTA, BED, VCF)

**Testing**: `cargo test` suite, `harness/compare_var.py` (Tier-1), `harness/compare.sh` (Tier-2)

**Target Platform**: Ubuntu 12c / 24 GB / amd64 (AVX2 SIMD), static musl target `x86_64-unknown-linux-musl`

**Project Type**: Library & CLI Tool

**Performance Goals**: Average execution speed $\ge$ 3x faster than VarDictJava v1.8.3, peak memory $<$ 4GB per thread

**Constraints**: Pure Rust (no C toolchain), no panics in `vardict-core`, Banker's rounding (`round_half_even`), exact intermediate columns (Simple 36/38, Amplicon 38/40, Somatic 55/61)

**Scale/Scope**: Bounded memory at any depth, Rayon work-stealing parallelism, tiled output identical to single-region runs

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Gate | Requirement | Status |
| :--- | :--- | :--- |
| **1. Output Parity** | Match VarDictJava v1.8.3 output exactly on valid inputs | ✅ Pass (All heuristics/formatting preserved) |
| **2. Two-Tier Validation** | Diff-testing of both VCF and intermediate `.var` outputs | ✅ Pass (Tier-1 and Tier-2 comparative scripts set up) |
| **3. Crate Architecture** | Core logic in `vardict-core` (#![forbid(unsafe_code)], no panic) | ✅ Pass (Ensured by library design and compiler flags) |
| **4. Scale-First Tiled Engine** | Genomic tiling with a halo for bounded memory | ✅ Pass (Region tile processing is the core architectural unit) |
| **5. Coordinate & Data Fidelity** | Typed coordinates (`Pos0`/`Pos1`), banker's rounding, IndexMap | ✅ Pass (Integrated into the core types and math libraries) |

All gates pass. No violations of constitution principles exist.

## Project Structure

### Documentation

```text
specs/001-vardict-rust-port/
├── plan.md              # This file
├── research.md          # Decisions, rationales, and alternatives
├── data-model.md        # Core entities and data representations
├── quickstart.md        # Runnable verification scenarios
└── contracts/
    └── cli.md           # CLI parameters and expected inputs/outputs
```

## Complexity Tracking

No violations of the Constitution exist. No additional complexity justifications are required.
