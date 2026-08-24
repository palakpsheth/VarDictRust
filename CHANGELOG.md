# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0] — unreleased

First release of VarDictRust: an idiomatic-Rust rewrite of
[VarDictJava](https://github.com/AstraZeneca-NGS/VarDictJava) v1.8.3
(single binary; no Java/R/Perl required).

### Added

- All four calling modes: Simple (germline), Somatic (paired tumor/normal),
  Amplicon, and Splicing/SV.
- Tier-2 VCF output plus the legacy intermediate `.var` table
  (`--emit-intermediate`).
- In-engine Fisher exact test (`--fisher`) — no R dependency.
- Bounded-memory region tiling + rayon parallelism for deep WGS.
- CRAM input support.

### Parity

Output validated against VarDictJava v1.8.3 (commit `0886261`), dataset =
the VarDictJava/testdata integration suite. Full results, per-slice
exceptions, and the divergence ledger live in
[README §Validation status](README.md#validation-status) and
[DIVERGENCES.md](DIVERGENCES.md).

### Fixed (robustness by design — no valid-input output change)

- Bounds-checked CIGAR processing fixes the entire VarDictJava `modifyCigar`
  crash family (#341/#243/#253/#283/#7/#38 and siblings).
- No panics in library code: NPE/AIOOBE/div-by-zero classes eliminated by the
  no-panic rule (`#![forbid(unsafe_code)]`, no `unwrap` outside tests).
- Bounded memory at any depth via region tiling addresses the OOM family
  (#197/#384).
- Deterministic output regardless of thread count (#282/#181/#356).
- Correct `##contig` lengths from the `.fai` (#320); fail-fast on missing BAM
  index (#187).
