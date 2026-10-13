# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Official Bioconda package distribution for `vardictrust` 1.0.0
  ([anaconda.org/bioconda/vardictrust](https://anaconda.org/bioconda/vardictrust))
  supporting `linux-64`, `linux-aarch64`, `osx-arm64`, and `osx-64`
  via `conda install -c bioconda -c conda-forge vardictrust`. Galaxy wrapper distribution
  remains an explicit follow-up.
- Advisory vardictcpp third arm everywhere VarDictRust compares against
  VarDictJava: pinned **vardictcpp 3 `hc52dbad_0` bioconda** (same
  `vardict-oracle` conda env as the Java oracle; pin mirrors the
  `vardict=2019.06.04` pattern in `benchmarks/env/vardict-oracle.yml`).
  Every gate arm records a `CPP: MATCH/MISMATCH/SKIP` line plus
  `compare_cpp.txt`; evidence `extra` gains `cpp_version`, `cpp_head`,
  `cpp_wall_s`, `cpp_maxrss_kb`, `cpp_f1`, `cpp_skip_reason`;
  `CPP_SKIP=1` escapes the arm. The release gate stays **java==rust** —
  cpp never masks a verdict, and `harness/acceptance.sh` fails loudly
  (not silently) when Tier-2 evidence lacks `cpp_version`/`cpp_skip_reason`.
  Known cpp divergences are normalized in `docs/upstream-bugs.md`.

- Three-way CRAM twin (`harness/cram_twin.sh`): rust and cpp each call the
  hg008 slice pair from BAM and CRAM; gate is rust-BAM == rust-CRAM
  (18,426 variants MATCH), cpp twin recorded advisory. The rust-CRAM arm
  runs with `RUST_MIN_STACK=64M`: the default 8M stack overflows on these
  high-depth slices (BAM unaffected) — a tracked engine frame-hog, not a
  parity issue. Rebase note (2026-10-12): pre-registered scale bar
  (`wall(th1)/wall(th8) >= 4.8`) measures **3.31×** on i7-13700H/
  `taskset 0-7` (single-thread got faster: 1984s vs 3436s on Ryzen) while
  the RSS leg passes (1.45× ≤ 1.5×); `scale.json` is recorded FAIL honestly
  pending threshold recalibration per host class. Java `-th 8` is
  OOM-killed on high-depth paired slices on 31G hosts (bisected 1✓/4✓/8✗);
  hg008/seqc2-exome evidences run at `-th 4` with cmds recorded verbatim.


## [1.0.0] - 2026-10-10

First release of VarDictRust: an idiomatic-Rust rewrite of
[VarDictJava](https://github.com/AstraZeneca-NGS/VarDictJava) v1.8.3
(single static binary; no Java/R/Perl required).

### Added

- All four calling modes: Simple (germline), Somatic (paired tumor/normal),
  Amplicon, and Splicing/SV.
- Complete VCFv4.2 output plus the legacy intermediate `.var` table
  (`--emit-intermediate`).
- In-engine Fisher exact test (`--fisher`) — no R dependency.
- Bounded-memory region tiling + rayon parallelism for deep WGS.
- CRAM input support.

### Parity

Output validated against VarDictJava v1.8.3 (commit `0886261`), across 92
integration fixtures, 17 differential benchmarks, and 5 ground-truth cohorts.
Full results, per-slice exceptions, and the divergence ledger live in
[README §Empirical Parity & Validation Status](README.md#empirical-parity--validation-status) and
[`DIVERGENCES.md`](DIVERGENCES.md).

### Fixed

Robustness improvements by design (no valid-input output changes relative to VarDictJava v1.8.3):

- Bounds-checked CIGAR processing fixes the entire VarDictJava `modifyCigar`
  crash family ([VarDictJava#341](https://github.com/AstraZeneca-NGS/VarDictJava/issues/341),
  [VarDictJava#243](https://github.com/AstraZeneca-NGS/VarDictJava/issues/243),
  [VarDictJava#253](https://github.com/AstraZeneca-NGS/VarDictJava/issues/253),
  [VarDictJava#283](https://github.com/AstraZeneca-NGS/VarDictJava/issues/283),
  [VarDictJava#7](https://github.com/AstraZeneca-NGS/VarDictJava/issues/7),
  [VarDictJava#38](https://github.com/AstraZeneca-NGS/VarDictJava/issues/38)).
- No panics in library code: NPE/AIOOBE/div-by-zero classes eliminated by the
  no-panic rule (`#![forbid(unsafe_code)]`, no `unwrap` outside tests).
- Bounded memory at any depth via region tiling addresses the OOM family
  ([VarDictJava#197](https://github.com/AstraZeneca-NGS/VarDictJava/issues/197),
  [VarDictJava#384](https://github.com/AstraZeneca-NGS/VarDictJava/issues/384)).
- Deterministic output regardless of thread count ([VarDictJava#282](https://github.com/AstraZeneca-NGS/VarDictJava/issues/282),
  [VarDictJava#181](https://github.com/AstraZeneca-NGS/VarDictJava/issues/181),
  [VarDictJava#356](https://github.com/AstraZeneca-NGS/VarDictJava/issues/356)).
- Correct `##contig` lengths from the `.fai` ([VarDictJava#320](https://github.com/AstraZeneca-NGS/VarDictJava/issues/320));
  fail-fast on missing BAM index ([VarDictJava#187](https://github.com/AstraZeneca-NGS/VarDictJava/issues/187)).

[Unreleased]: https://github.com/palakpsheth/VarDictRust/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/palakpsheth/VarDictRust/releases/tag/v1.0.0
