# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.0.2] - 2026-10-17

### Fixed

- Bounded-memory amplicon calling: per-base insert expansion replaced with an
  interval-overlap test plus a streaming `-p` path (active-interval sweep),
  eliminating OOMs on huge BED intervals with byte-identical output.

### Added

- Harness gates: pinned cpp `--chunk` (default 50000) recorded to evidence,
  `ORACLE_XMX_GB` heap sizing for the Java oracle, and a gated truvari SV
  wrapper with loud SKIP when truvari is absent.
- All 16 evidence artifacts recertified under the new engine fingerprint
  (acceptance 8/8 PASS).

## [1.0.1] - 2026-10-15

### Added

- Official Bioconda package distribution for `vardictrust`
  ([anaconda.org/bioconda/vardictrust](https://anaconda.org/bioconda/vardictrust))
  supporting `linux-64`, `linux-aarch64`, `osx-arm64`, and `osx-64`
  via `conda install -c bioconda -c conda-forge vardictrust`. Galaxy wrapper distribution
  remains an explicit follow-up.
- Advisory `vardictcpp` third arm everywhere VarDictRust compares against
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
  runs with `RUST_MIN_STACK=64M`: the default 2 MiB rayon worker stack overflows
  inside `noodles-cram 0.100.0` rANS/name-tokenizer decoding (~97 frames; BAM
  unaffected; pinned by `crates/vardictrust-core/tests/cram_stack.rs`), not a
  parity issue. Rebase note (2026-10-12; finalized 2026-10-14 per docs/validation/scale-recalibration-proposal.md sign-off): scale PASSES under the signed host-class bar `bar(P) = 0.6 × min(8, P)` (P=4 on i7-13700H `taskset 0-7` → bar 2.4; RSS leg unchanged ≤ 1.5×); medians across 3 full runs give **3.72×** (th1 1769s/th8 476s) with RSS 1.448× ≤ 1.5×, so `scale.json` records PASS (median 3.72 >= 2.4, RSS 1.448 <= 1.5). RSS headroom is thin (one run reached 1.489×) — variance watch only, no bar change. Java `-th 8` is
  OOM-killed on high-depth paired slices on 31G hosts (bisected 1✓/4✓/8✗);
  hg008/seqc2-exome evidences run at `-th 4` with cmds recorded verbatim.
- Expanded `harness/run_extended_coverage.sh` from 9 to 17 clinical calling
  modalities (`--crispr`, `-u`, `-A`, `--chimeric`, `-k 0`, `--fisher`, `-3`,
  `-O`, `-X`, `-c/-S/-E/-g`, `-M`, `--deldupvar`, `--DP ERR`) and added 19 CLI
  flag integration tests in `crates/vardictrust/tests/flags.rs` plus 2 CRAM
  decoder stack tests in `crates/vardictrust-core/tests/cram_stack.rs`
  (230 total workspace tests).
- Documentation and governance infrastructure: MkDocs Material site
  (`mkdocs.yml`, `harness/prepare_docs_site.py`, `.github/workflows/docs.yml`),
  automated GitHub Wiki sync (`wiki/`, `.github/workflows/wiki.yml`),
  ReadTheDocs API trigger (`harness/trigger_rtd.sh`), Playwright/GFM
  documentation UI test suite (`harness/tests/test_docs_ui.py`), and the
  end-to-end release engineering runbook in `CONTRIBUTING.md §9`.

### Changed

- Canonicalized workspace-only SemVer fields (`Cargo.toml`,
  `crates/vardictrust/Cargo.toml`, `Cargo.lock`) inside
  `harness/engine_fingerprint.sh` (`normalize_gate_blob()`) so pure release
  version bumps preserve the certified engine fingerprint (`388f4a545470d12a`)
  while external dependencies, profiles, features, and sources remain hashed.
- Re-certified all 16 evidence artifacts in `benchmarks/evidence/` under
  fingerprint `388f4a545470d12a` on the 14-core i7-13700H workstation (chip; benchmarks at -th 8, scale pin P=4 physical).

### Fixed

- Fixed multi-sample BAM header extraction on `-b "b1|b2"` in
  `crates/vardictrust/src/main.rs`, added warning when `-N` omits the second
  sample name, and guarded `end == 0` in `StructuralVariantsProcessor::find_inv`.
- Fixed `harness/record_evidence.py` stale `compare_cpp.txt` / `cpp.time` leak on
  `KEEP_DIR` reuse when `extra.json` records a pre-execution skip.

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

[Unreleased]: https://github.com/palakpsheth/VarDictRust/compare/v1.0.1...HEAD
[1.0.1]: https://github.com/palakpsheth/VarDictRust/compare/v1.0.0...v1.0.1
[1.0.0]: https://github.com/palakpsheth/VarDictRust/releases/tag/v1.0.0

