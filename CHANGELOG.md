# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.0.4] - 2026-10-25

### Fixed

- Achieved **89 / 89 MATCH** (`0` variant diffs across `11,227,154` called variants)
  against **VarDictJava v1.8.3** across the 12-cell full-chromosome, full-exome,
  ultra-deep MRD (`7,000x`–`173,000x`), RNA-seq (`--splice` & `-i`), CRAM 3.0 twin,
  and amplicon panel validation corpus ([Issue #42](https://github.com/palakpsheth/VarDictRust/issues/42),
  [PR #43](https://github.com/palakpsheth/VarDictRust/pull/43)):
  - **Bug A** (`crates/vardictrust-core/src/cigar/parse.rs`, `src/reference.rs`):
    Select the covering secondary `LoadedRegion` (`region_slice_for_span`) on
    `self.region.chr` during structural-variant `partial_pipeline` re-parsing and
    spliced alignments (`CigarParser.java:181–193`, `257`).
  - **Bug B** (`crates/vardictrust-core/src/reference.rs`): Truncate loaded
    reference sequences to `site_end = bases.len() - 17` after seeding
    `PendingSeedChunk` (`ReferenceResource.java:106–112`) so the trailing
    `SEED_1` (`17 bp`) lookahead is excluded from `region_sequences` and
    `reference_sequences` while remaining indexable by overlapping subsequent
    regions.
  - **Bug C** (`crates/vardictrust-core/src/bam.rs`): Normalize CRAM
    reference-matching bases decoded against soft-masked lowercase FASTA contigs
    (`a/c/g/t/n`) to uppercase (`A/C/G/T/N`) in `SamRecord::fill_from_alignment_record`,
    and reconstruct omitted CRAM 3.0 (`store_nm=0`) `NM` edit distance tags via
    `compute_cram_nm` per SAMtags §1.5.
  - **Bug D** (`crates/vardictrust/src/cli.rs`, `crates/vardictrust-core/src/config.rs`,
    `src/output/vcf.rs`): Replicate Perl falsy `$opt_f ? $opt_f : 0.02` fallback
    (`var2vcf_valid.pl:21`) via `valid_freq_threshold(config)` when `-f 0` is
    passed in single-sample mode, while preserving `0.0` for truthy decimal/scientific
    zero strings (`-f 0.0`, `-f 0e0`).
  - **Bug E** (`crates/vardictrust-core/src/postprocess/somatic.rs`): Propagate
    in-place `v2nt.is_noise()` coverage decrements (`total_pos_cov -= pos_coverage`)
    back into `v2.var_desc_to_variant` and `v2.variants` so subsequent tumor
    alleles at multi-allelic somatic loci observe the mutated `v2.variants[0]`
    (`SomaticPostProcessModule.java:172–186`).
  - **Bug F** (`crates/vardictrust-core/src/realign.rs`): Fix `<dupN>` $3'$ flank
    substring offset in `realignins` (`start_idx = sv_flanking`, matching
    `VariationRealigner.java:744–746` `substr(ins3, SVFLANK - ins3.length())`).
  - **Bug G** (`crates/vardictrust-core/src/sv.rs`): Query
    `var_data.ref_coverage.get(&(p5_adj + 1))` after `p5--` in `find_sv` $5'$
    candidate deletion scan (`StructuralVariantsProcessor.java:732, 748`).
  - **Bug H** (`crates/vardictrust-core/src/stats/fisher.rs`, `src/output/table.rs`,
    `src/output/vcf.rs`, `src/lib.rs`): Add `FisherExact::new_java` /
    `new_with_mode(..., config.fisher)` emulating `FisherExact.java:72, 112, 129`
    with 7-decimal `logdc` banker's rounding (`round_half_even(value, 7)`) and
    OpenJDK `DoubleStream.sum()` Kahan-Babuška compensated summation
    (`java_stream_sum`) for ultra-deep `--fisher` parity.
  - **Adaptive Region Tiling** (`crates/vardictrust-core/src/lib.rs`,
    `crates/vardictrust/src/cli.rs`): Execute BED intervals `<= 1 Mbp`
    (`DEFAULT_UNTILED_SPAN_MAX`) as single untiled intervals matching
    VarDictJava while tiling multi-megabase intervals `> 1 Mbp` into `100 kbp`
    (`DEFAULT_TILE_SIZE`) windows to preserve `< 1.5 GiB` bounded WGS memory;
    raise `MAX_TILE_SIZE` to `500 Mbp`.
- Resolved mate-window BAM index path fallback (`crates/vardictrust-core/src/lib.rs`)
  so tiled structural-variant mate queries accept both `foo.bam.bai` and
  `foo.bai` index naming conventions (`resolve_bam_index_path()`), matching the
  primary BAM reader (`bam.rs`).
- Wired `harness/scale_bench.sh` and `harness/run_integration.sh` to cache and
  replay `vardictcpp` invocations via `harness/golden_store.sh`, eliminating
  redundant multi-minute `vardictcpp` re-executions across test and scale loops.
- Added unit and real-data end-to-end regression tests across `crates/vardictrust-core`
  and `crates/vardictrust/tests/flags.rs` (`293` total workspace tests: `241` core
  lib + `4` core integration + `15` CLI unit + `33` CLI integration tests).
- Recertified all 16 release evidence artifacts in `benchmarks/evidence/*.json`
  under engine fingerprint `552856835a1eec01` (`perf_fingerprint`
  `9346a6e5419cc5e2`; `harness/acceptance.sh` 8/8 PASS).

## [1.0.3] - 2026-10-24

### Added

- Opt-in zero-default-overhead profiling subsystem (`crates/vardictrust/src/profile.rs`,
  `crates/vardictrust/tests/profile_smoke.rs`) with Cargo features `profile-cpu`
  (`pprof` + `inferno` SVG flamegraphs), `profile-heap` (`tikv-jemallocator` +
  `jemalloc_pprof`), and `profile-dhat` (`dhat::Alloc` allocation-site churn),
  driven by `harness/profile.sh` and stdlib-only `harness/profile_catalog.py`
  (`PR #35`).
- Golden Baseline Store (`harness/golden_store.sh`, `harness/golden_store.py`,
  `harness/tests/test_golden_store.py`) and turnkey G4 cohort stager
  (`harness/stage_g4.sh`) for hash-verified, host-scoped caching of pinned
  `VarDictJava v1.8.3` and `vardictcpp 3` oracle outputs (`PR #32`, `PR #33`).
- Complete 10-class `vardictcpp` divergence forensic catalog (`CPP-1` through
  `CPP-10`) in `docs/upstream-bugs.md`, cross-linked across `README.md`,
  `benchmarks/README.md`, and `docs/validation/omnibus-empirical-validation.md`
  (`PR #37`).
- In-tree Nextflow `nf-core/modules` packaging (`packaging/nf-core/modules/vardictrust/`,
  `packaging/nf-core/prepare_nfcore_pr.sh`) with multi-arch Seqera Wave
  containers (`linux/amd64`, `linux/arm64`) and `nf-test` snapshots.
- Four repository-specific Agent Skills under `.agents/skills/`
  (`vardictrust-issue`, `vardictrust-pr`, `vardictrust-profiling`, and
  `vardictrust-release` with stdlib-only `check_release_consistency.py`),
  validated in `harness/tests/test_bioconda_recipe.py`.

### Changed

- Zero-allocation `get_variation` / `get_variation_maybe` (`IndexMap`
  `raw_entry_mut_v1` with `FxBuildHasher`, stack UTF-8 key buffer) and
  `FxHashMap` `ref_coverage` + borrowed-key `increment` in `CigarParser`
  (`PR #35`).
- Replaced per-bucket `Vec` allocations in `JavaHashMap` with an OpenJDK-isomorphic flat intrusive node arena (`buckets: Vec<u32>`, `nodes: Vec<Node<K, V>>`, free-list reuse), preserving exact bucket/chain iteration and `treeify_bin` semantics while eliminating inner-bucket heap allocations and bare slice indexing (`PR #35`).
- BAM I/O and record decoding optimization (O-1 + O-3): share parsed `.bai`/`.csi`
  index and SAM header metadata across worker threads via `META_CACHE`, pool up
  to 4 BGZF readers per thread/path via `BAM_READER_POOL` (cleared via RAII
  guard around `run()`), and decode concrete `noodles::bam::Record` fields and
  auxiliary tags (`nM`/`NM`/`SA`/`MC`/`MQ`) in a single pass without per-read
  trait-object allocations (`PR #35`).
- `CigarModifier::modify_cigar`: added dynamic `(has_d, has_i, has_s)` byte
  pre-filters across pre-loop and fixpoint-loop regexes (recomputed on every
  `cigar_str` mutation) and replaced the 5 unconditional post-loop regexes with
  bounds-checked, allocation-free byte-slice matchers (`PR #35`).
- `VariationRealigner::realignlgins30` / `find_conseq` / `find_35_match` (NEW-1):
  hoisted `sc5_used` and integer distance bounds ahead of `BTreeMap` lookups in
  `realignlgins30`, cached empty consensus `Some(Vec::new())` in `find_conseq`,
  and pre-reversed `seq3` once for forward slice comparison in `find_35_match`
  (`PR #35`).
- `ToVarsBuilder` and `realignlgdel` optimization (O-9): replaced upfront
  `non_insertion_variants.clone()` and per-position `VariationMap::clone()` with
  a sparse copy-on-write `FxHashMap` overlay, eliminated per-position `Vec` and
  `Variant` clones and `HashSet<String>` allocations in
  `is_the_same_variation_on_ref`, and guarded `realignlgdel` soft-clip loops on
  live `cnt`/`used`/consensus length before cloning `Sclip` (`PR #35`).
- Reference k-mer seed indexing (`O-4 + O-5`): replaced `HashMap<Vec<u8>, Vec<u32>>`
  with `SeedMap` (`SeedKey` 17-byte inline key + `SeedPositions::One/Many` +
  `FxBuildHasher`) and eliminated redundant `bases.clone()` and `LoadedRegion`
  `String` clones in `get_reference` (`PR #35`).
- Stack-buffered `round_half_even` (`StackBuf<64>` in `types::round_half_even`
  with heap fallback for oversized values, O-12a) and allocation-free integer
  CIGAR deletion suffix check `cigar_contains_dlen_d` in
  `VariationRealigner::no_passing_reads` (O-7, `PR #35`).
- Enabled `noodles-bgzf/libdeflate` for whole-block BGZF deflate decompression
  in BAM/BAI/CSI/VCF.gz readers (O-2, +3.74% median wall-time improvement on
  `hg002` chr20, `PR #35`).
- Configured `tikv-jemallocator` as the default Linux global allocator with
  static `_rjem_malloc_conf` (`background_thread:true,dirty_decay_ms:250,muzzy_decay_ms:0,narenas:8`)
  and guarded `StructuralVariantsProcessor::find_sv` to filter out sub-`minr`
  soft-clip clusters before cloning `Sclip`, cutting peak RSS by 1.7×–4.2×
  (`hg002` `497.1 MiB`, `seqc2` `467.9 MiB`, `hg008` `313.1 MiB`, `chr1`
  249 Mbp `1.42 GiB` with `1.22×` flat RSS scaling, `PR #35`, `PR #37`).
- High-region-count Exome BED optimization (`EBO-1..EBO-3`, Issue #36, `PR #37`):
  added `CachedBgzfReader` with a 1-block thread-local cursor and 8-shard × 16-slot
  LRU process-wide `BGZF_BLOCK_CACHE` in `bam.rs` (`EBO-1`), deferred
  `Reference::seed` 12-mer/17-mer hashing until the first soft-clip realignment
  query via `OnceCell<RawSeedMap>` in `reference.rs` (`EBO-2`), and batched up
  to 16 consecutive small same-chromosome tiles in `regions::batch_tiles` and
  `runner::run_parallel_batched` (`EBO-3`), cutting `seqc2-exome-tier2` wall
  time by **2.51×** (`19m 28.78s` $\to$ `7m 46.34s`, `6.9×` faster than
  `VarDictJava`) and `scale.json` exome BED wall time across 1–8 threads by
  **1.38×–2.38×** (`3m 06.43s` at `-th 8`).
- Recertified all 16 release evidence artifacts in `benchmarks/evidence/*.json`
  under engine fingerprint `27f3667b32ac14f5` (`perf_fingerprint`
  `9f75643cf12aa988`; `harness/acceptance.sh` 8/8 PASS, `PR #37`).

### Fixed

- Per-tile BAM streaming for `chrM` RNA splicing (`-i`), ungated parity-restoring
  `splice-tile-ownership` and `mate-window-ownership` tile ownership rules, and
  lowercase-`n` `allele_key` passthrough (`PR #29`).

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

[Unreleased]: https://github.com/palakpsheth/VarDictRust/compare/v1.0.4...HEAD
[1.0.4]: https://github.com/palakpsheth/VarDictRust/compare/v1.0.3...v1.0.4
[1.0.3]: https://github.com/palakpsheth/VarDictRust/compare/v1.0.2...v1.0.3
[1.0.2]: https://github.com/palakpsheth/VarDictRust/compare/v1.0.1...v1.0.2
[1.0.1]: https://github.com/palakpsheth/VarDictRust/compare/v1.0.0...v1.0.1
[1.0.0]: https://github.com/palakpsheth/VarDictRust/releases/tag/v1.0.0

