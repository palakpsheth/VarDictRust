# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.1.0] - 2026-11-25

### Added

- **Clinical Short-Read Workflow Presets (`--preset <wgs|wes|amplicon|ctdna>`), Single-Primer Amplicon (`-sp`), Multi-Allelic VCF (`--vcf-all-alleles`), & `92`-Flag CLI Surface** ([`palakpsheth/VarDictRust#94`](https://github.com/palakpsheth/VarDictRust/issues/94), [PR #152](https://github.com/palakpsheth/VarDictRust/pull/152), [PR #154](https://github.com/palakpsheth/VarDictRust/pull/154), [PR #155](https://github.com/palakpsheth/VarDictRust/pull/155), [PR #157](https://github.com/palakpsheth/VarDictRust/pull/157), [PR #161](https://github.com/palakpsheth/VarDictRust/pull/161), [PR #162](https://github.com/palakpsheth/VarDictRust/pull/162)):
  - Added `--preset <wgs|wes|amplicon|ctdna>` (`WorkflowPreset` and `Config::from_preset` in `crates/vardictrust-core/src/config.rs`), keeping default invocations 100% byte-identical to `VarDictJava v1.8.3` (`--compat=vardictjava` / `--parity`) while stamping `##VarDictRustCompat=off` and `##VarDictRustPreset=<preset>` in VCF headers, setting `0`-based BED coordinates (`-c 1 -S 2 -E 3 -g 4 -z` on non-amplicon BEDs while keeping `-R` 1-based unless `-z` is passed explicitly), preserving UCSC 8-column CDS BED8 `SimpleMode` under `wgs`/`wes`/`ctdna` vs `AmpliconMode` (`-a 10:0.95`) under `--preset amplicon`, deduplicating paired somatic `Q20` `FILTER` tags, and giving explicit CLI flags full override precedence (`harness/verify_presets_public_data.sh`).
  - Added `--vcf-all-alleles` (`--all-alleles`, hidden `vardictcpp` alias `--vcfall`; [`AstraZeneca-NGS/VarDictJava#65`](https://github.com/AstraZeneca-NGS/VarDictJava/issues/65)), completing 18/18 `var2vcf_valid.pl` / `var2vcf_paired.pl` post-filter parity by porting `-A` multi-allelic VCF emission across Simple, Amplicon, and Paired Somatic modes via generic `SplitVarDeduper<P>` cross-tile deduplication and associative disk spill/merge (`SpillableVcfAccumulator`) while keeping `-A <INT>` as `VarDictJava`'s `sv_min_af`.
  - Added `--amplicon-single-primer` / `-sp` for single-primer extension (SPE) panels such as QIAseq V3 and ArcherDX ([`AstraZeneca-NGS/VarDictJava#301`](https://github.com/AstraZeneca-NGS/VarDictJava/issues/301), [`AstraZeneca-NGS/VarDictJava#303`](https://github.com/AstraZeneca-NGS/VarDictJava/pull/303), [`AstraZeneca-NGS/VarDictJava#326`](https://github.com/AstraZeneca-NGS/VarDictJava/issues/326), [`AstraZeneca-NGS/VarDict#191`](https://github.com/AstraZeneca-NGS/VarDict/issues/191)) and hidden double-dash `vardictcpp v2..v4` drop-in CLI aliases (`--vcfall`, `--noend`, `--vcfpass`, `--vcf`; [`MHH-Bioinformatics-Hematology/vardictcpp#4`](https://github.com/MHH-Bioinformatics-Hematology/vardictcpp/pull/4), [`bioconda/bioconda-recipes#68008`](https://github.com/bioconda/bioconda-recipes/pull/68008), [`bioconda/bioconda-recipes#71981`](https://github.com/bioconda/bioconda-recipes/pull/71981)), completing the `92`-flag CLI surface across `7` functional groups.
- **Dual-Mode `--compat=vardictjava|off` Governance & All `32` `CompatProfile` Correctness Switches** ([`palakpsheth/VarDictRust#34`](https://github.com/palakpsheth/VarDictRust/issues/34), [`palakpsheth/VarDictRust#94`](https://github.com/palakpsheth/VarDictRust/issues/94), [`palakpsheth/VarDictRust#95`](https://github.com/palakpsheth/VarDictRust/issues/95), [`palakpsheth/VarDictRust#158`](https://github.com/palakpsheth/VarDictRust/issues/158), [PR #152](https://github.com/palakpsheth/VarDictRust/pull/152), [PR #153](https://github.com/palakpsheth/VarDictRust/pull/153), [PR #156](https://github.com/palakpsheth/VarDictRust/pull/156), [PR #157](https://github.com/palakpsheth/VarDictRust/pull/157), [PR #159](https://github.com/palakpsheth/VarDictRust/pull/159)) — defaults to `COMPAT_VARDICTJAVA_1_8_3` (aliases `vardictjava-1.8.3`, `1.8.3`, `--parity`) for 100% byte-identical `VarDictJava v1.8.3` output while activating all `32` switches together under `--compat=off` (`##VarDictRustCompat=off`, zero per-bug `--fix` flags, validated across `harness/get_wave4_public_slices.sh`, `harness/verify_wave4_public_data.sh`, and `harness/verify_wave5_public_data.sh`):
  - *VCFv4.2 / VCFv4.3 Specification & Validator Compliance (8 switches):* `strict_vcf_ad_number_r` ([`clintval/vartovcf#101`](https://github.com/clintval/vartovcf/pull/101), [`AstraZeneca-NGS/VarDict#75`](https://github.com/AstraZeneca-NGS/VarDict/issues/75), [`AstraZeneca-NGS/VarDict#76`](https://github.com/AstraZeneca-NGS/VarDict/pull/76), [`EBIvariation/vcf-validator#99`](https://github.com/EBIvariation/vcf-validator/issues/99); enforces `Number=R` in single-sample `FORMAT/AD` at `VD == 0` sites by emitting `RD,0` when `ALT != "."` and scalar `RD` when `ALT == "."`), `canonical_vcf_gt` ([`clintval/vartovcf#99`](https://github.com/clintval/vartovcf/pull/99), [`clintval/vartovcf#54`](https://github.com/clintval/vartovcf/issues/54), [`AstraZeneca-NGS/VarDictJava#344`](https://github.com/AstraZeneca-NGS/VarDictJava/issues/344); canonicalizes unphased `1/0` to `0/1` and promotes rescued carrier `0/0` with `VD > 0` to `0/1`), `omit_nonsv_vcf_end` ([`clintval/vartovcf#105`](https://github.com/clintval/vartovcf/pull/105); restricts `INFO/END` to structural variants), `recompute_vartype_after_adj_complex` ([`clintval/vartovcf#109`](https://github.com/clintval/vartovcf/pull/109), [`clintval/vartovcf#147`](https://github.com/clintval/vartovcf/pull/147); recomputes `INFO/TYPE` after 5'-prefix stripping on complex alleles), `sanitize_iupac_alt_alleles` ([`AstraZeneca-NGS/VarDictJava#267`](https://github.com/AstraZeneca-NGS/VarDictJava/issues/267), [`AstraZeneca-NGS/VarDictJava#287`](https://github.com/AstraZeneca-NGS/VarDictJava/issues/287), [`AstraZeneca-NGS/VarDict#137`](https://github.com/AstraZeneca-NGS/VarDict/issues/137), [`bcbio/bcbio-nextgen#3321`](https://github.com/bcbio/bcbio-nextgen/issues/3321); reduces non-`ACGTN` IUPAC ambiguity codes in `ALT` outside symbolic `<...>` tags and suppresses non-pileup `REF == ALT` rows), `sanitize_sample_info` ([`AstraZeneca-NGS/VarDictJava#116`](https://github.com/AstraZeneca-NGS/VarDictJava/issues/116); bijective `%XX` encoding for `INFO/SAMPLE`), `canonical_sv_vcf_anchoring` ([`AstraZeneca-NGS/VarDictJava#398`](https://github.com/AstraZeneca-NGS/VarDictJava/issues/398), [`AstraZeneca-NGS/VarDictJava#225`](https://github.com/AstraZeneca-NGS/VarDictJava/issues/225), [`AstraZeneca-NGS/VarDictJava#229`](https://github.com/AstraZeneca-NGS/VarDictJava/pull/229), [`AstraZeneca-NGS/VarDict#111`](https://github.com/AstraZeneca-NGS/VarDict/pull/111), [`AstraZeneca-NGS/VarDictJava#382`](https://github.com/AstraZeneca-NGS/VarDictJava/issues/382), [`clintval/vartovcf#192`](https://github.com/clintval/vartovcf/pull/192), [`clintval/vartovcf#107`](https://github.com/clintval/vartovcf/pull/107), [`clintval/vartovcf#108`](https://github.com/clintval/vartovcf/pull/108), [`samtools/hts-specs`](https://github.com/samtools/hts-specs); VCFv4.2 1-bp upstream reference anchors, `##ALT=<ID=DEL|DUP|INV>`, signed `SVLEN`, and `END;SVTYPE;SVLEN`), and `standard_dup_symbolic_alt` ([`AstraZeneca-NGS/VarDictJava#342`](https://github.com/AstraZeneca-NGS/VarDictJava/issues/342); normalizes `< 100 bp` negative `<dup-N>` alleles to `<DUP>` with unified homozygous span accounting).
  - *Depth, Allele Count, Strand & Subtracted-Mean Invariants (5 switches):* `reconcile_vd_dp` ([`AstraZeneca-NGS/VarDictJava#92`](https://github.com/AstraZeneca-NGS/VarDictJava/issues/92), [`AstraZeneca-NGS/VarDictJava#386`](https://github.com/AstraZeneca-NGS/VarDictJava/issues/386), [`AstraZeneca-NGS/VarDictJava#400`](https://github.com/AstraZeneca-NGS/VarDictJava/issues/400); enforces `0 <= VD <= DP`, `0.0 <= AF <= 1.0`, `0.0 <= ADJAF <= 1.0`, `0 <= HICNT <= HICOV`, and `0.0 <= HIAF <= 1.0`), `reconcile_rd_vd_dp` ([`AstraZeneca-NGS/VarDict#168`](https://github.com/AstraZeneca-NGS/VarDict/issues/168), [`AstraZeneca-NGS/VarDictJava#323`](https://github.com/AstraZeneca-NGS/VarDictJava/issues/323); enforces `DP >= RD + VD`), `reconcile_strand_and_hq_sums` ([`AstraZeneca-NGS/VarDictJava#316`](https://github.com/AstraZeneca-NGS/VarDictJava/issues/316), [`clintval/vartovcf#120`](https://github.com/clintval/vartovcf/pull/120), [`clintval/vartovcf#142`](https://github.com/clintval/vartovcf/pull/142); enforces `ALD[0] + ALD[1] == VD` and `0 <= HICNT <= VD`), `nonnegative_subtracted_means` ([`bcbio/bcbio-nextgen#1704`](https://github.com/bcbio/bcbio-nextgen/issues/1704); clamps `PMEAN`/`QMEAN`/`MQ`/`NM`/`QUAL >= 0`), and `live_duprate` (`u64`-widened duplicate-read accounting across untiled and split-tile runs).
  - *Somatic, Amplicon & Multi-Contig Calling (4 switches):* `somatic_continue_after_bad_var` ([`joachimwolff/VarDictJava#4`](https://github.com/joachimwolff/VarDictJava/pull/4), [`clintval/vartovcf#139`](https://github.com/clintval/vartovcf/pull/139), [`clintval/vartovcf#149`](https://github.com/clintval/vartovcf/pull/149), [`AstraZeneca-NGS/VarDictJava#376`](https://github.com/AstraZeneca-NGS/VarDictJava/issues/376); continues evaluating lower-AF somatic candidates at multi-allelic loci when a top-ranked artifact fails `is_good_var`), `paired_somatic_amplicon` ([`AstraZeneca-NGS/VarDictJava#378`](https://github.com/AstraZeneca-NGS/VarDictJava/issues/378), [`AstraZeneca-NGS/VarDictJava#155`](https://github.com/AstraZeneca-NGS/VarDictJava/pull/155), [`AstraZeneca-NGS/VarDict#109`](https://github.com/AstraZeneca-NGS/VarDict/pull/109), [`MHH-Bioinformatics-Hematology/vardictcpp#1`](https://github.com/MHH-Bioinformatics-Hematology/vardictcpp/pull/1), [`MHH-Bioinformatics-Hematology/vardictcpp#3`](https://github.com/MHH-Bioinformatics-Hematology/vardictcpp/pull/3); paired `-b "tumor|normal" -a` calling with `GDAMP;TLAMP;NCAMP;AMPFLAG` and `FILTER=AMPBIAS`), `preserve_mixed_sex_contigs` (preserves distinct `X`/`chrX`, `Y`/`chrY`, and `MT`/`chrM` contigs), and `coordinate_sorted_vars` (`UB-3` / [`AstraZeneca-NGS/VarDictJava#282`](https://github.com/AstraZeneca-NGS/VarDictJava/issues/282); deterministic ascending genomic coordinate iteration).
  - *CIGAR, Indel Realignment & Adapter Chimera Hardening (7 switches):* `cigar_strip_outer_hard_clips` ([`AstraZeneca-NGS/VarDictJava#341`](https://github.com/AstraZeneca-NGS/VarDictJava/issues/341), [`AstraZeneca-NGS/VarDictJava#103`](https://github.com/AstraZeneca-NGS/VarDictJava/issues/103), [`Kennedy-Lab-UW/VarDictJava@9a0beb0d`](https://github.com/Kennedy-Lab-UW/VarDictJava/commit/9a0beb0d34), [`hedera-dx/VarDictJava@17d1cb51`](https://github.com/hedera-dx/VarDictJava/commit/17d1cb5197); strips outer `H`/`P` CIGAR operators without advancing `alignment_start` and treats interior `P` as non-consuming), `cigar_homopolymer_merge` ([`AstraZeneca-NGS/VarDictJava#402`](https://github.com/AstraZeneca-NGS/VarDictJava/pull/402); conserves read length across homopolymer multi-indels), `eqx_cigar_read_length` (full SAM `=`/`X` equivalence with `M` for `minimap2 --eqx` / `pbmm2`), `adaptor_all_consensus_dirs` ([`AstraZeneca-NGS/VarDictJava#240`](https://github.com/AstraZeneca-NGS/VarDictJava/issues/240), [`AstraZeneca-NGS/VarDictJava#241`](https://github.com/AstraZeneca-NGS/VarDictJava/pull/241), [`AstraZeneca-NGS/VarDict#119`](https://github.com/AstraZeneca-NGS/VarDict/pull/119), [`PolinaBevad/VarDictJava@92c62cf6`](https://github.com/PolinaBevad/VarDictJava/commit/92c62cf6c3520b877c394a4c9fd5fe0d371cd972); screens 5' and 3' Illumina/TruSeq adapters across all 20 soft-clip consensus sites), `shift3_3prime_alleles` ([`AstraZeneca-NGS/VarDictJava#292`](https://github.com/AstraZeneca-NGS/VarDictJava/issues/292); FASTA-matched `REF`/`ALT` sequences after `-3` 3' shifting), `exact_35_overlap`, and `exact_cigar_dlen_match`.
  - *Structural Variants, Region Bounds & IEEE Half-Even Rounding (8 switches):* `sv_respect_min_len` ([`AstraZeneca-NGS/VarDictJava#328`](https://github.com/AstraZeneca-NGS/VarDictJava/issues/328), [`AstraZeneca-NGS/VarDict#157`](https://github.com/AstraZeneca-NGS/VarDict/issues/157); enforces `-L <INT>` by materializing bounded `<= 32,768 bp` explicit sequence alleles for split-refined indels `< -L` and filtering discordant-only SVs `< -L`), `preserve_discordant_svs` (`UB-5`; retains discordant-only single-sample `<DUP>`/`<INV>`/`<DEL>` SVs unless `--min-splitreads` is explicitly passed), `mark_dup_sv_inclusive_end` (`UB-6`), `sv_lazy_mate_pileup`, `strict_region_bounds` ([`AstraZeneca-NGS/VarDictJava#369`](https://github.com/AstraZeneca-NGS/VarDictJava/issues/369)), `suppress_pileup_sentinels` (`UB-1`/`UB-2`), `pileup_preserve_shifted_ref` ([`AstraZeneca-NGS/VarDictJava#330`](https://github.com/AstraZeneca-NGS/VarDictJava/issues/330)), and `exact_binary_half_even`.
- **Cleanroom Pure-`MIT` Hypergeometric Engine, 8-Target Fuzzing, Zero-Allow Lint Ratchet, & `807`-Test Suite** ([`palakpsheth/VarDictRust#92`](https://github.com/palakpsheth/VarDictRust/issues/92), [PR #144](https://github.com/palakpsheth/VarDictRust/pull/144), [PR #149](https://github.com/palakpsheth/VarDictRust/pull/149), [PR #150](https://github.com/palakpsheth/VarDictRust/pull/150), [PR #156](https://github.com/palakpsheth/VarDictRust/pull/156), [PR #161](https://github.com/palakpsheth/VarDictRust/pull/161)):
  - Added cleanroom `MIT` `stats::hypergeom` saddlepoint binomial and hypergeometric routines (Catherine Loader 2000, NIST DLMF §5.11.1, and extended-precision `f80` summation; `0` `.var`/VCF diffs across `15,747` locked Fisher fixture rows and the full corpus), restoring pure `MIT` SPDX licensing across `vardictrust-core`, `vardictrust`, `CITATION.cff`, Bioconda (`linux-64`, `linux-aarch64`, `osx-64`, `osx-arm64` via `packaging/bioconda/prepare_bioconda_pr.sh`), and `nf-core/modules` (`packaging/nf-core/prepare_nfcore_pr.sh`).
  - Enforced the zero-allow 11-lint (`unwrap_used`, `expect_used`, `panic`, `unreachable`, `todo`, `unimplemented`, `dbg_macro`, `print_stdout`, `print_stderr`, `string_slice`, `indexing_slicing`) and `unsafe_code` ratchet across both crates (`harness/lint_ratchet.sh --scan-only`), expanded `cargo-fuzz` to `8` targets orchestrated by `harness/run_deep_fuzz.sh` (`34,396,380` ASan executions / `9,900` CPU-seconds with `0` crashes), and grew the automated Cargo workspace suite from `691` to `807` tests (`580` `vardictrust-core` lib + `110` `corpus_parity` + `42` `engine_behavior` + `2` `amplicon_stream_mem` + `3` `cram_stack` + `35` CLI bin + `35` CLI integration) and Python harness tests from `150` to `156` across `37` diagnostic scripts.

### Changed

- **Multi-Thread Skewed-Panel Scheduling & Bounded CRAM/BAM Memory Optimizations** ([PR #149](https://github.com/palakpsheth/VarDictRust/pull/149), [PR #150](https://github.com/palakpsheth/VarDictRust/pull/150), [PR #151](https://github.com/palakpsheth/VarDictRust/pull/151)), preserving 100% byte-identical `VarDictJava v1.8.3` output (`24/24` profiled VCFs byte-identical):
  - Decoupled weighted worker concurrency (`max_in_flight = 2 * num_threads`) from a `64 MiB` byte-weighted pending-output buffer (`FlowState`) and `32,768 bp` batch-span caps in `regions::batch_tiles`, eliminating head-of-line thread starvation on skewed targeted panels; added tile-scoped `CramContigPinGuard`, `1 MiB` chunked FASTA loading, and a `160 MiB` (`4 × 40 MiB`) byte-tracked LRU `CRAM_CONTAINER_CACHE`, cutting `HG008` somatic CRAM peak RSS by `-471.5 MB` and `hg008-sub` DHAT heap bytes by an additional `-9.64%` (`-289.9 MB`).
  - Added owned `Cow<'_, Vars>` post-processing fast-paths, streaming `AmpliconCoverSweep`, pre-prune `ToVarsBuilder::drain_vars_into` elision, watermarked 2-bit `SeedTable12`/`SeedTable17` k-mer indices (`-49.0%` index bytes), and regex-free `CigarModifier` cursor scanners (`-100%` runtime regex allocations), cutting slice DHAT byte churn by `-23.4%` to `-86.2%` and jemalloc exit `inuse_space` by `-51.9%` to `-95.1%`.
- **Packaging, Provenance, Documentation & Manuscript Synchronization** ([PR #142](https://github.com/palakpsheth/VarDictRust/pull/142), [PR #143](https://github.com/palakpsheth/VarDictRust/pull/143), [PR #145](https://github.com/palakpsheth/VarDictRust/pull/145), [PR #146](https://github.com/palakpsheth/VarDictRust/pull/146), [PR #147](https://github.com/palakpsheth/VarDictRust/pull/147), [PR #154](https://github.com/palakpsheth/VarDictRust/pull/154), [PR #155](https://github.com/palakpsheth/VarDictRust/pull/155), [PR #161](https://github.com/palakpsheth/VarDictRust/pull/161), [PR #162](https://github.com/palakpsheth/VarDictRust/pull/162)): synced prior `v1.0.10` Bioconda/[`bioconda/bioconda-recipes#71875`](https://github.com/bioconda/bioconda-recipes/pull/71875) and Nextflow/[`nf-core/modules#13830`](https://github.com/nf-core/modules/pull/13830) baselines, expanded `UPSTREAM.toml` cross-ecosystem lineage, recertified all `16` `benchmarks/evidence/*.json` manifests (`FP = b18477ae0a7b07c1`, `PERF_FP = 5c41be8047133842`, `3.24x` 8-thread `chr1` scaling in `1.26 GiB` flat RSS), and synchronized `README.md`, `ROADMAP.md`, `DIVERGENCES.md`, `docs/`, `zensical.toml`, `wiki/`, and `paper/{manuscript,supplement}.{md,pdf}`.

### Fixed

- **Spanning-SV `.csi` Indexing & Perl `$FREQ` / `--no-end` VCF Parity** ([PR #156](https://github.com/palakpsheth/VarDictRust/pull/156), [PR #159](https://github.com/palakpsheth/VarDictRust/pull/159)): propagated 16 kbp `Indexer::<LinearIndex>::new(14, 6)` `loffset` entries into `Index<BinnedIndex>` on `--csi` outputs so structural variants spanning multiple 64 KiB BGZF blocks are never missed by `tabix`/`bcftools view -r` ([`clintval/vartovcf#161`](https://github.com/clintval/vartovcf/pull/161)); restored verbatim Perl `$FREQ` string parity in `##FILTER=<ID=f...>` headers and `FILTER` tags (`-f 0.0` → `f0.0`, `-f 1e-5` → `f1e-5`; [`AstraZeneca-NGS/VarDict@009e017`](https://github.com/AstraZeneca-NGS/VarDict/commit/009e017), [`MHH-Bioinformatics-Hematology/vardictcpp#4`](https://github.com/MHH-Bioinformatics-Hematology/vardictcpp/pull/4)); decoupled explicit `--no-end`/`--noend` from `-E <INT>`; and keyed `check_release_consistency.py` Bioconda/`nf-core` license checks on each recipe's own pinned version.
- **Additional `VarDictJava v1.8.3` & `var2vcf` Edge-Case Parity + Hostile-Input & CI Hardening** ([PR #148](https://github.com/palakpsheth/VarDictRust/pull/148), [PR #150](https://github.com/palakpsheth/VarDictRust/pull/150), [PR #153](https://github.com/palakpsheth/VarDictRust/pull/153), [PR #156](https://github.com/palakpsheth/VarDictRust/pull/156)): preserved split-tile multi-BAM (`-b b1:b2:b3`) continuity while retaining upfront [`AstraZeneca-NGS/VarDictJava#187`](https://github.com/AstraZeneca-NGS/VarDictJava/issues/187) index checks; removed the non-Java `500 bp` `BamFeed` halo in `process_amplicon_region`; included SAM `=`/`X` in CRISPR `-J`/`-j` span accounting; enabled multi-tile `-i` splice trajectory replay and `-t` boundary dedup; preserved ASCII spaces in `#CHROM` sample headers while neutralizing Unicode line separators; formatted `-D` debug strings via `round_half_even` and propagated somatic `v2.debug`; seeded `realign::partial_pipeline` `max_read_length` across all four call sites while removing the non-Java SV `partial_pipeline_paths` writeback; added fail-fast `.fai` contig/bounds validation ([`AstraZeneca-NGS/VarDictJava#393`](https://github.com/AstraZeneca-NGS/VarDictJava/issues/393)), BED4 `3'_UTR` apostrophe verification ([`AstraZeneca-NGS/VarDict#172`](https://github.com/AstraZeneca-NGS/VarDict/issues/172)), CRAM `SliceHeader` / `InflightCramGuard` epoch bounds, `SIGPIPE`-safe harness matchers, and workflow `timeout-minutes` caps.

## [1.0.10] - 2026-11-18

### Added

- Expanded `cargo-fuzz` coverage (`SeedMap`-populated `fuzz_modify_cigar` and new `fuzz_parse_and_realign` target), added `37` harness diagnostic self-tests (`25` shell tests + `12` Python test suites, including single-pass FIFO-safe `compare_vcf.py` sample-column and `FORMAT` subfield count parity checks), and expanded the automated Cargo workspace test suite from `664` to `691` tests (`491` `vardictrust-core` library unit tests + `110` `corpus_parity` + `29` `engine_behavior` + `2` `amplicon_stream_mem` + `3` `cram_stack` + `29` CLI bin + `27` CLI integration binary tests) (`crates/vardictrust-core/fuzz/`, `harness/compare_vcf.py`, `harness/tests/`) ([PR #135](https://github.com/palakpsheth/VarDictRust/pull/135), [PR #136](https://github.com/palakpsheth/VarDictRust/pull/136), [PR #137](https://github.com/palakpsheth/VarDictRust/pull/137), [PR #138](https://github.com/palakpsheth/VarDictRust/pull/138), [PR #139](https://github.com/palakpsheth/VarDictRust/pull/139), [PR #140](https://github.com/palakpsheth/VarDictRust/pull/140)).

### Changed

- Landed SemVer-compatible zero-diff CPU, heap-allocation, and code-consolidation optimizations across `vardictrust-core` while preserving 100% `VarDictJava v1.8.3` (`0886261`) byte-for-byte parity ([PR #137](https://github.com/palakpsheth/VarDictRust/pull/137), [PR #138](https://github.com/palakpsheth/VarDictRust/pull/138), [PR #140](https://github.com/palakpsheth/VarDictRust/pull/140)):
  - Clipped CIGAR deletion coverage updates (`parse.rs`) to sorted merged `±5,000 bp` loaded reference windows with idempotent `i64` coordinate-compressed `finish_ref_coverage` sweeps for `> 64 bp` spans; flattened multi-tile `max_read_length` pre-scanning (`lib.rs`) into a non-nested `250 kbp` band pool with disjoint `alignment_start` band ownership, `4 MiB` per-band `FallbackSequential` cap, and a lazy `4 MiB` 2-bit `ACGT` + `'N'` side-map 12-mer index on 2nd+ parent-wide misses; and added `EXP_UNDERFLOW_LN_F64 = -746.0` 0-ULP `exp()` short-circuiting and loop-invariant `p3` hoisting in `fisher.rs`.
  - Replaced `swar_mismatch_count_u64` software `count_ones()` with a 3-op branchless multiply-shift horizontal byte sum, tightened `find_35_match` (`Iterator::position` candidate scanning + exact O(1) `mm0` IEEE-754 `0.1` ratio pre-rejects, `2.60×` faster), added `find_conseq_slice` eliminating all `Sclip` deep `BTreeMap`/`IndexMap` clones across `realignlgdel`/`find_del_paths`/`find_dup_disc_paths`, replaced quadratic `tmp3` scans in `realignlgins30` with binary-searched window slices, indexed `SOFTP2SV` upfront via `SclipSlot` in `find_sv`, and unified `find_mm5`/`find_mm3` and `.var`/genotype/regex formatting over stack-allocated buffers (`StackBuf`).
  - Reduced `JavaHashMap<i32, Vars>` to `96 B` (`Option<Box<V>>` cold fallback, `Option<Box<Vec<K>>>` `NewKeyLog`, `active_tree_bins: u32`, ASCII `java_string_hash`/`java_compare` fast-paths, `into_entries`), deleted the private `StrNode` Red-Black tree in `java_hashmap_key_order` in favor of `JavaHashMap<&str, usize>`, added `build_aligned_variants_pruned` to skip Phase-2 `collect_reference_variants` on doomed sub-threshold sites, consolidated `.var` serializers via `TableModeStyle` across Simple `36/38`, Amplicon `38/40`, and Somatic `55/61`, removed duplicate per-tile `BTreeMap` passes and unused private `SAM_READER_CACHE`/`ADAPTOR_MAP_CACHE`, added `i64` $O(1)$ pre-region `Match` skipping / single-base `Match` specialization / fused paired CIGAR scans / `>= 10S` and `>= 12S` `modify_pure_sm_in_place` fast-paths, and added 2-slot `BAM_META_TLS` (`BAM_CACHE_EPOCH`), `<= 4096 B` `BamFeedInner` stashes, `ACGTacgt=` fast-skip in `compute_cram_nm`, `SeedKey::pair_17_12`, and RAII `WorkerRegionGuard`.
- Synced in-tree `packaging/bioconda/recipes/vardictrust/meta.yaml` to `v1.0.9` verified GitHub release tarball SHA-256 `c40c3bf6180144d87ae39a52470b07e68b61f2ea41ef8c990deaf24529025284`, trimmed unused `regex` and `noodles` `bed`/`vcf` features from binary manifests, and retired `env RUST_MIN_STACK=67108864` and the last `# shortcut:` comment in `harness/cram_twin.sh` ([PR #134](https://github.com/palakpsheth/VarDictRust/pull/134), [PR #136](https://github.com/palakpsheth/VarDictRust/pull/136), [PR #140](https://github.com/palakpsheth/VarDictRust/pull/140)).

### Fixed

- Raised `MAX_CACHED_CRAM_CONTAINER_BYTES` in `crates/vardictrust-core/src/bam.rs` from `4 MiB` to `16 MiB` (`should_cache_cram_container` / `fetch_cram_container_in`), added per-key singleflight CRAM container decode coalescing (`InflightCramDecode` / `InflightCramGuard` on `CramContainerCacheShard`), upgraded `CRAM_META_TLS` from 1 slot to 2 slots (`lookup_cram_meta_tls` / `insert_cram_meta_tls`) eliminating `T <-> N` paired somatic TLS thrash, and buffered `prevalidate_cram_container_slices` through a `64 KiB` seek-retaining window (`SeekWindowReader`) so standard `10,000`–`10,240`-read `150–250 bp` Illumina CRAM containers (`~6.3–9.8 MiB`) decode once across concurrent workers and hit `CRAM_CONTAINER_CACHE` while keeping `8 × 4 = 32` decoded-record slots bounded to `<= 512 MiB` and bypassing oversized long-read containers ([PR #136](https://github.com/palakpsheth/VarDictRust/pull/136), [PR #139](https://github.com/palakpsheth/VarDictRust/pull/139)).
- Hardened hostile-input `.bai`/`.csi` index loading and text SAM CIGAR decoding uncovered by pre-release libFuzzer sweeps: added zero-allocation structural header and remaining-byte count pre-validation (`MAX_BAI_BYTES = 512 MiB`, `MAX_CSI_BYTES = 256 MiB`, `MAX_CSI_DECOMPRESSED_BYTES = 256 MiB`, `n_ref`/`n_bin`/`n_chunk`/`n_intv`/`l_aux` bounds) in `read_validated_bai_index` and `read_validated_csi_index` before delegating to `noodles`, and terminated non-advancing `cigar.iter()` `Err` loops immediately in `SamRecord::fill_from_alignment_record` and `intersects_record` (`crates/vardictrust-core/src/bam.rs`).
- Hardened CLI `-th`/`--threads` positional-BED normalization (`require_equals = true` + free-positional counting preserving Commons-CLI first-wins duplicate `-th`), stale `.tbi`/`.csi` sidecar cleanup and `0555` read-only preflight validation when overwriting BGZF `--output` paths without `--tbi`/`--csi`, `--output` parent-directory spill staging in `SpillableVcfAccumulator` with `std::env::temp_dir()` fallback, and `FileStamp`-validated process-wide `OnceLock` `.fai`/`.gzi` metadata sharing across worker threads (`crates/vardictrust/src/cli.rs`, `src/main.rs`, `crates/vardictrust-core/src/lib.rs`, `src/bam.rs`) ([PR #136](https://github.com/palakpsheth/VarDictRust/pull/136), [PR #138](https://github.com/palakpsheth/VarDictRust/pull/138)).
- Hardened harness BSD/macOS portability, comparator validation, and CI supply-chain pins: added `execvp`-safe `portable.sh --timeout`, `resolve_data_root`, and 16-path `find_gitleaks` while removing GNU-isms across 11 scripts; hardened `.github/workflows/release.yml` (`!announcement_is_prerelease`, step-scoped `CARGO_REGISTRY_TOKEN`, idempotent `vardictrust-core` + `vardictrust` publish, removed empty `announce`); bumped `setup-conda-env` to `conda-env-v3-` with portable CLI pins and pinned `playwright==1.62.0` in `docs.yml`; hardened `gcp_paper_runner.sh` (`--stage` allowlist + immutable SHA-256-verified installers); and added SHA-256 verification on `concordance_multimodal.sh` truth slices and `run_stage_6.sh` downloads ([PR #135](https://github.com/palakpsheth/VarDictRust/pull/135)).

## [1.0.9] - 2026-11-16

### Added

- Ported OpenJDK 8 `java.util.HashMap.TreeNode` Red-Black treeification (`treeify`, `balanceInsertion`, `rotateLeft`, `rotateRight`, `moveRootToFront`, `getTreeNode`, single-pass `probe_tree_slot` / `putTreeVal`), deletion rebalancing (`removeTreeNode` entry-topology untreeification, `balanceDeletion`, and 2-child interior link swapping), and `resize()` `TreeNode.split` (`UNTREEIFY_THRESHOLD = 6`) into `JavaHashMap<K, V>` behind a cold `Option<Box<Vec<TreeLink>>>` sidecar (`size_of::<Node<i32, i32>>() == 12`, `0` extra per-node bytes across `< 8`-collision biological cohorts), sealing `JavaHash` with UTF-16 `String` ordering (`encode_utf16().cmp(...)`) and signed `i32` spread-hash comparisons (`CWE-407` `2^16`-stride collision DoS closure; expanded the automated Cargo workspace test suite from `657` to `664` tests: `467` `vardictrust-core` library unit tests + `110` `corpus_parity` + `29` `engine_behavior` + `2` `amplicon_stream_mem` + `3` `cram_stack` + `28` CLI bin + `25` CLI integration binary tests) (`crates/vardictrust-core/src/javamap.rs`, `src/tovars.rs`, `benches/engine.rs`, `ROADMAP.md`) ([Issue #116](https://github.com/palakpsheth/VarDictRust/issues/116) `SEC-17b`, [PR #131](https://github.com/palakpsheth/VarDictRust/pull/131)).

### Changed

- Migrated the ReadTheDocs (`vardictrust.readthedocs.io`) and strict CI documentation platform from `mkdocs-material` (`mkdocs.yml`) to **Zensical** (`zensical.toml`, `zensical>=0.1.1,<0.2.0`), adding MathJax v3 configuration (`docs/javascripts/mathjax.js`), fence-aware Mermaid and GFM staging normalization in `harness/prepare_docs_site.py`, space-safe Zensical pre-commit/pre-push strict gates, pure-stdlib TOML `nav` validation in `check_release_consistency.py`, and 4-page Playwright Chromium DOM/MathJax/Mermaid/SVG UI verification (`zensical.toml`, `.readthedocs.yaml`, `.github/workflows/docs.yml`, `.githooks/{pre-commit,pre-push}`, `harness/prepare_docs_site.py`, `harness/tests/test_docs_ui.py`) ([Issue #123](https://github.com/palakpsheth/VarDictRust/issues/123), [PR #122](https://github.com/palakpsheth/VarDictRust/pull/122), [PR #129](https://github.com/palakpsheth/VarDictRust/pull/129)).
- Updated workspace and CI dependencies: `noodles` (`0.118 -> 0.119`, including `noodles-bam 0.98.0`, `noodles-cram 0.102.0`, `noodles-sam 0.93.0`, and `noodles-vcf 0.96.0`), `flate2` (`1.1.10 -> 1.1.11`), `lzma-rust2` (`0.22.0 -> 0.24.0`), `log` (`0.4.34 -> 0.4.35`), `pymdown-extensions` (`>=13.0`), transitive Cargo lockfile dependencies, and pinned GitHub Actions (`step-security/harden-runner`, `taiki-e/install-action`, `github/codeql-action/upload-sarif`) (`Cargo.toml`, `Cargo.lock`, `crates/vardictrust-core/fuzz/Cargo.toml`, `.github/workflows/*.yml`) ([PR #110](https://github.com/palakpsheth/VarDictRust/pull/110), [PR #111](https://github.com/palakpsheth/VarDictRust/pull/111), [PR #112](https://github.com/palakpsheth/VarDictRust/pull/112), [PR #115](https://github.com/palakpsheth/VarDictRust/pull/115), [PR #117](https://github.com/palakpsheth/VarDictRust/pull/117), [PR #124](https://github.com/palakpsheth/VarDictRust/pull/124)).

### Fixed

- Hardened release supply-chain bootstrap, dataset SHA-256 verification, tag ruleset governance, and deterministic complexity tests: replaced pipe-to-shell and unverified third-party toolchain installers in `.github/workflows/release.yml` with cached 4-architecture `cargo-dist v0.33.0` SHA-256 archive verification, `macos-15` `aarch64-apple-darwin` runner migration, SLSA v1.0 `resolvedDependencies[0]` alignment, `workflow_dispatch` recovery, and 6-job `step-security/harden-runner` audit; upgraded `harness/assets.tsv` to 6 columns with 40 verified SHA-256 digests, atomic `stage_verified_sha256`, and unconditional disk rehashing across `harness/fetch_data.sh` and slice stagers; enforced lockout-safe `refs/tags/v*` ruleset `24867690` in `harness/configure_github_repo.py`; and replaced sub-second wall-clock test assertions with deterministic step and probe bounds (`.github/workflows/release.yml`, `harness/assets.tsv`, `harness/fetch_data.sh`, `harness/lib/portable.sh`, `harness/configure_github_repo.py`, `crates/vardictrust-core/src/`) ([Issue #116](https://github.com/palakpsheth/VarDictRust/issues/116) `SEC-05` / `SEC-14` / `SEC-12b`, [PR #119](https://github.com/palakpsheth/VarDictRust/pull/119), [PR #120](https://github.com/palakpsheth/VarDictRust/pull/120), [PR #121](https://github.com/palakpsheth/VarDictRust/pull/121), [PR #130](https://github.com/palakpsheth/VarDictRust/pull/130), [PR #131](https://github.com/palakpsheth/VarDictRust/pull/131)).

## [1.0.8] - 2026-11-14

### Added

- Native direct file output, BGZF compression, Tabix/CSI indexing, run telemetry, and explicit parity pinning: `--output <PATH>` (streaming `.vcf`, `.var`, and `-i` splice counts directly to plain files or multithreaded BGZF `.gz`/`.bgz` via `noodles::bgzf::io::MultithreadedWriter` with `File::sync_all` verification), `--tbi` / `--index` pure-Rust Tabix (`.vcf.gz.tbi`) and `--csi` pure-Rust CSI (`.vcf.gz.csi`) index generation (`noodles::tabix` and `noodles::csi`, honoring `INFO/END` coordinate spans for symbolic `<DEL>`, `<DUP>`, `<INV>`, and `<INS>` records, with an actionable `--csi` hint when `--tbi` encounters coordinates exceeding Tabix's `512 Mbp` limit), `--no-end` (`var2vcf_valid.pl -E` `INFO/END` suppression), `--stats-json <PATH>` pipeline telemetry, `--parity` and `--compat vardictjava-1.8.3` CLI forward-compatibility pinning, and `UPSTREAM.toml` reference provenance manifest (`crates/vardictrust/src/cli.rs`, `src/main.rs`, `crates/vardictrust-core/src/lib.rs`, `src/output/vcf.rs`) ([Issue #93](https://github.com/palakpsheth/VarDictRust/issues/93), [Issue #94](https://github.com/palakpsheth/VarDictRust/issues/94), [PR #100](https://github.com/palakpsheth/VarDictRust/pull/100), [PR #102](https://github.com/palakpsheth/VarDictRust/pull/102), [PR #104](https://github.com/palakpsheth/VarDictRust/pull/104)).
- Commons-CLI single-dash long-option compatibility (`-fisher`, `-adaptor`, `-chimeric`, `-deldupvar`, `-nosv`, `-splice`, `-header`, `-dedup`, `-debug`, `-verbose`, `-amplicon`, `-DP`, `-default-printer`, `-defaultPrinter`, `-mfreq`, `-nmfreq`, `-insert-size`, `-insert-std`, `-ref-extension`, `-crispr`, `-?`/`--?`), `--default-printer`/`--insert-size`/`--insert-std`/`--ref-extension`/`--crispr` long aliases, and `-C` and `-Z` / `-downsample` / `--downsample` explicit `VardictError::Unsupported` reporting (`crates/vardictrust/src/cli.rs`) ([PR #104](https://github.com/palakpsheth/VarDictRust/pull/104)).
- Additional upstream VarDictJava CLI parity, short-option cluster normalization, thread-pool environment fallback, and workspace package metadata: `-T` / `--trim` (`-trim`) per-base pileup read trimming (`trim_bases_after` / `is_trim_at_opt_t_bases`), `-n <REGEXP>` sample-name regex extraction, `-VS` / `--VS` `SILENT`/`LENIENT` validation stringency acceptance (`-VS STRICT` fails fast with `VardictError::Unsupported`), clustered boolean-prefix short options (`-pz`, `-pz0`, `-uk`, `-ua`, `-pR`), `RAYON_NUM_THREADS` fallback when `-th` / `--threads` is omitted, and MSRV `rust-version = "1.91.0"` + `crates.io` metadata inheritance across `Cargo.toml`, `crates/vardictrust/Cargo.toml`, and `crates/vardictrust-core/Cargo.toml` (`crates/vardictrust/src/cli.rs`, `crates/vardictrust-core/src/config.rs`, `src/cigar/parse.rs`) ([PR #105](https://github.com/palakpsheth/VarDictRust/pull/105)).
- Cross-tile `.var` intermediate boundary deduplication via `SplitVarDeduper` (`SPLIT_VAR_RETENTION_WINDOW_BP = 200_000`) across split tiles in `SimpleMode` and `SomaticMode`, preserving within-tile multiplicity and `-p` shifted-singleton position-0 companion triples, verified by multi-tile `.var`/`.vcf`/SV/splice conservation suites (`crates/vardictrust-core/src/lib.rs`, `tests/engine_behavior.rs`) ([Issue #91](https://github.com/palakpsheth/VarDictRust/issues/91), [PR #99](https://github.com/palakpsheth/VarDictRust/pull/99), [PR #100](https://github.com/palakpsheth/VarDictRust/pull/100), [PR #101](https://github.com/palakpsheth/VarDictRust/pull/101), [PR #104](https://github.com/palakpsheth/VarDictRust/pull/104)).
- Certified 92-chromosome per-contig parity suite (`docs/validation/chr-parity-v107.md`, `benchmarks/chr-parity-v107.json`: `92 / 92 MATCH`, `0` diffs across 5 cohorts including full 22-autosome `GIAB HG002` WGS `8,111,171` variants), strict Markdown relative-link and MkDocs build audits, clinical/engineering SWOT ecosystem review (`docs/ecosystem-survey.md`), Bioconda and Nextflow `nf-core/modules` `v1.0.7` packaging sync (with `MIT AND GPL-2.0-or-later` block-sequence `meta.yml`), 29 zero-conda / 7 oracle CI self-test wiring, 5-target release performance budget guard, and live `--verify-cargo-list` release-consistency validation ([PR #86](https://github.com/palakpsheth/VarDictRust/pull/86), [PR #87](https://github.com/palakpsheth/VarDictRust/pull/87), [PR #88](https://github.com/palakpsheth/VarDictRust/pull/88), [PR #89](https://github.com/palakpsheth/VarDictRust/pull/89), [PR #97](https://github.com/palakpsheth/VarDictRust/pull/97)).
- Expanded the automated Cargo workspace test suite from `519` to `657` tests (`460` `vardictrust-core` library unit tests + `110` `corpus_parity` + `29` `engine_behavior` + `2` `amplicon_stream_mem` + `3` `cram_stack` + `28` CLI bin + `25` CLI integration binary tests) ([PR #104](https://github.com/palakpsheth/VarDictRust/pull/104), [PR #105](https://github.com/palakpsheth/VarDictRust/pull/105), [PR #106](https://github.com/palakpsheth/VarDictRust/pull/106), [PR #107](https://github.com/palakpsheth/VarDictRust/pull/107)).
- Consolidated `var2vcf_valid.pl` / `var2vcf_paired.pl` VCF post-filter CLI flags (`--pass-only`, `--paired-strict`, `--vcf-pval <FLOAT>`, `--min-splitreads <INT>`, `--no-pstd-filter`; porting `var2vcf_valid.pl` `-S`/`-T`/`-P 0` and `var2vcf_paired.pl` `-S`/`-M`/`-P` with verbatim `P$PVAL` token preservation and programmatic `resolve_effective_config` bounds validation) plus GNU-only single-dash typo guards failing fast with actionable `"did you mean '--...'?"` errors and tightened `split_bool_prefix_cluster` tail checks (`crates/vardictrust-core/src/config.rs`, `src/lib.rs`, `src/output/vcf.rs`, `crates/vardictrust/src/cli.rs`) ([PR #106](https://github.com/palakpsheth/VarDictRust/pull/106)).
- Promoted `-n <REGEXP>` and `-T` / `--trim` to `corpus` parity status (`0` `unexercised` CLI flags remaining in `README.md`) via zero-subprocess in-process `corpus_parity` coverage alongside `--fisher` VCF raw-cell formatting assertions, and closed the resolved `--fisher` VCF raw-cell divergence in `DIVERGENCES.md` (`crates/vardictrust-core/tests/corpus_parity.rs`, `golden/PROVENANCE`, `README.md`, `DIVERGENCES.md`) ([PR #106](https://github.com/palakpsheth/VarDictRust/pull/106)).
- Completed all 11 remaining live `var2vcf_valid.pl` / `var2vcf_paired.pl` VCF threshold, genotype, and cluster-filter CLI flags (`89` total CLI flags across 7 functional groups together with `--completions` and `--manpage`): `--vcf-af <FLOAT>` (`-f`), `--vcf-min-depth <INT>` (`-d`), `--vcf-min-vd <INT>` (`-v`), `--vcf-qmean <FLOAT>` (`-q`), `--vcf-mapq <FLOAT>` (`-Q`), `--vcf-pmean <FLOAT>` (`-p`), `--vcf-sn <FLOAT>` (`-o`), `--vcf-max-nm <FLOAT>` (`-m`), `--vcf-gt-freq <FLOAT>` (`-F`, validated in `[0, 1]`), `--vcf-msi <FLOAT>` (`-I`), and `--vcf-cluster-dist <INT>` (`-c`), modeled via `VcfThreshold<T>` with explicit GNU `--vcf-*` `defined()` (`Option::Some`) zero semantics (honoring `0` on explicit overrides while preserving legacy bare `-f 0` `0.02` single-sample fallback) and `SingleRecordStreamer` + 4-state `SingleVcfEvent` bounded-memory `Cluster<N>bp` sliding-window buffering across tiled and disk-spilled single-sample/amplicon streams (`crates/vardictrust-core/src/config.rs`, `src/lib.rs`, `src/output/vcf.rs`, `crates/vardictrust/src/cli.rs`) ([PR #107](https://github.com/palakpsheth/VarDictRust/pull/107)).
- Added zero-dependency `--completions <SHELL>` (`bash`, `zsh`, `fish`) and `--manpage` (`vardictrust(1)` `roff`) generators backed by `clap::CommandFactory` introspection, reorganized CLI `--help` (`-H` / `--help`) into 7 `help_heading` sections matching `README.md`, added cross-compile-safe completion and manpage installation in `packaging/bioconda/recipes/vardictrust/build.sh`, authored the [`docs/migration.md`](docs/migration.md) Rosetta Stone guide from `VarDictJava` + R + Perl (`var2vcf_valid.pl` / `var2vcf_paired.pl`) pipelines, recorded `[advisory_comparators.vardictcpp]` (`version = "3"`, `bioconda_build = "hc52dbad_0"`) in `UPSTREAM.toml`, and expanded `check_release_consistency.py` structural guards (`crates/vardictrust/src/completions_man.rs`, `src/cli.rs`, `docs/migration.md`, `packaging/bioconda/recipes/vardictrust/build.sh`, `UPSTREAM.toml`, `.agents/skills/vardictrust-release/scripts/check_release_consistency.py`) ([PR #108](https://github.com/palakpsheth/VarDictRust/pull/108)).
- Added end-to-end GRCh38 HLA/alt-contig colon-name support (e.g., `HLA-A*01:01:01:01:100-200[:GENE]`) across `-R` / `--region` and BED calling via right-to-left known-contig prefix matching in `regions::parse_region_flag` and zero-allocation `pub(crate) bam::build_noodles_region` (replacing string-formatted `"{}:{}-{}"` region reparsing), alongside MkDocs 1.6 top-level `not_in_nav` validation, Nextflow `nf-core` CRAM/CSI/ReadTheDocs packaging metadata, and MSRV `1.91.0` / `DIVERGENCES.md` release-consistency guards (`crates/vardictrust-core/src/regions.rs`, `src/bam.rs`, `src/reference.rs`, `mkdocs.yml`, `packaging/`, `harness/tests/`, `.agents/skills/vardictrust-release/`) ([PR #107](https://github.com/palakpsheth/VarDictRust/pull/107)).

### Changed

- Landed zero-output-diff Tier-1 performance and memory optimizations across CIGAR parsing, BAM/CRAM I/O, local realignment, Fisher exact statistics, and VCF/table serialization ([Issue #90](https://github.com/palakpsheth/VarDictRust/issues/90), [PR #98](https://github.com/palakpsheth/VarDictRust/pull/98), [PR #99](https://github.com/palakpsheth/VarDictRust/pull/99), [PR #100](https://github.com/palakpsheth/VarDictRust/pull/100), [PR #101](https://github.com/palakpsheth/VarDictRust/pull/101), [PR #102](https://github.com/palakpsheth/VarDictRust/pull/102), [PR #103](https://github.com/palakpsheth/VarDictRust/pull/103), [PR #104](https://github.com/palakpsheth/VarDictRust/pull/104), [PR #105](https://github.com/palakpsheth/VarDictRust/pull/105), [PR #106](https://github.com/palakpsheth/VarDictRust/pull/106)):
  - Added single-`Match` (`150M`) and pure `S`/`M` in-place `modify_cigar` fast-paths, zero-copy `seq`/`qual` read slices, `modify_cigar_preserving_buffers` scratch buffer capacity restoration, hoisted `adjusted_nm`, stack `ByteSet` (`[u64; 4]`), and direct byte-slice CIGAR scanners in `crates/vardictrust-core/src/cigar/modify.rs` and `src/cigar/parse.rs`.
  - Parallelized multi-tile `scan_multi_tile_max_read_lengths` pre-scans across `run_parallel_for_bam` and added an automatic 64 MiB worker stack floor for `.cram` inputs (`CRAM_WORKER_THREAD_STACK_FLOOR_BYTES`) (`crates/vardictrust-core/src/lib.rs`, `src/runner.rs`).
  - Expanded `BGZF_BLOCK_CACHE` from 8 to 32 shards (`BGZF_CACHE_SHARDS = 32`), added LUT/vectorized BAM sequence decode, switched `-t` duplicate tracking to typed `FxHashSet<DupKey>`, cached per-handle `SharedBamMeta` contig tables, decoded `BamFeed` CIGAR operations in a single pass, validated CRAM container spans lock-free via `Header::reference_sequences()`, shared `Arc<[u8]>` `Reference` spans (`RegionSpanSlot`), inlined `Reference::get_base` primary-region lookup, packed `[u64; 3]` `SeedKey` entries, and pre-sized `SeedMap` hash tables (`crates/vardictrust-core/src/bam.rs`, `src/reference.rs`).
  - Added inline single-allele `VariationEntries` storage with `JavaHashMap` new-key delta logs, Somatic `RefSiteMode` / `CompactRefSite` elision, and `O(W_mate log W_mate)` `sv::partial_pipeline` delta-log merging; eliminated `FisherExact::support` and `sum_elements` heap allocations via algebraic `support_len()`, 1-point fast return, `f80_sum_iter`, 4,096-entry thread-local LRU memoization, and stack `odds_ratio_buf` / `write_odds_ratio_into` formatting; removed redundant `Variant` and `Vars::var_desc_to_variant` clones via `Vars::variant_by_desc`, `ZERO_VARIANT`, `StackBuf`, `Cow<[u8]>` splice/marker helpers, single-winner `BTreeMap` bucketing, and single-buffer VCF formatting (`crates/vardictrust-core/src/stats/fisher.rs`, `src/stats/f80.rs`, `src/tovars.rs`, `src/variations/mod.rs`, `src/javamap.rs`, `src/realign.rs`, `src/sv.rs`, `src/postprocess/`, `src/output/`).
  - Refactored `PairedVcfRecord` to a single contiguous `String` buffer with `usize` byte-span `FILTER` offsets (`3 -> 1` heap `String` allocation per somatic VCF candidate, `u64` spill serialization in `SomaticVcfCandidate`, and zero-copy `PairedRecordStreamer` emission) (`crates/vardictrust-core/src/output/vcf.rs`, `src/lib.rs`) ([PR #106](https://github.com/palakpsheth/VarDictRust/pull/106)).
  - Replaced `perl_num_str` heap `String` allocations with stack-formatted `write_perl_num_str` across single-sample and paired VCF `FILTER`/`INFO` emission, short-circuited single-sample/amplicon `FisherExact` and `uniroot` strand-bias calculation in `build_single_vcf_event` under `--pass-only` when non-Fisher filters already fail, formatted worker candidates through a per-tile `String::with_capacity(512)` scratch buffer retaining exact-sized `String::from(scratch.as_str())` lines in `SingleVcfRecord`, added an `IndexMap::get_mut(chr)` fast-path in `SpillableVcfAccumulator::insert` avoiding per-candidate `chr.to_string()` allocations, reused pre-allocated `Variant.bias` capacity across `postprocess/{simple,amplicon,somatic}.rs`, and eliminated avoidable per-call heap allocations in `cigar/parse.rs`, `realign.rs`, `sv.rs`, and `tovars.rs` (`crates/vardictrust-core/src/output/vcf.rs`, `src/lib.rs`, `src/postprocess/`, `src/cigar/parse.rs`, `src/realign.rs`, `src/sv.rs`, `src/tovars.rs`) ([PR #107](https://github.com/palakpsheth/VarDictRust/pull/107)).

### Fixed

- Added `StagedFileTarget` atomic `.tmp.<pid>.<seq>` staging and `Drop` residue cleanup across `--output`, `--tbi`, `--csi`, and `--stats-json` (preserving live existing targets on early abort while falling back directly for character devices/FIFOs such as `/dev/null` or read-only parent directories), alongside input/output file collision guards refusing output paths that resolve to any input reference, `.fai`, BED, BAM/CRAM, or `.bai`/`.csi`/`.crai` index (`crates/vardictrust-core/src/lib.rs`, `crates/vardictrust/src/main.rs`) ([Issue #93](https://github.com/palakpsheth/VarDictRust/issues/93), [PR #104](https://github.com/palakpsheth/VarDictRust/pull/104), [PR #105](https://github.com/palakpsheth/VarDictRust/pull/105)).
- Eliminated `stderr` `BrokenPipe` (`EPIPE`) panics across `-y` verbose logging, `--DP ERR` debug printing, and profiling warnings via `verbose_eprintln!`, hardened hostile CIGAR/alignment/CRISPR arithmetic (`TLEN = i32::MIN`, `268435455D`, `i32::MAX` & `i64` CRISPR/SV/realignment coordinate saturation, missing `*` qualities, `-x` extension saturation, and `cov_cap` bounds) against signed overflow or OOM allocation, fixed paired-VCF omitted-`-N` sample extraction (`bam1|bam2`) and `#CHROM` sample-name whitespace normalization, enforced value-aware `UNSUPPORTED_FLAGS` and 1-based `-c/-S/-E/-s/-e` BED column validation without `-R`, emitted `-i` splice junction rows in deterministic sorted key order, added actionable 5/6-column BED column hints (`-c 1 -S 2 -E 3 -z`), and synchronized `vardictrust --version` license text to `MIT AND GPL-2.0-or-later` (`crates/vardictrust/src/cli.rs`, `src/main.rs`, `crates/vardictrust-core/src/config.rs`, `src/output/vcf.rs`, `src/cigar/`, `src/sv.rs`, `src/realign.rs`, `src/tovars.rs`, `src/bam.rs`, `src/regions.rs`, `src/lib.rs`) ([Issue #91](https://github.com/palakpsheth/VarDictRust/issues/91), [Issue #93](https://github.com/palakpsheth/VarDictRust/issues/93), [PR #98](https://github.com/palakpsheth/VarDictRust/pull/98), [PR #104](https://github.com/palakpsheth/VarDictRust/pull/104), [PR #105](https://github.com/palakpsheth/VarDictRust/pull/105)).
- Hardened remaining production arithmetic and hostile-input bounds (`AGENTS.md §9`): rejected hostile BAM/CRAM CIGAR reference-span overflows with `io::ErrorKind::InvalidData` across `bam_alignment_end`, `BamFeed`, and `validate_cram_span` and enforced unconditional `i32::MAX` coordinate/tag bounds in `RecordPreprocessor` (`bam.rs`), guarded slice-end wrap in `util::subseq` and sequence offset/length casts in `reference.rs`, widened iteration/exponent bounds and negative-input `NaN` guards in `stats/{zeroin,f80}.rs` and `stats/rmath/stirlerr.rs`, replaced unchecked arithmetic and `u32 -> i32` threshold casts with checked/saturating and `i64`-widened math across `cigar/parse.rs`, `realign.rs`, `sv.rs`, `tovars.rs`, `variations/`, `postprocess/somatic.rs`, and `output/table.rs`, and aligned `DIVERGENCES.md` with `##source=VarDict_v1.8.2` byte-parity and certified `3.14x` `scale.json` metrics (`crates/vardictrust-core/src/bam.rs`, `src/util.rs`, `src/reference.rs`, `src/stats/zeroin.rs`, `src/stats/f80.rs`, `src/stats/rmath/stirlerr.rs`, `src/cigar/parse.rs`, `src/realign.rs`, `src/sv.rs`, `src/tovars.rs`, `src/variations/`, `src/postprocess/somatic.rs`, `src/output/table.rs`, `DIVERGENCES.md`) ([PR #107](https://github.com/palakpsheth/VarDictRust/pull/107)).
- Completed comprehensive pre-`v1.0.8` adversarial security hardening across input parsing, VCF/stderr output sanitization, temp-file/symlink safety, and CI/supply-chain gates: enforced `.fai`/`.gzi`/`.crai`/CRAM container & slice-header allocation caps and out-of-bounds landmark guards before `noodles` container/slice reads (`reference.rs`, `bam.rs`, `output/vcf.rs`), rejected control/newline/NUL injection across `-G`, `-R`, `-N`, `.fai` contigs, and BED records before emitting `##fileformat=VCFv4.2` while preserving SAM `=` (`Seq::Eq`) bases and single-BAM colon paths (`chr17:100-200.bam`) during SV calling (`regions.rs`, `cigar/parse.rs`, `output/vcf.rs`, `sv.rs`, `lib.rs`), neutralized C0/C1/OSC terminal escapes on `stderr` CLI diagnostics and closed profiling/temp-file TOCTOU and `0555` read-only-parent symlink overwrite vectors via `O_CREAT | O_EXCL` (`0600`/`0700`), `O_NOFOLLOW`, and `DirectTempStaged` / `DirectPinned` staging (`main.rs`, `profile.rs`, `lib.rs`), and hardened `wiki.yml` (`persist-credentials: false`, `::add-mask::`, job-scoped `contents: write`), `release.yml` (`origin/main` ancestor tag guard), `deny.toml` (`yanked = "deny"`, `wildcards = "deny"`, `required-git-spec = "rev"`), `dependabot.yml` (`cooldown` `default-days: 7`), and `packaging/nf-core/modules/vardictrust/main.nf` shell quoting ([PR #114](https://github.com/palakpsheth/VarDictRust/pull/114)).

## [1.0.7] - 2026-11-01

### Added

- Bounded 512 MiB reference sequence cache (`DEFAULT_FASTA_CACHE_BUDGET_BYTES`) with active-`CramFeed` `Arc<ContigSlot>` pinning (`Arc::strong_count == 1` eviction), `RegistryDelegatingAdapter` handle pooling (`MAX_POOLED_FASTA_ADAPTERS = 8`), and `IndexedCrai` binary-searched slice seek (`crates/vardictrust-core/src/bam.rs`) ([PR #81](https://github.com/palakpsheth/VarDictRust/pull/81)).
- Multi-threaded amplicon execution across `(chrom, seg_idx)` groups via `runner::run_parallel` with deterministic input-order reassembly (`crates/vardictrust-core/src/lib.rs`) ([PR #82](https://github.com/palakpsheth/VarDictRust/pull/82)).
- Production no-panic lint ratchet (`harness/lint_ratchet.sh`) wired into `harness/pr_check.sh` and `harness/rgate.sh`, plus `x86_64` 80-bit `long double` oracle platform guard (`harness/run_oracle.sh`, `harness/compare.sh`, `harness/golden_store.py`) and POSIX-portable shell helpers (`harness/lib/portable.sh`) ([PR #71](https://github.com/palakpsheth/VarDictRust/pull/71), [PR #72](https://github.com/palakpsheth/VarDictRust/pull/72), [PR #80](https://github.com/palakpsheth/VarDictRust/pull/80)).
- Expanded the automated Cargo workspace test suite from `427` to `519` tests (`350` `vardictrust-core` library unit tests + `110` `corpus_parity` + `26` `engine_behavior` + `2` `amplicon_stream_mem` + `3` `cram_stack` + `19` CLI bin + `9` CLI integration binary tests) and expanded the harness diagnostic self-test suite to 36 scripts (25 shell tests and 11 Python test suites).

### Changed

- Replaced hot-path `VariationMap::get_variation` string formatting/cloning with `indexmap::map::RawEntryMut`, switched `Vars::var_desc_to_variant` to `FxBuildHasher`, boxed `VariationMap::sv: Option<Box<Sv>>`, added an integer `u128` fast path in `types::round_half_even`, and eliminated redundant per-variant heap allocations across `tovars.rs`, `sv.rs`, `util.rs`, and `variations/util.rs` (`crates/vardictrust-core/src/variations/mod.rs`, `src/types.rs`, `src/tovars.rs`, `src/sv.rs`, `src/util.rs`, `src/variations/util.rs`) ([Issue #63](https://github.com/palakpsheth/VarDictRust/issues/63), [PR #82](https://github.com/palakpsheth/VarDictRust/pull/82)).
- Enforced `hg008-cram-twin` (`cram_to_bam_wall_ratio <= 2.50`) and `cram_parity` (`cram_to_bam_wall_ratio <= 1.50`) in `harness/perf_budget.py`, `harness/pr_check.sh`, `harness/recertify.sh`, and `harness/gcp_paper_runner.sh` ([PR #73](https://github.com/palakpsheth/VarDictRust/pull/73), [PR #81](https://github.com/palakpsheth/VarDictRust/pull/81)).

### Fixed

- Resolved comprehensive A–N **VarDictJava v1.8.3** parity and robustness findings across [PR #67](https://github.com/palakpsheth/VarDictRust/pull/67)–[PR #82](https://github.com/palakpsheth/VarDictRust/pull/82):
  - Split `_rjem_malloc_conf` by target OS (`background_thread:true` on Linux vs `background_thread:false,dirty_decay_ms:250,muzzy_decay_ms:0,narenas:8` on non-Linux via `MALLOC_CONF_LINUX_DEFAULT`, `MALLOC_CONF_NON_LINUX_DEFAULT`, `malloc_conf_bytes_for_target` in `crates/vardictrust/src/main.rs`) ([PR #70](https://github.com/palakpsheth/VarDictRust/pull/70)).
  - Replicated `HashMap.putVal` 9-node treeify trigger (`bin_len >= TREEIFY_THRESHOLD` in `JavaHashMap::insert`), `JavaHashMap::get_or_put_with` / `variations::get_or_put_vars`, `ToVarsBuilder::build_aligned_variants` empty-`Vars` insertion parity before `collect_reference_variants`, and `STRICT_VAR_ORDER=1` unsorted `.var` row-order verification (`crates/vardictrust-core/src/javamap.rs`, `src/variations/mod.rs`, `src/tovars.rs`, `src/realign.rs`) ([PR #74](https://github.com/palakpsheth/VarDictRust/pull/74)).
  - Replicated `inc_cnt` / insertion `tp` signed `i32` arithmetic and `realignlgins` signed `max_read_length - 10` bounds, excluded `=`/`X` from `get_soft_clipped_length`/`get_match_insertion_length`/`prepare_sv_structures_for_analysis` (`query_consumed_len` in `bam.rs`), scaled `CigarModifier::modify_cigar`'s iteration bound with CIGAR operator count, and traversed `svffus`/`svrfus` in `is_best_sv_used` in OpenJDK `HashMap<String, List<Sclip>>` bucket order (`java_string_hash`, `java_hashmap_key_order` in `crates/vardictrust-core/src/cigar/modify.rs`, `src/cigar/parse.rs`, `src/realign.rs`, `src/javamap.rs`, `src/sv.rs`, `src/variations/mod.rs`) ([PR #75](https://github.com/palakpsheth/VarDictRust/pull/75)).
  - Replicated `DUPRATE=0` snapshot parity, `Utils.getRoundedValueToPrint` `"N."`/`"-0."` formatting (`push_java_rounded`), `FisherExact::get_odds_ratio` `java_double_to_string` scientific notation under `--fisher`, Apache Commons Math 3 `SaddlePointExpansion` / `HypergeometricDistribution` (`--fisher`), and R 4.x `fround` 5-decimal even-mantissa tie-break (`round_as_r`) (`crates/vardictrust-core/src/output/table.rs`, `src/output/vcf.rs`, `src/stats/fisher.rs`, `src/types.rs`, `src/variations/mod.rs`) ([PR #76](https://github.com/palakpsheth/VarDictRust/pull/76), [PR #77](https://github.com/palakpsheth/VarDictRust/pull/77)).
  - Replicated whole-stream VCF `$hash{chr}{start}` candidate accumulation and 4-decimal `round_half_even` AF winner selection across BED regions, Perl `sub reorder` chromosome ordering (`perl_reorder_chromosomes`), and `> 1 Mbp` multi-tile `max_read_length` pre-scan (`scan_multi_tile_max_read_lengths`) plus `-t` (`config.dedup`) boundary duplicate continuity (`crates/vardictrust-core/src/lib.rs`, `src/output/vcf.rs`, `src/cigar/modify.rs`) ([PR #78](https://github.com/palakpsheth/VarDictRust/pull/78)).
  - Replicated Commons-CLI `createNumber` (`parse_java_number_*`) & first-wins scalar options, `-F` `Integer.decode` & `samfilter_raw != "0"`, adaptor 6-mer `i + 6 < sequence.len()` & `ACGTacgt`-only complement (`realign::fill_adaptor_maps`), BED header/delimiter/4-col auto-detect (`java_split`, `detect_amplicon_bed`, `parse_amp_bed_from_reader` `java_hashmap_key_order`), `resolve_mode` precedence (`#378` / `--DP ERR`), amplicon segment `splice` chaining, and zero-coverage `-p` amplicon null-row formatting (`amplicon_38_into` / `amplicon_40_into`) (`crates/vardictrust/src/cli.rs`, `crates/vardictrust-core/src/realign.rs`, `src/regions.rs`, `src/postprocess/amplicon.rs`, `src/output/table.rs`) ([PR #69](https://github.com/palakpsheth/VarDictRust/pull/69), [PR #79](https://github.com/palakpsheth/VarDictRust/pull/79)).
  - Added truncated BGZF EOF (`UnexpectedEof`) detection, CRAM/FASTA contig span validation (`validate_cram_span`, `fetch_subsequence` past-end error) and 3' chromosome-end overhang padding (`CRAM_REF_OVERHANG_PAD = 65_536`), raw `Vec<u8>` quality byte preservation for `Phred >= 95` (`normalize_cram_base`), `i32` region coordinate bounds (`parse_int_with_commas` `<= i32::MAX`), checked integer conversions across `crates/vardictrust-core`, and the 16 MiB worker stack floor (`WORKER_THREAD_STACK_FLOOR_BYTES` / `run_with_worker_stack` in `crates/vardictrust-core/src/runner.rs`) ([PR #80](https://github.com/palakpsheth/VarDictRust/pull/80)).

## [1.0.6] - 2026-10-29

### Changed

- Increased Rayon worker thread pool stack size to `16 MiB` (`WORKER_THREAD_STACK_BYTES = 16 * 1024 * 1024` in `crates/vardictrust-core/src/runner.rs`) to prevent stack exhaustion on deep CRAM rANS/LZMA container decoding and recursive structural variant realignment paths, and eliminated per-region `HashMap` allocation in `reference::get_reference_with_chr_len` (`crates/vardictrust-core/src/reference.rs`, `src/realign.rs`) ([Issue #60](https://github.com/palakpsheth/VarDictRust/issues/60), [PR #61](https://github.com/palakpsheth/VarDictRust/pull/61), [PR #62](https://github.com/palakpsheth/VarDictRust/pull/62)).
- Expanded the in-process workspace test suite from `416` to `427` tests (`263` `vardictrust-core` library unit tests + `110` `corpus_parity` + `25` `engine_behavior` + `2` `amplicon_stream_mem` + `2` `cram_stack` + `17` CLI bin + `8` CLI integration binary tests) and codified the mandatory multi-skill PR review and inline review thread resolution gate (`AGENTS.md §15.8`, `CONTRIBUTING.md §6.1`, `.github/pull_request_template.md`, `.agents/skills/vardictrust-pr/SKILL.md`, `harness/tests/test_bioconda_recipe.py`) ([PR #61](https://github.com/palakpsheth/VarDictRust/pull/61), [PR #62](https://github.com/palakpsheth/VarDictRust/pull/62)).
- Synchronized in-tree Bioconda recipe metadata (`packaging/bioconda/recipes/vardictrust/meta.yaml`) with the `v1.0.5` source tarball digest ([PR #59](https://github.com/palakpsheth/VarDictRust/pull/59)).
- Simplified `harness/engine_fingerprint.sh` by removing the legacy `normalize_gate_blob()` / `blob_oid()` SemVer normalization workaround so `FP` and `PERF_FP` hash raw `HEAD:<path>` git blob OIDs directly under the Two-Tier Release Recertification Policy (`AGENTS.md §15.4`).
- Recertified all 16 release evidence artifacts in `benchmarks/evidence/*.json` under engine fingerprint `5865a2a0b24a8c17` (`perf_fingerprint` `c490c36317b8f25b`; `harness/acceptance.sh` 8/8 PASS) against live pinned `VarDictJava v1.8.3` and `vardictcpp` oracle executions on an isolated 32-vCPU GCP `n2-standard-32` host (`2 hit / 90 live` golden-store executions).

### Fixed

- Resolved 14 deep-coverage **VarDictJava v1.8.3** parity edge cases across `Simple`, `Somatic`, `Amplicon`, and `Splicing` modes and 7 sequencing chemistries (`Illumina`, `PacBio HiFi`, `ONT`, `Element`, `Ultima`, `MGI`, `Roche SBX`) at both `-f 0.01` and `-f 0` (`-r 2` and `-r 1` noise floor) ([Issue #60](https://github.com/palakpsheth/VarDictRust/issues/60), [PR #61](https://github.com/palakpsheth/VarDictRust/pull/61), [PR #62](https://github.com/palakpsheth/VarDictRust/pull/62)):
  - Cleared per-BAM duplicate tracking while preserving cumulative `total_reads`/`duplicate_reads` (`RecordPreprocessor::reset_for_next_bam`, `Scope.resetDuplicates`), replicated `BAMReader.hasMoreBams` (`if idx > 0 && !any_passed { break; }`) across colon-delimited multi-BAM inputs (`bam1:bam2`), and validated both tumor and normal BAM specifications upfront (`parse_colon_bams` in `crates/vardictrust-core/src/lib.rs`, `src/bam.rs`, `src/realign.rs`, `src/sv.rs`).
  - Seeded BAM2's `CigarParser` and `sv::find_all_sv` (`find_del_disc`) with BAM1's discovered `splice` set and propagated updated `max_read_length` and `splice` state across subsequent `combine_analysis` windows in somatic mode (`crates/vardictrust-core/src/lib.rs`, `src/postprocess/somatic.rs`, matching `SomaticMode.java:104-122` and `SomaticPostProcessModule.java:198-200, 303-311, 398-406`).
  - Updated `max_read_length` prior to the `0x800` supplementary-read early return (`CigarParser.java:181-190`), counted reference `N` bases as mismatches (`nmm++`, `m++`) with unconditional `offset` increment in `M`-after-`I` lookahead (`CigarParser.java:460-476`), and used signed `i32` loop bounds (`for k in 0..element_len`) when rescuing trailing mismatches into a $3'$ soft-clip (`CigarParser.java:1322-1457`) (`crates/vardictrust-core/src/cigar/parse.rs`).
  - Unconditionally cached adaptor-cleared empty consensus sequences (`soft_clip.sequence = Some(final_seq.clone())`) in `find_conseq` (`VariationUtils.java:136-148`, `crates/vardictrust-core/src/variations/util.rs`), delegated `ensure_reference_loaded` to `reference::get_reference_with_chr_len` (`VariationRealigner.java:149, 250, 916, 1020`, `crates/vardictrust-core/src/realign.rs`), preserved `vref.cnt` in $5'$/$3'$ SV inversion branches (`StructuralVariantsProcessor.java:781-795, 922-936`, `crates/vardictrust-core/src/sv.rs`), and allowed signed `next_var.fwd`/`rev` decrements during DUP/INV reference adjustment (`ToVarsBuilder.java:337-340`, `crates/vardictrust-core/src/tovars.rs`).
  - Selected somatic `max_af_variant` on 4-decimal banker-rounded `freq` (`SomaticPostProcessModule.java:118-127, 219-227`), emitted `tlamp = 0, ncamp = 0` on amplicon null pileup rows while persisting `good_variants` across amplicons (`AmpliconPostProcessModule.java:152-246`), canonicalized negative-zero `round_half_even` formatting (`-0` $	o$ `"0"`), applied two-stage `1dp -> 0dp` `HALF_EVEN` rounding to paired `MQ` (`var2vcf_paired.pl`), evaluated VCF filter thresholds on rounded `HIAF`/`SN`/`MSI`/`SBF`, and emitted `;CRISPR=` in single-sample VCF INFO when `-J` is enabled (`crates/vardictrust-core/src/lib.rs`, `src/postprocess/amplicon.rs`, `src/output/table.rs`, `src/output/vcf.rs`).

## [1.0.5] - 2026-10-27

### Added

- First-class in-process `vardictrust-core` test suites (`crates/vardictrust-core/tests/corpus_parity.rs` [110 tests] and `crates/vardictrust-core/tests/engine_behavior.rs` [25 tests], `416` total workspace tests: `252` core lib + `110` `corpus_parity` + `25` `engine_behavior` + `2` `amplicon_stream_mem` + `2` `cram_stack` + `17` CLI bin + `8` CLI integration binary tests) enforced with a `>= 90%` `cargo-llvm-cov` line coverage gate in CI ([PR #49](https://github.com/palakpsheth/VarDictRust/pull/49)).
- Public `vardictrust_core::run_with_writer` API for direct in-process streaming and testing without subprocess spawn overhead ([PR #49](https://github.com/palakpsheth/VarDictRust/pull/49)).
- Supply-chain, CI workflow, and git hook hardening: `.cargo/audit.toml` (`cargo-audit` 0.22+ schema with `deny = ["warnings"]`), `deny.toml` (`[graph] all-features = true` and `[sources]` registry lockdowns), reusable `.github/actions/setup-conda-env` composite action, expanded `.github/dependabot.yml` across all 19 direct crate dependencies plus `/docs` pip dependencies, and offline-safe `cargo deny` / `cargo audit` checks in `.githooks/pre-commit`, `harness/pr_check.sh`, and `harness/acceptance.sh` ([PR #49](https://github.com/palakpsheth/VarDictRust/pull/49), [PR #51](https://github.com/palakpsheth/VarDictRust/pull/51)).
- Automated documentation fingerprint and 16-file SHA-256 manifest synchronizer (`harness/sync_evidence_docs.py`) integrated into `harness/recertify.sh` ([PR #49](https://github.com/palakpsheth/VarDictRust/pull/49)).

### Changed

- Accelerated `VariationRealigner::find_35_match` via SWAR `u64` 8-byte word comparisons (`crates/vardictrust-core/src/realign.rs`), replaced eager CRAM slice vector buffering with streaming cached CRAM container decoding in `crates/vardictrust-core/src/bam.rs`, and borrowed single-allele somatic `Vars` in `crates/vardictrust-core/src/postprocess/somatic.rs` ([Issue #45](https://github.com/palakpsheth/VarDictRust/issues/45), [PR #46](https://github.com/palakpsheth/VarDictRust/pull/46)).
- Added a synchronous single-thread runner fast-path (`num_threads == 1`) in `crates/vardictrust-core/src/runner.rs` and eliminated ~582,600 redundant per-region `stat()` syscalls on BAM metadata cache hits in `crates/vardictrust-core/src/bam.rs` ([PR #49](https://github.com/palakpsheth/VarDictRust/pull/49)).
- Consolidated October 2026 dependency and GitHub Actions updates: bumped `noodles` (`0.117` $\to$ `0.118`), `noodles-bam` (`0.96` $\to$ `0.97`), `noodles-cram` (`0.100` $\to$ `0.101`), `noodles-fasta` (`0.67` $\to$ `0.68`), `noodles-sam` (`0.91` $\to$ `0.92`), `noodles-vcf` (`0.94` $\to$ `0.95`), `rustc-hash` (`1.1.0` $\to$ `2.1.3`), `lzma-rust2` (`0.21.0` $\to$ `0.22.0`), `clap_lex` (`1.1.0` $\to$ `1.1.1`), `libc` (`0.2.186` $\to$ `0.2.191`), `mkdocs` (`1.6.1`), `mkdocs-material` (`9.8.0`), `pymdown-extensions` (`12.3`), and SHA-pinned GitHub Actions (`step-security/harden-runner` `v2.22.2`, `taiki-e/install-action` `v2.87.31`, `actions/cache` `v6.1.0`, `actions/setup-python` `v7.1.0`, `github/codeql-action` `v4.39.1`) ([PR #52](https://github.com/palakpsheth/VarDictRust/pull/52), [PR #57](https://github.com/palakpsheth/VarDictRust/pull/57)).
- Audited and compressed `README.md` with Progressive-JPEG layout, 5-channel installation instructions (including Seqera Wave container digests and Nextflow `nf-core/modules`), and verified CLI flag documentation ([PR #50](https://github.com/palakpsheth/VarDictRust/pull/50)).
- Recertified all 16 release evidence artifacts in `benchmarks/evidence/*.json` under engine fingerprint `d925892846918c2a` (`perf_fingerprint` `b7decb95be6077f1`; `harness/acceptance.sh` 8/8 PASS) with 100% pinned `VarDictJava`/`vardictcpp` golden-store baseline reuse.

### Fixed

- Achieved **100% full-genome parity (`56,642 / 56,642` intervals, `8,111,171` variants, 0 diffs)** against **VarDictJava v1.8.3** across GIAB HG002 `chr1..chr22` ([Issue #47](https://github.com/palakpsheth/VarDictRust/issues/47), [PR #48](https://github.com/palakpsheth/VarDictRust/pull/48)):
  - Preserved soft-clip `unused` state across `find_sv` (`<DEL>`, `<DUP>`, `<INV>`) and synchronized `find_conseq` soft-clip mutations back to the live `Sclip` entry via `sync_sclip_after_conseq` in `crates/vardictrust-core/src/sv.rs`, matching `StructuralVariantsProcessor.java:687-965`.
  - Used `config.extension` (`-x`, default `0`) instead of `config.extension_nuc` (`-X`, default `3`) when computing reference seed fetch bounds in `crates/vardictrust-core/src/reference.rs` (`get_reference`), matching `ReferenceResource.java:85-89`.
- Removed the `paths:` filter from `.github/workflows/differential.yml` so the required `Dual-Engine Differential Parity (Tier-2)` status check executes on every pull request ([PR #50](https://github.com/palakpsheth/VarDictRust/pull/50)).
- Replaced hardcoded `noodles` version string matching in `harness/tests/test_evidence.sh` with a version-agnostic `sed -E` substitution ([PR #57](https://github.com/palakpsheth/VarDictRust/pull/57)).

## [1.0.4] - 2026-10-25

### Fixed

- Achieved **89 / 89 MATCH** (`fp=0, fn=0` across `1,202,444` compared oracle
  `.var`/`.vcf`/splice records and `36` multi-slice/mode suites) against
  **VarDictJava v1.8.3** across the 12-cell full-corpus validation sweep (Tier-P
  `30,267x`–`45,707x` AML panels `SRR15006540` / `SRR8657348` / `SRR15006375` /
  `SRR15006376` / `SRR15006386`, `51.8 GiB` HG002 STAR RNA-seq `-i` splice &
  cDNA, direct whole-chromosome `-R chr1:1-249250621`, full-cohort Amplicon /
  HG002 Deep WGS / HG008 Somatic WGS / SEQC2 WES `690,800` variants, `24`-arm
  WGS/WES MRD `-f 0.001` / `-f 0` / `--fisher` matrix, `30`-arm ultra-deep
  `7,361x`–`173,132x` MRD stress matrix, and `G4`/`E2` CRAM 3.0 twins;
  [Issue #42](https://github.com/palakpsheth/VarDictRust/issues/42),
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

[Unreleased]: https://github.com/palakpsheth/VarDictRust/compare/v1.1.0...HEAD
[1.1.0]: https://github.com/palakpsheth/VarDictRust/compare/v1.0.10...v1.1.0
[1.0.10]: https://github.com/palakpsheth/VarDictRust/compare/v1.0.9...v1.0.10
[1.0.9]: https://github.com/palakpsheth/VarDictRust/compare/v1.0.8...v1.0.9
[1.0.8]: https://github.com/palakpsheth/VarDictRust/compare/v1.0.7...v1.0.8
[1.0.7]: https://github.com/palakpsheth/VarDictRust/compare/v1.0.6...v1.0.7
[1.0.6]: https://github.com/palakpsheth/VarDictRust/compare/v1.0.5...v1.0.6
[1.0.5]: https://github.com/palakpsheth/VarDictRust/compare/v1.0.4...v1.0.5
[1.0.4]: https://github.com/palakpsheth/VarDictRust/compare/v1.0.3...v1.0.4
[1.0.3]: https://github.com/palakpsheth/VarDictRust/compare/v1.0.2...v1.0.3
[1.0.2]: https://github.com/palakpsheth/VarDictRust/compare/v1.0.1...v1.0.2
[1.0.1]: https://github.com/palakpsheth/VarDictRust/compare/v1.0.0...v1.0.1
[1.0.0]: https://github.com/palakpsheth/VarDictRust/releases/tag/v1.0.0

