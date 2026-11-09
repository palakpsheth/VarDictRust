# VarDictRust Roadmap

VarDictRust v1.0.0 replicates VarDictJava v1.8.3 output byte-for-byte on valid inputs, including the upstream behaviors documented in `DIVERGENCES.md` and [`docs/upstream-bugs.md`](docs/upstream-bugs.md).

To keep end-user CLI ergonomics simple while guaranteeing that rigorously tested VarDictJava v1.8.3 parity can always be reproduced, post-1.0 evolution follows a **Two-Phase `--parity` / `--modern` SemVer Model** (tracked in `VarDictRust#94`):
- **`v1.0.x` (Current Baseline):** Parity-only by default (no `--modern` or `--parity` CLI switches in `v1.0.7`); all default invocations are **100% byte-identical to VarDictJava v1.8.3**.
- **`v1.1.x`:** Default invocations remain **100% byte-identical to VarDictJava v1.8.3** (`--parity` implied by default). Introduces the **`--modern`** CLI switch (opting into the curated suite of clinical/spec corrections and parity-breaking speedups) and explicit **`--parity`** (`--compat vardictjava-1.8.3`) so pipeline scripts can forward-compatibly pin v1.8.3 behavior before upgrading to v2.0.
- **`v2.0.0+`:** **`--modern`** becomes the default out-of-the-box behavior, while **`--parity`** permanently reproduces exact VarDictJava v1.8.3 byte-for-byte behavior with a single flag.

Under the hood, each individual correction is still gated by an isolated `--fix <id>` toggle (used for bisection and CI verification) and lands in `--modern` only with:
1. A reproducible testcase demonstrating VarDictJava v1.8.3 behavior.
2. A dedicated unit test in `crates/vardictrust-core`.
3. A [Tier-1](docs/testing/test-matrix-and-ci-audit.md#tier-1) intermediate table and [Tier-2](docs/testing/test-matrix-and-ci-audit.md#tier-2) differential test run proving explicit `--parity` output (and default CLI output across all `v1.x` releases) remains byte-identical to v1.0.0 across all regression suites.
4. A diff-confinement proof verifying that `--fix <id>` alters only targeted records, produces zero off-target diffs, and maintains [Tier-3](docs/testing/test-matrix-and-ci-audit.md#tier-3) truth concordance (see [Master Test Matrix](docs/testing/test-matrix-and-ci-audit.md)).

---

## v1.1.0 Candidates — Unambiguous Bugs with Clear Correct Behavior

| ID | Upstream Reference | v1.8.3 Behavior (Replicated in v1.0.0) | Proposed `--fix` Behavior |
|:---|:-------------------|:--------------------------------------|:--------------------------|
| `402` | [PR #402](https://github.com/AstraZeneca-NGS/VarDictJava/pull/402) | `modifyCigar` yields `…D-1M…` for reads like `104M1I2M2D1M1D5M1D2M9D27M`; Java throws, read skipped | Repair CIGAR to `104M14D11I9D27M` per PR #402 |
| `292` | [#292](https://github.com/AstraZeneca-NGS/VarDictJava/issues/292) | `SHIFT3` calculation can yield `end < start` with `-3` | Bound `SHIFT3` offset to ensure non-negative coordinate spans |
| `342` | [#342](https://github.com/AstraZeneca-NGS/VarDictJava/issues/342) | `<dup-N>` emitted in `ALT` (non-standard VCF token) | Emit standard symbolic `<DUP>` |
| `369` | [#369](https://github.com/AstraZeneca-NGS/VarDictJava/issues/369) | Variants emitted outside target region `-R` | Restrict variant emission strictly to requested genomic intervals |
| `378` | [#378](https://github.com/AstraZeneca-NGS/VarDictJava/issues/378) | Paired + amplicon mode runs `AmpliconMode` on `bam1` in `.var` stream (`VarDictLauncher.java:56-61`), breaking `var2vcf_paired.pl` column alignment; VCF mode rejects paired + `-a` with `VardictError::Unsupported` | Support paired somatic amplicon calling or align tabular columns for `var2vcf_paired` |
| `ub1` | [UB-1](docs/upstream-bugs.md#ub-1--p-spurious-all-zero-record-at-position-00) | Position-0 row emitted in `.var` intermediate stream | Suppress spurious position-0 sentinel row |
| `ub2` | [UB-2](docs/upstream-bugs.md#ub-2--p-null-record-row-cell-formatting-empty-geno-00-bias-empty-flank-seqs) | `-p` null-record row has empty `GENO` and `"0;0"` bias | Emit standard null placeholders or omit row per VCF spec |
| `ub5` | [UB-5](docs/upstream-bugs.md#ub-5-simplemode-structural-variant-split-read-requirement-var2vcf_validpl197-205) | `var2vcf_valid.pl` drops `<DUP>`/`<INV>` when `splitreads < 1` | Preserve discordant-pair SVs lacking split reads |
| `ub6` | [UB-6](docs/upstream-bugs.md#ub-6-mate-pair-boundary-parameter-in-mark_dup_sv-structuralvariantsprocessorjava1468) | `markDUPSV` pre-decrements `pe`, truncating window boundary by 1 bp | Pass pre-decrement `pe + 1` window boundary |
| `ub7` | [UB-7](docs/upstream-bugs.md#ub-7-multi-allelic-site-single-winner-collapsing-var2vcf_validpl135-205) | Single-winner multi-allelic collapsing drops secondary alleles without `-A` | Retain secondary qualifying alleles in multi-allelic VCF records |
| `duprate-live` | [`CigarParser.java:78,162-163`](https://github.com/AstraZeneca-NGS/VarDictJava/blob/master/src/main/java/com/astrazeneca/vardict/modules/CigarParser.java#L78) | `initFromScope` copies `totalReads = 0` and `duplicateReads = 0` before `processor.nextRecord()` loop, so `duprate` is always `0` even with `-t` | Propagate live `RecordPreprocessor` duplicate/total read counts to `CigarParser` after the read loop |
| `rounding-underflow-tie` | `Utils.java:101-103` (OpenJDK `DigitList` fixed-point underflow branch) | `roundHalfEven` rounds `±5e-(d+1)` (`d in {3, 4, 5, 8}`) to `±0.0` despite exact binary `f64` exceeding `5×10^-(d+1)` | Round exact binary `f64` values exceeding `5×10^-(d+1)` up to `±10^-d` |
| `eqx-read-length` | `CigarParser.java:1605-1629,1990-1991` | `getSoftClippedLength`, `getMatchInsertionLength`, and SV `ALIGNED_LENGTH_MND` exclude `=`/`X` (`SequenceMatch`/`SequenceMismatch`), yielding `rlen1 = 0` and negative `PMEAN` (0 calls) on `=`/`X` BAMs | Count `=` and `X` alongside `M` in read-length and SV span calculations |
| `var2vcf-chrx-drop` | `var2vcf_valid.pl:265-285`, `var2vcf_paired.pl:271-291` (`sub reorder`) | When both `"X"` and `"chrX"` (or `"Y"`+`"chrY"`, `"MT"`+`"chrM"`, `"M"`) appear in the same stream, `if/elsif` plus the final `next if` skip drops `"chrX"`/`"chrY"`/`"chrM"`/`"M"` | Emit both unprefixed and `chr`-prefixed sex/mitochondrial contigs in `sub reorder` |

---

## Landed as Default (parity-restoring; never gated)

Owner-ruled 2026-10-22: these two Rust-only tiling artifacts both RESTORE VarDictJava parity (flag-on MATCHes java everywhere measured); the v1.0.2 defaults were never oracle-true here. Both shipped UNGATED as default behavior; the `--fix` ids were unreleased (PR #29 unmerged, no users), so the plumbing was removed, not kept as no-op aliases. (Governance exception to the §5 `--fix` rule above; AGENTS.md §5 amendment proposed as follow-up.)

| ID | Upstream Reference | v1.0.2 Behavior (Not True Parity) | Default Behavior Now |
|:---|:-------------------|:----------------------------------|:---------------------|
| `splice-tile-ownership` | Rust-only tiling artifact (no upstream issue — Java counts each junction once over the whole region, unconditionally: `CigarParser.java:86-110`, key `CigarParser.java:1289-1301`) | `-i` splice junctions double-counted under default 100 kb tiling: every tile whose halo-reaching reads touch a junction emits it | Emit each splice junction from exactly one tile (anchor-ownership by key `A`, region-edge clamped) |
| `mate-window-ownership` | Rust-only tiling artifact, sibling-audit G3 (no upstream issue — Java piles the whole region once, so a mate-window merge can never duplicate a row) | Tiled runs duplicate rows whose mate breakpoint falls in the same BED region but a different tile: the merging tile's `partial_pipeline` plants out-of-active-range keys (`realign.rs`, `sv.rs`) that pass the whole-region builder gate — intermediate tables (#1,#3) and somatic VCF (#4); simple VCF masks exact-`start_pos` cases via `pending` | Gate merged keys to the owning tile's active range (merge-time gating in `realign.rs`, planted-key attribution + builder ownership gate with Task-3-style edge clamping for `sv.rs`); single-region runs own everything so the gate never fires there. Accepted corner: realign-planted SV key outside the whole BED with no sv re-plant is dropped (see task-7-report F2) |
| `decimal-format-underflow-tie` (`M-D1`) | `Utils.java:101-103` (`DecimalFormat` `DigitList` fixed-point underflow tie) | Rust `round_half_even` rounded `±5e-(d+1)` (`d in {3, 4, 5, 8}`) to `±10^-d` using exact binary `f64` value | Match OpenJDK `DecimalFormat` (`±5e-(d+1) -> ±0.0` for `d in {3, 4, 5, 8}`), verified across `testdata/format/java_fmt_oracle.tsv` |
| `fisher-java-commons-math3` (`M-D2`) | `FisherExact.java:67-68, 186-197` (`org.apache.commons.math3.distribution.HypergeometricDistribution`) | `--fisher` mode called R `dhyper`/`phyper` in `logdc_dhyper` and `pnhyper`, diverging by 1 ulp on 5-decimal ties | Port Commons Math 3 `SaddlePointExpansion.logBinomialProbability` and `innerCumulativeProbability` in `java_fisher` mode, verified across `testdata/fisher/java_fisher_oracle.tsv` |
| `fisher-r-fround5` (`M-D3`) | `teststrandbias.R:24-26`, `testsomatic.R:31-35` (R 4.x `round(x, 5)` / `nmath/fround.c`) | R-mode `round_as_r` multiplied `x * 1e5` before `HALF_EVEN`, losing sub-ulp direction on exact dyadic `.5` ties | Port R 4.x 5-decimal candidate distance comparison (`floor(x*1e5)/1e5` vs `ceil(x*1e5)/1e5`), verified across `testdata/fisher/r_fisher_oracle.tsv` |
| `split-var-dedup` (`VarDictRust#91`) | Rust-only tiling artifact (no upstream issue — `VarDictJava v1.8.3` runs untiled per BED region, emitting each `.var` intermediate row once) | Split tiles emitted overlapping `.var` rows at tile halo boundaries before `var2vcf` post-processing deduplicated VCF records | `SplitVarDeduper` in `crates/vardictrust-core/src/lib.rs` enforces single-tile `.var` intermediate ownership (`tile.start..=tile.end`) across split tiles in `SimpleMode` and `SomaticMode` with a sliding `O(halo)` (`200 kbp`) boundary-retention window (`BTreeMap::split_off`), preserving within-tile multiplicity and `-p` shifted-singleton position-0 companion triples to match untiled `VarDictJava v1.8.3` output |

---

## Needs a Specification First (Debatable Heuristics; No Flag Until Specified)

Graduation from this table to a candidate `--fix` flag requires a community RFC and empirical concordance benchmarking against NIST GIAB and SEQC2 truth sets.

| Area | Upstream Reference | Current Behavior & Scope | Specification Requirement |
|:-----|:-------------------|:-------------------------|:--------------------------|
| **Depth Discrepancy (`VD > DP`)** | [#92](https://github.com/AstraZeneca-NGS/VarDictJava/issues/92), [#316](https://github.com/AstraZeneca-NGS/VarDictJava/issues/316), [#323](https://github.com/AstraZeneca-NGS/VarDictJava/issues/323), [#386](https://github.com/AstraZeneca-NGS/VarDictJava/issues/386), [#400](https://github.com/AstraZeneca-NGS/VarDictJava/issues/400) | Variant depth (`VD`) reported greater than total depth (`DP`) due to divergent filter accounting between pileup and variant extraction | Reconcile `DP` and `VD` accounting rules across filter boundaries and multi-allelic sites |
| **Filter Sensitivity** | [#62](https://github.com/AstraZeneca-NGS/VarDictJava/issues/62) | Heuristic over-filtering of candidate variants in low-frequency somatic calling | Calibrate filtering thresholds to balance sensitivity and specificity across depths |
| **Realignment Drops** | [#363](https://github.com/AstraZeneca-NGS/VarDictJava/issues/363) | Candidate variants dropped during local realignment around complex indels | Specify unified realignment scoring and fallback behavior for complex haplotypes |
| **Low-AF Noise Filtering** | [#381](https://github.com/AstraZeneca-NGS/VarDictJava/issues/381) | Ambiguous thresholds for low allele-frequency variant calls (< 1%) | Establish clear empirical noise model for deep targeted sequencing panels |
| **Soft-Clip Realignment** | [#282](https://github.com/AstraZeneca-NGS/VarDictJava/issues/282), [UB-3](docs/upstream-bugs.md#ub-3-var-row-order-across-positions-follows-hashmap-iteration-order-282-family) | Output ordering sensitivity and soft-clip realignment differences | Formalize deterministic soft-clip extension and trimming rules |
| **Indel Boundaries** | [#264](https://github.com/AstraZeneca-NGS/VarDictJava/issues/264) | Boundary alignment and coordinate shifting for adjacent insertions/deletions | Standardize 3' left/right alignment conventions |
| **Insertion Coverage** | [`ToVarsBuilder.java:523-538`](https://github.com/AstraZeneca-NGS/VarDictJava/blob/master/src/main/java/com/astrazeneca/vardict/modules/ToVarsBuilder.java#L523-L538) | `HICOV` calculation excludes insertions from local depth denominator | Clarify whether insertion reads contribute to local coverage metrics |
| **OpenJDK Tree-Bin Ordering** | `javamap.rs:211-217`, `tovars.rs:198-205` | When collision count reaches 8 in capacity >= 64, Java converts bucket to `TreeNode` (Red-Black tree) with `moveRootToFront`. Currently sets observability warning flag. 0 treeifications observed across all 92 fixtures and 10 biological cohorts ([Validation Report](docs/validation/omnibus-empirical-validation.md)). | Port exact OpenJDK `TreeNode` rotation, untreeification, and `moveRootToFront` pointer updates behind `--fix javamap-tree-bins` if synthetic multi-collision inputs require exact order preservation. |
| **Cleanroom Rmath Replacement** | `crates/vardictrust-core/src/stats/rmath/` | R `nmath` hypergeometric distribution routines are currently used for float-exact Fisher test calculations. While documented in `LICENSE-3RD-PARTY`, cleanroom reimplementation (e.g. via Apache-2.0 Commons Math 3 algorithms) will restore pure MIT licensing across the entire workspace. | Implement cleanroom hypergeometric distribution under Apache-2.0/MIT and verify floating-point equivalence ($\varepsilon \le 10^{-6}$) across all 15 benchmark cohorts ([Test Matrix](docs/testing/test-matrix-and-ci-audit.md)). |

---

## Open Upstream PRs Reviewed (v1.0.0 Baseline)

| Upstream PR | Title / Focus Area | Disposition in VarDictRust |
|:------------|:-------------------|:---------------------------|
| [PR #402](https://github.com/AstraZeneca-NGS/VarDictJava/pull/402) | Fix negative CIGAR deletion lengths in `modifyCigar` | Replicated crash-safety in v1.0.0; slated for `--fix 402` in v1.1.0 |
| [PR #371](https://github.com/AstraZeneca-NGS/VarDictJava/pull/371) | Logging and CLI argument handling | Evaluated: no output impact; Rust CLI handles argument parsing natively |
| [PR #372](https://github.com/AstraZeneca-NGS/VarDictJava/pull/372) | Build script and dependency updates | Evaluated: upstream Java/Gradle build maintenance; not applicable to Rust |
| [PR #335](https://github.com/AstraZeneca-NGS/VarDictJava/pull/335) | Exception handling in BAM reading | Evaluated: error handling managed via noodles and Rust `Result` types |

---

## Post-1.0 Features

All post-1.0 features will be introduced as opt-in modes, ensuring that default invocations remain byte-identical to the v1.0.0 baseline and maintain bounded peak memory ($\le 3\text{ GB}$):

- **Unique Molecular Identifiers (UMI):** Integrated duplex and single-index consensus calling (`--umi-tag RX --duplex`) within coordinate-sorted 100 kb streaming tiles (modeled on Fulcrum Genomics `fgumi`).
- **Read-Pair Orientation Bias Modeling (`F1R2`/`F2R1`):** Empirical Beta-Binomial orientation bias filtering (`--orientation-bias-model`) to suppress single-strand FFPE cytosine deamination (`C>T`/`G>A`) and acoustic-shearing 8-oxoG (`C>A`/`G>T`) artifacts at `<2%` VAF.
- **Long-Read NGS Support (`--preset hifi|ont`):** Calibrated homopolymer/indel error profiles, >125 bp SV clip thresholds, supplementary alignment (`SA:Z`) split-read chaining, `MM`/`ML` and `HP`/`PS` tag preservation, and step-bounded 2-piece affine Wavefront Alignment (`WFA2` / `block-aligner`, modeled on PacBio `sawfish`, `trgt`, and `aardvark`).
- **Native Compressed VCF & Run Telemetry:** `--output <PATH>` (plain `.vcf`/`.var` and multithreaded BGZF `.gz`/`.bgz` via `noodles::bgzf::io::MultithreadedWriter`) and `--stats-json <PATH>` pipeline telemetry landed in `v1.0.8` (`VarDictRust#93`); direct `.tbi`/`.csi` index generation (`noodles::tabix`), `--all-alleles` (resolving the legacy `-A` flag collision with `var2vcf_valid.pl`), and `--preset wgs|wes|amplicon|ctdna` are slated for `v1.1.0`.

---

## Audit Follow-ups — 2026-10-16 Full Rust Audit (parity-sensitive, needs oracle proof)

Deferred items from the full production-Rust audit. Each needs `harness/compare.sh` MATCH + diff-confinement before landing; default output stays v1.0 byte-identical.

| ID | Location | Behavior (v1.0 kept) | Proposed `--fix`/opt-in |
|:---|:---------|:---------------------|:-----------------------|
| `audit-somatic-dedup` | `lib.rs:797-993` | Somatic tiles emit directly with no cross-tile `pending` dedup (simple path has `BTreeMap` at `:351`). Tiling is Rust-only; `§10` requires tiled == single-region. | Reuse simple-path position-bucket dedup for somatic intermediate + paired VCF. Gate by tiled-vs-single `compare.sh`. |
| `audit-fisher-stream` | `stats/fisher.rs:40-43` | `support`/`logdc` Vecs materialize hypergeometric support (OOM at 1M-depth). Kept saturating counts only. | Streaming sum without vecs behind flag; prove `eps<=1e-6` on deep cohorts. |
| `audit-stream-pileup` | `lib.rs:84-93`, `bam.rs:347` | Per-tile `Vec<SamRecord>` unbounded (violates `§10` bounded-mem at WGS depth). | Iterator/streaming pileup or depth cap; prove tiled == single. |
| `audit-ref-multiregion` | `reference.rs:101-105` | `get_base` fallback reuses current-region offset for `reference_sequences` (wrong base after multi-region loads; latent single-region). Tests rely on fallback. | LANDED (#47 / PR #48): interval-stack `loaded_chunks` across auxiliary `get_reference` calls, `0` SV diffs across full-genome `HG002` (`4.7M` variants). |
| ` audit-amplicon-expand` | `lib.rs:657` | `for p in insert_start..=insert_end` per-base `BTreeMap` expansion OOMs on huge BED interval. | LANDED (Tasks 4b/4d): interval-overlap test + `-p` streaming, oracle-MATCHED, capped by in-tree mem test. |
| `sv-lazy-mate-pileup` | `sv.rs:874-1017` | `find_inv_sub` calls `partial_pipeline(mstart - 200, mend + 200)` at loop top before checking if any soft-clip exists at `start..end` (`softp`) or matches `reference` (`find_match_rev`). On NovaSeq 2-color WGS (`HG008-T\|N`), dark-cycle poly-G sinks (`chr2:32,916,196`, `25M` reads) trigger $O(N_{\text{windows}} \times 25\text{M})$ wasted mate-region BAM re-scans across `chr2`. | Gate `partial_pipeline` in `find_inv_sub` on first finding a candidate soft-clip (`seq.len() >= 12`) that matches `reference` via `find_match_rev` (plus optional `--max-sv-mate-depth`). |

---

## 2026-11-08 Critical Engineering & Ecosystem Review Priorities

Prioritized engineering and clinical follow-ups from the comprehensive survey of Fulcrum Genomics, PacBio Rust tools, clinical literature, and `VarDictRust` v1.0.7 hotspots (see [`docs/ecosystem-survey.md`](docs/ecosystem-survey.md#critical-review-v107); tracked in GitHub issues `VarDictRust#90`–`VarDictRust#96`):

| Tier | ID | Location | Action & Verification Gate |
|:-----|:---|:---------|:---------------------------|
| **Tier 1 (`v1.0.8` Zero-Diff, `#90`)** | `fastpath-modify-cigar` | `cigar/modify.rs`, `cigar/parse.rs`, `sv.rs`, `tovars.rs`, `variations/mod.rs` | **LANDED (`PR #98`–`#101`):** Single-`Match` (`150M`) in-place CIGAR fast-path and shared endpoint mismatch/clipping helpers (`cigar/modify.rs`); borrowed `&[u8]`/`&str` slices, stack key formatting, and zero-alloc split iterators across `sv.rs`, `tovars.rs`, `variations/mod.rs`, `cigar/parse.rs`, and `output/vcf.rs`; pre-allocated `VariationMap` capacities, `RefCell` scratch reuse (`parse.rs`, `tovars.rs`), and `write!`-based `.var` row formatting (`output/table.rs`). Gate: endpoint-mismatch parity fixtures + `pr_check.sh` 0 diffs. |
| **Tier 1 (`v1.0.x` Zero-Diff, `#90`)** | `parallel-max-rl-scan` | `lib.rs:69, 412-893, 1677, 2368` | Parallelize `scan_multi_tile_max_read_lengths` across `rayon` tiles via prefix-max fold and lower `DEFAULT_UNTILED_SPAN_MAX` from `1_000_000` to `100_000` (25–35% faster unchunked WGS). |
| **Tier 1 (`v1.0.8` Zero-Diff, `#90`)** | `direct-loaded-region` | `reference.rs`, `realign.rs`, `stats/fisher.rs` | **LANDED (`PR #98`, `PR #100`):** Inlined `Reference::get_base` primary-region fast-path with `SEED1_LEN` fallthrough (`reference.rs`); `Cow<[u8]>` `strip_markers` and stack mismatch buffers in realignment (`realign.rs`); eliminated `FisherExact::support: Vec<i32>` heap allocation via `support_len()` algebraic indexing and added `thread_local!` LRU memoization cache (`4,096` entries) for `fisher_exact` (`stats/fisher.rs`). |
| **Tier 1 (`v1.0.8` Zero-Diff, `#90`)** | `bgzf-cache-scale` | `bam.rs:743-744, 886-942, 3186` | **LANDED (`PR #98`):** Expanded `BGZF_BLOCK_CACHE` ($8 \to 32$ shards, `BGZF_CACHE_SHARDS = 32`) to eliminate 32-thread `fetch_block Arc` eviction churn, constructed `build_noodles_region` directly, and replaced `-t` `HashSet<String>` with typed `FxHashSet<DupKey>` (`MateRefKey`). |
| **Tier 1 (`v1.0.8` Zero-Diff, `#91`, `#93`)** | `cli-bed-cram-hardening` | `cli.rs`, `regions.rs`, `runner.rs`, `lib.rs`, `UPSTREAM.toml` | **LANDED (`PR #97`–`#101`):** Cross-tile `.var` dedup via `SplitVarDeduper` + tiled `.var`/`.vcf` conservation tests (`#91`, `lib.rs`); `--output <PATH>` (plain + multithreaded BGZF `.gz`/`.bgz` with `File::sync_all` verification), `--stats-json <PATH>` run telemetry, automatic 64 MiB worker stack floor for `.cram` inputs, `-c 1 -S 2 -E 3 -z` BED column hint (`#54`), `cli.rs` `--version` license sync, and `UPSTREAM.toml` provenance manifest (`#93`). |
| **Tier 2 (`v1.1.0` Cleanroom & Opt-In, `#92`)** | `cleanroom-rmath-mit` | `stats/rmath/`, `stats/fisher.rs:343-478`, `Cargo.toml:9` | Replace the 5 GPL-2.0 R `nmath` files (~350 LoC) with cleanroom Loader (2000) / Apache-2.0 saddle-point functions verified to $\varepsilon \le 10^{-6}$, restoring pure `MIT` workspace licensing. |
| **Tier 2 (`v1.1.0` Cleanroom & Opt-In, `#93`, `#94`, `#95`)** | `parity-modern-mode-ux` | `crates/vardictrust/src/cli.rs`, `main.rs` | Add `--parity` / `--modern` mode switches (`#94`), clinical/spec fixes (`#95`), `.tbi`/`.csi` index generation (`noodles::tabix`), `--all-alleles` ([`#65`](https://github.com/AstraZeneca-NGS/VarDictJava/issues/65)), and `--preset wgs\|wes\|amplicon\|ctdna` (`--output` and `--stats-json` landed in `v1.0.8`, `#93`). |
| **Tier 3 (`v2.0` Next-Gen, `#96`)** | `umi-f1r2-wfa2-longread` | `bam.rs`, `cigar/parse.rs`, `realign.rs`, `sv.rs` | Default to `--modern` (`--parity` preserves v1.8.3 forever); streaming tile-local UMI/duplex consensus (`--umi-tag RX --duplex`), `F1R2`/`F2R1` orientation bias filter, 2-bit DNA `XOR`+`popcnt`, step-bounded WFA2 (`60`-step cap + 6 GB global alignment budget), rolling ring-buffer pileup, and `--preset hifi\|ont`. |

