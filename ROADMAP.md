# VarDictRust Roadmap

VarDictRust v1.0.0 replicates VarDictJava v1.8.3 output byte-for-byte on valid inputs, including the upstream behaviors documented in `DIVERGENCES.md` and [`docs/upstream-bugs.md`](docs/upstream-bugs.md).

From v1.1.0, each correction is opt-in via `--fix <id>` (default off) and lands only with:
1. A reproducible testcase demonstrating VarDictJava v1.8.3 behavior.
2. A dedicated unit test in `crates/vardictrust-core`.
3. A [Tier-1](docs/testing/test-matrix-and-ci-audit.md#tier-1) intermediate table and [Tier-2](docs/testing/test-matrix-and-ci-audit.md#tier-2) differential test run proving default output remains byte-identical to v1.0.0 across all regression suites.
4. A diff-confinement proof verifying that `--fix <id>` alters only targeted records, produces zero off-target diffs, and maintains [Tier-3](docs/testing/test-matrix-and-ci-audit.md#tier-3) truth concordance (see [Master Test Matrix](docs/testing/test-matrix-and-ci-audit.md)).

---

## v1.1.0 Candidates — Unambiguous Bugs with Clear Correct Behavior

| ID | Upstream Reference | v1.8.3 Behavior (Replicated in v1.0.0) | Proposed `--fix` Behavior |
|:---|:-------------------|:--------------------------------------|:--------------------------|
| `402` | [PR #402](https://github.com/AstraZeneca-NGS/VarDictJava/pull/402) | `modifyCigar` yields `…D-1M…` for reads like `104M1I2M2D1M1D5M1D2M9D27M`; Java throws, read skipped | Repair CIGAR to `104M14D11I9D27M` per PR #402 |
| `292` | [#292](https://github.com/AstraZeneca-NGS/VarDictJava/issues/292) | `SHIFT3` calculation can yield `end < start` with `-3` | Bound `SHIFT3` offset to ensure non-negative coordinate spans |
| `342` | [#342](https://github.com/AstraZeneca-NGS/VarDictJava/issues/342) | `<dup-N>` emitted in `ALT` (non-standard VCF token) | Emit standard symbolic `<DUP>` |
| `369` | [#369](https://github.com/AstraZeneca-NGS/VarDictJava/issues/369) | Variants emitted outside target region `-R` | Restrict variant emission strictly to requested genomic intervals |
| `378` | [#378](https://github.com/AstraZeneca-NGS/VarDictJava/issues/378) | Paired + amplicon mode column shift in tabular stream | Correct column alignment (or emit explicit unsupported diagnostic) |
| `ub1` | [UB-1](docs/upstream-bugs.md#ub-1--p-spurious-all-zero-record-at-position-00) | Position-0 row emitted in `.var` intermediate stream | Suppress spurious position-0 sentinel row |
| `ub2` | [UB-2](docs/upstream-bugs.md#ub-2--p-null-record-row-cell-formatting-empty-geno-00-bias-empty-flank-seqs) | `-p` null-record row has empty `GENO` and `"0;0"` bias | Emit standard null placeholders or omit row per VCF spec |
| `ub5` | [UB-5](docs/upstream-bugs.md#ub-5-simplemode-structural-variant-split-read-requirement-var2vcf_validpl197-205) | `var2vcf_valid.pl` drops `<DUP>`/`<INV>` when `splitreads < 1` | Preserve discordant-pair SVs lacking split reads |
| `ub6` | [UB-6](docs/upstream-bugs.md#ub-6-mate-pair-boundary-parameter-in-mark_dup_sv-structuralvariantsprocessorjava1468) | `markDUPSV` pre-decrements `pe`, truncating window boundary by 1 bp | Pass pre-decrement `pe + 1` window boundary |
| `ub7` | [UB-7](docs/upstream-bugs.md#ub-7-multi-allelic-site-single-winner-collapsing-var2vcf_validpl135-205) | Single-winner multi-allelic collapsing drops secondary alleles without `-A` | Retain secondary qualifying alleles in multi-allelic VCF records |
| `splice-tile-ownership` | Rust-only tiling artifact (no upstream issue — Java counts each junction once over the whole region, unconditionally: `CigarParser.java:86-110`, key `CigarParser.java:1289-1301`) | `-i` splice junctions double-counted under default 100 kb tiling: every tile whose halo-reaching reads touch a junction emits it | Emit each splice junction from exactly one tile (anchor-ownership by key `A`, region-edge clamped). Graduation: RFC + GIAB/SEQC2 truth concordance pending (Tier-3 rule above) |
| `mate-window-ownership` | Rust-only tiling artifact, sibling-audit G3 (no upstream issue — Java piles the whole region once, so a mate-window merge can never duplicate a row) | Tiled runs duplicate rows whose mate breakpoint falls in the same BED region but a different tile: the merging tile's `partial_pipeline` plants out-of-active-range keys (`realign.rs`, `sv.rs`) that pass the whole-region builder gate — intermediate tables (#1,#3) and somatic VCF (#4); simple VCF masks exact-`start_pos` cases via `pending` | Gate merged keys to the owning tile's active range (merge-time gating in `realign.rs`, planted-key attribution + builder ownership gate with Task-3-style edge clamping for `sv.rs`); single-region runs own everything so the gate never fires there. Library default on (§10 tiled==single-region); CLI default off (v1.0 parity) |

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

- **Unique Molecular Identifiers (UMI):** Integrated duplex and single-index consensus calling within coordinate-sorted streaming windows.
- **Long-Read NGS Support:** Calibrated error profiles and alignment handling for PacBio HiFi and Oxford Nanopore.

---

## Audit Follow-ups — 2026-10-16 Full Rust Audit (parity-sensitive, needs oracle proof)

Deferred items from the full production-Rust audit. Each needs `harness/compare.sh` MATCH + diff-confinement before landing; default output stays v1.0 byte-identical.

| ID | Location | Behavior (v1.0 kept) | Proposed `--fix`/opt-in |
|:---|:---------|:---------------------|:-----------------------|
| `audit-somatic-dedup` | `lib.rs:797-993` | Somatic tiles emit directly with no cross-tile `pending` dedup (simple path has `BTreeMap` at `:351`). Tiling is Rust-only; `§10` requires tiled == single-region. | Reuse simple-path position-bucket dedup for somatic intermediate + paired VCF. Gate by tiled-vs-single `compare.sh`. |
| `audit-fisher-stream` | `stats/fisher.rs:40-43` | `support`/`logdc` Vecs materialize hypergeometric support (OOM at 1M-depth). Kept saturating counts only. | Streaming sum without vecs behind flag; prove `eps<=1e-6` on deep cohorts. |
| `audit-stream-pileup` | `lib.rs:84-93`, `bam.rs:347` | Per-tile `Vec<SamRecord>` unbounded (violates `§10` bounded-mem at WGS depth). | Iterator/streaming pileup or depth cap; prove tiled == single. |
| `audit-ref-multiregion` | `reference.rs:101-105` | `get_base` fallback reuses current-region offset for `reference_sequences` (wrong base after multi-region loads; latent single-region). Tests rely on fallback. | Track per-buffer start offsets; prove via multi-region oracle. |
| `audit-amplicon-expand` | `lib.rs:657` | `for p in insert_start..=insert_end` per-base `BTreeMap` expansion OOMs on huge BED interval. | LANDED (Tasks 4b/4d): interval-overlap test + `-p` streaming, oracle-MATCHED, capped by in-tree mem test. |
