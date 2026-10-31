# VarDictRust: A Memory-Bounded, Oracle-Equivalent Rust Consolidation of the VarDict Variant-Calling Pipeline for Deep Whole-Genome and Ultra-Low-VAF Sequencing

**Authors:** [Palak Sheth](https://www.linkedin.com/in/palak-p-sheth/)<sup>1,*</sup>

<sup>1</sup>Independent Researcher  
<sup>*</sup>Correspondence: [palakpsheth@gmail.com](mailto:palakpsheth@gmail.com) | [https://www.linkedin.com/in/palak-p-sheth/](https://www.linkedin.com/in/palak-p-sheth/)

---

## Abstract

### Background
`VarDict` (`VarDictJava` v1.8.3) is a heuristic-aware local realigner and variant caller. It serves clinical oncology, rare-disease, and circulating tumor DNA (ctDNA) minimal residual disease (MRD) pipelines. Its three-language execution pipeline streams binary alignment maps (BAMs) through a Java Virtual Machine (JVM) pileup engine. External R (`teststrandbias.R` / `testsomatic.R`) and Perl (`var2vcf_valid.pl` / `var2vcf_paired.pl`) scripts add inter-process serialization overhead. Together they cause unbounded per-interval memory growth and runtime crashes on deep sequencing datasets. The C++ reimplementation (`vardictcpp`) accelerates pre-chunked exome panels. It omits splicing, amplicon VCF, structural variant, and CRAM modes. It diverges statistically and retains unbounded per-interval memory on full chromosomes.

### Results
We built `VarDictRust` (`v{{m:version}}`), a single-binary, memory-safe Rust consolidation of the `VarDictJava` v1.8.3 pipeline. The core engine compiles under `#![forbid(unsafe_code)]` and follows the [rewrites.bio](https://rewrites.bio) specification. After replicating R's 5-decimal round-half-to-even (`HALF_EVEN`) Fisher exact test, `VarDictRust` matches `{{m:oracle_version}}` (`ε = {{m:eps_tol}}`) across `{{m:s2_modes_rust_match}}/{{m:s2_modes_total}}` operational modes (`95%` Clopper-Pearson confidence interval [CI] `{{m:cp95_13_13}}`). Parity is Tier-1/2 throughout: intermediate `.var` plus VCF 4.3, record-aligned on CHROM/POS/REF/ALT. Coverage extends to `{{m:s4_ext_cov_modalities}}` functional modalities (`95%` CI `{{m:cp95_17_17}}`) and `{{m:s4_multimodal_slices}}` multimodal slices (`95%` CI `{{m:cp95_18_18}}`). It also covers `{{m:s5_f0_noise_floor_strict_match}}` unfiltered (`-f {{m:f_zero}}`) somatic noise-floor sweeps (`95%` CI `{{m:cp95_noise_18_18}}`). It further spans `{{m:s5_cross_platform_arms}}` short-read chemistries and completes `{{m:s6_ga4gh_cells}}/{{m:s6_ga4gh_cells}}` Rust truth evaluations across GA4GH `v3.5` strata (`95%` CI `{{m:cp95_70_70}}`; Java–Rust parity holds on all completing paired runs — see `Table 6`).

On contiguous human chromosome 1 (`-R chr1:1-249250621`, `{{m:s2_rchr1_mbp}}` Mbp), `VarDictRust` finishes in `{{m:s2_rchr1_rust_wall_s}}` s. Peak resident set size (RSS) is `{{m:s2_rchr1_rust_rss_mb}}` MB. `VarDictJava` (`-Xmx200g`) needs `{{m:s2_rchr1_java200g_wall_s}}` s and `{{m:s2_rchr1_java200g_rss_gib}}` GiB. That is a `{{m:s2_rchr1_speedup_vs_java}}`× speedup at `{{m:s2_rchr1_rss_ratio_vs_java}}`× lower RAM. Against `{{m:cpp_version}}`, it is `{{m:s2_rchr1_speedup_vs_cpp}}`× faster at `{{m:s2_rchr1_rss_ratio_vs_cpp}}`× lower RAM. `VarDictRust` also reaches `{{m:s4_exome_32t_speedup_vs_1t}}`× `32`-thread whole-exome scaling. That is `{{m:s4_exome_th32_speedup_vs_java}}`× faster at `{{m:s4_exome_32t_rss_ratio_vs_java}}`× lower RAM than `VarDictJava`.

On a `128` GiB node, `VarDictRust` completes unchunked `{{m:s5_hg002_bam_gb}}` GB HG002 germline whole-genome sequencing (WGS) (`{{m:s5_hg002_contigs_rust_rss_gb}}` GB RAM). It also completes `{{m:s5_hg008_bam_gb}}` GB `364`× HG008 somatic WGS (`{{m:s5_hg008_contigs_rust_rss_gb}}` GB RAM). `VarDictJava` runs out of memory (`>{{m:s5_hg002_contigs_java_rss_gb}}` GB) and `vardictcpp` times out (`{{m:s5_hg008_contigs_cpp_rss_gb}}` GB). These runs are completion results, not parity: no Java oracle finished them. Removing `≤{{m:giab_arm_b_max_kbp}}` kbp BED pre-chunking recovers `+{{m:s6_v50q_gain_gbp}}` Gbp of `NIST GIAB v5.0q T2T-Q100` sequence. Recall rises by `+{{m:s6_v50q_snv_r_gain_pct}}` single-nucleotide variant [SNV] and `+{{m:s6_v50q_indel_r_gain_pct}}` insertion/deletion [indel] recall. Gains reach `+{{m:s6_cmrg_snv_r_gain_pct}}` SNV and `+{{m:s6_cmrg_indel_r_gain_pct}}` indel recall across `273` `NIST CMRG v1.00` genes. Mosaic recall is `{{m:s6_mosaic_rust_tp}}/{{m:s6_mosaic_truth_total}}` (`{{m:s6_mosaic_rust_r_pct}}`) on `NIST Mosaic SNV v1.10`.

### Conclusions
`VarDictRust` consolidates a three-language variant-calling pipeline into a `{{m:s5_h4_musl_mb}}` MB static Linux binary with exact numerical and allele equivalence (`ε = {{m:eps_tol}}`). Bounded tile-plus-halo streaming (`O(tile + halo)` memory) and zero runtime JVM/R/Perl dependencies (`{{m:s5_h4_conda_pkgs}}` Conda packages eliminated) enable deterministic calling across targeted panels, RNA-seq splice junctions, unchunked deep WGS, and `{{m:s6_tit_vaf_1_pct}}` variant allele frequency (VAF) titrations.

**Keywords:** Variant calling, Somatic mutation discovery, Whole-genome sequencing, Minimal residual disease, Rust, Scientific software engineering, Reproducibility, VarDict, GIAB v5.0q, HG008

> **Key Points**
>
>
> - **Single-Binary Oracle Equivalence:** `VarDictRust` (`v{{m:version}}`, `{{m:s5_h4_musl_mb}}` MB static binary) replaces `VarDictJava -> R -> Perl` (`{{m:s5_h4_conda_pkgs}}` Conda packages). It holds `100%` Tier-1/2 parity (intermediate `.var` plus VCF, `ε = {{m:eps_tol}}`) across `{{m:s2_modes_rust_match}}/{{m:s2_modes_total}}` modes, `{{m:optc_match_arms}}/{{m:optc_total_arms}}` sweep sub-arms, and `{{m:s6_ga4gh_cells}}/{{m:s6_ga4gh_cells}}` completed Rust GA4GH `v3.5` truth evaluations (Java–Rust parity on all completing paired runs).
>
> - **Bounded $O(\text{tile} + \text{halo})$ Streaming Memory:** Spans $> {{m:s1_tile_1m}}$ bp split into `{{m:rust_default_tile_kbp}}` kbp tiles with `{{m:rust_halo_kbp}}` kbp halos (`{{m:rust_splice_halo_kbp}}` kbp in `-i` splice mode). Chromosome 1 (`-R chr1:1-249250621`) then peaks at `{{m:s2_rchr1_rust_rss_mb}}` MB RSS in `{{m:s2_rchr1_rust_wall_s}}` s. That is `{{m:s2_rchr1_speedup_vs_java}}`× faster and `{{m:s2_rchr1_rss_ratio_vs_java}}`× leaner than `VarDictJava`.
>
> - **Unchunked Full-Genome Discovery:** Calling unchunked chromosomes (`chr1`–`chr22`) in `{{m:s5_hg002_contigs_rust_rss_gb}}`–`{{m:s5_hg008_contigs_rust_rss_gb}}` GB RAM avoids `VarDictJava` out-of-memory failures. It gains `+{{m:s6_v50q_snv_r_gain_pct}}` SNV / `+{{m:s6_v50q_indel_r_gain_pct}}` indel recall on `NIST GIAB v5.0q T2T-Q100`. It gains `+{{m:s6_cmrg_snv_r_gain_pct}}` SNV / `+{{m:s6_cmrg_indel_r_gain_pct}}` indel recall on `CMRG v1.00`.
>
> - **Somatic Indel & Low-VAF Sensitivity on CPUs:** `VarDictRust` exceeds `ClairS`, `DeepSomatic`, `Lancet2`, and `Ultima DeepVariant Somatic` in `HG008 v0.4` somatic indel recall (`{{m:s6_HG008_FullWGS_v04_All_VarDictRust_v106_indel_r_pct}}`). It achieves `{{m:s6_titration_0p5vaf_recall_ratio_vs_m2}}`× higher SNV recall than `Mutect2` and `Strelka2` at `{{m:s6_tit_vaf_1_pct}}` VAF (`1%` purity).

---

## Background

Germline variant callers such as `FreeBayes` [3], `Strelka2` [4], `Octopus` [5], and `DeepVariant` [6] model diploid genotype likelihoods. Clinical oncology assays must instead detect somatic SNVs, multinucleotide variants (MNVs), indels, and structural variants (SVs). They operate across variable tumor purities and cell-free DNA (cfDNA) fractions down to `{{m:s6_tit_vaf_1_pct}}`–`{{m:s6_tit_vaf_2_pct}}` VAF [1, 2, 7, 8]. Introduced by Lai et al. in 2016 [1], `VarDict` (`VarDictJava` v1.8.3) integrates local soft-clip and indel realignment (`CigarModifier` and `VariationRealigner`). It adds composite MNP phasing, amplicon primer-edge filtering (`-a`), RNA-seq splice-junction accounting (`-i`), and inline SV breakpoint assembly. `VarDict` is deployed across clinical panel pipelines, `bcbio-nextgen`, `nf-core/sarek`, and multi-caller somatic consensus ensembles [9, 10].

Reference `VarDictJava` v1.8.3 has four structural limits on deep WGS and cloud infrastructure:

1. **Three-Language Pipeline and IPC Overhead.** `VarDictJava` emits an intermediate `36`- to `61`-column tab-delimited `.var` stream. Column counts are `Simple` `36/38`, `Amplicon` `38/40`, and `Somatic` `55/61` without/with `--fisher`. That stream pipes to an R script (`teststrandbias.R` or `testsomatic.R`) for Fisher exact testing. It then pipes to a Perl script (`var2vcf_valid.pl` or `var2vcf_paired.pl`) for filtering and VCF 4.3 formatting (`Figure 1`). Even when `--fisher` runs inside Java, the Perl post-processor remains mandatory, requiring `{{m:s5_h4_conda_pkgs}}` Conda packages across OpenJDK, R (`statmod`), and Perl.

2. **Unbounded JVM Heap Retention on Contiguous Spans.** `VarDictJava` allocates per-position pileup maps (`Map<Integer, VariationMap<String, Variation>>` and reference base arrays) across the full BED interval or `-R` region. Variant finalization (`ToVarsBuilder`) and realignment (`VariationRealigner`) only run afterwards. On megabase or chromosome-scale intervals, RSS grows linearly with interval length (`O(L × depth)`). That causes `java.lang.OutOfMemoryError` crashes or garbage-collection stalls above `100`–`200` GiB RAM (`Table 2`, `Table 5`).

3. **Runtime Exceptions and Non-Determinism on Edge-Case Reads.** Over `30` upstream `VarDictJava` issues report `CigarModifier.modifyCigar` crashes. Affected issues are `#7`, `#38`, `#187`, `#197`, `#243`, `#249`, `#253`, `#282`, `#283`, `#320`, `#341`, `#380`, and `#384`. The crashes hit supplementary or heavily soft-clipped reads. The crash types are `NullPointerException`, `StringIndexOutOfBoundsException`, and `ArrayIndexOutOfBoundsException`. Other issues cover non-deterministic row ordering from unseeded `java.util.HashMap` iteration (`#282`), broken CRAM decoding (`#249`, `#380`), and truncated `##contig` VCF headers (`#320`).

4. **Incomplete C++ Port (`vardictcpp`).** Wolff (2026) introduced `{{m:cpp_version}}`, a C++/htslib reimplementation targeting single-sample and paired somatic speedups [11]. Our 3-way audit identifies ten functional and statistical divergences in `vardictcpp` (`Table 8`; `Supplementary Table S6`). Specifically, `vardictcpp` omits `##contig` and `INFO/END` VCF headers. It omits CRISPR (`-J`/`-j`), structural variant (`--nosv` forced by default), and RNA splice (`-i`) modes. It emits raw tab-separated text instead of VCF in amplicon (`-a`) mode. It rejects concatenated CLI flags. It diverges on CRAM decoding (`-G`), `b37` numeric contigs, and MNP phasing/tie-breaking. It consumes unbounded memory on large intervals (`{{m:s5_hg002_contigs_cpp_rss_gb}}` GB on HG002 WGS; `{{m:s5_hg008_contigs_cpp_rss_gb}}` GB with a `3`-hour timeout on HG008 WGS). Pre-chunked regions gain `+{{m:s5_hg002_contigs_cpp_excess}}` duplicate boundary calls (`--chunk {{m:cpp_chunk_bp}}`).

We built `VarDictRust` (`v{{m:version}}`) under the [rewrites.bio](https://rewrites.bio) prime directive. It produces a memory-safe, single-binary Rust implementation. Its output is byte- and float-equivalent (`ε = {{m:eps_tol}}`) to `VarDictJava` v1.8.3 on all valid inputs. It eliminates runtime panics and bounds peak memory to `O(tile + halo)` regardless of chromosome length or sequencing depth.

<a id="fig1_architecture_pipeline"></a>

![Figure 1: Multi-suite cloud benchmark architecture of VarDictRust.](figures/fig1.svg)

**Figure 1. Multi-suite cloud benchmark architecture of `VarDictRust`.** **(a)** Six evaluation suites plus the multi-regime parity sweep (`{{m:total_cells}}` experiments total), with per-suite experiment and archived-artifact counts; the hardware lock and comparator set are bannered in the panel. **(b)** Completed benchmark experiments per evaluation suite. **(c)** Archived verification artifacts per suite (`Supplementary Table S1`).

---

## Results

### Tier-1/2 Oracle Equivalence (Exact at `ε = {{m:eps_tol}}`, `HALF_EVEN`, Record-Aligned) Across Operational Modes and the 87-Configuration Sweep

Across our `{{m:total_cells}}`-experiment campaign (`Table 1`, `Figure 2a`, `Supplementary Table S1`), we tested `VarDictRust` (`v{{m:version}}`) against `{{m:oracle_version}}`. We checked both Tier-1 intermediate `.var` parity and Tier-2 VCF 4.3 parity. Tier-1 covers `Simple` `36/38`, `Amplicon` `38/40`, and `Somatic` `55/61` columns without/with `--fisher`. Tier-2 covers `chrom`, `pos`, `ref`, `alt`, `QUAL`, `FILTER`, `INFO` tags `SBF`, `ODDRATIO`, `SSF`, `MSI`, `MSILEN`, `SHIFT3`, `TYPE`, `STATUS`, and `FORMAT` genotypes. All comparisons use `ε = {{m:eps_tol}}` with `HALF_EVEN` rounding, record-aligned on CHROM/POS/REF/ALT.

Across the baseline ceiling, multi-region boundary deduplication, and tile-size invariance benchmarks, `VarDictRust` matched `VarDictJava` with zero discrepancies (`Table 1`; `Supplementary Table S1`):

- **EGFR/MRE11A Amplicon (`-a {{m:amp_param}}`):** `VarDictRust` emitted all `{{m:s1_amplicon_vars}}` variants in `{{m:s1_amplicon_rust_wall_s}}` s (`{{m:s1_amplicon_rust_rss_mb}}` MB RSS) vs `{{m:s1_amplicon_java_wall_s}}` s (`{{m:s1_amplicon_java_rss_mb}}` MB) for `VarDictJava` (`{{m:s1_amplicon_speedup_vs_java}}`× faster, `{{m:s1_amplicon_rss_ratio_vs_java}}`× lower RSS). `vardictcpp` (`{{m:s1_amplicon_cpp_wall_s}}` s, `{{m:s1_amplicon_cpp_rss_mb}}` MB) emits raw `.var` text instead of VCF in amplicon mode (`Table 8`; `Supplementary Table S6`).

- **NIST GIAB HG008 Somatic Slice:** `VarDictRust` matched all `{{m:s1_hg008_vars}}` candidate variants in `{{m:s1_hg008_rust_wall_s}}` s (`{{m:s1_hg008_rust_rss_mb}}` MB RSS). `VarDictJava` took `{{m:s1_hg008_java_wall_s}}` s (`{{m:s1_hg008_java_rss_gib}}` GiB RSS). That is `{{m:s1_hg008_speedup_vs_java}}`× faster at `{{m:s1_hg008_rss_ratio_vs_java}}`× lower RSS (`vardictcpp`: `{{m:s1_hg008_cpp_wall_s}}` s, `{{m:s1_hg008_cpp_rss_mb}}` MB).

- **HG002 Chromosome 20 Multi-Region Boundary:** `VarDictRust` matched all `{{m:s1_gate7_vars}}` variants in `{{m:s1_gate7_rust_wall_s}}` s (`{{m:s1_gate7_rust_rss_mb}}` MB RSS). `VarDictJava` took `{{m:s1_gate7_java_wall_s}}` s (`{{m:s1_gate7_java_rss_gib}}` GiB RSS). That is `{{m:s1_gate7_speedup_vs_java}}`× faster at `{{m:s1_gate7_rss_ratio_vs_java}}`× lower RSS.

- **Tile-Size Invariance:** We swept four tile sizes (`{{m:s1_tile_sizes_bp_list}}` bp, i.e., `{{m:s1_tile_sizes_mbp_list}}` Mbp; Methods) on SEQC2 HCC1395. Every size produced the exact same `{{m:s1_tile_vars}}` variants as single-region `VarDictJava` (`4/4` `MATCH`).

The `13`-mode functional sweep (`Table 1`, `Supplementary Table S6`) covers seven call families. These are Simple (`--fisher` and R pipe), Somatic (`--fisher` and R pipe), Amplicon (`-a`), and Splice (`-i`). They include Structural Variant (`SV` and `--nosv`), CRAM (`-G`), and all-positions (`-f {{m:f_zero}}`) calling. `VarDictRust` matched `VarDictJava` across `{{m:s2_modes_rust_match}}/{{m:s2_modes_total}}` modes (`100%` Tier-1/2 parity, exact Clopper-Pearson `95%` CI `{{m:cp95_13_13}}`). On the `{{m:s3_rna_bam_gib}}` GiB HG002 GM24385 RNA-seq alignment (`Table 1`, `Supplementary Table S1`), `VarDictRust` `-i` splice mode matched all `{{m:s3_rna_splice_rows}}` `.var` and `.splice.tsv` exon-junction rows. It took `{{m:s3_rna_rust_wall_s}}` s (`n = 3` technical replicates, range `{{m:s3_rna_rust_wall_range}}` s) at `{{m:s3_rna_rust_rss_mb}}` MB RSS. `VarDictJava` took `{{m:s3_rna_java_wall_s}}` s (`{{m:s3_rna_java_wall_range}}` s) at `{{m:s3_rna_java_rss_gib}}` GiB. That is `{{m:s3_rna_speedup_vs_java}}`× faster at `{{m:s3_rna_rss_ratio_vs_java}}`× lower RSS. `vardictcpp` does not implement `-i` (`Supplementary Table S6`). `VarDictRust` also passed `{{m:s4_ext_cov_modalities}}` extended coverage modalities (`95%` CI `{{m:cp95_17_17}}`). It passed `{{m:s4_multimodal_slices}}` multimodal slices (`{{m:s4_multimodal_vars}}` variants, `95%` CI `{{m:cp95_18_18}}`). Lossless BAM-vs-CRAM twin calling showed zero diff (`{{m:s4_cram_twin_wall_s}}` s, `{{m:s4_cram_twin_rss_mb}}` MB RSS).

Our `{{m:optc_total_arms}}`-configuration parameter sweep (`Table 1`, `Supplementary Table S2`) spans four calling regimes across `12` suites. These are standard `-f {{m:f_std}}`, WGS MRD `-f {{m:f_mrd_wgs}}`, ultra-deep `45000`× MRD `-f {{m:f_mrd_ud}}`, and unfiltered `-f {{m:f_zero}}` calling. `VarDictRust` achieved `{{m:optc_match_arms}}/{{m:optc_total_arms}}` (`{{m:optc_match_pct}}%`, `95%` CI `{{m:cp95_optc_sweep}}`) exact Tier-1 `.var` and Tier-2 `.vcf` parity with `VarDictJava`. Standard BED intervals ($\le {{m:s1_tile_1m}}$ bp) run untiled by default. The two multi-megabase single-region stress suites used single-window scope matching instead. We passed `--tile-size 10000000` on `10` Mbp `chr1` `-f {{m:f_zero}}` (`{{m:sw_span_f0_rows}}` `.var` rows / `{{m:s4_span_f0_vars}}` VCF variants). We passed `--tile-size 250000000` on `{{m:s2_rchr1_mbp}}` Mbp `-R chr1:1-249250621` (`{{m:sw_rchr1_rows}}` `.var` rows). Both produced `100%` Tier-1/2 sorted `.var` and `.vcf` parity.

{{table:table1}}

---

### Replication and Extension of the Wolff Acute Myeloid Leukemia Benchmark

We replicated the Wolff [11] acute myeloid leukemia (AML) targeted panel cohort (`SRR15006375`, `SRR15006376`, `SRR15006386`, `SRR15006540`, `SRR8657348`, and paired somatic `SRR15006376|SRR15006386`, spanning `{{m:s2_SRR15006386_intervals}}` to `{{m:s2_SRR8657348_intervals}}` BED intervals; `Table 2`, `Figure 2`, `Supplementary Table S3`) at `1` and `8` threads.

Across all `12` sample-by-thread configurations (`n = 3` technical replicates per engine), `VarDictRust` matched `VarDictJava` calls. Coverage was `100%` of `.var` (`{{m:s2_SRR15006375_var_rows}}`–`{{m:s2_SRR8657348_var_rows}}` rows) and `.vcf` (`{{m:s2_SRR15006375_vcf_vars}}`–`{{m:s2_SRR8657348_vcf_vars}}` variants) calls (Tier-1/2 parity). `VarDictRust` achieved a `{{m:s2_tierp_speedup_min}}`×–`{{m:s2_tierp_speedup_max}}`× speedup and `{{m:s2_tierp_rss_ratio_min}}`×–`{{m:s2_tierp_rss_ratio_max}}`× lower peak RSS (`{{m:s2_tierp_rust_rss_mb_min}}`–`{{m:s2_tierp_rust_rss_mb_max}}` MB vs `{{m:s2_tierp_java_rss_gb_min}}`–`{{m:s2_tierp_java_rss_gb_max}}` GB for `VarDictJava`). Two-sided Mann-Whitney $U$ tests compared runtime ($U = {{m:s2_tierp_mwu_wall_u}}$, $p = {{m:s2_tierp_mwu_wall_p}}$) and RSS ($U = {{m:s2_tierp_mwu_rss_u}}$, $p = {{m:s2_tierp_mwu_rss_p}}$). An exact paired sign-flip permutation gave $p = {{m:s2_tierp_sign_perm_p}}$ across $N = 12$. On the largest `{{m:s2_SRR8657348_intervals}}`-interval panel (`SRR8657348`), `1`-thread `VarDictRust` (`{{m:s2_SRR8657348_th1_rust_wall_s}}` s, range `{{m:s2_SRR8657348_th1_rust_wall_range}}` s) beat both `VarDictJava` (`{{m:s2_SRR8657348_th1_java_wall_s}}` s, `{{m:s2_SRR8657348_th1_speedup_vs_java}}`×) and `vardictcpp` (`{{m:s2_SRR8657348_th1_cpp_wall_s}}` s). At `8` threads (`-th 8`), `vardictcpp` finished the pre-chunked AML panels in `{{m:s2_tierp_cpp_th8_wall_min}}`–`{{m:s2_tierp_cpp_th8_wall_max}}` s vs `{{m:s2_tierp_rust_th8_wall_min}}`–`{{m:s2_tierp_rust_th8_wall_max}}` s for `VarDictRust`. That gap reflects `vardictcpp` hardcoding `--nosv` and skipping `##contig` and `INFO/END` VCF schema validation (`Table 8`; `Supplementary Table S6`).

Replicating Table 3 of Wolff [11] on contiguous human chromosome 1 (`-R chr1:1-249250621`, `{{m:s2_rchr1_mbp}}` Mbp on `{{m:s2_rchr1_sample}}`; `Table 2`, `Figure 2b`, `Supplementary Table S1`) separates tile-bounded streaming (`VarDictRust`) from whole-interval allocation (`VarDictJava` and `vardictcpp`):

- **`VarDictJava`:** Crashed with `java.lang.OutOfMemoryError` at `{{m:s2_rchr1_java4g_oom_wall_s}}` s under a `-Xmx4g` heap (`{{m:s2_rchr1_java4g_oom_rss_gib}}` GiB RSS). Under `-Xmx200g` it required `{{m:s2_rchr1_java200g_wall_s}}` s (`n = 3` range `{{m:s2_rchr1_java200g_wall_range}}` s) and `{{m:s2_rchr1_java200g_rss_gib}}` GiB peak RSS.

- **`vardictcpp`:** Required `{{m:s2_rchr1_cpp_wall_s}}` s (`n = 3` range `{{m:s2_rchr1_cpp_wall_range}}` s) and `{{m:s2_rchr1_cpp_rss_gib}}` GiB peak RSS because a single `-R` region bypasses external BED chunking.

- **`VarDictRust`:** Auto-tiled `-R chr1:1-249250621` into `{{m:rust_default_tile_kbp}}` kbp windows with `{{m:rust_halo_kbp}}` kbp halos. It finished in `{{m:s2_rchr1_rust_wall_s}}` s (`n = 3` range `{{m:s2_rchr1_rust_wall_range}}` s) at `{{m:s2_rchr1_rust_rss_mb}}` MB RSS. That is `{{m:s2_rchr1_speedup_vs_java}}`× faster and `{{m:s2_rchr1_rss_ratio_vs_java}}`× leaner than `VarDictJava` (`-Xmx200g`). It is `{{m:s2_rchr1_speedup_vs_cpp}}`× faster and `{{m:s2_rchr1_rss_ratio_vs_cpp}}`× leaner than `vardictcpp`.

Against external callers (`FreeBayes`, `LoFreq`, `Mutect2`, `Octopus`, `Strelka2`; `Table 2`, `Supplementary Table S1`), `8`-thread `VarDictRust` on `SRR8657348` took `{{m:s2_ss_VarDictRust_th8_wall_s}}` s (`{{m:s2_ss_VarDictRust_th8_rss_mb}}` MB). That beat `FreeBayes` (`{{m:s2_ss_FreeBayes_th8_wall_s}}` s), `LoFreq` (`{{m:s2_ss_LoFreq_th8_wall_s}}` s), `Mutect2` (`{{m:s2_ss_Mutect2_th8_wall_s}}` s), and `Octopus` (`{{m:s2_ss_octopus_th8_wall_s}}` s). On paired somatic `SRR15006376|SRR15006386` (`8` threads), `VarDictRust` finished in `{{m:s2_som_VarDictRust_th8_wall_s}}` s (`{{m:s2_som_VarDictRust_th8_rss_mb}}` MB) vs `{{m:s2_som_Strelka2_th8_wall_s}}` s (`Strelka2`) to `{{m:s2_som_Mutect2_th8_wall_s}}` s (`Mutect2`). It matched `VarDictJava` on all `{{m:s2_venn_VarDict}}` filtered somatic SNVs in the `{{m:s2_venn_union}}`-variant 4-caller Venn union (`Mutect2`: `{{m:s2_venn_Mutect2}}`; `LoFreq`: `{{m:s2_venn_LoFreq}}`; `Octopus`: `{{m:s2_venn_octopus}}`).

<a id="fig2_parity_and_tierp_replication"></a>

![Figure 2: Multi-regime oracle equivalence and AML targeted panel replication.](figures/fig2.svg)

**Figure 2. Multi-regime oracle equivalence and AML targeted panel replication.** **(a)** Exact Tier-1/2 parity sub-arms across the multi-regime sweep cells (`{{m:sweep_match_arms}}/{{m:sweep_total_arms}}` exact match), including the `{{m:s2_modes_total}}`-mode CLI matrix and ultra-deep MRD stress (`Supplementary Table S2`). **(b)** Tile-size invariance on SEQC2 HCC1395: identical calls across all tile sizes (`{{m:s1_tile_vars}}` variants; `Supplementary Table S1`). **(c)** Wall-clock runtime across the Wolff [11] AML panel cohort (`n = 3` technical replicates). **(d)** Somatic tumor-normal runtime and concordant SNV summary across field callers (`FreeBayes`, `LoFreq`, `Mutect2`, `Octopus`, `Strelka2`).

{{table:table2}}

---

### Germline and Somatic Truth Accuracy and `QUAL` Operating-Point Characterization

We evaluated `rtg vcfeval` [12] genotype-level precision, recall, and $F_1$ across germline and somatic truth sets (`Table 3`, `Figure 3`, `Supplementary Table S4`). On the hg19 HG002 chromosome 20 exome benchmark (`Table 3`), `VarDictRust` and `VarDictJava` produced identical SNV (`{{m:s2_fig2_VarDictRust_snv_p}}` precision, `{{m:s2_fig2_VarDictRust_snv_r}}` recall, $F_1 = {{m:s2_fig2_VarDictRust_snv_f1}}$) accuracy. The benchmark spans `{{m:s2_fig2_region_mbp}}` Mbp with `{{m:s2_fig2_parity_vars}}` candidate calls and `{{m:s2_fig2_VarDictRust_calls}}` filtered calls. Indel accuracy was identical too (`{{m:s2_fig2_VarDictRust_indel_p}}` precision, `{{m:s2_fig2_VarDictRust_indel_r}}` recall, $F_1 = {{m:s2_fig2_VarDictRust_indel_f1}}$). They exceeded `Mutect2` (`SNV` $F_1 = {{m:s2_fig2_Mutect2_snv_f1}}$, `Indel` $F_1 = {{m:s2_fig2_Mutect2_indel_f1}}$). They matched `Octopus` (`SNV` $F_1 = {{m:s2_fig2_octopus_snv_f1}}$), `FreeBayes` (`{{m:s2_fig2_FreeBayes_snv_f1}}`), and `LoFreq` (`{{m:s2_fig2_LoFreq_snv_f1}}`).

`VarDictJava` defaults to a low-VAF threshold (`-f {{m:f_std}}`). Unfiltered `PASS` callsets on `35`×–`60`× WGS/exome data therefore retain low-VAF candidates. Stratifying by `QUAL` separates them (`Table 3`, `Supplementary Table S4`, `Figure 3b`). Across all five GIAB and SEQC2 truth cohorts, `VarDictRust` and `VarDictJava` yielded identical `rtg vcfeval` ROC curves:

- **HG002 GRCh38 Slice (`NIST GIAB v4.2.1`):** Optimal `QUAL >= {{m:s4_t3_hg002_best_qual}}` filtering removes `{{m:s4_t3_hg002_fp_reduction_pct}}` of unfiltered false positives (`{{m:s4_t3_hg002_none_fp}}` to `{{m:s4_t3_hg002_best_fp}}` FP). It keeps `{{m:s4_t3_hg002_best_r}}` recall (`{{m:s4_t3_hg002_best_tp}}` TP). Combined $F_1$ rises from `{{m:s4_t3_hg002_none_f1}}` to `{{m:s4_t3_hg002_best_f1}}` (`HG005` $F_1$: `{{m:ev_hg005_none_f1}}` to `{{m:ev_hg005_best_f1}}`; `NA24385-b37` $F_1$: `{{m:ev_b37_none_f1}}` to `{{m:ev_b37_best_f1}}`; `Table 3`, `Supplementary Table S4`).

- **SEQC2 HCC1395/HCC1395BL Somatic Exome (`v1.2.1`):** Across `{{m:s4_strata_seqc2_rust_total}}` candidate calls (`{{m:s4_strata_seqc2_rust_pass}}` `PASS`), `QUAL >= {{m:s4_t3_seqc2_best_qual}}` raises combined somatic $F_1$ from `{{m:s4_t3_seqc2_none_f1}}` to `{{m:s4_t3_seqc2_best_f1}}` (`{{m:s4_t3_seqc2_best_tp}}` TP, `{{m:s4_t3_seqc2_best_fp}}` FP, precision `{{m:s4_t3_seqc2_best_p}}`, recall `{{m:s4_t3_seqc2_best_r}}`). In the native GRCh38 caller benchmark (`Table 3`), `VarDictRust` and `VarDictJava` reached `SNV` $F_1 = {{m:s5_fc_SEQC2_v121_VarDictRust_snv_f1}}$ and `Indel` $F_1 = {{m:s5_fc_SEQC2_v121_VarDictRust_indel_f1}}$. That exceeds `Mutect2` (`SNV` $F_1 = {{m:s5_fc_SEQC2_v121_Mutect2_snv_f1}}$) and `Strelka2` (`{{m:s5_fc_SEQC2_v121_Strelka2_snv_f1}}`) in SNV $F_1$. `VarDictRust` ran in `{{m:s5_fc_SEQC2_v121_VarDictRust_wall_s}}` s (`{{m:s5_fc_SEQC2_v121_VarDictRust_rss_mb}}` MB RSS) vs `{{m:s5_fc_SEQC2_v121_VarDictJava_wall_s}}` s (`VarDictJava`) and `{{m:s5_fc_SEQC2_v121_Mutect2_wall_s}}` s (`Mutect2`).

- **NIST GIAB HG008 Somatic Slice (`v0.3`):** Across `{{m:s4_strata_hg008_rust_total}}` candidates (`{{m:s4_strata_hg008_rust_pass}}` `PASS`), `QUAL >= {{m:s4_t3_hg008_best_qual}}` raises combined $F_1$ from `{{m:s4_t3_hg008_none_f1}}` to `{{m:s4_t3_hg008_best_f1}}` (`{{m:s4_t3_hg008_best_tp}}` TP, `{{m:s4_t3_hg008_best_fp}}` FP). Under standard `StrongSomatic`/`LikelySomatic` `PASS` filtering (`Table 3`), `VarDictRust` and `VarDictJava` achieved identical `SNV` $F_1 = {{m:s5_fc_HG008_v03_VarDictRust_snv_f1}}$ and top `Indel` $F_1 = {{m:s5_fc_HG008_v03_VarDictRust_indel_f1}}$. That outperforms `Mutect2` (`Indel` $F_1 = {{m:s5_fc_HG008_v03_Mutect2_indel_f1}}$), `Strelka2` (`{{m:s5_fc_HG008_v03_Strelka2_indel_f1}}`), and `LoFreq` (`{{m:s5_fc_HG008_v03_LoFreq_indel_f1}}`). `VarDictRust` ran in `{{m:s5_fc_HG008_v03_VarDictRust_wall_s}}` s (`{{m:s5_fc_HG008_v03_VarDictRust_rss_mb}}` MB RSS) vs `{{m:s5_fc_HG008_v03_VarDictJava_wall_s}}` s (`VarDictJava`) and `{{m:s5_fc_HG008_v03_Mutect2_wall_s}}` s (`Mutect2`).

<a id="fig3_truth_qual_and_titration"></a>

![Figure 3: Germline and somatic truth concordance, QUAL operating points, and runtime-memory efficiency.](figures/fig3.svg)

**Figure 3. Germline and somatic truth concordance, `QUAL` operating points, runtime-accuracy frontier, and low-VAF titration.** **(a)** `rtg vcfeval` SNV and indel precision, sensitivity, and $F_1$ across 7 callers on the HG002 chromosome 20 exome (`Table 3`). **(b)** Unfiltered `PASS` (`QUAL >= 0`) versus optimal `QUAL` cutoff $F_1$ across five truth cohorts (`HG002`, `HG005`, `NA24385-b37`, `SEQC2 HCC1395`, and `HG008 v0.3`; `Supplementary Table S4`). **(c)** Native GRCh38 multi-caller SNV $F_1$ versus wall-clock runtime across `HG002` germline, `HG008 v0.3` somatic, and `SEQC2` somatic benchmarks (`Table 3`). **(d)** Controlled *in silico* tumor-purity titration (`20%` down to `1%` purity / `{{m:s6_tit_vaf_20_pct}}` down to `{{m:s6_tit_vaf_1_pct}}` VAF) demonstrating `VarDictRust`'s `{{m:s6_titration_0p5vaf_recall_ratio_vs_m2}}`× recall and `{{m:s6_titration_0p5vaf_f1_ratio_vs_m2}}`× $F_1$ advantage at `{{m:s6_tit_vaf_1_pct}}` VAF (`Supplementary Table S4`).

{{table:table3}}

---

### Bounded Memory Scaling, Multicore Work-Stealing, and Multi-Tenant Node Packing

We tested `VarDictRust`'s `O(threads × (tile + halo))` memory bound across four scaling benchmarks (`Table 4`, `Figure 4`, `Supplementary Table S3`):

1. **Four-Assay Performance Gate (`Table 4`, `Supplementary Table S3`):** Across `n = 3` technical replicates, `VarDictRust` held wall-time ratios to `{{m:s4_pb_amplicon_tier2_wall_ratio}}`–`{{m:s4_pb_hg008_tier2_wall_ratio}}`. RSS ratios held to `{{m:s4_pb_hg008_tier2_rss_ratio}}`–`{{m:s4_pb_amplicon_tier2_rss_ratio}}` relative to `VarDictJava`. Targets were EGFR/MRE11A amplicon (`{{m:s4_pb_amplicon_tier2_rust_wall_s}}` s vs `{{m:s4_pb_amplicon_tier2_java_wall_s}}` s Java) and HG002 chr20 germline (`{{m:s4_pb_hg002_tier2_rust_wall_s}}` s vs `{{m:s4_pb_hg002_tier2_java_wall_s}}` s). They included HG008 somatic (`{{m:s4_pb_hg008_tier2_rust_wall_s}}` s vs `{{m:s4_pb_hg008_tier2_java_wall_s}}` s) and SEQC2 somatic exome (`{{m:s4_pb_seqc2_exome_tier2_rust_wall_s}}` s vs `{{m:s4_pb_seqc2_exome_tier2_java_wall_s}}` s).

2. **Whole-Exome Multicore Scaling (`1` to `32` Threads, `Table 4`, `Figure 4c`):** We used the SEQC2 HCC1395 whole exome (`{{m:s4_exome_total_vars}}` total variants, `{{m:s4_exome_vars}}` `PASS`, `n = 3` technical replicates). `VarDictRust` scaled from `{{m:s4_exome_th1_rust_wall_s}}` s (`{{m:s4_exome_th1_rust_rss_mb}}` MB RSS) at `1` thread. At `8` threads it reached `{{m:s4_exome_th8_rust_wall_s}}` s (`{{m:s4_exome_th8_rust_sp}}`× self-speedup). At `32` threads it reached `{{m:s4_exome_th32_rust_wall_s}}` s (`{{m:s4_exome_th32_rust_sp}}`× self-speedup, `{{m:s4_exome_th32_rust_rss_mb}}` MB RSS). `VarDictJava` plateaued above `8` threads (`{{m:s4_exome_th8_java_wall_s}}` s at `8` threads; `{{m:s4_exome_th32_java_wall_s}}` s at `32` threads). Its heap grew from `{{m:s4_exome_th1_java_rss_gib}}` GiB to `{{m:s4_exome_th32_java_rss_gib}}` GiB. At `32` threads `VarDictRust` is `{{m:s4_exome_th32_speedup_vs_java}}`× faster and `{{m:s4_exome_th32_rss_ratio_vs_java}}`× leaner than `VarDictJava`. It is `{{m:s4_exome_th32_rss_ratio_vs_cpp}}`× leaner than `vardictcpp` at `{{m:s4_exome_th32_cpp_rss_gib}}` GiB.

3. **Chromosome 1 Span Invariance (`Table 4`, `Figure 4d`):** Contiguous chromosome 1 span grew `250`× across `{{m:s4_span_1m}}`, `{{m:s4_span_10m}}`, `{{m:s4_span_100m}}`, and `{{m:s4_span_249m}}` bp (`1` to `249` Mbp). `VarDictRust`'s `249` Mbp-to-`10` Mbp RSS ratio stayed below `{{m:s4_span_h2_bound}}`× at every thread count. Ratios were `{{m:s4_span_th8_ratio_249m_10m}}`× (`{{m:s4_span_th8_249m_rust_gib}}` GiB) at `8` threads, `{{m:s4_span_th16_ratio_249m_10m}}`× (`{{m:s4_span_th16_249m_rust_gib}}` GiB) at `16` threads, and `{{m:s4_span_th32_ratio_249m_10m}}`× (`{{m:s4_span_th32_249m_rust_gib}}` GiB) at `32` threads. On `10` Mbp all-positions `-f {{m:f_zero}}` pileup (`{{m:s4_span_f0_vars}}` variants), `VarDictRust` matched `VarDictJava` in `{{m:s4_span_f0_rust_wall_s}}` s (`{{m:s4_span_f0_rust_gib}}` GiB RSS).

4. **Concurrent Node Packing (`Table 4`, `Figure 4e`):** Two concurrent `4`-core SEQC2 exome jobs finished in `{{m:s4_jpn_job1_wall_s}}` s (`{{m:s4_jpn_job1_rss_mb}}` MB) and `{{m:s4_jpn_job2_wall_s}}` s (`{{m:s4_jpn_job2_rss_mb}}` MB). That matches an isolated `4`-core run (`{{m:s4_exome_th4_rust_wall_s}}` s). Combined RSS was `{{m:s4_jpn_combined_gib}}` GiB. A `128` GiB node can pack `8` concurrent `4`-core `VarDictRust` exome jobs in `{{m:s4_jpn_8jobs_rust_gib}}` GiB total RAM versus `{{m:s4_jpn_8jobs_java_gib}}` GiB for `VarDictJava`.

<a id="fig4_memory_and_multicore_scaling"></a>

![Figure 4: Whole-exome multicore scaling, chromosome 1 span-to-RSS invariance, and node packing.](figures/fig4.svg)

**Figure 4. Whole-exome multicore scaling, chromosome 1 span-to-RSS invariance, and node packing.** **(a, b)** Peak memory across germline, somatic, `{{m:s3_rna_bam_gib}}` GiB RNA-seq, and contiguous chromosome 1 (`{{m:s2_rchr1_mbp}}` Mbp) workloads (`Table 4`). **(c)** SEQC2 HCC1395 whole-exome (`{{m:s4_exome_total_vars}}` variants) wall-clock time and peak RSS from `1` to `32` threads (`n = 3` technical replicates; `Supplementary Table S3`). **(d)** Chromosome 1 contiguous span scaling (`1` Mbp to `249` Mbp) demonstrating `VarDictRust`'s `O(tile + halo)` memory bound (`{{m:s4_span_th8_ratio_249m_10m}}`×–`{{m:s4_span_th32_ratio_249m_10m}}`× ratio over a `250`× span increase). **(e)** Concurrent multi-tenant node-packing efficiency (`2 × 4` cores, `{{m:s4_jpn_combined_gib}}` GiB combined peak RSS) and chromosome 1 `-f {{m:f_zero}}` all-positions stress.

{{table:table4}}

---

### Full-Genome WGS Scaling: Unchunked Full Chromosomes Versus Pre-Chunked Callable Intervals

We compared all three engines at `32` threads across two GRCh38 whole genomes (`Table 5`, `Figure 5`). They were `{{m:s5_hg002_bam_gb}}` GB NIST HG002 `2×250` bp germline WGS and `{{m:s5_hg008_bam_gb}}` GB NIST HG008 `364`× (`195`× tumor / `169`× normal) paired somatic WGS. Each genome was run under two region protocols:

- **Unchunked Full Chromosomes (`chr1`–`chr22`, `22` autosomes):** On HG002, `VarDictJava` (`-Xmx200g`) aborted with `java.lang.OutOfMemoryError` at `{{m:s5_hg002_contigs_java_rss_gb}}` GB RAM. `vardictcpp` finished in `{{m:s5_hg002_contigs_cpp_wall_str}}` at `{{m:s5_hg002_contigs_cpp_rss_gb}}` GB RSS (`{{m:s5_hg002_contigs_cpp_rss_ratio_vs_rust}}`× higher than `VarDictRust`). It emitted `{{m:s5_hg002_contigs_cpp_vars}}` variants (`+{{m:s5_hg002_contigs_cpp_excess}}` duplicate chunk-boundary calls relative to `VarDictRust`'s `{{m:s5_hg002_contigs_rust_vars}}`). `VarDictRust` finished all `22` autosomes in `{{m:s5_hg002_contigs_rust_wall_str}}` at `{{m:s5_hg002_contigs_rust_rss_gb}}` GB RSS. On `364`× HG008 paired somatic WGS, `VarDictJava` crashed out of memory at `{{m:s5_hg008_contigs_java_rss_gb}}` GB. `vardictcpp` reached `{{m:s5_hg008_contigs_cpp_rss_gb}}` GB RSS (`{{m:s5_hg008_contigs_cpp_rss_ratio_vs_rust}}`× higher than `VarDictRust`) before timing out at `3` hours. `VarDictRust` finished in `{{m:s5_hg008_contigs_rust_wall_str}}` at `{{m:s5_hg008_contigs_rust_rss_gb}}` GB RSS (`{{m:s5_hg008_contigs_rust_vars}}` variants). These unchunked runs are completion results: neither Java nor C++ finished them.

- **Pre-Chunked Callable Intervals (`≤{{m:giab_arm_b_max_kbp}}` kbp GIAB BED):** On pre-split HG002 intervals (`{{m:s5_hg002_giab_gbp}}` Gbp), `VarDictRust` finished in `{{m:s5_hg002_giab_rust_wall_str}}` at `{{m:s5_hg002_giab_rust_rss_gb}}` GB RSS. That is `{{m:s5_hg002_giab_speedup_vs_java}}`× faster and `{{m:s5_hg002_giab_rss_ratio_vs_java}}`× leaner than `VarDictJava` (`{{m:s5_hg002_giab_java_wall_str}}`, `{{m:s5_hg002_giab_java_rss_gb}}` GB). It is `{{m:s5_hg002_giab_rss_ratio_vs_cpp}}`× leaner than `vardictcpp` at `{{m:s5_hg002_giab_cpp_rss_gb}}` GB in `{{m:s5_hg002_giab_cpp_wall_str}}`. `VarDictRust` matched `VarDictJava` across all `{{m:s5_hg002_giab_vars}}` variants (`vardictcpp` differed by `{{m:s5_hg002_giab_cpp_diff_vars}}` at `{{m:s5_hg002_giab_cpp_vars}}`). It reached `rtg vcfeval` $F_1 = {{m:s5_hg002_wgs_f1_best}}$ (`{{m:s5_hg002_wgs_p_best}}` precision, `{{m:s5_hg002_wgs_r_best}}` recall; unfiltered $F_1 = {{m:s5_hg002_wgs_f1_none}}$). On pre-split HG008 intervals (`{{m:s5_hg008_giab_gbp}}` Gbp), `VarDictJava` (`{{m:s5_hg008_giab_java_rss_gb}}` GB) timed out at `4` hours (`{{m:s5_hg008_giab_java_wall_str}}`). `vardictcpp` (`{{m:s5_hg008_giab_cpp_rss_gb}}` GB) timed out at `3` hours (`{{m:s5_hg008_giab_cpp_wall_str}}`). `VarDictRust` finished in `{{m:s5_hg008_giab_rust_wall_str}}` at `{{m:s5_hg008_giab_rust_rss_gb}}` GB RSS (`{{m:s5_hg008_giab_rust_vars}}` variants; somatic $F_1 = {{m:s5_hg008_wgs_v03_f1_best}}$ on `GIAB v0.3` and `{{m:s5_hg008_wgs_v04_f1_best}}` on `GIAB v0.4`). This HG008 run is completion-only: no oracle finished it.

<a id="fig5_full_wgs_arm_a_vs_arm_b"></a>

![Figure 5: Full-genome WGS scaling across unchunked full chromosomes and pre-chunked callable intervals.](figures/fig5.svg)

**Figure 5. Full-genome WGS scaling across unchunked full chromosomes (`chr1`–`chr22`) and pre-chunked callable intervals (`≤{{m:giab_arm_b_max_kbp}}` kbp GIAB BED).** **(a)** Peak RSS (GB) and completion status on `{{m:s5_hg002_bam_gb}}` GB HG002 germline WGS and `{{m:s5_hg008_bam_gb}}` GB `364`× HG008 paired somatic WGS across unchunked chromosomes and pre-chunked intervals (`32` threads; `Table 5`). **(b)** Wall-clock runtime comparison showing `VarDictJava` out-of-memory and wall-clock timeout failures alongside `VarDictRust` completion (`{{m:s5_hg002_contigs_rust_wall_str}}` and `{{m:s5_hg008_contigs_rust_wall_str}}`). **(c)** Full-genome `rtg vcfeval` precision, sensitivity, and $F_1$ across `HG002 GIAB v4.2.1` (optimal `QUAL` and all-`PASS`) and `HG008 GIAB v0.3`/`v0.4` (`Table 5`).

{{table:table5}}

---

### Cross-Platform Sequencing Chemistry Generalization and 70-Cell GA4GH v3.5 Stratification

We evaluated `{{m:s5_cross_platform_arms}}` short-read sequencing chemistries across the full `{{m:s5_chr20_mbp}}` Mbp human chromosome 20 (`Table 6`, `Figure 6`, `Supplementary Table S4`, `Supplementary Table S5`). The panel holds `4` HG002 germline chemistries (`Illumina`, `Element AVITI`, `MGI DNBSEQ`, `Ultima Genomics UG100`). They were evaluated on standard GIAB and Ultima High-Confidence Region [`UG-HCR`] masks. It also holds `6` HG008 paired somatic chemistries (`Illumina`, `Element AVITI`, `PacBio Onso`, `Roche SBX`, `Ultima Genomics ppmSeq`, `Ultima Genomics UG100`).

`VarDictRust` completed all `{{m:s5_chr20_total_runs}}` standard (`-f {{m:f_std}}`) and all-frequency (`-f {{m:f_zero}}`) chromosome 20 runs (`100%` completion). `VarDictJava` completed `{{m:s5_chr20_java_completed_runs}}/{{m:s5_chr20_total_runs}}` runs and timed out on `{{m:s5_chr20_java_timeout_runs}}` deep runs (`Table 6`). On all `{{m:s5_chr20_java_completed_runs}}` runs where `VarDictJava` finished, `VarDictRust` achieved `100%` (`95%` CI `{{m:cp95_chr20_java_completed}}`) exact Tier-1/2 parity. SNV $F_1$ was identical across HG002 germline (`{{m:s5_chr20_hg2_std_snv_f1_range}}` standard; `{{m:s5_chr20_hg2_ug_snv_f1_range}}` `UG-HCR`). It was identical across HG008 somatic (`{{m:s5_chr20_hg8_f01_std_snv_f1_range}}` standard and `{{m:s5_chr20_hg8_f01_ug_snv_f1_range}}` `UG-HCR` at `-f {{m:f_std}}`; `{{m:s5_chr20_hg8_all_std_snv_f1_range}}` and `{{m:s5_chr20_hg8_all_ug_snv_f1_range}}` across all frequencies). The unfiltered somatic noise-floor sweep ran `18` runs (`-f {{m:f_zero}}`, `{{m:s5_nf_total_candidates}}` candidates, `{{m:s5_nf_total_pass}}` `PASS`; `Supplementary Table S4`). `VarDictRust` matched `VarDictJava` byte-for-byte on `{{m:s5_f0_noise_floor_strict_match}}` runs (`100%`, `95%` CI `{{m:cp95_noise_18_18}}`).

Crossing all `{{m:s5_cross_platform_arms}}` chemistries with `7` GA4GH `v3.5` genomic context strata (`4–6 bp homopolymers`, `7–11 bp homopolymers`, `≥12 bp homopolymers`, `tandem repeats`, `segmental duplications`, `low-mappability regions`, and `RefSeq coding exons (CDS)`) yielded a `{{m:s6_ga4gh_cells}}`-cell stratification matrix. All `{{m:s6_ga4gh_cells}}` Rust truth evaluations completed (`95%` CI `{{m:cp95_70_70}}`, mean SNV $F_1 = {{m:s6_ga4gh_mean_snv_f1}}$; `Supplementary Table S5`, `Figure 6c`); Java–Rust parity held at `100%` on the `{{m:s5_chr20_java_completed_runs}}/{{m:s5_chr20_java_completed_runs}}` completing paired chromosome 20 arms (per-arm rows in `Table 6`).

<a id="fig6_cross_platform_and_ga4gh"></a>

![Figure 6: Cross-platform short-read sequencing chemistry panel and 70-cell GA4GH v3.5 stratification heatmap.](figures/fig6.svg)

**Figure 6. Cross-platform short-read sequencing chemistry panel and `{{m:s6_ga4gh_cells}}`-cell GA4GH `v3.5` stratification heatmap.** **(a)** Completion rate (`VarDictRust` `{{m:s5_chr20_total_runs}}/{{m:s5_chr20_total_runs}}` vs `VarDictJava` `{{m:s5_chr20_java_completed_runs}}/{{m:s5_chr20_total_runs}}`) and accuracy across `{{m:s5_cross_platform_arms}}` germline and somatic chemistries on chromosome 20 (`{{m:s5_chr20_mbp}}` Mbp; `Table 6`). **(b)** Byte-equivalent candidate and `PASS` somatic variant counts (`{{m:s5_f0_noise_floor_strict_match}}` exact match) across the unfiltered somatic noise-floor sweep (`-f {{m:f_zero}}`; `Supplementary Table S4`). **(c)** `{{m:s6_ga4gh_cells}}`-cell (`10` chemistries × `7` GA4GH `v3.5` strata) SNV $F_1$ heatmap (`Supplementary Table S5`).

{{table:table6}}

---

### Unlocking GIAB v5.0q T2T-Q100, CMRG v1.00, and Mosaic v1.10 Discovery via Unchunked Calling

We compared full-genome HG002 callsets from unchunked chromosomes (`chr1`–`chr22`) versus pre-chunked `≤{{m:giab_arm_b_max_kbp}}` kbp `GIAB v4.2.1` intervals across three NIST truth benchmarks (`Table 7`, `Figure 7`):

1. **`NIST GIAB v5.0q T2T-Q100` (`{{m:s6_v50q_gbp}}` Gbp, diploid Telomere-to-Telomere Q100 assembly truth):** Because `VarDictJava` runs out of memory on unchunked chromosomes, legacy workflows restrict calling to pre-chunked `GIAB v4.2.1` intervals. Those cover `{{m:s5_hg002_giab_gbp}}` Gbp and miss `+{{m:s6_v50q_gain_gbp}}` Gbp of segmental duplications, pericentromeric satellites, and complex loci added in `v5.0q`. On `GIAB v4.2.1`, unchunked and pre-chunked calling matched (`SNV` recall `{{m:s6_NIST_v4.2.1_All_VarDictRust_FullContigs_snv_r}}` vs `{{m:s6_NIST_v4.2.1_All_VarDictRust_GIAB50kb_snv_r}}`, `SNV` $F_1 = {{m:s6_NIST_v4.2.1_All_VarDictRust_FullContigs_snv_f1}}$, `Indel` $F_1 = {{m:s6_NIST_v4.2.1_All_VarDictRust_FullContigs_indel_f1}}$). On `GIAB v5.0q T2T-Q100`, streaming unchunked full chromosomes with `VarDictRust` raised genome-wide SNV recall by `+{{m:s6_v50q_snv_r_gain_pct}}` (`{{m:s6_NIST_v5.0q_T2T_Q100_VarDictJava_GIAB50kb_snv_r}}` to `{{m:s6_v50q_snv_r_fullcontigs}}`, $F_1$ from `{{m:s6_NIST_v5.0q_T2T_Q100_VarDictJava_GIAB50kb_snv_f1}}` to `{{m:s6_v50q_snv_f1_fullcontigs}}`). Indel recall rose by `+{{m:s6_v50q_indel_r_gain_pct}}` (`{{m:s6_NIST_v5.0q_T2T_Q100_VarDictJava_GIAB50kb_indel_r}}` to `{{m:s6_v50q_indel_r_fullcontigs}}`, $F_1$ from `{{m:s6_NIST_v5.0q_T2T_Q100_VarDictJava_GIAB50kb_indel_f1}}` to `{{m:s6_v50q_indel_f1_fullcontigs}}`).

2. **`NIST Challenging Medically Relevant Genes` (`CMRG v1.00`, `273` clinical genes):** Many of the `273` CMRG genes lie outside the `v4.2.1` `50` kbp mask. Unchunked `VarDictRust` calling increased CMRG SNV recall by `+{{m:s6_cmrg_snv_r_gain_pct}}` (`{{m:s6_NIST_CMRG_v1.00_273Genes_VarDictJava_GIAB50kb_snv_r}}` to `{{m:s6_cmrg_snv_r_fullcontigs}}`, $F_1$ from `{{m:s6_NIST_CMRG_v1.00_273Genes_VarDictJava_GIAB50kb_snv_f1}}` to `{{m:s6_cmrg_snv_f1_fullcontigs}}`). Indel recall rose by `+{{m:s6_cmrg_indel_r_gain_pct}}` (`{{m:s6_NIST_CMRG_v1.00_273Genes_VarDictJava_GIAB50kb_indel_r}}` to `{{m:s6_cmrg_indel_r_fullcontigs}}`, $F_1$ from `{{m:s6_NIST_CMRG_v1.00_273Genes_VarDictJava_GIAB50kb_indel_f1}}` to `{{m:s6_cmrg_indel_f1_fullcontigs}}`).

3. **`NIST Mosaic SNV v1.10` (`{{m:s6_mosaic_truth_total}}` validated low-VAF mosaic SNVs in HG002, VAF `{{m:s6_mosaic_vaf_range_pct}}`):** `VarDictRust` and `VarDictJava` detected `{{m:s6_mosaic_rust_tp}}/{{m:s6_mosaic_truth_total}}` (`{{m:s6_mosaic_rust_r_pct}}`, `95%` CI `{{m:cp95_84_85}}`) validated mosaic variants (`Table 7`).

<a id="fig7_nextgen_truth_v50q_cmrg_mosaic"></a>

![Figure 7: Next-generation GIAB v5.0q T2T-Q100, CMRG v1.00, Mosaic v1.10, and HG008 v0.4 truth evaluation.](figures/fig7.svg)

**Figure 7. Next-generation `GIAB v5.0q T2T-Q100`, `CMRG v1.00`, `Mosaic v1.10`, and `HG008 v0.4` truth evaluation.** **(a)** Genome-wide and medically relevant gene SNV sensitivity unlocked by unchunked full-chromosome calling (`VarDictRust`) versus legacy `≤{{m:giab_arm_b_max_kbp}}` kbp `GIAB v4.2.1` interval calling on `NIST v4.2.1`, `NIST v5.0q T2T-Q100` (`+{{m:s6_v50q_snv_r_gain_pct}}` SNV recall), and `CMRG v1.00` (`+{{m:s6_cmrg_snv_r_gain_pct}}` SNV recall; `Table 7`). **(b)** `NIST Mosaic v1.10` sensitivity (`{{m:s6_mosaic_rust_tp}}/{{m:s6_mosaic_truth_total}}`, `{{m:s6_mosaic_rust_r_pct}}`) and `CMRG v1.00` indel recovery (`+{{m:s6_cmrg_indel_r_gain_pct}}` indel recall; `Table 7`). **(c)** `NIST GIAB HG008 v0.4` somatic SNV and indel $F_1$ across `{{m:s6_v04_total_callers}}` callers (`Table 7`).

{{table:table7}}

---

### Deep-Learning and Commercial Caller Comparison, Low-VAF Titration, and Supply-Chain Hardening

We compared `VarDictRust` against contemporary deep-learning, graph, and haplotype somatic callers. Three benchmarks cover them (`Table 3`, `Table 7`, `Table 8`, `Figure 3d`, `Figure 7c`, `Figure 8`, `Supplementary Table S6`):

- **`{{m:s6_v04_total_callers}}`-Caller Benchmark on `NIST GIAB HG008 v0.4` (`{{m:s6_v04_total_rows}}` evaluations; `Table 7`, `Figure 7c`):** We compared `VarDictRust` (`v{{m:version}}`) against official NIST GIAB and local callsets. Callers were `DRAGEN v4.2.4` (`Graph` and `Linear`), `DeepSomatic v1.6.0`, `ClairS v0.4.1`, `Lancet2` (`NYGC`), `Hartwig SAGE` [25], `Ultima DeepVariant Somatic`, `Mutect2`, `Strelka2`, and `LoFreq`. Scopes were full-genome (`All` and `NoGI`) and slice. On full-genome HG008 `v0.4` (`All`), `VarDictRust` reached `{{m:s6_HG008_FullWGS_v04_All_VarDictRust_v106_indel_r_pct}}` somatic indel recall (`Indel` $F_1 = {{m:s6_HG008_FullWGS_v04_All_VarDictRust_v106_indel_f1}}$, `SNV` $F_1 = {{m:s6_HG008_FullWGS_v04_All_VarDictRust_v106_snv_f1}}$). On `NoGI` it reached `{{m:s6_HG008_FullWGS_v04_NoGI_VarDictRust_v106_indel_r_pct}}` indel recall (`Indel` $F_1 = {{m:s6_HG008_FullWGS_v04_NoGI_VarDictRust_v106_indel_f1}}$, `SNV` $F_1 = {{m:s6_HG008_FullWGS_v04_NoGI_VarDictRust_v106_snv_f1}}$). This exceeded `ClairS v0.4.1` (`{{m:s6_HG008_FullWGS_v04_All_ClairS_v041_indel_r_pct}}` indel recall, $F_1 = {{m:s6_HG008_FullWGS_v04_All_ClairS_v041_indel_f1}}$) and `DeepSomatic v1.6.0` (`{{m:s6_HG008_FullWGS_v04_All_DeepSomatic_v160_indel_r_pct}}`, $F_1 = {{m:s6_HG008_FullWGS_v04_All_DeepSomatic_v160_indel_f1}}$). It also exceeded `Lancet2` (`{{m:s6_HG008_FullWGS_v04_All_Lancet2_NYGC_indel_r_pct}}`, $F_1 = {{m:s6_HG008_FullWGS_v04_All_Lancet2_NYGC_indel_f1}}$) and `Ultima DeepVariant Somatic` (`{{m:s6_HG008_FullWGS_v04_All_Ultima_DeepVariant_Somatic_indel_r_pct}}`, $F_1 = {{m:s6_HG008_FullWGS_v04_All_Ultima_DeepVariant_Somatic_indel_f1}}$). It approached FPGA-accelerated `DRAGEN v4.2.4 Graph` (`{{m:s6_HG008_FullWGS_v04_All_DRAGEN_v424_Graph_indel_r_pct}}`) on standard CPUs.

- **Controlled *in silico* Tumor-Purity Titration (`20%` to `1%` Purity / `{{m:s6_tit_vaf_range_pct}}` VAF; `Supplementary Table S4`, `Figure 3d`):** We built five HG008 tumor-in-normal admixtures. Purities were `20%`, `10%`, `5%`, `2%`, and `1%` with expected clonal VAF `{{m:s6_tit_vaf_20_pct}}`, `{{m:s6_tit_vaf_10_pct}}`, `{{m:s6_tit_vaf_5_pct}}`, `{{m:s6_tit_vaf_2_pct}}`, and `{{m:s6_tit_vaf_1_pct}}`. Each was called at `-f {{m:f_titration}}`. `VarDictRust` matched `VarDictJava` at every purity. It ran `{{m:s6_tit_speedup_range}}` faster (`{{m:s6_tit_rust_wall_range}}` s vs `{{m:s6_tit_java_wall_range}}` s) at `{{m:s6_tit_rss_ratio_range}}` lower RSS (`{{m:s6_tit_rust_rss_kb_range}}` KB vs `{{m:s6_tit_java_rss_kb_range}}` KB). At `5%` purity (`{{m:s6_tit_vaf_5_pct}}` VAF), `VarDictRust` recalled `{{m:s6_tit_p5_VarDictRust_snv_r_pct}}` of SNVs ($F_1 = {{m:s6_tit_p5_VarDictRust_snv_f1}}$). Indel recall was `{{m:s6_tit_p5_VarDictRust_indel_r_pct}}` ($F_1 = {{m:s6_tit_p5_VarDictRust_indel_f1}}$) vs `{{m:s6_tit_p5_Strelka2_snv_r_pct}}` (`Strelka2`), `{{m:s6_tit_p5_Mutect2_snv_r_pct}}` (`Mutect2`), and `{{m:s6_tit_p5_LoFreq_snv_r_pct}}` (`LoFreq`). At `1%` purity (`{{m:s6_tit_vaf_1_pct}}` VAF), `VarDictRust` retained `{{m:s6_tit_p1_VarDictRust_snv_r_pct}}` SNV recall ($F_1 = {{m:s6_titration_0p5vaf_rust_snv_f1}}$). That is `{{m:s6_titration_0p5vaf_recall_ratio_vs_m2}}`× higher recall and `{{m:s6_titration_0p5vaf_f1_ratio_vs_m2}}`× higher $F_1$ than `Mutect2` (`{{m:s6_tit_p1_Mutect2_snv_r_pct}}`, $F_1 = {{m:s6_titration_0p5vaf_m2_snv_f1}}$) and `Strelka2` (`{{m:s6_tit_p1_Strelka2_snv_r_pct}}`).

- **Supply-Chain Hardening and `vardictcpp` Audit (`Table 8`, `Figure 8`, `Supplementary Table S6`):** Compiling `VarDictRust` to a `{{m:s5_h4_musl_mb}}` MB static Linux binary (`{{m:s5_h4_gnu_mb}}` MB dynamic `glibc`) removes `{{m:s5_h4_conda_pkgs}}` runtime Conda packages (`OpenJDK`, `R`, `Perl`). That resolves all ten `vardictcpp` functional and statistical divergences (`Table 8`).

<a id="fig8_modern_callers_and_low_vaf"></a>

![Figure 8: Software supply-chain footprint, multicore scaling, and functional divergences in vardictcpp.](figures/fig8.svg)

**Figure 8. Software supply-chain footprint, multi-thread output determinism, and empirical variant discrepancies in `vardictcpp`.** **(a)** Deployment artifact footprint (`{{m:s5_h4_gnu_mb}}` MiB `glibc`, `{{m:s5_h4_musl_mb}}` MiB static `musl`, `0` runtime packages, `0` CVEs, SLSA Level 3) versus `VarDictJava` (`{{m:s5_h4_conda_pkgs}}` Conda packages) and `vardictcpp` (`Table 8`). **(b)** Multi-thread output determinism across `1`, `2`, `4`, `8`, `16`, and `32` threads on the SEQC2 HCC1395 whole exome (`{{m:s4_exome_vars}}` emitted variants at every thread count, `0` diffs vs `VarDictJava`; `Table 4`). **(c)** Empirical variant-count discrepancies in `vardictcpp` relative to the `VarDictJava` oracle across unchunked (`+{{m:s5_hg002_contigs_cpp_excess}}` `50` kbp chunk-boundary duplicates) and pre-chunked (`{{m:s5_hg002_giab_cpp_diff_vars}}` MNP/tie-break divergences) HG002 WGS, HG002 chromosome 20, SEQC2 somatic exome, and HG008 somatic slice (`Table 5`, `Table 8`, `Supplementary Table S6`).

{{table:table8}}

---

## Discussion

### Scientific Software Consolidation Under the `rewrites.bio` Specification
Rewriting a clinical variant caller for speed risks changing Smith-Waterman tie-breaking, Phred quality offsets, floating-point rounding, or strand-bias p-values, which forces clinical laboratories to re-validate pipelines. Conversely, line-by-line transliteration preserves memory growth and crash bugs. `VarDictRust` separates **v1.0.x oracle parity and crash/memory robustness** from **post-1.0 opt-in algorithmic corrections (`--fix <id>`)** under the [rewrites.bio](https://rewrites.bio) specification. In default mode (`v{{m:version}}`), `VarDictRust` reproduces the exact allele calls, `FILTER` flags, somatic `STATUS` classifications, and `HALF_EVEN`-rounded Fisher exact statistics of `{{m:oracle_version}}` (`Table 1`, `Supplementary Table S2`), while eliminating JVM `OutOfMemoryError` failures, `CigarModifier` bounds panics (`#![forbid(unsafe_code)]`), and `java.util.HashMap` row-order non-determinism (`#282`). Upstream output-affecting quirks in `VarDictJava`—including `VD > DP` under complex realignment (`#92`), aggressive `NM` filtering (`#62`), `SHIFT3` coordinate shifts (`#292`), and `<dup-N>` symbolic allele naming (`#342`)—are preserved in `v1.0.x` and gated behind opt-in `--fix` flags for `v1.1+` (`Table 8`).

### Impact of Bounded Streaming on T2T-Q100 and CMRG Discovery
Whole-interval heap allocation in `VarDictJava` and `vardictcpp` has historically restricted whole-genome `VarDict` calling to pre-chunked `≤{{m:giab_arm_b_max_kbp}}` kbp BED masks (`GIAB v4.2.1`) or external `{{m:cpp_chunk_kbp}}` kbp BED-splitting scripts (`Table 5`). Keeping intervals $\le {{m:s1_tile_1m}}$ bp untiled and partitioning spans $> {{m:s1_tile_1m}}$ bp into `{{m:rust_default_tile_kbp}}` kbp tiles with `{{m:rust_halo_kbp}}` kbp realignment halos (`Table 1`) allows `VarDictRust` to call unchunked full chromosomes in `{{m:s5_hg002_contigs_rust_rss_gb}}`–`{{m:s5_hg008_contigs_rust_rss_gb}}` GB RAM. Removing the `50` kbp BED mask recovers `+{{m:s6_v50q_gain_gbp}}` Gbp of segmental duplications and complex loci in `NIST GIAB v5.0q T2T-Q100` (`+{{m:s6_v50q_snv_r_gain_pct}}` SNV recall, `+{{m:s6_v50q_indel_r_gain_pct}}` indel recall) and raises recall across the `273` `NIST CMRG v1.00` genes (`+{{m:s6_cmrg_snv_r_gain_pct}}` SNV, `+{{m:s6_cmrg_indel_r_gain_pct}}` indel; `Table 7`).

### Trade-Offs Against Deep-Learning, Commercial, and C++ Callers
On `364`× `HG008 v0.4` somatic WGS (`Table 7`, `Figure 7c`), deep-learning (`DeepSomatic`, `ClairS`, `Ultima DeepVariant Somatic`) and graph/haplotype (`DRAGEN`, `Lancet2`, `Hartwig SAGE`, `Mutect2`, `Strelka2`) callers achieve higher SNV precision in repetitive sequence. `VarDictRust`, however, achieves higher somatic indel recall (`{{m:s6_HG008_FullWGS_v04_All_VarDictRust_v106_indel_r_pct}}`) than `ClairS`, `DeepSomatic`, `Lancet2`, and `Ultima DeepVariant Somatic`. It retains `{{m:s6_titration_0p5vaf_recall_ratio_vs_m2}}`× higher SNV recall than `Mutect2` and `Strelka2` at `{{m:s6_tit_vaf_1_pct}}` VAF (`1%` purity; `Supplementary Table S4`) on commodity CPUs without GPUs or FPGAs. On pre-chunked `{{m:s2_fig1_fig3_wgs_bed_intervals}}`-interval panels (`SRR8657348`, `8` threads), `vardictcpp` finishes raw `.var` calling in `{{m:s2_SRR8657348_th8_cpp_wall_s}}` s vs `{{m:s2_SRR8657348_th8_rust_wall_s}}` s for `VarDictRust`. That gap reflects `vardictcpp` hardcoding `--nosv` and skipping `##contig` and `INFO/END` VCF validation (`Table 8`). `VarDictRust` instead completes every unchunked chromosome run and all `{{m:s5_chr20_total_runs}}` deep chromosome 20 runs with full SV and VCF output.

### Operational Boundary Conditions
Three settings govern production deployment:

1. **Multi-Megabase Single-Region Calling (`--tile-size` vs. `{{m:rust_default_tile_kbp}}` kbp Auto-Tiling):** Intervals $\le {{m:s1_tile_1m}}$ bp run untiled by default. On unchunked $> {{m:s1_tile_1m}}$ bp regions, default `{{m:rust_default_tile_kbp}}` kbp auto-tiling bounds `-R chr1:1-249250621` memory to `{{m:s2_rchr1_rust_rss_mb}}` MB (`Table 2`). It also prevents `VarDictJava`'s cross-centromere `36`–`81` Mbp SV pairing and `VD = 74 > DP` state leak (`#92`). Passing `--tile-size` equal to the region span reproduces `VarDictJava`'s single-window output byte-for-byte. Use `--tile-size 10000000` on `10` Mbp `chr1` `-f {{m:f_zero}}`. Use `--tile-size 250000000` on `{{m:s2_rchr1_mbp}}` Mbp `-R chr1` (`{{m:sw_span_f0_rows}}` `.var` / `{{m:s4_span_f0_vars}}` `.vcf` and `{{m:sw_rchr1_rows}}` `.var` rows; `Table 1`, `Supplementary Table S2`).

2. **`QUAL` Filtering on Standard-Depth Germline WGS:** Because `-f {{m:f_std}}` targets low-VAF somatic discovery, `35`×–`60`× germline WGS pipelines need filtering. Filter at `QUAL >= {{m:s4_t3_hg002_best_qual}}` ($F_1 = {{m:s4_t3_hg002_best_f1}}$ vs `{{m:s4_t3_hg002_none_f1}}` unfiltered, removing `{{m:s4_t3_hg002_fp_reduction_pct}}` of false positives; `Table 3`, `Supplementary Table S4`). Alternatively set `-f` to `0.10`–`0.20`.

3. **Strict `v1.0.x` Parity Policy:** All `v1.0.x` releases preserve upstream `VarDictJava` v1.8.3 output quirks (`#62`, `#92`, `#292`, `#342`) for zero-diff drop-in migration before opt-in `--fix <id>` flags in `v1.1+`.

---

## Conclusions

`VarDictRust` (`v{{m:version}}`) consolidates the `VarDictJava -> R -> Perl` pipeline into a single static Rust binary. It holds exact Tier-1/2 oracle parity (`ε = {{m:eps_tol}}`, `HALF_EVEN`, record-aligned) across `{{m:s2_modes_total}}` operational modes, `{{m:s5_cross_platform_arms}}` sequencing chemistries (parity on all completing paired runs), and `{{m:s6_ga4gh_cells}}` completed Rust GA4GH `v3.5` truth evaluations. Bounding memory to `O(tile + halo)` delivers up to `{{m:s2_rchr1_speedup_vs_java}}`× faster chromosome-scale calling (`{{m:s2_rchr1_rss_ratio_vs_java}}`× lower RSS). It delivers `{{m:s4_exome_th32_speedup_vs_java}}`× faster `32`-thread whole-exome calling (`{{m:s4_exome_th32_rss_ratio_vs_java}}`× lower RSS). This enables unchunked whole-genome discovery across `GIAB v5.0q T2T-Q100` and `CMRG v1.00`.

---

## Methods

### Consolidated Rust Architecture and Coordinate Type Safety
`VarDictRust` is a two-crate Rust workspace: `vardictrust-core` (library engine, compiled with `#![forbid(unsafe_code)]` and zero `unwrap()`/`expect()`/`panic!()` calls outside unit tests) and `vardictrust` (`clap` CLI supporting all `VarDictJava` v1.8.3 flags). Distinct `Pos0` (`0`-based half-open BAM/CRAM coordinates from `noodles`) and `Pos1` (`1`-based inclusive VCF/VarDict coordinates) newtypes enforce explicit `.to_one_based()` and `.to_zero_based()` conversions. Sequences and qualities are processed as byte slices (`&[u8]`) with Phred+33 saturation arithmetic (`Phred(byte.saturating_sub(33))`).

### Bounded Tile-Plus-Halo Streaming and Deterministic Realignment
Intervals $\le {{m:s1_tile_1m}}$ bp execute as a single untiled window. Intervals $> {{m:s1_tile_1m}}$ bp (or with explicit `--tile-size`) are partitioned into non-overlapping ownership tiles $[T_k^{\text{start} }, T_k^{\text{end} })$ of default width `{{m:rust_default_tile_kbp}}` kbp. Each tile is padded on both sides by a `{{m:rust_halo_kbp}}` kbp realignment halo (`{{m:rust_splice_halo_kbp}}` kbp in `-i` splice mode; `Figure 1b`). A `Rayon` work-stealing pool streams alignments overlapping $[T_k^{\text{start} } - H, T_k^{\text{end} } + H)$ via `noodles` (multi-threaded BGZF BAM and reference-decoded CRAM `-G`). It runs CIGAR modification (`cigar/modify.rs`), pileup accumulation (`cigar/parse.rs`), local realignment (`realign.rs`), SV assembly (`sv.rs`), and variant finalization (`tovars.rs`). It emits only variants whose canonical `Pos1` start falls in $[T_k^{\text{start} }, T_k^{\text{end} })$. Per-tile pileup maps drop immediately upon tile completion (RAII), bounding heap memory to active tiles in flight. Where `VarDictJava` output order depends on `java.util.HashMap` bucket iteration (`#282`), `VarDictRust` emulates JDK8 `HashMap` spread, bin, resize, and treeify semantics via `JavaHashMap`, and uses `IndexMap` (`LinkedHashMap`) and `BTreeMap` (`TreeMap`) elsewhere.

### In-Process Statistical Post-Processing and Banker's Rounding (`HALF_EVEN`)
`VarDictRust` computes hypergeometric two-sided Fisher exact p-values and conditional maximum-likelihood odds ratios (`uniroot` zero-in with tolerance $\text{ulp}(1)^{0.25}$) natively in Rust (`stats/fisher.rs`), followed by VCF 4.3 formatting (`output/`), replacing R (`teststrandbias.R` / `testsomatic.R`) and Perl (`var2vcf_valid.pl` / `var2vcf_paired.pl`). All floating-point fields (`AF`, `AdjAF`, `HIAF`, `SBF`, `ODDRATIO`, `SSF`, `MSI`, `QUAL`) use IEEE 754 5-decimal round-half-to-even (`round_half_even`), matching R `round(x, 5)` and Java `DecimalFormat` (`HALF_EVEN`) in both no-`--fisher` (trailing zeros kept) and `--fisher` (trailing zeros stripped) modes within `ε = {{m:eps_tol}}` (`10^{-6}`).

### Cloud Benchmarking, Reproducibility Protocol, and Statistical Analysis
All `{{m:total_cells}}` benchmark cells ran on Google Cloud Platform (`{{m:host}}`; `{{m:nvme_raid0_tb}}` TB local NVMe RAID0 scratch; exact commands and SHA-256 manifests in `Supplementary Methods` and `Supplementary Table S1`). Timed cells ran `n = 3` independent technical replicates under `/usr/bin/time -v` to record wall-clock seconds and `Maximum resident set size (kbytes)`. Wall-clock timeouts surface as exit code `rc=124` (`3`-hour limit); JVM heap exhaustion surfaces as `rc=3` (`OutOfMemoryError`). Binary concordance proportions use exact two-sided Clopper-Pearson `95%` binomial CIs ($\alpha = 0.05$); see `Supplementary Note 3` for the Beta-inverse derivation and `v1.0.x` versus `v1.1+` `--fix` governance. Runtime and RSS across the `12` AML panel configurations were compared using two-sided Mann-Whitney $U$ tests and exact paired sign-flip permutation tests ($2^{-12}$). Truth concordance was evaluated with `rtg vcfeval` (`v3.12.1`) against official NIST GIAB and SEQC2 high-confidence BED masks (`HG002 v4.2.1`, `HG002 v5.0q T2T-Q100`, `CMRG v1.00`, `Mosaic v1.10`, `HG005 v4.2.1`, `SEQC2 v1.2.1`, and `HG008 v0.3`/`v0.4`). *In silico* tumor-purity titrations (`20%`, `10%`, `5%`, `2%`, `1%`) were built via seed-controlled `samtools view -s` subsampling and merging of HG008 tumor and normal BAMs (`Supplementary Methods`).

A single-command `VarDictRust` invocation replacing `VarDictJava | teststrandbias.R | var2vcf_valid.pl` is:
```bash
# Single-binary VarDictRust GERMLINE/SOMATIC calling with native Fisher exact test & VCF emission:
vardictrust -G GRCh38.fa -f {{m:f_std}} -N "TUMOR|NORMAL" \
  -b "tumor.bam|normal.bam" -th 32 --fisher -z -c 1 -S 2 -E 3 -g 4 \
  --vcf -A -E-filter 0.05 -P 0 target.bed > calls.vcf
```

---

## Declarations

### Ethics approval and consent to participate
Not applicable. All sequencing datasets analyzed in this study are publicly available, open-access reference materials from the NIST Genome in a Bottle (GIAB) Consortium (`HG002`, `HG005`, `HG008`), the FDA SEQC2 Consortium (`HCC1395`/`HCC1395BL`), and the NCBI Sequence Read Archive (`PRJNA744870` / `PRJNA523380`).

### Consent for publication
Not applicable.

### Availability of data and materials

- **Source Code and Binary Releases:** `VarDictRust` (`v{{m:version}}`, MIT License) is available at [https://github.com/palakpsheth/VarDictRust](https://github.com/palakpsheth/VarDictRust), with SLSA Level 3 static Linux binaries, Crates.io packages (`vardictrust` and `vardictrust-core`), Bioconda recipes, and `nf-core` modules.

- **Reproducible Benchmark Evidence and Manuscript Compiler:** All `{{m:total_cells}}` benchmark summaries, SHA-256 manifests, CLI scripts, `benchmarks/gcp-study-v106/master_gcp_metrics.json`, `benchmarks/gcp-study-v106/cells_60_inventory.tsv`, `paper/build_figures_and_tables.py`, and `paper/build_manuscript.py` are tracked in the repository. They are archived in `Supplementary Tables S1–S6`.

- **Public Accession Numbers:** Wolff [11] AML panel BAMs: NCBI SRA `SRR15006375`, `SRR15006376`, `SRR15006386`, `SRR15006540` (`PRJNA744870`), and `SRR8657348` (`PRJNA523380`). NIST GIAB reference datasets (`HG002` `v4.2.1`/`v5.0q T2T-Q100`/`CMRG v1.00`/`Mosaic v1.10`, `HG005`, and `HG008` `v0.3`/`v0.4` somatic callsets): `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/` (`HG008` somatic: `data_somatic/HG008/`; draft small-variant benchmarks `Liss_lab/analysis/NIST_HG008-T_somatic-smvar_DraftBenchmark_V0.3-20260425/` and `..._V0.4-20261006/`). SEQC2 `HCC1395`/`HCC1395BL` (`v1.2.1`): `PRJNA506071`.

### Competing interests
The author declares no competing financial or non-financial interests.

### Funding
Cloud compute on Google Cloud Platform was self-funded by the author; no external funding was received.

### Authors' contributions (CRediT)
**Palak Sheth:** Conceptualization, Methodology, Software, Validation, Formal Analysis, Investigation, Resources, Data Curation, Writing – Original Draft, Writing – Review & Editing, Visualization, Project Administration. AI pair-programming assistance (`Claude` / `Gemini`) was used under human verification and automated differential oracle-testing gates per the [rewrites.bio](https://rewrites.bio) transparency commitment.

### Acknowledgements
We thank Zhongwu Lai, Aleksandra Markovets, and the AstraZeneca-NGS team for creating and open-sourcing `VarDict` and `VarDictJava`, the NIST Genome in a Bottle Consortium and FDA SEQC2 Consortium for public truth benchmarks, and Lukas Wolff for open-sourcing `vardictcpp`.

---

## References

1. Lai Z, Markovets A, Ahdesmaki M, Chapman B, Hofmann O, McEwen R, et al. VarDict: a novel and versatile variant caller for next-generation sequencing in cancer research. *Nucleic Acids Res.* 2016;44(11):e108. doi:10.1093/nar/gkw227.
2. Koboldt DC. Best practices for variant calling in clinical sequencing. *Genome Med.* 2020;12(1):91. doi:10.1186/s13073-020-00791-w.
3. Garrison E, Marth G. Haplotype-based variant detection from short-read sequencing. *arXiv* [Preprint]. 2012;arXiv:1207.3907.
4. Kim S, Scheffler K, Halpern AL, Bekritsky MA, Noh E, Källberg M, et al. Strelka2: fast and accurate calling of germline and somatic variants. *Nat Methods.* 2018;15(8):591–4. doi:10.1038/s41592-018-0051-x.
5. Cooke DP, Wedge DC, Lunter G. A unified haplotype-based method for accurate and comprehensive variant calling. *Nat Biotechnol.* 2021;39(7):885–92. doi:10.1038/s41587-021-00861-3.
6. Poplin R, Chang PC, Alexander D, Schwartz S, Colthurst T, Ku A, et al. A universal SNP and small-indel variant caller using deep neural networks. *Nat Biotechnol.* 2018;36(10):983–7. doi:10.1038/nbt.4235.
7. Wilm A, Aw PPK, Bertrand D, Yeo GHT, Ong SH, Wong CH, et al. LoFreq: a sequence-quality aware, ultra-sensitive variant caller for uncovering cell-population heterogeneity from high-throughput sequencing datasets. *Nucleic Acids Res.* 2012;40(22):11189–201. doi:10.1093/nar/gks918.
8. Benjamin D, Sato T, Cibulskis K, Getz G, Stewart C, Lichtenstein L. Calling somatic SNVs and indels with Mutect2. *bioRxiv* [Preprint]. 2019;doi:10.1101/861054.
9. Garcia M, Juhos S, Larsson M, Olason PI, Martin M, Eisfeldt J, et al. Sarek: a portable workflow for whole-genome sequencing analysis of germline and somatic variants. *F1000Research.* 2020;9:63. doi:10.12688/f1000research.16665.2.
10. Fang LT, Zhu B, Zhao Y, Chen W, Yang Z, Kerrigan L, et al. Establishing community reference samples, data and call sets for benchmarking cancer mutation detection using whole-genome sequencing. *Nat Biotechnol.* 2021;39(9):1151–60. doi:10.1038/s41587-021-00993-6.
11. Wolff J. vardictcpp: a C++ reimplementation of VarDict variant caller. *GigaByte.* 2026. doi:10.46471/gigabyte.194.
12. Cleary JG, Braithwaite R, Gaastra K, Hilbush BS, Inglis S, Irvine SA, et al. Comparing variant call files for performance benchmarking of next-generation sequencing variant calling pipelines. *bioRxiv* [Preprint]. 2015;doi:10.1101/023754.
13. Wagner J, Olson ND, Harris L, Khan Z, Farek J, Mahmoud M, et al. Benchmarking challenging small variants with linked and long reads. *Cell Genomics.* 2022;2(5):100128. doi:10.1016/j.xgen.2022.100128.
14. Rhie A, Nurk S, Cechova M, Hoyt SJ, Taylor DJ, Altemose N, et al. The complete sequence of a human Y chromosome. *Nature.* 2023;621(7978):344–54. doi:10.1038/s41586-023-06457-y.
15. Wagner J, Olson ND, Harris L, Khan Z, Farek J, Mahmoud M, et al. Curated variation benchmarks for challenging medically relevant autosomal genes. *Nat Biotechnol.* 2022;40(5):672–80. doi:10.1038/s41587-021-01158-1.
16. Nurk S, Koren S, Rhie A, Rautiainen M, Bzikadze AV, Mikheenko A, et al. The complete sequence of a human genome. *Science.* 2022;376(6588):44–53. doi:10.1126/science.abj6987.
17. Krusche P, Trigg L, Boutros PC, Mason CE, De La Vega FM, Moore BL, et al. Best practices for benchmarking germline small-variant calls in human genomes. *Nat Biotechnol.* 2019;37(5):555–60. doi:10.1038/s41587-019-0054-x.
18. Zook JM, Catoe D, McDaniel J, Vang L, Spies N, Sidow A, et al. Extensive sequencing of seven human genomes to characterize benchmark reference materials. *Sci Data.* 2016;3:160025. doi:10.1038/sdata.2016.25.
19. McDaniel JH, Patel V, Olson ND, et al. Development and extensive sequencing of a broadly-consented Genome in a Bottle matched tumor-normal pair. *Sci Data.* 2025;12:1195. doi:10.1038/s41597-025-05438-2.
20. Daniels CA, Abdulkadir AA, Cleveland MH, McDaniel JH, et al. Characterization of subclonal variants in HG002 Genome in a Bottle reference material as a resource for benchmarking variant callers. *Cell Genom.* 2026;6(4):101104. doi:10.1016/j.xgen.2025.101104.
21. Behera S, Catreux S, Rossi M, Truong S, Huang Z, Ruehle M, et al. Comprehensive genome analysis and variant detection at scale using DRAGEN. *Nat Biotechnol.* 2024;43:1177–91. doi:10.1038/s41587-024-02382-1.
22. Park J, Cook DE, Chang PC, Kolesnikov A, Brambrink L, Mier JC, et al. Accurate somatic small variant discovery for multiple sequencing technologies with DeepSomatic. *Nat Biotechnol.* 2025;44(9):1569–78. doi:10.1038/s41587-025-02839-x.
23. Zheng Z, Chen L, Su J, Yu X, He M, Lee Y-L, et al. ClairS: a deep-learning method for long-read tumor–normal pair somatic small variant calling. *Nat Methods.* 2026;23(8):1507–20. doi:10.1038/s41592-026-03152-4.
24. Narzisi G, Corvelo A, Arora K, Bergmann EA, Shah M, Musunuri R, et al. Genome-wide somatic variant calling using localized colored de Bruijn graphs. *Commun Biol.* 2018;1:20. doi:10.1038/s42003-018-0023-9.
25. Hartwig Medical Foundation. SAGE somatic variant caller. HMFTools. https://github.com/hartwigmedical/hmftools/tree/master/sage.
26. Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, et al. Twelve years of SAMtools and BCFtools. *GigaScience.* 2021;10(2):giab008. doi:10.1093/gigascience/giab008.
