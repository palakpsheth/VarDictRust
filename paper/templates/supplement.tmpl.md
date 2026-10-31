# Supplementary Material: VarDictRust: A Memory-Bounded, Oracle-Equivalent Rust Consolidation of the VarDict Variant-Calling Pipeline for Deep Whole-Genome and Ultra-Low-VAF Sequencing

**Authors:** Palak Sheth<sup>1,*</sup>  
<sup>1</sup>Independent Researcher, [rewrites.bio](https://rewrites.bio) Initiative (`palak@rewrites.bio`)  
**Engine Release:** `VarDictRust` `v{{m:version}}` (commit `2a8e9a3`, engine fingerprint `{{m:fingerprint}}`, performance fingerprint `{{m:perf_fingerprint}}`)  
**Reference Oracles:** `{{m:oracle_version}}` and `{{m:cpp_version}}`  
**Cloud Benchmark Host:** `{{m:host}}` (`{{m:nvme_raid0_tb}}` TB NVMe RAID0 scratch; exploratory staging also verified on `c3d-standard-32`)

---

## Table of Contents

1. **Supplementary Methods: Exact Command-Line Invocations Across All 60 Benchmark Cells**
   - 1.1 Stage 1: Ceiling, Multi-Region Boundary, and Tile-Invariance Cells (`4` cells)
   - 1.2 Stage 2: Tier-P AML Replication, 13-Mode Sweep, Field Callers, and `-R chr1` (`17` cells)
   - 1.3 Stage 3: Full-Sample {{m:s3_rna_bam_gib}} GiB RNA-Seq Splice Calling (`1` cell)
   - 1.4 Stage 4: Whole-Exome Scaling, Chromosome 1 Span-to-RSS, Node Packing, and `QUAL` ROC (`15` cells)
   - 1.5 Option C: Comprehensive 87-Arm `.var` and `.vcf` Parity Sweep (`12` cells)
   - 1.6 Stage 5: Full-Genome WGS (Arm A vs Arm B), 10-Platform Chr20 Panel, and Supply Chain (`7` cells)
   - 1.7 Stage 6: GA4GH v3.5 Stratification, GIAB v5.0Q/CMRG/Mosaic, HG008 v0.4 11-Caller Matrix, and Low-VAF Titration (`4` cells)
2. **Supplementary Note 1: Analytical Derivation of JVM (`VarDictJava`) and `vardictcpp` Heap Retention (`Arm A` vs `Arm B`)**
3. **Supplementary Note 2: Ten-Point Technical Audit of `vardictcpp` Divergences (`CPP-1` through `CPP-10`)**
4. **Supplementary Note 3: Biostatistical Methodology, Exact Clopper-Pearson Intervals, and `--fix` Governance**
5. **Supplementary Tables (`Table S1` – `Table S6`) and Complete 60-Cell Verification Ledger**

---

## 1. Supplementary Methods: Exact Command-Line Invocations Across All 60 Benchmark Cells

All `{{m:total_cells}}` verification cells (`{{m:compute_cells}}` cloud compute cells plus `1` certified local-to-cloud evidence bundle `evidence_certified`) were executed on `{{m:host}}` (`Ubuntu 24.04 LTS`, Linux kernel `6.8`, `{{m:nvme_raid0_tb}}` TB 4-way NVMe RAID0 mounted at `/var/tmp/vardict-paper`). Every timed execution was wrapped in `/usr/bin/time -v` across `n = 3` independent technical replicates (`Stage 2`, `Stage 3`, `Stage 4`) or guarded by `timeout 10800` (`3` hours, exit code `rc=124` on timeout; `rc=3` on JVM `OutOfMemoryError`).

### 1.1 Stage 1: Ceiling, Multi-Region Boundary, and Tile-Invariance Cells (`s1_*`)
Stage 1 comprises `4` cells (`s1_ceiling_amplicon`, `s1_ceiling_hg008`, `s1_gate7_ref_multiregion`, and `s1_h1a_simple_tiles`):
```bash
# Cell 1: s1_ceiling_amplicon (EGFR/MRE11a amplicon -a {{m:amp_param}}, {{m:s1_amplicon_vars}} variants)
vardictrust -G hg19.fa -f {{m:f_std}} -N SAMPLE -b amplicon.bam \
  -a {{m:amp_param}} --fisher -c 1 -S 2 -E 3 -g 4 --vcf amplicon.bed > rust_amp.vcf

# Cell 2: s1_ceiling_hg008 (NIST GIAB HG008 somatic slice, {{m:s1_hg008_vars}} variants)
vardictrust -G GRCh38_GIABv3.fasta -f {{m:f_std}} -N "HG008-T|HG008-N" \
  -b "HG008-T.bam|HG008-N.bam" -th 8 --fisher -z -c 1 -S 2 -E 3 -g 4 \
  --vcf hg008_slice.bed > rust_hg008.vcf

# Cell 3: s1_gate7_ref_multiregion (HG002 chr20 multi-region boundary gate, {{m:s1_gate7_vars}} variants)
bash harness/compare.sh GRCh38.fa hg002_chr20.bam chr20_multiregion.bed HG002

# Cell 4: s1_h1a_simple_tiles (Tile-size invariance: {{m:s1_tile_sizes_bp_list}} bp -> {{m:s1_tile_vars}} variants)
for TILE in 500000 1000000 2000000 10000000; do
  VARDICTRUST_TILE_BP="${TILE}" vardictrust -G GRCh38.fa -f {{m:f_std}} -N HCC1395 \
    -b HCC1395.bam -th 8 --fisher -c 1 -S 2 -E 3 -g 4 --vcf seqc2_slice.bed > "tile_${TILE}.vcf"
done
```

### 1.2 Stage 2: Tier-P AML Replication, 13-Mode Sweep, Field Callers, and `-R chr1` (`s2_*`)
Stage 2 comprises `17` cells replicating and extending Wolff et al. (2026) [11]:
- **12 Tier-P v2 AML Panel Cells (`n = 3` technical replicates each at `-th 1` and `-th 8`):** `s2_tierp_v2_SRR15006375_th1`, `s2_tierp_v2_SRR15006375_th8`, `s2_tierp_v2_SRR15006376_th1`, `s2_tierp_v2_SRR15006376_th8`, `s2_tierp_v2_SRR15006386_th1`, `s2_tierp_v2_SRR15006386_th8`, `s2_tierp_v2_SRR15006540_th1`, `s2_tierp_v2_SRR15006540_th8`, `s2_tierp_v2_SRR8657348_th1`, `s2_tierp_v2_SRR8657348_th8`, `s2_tierp_v2_somatic_76_86_th1`, and `s2_tierp_v2_somatic_76_86_th8`.
- **Contiguous Chromosome 1 (`-R chr1:1-249250621`, `{{m:s2_rchr1_bp}}` bp):** `s2_table3_R_chr1_java_200g_n3` (`JAVA_OPTS="-Xmx4g"` OOM `rc=3` at `{{m:s2_rchr1_java4g_oom_wall_s}}` s; `JAVA_OPTS="-Xmx200g"` completes in `{{m:s2_rchr1_java200g_wall_s}}` s at `{{m:s2_rchr1_java200g_rss_gib}}` GiB RSS) and `s2_table3_R_chr1_rust_and_cpp` (`VarDictRust` `{{m:s2_rchr1_rust_wall_s}}` s, `{{m:s2_rchr1_rust_rss_mb}}` MB RSS vs `vardictcpp` `{{m:s2_rchr1_cpp_wall_s}}` s, `{{m:s2_rchr1_cpp_rss_gib}}` GiB RSS).
- **13-Mode Functional Sweep & Field Callers:** `s2_tierp_mode_parity_sweep` (`{{m:s2_modes_rust_match}}/{{m:s2_modes_total}}` modes), `s2_fig1_fig3_field_callers` (`FreeBayes`, `LoFreq`, `Mutect2`, `Octopus`, `Strelka2`), and `s2_fig2_hg002_chr20_truth` (`{{m:s2_fig2_parity_vars}}` candidate variants).

```bash
# Tier-P single-sample & paired somatic replication (s2_tierp_v2_*):
/usr/bin/time -v vardictrust -G hg19.fa -f {{m:f_std}} -N SRR8657348 \
  -b SRR8657348.bam -th 8 --fisher -z -c 1 -S 2 -E 3 -g 4 --vcf SRR8657348.bed > SRR8657348.vcf

# Table 3 contiguous chromosome 1 (-R chr1:1-249250621, s2_table3_R_chr1_*):
/usr/bin/time -v vardictrust -G hg19.fa -f {{m:f_std}} -N SRR15006375 \
  -b SRR15006375.bam -th 8 -R chr1:1-249250621 > rust_chr1.var
JAVA_OPTS="-Xmx200g" /usr/bin/time -v VarDict -G hg19.fa -f {{m:f_std}} -N SRR15006375 \
  -b SRR15006375.bam -th 8 -R chr1:1-249250621 > java_chr1.var
```

### 1.3 Stage 3: Full-Sample {{m:s3_rna_bam_gib}} GiB RNA-Seq Splice Calling (`s3_rna_splice_parity_n3`)
Cell `s3_rna_splice_parity_n3` (`n = 3` technical replicates) evaluates `-i` RNA splice-junction mode on the `{{m:s3_rna_bam_gib}}` GiB HG002 GM24385 RNA-seq alignment (`{{m:s3_rna_splice_rows}}` `.var` and `.splice.tsv` records, `{{m:s3_rna_rust_wall_s}}` s vs `{{m:s3_rna_java_wall_s}}` s, `{{m:s3_rna_speedup_vs_java}}`× speedup, `{{m:s3_rna_rss_ratio_vs_java}}`× RSS reduction):
```bash
/usr/bin/time -v vardictrust -G GRCh38.fa -f {{m:f_std}} -N HG002_RNA \
  -b HG002_GM24385_RNA.bam -i -th 8 -c 1 -S 2 -E 3 -g 4 gencode_exons.bed > rust_rna.var
```

### 1.4 Stage 4: Whole-Exome Scaling, Chromosome 1 Span-to-RSS, Node Packing, and `QUAL` ROC (`s4_*`)
Stage 4 comprises `15` cells (`s4_perf_budget_4targets_n3`, `s4_scale_exome_3way_th1`, `s4_scale_exome_3way_th2`, `s4_scale_exome_3way_th4`, `s4_scale_exome_3way_th8`, `s4_scale_exome_3way_th16`, `s4_scale_exome_3way_th32`, `s4_span_chr1_3way_th8`, `s4_span_chr1_3way_th16`, `s4_span_chr1_3way_th32`, `s4_span_chr1_f0`, `s4_jobs_per_node_2x4c_n3`, `s4_multimodal_and_cram_twin`, `s4_tier3_truth_and_qual_sweep`, and `evidence_certified`):
```bash
# Whole-exome scaling across 1, 2, 4, 8, 16, 32 threads (s4_scale_exome_3way_th1..th32):
for TH in 1 2 4 8 16 32; do
  /usr/bin/time -v vardictrust -G GRCh38.fa -f {{m:f_std}} -N "HCC1395|HCC1395BL" \
    -b "HCC1395.bam|HCC1395BL.bam" -th "${TH}" --fisher -z -c 1 -S 2 -E 3 -g 4 \
    --vcf seqc2_exome.bed > "seqc2_th${TH}.vcf"
done

# Chromosome 1 span-to-RSS invariance across 1M, 10M, 50M, 100M, 249M bp (s4_span_chr1_3way_th8..th32):
for SPAN in 1000000 10000000 50000000 100000000 248956422; do
  /usr/bin/time -v vardictrust -G GRCh38.fa -f {{m:f_std}} -N HG002 \
    -b HG002_chr1.bam -th 32 -R "chr1:1-${SPAN}" > "chr1_${SPAN}.var"
done

# All-positions (-f {{m:f_zero}}) pileup stress on chr1 10 Mbp (s4_span_chr1_f0) and 2x4c node packing (s4_jobs_per_node_2x4c_n3):
/usr/bin/time -v vardictrust -G GRCh38.fa -f {{m:f_zero}} -N HG002 -b HG002_chr1.bam -th 16 -R chr1:10000001-20000000 > chr1_f0.var
```

### 1.5 Option C: Comprehensive 87-Arm `.var` and `.vcf` Parity Sweep (`sweep_*`)
The `12` Option C cells (`sweep_4targets_full`, `sweep_R_chr1_full`, `sweep_hg002_chr20_exome`, `sweep_multimodal_and_cram_twin`, `sweep_rna_hg002_52gb`, `sweep_span_chr1_f0`, `sweep_tierp_13modes`, `sweep_tierp_SRR15006540`, `sweep_tierp_SRR8657348`, `sweep_tierp_somatic_76_86`, `sweep_ultradeep_mrd_stress`, and `sweep_wgs_mrd_matrix`) validate `.var` and `.vcf` equivalence across `{{m:sweep_total_arms}}` comparative arms (`{{m:sweep_match_arms}}/{{m:sweep_total_arms}}` `MATCH`, `{{m:sweep_match_pct}}`, `Table S2`):
```bash
# Ultra-deep 45000x ctDNA MRD stress (-f {{m:f_mrd_ud}} and -f {{m:f_zero}} with --fisher, sweep_ultradeep_mrd_stress):
python3 harness/compare_var.py golden/ultradeep_f0005.var rust_ultradeep_f0005.var --eps {{m:eps_tol}}
# WGS & exome MRD matrix (-f {{m:f_mrd_wgs}} and -f {{m:f_zero}}, sweep_wgs_mrd_matrix):
python3 harness/compare_vcf.py golden/hg005_300x_mrd.vcf rust_hg005_300x_mrd.vcf --eps {{m:eps_tol}}
```

### 1.6 Stage 5: Full-Genome WGS (Arm A vs Arm B), 10-Platform Chr20 Panel, and Supply Chain (`s5_*`)
Stage 5 comprises `7` cells (`s5_wgs_hg002_full_contigs_3way_th32`, `s5_wgs_hg002_full_giab_3way_th32`, `s5_wgs_hg008_full_contigs_3way_th32`, `s5_wgs_hg008_full_giab_3way_th32`, `s5_grch38_field_callers_truth`, `s5_cross_platform_short_read_chr20_panel`, and `s5_h4_maintenance_security_matrix`):
```bash
# Arm A: Un-chunked full chromosomes (chr1..chr22, s5_wgs_hg002_full_contigs_3way_th32 & s5_wgs_hg008_full_contigs_3way_th32):
timeout 10800 /usr/bin/time -v vardictrust -G GRCh38_GIABv3.fasta -f {{m:f_std}} -N HG002 \
  -b HG002_2x250_131GB.bam -th 32 --fisher -z -c 1 -S 2 -E 3 -g 4 --vcf grch38_22autosomes_full.bed > hg002_arm_a.vcf

timeout 10800 /usr/bin/time -v vardictrust -G GRCh38_GIABv3.fasta -f {{m:f_std}} -N "HG008-T|HG008-N" \
  -b "HG008-T_195x.bam|HG008-N_169x.bam" -th 32 --fisher -z -c 1 -S 2 -E 3 -g 4 --vcf grch38_22autosomes_full.bed > hg008_arm_a.vcf

# Arm B: Pre-chunked <=50 kbp GIAB callable intervals (s5_wgs_hg002_full_giab_3way_th32 & s5_wgs_hg008_full_giab_3way_th32):
timeout 10800 /usr/bin/time -v vardictrust -G GRCh38_GIABv3.fasta -f {{m:f_std}} -N HG002 \
  -b HG002_2x250_131GB.bam -th 32 --fisher -z -c 1 -S 2 -E 3 -g 4 --vcf HG002_GRCh38_v421_50kb.bed > hg002_arm_b.vcf
```

### 1.7 Stage 6: GA4GH v3.5 Stratification, GIAB v5.0Q/CMRG/Mosaic, HG008 v0.4 11-Caller Matrix, and Low-VAF Titration (`s6_*`)
Stage 6 comprises `4` cells (`s6_ga4gh_v35_stratification_matrix`, `s6_hg002_v50q_cmrg_mosaic_eval`, `s6_hg008_modern_ml_and_commercial_callers`, and `s6_low_vaf_titration_truth_curve`):
```bash
# Controlled in silico tumor-purity titration (20%, 10%, 5%, 2%, 1% purity at -f {{m:f_titration}}, s6_low_vaf_titration_truth_curve):
for PURITY in 20 10 5 2 1; do
  vardictrust -G GRCh38_GIABv3.fasta -f {{m:f_titration}} -N "HG008_P${PURITY}|HG008-N" \
    -b "hg008_titration_p${PURITY}.bam|HG008-N.bam" -th 16 --fisher -z -c 1 -S 2 -E 3 -g 4 \
    --vcf hg008_slice.bed > "rust_titration_p${PURITY}.vcf"
done

# Next-gen GIAB v5.0Q T2T-Q100, CMRG v1.00, and Mosaic v1.10 evaluation (s6_hg002_v50q_cmrg_mosaic_eval):
rtg vcfeval -b HG002_GRCh38_v5.0q_T2T.vcf.gz -e HG002_GRCh38_v5.0q_T2T.bed \
  -c hg002_arm_a.vcf.gz -t GRCh38.sdf -o eval_v50q_arm_a
```

---

## 2. Supplementary Note 1: Analytical Derivation of JVM (`VarDictJava`) and `vardictcpp` Heap Retention (`Arm A` vs `Arm B`)

In `VarDictJava` v1.8.3 (`CigarParser.java` and `VariationRealigner.java`), the pileup engine processes each BED region $[S, E]$ of length $L = E - S + 1$ by instantiating two region-wide heap structures before emitting a single variant row:
1. A reference base map `Reference.referenceSequences` (`Map<Integer, Character>`) containing $O(L)$ boxed `java.lang.Integer` (`16` bytes) and `java.lang.Character` (`16` bytes) objects plus `HashMap$Node` headers (`32` bytes), consuming $\approx 64\text{–}96$ bytes per reference base ($\approx 16\text{–}24$ GiB for human `chr1` of length $L = {{m:s2_rchr1_bp}}$ bp per worker thread).
2. A position-indexed variation map `Map<Integer, VariationMap<String, Variation>> nonInsertionVariants` and `insertionVariants`, where every covered genomic coordinate $p \in [S, E]$ allocates a `VariationMap` and one `Variation` object per observed allele (including the reference allele). At sequencing depth $D$, the active heap footprint per worker thread scales as:
   $$\text{RSS}_{\text{Java / C++} }(L, T) = \Theta\!\left(\sum_{t=1}^{T} L_t \cdot (c_{\text{ref} } + c_{\text{var} } \cdot \bar{a}(D))\right)$$
   where $L_t$ is the length of the region assigned to thread $t$ and $\bar{a}(D)$ is the mean number of candidate alleles per locus. When $L_t \le {{m:giab_arm_b_max_kbp}}$ kbp (Arm B), $L_t$ is small enough that GC can reclaim objects between intervals (`{{m:s5_hg002_giab_java_rss_gb}}` GB on HG002 Arm B; `{{m:s5_hg008_giab_java_rss_gb}}` GB with timeout `rc=124` on HG008 Arm B). When $L_t \in [50, 249]$ Mbp (Arm A un-chunked chromosomes or `-R chr1:1-249250621`), live objects cannot be garbage-collected until the entire chromosome completes, driving `VarDictJava` to `{{m:s2_rchr1_java200g_rss_gib}}` GiB on `chr1` alone and `>{{m:s5_hg002_contigs_java_rss_gb}}` GB (`rc=3` `OutOfMemoryError`) on multi-threaded full genomes (`Table 2`, `Table 5`). Similarly, `vardictcpp` retains per-region C++ `robin_hood` / `std::unordered_map` structures across each interval, consuming `{{m:s2_rchr1_cpp_rss_gib}}` GiB on `chr1`, `{{m:s5_hg002_contigs_cpp_rss_gb}}` GB on HG002 Arm A, and `{{m:s5_hg008_contigs_cpp_rss_gb}}` GB (`rc=124` timeout) on HG008 Arm A unless an external shell script pre-splits the BED file (`CPP_CHUNK={{m:cpp_chunk_bp}}`).

By contrast, `VarDictRust` partitions any region of arbitrary length $L$ into fixed tiles of width $W = {{m:s1_tile_1m}}$ bp with halo $H = {{m:rust_halo_kbp}}$ kbp (`Figure 1b`), bounding memory strictly to:
$$\text{RSS}_{\text{VarDictRust} }(L, T) = O\!\left(T \cdot (W + 2H)\right) \perp L$$
Empirically, increasing $L$ by `250`× on chromosome 1 (`1` Mbp to `249` Mbp, `s4_span_chr1_3way_th8`..`th32`) changes `VarDictRust` RSS by only `{{m:s4_span_th8_ratio_249m_10m}}`×–`{{m:s4_span_th32_ratio_249m_10m}}`× (`{{m:s4_span_th8_249m_rust_gib}}`–`{{m:s4_span_th32_249m_rust_gib}}` GiB), enabling un-chunked full-genome calling in `{{m:s5_hg002_contigs_rust_rss_gb}}` GB (HG002) and `{{m:s5_hg008_contigs_rust_rss_gb}}` GB (HG008).

---

## 3. Supplementary Note 2: Ten-Point Technical Audit of `vardictcpp` Divergences (`CPP-1` through `CPP-10`)

Across cells `s2_tierp_mode_parity_sweep`, `s5_h4_maintenance_security_matrix`, and `s5_wgs_hg002_full_contigs_3way_th32`, we cataloged ten reproducibility and capability gaps in `{{m:cpp_version}}` (`Table 8`, `Table S6`):
- **`CPP-1` (No Native VCF Output):** `vardictcpp` only outputs raw intermediate `.var` text and still requires piping into external Perl scripts (`var2vcf_valid.pl` / `var2vcf_paired.pl`) to generate VCF files.
- **`CPP-2` (Missing Modes: Amplicon `-a`, Splice `-i`, and Structural Variants):** `vardictcpp` hardcodes `--nosv` and does not implement structural variant calling (`StructuralVariantsProcessor`), amplicon primer clipping (`-a`), or RNA-seq splice junction reporting (`-i`).
- **`CPP-3` (CRAM Incompatibility):** Passing a CRAM alignment with `-G ref.fa` fails or emits empty output in `vardictcpp`, whereas `VarDictRust` natively decodes CRAM via `noodles`.
- **`CPP-4` (`-R` Region Parser Crash):** Passing a chromosome name without colon coordinates (e.g., `-R chr20`) triggers an unhandled `std::out_of_range` / `SIGABRT` crash in `vardictcpp`.
- **`CPP-5` (Unbounded Interval Memory Without External `CPP_CHUNK`):** Without an external BED-splitting wrapper (`CPP_CHUNK={{m:cpp_chunk_bp}}`), `vardictcpp` consumes `{{m:s2_rchr1_cpp_rss_gib}}` GiB on `chr1`, `{{m:s5_hg002_contigs_cpp_rss_gb}}` GB on HG002 WGS, and `{{m:s5_hg008_contigs_cpp_rss_gb}}` GB (`rc=124` timeout) on HG008 WGS.
- **`CPP-6` (Chunk-Boundary Duplicate Emission):** When `vardictcpp` runs on chunked or contiguous whole genomes, it emits `{{m:s5_hg002_contigs_cpp_excess}}` excess duplicate/boundary-shifted calls on HG002 Arm A (`{{m:s5_hg002_contigs_cpp_vars}}` vs `{{m:s5_hg002_contigs_rust_vars}}`) and `+{{m:s5_hg002_giab_cpp_diff_vars}}` on Arm B (`{{m:s5_hg002_giab_cpp_vars}}` vs `{{m:s5_hg002_giab_vars}}`).
- **`CPP-7` (Fisher Exact Odds-Ratio & Rounding Drift):** `vardictcpp` diverges from R's `fisher.test` conditional MLE root-finder and `HALF_EVEN` 5-decimal rounding on edge contingency tables.
- **`CPP-8` (Single-End Realignment Drop):** In `s2_tierp_mode_parity_sweep`, `vardictcpp` diverges on `1/13` modes (`{{m:s2_modes_cpp_exact}}/{{m:s2_modes_total}}` exact) due to indel realignment differences on single-end/unpaired reads.
- **`CPP-9` (High-Purity Titration Indel Recall Drop):** At `20%` tumor purity (`s6_low_vaf_titration_truth_curve`), `vardictcpp` drops somatic indel $F_1$ to `{{m:s6_tit_p20_vardictcpp_indel_f1}}` (recall `{{m:s6_tit_p20_vardictcpp_indel_r_pct}}`) and SNV $F_1$ to `{{m:s6_tit_p20_vardictcpp_snv_f1}}` (recall `{{m:s6_tit_p20_vardictcpp_snv_r_pct}}`), whereas `VarDictRust` matches `VarDictJava` exactly (`Indel` $F_1 = {{m:s6_tit_p20_VarDictRust_indel_f1}}$, recall `{{m:s6_tit_p20_VarDictRust_indel_r_pct}}`; `SNV` $F_1 = {{m:s6_tit_p20_VarDictRust_snv_f1}}$, recall `{{m:s6_tit_p20_VarDictRust_snv_r_pct}}`).
- **`CPP-10` (Dynamic Shared-Library Fragility):** `vardictcpp` requires dynamic linkage against system `libhts`, `libcurl`, `liblzma`, `libbz2`, and `libdeflate`, plus a Perl runtime for VCF generation, whereas `VarDictRust` ships as a `{{m:s5_h4_musl_mb}}` MB zero-dependency static `musl` binary.

---

## 4. Supplementary Note 3: Biostatistical Methodology, Exact Clopper-Pearson Intervals, and `--fix` Governance

### 4.1 Exact Clopper-Pearson Binomial Confidence Intervals
For all binary oracle-concordance and completion rates ($k$ successes out of $n$ independent evaluation cells/arms), we compute exact two-sided $1 - \alpha = 95\%$ Clopper-Pearson binomial confidence intervals $[p_{\text{lo} }, p_{\text{hi} }]$ via the inverse Beta distribution:
$$p_{\text{lo} } = \text{Beta}^{-1}\!\left(\frac{\alpha}{2};\, k,\, n - k + 1\right), \quad p_{\text{hi} } = \text{Beta}^{-1}\!\left(1 - \frac{\alpha}{2};\, k + 1,\, n - k\right)$$
with $p_{\text{lo} } = 0$ when $k = 0$ and $p_{\text{hi} } = 1$ when $k = n$:
- Complete `{{m:total_cells}}`-cell verification campaign (`{{m:total_cells}}/{{m:total_cells}}` exit `0`): `100%` (`95%` CI `{{m:cp95_60_60}}`)
- Stage 2 `13`-mode functional parity sweep (`{{m:s2_modes_rust_match}}/{{m:s2_modes_total}}`): `100%` (`95%` CI `{{m:cp95_13_13}}`)
- Extended coverage Tier-2 modalities (`{{m:s4_ext_cov_modalities}}`): `100%` (`95%` CI `{{m:cp95_17_17}}`)
- Multimodal Tier-2 slices (`{{m:s4_multimodal_slices}}`): `100%` (`95%` CI `{{m:cp95_18_18}}`)
- Cross-platform chr20 runs where `VarDictJava` completed (`{{m:s5_chr20_java_completed_runs}}/{{m:s5_chr20_java_completed_runs}}`): `100%` (`95%` CI `{{m:cp95_chr20_java_completed}}`)
- Cross-platform `-f {{m:f_zero}}` somatic noise-floor `STRICT-BYTES` sweep (`{{m:s5_f0_noise_floor_strict_match}}`): `100%` (`95%` CI `{{m:cp95_noise_18_18}}`)
- GA4GH `v3.5` `{{m:s6_ga4gh_cells}}`-cell stratification matrix (`{{m:s6_ga4gh_cells}}/{{m:s6_ga4gh_cells}}`): `100%` (`95%` CI `{{m:cp95_70_70}}`)
- Option C parity sweep across all arms (`{{m:sweep_match_arms}}/{{m:sweep_total_arms}}`): `{{m:sweep_match_pct}}` (`95%` CI `{{m:cp95_78_87}}`), and `{{m:sweep_match_arms}}/{{m:sweep_match_arms}}` (`100%`, `95%` CI `{{m:cp95_78_78_single_region}}`) on single-region runs
- NIST Mosaic SNV `v1.10` sensitivity (`{{m:s6_mosaic_tp_best}}`): `{{m:s6_mosaic_recall_pct}}` (`95%` CI `{{m:cp95_84_85}}`)

### 4.2 Non-Parametric Runtime and Memory Hypothesis Testing
Across the $N = 12$ Tier-P sample-by-thread cells (`s2_tierp_v2_*`, each measured across `n = 3` technical replicates), we report both the cell mean and full replicate range (`min–max`, `Table S3`) rather than normal-theory standard errors (which are unstable at $n = 3$). Comparing `VarDictRust` against `VarDictJava` across the $N = 12$ cells yields two-sided Mann-Whitney $U = {{m:s2_tierp_mwu_wall_u}}$ ($p = {{m:s2_tierp_mwu_wall_p}}$) for wall-clock time and $U = {{m:s2_tierp_mwu_rss_u}}$ ($p = {{m:s2_tierp_mwu_rss_p}}$) for peak RSS, and an exact paired sign-flip permutation test p-value of $p = 2^{-12} = {{m:s2_tierp_sign_perm_p}}$.

### 4.3 Release Governance: `v1.0.x` Parity Versus `v1.1+` Opt-In `--fix` Flags
Under `ROADMAP.md`, all `v1.0.x` releases (`v{{m:version}}`) enforce strict oracle parity on valid inputs. True output-changing upstream bug fixes (e.g., `#92` `VD > DP`, `#62` over-filtering, `#292` `SHIFT3`, `#342` `<dup-N>` symbolic alleles) are gated behind explicit `--fix <id>` flags starting in `v1.1.0`, while parity-restoring fixes (`splice-tile-ownership`, `mate-window-ownership`) ship as verified defaults.

---

## 5. Supplementary Tables (`Table S1` – `Table S6`) and Complete 60-Cell Verification Ledger

### Table S1: Complete 60-Cell Cloud Benchmark Inventory and Provenance Summary

{{table:supp_table1}}

#### Machine-Verified 60-Cell Ledger (`benchmarks/gcp-study-v106/cells_60_inventory.tsv`)

{{cells_60_ledger}}

---

### Table S2: Option C 87-Arm `.var` and `.vcf` Differential Parity Matrix

{{table:supp_table2}}

---

### Table S3: Stage 2 Tier-P v2 Technical Replicate (`n = 3`) Distributions (`Mean` and `Min–Max`)

{{table:supp_table3}}

---

### Table S4: Stage 4 `QUAL` Threshold ROC Operating Points and Stage 5 `-f 0` Somatic Noise-Floor Sweep

{{table:supp_table4}}

---

### Table S5: Stage 6 70-Cell GA4GH `v3.5` Genomic Context Stratification Matrix (`10` Platforms × `7` Strata)

{{table:supp_table5}}

---

### Table S6: Stage 2 13-Mode Functional Sweep, Structural Variant/Splice/Amplicon Parity, and `CPP-1..10` Matrix

{{table:supp_table6}}
