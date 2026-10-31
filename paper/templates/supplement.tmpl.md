# Supplementary Material: VarDictRust: A Memory-Bounded, Oracle-Equivalent Rust Consolidation of the VarDict Variant-Calling Pipeline for Deep Whole-Genome and Ultra-Low-VAF Sequencing

**Authors:** [Palak Sheth](https://www.linkedin.com/in/palak-p-sheth/)<sup>1,*</sup>  
<sup>1</sup>Independent Researcher ([palakpsheth@gmail.com](mailto:palakpsheth@gmail.com) | [https://www.linkedin.com/in/palak-p-sheth/](https://www.linkedin.com/in/palak-p-sheth/))  
**Engine Release:** `VarDictRust` `v{{m:version}}` (commit `{{m:commit}}`, engine fingerprint `{{m:fingerprint}}`, performance fingerprint `{{m:perf_fingerprint}}`)  
**Reference Oracles:** `{{m:oracle_version}}` and `{{m:cpp_version}}`  
**Cloud Benchmark Host:** `{{m:host}}` (`{{m:nvme_raid0_tb}}` TB NVMe RAID0 scratch; exploratory staging also verified on `c3d-standard-32`)

---

## Table of Contents

1. **Supplementary Methods: Exact Command-Line Invocations Across All 60 Benchmark Cells**
    - 1.0 Manuscript Figure and Table to 60-Cell Provenance Crosswalk
    - 1.1 Domain 1: Ceiling, Multi-Region Boundary, and Tile-Invariance Cells (`4` cells)
    - 1.2 Domain 2: Wolff et al. AML Panel Replication, 13-Mode Sweep, Field Callers, and `-R chr1` (`17` cells)
    - 1.3 Domain 3: Full-Sample {{m:s3_rna_bam_gib}} GiB RNA-Seq Splice Calling (`1` cell)
    - 1.4 Domain 4: Whole-Exome Scaling, Chromosome 1 Span-to-RSS, Node Packing, and `QUAL` ROC (`15` cells)
    - 1.5 Multi-Regime Parameter Sweep: Comprehensive 87-Arm `.var` and `.vcf` Parity Sweep (`12` cells)
    - 1.6 Domain 5: Full-Genome WGS (Unchunked Full Chromosomes vs. Pre-Chunked Callable Intervals), 10-Platform Chr20 Panel, and Supply Chain (`7` cells)
    - 1.7 Domain 6: GA4GH v3.5 Stratification, GIAB v5.0q/CMRG/Mosaic, HG008 v0.4 11-Caller Matrix, and Low-VAF Titration (`4` cells)
2. **Supplementary Note 1: Analytical Derivation of JVM (`VarDictJava`) and `vardictcpp` Heap Retention (Unchunked Chromosomes vs. Pre-Chunked Callable Intervals)**
3. **Supplementary Note 2: Ten-Point Technical Audit of `vardictcpp` Divergences (`CPP-1` through `CPP-10`)**
4. **Supplementary Note 3: Biostatistical Methodology, Exact Clopper-Pearson Intervals, and `--fix` Governance**
5. **Supplementary Tables (`Supplementary Table S1` – `Supplementary Table S6`) and Complete 60-Cell Verification Ledger**
6. **References**

---

## 1. Supplementary Methods: Exact Command-Line Invocations Across All 60 Benchmark Cells

All `{{m:total_cells}}` verification cells (`{{m:compute_cells}}` cloud compute cells plus `1` certified local-to-cloud evidence bundle `evidence_certified`) ran on `{{m:host}}` (`Ubuntu 24.04 LTS`, Linux kernel `6.8`, `{{m:nvme_raid0_tb}}` TB 4-way NVMe RAID0 mounted at `/var/tmp/vardict-paper`). Timed runs used `/usr/bin/time -v` across `n = 3` independent technical replicates or `timeout 10800` (`3` hours, exit code `rc=124` on timeout; `rc=3` on JVM `OutOfMemoryError`).

### 1.0 Manuscript Figure and Table to 60-Cell Provenance Crosswalk
The crosswalk below maps each main-text Figure (`Figure 1`–`Figure 8`) and Table (`Table 1`–`Table 8`) to its benchmark directory in `benchmarks/gcp-study-v106/` (`Supplementary Table S1`):
- **`Figure 1` & `Table 1` (Architecture, Operational Modes, and 87-Arm Equivalence Sweep):** `s1_ceiling_amplicon`, `s1_ceiling_hg008`, `s1_gate7_ref_multiregion`, `s1_h1a_simple_tiles`, `s2_tierp_mode_parity_sweep`, `s3_rna_splice_parity_n3`, `s4_multimodal_and_cram_twin`, `evidence_certified`, and the `12` `sweep_*` cells (`sweep_4targets_full`, `sweep_R_chr1_full`, `sweep_hg002_chr20_exome`, `sweep_multimodal_and_cram_twin`, `sweep_rna_hg002_52gb`, `sweep_span_chr1_f0`, `sweep_tierp_13modes`, `sweep_tierp_SRR15006540`, `sweep_tierp_SRR8657348`, `sweep_tierp_somatic_76_86`, `sweep_ultradeep_mrd_stress`, `sweep_wgs_mrd_matrix`).
- **`Figure 2` & `Table 2` (Wolff et al. AML Cohort, Field Callers, and Contiguous `-R chr1`):** `s2_tierp_v2_SRR15006375_th1`, `s2_tierp_v2_SRR15006375_th8`, `s2_tierp_v2_SRR15006376_th1`, `s2_tierp_v2_SRR15006376_th8`, `s2_tierp_v2_SRR15006386_th1`, `s2_tierp_v2_SRR15006386_th8`, `s2_tierp_v2_SRR15006540_th1`, `s2_tierp_v2_SRR15006540_th8`, `s2_tierp_v2_SRR8657348_th1`, `s2_tierp_v2_SRR8657348_th8`, `s2_tierp_v2_somatic_76_86_th1`, `s2_tierp_v2_somatic_76_86_th8`, `s2_table3_R_chr1_java_200g_n3`, `s2_table3_R_chr1_rust_and_cpp`, and `s2_fig1_fig3_field_callers`.
- **`Figure 3` & `Table 3` (Germline/Somatic Truth Concordance, `QUAL` Operating Points, and Low-VAF Titration):** `s2_fig2_hg002_chr20_truth`, `s4_tier3_truth_and_qual_sweep`, `s5_grch38_field_callers_truth`, `s6_low_vaf_titration_truth_curve`, and `evidence_certified`.
- **`Figure 4` & `Table 4` (Multicore Scaling, Chromosome 1 Span-to-RSS Invariance, and Node Packing):** `s4_perf_budget_4targets_n3`, `s4_scale_exome_3way_th1`, `s4_scale_exome_3way_th2`, `s4_scale_exome_3way_th4`, `s4_scale_exome_3way_th8`, `s4_scale_exome_3way_th16`, `s4_scale_exome_3way_th32`, `s4_span_chr1_3way_th8`, `s4_span_chr1_3way_th16`, `s4_span_chr1_3way_th32`, `s4_span_chr1_f0`, and `s4_jobs_per_node_2x4c_n3`.
- **`Figure 5` & `Table 5` (Full-Genome WGS Unchunked Contigs vs. Pre-Chunked Callable Intervals):** `s5_wgs_hg002_full_contigs_3way_th32`, `s5_wgs_hg002_full_giab_3way_th32`, `s5_wgs_hg008_full_contigs_3way_th32`, and `s5_wgs_hg008_full_giab_3way_th32`.
- **`Figure 6` & `Table 6` (10-Chemistry Chr20 Panel, `-f 0.0` Noise Floor, and 70-Cell GA4GH `v3.5` Matrix):** `s5_cross_platform_short_read_chr20_panel` and `s6_ga4gh_v35_stratification_matrix`.
- **`Figure 7` & `Table 7` (Next-Generation `GIAB v5.0q T2T-Q100`, `CMRG v1.00`, `Mosaic v1.10`, and 11-Caller `HG008 v0.4` Benchmark):** `s6_hg002_v50q_cmrg_mosaic_eval` and `s6_hg008_modern_ml_and_commercial_callers`.
- **`Figure 8` & `Table 8` (Software Supply Chain and 10-Point `vardictcpp` Divergence Audit):** `s5_h4_maintenance_security_matrix` and `s2_tierp_mode_parity_sweep`.

### 1.1 Domain 1: Ceiling, Multi-Region Boundary, and Tile-Invariance Cells (`s1_*`)
Domain 1 comprises `4` cells (`s1_ceiling_amplicon`, `s1_ceiling_hg008`, `s1_gate7_ref_multiregion`, and `s1_h1a_simple_tiles`):
```bash
# Cell 1: s1_ceiling_amplicon (EGFR/MRE11a amplicon -a {{m:amp_param}}, {{m:s1_amplicon_vars}} variants)
vardictrust -G hg19.fa -f {{m:f_std}} -N SAMPLE -b amplicon.bam \
  -a {{m:amp_param}} --fisher -c 1 -S 2 -E 3 -g 4 --vcf amplicon.bed > rust_amp.vcf

# Cell 2: s1_ceiling_hg008 (NIST GIAB HG008 somatic slice, {{m:s1_hg008_vars}} variants)
vardictrust -G GRCh38_GIABv3.fasta -f {{m:f_std}} -N "HG008-T|HG008-N" \
  -b "HG008-T.bam|HG008-N.bam" -th 8 --fisher -z -c 1 -S 2 -E 3 -g 4 \
  --vcf hg008_slice.bed > rust_hg008.vcf

# Cell 3: s1_gate7_ref_multiregion (HG002 chr20 multi-region boundary gate, {{m:s1_gate7_vars}} variants)
vardictrust -G GRCh38.fa -f {{m:f_std}} -N HG002 -b hg002_chr20.bam \
  -th 8 --fisher -c 1 -S 2 -E 3 -g 4 --vcf chr20_multiregion.bed > rust_chr20_multiregion.vcf

# Cell 4: s1_h1a_simple_tiles (Tile-size invariance: {{m:s1_tile_sizes_bp_list}} bp -> {{m:s1_tile_vars}} variants)
for TILE in 500000 1000000 2000000 10000000; do
  VARDICTRUST_TILE_BP="${TILE}" vardictrust -G GRCh38.fa -f {{m:f_std}} -N HCC1395 \
    -b HCC1395.bam -th 8 --fisher -c 1 -S 2 -E 3 -g 4 --vcf seqc2_slice.bed > "tile_${TILE}.vcf"
done
```

### 1.2 Domain 2: Wolff et al. AML Panel Replication, 13-Mode Sweep, Field Callers, and `-R chr1` (`s2_*`)
Domain 2 comprises `17` cells replicating and extending Wolff et al. (2026) [11]:
- **12 AML Panel Cells (`n = 3` technical replicates each at `-th 1` and `-th 8`):** `s2_tierp_v2_SRR15006375_th1`, `s2_tierp_v2_SRR15006375_th8`, `s2_tierp_v2_SRR15006376_th1`, `s2_tierp_v2_SRR15006376_th8`, `s2_tierp_v2_SRR15006386_th1`, `s2_tierp_v2_SRR15006386_th8`, `s2_tierp_v2_SRR15006540_th1`, `s2_tierp_v2_SRR15006540_th8`, `s2_tierp_v2_SRR8657348_th1`, `s2_tierp_v2_SRR8657348_th8`, `s2_tierp_v2_somatic_76_86_th1`, and `s2_tierp_v2_somatic_76_86_th8`.
- **Contiguous Chromosome 1 (`-R chr1:1-249250621`, `{{m:s2_rchr1_bp}}` bp on `{{m:s2_rchr1_sample}}`):** `s2_table3_R_chr1_java_200g_n3` (`JAVA_OPTS="-Xmx4g"` OOM `rc=3` at `{{m:s2_rchr1_java4g_oom_wall_s}}` s; `JAVA_OPTS="-Xmx200g"` completes in `{{m:s2_rchr1_java200g_wall_s}}` s at `{{m:s2_rchr1_java200g_rss_gib}}` GiB RSS) and `s2_table3_R_chr1_rust_and_cpp` (`VarDictRust` `{{m:s2_rchr1_rust_wall_s}}` s, `{{m:s2_rchr1_rust_rss_mb}}` MB RSS vs `vardictcpp` `{{m:s2_rchr1_cpp_wall_s}}` s, `{{m:s2_rchr1_cpp_rss_gib}}` GiB RSS).
- **13-Mode Functional Sweep & Field Callers:** `s2_tierp_mode_parity_sweep` (`{{m:s2_modes_rust_match}}/{{m:s2_modes_total}}` modes), `s2_fig1_fig3_field_callers` (`FreeBayes`, `LoFreq`, `Mutect2`, `Octopus`, `Strelka2`), and `s2_fig2_hg002_chr20_truth` (`{{m:s2_fig2_parity_vars}}` candidate variants).

```bash
# AML single-sample & paired somatic replication (s2_tierp_v2_*):
/usr/bin/time -v vardictrust -G hg19.fa -f {{m:f_std}} -N {{m:s2_rchr1_sample}} \
  -b {{m:s2_rchr1_sample}}.bam -th 8 --fisher -z -c 1 -S 2 -E 3 -g 4 --vcf {{m:s2_rchr1_sample}}.bed > {{m:s2_rchr1_sample}}.vcf

# Table 3 contiguous chromosome 1 (-R chr1:1-249250621 on {{m:s2_rchr1_sample}}, s2_table3_R_chr1_*):
/usr/bin/time -v vardictrust -G hg19.fa -f {{m:f_std}} -N {{m:s2_rchr1_sample}} \
  -b {{m:s2_rchr1_sample}}.bam -th 8 -R chr1:1-249250621 > rust_chr1.var
JAVA_OPTS="-Xmx200g" /usr/bin/time -v VarDict -G hg19.fa -f {{m:f_std}} -N {{m:s2_rchr1_sample}} \
  -b {{m:s2_rchr1_sample}}.bam -th 8 -R chr1:1-249250621 > java_chr1.var
```

### 1.3 Domain 3: Full-Sample {{m:s3_rna_bam_gib}} GiB RNA-Seq Splice Calling (`s3_rna_splice_parity_n3`)
Cell `s3_rna_splice_parity_n3` (`n = 3` technical replicates) evaluates `-i` RNA splice-junction mode on the `{{m:s3_rna_bam_gib}}` GiB HG002 GM24385 RNA-seq alignment (`{{m:s3_rna_splice_rows}}` `.var` and `.splice.tsv` records, `{{m:s3_rna_rust_wall_s}}` s vs `{{m:s3_rna_java_wall_s}}` s, `{{m:s3_rna_speedup_vs_java}}`× speedup, `{{m:s3_rna_rss_ratio_vs_java}}`× RSS reduction):
```bash
/usr/bin/time -v vardictrust -G GRCh38.fa -f {{m:f_std}} -N HG002_RNA \
  -b HG002_GM24385_RNA.bam -i -th 8 -c 1 -S 2 -E 3 -g 4 gencode_exons.bed > rust_rna.var
```

### 1.4 Domain 4: Whole-Exome Scaling, Chromosome 1 Span-to-RSS, Node Packing, and `QUAL` ROC (`s4_*`)
Domain 4 comprises `15` cells (`s4_perf_budget_4targets_n3`, `s4_scale_exome_3way_th1`, `s4_scale_exome_3way_th2`, `s4_scale_exome_3way_th4`, `s4_scale_exome_3way_th8`, `s4_scale_exome_3way_th16`, `s4_scale_exome_3way_th32`, `s4_span_chr1_3way_th8`, `s4_span_chr1_3way_th16`, `s4_span_chr1_3way_th32`, `s4_span_chr1_f0`, `s4_jobs_per_node_2x4c_n3`, `s4_multimodal_and_cram_twin`, `s4_tier3_truth_and_qual_sweep`, and `evidence_certified`):
```bash
# Whole-exome scaling across 1, 2, 4, 8, 16, 32 threads on SEQC2 HCC1395 (s4_scale_exome_3way_th1..th32):
for TH in 1 2 4 8 16 32; do
  /usr/bin/time -v vardictrust -G GRCh38.fa -f {{m:f_std}} -N HCC1395 \
    -b HCC1395.bam -th "${TH}" --fisher -z -c 1 -S 2 -E 3 -g 4 \
    --vcf seqc2_exome.bed > "seqc2_th${TH}.vcf"
done

# Chromosome 1 span-to-RSS invariance across 1M, 10M, 100M, 249M bp on SEQC2 HCC1395 (s4_span_chr1_3way_th8..th32):
for SPAN in 1000000 10000000 100000000 248956422; do
  /usr/bin/time -v vardictrust -G GRCh38.fa -f {{m:f_std}} -N HCC1395 \
    -b HCC1395.bam -th 32 -R "chr1:1-${SPAN}" > "chr1_${SPAN}.var"
done

# All-positions (-f {{m:f_zero}}) pileup stress on chr1 10 Mbp (s4_span_chr1_f0) and 2x4c node packing (s4_jobs_per_node_2x4c_n3):
/usr/bin/time -v vardictrust -G GRCh38.fa -f {{m:f_zero}} -N HCC1395 -b HCC1395.bam -th 8 -R chr1:10000000-20000000 > chr1_f0.var
```

### 1.5 Multi-Regime Parameter Sweep: Comprehensive 87-Arm `.var` and `.vcf` Parity Sweep (`sweep_*`)
The `12` multi-regime sweep cells (`sweep_4targets_full`, `sweep_R_chr1_full`, `sweep_hg002_chr20_exome`, `sweep_multimodal_and_cram_twin`, `sweep_rna_hg002_52gb`, `sweep_span_chr1_f0`, `sweep_tierp_13modes`, `sweep_tierp_SRR15006540`, `sweep_tierp_SRR8657348`, `sweep_tierp_somatic_76_86`, `sweep_ultradeep_mrd_stress`, and `sweep_wgs_mrd_matrix`) validate `.var` and `.vcf` equivalence across `{{m:sweep_total_arms}}` comparative arms (`{{m:sweep_match_arms}}/{{m:sweep_total_arms}}` `MATCH`, `{{m:sweep_match_pct}}`, `Supplementary Table S2`):
```bash
# Ultra-deep 45000x ctDNA MRD stress (-f {{m:f_mrd_ud}} and -f {{m:f_zero}} with --fisher, sweep_ultradeep_mrd_stress):
vardictrust -G GRCh38.fa -f {{m:f_mrd_ud}} -N SAMPLE -b ultradeep_45000x.bam --fisher -c 1 -S 2 -E 3 -g 4 panel.bed > rust_ultradeep_f0005.var
# WGS & exome MRD matrix (-f {{m:f_mrd_wgs}} and -f {{m:f_zero}}, sweep_wgs_mrd_matrix):
vardictrust -G GRCh38.fa -f {{m:f_mrd_wgs}} -N HG005 -b hg005_300x.bam --fisher -c 1 -S 2 -E 3 -g 4 --vcf mrd_targets.bed > rust_hg005_300x_mrd.vcf
```

### 1.6 Domain 5: Full-Genome WGS (Unchunked Full Chromosomes vs. Pre-Chunked Callable Intervals), 10-Platform Chr20 Panel, and Supply Chain (`s5_*`)
Domain 5 comprises `7` cells (`s5_wgs_hg002_full_contigs_3way_th32`, `s5_wgs_hg002_full_giab_3way_th32`, `s5_wgs_hg008_full_contigs_3way_th32`, `s5_wgs_hg008_full_giab_3way_th32`, `s5_grch38_field_callers_truth`, `s5_cross_platform_short_read_chr20_panel`, and `s5_h4_maintenance_security_matrix`):
```bash
# Unchunked full chromosomes (chr1..chr22, s5_wgs_hg002_full_contigs_3way_th32 & s5_wgs_hg008_full_contigs_3way_th32):
timeout 10800 /usr/bin/time -v vardictrust -G GRCh38_GIABv3.fasta -f {{m:f_std}} -N HG002 \
  -b HG002_2x250_131GB.bam -th 32 --fisher -z -c 1 -S 2 -E 3 -g 4 --vcf grch38_22autosomes_full.bed > hg002_unchunked.vcf

timeout 10800 /usr/bin/time -v vardictrust -G GRCh38_GIABv3.fasta -f {{m:f_std}} -N "HG008-T|HG008-N" \
  -b "HG008-T_195x.bam|HG008-N_169x.bam" -th 32 -U --fisher -z -c 1 -S 2 -E 3 -g 4 --vcf grch38_22autosomes_full.bed > hg008_unchunked.vcf

# Pre-chunked <=50 kbp GIAB callable intervals (s5_wgs_hg002_full_giab_3way_th32 & s5_wgs_hg008_full_giab_3way_th32):
timeout 10800 /usr/bin/time -v vardictrust -G GRCh38_GIABv3.fasta -f {{m:f_std}} -N HG002 \
  -b HG002_2x250_131GB.bam -th 32 --fisher -z -c 1 -S 2 -E 3 -g 4 --vcf HG002_GRCh38_v421_50kb.bed > hg002_giab50kb.vcf

timeout 10800 /usr/bin/time -v vardictrust -G GRCh38_GIABv3.fasta -f {{m:f_std}} -N "HG008-T|HG008-N" \
  -b "HG008-T_195x.bam|HG008-N_169x.bam" -th 32 -U --fisher -z -c 1 -S 2 -E 3 -g 4 --vcf HG008_GRCh38_v0.3_50kb.bed > hg008_giab50kb.vcf
```

### 1.7 Domain 6: GA4GH v3.5 Stratification, GIAB v5.0q/CMRG/Mosaic, HG008 v0.4 11-Caller Matrix, and Low-VAF Titration (`s6_*`)
Domain 6 comprises `4` cells (`s6_ga4gh_v35_stratification_matrix`, `s6_hg002_v50q_cmrg_mosaic_eval`, `s6_hg008_modern_ml_and_commercial_callers`—including the official NIST GIAB `Liss_lab` somatic callsets for `DRAGEN v4.2.4`, `DeepSomatic v1.6.0`, and `ClairS v0.4.1`, alongside `Lancet2 (NYGC)`, `Hartwig SAGE`, and `Ultima DeepVariant Somatic`—and `s6_low_vaf_titration_truth_curve`):
```bash
# Controlled in silico tumor-purity titration (20%, 10%, 5%, 2%, 1% purity at -f {{m:f_titration}}, s6_low_vaf_titration_truth_curve):
for PURITY in 20 10 5 2 1; do
  vardictrust -G GRCh38_GIABv3.fasta -f {{m:f_titration}} -N "HG008_P${PURITY}|HG008-N" \
    -b "hg008_titration_p${PURITY}.bam|HG008-N.bam" -th 16 --fisher -z -c 1 -S 2 -E 3 -g 4 \
    --vcf hg008_slice.bed > "rust_titration_p${PURITY}.vcf"
done

# Next-gen GIAB v5.0q T2T-Q100, CMRG v1.00, and Mosaic v1.10 evaluation (s6_hg002_v50q_cmrg_mosaic_eval):
rtg vcfeval -b HG002_GRCh38_v5.0q_T2T.vcf.gz -e HG002_GRCh38_v5.0q_T2T.bed \
  -c hg002_unchunked.vcf.gz -t GRCh38.sdf -o eval_v50q_unchunked
```

---

## 2. Supplementary Note 1: Analytical Derivation of JVM (`VarDictJava`) and `vardictcpp` Heap Retention (Unchunked Chromosomes vs. Pre-Chunked Callable Intervals)

In `VarDictJava` v1.8.3 (`CigarParser.java` and `VariationRealigner.java`), the pileup engine processes each BED region $[S, E]$ of length $L = E - S + 1$ by instantiating two region-wide heap structures before emitting any variant row:
1. A reference base map `Reference.referenceSequences` (`Map<Integer, Character>`) holding $O(L)$ boxed `java.lang.Integer` (`16` bytes) and `java.lang.Character` (`16` bytes) objects plus `HashMap$Node` headers (`32` bytes), consuming $\approx 64\text{–}96$ bytes per reference base ($\approx 16\text{–}24$ GiB for human `chr1` of length $L = {{m:s2_rchr1_bp}}$ bp per worker thread).
2. A position-indexed variation map `Map<Integer, VariationMap<String, Variation>> nonInsertionVariants` and `insertionVariants`, where each covered coordinate $p \in [S, E]$ allocates a `VariationMap` and one `Variation` object per observed allele. At sequencing depth $D$, the active heap footprint across $T$ worker threads scales as:
   $$\text{RSS}_{\text{Java / C++} }(L, T) = \Theta\!\left(\sum_{t=1}^{T} L_t \cdot (c_{\text{ref} } + c_{\text{var} } \cdot \bar{a}(D))\right)$$
   where $L_t$ is the region length assigned to thread $t$ and $\bar{a}(D)$ is the mean allele count per locus. When $L_t \le {{m:giab_arm_b_max_kbp}}$ kbp (pre-chunked intervals), garbage collection reclaims objects between intervals (`{{m:s5_hg002_giab_java_rss_gb}}` GB on HG002; `{{m:s5_hg008_giab_java_rss_gb}}` GB with timeout `rc=124` on HG008). When $L_t \in [50, 249]$ Mbp (unchunked chromosomes or `-R chr1:1-249250621`), live objects remain reachable until the chromosome finishes, driving `VarDictJava` to `{{m:s2_rchr1_java200g_rss_gib}}` GiB on `chr1` and `>{{m:s5_hg002_contigs_java_rss_gb}}` GB (`rc=3` `OutOfMemoryError`) on multi-threaded full genomes (`Table 2`, `Table 5`). `vardictcpp` similarly retains per-region C++ hash tables across each interval, consuming `{{m:s2_rchr1_cpp_rss_gib}}` GiB on `chr1`, `{{m:s5_hg002_contigs_cpp_rss_gb}}` GB on HG002 unchunked chromosomes, and `{{m:s5_hg008_contigs_cpp_rss_gb}}` GB (`rc=124` timeout) on HG008 unchunked chromosomes unless an external script pre-splits the BED file (`CPP_CHUNK={{m:cpp_chunk_bp}}`).

`VarDictRust` leaves intervals $\le {{m:s1_tile_1m}}$ bp untiled and partitions any region $L > {{m:s1_tile_1m}}$ bp into tiles of default width $W = {{m:rust_default_tile_kbp}}$ kbp with halo $H = {{m:rust_halo_kbp}}$ kbp (`{{m:rust_splice_halo_kbp}}` kbp in `-i` splice mode; `Figure 1b`), bounding memory to:
$$\text{RSS}_{\text{VarDictRust} }(L, T) = O\!\left(T \cdot (W + 2H)\right) \perp L$$
Increasing $L$ by `250`× on chromosome 1 (`1` Mbp to `249` Mbp, `s4_span_chr1_3way_th8`..`th32`) changes `VarDictRust` RSS by `{{m:s4_span_th8_ratio_249m_10m}}`×–`{{m:s4_span_th32_ratio_249m_10m}}`× (`{{m:s4_span_th8_249m_rust_gib}}`–`{{m:s4_span_th32_249m_rust_gib}}` GiB), allowing unchunked full-genome calling in `{{m:s5_hg002_contigs_rust_rss_gb}}` GB (HG002) and `{{m:s5_hg008_contigs_rust_rss_gb}}` GB (HG008).

---

## 3. Supplementary Note 2: Ten-Point Technical Audit of `vardictcpp` Divergences (`CPP-1` through `CPP-10`)

Across cells `s2_tierp_mode_parity_sweep`, `s5_h4_maintenance_security_matrix`, `s5_wgs_hg002_full_contigs_3way_th32`, and `s6_low_vaf_titration_truth_curve`, we cataloged ten reproducibility and capability gaps in `{{m:cpp_version}}` (`Table 8`, `Supplementary Table S6`):
- **`CPP-1` (VCF Header & `INFO/END` Schema):** `vardictcpp` omits the `##INFO=<ID=END>` header definition and diverges from `var2vcf_valid.pl` / `var2vcf_paired.pl` VCF 4.3 schema formatting.
- **`CPP-2` (CRISPR (`-J`/`-j`) & Structural Variant Modes):** `vardictcpp` sets `--nosv` (`-U`) by default and omits CRISPR cleavage-site (`-J`/`-j`) and structural variant (`StructuralVariantsProcessor`) assembly (`1/13` modes failed in `s2_tierp_mode_parity_sweep`, `{{m:s2_modes_cpp_exact}}/{{m:s2_modes_total}}` exact).
- **`CPP-3` (Amplicon (`-a`) Output Format):** In amplicon primer-clipping mode (`-a {{m:amp_param}}`, `s1_ceiling_amplicon`), `vardictcpp` emits raw tab-separated `.var` text instead of VCF (`0` VCF lines), breaking downstream `bcftools` and `rtg vcfeval` consumers.
- **`CPP-4` (Concatenated CLI Flag Parsing):** Passing bundled single-character flags (e.g., `-z -c 1` combined or legacy `VarDictJava` bundled options) causes `vardictcpp` to fail CLI parsing with exit code `120`.
- **`CPP-5` (50 kbp Chunk Boundary Deduplication):** Without external BED splitting (`--chunk {{m:cpp_chunk_bp}}` / `CPP_CHUNK={{m:cpp_chunk_bp}}`), `vardictcpp` consumes `{{m:s2_rchr1_cpp_rss_gib}}` GiB on `chr1`, `{{m:s5_hg002_contigs_cpp_rss_gb}}` GB on HG002 unchunked chromosomes, and `{{m:s5_hg008_contigs_cpp_rss_gb}}` GB (`rc=124` timeout) on HG008 unchunked chromosomes. When intervals are chunked, `vardictcpp` lacks boundary deduplication, emitting `+{{m:s5_hg002_contigs_cpp_excess}}` duplicate/shifted calls on HG002 unchunked chromosomes (`{{m:s5_hg002_contigs_cpp_vars}}` vs `{{m:s5_hg002_contigs_rust_vars}}`) and `+{{m:s5_hg002_giab_cpp_diff_vars}}` on pre-chunked intervals (`{{m:s5_hg002_giab_cpp_vars}}` vs `{{m:s5_hg002_giab_vars}}`).
- **`CPP-6` (Missing `##contig` Header Lines):** `vardictcpp` omits `##contig=<ID=...,length=...>` lines from VCF headers (`VarDictJava` issue `#320`), causing strict VCF validators and `rtg vcfeval` header checks to warn or fail.
- **`CPP-7` (GRCh37 / b37 Chromosome Naming):** `vardictcpp` hardcodes `chr`-prefix assumptions in internal chromosome handling, breaking or diverging on `b37` / `GRCh37` numeric contigs (`1`..`22`, `X`, `Y`, `MT`).
- **`CPP-8` (BAM/CRAM Reference & Decode):** `vardictcpp` diverges on CRAM reference decoding (`-G ref.fa`) and requires dynamic linkage against external `libhts`, `libcurl`, `liblzma`, `libbz2`, and `libdeflate`, whereas `VarDictRust` natively decodes BAM and CRAM in a `{{m:s5_h4_musl_mb}}` MB zero-dependency static `musl` binary (`s4_multimodal_and_cram_twin`).
- **`CPP-9` (MNP Phasing, Tie-Break & Somatic SV Drops):** `vardictcpp` diverges from `VarDictJava` on complex MNP phasing, realignment tie-breaking, and somatic variant retention. At `20%` tumor purity (`s6_low_vaf_titration_truth_curve`), `vardictcpp` drops somatic indel $F_1$ to `{{m:s6_tit_p20_vardictcpp_indel_f1}}` (recall `{{m:s6_tit_p20_vardictcpp_indel_r_pct}}`) and SNV $F_1$ to `{{m:s6_tit_p20_vardictcpp_snv_f1}}` (recall `{{m:s6_tit_p20_vardictcpp_snv_r_pct}}`), whereas `VarDictRust` matches `VarDictJava` (`Indel` $F_1 = {{m:s6_tit_p20_VarDictRust_indel_f1}}$, recall `{{m:s6_tit_p20_VarDictRust_indel_r_pct}}`; `SNV` $F_1 = {{m:s6_tit_p20_VarDictRust_snv_f1}}$, recall `{{m:s6_tit_p20_VarDictRust_snv_r_pct}}`).
- **`CPP-10` (RNA-Seq Splice Mode (`-i`.var / `.splice.tsv`):** `vardictcpp` does not implement `-i` exon-junction splice-site accounting (`s3_rna_splice_parity_n3`, `{{m:s3_rna_bam_gib}}` GiB RNA-seq BAM, `{{m:s3_rna_splice_rows}}` records).

---

## 4. Supplementary Note 3: Biostatistical Methodology, Exact Clopper-Pearson Intervals, and `--fix` Governance

### 4.1 Exact Clopper-Pearson Binomial Confidence Intervals
For binary oracle-concordance and completion rates ($k$ successes out of $n$ independent cells/arms), we compute exact two-sided $1 - \alpha = 95\%$ Clopper-Pearson binomial CIs $[p_{\text{lo} }, p_{\text{hi} }]$ via the inverse Beta distribution:
$$p_{\text{lo} } = \text{Beta}^{-1}\!\left(\frac{\alpha}{2};\, k,\, n - k + 1\right), \quad p_{\text{hi} } = \text{Beta}^{-1}\!\left(1 - \frac{\alpha}{2};\, k + 1,\, n - k\right)$$
with $p_{\text{lo} } = 0$ when $k = 0$ and $p_{\text{hi} } = 1$ when $k = n$:
- Complete `{{m:total_cells}}`-cell verification campaign (`{{m:total_cells}}/{{m:total_cells}}` exit `0`): `100%` (`95%` CI `{{m:cp95_60_60}}`)
- `13`-mode functional parity sweep (`{{m:s2_modes_rust_match}}/{{m:s2_modes_total}}`): `100%` (`95%` CI `{{m:cp95_13_13}}`)
- Extended coverage functional modalities (`{{m:s4_ext_cov_modalities}}`): `100%` (`95%` CI `{{m:cp95_17_17}}`)
- Multimodal clinical and edge-case slices (`{{m:s4_multimodal_slices}}`): `100%` (`95%` CI `{{m:cp95_18_18}}`)
- Cross-platform chr20 runs where `VarDictJava` completed (`{{m:s5_chr20_java_completed_runs}}/{{m:s5_chr20_java_completed_runs}}`): `100%` (`95%` CI `{{m:cp95_chr20_java_completed}}`)
- Cross-platform `-f {{m:f_zero}}` somatic noise-floor `STRICT-BYTES` sweep (`{{m:s5_f0_noise_floor_strict_match}}`): `100%` (`95%` CI `{{m:cp95_noise_18_18}}`)
- GA4GH `v3.5` `{{m:s6_ga4gh_cells}}`-cell stratification matrix (`{{m:s6_ga4gh_cells}}/{{m:s6_ga4gh_cells}}`): `100%` (`95%` CI `{{m:cp95_70_70}}`)
- Multi-regime parity sweep across `12` suites (`{{m:optc_match_arms}}/{{m:optc_total_arms}}`): `{{m:optc_match_pct}}%` (`95%` CI `{{m:cp95_optc_sweep}}`; includes untiled $\le {{m:s1_tile_1m}}$ bp BED suites and single-window `--tile-size` runs on `sweep_span_chr1_f0` [`{{m:sw_span_f0_rows}}` `.var` / `{{m:s4_span_f0_vars}}` `.vcf`] and `sweep_R_chr1_full` [`{{m:sw_rchr1_rows}}` `.var`]; `Supplementary Table S2`)
- NIST Mosaic SNV `v1.10` sensitivity (`{{m:s6_mosaic_rust_tp}}/{{m:s6_mosaic_truth_total}}`): `{{m:s6_mosaic_rust_r_pct}}` (`95%` CI `{{m:cp95_84_85}}`)

### 4.2 Non-Parametric Runtime and Memory Hypothesis Testing
Across the $N = 12$ AML sample-by-thread cells (`s2_tierp_v2_*`, `n = 3` technical replicates each), we report the mean and replicate range (`min–max`, `Supplementary Table S3`). Comparing `VarDictRust` against `VarDictJava` across $N = 12$ yields two-sided Mann-Whitney $U = {{m:s2_tierp_mwu_wall_u}}$ ($p = {{m:s2_tierp_mwu_wall_p}}$) for wall time, $U = {{m:s2_tierp_mwu_rss_u}}$ ($p = {{m:s2_tierp_mwu_rss_p}}$) for peak RSS, and an exact paired sign-flip permutation p-value of $p = 2^{-12} = {{m:s2_tierp_sign_perm_p}}$.

### 4.3 Release Governance: `v1.0.x` Parity Versus `v1.1+` Opt-In `--fix` Flags
All `v1.0.x` releases (`v{{m:version}}`) enforce strict oracle parity on valid inputs. Output-changing upstream fixes (`#92` `VD > DP`, `#62` over-filtering, `#292` `SHIFT3`, `#342` `<dup-N>` symbolic alleles) are gated behind opt-in `--fix <id>` flags starting in `v1.1.0`, while parity-restoring fixes (`splice-tile-ownership`, `mate-window-ownership`) ship as verified defaults.

---

## 5. Supplementary Tables (`Supplementary Table S1` – `Supplementary Table S6`) and Complete 60-Cell Verification Ledger

{{table:supp_table1}}

---

{{table:supp_table2}}

---

{{table:supp_table3}}

---

{{table:supp_table4}}

---

{{table:supp_table5}}

---

{{table:supp_table6}}

---

## 6. References

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
11. Wolff L, et al. vardictcpp: a high-performance C++ reimplementation of VarDict. *Bioinformatics.* 2026;42(2):btag042. doi:10.1093/bioinformatics/btag042.
12. Cleary JG, Braithwaite R, Gaastra K, Hilbush BS, Inglis S, Irvine SA, et al. Comparing variant call files for performance benchmarking of next-generation sequencing variant calling pipelines. *bioRxiv* [Preprint]. 2015;doi:10.1101/023754.
13. Wagner J, Olson ND, Harris L, Khan Z, Farek J, Mahmoud M, et al. Benchmarking challenging small variants with linked and long reads. *Cell Genomics.* 2022;2(5):100128. doi:10.1016/j.xgen.2022.100128.
14. Rhie A, Nurk S, Cechova M, Hoyt SJ, Taylor DJ, Altemose N, et al. The complete sequence of a human Y chromosome. *Nature.* 2023;621(7978):344–54. doi:10.1038/s41586-023-06457-y.
15. Wagner J, Olson ND, Harris L, Khan Z, Farek J, Mahmoud M, et al. Curated variation benchmarks for challenging medically relevant autosomal genes. *Nat Biotechnol.* 2022;40(5):672–80. doi:10.1038/s41587-021-01158-1.
16. Nurk S, Koren S, Rhie A, Rautiainen M, Bzikadze AV, Mikheenko A, et al. The complete sequence of a human genome. *Science.* 2022;376(6588):44–53. doi:10.1126/science.abj6987.
17. Krusche P, Trigg L, Boutros PC, Mason CE, De La Vega FM, Moore BL, et al. Best practices for benchmarking germline small-variant calls in human genomes. *Nat Biotechnol.* 2019;37(5):555–60. doi:10.1038/s41587-019-0054-x.
18. Zook JM, Catoe D, McDaniel J, Vang L, Spies N, Sidow A, et al. Extensive sequencing of seven human genomes to characterize benchmark reference materials. *Sci Data.* 2016;3:160025. doi:10.1038/sdata.2016.25.
19. McDaniel J, Olson ND, Wagner J, et al. Authoritative characterization of somatic variation in the pancreatic ductal adenocarcinoma reference cell line HG008. *NIST Genome in a Bottle Consortium Technical Report.* 2025;v0.4.
20. Daniels C, et al. Comprehensive evaluation of germline variant callers on human chromosome 20 exome sequencing benchmarks. *Brief Bioinform.* 2024;25(3):bbae142.
21. Behera S, Catreux S, Rossi M, Truong S, Huang Z, Ruehle M, et al. Comprehensive genome analysis and variant detection at scale using DRAGEN. *Nat Biotechnol.* 2024;43:1177–89. doi:10.1038/s41587-024-02382-1.
22. Park J, Cook DE, Chang PC, Kolesnikov A, Brambrink L, Miga KH, et al. DeepSomatic: accurate somatic small variant discovery for multiple sequencing technologies. *bioRxiv* [Preprint]. 2024;doi:10.1101/2024.08.16.608331.
23. Zheng Z, Su J, Chen L, Lee Yan V, Lam TW, Luo R. ClairS: a deep-learning method for long-read somatic small variant calling. *bioRxiv* [Preprint]. 2023;doi:10.1101/2023.08.17.553778.
24. Narzisi G, Corvelo A, Arora K, Bergmann EA, Shah M, Musunuri R, et al. Genome-wide somatic variant calling using localized colored de Bruijn graphs. *Commun Biol.* 2018;1:20. doi:10.1038/s42003-018-0023-9.
25. Priestley P, Baber J, Lolkema MP, Steeghs N, de Bruijn E, Shale C, et al. Pan-cancer whole-genome analyses of metastatic solid tumours. *Nature.* 2019;575(7781):210–6. doi:10.1038/s41586-019-1689-y.
26. Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, et al. Twelve years of SAMtools and BCFtools. *GigaScience.* 2021;10(2):giab008. doi:10.1093/gigascience/giab008.
