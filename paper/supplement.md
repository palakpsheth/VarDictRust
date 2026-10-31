# Supplementary Material: VarDictRust: A Memory-Bounded, Oracle-Equivalent Rust Consolidation of the VarDict Variant-Calling Pipeline for Deep Whole-Genome and Ultra-Low-VAF Sequencing

**Authors:** [Palak Sheth](https://www.linkedin.com/in/palak-p-sheth/)<sup>1,*</sup>  
<sup>1</sup>Independent Researcher ([palakpsheth@gmail.com](mailto:palakpsheth@gmail.com) | [https://www.linkedin.com/in/palak-p-sheth/](https://www.linkedin.com/in/palak-p-sheth/))  
**Engine Release:** `VarDictRust` `v1.0.6` (commit `2a8e9a3`, engine fingerprint `5865a2a0b24a8c17`, performance fingerprint `c490c36317b8f25b`)  
**Reference Oracles:** `VarDictJava v1.8.3 (0886261)` and `vardictcpp ae93fe39dc (Wolff 2026)`  
**Cloud Benchmark Host:** `GCP n2-standard-32 (Intel Xeon Platinum 8373C Ice Lake, 32 vCPUs, 128 GiB RAM, 4x375 GB NVMe RAID0)` (`1.5` TB NVMe RAID0 scratch; exploratory staging also verified on `c3d-standard-32`)

---

## Table of Contents

1. **Supplementary Methods: Exact Command-Line Invocations Across All 60 Benchmark Cells**
    - 1.0 Manuscript Figure and Table to 60-Cell Provenance Crosswalk
    - 1.1 Domain 1: Ceiling, Multi-Region Boundary, and Tile-Invariance Cells (`4` cells)
    - 1.2 Domain 2: Wolff AML Panel Replication, 13-Mode Sweep, Field Callers, and `-R chr1` (`17` cells)
    - 1.3 Domain 3: Full-Sample 51.8 GiB RNA-Seq Splice Calling (`1` cell)
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

All `60` verification cells ran on `GCP n2-standard-32 (Intel Xeon Platinum 8373C Ice Lake, 32 vCPUs, 128 GiB RAM, 4x375 GB NVMe RAID0)` (`Ubuntu 24.04 LTS`, Linux kernel `6.8`, `1.5` TB 4-way NVMe RAID0 scratch mounted at `/var/tmp/vardict-paper`). They comprise `59` cloud compute cells plus `1` certified local-to-cloud evidence bundle `evidence_certified`. Timed runs used `/usr/bin/time -v` across `n = 3` independent technical replicates. The cap was `timeout 10800` (`3` hours, exit code `rc=124` on timeout; `rc=3` on JVM `OutOfMemoryError`).

### 1.0 Manuscript Figure and Table to 60-Cell Provenance Crosswalk

The crosswalk below maps each main-text Figure (`Figure 1`–`Figure 8`) and Table (`Table 1`–`Table 8`) to its benchmark directory in `benchmarks/gcp-study-v106/`. `Supplementary Table S1` holds the full ledger.

- **`Figure 1` & `Table 1` (Architecture, Operational Modes, and 87-Arm Equivalence Sweep):** `s1_ceiling_amplicon`, `s1_ceiling_hg008`, `s1_gate7_ref_multiregion`, `s1_h1a_simple_tiles`, `s2_tierp_mode_parity_sweep`, `s3_rna_splice_parity_n3`, `s4_multimodal_and_cram_twin`, and `evidence_certified`. The `12` `sweep_*` cells add `sweep_4targets_full`, `sweep_R_chr1_full`, `sweep_hg002_chr20_exome`, `sweep_multimodal_and_cram_twin`, `sweep_rna_hg002_52gb`, `sweep_span_chr1_f0`, `sweep_tierp_13modes`, `sweep_tierp_SRR15006540`, `sweep_tierp_SRR8657348`, `sweep_tierp_somatic_76_86`, `sweep_ultradeep_mrd_stress`, and `sweep_wgs_mrd_matrix`.

- **`Figure 2` & `Table 2` (Wolff AML Cohort, Field Callers, and Contiguous `-R chr1`):** The `12` AML panel cells are `s2_tierp_v2_SRR15006375_th1`, `s2_tierp_v2_SRR15006375_th8`, `s2_tierp_v2_SRR15006376_th1`, `s2_tierp_v2_SRR15006376_th8`, `s2_tierp_v2_SRR15006386_th1`, and `s2_tierp_v2_SRR15006386_th8`. They continue with `s2_tierp_v2_SRR15006540_th1`, `s2_tierp_v2_SRR15006540_th8`, `s2_tierp_v2_SRR8657348_th1`, `s2_tierp_v2_SRR8657348_th8`, `s2_tierp_v2_somatic_76_86_th1`, and `s2_tierp_v2_somatic_76_86_th8`. Contiguous `-R chr1` adds `s2_table3_R_chr1_java_200g_n3` and `s2_table3_R_chr1_rust_and_cpp`. Field callers add `s2_fig1_fig3_field_callers`.

- **`Figure 3` & `Table 3` (Germline/Somatic Truth Concordance and `QUAL` Operating Points):** `s2_fig2_hg002_chr20_truth`, `s4_tier3_truth_and_qual_sweep`, `s5_grch38_field_callers_truth`, and `evidence_certified`. The full unfiltered-vs-optimal `QUAL` ROC pairs and the low-VAF titration detail matrix (`s6_low_vaf_titration_truth_curve`) live in `Supplementary Table S4`.

- **`Figure 4` & `Table 4` (Multicore Scaling, Chromosome 1 Span-to-RSS Invariance, and Node Packing):** `s4_perf_budget_4targets_n3`, `s4_scale_exome_3way_th1`, `s4_scale_exome_3way_th2`, `s4_scale_exome_3way_th4`, `s4_scale_exome_3way_th8`, `s4_scale_exome_3way_th16`, and `s4_scale_exome_3way_th32`. Span invariance adds `s4_span_chr1_3way_th8`, `s4_span_chr1_3way_th16`, `s4_span_chr1_3way_th32`, and `s4_span_chr1_f0`. Node packing adds `s4_jobs_per_node_2x4c_n3`.

- **`Figure 5` & `Table 5` (Full-Genome WGS Unchunked Contigs vs. Pre-Chunked Callable Intervals):** `s5_wgs_hg002_full_contigs_3way_th32`, `s5_wgs_hg002_full_giab_3way_th32`, `s5_wgs_hg008_full_contigs_3way_th32`, and `s5_wgs_hg008_full_giab_3way_th32`.

- **`Figure 6` & `Table 6` (10-Chemistry Chr20 Panel, `-f 0.0` Noise Floor, and 70-Cell GA4GH `v3.5` Matrix):** `s5_cross_platform_short_read_chr20_panel` and `s6_ga4gh_v35_stratification_matrix`.

- **`Figure 7` & `Table 7` (Next-Generation `GIAB v5.0q T2T-Q100`, `CMRG v1.00`, `Mosaic v1.10`, and 11-Caller `HG008 v0.4` Benchmark):** `s6_hg002_v50q_cmrg_mosaic_eval` and `s6_hg008_modern_ml_and_commercial_callers`.

- **`Figure 8` & `Table 8` (Software Supply Chain and 10-Point `vardictcpp` Divergence Audit):** `s5_h4_maintenance_security_matrix` and `s2_tierp_mode_parity_sweep`.

### 1.1 Domain 1: Ceiling, Multi-Region Boundary, and Tile-Invariance Cells (`s1_*`)
Domain 1 comprises `4` cells (`s1_ceiling_amplicon`, `s1_ceiling_hg008`, `s1_gate7_ref_multiregion`, and `s1_h1a_simple_tiles`):
```bash
# Cell 1: s1_ceiling_amplicon (EGFR/MRE11a amplicon -a 10:0.95, 2 variants)
vardictrust -G hg19.fa -f 0.01 -N SAMPLE -b amplicon.bam \
  -a 10:0.95 --fisher -c 1 -S 2 -E 3 -g 4 --vcf amplicon.bed > rust_amp.vcf

# Cell 2: s1_ceiling_hg008 (NIST GIAB HG008 somatic slice, 18,426 variants)
vardictrust -G GRCh38_GIABv3.fasta -f 0.01 -N "HG008-T|HG008-N" \
  -b "HG008-T.bam|HG008-N.bam" -th 8 --fisher -z -c 1 -S 2 -E 3 -g 4 \
  --vcf hg008_slice.bed > rust_hg008.vcf

# Cell 3: s1_gate7_ref_multiregion (HG002 chr20 multi-region boundary gate, 2,572 variants)
vardictrust -G GRCh38.fa -f 0.01 -N HG002 -b hg002_chr20.bam \
  -th 8 --fisher -c 1 -S 2 -E 3 -g 4 --vcf chr20_multiregion.bed > rust_chr20_multiregion.vcf

# Cell 4: s1_h1a_simple_tiles (Tile-size invariance: 500,000, 1,000,000, 2,000,000, 10,000,000 bp -> 5,090 variants)
for TILE in 500000 1000000 2000000 10000000; do
  VARDICTRUST_TILE_BP="${TILE}" vardictrust -G GRCh38.fa -f 0.01 -N HCC1395 \
    -b HCC1395.bam -th 8 --fisher -c 1 -S 2 -E 3 -g 4 --vcf seqc2_slice.bed > "tile_${TILE}.vcf"
done
```

### 1.2 Domain 2: Wolff AML Panel Replication, 13-Mode Sweep, Field Callers, and `-R chr1` (`s2_*`)
Domain 2 comprises `17` cells replicating and extending Wolff (2026) [11]:

- **12 AML Panel Cells (`n = 3` technical replicates each at `-th 1` and `-th 8`):** `s2_tierp_v2_SRR15006375_th1`, `s2_tierp_v2_SRR15006375_th8`, `s2_tierp_v2_SRR15006376_th1`, `s2_tierp_v2_SRR15006376_th8`, `s2_tierp_v2_SRR15006386_th1`, `s2_tierp_v2_SRR15006386_th8`, `s2_tierp_v2_SRR15006540_th1`, `s2_tierp_v2_SRR15006540_th8`, `s2_tierp_v2_SRR8657348_th1`, `s2_tierp_v2_SRR8657348_th8`, `s2_tierp_v2_somatic_76_86_th1`, and `s2_tierp_v2_somatic_76_86_th8`.

- **Contiguous Chromosome 1 (`-R chr1:1-249250621`, `249,250,621` bp on `SRR8657348`):** `s2_table3_R_chr1_java_200g_n3` covers `JAVA_OPTS="-Xmx4g"` OOM `rc=3` at `31.08` s. It also covers `JAVA_OPTS="-Xmx200g"` completing in `501.83` s at `106.47` GiB RSS. `s2_table3_R_chr1_rust_and_cpp` compares `VarDictRust` (`9.75` s, `318.6` MB RSS) against `vardictcpp` (`105.89` s, `10.58` GiB RSS).

- **13-Mode Functional Sweep & Field Callers:** `s2_tierp_mode_parity_sweep` (`13/13` modes), `s2_fig1_fig3_field_callers` (`FreeBayes`, `LoFreq`, `Mutect2`, `Octopus`, `Strelka2`), and `s2_fig2_hg002_chr20_truth` (`8,233` candidate variants).

```bash
# AML single-sample & paired somatic replication (s2_tierp_v2_*):
/usr/bin/time -v vardictrust -G hg19.fa -f 0.01 -N SRR8657348 \
  -b SRR8657348.bam -th 8 --fisher -z -c 1 -S 2 -E 3 -g 4 --vcf SRR8657348.bed > SRR8657348.vcf

# Table 3 contiguous chromosome 1 (-R chr1:1-249250621 on SRR8657348, s2_table3_R_chr1_*):
/usr/bin/time -v vardictrust -G hg19.fa -f 0.01 -N SRR8657348 \
  -b SRR8657348.bam -th 8 -R chr1:1-249250621 > rust_chr1.var
JAVA_OPTS="-Xmx200g" /usr/bin/time -v VarDict -G hg19.fa -f 0.01 -N SRR8657348 \
  -b SRR8657348.bam -th 8 -R chr1:1-249250621 > java_chr1.var
```

### 1.3 Domain 3: Full-Sample 51.8 GiB RNA-Seq Splice Calling (`s3_rna_splice_parity_n3`)

Cell `s3_rna_splice_parity_n3` (`n = 3` technical replicates) evaluates `-i` RNA splice-junction mode on the `51.8` GiB HG002 GM24385 RNA-seq alignment. It produced `12,668` `.var` and `.splice.tsv` records in `4.94` s vs `23.31` s (`4.72`× speedup, `13.6`× RSS reduction):
```bash
/usr/bin/time -v vardictrust -G GRCh38.fa -f 0.01 -N HG002_RNA \
  -b HG002_GM24385_RNA.bam -i -th 8 -c 1 -S 2 -E 3 -g 4 gencode_exons.bed > rust_rna.var
```

### 1.4 Domain 4: Whole-Exome Scaling, Chromosome 1 Span-to-RSS, Node Packing, and `QUAL` ROC (`s4_*`)
Domain 4 comprises `15` cells. The scaling series runs `s4_perf_budget_4targets_n3`, `s4_scale_exome_3way_th1`, `s4_scale_exome_3way_th2`, `s4_scale_exome_3way_th4`, `s4_scale_exome_3way_th8`, `s4_scale_exome_3way_th16`, and `s4_scale_exome_3way_th32`. Span invariance adds `s4_span_chr1_3way_th8`, `s4_span_chr1_3way_th16`, `s4_span_chr1_3way_th32`, and `s4_span_chr1_f0`. Node packing, twin calling, truth sweeps, and certification add `s4_jobs_per_node_2x4c_n3`, `s4_multimodal_and_cram_twin`, `s4_tier3_truth_and_qual_sweep`, and `evidence_certified`:
```bash
# Whole-exome scaling across 1, 2, 4, 8, 16, 32 threads on SEQC2 HCC1395 (s4_scale_exome_3way_th1..th32):
for TH in 1 2 4 8 16 32; do
  /usr/bin/time -v vardictrust -G GRCh38.fa -f 0.01 -N HCC1395 \
    -b HCC1395.bam -th "${TH}" --fisher -z -c 1 -S 2 -E 3 -g 4 \
    --vcf seqc2_exome.bed > "seqc2_th${TH}.vcf"
done

# Chromosome 1 span-to-RSS invariance across 1M, 10M, 100M, 249M bp on SEQC2 HCC1395 (s4_span_chr1_3way_th8..th32):
for SPAN in 1000000 10000000 100000000 248956422; do
  /usr/bin/time -v vardictrust -G GRCh38.fa -f 0.01 -N HCC1395 \
    -b HCC1395.bam -th 32 -R "chr1:1-${SPAN}" > "chr1_${SPAN}.var"
done

# All-positions (-f 0) pileup stress on chr1 10 Mbp (s4_span_chr1_f0) and 2x4c node packing (s4_jobs_per_node_2x4c_n3):
/usr/bin/time -v vardictrust -G GRCh38.fa -f 0 -N HCC1395 -b HCC1395.bam -th 8 -R chr1:10000000-20000000 > chr1_f0.var
```

### 1.5 Multi-Regime Parameter Sweep: Comprehensive 87-Arm `.var` and `.vcf` Parity Sweep (`sweep_*`)

The `12` multi-regime sweep cells are `sweep_4targets_full`, `sweep_R_chr1_full`, `sweep_hg002_chr20_exome`, `sweep_multimodal_and_cram_twin`, `sweep_rna_hg002_52gb`, and `sweep_span_chr1_f0`. They continue with `sweep_tierp_13modes`, `sweep_tierp_SRR15006540`, `sweep_tierp_SRR8657348`, `sweep_tierp_somatic_76_86`, `sweep_ultradeep_mrd_stress`, and `sweep_wgs_mrd_matrix`. Together they validate `.var` and `.vcf` equivalence across `87` comparative arms (`87/87` `MATCH`, `100.00%`, `Supplementary Table S2`):
```bash
# Ultra-deep 45000x ctDNA MRD stress (-f 0.0005 and -f 0 with --fisher, sweep_ultradeep_mrd_stress):
vardictrust -G GRCh38.fa -f 0.0005 -N SAMPLE -b ultradeep_45000x.bam --fisher -c 1 -S 2 -E 3 -g 4 panel.bed > rust_ultradeep_f0005.var
# WGS & exome MRD matrix (-f 0.001 and -f 0, sweep_wgs_mrd_matrix):
vardictrust -G GRCh38.fa -f 0.001 -N HG005 -b hg005_300x.bam --fisher -c 1 -S 2 -E 3 -g 4 --vcf mrd_targets.bed > rust_hg005_300x_mrd.vcf
```

### 1.6 Domain 5: Full-Genome WGS (Unchunked Full Chromosomes vs. Pre-Chunked Callable Intervals), 10-Platform Chr20 Panel, and Supply Chain (`s5_*`)
Domain 5 comprises `7` cells (`s5_wgs_hg002_full_contigs_3way_th32`, `s5_wgs_hg002_full_giab_3way_th32`, `s5_wgs_hg008_full_contigs_3way_th32`, `s5_wgs_hg008_full_giab_3way_th32`, `s5_grch38_field_callers_truth`, `s5_cross_platform_short_read_chr20_panel`, and `s5_h4_maintenance_security_matrix`):
```bash
# Unchunked full chromosomes (chr1..chr22, s5_wgs_hg002_full_contigs_3way_th32 & s5_wgs_hg008_full_contigs_3way_th32):
timeout 10800 /usr/bin/time -v vardictrust -G GRCh38_GIABv3.fasta -f 0.01 -N HG002 \
  -b HG002_2x250_131GB.bam -th 32 --fisher -z -c 1 -S 2 -E 3 -g 4 --vcf grch38_22autosomes_full.bed > hg002_unchunked.vcf

timeout 10800 /usr/bin/time -v vardictrust -G GRCh38_GIABv3.fasta -f 0.01 -N "HG008-T|HG008-N" \
  -b "HG008-T_195x.bam|HG008-N_169x.bam" -th 32 -U --fisher -z -c 1 -S 2 -E 3 -g 4 --vcf grch38_22autosomes_full.bed > hg008_unchunked.vcf

# Pre-chunked <=50 kbp GIAB callable intervals (s5_wgs_hg002_full_giab_3way_th32 & s5_wgs_hg008_full_giab_3way_th32):
timeout 10800 /usr/bin/time -v vardictrust -G GRCh38_GIABv3.fasta -f 0.01 -N HG002 \
  -b HG002_2x250_131GB.bam -th 32 --fisher -z -c 1 -S 2 -E 3 -g 4 --vcf HG002_GRCh38_v421_50kb.bed > hg002_giab50kb.vcf

timeout 10800 /usr/bin/time -v vardictrust -G GRCh38_GIABv3.fasta -f 0.01 -N "HG008-T|HG008-N" \
  -b "HG008-T_195x.bam|HG008-N_169x.bam" -th 32 -U --fisher -z -c 1 -S 2 -E 3 -g 4 --vcf HG008_GRCh38_v0.3_50kb.bed > hg008_giab50kb.vcf
```

### 1.7 Domain 6: GA4GH v3.5 Stratification, GIAB v5.0q/CMRG/Mosaic, HG008 v0.4 11-Caller Matrix, and Low-VAF Titration (`s6_*`)
Domain 6 comprises `4` cells: `s6_ga4gh_v35_stratification_matrix`, `s6_hg002_v50q_cmrg_mosaic_eval`, `s6_hg008_modern_ml_and_commercial_callers`, and `s6_low_vaf_titration_truth_curve`. The caller-matrix cell includes the official NIST GIAB `Liss_lab` somatic callsets for `DRAGEN v4.2.4`, `DeepSomatic v1.6.0`, and `ClairS v0.4.1`. It also includes `Lancet2 (NYGC)`, `Hartwig SAGE`, and `Ultima DeepVariant Somatic` alongside local callsets:
```bash
# Controlled in silico tumor-purity titration (20%, 10%, 5%, 2%, 1% purity at -f 0.002, s6_low_vaf_titration_truth_curve):
for PURITY in 20 10 5 2 1; do
  vardictrust -G GRCh38_GIABv3.fasta -f 0.002 -N "HG008_P${PURITY}|HG008-N" \
    -b "hg008_titration_p${PURITY}.bam|HG008-N.bam" -th 16 --fisher -z -c 1 -S 2 -E 3 -g 4 \
    --vcf hg008_slice.bed > "rust_titration_p${PURITY}.vcf"
done

# Next-gen GIAB v5.0q T2T-Q100, CMRG v1.00, and Mosaic v1.10 evaluation (s6_hg002_v50q_cmrg_mosaic_eval):
rtg vcfeval -b HG002_GRCh38_v5.0q_T2T.vcf.gz -e HG002_GRCh38_v5.0q_T2T.bed \
  -c hg002_unchunked.vcf.gz -t GRCh38.sdf -o eval_v50q_unchunked
```

---

## 2. Supplementary Note 1: Analytical Derivation of JVM (`VarDictJava`) and `vardictcpp` Heap Retention (Unchunked Chromosomes vs. Pre-Chunked Callable Intervals)

In `VarDictJava` v1.8.3 (`CigarParser.java` and `VariationRealigner.java`), the pileup engine processes each BED region $[S, E]$ of length $L = E - S + 1$. It instantiates two region-wide heap structures before emitting any variant row:

1. A reference base map `Reference.referenceSequences` (`Map<Integer, Character>`) holding $O(L)$ boxed `java.lang.Integer` (`16` bytes) and `java.lang.Character` (`16` bytes) objects plus `HashMap$Node` headers (`32` bytes). That costs $\approx 64\text{–}96$ bytes per reference base ($\approx 16\text{–}24$ GiB for human `chr1` of length $L = 249,250,621$ bp per worker thread).

2. A position-indexed variation map `Map<Integer, VariationMap<String, Variation>> nonInsertionVariants` and `insertionVariants`. Each covered coordinate $p \in [S, E]$ allocates a `VariationMap` and one `Variation` object per observed allele. At sequencing depth $D$, the active heap footprint across $T$ worker threads scales as:
   $$\text{RSS}_{\text{Java / C++} }(L, T) = \Theta\!\left(\sum_{t=1}^{T} L_t \cdot (c_{\text{ref} } + c_{\text{var} } \cdot \bar{a}(D))\right)$$
   where $L_t$ is the region length assigned to thread $t$ and $\bar{a}(D)$ is the mean allele count per locus. When $L_t \le 50$ kbp (pre-chunked intervals), garbage collection reclaims objects between intervals (`34.14` GB on HG002; `62.51` GB with timeout `rc=124` on HG008). When $L_t \in [50, 249]$ Mbp (unchunked chromosomes or `-R chr1:1-249250621`), live objects remain reachable until the chromosome finishes. That drives `VarDictJava` to `106.47` GiB on `chr1` and `>213.7` GB (`rc=3` `OutOfMemoryError`) on multi-threaded full genomes (`Table 2`, `Table 5`). `vardictcpp` similarly retains per-region C++ hash tables across each interval. It consumes `10.58` GiB on `chr1`, `31.91` GB on HG002 unchunked chromosomes, and `91.85` GB (`rc=124` timeout) on HG008 unchunked chromosomes. An external script must pre-split the BED file (`CPP_CHUNK=50,000`) to avoid this.

`VarDictRust` leaves intervals $\le 1,000,000$ bp untiled. Any region $L > 1,000,000$ bp splits into tiles. Default width is $W = 100$ kbp with halo $H = 10$ kbp (`200` kbp in `-i` splice mode; `Figure 1b`). Memory is then bounded by:
$$\text{RSS}_{\text{VarDictRust} }(L, T) = O\!\left(T \cdot (W + 2H)\right) \perp L$$
Increasing $L$ by `250`× on chromosome 1 (`1` Mbp to `249` Mbp, `s4_span_chr1_3way_th8`..`th32`) changes `VarDictRust` RSS by `1.18`×–`1.41`× (`1.38`–`3.69` GiB), allowing unchunked full-genome calling in `4.93` GB (HG002) and `12.14` GB (HG008).

---

## 3. Supplementary Note 2: Ten-Point Technical Audit of `vardictcpp` Divergences (`CPP-1` through `CPP-10`)

Across cells `s2_tierp_mode_parity_sweep`, `s5_h4_maintenance_security_matrix`, `s5_wgs_hg002_full_contigs_3way_th32`, and `s6_low_vaf_titration_truth_curve`, we cataloged ten reproducibility and capability gaps in `vardictcpp ae93fe39dc (Wolff 2026)` (`Table 8`, `Supplementary Table S6`):

- **`CPP-1` (VCF Header & `INFO/END` Schema):** `vardictcpp` omits the `##INFO=<ID=END>` header definition and diverges from `var2vcf_valid.pl` / `var2vcf_paired.pl` VCF 4.3 schema formatting.

- **`CPP-2` (CRISPR (`-J`/`-j`) & Structural Variant Modes):** `vardictcpp` sets `--nosv` (`-U`) by default. It omits CRISPR cleavage-site (`-J`/`-j`) and structural variant (`StructuralVariantsProcessor`) assembly. That failed `1/13` modes in `s2_tierp_mode_parity_sweep` (`12/13` exact).

- **`CPP-3` (Amplicon (`-a`) Output Format):** In amplicon primer-clipping mode (`-a 10:0.95`, `s1_ceiling_amplicon`), `vardictcpp` emits raw tab-separated `.var` text instead of VCF (`0` VCF lines), breaking downstream `bcftools` and `rtg vcfeval` consumers.

- **`CPP-4` (Concatenated CLI Flag Parsing):** Bundled single-character flags (e.g., `-z -c 1` combined or legacy `VarDictJava` bundled options) break `vardictcpp` CLI parsing. It exits with code `120`.

- **`CPP-5` (50 kbp Chunk Boundary Deduplication):** Without external BED splitting (`--chunk 50,000` / `CPP_CHUNK=50,000`), `vardictcpp` consumes `10.58` GiB on `chr1`. It consumes `31.91` GB on HG002 unchunked chromosomes and `91.85` GB (`rc=124` timeout) on HG008 unchunked chromosomes. When intervals are chunked, `vardictcpp` lacks boundary deduplication. It emits `+309,569` duplicate/shifted calls on HG002 unchunked chromosomes (`9,953,338` vs `9,643,769`) and `+394` on pre-chunked intervals (`8,111,565` vs `8,111,171`).

- **`CPP-6` (Missing `##contig` Header Lines):** `vardictcpp` omits `##contig=<ID=...,length=...>` lines from VCF headers (`VarDictJava` issue `#320`). Strict VCF validators and `rtg vcfeval` header checks then warn or fail.

- **`CPP-7` (GRCh37 / b37 Chromosome Naming):** `vardictcpp` hardcodes `chr`-prefix assumptions in internal chromosome handling. That breaks or diverges on `b37` / `GRCh37` numeric contigs (`1`..`22`, `X`, `Y`, `MT`).

- **`CPP-8` (BAM/CRAM Reference & Decode):** `vardictcpp` diverges on CRAM reference decoding (`-G ref.fa`). It requires dynamic linkage against external `libhts`, `libcurl`, `liblzma`, `libbz2`, and `libdeflate`. `VarDictRust` natively decodes BAM and CRAM in a `4.56` MB zero-dependency static `musl` binary (`s4_multimodal_and_cram_twin`).

- **`CPP-9` (MNP Phasing, Tie-Break & Somatic SV Drops):** `vardictcpp` diverges from `VarDictJava` on complex MNP phasing, realignment tie-breaking, and somatic variant retention. At `20%` tumor purity (`s6_low_vaf_titration_truth_curve`), `vardictcpp` drops somatic indel $F_1$ to `0.2466` (recall `15.65%`). SNV $F_1$ drops to `0.6914` (recall `57.07%`). `VarDictRust` matches `VarDictJava` instead (`Indel` $F_1 = 0.3104$, recall `20.75%`; `SNV` $F_1 = 0.7512$, recall `65.59%`).

- **`CPP-10` (RNA-Seq Splice Mode (`-i`.var / `.splice.tsv`):** `vardictcpp` does not implement `-i` exon-junction splice-site accounting (`s3_rna_splice_parity_n3`, `51.8` GiB RNA-seq BAM, `12,668` records).

---

## 4. Supplementary Note 3: Biostatistical Methodology, Exact Clopper-Pearson Intervals, and `--fix` Governance

### 4.1 Exact Clopper-Pearson Binomial Confidence Intervals

For binary oracle-concordance (Tier-1 `.var` and Tier-2 `.vcf`) and completion rates, let $k$ be successes out of $n$ independent cells or arms. We compute exact two-sided $1 - \alpha = 95\%$ Clopper-Pearson binomial CIs $[p_{\text{lo} }, p_{\text{hi} }]$ via the inverse Beta distribution:

$$p_{\text{lo} } = \text{Beta}^{-1}\!\left(\frac{\alpha}{2};\, k,\, n - k + 1\right), \quad p_{\text{hi} } = \text{Beta}^{-1}\!\left(1 - \frac{\alpha}{2};\, k + 1,\, n - k\right)$$

with $p_{\text{lo} } = 0$ when $k = 0$ and $p_{\text{hi} } = 1$ when $k = n$:
- Complete `60`-cell verification campaign (`60/60` exit `0`): `100%` (`95%` CI `94.04%–100.00%`)
- `13`-mode functional parity sweep (`13/13`): `100%` (`95%` CI `75.29%–100.00%`)
- Extended coverage functional modalities (`17/17`): `100%` (`95%` CI `80.49%–100.00%`)
- Multimodal clinical and edge-case slices (`18/18`): `100%` (`95%` CI `81.47%–100.00%`)
- Cross-platform chr20 runs where `VarDictJava` completed (`11/11`): `100%` (`95%` CI `71.51%–100.00%`)
- Cross-platform `-f 0` somatic noise-floor `STRICT-BYTES` sweep (`18/18`): `100%` (`95%` CI `81.47%–100.00%`)
- GA4GH `v3.5` `70`-cell stratification matrix (`70/70`): `100%` (`95%` CI `94.87%–100.00%`)
- Multi-regime parity sweep across `12` suites (`87/87`): `100.00%` (`95%` CI `95.85%–100.00%`). Scope covers untiled $\le 1,000,000$ bp BED suites plus single-window `--tile-size` runs on `sweep_span_chr1_f0` (`7,077` `.var` / `7,032` `.vcf`) and `sweep_R_chr1_full` (`4,705` `.var`; `Supplementary Table S2`)
- NIST Mosaic SNV `v1.10` sensitivity (`84/85`): `98.82%` (`95%` CI `93.62%–99.97%`)

### 4.2 Non-Parametric Runtime and Memory Hypothesis Testing

Across the $N = 12$ AML sample-by-thread cells (`s2_tierp_v2_*`, `n = 3` technical replicates each), we report the mean. Replicate range is min–max (`Supplementary Table S3`). Comparing `VarDictRust` against `VarDictJava` across $N = 12$ gives two-sided Mann-Whitney $U = 4.0$ ($p = 9.735e-05$) for wall time. For peak RSS, $U = 0.0$ ($p = 3.658e-05$). An exact paired sign-flip permutation test gives $p = 2^{-12} = 0.00024$.

### 4.3 Release Governance: `v1.0.x` Parity Versus `v1.1+` Opt-In `--fix` Flags
All `v1.0.x` releases (`v1.0.6`) enforce strict oracle parity on valid inputs. Output-changing upstream fixes (`#92` `VD > DP`, `#62` over-filtering, `#292` `SHIFT3`, `#342` `<dup-N>` symbolic alleles) are gated behind opt-in `--fix <id>` flags starting in `v1.1.0`, while parity-restoring fixes (`splice-tile-ownership`, `mate-window-ownership`) ship as verified defaults.

---

## 5. Supplementary Tables (`Supplementary Table S1` – `Supplementary Table S6`) and Complete 60-Cell Verification Ledger

<a id="supp_table1_60cell_inventory"></a>

### Supplementary Table S1. Complete 60-Experiment Benchmark Provenance Inventory, Execution Status & Manuscript Crosswalk

| Archive Directory ID | Evaluation Suite | Assay Class<sup>a</sup> | Exit Code<sup>b</sup> | Artifacts | Manuscript Targets | Parity / Execution Verdict | Experiment Description |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| s1_ceiling_amplicon | Baseline Equivalence & Tiling | Targeted Panel | 0 | 6 | Figure 1, Figure 4, Table 1, Supplementary Table S1 | 100% Exact Match (2 variants) | EGFR/MRE11a amplicon (-a) 3-way ceiling & exact parity |
| s1_ceiling_hg008 | Baseline Equivalence & Tiling | Whole Genome | 0 | 6 | Figure 1, Figure 4, Table 1, Supplementary Table S1 | 100% Exact Match (18,426 variants) | NIST GIAB HG008 somatic slice 3-way ceiling & exact parity |
| s1_gate7_ref_multiregion | Baseline Equivalence & Tiling | Whole Genome | 0 | 5 | Figure 1, Figure 2, Table 1, Supplementary Table S1 | 100% Exact Match (2,572 variants) | HG002 chr20 multi-region reference boundary parity gate |
| s1_h1a_simple_tiles | Baseline Equivalence & Tiling | Whole Exome | 0 | 5 | Figure 1, Figure 2, Table 1, Supplementary Table S1 | 100% Exact Match across 4 tile sizes (5,090 variants) | Tile-size invariance (0.5, 1, 2, 10 Mbp) on SEQC2 HCC1395 |
| s2_fig1_fig3_field_callers | AML Panel & Field-Caller Replication [11] | Targeted Panel | 0 | 73 | Figure 1, Figure 2, Table 2, Supplementary Table S1 | 16/16 field-caller runs completed (VarDictRust == VarDictJava) | Wolff [11] Fig 1 & Fig 3 field callers (FreeBayes, LoFreq, Mutect2, Octopus, Strelka2) at 1c/8c |
| s2_fig2_hg002_chr20_truth | AML Panel & Field-Caller Replication [11] | Whole Genome | 0 | 20 | Figure 1, Figure 3, Table 2, Table 3, Supplementary Table S1 | 100% Exact Match (7,422 variants; 7 callers evaluated) | Wolff [11] Fig 2 HG002 chr20 exome 7-caller truth accuracy |
| s2_table3_R_chr1_java_200g_n3 | AML Panel & Field-Caller Replication [11] | Targeted Panel | 0 | 6 | Figure 1, Figure 4, Table 2, Table 4, Supplementary Table S1 | VarDictJava -Xmx4g OOM (rc=3); completes at -Xmx200g (n=3) | Wolff [11] Table 3 -R chr1:1-249250621 Java -Xmx4g OOM & -Xmx200g completion (n=3) |
| s2_table3_R_chr1_rust_and_cpp | AML Panel & Field-Caller Replication [11] | Targeted Panel | 0 | 10 | Figure 1, Figure 4, Table 2, Table 4, Supplementary Table S1 | VarDictRust & vardictcpp complete in bounded memory (n=3) | Wolff [11] Table 3 -R chr1:1-249250621 VarDictRust vs vardictcpp (n=3) |
| s2_tierp_mode_parity_sweep | AML Panel & Field-Caller Replication [11] | Targeted Panel | 0 | 14 | Figure 1, Figure 2, Figure 8, Table 1, Table 2, Table 8, Supplementary Table S1, Supplementary Table S6 | 100% Exact Match across all 13 functional CLI modes | 13-mode 3-way functional parity sweep (Simple, Somatic, Amplicon, Splice, SV, CRAM, Fisher, -f 0) |
| s2_tierp_v2_SRR15006375_th1 | AML Panel & Field-Caller Replication [11] | Targeted Panel | 0 | 19 | Figure 1, Figure 2, Table 2, Supplementary Table S1, Supplementary Table S3 | 100% Exact Match (4,056 .var rows; 3,966 VCF variants) | Clinical AML AML panel SRR15006375 3-way parity & timing at 1 thread (n=3) |
| s2_tierp_v2_SRR15006375_th8 | AML Panel & Field-Caller Replication [11] | Targeted Panel | 0 | 19 | Figure 1, Figure 2, Table 2, Supplementary Table S1, Supplementary Table S3 | 100% Exact Match (4,056 .var rows; 3,966 VCF variants) | Clinical AML AML panel SRR15006375 3-way parity & timing at 8 threads (n=3) |
| s2_tierp_v2_SRR15006376_th1 | AML Panel & Field-Caller Replication [11] | Targeted Panel | 0 | 19 | Figure 1, Figure 2, Table 2, Supplementary Table S1, Supplementary Table S3 | 100% Exact Match (12,739 .var rows; 12,655 VCF variants) | Clinical AML AML panel SRR15006376 3-way parity & timing at 1 thread (n=3) |
| s2_tierp_v2_SRR15006376_th8 | AML Panel & Field-Caller Replication [11] | Targeted Panel | 0 | 19 | Figure 1, Figure 2, Table 2, Supplementary Table S1, Supplementary Table S3 | 100% Exact Match (12,739 .var rows; 12,655 VCF variants) | Clinical AML AML panel SRR15006376 3-way parity & timing at 8 threads (n=3) |
| s2_tierp_v2_SRR15006386_th1 | AML Panel & Field-Caller Replication [11] | Targeted Panel | 0 | 19 | Figure 1, Figure 2, Table 2, Supplementary Table S1, Supplementary Table S3 | 100% Exact Match (4,014 .var rows; 3,925 VCF variants) | Clinical AML AML panel SRR15006386 3-way parity & timing at 1 thread (n=3) |
| s2_tierp_v2_SRR15006386_th8 | AML Panel & Field-Caller Replication [11] | Targeted Panel | 0 | 19 | Figure 1, Figure 2, Table 2, Supplementary Table S1, Supplementary Table S3 | 100% Exact Match (4,014 .var rows; 3,925 VCF variants) | Clinical AML AML panel SRR15006386 3-way parity & timing at 8 threads (n=3) |
| s2_tierp_v2_SRR15006540_th1 | AML Panel & Field-Caller Replication [11] | Targeted Panel | 0 | 19 | Figure 1, Figure 2, Table 2, Supplementary Table S1, Supplementary Table S3 | 100% Exact Match (16,031 .var rows; 15,912 VCF variants) | Clinical AML AML panel SRR15006540 3-way parity & timing at 1 thread (n=3) |
| s2_tierp_v2_SRR15006540_th8 | AML Panel & Field-Caller Replication [11] | Targeted Panel | 0 | 19 | Figure 1, Figure 2, Table 2, Supplementary Table S1, Supplementary Table S3 | 100% Exact Match (16,031 .var rows; 15,912 VCF variants) | Clinical AML AML panel SRR15006540 3-way parity & timing at 8 threads (n=3) |
| s2_tierp_v2_SRR8657348_th1 | AML Panel & Field-Caller Replication [11] | Targeted Panel | 0 | 19 | Figure 1, Figure 2, Table 2, Supplementary Table S1, Supplementary Table S3 | 100% Exact Match (54,448 .var rows; 53,224 VCF variants) | Clinical AML AML panel SRR8657348 (773,667 intervals) 3-way parity & timing at 1 thread (n=3) |
| s2_tierp_v2_SRR8657348_th8 | AML Panel & Field-Caller Replication [11] | Targeted Panel | 0 | 19 | Figure 1, Figure 2, Table 2, Supplementary Table S1, Supplementary Table S3 | 100% Exact Match (54,448 .var rows; 53,224 VCF variants) | Clinical AML AML panel SRR8657348 (773,667 intervals) 3-way parity & timing at 8 threads (n=3) |
| s2_tierp_v2_somatic_76_86_th1 | AML Panel & Field-Caller Replication [11] | Targeted Panel | 0 | 17 | Figure 1, Figure 2, Table 2, Supplementary Table S1, Supplementary Table S3 | 100% Exact Match (14,126 .var rows; 14,000 VCF variants) | Clinical AML somatic pair SRR15006376 / SRR15006386 3-way parity & timing at 1 thread (n=3) |
| s2_tierp_v2_somatic_76_86_th8 | AML Panel & Field-Caller Replication [11] | Targeted Panel | 0 | 17 | Figure 1, Figure 2, Table 2, Supplementary Table S1, Supplementary Table S3 | 100% Exact Match (14,126 .var rows; 14,000 VCF variants) | Clinical AML somatic pair SRR15006376 / SRR15006386 3-way parity & timing at 8 threads (n=3) |
| s3_rna_splice_parity_n3 | Transcriptome RNA-Seq Splice Parity | Transcriptome RNA-Seq | 0 | 10 | Figure 1, Figure 4, Table 1, Table 4, Supplementary Table S1 | 100% Exact Match (12,668 rows) | Full-sample 51.8 GiB HG002 GM24385 RNA-seq (-i splice mode) parity & timing (n=3) |
| s4_jobs_per_node_2x4c_n3 | Multi-Thread Scaling, Span & ROC | Whole Exome | 0 | 7 | Figure 1, Figure 4, Table 4, Supplementary Table S1 | 2x4-thread node packing completed (0.859 GiB combined RSS, n=3) | Concurrent 2x4-core node-packing throughput & RSS on SEQC2 HCC1395 exome (n=3) |
| s4_multimodal_and_cram_twin | Multi-Thread Scaling, Span & ROC | Whole Genome | 0 | 3 | Figure 1, Figure 2, Table 1, Supplementary Table S1 | 100% Exact Match (0 textual differences within eps=1e-6) | Lossless BAM-vs-CRAM twin zero-diff check & version verification |
| s4_perf_budget_4targets_n3 | Multi-Thread Scaling, Span & ROC | Panel, Exome & WGS | 0 | 54 | Figure 1, Figure 4, Table 4, Supplementary Table S1 | 3-way timing & RSS verified across all 4 benchmark targets (n=3) | 3-way performance budget across amplicon, HG002 germline, SEQC2 somatic exome, and HG008 somatic (n=3) |
| s4_scale_exome_3way_th1 | Multi-Thread Scaling, Span & ROC | Whole Exome | 0 | 11 | Figure 1, Figure 4, Figure 8, Table 4, Supplementary Table S1, Supplementary Table S3 | 100% Exact Match (339,698 variants) | SEQC2 HCC1395 whole-exome 3-way scaling at 1 thread (n=3) |
| s4_scale_exome_3way_th2 | Multi-Thread Scaling, Span & ROC | Whole Exome | 0 | 12 | Figure 1, Figure 4, Figure 8, Table 4, Supplementary Table S1, Supplementary Table S3 | 100% Exact Match (339,698 variants) | SEQC2 HCC1395 whole-exome 3-way scaling at 2 threads (n=3) |
| s4_scale_exome_3way_th4 | Multi-Thread Scaling, Span & ROC | Whole Exome | 0 | 12 | Figure 1, Figure 4, Figure 8, Table 4, Supplementary Table S1, Supplementary Table S3 | 100% Exact Match (339,698 variants) | SEQC2 HCC1395 whole-exome 3-way scaling at 4 threads (n=3) |
| s4_scale_exome_3way_th8 | Multi-Thread Scaling, Span & ROC | Whole Exome | 0 | 12 | Figure 1, Figure 4, Figure 8, Table 4, Supplementary Table S1, Supplementary Table S3 | 100% Exact Match (339,698 variants) | SEQC2 HCC1395 whole-exome 3-way scaling at 8 threads (n=3) |
| s4_scale_exome_3way_th16 | Multi-Thread Scaling, Span & ROC | Whole Exome | 0 | 12 | Figure 1, Figure 4, Figure 8, Table 4, Supplementary Table S1, Supplementary Table S3 | 100% Exact Match (339,698 variants) | SEQC2 HCC1395 whole-exome 3-way scaling at 16 threads (n=3) |
| s4_scale_exome_3way_th32 | Multi-Thread Scaling, Span & ROC | Whole Exome | 0 | 12 | Figure 1, Figure 4, Figure 8, Table 4, Supplementary Table S1, Supplementary Table S3 | 100% Exact Match (339,698 variants) | SEQC2 HCC1395 whole-exome 3-way scaling at 32 threads (n=3) |
| s4_span_chr1_3way_th8 | Multi-Thread Scaling, Span & ROC | Whole Exome | 0 | 37 | Figure 1, Figure 4, Table 4, Supplementary Table S1, Supplementary Table S3 | Bounded memory verified across 1–249 Mbp spans (<1.50x RSS ratio, n=3) | Chromosome 1 span-to-RSS scaling (1, 10, 100, 249 Mbp) at 8 threads (n=3) |
| s4_span_chr1_3way_th16 | Multi-Thread Scaling, Span & ROC | Whole Exome | 0 | 37 | Figure 1, Figure 4, Table 4, Supplementary Table S1, Supplementary Table S3 | Bounded memory verified across 1–249 Mbp spans (<1.50x RSS ratio, n=3) | Chromosome 1 span-to-RSS scaling (1, 10, 100, 249 Mbp) at 16 threads (n=3) |
| s4_span_chr1_3way_th32 | Multi-Thread Scaling, Span & ROC | Whole Exome | 0 | 37 | Figure 1, Figure 4, Table 4, Supplementary Table S1, Supplementary Table S3 | Bounded memory verified across 1–249 Mbp spans (<1.50x RSS ratio, n=3) | Chromosome 1 span-to-RSS scaling (1, 10, 100, 249 Mbp) at 32 threads (n=3) |
| s4_span_chr1_f0 | Multi-Thread Scaling, Span & ROC | Whole Exome | 0 | 9 | Figure 1, Figure 4, Table 1, Table 4, Supplementary Table S1 | 100% Exact Match (7,032 variants) | Chromosome 1 10 Mbp all-positions (-f 0) parity & memory stress (n=3) |
| s4_tier3_truth_and_qual_sweep | Multi-Thread Scaling, Span & ROC | Exome & WGS | 0 | 24 | Figure 1, Figure 3, Figure 7, Table 3, Table 7, Supplementary Table S1, Supplementary Table S4 | 100% Exact Match across HG002, SEQC2 & HG008 QUAL ROC sweeps | rtg vcfeval truth concordance & QUAL ROC sweep (0..200) on HG002, SEQC2, HG008 |
| evidence_certified | Multi-Thread Scaling, Span & ROC | Panel, Exome & WGS | 0 | 16 | Figure 1, Figure 2, Figure 3, Table 1, Table 3, Supplementary Table S1 | 100% Exact Match (17/17 extended modes, 18/18 slices, HG005 & NA24385) | GCP-certified evidence bundle (17-modality extended coverage, 18-slice multimodal, HG005, NA24385-b37) |
| sweep_4targets_full | Full-Core Multi-Regime Parity Sweep | Panel, Exome & WGS | 0 | 15 | Figure 1, Figure 2, Table 1, Supplementary Table S1, Supplementary Table S2 | 100% Exact Match (4/4 sub-arms) | Full-core multi-regime .var & .vcf parity sweep across amplicon, hg002, seqc2, hg008 |
| sweep_R_chr1_full | Full-Core Multi-Regime Parity Sweep | Targeted Panel | 0 | 6 | Figure 1, Figure 2, Table 1, Supplementary Table S1, Supplementary Table S2 | 100% Exact Match (1/1 sub-arm) | Full-core multi-regime .var parity sweep on -R chr1:1-249250621 (249 Mbp single region) |
| sweep_hg002_chr20_exome | Full-Core Multi-Regime Parity Sweep | Whole Genome | 0 | 9 | Figure 1, Figure 2, Table 1, Supplementary Table S1, Supplementary Table S2 | 100% Exact Match (2/2 sub-arms) | Full-core multi-regime .var & .vcf parity sweep on HG002 chr20 exome |
| sweep_multimodal_and_cram_twin | Full-Core Multi-Regime Parity Sweep | Whole Genome | 0 | 11 | Figure 1, Figure 2, Table 1, Supplementary Table S1, Supplementary Table S2 | 100% Exact Match (4/4 sub-arms) | Full-core multi-regime multimodal slices, extended coverage, and HG008 CRAM twin timing & parity |
| sweep_rna_hg002_52gb | Full-Core Multi-Regime Parity Sweep | Transcriptome RNA-Seq | 0 | 10 | Figure 1, Figure 2, Table 1, Supplementary Table S1, Supplementary Table S2 | 100% Exact Match (2/2 sub-arms) | Full-core multi-regime .var & .splice.tsv parity sweep on 51.8 GiB HG002 RNA-seq |
| sweep_span_chr1_f0 | Full-Core Multi-Regime Parity Sweep | Whole Exome | 0 | 9 | Figure 1, Figure 2, Table 1, Supplementary Table S1, Supplementary Table S2 | 100% Exact Match (2/2 sub-arms) | Full-core multi-regime .var & .vcf parity sweep on chr1 10 Mbp at -f 0 |
| sweep_tierp_13modes | Full-Core Multi-Regime Parity Sweep | Targeted Panel | 0 | 49 | Figure 1, Figure 2, Table 1, Supplementary Table S1, Supplementary Table S2 | 100% Exact Match (12/12 sub-arms) | Full-core multi-regime .var parity sweep across 12 Clinical AML operational mode sub-arms |
| sweep_tierp_SRR15006540 | Full-Core Multi-Regime Parity Sweep | Targeted Panel | 0 | 9 | Figure 1, Figure 2, Table 1, Supplementary Table S1, Supplementary Table S2 | 100% Exact Match (2/2 sub-arms) | Full-core multi-regime .var & .vcf parity sweep on Clinical AML SRR15006540 |
| sweep_tierp_SRR8657348 | Full-Core Multi-Regime Parity Sweep | Targeted Panel | 0 | 9 | Figure 1, Figure 2, Table 1, Supplementary Table S1, Supplementary Table S2 | 100% Exact Match (2/2 sub-arms) | Full-core multi-regime .var & .vcf parity sweep on Clinical AML SRR8657348 |
| sweep_tierp_somatic_76_86 | Full-Core Multi-Regime Parity Sweep | Targeted Panel | 0 | 9 | Figure 1, Figure 2, Table 1, Supplementary Table S1, Supplementary Table S2 | 100% Exact Match (2/2 sub-arms) | Full-core multi-regime .var & .vcf parity sweep on Clinical AML somatic pair SRR15006376 / SRR15006386 |
| sweep_ultradeep_mrd_stress | Full-Core Multi-Regime Parity Sweep | Panel & WGS | 0 | 107 | Figure 1, Figure 2, Table 1, Supplementary Table S1, Supplementary Table S2 | 100% Exact Match (30/30 sub-arms) | Full-core multi-regime ultra-deep MRD stress test (-f 0.0005 and -f 0 with --fisher up to 45,000x) |
| sweep_wgs_mrd_matrix | Full-Core Multi-Regime Parity Sweep | Exome & WGS | 0 | 86 | Figure 1, Figure 2, Table 1, Supplementary Table S1, Supplementary Table S2 | 100% Exact Match (24/24 sub-arms) | Full-core multi-regime WGS & exome MRD matrix (-f 0.001 and -f 0 across HG002, HG005 300x, SEQC2, HG008) |
| s5_cross_platform_short_read_chr20_panel | Full-Genome WGS, 10-Platform & Security | 10-Platform Chr20 WGS | 0 | 420 | Figure 1, Figure 6, Table 6, Supplementary Table S1, Supplementary Table S4 | 100% Exact Match (11/11 completing Chr20 runs & 18/18 noise-floor runs) | 10-platform sequencing chemistry chr20 panel (Part A -f 0.01/-f 0 std/UG-HCR; Part B 18-arm noise floor; Ultima CRAM twin) |
| s5_grch38_field_callers_truth | Full-Genome WGS, 10-Platform & Security | Whole Genome | 0 | 78 | Figure 1, Figure 3, Table 3, Supplementary Table S1 | 19/19 caller-cohort evaluations completed (VarDictRust == VarDictJava) | Native GRCh38 germline (HG002 chr20) & somatic (HG008 slice, SEQC2) field-caller truth evaluation |
| s5_h4_maintenance_security_matrix | Full-Genome WGS, 10-Platform & Security | Cross-Modal Security & Parity | 0 | 13 | Figure 1, Figure 8, Table 8, Supplementary Table S1, Supplementary Table S6 | Static musl equivalence verified & CPP-1..10 cataloged | supply-chain security, static musl binary, CLI flag parity, and CPP-1..10 divergence matrix |
| s5_wgs_hg002_full_contigs_3way_th32 | Full-Genome WGS, 10-Platform & Security | Full-Genome WGS (chr1–chr22) | 0 | 14 | Figure 1, Figure 5, Table 5, Supplementary Table S1 | VarDictJava OOM (rc=3); VarDictRust completes (9,643,769 variants) | Full-genome HG002 (130.8 GB 2x250 bp) (22 unchunked chromosomes) 3-way scaling |
| s5_wgs_hg002_full_giab_3way_th32 | Full-Genome WGS, 10-Platform & Security | Full-Genome WGS (chr1–chr22) | 0 | 27 | Figure 1, Figure 5, Table 1, Table 5, Supplementary Table S1 | 100% Exact Match (8,111,171 variants) | Full-genome HG002 (130.8 GB 2x250 bp) (<=50 kbp GIAB callable intervals) 3-way parity & truth |
| s5_wgs_hg008_full_contigs_3way_th32 | Full-Genome WGS, 10-Platform & Security | Full-Genome WGS (chr1–chr22) | 0 | 12 | Figure 1, Figure 5, Table 5, Supplementary Table S1 | VarDictJava OOM (rc=3) & C++ timeout; VarDictRust completes (13,067,287 variants) | Full-genome HG008 (329.0 GB 364x T/N) (22 unchunked chromosomes) 3-way scaling |
| s5_wgs_hg008_full_giab_3way_th32 | Full-Genome WGS, 10-Platform & Security | Full-Genome WGS (chr1–chr22) | 0 | 20 | Figure 1, Figure 5, Table 5, Supplementary Table S1 | VarDictJava & C++ timeout (rc=124); VarDictRust completes (10,195,118 variants) | Full-genome HG008 (329.0 GB 364x T/N) (<=50 kbp GIAB v0.4 callable intervals) 3-way scaling & v0.3/v0.4 truth |
| s6_ga4gh_v35_stratification_matrix | Next-Gen GIAB, GA4GH v3.5 & Low-VAF Titration | GA4GH v3.5 Stratification [17] | 0 | 2 | Figure 1, Figure 6, Table 6, Supplementary Table S1, Supplementary Table S5 | 100% Exact Match across all 70 GA4GH v3.5 stratification cells | 70-cell GA4GH v3.5 genomic stratification matrix (7 strata x 10 sequencing platforms) |
| s6_hg002_v50q_cmrg_mosaic_eval | Next-Gen GIAB, GA4GH v3.5 & Low-VAF Titration | Next-Gen GIAB Truth [14–16, 20] | 0 | 2 | Figure 1, Figure 7, Table 7, Supplementary Table S1 | 12/12 evaluations completed (100% Match on callable BED) | HG002 GIAB v5.0q T2T-Q100, CMRG v1.00, and Mosaic SNV v1.10 evaluation |
| s6_hg008_modern_ml_and_commercial_callers | Next-Gen GIAB, GA4GH v3.5 & Low-VAF Titration | GIAB HG008 v0.4 Somatic [19] | 0 | 2 | Figure 1, Figure 7, Table 7, Supplementary Table S1 | 27/27 caller-benchmark evaluations completed across 11 callers | HG008 GIAB v0.4 benchmark across VarDictRust, 7 NIST Liss_lab callers (DRAGEN, DeepSomatic, ClairS, Lancet2, SAGE, Ultima DV), Mutect2, Strelka2, LoFreq |
| s6_low_vaf_titration_truth_curve | Next-Gen GIAB, GA4GH v3.5 & Low-VAF Titration | Somatic Low-VAF Titration [19] | 0 | 42 | Figure 1, Figure 3, Table 3, Supplementary Table S1 | 100% Exact Match across all 5 purities (30/30 caller-purity runs) | HG008 somatic tumor-purity titration (20%->1% purity / 10%->0.5% VAF) at -f 0.002 across 6 callers |

*<sup>a</sup>Assay classes span targeted clinical panels and amplicons [11], whole-exome sequencing [10], 51.8 GiB transcriptome RNA-seq [13], 10-platform Chromosome 20 sequencing chemistries [6, 13, 19], and full-genome germline/somatic WGS [13–19, 21–25]. <sup>b</sup>Every experiment completed with process exit code 0 across 60/60 benchmark archive directories (exact Clopper-Pearson 95% CI 94.04%–100.00%) and maps directly to the figures and tables listed in Column 6. Canonical ledger: benchmarks/gcp-study-v106/cells_60_inventory.tsv (with SHA-256 manifests) is the single source of truth; this table is a generated formatted view of that file.*

---

<a id="supp_table2_option_c_87arms"></a>

### Supplementary Table S2. Baseline Equivalence Suites and Complete 87-Arm Full-Core Intermediate (.var) & VCF Parity Sweep

| Evaluation Suite | Sub-Arm / Configuration & Format<sup>a</sup> | VarDictJava [1] Rows / Variants | VarDictRust Rows / Variants | FP vs. Java | FN vs. Java | Parity Verdict<sup>a</sup> |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| Amplicon Baseline | EGFR/MRE11A targeted amplicon calling (-a 10:0.95) (VCF) | 2 | 2 | 0 | 0 | 100% Exact Match |
| Somatic Benchmark Slice | NIST GIAB HG008 [19] tumor-normal somatic slice (VCF) | 18,426 | 18,426 | 0 | 0 | 100% Exact Match |
| Multi-Region Reference Offset | GIAB HG002 [13, 18] Chr20 multi-region coordinate offset verification (VCF) | 2,572 | 2,572 | 0 | 0 | 100% Exact Match |
| Tile-Size Invariance<sup>b</sup> | SEQC2 HCC1395 [10] halo-tiling invariance (0.5, 1, 2, 10 Mbp tiles) (VCF) | 5,090 | 5,090 | 0 | 0 | 100% Exact Match (all 4 tile sizes)<sup>b</sup> |
| Lossless BAM-vs-CRAM Twin | GIAB HG008 [19] lossless BAM-vs-CRAM [26] twin & 18 multimodal slices (VCF) | 1,129 | 1,129 | 0 | 0 | 100% Exact Match (0 byte differences) |
| 4-Assay Regression Suite | Targeted EGFR/MRE11A amplicon (-a) (VCF) | 2 | 2 | 0 | 0 | 100% Exact Match |
| 4-Assay Regression Suite | GIAB HG002 [13, 18] deep WGS slice (VCF) | 227,562 | 227,562 | 0 | 0 | 100% Exact Match |
| 4-Assay Regression Suite | GIAB HG008 [19] tumor-normal somatic slice (VCF) | 18,426 | 18,426 | 0 | 0 | 100% Exact Match |
| 4-Assay Regression Suite | SEQC2 HCC1395 [10] paired somatic exome slice (VCF) | 623,444 | 623,444 | 0 | 0 | 100% Exact Match |
| Contiguous 249 Mbp -R chr1 [11]<sup>c</sup> | Unchunked -R chr1:1-249250621 single window (--tile-size 250000000) (.var)<sup>c</sup> | 4,705 | 4,705 | 0 | 0 | 100% Exact Match<sup>c</sup> |
| GIAB HG002 Chr20 Exome [11] | GIAB HG002 Chr20 exome [11] (.var) | 8,330 | 8,330 | 0 | 0 | 100% Exact Match |
| GIAB HG002 Chr20 Exome [11] | GIAB HG002 Chr20 exome [11] (VCF) | 8,233 | 8,233 | 0 | 0 | 100% Exact Match |
| Multimodal & CRAM Twin Suite [26] | 17-modality extended CLI coverage suite (.var & VCF) | 17 | 17 | 0 | 0 | 100% Exact Match |
| Multimodal & CRAM Twin Suite [26] | 18-slice multimodal regression suite (VCF) | 18 | 18 | 0 | 0 | 100% Exact Match |
| Multimodal & CRAM Twin Suite [26] | GIAB HG008 [19] CRAM-vs-BAM [26] twin slice (VCF) | 1 | 1 | 0 | 0 | 100% Exact Match |
| Multimodal & CRAM Twin Suite [26] | Lossless CRAM-vs-BAM [26] bitwise twin verification (VCF) | 40 | 40 | 0 | 0 | 100% Exact Match |
| 51.8 GiB HG002 RNA-Seq [13] | 51.8 GiB HG002 GM24385 [13] splice junction table (-i) (.splice.tsv) | 12,668 | 12,668 | 0 | 0 | 100% Exact Match |
| 51.8 GiB HG002 RNA-Seq [13] | 51.8 GiB HG002 GM24385 [13] cDNA variant calls (-i) (.var) | 1,230 | 1,230 | 0 | 0 | 100% Exact Match |
| 10 Mbp Chr1 All-Positions (-f 0)<sup>c</sup> | Chromosome 1 10 Mbp at -f 0 (--tile-size 10000000) (.var)<sup>c</sup> | 7,077 | 7,077 | 0 | 0 | 100% Exact Match<sup>c</sup> |
| 10 Mbp Chr1 All-Positions (-f 0)<sup>c</sup> | Chromosome 1 10 Mbp at -f 0 (--tile-size 10000000) (VCF)<sup>c</sup> | 7,032 | 7,032 | 0 | 0 | 100% Exact Match<sup>c</sup> |
| 12-Sub-Arm CLI Mode Matrix [11] | Simple germline calling (-f 0.01) (.var) | 228 | 228 | 0 | 0 | 100% Exact Match |
| 12-Sub-Arm CLI Mode Matrix [11] | Simple zero-threshold MRD calling (-f 0) (.var) | 248 | 248 | 0 | 0 | 100% Exact Match |
| 12-Sub-Arm CLI Mode Matrix [11] | Simple germline calling with native Fisher exact test (--fisher) (.var) | 650 | 650 | 0 | 0 | 100% Exact Match |
| 12-Sub-Arm CLI Mode Matrix [11] | Chromosome 17 targeted clinical panel (.var) | 83 | 83 | 0 | 0 | 100% Exact Match |
| 12-Sub-Arm CLI Mode Matrix [11] | Structural variant deletion regions (.var) | 3 | 3 | 0 | 0 | 100% Exact Match |
| 12-Sub-Arm CLI Mode Matrix [11] | Low-coverage germline boundary calling (.var) | 146 | 146 | 0 | 0 | 100% Exact Match |
| 12-Sub-Arm CLI Mode Matrix [11] | Chromosome 17 amplicon awareness (-a) (.var) | 145 | 145 | 0 | 0 | 100% Exact Match |
| 12-Sub-Arm CLI Mode Matrix [11] | Amplicon calling with native Fisher exact test (-a --fisher) (.var) | 145 | 145 | 0 | 0 | 100% Exact Match |
| 12-Sub-Arm CLI Mode Matrix [11] | Paired tumor-normal somatic calling (-f 0.01) (.var) | 664 | 664 | 0 | 0 | 100% Exact Match |
| 12-Sub-Arm CLI Mode Matrix [11] | Paired somatic zero-threshold MRD calling (-f 0) (.var) | 1,396 | 1,396 | 0 | 0 | 100% Exact Match |
| 12-Sub-Arm CLI Mode Matrix [11] | Paired somatic calling with native Fisher exact test (--fisher) (.var) | 527 | 527 | 0 | 0 | 100% Exact Match |
| 12-Sub-Arm CLI Mode Matrix [11] | CRAM-vs-BAM [26] lossless twin intermediate table (-G) (.var) | 228 | 228 | 0 | 0 | 100% Exact Match |
| AML Panel SRR15006540 [11] | Clinical AML panel SRR15006540 [11] (.var) | 16,031 | 16,031 | 0 | 0 | 100% Exact Match |
| AML Panel SRR15006540 [11] | Clinical AML panel SRR15006540 [11] (VCF) | 15,912 | 15,912 | 0 | 0 | 100% Exact Match |
| AML Panel SRR8657348 [11] | Clinical AML panel SRR8657348 [11] (.var) | 54,448 | 54,448 | 0 | 0 | 100% Exact Match |
| AML Panel SRR8657348 [11] | Clinical AML panel SRR8657348 [11] (VCF) | 53,224 | 53,224 | 0 | 0 | 100% Exact Match |
| AML Somatic Pair (SRR15006376 / SRR15006386) [11] | Clinical AML somatic pair SRR15006376 / SRR15006386 [11] (.var) | 14,126 | 14,126 | 0 | 0 | 100% Exact Match |
| AML Somatic Pair (SRR15006376 / SRR15006386) [11] | Clinical AML somatic pair SRR15006376 / SRR15006386 [11] (VCF) | 14,000 | 14,000 | 0 | 0 | 100% Exact Match |
| Ultra-Deep MRD Stress Suite | SRR15006375 (-f 0.0005 --fisher) (.var) | 2,780 | 2,780 | 0 | 0 | 100% Exact Match |
| Ultra-Deep MRD Stress Suite | SRR15006375 (-f 0.0005 --fisher) (VCF) | 2,645 | 2,645 | 0 | 0 | 100% Exact Match |
| Ultra-Deep MRD Stress Suite | SRR15006375 (-f 0 --fisher) (.var) | 6,182 | 6,182 | 0 | 0 | 100% Exact Match |
| Ultra-Deep MRD Stress Suite | SRR15006375 (-f 0 --fisher) (VCF) | 5,685 | 5,685 | 0 | 0 | 100% Exact Match |
| Ultra-Deep MRD Stress Suite | SRR15006540 (-f 0.0005 --fisher) (.var) | 5,046 | 5,046 | 0 | 0 | 100% Exact Match |
| Ultra-Deep MRD Stress Suite | SRR15006540 (-f 0.0005 --fisher) (VCF) | 4,711 | 4,711 | 0 | 0 | 100% Exact Match |
| Ultra-Deep MRD Stress Suite | SRR15006540 (-f 0 --fisher) (.var) | 7,044 | 7,044 | 0 | 0 | 100% Exact Match |
| Ultra-Deep MRD Stress Suite | SRR15006540 (-f 0 --fisher) (VCF) | 6,313 | 6,313 | 0 | 0 | 100% Exact Match |
| Ultra-Deep MRD Stress Suite | SRR15006386 (-f 0.0005 --fisher) (.var) | 3,028 | 3,028 | 0 | 0 | 100% Exact Match |
| Ultra-Deep MRD Stress Suite | SRR15006386 (-f 0.0005 --fisher) (VCF) | 2,888 | 2,888 | 0 | 0 | 100% Exact Match |
| Ultra-Deep MRD Stress Suite | SRR15006386 (-f 0 --fisher) (.var) | 6,153 | 6,153 | 0 | 0 | 100% Exact Match |
| Ultra-Deep MRD Stress Suite | SRR15006386 (-f 0 --fisher) (VCF) | 5,608 | 5,608 | 0 | 0 | 100% Exact Match |
| Ultra-Deep MRD Stress Suite | SRR15006376 (-f 0.0005 --fisher) (.var) | 2,844 | 2,844 | 0 | 0 | 100% Exact Match |
| Ultra-Deep MRD Stress Suite | SRR15006376 (-f 0.0005 --fisher) (VCF) | 2,732 | 2,732 | 0 | 0 | 100% Exact Match |
| Ultra-Deep MRD Stress Suite | SRR15006376 (-f 0 --fisher) (.var) | 3,841 | 3,841 | 0 | 0 | 100% Exact Match |
| Ultra-Deep MRD Stress Suite | SRR15006376 (-f 0 --fisher) (VCF) | 3,630 | 3,630 | 0 | 0 | 100% Exact Match |
| Ultra-Deep MRD Stress Suite | AML Somatic SRR15006376 / SRR15006386 (29kx depth) (-f 0.0005 --fisher) (.var) | 4,045 | 4,045 | 0 | 0 | 100% Exact Match |
| Ultra-Deep MRD Stress Suite | AML Somatic SRR15006376 / SRR15006386 (29kx depth) (-f 0.0005 --fisher) (VCF) | 3,855 | 3,855 | 0 | 0 | 100% Exact Match |
| Ultra-Deep MRD Stress Suite | AML Somatic SRR15006376 / SRR15006386 (29kx depth) (-f 0 --fisher) (.var) | 6,890 | 6,890 | 0 | 0 | 100% Exact Match |
| Ultra-Deep MRD Stress Suite | AML Somatic SRR15006376 / SRR15006386 (29kx depth) (-f 0 --fisher) (VCF) | 6,309 | 6,309 | 0 | 0 | 100% Exact Match |
| Ultra-Deep MRD Stress Suite | AML Somatic SRR15006375 / SRR15006386 (45kx depth) (-f 0.0005 --fisher) (.var) | 4,481 | 4,481 | 0 | 0 | 100% Exact Match |
| Ultra-Deep MRD Stress Suite | AML Somatic SRR15006375 / SRR15006386 (45kx depth) (-f 0.0005 --fisher) (VCF) | 4,262 | 4,262 | 0 | 0 | 100% Exact Match |
| Ultra-Deep MRD Stress Suite | AML Somatic SRR15006375 / SRR15006386 (45kx depth) (-f 0 --fisher) (.var) | 7,474 | 7,474 | 0 | 0 | 100% Exact Match |
| Ultra-Deep MRD Stress Suite | AML Somatic SRR15006375 / SRR15006386 (45kx depth) (-f 0 --fisher) (VCF) | 6,828 | 6,828 | 0 | 0 | 100% Exact Match |
| Ultra-Deep MRD Stress Suite | GIAB HG008 [19] Somatic (15kx depth) (-f 0.0005 --fisher) (.var) | 19 | 19 | 0 | 0 | 100% Exact Match |
| Ultra-Deep MRD Stress Suite | GIAB HG008 [19] Somatic (15kx depth) (-f 0.0005 --fisher) (VCF) | 16 | 16 | 0 | 0 | 100% Exact Match |
| Ultra-Deep MRD Stress Suite | GIAB HG002 [13] Germline (7kx depth) (-f 0.0005 --fisher) (.var) | 1,008 | 1,008 | 0 | 0 | 100% Exact Match |
| Ultra-Deep MRD Stress Suite | GIAB HG002 [13] Germline (7kx depth) (-f 0.0005 --fisher) (VCF) | 640 | 640 | 0 | 0 | 100% Exact Match |
| Ultra-Deep MRD Stress Suite | HG002 GM24385 [13] RNA-Seq (173kx depth) (-f 0.001 --fisher) (.var) | 139 | 139 | 0 | 0 | 100% Exact Match |
| Ultra-Deep MRD Stress Suite | HG002 GM24385 [13] RNA-Seq (173kx depth) (-f 0.001 --fisher) (VCF) | 120 | 120 | 0 | 0 | 100% Exact Match |
| WGS & WES Low-VAF MRD Matrix | GIAB HG002 [13, 18] WGS MRD (-f 0.001) (.var) | 15,083 | 15,083 | 0 | 0 | 100% Exact Match |
| WGS & WES Low-VAF MRD Matrix | GIAB HG002 [13, 18] WGS MRD (-f 0.001) (VCF) | 14,769 | 14,769 | 0 | 0 | 100% Exact Match |
| WGS & WES Low-VAF MRD Matrix | GIAB HG002 [13, 18] WGS MRD (-f 0) (.var) | 15,083 | 15,083 | 0 | 0 | 100% Exact Match |
| WGS & WES Low-VAF MRD Matrix | GIAB HG002 [13, 18] WGS MRD (-f 0) (VCF) | 14,769 | 14,769 | 0 | 0 | 100% Exact Match |
| WGS & WES Low-VAF MRD Matrix | GIAB HG008 [19] Somatic WGS MRD (-f 0.001) (.var) | 21,379 | 21,379 | 0 | 0 | 100% Exact Match |
| WGS & WES Low-VAF MRD Matrix | GIAB HG008 [19] Somatic WGS MRD (-f 0.001) (VCF) | 20,293 | 20,293 | 0 | 0 | 100% Exact Match |
| WGS & WES Low-VAF MRD Matrix | GIAB HG008 [19] Somatic WGS MRD (-f 0.001 --fisher) (.var) | 21,379 | 21,379 | 0 | 0 | 100% Exact Match |
| WGS & WES Low-VAF MRD Matrix | GIAB HG008 [19] Somatic WGS MRD (-f 0.001 --fisher) (VCF) | 20,293 | 20,293 | 0 | 0 | 100% Exact Match |
| WGS & WES Low-VAF MRD Matrix | GIAB HG008 [19] Somatic WGS MRD (-f 0) (.var) | 21,379 | 21,379 | 0 | 0 | 100% Exact Match |
| WGS & WES Low-VAF MRD Matrix | GIAB HG008 [19] Somatic WGS MRD (-f 0) (VCF) | 20,293 | 20,293 | 0 | 0 | 100% Exact Match |
| WGS & WES Low-VAF MRD Matrix | GIAB HG005 [13] 300x WGS chr3 (-f 0.001) (.var) | 339 | 339 | 0 | 0 | 100% Exact Match |
| WGS & WES Low-VAF MRD Matrix | GIAB HG005 [13] 300x WGS chr3 (-f 0.001) (VCF) | 325 | 325 | 0 | 0 | 100% Exact Match |
| WGS & WES Low-VAF MRD Matrix | GIAB HG005 [13] 300x WGS chr3 (-f 0) (.var) | 339 | 339 | 0 | 0 | 100% Exact Match |
| WGS & WES Low-VAF MRD Matrix | GIAB HG005 [13] 300x WGS chr3 (-f 0) (VCF) | 325 | 325 | 0 | 0 | 100% Exact Match |
| WGS & WES Low-VAF MRD Matrix | GIAB HG005 [13] 300x WGS chr6 (-f 0.001) (.var) | 305 | 305 | 0 | 0 | 100% Exact Match |
| WGS & WES Low-VAF MRD Matrix | GIAB HG005 [13] 300x WGS chr6 (-f 0.001) (VCF) | 299 | 299 | 0 | 0 | 100% Exact Match |
| WGS & WES Low-VAF MRD Matrix | GIAB HG005 [13] 300x WGS chr21 (-f 0.001) (.var) | 286 | 286 | 0 | 0 | 100% Exact Match |
| WGS & WES Low-VAF MRD Matrix | GIAB HG005 [13] 300x WGS chr21 (-f 0.001) (VCF) | 284 | 284 | 0 | 0 | 100% Exact Match |
| WGS & WES Low-VAF MRD Matrix | SEQC2 HCC1395 [10] Somatic Exome MRD (-f 0.001) (.var) | 5,261 | 5,261 | 0 | 0 | 100% Exact Match |
| WGS & WES Low-VAF MRD Matrix | SEQC2 HCC1395 [10] Somatic Exome MRD (-f 0.001) (VCF) | 5,206 | 5,206 | 0 | 0 | 100% Exact Match |
| WGS & WES Low-VAF MRD Matrix | SEQC2 HCC1395 [10] Somatic Exome MRD (-f 0.001 --fisher) (.var) | 5,261 | 5,261 | 0 | 0 | 100% Exact Match |
| WGS & WES Low-VAF MRD Matrix | SEQC2 HCC1395 [10] Somatic Exome MRD (-f 0.001 --fisher) (VCF) | 5,206 | 5,206 | 0 | 0 | 100% Exact Match |
| WGS & WES Low-VAF MRD Matrix | SEQC2 HCC1395 [10] Somatic Exome MRD (-f 0) (.var) | 5,261 | 5,261 | 0 | 0 | 100% Exact Match |
| WGS & WES Low-VAF MRD Matrix | SEQC2 HCC1395 [10] Somatic Exome MRD (-f 0) (VCF) | 5,206 | 5,206 | 0 | 0 | 100% Exact Match |

*<sup>a</sup>All 87/87 (100.00%) comparative sub-arms across all 12 multi-regime sweep suites and all 5 baseline equivalence suites achieve 100% exact intermediate (.var, .splice.tsv) and VCF parity against VarDictJava v1.8.3 [1] (0 false positives and 0 false negatives at eps=1e-6). <sup>b</sup>Halo-tiling invariance confirms identical variant output across 0.5, 1.0, 2.0, and 10.0 Mbp tile sizes on SEQC2 HCC1395 [10]. <sup>c</sup>Standard BED intervals (<=1 Mbp) run untiled by default; the two multi-megabase single-interval stress suites (10 Mbp chr1 at -f 0 and 249.25 Mbp -R chr1) use --tile-size matching the region span so that single-window VarDictRust is compared directly against single-window VarDictJava [1]. Source ledger: sweep arm counts and verdicts are generated from master_gcp_metrics.json; this table is the formatted view.*

---

<a id="supp_table3_n3_technical_replicates"></a>

### Supplementary Table S3. Complete Technical Replicate (n=3) and Full-WGS Wall-Clock and Peak RSS Ledger

| Benchmark Suite | Workload / Cohort | Threads | Caller / Engine | Replicate Wall Times (s)<sup>a</sup> | Wall Mean (Min–Max) (s)<sup>a</sup> | Replicate Peak RSS (KB)<sup>a</sup> | Peak RSS Mean (MB / GiB) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| AML Panel Replication [11] | AML Panel SRR15006375 | 1 | VarDictRust | 121.14, 122.26, 122.38 | 121.93 (121.14–122.38) | 215,108, 215,236, 215,736 | 210.3 MB (0.205 GiB) |
| AML Panel Replication [11] | AML Panel SRR15006375 | 1 | VarDictJava [1] | 404.89, 418.94, 400.74 | 408.19 (400.74–418.94) | 2,590,208, 2,591,732, 2,612,224 | 2537.2 MB (2.478 GiB) |
| AML Panel Replication [11] | AML Panel SRR15006375 | 1 | vardictcpp [11] | 109.90, 109.88, 109.19 | 109.66 (109.19–109.90) | 65,292, 65,280, 65,312 | 63.8 MB (0.062 GiB) |
| AML Panel Replication [11] | AML Panel SRR15006375 | 8 | VarDictRust | 49.39, 49.59, 49.44 | 49.47 (49.39–49.59) | 278,052, 274,168, 277,424 | 270.1 MB (0.264 GiB) |
| AML Panel Replication [11] | AML Panel SRR15006375 | 8 | VarDictJava [1] | 229.14, 235.27, 233.69 | 232.70 (229.14–235.27) | 4,665,344, 4,511,724, 4,371,936 | 4410.5 MB (4.307 GiB) |
| AML Panel Replication [11] | AML Panel SRR15006375 | 8 | vardictcpp [11] | 15.42, 15.32, 15.30 | 15.35 (15.30–15.42) | 342,664, 334,764, 331,092 | 328.3 MB (0.321 GiB) |
| AML Panel Replication [11] | AML Panel SRR15006376 | 1 | VarDictRust | 80.84, 80.32, 80.08 | 80.41 (80.08–80.84) | 231,180, 231,156, 230,976 | 225.7 MB (0.220 GiB) |
| AML Panel Replication [11] | AML Panel SRR15006376 | 1 | VarDictJava [1] | 326.10, 340.26, 336.12 | 334.16 (326.10–340.26) | 4,428,032, 4,031,016, 4,579,804 | 4244.4 MB (4.145 GiB) |
| AML Panel Replication [11] | AML Panel SRR15006376 | 1 | vardictcpp [11] | 83.78, 84.32, 83.03 | 83.71 (83.03–84.32) | 75,804, 75,808, 75,804 | 74.0 MB (0.072 GiB) |
| AML Panel Replication [11] | AML Panel SRR15006376 | 8 | VarDictRust | 48.86, 48.70, 48.65 | 48.74 (48.65–48.86) | 287,200, 288,108, 287,872 | 281.0 MB (0.274 GiB) |
| AML Panel Replication [11] | AML Panel SRR15006376 | 8 | VarDictJava [1] | 174.29, 150.76, 168.04 | 164.36 (150.76–174.29) | 7,744,512, 8,044,032, 7,073,792 | 7442.2 MB (7.268 GiB) |
| AML Panel Replication [11] | AML Panel SRR15006376 | 8 | vardictcpp [11] | 13.09, 12.50, 12.96 | 12.85 (12.50–13.09) | 351,772, 364,060, 371,688 | 354.0 MB (0.346 GiB) |
| AML Panel Replication [11] | AML Panel SRR15006386 | 1 | VarDictRust | 114.04, 113.35, 113.17 | 113.52 (113.17–114.04) | 191,748, 191,976, 192,300 | 187.5 MB (0.183 GiB) |
| AML Panel Replication [11] | AML Panel SRR15006386 | 1 | VarDictJava [1] | 334.16, 338.80, 331.87 | 334.94 (331.87–338.80) | 2,596,316, 2,423,020, 2,425,344 | 2423.4 MB (2.367 GiB) |
| AML Panel Replication [11] | AML Panel SRR15006386 | 1 | vardictcpp [11] | 98.37, 98.14, 98.09 | 98.20 (98.09–98.37) | 57,516, 57,516, 57,520 | 56.2 MB (0.055 GiB) |
| AML Panel Replication [11] | AML Panel SRR15006386 | 8 | VarDictRust | 43.81, 43.88, 43.64 | 43.78 (43.64–43.88) | 312,272, 312,860, 313,136 | 305.4 MB (0.298 GiB) |
| AML Panel Replication [11] | AML Panel SRR15006386 | 8 | VarDictJava [1] | 178.99, 182.17, 180.70 | 180.62 (178.99–182.17) | 4,535,524, 4,536,296, 4,551,680 | 4434.7 MB (4.331 GiB) |
| AML Panel Replication [11] | AML Panel SRR15006386 | 8 | vardictcpp [11] | 16.69, 15.01, 14.95 | 15.55 (14.95–16.69) | 300,084, 309,872, 307,332 | 298.6 MB (0.292 GiB) |
| AML Panel Replication [11] | AML Panel SRR15006540 | 1 | VarDictRust | 125.49, 125.61, 126.91 | 126.00 (125.49–126.91) | 199,412, 199,396, 199,072 | 194.6 MB (0.190 GiB) |
| AML Panel Replication [11] | AML Panel SRR15006540 | 1 | VarDictJava [1] | 406.72, 412.42, 402.18 | 407.11 (402.18–412.42) | 4,630,248, 4,671,184, 4,745,216 | 4572.5 MB (4.465 GiB) |
| AML Panel Replication [11] | AML Panel SRR15006540 | 1 | vardictcpp [11] | 111.55, 111.32, 111.66 | 111.51 (111.32–111.66) | 78,084, 78,088, 78,084 | 76.3 MB (0.074 GiB) |
| AML Panel Replication [11] | AML Panel SRR15006540 | 8 | VarDictRust | 77.89, 77.99, 78.24 | 78.04 (77.89–78.24) | 233,572, 227,592, 229,364 | 224.8 MB (0.220 GiB) |
| AML Panel Replication [11] | AML Panel SRR15006540 | 8 | VarDictJava [1] | 179.54, 186.40, 187.21 | 184.38 (179.54–187.21) | 3,863,552, 4,506,172, 4,187,860 | 4087.8 MB (3.992 GiB) |
| AML Panel Replication [11] | AML Panel SRR15006540 | 8 | vardictcpp [11] | 16.77, 16.31, 17.00 | 16.69 (16.31–17.00) | 342,792, 350,980, 347,656 | 339.0 MB (0.331 GiB) |
| AML Panel Replication [11] | AML Panel SRR8657348 | 1 | VarDictRust | 155.05, 155.12, 154.77 | 154.98 (154.77–155.12) | 486,840, 486,476, 483,624 | 474.3 MB (0.463 GiB) |
| AML Panel Replication [11] | AML Panel SRR8657348 | 1 | VarDictJava [1] | 862.51, 857.96, 863.07 | 861.18 (857.96–863.07) | 2,079,488, 2,078,976, 2,073,724 | 2028.7 MB (1.981 GiB) |
| AML Panel Replication [11] | AML Panel SRR8657348 | 1 | vardictcpp [11] | 175.43, 173.98, 174.02 | 174.48 (173.98–175.43) | 146,848, 146,932, 146,936 | 143.5 MB (0.140 GiB) |
| AML Panel Replication [11] | AML Panel SRR8657348 | 8 | VarDictRust | 39.76, 39.62, 39.50 | 39.63 (39.50–39.76) | 492,940, 485,740, 536,324 | 493.2 MB (0.482 GiB) |
| AML Panel Replication [11] | AML Panel SRR8657348 | 8 | VarDictJava [1] | 212.70, 217.72, 217.12 | 215.85 (212.70–217.72) | 3,713,444, 3,823,360, 3,822,592 | 3697.7 MB (3.611 GiB) |
| AML Panel Replication [11] | AML Panel SRR8657348 | 8 | vardictcpp [11] | 24.71, 24.68, 24.64 | 24.68 (24.64–24.71) | 255,472, 259,568, 258,544 | 251.8 MB (0.246 GiB) |
| AML Panel Replication [11] | AML Panel Somatic Pair (SRR15006376 / SRR15006386) | 1 | VarDictRust | 230.96, 230.70, 229.92 | 230.53 (229.92–230.96) | 255,464, 252,304, 253,316 | 247.8 MB (0.242 GiB) |
| AML Panel Replication [11] | AML Panel Somatic Pair (SRR15006376 / SRR15006386) | 1 | VarDictJava [1] | 720.47, 728.53, 729.27 | 726.09 (720.47–729.27) | 4,415,396, 4,396,032, 4,420,820 | 4307.4 MB (4.206 GiB) |
| AML Panel Replication [11] | AML Panel Somatic Pair (SRR15006376 / SRR15006386) | 1 | vardictcpp [11] | 215.50, 215.04, 214.94 | 215.16 (214.94–215.50) | 148,252, 148,252, 148,504 | 144.9 MB (0.141 GiB) |
| AML Panel Replication [11] | AML Panel Somatic Pair (SRR15006376 / SRR15006386) | 8 | VarDictRust | 134.96, 133.83, 137.51 | 135.43 (133.83–137.51) | 383,724, 396,596, 393,488 | 382.1 MB (0.373 GiB) |
| AML Panel Replication [11] | AML Panel Somatic Pair (SRR15006376 / SRR15006386) | 8 | VarDictJava [1] | 373.79, 371.55, 374.68 | 373.34 (371.55–374.68) | 5,850,600, 5,814,784, 5,783,296 | 5679.9 MB (5.547 GiB) |
| AML Panel Replication [11] | AML Panel Somatic Pair (SRR15006376 / SRR15006386) | 8 | vardictcpp [11] | 31.60, 31.64, 31.75 | 31.66 (31.60–31.75) | 824,788, 803,560, 786,332 | 786.0 MB (0.768 GiB) |
| Unchunked Chromosome 1 (-R chr1) [11] | SRR8657348 Unchunked -R chr1 (249.25 Mbp) | 8 | VarDictJava [1] (-Xmx200g) | 509.65, 492.93, 502.90 | 501.83 (492.93–509.65) | 112,714,496, 111,591,168, 110,633,216 | 109029.6 MB (106.474 GiB) |
| Unchunked Chromosome 1 (-R chr1) [11] | SRR8657348 Unchunked -R chr1 (249.25 Mbp) | 8 | VarDictRust | 9.80, 9.81, 9.64 | 9.75 (9.64–9.81) | 320,192, 319,576, 339,012 | 318.6 MB (0.311 GiB) |
| Unchunked Chromosome 1 (-R chr1) [11] | SRR8657348 Unchunked -R chr1 (249.25 Mbp) | 8 | vardictcpp [11] | 106.58, 105.07, 106.03 | 105.89 (105.07–106.58) | 11,098,040, 11,098,040, 11,098,040 | 10837.9 MB (10.584 GiB) |
| Somatic Field Callers [11]<sup>b</sup> | AML Somatic Pair (SRR15006376 / SRR15006386) | 8 | VarDictRust | 135.43 | 135.43 (single run)<sup>b</sup> | 391,269 | 382.1 MB (0.373 GiB) |
| Somatic Field Callers [11]<sup>b</sup> | AML Somatic Pair (SRR15006376 / SRR15006386) | 8 | VarDictJava [1] | 373.34 | 373.34 (single run)<sup>b</sup> | 5,816,227 | 5679.9 MB (5.547 GiB) |
| Somatic Field Callers [11]<sup>b</sup> | AML Somatic Pair (SRR15006376 / SRR15006386) | 8 | vardictcpp [11] | 31.66 | 31.66 (single run)<sup>b</sup> | 804,893 | 786.0 MB (0.768 GiB) |
| Somatic Field Callers [11]<sup>b</sup> | AML Somatic Pair (SRR15006376 / SRR15006386) | 8 | Strelka2 [4] | 199.96 | 199.96 (single run)<sup>b</sup> | 1,184,452 | 1156.7 MB (1.130 GiB) |
| Somatic Field Callers [11]<sup>b</sup> | AML Somatic Pair (SRR15006376 / SRR15006386) | 8 | LoFreq [7] | 521.82 | 521.82 (single run)<sup>b</sup> | 2,146,108 | 2095.8 MB (2.047 GiB) |
| Somatic Field Callers [11]<sup>b</sup> | AML Somatic Pair (SRR15006376 / SRR15006386) | 8 | Octopus [5] | 450.87 | 450.87 (single run)<sup>b</sup> | 2,521,004 | 2461.9 MB (2.404 GiB) |
| Somatic Field Callers [11]<sup>b</sup> | AML Somatic Pair (SRR15006376 / SRR15006386) | 8 | Mutect2 [8] | 969.01 | 969.01 (single run)<sup>b</sup> | 2,355,632 | 2300.4 MB (2.247 GiB) |
| Transcriptome RNA-Seq Splice (-i) | 51.8 GiB GIAB HG002 GM24385 [13] RNA-Seq (-i) | 8 | VarDictRust | 5.01, 4.90, 4.90 | 4.94 (4.90–5.01) | 361,360, 368,040, 376,072 | 359.9 MB (0.351 GiB) |
| Transcriptome RNA-Seq Splice (-i) | 51.8 GiB GIAB HG002 GM24385 [13] RNA-Seq (-i) | 8 | VarDictJava [1] | 23.21, 23.03, 23.70 | 23.31 (23.03–23.70) | 5,324,544, 4,538,880, 5,217,280 | 4909.1 MB (4.794 GiB) |
| 4-Target Performance Budget | Targeted EGFR/MRE11A Amplicon Panel | 8 | VarDictRust | 0.01, 0.01, 0.01 | 0.01 (0.01–0.01) | 11,776, 11,776, 12,032 | 11.6 MB (0.011 GiB) |
| 4-Target Performance Budget | Targeted EGFR/MRE11A Amplicon Panel | 8 | VarDictJava [1] | 3.06, 3.12, 3.07 | 3.08 (3.06–3.12) | 122,744, 122,968, 122,964 | 120.0 MB (0.117 GiB) |
| 4-Target Performance Budget | Targeted EGFR/MRE11A Amplicon Panel | 8 | vardictcpp [11] | 5.37, 5.22, 5.24 | 5.28 (5.22–5.37) | 122,600, 122,576, 122,572 | 119.7 MB (0.117 GiB) |
| 4-Target Performance Budget | GIAB HG002 [13, 18] Germline Benchmark Slice | 8 | VarDictRust | 21.30, 21.33, 21.29 | 21.31 (21.29–21.33) | 486,236, 519,872, 487,080 | 486.1 MB (0.475 GiB) |
| 4-Target Performance Budget | GIAB HG002 [13, 18] Germline Benchmark Slice | 8 | VarDictJava [1] | 159.59, 152.13, 153.25 | 154.99 (152.13–159.59) | 10,989,804, 11,341,096, 10,761,680 | 10772.3 MB (10.520 GiB) |
| 4-Target Performance Budget | GIAB HG002 [13, 18] Germline Benchmark Slice | 8 | vardictcpp [11] | 21.14, 20.69, 20.10 | 20.64 (20.10–21.14) | 879,700, 880,956, 881,460 | 860.1 MB (0.840 GiB) |
| 4-Target Performance Budget | GIAB HG008 [19] Somatic Benchmark Slice | 8 | VarDictRust | 17.14, 17.19, 17.11 | 17.15 (17.11–17.19) | 522,756, 518,000, 528,268 | 510.8 MB (0.499 GiB) |
| 4-Target Performance Budget | GIAB HG008 [19] Somatic Benchmark Slice | 8 | VarDictJava [1] | 79.15, 78.54, 77.28 | 78.32 (77.28–79.15) | 17,866,968, 18,674,616, 17,734,080 | 17667.9 MB (17.254 GiB) |
| 4-Target Performance Budget | GIAB HG008 [19] Somatic Benchmark Slice | 8 | vardictcpp [11] | 17.30, 17.17, 16.91 | 17.13 (16.91–17.30) | 578,696, 584,836, 582,644 | 568.4 MB (0.555 GiB) |
| 4-Target Performance Budget | SEQC2 HCC1395 [10] Paired Somatic Exome Slice | 8 | VarDictRust | 279.29, 277.68, 277.52 | 278.16 (277.52–279.29) | 526,640, 526,944, 529,060 | 515.2 MB (0.503 GiB) |
| 4-Target Performance Budget | SEQC2 HCC1395 [10] Paired Somatic Exome Slice | 8 | VarDictJava [1] | 3967.00, 4030.00, 4045.00 | 4014.00 (3967.00–4045.00) | 14,618,112, 14,880,508, 16,471,012 | 14964.1 MB (14.613 GiB) |
| 4-Target Performance Budget | SEQC2 HCC1395 [10] Paired Somatic Exome Slice | 8 | vardictcpp [11] | 293.71, 294.67, 287.70 | 292.03 (287.70–294.67) | 2,314,000, 2,274,252, 2,299,560 | 2242.1 MB (2.190 GiB) |
| Whole-Exome Thread Scaling | SEQC2 HCC1395 [10] Single-Sample Whole Exome | 1 | VarDictRust | 805.74, 806.69, 825.81 | 812.75 (805.74–825.81) | 351,724, 351,844, 351,652 | 343.5 MB (0.335 GiB) |
| Whole-Exome Thread Scaling | SEQC2 HCC1395 [10] Single-Sample Whole Exome | 1 | VarDictJava [1] | 5767.00, 5789.00, 5717.00 | 5757.67 (5717.00–5789.00) | 1,857,212, 1,887,008, 1,903,816 | 1838.5 MB (1.795 GiB) |
| Whole-Exome Thread Scaling | SEQC2 HCC1395 [10] Single-Sample Whole Exome | 1 | vardictcpp [11] | 993.97, 1010.80, 1012.37 | 1005.71 (993.97–1012.37) | 1,264,172, 1,263,912, 1,264,168 | 1234.5 MB (1.206 GiB) |
| Whole-Exome Thread Scaling | SEQC2 HCC1395 [10] Single-Sample Whole Exome | 2 | VarDictRust | 405.60, 411.13, 410.70 | 409.14 (405.60–411.13) | 428,676, 428,772, 417,680 | 415.1 MB (0.405 GiB) |
| Whole-Exome Thread Scaling | SEQC2 HCC1395 [10] Single-Sample Whole Exome | 2 | VarDictJava [1] | 5762.00, 5839.00, 5806.00 | 5802.33 (5762.00–5839.00) | 5,351,168, 4,222,464, 5,321,168 | 4848.6 MB (4.735 GiB) |
| Whole-Exome Thread Scaling | SEQC2 HCC1395 [10] Single-Sample Whole Exome | 2 | vardictcpp [11] | 529.51, 530.65, 531.93 | 530.70 (529.51–531.93) | 1,357,128, 1,357,144, 1,357,020 | 1325.3 MB (1.294 GiB) |
| Whole-Exome Thread Scaling | SEQC2 HCC1395 [10] Single-Sample Whole Exome | 4 | VarDictRust | 211.15, 212.32, 213.44 | 212.30 (211.15–213.44) | 462,724, 454,340, 444,064 | 443.1 MB (0.433 GiB) |
| Whole-Exome Thread Scaling | SEQC2 HCC1395 [10] Single-Sample Whole Exome | 4 | VarDictJava [1] | 2254.41, 2276.68, 2280.43 | 2270.51 (2254.41–2280.43) | 6,932,480, 7,551,728, 6,877,428 | 6953.7 MB (6.791 GiB) |
| Whole-Exome Thread Scaling | SEQC2 HCC1395 [10] Single-Sample Whole Exome | 4 | vardictcpp [11] | 267.51, 267.11, 267.34 | 267.32 (267.11–267.51) | 1,438,412, 1,434,796, 1,437,096 | 1403.1 MB (1.370 GiB) |
| Whole-Exome Thread Scaling | SEQC2 HCC1395 [10] Single-Sample Whole Exome | 8 | VarDictRust | 117.36, 119.47, 119.88 | 118.90 (117.36–119.88) | 546,844, 553,576, 548,892 | 536.9 MB (0.524 GiB) |
| Whole-Exome Thread Scaling | SEQC2 HCC1395 [10] Single-Sample Whole Exome | 8 | VarDictJava [1] | 1752.96, 1911.57, 1779.37 | 1814.63 (1752.96–1911.57) | 9,542,144, 9,317,612, 9,601,380 | 9264.7 MB (9.048 GiB) |
| Whole-Exome Thread Scaling | SEQC2 HCC1395 [10] Single-Sample Whole Exome | 8 | vardictcpp [11] | 138.32, 138.71, 140.40 | 139.14 (138.32–140.40) | 1,571,092, 1,571,564, 1,557,576 | 1530.0 MB (1.494 GiB) |
| Whole-Exome Thread Scaling | SEQC2 HCC1395 [10] Single-Sample Whole Exome | 16 | VarDictRust | 75.89, 76.76, 74.66 | 75.77 (74.66–76.76) | 613,520, 618,784, 633,916 | 607.5 MB (0.593 GiB) |
| Whole-Exome Thread Scaling | SEQC2 HCC1395 [10] Single-Sample Whole Exome | 16 | VarDictJava [1] | 1885.00, 1929.42, 1866.57 | 1893.66 (1866.57–1929.42) | 18,610,904, 17,986,784, 18,267,580 | 17859.8 MB (17.441 GiB) |
| Whole-Exome Thread Scaling | SEQC2 HCC1395 [10] Single-Sample Whole Exome | 16 | vardictcpp [11] | 73.64, 73.74, 75.12 | 74.17 (73.64–75.12) | 1,760,284, 1,774,788, 1,777,884 | 1729.5 MB (1.689 GiB) |
| Whole-Exome Thread Scaling | SEQC2 HCC1395 [10] Single-Sample Whole Exome | 32 | VarDictRust | 67.27, 68.20, 68.10 | 67.86 (67.27–68.20) | 763,156, 761,796, 770,404 | 747.2 MB (0.730 GiB) |
| Whole-Exome Thread Scaling | SEQC2 HCC1395 [10] Single-Sample Whole Exome | 32 | VarDictJava [1] | 1969.99, 2001.39, 2003.54 | 1991.64 (1969.99–2003.54) | 23,221,636, 23,172,288, 22,081,264 | 22290.1 MB (21.768 GiB) |
| Whole-Exome Thread Scaling | SEQC2 HCC1395 [10] Single-Sample Whole Exome | 32 | vardictcpp [11] | 63.41, 64.15, 63.49 | 63.68 (63.41–64.15) | 2,151,852, 2,164,356, 2,164,628 | 2109.7 MB (2.060 GiB) |
| Chromosome 1 Contiguous Span Scaling | Chromosome 1 Contiguous Span (1 Mbp) | 8 | VarDictRust | 1.16, 1.17, 1.20 | 1.18 (1.16–1.20) | 532,652, 531,732, 531,972 | 519.6 MB (0.507 GiB) |
| Chromosome 1 Contiguous Span Scaling | Chromosome 1 Contiguous Span (1 Mbp) | 8 | vardictcpp [11] | 5.50, 5.45, 5.48 | 5.48 (5.45–5.50) | 195,336, 195,592, 195,080 | 190.8 MB (0.186 GiB) |
| Chromosome 1 Contiguous Span Scaling | Chromosome 1 Contiguous Span (10 Mbp) | 8 | VarDictRust | 1.60, 1.60, 1.60 | 1.60 (1.60–1.60) | 1,243,864, 1,184,828, 1,264,216 | 1202.1 MB (1.174 GiB) |
| Chromosome 1 Contiguous Span Scaling | Chromosome 1 Contiguous Span (10 Mbp) | 8 | vardictcpp [11] | 6.73, 6.67, 6.66 | 6.69 (6.66–6.73) | 384,668, 386,380, 384,580 | 376.2 MB (0.367 GiB) |
| Chromosome 1 Contiguous Span Scaling | Chromosome 1 Contiguous Span (100 Mbp) | 8 | VarDictRust | 12.20, 11.98, 11.80 | 11.99 (11.80–12.20) | 1,277,644, 1,245,396, 1,311,332 | 1248.2 MB (1.219 GiB) |
| Chromosome 1 Contiguous Span Scaling | Chromosome 1 Contiguous Span (100 Mbp) | 8 | vardictcpp [11] | 16.45, 16.32, 16.36 | 16.38 (16.32–16.45) | 489,044, 456,176, 494,436 | 468.6 MB (0.458 GiB) |
| Chromosome 1 Contiguous Span Scaling | Chromosome 1 Contiguous Span (249.25 Mbp) | 8 | VarDictRust | 27.18, 27.21, 27.77 | 27.39 (27.18–27.77) | 1,455,976, 1,458,324, 1,428,280 | 1413.6 MB (1.380 GiB) |
| Chromosome 1 Contiguous Span Scaling | Chromosome 1 Contiguous Span (249.25 Mbp) | 8 | vardictcpp [11] | 30.39, 30.33, 30.27 | 30.33 (30.27–30.39) | 857,536, 814,916, 821,808 | 811.9 MB (0.793 GiB) |
| Chromosome 1 Contiguous Span Scaling | Chromosome 1 Contiguous Span (1 Mbp) | 16 | VarDictRust | 1.14, 1.12, 1.12 | 1.13 (1.12–1.14) | 532,756, 533,136, 532,956 | 520.5 MB (0.508 GiB) |
| Chromosome 1 Contiguous Span Scaling | Chromosome 1 Contiguous Span (1 Mbp) | 16 | vardictcpp [11] | 5.13, 5.09, 5.16 | 5.13 (5.09–5.16) | 354,976, 352,696, 354,120 | 345.6 MB (0.338 GiB) |
| Chromosome 1 Contiguous Span Scaling | Chromosome 1 Contiguous Span (10 Mbp) | 16 | VarDictRust | 0.94, 0.94, 0.96 | 0.95 (0.94–0.96) | 1,799,560, 1,909,712, 1,909,552 | 1829.0 MB (1.786 GiB) |
| Chromosome 1 Contiguous Span Scaling | Chromosome 1 Contiguous Span (10 Mbp) | 16 | vardictcpp [11] | 5.50, 5.76, 5.83 | 5.70 (5.50–5.83) | 690,112, 690,496, 694,340 | 675.4 MB (0.660 GiB) |
| Chromosome 1 Contiguous Span Scaling | Chromosome 1 Contiguous Span (100 Mbp) | 16 | VarDictRust | 6.23, 6.23, 6.38 | 6.28 (6.23–6.38) | 2,036,552, 2,025,764, 2,006,076 | 1975.4 MB (1.929 GiB) |
| Chromosome 1 Contiguous Span Scaling | Chromosome 1 Contiguous Span (100 Mbp) | 16 | vardictcpp [11] | 10.75, 10.72, 10.56 | 10.68 (10.56–10.75) | 836,088, 828,196, 814,352 | 806.9 MB (0.788 GiB) |
| Chromosome 1 Contiguous Span Scaling | Chromosome 1 Contiguous Span (249.25 Mbp) | 16 | VarDictRust | 14.60, 15.14, 14.81 | 14.85 (14.60–15.14) | 2,439,756, 2,382,796, 2,433,268 | 2361.9 MB (2.307 GiB) |
| Chromosome 1 Contiguous Span Scaling | Chromosome 1 Contiguous Span (249.25 Mbp) | 16 | vardictcpp [11] | 17.70, 18.57, 18.90 | 18.39 (17.70–18.90) | 1,209,508, 1,203,464, 1,169,436 | 1166.2 MB (1.139 GiB) |
| Chromosome 1 Contiguous Span Scaling | Chromosome 1 Contiguous Span (1 Mbp) | 32 | VarDictRust | 1.17, 1.17, 1.16 | 1.17 (1.16–1.17) | 535,504, 533,984, 535,260 | 522.4 MB (0.510 GiB) |
| Chromosome 1 Contiguous Span Scaling | Chromosome 1 Contiguous Span (1 Mbp) | 32 | vardictcpp [11] | 5.39, 5.42, 5.55 | 5.45 (5.39–5.55) | 399,284, 397,120, 400,828 | 389.7 MB (0.381 GiB) |
| Chromosome 1 Contiguous Span Scaling | Chromosome 1 Contiguous Span (10 Mbp) | 32 | VarDictRust | 0.91, 0.92, 0.91 | 0.91 (0.91–0.92) | 2,788,972, 2,733,968, 2,726,516 | 2685.4 MB (2.622 GiB) |
| Chromosome 1 Contiguous Span Scaling | Chromosome 1 Contiguous Span (10 Mbp) | 32 | vardictcpp [11] | 6.04, 6.10, 6.16 | 6.10 (6.04–6.16) | 1,262,680, 1,249,556, 1,240,044 | 1221.5 MB (1.193 GiB) |
| Chromosome 1 Contiguous Span Scaling | Chromosome 1 Contiguous Span (100 Mbp) | 32 | VarDictRust | 5.34, 5.38, 5.16 | 5.29 (5.16–5.38) | 2,946,828, 3,016,304, 3,043,928 | 2932.0 MB (2.863 GiB) |
| Chromosome 1 Contiguous Span Scaling | Chromosome 1 Contiguous Span (100 Mbp) | 32 | vardictcpp [11] | 10.79, 10.30, 10.00 | 10.36 (10.00–10.79) | 1,443,032, 1,493,356, 1,510,604 | 1447.6 MB (1.414 GiB) |
| Chromosome 1 Contiguous Span Scaling | Chromosome 1 Contiguous Span (249.25 Mbp) | 32 | VarDictRust | 11.90, 11.96, 11.76 | 11.87 (11.76–11.96) | 3,828,452, 3,939,028, 3,823,184 | 3773.0 MB (3.685 GiB) |
| Chromosome 1 Contiguous Span Scaling | Chromosome 1 Contiguous Span (249.25 Mbp) | 32 | vardictcpp [11] | 16.17, 16.25, 16.26 | 16.23 (16.17–16.26) | 1,821,828, 1,762,144, 1,756,708 | 1738.5 MB (1.698 GiB) |
| All-Positions Pileup (-f 0) | Chromosome 1 10 Mbp All-Positions (-f 0) | 8 | VarDictRust | 20.09, 19.41, 19.24 | 19.58 (19.24–20.09) | 6,185,864, 6,185,476, 6,185,552 | 6040.6 MB (5.899 GiB) |
| All-Positions Pileup (-f 0) | Chromosome 1 10 Mbp All-Positions (-f 0) | 8 | vardictcpp [11] | 6.04, 5.99, 5.99 | 6.01 (5.99–6.04) | 386,612, 380,596, 386,420 | 375.5 MB (0.367 GiB) |
| Concurrent 2x4-Thread Node Packing | SEQC2 [10] Whole Exome (Concurrent Job 1) | 4 | VarDictRust | 212.56, 209.64, 211.94 | 211.38 (209.64–212.56) | 449,244, 446,504, 448,288 | 437.5 MB (0.427 GiB) |
| Concurrent 2x4-Thread Node Packing | SEQC2 [10] Whole Exome (Concurrent Job 2) | 4 | VarDictRust | 212.43, 209.86, 210.85 | 211.05 (209.86–212.43) | 453,764, 451,428, 452,540 | 442.0 MB (0.432 GiB) |
| Full-Genome WGS Scaling (chr1–chr22)<sup>c</sup> | GIAB HG002 [13, 18] Full WGS (22 Unchunked Chromosomes) | 32 | VarDictRust | 1516.61 | 25:16.61 (1,516.61 s)<sup>c</sup> | 4,928,384 | 4.93 GB (4.70 GiB) |
| Full-Genome WGS Scaling (chr1–chr22)<sup>c</sup> | GIAB HG002 [13, 18] Full WGS (<=50 kbp GIAB Callable BED) | 32 | VarDictRust | 1221.21 | 20:21.21 (1,221.21 s)<sup>c</sup> | 1,433,128 | 1.43 GB (1.37 GiB) |
| Full-Genome WGS Scaling (chr1–chr22)<sup>c</sup> | GIAB HG008 [19] 364x Full WGS (22 Unchunked Chromosomes, -U) | 32 | VarDictRust | 10192.00 | 2:49:52 (10,192.00 s)<sup>c</sup> | 12,135,724 | 12.14 GB (11.57 GiB) |
| Full-Genome WGS Scaling (chr1–chr22)<sup>c</sup> | GIAB HG008 [19] 364x Full WGS (<=50 kbp GIAB Callable BED, -U) | 32 | VarDictRust | 6355.00 | 1:45:55 (6,355.00 s)<sup>c</sup> | 4,877,136 | 4.88 GB (4.65 GiB) |

*<sup>a</sup>Unless otherwise noted, rows report all n=3 technical replicate measurements along with the replicate mean (min–max) on GCP n2-standard-32 (Intel Xeon Platinum 8373C Ice Lake) instances, measured via /usr/bin/time -v. <sup>b</sup>The 7-caller somatic comparison on SRR15006376 / SRR15006386 [11] profiles single-run 8-thread executions across VarDictRust, VarDictJava [1], vardictcpp [11], Strelka2 [4], LoFreq [7], Octopus [5], and Mutect2 [8]. <sup>c</sup>Full-genome 32-thread WGS runs on GIAB HG002 (130.8 GB) [13, 18] and GIAB HG008 (329.0 GB, 364x combined depth) [19] report single-run wall-clock times as MM:SS.ss or H:MM:SS (with total seconds in parentheses). Source ledger: replicate measurements aggregate the stage5 and stage6 TSVs via master_gcp_metrics.json; this table is the formatted view.*

---

<a id="supp_table4_qual_roc_and_noise_floor"></a>

### Supplementary Table S4. Verified Truth-Set Operating Points, Variant-Class Stratification, 18-Arm Somatic Noise-Floor Matrix & Low-VAF Titration Detail

| Benchmark Suite | Cohort / Platform | Stratum / Regime / Operating Point | Caller / Parity Status | Primary Metric 1 (TP / Candidates / SNV) | Primary Metric 2 (FP / PASS / InDel) | Summary Score / F1 [12] |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| QUAL ROC Operating Points<sup>a</sup> | GIAB HG002 GRCh38 [13, 18] | All-PASS (Unfiltered, -f 0.01) | VarDictRust == VarDictJava [1] | TP=23,062 | FP=22,458, FN=813 | P=0.5038, R=0.9659, F1=0.6622 |
| QUAL ROC Operating Points<sup>a</sup> | GIAB HG002 GRCh38 [13, 18] | Optimal QUAL (>= 128.0) | VarDictRust == VarDictJava [1] | TP=23,004 | FP=553, FN=871 | P=0.9763, R=0.9635, F1=0.9699 |
| QUAL ROC Operating Points<sup>a</sup> | SEQC2 HCC1395 [10] | All-PASS (Unfiltered, -f 0.01) | VarDictRust == VarDictJava [1] | TP=1,069 | FP=11,334, FN=140 | P=0.0852, R=0.8842, F1=0.1555 |
| QUAL ROC Operating Points<sup>a</sup> | SEQC2 HCC1395 [10] | Optimal QUAL (>= 126.0) | VarDictRust == VarDictJava [1] | TP=794 | FP=293, FN=415 | P=0.7279, R=0.6567, F1=0.6905 |
| QUAL ROC Operating Points<sup>a</sup> | GIAB HG008 v0.3 [19] | All-PASS (Unfiltered, -f 0.01) | VarDictRust == VarDictJava [1] | TP=619 | FP=1,565, FN=288 | P=0.2687, R=0.6825, F1=0.3856 |
| QUAL ROC Operating Points<sup>a</sup> | GIAB HG008 v0.3 [19] | Optimal QUAL (>= 176.0) | VarDictRust == VarDictJava [1] | TP=593 | FP=38, FN=314 | P=0.9354, R=0.6538, F1=0.7696 |
| QUAL ROC Operating Points<sup>a</sup> | GIAB HG005 GRCh38 [13] | All-PASS (Unfiltered, -f 0.01) | VarDictRust == VarDictJava [1] | TP=29 | FP=2, FN=1 | P=0.9355, R=0.9667, F1=0.9508 |
| QUAL ROC Operating Points<sup>a</sup> | GIAB HG005 GRCh38 [13] | Optimal QUAL (>= 226.0) | VarDictRust == VarDictJava [1] | TP=29 | FP=1, FN=1 | P=0.9667, R=0.9667, F1=0.9667 |
| QUAL ROC Operating Points<sup>a</sup> | GIAB NA24385 GRCh37 [13] | All-PASS (Unfiltered, -f 0.01) | VarDictRust == VarDictJava [1] | TP=5 | FP=81, FN=0 | P=0.0581, R=1.0000, F1=0.1099 |
| QUAL ROC Operating Points<sup>a</sup> | GIAB NA24385 GRCh37 [13] | Optimal QUAL (>= 160.0) | VarDictRust == VarDictJava [1] | TP=5 | FP=0, FN=0 | P=1.0000, R=1.0000, F1=1.0000 |
| Variant-Class Stratification<sup>b</sup> | GIAB HG002 GRCh38 [13, 18] | Variant-Class & Allele-Frequency Breakdown | VarDictRust | SNV=209,030, Ins=2,305, Del=10,975 | MNP=2,065, SV=3,187 | PASS=48,447 / Total=227,562 |
| Variant-Class Stratification<sup>b</sup> | GIAB HG002 GRCh38 [13, 18] | Variant-Class & Allele-Frequency Breakdown | VarDictJava [1] | SNV=209,030, Ins=2,305, Del=10,975 | MNP=2,065, SV=3,187 | PASS=48,447 / Total=227,562 |
| Variant-Class Stratification<sup>b</sup> | GIAB HG002 GRCh38 [13, 18] | Variant-Class & Allele-Frequency Breakdown | vardictcpp [11] | SNV=209,066, Ins=2,309, Del=10,949 | MNP=2,049, SV=3,187 | PASS=48,431 / Total=227,560 |
| Variant-Class Stratification<sup>b</sup> | SEQC2 HCC1395 [10] | Variant-Class & Allele-Frequency Breakdown | VarDictRust | SNV=594,344, Ins=8,966, Del=15,716 | MNP=4,383, SV=35 | PASS=120,307 / Total=623,444 |
| Variant-Class Stratification<sup>b</sup> | SEQC2 HCC1395 [10] | Variant-Class & Allele-Frequency Breakdown | VarDictJava [1] | SNV=594,344, Ins=8,966, Del=15,716 | MNP=4,383, SV=35 | PASS=120,307 / Total=623,444 |
| Variant-Class Stratification<sup>b</sup> | SEQC2 HCC1395 [10] | Variant-Class & Allele-Frequency Breakdown | vardictcpp [11] | SNV=594,344, Ins=8,966, Del=15,715 | MNP=4,383, SV=34 | PASS=120,305 / Total=623,442 |
| Variant-Class Stratification<sup>b</sup> | GIAB HG008 v0.3 [19] | Variant-Class & Allele-Frequency Breakdown | VarDictRust | SNV=15,213, Ins=1,161, Del=1,410 | MNP=614, SV=28 | PASS=11,358 / Total=18,426 |
| Variant-Class Stratification<sup>b</sup> | GIAB HG008 v0.3 [19] | Variant-Class & Allele-Frequency Breakdown | VarDictJava [1] | SNV=15,213, Ins=1,161, Del=1,410 | MNP=614, SV=28 | PASS=11,358 / Total=18,426 |
| Variant-Class Stratification<sup>b</sup> | GIAB HG008 v0.3 [19] | Variant-Class & Allele-Frequency Breakdown | vardictcpp [11] | SNV=15,214, Ins=1,161, Del=1,408 | MNP=614, SV=28 | PASS=11,356 / Total=18,425 |
| 18-Arm Somatic Noise Floor<sup>c</sup> | Element AVITI 111x | Min VAF 1%, >=2 Reads (-f 0.01 -r 2) | 100% Exact Match (VarDictRust == VarDictJava [1]) | Total Candidates=61 | PASS Somatic=13 | 0 differences vs. VarDictJava (eps=1e-6)<sup>c</sup> |
| 18-Arm Somatic Noise Floor<sup>c</sup> | Element AVITI 111x | All VAFs, >=2 Reads (-f 0 -r 2) | 100% Exact Match (VarDictRust == VarDictJava [1]) | Total Candidates=61 | PASS Somatic=13 | 0 differences vs. VarDictJava (eps=1e-6)<sup>c</sup> |
| 18-Arm Somatic Noise Floor<sup>c</sup> | Element AVITI 111x | Singleton Noise Floor (-f 0 -r 1) | 100% Exact Match (VarDictRust == VarDictJava [1]) | Total Candidates=186 | PASS Somatic=13 | 0 differences vs. VarDictJava (eps=1e-6)<sup>c</sup> |
| 18-Arm Somatic Noise Floor<sup>c</sup> | PacBio Onso 136x | Min VAF 1%, >=2 Reads (-f 0.01 -r 2) | 100% Exact Match (VarDictRust == VarDictJava [1]) | Total Candidates=76 | PASS Somatic=19 | 0 differences vs. VarDictJava (eps=1e-6)<sup>c</sup> |
| 18-Arm Somatic Noise Floor<sup>c</sup> | PacBio Onso 136x | All VAFs, >=2 Reads (-f 0 -r 2) | 100% Exact Match (VarDictRust == VarDictJava [1]) | Total Candidates=78 | PASS Somatic=19 | 0 differences vs. VarDictJava (eps=1e-6)<sup>c</sup> |
| 18-Arm Somatic Noise Floor<sup>c</sup> | PacBio Onso 136x | Singleton Noise Floor (-f 0 -r 1) | 100% Exact Match (VarDictRust == VarDictJava [1]) | Total Candidates=185 | PASS Somatic=19 | 0 differences vs. VarDictJava (eps=1e-6)<sup>c</sup> |
| 18-Arm Somatic Noise Floor<sup>c</sup> | Illumina 195x | Min VAF 1%, >=2 Reads (-f 0.01 -r 2) | 100% Exact Match (VarDictRust == VarDictJava [1]) | Total Candidates=115 | PASS Somatic=24 | 0 differences vs. VarDictJava (eps=1e-6)<sup>c</sup> |
| 18-Arm Somatic Noise Floor<sup>c</sup> | Illumina 195x | All VAFs, >=2 Reads (-f 0 -r 2) | 100% Exact Match (VarDictRust == VarDictJava [1]) | Total Candidates=131 | PASS Somatic=24 | 0 differences vs. VarDictJava (eps=1e-6)<sup>c</sup> |
| 18-Arm Somatic Noise Floor<sup>c</sup> | Illumina 195x | Singleton Noise Floor (-f 0 -r 1) | 100% Exact Match (VarDictRust == VarDictJava [1]) | Total Candidates=1,199 | PASS Somatic=23 | 0 differences vs. VarDictJava (eps=1e-6)<sup>c</sup> |
| 18-Arm Somatic Noise Floor<sup>c</sup> | Ultima ppmSeq 157x | Min VAF 1%, >=2 Reads (-f 0.01 -r 2) | 100% Exact Match (VarDictRust == VarDictJava [1]) | Total Candidates=574 | PASS Somatic=133 | 0 differences vs. VarDictJava (eps=1e-6)<sup>c</sup> |
| 18-Arm Somatic Noise Floor<sup>c</sup> | Ultima ppmSeq 157x | All VAFs, >=2 Reads (-f 0 -r 2) | 100% Exact Match (VarDictRust == VarDictJava [1]) | Total Candidates=581 | PASS Somatic=136 | 0 differences vs. VarDictJava (eps=1e-6)<sup>c</sup> |
| 18-Arm Somatic Noise Floor<sup>c</sup> | Ultima ppmSeq 157x | Singleton Noise Floor (-f 0 -r 1) | 100% Exact Match (VarDictRust == VarDictJava [1]) | Total Candidates=2,791 | PASS Somatic=159 | 0 differences vs. VarDictJava (eps=1e-6)<sup>c</sup> |
| 18-Arm Somatic Noise Floor<sup>c</sup> | Roche SBX 119x | Min VAF 1%, >=2 Reads (-f 0.01 -r 2) | 100% Exact Match (VarDictRust == VarDictJava [1]) | Total Candidates=770 | PASS Somatic=293 | 0 differences vs. VarDictJava (eps=1e-6)<sup>c</sup> |
| 18-Arm Somatic Noise Floor<sup>c</sup> | Roche SBX 119x | All VAFs, >=2 Reads (-f 0 -r 2) | 100% Exact Match (VarDictRust == VarDictJava [1]) | Total Candidates=788 | PASS Somatic=295 | 0 differences vs. VarDictJava (eps=1e-6)<sup>c</sup> |
| 18-Arm Somatic Noise Floor<sup>c</sup> | Roche SBX 119x | Singleton Noise Floor (-f 0 -r 1) | 100% Exact Match (VarDictRust == VarDictJava [1]) | Total Candidates=3,354 | PASS Somatic=296 | 0 differences vs. VarDictJava (eps=1e-6)<sup>c</sup> |
| 18-Arm Somatic Noise Floor<sup>c</sup> | Ultima UG100 306x | Min VAF 1%, >=2 Reads (-f 0.01 -r 2) | 100% Exact Match (VarDictRust == VarDictJava [1]) | Total Candidates=944 | PASS Somatic=376 | 0 differences vs. VarDictJava (eps=1e-6)<sup>c</sup> |
| 18-Arm Somatic Noise Floor<sup>c</sup> | Ultima UG100 306x | All VAFs, >=2 Reads (-f 0 -r 2) | 100% Exact Match (VarDictRust == VarDictJava [1]) | Total Candidates=1,519 | PASS Somatic=472 | 0 differences vs. VarDictJava (eps=1e-6)<sup>c</sup> |
| 18-Arm Somatic Noise Floor<sup>c</sup> | Ultima UG100 306x | Singleton Noise Floor (-f 0 -r 1) | 100% Exact Match (VarDictRust == VarDictJava [1]) | Total Candidates=4,994 | PASS Somatic=538 | 0 differences vs. VarDictJava (eps=1e-6)<sup>c</sup> |
| Low-VAF Tumor-Purity Titration (Detail)<sup>d</sup> | GIAB HG008 20% Purity (10.0% VAF) [19] | -f 0.002 Somatic Calling Regime | VarDictRust; 100% Exact Match (VarDictRust == VarDictJava) | 0:08.72 (8.72 s) / 762.3 MB | 0.8790 / 0.6559 / 0.7512 | 0.6162 / 0.2075 / 0.3104 |
| Low-VAF Tumor-Purity Titration (Detail)<sup>d</sup> | GIAB HG008 20% Purity (10.0% VAF) [19] | -f 0.002 Somatic Calling Regime | VarDictJava [1]; 100% Exact Match (VarDictRust == VarDictJava) | 1:06.76 (66.76 s) / 12027.7 MB | 0.8790 / 0.6559 / 0.7512 | 0.6162 / 0.2075 / 0.3104 |
| Low-VAF Tumor-Purity Titration (Detail)<sup>d</sup> | GIAB HG008 20% Purity (10.0% VAF) [19] | -f 0.002 Somatic Calling Regime | vardictcpp [11]; MNP / Tie-Break Divergence vs. Java<sup>e</sup> | 0:10.94 (10.94 s) / 1048.3 MB | 0.8770 / 0.5707 / 0.6914 | 0.5823 / 0.1565 / 0.2466 |
| Low-VAF Tumor-Purity Titration (Detail)<sup>d</sup> | GIAB HG008 20% Purity (10.0% VAF) [19] | -f 0.002 Somatic Calling Regime | Strelka2 [4]; External Somatic Caller | 0:36.62 (36.62 s) / 118.8 MB | 1.0000 / 0.7070 / 0.8283 | 1.0000 / 0.1156 / 0.2073 |
| Low-VAF Tumor-Purity Titration (Detail)<sup>d</sup> | GIAB HG008 20% Purity (10.0% VAF) [19] | -f 0.002 Somatic Calling Regime | LoFreq [7]; External Somatic Caller | 1:18.14 (78.14 s) / 252.2 MB | 1.0000 / 0.4072 / 0.5787 | 0.0000 / 0.0000 / 0.0000 |
| Low-VAF Tumor-Purity Titration (Detail)<sup>d</sup> | GIAB HG008 20% Purity (10.0% VAF) [19] | -f 0.002 Somatic Calling Regime | Mutect2 [8]; External Somatic Caller | 7:10.48 (430.48 s) / 1130.0 MB | 1.0000 / 0.5826 / 0.7363 | 0.9032 / 0.0952 / 0.1723 |
| Low-VAF Tumor-Purity Titration (Detail)<sup>d</sup> | GIAB HG008 10% Purity (5.0% VAF) [19] | -f 0.002 Somatic Calling Regime | VarDictRust; 100% Exact Match (VarDictRust == VarDictJava) | 0:08.76 (8.76 s) / 839.4 MB | 0.6726 / 0.5179 / 0.5852 | 0.4082 / 0.1361 / 0.2041 |
| Low-VAF Tumor-Purity Titration (Detail)<sup>d</sup> | GIAB HG008 10% Purity (5.0% VAF) [19] | -f 0.002 Somatic Calling Regime | VarDictJava [1]; 100% Exact Match (VarDictRust == VarDictJava) | 1:06.58 (66.58 s) / 13019.7 MB | 0.6726 / 0.5179 / 0.5852 | 0.4082 / 0.1361 / 0.2041 |
| Low-VAF Tumor-Purity Titration (Detail)<sup>d</sup> | GIAB HG008 10% Purity (5.0% VAF) [19] | -f 0.002 Somatic Calling Regime | vardictcpp [11]; VCF Header / Non-Benchmark Indel Divergence<sup>e</sup> | 0:10.91 (10.91 s) / 1020.3 MB | 0.6726 / 0.5179 / 0.5852 | 0.4082 / 0.1361 / 0.2041 |
| Low-VAF Tumor-Purity Titration (Detail)<sup>d</sup> | GIAB HG008 10% Purity (5.0% VAF) [19] | -f 0.002 Somatic Calling Regime | Strelka2 [4]; External Somatic Caller | 0:36.61 (36.61 s) / 118.9 MB | 1.0000 / 0.4719 / 0.6412 | 1.0000 / 0.0476 / 0.0909 |
| Low-VAF Tumor-Purity Titration (Detail)<sup>d</sup> | GIAB HG008 10% Purity (5.0% VAF) [19] | -f 0.002 Somatic Calling Regime | LoFreq [7]; External Somatic Caller | 1:14.64 (74.64 s) / 251.8 MB | 1.0000 / 0.0716 / 0.1335 | 0.0000 / 0.0000 / 0.0000 |
| Low-VAF Tumor-Purity Titration (Detail)<sup>d</sup> | GIAB HG008 10% Purity (5.0% VAF) [19] | -f 0.002 Somatic Calling Regime | Mutect2 [8]; External Somatic Caller | 7:16.35 (436.35 s) / 1171.3 MB | 1.0000 / 0.3237 / 0.4891 | 0.9375 / 0.0510 / 0.0968 |
| Low-VAF Tumor-Purity Titration (Detail)<sup>d</sup> | GIAB HG008 5% Purity (2.5% VAF) [19] | -f 0.002 Somatic Calling Regime | VarDictRust; 100% Exact Match (VarDictRust == VarDictJava) | 0:08.65 (8.65 s) / 822.7 MB | 0.4932 / 0.3680 / 0.4215 | 0.2183 / 0.1054 / 0.1422 |
| Low-VAF Tumor-Purity Titration (Detail)<sup>d</sup> | GIAB HG008 5% Purity (2.5% VAF) [19] | -f 0.002 Somatic Calling Regime | VarDictJava [1]; 100% Exact Match (VarDictRust == VarDictJava) | 1:06.17 (66.17 s) / 12117.8 MB | 0.4932 / 0.3680 / 0.4215 | 0.2183 / 0.1054 / 0.1422 |
| Low-VAF Tumor-Purity Titration (Detail)<sup>d</sup> | GIAB HG008 5% Purity (2.5% VAF) [19] | -f 0.002 Somatic Calling Regime | vardictcpp [11]; VCF Header / Non-Benchmark Indel Divergence<sup>e</sup> | 0:10.90 (10.90 s) / 1077.7 MB | 0.4932 / 0.3680 / 0.4215 | 0.2183 / 0.1054 / 0.1422 |
| Low-VAF Tumor-Purity Titration (Detail)<sup>d</sup> | GIAB HG008 5% Purity (2.5% VAF) [19] | -f 0.002 Somatic Calling Regime | Strelka2 [4]; External Somatic Caller | 0:36.56 (36.56 s) / 118.8 MB | 1.0000 / 0.1891 / 0.3181 | 1.0000 / 0.0272 / 0.0530 |
| Low-VAF Tumor-Purity Titration (Detail)<sup>d</sup> | GIAB HG008 5% Purity (2.5% VAF) [19] | -f 0.002 Somatic Calling Regime | LoFreq [7]; External Somatic Caller | 1:14.67 (74.67 s) / 251.5 MB | 1.0000 / 0.0051 / 0.0102 | 0.0000 / 0.0000 / 0.0000 |
| Low-VAF Tumor-Purity Titration (Detail)<sup>d</sup> | GIAB HG008 5% Purity (2.5% VAF) [19] | -f 0.002 Somatic Calling Regime | Mutect2 [8]; External Somatic Caller | 7:14.72 (434.72 s) / 1131.1 MB | 1.0000 / 0.0988 / 0.1798 | 0.8889 / 0.0272 / 0.0528 |
| Low-VAF Tumor-Purity Titration (Detail)<sup>d</sup> | GIAB HG008 2% Purity (1.0% VAF) [19] | -f 0.002 Somatic Calling Regime | VarDictRust; 100% Exact Match (VarDictRust == VarDictJava) | 0:08.60 (8.60 s) / 813.6 MB | 0.2430 / 0.1175 / 0.1584 | 0.1391 / 0.0544 / 0.0782 |
| Low-VAF Tumor-Purity Titration (Detail)<sup>d</sup> | GIAB HG008 2% Purity (1.0% VAF) [19] | -f 0.002 Somatic Calling Regime | VarDictJava [1]; 100% Exact Match (VarDictRust == VarDictJava) | 1:05.16 (65.16 s) / 13443.3 MB | 0.2430 / 0.1175 / 0.1584 | 0.1391 / 0.0544 / 0.0782 |
| Low-VAF Tumor-Purity Titration (Detail)<sup>d</sup> | GIAB HG008 2% Purity (1.0% VAF) [19] | -f 0.002 Somatic Calling Regime | vardictcpp [11]; VCF Header / Non-Benchmark Indel Divergence<sup>e</sup> | 0:11.00 (11.00 s) / 1088.3 MB | 0.2430 / 0.1175 / 0.1584 | 0.1391 / 0.0544 / 0.0782 |
| Low-VAF Tumor-Purity Titration (Detail)<sup>d</sup> | GIAB HG008 2% Purity (1.0% VAF) [19] | -f 0.002 Somatic Calling Regime | Strelka2 [4]; External Somatic Caller | 0:36.70 (36.70 s) / 118.8 MB | 1.0000 / 0.0187 / 0.0368 | 0.0000 / 0.0000 / 0.0000 |
| Low-VAF Tumor-Purity Titration (Detail)<sup>d</sup> | GIAB HG008 2% Purity (1.0% VAF) [19] | -f 0.002 Somatic Calling Regime | LoFreq [7]; External Somatic Caller | 1:14.80 (74.80 s) / 251.5 MB | 0.0000 / 0.0000 / 0.0000 | 0.0000 / 0.0000 / 0.0000 |
| Low-VAF Tumor-Purity Titration (Detail)<sup>d</sup> | GIAB HG008 2% Purity (1.0% VAF) [19] | -f 0.002 Somatic Calling Regime | Mutect2 [8]; External Somatic Caller | 6:59.88 (419.88 s) / 1178.6 MB | 1.0000 / 0.0119 / 0.0236 | 0.0000 / 0.0000 / 0.0000 |
| Low-VAF Tumor-Purity Titration (Detail)<sup>d</sup> | GIAB HG008 1% Purity (0.5% VAF) [19] | -f 0.002 Somatic Calling Regime | VarDictRust; 100% Exact Match (VarDictRust == VarDictJava) | 0:08.53 (8.53 s) / 803.0 MB | 0.0856 / 0.0426 / 0.0569 | 0.0660 / 0.0238 / 0.0350 |
| Low-VAF Tumor-Purity Titration (Detail)<sup>d</sup> | GIAB HG008 1% Purity (0.5% VAF) [19] | -f 0.002 Somatic Calling Regime | VarDictJava [1]; 100% Exact Match (VarDictRust == VarDictJava) | 1:05.41 (65.41 s) / 12946.3 MB | 0.0856 / 0.0426 / 0.0569 | 0.0660 / 0.0238 / 0.0350 |
| Low-VAF Tumor-Purity Titration (Detail)<sup>d</sup> | GIAB HG008 1% Purity (0.5% VAF) [19] | -f 0.002 Somatic Calling Regime | vardictcpp [11]; VCF Header / Non-Benchmark Indel Divergence<sup>e</sup> | 0:10.95 (10.95 s) / 1020.9 MB | 0.0856 / 0.0426 / 0.0569 | 0.0660 / 0.0238 / 0.0350 |
| Low-VAF Tumor-Purity Titration (Detail)<sup>d</sup> | GIAB HG008 1% Purity (0.5% VAF) [19] | -f 0.002 Somatic Calling Regime | Strelka2 [4]; External Somatic Caller | 0:36.58 (36.58 s) / 119.0 MB | 1.0000 / 0.0017 / 0.0034 | 0.0000 / 0.0000 / 0.0000 |
| Low-VAF Tumor-Purity Titration (Detail)<sup>d</sup> | GIAB HG008 1% Purity (0.5% VAF) [19] | -f 0.002 Somatic Calling Regime | LoFreq [7]; External Somatic Caller | 1:14.51 (74.51 s) / 252.0 MB | 0.0000 / 0.0000 / 0.0000 | 0.0000 / 0.0000 / 0.0000 |
| Low-VAF Tumor-Purity Titration (Detail)<sup>d</sup> | GIAB HG008 1% Purity (0.5% VAF) [19] | -f 0.002 Somatic Calling Regime | Mutect2 [8]; External Somatic Caller | 6:47.62 (407.62 s) / 1180.5 MB | 1.0000 / 0.0017 / 0.0034 | 0.0000 / 0.0000 / 0.0000 |
| Chr20 Exome 7-Caller Benchmark | GIAB HG002 Chr20 Exome [11, 13] | 7-Caller Truth Benchmark | VarDictRust == VarDictJava [1] | Calls=7,422 | SNV F1=0.9858 | InDel F1=0.8956 (see Table 3) |
| Native GRCh38 Caller Matrix | GIAB HG002 [13], SEQC2 [10] & HG008 [19] GRCh38 | 19-Run Native GRCh38 Caller Matrix | VarDictRust == VarDictJava [1] | HG002 SNV F1=0.6619 | HG002 InDel F1=0.7177 | 100% Exact Match (see Table 3) |
| Next-Gen Germline & Mosaic Suite | GIAB HG002 v4.2.1 [13], v5.0q T2T-Q100 [14, 16], CMRG [15] & Mosaic [20] | 12-Evaluation Next-Gen Germline & Mosaic Suite | Full Contigs (Rust) & <=50 kbp BED (VarDictRust == VarDictJava [1]) | v5.0q SNV F1=0.6650 (Full) / 0.6461 (BED) | CMRG SNV F1=0.6791 (Full) / 0.6132 (BED) | Mosaic Recall=98.82% (see Table 7) |
| HG008 v0.4 Somatic Landscape | GIAB HG008 v0.4 Somatic Benchmark [19] | 11-Caller Modern ML & Commercial Landscape | VarDictRust v1.0.6 Full WGS & Slice | Full WGS SNV F1=0.7886 | Full WGS InDel F1=0.4898 | 27 evaluations across 11 callers (see Table 7) |

*<sup>a</sup>VarDictRust and VarDictJava v1.8.3 [1] produce identical rtg vcfeval [12] TP, FP, FN, precision, recall, and F1 scores at the unfiltered (All-PASS) and optimal QUAL operating points of all 5 truth cohorts (displayed rows; per-cohort optima: HG002 >= 128, SEQC2 >= 126, HG008 >= 176). <sup>b</sup>Variant-class stratification compares VarDictRust, VarDictJava [1], and vardictcpp [11] across SNVs, insertions, deletions, complex MNPs, and structural variants (SVs), showing that vardictcpp splits complex MNPs on HG002 and omits SVs on SEQC2 and HG008 (CPP-2, CPP-9). <sup>c</sup>Across all 18/18 targeted Chromosome 20 somatic noise-floor runs (6 HG008 platforms x 3 low-VAF regimes: -f 0.01 -r 2, -f 0 -r 2, and singleton -f 0 -r 1), VarDictRust and VarDictJava emit identical candidate and PASS variant sets. <sup>d</sup>In silico tumor-purity titration admixtures were generated from GIAB HG008 [19] tumor and matched normal BAMs via samtools view -s [26] at 20%, 10%, 5%, 2%, and 1% tumor purity (10.0% to 0.5% expected clonal VAF) and called at -f 0.002; undefined 0/0 precision and F1 values at 0.0000 recall are normalized to 0.0000. <sup>e</sup>While VarDictRust matches VarDictJava [1] across all 5 titration purities, vardictcpp [11] omits ##contig headers (CPP-6) and exhibits MNP/tie-break divergence at 20% purity (SNV F1 0.7512 vs. 0.6914) plus non-benchmark indel differences at 10%–1% purity. Source ledger: truth, stratification, and titration values are generated from the stage5 and stage6 TSVs via master_gcp_metrics.json; this table is the formatted view.*

---

<a id="supp_table5_ga4gh_70cell_matrix"></a>

### Supplementary Table S5. Complete 70-Cell GA4GH v3.5 [17] Genomic Stratification Matrix (7 Strata x 10 Sequencing Platforms)

| Cohort | Sequencing Platform<sup>a</sup> | GA4GH v3.5 Stratum [17]<sup>b</sup> | SNV Precision [12] | SNV Recall [12] | SNV F1 [12] | InDel Precision / Recall / F1 [12]<sup>c</sup> |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| GIAB HG002 Germline [13, 18] | Illumina | Homopolymer (4–6 bp) | 0.4869 | 0.9645 | 0.6471 | 0.4874 / 0.9134 / 0.6356 |
| GIAB HG002 Germline [13, 18] | Illumina | Homopolymer (7–11 bp) | 0.6056 | 0.9278 | 0.7328 | 0.6079 / 0.9474 / 0.7406 |
| GIAB HG002 Germline [13, 18] | Illumina | Homopolymer (>=12 bp) | 0.6828 | 0.8704 | 0.7653 | 0.8869 / 0.7491 / 0.8122 |
| GIAB HG002 Germline [13, 18] | Illumina | Tandem Repeats | 0.5462 | 0.9126 | 0.6833 | 0.8782 / 0.8357 / 0.8564 |
| GIAB HG002 Germline [13, 18] | Illumina | Segmental Duplications | 0.5343 | 0.8980 | 0.6700 | 0.5650 / 0.8514 / 0.6792 |
| GIAB HG002 Germline [13, 18] | Illumina | Low Mappability | 0.5202 | 0.8309 | 0.6398 | 0.4419 / 0.8199 / 0.5742 |
| GIAB HG002 Germline [13, 18] | Illumina | RefSeq CDS | 0.3644 | 0.9831 | 0.5317 | 0.1209 / 1.0000 / 0.2157 |
| GIAB HG002 Germline [13, 18] | Element AVITI | Homopolymer (4–6 bp) | 0.9491 | 0.9685 | 0.9587 | 0.9371 / 0.9308 / 0.9339 |
| GIAB HG002 Germline [13, 18] | Element AVITI | Homopolymer (7–11 bp) | 0.9657 | 0.9325 | 0.9488 | 0.9626 / 0.9533 / 0.9579 |
| GIAB HG002 Germline [13, 18] | Element AVITI | Homopolymer (>=12 bp) | 0.9793 | 0.8765 | 0.9251 | 0.9882 / 0.8968 / 0.9403 |
| GIAB HG002 Germline [13, 18] | Element AVITI | Tandem Repeats | 0.8662 | 0.9383 | 0.9008 | 0.9430 / 0.8396 / 0.8883 |
| GIAB HG002 Germline [13, 18] | Element AVITI | Segmental Duplications | 0.7358 | 0.8980 | 0.8088 | 0.7485 / 0.8446 / 0.7937 |
| GIAB HG002 Germline [13, 18] | Element AVITI | Low Mappability | 0.6786 | 0.8419 | 0.7515 | 0.6406 / 0.7640 / 0.6969 |
| GIAB HG002 Germline [13, 18] | Element AVITI | RefSeq CDS | 0.9311 | 0.9887 | 0.9591 | 0.7857 / 1.0000 / 0.8800 |
| GIAB HG002 Germline [13, 18] | MGI DNBSEQ | Homopolymer (4–6 bp) | 0.5416 | 0.9632 | 0.6934 | 0.4645 / 0.9279 / 0.6191 |
| GIAB HG002 Germline [13, 18] | MGI DNBSEQ | Homopolymer (7–11 bp) | 0.5844 | 0.9284 | 0.7173 | 0.2319 / 0.9519 / 0.3729 |
| GIAB HG002 Germline [13, 18] | MGI DNBSEQ | Homopolymer (>=12 bp) | 0.4531 | 0.8796 | 0.5981 | 0.7337 / 0.8764 / 0.7987 |
| GIAB HG002 Germline [13, 18] | MGI DNBSEQ | Tandem Repeats | 0.4990 | 0.9326 | 0.6502 | 0.6606 / 0.8336 / 0.7371 |
| GIAB HG002 Germline [13, 18] | MGI DNBSEQ | Segmental Duplications | 0.5267 | 0.8960 | 0.6634 | 0.4792 / 0.8581 / 0.6150 |
| GIAB HG002 Germline [13, 18] | MGI DNBSEQ | Low Mappability | 0.3691 | 0.8269 | 0.5104 | 0.3342 / 0.7888 / 0.4695 |
| GIAB HG002 Germline [13, 18] | MGI DNBSEQ | RefSeq CDS | 0.4110 | 0.9831 | 0.5796 | 0.2703 / 0.9091 / 0.4167 |
| GIAB HG002 Germline [13, 18] | Ultima UG100 | Homopolymer (4–6 bp) | 0.8208 | 0.9409 | 0.8768 | 0.0195 / 0.7546 / 0.0381 |
| GIAB HG002 Germline [13, 18] | Ultima UG100 | Homopolymer (7–11 bp) | 0.5572 | 0.7097 | 0.6243 | 0.0514 / 0.4792 / 0.0928 |
| GIAB HG002 Germline [13, 18] | Ultima UG100 | Homopolymer (>=12 bp) | 0.2666 | 0.6944 | 0.3853 | 0.2957 / 0.4573 / 0.3591 |
| GIAB HG002 Germline [13, 18] | Ultima UG100 | Tandem Repeats | 0.7497 | 0.9100 | 0.8221 | 0.2424 / 0.7895 / 0.3709 |
| GIAB HG002 Germline [13, 18] | Ultima UG100 | Segmental Duplications | 0.8726 | 0.8755 | 0.8740 | 0.0277 / 0.7297 / 0.0535 |
| GIAB HG002 Germline [13, 18] | Ultima UG100 | Low Mappability | 0.8497 | 0.7799 | 0.8133 | 0.0359 / 0.6708 / 0.0682 |
| GIAB HG002 Germline [13, 18] | Ultima UG100 | RefSeq CDS | 0.7489 | 0.9568 | 0.8402 | 0.0016 / 0.7273 / 0.0033 |
| GIAB HG008 Somatic [19] | Illumina | Homopolymer (4–6 bp) | 0.9545 | 0.8750 | 0.9130 | 0.7778 / 0.4118 / 0.5385 |
| GIAB HG008 Somatic [19] | Illumina | Homopolymer (7–11 bp) | 1.0000 | 1.0000 | 1.0000 | 1.0000 / 0.7647 / 0.8667 |
| GIAB HG008 Somatic [19] | Illumina | Homopolymer (>=12 bp) | 0.8000 | 0.5000 | 0.6154 | 0.8235 / 0.3373 / 0.4786 |
| GIAB HG008 Somatic [19] | Illumina | Tandem Repeats | 0.1296 | 0.6364 | 0.2154 | 0.0563 / 0.2105 / 0.0889 |
| GIAB HG008 Somatic [19] | Illumina | Segmental Duplications | 0.9000 | 0.7500 | 0.8182 | 0.8333 / 0.5556 / 0.6667 |
| GIAB HG008 Somatic [19] | Illumina | Low Mappability | 0.3158 | 1.0000 | 0.4800 | 0.3000 / 1.0000 / 0.4615 |
| GIAB HG008 Somatic [19] | Illumina | RefSeq CDS | 1.0000 | 1.0000 | 1.0000 | No Truth InDels in Stratum<sup>c</sup> |
| GIAB HG008 Somatic [19] | Element AVITI | Homopolymer (4–6 bp) | 0.8958 | 0.8958 | 0.8958 | 0.9375 / 0.8824 / 0.9091 |
| GIAB HG008 Somatic [19] | Element AVITI | Homopolymer (7–11 bp) | 1.0000 | 1.0000 | 1.0000 | 1.0000 / 0.8824 / 0.9375 |
| GIAB HG008 Somatic [19] | Element AVITI | Homopolymer (>=12 bp) | 0.8000 | 0.5000 | 0.6154 | 0.9500 / 0.6867 / 0.7972 |
| GIAB HG008 Somatic [19] | Element AVITI | Tandem Repeats | 0.0909 | 0.6364 | 0.1591 | 0.1154 / 0.3158 / 0.1690 |
| GIAB HG008 Somatic [19] | Element AVITI | Segmental Duplications | 0.7143 | 0.8333 | 0.7692 | 0.8333 / 0.5556 / 0.6667 |
| GIAB HG008 Somatic [19] | Element AVITI | Low Mappability | 0.2683 | 0.9167 | 0.4151 | 0.2727 / 1.0000 / 0.4286 |
| GIAB HG008 Somatic [19] | Element AVITI | RefSeq CDS | 1.0000 | 1.0000 | 1.0000 | No Truth InDels in Stratum<sup>c</sup> |
| GIAB HG008 Somatic [19] | PacBio Onso | Homopolymer (4–6 bp) | 0.9286 | 0.8125 | 0.8667 | 0.6667 / 0.7059 / 0.6857 |
| GIAB HG008 Somatic [19] | PacBio Onso | Homopolymer (7–11 bp) | 0.7500 | 0.7500 | 0.7500 | 0.8125 / 0.7647 / 0.7879 |
| GIAB HG008 Somatic [19] | PacBio Onso | Homopolymer (>=12 bp) | 0.8000 | 0.5000 | 0.6154 | 0.9231 / 0.5783 / 0.7111 |
| GIAB HG008 Somatic [19] | PacBio Onso | Tandem Repeats | 0.0976 | 0.7273 | 0.1720 | 0.0625 / 0.1579 / 0.0896 |
| GIAB HG008 Somatic [19] | PacBio Onso | Segmental Duplications | 0.6000 | 0.7500 | 0.6667 | 0.6667 / 0.4444 / 0.5333 |
| GIAB HG008 Somatic [19] | PacBio Onso | Low Mappability | 0.2353 | 1.0000 | 0.3810 | 0.2143 / 1.0000 / 0.3529 |
| GIAB HG008 Somatic [19] | PacBio Onso | RefSeq CDS | 1.0000 | 1.0000 | 1.0000 | No Truth InDels in Stratum<sup>c</sup> |
| GIAB HG008 Somatic [19] | Roche SBX | Homopolymer (4–6 bp) | 0.8571 | 0.8750 | 0.8660 | 0.4000 / 0.2353 / 0.2963 |
| GIAB HG008 Somatic [19] | Roche SBX | Homopolymer (7–11 bp) | 0.8000 | 1.0000 | 0.8889 | 0.7000 / 0.4118 / 0.5185 |
| GIAB HG008 Somatic [19] | Roche SBX | Homopolymer (>=12 bp) | 1.0000 | 0.5000 | 0.6667 | 0.3913 / 0.1084 / 0.1698 |
| GIAB HG008 Somatic [19] | Roche SBX | Tandem Repeats | 0.7143 | 0.4545 | 0.5556 | 0.5000 / 0.1053 / 0.1739 |
| GIAB HG008 Somatic [19] | Roche SBX | Segmental Duplications | 0.8889 | 0.6667 | 0.7619 | 0.7500 / 0.3333 / 0.4615 |
| GIAB HG008 Somatic [19] | Roche SBX | Low Mappability | 0.7857 | 0.9167 | 0.8462 | 0.5000 / 1.0000 / 0.6667 |
| GIAB HG008 Somatic [19] | Roche SBX | RefSeq CDS | 1.0000 | 1.0000 | 1.0000 | No Truth InDels in Stratum<sup>c</sup> |
| GIAB HG008 Somatic [19] | Ultima ppmSeq | Homopolymer (4–6 bp) | 0.8667 | 0.8125 | 0.8387 | 0.2857 / 0.1176 / 0.1667 |
| GIAB HG008 Somatic [19] | Ultima ppmSeq | Homopolymer (7–11 bp) | 1.0000 | 0.2500 | 0.4000 | 0.0000 / 0.0000 / 0.0000 |
| GIAB HG008 Somatic [19] | Ultima ppmSeq | Homopolymer (>=12 bp) | 1.0000 | 0.1250 | 0.2222 | 1.0000 / 0.0241 / 0.0471 |
| GIAB HG008 Somatic [19] | Ultima ppmSeq | Tandem Repeats | 0.2500 | 0.6364 | 0.3590 | 0.0435 / 0.0526 / 0.0476 |
| GIAB HG008 Somatic [19] | Ultima ppmSeq | Segmental Duplications | 0.7273 | 0.6667 | 0.6957 | 0.5000 / 0.2222 / 0.3077 |
| GIAB HG008 Somatic [19] | Ultima ppmSeq | Low Mappability | 0.6154 | 0.6667 | 0.6400 | 1.0000 / 0.3333 / 0.5000 |
| GIAB HG008 Somatic [19] | Ultima ppmSeq | RefSeq CDS | 1.0000 | 1.0000 | 1.0000 | No Truth InDels in Stratum<sup>c</sup> |
| GIAB HG008 Somatic [19] | Ultima UG100 | Homopolymer (4–6 bp) | 0.9512 | 0.8125 | 0.8764 | 0.5000 / 0.1176 / 0.1905 |
| GIAB HG008 Somatic [19] | Ultima UG100 | Homopolymer (7–11 bp) | 1.0000 | 0.2500 | 0.4000 | 0.0000 / 0.0000 / 0.0000 |
| GIAB HG008 Somatic [19] | Ultima UG100 | Homopolymer (>=12 bp) | 1.0000 | 0.1250 | 0.2222 | 0.2308 / 0.0361 / 0.0625 |
| GIAB HG008 Somatic [19] | Ultima UG100 | Tandem Repeats | 0.6250 | 0.4545 | 0.5263 | 0.1429 / 0.0526 / 0.0769 |
| GIAB HG008 Somatic [19] | Ultima UG100 | Segmental Duplications | 1.0000 | 0.8333 | 0.9091 | 1.0000 / 0.2222 / 0.3636 |
| GIAB HG008 Somatic [19] | Ultima UG100 | Low Mappability | 0.9167 | 0.9167 | 0.9167 | 0.5000 / 0.3333 / 0.4000 |
| GIAB HG008 Somatic [19] | Ultima UG100 | RefSeq CDS | 1.0000 | 1.0000 | 1.0000 | No Truth InDels in Stratum<sup>c</sup> |

*<sup>a</sup>Evaluates Chromosome 20 callsets at -f 0.01 across 4 GIAB HG002 [13, 18] germline platforms and 6 GIAB HG008 [19] somatic platforms; across all 70/70 stratification cells (Clopper-Pearson 95% CI 94.87%–100.00%), VarDictRust v1.0.6 and VarDictJava v1.8.3 [1] exhibit 0 site differences. <sup>b</sup>Genomic strata follow the NIST GIAB / GA4GH v3.5 [17] stratification BED definitions on GRCh38 (homopolymers of 4–6 bp, 7–11 bp, and >=12 bp, all tandem repeats, segmental duplications, low-mappability regions, and RefSeq coding sequences). <sup>c</sup>Within the Chromosome 20 RefSeq CDS stratum of GIAB HG008 v0.4 [19], the somatic truth set contains SNVs only and zero truth InDels. Source ledger: stratification values are generated from the stage6 TSVs via master_gcp_metrics.json; this table is the formatted view.*

---

<a id="supp_table6_sv_and_amplicon_splice"></a>

### Supplementary Table S6. 13-Mode Functional Parity Matrix and CPP-1 to CPP-10 vardictcpp [11] Divergence Taxonomy

| Evaluation Category | CLI Mode / Divergence Class | VarDictJava v1.8.3 [1] Rows / Behavior | VarDictRust v1.0.6 Rows / Behavior | vardictcpp [11] Rows / Behavior | Verification Verdict & Root Cause |
| :--- | :--- | :--- | :--- | :--- | :--- |
| 13-Mode Functional Parity Suite<sup>a</sup> | Simple Germline (-f 0.01) | 228 | 228 (FP=0, FN=0) | 228 (FP=0, FN=0) | 100% Exact Match (VarDictRust == VarDictJava) |
| 13-Mode Functional Parity Suite<sup>a</sup> | Simple Zero-Threshold MRD (-f 0) | 248 | 248 (FP=0, FN=0) | 248 (FP=0, FN=0) | 100% Exact Match (VarDictRust == VarDictJava) |
| 13-Mode Functional Parity Suite<sup>a</sup> | Simple Chunked Mode (--chunk 50,000) | 228 | 228 (FP=0, FN=0) | 228 (FP=0, FN=0) | 100% Exact Match (VarDictRust == VarDictJava) |
| 13-Mode Functional Parity Suite<sup>a</sup> | Simple Germline with Native Fisher (--fisher)<sup>a</sup> | 650 | 650 (FP=0, FN=0) | 650 (FP=1, FN=1) | 100% Exact Match (VarDictRust == VarDictJava) |
| 13-Mode Functional Parity Suite<sup>a</sup> | Targeted Clinical Panel (Chr17) | 83 | 83 (FP=0, FN=0) | 83 (FP=0, FN=0) | 100% Exact Match (VarDictRust == VarDictJava) |
| 13-Mode Functional Parity Suite<sup>a</sup> | Structural Variant Deletion Regions | 3 | 3 (FP=0, FN=0) | 3 (FP=0, FN=0) | 100% Exact Match (VarDictRust == VarDictJava) |
| 13-Mode Functional Parity Suite<sup>a</sup> | Low-Coverage Boundary Calling | 146 | 146 (FP=0, FN=0) | 146 (FP=0, FN=0) | 100% Exact Match (VarDictRust == VarDictJava) |
| 13-Mode Functional Parity Suite<sup>a</sup> | Amplicon Calling (-a, Chr17) | 145 | 145 (FP=0, FN=0) | 145 (FP=0, FN=0) | 100% Exact Match (VarDictRust == VarDictJava) |
| 13-Mode Functional Parity Suite<sup>a</sup> | Amplicon Calling with Native Fisher (-a --fisher) | 145 | 145 (FP=0, FN=0) | 145 (FP=0, FN=0) | 100% Exact Match (VarDictRust == VarDictJava) |
| 13-Mode Functional Parity Suite<sup>a</sup> | Paired Tumor-Normal Somatic (-f 0.01) | 664 | 664 (FP=0, FN=0) | 664 (FP=0, FN=0) | 100% Exact Match (VarDictRust == VarDictJava) |
| 13-Mode Functional Parity Suite<sup>a</sup> | Paired Somatic Zero-Threshold MRD (-f 0) | 1,396 | 1,396 (FP=0, FN=0) | 1,396 (FP=0, FN=0) | 100% Exact Match (VarDictRust == VarDictJava) |
| 13-Mode Functional Parity Suite<sup>a</sup> | Paired Somatic with Native Fisher (--fisher) | 527 | 527 (FP=0, FN=0) | 527 (FP=0, FN=0) | 100% Exact Match (VarDictRust == VarDictJava) |
| 13-Mode Functional Parity Suite<sup>a</sup> | Lossless CRAM-vs-BAM Twin (-G) [26]<sup>b</sup> | Fails (NoClassDefFoundError on CRAM)<sup>b</sup> | 228 (FP=0, FN=0) | 228 (FP=0, FN=0) | 100% Exact Match vs. BAM Twin<sup>b</sup> |
| vardictcpp Divergence Taxonomy [11]<sup>c</sup> | CPP-1: VCF Header & INFO/END Schema | Standard VCF 4.2 + INFO/END tag | Exact v1.8.3 header & INFO/END parity | Malformed or omitted INFO/END & header lines | Schema mismatch in downstream VCF parsers |
| vardictcpp Divergence Taxonomy [11]<sup>c</sup> | CPP-2: CRISPR (-J/-j) & Structural Variant Modes | Supported | Supported (100% Exact Match) | Unimplemented / drops SV calls | Verified across 17-modality extended coverage suite |
| vardictcpp Divergence Taxonomy [11]<sup>c</sup> | CPP-3: Amplicon (-a) Output Format | 38/40-col .var + VCF via var2vcf_valid.pl | Exact 38/40-col .var + VCF parity | Emits raw intermediate TSV instead of VCF | Incompatible with standard VCF comparators & bcftools [26] |
| vardictcpp Divergence Taxonomy [11]<sup>c</sup> | CPP-4: Concatenated CLI Flag Parsing | Supported | Supported (60-flag CLI parity) | Exits with code 120 on bundled CLI flags | Breaks existing VarDict pipeline CLI wrappers |
| vardictcpp Divergence Taxonomy [11]<sup>c</sup> | CPP-5: Chunk Boundary Deduplication (--chunk 50,000) | Single-region calling (or OOM on contigs) | 100 kbp halo tiling + cross-tile dedup | Emits +309,569 duplicate boundary calls on unchunked HG002 WGS | Uncoordinated 50 kbp window splitting without dedup |
| vardictcpp Divergence Taxonomy [11]<sup>c</sup> | CPP-6: Missing ##contig VCF Header Lines | Emits ##contig lines from reference .fai | Emits exact ##contig lines from reference .fai | Omits ##contig header lines | Requires manual header injection prior to rtg vcfeval [12] |
| vardictcpp Divergence Taxonomy [11]<sup>c</sup> | CPP-7: GRCh37 / b37 Numeric Chromosome Naming | Supported | Supported (100% Match on GIAB NA24385 GRCh37 [13]) | Hardcoded chr-prefix assumptions | Verified on GIAB NA24385 GRCh37 v4.2.1 benchmark [13] |
| vardictcpp Divergence Taxonomy [11]<sup>c</sup> | CPP-8: BAM/CRAM [26] Reference & Decode | BAM (CRAM requires external htsjdk classpath) | Native BAM + CRAM (-G) zero-difference parity | Fails or diverges on CRAM twin alignments | Verified on GIAB HG008 [19] & Ultima CRAM twins |
| vardictcpp Divergence Taxonomy [11]<sup>c</sup> | CPP-9: Complex MNP Splitting & Somatic SV Drops | Reference pileup & local realigner | Exact parity with VarDictJava [1] | Differs on 394 callable HG002 WGS calls & drops SVs | Local realignment & tie-break divergence vs. oracle |
| vardictcpp Divergence Taxonomy [11]<sup>c</sup> | CPP-10: Transcriptome RNA-Seq Splice Mode (-i) | Supported (12,668 rows on HG002 mRNA [13]) | Supported (12,668 rows, 100% Exact Match) | Unimplemented (unsupported flag) | Verified on 51.8 GiB GIAB HG002 GM24385 STAR BAM [13] |

*<sup>a</sup>Evaluates all 13 functional CLI modes across VarDictJava v1.8.3 [1], VarDictRust v1.0.6, and vardictcpp [11]; VarDictRust achieves 0 false positives (FP=0) and 0 false negatives (FN=0) across all 13/13 modes, whereas vardictcpp exhibits FP=1, FN=1 on simple_fisher. <sup>b</sup>Stock Bioconda VarDictJava v1.8.3 [1] fails with java.lang.NoClassDefFoundError on CRAM inputs; VarDictRust decodes CRAM natively via pure-Rust noodles and matches the lossless BAM twin oracle with 0 differences. <sup>c</sup>Catalogs the 10 functional, schema, and algorithmic divergence classes (CPP-1 through CPP-10) empirically observed in vardictcpp (commit ae93fe39dc) [11] during the 60-experiment benchmark study. Source ledger: mode-parity rows are generated from the sweep matrix via master_gcp_metrics.json; this table is the formatted view.*

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
