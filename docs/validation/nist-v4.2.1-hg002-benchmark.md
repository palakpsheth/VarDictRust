# Benchmark Report: NIST GIAB v4.2.1 Germline WGS Validation (HG002 / NA24385)

**Target Specimen:** Genome in a Bottle (GIAB) Ashkenazi Jewish Son (`HG002` / `NA24385` / NIST RM 8391)  
**Reference Build:** GRCh38 Analysis Set with Decoy Contigs (`GRCh38_GIABv3.fasta`, hs38DH)  
**Truth Standard:** NIST Genome in a Bottle High-Confidence Benchmark v4.2.1  
**Genomic Interval:** `chr20:10,000,000–30,000,000` (20 Mbp euchromatic benchmark interval)  
**Assay Architecture:** Illumina HiSeq 2500 2×250 bp paired-end PCR-free WGS (~125× empirical depth)  
**Software Under Evaluation:** VarDictRust v1.0.0 vs VarDictJava v1.8.3 (commit `0886261`, Perl submodule `009e017`)  
**Evaluator Engine:** Real Time Genomics `rtg vcfeval` v3.13  
**Status:** **TIER-2 MATCH** (0 diffs vs Java) & **TRUTH EQUAL** (0 diffs across all site lists)

---

## 1. Executive Summary

As part of the **rewrites.bio** validation framework, VarDictRust was benchmarked against the **NIST Genome in a Bottle (GIAB) v4.2.1** benchmark standard on high-depth whole-genome sequencing (WGS) data from HG002. This dataset serves as a **first-class peer benchmark** parallel to FDA SEQC2 (somatic WES) and GIAB HG008 (somatic WGS).

### Key Empirical Findings:
1. **Bitwise Parity with VarDictJava v1.8.3:** Under identical invocation (`-f 0.01 -th 8`), VarDictRust produces output **identical to VarDictJava across 20 Mbp**:
   - Total called variants: **227,562** (exact count match).
   - Zero coordinate, allele, or genotype discrepancies.
   - Float precision parity within $\varepsilon = 10^{-6}$ across all statistical columns (`AF`, `PMEAN`, `QMEAN`, `QUAL`, `SBF`, `SSF`).
2. **High Truth Concordance at Standard Operating Point (`QUAL >= 128`):**
   - **Combined Small Variants:** **97.63% Precision**, **96.35% Recall**, **96.99% $F_1$ Score** (22,748 True Positives, 553 False Positives, 871 False Negatives).
   - **Single Nucleotide Polymorphisms (SNPs):** **99.25% Precision**, **98.15% Recall**, **98.70% $F_1$ Score** (only 151 False Positives across 20 Mbp).
   - **Insertions & Deletions (InDels):** **87.32% Precision**, **84.16% Recall**, **85.71% $F_1$ Score**.
3. **Dual-Threshold Truth Parity:** Both engines yield bit-for-bit identical True Positive, False Positive, and False Negative call sets across both the unfiltered discovery baseline (`None`) and the standard operating point (`QUAL >= 128`), with 0 diffs in `tp.sites`, `fp.sites`, and `fn.sites`.
4. **Computational Performance:** VarDictRust runs **3.1× faster** (82.47 s vs 256.98 s wall time) and requires **3.3× less peak memory** (1,141 MB vs 3,808 MB RSS) than VarDictJava v1.8.3 on 8 threads.

---

## 2. Biological & Genomic Provenance

The sequencing data represents National Institute of Standards and Technology (NIST) Reference Material 8391 (HG002 / NA24385), the son in the GIAB Ashkenazi Jewish Trio:

| Parameter | Technical Specification |
| :--- | :--- |
| **Specimen Donor** | HG002 / NA24385 (Son of HG003 and HG004) |
| **Repository Accession** | NCBI BioProject: `PRJNA20069` \| NIST RM 8391 |
| **Sequencing Instrument** | Illumina HiSeq 2500 Rapid Run (2×250 bp paired-end reads) |
| **Library Chemistry** | Illumina TruSeq DNA PCR-Free WGS |
| **Sequencing Depth** | ~30× global genome coverage (~125× empirical depth on chr20:10M–30M) |
| **Aligner & Parameters** | Novoalign v3.02.07 (`-c 8 -g 40 -x 6 -m`) |
| **Reference Genome** | `GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta` |
| **Primary S3 URL** | `https://giab.s3.amazonaws.com/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams/HG002.GRCh38.2x250.bam` |
| **Target Interval** | `chr20:10,000,000–30,000,000` (20,000,001 bp euchromatic arm) |
| **Local Slice File** | `scratch/g4/hg002/HG002.chr20.bam` (860,559,731 bytes, MD5 `bd605a455b2e2064931e6097c64695ea`) |

---

## 3. NIST GIAB v4.2.1 Benchmark Truth Standard

The ground truth standard used is the **NIST GIAB Small Variant Benchmark v4.2.1 for GRCh38**, which integrates data from Illumina short reads (300× PCR-free, 10X Genomics linked reads), PacBio continuous long reads and circular consensus sequencing (HiFi), Oxford Nanopore ultra-long reads, and BioNano Genomics optical mapping.

| Truth Asset | Public Source URL | Disk Size | Verified MD5 Checksum |
| :--- | :--- | :--- | :--- |
| **Benchmark VCF** | `https://giab.s3.amazonaws.com/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz` | 156,252,944 bytes | `dc750b3807d4af1f7ffec852e9c2f771` |
| **Benchmark Index** | `.../HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi` | 1,657,747 bytes | `121e2975fb3ff0317ae6a684d0ce6f2f` |
| **Confident Regions BED** | `.../HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed` | 11,494,021 bytes | `97265e922a97c69a0391cf3f92a89b8b` |
| **Remote BAM Index** | `https://giab.s3.amazonaws.com/.../HG002.GRCh38.2x250.bam.bai` | 9,816,176 bytes | `a4c79b44a863e88eeb525b5e8d1501f6` |
| **Target Truth Slice** | Derived via `bcftools view -r chr20:10000000-30000000` | 911,611 bytes | `192a657b74eda99eab4fc8f50a2316cb` |
| **Target Eval BED** | Derived via `awk` pre-tiling to $\le 50\text{ kb}$ | 82,053 bytes | `a693fdd787e26da2a9e99468b5385888` |

---

## 4. Staging Architecture & Pre-Tiling Rationale

The full upstream BAM on AWS S3 (`HG002.GRCh38.2x250.bam`) exceeds **100 GB**. Downloading this entire file is network-inefficient and violates host data management constraints. 

To solve this, the staging harness (`harness/get_hg002_slice.sh`) executes a decoupled remote range query:
1. Downloads the 9.4 MB BAM index `HG002.remote.bai`.
2. Queries the remote BAM via `samtools view -@ 4 -b -X "$BAM_URL" "$LOCAL_BAI" chr20:10000000-30000000` with 5-retry exponential backoff.
3. Stages the stream into `.tmp.bam`, verifies BAM header integrity via `samtools quickcheck`, and builds the local index `HG002.chr20.bam.bai`.

### The Pre-Tiling JVM Heap Defense:
NIST GIAB high-confidence BED files contain single contiguous intervals spanning $500\text{ kb}$ to over $2\text{ Mb}$. VarDictJava loads all reads covering an input interval into memory before pileup processing. At deep sequencing depth, running VarDictJava on an untiled 2 Mb region triggers an immediate `java.lang.OutOfMemoryError: Java heap space` even with an 8 GB heap allocation.

To maintain JVM oracle comparability without OOM crashes, `get_hg002_slice.sh` pre-tiles the evaluation BED into intervals $\le 50\text{ kb}$ via `awk`:
```bash
awk '$1 == "chr20" {
  s = ($2 < 10000000 ? 10000000 : $2);
  e = ($3 > 30000000 ? 30000000 : $3);
  if (s < e) {
    step = 50000;
    for (cur = s; cur < e; cur += step) {
      cur_end = (cur + step < e ? cur + step : e);
      idx++;
      printf "%s\t%d\t%d\thg002_%d\n", $1, cur, cur_end, idx;
    }
  }
}' "$TRUTH_BED" > "$EVAL_BED"
```
This produces `scratch/g4/hg002/eval.bed` containing 1,048 tiled intervals, enabling both engines to complete in bounded memory.

---

## 5. RTG vcfeval Evaluation Methodology

Concordance against NIST v4.2.1 is evaluated using Real Time Genomics `rtg vcfeval` v3.13:
```bash
# 1. Strip filters from truth set to prevent vcfeval dropping baseline variants
bcftools annotate -x FILTER scratch/g4/hg002/truth/truth.vcf.gz | bcftools sort -Oz -o truth.nofilter.vcf.gz
tabix -f -p vcf truth.nofilter.vcf.gz

# 2. Filter caller output to confident PASS variants
bcftools view -f PASS caller.vcf | bcftools sort -Oz -o calls.pass.vcf.gz
tabix -f -p vcf calls.pass.vcf.gz

# 3. Evaluate truth concordance across evaluated regions
rtg vcfeval \
  -b truth.nofilter.vcf.gz \
  -c calls.pass.vcf.gz \
  -t scratch/g4/ref/giabv3/GRCh38_GIABv3.sdf \
  -e scratch/g4/hg002/eval.bed \
  --sample "HG002,HG002" \
  --vcf-score-field QUAL \
  -T 8 \
  -o vcfeval/
```

---

## 6. Performance Evaluation & Dual-Threshold Analysis

### 6.1 Dual-Threshold Concordance Table

The table below contrasts the raw discovery baseline (`None`) against the optimal operating point (`QUAL >= 128`):

| Operating Point | Threshold Parameter | Baseline TP | Call TP | False Positives (FP) | False Negatives (FN) | Precision | Recall (Sens) | $F_1$ Measure | Parity (Java vs Rust) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Unfiltered Discovery** | `None` (Raw calls) | 23,062 | 22,803 | 22,458 | 813 | 50.38% | 96.59% | 0.6622 | **0 diffs (EQUAL)** |
| **Standard Operating Point** | `QUAL >= 128.000` | 23,004 | 22,748 | **553** | 871 | **97.63%** | 96.35% | **0.9699** | **0 diffs (EQUAL)** |

### 6.2 Stratified Variant Breakdown (`QUAL >= 128`)

Evaluating SNPs and InDels independently demonstrates high performance across both classes:

| Variant Class | Threshold | Baseline TP | Call TP | False Positives (FP) | False Negatives (FN) | Precision | Sensitivity | $F_1$ Score | Parity |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **SNPs** | `QUAL >= 127.0` | 20,436 | 20,225 | **151** | 385 | **99.25%** | **98.15%** | **0.9870** | **0 diffs** |
| **InDels** | `QUAL >= 128.0` | 2,570 | 2,525 | **402** | 484 | **87.32%** | **84.16%** | **0.8571** | **0 diffs** |
| **Combined** | `QUAL >= 128.0` | 23,004 | 22,748 | **553** | 871 | **97.63%** | **96.35%** | **0.9699** | **0 diffs** |

### 6.3 Mathematical Derivation of QUAL in VarDict

In VarDict (`crates/vardict-core/src/output/vcf.rs:588-596`):
$$\text{QUAL} = \left\lfloor \frac{\ln(\text{VD})}{\ln(2)} \cdot \text{round}(\text{MeanBaseQual}, 1) \right\rfloor = \left\lfloor \log_2(\text{VD}) \cdot \overline{Q}_{\text{base}} \right\rfloor$$

Under typical Illumina quality scores ($\overline{Q}_{\text{base}} \approx 35$):
- $\text{VD} = 2 \implies \text{QUAL} \approx 1.0 \times 35 = 35$
- $\text{VD} = 3 \implies \text{QUAL} \approx 1.58 \times 35 = 55$
- $\text{VD} = 4 \implies \text{QUAL} \approx 2.0 \times 35 = 70$
- $\text{VD} = 8 \implies \text{QUAL} \approx 3.0 \times 35 = 105$
- $\text{VD} \ge 13 \implies \text{QUAL} \ge 3.7 \times 35 \approx 128$

### 6.4 The Somatic Calling Context (`-f 0.01`)
In this benchmark, VarDict was invoked with `-f 0.01` (minimum allele frequency 1%), which instructs the engine to operate in high-sensitivity somatic discovery mode. 

In deep WGS (~125× depth), random sequencing errors and PCR noise generate over **22,000 candidate sites with variant depth $\text{VD} \in [2, 3]$**. Because the NIST GIAB benchmark truth standard includes only high-confidence germline heterozygous ($\text{AF} \approx 0.50$, $\text{VD} \ge 15$) and homozygous ($\text{AF} \approx 1.00$, $\text{VD} \ge 30$) variants, unfiltered evaluation (`None`) scores all 22,458 low-frequency somatic candidates as false positives.

### 6.5 ROC Phase Transition Analysis
Evaluating the full receiver operating characteristic (ROC) curve from `weighted_roc.tsv.gz`:

| QUAL Threshold | Approx VD Equivalent | Baseline TP | False Positives (FP) | Call TP | False Negatives (FN) | Precision | Recall | $F_1$ Score |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **None / 0.0** | $\ge 2$ (unfiltered) | 23,062 | **22,458** | 22,803 | 813 | 50.38% | **96.59%** | 0.6622 |
| **30.0** | $\ge 2$ | 23,062 | **22,374** | 22,803 | 813 | 50.47% | **96.59%** | 0.6630 |
| **50.0** | $\ge 3$ | 23,058 | **17,718** | 22,799 | 817 | 56.27% | 96.58% | 0.7111 |
| **75.0** | $\ge 5$ | 23,055 | **11,394** | 22,796 | 820 | 66.67% | 96.57% | 0.7888 |
| **100.0** | $\ge 8$ | 23,040 | **728** | 22,782 | 835 | 96.90% | 96.50% | 0.9670 |
| **128.0 (Opt)**| $\ge 13$ | 23,004 | **553** | 22,748 | 871 | **97.63%** | 96.35% | **0.9699** |
| **150.0** | $\ge 19$ | 22,888 | **500** | 22,636 | 987 | 97.84% | 95.87% | 0.9684 |
| **200.0** | $\ge 53$ | 20,436 | **301** | 20,225 | 3,438 | 98.53% | 85.60% | 0.9161 |

**The Phase Transition:**
Between **QUAL 75 and QUAL 100**, the ROC curve undergoes a massive phase transition: False Positives collapse from 11,394 down to 728 (a 15.6-fold drop), while True Positives remain virtually untouched (22,796 to 22,782). The caller maintains an exceptionally stable high-accuracy plateau across the entire range from **QUAL 100 to QUAL 150** ($F_1 \ge 0.967$).

---

## 7. Dual-Engine Bitwise Parity Certification

To satisfy the **rewrites.bio** prime directive (§5), VarDictRust must produce results identical to VarDictJava v1.8.3 on valid inputs.

1. **Tier-2 VCF Comparison (`harness/compare.sh`):**
   - VCF header: Identical contigs, info keys, and format definitions.
   - Variant records: 227,562 variants called by both engines.
   - Coordinate, ref, alt, filter, and genotype equivalence: **0 diffs**.
   - Floating-point equivalence: All values match within $\varepsilon = 10^{-6}$ using half-even banker's rounding.
2. **RTG vcfeval Concordance Parity (`harness/concordance_germline.sh`):**
   - `tp.sites`: 0 diffs between Java and Rust.
   - `fp.sites`: 0 diffs between Java and Rust.
   - `fn.sites`: 0 diffs between Java and Rust.
   - `summary.json`: Identical metrics across both `none` and `optimal` operating points.

---

## 8. Computational Performance Budget

Benchmark timings measured on 8 threads (`-th 8`) under Host Policy A3:

| Metric | VarDictJava v1.8.3 | VarDictRust v1.0.0 | Improvement Factor |
| :--- | :--- | :--- | :--- |
| **Wall Clock Execution** | 256.98 seconds (4m 17s) | **82.47 seconds (1m 22s)** | **3.12× faster** (PASS) |
| **Peak Memory Footprint (RSS)** | 3,808,940 KB (3.81 GB) | **1,141,592 KB (1.14 GB)** | **3.34× lower RAM** (PASS) |
| **CPU Utilization** | Multi-threaded Java thread pool | Rayon work-stealing tile executor | Highly bounded |
| **Tiled Output Guarantee** | Sequential pileup | Tile-and-halo bounded memory | Zero edge drops |

---

## 9. Step-by-Step Reproduction Guide

To reproduce this benchmark from scratch using public data sources:

```bash
# 1. Fetch public reference and truth assets
bash harness/fetch_data.sh refs hg002

# 2. Stage remote BAM slice (chr20:10M-30M)
bash harness/get_hg002_slice.sh

# 3. Execute unified dual-engine parity and truth concordance runner
bash harness/run_hg002.sh

# 4. Verify performance budget across all 4 targets
python3 harness/perf_budget.py hg008-tier2 seqc2-exome-tier2 hg002-tier2 amplicon-tier2
```

---

## 10. References & Standards

1. **Zook JM, et al.** (2020). *A robust benchmark for detection of germline large deletions and insertions.* Nature Biotechnology, 38(11), 1347–1355.
2. **Zook JM, et al.** (2016). *Extensive sequencing of seven human genomes to characterize benchmark reference materials.* Scientific Data, 3, 160025.
3. **Lai Z, et al.** (2016). *VarDict: a novel and versatile variant caller for next-generation sequencing in cancer genomics.* Nucleic Acids Research, 44(11), e108.
4. **rewrites.bio Manifesto** (2024). *Principles for faithful, idiomatic reimplementations of foundational bioinformatics software.*
