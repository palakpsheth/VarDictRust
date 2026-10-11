# Multi-Slice Multi-Modal Public Dataset Expansion & Complete Documentation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Ingest, stream, and evaluate 18 non-overlapping biological slices across 6 diverse public datasets spanning 12 distinct human chromosomes to prove 100% byte-for-byte and strict float parity against VarDictJava v1.8.3 on valid inputs, author exhaustive publication-grade documentation, and certify Task 24h for the v1.0.0 release.

**Architecture:** Bounded-interval HTTP range request streaming via `samtools view -b -X` for remote WGS datasets (Hartwig COLO829, Illumina Platinum Genomes NA12878, GIAB HG005) coupled with local high-coverage slices for somatic exome (SEQC2 HCC1395), targeted amplicon (Horizon/L861Q), and RNA-seq transcriptome (STAR junctions). Standalone dual-engine runner executing VarDictJava v1.8.3 oracle against VarDictRust v1.0.0 with automated bitwise diffs and strict VCF comparison. Comprehensive publication-grade documentation report detailing genomic rationale, mutational spectra, code paths exercised, and reproducibility commands.

**Tech Stack:** Rust (release binary `target/release/vardictrust`), Java 8 / Conda (`vardict-oracle`), SAMtools 1.18+, Python 3 (NumPy, SciPy), Bash, Git.

**Spec:** `multislice_multimodal_dataset_expansion_plan.md`

## Global Constraints

- rewrites.bio Manifesto: §0 evidence before completion; §0b integrity; §5 exact parity on valid inputs.
- Certified Release Engine Fingerprint: Must remain strictly invariant at `7e9b1be766f73b17` (zero `crates/` modifications; streaming harness, slices, evidence, and docs expanded).
- Host Policy A3: Thread concurrency <= 8 (`-th 8`), memory wrapping under `courtesy.sh`, respectful bandwidth, zero whole-genome downloads (stream 20 kb to 500 kb slices via O(log N) BAM index range requests).
- Host Policy A10: Zero `Co-Authored-By` commit trailers.
- Success Criterion: 18/18 MATCH (0 variant diffs, 0 filter diffs, float eps <= 1e-6) against VarDictJava v1.8.3 live oracle.

---

### Task 1: Create Multi-Slice BED Intervals

**Files:**
- Create: `harness/slices/multimodal/colo829.bed`
- Create: `harness/slices/multimodal/platinum.bed`
- Create: `harness/slices/multimodal/seqc2.bed`
- Create: `harness/slices/multimodal/hg005.bed`
- Create: `harness/slices/multimodal/amplicon.bed`
- Create: `harness/slices/multimodal/rna.bed`

**Interfaces:**
- Consumes: Standard 0-based half-open BED3/BED4 format with chromosome, start, end, and optional gene/target name.
- Produces: 6 version-controlled BED files defining exactly 3 non-overlapping evaluation intervals per dataset across 12 distinct chromosomes.

- [ ] **Step 1: Write the BED files for all 6 datasets**

```bash
mkdir -p harness/slices/multimodal
```

Write `harness/slices/multimodal/colo829.bed` (b37 coordinates: 1, 3, 20):
```tsv
1	10000000	10200000	UV_Dinucleotide_Region
3	10000000	10200000	Somatic_LOH_Region
20	45000000	45200000	SV_Breakpoint_Cluster
```

Write `harness/slices/multimodal/platinum.bed` (GRCh38 coordinates: chr2, chr17, chrX):
```tsv
chr2	50000000	50200000	Large_Autosome_Baseline
chr17	40000000	40200000	High_GC_BRCA1_Locus
chrX	15000000	15200000	Sex_Chromosome_NonAutosomal
```

Write `harness/slices/multimodal/seqc2.bed` (GRCh38 coordinates: chr7, chr11, chr12):
```tsv
chr7	140400000	140600000	BRAF_Oncogenic_Hotspot
chr11	530000	550000	HRAS_Exon_Hotspot
chr12	25200000	25400000	KRAS_Target_Region
```

Write `harness/slices/multimodal/hg005.bed` (GRCh38 coordinates: chr3, chr6, chr21):
```tsv
chr3	40000000	40100000	Han_Chinese_Ancestral_Haplotypes
chr6	30000000	30100000	MHC_HLA_Polymorphism_Region
chr21	20000000	20100000	Acrocentric_SegDup_Boundary
```

Write `harness/slices/multimodal/amplicon.bed` (hg19 coordinates: chr7 EGFR hotspots):
```tsv
chr7	55241000	55243000	EGFR_Exon19_Deletions
chr7	55248000	55250000	EGFR_L858R_Hotspot
chr7	55259000	55260500	EGFR_T790M_Resistance
```

Write `harness/slices/multimodal/rna.bed` (hg38 coordinates: chr8 splicing junctions):
```tsv
chr8	90534611	90550000	Exon_Intron_Junction_1
chr8	90550000	90580000	Multi_Exon_Internal_Splicing
chr8	90580000	90620070	Exon_Intron_Junction_2
```

- [ ] **Step 2: Verify BED syntax and chromosome coordinate compatibility**

Run: `awk '{print $1, $2, $3, $3-$2}' harness/slices/multimodal/*.bed`
Expected: 18 non-overlapping intervals across 6 files, each interval between 1.5 kb and 200 kb.

---

### Task 2: Implement Slicing and Ingestion Script

**Files:**
- Create: `harness/get_multimodal_slices.sh`

**Interfaces:**
- Consumes: Remote BAM/BAI URLs (ENA, Google Cloud, NCBI GIAB) and local repository assets (`scratch/g4/seqc2/`, `VarDictJava/`).
- Produces: 18 BAM slices and index files in `scratch/g4/multimodal/`, each passing `samtools quickcheck`.

- [ ] **Step 1: Write `harness/get_multimodal_slices.sh`**

Script must support:
- Creating `scratch/g4/multimodal/`
- Streaming remote slices via `samtools view -b -X <bam_url> <bai_url> <region> > <out_bam>`
- Extracting local slices for SEQC2, amplicon, and RNA-seq
- Generating BAM index `.bai` via `samtools index`
- Validating every BAM with `samtools quickcheck`

- [ ] **Step 2: Execute `harness/get_multimodal_slices.sh` and verify all slices**

Run: `bash harness/get_multimodal_slices.sh`
Expected: 18 BAM files generated, total disk size < 150 MB, all 18 pass `samtools quickcheck`.

---

### Task 3: Implement Automated Dual-Engine Parity Runner

**Files:**
- Create: `harness/run_multimodal_slices.sh`

**Interfaces:**
- Consumes: 18 BAM slices, reference FASTAs, and BED intervals.
- Produces: Dual-engine output files in `/tmp/vardict_multimodal_*/`, differential comparison logs, summary scoreboard, and machine-readable `benchmarks/evidence/multimodal-slices-tier2.json`.

- [ ] **Step 1: Write `harness/run_multimodal_slices.sh`**

Implement the 18-slice evaluation matrix:
1. Slices 1A-1C: COLO829 tumor-normal somatic mode (`-b tumor -B normal`, b37, `--fisher`).
2. Slices 2A-2C: NA12878 germline mode (GRCh38, SimpleMode).
3. Slices 3A-3C: SEQC2 somatic mode (HCC1395 tumor, HCC1395BL normal, GRCh38, `--fisher`).
4. Slices 4A-4C: HG005 germline mode (GRCh38, SimpleMode).
5. Slices 5A-5C: Amplicon targeted mode (hg19, `-a`, amplicon 38-column / 40-column).
6. Slices 6A-6C: RNA-seq splicing mode (hg38, `-i --emit-intermediate`, 4-column intron table).

Comparator rules:
- Table outputs (RNA splicing): `diff -u <(sort java) <(sort rust)`.
- VCF outputs: `python3 harness/compare_vcf.py --strict-bytes java.vcf rust.vcf --eps 1e-6`.
- Assert: 0 diff lines, 0 variant mismatches, 0 filter mismatches.
- Generate `benchmarks/evidence/multimodal-slices-tier2.json`.

- [ ] **Step 2: Run `harness/run_multimodal_slices.sh` under Host Policy A3**

Run: `EST_MIN=15 bash harness/courtesy.sh bash harness/run_multimodal_slices.sh`
Expected: 18/18 MATCH (0 diffs across all 18 slices).

---

### Task 4: Author Publication-Grade Validation Documentation

**Files:**
- Create: `docs/validation/multimodal-multislice-validation.md`
- Create: `benchmarks/evidence/multimodal-slices-tier2.json`

**Interfaces:**
- Consumes: Live run metrics, slice coordinates, mutational spectra data, and comparison logs from Task 3.
- Produces: Complete, exhaustive documentation report covering biological rationale, genomic architecture, mutational profiles, code path mapping, and 0-diff parity scoreboards.

- [ ] **Step 1: Write `docs/validation/multimodal-multislice-validation.md`**

Sections required:
1. Executive Summary & Prime Directive Parity (18/18 MATCH, 0 diffs).
2. Biological & Mutational Diversity Architecture (UV signatures, somatic LOH, amplicon primer trimming, deep ctDNA, RNA splicing, ancestral polymorphisms).
3. Complete 18-Slice Genomic Table (Dataset, Slice ID, Chromosome, Coordinates, Span, Reference Assembly, Targets, Code Paths Stressed).
4. Code Path Mapping in `vardict-core` (modules, functions, lines exercised).
5. Dual-Engine Parity Scoreboards (Java vs Rust record counts, variant counts, filter counts, diffs).
6. Multi-Interval BED Scheduling & Scalability Analysis ($O(1)$ memory consumption, thread parallelism).
7. Reproducibility & Audit Trail (Exact streaming and running commands).

- [ ] **Step 2: Verify `benchmarks/evidence/multimodal-slices-tier2.json` syntax and release fingerprint**

Run: `python3 -m json.tool benchmarks/evidence/multimodal-slices-tier2.json > /dev/null`
Expected: Valid JSON with `fingerprint: "7e9b1be766f73b17"` and `slices_passed: 18`.

---

### Task 5: Document Integration, Gate Verification, and Commit

**Files:**
- Modify: `docs/validation/omnibus-empirical-validation.md`
- Modify: `README.md`
- Modify: `docs/release/v1.0.0-pr.md`
- Modify: `release_plan.md`
- Modify: `PLAN.md`

**Interfaces:**
- Consumes: Task 1-4 deliverables.
- Produces: Integrated documentation suite, passed acceptance gate (8/8 PASS), passed PR checks (6/6 PASS), invariant release engine fingerprint `7e9b1be766f73b17`.

- [ ] **Step 1: Update documentation files to incorporate Task 24h results**
- [ ] **Step 2: Verify release engine fingerprint invariance**

Run: `bash harness/engine_fingerprint.sh`
Expected: `7e9b1be766f73b17`

- [ ] **Step 3: Run master acceptance gate and PR checks**

Run: `bash harness/acceptance.sh && bash harness/pr_check.sh`
Expected: 8/8 PASS and 6/6 PR stages PASS.

- [ ] **Step 4: Commit changes cleanly**

```bash
git add harness/slices/multimodal/ harness/get_multimodal_slices.sh harness/run_multimodal_slices.sh benchmarks/evidence/multimodal-slices-tier2.json docs/validation/multimodal-multislice-validation.md docs/validation/omnibus-empirical-validation.md README.md docs/release/v1.0.0-pr.md release_plan.md docs/superpowers/plans/
git commit -m "T024h: verify and document multi-slice multi-modal public dataset expansion for v1.0.0"
```
