<a id="readme"></a>
# VarDictRust (`vardictrust`)

[![CI](https://github.com/palakpsheth/VarDictRust/actions/workflows/ci.yml/badge.svg)](https://github.com/palakpsheth/VarDictRust/actions)
[![Documentation Status](https://img.shields.io/readthedocs/vardictrust/latest?label=docs)](https://vardictrust.readthedocs.io/en/latest/)
[![Release](https://img.shields.io/github/v/release/palakpsheth/VarDictRust?color=blue)](https://github.com/palakpsheth/VarDictRust/releases)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/vardictrust.svg?color=green)](https://anaconda.org/bioconda/vardictrust)
[![rewrites.bio](https://img.shields.io/badge/rewrites.bio-manifesto-blueviolet.svg)](https://rewrites.bio)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Security Policy](https://img.shields.io/badge/security-policy-blue.svg)](.github/SECURITY.md)

**VarDictRust** is an ultra-fast, memory-safe variant caller for next-generation sequencing (NGS) that delivers byte-identical results to [VarDictJava](https://github.com/AstraZeneca-NGS/VarDictJava) v1.8.3 while running **4.3×–6.9× faster on WGS/WES** (up to **36.4× faster** on amplicons) in **7.7×–10.8× less memory**. It consolidates the Java → R → Perl pipeline into a **single static binary**, eliminating JVM crashes, WGS out-of-memory errors, and runtime dependency friction.

- **Exact Parity:** Bitwise identical variant calls, `FILTER` tags, and somatic `STATUS` decisions to VarDictJava v1.8.3 across all cohorts ($\varepsilon = 10^{-6}$, 0 diffs).
- **High-Throughput Speed:** **4.3×–6.9× faster** on WGS/WES, **36.4× faster** on amplicon panels, and **3.58× parallel scaling** on 8 threads across full chromosomes.
- **Bounded $O(1)$ Memory:** Flat RSS scaling ($1.20\times$ from 10 Mbp to 249 Mbp on Chr1, `1.16` $\to$ `1.39 GiB`), using **7.7×–10.8× less RAM** than VarDictJava on WGS/WES.
- **Zero Runtime Dependencies:** Single static binary (`x86_64-unknown-linux-musl`). No JVM, R, Perl, or Conda needed.
- **Robust by Construction:** `#![forbid(unsafe_code)]`. Malformed CIGARs, missing indexes, and corrupt records return actionable errors instead of crashing.

```bash
# Somatic paired calling (tumor vs normal) — emits standard VCFv4.2 to stdout
vardictrust -G GRCh38.fa -f 0.01 -N "TUMOR|NORMAL" -b "tumor.bam|normal.bam" -R chr7:55000000-55300000 > somatic.vcf

# Output format: standard VCFv4.2 with VarDict somatic annotations
# CHROM  POS       ID  REF  ALT  QUAL  FILTER  INFO                                     FORMAT                 TUMOR                NORMAL
# chr7   55152340  .   G    A    184   PASS    TYPE=SNV;AF=0.245;STATUS=StrongSomatic  GT:DP:VD:AD:AF:RD:ALD  0/1:120:29:91,29...  0/0:115:0:115,0...
```

Built under the [rewrites.bio](https://rewrites.bio) manifesto: *a faster, more robust tool that produces the same results*. For biological provenance, mathematical proofs, and cohort datasets, see the [Omnibus Empirical Validation Report](docs/validation/omnibus-empirical-validation.md).

---

## Table of Contents

- [Installation](#installation)
- [Usage & Pipeline Recipes](#usage--pipeline-recipes)
- [Operational Modalities](#operational-modalities)
- [Performance & Scaling Scoreboard](#performance--scaling-scoreboard)
- [Empirical Parity & Validation Status](#empirical-parity--validation-status)
- [Architectural Parity & Upstream Policy](#architectural-parity--upstream-policy)
- [CLI Parameter Reference](#cli-parameter-reference)
- [AI Tooling & Verification Transparency](#ai-tooling--verification-transparency)
- [Attribution & Citation](#attribution--citation)
- [License](#license)

---

## Installation

VarDictRust is distributed as a standalone static executable with zero shared-library dependencies, as a Cargo crate, as a multi-arch container image, and via Bioconda.

### 1. Pre-Compiled Static Binary (Recommended for Linux x86_64)

Download the standalone `musl` binary from [GitHub Releases](https://github.com/palakpsheth/VarDictRust/releases/latest):

```bash
# Download release archive and SHA-256 checksum
curl -LO https://github.com/palakpsheth/VarDictRust/releases/latest/download/vardictrust-x86_64-unknown-linux-musl.tar.xz
curl -LO https://github.com/palakpsheth/VarDictRust/releases/latest/download/vardictrust-x86_64-unknown-linux-musl.tar.xz.sha256

# Verify SHA256 integrity
sha256sum --check vardictrust-x86_64-unknown-linux-musl.tar.xz.sha256

# Extract and install
tar -xf vardictrust-x86_64-unknown-linux-musl.tar.xz
sudo mv vardictrust /usr/local/bin/
```

### 2. Install via Cargo

```bash
# Install from crates.io
cargo install vardictrust

# Or build the latest release directly from source
cargo install --git https://github.com/palakpsheth/VarDictRust.git --tag v1.0.1
```

### 3. Containerized Execution (Docker, Podman, Apptainer)

```bash
# Docker / Podman
docker run --rm -v $(pwd):/data ghcr.io/palakpsheth/vardictrust:1.0.1 \
  vardictrust -G /data/ref.fa -b /data/sample.bam -R chr1:1000-2000 > output.vcf

# Apptainer / Singularity (HPC Clusters)
apptainer exec docker://ghcr.io/palakpsheth/vardictrust:1.0.1 vardictrust --help
```

### 4. Conda / Bioconda

Install via Bioconda (`linux-64`, `linux-aarch64`, `osx-arm64`, `osx-64`):

```bash
conda install -c bioconda -c conda-forge vardictrust
```

### 10-Second Installation Verification

Verify that BAM decoding, CIGAR realignment, and VCF emission function correctly:

```bash
# 1. Verify build and version metadata
vardictrust --version

# 2. Functional smoke test (from within a repository checkout; emits single synthetic SNV in <0.1s)
vardictrust -G fixtures/micro/micro.ref.fa -b fixtures/micro/micro.bam -R chr1:40-260
```

---

<a id="usage-recipes"></a>
<a id="usage--pipeline-recipes"></a>
## Usage & Pipeline Recipes

> [!WARNING]
> **The Silent BED Coordinate Trap (`-c 1 -S 2 -E 3`):**
> - **Standard 3/4-Column BED Files (`chr start end [name]`):** VarDict's historical default column indices are `-c 3 -S 7 -E 8` (1-based), tailored for legacy 8-column amplicon BED files. When supplying standard 3-column BED files, you **must explicitly pass `-c 1 -S 2 -E 3`** (and `-z` if coordinates are 0-based half-open) to avoid parser aborts.
> - **Direct Region Specification:** Alternatively, bypass BED parsing entirely with `-R chr:start-end[:gene]`.

### 1. Germline Single-Sample Calling

Constitutional variant discovery from single-sample WGS, WES, or targeted panels:

```bash
vardictrust -G GRCh38.fa -f 0.01 -N SAMPLE -b sample.bam -c 1 -S 2 -E 3 -z targets.bed > germline.vcf
```

### 2. Somatic Paired Calling (Tumor vs Normal)

Pair-matched somatic mutation discovery with loss-of-heterozygosity (LOH) detection and somatic classification (`StrongSomatic`, `LikelySomatic`, `Germline`, `SampleSpecific`):

```bash
vardictrust -G GRCh38.fa -f 0.01 -N "TUMOR|NORMAL" -b "tumor.bam|normal.bam" -c 1 -S 2 -E 3 -z target.bed > somatic.vcf
```

### 3. Targeted Amplicon Panels (8-Column BED & Primer Trimming)

Trims primer bias artifacts (`AMPBIAS`) using overlapping primer coordinates from an 8-column BED file:

```bash
vardictrust -G GRCh38.fa -f 0.01 -N SAMPLE -a 10:0.95 -b amplicon.bam amplicon_8col.bed > amplicon.vcf
```

### 4. RNA-Seq Splicing Analysis

Quantifies splice junction coverage and outputs a 4-column intron count matrix (`-i` / `--splice`):

```bash
vardictrust -G GRCh38.fa -i -b rna.bam -c 1 -S 2 -E 3 -z junctions.bed > splice_counts.txt
```

### 5. High-Throughput Streaming UNIX Pipeline

Stream directly from indexed BAM/CRAM files into compressed and indexed VCF archives:

```bash
vardictrust -G GRCh38.fa -f 0.01 -N SAMPLE -b sample.bam -R chr1:1000000-2000000 -th 8 \
  | bcftools filter -i 'FILTER="PASS"' -Oz -o filtered.vcf.gz \
  && tabix -p vcf filtered.vcf.gz
```

---

## Operational Modalities

| Modality | Engine Mode & Flag Signature | Input Assays | Output Schema | Typical Use Case |
| :--- | :--- | :--- | :--- | :--- |
| **Germline Calling** | `SimpleMode` (`-b sample.bam`) | WGS, WES, Targeted Panels | `VCFv4.2` (`10` cols) or `.var` (`36` / `38` cols with `--fisher`) | Constitutional SNV, InDel, and MNV discovery; family trios; population cohorts. |
| **Somatic Paired Calling** | `SomaticMode` (`-b "tumor\|normal"`) | Paired Tumor / Normal WGS & WES | `VCFv4.2` (`11` cols) or `.var` (`55` / `61` cols with `--fisher`) | Cancer somatic mutations down to 1% AF, LOH detection, `STATUS` classification. |
| **Targeted Amplicon Panel** | `AmpliconMode` (`-a 10:0.95`) | AmpliSeq, HaloPlex, PCR Panels | `VCFv4.2` (`10` cols) or `.var` (`38` / `40` cols with `--fisher`) | Overlapping amplicon panels; trims primer bias ([`AMPBIAS`](DIVERGENCES.md#amplicon-primer-trimming)) via 8-column BED. |
| **RNA-Seq Splicing** | `SplicingMode` (`-i` / `--splice`) | Spliced RNA-Seq BAM / CRAM | `4` columns (`qstart-qend`, `cnt`, `mean_pos`, `strand`) | Exon-intron splice junction quantification and split-read intron coverage matrix. |
| **Structural Variants & CRISPR** | `SimpleMode` / `SomaticMode` (`-L`, `-J`, `-j`) | Short-Read WGS, WES, CRISPR Assays | `VCFv4.2` (`<DEL>`, `<DUP>`, `<INV>`, `<INS>`) + `CRISPR` tag | Large split-read/discordant-mate SV discovery (disable via `-U`) and CRISPR cut-site (`-J`/`-j`) scoring. |
| **Pileup-Only Inspection** | `SimpleMode` / `AmpliconMode` (`-p` / `--pileup`) | Any Indexed BAM / CRAM | `VCFv4.2` or `.var` across every reference locus | Unfiltered per-base depth, reference quality, and non-variant pileup auditing. |

> [!NOTE]
> **Assumptions & Coordinate Conventions:**
> - **Coordinate Contracts:** BAM/CRAM alignment coordinates are parsed as 0-based half-open intervals (`[start, end)`), while VCF output records are strictly 1-based closed intervals (`[start, end]`). Passing `-z` (or `--zero-based`) shifts 0-based BED intervals by `+1`.
> - **Statistical Rounding:** All variant frequencies (`AF`), quality scores (`QUAL`), and p-values adhere to IEEE 754 round-to-nearest, ties-to-even (**Banker's rounding**), replicating Java's `RoundingMode.HALF_EVEN` and R's `round()` ([`round_half_even`](DIVERGENCES.md#statistical-rounding--floating-point-half-even)).
> - **Multi-Allelic Collapsing:** Standard single-sample mode defaults to single-winner top-AF candidate selection per locus matching [`var2vcf_valid.pl`](DIVERGENCES.md#consolidation-of-the-three-language-pipeline).

---

<a id="performance--scaling-scoreboard"></a>
<a id="performance-scaling-scoreboard"></a>
## Performance & Scaling Scoreboard

- **Exact Parity:** 0 variant discrepancies and identical `FILTER`/`STATUS` decisions against VarDictJava v1.8.3 across all 17 dual-engine differential benchmark runs and 92 integration fixtures ($\varepsilon = 10^{-6}$).
- **High-Throughput Speedup:** **4.3×–6.9× faster** on deep WGS/WES cohorts, **36.4× faster** on targeted amplicon assays, and **3.58× parallel speedup** on 8 threads across whole human chromosomes.
- **Bounded $O(1)$ Memory:** Flat resident memory scaling ($1.20\times$ peak RSS ratio across a $25\times$ increase in region span on Chromosome 1, `1.16 GiB` $\to$ `1.39 GiB`), using **7.7×–10.8× less memory** than VarDictJava on WGS/WES.
- **Self-Contained Deployment:** Single statically-linked binary (`x86_64-unknown-linux-musl`) with zero runtime dependencies.

Evaluated on Intel i7-13700H, Linux x86_64, 8 worker threads (`-th 8`). Per-dataset rows report current single-run certified evidence (`hg008-tier2` and `seqc2-exome-tier2` at `-th 4`²; exact commands in evidence). The Whole Chromosome 1 Scale row reports the 3-run median certified [`scale.json`](benchmarks/evidence/scale.json) evidence:

| Workload / Benchmark | Threads | VarDictJava Wall Time | vardictcpp Wall Time¹ | VarDictRust Wall Time | Speedup (vs Java) | VarDictJava Peak RSS | vardictcpp Peak RSS¹ | VarDictRust Peak RSS | Memory Use (vs Java) | Evidence Artifact |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :---: |
| **Targeted Amplicon Panel** | `8` (`-th 8`) | `0:02.55` | `0:03.75` ([`CPP-3`](docs/upstream-bugs.md#cpp-3)) | **`0:00.07`** | **36.4× faster** | `122.2 MiB` | `88.9 MiB` | **`58.7 MiB`** | **2.1× less** | [`amplicon-tier2.json`](benchmarks/evidence/amplicon-tier2.json) |
| **GIAB HG002 (chr20:10M–30M, 30×)** | `8` (`-th 8`) | `1:31.86` | `0:18.87` ([`CPP-9`](docs/upstream-bugs.md#cpp-9)) | **`0:20.82`** | **4.4× faster** | `4.20 GiB` | `853.1 MiB` | **`456.8 MiB`** | **9.4× less** | [`hg002-tier2.json`](benchmarks/evidence/hg002-tier2.json) |
| **SEQC2 HCC1395 Somatic Exome** | `4` (`-th 4`²) | `53:53.76` | `7:11.99` ([`CPP-9`](docs/upstream-bugs.md#cpp-9)) | **`7:46.34`** | **6.9× faster** | `3.52 GiB` | `2.15 GiB` | **`467.9 MiB`** | **7.7× less** | [`seqc2-exome-tier2.json`](benchmarks/evidence/seqc2-exome-tier2.json) |
| **GIAB HG008 Somatic PDAC** | `4` (`-th 4`²) | `1:27.68` | `0:21.96` ([`CPP-9`](docs/upstream-bugs.md#cpp-9)) | **`0:20.44`** | **4.3× faster** | `3.29 GiB` | `303.6 MiB` | **`313.1 MiB`** | **10.8× less** | [`hg008-tier2.json`](benchmarks/evidence/hg008-tier2.json) |
| **Whole Chromosome 1 Scale (249 Mbp)** | `1 → 8` (`-th 8`³) | `15:24.38` (`-th 1`³) | `2:48.27` ([`CPP-5`](docs/upstream-bugs.md#cpp-5)) | **`4:17.87`** | **3.58× parallel** | `>31.0 GiB` (`OOM`³) | `778.6 MiB` ([`CPP-5`](docs/upstream-bugs.md#cpp-5)) | **`1.39 GiB`** | **1.20× flat RSS**⁴ | [`scale.json`](benchmarks/evidence/scale.json) |

¹ **`vardictcpp` comparison (`vardictcpp 3 hc52dbad_0`, bioconda).** Runs alongside each benchmark under identical thread and timer constraints (`LC_ALL=C /usr/bin/time -v`). Release acceptance gates strictly on `VarDictRust == VarDictJava v1.8.3`; linked `CPP-*` codes document where `vardictcpp` diverges from `VarDictJava` (e.g., [`CPP-3`](docs/upstream-bugs.md#cpp-3) TSV-only `-a`, [`CPP-5`](docs/upstream-bugs.md#cpp-5) `22.5 GiB` unchunked `-R` RSS vs `778.6 MiB` with `--chunk 50000`, [`CPP-9`](docs/upstream-bugs.md#cpp-9) caller variant deltas).
² `hg008-tier2` and `seqc2-exome-tier2` run at `-th 4` because `VarDictJava` at `-th 8` is OOM-killed on 31 GiB hosts (bisected `1✓ / 4✓ / 8✗`); exact commands in evidence.
³ On `Whole Chromosome 1 Scale`, `VarDictJava` at `-th 8` OOMs on the 31 GiB i7 host (completed in `57:16.00` on an 8-core Ryzen host; see [`docs/release/v1.0.0-pr.md`](docs/release/v1.0.0-pr.md)); the `15:24.38` (`924.38 s`) → **`4:17.87`** (`257.87 s`) wall-time columns report the 3-run median 1-core vs 8-core `VarDictRust` parallel scaling on `chr1`.
⁴ **Flat $O(1)$ memory ratio:** `1.39 GiB` (`1,459,504 KB`) across all 249 Mbp of `chr1` vs `1.16 GiB` (`1,219,156 KB`) on a 10 Mbp slice (**1.20× RSS ratio** across a 25× span increase; pre-registered gate $\le 1.50\times$).

> [!NOTE]
> **Bounded $O(1)$ Memory at Any Depth:** VarDictJava holds entire regional chromosome structures in JVM heap collections, leading to JVM heap exhaustion and Out-Of-Memory errors on deep WGS ([Issues #197 and #384](https://github.com/AstraZeneca-NGS/VarDictJava/issues/197)). VarDictRust processes genomic regions via **bounded tiles with halo intervals** using RAII streaming. Across a 25-fold increase in genomic span (10 Mbp to all 249 Mbp of chromosome 1), VarDictRust peak RSS scales by only **$1.20\times$** (`1.39 GiB` vs `1.16 GiB` = **1.20× RSS ratio**, certified gate $\le 1.50\times$).

---

<a id="empirical-parity--validation-status"></a>
<a id="empirical-parity-validation-status"></a>
## Empirical Parity & Validation Status

Every release candidate must pass dual-engine differential testing against **VarDictJava v1.8.3** (commit `0886261`, Perl submodule `009e017`) and match physical biological truth across gold-standard reference benchmarks:

### 1. Three-Way Parity Matrix (VarDictJava vs VarDictRust vs vardictcpp)

All comparisons evaluate the full end-to-end VCF output under float tolerance $\varepsilon = 10^{-6}$:

| Dataset / Benchmark | Modality | Target Region / Assay | Verification Tier | VarDictJava Variants | VarDictRust Variants | vardictcpp Variants¹ | Rust vs Java Parity | vardictcpp vs Java Parity¹ |
| :--- | :--- | :--- | :---: | :--- | :--- | :--- | :--- | :--- |
| **Upstream Integration Corpus** | `Simple`, `Somatic`, `Amplicon` | 92 micro-fixtures | **[Tier-1](docs/testing/test-matrix-and-ci-audit.md#tier-1)** | `91 PASS / 1 XFAIL*` | **`91 PASS / 1 XFAIL*`** | `34 PASS / 32 MISMATCH / 26 SKIP` | **MATCH** (0 live diffs) | [`MISMATCH (CPP-1)`](docs/upstream-bugs.md#cpp-1), [`SKIP (CPP-4)`](docs/upstream-bugs.md#cpp-4) |
| **GIAB HG002 Germline WGS** | `SimpleMode` (Germline) | `chr20:10M–30M` (30×) | **[Tier-2](docs/testing/test-matrix-and-ci-audit.md#tier-2)** | `227,562` | **`227,562`** | `227,560` | **MATCH** (0 diffs, $\varepsilon = 10^{-6}$) | [`MISMATCH (CPP-9)`](docs/upstream-bugs.md#cpp-9), [`CPP-1`](docs/upstream-bugs.md#cpp-1) |
| **SEQC2 HCC1395/BL Exome** | `SomaticMode` (`tumor\|normal`) | Whole Exome (WES) | **[Tier-2](docs/testing/test-matrix-and-ci-audit.md#tier-2)** | `623,444` | **`623,444`** | `623,442` | **MATCH** (0 diffs, $\varepsilon = 10^{-6}$) | [`MISMATCH (CPP-9)`](docs/upstream-bugs.md#cpp-9), [`CPP-1`](docs/upstream-bugs.md#cpp-1) |
| **SEQC2 Somatic Slices** | `SomaticMode` (`tumor\|normal`) | `chr1, 2, 3, 17, 19, 21` slices | **[Tier-2](docs/testing/test-matrix-and-ci-audit.md#tier-2)** | `4,151` | **`4,151`** | `4,151` | **MATCH** (0 diffs in SNV/InDel/SV) | [`MISMATCH (CPP-1)`](docs/upstream-bugs.md#cpp-1), [`CPP-9`](docs/upstream-bugs.md#cpp-9) |
| **GIAB HG008 Somatic PDAC** | `SomaticMode` (`tumor\|normal`) | PDAC high-depth slices | **[Tier-2](docs/testing/test-matrix-and-ci-audit.md#tier-2)** | `18,426` | **`18,426`** | `18,425` | **MATCH** (0 diffs in somatic status) | [`MISMATCH (CPP-9)`](docs/upstream-bugs.md#cpp-9), [`CPP-1`](docs/upstream-bugs.md#cpp-1) |
| **GIAB HG008 CRAM Twin** | `SomaticMode` (CRAM vs BAM) | PDAC CRAM reference decode | **[Tier-2](docs/testing/test-matrix-and-ci-audit.md#tier-2)** | `18,426` | **`18,426`** | `18,425` | **MATCH** (0 diffs, BAM == CRAM) | [`MISMATCH (CPP-9)`](docs/upstream-bugs.md#cpp-9), [`CPP-1`](docs/upstream-bugs.md#cpp-1) |
| **Targeted Amplicon Panel** | `AmpliconMode` (`-a 10:0.95`) | Overlapping primer assay | **[Tier-2](docs/testing/test-matrix-and-ci-audit.md#tier-2)** | `2` | **`2`** | `N/A` (TSV-only) | **MATCH** (0 diffs, `AMPBIAS`) | [`SKIP (CPP-3)`](docs/upstream-bugs.md#cpp-3) |
| **Chr1 Scale Benchmark** | `SimpleMode` (Germline) | `chr1:1–248,956,422` (249 Mbp) | **[Level 5 Scale](docs/testing/test-matrix-and-ci-audit.md#level-5)** | `2,840,119` | **`2,840,119`** | `2,840,119` | **MATCH** (0 diffs across `chr1`) | [`MISMATCH (CPP-5)`](docs/upstream-bugs.md#cpp-5), [`CPP-1`](docs/upstream-bugs.md#cpp-1) |
| **Extended Modalities Suite** | `Splicing`, `<INV>`, `<DUP>`, CRISPR, Flags | 17 biological assays | **[Tier-2](docs/testing/test-matrix-and-ci-audit.md#tier-2)** | `4,993` (`17/17 PASS`) | **`4,993` (`17/17 PASS`)** | `2,304` (`0/17 PASS`) | **MATCH** (0 diffs in `.var` & VCF) | [`SKIP (CPP-3)`](docs/upstream-bugs.md#cpp-3), [`MISMATCH (CPP-1)`](docs/upstream-bugs.md#cpp-1), [`CPP-2`](docs/upstream-bugs.md#cpp-2), [`CPP-7`](docs/upstream-bugs.md#cpp-7), [`CPP-8`](docs/upstream-bugs.md#cpp-8), [`CPP-10`](docs/upstream-bugs.md#cpp-10) |
| **Multi-Modal Multi-Slice Suite** | `Germline`, `Somatic`, `RNA` | 18 slices across 6 cohorts & 9 chrs | **[Tier-2](docs/testing/test-matrix-and-ci-audit.md#tier-2)** | `1,129` (`18/18 PASS`) | **`1,129` (`18/18 PASS`)** | `1,103` (`0/18 PASS`) | **MATCH** (0 diffs, GRCh38 & b37) | `SKIP (-i)`, [`MISMATCH (CPP-1)`](docs/upstream-bugs.md#cpp-1), [`CPP-9`](docs/upstream-bugs.md#cpp-9) |

\* *The 91 PASS / 1 XFAIL metric reflects comparison against stale pre-computed static files in the upstream repository (`VarDictJava/testdata/integrationtestcases/`). For the single static mismatch (`hard_clip_case.fa`), the upstream fixture encodes conflicting reference chunks ([Class D / UB-4](docs/upstream-bugs.md#ub-4)). When evaluated against the live `VarDictJava v1.8.3` oracle, VarDictRust produces 100% byte-identical output (7/7 table rows, 5/5 VCF variants, 0 diffs). Thus, live oracle parity is 100% across all 92 fixtures.*
¹ **`vardictcpp` third-arm comparison (`vardictcpp 3 hc52dbad_0`, bioconda).** Executed alongside each Tier-1/Tier-2 comparator arm; the release gate remains strictly `VarDictRust == VarDictJava v1.8.3`. Linked `CPP-*` identifiers reference the root-cause divergence catalog in [`docs/upstream-bugs.md`](docs/upstream-bugs.md#cpp-divergences).

- **Intermediate Table Parity ([Tier-1](docs/testing/test-matrix-and-ci-audit.md#tier-1)):** Emitted intermediate `.var` tables match 100% of upstream columns and ordering. Slashes denote standard mode vs. statistical `--fisher` mode (Simple: 36/38 columns; Amplicon: 38/40 columns; Somatic: 55/61 columns). Standard CLI output is standard VCFv4.2 (10 or 11 columns); Splicing Mode (`-i`) outputs a 4-column intron count table.

<a id="2-ground-truth-concordance-tier-3-via-rtg-vcfeval-vardictcpp-advisory"></a>
### 2. Ground Truth Concordance ([Tier-3](docs/testing/test-matrix-and-ci-audit.md#tier-3) via RTG vcfeval)

Evaluated against physical ground truth using Real Time Genomics [`rtg vcfeval`](docs/validation/omnibus-empirical-validation.md#pillar-1):

| Benchmark Dataset | Truth Standard | Assay / Threshold | Java & Rust Sites (TP / FP / FN) | vardictcpp Sites (TP / FP / FN)¹ | VarDictJava vcfeval | VarDictRust vcfeval | vardictcpp vcfeval¹ | Rust vs Java | vardictcpp Notes¹ |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :---: | :--- |
| **GIAB HG002** (Operating Point) | NIST v4.2.1 High-Confidence | Germline (`QUAL >= 128`) | `22,748 / 553 / 871` | `22,748 / 553 / 871` | `F1: 0.9699 (P: 0.9763, R: 0.9635)` | **`F1: 0.9699 (P: 0.9763, R: 0.9635)`** | `F1: 0.9699 (P: 0.9763, R: 0.9635)` | **EQUAL** | [`CPP-6`](docs/upstream-bugs.md#cpp-6) |
| **GIAB HG002** (Raw Baseline) | NIST v4.2.1 High-Confidence | Germline (`None` / raw) | `22,803 / 22,458 / 813` | `22,803 / 22,442 / 813` | `F1: 0.6622 (P: 0.5038, R: 0.9659)` | **`F1: 0.6622 (P: 0.5038, R: 0.9659)`** | `F1: 0.6624 (P: 0.5040, R: 0.9659)` | **EQUAL** | [`CPP-6`](docs/upstream-bugs.md#cpp-6), [`CPP-9`](docs/upstream-bugs.md#cpp-9) |
| **GIAB HG005** (Operating Point) | NIST v4.2.1 Han Chinese | Germline (`QUAL >= 226`) | `29 / 1 / 1` | `29 / 1 / 1` | `F1: 0.9667 (P: 0.9667, R: 0.9667)` | **`F1: 0.9667 (P: 0.9667, R: 0.9667)`** | `F1: 0.9667 (P: 0.9667, R: 0.9667)` | **EQUAL** | [`CPP-6`](docs/upstream-bugs.md#cpp-6) |
| **GIAB NA24385 b37** (Operating Point) | NIST v4.2.1 Ashkenazi Son | Germline b37 (`QUAL >= 160`) | `5 / 0 / 0` | `5 / 0 / 0` | `F1: 1.0000 (P: 1.0000, R: 1.0000)` | **`F1: 1.0000 (P: 1.0000, R: 1.0000)`** | `F1: 1.0000 (P: 1.0000, R: 1.0000)` | **EQUAL** | [`CPP-6`](docs/upstream-bugs.md#cpp-6), [`CPP-7`](docs/upstream-bugs.md#cpp-7) |
| **SEQC2 HCC1395** (Raw Baseline) | SEQC2 v1.2.1 Somatic Truth | Somatic (`None` / raw) | `1,056 / 11,334 / 140` | `1,056 / 11,333 / 140` | `F1: 0.1555 (P: 0.0852, R: 0.8842)` | **`F1: 0.1555 (P: 0.0852, R: 0.8842)`** | `F1: 0.1555 (P: 0.0852, R: 0.8842)` | **EQUAL** | [`CPP-6`](docs/upstream-bugs.md#cpp-6), [`CPP-9`](docs/upstream-bugs.md#cpp-9) |
| **GIAB HG008** (Raw Baseline) | NIST Draft Benchmark v0.3 | Somatic (`None` / raw) | `575 / 1,565 / 288` | `575 / 1,563 / 288` | `F1: 0.3856 (P: 0.2687, R: 0.6825)` | **`F1: 0.3856 (P: 0.2687, R: 0.6825)`** | `F1: 0.3858 (P: 0.2689, R: 0.6825)` | **EQUAL** | [`CPP-6`](docs/upstream-bugs.md#cpp-6), [`CPP-9`](docs/upstream-bugs.md#cpp-9) |

¹ **`vardictcpp` vcfeval normalization (`vardictcpp 3 hc52dbad_0`).** Because native `vardictcpp -vcf` omits `##contig` header lines ([`CPP-6`](docs/upstream-bugs.md#cpp-6)) and strips `chr` on ambiguous b37 FASTAs ([`CPP-7`](docs/upstream-bugs.md#cpp-7)), the concordance harness injects reference contigs before running `rtg vcfeval`. Residual 1–16 FP count differences on raw HG002, SEQC2, and HG008 reflect upstream [`CPP-9`](docs/upstream-bugs.md#cpp-9) MNP/SV caller deltas; `VarDictRust` and `VarDictJava` produce identical TP, FP, FN, Precision, Recall, and $F_1$ across every threshold.

*For full scientific details, reference checksums, and exact reproduction commands across all 10 public cohorts and 17 operational modalities, see the [`Omnibus Empirical Validation Report`](docs/validation/omnibus-empirical-validation.md). For the complete Seven-Level Testing Matrix and CI/CD audit, see [`docs/testing/test-matrix-and-ci-audit.md`](docs/testing/test-matrix-and-ci-audit.md).*

---

<a id="architectural-parity--upstream-policy"></a>
<a id="architectural-parity-upstream-policy"></a>
## Architectural Parity & Upstream Policy

VarDictRust is built under strict compliance with the **[rewrites.bio](https://rewrites.bio) Prime Directive**:

1. **Fidelity over Innovation:** VarDictRust v1.0.x matches **every** output decision of VarDictJava v1.8.3 on valid inputs. This includes IEEE 754 half-even banker's rounding ([`round_half_even`](DIVERGENCES.md#statistical-rounding--floating-point-half-even)), exact odds ratios, OpenJDK hash iteration layout via the deterministic [`JavaHashMap`](DIVERGENCES.md#ub-3-var-row-order-across-positions-follows-hashmap-iteration-order-282-family) engine, and in-engine somatic status determination ([`SomaticPostProcessModule`](DIVERGENCES.md#somatic-status-determination)).
2. **Robustness by Construction:** While VarDictJava crashes on malformed CIGAR strings ([`modifyCigar` issue #341 family](docs/upstream-bugs.md) and [PR #402](DIVERGENCES.md#category-1-cigar-and-realignment-resilience)), null pointers, or missing indices ([Issue #187](DIVERGENCES.md#category-7-cli-diagnostics-and-error-messages)), VarDictRust enforces `#![forbid(unsafe_code)]` and returns actionable `Result<_, VardictError>` errors without panicking.
3. **Bounded Memory at Any Depth:** While VarDictJava exhausts JVM heap memory on deep WGS ([Issues #197 and #384](https://github.com/AstraZeneca-NGS/VarDictJava/issues/197)), VarDictRust processes reads using bounded genomic tiles with halo intervals, maintaining flat resident memory ($O(1)$ RSS scaling).
4. **Bug Release Policy:** Upstream VarDictJava output quirks are faithfully replicated in v1.0.x and tracked in [`docs/upstream-bugs.md`](docs/upstream-bugs.md) (including [UB-1 through UB-7](docs/upstream-bugs.md)). Corrective fixes will ship in v1.1.0+ behind explicit opt-in flags (`--fix <id>`) per [`ROADMAP.md`](ROADMAP.md).
5. **Documented Quirks & Testing Matrix:** All design choices and parity validations are detailed in [`DIVERGENCES.md`](DIVERGENCES.md) and [`docs/testing/test-matrix-and-ci-audit.md`](docs/testing/test-matrix-and-ci-audit.md).

---

## CLI Parameter Reference

Flag parity with VarDictJava v1.8.3 is defined against [`CmdParser.java`](VarDictJava/src/main/java/com/astrazeneca/vardict/CmdParser.java) and verified by differential testing:

*Evidence Status Vocabulary:*
- `corpus` — Verified against integration corpus test cases or [Tier-2](docs/testing/test-matrix-and-ci-audit.md#tier-2) oracle sweeps.
- `unexercised` — Implemented to Java source specification, but lacks explicit differential corpus fixtures.
- `no-output` — Performance or concurrency flag; does not alter output bytes.
- `parity-exempt` — Accepted for CLI compatibility matching upstream semantics.

### 1. Core Inputs & Outputs

| Short | Long | Description | Default | Status |
| :--- | :--- | :--- | :--- | :--- |
| `-G` | `--reference` | Reference FASTA (`.fai` required; hg19/hg38/mm10 shortcuts accepted) | required | `corpus` |
| `-b` | `--bam` | Indexed BAM/CRAM; somatic: `"tumor.bam\|normal.bam"` | required | `corpus` |
| `-N` | `--name` | Sample name (somatic: `"TUMOR\|NORMAL"`) | derived from BAM | `corpus` |
| `-R` | `--region` | Region `chr:start-end[:gene]` instead of BED | — | `corpus` |
| positional | — | BED file (used when `-R` absent) | — | `corpus` |
| `-h` | `--header` | Print `.var` column header before results | off | `corpus` |
| — | `--emit-intermediate` | Emit legacy `.var` table instead of VCF | off | `corpus` |
| `-d` | — | Delimiter for splitting region info | `\t` | `corpus` |
| `-v` | — | Accepted and ignored (never read by VarDictJava v1.8.3) | — | `parity-exempt` |
| `-y` | `--verbose` | Verbose mode | off | `no-output` |
| `-D` | `--debug` | Debug mode (append per-read debug columns to `.var`) | off | `corpus` |
| — | `--DP` | Debug print destination (`OUT` or `ERR`; alias: `--defaultPrinter`) | — | `corpus` |

### 2. Coordinate & BED Column Mapping

| Short | Long | Description | Default | Status |
| :--- | :--- | :--- | :--- | :--- |
| `-c` | — | BED chromosome column (1-based; default 3 for 8-col amplicon BED) | `3` | `corpus` |
| `-S` | — | BED start column (default 7 for 8-col amplicon BED) | `7` | `corpus` |
| `-E` | — | BED end column (default 8 for 8-col amplicon BED) | `8` | `corpus` |
| `-s` | — | BED segment-start column | `10` (= `-S` if set) | `corpus` |
| `-e` | — | BED segment-end column | `11` (= `-E` if set) | `corpus` |
| `-g` | — | BED gene column | `13` | `corpus` |
| `-z` | — | BED coordinates are 0-based half-open | off (1-based) | `corpus` |
| `-x` | — | BED region extension (bp) | `0` | `corpus` |
| `-X` | — | Nucleotide extension | `2` | `corpus` |
| `-Y` | — | Max region length | `1200` | `corpus` |

### 3. Variant Calling & Quality Thresholds

| Short | Long | Description | Default | Status |
| :--- | :--- | :--- | :--- | :--- |
| `-f` | `--af-threshold` | Allele-frequency threshold | `0.01` | `corpus` |
| `-r` | — | Minimum variant reads | `2` | `corpus` |
| `-q` | — | Minimum base quality | `22.5` | `corpus` |
| `-Q` | — | Minimum mapping quality | — | `corpus` |
| `-P` | — | Mean read position filter | `5` | `corpus` |
| `-m` | — | Maximum mismatches | `8` | `corpus` |
| `-M` | — | Min match+insertion read length (shorter reads ignored) | `0` | `corpus` |
| `-B` | — | Minimum reads for bias | `2` | `corpus` |
| `-o` | — | qratio filter | `1.5` | `corpus` |
| `-O` | — | Output MAPQ threshold | `0` | `corpus` |
| `-F` | — | SAM flag filter | `0x504` | `corpus` |
| `-k` | — | Local realignment 0/1 | `1` | `corpus` |
| `-3` | — | Move indels to 3′ ([`moveIndelsTo3`](DIVERGENCES.md)) | off | `corpus` |
| `-K` | — | Include Ns in depth calculation | off | `corpus` |
| `-p` | — | Pileup-only mode | off | `corpus` |
| `-t` | `--dedup` | Remove duplicated reads (one pair per start position) | off | `corpus` |
| `-u` | — | Unique mode: count overlapping mates once (forward read) | off | `corpus` |
| — | `--UN` | Unique mode: count once using first read only | off | `corpus` |
| — | `--adaptor` | Adaptor sequence(s), comma-separated, excluded from realignment | — | `corpus` |

### 4. Modality Flags: Somatic & Amplicon

| Short | Long | Description | Default | Status |
| :--- | :--- | :--- | :--- | :--- |
| `-V` | — | Min normal AF (somatic mode) | `0.05` | `corpus` |
| — | `--fisher` | In-engine Fisher exact test (p-value / odds ratio) | off | `corpus` |
| `-a` | `--amplicon` | Amplicon mode `"int:float"` (distance:fraction) | auto `10:0.95` on 8-col BED | `corpus` |
| `-i` | `--splice` | RNA splicing mode (emits 4-column intron count table) | off | `corpus` |
| — | `--mfreq` | Mono frequency threshold | `0.25` | `corpus` |
| — | `--nmfreq` | Normal mono frequency threshold | `0.1` | `corpus` |
| `-J` | — | CRISPR cutting site coordinate | — | `corpus` |
| `-j` | — | CRISPR flanking size | — | `corpus` |

### 5. Structural Variants & Complex Indels

| Short | Long | Description | Default | Status |
| :--- | :--- | :--- | :--- | :--- |
| `-U` | `--nosv` | Disable structural variant calling | off | `corpus` |
| `-L` | — | SV minimum length in bp | `1000` | `corpus` |
| `-w` | — | SV read length | `300` | `corpus` |
| `-W` | — | SV flanking length | `100` | `corpus` |
| `-A` | — | SV minimum AF | `4` | `corpus` |
| `-I` | — | Indel size | `50` | `corpus` |
| — | `--chimeric` | Chimeric read handling | off | `corpus` |
| — | `--deldupvar` | Delete/dup variant handling | off | `corpus` |

### 6. Performance & Concurrency

| Short | Long | Description | Default | Status |
| :--- | :--- | :--- | :--- | :--- |
| — | `--threads` | Thread count (`-th`/`--th` aliases; bare `-th` = all cores) | all cores | `no-output` |
| — | `--tile-size` | Region tile size in bp (bounded-memory WGS) | auto | `no-output` |

*Reproducible Profiling Hooks:* When compiled with opt-in Cargo features (`profile-cpu`, `profile-heap`, or `profile-dhat`), `vardictrust` honors `VARDICT_PROFILE_CPU`, `VARDICT_PROFILE_HEAP`, and `VARDICT_PROFILE_DHAT` output directory environment variables (orchestrated via [`harness/profile.sh`](harness/README.md)).

---

<a id="ai-tooling--verification-transparency"></a>
<a id="ai-tooling-verification-transparency"></a>
## AI Tooling & Verification Transparency

Per the [rewrites.bio](https://rewrites.bio) manifesto, we disclose all AI tools, their operational roles, and human verification gates:

- **AI Systems Used:** Google DeepMind Gemini / Antigravity CLI, Anthropic Claude Code, and Qoder.
- **Scope of AI Assistance:** Implementation of core Rust modules under direct human architectural instruction, transliteration of legacy R/Perl statistical functions into pure Rust, test harness development, and adversarial differential verification.
- **Verification Guarantee (Evidence Before Completion):** **Zero** AI-generated code was accepted based on model assertion. Every single algorithmic component was validated against VarDictJava v1.8.3 across 275 automated Cargo workspace tests (227 core library tests, 2 amplicon streaming memory tests, 2 CRAM stack tests, 13 CLI tests, and 31 integration binary tests), the 92-case integration corpus, 10 slice fixtures, 18 multi-modal slices across 6 cohorts, and multi-gigabyte public benchmark datasets (SEQC2, GIAB HG008, GIAB HG002, GIAB HG005, GIAB NA24385 b37) under strict float tolerance $\varepsilon = 10^{-6}$, with complete push-button reproducibility documented in [`docs/testing/test-matrix-and-ci-audit.md`](docs/testing/test-matrix-and-ci-audit.md).

---

<a id="attribution--citation"></a>
<a id="attribution-citation"></a>
## Attribution & Citation

VarDictRust is authored and maintained by **Palak Sheth** (<palakpsheth@gmail.com>), developed as an idiomatic Rust rewrite of **VarDict** and **VarDictJava** originally created by Zhongwu Lai and the AstraZeneca-NGS team:

> Lai Z, Markovets A, Ahdesmaki M, et al. VarDict: a novel and versatile variant caller for next-generation sequencing in cancer research. *Nucleic Acids Research*, 2016;44(11):e108. doi:[10.1093/nar/gkw227](https://doi.org/10.1093/nar/gkw227)

```bibtex
@article{lai2016vardict,
  author    = {Lai, Zhongwu and Markovets, Alena and Ahdesmaki, Miika and Chapman, Brad and Hofmann, Oliver and McEwen, Robert and Johnson, Justin and Dougherty, Brian and Barrett, J. Carl and Dry, Jonathan R.},
  title     = {{VarDict: a novel and versatile variant caller for next-generation sequencing in cancer research}},
  journal   = {Nucleic Acids Research},
  volume    = {44},
  number    = {11},
  pages     = {e108},
  year      = {2016},
  doi       = {10.1093/nar/gkw227}
}

@article{seqc2_2021_reference,
  author    = {Fang, Li Tai and Zhu, Baolin and Zhao, Yongmei and Chen, Wanqiu and Yang, Zhaowei and Xu, Leming and others},
  title     = {{Establishing community reference samples, data and call sets for cancer genomics}},
  journal   = {Nature Biotechnology},
  volume    = {39},
  number    = {9},
  pages     = {1151--1160},
  year      = {2021},
  doi       = {10.1038/s41587-021-00993-6}
}

@article{giab_2020_benchmark,
  author    = {Zook, Justin M. and Hansen, Nancy F. and Olson, Nathan D. and Chapman, Lindsay and Mullikin, James C. and Xiao, Chunlin and others},
  title     = {{A robust benchmark for detection of germline large deletions and insertions}},
  journal   = {Nature Biotechnology},
  volume    = {38},
  number    = {11},
  pages     = {1347--1355},
  year      = {2020},
  doi       = {10.1038/s41587-020-0538-8}
}
```

If you use VarDictRust in published research, please cite both the original Lai et al. (2016) paper and this software:

```bibtex
@software{sheth2026vardictrust,
  author    = {Sheth, Palak},
  title     = {{VarDictRust: Memory-Safe, High-Performance Variant Calling for NGS}},
  year      = {2026},
  version   = {1.0.1},
  url       = {https://github.com/palakpsheth/VarDictRust}
}
```
Full citation metadata is available in [`CITATION.cff`](CITATION.cff).

---

## License

VarDictRust is open-source software licensed under the [MIT License](LICENSE). Third-party attribution, component licenses (including R Project `nmath` distribution notices), and copyright disclosures are cataloged in [`LICENSE-3RD-PARTY`](LICENSE-3RD-PARTY).
