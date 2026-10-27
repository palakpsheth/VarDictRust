<a id="readme"></a>
# VarDictRust (`vardictrust`)

[![CI](https://github.com/palakpsheth/VarDictRust/actions/workflows/ci.yml/badge.svg)](https://github.com/palakpsheth/VarDictRust/actions/workflows/ci.yml)
[![Differential Parity](https://github.com/palakpsheth/VarDictRust/actions/workflows/differential.yml/badge.svg)](https://github.com/palakpsheth/VarDictRust/actions/workflows/differential.yml)
[![Documentation Status](https://img.shields.io/readthedocs/vardictrust/latest?label=docs)](https://vardictrust.readthedocs.io/en/latest/)
[![Release](https://img.shields.io/github/v/release/palakpsheth/VarDictRust?color=blue)](https://github.com/palakpsheth/VarDictRust/releases)
[![Crates.io](https://img.shields.io/crates/v/vardictrust.svg?logo=rust)](https://crates.io/crates/vardictrust)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/vardictrust.svg?color=green)](https://anaconda.org/bioconda/vardictrust)
[![Nextflow nf-core](https://img.shields.io/badge/nf--core-vardictrust-23aa62?logo=nextflow&logoColor=white)](packaging/nf-core/modules/vardictrust/main.nf)
[![rewrites.bio](https://img.shields.io/badge/rewrites.bio-manifesto-blueviolet.svg)](https://rewrites.bio)
[![License: MIT AND GPL-2.0-or-later](https://img.shields.io/badge/License-MIT%20AND%20GPL--2.0--or--later-yellow.svg)](LICENSE)
[![Security Policy](https://img.shields.io/badge/security-policy-blue.svg)](.github/SECURITY.md)

**VarDictRust** is a fast, memory-safe NGS variant caller that delivers byte-identical results to [VarDictJava](https://github.com/AstraZeneca-NGS/VarDictJava) v1.8.3 while running **4.3×–6.9× faster on WGS/WES** (up to **36.4× faster** on amplicons) in **7.7×–10.8× less memory**. It consolidates the Java → R → Perl pipeline into a **single binary**, eliminating JVM crashes, WGS out-of-memory errors, and runtime dependency friction.

- **Exact Parity:** Bitwise identical variant calls, `FILTER` tags, and somatic `STATUS` decisions to VarDictJava v1.8.3 across all cohorts ($\varepsilon = 10^{-6}$, 0 diffs).
- **High-Throughput Speed:** **4.3×–6.9× faster** on WGS/WES, **36.4× faster** on amplicon panels, and **3.48× parallel scaling** on 8 threads across full chromosomes.
- **Bounded $O(1)$ Memory:** Flat RSS scaling ($1.22\times$ from 10 Mbp to 249 Mbp on Chr1, `1.16` $\to$ `1.42 GiB`), using **7.7×–10.8× less RAM** than VarDictJava on WGS/WES.
- **Zero External Runtimes:** Prebuilt static Linux (`x86_64`/`aarch64` `musl`) and native Linux (`x86_64` `gnu`) / macOS (`arm64`/`x86_64`) binaries—no JVM, R, or Perl.
- **Robust by Construction:** `#![forbid(unsafe_code)]`. Malformed CIGARs, missing indexes, and corrupt records return actionable errors instead of crashing.

```bash
# Install via Bioconda (or see static binaries, Cargo, Seqera Wave containers & Nextflow below)
conda install -c bioconda -c conda-forge vardictrust

# Somatic paired calling (tumor vs normal) — emits standard VCFv4.2 to stdout
vardictrust -G GRCh38.fa -f 0.01 -N "TUMOR|NORMAL" -b "tumor.bam|normal.bam" -R chr7:55000000-55300000 > somatic.vcf

# Output format: standard VCFv4.2 with VarDict somatic annotations
# CHROM  POS       ID  REF  ALT  QUAL  FILTER  INFO                                     FORMAT                 TUMOR                NORMAL
# chr7   55152340  .   G    A    184   PASS    TYPE=SNV;AF=0.245;STATUS=StrongSomatic  GT:DP:VD:AD:AF:RD:ALD  0/1:120:29:91,29...  0/0:115:0:115,0...
```

Built under the [rewrites.bio](https://rewrites.bio) manifesto: *a faster, more robust tool that produces the same results*. See the [Omnibus Empirical Validation Report](docs/validation/omnibus-empirical-validation.md) for biological provenance, mathematical proofs, and cohort datasets.

---

## Table of Contents

- [Performance & Scaling Scoreboard](#performance--scaling-scoreboard)
- [Installation](#installation)
- [Usage & Pipeline Recipes](#usage--pipeline-recipes)
- [Operational Modalities](#operational-modalities)
- [Empirical Parity & Validation Status](#empirical-parity--validation-status)
- [Architectural Parity & Upstream Policy](#architectural-parity--upstream-policy)
- [CLI Parameter Reference](#cli-parameter-reference)
- [AI Tooling & Verification Transparency](#ai-tooling--verification-transparency)
- [Attribution & Citation](#attribution--citation)
- [License](#license)

---

<a id="performance--scaling-scoreboard"></a>
<a id="performance-scaling-scoreboard"></a>
## Performance & Scaling Scoreboard

Evaluated on Intel i7-13700H, Linux x86_64, 8 worker threads (`-th 8`; `hg008-tier2` and `seqc2-exome-tier2` at `-th 4`²), comparing `vardictrust` against the `VarDictJava v1.8.3` pipeline (Java + R + Perl) and `vardictcpp`¹:

| Workload / Benchmark | Threads | VarDictJava Wall Time | vardictcpp Wall Time¹ | VarDictRust Wall Time | Speedup (vs Java) | VarDictJava Peak RSS | vardictcpp Peak RSS¹ | VarDictRust Peak RSS | Memory Use (vs Java) | Evidence Artifact |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :---: |
| **Targeted Amplicon Panel** | `8` (`-th 8`) | `0:02.55` | `0:03.75` ([`CPP-3`](docs/upstream-bugs.md#cpp-3)) | **`0:00.07`** | **36.4× faster** | `122.2 MiB` | `88.9 MiB` | **`58.7 MiB`** | **2.1× less** | [`amplicon-tier2.json`](benchmarks/evidence/amplicon-tier2.json) |
| **GIAB HG002 (chr20:10M–30M, 30×)** | `8` (`-th 8`) | `1:31.86` | `0:18.87` ([`CPP-9`](docs/upstream-bugs.md#cpp-9)) | **`0:14.75`** | **6.2× faster** | `4.20 GiB` | `853.1 MiB` | **`497.1 MiB`** | **8.7× less** | [`hg002-tier2.json`](benchmarks/evidence/hg002-tier2.json) |
| **SEQC2 HCC1395 Somatic Exome** | `4` (`-th 4`²) | `53:53.76` | `7:11.99` ([`CPP-9`](docs/upstream-bugs.md#cpp-9)) | **`7:46.34`** | **6.9× faster** | `3.52 GiB` | `2.15 GiB` | **`467.9 MiB`** | **7.7× less** | [`seqc2-exome-tier2.json`](benchmarks/evidence/seqc2-exome-tier2.json) |
| **GIAB HG008 Somatic PDAC** | `4` (`-th 4`²) | `1:27.68` | `0:21.96` ([`CPP-9`](docs/upstream-bugs.md#cpp-9)) | **`0:20.44`** | **4.3× faster** | `3.29 GiB` | `303.6 MiB` | **`313.1 MiB`** | **10.8× less** | [`hg008-tier2.json`](benchmarks/evidence/hg008-tier2.json) |
| **Whole Chromosome 1 Scale (249 Mbp)** | `1 → 8` (`-th 8`³) | `10:48.98` (Rust `-th 1`³; Java `OOM`) | `2:48.27` ([`CPP-5`](docs/upstream-bugs.md#cpp-5)) | **`3:06.43`** | **3.48× parallel** | `>31.0 GiB` (`OOM`³) | `778.6 MiB` ([`CPP-5`](docs/upstream-bugs.md#cpp-5)) | **`1.42 GiB`** | **1.22× flat RSS**⁴ | [`scale.json`](benchmarks/evidence/scale.json) |

¹ **`vardictcpp` comparison (`vardictcpp 3 hc52dbad_0`, bioconda).** Runs beside each benchmark under identical thread and timer constraints (`LC_ALL=C /usr/bin/time -v`). Release acceptance gates strictly on `VarDictRust == VarDictJava v1.8.3`; linked `CPP-*` codes document where `vardictcpp` diverges from `VarDictJava` (e.g., [`CPP-3`](docs/upstream-bugs.md#cpp-3) TSV-only `-a`, [`CPP-5`](docs/upstream-bugs.md#cpp-5) `22.5 GiB` unchunked `-R` RSS vs `778.6 MiB` with `--chunk 50000`, [`CPP-9`](docs/upstream-bugs.md#cpp-9) caller variant deltas).
² `hg008-tier2` and `seqc2-exome-tier2` run at `-th 4` because `VarDictJava` at `-th 8` is OOM-killed on 31 GiB hosts (bisected `1✓ / 4✓ / 8✗`); exact commands in evidence.
³ On `Whole Chromosome 1 Scale`, `VarDictJava` at `-th 8` OOMs on the 31 GiB i7 host (completed in `57:16.00` on an 8-core Ryzen host; see the [`Omnibus Empirical Validation Report`](docs/validation/omnibus-empirical-validation.md)); `10:48.98` (`648.98 s`) → **`3:06.43`** (`186.43 s`) reports the 3-run median 1-core vs 8-core `VarDictRust` parallel scaling on `chr1`.
⁴ **Flat $O(1)$ memory ratio ([Issues #197 and #384](https://github.com/AstraZeneca-NGS/VarDictJava/issues/197)):** Bounded halo-tile RAII streaming keeps peak RSS at `1.42 GiB` (`1,493,656 KB`) across all 249 Mbp of `chr1` vs `1.16 GiB` (`1,221,380 KB`) on a 10 Mbp slice (**1.22× RSS ratio** across a 25× span increase; pre-registered gate $\le 1.50\times$).

---

## Installation

VarDictRust ships across five channels: pre-compiled static/native binaries (`cargo-dist` with SLSA Level 3 provenance), Cargo (`crates.io`), multi-arch Seqera Wave containers, Bioconda (`conda`/`mamba`/`pixi`), and Nextflow (`nf-core/modules`).

### 1. Pre-Compiled Binaries (`cargo-dist` & SLSA Level 3 Provenance)

Prebuilt archives on [GitHub Releases](https://github.com/palakpsheth/VarDictRust/releases/latest) cover 5 targets: static Linux (`x86_64-unknown-linux-musl`, `aarch64-unknown-linux-musl`), native Linux (`x86_64-unknown-linux-gnu`), and native macOS (`aarch64-apple-darwin`, `x86_64-apple-darwin`):

```bash
# Option A: One-line shell installer (auto-detects Linux/macOS architecture)
curl --proto '=https' --tlsv1.2 -LsSf https://github.com/palakpsheth/VarDictRust/releases/latest/download/vardictrust-installer.sh | sh

# Option B: Manual download of standalone Linux x86_64 musl archive + SHA-256 & SLSA verification
curl -LO https://github.com/palakpsheth/VarDictRust/releases/latest/download/vardictrust-x86_64-unknown-linux-musl.tar.xz
curl -LO https://github.com/palakpsheth/VarDictRust/releases/latest/download/vardictrust-x86_64-unknown-linux-musl.tar.xz.sha256
sha256sum --check vardictrust-x86_64-unknown-linux-musl.tar.xz.sha256
gh attestation verify --owner palakpsheth vardictrust-x86_64-unknown-linux-musl.tar.xz
tar -xf vardictrust-x86_64-unknown-linux-musl.tar.xz && sudo mv vardictrust /usr/local/bin/
```

### 2. Install via Cargo

```bash
# Install from crates.io
cargo install vardictrust

# Or build the latest tagged release from source
cargo install --git https://github.com/palakpsheth/VarDictRust.git --tag v1.0.4
```

### 3. Containerized Execution (Docker, Podman, Apptainer via Seqera Wave)

Multi-arch Seqera Wave containers bundle `vardictrust` and `htslib` (`bgzip`/`tabix`) from the published Bioconda package (pinned digests and exact package version tracked in [`packaging/nf-core/modules/vardictrust/meta.yml`](packaging/nf-core/modules/vardictrust/meta.yml) and [`environment.yml`](packaging/nf-core/modules/vardictrust/environment.yml)):

```bash
# Docker / Podman (linux/amd64: 1920a853de03683c | linux/arm64: 88ae818a180b4691)
docker run --rm -v $(pwd):/data community.wave.seqera.io/library/htslib_vardictrust:1920a853de03683c \
  vardictrust -G /data/ref.fa -b /data/sample.bam -R chr1:1000-2000 > output.vcf

# Apptainer / Singularity ORAS (linux/amd64: 8601ae9fc8d7b8d5 | linux/arm64: 01fadb4d36749bdf)
apptainer exec oras://community.wave.seqera.io/library/htslib_vardictrust:8601ae9fc8d7b8d5 vardictrust --help
```

### 4. Conda / Bioconda (`conda`, `mamba`, `pixi`)

Pre-compiled packages (`linux-64`, `linux-aarch64`, `osx-arm64`, `osx-64`) on [Bioconda](https://anaconda.org/bioconda/vardictrust) ([`meta.yaml`](packaging/bioconda/recipes/vardictrust/meta.yaml)) install a standalone binary with no JVM, R, or Perl dependencies:

```bash
# Install into active environment (conda or mamba)
conda install -c bioconda -c conda-forge vardictrust
mamba install -c bioconda -c conda-forge vardictrust

# Create isolated environment
conda create -n vardictrust -c bioconda -c conda-forge vardictrust && conda activate vardictrust

# Install globally via pixi
pixi global install --channel conda-forge --channel bioconda vardictrust
```

### 5. Nextflow (`nf-core/modules`)

Install the [`vardictrust` module](packaging/nf-core/modules/vardictrust/main.nf) ([`meta.yml`](packaging/nf-core/modules/vardictrust/meta.yml), multi-arch Wave containers + `bgzip` output) into any Nextflow DSL2 pipeline:

```bash
nf-core modules install vardictrust
```

```groovy
include { VARDICTRUST } from './modules/nf-core/vardictrust/main'

workflow {
    // Germline: bams = sample.bam | Somatic: bams = [tumor.bam, normal.bam]
    // Default task.ext.args = '-c 1 -S 2 -E 3 -z'
    VARDICTRUST (
        ch_bam_bai_bed,   // tuple val(meta), path(bams), path(bais), path(bed)
        ch_fasta,         // tuple val(meta2), path(fasta)
        ch_fasta_fai      // tuple val(meta3), path(fasta_fai)
    )
}
```

### 10-Second Installation Verification

Verify BAM decoding, CIGAR realignment, and VCF emission:

```bash
# 1. Verify build and version metadata
vardictrust --version

# 2. Functional smoke test (from a repository checkout; emits single synthetic SNV in <0.1s)
vardictrust -G fixtures/micro/micro.ref.fa -b fixtures/micro/micro.bam -R chr1:40-260
```

---

<a id="usage-recipes"></a>
<a id="usage--pipeline-recipes"></a>
## Usage & Pipeline Recipes

> [!WARNING]
> **The Silent BED Coordinate Trap (`-c 1 -S 2 -E 3`):**
> - **Standard 3/4-Column BED Files (`chr start end [name]`):** VarDict defaults to `-c 3 -S 7 -E 8` (1-based) for legacy 8-column amplicon BED files. Pass **`-c 1 -S 2 -E 3`** (plus **`-z`** for 0-based half-open BED files) on 3-column BEDs to prevent parser aborts.
> - **Direct Region Specification:** Pass `-R chr:start-end[:gene]` to bypass BED parsing.

### 1. Germline Single-Sample Calling

Constitutional variant discovery from single-sample WGS, WES, or targeted panels:

```bash
vardictrust -G GRCh38.fa -f 0.01 -N SAMPLE -b sample.bam -c 1 -S 2 -E 3 -z targets.bed > germline.vcf
```

### 2. Somatic Paired Calling (Tumor vs Normal)

Pair-matched somatic calling with LOH detection and somatic `STATUS` classification (`StrongSomatic`, `LikelySomatic`, `Germline`, `SampleSpecific`):

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

Stream from indexed BAM/CRAM files into compressed, tabix-indexed VCF archives:

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
| **Pileup-Only Inspection** | `SimpleMode` / `AmpliconMode` (`-p`) | Any Indexed BAM / CRAM | `VCFv4.2` or `.var` across every reference locus | Unfiltered per-base depth, reference quality, and non-variant pileup auditing. |

> [!NOTE]
> **Assumptions & Coordinate Conventions:**
> - **Coordinate Contracts:** BAM/CRAM alignments use 0-based half-open intervals (`[start, end)`); VCF records use 1-based closed intervals (`[start, end]`). Passing `-z` shifts 0-based BED intervals by `+1`.
> - **Statistical Rounding:** Variant frequencies (`AF`), quality scores (`QUAL`), and p-values use IEEE 754 round-to-nearest-even (**Banker's rounding**), matching Java `RoundingMode.HALF_EVEN` and R `round()` ([`round_half_even`](DIVERGENCES.md#statistical-rounding--floating-point-half-even)).
> - **Multi-Allelic Collapsing:** Single-sample mode selects the top-AF candidate per locus matching [`var2vcf_valid.pl`](DIVERGENCES.md#consolidation-of-the-three-language-pipeline).

---

<a id="validation-status"></a>
<a id="empirical-parity--validation-status"></a>
<a id="empirical-parity-validation-status"></a>
## Empirical Parity & Validation Status

Every release passes dual-engine differential testing against **VarDictJava v1.8.3** (commit `0886261`, Perl submodule `009e017`) and matches biological ground truth across reference benchmarks:

### 1. Three-Way Parity Matrix (VarDictJava vs VarDictRust vs vardictcpp)

All comparisons evaluate end-to-end VCF output under float tolerance $\varepsilon = 10^{-6}$:

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

\* *The 91 PASS / 1 XFAIL metric compares against static files in `VarDictJava/testdata/integrationtestcases/`. For `hard_clip_case.fa`, the upstream fixture encodes conflicting reference chunks ([Class D / UB-4](docs/upstream-bugs.md#ub-4)); against live `VarDictJava v1.8.3`, VarDictRust matches 100% (7/7 `.var` rows, 5/5 VCF variants, 0 diffs).*
¹ **`vardictcpp` third-arm comparison (`vardictcpp 3 hc52dbad_0`, bioconda).** Executed beside each Tier-1/Tier-2 comparator arm; the release gate remains strictly `VarDictRust == VarDictJava v1.8.3`. Linked `CPP-*` identifiers reference [`docs/upstream-bugs.md`](docs/upstream-bugs.md#cpp-divergences).

- **Intermediate Table Parity ([Tier-1](docs/testing/test-matrix-and-ci-audit.md#tier-1)):** Intermediate `.var` tables match 100% of upstream columns and ordering (Simple: 36/38; Amplicon: 38/40; Somatic: 55/61 without/with `--fisher`). Standard CLI output is VCFv4.2 (10 or 11 columns); `-i` emits a 4-column intron count table.
- **Full-Corpus 12-Cell Sweep ([Omnibus Empirical Validation Report](docs/validation/omnibus-empirical-validation.md)):** Certified **89 / 89 MATCH** (`fp=0, fn=0` across `1,202,444` oracle `.var`/`.vcf`/splice records) spanning `7,361×–173,132×` ultra-deep AML MRD panels, `51.8 GiB` HG002 STAR RNA-Seq (`-i` splice & cDNA), whole-chromosome `-R chr1:1-249250621`, `--fisher` stress matrices, and CRAM 3.0 twins.

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

¹ **`vardictcpp` vcfeval normalization (`vardictcpp 3 hc52dbad_0`).** Native `vardictcpp -vcf` omits `##contig` headers ([`CPP-6`](docs/upstream-bugs.md#cpp-6)) and strips `chr` on ambiguous b37 FASTAs ([`CPP-7`](docs/upstream-bugs.md#cpp-7)); the concordance harness injects contigs before `rtg vcfeval`. Residual 1–16 FP deltas on raw HG002, SEQC2, and HG008 reflect [`CPP-9`](docs/upstream-bugs.md#cpp-9); `VarDictRust` and `VarDictJava` produce identical TP, FP, FN, Precision, Recall, and $F_1$ across every threshold.

*See the [`Omnibus Empirical Validation Report`](docs/validation/omnibus-empirical-validation.md) for reference checksums and reproduction commands across all 10 cohorts and 17 modalities, and [`docs/testing/test-matrix-and-ci-audit.md`](docs/testing/test-matrix-and-ci-audit.md) for the Seven-Level Testing Matrix.*

---

<a id="architectural-parity--upstream-policy"></a>
<a id="architectural-parity-upstream-policy"></a>
## Architectural Parity & Upstream Policy

VarDictRust enforces the **[rewrites.bio](https://rewrites.bio) Prime Directive**:

1. **Fidelity over Innovation:** Matches every VarDictJava v1.8.3 output decision on valid inputs, including IEEE 754 banker's rounding ([`round_half_even`](DIVERGENCES.md#statistical-rounding--floating-point-half-even)), exact odds ratios, OpenJDK [`JavaHashMap`](DIVERGENCES.md#ub-3-var-row-order-across-positions-follows-hashmap-iteration-order-282-family) iteration order, and in-engine somatic classification ([`SomaticPostProcessModule`](DIVERGENCES.md#somatic-status-determination)).
2. **Robustness by Construction:** Enforces `#![forbid(unsafe_code)]` and returns `Result<_, VardictError>` without panicking where VarDictJava crashes on malformed CIGARs ([`modifyCigar` #341 family](docs/upstream-bugs.md), [PR #402](DIVERGENCES.md#category-1-cigar-and-realignment-resilience)), null pointers, or missing indexes ([Issue #187](DIVERGENCES.md#category-7-cli-diagnostics-and-error-messages)).
3. **Bounded Memory at Any Depth:** Processes reads in bounded halo tiles, eliminating whole-chromosome JVM heap exhaustion ([Issues #197 and #384](https://github.com/AstraZeneca-NGS/VarDictJava/issues/197)).
4. **Bug Release Policy:** Replicates upstream output quirks in v1.0.x ([UB-1 through UB-7](docs/upstream-bugs.md)); output-changing fixes ship in v1.1.0+ behind opt-in `--fix <id>` flags per [`ROADMAP.md`](ROADMAP.md).
5. **Documented Quirks & Testing Matrix:** Cataloged in [`DIVERGENCES.md`](DIVERGENCES.md) and [`docs/testing/test-matrix-and-ci-audit.md`](docs/testing/test-matrix-and-ci-audit.md).

---

## CLI Parameter Reference

Flag parity with VarDictJava v1.8.3 is defined against [`CmdParser.java`](https://github.com/AstraZeneca-NGS/VarDictJava/blob/0886261/src/main/java/com/astrazeneca/vardict/CmdParser.java) and verified by differential testing:

*Evidence Status Vocabulary:*
- `corpus` — Verified against integration corpus test cases or [Tier-2](docs/testing/test-matrix-and-ci-audit.md#tier-2) oracle sweeps.
- `unexercised` — Implemented to Java source specification without dedicated corpus fixtures (`0` flags).
- `no-output` — Performance or concurrency flag; does not alter output bytes.
- `parity-exempt` — Accepted for CLI compatibility matching upstream semantics.

### 1. Core Inputs & Outputs

| Short | Long | Description | Default | Status |
| :--- | :--- | :--- | :--- | :--- |
| `-G` | `--reference` | Reference FASTA file path (`.fai` index required) | required | `corpus` |
| `-b` | `--bam` | Indexed BAM/CRAM; somatic: `"tumor.bam\|normal.bam"` | required | `corpus` |
| `-N` | `--name` | Sample name (somatic: `"TUMOR\|NORMAL"`) | derived from BAM | `corpus` |
| `-R` | `--region` | Region `chr:start-end[:gene]` instead of BED | — | `corpus` |
| positional | — | BED file (used when `-R` absent) | — | `corpus` |
| `-h` | `--header` | Print `.var` column header before results (not help; use `-H`/`--help`) | off | `corpus` |
| `-H` | `--help` | Print CLI help message and exit | — | `parity-exempt` |
| — | `--version` | Print version, upstream provenance, and AI disclosure metadata | — | `parity-exempt` |
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
| `-3` | — | Move indels to 3′ ([`moveIndelsTo3`](ROADMAP.md)) | off | `corpus` |
| `-K` | — | Include Ns in depth calculation | off | `corpus` |
| `-p` | — | Pileup-only mode | off | `corpus` |
| `-t` | `--dedup` | Remove duplicated reads (one pair per start position) | off | `corpus` |
| `-u` | — | Unique mode: count overlapping mates once (forward read) | off | `corpus` |
| `-UN` | `--UN` | Unique mode: count overlapping mates once (first read in pair) | off | `corpus` |
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
| — | `--tile-size` | Region tile size in bp (bounded-memory WGS) | auto (`<= 1 Mbp` untiled; `> 1 Mbp` `100 kbp` tiles; max `500 Mbp`) | `no-output` |

*Loud-Error Legacy Flags (`AGENTS.md §12`):* Passing any of the 5 unimplemented VarDictJava options (7 flag spellings: `-n`, `-C`, `-T` / `--trim`, `--VS`, `-Z` / `--downsample`) fails fast with `VardictError::Unsupported("<flag> not implemented in v1.0; use VarDictJava")` rather than ignoring the flag.

*Reproducible Profiling Hooks:* When compiled with opt-in Cargo features (`profile-cpu`, `profile-heap`, or `profile-dhat`), `vardictrust` honors `VARDICT_PROFILE_CPU`, `VARDICT_PROFILE_HEAP`, and `VARDICT_PROFILE_DHAT` output directory environment variables (driven via [`harness/profile.sh`](harness/README.md)).

---

<a id="ai-disclosure"></a>
<a id="ai-tooling--verification-transparency"></a>
<a id="ai-tooling-verification-transparency"></a>
## AI Tooling & Verification Transparency

Per the [rewrites.bio](https://rewrites.bio) manifesto, we disclose all AI tools, their operational roles, and human verification gates:

- **AI Systems Used:** Google DeepMind Gemini / Antigravity CLI, Anthropic Claude Code, and Qoder.
- **Scope of AI Assistance:** Implementation of core Rust modules under human architectural direction, transliteration of R/Perl statistical routines into pure Rust, test harness engineering, and adversarial differential verification.
- **Verification Guarantee (Evidence Before Completion):** **Zero** AI-generated code was accepted on model assertion. Every algorithmic component was validated against VarDictJava v1.8.3 across 416 automated Cargo workspace tests (252 core library tests, 110 corpus and modality parity tests, 25 engine behavior tests, 2 amplicon streaming memory tests, 2 CRAM stack tests, 17 CLI tests, and 8 CLI integration binary tests) enforcing a `>= 90%` line-coverage gate (`cargo-llvm-cov`), 25 harness diagnostic self-tests, the 92-case integration corpus, 10 slice fixtures, 18 multi-modal slices across 6 cohorts, the 89-job 12-cell full-corpus sweep, and multi-gigabyte public benchmarks (SEQC2, GIAB HG008, GIAB HG002, GIAB HG005, GIAB NA24385 b37) under float tolerance $\varepsilon = 10^{-6}$, with reproducibility documented in [`docs/testing/test-matrix-and-ci-audit.md`](docs/testing/test-matrix-and-ci-audit.md) and workflow skills in [`.agents/skills/`](AGENTS.md#19-repository-specific-agent-skills).

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

If you use VarDictRust in published research, cite both Lai et al. (2016) and this software:

```bibtex
@software{sheth2026vardictrust,
  author    = {Sheth, Palak},
  title     = {{VarDictRust: Memory-Safe, High-Performance Variant Calling for NGS}},
  year      = {2026},
  version   = {1.0.4},
  url       = {https://github.com/palakpsheth/VarDictRust}
}
```
Full citation metadata is available in [`CITATION.cff`](CITATION.cff).

---

## License

Original VarDictRust source code is licensed under the [MIT License](LICENSE). Combined source and binary distributions carry the SPDX license expression `MIT AND GPL-2.0-or-later`, covering translated R Project `nmath` statistical routines (`crates/vardictrust-core/src/stats/rmath/`) and third-party dependency notices cataloged in [`LICENSE-3RD-PARTY`](LICENSE-3RD-PARTY) and [`deny.toml`](deny.toml).
