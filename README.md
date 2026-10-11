# VarDictRust (`vardictrust`)

[![CI](https://github.com/palakpsheth/VarDictRust/actions/workflows/ci.yml/badge.svg)](https://github.com/palakpsheth/VarDictRust/actions)
[![Security Policy](https://img.shields.io/badge/security-policy-blue.svg)](.github/SECURITY.md)
[![SLSA 3](https://img.shields.io/badge/SLSA-Level%203-green.svg)](https://slsa.dev)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![rewrites.bio](https://img.shields.io/badge/rewrites.bio-manifesto-blueviolet.svg)](https://rewrites.bio)

An idiomatic, memory-safe **Rust** rewrite of [VarDictJava](https://github.com/AstraZeneca-NGS/VarDictJava) v1.8.3 — a sensitive somatic, germline, and amplicon variant caller for next-generation sequencing (NGS) data.

VarDictRust consolidates the original 3-language pipeline (`Java` → `R` → `Perl`) into a **single, statically-compiled binary** that produces output **identical to VarDictJava v1.8.3 on valid inputs ($\varepsilon = 10^{-6}$, 0 variant diffs)**, eliminates the JVM/Perl/R dependency stack, and **scales to deep whole-genome sequencing (WGS) with bounded $O(1)$ memory**.

Built under the [rewrites.bio](https://rewrites.bio) manifesto: *a faster, more robust tool that produces the same results*. See the comprehensive [Omnibus Master Empirical Validation Report](docs/validation/omnibus-empirical-validation.md) for full biological provenance, algorithmic proofs, and benchmark datasets.

```bash
# Quick Example: Somatic paired calling (tumor vs normal)
vardictrust -G GRCh38.fa -f 0.01 -N "TUMOR|NORMAL" -b "tumor.bam|normal.bam" target.bed > somatic.vcf
```

---

## Performance & Scaling Scoreboard

- **Exact Parity:** 0 variant discrepancies and identical `FILTER`/`STATUS` decisions against VarDictJava v1.8.3 across all 15 dual-engine differential benchmark runs and 92 integration fixtures ($\varepsilon = 10^{-6}$).
- **High-Throughput Speedup:** Up to **30.0× faster** on targeted amplicon assays and **4.85× parallel speedup** on 8 cores across whole human chromosomes.
- **Bounded $O(1)$ Memory:** Flat resident memory scaling ($1.39\times$ peak RSS ratio across a 25-fold chromosome span increase), eliminating JVM heap exhaustion on deep WGS.
- **Self-Contained Deployment:** Single statically-linked binary (`x86_64-unknown-linux-musl`) with zero runtime dependencies (no JVM, R, Perl, or conda environment needed).

Recorded on identical hardware (AMD Ryzen / Ubuntu 24.04 LTS, NVMe SSD) pinned to 8 CPU threads:

| Workload / Benchmark | Threads | VarDictJava Wall Time | VarDictRust Wall Time | Speedup Ratio | VarDictJava Peak RSS | VarDictRust Peak RSS | Memory Scaling Ratio |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Targeted Amplicon Panel** | 8 cores | 0:02.10 | 0:00.07 | **30.0× faster** | 120.8 MB | 50.8 MB | 0.42× (2.4× lower) |
| **GIAB HG002 (chr20:10M–30M, 30×)** | 8 cores | 4:16.98 | 1:22.47 | **3.1× faster** | 3.6 GB | 1.1 GB | 0.30× (3.3× lower) |
| **SEQC2 HCC1395 Somatic Exome** | 8 cores | 54:23.00 | 19:45.73 | **2.8× faster** | 4.1 GB | 1.4 GB | 0.36× (2.8× lower) |
| **GIAB HG008 Somatic PDAC** | 8 cores | 1:10.99 | 0:49.48 | **1.4× faster** | 5.7 GB | 1.2 GB | 0.21× (4.8× lower) |
| **Whole Chromosome 1 Scale (249 Mbp)** | 1 → 8 cores | 57m 16s (1 core) | 11m 49s (8 cores) | **4.85× parallel** | OOM / Unbounded | 2.9 GB (249 Mb) vs 2.1 GB (10 Mb) | **1.39× flat RSS scaling** |

> [!NOTE]
> **Bounded $O(1)$ Memory at Any Depth:** VarDictJava holds entire regional chromosome structures in JVM heap collections, leading to JVM heap exhaustion and Out-Of-Memory errors on deep WGS (issues #197, #384). VarDictRust processes genomic regions via **bounded tiles with halo intervals** using RAII streaming. Across a 25-fold increase in genomic span (10 Mbp to all 249 Mbp of chromosome 1), VarDictRust peak RSS scales by only **$1.39\times$** (2.11 GB → 2.94 GB).

---

## Three Modalities to Run It

| Modality | Flag Signature | Input Assays | Output Columns | Typical Use Case |
| :--- | :--- | :--- | :--- | :--- |
| **Germline Calling** | `SimpleMode` (`-b sample.bam`) | WGS, WES, Panels | 36 (default) / 38 (`--fisher`) | Constitutional variant discovery, family trios, population genomics. |
| **Somatic Paired Calling** | `SomaticMode` (`-b "tumor\|normal"`) | Tumor / Normal WGS / WES | 55 (default) / 61 (`--fisher`) | Cancer somatic mutations, subclonal alleles down to 1% AF, LOH detection. |
| **Targeted Amplicon Panel** | `AmpliconMode` (`-a 10:0.95`) | AmpliSeq, HaloPlex, PCR | 38 (default) / 40 (`--fisher`) | Overlapping amplicon panels; trims primer bias (`AMPBIAS`) using 8-column BED. |

> [!NOTE]
> **Assumptions & Coordinate Conventions:**
> - **Coordinate Contracts:** BAM/CRAM alignment coordinates are parsed as 0-based half-open intervals (`[start, end)`), while VCF output records are strictly 1-based closed intervals (`[start, end]`). Passing `-z 1` (or `--zero-based`) shifts 0-based BED intervals by `+1`.
> - **Statistical Rounding:** All variant frequencies (`AF`), quality scores (`QUAL`), and p-values adhere to IEEE 754 round-to-nearest, ties-to-even (**Banker's rounding**), replicating Java's `RoundingMode.HALF_EVEN` and R's `round()`.
> - **Multi-Allelic Collapsing:** Standard single-sample mode defaults to single-winner top-AF candidate selection per locus matching `var2vcf_valid.pl`.

---

## Quick Start & Pipeline Recipes

### Installation
VarDictRust is distributed as a single standalone, statically-linked executable (`x86_64-unknown-linux-musl`) with zero external dependencies:

```bash
# Option 1: Cargo install
cargo install vardictrust

# Option 2: Pre-built static musl binary (no glibc dependency)
curl -LsSf https://github.com/palakpsheth/VarDictRust/releases/latest/download/vardictrust-x86_64-unknown-linux-musl.tar.gz | tar xz
sudo mv vardictrust /usr/local/bin/

# 10-second verification (confirms installation & prints build metadata)
vardictrust --version
```

### 1. Germline Single-Sample Calling
```bash
vardictrust -G ref.fa -f 0.01 -N SAMPLE -b sample.bam regions.bed > germline.vcf
```

### 2. Somatic Paired Tumor / Normal Calling
```bash
vardictrust -G ref.fa -f 0.01 -N "TUMOR|NORMAL" -b "tumor.bam|normal.bam" target.bed > somatic.vcf
```

### 3. Targeted Amplicon Panel Calling (8-Column BED)
```bash
vardictrust -G ref.fa -f 0.01 -N SAMPLE -a 10:0.95 -b amplicon.bam amplicons.bed > amplicon.vcf
```

### 4. High-Throughput UNIX Streaming Pipeline
VarDictRust streams standards-compliant VCF directly to stdout, enabling pipe integration with `bcftools` and `bgzip`:
```bash
vardictrust -G ref.fa -f 0.01 -N "TUMOR|NORMAL" -b "tumor.bam|normal.bam" target.bed -th 8 \
  | bcftools filter -i 'FILTER="PASS" && STATUS="StrongSomatic"' -Oz -o somatic.filtered.vcf.gz \
  && tabix -p vcf somatic.filtered.vcf.gz
```

---

## Empirical Parity & Validation Status

Every release candidate must pass dual-engine differential testing against **VarDictJava v1.8.3** (commit `0886261`, Perl submodule `009e017`) and match physical biological truth across gold-standard reference benchmarks:

### 1. Dual-Engine Parity Matrix (VarDictJava vs VarDictRust)
All comparisons evaluate the full end-to-end VCF output under float tolerance $\varepsilon = 10^{-6}$:

| Dataset / Benchmark | Modality | Target Region / Assay | VarDictJava Variants | VarDictRust Variants | Tier-2 Parity | Parity Conformance Details |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Upstream Integration Corpus** | Simple, Somatic, Amplicon | 92 micro-fixtures | 91 PASS / 1 FAIL* | 91 PASS / 1 FAIL* | **100% MATCH*** | Identical intermediate `.var` tables |
| **GIAB HG002 Germline WGS** | SimpleMode (Germline) | chr20:10M–30M (30×) | 227,562 | 227,562 | **MATCH** | 0 diffs across 20 Mbp ($\varepsilon = 10^{-6}$) |
| **SEQC2 HCC1395/BL** | SomaticMode (Tumor/Normal) | Whole Exome (WES) | 623,444 | 623,444 | **MATCH** | 0 diffs across all exons ($\varepsilon = 10^{-6}$) |
| **SEQC2 Somatic Slices** | SomaticMode (Tumor/Normal) | chr1, 2, 3, 17, 19, 21 slices | 4,151 | 4,151 | **MATCH** | 0 diffs in SNV, InDel, and SV |
| **GIAB HG008 Somatic** | SomaticMode (Tumor/Normal) | PDAC high-depth slices | 18,426 | 18,426 | **MATCH** | 0 diffs in somatic status |
| **Targeted Amplicon Panel** | AmpliconMode (`-a 10:0.95`) | Overlapping primer assay | 2 | 2 | **MATCH** | 0 diffs with primer trimming & AMPBIAS |
| **Chr1 Scale Benchmark** | SimpleMode (Germline) | chr1: 1 Mbp to 249 Mbp | 2,840,119 | 2,840,119 | **MATCH** | 0 diffs across whole chromosome |
| **Extended Modalities Suite** | RNA Splicing, `<INV>`, `<DUP>`, CRISPR, Flags | 9 real biological assays | 9 PASS | 9 PASS | **100% MATCH** | 0 diffs in table & VCF strict-bytes |
| **Multi-Modal Multi-Slice Suite** | Germline, Somatic, RNA | 18 slices across 6 datasets & 9 chromosomes | 1,129 | 1,129 | **100% MATCH** | 0 diffs across GRCh38 & GRCh37 |

\* *The 91 PASS / 1 FAIL metric reflects comparison against stale pre-computed static files in the upstream repository (`VarDictJava/testdata/integrationtestcases/`). For the single static mismatch (`hard_clip_case.fa`), the upstream fixture encodes conflicting reference chunks (Class D / UB-4 in [`docs/upstream-bugs.md`](docs/upstream-bugs.md)). When evaluated against the live `VarDictJava v1.8.3` oracle, VarDictRust produces 100% byte-identical output (7/7 table rows, 5/5 VCF variants, 0 diffs). Thus, live oracle parity is 100% across all 92 fixtures.*

- **Intermediate Table Parity (100% Column Parity):** Emitted intermediate `.var` tables match 100% of upstream columns and ordering. Slashes denote standard mode vs. statistical `--fisher` mode (Simple: 36/38 columns; Amplicon: 38/40 columns; Somatic: 55/61 columns). Standard CLI output is standard VCFv4.2 (10 or 11 columns); Splicing Mode (`-i`) outputs a 4-column intron count table.

### 2. Ground Truth Concordance (RTG vcfeval)
Evaluated against physical ground truth using Real Time Genomics `rtg vcfeval`:

| Benchmark Dataset | Truth Standard | Variant Assay / Threshold | Evaluated Sites (TP / FP / FN) | VarDictJava vcfeval | VarDictRust vcfeval | Truth Concordance |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **GIAB HG002** (Operating Point) | NIST v4.2.1 High-Confidence | Germline Small (`QUAL >= 128`) | 22,748 TP / **553 FP** / 871 FN | **F1: 0.9699 (P: 0.9763, R: 0.9635)** | **F1: 0.9699 (P: 0.9763, R: 0.9635)** | **EQUAL** |
| **GIAB HG002** (Raw Baseline) | NIST v4.2.1 High-Confidence | Germline Small (`None` / raw) | 22,803 TP / 22,458 FP / 813 FN | F1: 0.6622 (P: 0.5038, R: 0.9659) | F1: 0.6622 (P: 0.5038, R: 0.9659) | **EQUAL** |
| **SEQC2 HCC1395** | SEQC2 v1.2.1 Somatic Truth | High-Confidence Somatic SNV/InDel | 1,056 TP / 11,334 FP / 140 FN | F1: 0.1555 (P: 0.0852, R: 0.8842) | F1: 0.1555 (P: 0.0852, R: 0.8842) | **EQUAL** |
| **GIAB HG008** | NIST Draft Benchmark v0.3 | Somatic Small Variants | 575 TP / 1,565 FP / 288 FN | F1: 0.3856 (P: 0.2687, R: 0.6825) | F1: 0.3856 (P: 0.2687, R: 0.6825) | **EQUAL** |

> [!NOTE]
> **Interpreting Raw Discovery Sensitivity vs. Consensus Truth Set Precision:**
> VarDict was designed as an ultra-sensitive discovery engine capable of calling subclonal cancer mutations down to 1% allele frequency (`-f 0.01`). Ground truth consensus sets (NIST v4.2.1, SEQC2 v1.2.1) are intentionally restricted to constitutional diploid or high-confidence somatic consensus variants, so raw discovery calls without post-filtering score thousands of subclonal candidates as technical false positives. At the standard operating point (`QUAL >= 128`), false positives collapse to **553** while preserving 96.35% recall (**96.99% $F_1$ score**). Crucially, **VarDictJava and VarDictRust produce bitwise identical site lists** across all benchmark sets. For full biological details, see [Section 5.1 of the Validation Report](docs/validation/omnibus-empirical-validation.md#51-biological--methodological-interpretation-of-raw-caller-precision) and [Pillar 1: Full NIST v4.2.1 HG002 Benchmark Evaluation](docs/validation/omnibus-empirical-validation.md#2-pillar-1-full-nist-v421-hg002-benchmark-evaluation).

*For full scientific details, reference checksums, and exact reproduction commands across all 10 public cohorts and 9 operational modalities, see the [`Omnibus Empirical Validation Report`](docs/validation/omnibus-empirical-validation.md). For the complete Seven-Level Testing Matrix and CI/CD audit, see [`docs/testing/test-matrix-and-ci-audit.md`](docs/testing/test-matrix-and-ci-audit.md).*


---

## Architectural Parity & Upstream Policy

VarDictRust is built under strict compliance with the **rewrites.bio Prime Directive**:
1. **Fidelity over Innovation:** VarDictRust v1.0.0 matches **every** output decision of VarDictJava v1.8.3 on valid inputs, including floating-point half-even banker's rounding (`round_half_even`), exact odds ratios, and in-engine somatic status determination (`SomaticPostProcessModule`).
2. **Robustness by Construction:** While VarDictJava crashes on malformed CIGAR strings (issue #341 family), NPEs, or missing indexes, VarDictRust enforces `#![forbid(unsafe_code)]` and returns actionable `Result<_, VardictError>` errors without panicking.
3. **Bug Release Policy:** Known upstream VarDictJava **output bugs** are faithfully replicated in v1.0.0 and tracked in [`docs/upstream-bugs.md`](docs/upstream-bugs.md). Corrective bug fixes will ship in v1.1.0+ behind explicit opt-in flags (`--fix <id>`) per [`ROADMAP.md`](ROADMAP.md).
4. **Documented Quirks:** All known design divergences are tracked in [`DIVERGENCES.md`](DIVERGENCES.md).

---

## Supported Flags

Flag parity with VarDictJava v1.8.3 is defined against `CmdParser.java` and verified by differential testing:

*Evidence Status Vocabulary:*
- `corpus` — Verified against integration corpus test cases or Tier-2 oracle sweeps.
- `unexercised` — Implemented to Java source specification, but lacks explicit differential corpus fixtures.
- `no-output` — Performance or concurrency flag; does not alter output bytes.
- `parity-exempt` — Accepted for CLI compatibility matching upstream semantics.

| Short | Long | Description | Default | Status |
| :--- | :--- | :--- | :--- | :--- |
| `-G` | `--reference` | Reference FASTA (`.fai` required; hg19/hg38/mm10 shortcuts accepted) | required | `corpus` |
| `-b` | `--bam` | Indexed BAM/CRAM; somatic: `"tumor.bam\|normal.bam"` | required | `corpus` |
| `-f` | `--af-threshold` | Allele-frequency threshold | `0.01` | `corpus` |
| `-N` | `--name` | Sample name (somatic: `"TUMOR\|NORMAL"`) | derived from BAM | `corpus` |
| `-R` | `--region` | Region `chr:start-end[:gene]` instead of BED | — | `corpus` |
| `-c` | — | BED chromosome column (1-based) | `3` | `corpus` |
| `-S` | — | BED start column | `7` | `corpus` |
| `-E` | — | BED end column | `8` | `corpus` |
| `-s` | — | BED segment-start column | `10` (= `-S` if set) | `unexercised` |
| `-e` | — | BED segment-end column | `11` (= `-E` if set) | `unexercised` |
| `-g` | — | BED gene column | `13` | `corpus` |
| — | `--threads` | Thread count (`-th`/`--th` aliases; bare `-th` = all cores) | all cores | `no-output` |
| `-k` | — | Local realignment 0/1 | `1` | `corpus` |
| `-a` | `--amplicon` | Amplicon mode `"int:float"` (distance:fraction) | auto `10:0.95` on 8-col BED | `corpus` |
| `-r` | — | Minimum variant reads | `2` | `corpus` |
| `-q` | — | Minimum base quality | `22.5` | `corpus` |
| `-Q` | — | Minimum mapping quality | — | `corpus` |
| `-P` | — | Mean read position filter | `5` | `corpus` |
| `-m` | — | Maximum mismatches | `8` | `corpus` |
| `-M` | — | Min match+insertion read length (shorter reads ignored) | `0` | `unexercised` |
| `-V` | — | Min normal AF (somatic) | `0.05` | `corpus` |
| `-U` | `--nosv` | Disable SV calling | off | `corpus` |
| `-L` | — | SV minimum length | `1000` | `corpus` |
| — | `--emit-intermediate` | Emit legacy `.var` table instead of VCF | off | `corpus` |
| — | `--fisher` | In-engine Fisher exact test (p-value / odds ratio) | off | `corpus` |
| — | `--tile-size` | Region tile size in bp (bounded-memory WGS) | auto | `no-output` |
| `-x` | — | BED region extension (bp) | `0` | `unexercised` |
| `-X` | — | Nucleotide extension | `2` | `corpus` |
| `-3` | — | Move indels to 3′ (Java `moveIndelsTo3`) | off | `unexercised` |
| `-B` | — | Minimum reads for bias | `2` | `corpus` |
| `-o` | — | qratio filter | `1.5` | `corpus` |
| `-O` | — | Output MAPQ threshold | `0` | `unexercised` |
| `-I` | — | Indel size | `50` | `corpus` |
| `-z` | — | BED coordinates are 0-based | off (1-based) | `corpus` |
| `-F` | — | SAM flag filter | `0x504` | `corpus` |
| `-w` | — | SV read length | `300` | `corpus` |
| `-W` | — | SV flanking | `100` | `corpus` |
| `-A` | — | SV minimum AF | `4` | `corpus` |
| `-Y` | — | Max region length | `1200` | `corpus` |
| `-h` | `--header` | Print `.var` column header before results | off | `corpus` |
| `-v` | — | Accepted and ignored (never read by VarDictJava v1.8.3) | — | `parity-exempt` |
| — | `--chimeric` | Chimeric read handling | off | `corpus` |
| — | `--deldupvar` | Delete/dup variant handling | off | `unexercised` |
| — | `--mfreq` | Mono frequency threshold | `0.25` | `corpus` |
| — | `--nmfreq` | Normal mono frequency threshold | `0.1` | `corpus` |
| `-J` | — | CRISPR cutting site | — | `corpus` |
| `-j` | — | CRISPR flanking size | — | `corpus` |
| `-i` | `--splice` | Splice mode | off | `corpus` |
| `-p` | — | Pileup-only mode | off | `corpus` |
| `-D` | `--debug` | Debug mode (append per-read debug columns to `.var`) | off | `corpus` |
| — | `--DP` | Debug print destination (`OUT` or `ERR`) | — | `unexercised` |
| `-y` | `--verbose` | Verbose mode | off | `no-output` |
| `-t` | `--dedup` | Remove duplicated reads (one pair per start position) | off | `unexercised` |
| — | `--adaptor` | Adaptor sequence(s), comma-separated, excluded from realignment | — | `corpus` |
| `-d` | — | Delimiter for splitting region info | `\t` | `corpus` |
| `-u` | — | Unique mode: count overlapping mates once (forward read) | off | `corpus` |
| — | `--UN` | Unique mode: count once using first read only | off | `corpus` |
| `-K` | — | Include Ns in depth calculation | off | `corpus` |
| positional | — | BED file (used when `-R` absent) | — | `corpus` |

---

## AI Tooling & Verification Transparency Disclosure

Per the [rewrites.bio](https://rewrites.bio) manifesto, we disclose all AI tools, their operational roles, and human verification gates:

- **AI Coding Systems Used:** Google DeepMind Gemini / Antigravity CLI, Anthropic Claude Code, Qoder, and spec-kit.
- **Scope of AI Assistance:** Implementation of core Rust modules under direct human architectural instruction, transliteration of legacy R/Perl statistical functions into pure Rust (`statrs`), test harness development, and adversarial failure analysis.
- **Verification Guarantee (Evidence Before Completion):** **Zero** AI-generated code was accepted based on model assertion. Every single algorithmic component was validated against VarDictJava v1.8.3 across 208 automated workspace tests, the 92-case integration corpus, 10 slice fixtures, 18 multi-modal slices across 6 cohorts, and multi-gigabyte public benchmark datasets (SEQC2, GIAB HG008, GIAB HG002, GIAB HG005, GIAB NA24385 b37) under strict float tolerance $\varepsilon = 10^{-6}$, with complete push-button reproducibility documented in [`docs/testing/test-matrix-and-ci-audit.md`](docs/testing/test-matrix-and-ci-audit.md).

---

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
  version   = {1.0.0},
  url       = {https://github.com/palakpsheth/VarDictRust}
}
```
Full citation metadata is available in [`CITATION.cff`](CITATION.cff).

---

## License

VarDictRust is open-source software licensed under the [MIT License](LICENSE). Third-party attribution, component licenses (including R Project `nmath` distribution notices), and copyright disclosures are cataloged in [`LICENSE-3RD-PARTY`](LICENSE-3RD-PARTY).
