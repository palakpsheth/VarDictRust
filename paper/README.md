# VarDictRust Paper Reproducibility Package (`v1.1.0` Clinical Study & `v1.0.6` Cloud Study)

This package provides a deterministic, single-command reproducibility pipeline for all figures, tables, text statistics, and LinkedIn dissemination assets accompanying the `VarDictRust` evaluation manuscript:

> **"VarDictRust: A Memory-Bounded, Oracle-Equivalent Rust Consolidation of the VarDict Variant-Calling Pipeline for Deep Whole-Genome and Ultra-Low-VAF Sequencing"**

---

## 1. Quickstart: Offline Reproduction (After One-Time Setup)

Two complementary, cryptographically checksummed empirical study bundles power the manuscript without mixing cross-host timings:
1. **`benchmarks/gcp-study-v106/` (`Figures 1–8`, `Tables 1–8`, `Supplementary Tables S1–S6`):** All aggregated metrics, per-cell runtime/RSS measurements, truth evaluations, and SHA-256 digests from the 60 Google Cloud Platform (GCP `n2-standard-32`) study cells (`v1.0.6` cloud benchmark baseline, engine fingerprint `5865a2a0b24a8c17`, performance fingerprint `c490c36317b8f25b`).
2. **`benchmarks/v110-compat-study/` (`Figure 9`, `Table 9`, `Supplementary Table S7`):** Host-independent deterministic callset, strict VCFv4.2/v4.3 schema-validator, `bcftools norm --check-ref e` FASTA concordance, `QIAseq V3` `-sp` single-primer extension, and diploid `rtg vcfeval` subclonal somatic metrics across 6 evaluation domains plus the 35-entry `CompatProfile` / CLI capability catalog for **VarDictRust `v1.1.0`** (`b18477ae0a7b07c1`, performance fingerprint `5c41be8047133842`; verified or re-collected by [`paper/collect_v110_compat_metrics.py`](collect_v110_compat_metrics.py)).

Alongside these two committed datasets, [`paper/build_manuscript.py`](build_manuscript.py) and [`paper/render_pdf.py`](render_pdf.py) dynamically read `Cargo.toml`, `.agents/skills/vardictrust-release/scripts/check_release_consistency.py`, `crates/vardictrust-core/src/`, and `benchmarks/evidence/*.json` to synchronize the current release architecture (`92` CLI flags across `7` functional groups including 18/18 `var2vcf_*.pl` post-filters via `--vcf-all-alleles` / `--all-alleles`, `--amplicon-single-primer` / `-sp`, `--preset` across `wgs`, `wes`, `amplicon`, and `ctdna`, dual-mode `--compat=vardictjava|off` `CompatProfile` governance with `32` opt-in switches, multithreaded BGZF `--output` + `--tbi`/`--csi` indexing, OpenJDK 8 `TreeNode` `JavaHashMap`, `8 MiB` `SpillableVcfAccumulator` VCF spill cap, singleflight `16 MiB` Illumina CRAM container cache + 2-slot `CRAM_META_TLS` + `SeekWindowReader`, `.bai`/`.csi` structural header pre-validation, automated Cargo workspace tests + `37` harness self-tests, and `cargo-fuzz` harnesses).

> **Prerequisites (one-time, needs network):** building the Python virtual environment (`.venv`) via `uv` downloads dependencies from PyPI, and PDF rendering needs a headless Chrome/Chromium binary. Do this once with network access; everything after it runs fully offline from committed data. To run the pipeline on a machine without Chrome, export `PAPER_ALLOW_STALE_PDF=1` to reuse the checked-in PDFs — without it, `reproduce_all.sh` fails fast instead of reporting stale PDFs as reproduced (see §5).

To reproduce every figure, table, compiled manuscript, LinkedIn carousel slide, and verification gate from scratch:

```bash
# From the repository root:
bash paper/reproduce_all.sh
```

This single command automatically:
1. Bootstraps a local Python virtual environment (`.venv`) using `uv` with dependencies from [`paper/requirements.txt`](requirements.txt).
2. Verifies the SHA-256 hashes and invariants of both committed study bundles against [`benchmarks/gcp-study-v106/SHA256SUMS.tsv`](../benchmarks/gcp-study-v106/SHA256SUMS.tsv) and [`benchmarks/v110-compat-study/SHA256SUMS.tsv`](../benchmarks/v110-compat-study/SHA256SUMS.tsv) (pass `--run-v110-eval` to re-execute all 6 deterministic `v1.1.0` domains locally).
3. Regenerates all **9 Main Figures** ([`paper/figures/fig1..fig9.{svg,pdf}`](figures/)) and all **16 Tables** ([`paper/tables/table1..table9.{md,tex}`](tables/) and [`paper/tables/supp_table1..supp_table7.{md,tex}`](tables/)).
4. Compiles and byte-verifies the manuscript and supplement ([`paper/manuscript.md`](manuscript.md), [`paper/supplement.md`](supplement.md)) and renders typeset PDFs ([`paper/manuscript.pdf`](manuscript.pdf), [`paper/supplement.pdf`](supplement.pdf)) with embedded vector figures.
5. Generates the 10-slide high-resolution LinkedIn carousel in `scratch/linkedin/`.
6. Executes all unit test suites, release consistency audits, GCP runner plan checks, and documentation evidence synchronizers.

---

## 2. Online Reproduction from Raw GCP Cell Artifacts

If you have cloud storage access or wish to re-aggregate directly from raw cell outputs (per-cell stdout logs, `.time` files, `rtg vcfeval` summaries, and `.var` diff receipts):

```bash
# Synchronize raw cell logs from Google Cloud Storage (~250 MB compressed metadata):
gcloud storage rsync -r \
  -x '.*\.(vcf|vcf\.gz|tbi|csi|var|bam|bai|cram|crai|bed|bed\.gz|stdout|err|fa|fai|sdf)$|(^|/).*\.(java|rust|cpp|splice|sorted)\.(tsv|var)$|(^|/)[^/]*\.var\.diff\.txt$|(^|/)(strelka_.*|fb_chunks|.*_roc\.tsv\.gz)(/|$)' \
  gs://vardict-paper-trial-results/cells-v106/ \
  benchmarks/data/gcp-study-v106/

# Or run reproduce_all.sh with the --sync-gcs flag:
bash paper/reproduce_all.sh --sync-gcs
```

When `benchmarks/data/gcp-study-v106/` contains the 60 raw cell directories, [`paper/aggregate_gcp_study.py`](aggregate_gcp_study.py) re-extracts all metrics from raw telemetry, audits SHA-256 checksums, and overwrites [`benchmarks/gcp-study-v106/master_gcp_metrics.json`](../benchmarks/gcp-study-v106/master_gcp_metrics.json).

---

## 3. Directory Layout

```
paper/
├── README.md                              # This reproduction & submission guide
├── requirements.txt                       # Pinned Python dependencies (matplotlib, pytest, etc.)
├── reproduce_all.sh                       # Master one-command offline & online runner
├── aggregate_gcp_study.py                 # 60-cell GCP JSON/TSV metric aggregator & validator
├── collect_v110_compat_metrics.py         # Deterministic v1.1.0 --compat=off & --preset collector/verifier
├── build_figures_and_tables.py            # Figure (1..9 SVG/PDF) & Table (1..9, S1..S7 MD/TeX) builder
├── build_manuscript.py                    # Template-interpolating manuscript/supplement compiler
├── render_pdf.py                          # Typeset HTML5 + headless Chrome + Ghostscript PDF renderer
├── templates/
│   ├── manuscript.tmpl.md                 # Full paper template with dynamic Jinja-style macros
│   └── supplement.tmpl.md                 # Supplementary note template with dynamic macros
├── figures/                               # Regenerated publication-ready vector figures
│   ├── fig{1..9}.svg                      # Vector graphics for web / digital viewing
│   └── fig{1..9}.pdf                      # Vector graphics for camera-ready print publication
├── tables/                                # Regenerated paired Markdown and LaTeX tables
│   ├── table{1..9}.{md,tex}               # Main Tables 1 to 9 (GFM Markdown and LaTeX)
│   └── supp_table{1..7}.{md,tex}          # Supplementary Tables S1 to S7 (GFM Markdown and LaTeX)
├── tests/                                 # Pytest test suite for end-to-end reproducibility
│   ├── test_aggregate_gcp_study.py
│   ├── test_collect_v110_compat_metrics.py
│   ├── test_build_figures_and_tables.py
│   └── test_paper_reproducibility.py
├── manuscript.md                          # Compiled ready-to-render main manuscript (Markdown)
├── manuscript.pdf                         # Compiled typeset main manuscript (PDF with vector figures)
├── supplement.md                          # Compiled ready-to-render supplementary material (Markdown)
└── supplement.pdf                         # Compiled typeset supplementary material (PDF)
```

---

## 4. Mapping of Figures, Tables, and Study Bundles

`Figures 1–8` and `Tables 1–8` (`Tables S1–S6`) derive from the 60 cloud study cells on identical GCP Intel Xeon Platinum 8373C Ice Lake instances (`n2-standard-32`, `benchmarks/gcp-study-v106/`), while `Figure 9`, `Table 9`, and `Table S7` derive from the host-independent deterministic `v1.1.0` clinical/validator bundle (`benchmarks/v110-compat-study/`):

| Target Artifact | Title / Scope | Primary Study Bundle / Cells | Experimental Stage / Domain |
| :--- | :--- | :--- | :--- |
| **Figure 1** | Study Architecture & 60-Cell Experimental Design | All 60 directories in `cells_60_inventory.tsv` | Stages 1–6 + Option C sweeps |
| **Figure 2** | Algorithmic Parity (H1/H1a) & vardictcpp Replication | `s1_*`, `s2_tierp_*`, `sweep_*` | Stage 1, Stage 2, Sweeps |
| **Figure 3** | Germline & Somatic Truth Accuracy (H3) & Titration | `s2_fig2_*`, `s4_tier3_*`, `s5_grch38_*`, `s6_low_vaf_*` | Stage 2, Stage 4, Stage 5, Stage 6 |
| **Figure 4** | Bounded Memory & Multi-Core Scaling (H2) | `s1_ceiling_*`, `s2_table3_*`, `s3_rna_*`, `s4_scale_*`, `s4_span_*`, `s4_perf_*`, `s4_jobs_*` | Stage 1, Stage 2, Stage 3, Stage 4 |
| **Figure 5** | Full-Genome chr1..chr22 Germline & Somatic WGS Scaling | `s5_wgs_hg002_*`, `s5_wgs_hg008_*` | Stage 5 (WGS Arms A & B) |
| **Figure 6** | 10-Platform Sequencing Chemistry & GA4GH v3.5 | `s5_cross_platform_*`, `s6_ga4gh_v35_*` | Stage 5, Stage 6 |
| **Figure 7** | Next-Gen Truth (v5.0Q, CMRG, Mosaic) & HG008 Landscape | `s4_tier3_*`, `s6_hg002_v50q_*`, `s6_hg008_modern_*` | Stage 4, Stage 6 |
| **Figure 8** | Supply Chain, Determinism, CLI Parity & CPP-1..12 | `s2_tierp_mode_*`, `s4_scale_*`, `s5_h4_*` | Stage 2, Stage 4, Stage 5 |
| **Figure 9** | `v1.1.0` Opt-In `--compat=off` Schema, Diploid `SEQC2` & `--preset` | `benchmarks/v110-compat-study/v110_compat_metrics.json` | Deterministic `v1.1.0` Domains 1–5 |
| **Table 1** | Algorithmic Parity (H1/H1a) & Option C Summary | `s1_*`, `s2_tierp_*`, `s3_rna_*`, `s4_*`, `sweep_*` | Stage 1, Stage 2, Stage 3, Stage 4, Sweeps |
| **Table 2** | vardictcpp Direct Replication & Table 3 -R chr1 | `s2_fig1_*`, `s2_fig2_*`, `s2_table3_*`, `s2_tierp_*` | Stage 2 |
| **Table 3** | Germline/Somatic Truth, Native GRCh38 & Titration | `s2_fig2_*`, `s4_tier3_*`, `s5_grch38_*`, `s6_low_vaf_*` | Stage 2, Stage 4, Stage 5, Stage 6 |
| **Table 4** | Bounded Memory (H2), Multi-Core Scaling & Span | `s2_table3_*`, `s3_rna_*`, `s4_jobs_*`, `s4_perf_*`, `s4_scale_*`, `s4_span_*` | Stage 2, Stage 3, Stage 4 |
| **Table 5** | Full-Genome Germline & 364x Somatic WGS Scaling | `s5_wgs_hg002_*`, `s5_wgs_hg008_*` | Stage 5 (WGS Arms A & B) |
| **Table 6** | 10-Platform Chemistry Generalization & GA4GH v3.5 | `s5_cross_platform_*`, `s6_ga4gh_v35_*` | Stage 5, Stage 6 |
| **Table 7** | Next-Gen Truth Benchmarks & HG008 Landscape | `s4_tier3_*`, `s6_hg002_v50q_*`, `s6_hg008_modern_*` | Stage 4, Stage 6 |
| **Table 8** | Supply-Chain Security (H4), Static musl & Divergences | `s2_tierp_mode_*`, `s5_h4_maintenance_security_matrix` | Stage 2, Stage 5 |
| **Table 9** | `v1.1.0` `--compat=off`, `--preset`, & Validator/Truth Summary | `benchmarks/v110-compat-study/v110_compat_metrics.json` | Deterministic `v1.1.0` Domains 1–6 |
| **Table S1** | Complete 60-Directory GCP Study Inventory | All 60 directories in `cells_60_inventory.tsv` | All 6 Stages + Sweeps |
| **Table S2** | Complete 87-Arm Option C Full-Core Parity Sweep | `sweep_*` (12 sweep directories, 87 evaluated arms) | Option C Full-Core Sweeps |
| **Table S3** | Complete Technical Replicate (n=3) Wall-Clock & RSS | `s2_tierp_*`, `s4_scale_*`, `s4_span_*`, `s4_perf_*`, etc. | Technical Replicates (n=3) |
| **Table S4** | QUAL ROC Sweep & 18-Arm Somatic Noise Floor | `s4_tier3_*`, `s5_cross_platform_*` (Part B) | Stage 4, Stage 5 |
| **Table S5** | Complete 70-Cell GA4GH v3.5 Genomic Stratification | `s6_ga4gh_v35_stratification_matrix` | Stage 6 |
| **Table S6** | 13-Mode Functional Parity & CPP-1..12 Catalog | `s2_tierp_mode_*`, `s5_h4_maintenance_security_matrix` | Stage 2, Stage 5 |
| **Table S7** | Complete 35-Entry `CompatProfile` (32 Switches) & CLI Catalog | `benchmarks/v110-compat-study/v110_compat_metrics.json` | Deterministic `v1.1.0` Domains 1–6 |

---

## 5. Verification Commands

To run individual test and verification suites without executing the entire pipeline:

```bash
# 1. Verify committed v1.1.0 deterministic study bundle & 32-switch config.rs alignment:
.venv/bin/python paper/collect_v110_compat_metrics.py --verify

# 2. Run pytest suite across all paper generators, carousels, and docs UI:
.venv/bin/python -m pytest paper/tests/ harness/tests/test_linkedin_carousel.py harness/tests/test_docs_ui.py -v

# 3. Check release version and manifest consistency:
python3 .agents/skills/vardictrust-release/scripts/check_release_consistency.py --phase pre-tag --require-head-fingerprint

# 4. Check dry-run GCP study runner plan and local caching:
bash harness/gcp_paper_runner.sh --plan-check

# 5. Verify release documentation surfaces match evidence manifests:
python3 harness/sync_evidence_docs.py --check
```

`reproduce_all.sh` requires a headless Chrome/Chromium binary to render PDFs. If none is found it exits before compiling (no stale PDFs reported as reproduced). For Chromeless environments (e.g. docs-only CI), `PAPER_ALLOW_STALE_PDF=1 bash paper/reproduce_all.sh` reuses the checked-in PDFs with a loud warning; Markdown/evidence gates still run.

---

## 6. Preprint (`bioRxiv`) & Journal (*Genome Biology* / *Bioinformatics*) Submission Guide

### 6.1 Uploadable Asset Bundle
All camera-ready and editorial-portal submission artifacts are generated deterministically by `bash paper/reproduce_all.sh`:
- **Typeset Main Manuscript & Supplement PDFs (`bioRxiv` & Initial Peer Review):** [`paper/manuscript.pdf`](manuscript.pdf) (single-file typeset PDF with inline vector figures `Figure 1`–`Figure 9` and `Table 1`–`Table 9`) and [`paper/supplement.pdf`](supplement.pdf) (`Supplementary Methods 1.0–1.8`, `Supplementary Notes 1–3`, and `Supplementary Tables S1–S7`).
- **Standalone Camera-Ready Vector Figures (`300+` DPI / TrueType Type-42 PDF & Searchable SVG):** [`paper/figures/fig1.pdf`](figures/fig1.pdf)–[`paper/figures/fig9.pdf`](figures/fig9.pdf) and [`paper/figures/fig1.svg`](figures/fig1.svg)–[`paper/figures/fig9.svg`](figures/fig9.svg) (Okabe-Ito colorblind-safe palette).
- **Standalone LaTeX & Markdown Tables:** [`paper/tables/table1.tex`](tables/table1.tex)–[`paper/tables/table9.tex`](tables/table9.tex) and [`paper/tables/supp_table1.tex`](tables/supp_table1.tex)–[`paper/tables/supp_table7.tex`](tables/supp_table7.tex) (`booktabs` + `\resizebox{\textwidth}{!}{...}` ready for Springer Nature `sn-jnl` or Oxford `bioinfo` LaTeX templates).

### 6.2 Manuscript Structure & Formatting Metrics
- **Structured Abstract:** Strictly `<= 350` words across `Background`, `Results`, and `Conclusions` (*Genome Biology* Software article limit `250–350` words), followed by `10` indexing keywords and a `5`-bullet **Key Points** box.
- **Display Items:** `9` Main Figures (`Figures 1–9`), `9` Main Tables (`Tables 1–9`), `7` Supplementary Tables (`Tables S1–S7`), and `27` numbered, DOI-verified references.
- **Preprint Licensing & Category (`bioRxiv`):** Subject area **Bioinformatics**; recommended license **CC BY 4.0** (matching permissive `MIT` software licensing in [`LICENSE`](../LICENSE)).

### 6.3 Portal Metadata, Disclosures, & Data Availability Checklist
- **Corresponding Author:** Palak Sheth (`palakpsheth@gmail.com`; [LinkedIn](https://www.linkedin.com/in/palak-p-sheth/); attach ORCID iD in Springer Nature Snapp / Oxford ScholarOne / `bioRxiv` portal metadata).
- **CRediT Taxonomy & Competing Interests:** Sole-author CRediT roles (`Conceptualization` through `Visualization`), self-funded cloud compute disclosure (`No external funding`), and `No competing interests` statement included in [`manuscript.md`](manuscript.md) under `## Declarations`.
- **AI-Assisted Development & Transparency Disclosure:** Included under `### AI-assisted development and transparency disclosure` in [`manuscript.md`](manuscript.md) per ICMJE, Springer Nature/*Genome Biology*, Oxford University Press, and [rewrites.bio](https://rewrites.bio) policies (disclosing `Google Gemini / Antigravity CLI` and `Anthropic Claude`, non-authorship of AI tools, full human responsibility, and automated differential/workspace test gates).
- **Software & Public Data Accessions:**
  - **Software (`MIT` License, SLSA Level 3):** `https://github.com/palakpsheth/VarDictRust` (tag `v1.1.0`; archive a permanent Zenodo software DOI snapshot upon preprint posting and paste the `10.5281/zenodo.*` DOI into portal metadata).
  - **Public Sequencing Benchmarks:** NCBI SRA `PRJNA744870` (`SRR15006375`, `SRR15006376`, `SRR15006386`, `SRR15006540`), `PRJNA523380` (`SRR8657348`), and `PRJNA506071` (`SEQC2 HCC1395/HCC1395BL`); NIST GIAB FTP (`HG002` `v4.2.1`, `v5.0q T2T-Q100`, `CMRG v1.00`, `Mosaic v1.10`; `HG005`; `HG008` `v0.3`/`v0.4`); and `ENA SRP295113` + Figshare (`doi:10.6084/m9.figshare.19128005.v2`, `doi:10.6084/m9.figshare.13511829.v4`) for `QIAseq V3` and Illumina paired-end oncopanel benchmarks.

### 6.4 Suggested Cover-Letter Highlights (*Genome Biology* / *Bioinformatics*)
1. **Drop-In Pipeline Consolidation Without Re-Validation Friction:** `VarDictRust` consolidates `VarDictJava -> R -> Perl` (`121` Conda packages) into a single `3.46 MB` `#![forbid(unsafe_code)]` static Rust binary that preserves `100%` Tier-1/2 numerical and allele equivalence (`ε = 1e-6`) across `13/13` operational modes and `10` short-read sequencing chemistries by default (`--compat=vardictjava`).
2. **Bounded $O(\text{tile} + \text{halo})$ Whole-Genome Scaling:** Eliminates whole-chromosome JVM out-of-memory crashes (`38.2×` faster and `110.2×` leaner than `VarDictJava -Xmx200g` on contiguous `chr1`), completing unchunked `130.8 GB` `HG002` and `329.0 GB` `364×` `HG008` somatic WGS on a `128 GiB` node and recovering `+2.46%` SNV / `+2.43%` indel recall on `NIST GIAB v5.0q T2T-Q100` (`+4.18%` / `+5.74%` on `CMRG v1.00`).
3. **Governed Opt-In Clinical Correctness (`v1.1.0`):** Isolates `32` upstream VCFv4.2/v4.3 schema, SV, somatic multi-allelic, adapter-chimera, and `--eqx` fixes behind `--compat=off` plus `--preset` (`wgs`, `wes`, `amplicon`, `ctdna`) and `-sp`, eliminating `11,612` strict VCF schema defects across `19` callsets and recovering `+31.7%` subclonal `SEQC2` somatic recall (`+0.1538` optimal $F_1$) in diploid `rtg vcfeval` without `--squash-ploidy`.
