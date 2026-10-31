# VarDictRust v1.0.6 Paper Reproducibility Package

This package provides a deterministic, single-command reproducibility pipeline for all figures, tables, text statistics, and LinkedIn dissemination assets accompanying the VarDictRust evaluation manuscript:

> **"Scaling Somatic and Germline Variant Calling with Byte Parity: An Empirical Study of VarDictRust Across 60 Cloud Workloads, 10 Sequencing Chemistries, and Full Human Genomes"**

---

## 1. Quickstart: 100% Offline Reproduction

All aggregated metrics, per-cell runtime/RSS measurements, truth evaluations, and cryptographic digests from the 60 Google Cloud Platform (GCP) study cells are committed in-tree under [`benchmarks/gcp-study-v106/`](../benchmarks/gcp-study-v106/).

To reproduce every figure, table, compiled manuscript, LinkedIn carousel slide, and verification gate from scratch:

```bash
# From the repository root:
bash paper/reproduce_all.sh
```

This single command automatically:
1. Bootstraps a local Python virtual environment (`.venv`) using `uv` with dependencies from [`paper/requirements.txt`](requirements.txt).
2. Verifies the SHA-256 hashes of all committed study synthesis tables against [`benchmarks/gcp-study-v106/SHA256SUMS.tsv`](../benchmarks/gcp-study-v106/SHA256SUMS.tsv).
3. Regenerates all **8 Main Figures** ([`paper/figures/fig1..fig8.{svg,pdf}`](figures/)) and all **14 Tables** ([`paper/tables/table1..table8.{md,tex}`](tables/) and [`paper/tables/supp_table1..supp_table6.{md,tex}`](tables/)).
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
├── README.md                      # This reproduction guide
├── requirements.txt               # Pinned Python dependencies (matplotlib, pytest, etc.)
├── reproduce_all.sh               # Master one-command offline & online runner
├── aggregate_gcp_study.py         # 60-cell JSON/TSV metric aggregator & validator
├── build_figures_and_tables.py    # Figure (1..8 SVG/PDF) & Table (1..8, S1..S6 MD/TeX) builder
├── build_manuscript.py            # Template-interpolating manuscript/supplement compiler
├── render_pdf.py                  # Typeset HTML5 + headless Chrome + Ghostscript PDF renderer
├── templates/
│   ├── manuscript.tmpl.md         # Full paper template with dynamic Jinja-style macros
│   └── supplement.tmpl.md         # Supplementary note template with dynamic macros
├── figures/                       # Regenerated publication-ready vector figures
│   ├── fig{1..8}.svg              # Vector graphics for web / digital viewing
│   └── fig{1..8}.pdf              # Vector graphics for camera-ready print publication
├── tables/                        # Regenerated paired Markdown and LaTeX tables
│   ├── table{1..8}.{md,tex}       # Main Tables 1 to 8 (GFM Markdown and LaTeX)
│   └── supp_table{1..6}.{md,tex}  # Supplementary Tables S1 to S6 (GFM Markdown and LaTeX)
├── tests/                         # Pytest test suite for end-to-end reproducibility
│   ├── test_aggregate_gcp_study.py
│   ├── test_build_figures_and_tables.py
│   └── test_paper_reproducibility.py
├── manuscript.md                  # Compiled ready-to-render main manuscript (Markdown)
├── manuscript.pdf                 # Compiled typeset main manuscript (PDF with vector figures)
├── supplement.md                  # Compiled ready-to-render supplementary material (Markdown)
└── supplement.pdf                 # Compiled typeset supplementary material (PDF)
```

---

## 4. Mapping of Figures, Tables, and GCP Study Cells

The 60 cloud study cells span 6 experimental stages plus multi-core Option C parameter sweeps on identical GCP Intel Xeon Platinum 8373C Ice Lake instances (`n2-standard-32`):

| Target Artifact | Title / Scope | Primary GCP Study Cells | Experimental Stage |
| :--- | :--- | :--- | :--- |
| **Figure 1** | Study Architecture & 60-Cell Experimental Design | All 60 directories in `cells_60_inventory.tsv` | Stages 1–6 + Option C sweeps |
| **Figure 2** | Algorithmic Parity (H1/H1a) & vardictcpp Replication | `s1_*`, `s2_tierp_*`, `sweep_*` | Stage 1, Stage 2, Sweeps |
| **Figure 3** | Germline & Somatic Truth Accuracy (H3) & Titration | `s2_fig2_*`, `s4_tier3_*`, `s5_grch38_*`, `s6_low_vaf_*` | Stage 2, Stage 4, Stage 5, Stage 6 |
| **Figure 4** | Bounded Memory & Multi-Core Scaling (H2) | `s1_ceiling_*`, `s2_table3_*`, `s3_rna_*`, `s4_scale_*`, `s4_span_*`, `s4_perf_*`, `s4_jobs_*` | Stage 1, Stage 2, Stage 3, Stage 4 |
| **Figure 5** | Full-Genome chr1..chr22 Germline & Somatic WGS Scaling | `s5_wgs_hg002_*`, `s5_wgs_hg008_*` | Stage 5 (WGS Arms A & B) |
| **Figure 6** | 10-Platform Sequencing Chemistry & GA4GH v3.5 | `s5_cross_platform_*`, `s6_ga4gh_v35_*` | Stage 5, Stage 6 |
| **Figure 7** | Next-Gen Truth (v5.0Q, CMRG, Mosaic) & HG008 Landscape | `s4_tier3_*`, `s6_hg002_v50q_*`, `s6_hg008_modern_*` | Stage 4, Stage 6 |
| **Figure 8** | Supply Chain, Determinism, CLI Parity & CPP-1..10 | `s2_tierp_mode_*`, `s4_scale_*`, `s5_h4_*` | Stage 2, Stage 4, Stage 5 |
| **Table 1** | Algorithmic Parity (H1/H1a) & Option C Summary | `s1_*`, `s2_tierp_*`, `s3_rna_*`, `s4_*`, `sweep_*` | Stage 1, Stage 2, Stage 3, Stage 4, Sweeps |
| **Table 2** | vardictcpp Direct Replication & Table 3 -R chr1 | `s2_fig1_*`, `s2_fig2_*`, `s2_table3_*`, `s2_tierp_*` | Stage 2 |
| **Table 3** | Germline/Somatic Truth, Native GRCh38 & Titration | `s2_fig2_*`, `s4_tier3_*`, `s5_grch38_*`, `s6_low_vaf_*` | Stage 2, Stage 4, Stage 5, Stage 6 |
| **Table 4** | Bounded Memory (H2), Multi-Core Scaling & Span | `s2_table3_*`, `s3_rna_*`, `s4_jobs_*`, `s4_perf_*`, `s4_scale_*`, `s4_span_*` | Stage 2, Stage 3, Stage 4 |
| **Table 5** | Full-Genome Germline & 364x Somatic WGS Scaling | `s5_wgs_hg002_*`, `s5_wgs_hg008_*` | Stage 5 (WGS Arms A & B) |
| **Table 6** | 10-Platform Chemistry Generalization & GA4GH v3.5 | `s5_cross_platform_*`, `s6_ga4gh_v35_*` | Stage 5, Stage 6 |
| **Table 7** | Next-Gen Truth Benchmarks & HG008 Landscape | `s4_tier3_*`, `s6_hg002_v50q_*`, `s6_hg008_modern_*` | Stage 4, Stage 6 |
| **Table 8** | Supply-Chain Security (H4), Static musl & Divergences | `s2_tierp_mode_*`, `s5_h4_maintenance_security_matrix` | Stage 2, Stage 5 |
| **Table S1** | Complete 60-Directory GCP Study Inventory | All 60 directories in `cells_60_inventory.tsv` | All 6 Stages + Sweeps |
| **Table S2** | Complete 87-Arm Option C Full-Core Parity Sweep | `sweep_*` (12 sweep directories, 87 evaluated arms) | Option C Full-Core Sweeps |
| **Table S3** | Complete Technical Replicate (n=3) Wall-Clock & RSS | `s2_tierp_*`, `s4_scale_*`, `s4_span_*`, `s4_perf_*`, etc. | Technical Replicates (n=3) |
| **Table S4** | QUAL ROC Sweep & 18-Arm Somatic Noise Floor | `s4_tier3_*`, `s5_cross_platform_*` (Part B) | Stage 4, Stage 5 |
| **Table S5** | Complete 70-Cell GA4GH v3.5 Genomic Stratification | `s6_ga4gh_v35_stratification_matrix` | Stage 6 |
| **Table S6** | 13-Mode Functional Parity & CPP-1..10 Catalog | `s2_tierp_mode_*`, `s5_h4_maintenance_security_matrix` | Stage 2, Stage 5 |

---

## 5. Verification Commands

To run individual test and verification suites without executing the entire pipeline:

```bash
# 1. Run pytest suite across all paper generators, carousels, and docs UI:
.venv/bin/pytest paper/tests/ harness/tests/test_linkedin_carousel.py harness/tests/test_docs_ui.py -v

# 2. Check release version and manifest consistency:
python3 .agents/skills/vardictrust-release/scripts/check_release_consistency.py --version 1.0.6 --phase post-bioconda --require-head-fingerprint

# 3. Check dry-run GCP study runner plan and local caching:
bash harness/gcp_paper_runner.sh --plan-check

# 4. Verify release documentation surfaces match evidence manifests:
python3 harness/sync_evidence_docs.py --check
```
