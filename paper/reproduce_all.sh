#!/usr/bin/env bash
# Reproduce all VarDictRust v1.0.6 paper figures, tables, manuscript text, and LinkedIn carousel.
# Default mode is 100% offline using git-tracked benchmarks/gcp-study-v106/.
# Pass --sync-gcs to re-pull raw cell outputs from gs://vardict-paper-trial-results/cells-v106/.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

export PATH="$HOME/.local/bin:$PATH"
UV_BIN="${UV_BIN:-$(command -v uv 2>/dev/null || echo /home/palak/.local/bin/uv)}"

if [[ ! -x "$ROOT/.venv/bin/pytest" || ! -x "$ROOT/.venv/bin/python" ]]; then
  if [[ ! -x "$UV_BIN" ]]; then
    echo "ERROR: uv binary not found at '$UV_BIN'. Please install uv or set UV_BIN." >&2
    exit 1
  fi
  echo "Bootstrapping Python virtual environment via uv..."
  "$UV_BIN" venv "$ROOT/.venv"
  "$UV_BIN" pip install --quiet --python "$ROOT/.venv/bin/python" -r "$ROOT/paper/requirements.txt"
fi

PY="$ROOT/.venv/bin/python"

if [[ "${1:-}" == "--sync-gcs" ]]; then
  echo "Syncing raw GCP study cell results from GCS..."
  mkdir -p benchmarks/data/gcp-study-v106
  if command -v gcloud >/dev/null 2>&1; then
    gcloud storage rsync -r \
      -x '.*\.(vcf|vcf\.gz|tbi|csi|var|bam|bai|cram|crai|bed|bed\.gz|stdout|err|fa|fai|sdf)$|(^|/).*\.(java|rust|cpp|splice|sorted)\.(tsv|var)$|(^|/)[^/]*\.var\.diff\.txt$|(^|/)(strelka_.*|fb_chunks|.*_roc\.tsv\.gz)(/|$)' \
      gs://vardict-paper-trial-results/cells-v106/ \
      benchmarks/data/gcp-study-v106/
  elif command -v gsutil >/dev/null 2>&1; then
    gsutil -m -o GSUtil:parallel_process_count=1 -o GSUtil:parallel_thread_count=16 rsync -r \
      -x '.*\.(vcf|vcf\.gz|tbi|csi|var|bam|bai|cram|crai|bed|bed\.gz|stdout|err|fa|fai|sdf)$|(^|/).*\.(java|rust|cpp|splice|sorted)\.(tsv|var)$|(^|/)[^/]*\.var\.diff\.txt$|(^|/)(strelka_.*|fb_chunks|.*_roc\.tsv\.gz)(/|$)' \
      gs://vardict-paper-trial-results/cells-v106/ \
      benchmarks/data/gcp-study-v106/
  else
    echo "ERROR: neither 'gcloud' nor 'gsutil' is available to perform --sync-gcs." >&2
    exit 1
  fi
fi

RAW_CELLS_DIR="$ROOT/benchmarks/data/gcp-study-v106"
if [[ -d "$RAW_CELLS_DIR" ]] && [[ $(find "$RAW_CELLS_DIR" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | wc -l) -ge 60 ]]; then
  echo "Aggregating 60 raw GCP study cells from $RAW_CELLS_DIR..."
  "$PY" paper/aggregate_gcp_study.py --raw-dir "$RAW_CELLS_DIR" --out-dir benchmarks/gcp-study-v106
else
  echo "Verifying committed study artifacts against benchmarks/gcp-study-v106/SHA256SUMS.tsv..."
  (cd benchmarks/gcp-study-v106 && awk 'NR>1 {print $1 "  " $2}' SHA256SUMS.tsv | sha256sum -c -)
  "$PY" paper/aggregate_gcp_study.py --verify
fi

echo "Regenerating 8 Main Figures, 8 Main Tables, and 6 Supplementary Tables..."
"$PY" paper/build_figures_and_tables.py --input benchmarks/gcp-study-v106/master_gcp_metrics.json
"$PY" paper/build_figures_and_tables.py --verify

echo "Compiling manuscript.md and supplement.md..."
"$PY" paper/build_manuscript.py
"$PY" paper/build_manuscript.py --check

echo "Generating LinkedIn carousel slides..."
"$PY" harness/generate_linkedin_carousel.py \
  --evidence-dir benchmarks/evidence \
  --gcp-metrics benchmarks/gcp-study-v106/master_gcp_metrics.json \
  --out-dir scratch/linkedin

echo "Running paper, carousel, release-consistency, docs-UI, and evidence-sync gates..."
"$PY" -m pytest paper/tests/ harness/tests/test_linkedin_carousel.py harness/tests/test_docs_ui.py -v
"$PY" .agents/skills/vardictrust-release/scripts/check_release_consistency.py \
  --version 1.0.6 --phase post-bioconda --require-head-fingerprint
bash harness/gcp_paper_runner.sh --plan-check
"$PY" harness/sync_evidence_docs.py --check

echo "PASS: All 60 cells, 8 figures, 14 tables, manuscript, supplement, and 10 carousel slides reproduced."
