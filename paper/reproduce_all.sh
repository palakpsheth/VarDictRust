#!/usr/bin/env bash
# Reproduce all VarDictRust v1.0.6 cloud & v1.1.0 clinical paper figures, tables, manuscript/supplement PDFs, and LinkedIn carousel.
# Default mode is 100% offline using git-tracked benchmarks/gcp-study-v106/ and benchmarks/v110-compat-study/.
# Usage: bash paper/reproduce_all.sh [--sync-gcs] [--reaggregate-gcp] [--run-v110-eval] [--require-head-fingerprint]
#   --sync-gcs                  Re-pull raw cell outputs from gs://vardict-paper-trial-results/cells-v106/.
#   --reaggregate-gcp           Re-extract benchmarks/gcp-study-v106/ from local benchmarks/data/gcp-study-v106/ raw cells.
#   --run-v110-eval             Re-execute all 6 deterministic v1.1.0 --compat=off & --preset domains locally.
#   --require-head-fingerprint  Enforce live HEAD engine fingerprint parity in check_release_consistency.py
#                               (also enabled automatically on release/* branches/PRs and refs/tags/v*).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

export PATH="$HOME/.local/bin:$PATH"
UV_BIN="${UV_BIN:-$(command -v uv 2>/dev/null || echo "$HOME/.local/bin/uv")}"

SYNC_GCS=0
REAGGREGATE_GCP=0
RUN_V110_EVAL=0
REQUIRE_HEAD_FP=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --sync-gcs)
      SYNC_GCS=1
      ;;
    --reaggregate-gcp)
      REAGGREGATE_GCP=1
      ;;
    --run-v110-eval)
      RUN_V110_EVAL=1
      ;;
    --require-head-fingerprint)
      REQUIRE_HEAD_FP=1
      ;;
    *)
      echo "ERROR: unknown flag '$1'" >&2
      exit 2
      ;;
  esac
  shift
done

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

# Fail fast when no Chrome/Chromium renderer is available: otherwise the
# pipeline would silently reuse checked-in PDFs and report stale artifacts
# as reproduced. Chromeless environments may opt into tracked-PDF reuse with
# PAPER_ALLOW_STALE_PDF=1 (Markdown/evidence gates still run).
if [[ "${PAPER_ALLOW_STALE_PDF:-0}" != "1" ]]; then
  if ! command -v google-chrome >/dev/null 2>&1 && ! command -v chromium >/dev/null 2>&1 && ! command -v chromium-browser >/dev/null 2>&1 && [[ ! -x /usr/bin/google-chrome ]]; then
    echo "ERROR: no headless Chrome/Chromium binary found; cannot render manuscript/supplement PDFs." >&2
    echo "Install google-chrome or chromium, or run with PAPER_ALLOW_STALE_PDF=1 to reuse the checked-in PDFs (figures/text changes will NOT be re-rendered)." >&2
    exit 1
  fi
fi

if [[ "$SYNC_GCS" -eq 1 ]]; then
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
RAW_CELL_COUNT=0
if [[ -d "$RAW_CELLS_DIR" ]]; then
  RAW_CELL_COUNT="$(find "$RAW_CELLS_DIR" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | wc -l | tr -d '[:space:]')"
fi
if [[ "$SYNC_GCS" -eq 1 || "$REAGGREGATE_GCP" -eq 1 ]]; then
  if [[ "$RAW_CELL_COUNT" -lt 60 ]]; then
    echo "ERROR: explicit GCP reaggregation requested (--reaggregate-gcp / --sync-gcs), but '$RAW_CELLS_DIR' contains $RAW_CELL_COUNT raw cell directories (need >= 60)." >&2
    exit 1
  fi
  echo "Aggregating 60 raw GCP study cells from $RAW_CELLS_DIR..."
  "$PY" paper/aggregate_gcp_study.py --raw-dir "$RAW_CELLS_DIR" --out-dir benchmarks/gcp-study-v106
elif [[ ! -f "$ROOT/benchmarks/gcp-study-v106/SHA256SUMS.tsv" && "$RAW_CELL_COUNT" -ge 60 ]]; then
  echo "Aggregating 60 raw GCP study cells from $RAW_CELLS_DIR..."
  "$PY" paper/aggregate_gcp_study.py --raw-dir "$RAW_CELLS_DIR" --out-dir benchmarks/gcp-study-v106
else
  echo "Verifying committed study artifacts against benchmarks/gcp-study-v106/SHA256SUMS.tsv..."
  (cd benchmarks/gcp-study-v106 && awk 'NR>1 {print $1 "  " $2}' SHA256SUMS.tsv | sha256sum -c -)
  "$PY" paper/aggregate_gcp_study.py --verify
fi

if [[ "$RUN_V110_EVAL" -eq 1 ]]; then
  echo "Collecting deterministic v1.1.0 --compat=off & --preset metrics..."
  "$PY" paper/collect_v110_compat_metrics.py --run
else
  echo "Verifying committed v1.1.0 study bundle (benchmarks/v110-compat-study/)..."
  "$PY" paper/collect_v110_compat_metrics.py --verify
fi

echo "Regenerating 9 Main Figures, 9 Main Tables, and 7 Supplementary Tables..."
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
cur_branch="$(git symbolic-ref -q --short HEAD 2>/dev/null || true)"
if [[ "$REQUIRE_HEAD_FP" -eq 1 || "${GITHUB_HEAD_REF:-}" == release/* || "${GITHUB_REF:-}" == refs/tags/v* || "$cur_branch" == release/* ]]; then
  "$PY" .agents/skills/vardictrust-release/scripts/check_release_consistency.py \
    --phase pre-tag --require-head-fingerprint
else
  "$PY" .agents/skills/vardictrust-release/scripts/check_release_consistency.py \
    --phase pre-tag
fi
bash harness/gcp_paper_runner.sh --plan-check
"$PY" harness/sync_evidence_docs.py --check

echo "PASS: All 60 cells + 6 v1.1.0 domains, 9 figures, 16 tables, manuscript, supplement, and 10 carousel slides reproduced."
