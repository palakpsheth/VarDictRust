#!/usr/bin/env bash
# Self-test for harness/gcp_paper_runner.sh: verifies 2-strike Spot -> On-Demand
# escalation policy, 2048GB pd-ssd disk configuration, GOLDEN_MODE=off guard,
# and per-cell checkpoint skip/resume semantics.
set -euo pipefail

bash harness/gcp_paper_runner.sh --plan-check

# Stage validation & supply-chain hardening checks (with mock gcloud/gsutil blocking any live GCP calls)
MOCK_GCP_BIN="$(mktemp -d)"
trap 'rm -rf "$MOCK_GCP_BIN"' EXIT
for tool in gcloud gsutil; do
  printf '#!/usr/bin/env bash\necho "MOCK %s invoked unexpectedly" >&2\nexit 99\n' "$tool" > "$MOCK_GCP_BIN/$tool"
  chmod +x "$MOCK_GCP_BIN/$tool"
done
if out=$(PATH="$MOCK_GCP_BIN:$PATH" bash harness/lib/portable.sh --timeout 2s bash harness/gcp_paper_runner.sh --controller --stage "1;echo_pwn" 2>&1); then
  echo "FAIL test_gcp_paper_runner: --controller accepted invalid --stage" >&2
  exit 1
elif grep -q "MOCK " <<<"$out"; then
  echo "FAIL test_gcp_paper_runner: --controller reached gcloud/gsutil before validating --stage" >&2
  exit 1
fi
if out=$(PATH="$MOCK_GCP_BIN:$PATH" WORKER_STAGE="bad_stage" bash harness/lib/portable.sh --timeout 2s bash harness/gcp_paper_runner.sh --controller 2>&1); then
  echo "FAIL test_gcp_paper_runner: --controller accepted invalid WORKER_STAGE env" >&2
  exit 1
elif grep -q "MOCK " <<<"$out"; then
  echo "FAIL test_gcp_paper_runner: --controller reached gcloud/gsutil before validating WORKER_STAGE" >&2
  exit 1
fi
if out=$(PATH="$MOCK_GCP_BIN:$PATH" bash harness/lib/portable.sh --timeout 2s bash harness/gcp_paper_runner.sh --worker --stage "bad_stage" 2>&1); then
  echo "FAIL test_gcp_paper_runner: --worker accepted invalid --stage" >&2
  exit 1
elif grep -q "MOCK " <<<"$out"; then
  echo "FAIL test_gcp_paper_runner: --worker reached preflight/gsutil before validating --stage" >&2
  exit 1
fi
for forbidden in "Miniconda3-latest-Linux-x86_64.sh" "sh.rustup.rs" "export WORK_DIR='\${WORK_DIR}'" "export GCS_CKPT='\${GCS_CKPT}'"; do
  if grep -Fq "$forbidden" harness/gcp_paper_runner.sh; then
    echo "FAIL test_gcp_paper_runner: forbidden unpinned/unescaped construct in gcp_paper_runner.sh: $forbidden" >&2
    exit 1
  fi
done
grep -Fq "printf -v q_work_dir '%q'" harness/gcp_paper_runner.sh || {
  echo "FAIL test_gcp_paper_runner: missing printf %q escaping for WORK_DIR" >&2
  exit 1
}
grep -Fq "76cf791877078e74410f64e5c2ef7da4e217fd7a433c2c3674c9d9695225c3bb" harness/gcp_paper_runner.sh || {
  echo "FAIL test_gcp_paper_runner: missing pinned Miniconda SHA-256" >&2
  exit 1
}
grep -Fq "b1c8cb14e30889e7e2999e6b69e62d9e9b9ffc0f4df8180564d306d9680e6c55" harness/gcp_paper_runner.sh || {
  echo "FAIL test_gcp_paper_runner: missing pinned rustup-init SHA-256" >&2
  exit 1
}
grep -Fq "f36059985ed4b4cc1b5809d11750102d73a5fc16dde5e81ae0eea738412a4a24" harness/concordance_multimodal.sh || {
  echo "FAIL test_gcp_paper_runner: missing hg005 truth body SHA-256 in concordance_multimodal.sh" >&2
  exit 1
}
grep -Fq "43174566bca0425be270d110dd9c9897b5508eb06644ba5198fd7a5f185b33ff" harness/concordance_multimodal.sh || {
  echo "FAIL test_gcp_paper_runner: missing na24385_b37 truth body SHA-256 in concordance_multimodal.sh" >&2
  exit 1
}
grep -Fq "stage_verified_sha256" harness/run_stage_6.sh || {
  echo "FAIL test_gcp_paper_runner: run_stage_6.sh must route downloads through stage_verified_sha256" >&2
  exit 1
}

echo "PASS test_gcp_paper_runner"
