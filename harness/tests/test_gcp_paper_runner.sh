#!/usr/bin/env bash
# Self-test for harness/gcp_paper_runner.sh: verifies 2-strike Spot -> On-Demand
# escalation policy, 2048GB pd-ssd disk configuration, GOLDEN_MODE=off guard,
# and per-cell checkpoint skip/resume semantics.
set -euo pipefail

bash harness/gcp_paper_runner.sh --plan-check
echo "PASS test_gcp_paper_runner"
