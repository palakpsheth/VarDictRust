#!/usr/bin/env bash
set -euo pipefail

# Case 1: Overlaps maintenance (03:00 + 60 min reaches 04:00, overlapping 03:40-05:15) -> exit 4
set +e
MIN_FREE_GB=0 MIN_AVAIL_GB=0 COURTESY_NOW=0300 EST_MIN=60 bash harness/courtesy.sh true >/dev/null 2>&1
rc=$?
set -e
[ "$rc" -eq 4 ] || { echo "FAIL test_courtesy: expected 4 for 0300+60m, got $rc" >&2; exit 1; }

# Case 2: Overlaps maintenance crossing midnight (23:00 + 300 min reaches 04:00) -> exit 4
set +e
MIN_FREE_GB=0 MIN_AVAIL_GB=0 COURTESY_NOW=2300 EST_MIN=300 bash harness/courtesy.sh true >/dev/null 2>&1
rc=$?
set -e
[ "$rc" -eq 4 ] || { echo "FAIL test_courtesy: expected 4 for 2300+300m, got $rc" >&2; exit 1; }

# Case 3: Outside maintenance window, normal execution -> exit 0
MIN_FREE_GB=0 MIN_AVAIL_GB=0 COURTESY_NOW=0600 EST_MIN=60 bash harness/courtesy.sh true

# Case 4: Insufficient disk space -> exit 3
set +e
MIN_FREE_GB=999999 COURTESY_NOW=0600 EST_MIN=60 bash harness/courtesy.sh true >/dev/null 2>&1
rc=$?
set -e
[ "$rc" -eq 3 ] || { echo "FAIL test_courtesy: expected 3 for MIN_FREE_GB, got $rc" >&2; exit 1; }

# Case 5: Insufficient available RAM -> exit 5
set +e
MIN_AVAIL_GB=999999 COURTESY_NOW=0600 EST_MIN=60 bash harness/courtesy.sh true >/dev/null 2>&1
rc=$?
set -e
[ "$rc" -eq 5 ] || { echo "FAIL test_courtesy: expected 5 for MIN_AVAIL_GB, got $rc" >&2; exit 1; }

# Case 6: SKIP_BLACKOUT=1 with EST_MIN overlapping blackout (inside window) skips exit 4
SKIP_BLACKOUT=1 MIN_FREE_GB=0 MIN_AVAIL_GB=0 COURTESY_NOW=0400 EST_MIN=60 bash harness/courtesy.sh true

# Case 7: SKIP_BLACKOUT=1 still enforces floors -> exit 3
set +e
SKIP_BLACKOUT=1 MIN_FREE_GB=999999 COURTESY_NOW=0400 EST_MIN=60 bash harness/courtesy.sh true >/dev/null 2>&1
rc=$?
set -e
[ "$rc" -eq 3 ] || { echo "FAIL test_courtesy: expected 3 for SKIP_BLACKOUT floors, got $rc" >&2; exit 1; }

echo "PASS test_courtesy"
