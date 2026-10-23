#!/usr/bin/env bash
# Usage: EST_MIN=<minutes> courtesy.sh <cmd...>   — decision A3: never touch services.
set -euo pipefail
: "${EST_MIN:?courtesy: set EST_MIN (estimated runtime, minutes)}"
free=$(df -BG --output=avail "${COURTESY_DIR:-.}" | tail -1 | tr -dc 0-9)
[ "$free" -ge "${MIN_FREE_GB:-150}" ] || { echo "courtesy: ${free}G free < ${MIN_FREE_GB:-150}G" >&2; exit 3; }
avail=$(( $(awk '/MemAvailable/{print $2}' /proc/meminfo) / 1048576 ))
[ "$avail" -ge "${MIN_AVAIL_GB:-12}" ] || { echo "courtesy: ${avail}G RAM available < ${MIN_AVAIL_GB:-12}G" >&2; exit 5; }
hhmm=${COURTESY_NOW:-$(date +%H%M)}
now=$((10#${hhmm:0:2} * 60 + 10#${hhmm:2:2})); end=$((now + EST_MIN)); ws=220; we=315
for off in 0 1440; do
  if [ "$now" -lt $((we + off)) ] && [ "$end" -gt $((ws + off)) ]; then
    echo "courtesy: [$hhmm +${EST_MIN}m] overlaps 03:40-05:15 maintenance; start after 05:15 or split" >&2; exit 4
  fi
done
export COURTESY_ACTIVE=1
exec nice -n 19 ionice -c 3 "$@"
