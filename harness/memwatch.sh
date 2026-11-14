#!/usr/bin/env bash
# Usage: memwatch.sh [max_gib=4] — safety net: SIGKILL any vardictrust/VarDict process that runs in the
# user session (user@UID.service, which systemd-oomd watches) and exceeds max_gib RSS. Runs launched
# via harness/sysscope.sh live in system.slice with their own MemoryMax and are left alone.
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib/portable.sh"
MAX_KB=$(( ${1:-4} * 1048576 ))
LOG="${MEMWATCH_LOG:-${XDG_RUNTIME_DIR:-$HOME/.local/state}/vardictrust/memwatch.log}"
mkdir -p -m 700 "$(dirname "$LOG")"
while sleep 3; do
  ps -u "$(id -u)" -o pid=,rss=,args= | awk -v m="$MAX_KB" '$2>m && ($0 ~ /target\/(debug|release)\/vardictrust/ || $0 ~ /VarDict.*-G /) {print $1, $2}' |
  while read -r pid rss; do
    is_user_session_pid "$pid" || continue
    mkdir -p -m 700 "$(dirname "$LOG")"
    echo "$(date -u '+%Y-%m-%dT%H:%M:%SZ') memwatch: killing pid $pid rss=$((rss/1024))MiB (user-session run over limit)" | tee -a "$LOG"
    kill -9 "$pid"
  done
done
