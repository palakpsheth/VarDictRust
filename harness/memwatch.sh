#!/usr/bin/env bash
# Usage: memwatch.sh [max_gib=4] — safety net: SIGKILL any vardictrust/VarDict process that runs in the
# user session (user@UID.service, which systemd-oomd watches) and exceeds max_gib RSS. Runs launched
# via harness/sysscope.sh live in system.slice with their own MemoryMax and are left alone.
MAX_KB=$(( ${1:-4} * 1048576 ))
while sleep 3; do
  ps -eo pid=,rss=,args= | awk -v m="$MAX_KB" '$2>m && ($0 ~ /target\/(debug|release)\/vardictrust/ || $0 ~ /VarDict.*-G /) {print $1, $2}' |
  while read -r pid rss; do
    grep -q 'user@' "/proc/$pid/cgroup" 2>/dev/null || continue
    echo "$(date -Is) memwatch: killing pid $pid rss=$((rss/1024))MiB (user-session run over limit)" | tee -a /tmp/memwatch.log
    kill -9 "$pid"
  done
done
