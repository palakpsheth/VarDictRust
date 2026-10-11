#!/usr/bin/env bash
# Usage: sysscope.sh <cmd...> — run a heavy job in a transient SYSTEM scope (system.slice) as the
# current user, outside user@UID.service, which systemd-oomd watches (it SIGKILLs the whole
# terminal scope, incl. agy, when that cgroup's memory pressure exceeds 50% for 20s).
# Transient only: no unit files or config are changed and the scope vanishes when the job ends.
set -euo pipefail
exec sudo -n --preserve-env=PATH,HOME,JAVA_OPTS,KEEP_DIR,MIN_VARIANTS,EST_MIN,CMP_STRICT,COURTESY_DIR,CONDA_PREFIX,CONDA_DEFAULT_ENV \
  systemd-run --scope --quiet --collect --slice=system.slice \
  -p MemoryMax="${SYSSCOPE_MEM_MAX:-14G}" -p MemorySwapMax="${SYSSCOPE_SWAP_MAX:-2G}" \
  --uid="$(id -u)" --gid="$(id -g)" \
  --setenv=PATH="$PATH" --setenv=HOME="$HOME" --setenv=JAVA_OPTS="${JAVA_OPTS:-}" \
  --setenv=KEEP_DIR="${KEEP_DIR:-}" --setenv=MIN_VARIANTS="${MIN_VARIANTS:-}" --setenv=EST_MIN="${EST_MIN:-}" \
  --setenv=CMP_STRICT="${CMP_STRICT:-}" \
  -- "$@"
