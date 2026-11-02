#!/usr/bin/env bash
# harness/lib/portable.sh — Cross-platform (GNU/Linux + macOS BSD) shell helpers.

PORTABLE_SH_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
export PORTABLE_SH_PATH

_CONDA_ROOT="${CONDA_ROOT:-${CONDA_PREFIX:-$HOME/miniconda3}}"
_CONDA_ROOT="${_CONDA_ROOT%/envs/*}"
G4_TOOLS_BIN="${G4_TOOLS_BIN:-$_CONDA_ROOT/envs/g4-tools/bin}"
unset _CONDA_ROOT
VARDICT_PROFILES_DIR="${VARDICT_PROFILES_DIR:-$HOME/Documents/vardictrust_profiles}"
export G4_TOOLS_BIN VARDICT_PROFILES_DIR

ensure_g4_tools_path() {
    if [ "${NO_G4_TOOLS:-0}" = 0 ] && [ -d "$G4_TOOLS_BIN" ]; then
        case ":$PATH:" in
            *":$G4_TOOLS_BIN:"*) ;;
            *) export PATH="$G4_TOOLS_BIN:$PATH" ;;
        esac
    fi
}

find_gitleaks() {
    if command -v gitleaks >/dev/null 2>&1; then
        command -v gitleaks
        return 0
    elif [ -x "$G4_TOOLS_BIN/gitleaks" ]; then
        printf '%s\n' "$G4_TOOLS_BIN/gitleaks"
        return 0
    fi
    return 1
}

now_ms() {
    if command -v python3 >/dev/null 2>&1; then
        python3 -c 'import time; print(int(time.time() * 1000))'
    elif command -v perl >/dev/null 2>&1; then
        perl -MTime::HiRes=time -e 'printf "%d\n", int(time() * 1000)'
    elif command -v gdate >/dev/null 2>&1; then
        gdate +%s%3N
    else
        printf '%d000\n' "$(date +%s)"
    fi
}

portable_time_v() {
    if [ -x /usr/bin/time ] && /usr/bin/time -v true >/dev/null 2>&1; then
        /usr/bin/time -v "$@"
        return $?
    elif command -v gtime >/dev/null 2>&1 && gtime -v true >/dev/null 2>&1; then
        gtime -v "$@"
        return $?
    fi

    local append=0
    local out_file=""
    while [ $# -gt 0 ]; do
        case "$1" in
            -a)
                append=1
                shift
                ;;
            -o)
                out_file="$2"
                shift 2
                ;;
            --)
                shift
                break
                ;;
            -*)
                shift
                ;;
            *)
                break
                ;;
        esac
    done

    local t0_ms t1_ms elapsed_ms elapsed_sec elapsed_min elapsed_rem
    t0_ms="$(now_ms)"
    local tmp_time
    tmp_time="$(mktemp "${TMPDIR:-/tmp}/portable_time.XXXXXX")"
    local rc=0
    if [ -x /usr/bin/time ] && /usr/bin/time -l true >/dev/null 2>&1; then
        { /usr/bin/time -l sh -c 'exec "$@" 2>&3 3>&-' sh "$@" 2>"$tmp_time"; } 3>&2 || rc=$?
    elif command -v python3 >/dev/null 2>&1; then
        python3 -c '
import resource, subprocess, sys
tmp_time = sys.argv[1]
proc = subprocess.run(sys.argv[2:])
ru = resource.getrusage(resource.RUSAGE_CHILDREN)
with open(tmp_time, "w", encoding="utf-8") as fh:
    fh.write(f"  {int(ru.ru_maxrss)}  maximum resident set size\n")
sys.exit(proc.returncode)
' "$tmp_time" "$@" || rc=$?
    else
        rm -f "$tmp_time"
        echo "portable_time_v: no supported timing backend (/usr/bin/time, gtime, or python3)" >&2
        return 127
    fi
    t1_ms="$(now_ms)"
    elapsed_ms=$(( t1_ms - t0_ms ))
    [ "$elapsed_ms" -lt 0 ] && elapsed_ms=0
    elapsed_sec=$(( elapsed_ms / 1000 ))
    elapsed_min=$(( elapsed_sec / 60 ))
    elapsed_rem=$(( elapsed_sec % 60 ))
    local centis=$(( (elapsed_ms % 1000) / 10 ))

    local rss_kb=0
    if [ -s "$tmp_time" ]; then
        local raw_rss
        raw_rss="$(awk '/maximum resident set size/{print $1}' "$tmp_time" | tail -n 1)"
        if [ -n "$raw_rss" ]; then
            if [ "$(uname -s)" = "Darwin" ]; then
                rss_kb=$(( raw_rss / 1024 ))
            else
                rss_kb="$raw_rss"
            fi
        fi
    fi
    rm -f "$tmp_time"

    local report
    report="$(printf '\tCommand being timed: "%s"\n\tUser time (seconds): 0.00\n\tSystem time (seconds): 0.00\n\tElapsed (wall clock) time (h:mm:ss or m:ss): %d:%02d.%02d\n\tMaximum resident set size (kbytes): %d\n\tExit status: %d\n' \
        "$*" "$elapsed_min" "$elapsed_rem" "$centis" "$rss_kb" "$rc")"

    if [ -n "$out_file" ]; then
        if [ "$append" -eq 1 ]; then
            printf '%s\n' "$report" >> "$out_file"
        else
            printf '%s\n' "$report" > "$out_file"
        fi
    else
        printf '%s\n' "$report" >&2
    fi
    return "$rc"
}

if [ -x /usr/bin/time ] && /usr/bin/time -v true >/dev/null 2>&1; then
    TIME_V_BIN="/usr/bin/time"
    TIME_V_ARGS=(-v)
    TIME_V="/usr/bin/time -v"
elif command -v gtime >/dev/null 2>&1 && gtime -v true >/dev/null 2>&1; then
    TIME_V_BIN="$(command -v gtime)"
    TIME_V_ARGS=(-v)
    TIME_V="$TIME_V_BIN -v"
else
    TIME_V_BIN="$PORTABLE_SH_PATH"
    TIME_V_ARGS=(--time-v)
    TIME_V="$PORTABLE_SH_PATH --time-v"
fi
TIME_V_CMD=("$TIME_V_BIN" "${TIME_V_ARGS[@]}")
export TIME_V TIME_V_BIN

mem_avail_gb() {
    if [ -r /proc/meminfo ]; then
        awk '/MemAvailable/{printf "%d\n", int($2 / 1048576)}' /proc/meminfo
    elif command -v vm_stat >/dev/null 2>&1 && command -v sysctl >/dev/null 2>&1; then
        local page_size free_pages
        page_size="$(sysctl -n hw.pagesize 2>/dev/null || echo 4096)"
        free_pages="$(vm_stat 2>/dev/null | awk '/Pages (free|inactive|speculative):/{gsub(/\./,"",$NF); sum+=$NF} END{print sum+0}')"
        if [ "${free_pages:-0}" -gt 0 ]; then
            printf '%d\n' "$(( (free_pages * page_size) / 1073741824 ))"
        else
            local mem_bytes
            mem_bytes="$(sysctl -n hw.memsize 2>/dev/null || echo 17179869184)"
            printf '%d\n' "$(( mem_bytes / 1073741824 ))"
        fi
    else
        echo 16
    fi
}

disk_free_gb() {
    local target="${1:-.}"
    df -Pk "$target" | awk 'NR==2 {printf "%d\n", int($4 / 1048576)}'
}

disk_total_gb() {
    local target="${1:-.}"
    df -Pk "$target" | awk 'NR==2 {printf "%d\n", int($2 / 1048576)}'
}

portable_loadavg() {
    if [ -r /proc/loadavg ]; then
        cut -d' ' -f1-3 /proc/loadavg
    elif command -v sysctl >/dev/null 2>&1 && sysctl -n vm.loadavg >/dev/null 2>&1; then
        sysctl -n vm.loadavg | tr -d '{}' | awk '{print $1, $2, $3}'
    else
        uptime | awk -F'load average[s]*: ' '{print $2}' | tr -d ',' | awk '{print $1, $2, $3}'
    fi
}

is_user_session_pid() {
    local pid="$1"
    if [ -r "/proc/$pid/cgroup" ]; then
        grep -q 'user@' "/proc/$pid/cgroup" 2>/dev/null
        return $?
    fi
    local p_uid
    p_uid="$(ps -o uid= -p "$pid" 2>/dev/null | tr -d '[:space:]')"
    [ -n "$p_uid" ] && [ "$p_uid" = "$(id -u)" ]
}

nice_idle() {
    if command -v ionice >/dev/null 2>&1; then
        nice -n 19 ionice -c 3 "$@"
    else
        nice -n 19 "$@"
    fi
}

portable_md5() {
    if [ $# -gt 0 ]; then
        if command -v md5sum >/dev/null 2>&1; then
            md5sum "$1" | awk '{print $1}'
        elif [ -x /sbin/md5 ]; then
            /sbin/md5 -q "$1"
        elif command -v openssl >/dev/null 2>&1; then
            openssl dgst -md5 "$1" | awk '{print $NF}'
        else
            python3 -c 'import hashlib,sys; h=hashlib.md5(); f=open(sys.argv[1],"rb"); [h.update(b) for b in iter(lambda: f.read(1<<20), b"")]; print(h.hexdigest())' "$1"
        fi
    else
        if command -v md5sum >/dev/null 2>&1; then
            md5sum | awk '{print $1}'
        elif [ -x /sbin/md5 ]; then
            /sbin/md5 -q
        elif command -v openssl >/dev/null 2>&1; then
            openssl dgst -md5 | awk '{print $NF}'
        else
            python3 -c 'import hashlib,sys; print(hashlib.md5(sys.stdin.buffer.read()).hexdigest())'
        fi
    fi
}

md5() {
    portable_md5 "$@"
}

portable_sha256() {
    if [ $# -gt 0 ]; then
        if command -v sha256sum >/dev/null 2>&1; then
            sha256sum "$1" | awk '{print $1}'
        elif command -v shasum >/dev/null 2>&1; then
            shasum -a 256 "$1" | awk '{print $1}'
        elif command -v openssl >/dev/null 2>&1; then
            openssl dgst -sha256 "$1" | awk '{print $NF}'
        else
            python3 -c 'import hashlib,sys; h=hashlib.sha256(); f=open(sys.argv[1],"rb"); [h.update(b) for b in iter(lambda: f.read(1<<20), b"")]; print(h.hexdigest())' "$1"
        fi
    else
        if command -v sha256sum >/dev/null 2>&1; then
            sha256sum | awk '{print $1}'
        elif command -v shasum >/dev/null 2>&1; then
            shasum -a 256 | awk '{print $1}'
        elif command -v openssl >/dev/null 2>&1; then
            openssl dgst -sha256 | awk '{print $NF}'
        else
            python3 -c 'import hashlib,sys; print(hashlib.sha256(sys.stdin.buffer.read()).hexdigest())'
        fi
    fi
}

sha256() {
    portable_sha256 "$@"
}

portable_filesize() {
    local f="$1"
    [ -f "$f" ] || return 1
    wc -c < "$f" | tr -d '[:space:]'
}

filesize() {
    portable_filesize "$@"
}

export -f ensure_g4_tools_path find_gitleaks now_ms portable_time_v mem_avail_gb \
    disk_free_gb disk_total_gb portable_loadavg is_user_session_pid nice_idle \
    portable_md5 md5 portable_sha256 sha256 portable_filesize filesize >/dev/null 2>&1 || true

if [ "${1:-}" = "--time-v" ]; then
    shift
    portable_time_v "$@"
    exit $?
fi
