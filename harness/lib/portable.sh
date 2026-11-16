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
    fi
    local conda_base=""
    if [ -n "${CONDA_PREFIX:-}" ]; then
        conda_base="${CONDA_PREFIX%/envs/*}"
    fi
    local cand
    for cand in \
        "${G4_TOOLS_BIN:+$G4_TOOLS_BIN/gitleaks}" \
        "${CONDA_PREFIX:+$CONDA_PREFIX/bin/gitleaks}" \
        "${conda_base:+$conda_base/envs/g4-tools/bin/gitleaks}" \
        "${conda_base:+$conda_base/envs/vardict-oracle/bin/gitleaks}" \
        "${HOME:-}/miniconda3/envs/g4-tools/bin/gitleaks" \
        "${HOME:-}/miniconda3/envs/vardict-oracle/bin/gitleaks" \
        "${HOME:-}/miniconda3/bin/gitleaks" \
        "${HOME:-}/miniforge3/envs/g4-tools/bin/gitleaks" \
        "${HOME:-}/miniforge3/envs/vardict-oracle/bin/gitleaks" \
        "${HOME:-}/miniforge3/bin/gitleaks" \
        "${HOME:-}/.local/bin/gitleaks" \
        "/opt/conda/envs/g4-tools/bin/gitleaks" \
        "/opt/conda/envs/vardict-oracle/bin/gitleaks" \
        "/opt/conda/bin/gitleaks" \
        "/usr/local/bin/gitleaks" \
        "/opt/homebrew/bin/gitleaks"
    do
        if [ -n "$cand" ] && [ -x "$cand" ]; then
            printf '%s\n' "$cand"
            return 0
        fi
    done
    return 1
}

resolve_data_root() {
    if [ -n "${DATA_ROOT:-}" ]; then
        printf '%s\n' "$DATA_ROOT"
        return 0
    fi
    local root_dir="${1:-${ROOT_DIR:-$(pwd)}}"
    root_dir="$(cd "$root_dir" 2>/dev/null && pwd || printf '%s' "$root_dir")"
    local common_git_dir=""
    common_git_dir="$(git -C "$root_dir" rev-parse --path-format=absolute --git-common-dir 2>/dev/null || true)"
    if [ -z "$common_git_dir" ]; then
        local rel_common
        rel_common="$(git -C "$root_dir" rev-parse --git-common-dir 2>/dev/null || true)"
        if [ -n "$rel_common" ]; then
            common_git_dir="$(cd "$root_dir" && cd "$rel_common" 2>/dev/null && pwd || true)"
        fi
    fi
    if [ -n "$common_git_dir" ]; then
        local main_root
        main_root="$(cd "$(dirname "$common_git_dir")" 2>/dev/null && pwd || true)"
        if [ -n "$main_root" ] && [ "$main_root" != "$root_dir" ] \
           && { [ -d "$main_root/scratch" ] || [ -d "$main_root/VarDictJava" ]; }; then
            if [ -d "$root_dir/scratch/g4/ref" ] && [ ! -L "$root_dir/scratch/g4/ref" ] \
               && [ ! -L "$root_dir/scratch/g4" ] && [ ! -L "$root_dir/scratch" ] \
               && { [ -d "$root_dir/VarDictJava/testdata" ] || [ ! -d "$main_root/VarDictJava/testdata" ]; }; then
                printf '%s\n' "$root_dir"
                return 0
            fi
            printf '%s\n' "$main_root"
            return 0
        fi
    fi
    printf '%s\n' "$root_dir"
}

portable_timeout() {
    if [ $# -lt 2 ]; then
        echo "portable_timeout: usage: portable_timeout <duration> <cmd> [args...]" >&2
        return 125
    fi
    local dur="$1"
    shift
    local backend="${PORTABLE_TIMEOUT_BACKEND:-auto}"
    if [ "$backend" = "auto" ] || [ "$backend" = "gnu" ]; then
        if command -v timeout >/dev/null 2>&1 && timeout 1s true >/dev/null 2>&1; then
            timeout "$dur" "$@"
            return $?
        elif command -v gtimeout >/dev/null 2>&1 && gtimeout 1s true >/dev/null 2>&1; then
            gtimeout "$dur" "$@"
            return $?
        elif [ "$backend" = "gnu" ]; then
            echo "portable_timeout: GNU timeout backend requested but unavailable" >&2
            return 127
        fi
    fi
    python3 -c '
import os, re, signal, subprocess, sys, time
raw = sys.argv[1]
cmd = sys.argv[2:]
m = re.fullmatch(r"([0-9]+(?:\.[0-9]+)?|\.[0-9]+)([smhd]?)", raw.strip())
if not m:
    sys.stderr.write(f"portable_timeout: invalid duration {raw!r}\n")
    sys.exit(125)
secs = float(m.group(1)) * {"": 1.0, "s": 1.0, "m": 60.0, "h": 3600.0, "d": 86400.0}[m.group(2)]
try:
    proc = subprocess.Popen(cmd, start_new_session=True)
except FileNotFoundError as exc:
    sys.stderr.write(f"portable_timeout: {cmd[0]}: {exc.strerror}\n")
    sys.exit(127)
except PermissionError as exc:
    sys.stderr.write(f"portable_timeout: {cmd[0]}: {exc.strerror}\n")
    sys.exit(126)
def _fwd(signum, _frame):
    try:
        os.killpg(proc.pid, signum)
    except ProcessLookupError:
        pass
for _sig in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
    signal.signal(_sig, _fwd)
try:
    rc = proc.wait(timeout=secs if secs > 0 else None)
except subprocess.TimeoutExpired:
    try:
        os.killpg(proc.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    deadline = time.monotonic() + 2.0
    while proc.poll() is None and time.monotonic() < deadline:
        time.sleep(0.02)
    try:
        os.killpg(proc.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    proc.wait()
    sys.exit(124)
sys.exit(128 + (-rc) if rc < 0 else rc)
' "$dur" "$@"
}

TIMEOUT_CMD=("$PORTABLE_SH_PATH" --timeout)
TIMEOUT="$PORTABLE_SH_PATH --timeout"
export TIMEOUT

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
    local backend="${PORTABLE_TIME_BACKEND:-auto}"
    local time_bin="${PORTABLE_TIME_BIN:-/usr/bin/time}"
    if [ "$backend" = "auto" ] || [ "$backend" = "gnu" ]; then
        if [ -x "$time_bin" ] && "$time_bin" -v true >/dev/null 2>&1; then
            "$time_bin" -v "$@"
            return $?
        elif command -v gtime >/dev/null 2>&1 && gtime -v true >/dev/null 2>&1; then
            gtime -v "$@"
            return $?
        elif [ "$backend" = "gnu" ]; then
            echo "portable_time_v: GNU time backend requested but unavailable" >&2
            return 127
        fi
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
    if { [ "$backend" = "auto" ] || [ "$backend" = "bsd" ]; } && [ -x "$time_bin" ] && "$time_bin" -l true >/dev/null 2>&1; then
        { "$time_bin" -l sh -c 'exec "$@" 2>&3 3>&-' sh "$@" 2>"$tmp_time"; } 3>&2 || rc=$?
    elif { [ "$backend" = "auto" ] || [ "$backend" = "python" ]; } && command -v python3 >/dev/null 2>&1; then
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
            if [ "${PORTABLE_UNAME_S:-$(uname -s)}" = "Darwin" ]; then
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
        grep -Eq "user@$(id -u)\\.service" "/proc/$pid/cgroup" 2>/dev/null
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

exec_nice_idle() {
    if command -v ionice >/dev/null 2>&1; then
        exec nice -n 19 ionice -c 3 "$@"
    else
        exec nice -n 19 "$@"
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

stage_verified_sha256() {
    local url="${1:-}"
    local expected_sha256="${2:-}"
    local dest="${3:-}"
    if [ -z "$url" ] || [ -z "$expected_sha256" ] || [ -z "$dest" ]; then
        echo "stage_verified_sha256: usage: stage_verified_sha256 <url> <sha256> <dest>" >&2
        return 1
    fi
    if [[ ! "$expected_sha256" =~ ^[0-9a-f]{64}$ ]]; then
        echo "stage_verified_sha256: invalid expected SHA-256 digest '$expected_sha256' for $dest" >&2
        return 1
    fi
    if [ -f "$dest" ]; then
        local existing_sha256
        existing_sha256="$(portable_sha256 "$dest")"
        if [ "$existing_sha256" = "$expected_sha256" ]; then
            return 0
        fi
    fi
    mkdir -p "$(dirname "$dest")"
    local part="${dest}.part"
    rm -f "$part"
    if ! curl -sSLf --retry 5 --retry-delay 5 -o "$part" "$url"; then
        rm -f "$part"
        echo "stage_verified_sha256: download failed for $dest ($url)" >&2
        return 1
    fi
    local got_sha256
    got_sha256="$(portable_sha256 "$part")"
    if [ "$got_sha256" != "$expected_sha256" ]; then
        rm -f "$part"
        echo "stage_verified_sha256: SHA-256 MISMATCH for $dest (got=$got_sha256 want=$expected_sha256)" >&2
        return 1
    fi
    mv -f "$part" "$dest"
}

stream_verified_bam_slice() {
    local bed_file=""
    local threads=""
    while [ $# -gt 0 ]; do
        case "$1" in
            --bed)
                bed_file="${2:-}"
                shift 2
                ;;
            --threads)
                threads="${2:-}"
                shift 2
                ;;
            --)
                shift
                break
                ;;
            *)
                break
                ;;
        esac
    done
    local label="${1:-}"
    local bam_url="${2:-}"
    local bai_url="${3:-}"
    local bai_sha256="${4:-}"
    local local_bai="${5:-}"
    local region="${6:-}"
    local out_bam="${7:-}"
    if [ -z "$label" ] || [ -z "$bam_url" ] || [ -z "$bai_url" ] || [ -z "$bai_sha256" ] || [ -z "$local_bai" ] || [ -z "$out_bam" ]; then
        echo "stream_verified_bam_slice: usage: stream_verified_bam_slice [--bed <bed>] [--threads <N>] <label> <bam_url> <bai_url> <bai_sha256> <local_bai> <region> <out_bam>" >&2
        return 1
    fi
    if [ -z "$bed_file" ] && [ -z "$region" ]; then
        echo "stream_verified_bam_slice: either --bed <bed> or <region> is required for $label" >&2
        return 1
    fi
    if [ -s "$out_bam" ] && [ -s "${out_bam}.bai" ] && samtools quickcheck "$out_bam" 2>/dev/null; then
        echo "  [cached] $label ($out_bam)"
        return 0
    fi
    mkdir -p "$(dirname "$out_bam")"
    stage_verified_sha256 "$bai_url" "$bai_sha256" "$local_bai"

    local max_retries="${BAM_SLICE_MAX_RETRIES:-5}"
    local delay="${BAM_SLICE_INITIAL_DELAY:-2}"
    local attempt=1
    local tmp_bam="${out_bam}.tmp"
    local tmp_unsorted="${out_bam}.unsorted.tmp"
    local th_args=()
    if [ -n "$threads" ]; then
        th_args=(-@ "$threads")
    fi

    while [ "$attempt" -le "$max_retries" ]; do
        echo "Streaming $label from $bam_url (attempt $attempt/$max_retries)..."
        rm -f "$tmp_bam" "$tmp_unsorted"
        local stream_ok=0
        if [ -n "$bed_file" ]; then
            if samtools view ${th_args[@]+"${th_args[@]}"} -b -M -L "$bed_file" -X "$bam_url" "$local_bai" > "$tmp_unsorted" \
               && samtools sort ${th_args[@]+"${th_args[@]}"} -o "$tmp_bam" "$tmp_unsorted" \
               && samtools quickcheck "$tmp_bam"; then
                stream_ok=1
            fi
        else
            if samtools view ${th_args[@]+"${th_args[@]}"} -b -X "$bam_url" "$local_bai" "$region" > "$tmp_bam" \
               && samtools quickcheck "$tmp_bam"; then
                stream_ok=1
            fi
        fi
        rm -f "$tmp_unsorted"
        if [ "$stream_ok" -eq 1 ]; then
            mv -f "$tmp_bam" "$out_bam"
            samtools index ${th_args[@]+"${th_args[@]}"} "$out_bam"
            samtools quickcheck "$out_bam"
            return 0
        fi
        rm -f "$tmp_bam"
        attempt=$((attempt + 1))
        if [ "$attempt" -le "$max_retries" ]; then
            echo "Sleeping $delay seconds before retry..." >&2
            sleep "$delay"
            delay=$((delay * 2))
        fi
    done
    echo "Failed to stream $label BAM slice after $max_retries attempts" >&2
    return 1
}

export -f ensure_g4_tools_path find_gitleaks resolve_data_root portable_timeout now_ms portable_time_v mem_avail_gb \
    disk_free_gb disk_total_gb portable_loadavg is_user_session_pid nice_idle exec_nice_idle \
    portable_md5 md5 portable_sha256 sha256 portable_filesize filesize stage_verified_sha256 stream_verified_bam_slice >/dev/null 2>&1 || true

if [ "${1:-}" = "--time-v" ]; then
    shift
    portable_time_v "$@"
    exit $?
elif [ "${1:-}" = "--timeout" ]; then
    shift
    portable_timeout "$@"
    exit $?
fi
