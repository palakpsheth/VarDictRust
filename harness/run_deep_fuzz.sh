#!/usr/bin/env bash
# harness/run_deep_fuzz.sh — Courteous, POSIX-portable multi-target cargo-fuzz runner.
#
# Pre-builds all 8 fuzz targets sequentially under nice_idle + CARGO_BUILD_JOBS=4
# to prevent concurrent ASan rustc lock/RSS spikes, then runs <=2 concurrent
# libFuzzer workers with isolated temp corpus/artifact directories.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck source=harness/lib/portable.sh
source "$SCRIPT_DIR/lib/portable.sh"

DRY_RUN=0
FAST_SECS=900
DEEP_SECS=1500
JOBS=2

usage() {
    cat <<'EOF'
Usage: bash harness/run_deep_fuzz.sh [--dry-run] [--fast-secs <N>] [--deep-secs <N>] [--jobs <1|2>]

Options:
  --dry-run         Validate arguments and target files, print schedule, and exit 0
  --fast-secs <N>   Per-target duration in seconds for fast in-memory targets (default: 900)
  --deep-secs <N>   Per-target duration in seconds for deep pipeline/I/O targets (default: 1500)
  --jobs <1|2>      Maximum concurrent libFuzzer workers after sequential build (default: 2)
  -h, --help        Show this help message and exit 0

Environment:
  FUZZ_OUT_DIR      Output directory for summary.json, logs/, and artifacts/ (default: scratch/fuzz-campaign)
  MIN_AVAIL_GB      Minimum available RAM in GiB required before live execution (default: 12)
  CARGO_BUILD_JOBS  Parallel rustc codegen jobs during sequential build (default: 4)
EOF
}

while [ $# -gt 0 ]; do
    case "$1" in
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        --fast-secs)
            [ $# -ge 2 ] || { echo "ERROR: --fast-secs requires an integer argument" >&2; exit 2; }
            FAST_SECS="$2"
            shift 2
            ;;
        --deep-secs)
            [ $# -ge 2 ] || { echo "ERROR: --deep-secs requires an integer argument" >&2; exit 2; }
            DEEP_SECS="$2"
            shift 2
            ;;
        --jobs)
            [ $# -ge 2 ] || { echo "ERROR: --jobs requires 1 or 2" >&2; exit 2; }
            JOBS="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "ERROR: unknown argument '$1'" >&2
            usage >&2
            exit 2
            ;;
    esac
done

if ! [[ "$FAST_SECS" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: --fast-secs must be a positive integer (got '$FAST_SECS')" >&2
    exit 2
fi
if ! [[ "$DEEP_SECS" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: --deep-secs must be a positive integer (got '$DEEP_SECS')" >&2
    exit 2
fi
if [ "$JOBS" != "1" ] && [ "$JOBS" != "2" ]; then
    echo "ERROR: --jobs must be 1 or 2 (got '$JOBS')" >&2
    exit 2
fi

FUZZ_CRATE_DIR="$ROOT_DIR/crates/vardictrust-core"
FUZZ_DIR="$FUZZ_CRATE_DIR/fuzz"
FUZZ_MANIFEST="$FUZZ_DIR/Cargo.toml"
FUZZ_TARGETS_DIR="$FUZZ_DIR/fuzz_targets"

if [ ! -f "$FUZZ_MANIFEST" ]; then
    echo "ERROR: missing fuzz manifest: $FUZZ_MANIFEST" >&2
    exit 1
fi

# Ordered schedule: 3 deep targets + 5 fast targets (paired so --jobs 2 runs in 4 balanced batches).
SCHEDULE_TARGETS=(
    "fuzz_cram_and_fai"
    "fuzz_parse_and_realign"
    "fuzz_sam_aux"
    "fuzz_somatic_fisher_and_amplicon"
    "fuzz_vcf_streamer"
    "fuzz_modify_cigar"
    "fuzz_bed_parser"
    "fuzz_realign"
)
SCHEDULE_TIERS=(
    "deep"
    "deep"
    "deep"
    "fast"
    "fast"
    "fast"
    "fast"
    "fast"
)

for tgt in "${SCHEDULE_TARGETS[@]}"; do
    tgt_file="$FUZZ_TARGETS_DIR/${tgt}.rs"
    if [ ! -f "$tgt_file" ]; then
        echo "ERROR: missing fuzz target source file: $tgt_file" >&2
        exit 1
    fi
done

FUZZ_OUT_RAW="${FUZZ_OUT_DIR:-scratch/fuzz-campaign}"
case "$FUZZ_OUT_RAW" in
    /*) OUT_DIR="$FUZZ_OUT_RAW" ;;
    *)  OUT_DIR="$ROOT_DIR/$FUZZ_OUT_RAW" ;;
esac

total_cpu_secs=$(( 3 * DEEP_SECS + 5 * FAST_SECS ))
echo "[deep-fuzz] Schedule: 8 targets | jobs=$JOBS | fast=${FAST_SECS}s | deep=${DEEP_SECS}s | total_cpu=${total_cpu_secs}s | out=$OUT_DIR"

idx=0
batch_num=1
total_targets="${#SCHEDULE_TARGETS[@]}"
while [ "$idx" -lt "$total_targets" ]; do
    b_end=$(( idx + JOBS ))
    [ "$b_end" -gt "$total_targets" ] && b_end="$total_targets"
    batch_desc=""
    k="$idx"
    while [ "$k" -lt "$b_end" ]; do
        t_name="${SCHEDULE_TARGETS[$k]}"
        t_tier="${SCHEDULE_TIERS[$k]}"
        if [ "$t_tier" = "deep" ]; then
            t_secs="$DEEP_SECS"
        else
            t_secs="$FAST_SECS"
        fi
        if [ -n "$batch_desc" ]; then
            batch_desc="${batch_desc}, "
        fi
        batch_desc="${batch_desc}${t_name} (${t_tier}, ${t_secs}s)"
        k=$(( k + 1 ))
    done
    echo "  Batch ${batch_num}: ${batch_desc}"
    batch_num=$(( batch_num + 1 ))
    idx="$b_end"
done

# Short-circuit --dry-run before probing cargo +nightly fuzz, checking mem_avail_gb, or touching scratch/.
if [ "$DRY_RUN" -eq 1 ]; then
    echo "[deep-fuzz] Dry-run validation succeeded (0 files written)."
    exit 0
fi

MIN_AVAIL="${MIN_AVAIL_GB:-12}"
AVAIL_GB="$(mem_avail_gb)"
if [ "${AVAIL_GB:-0}" -lt "$MIN_AVAIL" ]; then
    echo "ERROR: insufficient available memory (${AVAIL_GB} GiB < MIN_AVAIL_GB=${MIN_AVAIL} GiB)" >&2
    exit 1
fi

if ! (cd "$FUZZ_CRATE_DIR" && cargo +nightly fuzz --version >/dev/null 2>&1); then
    echo "ERROR: 'cargo +nightly fuzz' is required for live deep fuzzing (install nightly toolchain + cargo-fuzz)" >&2
    exit 1
fi

BUILD_JOBS="${CARGO_BUILD_JOBS:-4}"
ART_OUT_DIR="$OUT_DIR/artifacts"
LOG_OUT_DIR="$OUT_DIR/logs"
SUMMARY_JSON="$OUT_DIR/summary.json"

HAD_IN_TREE_CORPUS=0
HAD_IN_TREE_ARTIFACTS=0
[ -d "$FUZZ_DIR/corpus" ] && HAD_IN_TREE_CORPUS=1
[ -d "$FUZZ_DIR/artifacts" ] && HAD_IN_TREE_ARTIFACTS=1

mkdir -p "$ART_OUT_DIR" "$LOG_OUT_DIR"
TMP_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/vdr-deep-fuzz.XXXXXX")"
mkdir -p "$TMP_ROOT/corpus" "$TMP_ROOT/artifacts" "$TMP_ROOT/meta"

ACTIVE_PIDS=()
CLEANED_UP=0

copy_staged_artifacts() {
    local src_root tgt_dir art_f base_f tgt_base
    for src_root in "${TMP_ROOT:+$TMP_ROOT/artifacts}" "$FUZZ_DIR/artifacts"; do
        [ -n "$src_root" ] && [ -d "$src_root" ] || continue
        mkdir -p "$ART_OUT_DIR" 2>/dev/null || true
        for tgt_dir in "$src_root"/*; do
            [ -d "$tgt_dir" ] || continue
            tgt_base="$(basename "$tgt_dir")"
            for art_f in "$tgt_dir"/*; do
                [ -f "$art_f" ] || continue
                base_f="$(basename "$art_f")"
                cp -f "$art_f" "$ART_OUT_DIR/${tgt_base}__${base_f}" 2>/dev/null || true
            done
        done
    done
}

cleanup_campaign() {
    local exit_rc=$?
    if [ "$CLEANED_UP" -eq 1 ]; then
        return "$exit_rc"
    fi
    CLEANED_UP=1
    trap - EXIT INT TERM

    local pid
    for pid in ${ACTIVE_PIDS[@]+"${ACTIVE_PIDS[@]}"}; do
        if kill -0 "$pid" 2>/dev/null; then
            python3 -c '
import os, signal, subprocess, sys
root_pid = int(sys.argv[1])
out = subprocess.run(["ps", "-eo", "pid=,ppid="], capture_output=True, text=True, check=False).stdout
children = {}
for line in out.splitlines():
    parts = line.split()
    if len(parts) == 2 and parts[0].isdigit() and parts[1].isdigit():
        children.setdefault(int(parts[1]), []).append(int(parts[0]))
stack = [root_pid]
order = []
while stack:
    cur = stack.pop()
    order.append(cur)
    stack.extend(children.get(cur, []))
for p in reversed(order):
    try:
        os.kill(p, signal.SIGTERM)
    except OSError:
        pass
' "$pid" 2>/dev/null || kill "$pid" 2>/dev/null || true
        fi
    done
    for pid in ${ACTIVE_PIDS[@]+"${ACTIVE_PIDS[@]}"}; do
        wait "$pid" 2>/dev/null || true
    done

    copy_staged_artifacts

    if [ "$HAD_IN_TREE_CORPUS" -eq 0 ] && [ -d "$FUZZ_DIR/corpus" ]; then
        rm -rf "$FUZZ_DIR/corpus"
    fi
    if [ "$HAD_IN_TREE_ARTIFACTS" -eq 0 ] && [ -d "$FUZZ_DIR/artifacts" ]; then
        rm -rf "$FUZZ_DIR/artifacts"
    fi

    if [ -n "${TMP_ROOT:-}" ] && [ -d "$TMP_ROOT" ]; then
        rm -rf "$TMP_ROOT"
    fi

    [ -d /dev/shm ] && find /dev/shm -maxdepth 1 -user "$(id -un)" -name 'vdr-fuzz-*' -exec rm -rf {} + 2>/dev/null || true
    return "$exit_rc"
}

on_signal() {
    cleanup_campaign
    exit 130
}

trap cleanup_campaign EXIT
trap on_signal INT TERM

echo "[deep-fuzz] Pre-building all 8 fuzz binaries sequentially (CARGO_BUILD_JOBS=$BUILD_JOBS, avail_ram=${AVAIL_GB}GiB)..."
(
    cd "$FUZZ_CRATE_DIR"
    CARGO_BUILD_JOBS="$BUILD_JOBS" nice_idle cargo +nightly fuzz build
)

CAMPAIGN_T0_MS="$(now_ms)"
OVERALL_FAIL=0

idx=0
batch_num=1
while [ "$idx" -lt "$total_targets" ]; do
    b_end=$(( idx + JOBS ))
    [ "$b_end" -gt "$total_targets" ] && b_end="$total_targets"
    echo "[deep-fuzz] Starting batch ${batch_num}..."

    ACTIVE_PIDS=()
    BATCH_TARGETS=()
    k="$idx"
    while [ "$k" -lt "$b_end" ]; do
        target="${SCHEDULE_TARGETS[$k]}"
        tier="${SCHEDULE_TIERS[$k]}"
        if [ "$tier" = "deep" ]; then
            secs="$DEEP_SECS"
        else
            secs="$FAST_SECS"
        fi
        corp_dir="$TMP_ROOT/corpus/$target"
        art_dir="$TMP_ROOT/artifacts/$target"
        log_file="$LOG_OUT_DIR/${target}.log"
        meta_file="$TMP_ROOT/meta/${target}.tsv"
        mkdir -p "$corp_dir" "$art_dir"
        rm -f "$log_file"

        (
            child_pid=""
            trap '[ -n "$child_pid" ] && kill "$child_pid" 2>/dev/null || true; exit 130' INT TERM
            t0_ms="$(now_ms)"
            (
                cd "$FUZZ_CRATE_DIR"
                CARGO_BUILD_JOBS="$BUILD_JOBS" nice_idle cargo +nightly fuzz run "$target" "$corp_dir" -- \
                    -max_total_time="$secs" \
                    -timeout=10 \
                    -rss_limit_mb=2048 \
                    -artifact_prefix="$art_dir/" \
                    >"$log_file" 2>&1
            ) &
            child_pid=$!
            rc=0
            wait "$child_pid" || rc=$?
            t1_ms="$(now_ms)"
            elapsed_ms=$(( t1_ms - t0_ms ))
            [ "$elapsed_ms" -lt 0 ] && elapsed_ms=0
            printf '%s\t%s\t%s\t%s\t%s\n' "$target" "$tier" "$secs" "$rc" "$elapsed_ms" > "$meta_file"
            exit "$rc"
        ) &
        wpid=$!
        ACTIVE_PIDS+=("$wpid")
        BATCH_TARGETS+=("$target")
        k=$(( k + 1 ))
    done

    w_idx=0
    while [ "$w_idx" -lt "${#ACTIVE_PIDS[@]}" ]; do
        wpid="${ACTIVE_PIDS[$w_idx]}"
        wtarget="${BATCH_TARGETS[$w_idx]}"
        wrc=0
        wait "$wpid" || wrc=$?
        if [ "$wrc" -ne 0 ]; then
            echo "[deep-fuzz] FAIL: $wtarget exited with status $wrc (see $LOG_OUT_DIR/${wtarget}.log)" >&2
            OVERALL_FAIL=1
        else
            echo "[deep-fuzz] PASS: $wtarget completed cleanly"
        fi
        w_idx=$(( w_idx + 1 ))
    done
    ACTIVE_PIDS=()

    # Reclaim any stray /dev/shm/vdr-fuzz-* directories between batches as well.
    [ -d /dev/shm ] && find /dev/shm -maxdepth 1 -user "$(id -un)" -name 'vdr-fuzz-*' -exec rm -rf {} + 2>/dev/null || true

    batch_num=$(( batch_num + 1 ))
    idx="$b_end"
done

CAMPAIGN_T1_MS="$(now_ms)"
WALL_ELAPSED_MS=$(( CAMPAIGN_T1_MS - CAMPAIGN_T0_MS ))
[ "$WALL_ELAPSED_MS" -lt 0 ] && WALL_ELAPSED_MS=0

copy_staged_artifacts

GIT_HEAD="$(git -C "$ROOT_DIR" rev-parse HEAD 2>/dev/null || echo unknown)"

python3 - "$TMP_ROOT" "$LOG_OUT_DIR" "$ART_OUT_DIR" "$SUMMARY_JSON" "$GIT_HEAD" "$FAST_SECS" "$DEEP_SECS" "$JOBS" "$WALL_ELAPSED_MS" "${SCHEDULE_TARGETS[@]}" <<'PY'
import json
import pathlib
import re
import sys

tmp_root = pathlib.Path(sys.argv[1])
log_dir = pathlib.Path(sys.argv[2])
art_out_dir = pathlib.Path(sys.argv[3])
summary_path = pathlib.Path(sys.argv[4])
git_head = sys.argv[5]
fast_secs = int(sys.argv[6])
deep_secs = int(sys.argv[7])
jobs = int(sys.argv[8])
wall_elapsed_ms = int(sys.argv[9])
targets_order = sys.argv[10:]

done_re = re.compile(r"Done\s+(\d+)\s+runs\s+in\s+(\d+)\s+second")
step_re = re.compile(r"^#(\d+)\s+")
cov_re = re.compile(r"\bcov:\s*(\d+)")
ft_re = re.compile(r"\bft:\s*(\d+)")
exec_re = re.compile(r"\bexec/s:\s*(\d+)")
rss_re = re.compile(r"\brss:\s*(\d+)Mb")

target_entries = []
total_execs = 0
passed = 0
failed = 0
all_artifacts = []

for tgt in targets_order:
    meta_path = tmp_root / "meta" / f"{tgt}.tsv"
    tier = "fast"
    budget_secs = fast_secs
    exit_code = 1
    elapsed_ms = 0
    if meta_path.is_file():
        parts = meta_path.read_text(encoding="utf-8").strip().split("\t")
        if len(parts) >= 5:
            tier = parts[1]
            budget_secs = int(parts[2])
            exit_code = int(parts[3])
            elapsed_ms = int(parts[4])

    log_path = log_dir / f"{tgt}.log"
    runs = 0
    cov = 0
    ft = 0
    exec_per_sec = 0
    rss_mb = 0
    if log_path.is_file():
        log_txt = log_path.read_text(encoding="utf-8", errors="replace")
        for line in log_txt.splitlines():
            sm = step_re.match(line)
            if not sm:
                continue
            runs = max(runs, int(sm.group(1)))
            cm = cov_re.search(line)
            if cm:
                cov = int(cm.group(1))
            fm = ft_re.search(line)
            if fm:
                ft = int(fm.group(1))
            em = exec_re.search(line)
            if em:
                exec_per_sec = int(em.group(1))
            rm = rss_re.search(line)
            if rm:
                rss_mb = max(rss_mb, int(rm.group(1)))
        dm = done_re.search(log_txt)
        if dm:
            runs = max(runs, int(dm.group(1)))

    tgt_art_dir = tmp_root / "artifacts" / tgt
    tgt_arts = []
    if tgt_art_dir.is_dir():
        for p in sorted(tgt_art_dir.iterdir()):
            if p.is_file():
                saved_name = f"{tgt}__{p.name}"
                tgt_arts.append(saved_name)
                all_artifacts.append(saved_name)

    ok = exit_code == 0 and len(tgt_arts) == 0
    if ok:
        passed += 1
    else:
        failed += 1
    total_execs += runs

    target_entries.append(
        {
            "target": tgt,
            "tier": tier,
            "budget_secs": budget_secs,
            "elapsed_ms": elapsed_ms,
            "exit_code": exit_code,
            "runs": runs,
            "cov": cov,
            "ft": ft,
            "exec_per_sec": exec_per_sec,
            "peak_rss_mb": rss_mb,
            "artifacts": tgt_arts,
            "status": "PASS" if ok else "FAIL",
        }
    )

summary = {
    "status": "PASS" if failed == 0 else "FAIL",
    "git_head": git_head,
    "fast_secs": fast_secs,
    "deep_secs": deep_secs,
    "jobs": jobs,
    "timeout_secs": 10,
    "rss_limit_mb": 2048,
    "total_targets": len(targets_order),
    "passed_targets": passed,
    "failed_targets": failed,
    "total_executions": total_execs,
    "total_artifacts": len(all_artifacts),
    "artifacts": all_artifacts,
    "wall_elapsed_ms": wall_elapsed_ms,
    "targets": target_entries,
}

summary_path.parent.mkdir(parents=True, exist_ok=True)
summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
print(
    f"[deep-fuzz] Summary written to {summary_path}: status={summary['status']} "
    f"targets={passed}/{len(targets_order)} total_executions={total_execs} artifacts={len(all_artifacts)}"
)
if failed > 0:
    sys.exit(1)
PY

if [ "$OVERALL_FAIL" -ne 0 ]; then
    exit 1
fi
