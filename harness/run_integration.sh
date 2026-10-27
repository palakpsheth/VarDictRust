#!/usr/bin/env bash
set -euo pipefail

export LC_ALL=C
export LANG=C

sha256_hash() {
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum | awk '{print $1}'
    elif command -v shasum >/dev/null 2>&1; then
        shasum -a 256 | awk '{print $1}'
    else
        python3 -c 'import hashlib,sys; print(hashlib.sha256(sys.stdin.buffer.read()).hexdigest())'
    fi
}

# Run VarDictRust against VarDictJava integration test cases.
# Usage: run_integration.sh [test_case_file ...]
# If no files given, runs all Simple mode tests.
#
# Each test case file has:
#   Line 1: CSV header "Mode,fasta,bam,chr,start,end[,ampStart,ampEnd][,extra_args]"
#   Lines 2+: expected tab-delimited output

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
COMMON_GIT_DIR="$(git -C "$ROOT_DIR" rev-parse --path-format=absolute --git-common-dir 2>/dev/null || true)"
if [[ -n "$COMMON_GIT_DIR" && -d "$(dirname "$COMMON_GIT_DIR")/VarDictJava" ]]; then
    DATA_ROOT="$(cd "$(dirname "$COMMON_GIT_DIR")" && pwd)"
else
    DATA_ROOT="$ROOT_DIR"
fi
TESTDATA_DIR="VarDictJava/testdata/integrationtestcases"
FASTA_DIR="VarDictJava/testdata/fastas"
BAM_DIR="VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests"
if [[ ! -d "$TESTDATA_DIR" && -d "$DATA_ROOT/$TESTDATA_DIR" ]]; then
    TESTDATA_DIR="$DATA_ROOT/$TESTDATA_DIR"
    FASTA_DIR="$DATA_ROOT/$FASTA_DIR"
    BAM_DIR="$DATA_ROOT/$BAM_DIR"
fi
if [[ -z "${BINARY:-}" ]]; then
    if [[ -n "${VARDICT_BIN:-}" && -x "${VARDICT_BIN}" ]]; then
        BINARY="\"$VARDICT_BIN\" -th 1"
    elif [[ -x "$ROOT_DIR/target/release/vardictrust" ]]; then
        BINARY="\"$ROOT_DIR/target/release/vardictrust\" -th 1"
    else
        cargo build -q --release -p vardictrust
        BINARY="\"$ROOT_DIR/target/release/vardictrust\" -th 1"
    fi
fi

# Auto-materialize test FASTAs from upstream CSVs if not already present
if [[ ! -f "$FASTA_DIR/GRCh37.fa" || ! -f "$FASTA_DIR/GRCh37.fa.fai" ]]; then
    (cd "$DATA_ROOT" && python3 "$ROOT_DIR/harness/convert_fastas.py")
fi

# Load XFAIL defect registry
XFAIL_REGISTRY="$SCRIPT_DIR/integration_xfail.tsv"
declare -A XFAIL_CASES=()
declare -A XFAIL_ISSUES=()
declare -A XFAIL_SHAS=()

if [[ -f "$XFAIL_REGISTRY" ]]; then
    while IFS=$'\t' read -r fixture classification issue_id target_fix expected_sha256 rationale || [[ -n "$fixture" ]]; do
        fixture="${fixture//$'\r'/}"
        fixture="${fixture#"${fixture%%[![:space:]]*}"}"
        fixture="${fixture%"${fixture##*[![:space:]]}"}"
        classification="${classification//$'\r'/}"
        issue_id="${issue_id//$'\r'/}"
        expected_sha256="${expected_sha256//$'\r'/}"
        expected_sha256="${expected_sha256#"${expected_sha256%%[![:space:]]*}"}"
        expected_sha256="${expected_sha256%"${expected_sha256##*[![:space:]]}"}"
        [[ -z "$fixture" || "$fixture" =~ ^# ]] && continue
        if [[ "$classification" == "XFAIL" ]]; then
            if [[ ! "$expected_sha256" =~ ^[0-9a-f]{64}$ ]]; then
                echo "FATAL: $XFAIL_REGISTRY: entry '$fixture' missing valid 64-char hex expected_sha256 (got '$expected_sha256')" >&2
                exit 1
            fi
            XFAIL_SHAS["$fixture"]="$expected_sha256"
        fi
        XFAIL_CASES["$fixture"]="$classification"
        XFAIL_ISSUES["$fixture"]="$issue_id"
    done < "$XFAIL_REGISTRY"
fi

PASS=0
XFAIL=0
FAIL=0
XPASS=0
SKIP=0
ERRORS=()

# Advisory cpp third arm (Tasks T3A-04): availability probed once; per-fixture
# CPP lines never touch the PASS/XFAIL/FAIL gate or exit code (Task 3 ruling).
CPP_MISMATCH=0
CPP_AVAILABLE=0
CPP_SKIP_REASON=""
if [[ "${CPP_SKIP:-0}" == "1" || "${CPP_SKIP:-0}" == "true" ]]; then
    CPP_SKIP_REASON="disabled via CPP_SKIP=1"
elif ! command -v conda >/dev/null 2>&1; then
    CPP_SKIP_REASON="conda not on PATH"
elif [[ ! -x "$(dirname "$(dirname "$(command -v conda)")")/envs/vardict-oracle/bin/vardictcpp" ]] && ! conda run -n vardict-oracle vardictcpp --version >/dev/null 2>&1; then
    CPP_SKIP_REASON="vardictcpp not available in conda env vardict-oracle"
else
    CPP_AVAILABLE=1
fi

run_test() {
    local testfile="$1"
    local basename_file
    basename_file=$(basename "$testfile")

    # Parse first line for the early skip checks (re-parsed carefully below)
    IFS=',' read -r mode fasta bam chr _ < <(head -1 "$testfile")

    # Skip unsupported modes
    if [[ "$mode" != "Simple" && "$mode" != "Somatic" && "$mode" != "Amplicon" ]]; then
        echo "SKIP (mode=$mode): $basename_file"
        SKIP=$((SKIP+1))
        return
    fi

    local fasta_path="$FASTA_DIR/$fasta"
    if [[ ! -f "$fasta_path" ]]; then
        echo "SKIP (no fasta): $basename_file"
        SKIP=$((SKIP+1))
        return
    fi

    # Parse BAM (handle somatic pipe)
    local bam_arg
    if [[ "$bam" == *"|"* ]]; then
        local bam1 bam2
        IFS='|' read -r bam1 bam2 <<< "$bam"
        bam_arg="$BAM_DIR/$bam1|$BAM_DIR/$bam2"
    else
        bam_arg="$BAM_DIR/$bam"
    fi

    if [[ ! -f "${bam_arg%%|*}" ]]; then
        echo "SKIP (no bam): $basename_file"
        SKIP=$((SKIP+1))
        return
    fi

    # Extra args - handle the last field which may have more mixed-in args
    # The CSV format is: Mode,fasta,bam,chr,start,end[,ampStart,ampEnd][,extra]
    # But sometimes start-end is a single field
    local full_line
    full_line=$(head -1 "$testfile")
    
    # Re-parse more carefully: fields are comma-separated
    IFS=',' read -ra fields_raw <<< "$full_line"
    local fields=()
    for fld in "${fields_raw[@]}"; do
        # Trim leading and trailing whitespace
        fld="${fld#"${fld%%[![:space:]]*}"}"
        fld="${fld%"${fld##*[![:space:]]}"}"
        fields+=("$fld")
    done
    mode="${fields[0]}"
    fasta="${fields[1]}"
    bam="${fields[2]}"
    chr="${fields[3]}"
    
    local start_val end_val
    local extra_flags=""
    
    # Field 4 could be "start-end" or just "start"
    if [[ "${fields[4]}" == *"-"* ]]; then
        IFS='-' read -r start_val end_val <<< "${fields[4]}"
    else
        start_val="${fields[4]}"
        end_val="${fields[5]:-}"
    fi

    # Remaining fields could be amp or extra args
    local field_idx=5
    if [[ "${fields[4]}" != *"-"* ]]; then
        field_idx=6
    fi
    
    # Check if remaining fields look like amplicon or extra flags
    while [[ $field_idx -lt ${#fields[@]} ]]; do
        local fld="${fields[$field_idx]}"
        if [[ "$fld" == -* ]] || [[ "$fld" == "" ]]; then
            extra_flags="$extra_flags $fld"
        fi
        field_idx=$((field_idx+1))
    done

    # Create temp BED
    local tmpbed
    tmpbed=$(mktemp --suffix=.bed)
    if [[ "$mode" == "Amplicon" ]]; then
        local amp_start="${fields[6]}"
        local amp_end="${fields[7]}"
        echo -e "${chr}\t${start_val}\t${end_val}\ttestbed\t0\t0\t${amp_start}\t${amp_end}" > "$tmpbed"
    else
        echo -e "${chr}\t${start_val}\t${end_val}\ttestbed" > "$tmpbed"
    fi

    # Get expected output (everything after line 1)
    local expected
    expected=$(tail -n +2 "$testfile")

    # Build command
    local cmd="$BINARY -z -c 1 -S 2 -E 3 -g 4 -G $fasta_path -b '$bam_arg' --emit-intermediate $extra_flags $tmpbed"

    # Run with exit code and error capture
    local tmperr
    tmperr=$(mktemp)
    local exit_code=0
    local actual
    actual=$(eval $cmd 2>"$tmperr") || exit_code=$?

    # Engine crash / non-zero exit code must NEVER be masked as XFAIL
    if [[ $exit_code -ne 0 ]]; then
        echo "CRASH/ERROR (exit $exit_code): $basename_file"
        cat "$tmperr" >&2
        rm -f "$tmperr" "$tmpbed"
        FAIL=$((FAIL+1))
        ERRORS+=("CRASH: $basename_file")
        return
    fi
    rm -f "$tmperr"

    local is_xfail=0
    local issue=""
    local expected_sha=""
    if [[ -v XFAIL_CASES["$basename_file"] ]]; then
        is_xfail=1
        issue="${XFAIL_ISSUES[$basename_file]:-defect}"
        expected_sha="${XFAIL_SHAS[$basename_file]:-}"
    fi

    # Compare
    local sorted_actual sorted_expected
    sorted_actual=$(printf '%s\n' "$actual" | LC_ALL=C sort)
    sorted_expected=$(printf '%s\n' "$expected" | LC_ALL=C sort)
    if [[ "$sorted_actual" == "$sorted_expected" ]]; then
        if [[ $is_xfail -eq 1 ]]; then
            echo "XPASS: $basename_file (unexpected pass for registered $issue)"
            XPASS=$((XPASS+1))
            ERRORS+=("XPASS: $basename_file")
        else
            echo "MATCH: $basename_file"
            PASS=$((PASS+1))
        fi
    else
        if [[ $is_xfail -eq 1 ]]; then
            local actual_sha256
            actual_sha256=$(printf '%s\n' "$sorted_actual" | sha256_hash)
            if [[ -n "$expected_sha" && "$actual_sha256" == "$expected_sha" ]]; then
                echo "XFAIL: $basename_file ($issue upstream fixture defect, verified against live-oracle signature)"
                XFAIL=$((XFAIL+1))
            else
                echo "FAIL/REGRESSION: $basename_file ($issue actual sha256 $actual_sha256 != expected $expected_sha)"
                FAIL=$((FAIL+1))
                ERRORS+=("REGRESSION: $basename_file (actual sha256 $actual_sha256 != expected $expected_sha)")
                
                # Show diff (first 20 lines)
                diff <(echo "$sorted_expected") <(echo "$sorted_actual") | head -40 || true
                echo "---"
            fi
        else
            echo "MISMATCH: $basename_file"
            FAIL=$((FAIL+1))
            ERRORS+=("MISMATCH: $basename_file")
            
            # Show diff (first 20 lines)
            diff <(echo "$sorted_expected") <(echo "$sorted_actual") | head -40 || true
            echo "---"
        fi
    fi

    # Advisory cpp third arm: vardictcpp TSV with flags identical to the gate
    # (-z boolean, positional $tmpbed; -N pinned to the oracle sample label
    # since cpp requires -N and never derives it from the BAM name). This cell
    # intentionally does NOT reuse run_cpp.sh: Tier-1 needs the TSV surface
    # (run_cpp.sh is VCF-only) plus per-fixture -N pinning and --cols 3-61
    # normalization below. Tier-1 compare_var.py --cols 3-61: col 1 is our own
    # -N echoed back, col 2 is the known-normalized BED-name echo.
    # -i fixtures: cpp ignores -i (no splice-table surface, probed live) -> SKIP.
    if [[ $CPP_AVAILABLE -eq 1 && " $extra_flags " == *" -i "* ]]; then
        echo "CPP: SKIP: $basename_file (-i splice-table surface unsupported by vardictcpp)"
    elif [[ $CPP_AVAILABLE -eq 1 ]]; then
        if [[ "$start_val" =~ ^[0-9]+$ && "$end_val" =~ ^[0-9]+$ ]]; then
            local cpp_sample=${expected%%$'\n'*}
            cpp_sample=${cpp_sample%%$'\t'*}
            if [[ -z "$cpp_sample" ]]; then cpp_sample=$(basename "${bam_arg%%|*}" .bam); fi
            local -a extra_argv=()
            if [[ -n "$extra_flags" ]]; then
                read -r -a extra_argv <<< "$extra_flags"
            fi
            for i in "${!extra_argv[@]}"; do [[ "${extra_argv[$i]}" == "-nosv" ]] && extra_argv[$i]="--nosv"; done
            local -a gs_inputs=(--input "$fasta_path")
            [[ -f "$fasta_path.fai" ]] && gs_inputs+=(--input "$fasta_path.fai")
            if [[ "$bam_arg" == *"|"* ]]; then
                gs_inputs+=(--input "${bam_arg%%|*}" --input "${bam_arg#*|}")
            else
                gs_inputs+=(--input "$bam_arg")
            fi
            gs_inputs+=(--input "$tmpbed")
            local tmpcpp tmperr
            tmpcpp=$(mktemp); tmperr=$(mktemp)
            local cpp_exit=0
            bash "$SCRIPT_DIR/golden_store.sh" raw --pin-cpp-exit --kind cpp --label "int-$basename_file" --out "$tmpcpp" --err-out "$tmperr" "${gs_inputs[@]}" -- conda run -n vardict-oracle vardictcpp -z -c 1 -S 2 -E 3 -g 4 -G "$fasta_path" -b "$bam_arg" -N "$cpp_sample" "${extra_argv[@]}" "$tmpbed" || cpp_exit=$?
            if [[ $cpp_exit -ne 0 ]]; then
                local first_err
                first_err=$(grep -vE '^(ERROR conda\.cli|cpp_chunk=|[[:space:]]*$)' "$tmperr" 2>/dev/null | head -n 1 | tr -d '\r' | cut -c 1-140 || true)
                echo "CPP: SKIP: $basename_file (cpp engine exit $cpp_exit${first_err:+: $first_err})"
                rm -f "$tmpcpp" "$tmperr"
            else
                rm -f "$tmperr"
                local tmpexp tmpgot
                tmpexp=$(mktemp); tmpgot=$(mktemp)
                printf '%s\n' "$expected" | LC_ALL=C sort > "$tmpexp"
                LC_ALL=C sort "$tmpcpp" > "$tmpgot"
                rm -f "$tmpcpp"
                local cpp_msg cpp_rc=0
                cpp_msg=$(python3 "$SCRIPT_DIR/compare_var.py" --eps 1e-6 --cols 3-61 "$tmpexp" "$tmpgot" 2>&1) || cpp_rc=$?
                rm -f "$tmpexp" "$tmpgot"
                if [[ $cpp_rc -eq 0 ]]; then
                    echo "CPP: MATCH ($basename_file)"
                else
                    echo "CPP-MISMATCH: $basename_file ($(head -n 1 <<< "$cpp_msg"))"
                    CPP_MISMATCH=$((CPP_MISMATCH+1))
                fi
            fi
        else
            echo "CPP: SKIP: $basename_file (non-numeric region $chr:$start_val-$end_val)"
        fi
    fi
    rm -f "$tmpbed"
}

# If specific files passed, test them; otherwise test all cases in parallel
if [[ $# -gt 0 ]]; then
    TEST_FILES=("$@")
else
    TEST_FILES=("$TESTDATA_DIR"/*.txt)
fi

DEFAULT_JOBS=8
if [[ "${CI:-}" == "true" || "${GOLDEN_MODE:-auto}" != "auto" || ! -d "$ROOT_DIR/scratch/golden-store/cells" || -z "$(ls -A "$ROOT_DIR/scratch/golden-store/cells" 2>/dev/null)" ]]; then
    DEFAULT_JOBS=1
fi
JOBS="${INTEGRATION_JOBS:-$DEFAULT_JOBS}"
SUITE_START_MS=$(date +%s%3N)

if [[ "$JOBS" -le 1 || "${#TEST_FILES[@]}" -le 1 ]]; then
    for idx in "${!TEST_FILES[@]}"; do
        run_test "${TEST_FILES[$idx]}"
        done_n=$((idx + 1))
        if [[ -n "${INTEGRATION_PROGRESS_FD:-}" ]] && (( done_n % 25 == 0 || done_n == ${#TEST_FILES[@]} )); then
            printf '[%s] [integration] [%d/%d] fixtures complete (%dms)\n' "$(date '+%H:%M:%S')" "$done_n" "${#TEST_FILES[@]}" "$(( $(date +%s%3N) - SUITE_START_MS ))" >&"$INTEGRATION_PROGRESS_FD"
        fi
    done
else
    WORK_TMP=$(mktemp -d)
    trap 'rm -rf "$WORK_TMP"' EXIT
    (
        for idx in "${!TEST_FILES[@]}"; do
            while [[ $(jobs -rp | wc -l) -ge $JOBS ]]; do
                wait -n 2>/dev/null || sleep 0.02
            done
            {
                PASS=0; XFAIL=0; FAIL=0; SKIP=0; XPASS=0; CPP_MISMATCH=0; ERRORS=()
                ( set +e; run_test "${TEST_FILES[$idx]}" > "$WORK_TMP/$idx.out" 2>"$WORK_TMP/$idx.err"; \
                  printf '%d\t%d\t%d\t%d\t%d\t%d\n' "$PASS" "$XFAIL" "$FAIL" "$SKIP" "$XPASS" "$CPP_MISMATCH" > "$WORK_TMP/$idx.meta.tmp" && mv "$WORK_TMP/$idx.meta.tmp" "$WORK_TMP/$idx.meta"; \
                  if [[ ${#ERRORS[@]} -gt 0 ]]; then printf '%s\n' "${ERRORS[@]}" > "$WORK_TMP/$idx.errs"; fi )
            } &
        done
        wait
    ) &
    POOL_PID=$!
    for idx in "${!TEST_FILES[@]}"; do
        while [[ ! -f "$WORK_TMP/$idx.meta" ]] && kill -0 "$POOL_PID" 2>/dev/null; do
            sleep 0.02
        done
        [[ -s "$WORK_TMP/$idx.err" ]] && cat "$WORK_TMP/$idx.err" >&2
        [[ -f "$WORK_TMP/$idx.out" ]] && cat "$WORK_TMP/$idx.out"
        if [[ -f "$WORK_TMP/$idx.meta" ]]; then
            IFS=$'\t' read -r p xf fl sk xp cm < "$WORK_TMP/$idx.meta"
            PASS=$((PASS + p))
            XFAIL=$((XFAIL + xf))
            FAIL=$((FAIL + fl))
            SKIP=$((SKIP + sk))
            XPASS=$((XPASS + xp))
            CPP_MISMATCH=$((CPP_MISMATCH + cm))
        else
            FAIL=$((FAIL + 1))
            ERRORS+=("MISSING_META: ${TEST_FILES[$idx]}")
        fi
        if [[ -f "$WORK_TMP/$idx.errs" ]]; then
            while IFS= read -r e; do
                ERRORS+=("$e")
            done < "$WORK_TMP/$idx.errs"
        fi
        done_n=$((idx + 1))
        if [[ -n "${INTEGRATION_PROGRESS_FD:-}" ]] && (( done_n % 25 == 0 || done_n == ${#TEST_FILES[@]} )); then
            printf '[%s] [integration] [%d/%d] fixtures complete (%dms)\n' "$(date '+%H:%M:%S')" "$done_n" "${#TEST_FILES[@]}" "$(( $(date +%s%3N) - SUITE_START_MS ))" >&"$INTEGRATION_PROGRESS_FD"
        fi
    done
    wait "$POOL_PID" 2>/dev/null || true
fi

echo ""
echo "===== Results ====="
echo "PASS: $PASS  XFAIL: $XFAIL  FAIL: $FAIL  SKIP: $SKIP  XPASS: $XPASS"
if [[ $CPP_AVAILABLE -eq 1 ]]; then
    echo "CPP-MISMATCH: $CPP_MISMATCH"
else
    echo "CPP: SKIP: $CPP_SKIP_REASON"
fi

exit_code=0
if [[ $# -eq 0 ]]; then
    if [[ $SKIP -gt 0 ]]; then
        echo "ERROR: Full integration suite had $SKIP skipped tests!" >&2
        exit_code=1
    fi
    if [[ $PASS -ne 91 || $XFAIL -ne 1 || $FAIL -ne 0 || $XPASS -ne 0 ]]; then
        echo "ERROR: Corpus count invariant violated! Expected 91 PASS / 1 XFAIL / 0 FAIL / 0 SKIP / 0 XPASS, got $PASS PASS / $XFAIL XFAIL / $FAIL FAIL / $SKIP SKIP / $XPASS XPASS." >&2
        exit_code=1
    fi
fi

if [[ ${#ERRORS[@]} -gt 0 || $FAIL -gt 0 || $XPASS -gt 0 ]]; then
    echo "Failed tests:"
    for e in "${ERRORS[@]}"; do
        echo "  - $e"
    done
    exit_code=1
fi

exit "$exit_code"
