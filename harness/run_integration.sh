#!/usr/bin/env bash
set -euo pipefail

# Run VarDictRust against VarDictJava integration test cases.
# Usage: run_integration.sh [test_case_file ...]
# If no files given, runs all Simple mode tests.
#
# Each test case file has:
#   Line 1: CSV header "Mode,fasta,bam,chr,start,end[,ampStart,ampEnd][,extra_args]"
#   Lines 2+: expected tab-delimited output

TESTDATA_DIR="VarDictJava/testdata/integrationtestcases"
FASTA_DIR="VarDictJava/testdata/fastas"
BAM_DIR="VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests"
BINARY="cargo run -q --release -p vardictrust --"

# Auto-materialize test FASTAs from upstream CSVs if not already present
if [[ ! -f "$FASTA_DIR/GRCh37.fa" || ! -f "$FASTA_DIR/GRCh37.fa.fai" ]]; then
    python3 harness/convert_fastas.py
fi

# Load XFAIL defect registry
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
XFAIL_REGISTRY="$SCRIPT_DIR/integration_xfail.tsv"
declare -A XFAIL_CASES=()
declare -A XFAIL_ISSUES=()

if [[ -f "$XFAIL_REGISTRY" ]]; then
    while IFS=$'\t' read -r fixture classification issue_id target_fix rationale || [[ -n "$fixture" ]]; do
        fixture="${fixture//$'\r'/}"
        fixture="${fixture#"${fixture%%[![:space:]]*}"}"
        fixture="${fixture%"${fixture##*[![:space:]]}"}"
        classification="${classification//$'\r'/}"
        issue_id="${issue_id//$'\r'/}"
        [[ -z "$fixture" || "$fixture" =~ ^# ]] && continue
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

    rm -f "$tmpbed"

    # Engine crash / non-zero exit code must NEVER be masked as XFAIL
    if [[ $exit_code -ne 0 ]]; then
        echo "CRASH/ERROR (exit $exit_code): $basename_file"
        cat "$tmperr" >&2
        rm -f "$tmperr"
        FAIL=$((FAIL+1))
        ERRORS+=("CRASH: $basename_file")
        return
    fi
    rm -f "$tmperr"

    local is_xfail=0
    local issue=""
    if [[ -v XFAIL_CASES["$basename_file"] ]]; then
        is_xfail=1
        issue="${XFAIL_ISSUES[$basename_file]:-defect}"
    fi

    # Compare
    local sorted_actual sorted_expected
    sorted_actual=$(echo "$actual" | sort)
    sorted_expected=$(echo "$expected" | sort)
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
            echo "XFAIL: $basename_file ($issue upstream fixture defect)"
            XFAIL=$((XFAIL+1))
        else
            echo "MISMATCH: $basename_file"
            FAIL=$((FAIL+1))
            ERRORS+=("MISMATCH: $basename_file")
            
            # Show diff (first 20 lines)
            diff <(echo "$sorted_expected") <(echo "$sorted_actual") | head -40 || true
            echo "---"
        fi
    fi
}

# If specific files passed, test them; otherwise test all Simple cases
if [[ $# -gt 0 ]]; then
    for f in "$@"; do
        run_test "$f"
    done
else
    for f in "$TESTDATA_DIR"/*.txt; do
        run_test "$f"
    done
fi

echo ""
echo "===== Results ====="
echo "PASS: $PASS  XFAIL: $XFAIL  FAIL: $FAIL  SKIP: $SKIP  XPASS: $XPASS"

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
