#!/usr/bin/env bash
# harness/lib/concordance_cpp.sh — Shared advisory vardictcpp third-arm vcfeval
# header scaffolding, N-REF/chr normalization, and extra.json writer for
# concordance_germline.sh and concordance_pair.sh.

concordance_cpp_run_third_arm() {
  local cpp_vcf="${1:-}"
  local truth_vcf="${2:-}"
  local java_vcf="${3:-}"
  local out_dir="${4:-}"
  local eval_fn="${5:-}"

  CPP_F1=""
  CPP_F1_NORM=""
  CPP_NOTE="SKIP: no cpp VCF supplied"

  if [ "${CPP_SKIP:-0}" = "1" ] || [ "${CPP_SKIP:-0}" = "true" ]; then
    CPP_NOTE="SKIP: disabled via CPP_SKIP=1"
    return 0
  fi
  if [ -z "$cpp_vcf" ]; then
    return 0
  fi
  if [ ! -f "$cpp_vcf" ]; then
    CPP_NOTE="SKIP: cpp VCF not found ($cpp_vcf)"
    return 0
  fi
  if ! command -v rtg &>/dev/null || ! command -v bcftools &>/dev/null; then
    CPP_NOTE="SKIP: rtg/bcftools missing for cpp arm"
    return 0
  fi

  # vardictcpp's native -vcf emitter carries no ##contig lines (Task-3 ruling),
  # which tabix requires for indexing. Complete missing ##contig lines from the
  # truth header while keeping ##fileformat on line 1 (scaffolding only).
  mkdir -p "$out_dir/cpp"
  if { { bcftools view -h "$truth_vcf" | grep '^##contig=<ID=' || true;
         grep '^##contig=<ID=' "$cpp_vcf" || true; } \
          | awk -F'[=,>]' '!seen[$3]++' > "$out_dir/cpp.contigs.hdr" \
      && { grep -m1 '^##fileformat=' "$cpp_vcf" || true; cat "$out_dir/cpp.contigs.hdr";
           grep -v '^##fileformat=' "$cpp_vcf"; } > "$out_dir/cpp.src.vcf" \
      && "$eval_fn" "cpp" "$out_dir/cpp.src.vcf" \
      && CPP_F1=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["f_measure"])' "$out_dir/cpp/summary.json"); } 2>"$out_dir/cpp.vcfeval.err"; then
    local first_chrom target_chr=0
    first_chrom=$(grep -v '^#' "$java_vcf" | head -n 1 | cut -f 1 || true)
    [[ "$first_chrom" =~ ^chr ]] && target_chr=1
    python3 - "$out_dir/cpp.src.vcf" "$out_dir/cpp.norm.src.vcf" "$target_chr" <<'PY'
import re, sys
src, dst, target_chr = sys.argv[1], sys.argv[2], sys.argv[3] == "1"
n_pattern = re.compile(r'^[Nn]+$')
with open(src, "r", encoding="utf-8") as fin, open(dst, "w", encoding="utf-8") as fout:
    for line in fin:
        if line.startswith("#"):
            fout.write(line)
            continue
        parts = line.split("\t")
        if len(parts) >= 5:
            if n_pattern.match(parts[3]):
                continue
            chrom = parts[0]
            if target_chr and not chrom.startswith("chr"):
                parts[0] = "chr" + chrom
            elif not target_chr and chrom.startswith("chr"):
                parts[0] = chrom[3:]
            fout.write("\t".join(parts))
        else:
            fout.write(line)
PY
    if cmp -s "$out_dir/cpp.src.vcf" "$out_dir/cpp.norm.src.vcf"; then
      CPP_F1_NORM="$CPP_F1"
    else
      mkdir -p "$out_dir/cpp_norm"
      if "$eval_fn" "cpp_norm" "$out_dir/cpp.norm.src.vcf" 2>"$out_dir/cpp_norm.vcfeval.err"; then
        if ! CPP_F1_NORM=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["f_measure"])' "$out_dir/cpp_norm/summary.json" 2>/dev/null); then
          CPP_F1_NORM=""
          echo "WARNING: failed to parse cpp_norm/summary.json; omitting cpp_f1_normalized" >&2
        fi
      else
        CPP_F1_NORM=""
        echo "WARNING: cpp_norm evaluation failed; omitting cpp_f1_normalized (see $out_dir/cpp_norm.vcfeval.err)" >&2
      fi
    fi
  else
    local first_err
    first_err=$(grep -vE '^(ERROR conda\.cli|cpp_chunk=|[[:space:]]*$)' "$out_dir/cpp.vcfeval.err" 2>/dev/null | head -n 1 | tr -d '\r' | cut -c 1-140 || true)
    CPP_NOTE="SKIP: cpp vcfeval failed (${first_err:-see $out_dir/cpp.vcfeval.err})"
    rm -rf "$out_dir/cpp" "$out_dir/cpp.stdout" "$out_dir/cpp.src.vcf" "$out_dir/cpp.contigs.hdr"
  fi
}

concordance_cpp_write_extra() {
  local out_dir="$1"
  local cpp_f1="${2:-}"
  local cpp_note="${3:-SKIP: no cpp VCF supplied}"
  local cpp_f1_norm="${4:-}"
  local cpp_ver
  cpp_ver=$(conda list -n vardict-oracle vardictcpp 2>/dev/null | awk '/^vardictcpp[[:space:]]/{print "vardictcpp " $2 " " $3 " bioconda"}' || true)
  [ -n "$cpp_ver" ] || cpp_ver="vardictcpp 3 hc52dbad_0 bioconda"
  python3 - "$out_dir" "$cpp_f1" "$cpp_note" "$cpp_ver" "$cpp_f1_norm" <<'PY'
import json, sys
o, f1, note, ver = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
f1_norm = sys.argv[5] if len(sys.argv) > 5 else ""
extra = {"java": json.load(open(o + "/java/summary.json")),
         "rust": json.load(open(o + "/rust/summary.json"))}
if f1:
    extra["cpp_f1"] = float(f1)
    if f1_norm:
        extra["cpp_f1_normalized"] = float(f1_norm)
    extra["cpp_version"] = ver
    j, r, c = extra["java"]["f_measure"], extra["rust"]["f_measure"], extra["cpp_f1"]
    print(f"CPP: f1={c} (java={j} rust={r}) rank {1 + (c < j) + (c < r)}/3 -- advisory, gate stays java==rust EQUAL")
else:
    extra["cpp_skip_reason"] = note
    print(f"CPP: {note}")
json.dump(extra, open(o + "/extra.json", "w"), indent=2, sort_keys=True)
PY
}
