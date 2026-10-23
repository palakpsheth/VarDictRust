#!/usr/bin/env bash
# Usage: concordance_pair.sh <name> <java.vcf> <rust.vcf> <truth.vcf.gz> <eval.bed> <ref.sdf> <outdir> <tumor> [cpp.vcf]
# Gate (D-27): identical 'None'-row counts AND identical tp/tp-baseline/fp/fn site lists. Not a dF1 tolerance.
# Optional 9th arg adds the advisory vardictcpp third arm (pinned vardictcpp 3
# hc52dbad_0 bioconda): a second concordance.sh vcfeval run with identical
# truth/BED/sdf. cpp F1 lands in extra.json (cpp_f1 + cpp_version) and stdout
# as a rank line; the gate and exit code stay java==rust only (Task-3 ruling:
# header/END-key deltas are known-normalized, cpp rank never gated).
set -euo pipefail

if ! command -v rtg &>/dev/null || ! command -v bcftools &>/dev/null; then
  export PATH="/home/palak/miniconda3/envs/g4-tools/bin:$PATH"
fi

NAME=$1; J=$2; R=$3; TRUTH=$4; BED=$5; SDF=$6; OUT=$7; TS=$8; CPP_VCF=${9:-}
bash harness/concordance.sh "$J" "$TRUTH" "$BED" "$SDF" "$OUT/java" "$TS" >/dev/null
bash harness/concordance.sh "$R" "$TRUTH" "$BED" "$SDF" "$OUT/rust" "$TS" >/dev/null
# Third arm: vardictcpp (advisory only; any failure degrades to a loud SKIP —
# allowed only for the cpp arm; java/rust behavior above is unchanged).
CPP_F1=""; CPP_F1_NORM=""; CPP_NOTE="SKIP: no cpp VCF supplied"
if [ "${CPP_SKIP:-0}" = "1" ] || [ "${CPP_SKIP:-0}" = "true" ]; then
  CPP_NOTE="SKIP: disabled via CPP_SKIP=1"
elif [ -n "$CPP_VCF" ]; then
  if [ ! -f "$CPP_VCF" ]; then
    CPP_NOTE="SKIP: cpp VCF not found ($CPP_VCF)"
  elif ! command -v rtg &>/dev/null || ! command -v bcftools &>/dev/null; then
    CPP_NOTE="SKIP: rtg/bcftools missing for cpp arm"
  else
    # vardictcpp's native -vcf emitter carries no ##contig lines (Task-3 ruling)
    # but tabix (inside concordance.sh) requires them for the calls index:
    # complete the missing ones from the truth header (truth first, so it wins
    # the dedup; ##fileformat stays line 1 as VCF requires). Scaffolding only —
    # no record, FILTER, QUAL or INFO is touched, so the vcfeval F1 is
    # unaffected. (Twin block in concordance_germline.sh.)
    mkdir -p "$OUT/cpp"
    if { { bcftools view -h "$TRUTH" | grep '^##contig=<ID=' || true;
           grep '^##contig=<ID=' "$CPP_VCF" || true; } \
            | awk -F'[=,>]' '!seen[$3]++' > "$OUT/cpp.contigs.hdr" \
        && { grep -m1 '^##fileformat=' "$CPP_VCF" || true; cat "$OUT/cpp.contigs.hdr";
             grep -v '^##fileformat=' "$CPP_VCF"; } > "$OUT/cpp.src.vcf" \
        && bash harness/concordance.sh "$OUT/cpp.src.vcf" "$TRUTH" "$BED" "$SDF" "$OUT/cpp" "$TS" >"$OUT/cpp.stdout" \
        && CPP_F1=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["f_measure"])' "$OUT/cpp/summary.json"); } 2>"$OUT/cpp.vcfeval.err"; then
      first_chrom=$(grep -v '^#' "$J" | head -n 1 | cut -f 1 || true)
      target_chr=0
      [[ "$first_chrom" =~ ^chr ]] && target_chr=1
      python3 - "$OUT/cpp.src.vcf" "$OUT/cpp.norm.src.vcf" "$target_chr" <<'PY'
import sys, re
src, dst, target_chr = sys.argv[1], sys.argv[2], sys.argv[3] == "1"
n_pattern = re.compile(r'^[Nn]+$')
with open(src, "r") as fin, open(dst, "w") as fout:
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
      if cmp -s "$OUT/cpp.src.vcf" "$OUT/cpp.norm.src.vcf"; then
        CPP_F1_NORM="$CPP_F1"
      else
        mkdir -p "$OUT/cpp_norm"
        if bash harness/concordance.sh "$OUT/cpp.norm.src.vcf" "$TRUTH" "$BED" "$SDF" "$OUT/cpp_norm" "$TS" >"$OUT/cpp_norm.stdout" 2>"$OUT/cpp_norm.vcfeval.err"; then
          if ! CPP_F1_NORM=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["f_measure"])' "$OUT/cpp_norm/summary.json" 2>/dev/null); then
            CPP_F1_NORM=""
            echo "WARNING: failed to parse cpp_norm/summary.json; omitting cpp_f1_normalized" >&2
          fi
        else
          CPP_F1_NORM=""
          echo "WARNING: cpp_norm concordance failed; omitting cpp_f1_normalized (see $OUT/cpp_norm.vcfeval.err)" >&2
        fi
      fi
    else
      first_err=$(grep -vE '^(ERROR conda\.cli|cpp_chunk=|[[:space:]]*$)' "$OUT/cpp.vcfeval.err" 2>/dev/null | head -n 1 | tr -d '\r' | cut -c 1-140 || true)
      CPP_NOTE="SKIP: cpp vcfeval failed (${first_err:-see $OUT/cpp.vcfeval.err})"
      rm -rf "$OUT/cpp" "$OUT/cpp.stdout" "$OUT/cpp.src.vcf" "$OUT/cpp.contigs.hdr"
    fi
  fi
fi
res=EQUAL; cmp -s "$OUT/java/summary.json" "$OUT/rust/summary.json" || res=DIFFERENT
for f in tp tp-baseline fp fn; do cmp -s "$OUT/java/$f.sites" "$OUT/rust/$f.sites" || res=DIFFERENT; done
# extra.json: java/rust summaries (unchanged) + cpp_f1/cpp_version on success,
# else loud cpp_skip_reason. The advisory CPP rank line is stdout only.
CPP_VER=$(conda list -n vardict-oracle vardictcpp 2>/dev/null | awk '/^vardictcpp[[:space:]]/{print "vardictcpp " $2 " " $3 " bioconda"}' || true)
[ -n "$CPP_VER" ] || CPP_VER="vardictcpp 3 hc52dbad_0 bioconda"
python3 - "$OUT" "$CPP_F1" "$CPP_NOTE" "$CPP_VER" "$CPP_F1_NORM" <<'PY'
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
python3 harness/record_evidence.py "$NAME" "$res" "$OUT" -- "$0" "$@"
[ "$res" = EQUAL ]
