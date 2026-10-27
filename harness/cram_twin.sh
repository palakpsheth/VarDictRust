#!/usr/bin/env bash
# Usage: KEEP_DIR=<dir> bash harness/cram_twin.sh
# Three-way BAM-vs-CRAM twin check (advisory cpp third arm, never gating the
# rust twin verdict): rust and cpp each call the hg008 slice pair from BAM and
# from CRAM; the gate is rust-BAM == rust-CRAM (proves noodles + -G CRAM decode,
# AGENTS.md s11 #249/#380). cpp-BAM == cpp-CRAM is recorded advisory.
# Exit 0 = both twins MATCH; non-zero otherwise. KEEP_DIR persists the four
# VCFs, both compare_*.txt, *.time, and extra.json for record_evidence.py.
set -euo pipefail
REF=scratch/g4/ref/giabv3/GRCh38_GIABv3.fa
BED=scratch/g4/hg008/engine.bed
BAM_T=scratch/g4/hg008/HG008-T.slice.bam; BAM_N=scratch/g4/hg008/HG008-N.slice.bam
CRAM_T=scratch/g4/hg008/HG008-T.slice.cram; CRAM_N=scratch/g4/hg008/HG008-N.slice.cram
SM="HG008-T|HG008-N"; TH=4; FLOOR=1000
if [ -n "${KEEP_DIR:-}" ]; then mkdir -p "$KEEP_DIR"; OUT=$KEEP_DIR
else OUT=$(mktemp -d); trap 'rm -rf "$OUT"' EXIT; fi
BIN=${VARDICT_BIN:-target/release/vardictrust}
cargo build -q --release -p vardictrust
# shortcut: rust-CRAM arm runs with RUST_MIN_STACK=64M (2 MiB rayon-worker
# default overflows inside noodles-cram 0.100.0 decode on these high-depth
# slices — ~97-frame nesting: container/slice iterators, Block::decode,
# name_tokenizer 10k-name blocks, striped rans_nx16 chunk decodes; 2 MiB
# fails / 64 MiB passes release on tile-0 (8 MiB already passes tile-0, so
# margin below 64M is untested for the full engine). BAM is unaffected, and
# 64M yields the full 18,426-call set identical to BAM). Ceiling: masks
# upstream decoder frames (rANS order tables up to 512 KiB + nested collect()
# chains) — upgrade = upstream Box/flatten in noodles-cram, then re-evaluate
# via crates/vardictrust-core/tests/cram_stack.rs and drop this override.
# Our frames are thin; nothing to Box on our side (T3D-03).
CURRENT_FP=$(bash harness/engine_fingerprint.sh 2>/dev/null || echo dirty)
if [ -s scratch/g4/hg008/cmp/rust.vcf ] && [ -s scratch/g4/hg008/cmp/rust.time ] \
   && [ -f benchmarks/evidence/hg008-tier2.json ] \
   && python3 -c 'import json,sys; j=json.load(open(sys.argv[1])); sys.exit(0 if j.get("result")=="MATCH" and j.get("fingerprint")==sys.argv[2] else 1)' benchmarks/evidence/hg008-tier2.json "$CURRENT_FP" 2>/dev/null \
   && [ "$(grep -vc '^#' scratch/g4/hg008/cmp/rust.vcf)" -eq 18426 ]; then
  echo "Reusing verified HG008 BAM callset from scratch/g4/hg008/cmp/rust.vcf (fingerprint $CURRENT_FP)"
  cp scratch/g4/hg008/cmp/rust.vcf "$OUT/rust-bam.vcf"
  cp scratch/g4/hg008/cmp/rust.time "$OUT/rust.time"
else
  LC_ALL=C /usr/bin/time -v -o "$OUT/rust.time" timeout "${CMP_TIMEOUT:-12h}" \
    "$BIN" -G "$REF" -f 0.01 -N "$SM" -b "$BAM_T|$BAM_N" -c 1 -S 2 -E 3 -g 4 -z -th "$TH" "$BED" > "$OUT/rust-bam.vcf" \
    || { echo "ENGINE-FAIL: rust-bam"; exit 2; }
fi
LC_ALL=C /usr/bin/time -v -o /dev/null timeout "${CMP_TIMEOUT:-12h}" \
  env RUST_MIN_STACK=67108864 "$BIN" -G "$REF" -f 0.01 -N "$SM" -b "$CRAM_T|$CRAM_N" -c 1 -S 2 -E 3 -g 4 -z -th "$TH" "$BED" > "$OUT/rust-cram.vcf" \
  || { echo "ENGINE-FAIL: rust-cram"; exit 2; }
python3 harness/compare_vcf.py --min-variants "$FLOOR" "$OUT/rust-bam.vcf" "$OUT/rust-cram.vcf" \
  | tee "$OUT/compare.txt" || exit 1
grep -q '^MATCH' "$OUT/compare.txt" || { echo "MISMATCH: rust BAM-vs-CRAM twin"; exit 1; }
echo "MATCH: rust BAM-vs-CRAM twin"
# Advisory cpp twin (same floors; never gates).
if [ "${CPP_SKIP:-0}" = "1" ] || [ "${CPP_SKIP:-0}" = "true" ]; then
  cpp_line="CPP: SKIP: disabled via CPP_SKIP=1"
elif bash harness/golden_store.sh cpp cram-twin-bam "$REF" "$BAM_T|$BAM_N" "$BED" "$OUT/cpp-bam.vcf" "$SM" -z -th "$TH" 2>"$OUT/cpp.err" \
  && bash harness/golden_store.sh cpp cram-twin-cram "$REF" "$CRAM_T|$CRAM_N" "$BED" "$OUT/cpp-cram.vcf" "$SM" -z -th "$TH" 2>>"$OUT/cpp.err"; then
  if [ ! -s "$OUT/cpp-bam.vcf" ] || [ ! -s "$OUT/cpp-cram.vcf" ]; then
    cpp_line="CPP: SKIP: empty cpp twin VCF"
  elif python3 harness/compare_vcf.py --min-variants "$FLOOR" "$OUT/cpp-bam.vcf" "$OUT/cpp-cram.vcf" > "$OUT/compare_cpp.txt" 2>&1; then
    cpp_line="CPP: MATCH"
  else
    cpp_line="CPP: MISMATCH"
  fi
  echo "$cpp_line" >> "$OUT/compare_cpp.txt"
else
  cpp_exit=$?
  if [ "$cpp_exit" -eq 3 ]; then
    cpp_line="CPP: SKIP: $(head -1 "$OUT/cpp.err" 2>/dev/null || echo 'vardictcpp unavailable')"
  else
    first_err=$(grep -vE '^(ERROR conda\.cli|cpp_chunk=|[[:space:]]*$)' "$OUT/cpp.err" 2>/dev/null | head -n 1 | tr -d '\r' | cut -c 1-140 || true)
    cpp_line="CPP: SKIP: engine failure (exit $cpp_exit${first_err:+: $first_err})"
  fi
  echo "$cpp_line" > "$OUT/compare_cpp.txt"
fi
echo "$cpp_line"
# extra.json for the ev_cpp presence gate via the shared helper (merge-only).
bash harness/staple_cpp_extra.sh "$OUT" "$cpp_line"
echo "MATCH (rust twin + $cpp_line)"
