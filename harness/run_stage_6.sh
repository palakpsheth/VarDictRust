#!/usr/bin/env bash
# =============================================================================
# Stage 6: Literature-Gap Closure Suite
#   1. s6_hg002_v50q_cmrg_mosaic_eval
#   2. s6_hg008_modern_ml_and_commercial_callers
#   3. s6_low_vaf_titration_truth_curve
#   4. s6_ga4gh_v35_stratification_matrix
# =============================================================================
set -euo pipefail

ensure_vcf_clean_helper() {
  local helper_py="${WORK_DIR:-/var/tmp/vardict-paper}/vcf_clean_chr1_22.py"
  mkdir -p "$(dirname "$helper_py")"
  if [ ! -s "$helper_py" ]; then
    cat > "$helper_py" <<'PY'
import sys

fai_path = sys.argv[1]
mode = sys.argv[2]  # germline | somatic | somatic_vardict
valid_chroms = {f"chr{i}" for i in range(1, 23)}
contig_lens = {}
with open(fai_path) as f:
    for line in f:
        parts = line.strip().split("\t")
        if parts[0] in valid_chroms:
            contig_lens[parts[0]] = parts[1]

wrote_contigs = False
has_gt_header = False
sample_name = "SAMPLE"

for raw in sys.stdin:
    if raw.startswith("##"):
        if raw.startswith("##contig="):
            continue
        if raw.startswith("##FORMAT=<ID=GT,"):
            has_gt_header = True
        sys.stdout.write(raw)
        continue
    if raw.startswith("#CHROM"):
        if not wrote_contigs:
            for c in sorted(contig_lens, key=lambda x: int(x[3:])):
                sys.stdout.write(f"##contig=<ID={c},length={contig_lens[c]}>\n")
            wrote_contigs = True
        if not has_gt_header:
            sys.stdout.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        cols = raw.rstrip("\n").split("\t")
        if len(cols) >= 10:
            sample_name = cols[9]
        sys.stdout.write("\t".join(cols[:8] + ["FORMAT", sample_name]) + "\n")
        continue
    cols = raw.rstrip("\n").split("\t")
    if len(cols) < 8:
        continue
    chrom, pos, vid, ref, alt, qual, flt, info = cols[:8]
    if chrom not in valid_chroms:
        continue
    if flt not in ("PASS", "."):
        continue
    if mode == "somatic_vardict":
        if "STATUS=StrongSomatic" not in info and "STATUS=LikelySomatic" not in info:
            continue
    if alt.startswith("<") or any(ch in alt for ch in "[]"):
        continue
    gt = "0/1"
    if mode == "germline" and len(cols) >= 10:
        fmt_keys = cols[8].split(":")
        smp_vals = cols[9].split(":")
        if "GT" in fmt_keys:
            idx = fmt_keys.index("GT")
            if idx < len(smp_vals) and smp_vals[idx] not in (".", "./.", ".|."):
                gt = smp_vals[idx]
    sys.stdout.write("\t".join([chrom, pos, vid, ref, alt, qual, "PASS", info, "GT", gt]) + "\n")
PY
  fi
}

eval_single_vcf_against_bed() {
  local vcf_in="$1"
  local ref="$2"
  local sdf="$3"
  local truth_vcf="$4"
  local eval_bed="$5"
  local mode="$6" # germline | somatic | somatic_vardict
  local work_dir="$7"
  local pre_normalized="${8:-0}"

  ensure_vcf_clean_helper
  local helper_py="${WORK_DIR:-/var/tmp/vardict-paper}/vcf_clean_chr1_22.py"

  if [ -f "$work_dir/rtg_snp/summary.txt" ] && [ -f "$work_dir/rtg_indel/summary.txt" ]; then
    local snp_p="NA" snp_r="NA" snp_f="NA" indel_p="NA" indel_r="NA" indel_f="NA"
    for cls in snp indel; do
      read -r p r f <<< "$(awk 'NR>2 && NF>=8' "$work_dir/rtg_${cls}/summary.txt" | sort -k8 -g | tail -1 | awk '{print $6,$7,$8}')"
      if [ -n "${f:-}" ]; then
        eval "${cls}_p=$p; ${cls}_r=$r; ${cls}_f=$f"
      fi
    done
    printf "%s\t%s\t%s\t%s\t%s\t%s" "$snp_p" "$snp_r" "$snp_f" "$indel_p" "$indel_r" "$indel_f"
    return 0
  fi

  rm -rf "$work_dir"
  mkdir -p "$work_dir"

  local calls_gz="$work_dir/calls.vcf.gz"
  if [ "$pre_normalized" = "1" ] && [ -s "$vcf_in" ]; then
    cp -f "$vcf_in" "$calls_gz"
    if [ -s "${vcf_in}.tbi" ]; then
      cp -f "${vcf_in}.tbi" "${calls_gz}.tbi"
    else
      tabix -f -p vcf "$calls_gz"
    fi
  else
    zcat -f "$vcf_in" 2>/dev/null \
      | python3 "$helper_py" "${ref}.fai" "$mode" \
      | bcftools norm -c x -f "$ref" -m -any 2>/dev/null \
      | bcftools sort -Oz -o "$calls_gz"
    tabix -f -p vcf "$calls_gz"
  fi

  local tr_snp="${truth_vcf}.norm_snp.vcf.gz"
  local tr_indel="${truth_vcf}.norm_indel.vcf.gz"
  if [ ! -s "$tr_snp" ] || [ ! -s "${tr_snp}.tbi" ] || [ ! -s "$tr_indel" ] || [ ! -s "${tr_indel}.tbi" ]; then
    local tr_gz="$work_dir/truth.vcf.gz"
    bcftools annotate -x FILTER "$truth_vcf" 2>/dev/null \
      | bcftools norm -c x -f "$ref" -m -any 2>/dev/null \
      | bcftools sort -Oz -o "$tr_gz"
    tabix -f -p vcf "$tr_gz"
    bcftools view -v snp "$tr_gz" -Oz -o "$tr_snp"
    tabix -f -p vcf "$tr_snp"
    bcftools view -v indel "$tr_gz" -Oz -o "$tr_indel"
    tabix -f -p vcf "$tr_indel"
  fi

  local sample_args=()
  if [ "$mode" != "germline" ]; then
    sample_args=(--sample "ALT,ALT" --vcf-score-field QUAL)
  fi

  local c_snp="$work_dir/c_snp.vcf.gz"
  local c_indel="$work_dir/c_indel.vcf.gz"
  if [ "$pre_normalized" = "1" ] && [ -s "${vcf_in}.snp.vcf.gz" ] && [ -s "${vcf_in}.indel.vcf.gz" ]; then
    c_snp="${vcf_in}.snp.vcf.gz"
    c_indel="${vcf_in}.indel.vcf.gz"
  else
    bcftools view -v snp "$calls_gz" -Oz -o "$c_snp"
    tabix -f -p vcf "$c_snp"
    bcftools view -v indel "$calls_gz" -Oz -o "$c_indel"
    tabix -f -p vcf "$c_indel"
    if [ "$pre_normalized" = "1" ]; then
      cp -f "$c_snp" "${vcf_in}.snp.vcf.gz" && cp -f "${c_snp}.tbi" "${vcf_in}.snp.vcf.gz.tbi"
      cp -f "$c_indel" "${vcf_in}.indel.vcf.gz" && cp -f "${c_indel}.tbi" "${vcf_in}.indel.vcf.gz.tbi"
    fi
  fi

  local snp_p="NA" snp_r="NA" snp_f="NA" indel_p="NA" indel_r="NA" indel_f="NA"
  for cls in snp indel; do
    local c_cls="$c_snp"
    local t_cls="$tr_snp"
    if [ "$cls" = "indel" ]; then
      c_cls="$c_indel"
      t_cls="$tr_indel"
    fi
    local n_tr
    n_tr=$(zgrep -vc '^#' "$t_cls" || true)
    if [ "${n_tr:-0}" -gt 0 ]; then
      rtg vcfeval -b "$t_cls" -c "$c_cls" \
        -t "$sdf" -e "$eval_bed" --squash-ploidy "${sample_args[@]}" -T 16 \
        -o "$work_dir/rtg_${cls}" >/dev/null 2>&1 || true
      if [ -f "$work_dir/rtg_${cls}/summary.txt" ]; then
        read -r p r f <<< "$(awk 'NR>2 && NF>=8' "$work_dir/rtg_${cls}/summary.txt" | sort -k8 -g | tail -1 | awk '{print $6,$7,$8}')"
        if [ -n "${f:-}" ]; then
          eval "${cls}_p=$p; ${cls}_r=$r; ${cls}_f=$f"
        fi
      fi
    fi
  done
  printf "%s\t%s\t%s\t%s\t%s\t%s" "$snp_p" "$snp_r" "$snp_f" "$indel_p" "$indel_r" "$indel_f"
}

run_s6_vm2() {
  # ---------------------------------------------------------------------------
  # Cell 1: s6_hg002_v50q_cmrg_mosaic_eval
  # ---------------------------------------------------------------------------
  run_checkpointed_cell "stage6" "s6_hg002_v50q_cmrg_mosaic_eval" 32 "115G" bash -c '
    set -euo pipefail
    source harness/gcp_paper_runner.sh --source-only
    source harness/run_stage_6.sh
    ensure_vcf_clean_helper

    REF="$DATA_DIR/ref/giabv3/GRCh38_GIABv3.fasta"
    SDF="$DATA_DIR/ref/giabv3/GRCh38_GIABv3.sdf"
    TR_DIR="$DATA_DIR/hg002/truth_s6"
    mkdir -p "$TR_DIR" "$DATA_DIR/hg002/truth"

    # Pull cached Stage-6 HG002 truth sets (including pre-normalized .norm_snp/.norm_indel VCFs) from GCS
    if [ -z "${LOCAL_CKPT_DIR:-}" ]; then
      setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -q rsync -r "${GCS_BAMS}/g4-cache/hg002/truth_s6" "$TR_DIR" 2>/dev/null || true
      setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -q rsync -r "${GCS_BAMS}/g4-cache/hg002/truth" "$DATA_DIR/hg002/truth" 2>/dev/null || true
    fi

    # 1. Download NIST GIAB HG002 v5.0q (2026-04-07), CMRG v1.00, and Mosaic SNV v1.10 if not cached
    BASE_HG2="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son"
    [ -s "$TR_DIR/HG002_GRCh38_v5.0q_smvar.vcf.gz" ] || {
      curl -sSLf --retry 5 "$BASE_HG2/v5.0q/HG002_GRCh38_v5.0q_smvar.vcf.gz" -o "$TR_DIR/HG002_GRCh38_v5.0q_smvar.vcf.gz"
      curl -sSLf --retry 5 "$BASE_HG2/v5.0q/HG002_GRCh38_v5.0q_smvar.vcf.gz.tbi" -o "$TR_DIR/HG002_GRCh38_v5.0q_smvar.vcf.gz.tbi"
      curl -sSLf --retry 5 "$BASE_HG2/v5.0q/HG002_GRCh38_v5.0q_smvar.benchmark.bed" -o "$TR_DIR/HG002_GRCh38_v5.0q_smvar.benchmark.bed"
    }
    [ -s "$TR_DIR/HG002_GRCh38_CMRG_smallvar_v1.00.vcf.gz" ] || {
      curl -sSLf --retry 5 "$BASE_HG2/CMRG_v1.00/GRCh38/SmallVariant/HG002_GRCh38_CMRG_smallvar_v1.00.vcf.gz" -o "$TR_DIR/HG002_GRCh38_CMRG_smallvar_v1.00.vcf.gz"
      curl -sSLf --retry 5 "$BASE_HG2/CMRG_v1.00/GRCh38/SmallVariant/HG002_GRCh38_CMRG_smallvar_v1.00.vcf.gz.tbi" -o "$TR_DIR/HG002_GRCh38_CMRG_smallvar_v1.00.vcf.gz.tbi"
      curl -sSLf --retry 5 "$BASE_HG2/CMRG_v1.00/GRCh38/SmallVariant/HG002_GRCh38_CMRG_smallvar_v1.00.bed" -o "$TR_DIR/HG002_GRCh38_CMRG_smallvar_v1.00.bed"
    }
    [ -s "$TR_DIR/HG002_GRCh38_MosaicSNVv1.1.vcf.gz" ] || {
      curl -sSLf --retry 5 "$BASE_HG2/mosaic_v1.10/GRCh38/SNV/HG002_GRCh38_MosaicSNVv1.1.vcf.gz" -o "$TR_DIR/HG002_GRCh38_MosaicSNVv1.1.vcf.gz"
      curl -sSLf --retry 5 "$BASE_HG2/mosaic_v1.10/GRCh38/SNV/HG002_GRCh38_MosaicSNVv1.1.vcf.gz.tbi" -o "$TR_DIR/HG002_GRCh38_MosaicSNVv1.1.vcf.gz.tbi"
      curl -sSLf --retry 5 "$BASE_HG2/mosaic_v1.10/GRCh38/SNV/HG002_GRCh38_MosaicSNVv1.1.bed" -o "$TR_DIR/HG002_GRCh38_MosaicSNVv1.1.bed"
    }
    [ -s "$DATA_DIR/hg002/truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" ] || {
      curl -sSLf --retry 5 "$BASE_HG2/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" -o "$DATA_DIR/hg002/truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
      curl -sSLf --retry 5 "$BASE_HG2/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi" -o "$DATA_DIR/hg002/truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi"
      curl -sSLf --retry 5 "$BASE_HG2/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" -o "$DATA_DIR/hg002/truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
    }

    # Filter BEDs to autosomes chr1..chr22 matching our full-WGS callsets
    awk -F"\t" '\''$1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/'\'' "$TR_DIR/HG002_GRCh38_v5.0q_smvar.benchmark.bed" > "$TR_DIR/v5.0q.chr1_22.bed"
    awk -F"\t" '\''$1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/'\'' "$TR_DIR/HG002_GRCh38_CMRG_smallvar_v1.00.bed" > "$TR_DIR/cmrg_v1.00.chr1_22.bed"
    awk -F"\t" '\''$1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/'\'' "$TR_DIR/HG002_GRCh38_MosaicSNVv1.1.bed" > "$TR_DIR/mosaic_v1.10.chr1_22.bed"

    # 2. Restore full-WGS HG002 VCFs from Stage 5 & cached VarDictJava normalized VCF / vcfeval from oracle-cache-v183
    WORK_S6="$WORK_DIR/tmp_s6_hg002"
    mkdir -p "$WORK_S6"
    if [ -z "${LOCAL_CKPT_DIR:-}" ]; then
      setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -q rsync -r "${GCS_ORACLE_CACHE}/stage6/hg002_v50q_cmrg_mosaic" "$WORK_S6" 2>/dev/null || true
    fi

    CELLS_SUB="${GCS_CELLS_SUBDIR:-cells-v106}"
    CONTIGS_VCF="$WORK_DIR/cells/s5_wgs_hg002_full_contigs_3way_th32/rust_contigs_th32.vcf.gz"
    GIAB_RUST_VCF="$WORK_DIR/cells/s5_wgs_hg002_full_giab_3way_th32/rust_giab_th32.vcf.gz"
    GIAB_JAVA_VCF="$WORK_DIR/cells/s5_wgs_hg002_full_giab_3way_th32/java_giab_th32.vcf.gz"
    [ -s "$CONTIGS_VCF" ] || setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -q cp "${GCS_RESULTS}/${CELLS_SUB}/s5_wgs_hg002_full_contigs_3way_th32/rust_contigs_th32.vcf.gz" "$WORK_S6/rust_contigs.vcf.gz"
    [ -s "$CONTIGS_VCF" ] && CONTIGS_VCF_USE="$CONTIGS_VCF" || CONTIGS_VCF_USE="$WORK_S6/rust_contigs.vcf.gz"
    [ -s "$GIAB_RUST_VCF" ] || setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -q cp "${GCS_RESULTS}/${CELLS_SUB}/s5_wgs_hg002_full_giab_3way_th32/rust_giab_th32.vcf.gz" "$WORK_S6/rust_giab.vcf.gz"
    [ -s "$GIAB_RUST_VCF" ] && GIAB_RUST_USE="$GIAB_RUST_VCF" || GIAB_RUST_USE="$WORK_S6/rust_giab.vcf.gz"
    if [ ! -s "$WORK_S6/VarDictJava_GIAB50kb.norm.vcf.gz" ]; then
      if [ ! -s "$GIAB_JAVA_VCF" ]; then
        setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -q cp "${GCS_ORACLE_CACHE}/hg002_full/java_giab_th32.vcf.gz" "$WORK_S6/java_giab.vcf.gz" || \
        setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -q cp "${GCS_RESULTS}/${CELLS_SUB}/s5_wgs_hg002_full_giab_3way_th32/java_giab_th32.vcf.gz" "$WORK_S6/java_giab.vcf.gz"
        GIAB_JAVA_USE="$WORK_S6/java_giab.vcf.gz"
      else
        GIAB_JAVA_USE="$GIAB_JAVA_VCF"
      fi
    else
      GIAB_JAVA_USE="$WORK_S6/VarDictJava_GIAB50kb.norm.vcf.gz"
    fi

    # Always refresh VarDictRust normalized VCFs for the current binary; reuse VarDictJava_GIAB50kb.norm.vcf.gz if cached
    if [ -x "target/release/vardictrust" ]; then
      for rlbl in VarDictRust_FullContigs VarDictRust_GIAB50kb; do
        if [ -f "$WORK_S6/${rlbl}.norm.vcf.gz" ] && [ ! "$WORK_S6/${rlbl}.norm.vcf.gz" -nt "target/release/vardictrust" ]; then
          rm -rf "$WORK_S6/${rlbl}.norm.vcf.gz"* "$WORK_S6"/ev_*_"${rlbl}"
        fi
      done
    fi

    # Pre-normalize each callset once in parallel to speed up multi-benchmark evaluation
    for pair in "VarDictRust_FullContigs:$CONTIGS_VCF_USE" "VarDictRust_GIAB50kb:$GIAB_RUST_USE" "VarDictJava_GIAB50kb:$GIAB_JAVA_USE"; do
      (
        lbl="${pair%%:*}"
        src="${pair#*:}"
        if [ ! -s "$WORK_S6/${lbl}.norm.vcf.gz" ] || [ ! -s "$WORK_S6/${lbl}.norm.vcf.gz.tbi" ]; then
          rm -f "$WORK_S6/${lbl}.norm.vcf.gz" "$WORK_S6/${lbl}.norm.vcf.gz.tbi" "$WORK_S6/${lbl}.norm.tmp.vcf.gz"
          zcat -f "$src" 2>/dev/null \
            | python3 "${WORK_DIR}/vcf_clean_chr1_22.py" "${REF}.fai" "germline" \
            | bcftools norm --threads 4 -c x -f "$REF" -m -any 2>/dev/null \
            | bcftools sort -m 8G -T "$WORK_S6/sort_${lbl}" -Oz -o "$WORK_S6/${lbl}.norm.tmp.vcf.gz"
          mv -f "$WORK_S6/${lbl}.norm.tmp.vcf.gz" "$WORK_S6/${lbl}.norm.vcf.gz"
          tabix -f -p vcf "$WORK_S6/${lbl}.norm.vcf.gz"
        fi
      ) &
    done
    wait

    OUT_TSV="$KEEP_DIR/hg002_v50q_cmrg_mosaic_summary.tsv"
    printf "benchmark\tcaller\tsnv_p\tsnv_r\tsnv_f1\tindel_p\tindel_r\tindel_f1\n" > "$OUT_TSV"

    for bm_spec in \
      "NIST_v4.2.1_All|$DATA_DIR/hg002/truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz|$DATA_DIR/hg002/truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" \
      "NIST_v5.0q_T2T_Q100|$TR_DIR/HG002_GRCh38_v5.0q_smvar.vcf.gz|$TR_DIR/v5.0q.chr1_22.bed" \
      "NIST_CMRG_v1.00_273Genes|$TR_DIR/HG002_GRCh38_CMRG_smallvar_v1.00.vcf.gz|$TR_DIR/cmrg_v1.00.chr1_22.bed" \
      "NIST_Mosaic_SNV_v1.10|$TR_DIR/HG002_GRCh38_MosaicSNVv1.1.vcf.gz|$TR_DIR/mosaic_v1.10.chr1_22.bed"; do
      IFS="|" read -r bm_name tr_vcf tr_bed <<< "$bm_spec"
      for lbl in VarDictRust_FullContigs VarDictRust_GIAB50kb VarDictJava_GIAB50kb; do
        metrics=$(eval_single_vcf_against_bed "$WORK_S6/${lbl}.norm.vcf.gz" "$REF" "$SDF" "$tr_vcf" "$tr_bed" "germline" "$WORK_S6/ev_${bm_name}_${lbl}" 1)
        printf "%s\t%s\t%s\n" "$bm_name" "$lbl" "$metrics" | tee -a "$OUT_TSV"
      done
    done

    # Persist truth_s6, normalized v4.2.1 truth, and VarDictJava_GIAB50kb normalized VCF + ev_* to GCS for future Rust-only runs
    if [ -z "${LOCAL_CKPT_DIR:-}" ]; then
      setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -o GSUtil:parallel_composite_upload_threshold=0 -q rsync -r "$TR_DIR" "${GCS_BAMS}/g4-cache/hg002/truth_s6" 2>/dev/null || true
      setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -o GSUtil:parallel_composite_upload_threshold=0 -q rsync -r "$DATA_DIR/hg002/truth" "${GCS_BAMS}/g4-cache/hg002/truth" 2>/dev/null || true
      for jf in "$WORK_S6"/VarDictJava_GIAB50kb.norm.vcf.gz*; do
        [ -s "$jf" ] && setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -o GSUtil:parallel_composite_upload_threshold=0 -q cp "$jf" "${GCS_ORACLE_CACHE}/stage6/hg002_v50q_cmrg_mosaic/" 2>/dev/null || true
      done
      for jdir in "$WORK_S6"/ev_*_VarDictJava_GIAB50kb; do
        [ -d "$jdir" ] && setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -o GSUtil:parallel_composite_upload_threshold=0 -q rsync -r "$jdir" "${GCS_ORACLE_CACHE}/stage6/hg002_v50q_cmrg_mosaic/$(basename "$jdir")" 2>/dev/null || true
      done
    fi
  '

  # ---------------------------------------------------------------------------
  # Cell 2: s6_hg008_modern_ml_and_commercial_callers
  # ---------------------------------------------------------------------------
  run_checkpointed_cell "stage6" "s6_hg008_modern_ml_and_commercial_callers" 32 "115G" bash -c '
    set -euo pipefail
    source harness/gcp_paper_runner.sh --source-only
    source harness/run_stage_6.sh
    ensure_vcf_clean_helper

    REF="$DATA_DIR/ref/giabv3/GRCh38_GIABv3.fasta"
    SDF="$DATA_DIR/ref/giabv3/GRCh38_GIABv3.sdf"
    TRUTH_V04="$DATA_DIR/hg008/truth/HG008-T_somatic_smvar_benchmark_v0.4_GRCh38_tumorvariants.vcf.gz"
    BED_V04_ALL="$DATA_DIR/hg008/truth/HG008-T_somatic_smvar_benchmark_v0.4_GRCh38_all.bed"
    BED_V04_NOGI="$DATA_DIR/hg008/truth/HG008-T_somatic_smvar_benchmark_v0.4_GRCh38_nogermlineinterference.bed"
    SLICE_BED="$DATA_DIR/hg008/eval.bed"

    WORK_ML="$WORK_DIR/tmp_s6_hg008_ml"
    EXT_DIR="$DATA_DIR/hg008/ext_callsets"
    mkdir -p "$WORK_ML" "$EXT_DIR"
    if [ -z "${LOCAL_CKPT_DIR:-}" ]; then
      setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -q rsync -r "${GCS_BAMS}/g4-cache/hg008/truth" "$DATA_DIR/hg008/truth" 2>/dev/null || true
      setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -q rsync -r "${GCS_BAMS}/g4-cache/hg008/ext_callsets" "$EXT_DIR" 2>/dev/null || true
      setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -q rsync -r "${GCS_ORACLE_CACHE}/stage6/hg008_ml_callers" "$WORK_ML" 2>/dev/null || true
    fi
    bedtools intersect -a "$SLICE_BED" -b "$BED_V04_ALL" > "$WORK_ML/slice_v04_all.bed"

    BASE_LISS="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/analysis"
    [ -s "$EXT_DIR/dragen_v4.2.4_ilmn_tumor.vcf.gz" ] || curl -sSLf --retry 5 "$BASE_LISS/DRAGEN-v4.2.4_ILMN-WGS_20240312/standard/dragen_4.2.4_HG008-mosaic_tumor.hard-filtered.vcf.gz" -o "$EXT_DIR/dragen_v4.2.4_ilmn_tumor.vcf.gz"
    [ -s "$EXT_DIR/dragen_v4.2.4_graph_ilmn.vcf.gz" ] || curl -sSLf --retry 5 "$BASE_LISS/DRAGEN-v4.2.4_ILMN-WGS_20240312/graph/dragen_4.2.4_HG008-mosaic_tumor.hard-filtered.vcf.gz" -o "$EXT_DIR/dragen_v4.2.4_graph_ilmn.vcf.gz"
    [ -s "$EXT_DIR/lancet2_ilmn.vcf.gz" ] || curl -sSLf --retry 5 "$BASE_LISS/NYGC_ILMN_Lancet2_20241016/HG008-T__vs__HG008-N-D.Lancet2.ml_scores.PASS.vcf.gz" -o "$EXT_DIR/lancet2_ilmn.vcf.gz"
    [ -s "$EXT_DIR/deepsomatic_v1.6.0_hifi.vcf.gz" ] || curl -sSLf --retry 5 "$BASE_LISS/UCSC_HiFi_DeepSomatic-SNV-INDEL_20240308/0201_HG008_HiFi_DeepSomatic_v1.6.0_model354_somaticOnly.vcf.gz" -o "$EXT_DIR/deepsomatic_v1.6.0_hifi.vcf.gz"
    if [ ! -s "$EXT_DIR/clairs_v0.4.1_hifi.vcf.gz" ]; then
      curl -sSLf --retry 5 "$BASE_LISS/HKU_HiFi_ClairS-v0.4.1-SNV-INDEL_20241201/HG008_T_N-P_HiFi_GRCh38-GIABv3_ClairS-SSRS_v0.4.1_snv.vcf.gz" -o "$WORK_ML/clairs_snv.vcf.gz"
      curl -sSLf --retry 5 "$BASE_LISS/HKU_HiFi_ClairS-v0.4.1-SNV-INDEL_20241201/HG008_T_N-P_HiFi_GRCh38-GIABv3_ClairS-SSRS_v0.4.1_indel.vcf.gz" -o "$WORK_ML/clairs_indel.vcf.gz"
      bcftools concat "$WORK_ML/clairs_snv.vcf.gz" "$WORK_ML/clairs_indel.vcf.gz" | bcftools sort -Oz -o "$EXT_DIR/clairs_v0.4.1_hifi.vcf.gz"
      tabix -f -p vcf "$EXT_DIR/clairs_v0.4.1_hifi.vcf.gz"
    fi
    [ -s "$EXT_DIR/ultima_deepvariant_somatic.vcf.gz" ] || curl -sSLf --retry 5 "$BASE_LISS/Ultima_DeepVariant-somatic-SNV-INDEL_20240822/HG008_edv_AF_recalibration.result.PASS.vcf.gz" -o "$EXT_DIR/ultima_deepvariant_somatic.vcf.gz"
    [ -s "$EXT_DIR/hartwig_sage_ilmn.vcf.gz" ] || curl -sSLf --retry 5 "$BASE_LISS/Hartwig_Oncoanalyzer_BCM-ILMN_20260204/sage/somatic/HG008-T.sage.somatic.vcf.gz" -o "$EXT_DIR/hartwig_sage_ilmn.vcf.gz"

    # Restore full-genome VarDictRust HG008 VCF from s5_wgs_hg008_full_giab_3way_th32 (always refresh for current release)
    CELLS_SUB="${GCS_CELLS_SUBDIR:-cells-v106}"
    RUST_FULL="$WORK_DIR/cells/s5_wgs_hg008_full_giab_3way_th32/rust_hg008_giab_th32.vcf.gz"
    if [ -s "$RUST_FULL" ]; then
      cp -f "$RUST_FULL" "$WORK_ML/VarDictRust_FullWGS.vcf.gz"
    else
      setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -q cp "${GCS_RESULTS}/${CELLS_SUB}/s5_wgs_hg008_full_giab_3way_th32/rust_hg008_giab_th32.vcf.gz" "$WORK_ML/VarDictRust_FullWGS.vcf.gz"
    fi
    rm -rf "$WORK_ML/norm_VarDictRust"* "$WORK_ML/ev_all_VarDictRust"* "$WORK_ML/ev_nogi_VarDictRust"* "$WORK_ML/ev_slice_VarDictRust"*

    OUT_TSV="$KEEP_DIR/hg008_modern_ml_commercial_callers_v04.tsv"
    printf "scope\tcaller\tchemistry\tsnv_p\tsnv_r\tsnv_f1\tindel_p\tindel_r\tindel_f1\n" > "$OUT_TSV"

    for spec in \
      "VarDictRust_v1.0.6|$WORK_ML/VarDictRust_FullWGS.vcf.gz|Illumina_195x_169x|somatic_vardict" \
      "DRAGEN_v4.2.4_Linear|$EXT_DIR/dragen_v4.2.4_ilmn_tumor.vcf.gz|Illumina_WGS|somatic" \
      "DRAGEN_v4.2.4_Graph|$EXT_DIR/dragen_v4.2.4_graph_ilmn.vcf.gz|Illumina_WGS|somatic" \
      "Lancet2_NYGC|$EXT_DIR/lancet2_ilmn.vcf.gz|Illumina_WGS|somatic" \
      "Hartwig_SAGE_2026|$EXT_DIR/hartwig_sage_ilmn.vcf.gz|Illumina_WGS|somatic" \
      "DeepSomatic_v1.6.0|$EXT_DIR/deepsomatic_v1.6.0_hifi.vcf.gz|PacBio_HiFi|somatic" \
      "ClairS_v0.4.1|$EXT_DIR/clairs_v0.4.1_hifi.vcf.gz|PacBio_HiFi|somatic" \
      "Ultima_DeepVariant_Somatic|$EXT_DIR/ultima_deepvariant_somatic.vcf.gz|Ultima_UG100|somatic"; do
      IFS="|" read -r cname cvcf chem cmode <<< "$spec"
      norm_vcf="$WORK_ML/norm_${cname}/calls.vcf.gz"
      if [ ! -s "$norm_vcf" ] || [ ! -s "${norm_vcf}.tbi" ]; then
        eval_single_vcf_against_bed "$cvcf" "$REF" "$SDF" "$TRUTH_V04" "$BED_V04_ALL" "$cmode" "$WORK_ML/norm_${cname}" 0 >/dev/null
      fi

      m_all=$(eval_single_vcf_against_bed "$norm_vcf" "$REF" "$SDF" "$TRUTH_V04" "$BED_V04_ALL" "somatic" "$WORK_ML/ev_all_${cname}" 1)
      printf "HG008_FullWGS_v0.4_All\t%s\t%s\t%s\n" "$cname" "$chem" "$m_all" | tee -a "$OUT_TSV"

      m_nogi=$(eval_single_vcf_against_bed "$norm_vcf" "$REF" "$SDF" "$TRUTH_V04" "$BED_V04_NOGI" "somatic" "$WORK_ML/ev_nogi_${cname}" 1)
      printf "HG008_FullWGS_v0.4_NoGI\t%s\t%s\t%s\n" "$cname" "$chem" "$m_nogi" | tee -a "$OUT_TSV"

      m_slice=$(eval_single_vcf_against_bed "$norm_vcf" "$REF" "$SDF" "$TRUTH_V04" "$WORK_ML/slice_v04_all.bed" "somatic" "$WORK_ML/ev_slice_${cname}" 1)
      printf "HG008_Slice_v0.4_All\t%s\t%s\t%s\n" "$cname" "$chem" "$m_slice" | tee -a "$OUT_TSV"
    done

    # Also score the CPU field callers from s5_grch38_field_callers_truth against HG008 v0.4 on the slice
    for fc in Mutect2 Strelka2 LoFreq; do
      fc_vcf="$WORK_DIR/cells/s5_grch38_field_callers_truth/eval_HG008_v0.3_${fc}/calls.vcf.gz"
      if [ ! -s "$fc_vcf" ]; then
        setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -q cp "${GCS_ORACLE_CACHE}/cells/s5_grch38_field_callers_truth/eval_HG008_v0.3_${fc}/calls.vcf.gz" "$WORK_ML/fc_${fc}.vcf.gz" 2>/dev/null || \
        setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -q cp "${GCS_RESULTS}/${CELLS_SUB}/s5_grch38_field_callers_truth/eval_HG008_v0.3_${fc}/calls.vcf.gz" "$WORK_ML/fc_${fc}.vcf.gz" 2>/dev/null || true
        fc_vcf="$WORK_ML/fc_${fc}.vcf.gz"
      fi
      if [ -s "$fc_vcf" ]; then
        m_fc=$(eval_single_vcf_against_bed "$fc_vcf" "$REF" "$SDF" "$TRUTH_V04" "$WORK_ML/slice_v04_all.bed" "somatic" "$WORK_ML/ev_slice_${fc}" 1)
        printf "HG008_Slice_v0.4_All\t%s\tIllumina_195x_169x\t%s\n" "$fc" "$m_fc" | tee -a "$OUT_TSV"
      fi
    done

    # Persist external callsets and non-Rust normalized VCFs + rtg vcfeval summaries to GCS
    if [ -z "${LOCAL_CKPT_DIR:-}" ]; then
      setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -o GSUtil:parallel_composite_upload_threshold=0 -q rsync -r "$EXT_DIR" "${GCS_BAMS}/g4-cache/hg008/ext_callsets" 2>/dev/null || true
      setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -o GSUtil:parallel_composite_upload_threshold=0 -q rsync -r "$DATA_DIR/hg008/truth" "${GCS_BAMS}/g4-cache/hg008/truth" 2>/dev/null || true
      setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -o GSUtil:parallel_composite_upload_threshold=0 -q rsync -r -x "(^|/).*VarDictRust.*(/|$)" "$WORK_ML" "${GCS_ORACLE_CACHE}/stage6/hg008_ml_callers" 2>/dev/null || true
    fi
  '

  # ---------------------------------------------------------------------------
  # Cell 3: s6_low_vaf_titration_truth_curve
  # ---------------------------------------------------------------------------
  run_checkpointed_cell "stage6" "s6_low_vaf_titration_truth_curve" 32 "115G" bash -c '
    set -euo pipefail
    source harness/gcp_paper_runner.sh --source-only
    source harness/run_stage_6.sh
    cell_s6_low_vaf_titration
  '
}

cell_s6_low_vaf_titration() {
  set -euo pipefail
  ensure_field_caller_envs
  ensure_vcf_clean_helper

  local REF="$DATA_DIR/ref/giabv3/GRCh38_GIABv3.fasta"
  local SDF="$DATA_DIR/ref/giabv3/GRCh38_GIABv3.sdf"
  [ -f "${REF%.*}.dict" ] || samtools dict "$REF" > "${REF%.*}.dict"
  [ -f "${REF}.dict" ]    || cp -f "${REF%.*}.dict" "${REF}.dict"

  local T_BAM="$DATA_DIR/hg008/HG008-T.slice.bam"
  local N_BAM="$DATA_DIR/hg008/HG008-N.slice.bam"
  local ENGINE_BED="$DATA_DIR/hg008/engine.bed"
  local EVAL_BED="$DATA_DIR/hg008/eval.bed"
  local TRUTH_V04="$DATA_DIR/hg008/truth/HG008-T_somatic_smvar_benchmark_v0.4_GRCh38_tumorvariants.vcf.gz"
  local BED_V04_ALL="$DATA_DIR/hg008/truth/HG008-T_somatic_smvar_benchmark_v0.4_GRCh38_all.bed"

  local WORK_TITR="$WORK_DIR/tmp_s6_titration"
  mkdir -p "$WORK_TITR"
  if [ -z "${LOCAL_CKPT_DIR:-}" ]; then
    setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -q rsync -r "${GCS_BAMS}/hg008_titration" "$WORK_TITR" 2>/dev/null || true
    setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -q rsync -r "${GCS_ORACLE_CACHE}/stage6/titration" "$WORK_TITR" 2>/dev/null || true
    cp -n "$WORK_TITR"/p*_*.time "$KEEP_DIR/" 2>/dev/null || true
  fi
  bedtools intersect -a "$EVAL_BED" -b "$BED_V04_ALL" > "$WORK_TITR/eval_v04.bed"
  bgzip -c "$EVAL_BED" > "$WORK_TITR/eval.bed.gz" && tabix -f -p bed "$WORK_TITR/eval.bed.gz"

  local OUT_TSV="$KEEP_DIR/low_vaf_titration_curve.tsv"
  printf "purity_pct\texpected_clonal_vaf_pct\tcaller\twall_s\tpeak_kb\tparity_vs_java\tsnv_p\tsnv_r\tsnv_f1\tindel_p\tindel_r\tindel_f1\n" > "$OUT_TSV"

  local spec pur vaf_exp ft fn ptag ADMIX_BAM parity_rust parity_cpp cspec cname cvcf ctime cpar cmode wall rss metrics
  for spec in \
    "20|10.0|1538|7101|p50" \
    "10|5.0|0769|7988|p25" \
    "5|2.5|0385|8432|p12" \
    "2|1.0|0154|8698|p06" \
    "1|0.5|0077|8787|p03"; do
    IFS="|" read -r pur vaf_exp ft fn ptag <<< "$spec"
    echo "=== [Titration] Purity=${pur}% (Expected Clonal VAF ~${vaf_exp}%) ==="
    ADMIX_BAM="$WORK_TITR/hg008_T_p${pur}.bam"
    if [ ! -s "$ADMIX_BAM" ] && [ -s "$WORK_TITR/admix_${ptag}.bam" ]; then
      ln -sf "admix_${ptag}.bam" "$ADMIX_BAM"
      [ -s "$WORK_TITR/admix_${ptag}.bam.bai" ] && ln -sf "admix_${ptag}.bam.bai" "${ADMIX_BAM}.bai"
    fi
    if [ ! -s "$ADMIX_BAM" ]; then
      samtools view -@ 16 -b -s "42.${ft}" "$T_BAM" > "$WORK_TITR/t_${pur}.bam"
      samtools view -@ 16 -b -s "84.${fn}" "$N_BAM" > "$WORK_TITR/n_${pur}.bam"
      samtools merge -@ 16 -f "$WORK_TITR/merged_${pur}.bam" "$WORK_TITR/t_${pur}.bam" "$WORK_TITR/n_${pur}.bam"
      samtools addreplacerg -r $'@RG\tID:HG008-T\tSM:HG008-T' -@ 16 -o "$ADMIX_BAM" "$WORK_TITR/merged_${pur}.bam"
      samtools index -@ 16 "$ADMIX_BAM"
      rm -f "$WORK_TITR/t_${pur}.bam" "$WORK_TITR/n_${pur}.bam" "$WORK_TITR/merged_${pur}.bam"
      if [ -z "${LOCAL_CKPT_DIR:-}" ]; then
        setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -o GSUtil:parallel_composite_upload_threshold=0 -q cp "$ADMIX_BAM" "${ADMIX_BAM}.bai" "${GCS_BAMS}/hg008_titration/" 2>/dev/null || true
      fi
    fi

    # Map ptag (p50/p25/p12/p06/p03) cached VCFs/times if restored under ptag naming
    for ext_c in java.vcf cpp.vcf strelka2.vcf.gz strelka2.vcf.gz.tbi lofreq.vcf.gz lofreq.vcf.gz.tbi mutect2.vcf.gz mutect2.vcf.gz.tbi; do
      if [ ! -s "$WORK_TITR/p${pur}_${ext_c}" ] && [ -s "$WORK_TITR/${ptag}_${ext_c}" ]; then
        cp -f "$WORK_TITR/${ptag}_${ext_c}" "$WORK_TITR/p${pur}_${ext_c}"
      fi
    done
    for ext_t in java.time cpp.time strelka2.time lofreq.time mutect2.time; do
      if [ ! -s "$KEEP_DIR/p${pur}_${ext_t}" ] && [ -s "$KEEP_DIR/${ptag}_${ext_t}" ]; then
        cp -f "$KEEP_DIR/${ptag}_${ext_t}" "$KEEP_DIR/p${pur}_${ext_t}"
      fi
    done

    # 1. VarDictRust (-f 0.002 -th 16; always run fresh for current VarDictRust binary)
    rm -rf "$WORK_TITR/ev_p${pur}_VarDictRust"
    /usr/bin/time -v -o "$KEEP_DIR/p${pur}_rust.time" \
      target/release/vardictrust -G "$REF" -f 0.002 -N "HG008-T|HG008-N-D" \
      -b "${ADMIX_BAM}|${N_BAM}" -c 1 -S 2 -E 3 -g 4 -th 16 "$ENGINE_BED" > "$WORK_TITR/p${pur}_rust.vcf"

    # 2. VarDictJava v1.8.3 (-f 0.002 -th 16)
    if [ ! -s "$WORK_TITR/p${pur}_java.vcf" ]; then
      /usr/bin/time -v -o "$KEEP_DIR/p${pur}_java.time" \
        env ORACLE_XMX_GB=32 bash harness/run_oracle.sh "$REF" "${ADMIX_BAM}|${N_BAM}" "$ENGINE_BED" "$WORK_TITR/p${pur}_java.vcf" "HG008-T|HG008-N-D" -f 0.002 -th 16 # golden:bypass gcp-oracle-cache-v183
    fi

    # 3. vardictcpp (-f 0.002 -th 16)
    if [ ! -s "$WORK_TITR/p${pur}_cpp.vcf" ]; then
      /usr/bin/time -v -o "$KEEP_DIR/p${pur}_cpp.time" \
        env CPP_CHUNK=50000 bash harness/run_cpp.sh "$REF" "${ADMIX_BAM}|${N_BAM}" "$ENGINE_BED" "$WORK_TITR/p${pur}_cpp.vcf" "HG008-T|HG008-N-D" -f 0.002 -th 16 # golden:bypass gcp-oracle-cache-v183
    fi

    parity_rust="MISMATCH"
    if python3 harness/compare_vcf.py "$WORK_TITR/p${pur}_java.vcf" "$WORK_TITR/p${pur}_rust.vcf" > "$KEEP_DIR/p${pur}_compare_rust_java.txt" 2>&1; then
      parity_rust="MATCH"
    fi
    parity_cpp="MISMATCH"
    if python3 harness/compare_vcf.py "$WORK_TITR/p${pur}_java.vcf" "$WORK_TITR/p${pur}_cpp.vcf" > "$KEEP_DIR/p${pur}_compare_cpp_java.txt" 2>&1; then
      parity_cpp="MATCH"
    fi

    # 4. Strelka2 2.9.10
    if [ ! -s "$WORK_TITR/p${pur}_strelka2.vcf.gz" ]; then
      rm -rf "$WORK_TITR/strelka_p${pur}"
      /usr/bin/time -v -o "$KEEP_DIR/p${pur}_strelka2.time" \
        bash -c 'conda run -n cl_strelka configureStrelkaSomaticWorkflow.py --tumorBam "$1" --normalBam "$2" --referenceFasta "$3" --callRegions "$4" --runDir "$5" && conda run -n cl_strelka "$5/runWorkflow.py" -m local -j 16' _ "$ADMIX_BAM" "$N_BAM" "$REF" "$WORK_TITR/eval.bed.gz" "$WORK_TITR/strelka_p${pur}" >/dev/null 2>&1
      bcftools concat -a "$WORK_TITR/strelka_p${pur}/results/variants/somatic.snvs.vcf.gz" "$WORK_TITR/strelka_p${pur}/results/variants/somatic.indels.vcf.gz" \
        | bcftools sort -Oz -o "$WORK_TITR/p${pur}_strelka2.vcf.gz"
      tabix -f -p vcf "$WORK_TITR/p${pur}_strelka2.vcf.gz"
    fi

    # 5. LoFreq 2.1.5
    if [ ! -s "$WORK_TITR/p${pur}_lofreq.vcf.gz" ]; then
      rm -f "$WORK_TITR/lf_p${pur}_"*
      /usr/bin/time -v -o "$KEEP_DIR/p${pur}_lofreq.time" \
        conda run -n cl_lofreq lofreq somatic -n "$N_BAM" -t "$ADMIX_BAM" -f "$REF" -l "$EVAL_BED" -o "$WORK_TITR/lf_p${pur}_" --threads 16 >/dev/null 2>&1 || true
      if [ -f "$WORK_TITR/lf_p${pur}_somatic_final.snvs.vcf.gz" ]; then
        cp -f "$WORK_TITR/lf_p${pur}_somatic_final.snvs.vcf.gz" "$WORK_TITR/p${pur}_lofreq.vcf.gz"
        tabix -f -p vcf "$WORK_TITR/p${pur}_lofreq.vcf.gz"
      fi
    fi

    # 6. Mutect2 4.6.2.0
    if [ ! -s "$WORK_TITR/p${pur}_mutect2.vcf.gz" ]; then
      /usr/bin/time -v -o "$KEEP_DIR/p${pur}_mutect2.time" \
        bash -c 'conda run -n cl_gatk4 gatk Mutect2 -R "$1" -I "$2" -I "$3" -normal HG008-N-D -L "$4" --native-pair-hmm-threads 16 -O "$5" && conda run -n cl_gatk4 gatk FilterMutectCalls -R "$1" -V "$5" -O "$6"' _ "$REF" "$ADMIX_BAM" "$N_BAM" "$EVAL_BED" "$WORK_TITR/m2_p${pur}.unfilt.vcf.gz" "$WORK_TITR/p${pur}_mutect2.vcf.gz" >/dev/null 2>&1
    fi

    for cspec in \
      "VarDictRust|$WORK_TITR/p${pur}_rust.vcf|$KEEP_DIR/p${pur}_rust.time|$parity_rust|somatic_vardict" \
      "VarDictJava|$WORK_TITR/p${pur}_java.vcf|$KEEP_DIR/p${pur}_java.time|ORACLE|somatic_vardict" \
      "vardictcpp|$WORK_TITR/p${pur}_cpp.vcf|$KEEP_DIR/p${pur}_cpp.time|$parity_cpp|somatic_vardict" \
      "Strelka2|$WORK_TITR/p${pur}_strelka2.vcf.gz|$KEEP_DIR/p${pur}_strelka2.time|NA|somatic" \
      "LoFreq|$WORK_TITR/p${pur}_lofreq.vcf.gz|$KEEP_DIR/p${pur}_lofreq.time|NA|somatic" \
      "Mutect2|$WORK_TITR/p${pur}_mutect2.vcf.gz|$KEEP_DIR/p${pur}_mutect2.time|NA|somatic"; do
      IFS="|" read -r cname cvcf ctime cpar cmode <<< "$cspec"
      [ -s "$cvcf" ] || continue
      wall=$(awk -F": " '/Elapsed/{print $2}' "$ctime" 2>/dev/null | tail -1 || true)
      rss=$(awk -F": " '/Maximum resident/{print $2}' "$ctime" 2>/dev/null | tr -dc "0-9" || true)
      metrics=$(eval_single_vcf_against_bed "$cvcf" "$REF" "$SDF" "$TRUTH_V04" "$WORK_TITR/eval_v04.bed" "$cmode" "$WORK_TITR/ev_p${pur}_${cname}" 0)
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$pur" "$vaf_exp" "$cname" "${wall:-NA}" "${rss:-NA}" "$cpar" "$metrics" | tee -a "$OUT_TSV"
    done

    # Persist non-Rust titration VCFs, .time files, and rtg vcfeval summaries to GCS oracle-cache-v183
    if [ -z "${LOCAL_CKPT_DIR:-}" ]; then
      for nrf in "$WORK_TITR/p${pur}_java.vcf" "$WORK_TITR/p${pur}_cpp.vcf" "$WORK_TITR/p${pur}_strelka2.vcf.gz"* "$WORK_TITR/p${pur}_lofreq.vcf.gz"* "$WORK_TITR/p${pur}_mutect2.vcf.gz"* "$KEEP_DIR/p${pur}_java.time" "$KEEP_DIR/p${pur}_cpp.time" "$KEEP_DIR/p${pur}_strelka2.time" "$KEEP_DIR/p${pur}_lofreq.time" "$KEEP_DIR/p${pur}_mutect2.time"; do
        [ -s "$nrf" ] && setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -o GSUtil:parallel_composite_upload_threshold=0 -q cp "$nrf" "${GCS_ORACLE_CACHE}/stage6/titration/" 2>/dev/null || true
      done
      for nr_ev in "$WORK_TITR/ev_p${pur}_VarDictJava" "$WORK_TITR/ev_p${pur}_vardictcpp" "$WORK_TITR/ev_p${pur}_Strelka2" "$WORK_TITR/ev_p${pur}_LoFreq" "$WORK_TITR/ev_p${pur}_Mutect2"; do
        [ -d "$nr_ev" ] && setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -o GSUtil:parallel_composite_upload_threshold=0 -q rsync -r "$nr_ev" "${GCS_ORACLE_CACHE}/stage6/titration/$(basename "$nr_ev")" 2>/dev/null || true
      done
    fi
  done
}

run_s6_vm1() {
  # ---------------------------------------------------------------------------
  # Cell 4: s6_ga4gh_v35_stratification_matrix
  # ---------------------------------------------------------------------------
  run_checkpointed_cell "stage6" "s6_ga4gh_v35_stratification_matrix" 32 "115G" bash -c '
    set -euo pipefail
    source harness/gcp_paper_runner.sh --source-only
    source harness/run_stage_6.sh
    cell_s6_ga4gh_v35_stratification
  '
}

cell_s6_ga4gh_v35_stratification() {
  set -euo pipefail
  ensure_vcf_clean_helper
  local REF_GIAB="$DATA_DIR/ref/giabv3/GRCh38_GIABv3.fasta"
  local SDF_GIAB="$DATA_DIR/ref/giabv3/GRCh38_GIABv3.sdf"
  local XP_DIR="$WORK_DIR/data/cross_platform_chr20"
  local XP_CELL_DIR="$WORK_DIR/cells/s5_cross_platform_short_read_chr20_panel"
  local CELLS_SUB="${GCS_CELLS_SUBDIR:-cells-v106}"

  local TR_HG2="$DATA_DIR/hg002/truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
  local BED_HG2_CHR20="$XP_DIR/hg002_chr20_v421_eval.bed"
  local TR_HG8="$DATA_DIR/hg008/truth/HG008-T_somatic_smvar_benchmark_v0.4_GRCh38_tumorvariants.vcf.gz"
  local BED_HG8_CHR20="$XP_DIR/hg008_chr20_v04_eval.bed"

  local STRAT_DIR="$DATA_DIR/stratifications_v35"
  mkdir -p "$STRAT_DIR" "$XP_DIR" "$XP_CELL_DIR"
  if [ -z "${LOCAL_CKPT_DIR:-}" ]; then
    setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -q rsync -r "${GCS_BAMS}/g4-cache/stratifications_v35" "$STRAT_DIR" 2>/dev/null || true
    [ -s "$BED_HG2_CHR20" ] || setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -q cp "${GCS_BAMS}/cross_platform_chr20/hg002_chr20_v421_eval.bed" "$BED_HG2_CHR20" 2>/dev/null || true
    [ -s "$BED_HG8_CHR20" ] || setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -q cp "${GCS_BAMS}/cross_platform_chr20/hg008_chr20_v04_eval.bed" "$BED_HG8_CHR20" 2>/dev/null || true
  fi
  local BASE_STRAT="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.5/GRCh38@all"

  local sspec sname spath arm_spec arm ref sdf norm_vcf metrics
  for sspec in \
    "Homopolymer_4to6bp|LowComplexity/GRCh38_SimpleRepeat_homopolymer_4to6_slop5.bed.gz" \
    "Homopolymer_7to11bp|LowComplexity/GRCh38_SimpleRepeat_homopolymer_7to11_slop5.bed.gz" \
    "Homopolymer_ge12bp|LowComplexity/GRCh38_SimpleRepeat_homopolymer_ge12_slop5.bed.gz" \
    "AllTandemRepeats|LowComplexity/GRCh38_AllTandemRepeats.bed.gz" \
    "SegmentalDuplications|SegmentalDuplications/GRCh38_segdups.bed.gz" \
    "LowMappability|Mappability/GRCh38_lowmappabilityall.bed.gz" \
    "RefSeq_CDS|Functional/GRCh38_refseq_cds.bed.gz"; do
    IFS="|" read -r sname spath <<< "$sspec"
    if [ ! -s "$STRAT_DIR/${sname}.chr20.bed" ]; then
      curl -sSLf --retry 5 "${BASE_STRAT}/${spath}" | zcat | awk -F"\t" '$1 == "chr20"' > "$STRAT_DIR/${sname}.chr20.bed"
    fi
    [ -s "$STRAT_DIR/hg002_${sname}.bed" ] || bedtools intersect -a "$BED_HG2_CHR20" -b "$STRAT_DIR/${sname}.chr20.bed" > "$STRAT_DIR/hg002_${sname}.bed"
    [ -s "$STRAT_DIR/hg008_${sname}.bed" ] || bedtools intersect -a "$BED_HG8_CHR20" -b "$STRAT_DIR/${sname}.chr20.bed" > "$STRAT_DIR/hg008_${sname}.bed"
  done
  if [ -z "${LOCAL_CKPT_DIR:-}" ]; then
    setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -o GSUtil:parallel_composite_upload_threshold=0 -q rsync -r "$STRAT_DIR" "${GCS_BAMS}/g4-cache/stratifications_v35" 2>/dev/null || true
  fi

  local WORK_ST="$WORK_DIR/tmp_s6_strat"
  rm -rf "$WORK_ST" && mkdir -p "$WORK_ST"
  local OUT_TSV="$KEEP_DIR/ga4gh_v35_cross_platform_stratification.tsv"
  printf "cohort\tarm\tstratum\tsnv_p\tsnv_r\tsnv_f1\tindel_p\tindel_r\tindel_f1\n" > "$OUT_TSV"

  for arm_spec in \
    "hg002_illumina|$REF_GIAB|$SDF_GIAB" \
    "hg002_element_aviti|$REF_GIAB|$SDF_GIAB" \
    "hg002_mgi_dnbseq|$REF_GIAB|$SDF_GIAB" \
    "hg002_ultima_ug100|$REF_GIAB|$SDF_GIAB"; do
    IFS="|" read -r arm ref sdf <<< "$arm_spec"
    if [ ! -s "$XP_CELL_DIR/${arm}_f0.01_rust.vcf" ] && [ -z "${LOCAL_CKPT_DIR:-}" ]; then
      setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -q cp "${GCS_RESULTS}/${CELLS_SUB}/s5_cross_platform_short_read_chr20_panel/${arm}_f0.01_rust.vcf" "$XP_CELL_DIR/${arm}_f0.01_rust.vcf" 2>/dev/null || true
    fi
    norm_vcf="$XP_CELL_DIR/ev_${arm}_f0.01_rust_std/calls.vcf.gz"
    if [ ! -s "$norm_vcf" ]; then
      eval_single_vcf_against_bed "$XP_CELL_DIR/${arm}_f0.01_rust.vcf" "$ref" "$sdf" "$TR_HG2" "$BED_HG2_CHR20" "germline" "$WORK_ST/norm_${arm}" 0 >/dev/null
      norm_vcf="$WORK_ST/norm_${arm}/calls.vcf.gz"
    fi
    for sname in Homopolymer_4to6bp Homopolymer_7to11bp Homopolymer_ge12bp AllTandemRepeats SegmentalDuplications LowMappability RefSeq_CDS; do
      metrics=$(eval_single_vcf_against_bed "$norm_vcf" "$ref" "$sdf" "$TR_HG2" "$STRAT_DIR/hg002_${sname}.bed" "germline" "$WORK_ST/ev_${arm}_${sname}" 1)
      printf "HG002_Germline\t%s\t%s\t%s\n" "$arm" "$sname" "$metrics" | tee -a "$OUT_TSV"
    done
  done

  for arm_spec in \
    "hg008_illumina|$REF_GIAB|$SDF_GIAB" \
    "hg008_element_aviti|$REF_GIAB|$SDF_GIAB" \
    "hg008_pacbio_onso|$REF_GIAB|$SDF_GIAB" \
    "hg008_roche_sbx|$REF_GIAB|$SDF_GIAB" \
    "hg008_ultima_ppmseq|$REF_GIAB|$SDF_GIAB" \
    "hg008_ultima_ug100|$REF_GIAB|$SDF_GIAB"; do
    IFS="|" read -r arm ref sdf <<< "$arm_spec"
    if [ ! -s "$XP_CELL_DIR/${arm}_f0.01_rust.vcf" ] && [ -z "${LOCAL_CKPT_DIR:-}" ]; then
      setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -q cp "${GCS_RESULTS}/${CELLS_SUB}/s5_cross_platform_short_read_chr20_panel/${arm}_f0.01_rust.vcf" "$XP_CELL_DIR/${arm}_f0.01_rust.vcf" 2>/dev/null || true
    fi
    norm_vcf="$XP_CELL_DIR/ev_${arm}_f0.01_rust_std/calls.vcf.gz"
    if [ ! -s "$norm_vcf" ]; then
      eval_single_vcf_against_bed "$XP_CELL_DIR/${arm}_f0.01_rust.vcf" "$ref" "$sdf" "$TR_HG8" "$BED_HG8_CHR20" "somatic_vardict" "$WORK_ST/norm_${arm}" 0 >/dev/null
      norm_vcf="$WORK_ST/norm_${arm}/calls.vcf.gz"
    fi
    for sname in Homopolymer_4to6bp Homopolymer_7to11bp Homopolymer_ge12bp AllTandemRepeats SegmentalDuplications LowMappability RefSeq_CDS; do
      metrics=$(eval_single_vcf_against_bed "$norm_vcf" "$ref" "$sdf" "$TR_HG8" "$STRAT_DIR/hg008_${sname}.bed" "somatic" "$WORK_ST/ev_${arm}_${sname}" 1)
      printf "HG008_Somatic\t%s\t%s\t%s\n" "$arm" "$sname" "$metrics" | tee -a "$OUT_TSV"
    done
  done
}

run_stage_6() {
  echo "=== [Stage 6: Literature-Gap Closure Suite (HG002 v5.0q/CMRG/Mosaic, HG008 Modern ML/Commercial Callers, Low-VAF Titration, GA4GH v3.5 Stratification)] ==="
  unset SOCKET_PIN_CPUS CELL_ALLOWLIST
  export DATA_DIR="scratch/g4"
  run_s6_vm2
  run_s6_vm1
  gcs_mark_done "STAGE6_DONE"
  gcs_write_heartbeat "stage6" "STAGE6_DONE" "COMPLETED"
}
