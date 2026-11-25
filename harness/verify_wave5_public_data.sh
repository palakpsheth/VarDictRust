#!/usr/bin/env bash
# =============================================================================
# harness/verify_wave5_public_data.sh
# Comprehensive 7-gate real-world clinical & validator stress-test suite for
# Wave-5 --csi spanning-SV indexing and CompatProfile switches 22..32 (bringing
# --compat=off to 32 switches total).
#
# Cross-Ecosystem Provenance & Citation Registry (owner/repo#N + DOIs):
# - Gate 1 (Ungated --csi LinearIndex loffset fix for spanning SVs):
#   clintval/vartovcf#161, samtools/hts-specs CSIv1 specification.
# - Gate 2 (Switch 22: cigar_strip_outer_hard_clips):
#   AstraZeneca-NGS/VarDictJava#341, AstraZeneca-NGS/VarDictJava#103,
#   Kennedy-Lab-UW/VarDictJava@9a0beb0d, hedera-dx/VarDictJava@17d1cb51.
# - Gate 3 (Switches 23..25: nonnegative_subtracted_means,
#   reconcile_strand_and_hq_sums, reconcile_rd_vd_dp):
#   bcbio/bcbio-nextgen#1704, AstraZeneca-NGS/VarDictJava#316,
#   AstraZeneca-NGS/VarDictJava#323, AstraZeneca-NGS/VarDict#168,
#   clintval/vartovcf#102, clintval/vartovcf#120, clintval/vartovcf#142.
# - Gate 4 (Switch 26: sanitize_iupac_alt_alleles):
#   AstraZeneca-NGS/VarDictJava#267, AstraZeneca-NGS/VarDictJava#287,
#   AstraZeneca-NGS/VarDict#137, bcbio/bcbio-nextgen#3321, bioinform/somaticseq.
# - Gate 5 (Switch 27: sv_respect_min_len):
#   AstraZeneca-NGS/VarDictJava#328, AstraZeneca-NGS/VarDict#157.
# - Gate 6 (Switch 28: pileup_preserve_shifted_ref):
#   AstraZeneca-NGS/VarDictJava#330.
# - Gate 7 (Switches 29..32: recompute_vartype_after_adj_complex,
#   strict_vcf_ad_number_r, omit_nonsv_vcf_end, canonical_vcf_gt):
#   clintval/vartovcf#54, clintval/vartovcf#99, clintval/vartovcf#101,
#   clintval/vartovcf#105, clintval/vartovcf#109, clintval/vartovcf#147,
#   AstraZeneca-NGS/VarDict#75, AstraZeneca-NGS/VarDict#76,
#   AstraZeneca-NGS/VarDictJava#344, AstraZeneca-NGS/VarDict#186,
#   EBIvariation/vcf-validator#99, sigven/pcgr#12;
#   Cleary et al., bioRxiv 2015 (DOI: 10.1101/023754, rtg vcfeval);
#   Fang et al., Nat. Biotechnol. 2021 (DOI: 10.1038/s41587-021-00993-6, SEQC2).
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
source "$SCRIPT_DIR/lib/portable.sh"

ensure_g4_tools_path
for _conda_env in "${CONDA_ROOT:-${HOME:-}/miniconda3}/envs/g4-tools/bin" \
                  "${CONDA_ROOT:-${HOME:-}/miniconda3}/envs/vardict-oracle/bin"; do
    if [ -d "$_conda_env" ]; then
        case ":$PATH:" in
            *":$_conda_env:"*) ;;
            *) export PATH="$_conda_env:$PATH" ;;
        esac
    fi
done
unset _conda_env

for _req_tool in samtools bcftools tabix bgzip rtg python3 gzip; do
    if ! command -v "$_req_tool" >/dev/null 2>&1; then
        echo "ERROR: required tool '$_req_tool' not found on PATH ($PATH)" >&2
        exit 1
    fi
done
unset _req_tool

VARDICTRUST_BIN="${VARDICTRUST_BIN:-${VARDICT_BIN:-$ROOT_DIR/target/release/vardictrust}}"
if [ ! -x "$VARDICTRUST_BIN" ]; then
    echo "ERROR: VARDICTRUST_BIN is not executable at '$VARDICTRUST_BIN' (run 'cargo build --release' first)" >&2
    exit 1
fi

DATA_ROOT="$(resolve_data_root "$ROOT_DIR")"
G4_DIR="${DATA_DIR:-$DATA_ROOT/scratch/g4}"
WAVE4_DIR="${WAVE4_DIR:-$G4_DIR/wave4}"

if [ ! -f "$WAVE4_DIR/.done" ]; then
    DATA_DIR="$G4_DIR" WAVE4_DIR="$WAVE4_DIR" bash "$SCRIPT_DIR/get_wave4_public_slices.sh"
fi

FA_HG19="$DATA_ROOT/VarDictJava/testdata/fastas/hg19.fa"
FA_B37="$DATA_ROOT/VarDictJava/testdata/fastas/GRCh37.fa"
FA_HS38DH="$DATA_ROOT/VarDictJava/testdata/fastas/hs38DH.fa"
FA_HG1K="$DATA_ROOT/VarDictJava/testdata/fastas/hg1k.fa"
if [ ! -f "$FA_HG19" ] || [ ! -f "$FA_B37" ] || [ ! -f "$FA_HS38DH" ] || [ ! -f "$FA_HG1K" ]; then
    python3 "$ROOT_DIR/harness/convert_fastas.py"
fi

GRCH38_REF="$G4_DIR/ref/d1vd1/GRCh38.d1.vd1.fa"
GIAB_REF="$G4_DIR/ref/giabv3/GRCh38_GIABv3.fa"
if [ ! -f "$GIAB_REF" ] && [ -f "$G4_DIR/ref/giabv3/GRCh38_GIABv3.fasta" ]; then
    GIAB_REF="$G4_DIR/ref/giabv3/GRCh38_GIABv3.fasta"
fi

VJ_BAMS="$DATA_ROOT/VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests"

for _req_file in "$GRCH38_REF" "${GRCH38_REF}.fai" "$GIAB_REF" "${GIAB_REF}.fai" \
                 "$G4_DIR/hg002/HG002.chr20.bam" \
                 "$G4_DIR/hg008/HG008-T.slice.bam" "$G4_DIR/hg008/HG008-N.slice.bam" \
                 "$G4_DIR/seqc2/HCC1395.bam" "$G4_DIR/seqc2/HCC1395BL.bam" \
                 "$G4_DIR/seqc2/truth/truth.vcf.gz"; do
    if [ ! -s "$_req_file" ]; then
        echo "ERROR: missing Wave-5 clinical dataset prerequisite: '$_req_file'" >&2
        exit 1
    fi
done
unset _req_file

WORK_TMP="$(mktemp -d "${TMPDIR:-/tmp}/vdr_wave5_verify.XXXXXX")"
trap 'rm -rf "$WORK_TMP"' EXIT

T0_MS="$(now_ms)"
echo "================================================================================"
echo " Wave-5 7-Gate Clinical & Validator Stress-Test Suite (32 CompatProfile Switches)"
echo " Binary: $VARDICTRUST_BIN"
echo " Data:   $G4_DIR"
echo "================================================================================"

# =============================================================================
# Gate 1: --csi Multi-BGZF-Block Spanning-SV tabix & bcftools view -r Gate
#         (clintval/vartovcf#161)
# =============================================================================
echo "[1/7] Verifying multi-BGZF-block --csi spanning-SV queries via tabix & bcftools view -r..."
python3 - "$VARDICTRUST_BIN" "$FA_HG19" "$GIAB_REF" "$VJ_BAMS" "$G4_DIR" "$WORK_TMP" <<'PY'
import gzip, os, struct, subprocess, sys

vdr, fa_hg19, giab_ref, vj_bams, g4_dir, tmp = sys.argv[1:7]

def count_bgzf_blocks(gz_path):
    """Count non-EOF BGZF blocks via RFC1952 extra subfield 'BC'."""
    blocks = 0
    size = os.path.getsize(gz_path)
    with open(gz_path, "rb") as fh:
        pos = 0
        while pos < size:
            hdr = fh.read(12)
            if not hdr:
                break
            if len(hdr) < 12 or hdr[0:2] != b"\x1f\x8b" or (hdr[3] & 4) == 0:
                raise AssertionError(f"{gz_path}: invalid BGZF block header at offset {pos}")
            xlen = struct.unpack("<H", hdr[10:12])[0]
            extra = fh.read(xlen)
            bsize = None
            idx = 0
            while idx + 4 <= len(extra):
                si1, si2, slen = extra[idx], extra[idx + 1], struct.unpack("<H", extra[idx + 2 : idx + 4])[0]
                if si1 == 66 and si2 == 67 and slen == 2:
                    bsize = struct.unpack("<H", extra[idx + 4 : idx + 6])[0] + 1
                    break
                idx += 4 + slen
            if bsize is None:
                raise AssertionError(f"{gz_path}: missing BC subfield at offset {pos}")
            if bsize > 28:
                blocks += 1
            pos += bsize
            fh.seek(pos)
    return blocks

# 1a) Colo829 DELDiscReverse with -p across chr10:89700380-89713976:
#     Crosses the 16,384-bp CSI leaf-bin boundary at 5475 * 16384 = 89,702,400
#     while emitting ~13,600 pileup rows across >= 8 BGZF blocks!
lgdel_gz = os.path.join(tmp, "colo829_deldisc_pileup.vcf.gz")
subprocess.run(
    [vdr, "-G", fa_hg19, "-b", os.path.join(vj_bams, "Colo829.chr10.DELDiscReverse.bam"),
     "-R", "chr10:89700380-89713976", "-f", "0.001", "-p", "--compat=off", "--vcfall",
     "--output", lgdel_gz, "--csi"],
    check=True,
)
lgdel_blocks = count_bgzf_blocks(lgdel_gz)
if lgdel_blocks < 8:
    raise AssertionError(f"Expected >= 8 BGZF blocks in {lgdel_gz}, got {lgdel_blocks}")

del_svs = []
with gzip.open(lgdel_gz, "rt", encoding="utf-8") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if p[4] == "<DEL>":
            kv = dict(x.split("=", 1) for x in p[7].split(";") if "=" in x)
            del_svs.append((p[0], int(p[1]), int(kv["END"]), p[4]))
if not del_svs:
    raise AssertionError("Expected <DEL> in Colo829 DELDiscReverse pileup --csi VCF")

# 1b) Multi-BGZF-block spanning <DEL>, <DUP>, and <INV> BAM over chr20:46840000-46910000
#     (each SV spans > 36 kbp across >= 3 distinct 16,384-bp CSI leaf bins:
#      2859*16384=46841856, 2860*16384=46858240, 2861*16384=46874624, 2862*16384=46891008)
#     paired with dense read coverage across 46841500..46896000 so -p writes >= 15 BGZF blocks!
fa_chr20 = subprocess.check_output(
    ["samtools", "faidx", fa_hg19, "chr20:46839001-46915000"], text=True
)
chr20_seq = "".join(l.strip() for l in fa_chr20.splitlines() if not l.startswith(">")).upper()
chr20_origin = 46839001  # 1-based

def ref_slice(p1_start, length):
    off = p1_start - chr20_origin
    return chr20_seq[off : off + length]

def revcomp(s):
    tr = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return s.translate(tr)[::-1]

sam_path = os.path.join(tmp, "span3_sv.sam")
bam_path = os.path.join(tmp, "span3_sv.bam")
chr20_len = None
with open(fa_hg19 + ".fai", "r", encoding="utf-8") as fh:
    for line in fh:
        p = line.split("\t")
        if p[0] == "chr20":
            chr20_len = int(p[1])
            break
if not chr20_len:
    raise AssertionError("chr20 missing from hg19.fa.fai")

records = []
# Dense wild-type tiling across 46841000..46896000 so -p generates ~55,000 rows (> 40 BGZF blocks)
r_id = 0
for pos in range(46841000, 46896000, 35):
    for strand_flag in (0, 16):
        r_id += 1
        seq = ref_slice(pos, 100)
        records.append((pos, f"wt_{r_id}\t{strand_flag}\tchr20\t{pos}\t60\t100M\t*\t0\t0\t{seq}\t{'I'*100}"))

# Spanning <DEL>: 46841600..46885000 (~43.4 kbp across 3 CSI leaf-bin boundaries) via 3' + 5' split reads + discordant pairs
del_l, del_r = 46841600, 46885000
for i in range(10):
    s1 = del_l - 60
    seq1 = ref_slice(s1, 60) + ref_slice(del_r, 30)
    records.append((s1, f"sdel_f_{i}\t99\tchr20\t{s1}\t60\t60M30S\t=\t{del_r}\t{del_r + 60 - s1}\t{seq1}\t{'I'*90}"))
    s2 = del_r
    seq2 = ref_slice(del_l - 30, 30) + ref_slice(s2, 60)
    records.append((s2, f"sdel_f_{i}\t147\tchr20\t{s2}\t60\t30S60M\t=\t{s1}\t{-(del_r + 60 - s1)}\t{seq2}\t{'I'*90}"))

# Spanning <DUP>: 46842200..46887500 (~45.3 kbp across 3 CSI leaf-bin boundaries)
dup_l, dup_r = 46842200, 46887500
for i in range(10):
    # Reverse read at dup_l with 5' soft-clip matching dup_r - 29 .. dup_r
    seq_r = ref_slice(dup_r - 29, 30) + ref_slice(dup_l, 70)
    records.append((dup_l, f"sdup_{i}\t81\tchr20\t{dup_l}\t60\t30S70M\t=\t{dup_r - 69}\t{dup_r - dup_l + 1}\t{seq_r}\t{'I'*100}"))
    # Forward read ending at dup_r with 3' soft-clip matching dup_l .. dup_l + 29
    sf = dup_r - 69
    seq_f = ref_slice(sf, 70) + ref_slice(dup_l, 30)
    records.append((sf, f"sdup_{i}\t161\tchr20\t{sf}\t60\t70M30S\t=\t{dup_l}\t{-(dup_r - dup_l + 1)}\t{seq_f}\t{'I'*100}"))

# Spanning <INV>: 46843000..46889500 (~46.5 kbp across 3 CSI leaf-bin boundaries, Side::Three 70M30S at inv_r + forward mate at inv_l)
inv_l, inv_r = 46843000, 46889500
for i in range(10):
    sm = inv_l
    seq_m = ref_slice(sm, 80)
    sf = inv_r - 69
    seq_f = ref_slice(sf, 70) + revcomp(ref_slice(inv_l - 19, 20))
    records.append((sm, f"sinv_m_{i}\t1\tchr20\t{sm}\t60\t80M\t=\t{sf}\t{sf + 70 - sm}\t{seq_m}\t{'I'*80}"))
    records.append((sf, f"sinv_s_{i}\t1\tchr20\t{sf}\t60\t70M20S\t=\t{sm}\t{-(sf + 70 - sm)}\t{seq_f}\t{'I'*90}"))

records.sort(key=lambda x: x[0])
with open(sam_path, "w", encoding="utf-8") as out:
    out.write(f"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr20\tLN:{chr20_len}\n")
    for _, line in records:
        out.write(line + "\n")
subprocess.run(["samtools", "view", "-b", "-o", bam_path, sam_path], check=True)
subprocess.run(["samtools", "index", bam_path], check=True)

span3_gz = os.path.join(tmp, "span3_sv_pileup.vcf.gz")
subprocess.run(
    [vdr, "-G", fa_hg19, "-b", bam_path, "-R", "chr20:46841200-46895000",
     "-f", "0.001", "-p", "--compat=off", "--vcfall", "--output", span3_gz, "--csi", "-th", "4"],
    check=True,
)
span3_blocks = count_bgzf_blocks(span3_gz)
if span3_blocks < 20:
    raise AssertionError(f"Expected >= 20 BGZF blocks in {span3_gz}, got {span3_blocks}")

found_types = {}
with gzip.open(span3_gz, "rt", encoding="utf-8") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if p[4] in ("<DEL>", "<DUP>", "<INV>"):
            kv = dict(x.split("=", 1) for x in p[7].split(";") if "=" in x)
            s_pos, e_pos = int(p[1]), int(kv["END"])
            if e_pos - s_pos >= 35000:
                found_types[p[4]] = (p[0], s_pos, e_pos, p[4])

for req_sv in ("<DEL>", "<DUP>", "<INV>"):
    if req_sv not in found_types:
        raise AssertionError(f"Missing >= 35 kbp spanning {req_sv} in {span3_gz}: got {found_types}")

# Query EVERY spanning SV (Colo829 <DEL> + chr20 <DEL>/<DUP>/<INV>) at interior leaf-bin coordinates
all_checks = [(lgdel_gz, del_svs[0])] + [(span3_gz, found_types[k]) for k in ("<DEL>", "<DUP>", "<INV>")]
queried_points = 0
for vcf_gz, (chrom, s_pos, e_pos, sv_alt) in all_checks:
    start_bin = (s_pos - 1) >> 14
    end_bin = (e_pos - 1) >> 14
    if end_bin <= start_bin:
        raise AssertionError(f"{sv_alt} ({chrom}:{s_pos}-{e_pos}) did not cross a 16,384-bp CSI leaf-bin boundary")
    for b in range(start_bin + 1, end_bin + 1):
        q_coord = min(e_pos - 1, (b << 14) + 500)
        if q_coord <= s_pos:
            continue
        reg = f"{chrom}:{q_coord}-{q_coord}"
        tbx_out = subprocess.check_output(["tabix", vcf_gz, reg], text=True)
        bcf_out = subprocess.check_output(["bcftools", "view", "-H", "-r", reg, vcf_gz], text=True)
        for tool_name, txt in (("tabix", tbx_out), ("bcftools view -r", bcf_out)):
            if not any(
                l.split("\t")[1] == str(s_pos) and l.split("\t")[4] == sv_alt
                for l in txt.splitlines() if l.strip()
            ):
                raise AssertionError(
                    f"{tool_name} dropped spanning {sv_alt} ({chrom}:{s_pos}-{e_pos}) at interior leaf-bin query {reg}"
                )
        queried_points += 1

# 1c) GIAB HG008 multi-BGZF-block somatic VCF (--csi) across chr1 slices
hg008_gz = os.path.join(tmp, "hg008_multi_csi.vcf.gz")
subprocess.run(
    [vdr, "-G", giab_ref,
     "-b", f"{g4_dir}/hg008/HG008-T.slice.bam|{g4_dir}/hg008/HG008-N.slice.bam",
     "-c", "1", "-S", "2", "-E", "3", "-f", "0.005", "--fisher", "--compat=off", "--vcfall",
     "--output", hg008_gz, "--csi", "-th", "4", os.path.join(g4_dir, "hg008", "slices.bed")],
    check=True,
)
hg008_blocks = count_bgzf_blocks(hg008_gz)
if hg008_blocks < 2:
    raise AssertionError(f"Expected >= 2 BGZF blocks on GIAB HG008 slices, got {hg008_blocks}")
tbx_hg8 = [tuple(l.split("\t")[:5]) for l in subprocess.check_output(["tabix", hg008_gz, "chr1:72271446-72303047"], text=True).strip().splitlines() if l.strip()]
bcf_hg8 = [tuple(l.split("\t")[:5]) for l in subprocess.check_output(["bcftools", "view", "-H", "-r", "chr1:72271446-72303047", hg008_gz], text=True).strip().splitlines() if l.strip()]
if not tbx_hg8 or tbx_hg8 != bcf_hg8:
    raise AssertionError(f"GIAB HG008 --csi tabix vs bcftools mismatch ({len(tbx_hg8)} vs {len(bcf_hg8)})")

print(
    f"  PASS Gate 1: Colo829 LGDEL3 ({lgdel_blocks} BGZF blocks), chr20 <DEL>/<DUP>/<INV> ({span3_blocks} BGZF blocks, "
    f"{queried_points} interior leaf-bin queries, 0 dropped SVs), & GIAB HG008 ({hg008_blocks} BGZF blocks)"
)
PY

# =============================================================================
# Gate 2: SEQ-Matched Hard-Clip Clinical BAM Twins (90M3D15M3I vs 90M3D15M3I20H)
#         (AstraZeneca-NGS/VarDictJava#341, AstraZeneca-NGS/VarDictJava#103,
#          Kennedy-Lab-UW/VarDictJava@9a0beb0d)
# =============================================================================
echo "[2/7] Verifying SEQ-matched hard-clipped clinical BAM twins on SEQC2 HCC1395 & QIAseq..."
python3 - "$VARDICTRUST_BIN" "$GRCH38_REF" "$G4_DIR" "$WAVE4_DIR" "$WORK_TMP" <<'PY'
import os, re, subprocess, sys

vdr, grch38_ref, g4_dir, w4_dir, tmp = sys.argv[1:6]
cigar_re = re.compile(r"(\d+)([MIDNSHP=X])")

def strip_outer_hp(cigar):
    ops = [(int(n), op) for n, op in cigar_re.findall(cigar)]
    while ops and ops[0][1] in ("H", "P"):
        ops.pop(0)
    while ops and ops[-1][1] in ("H", "P"):
        ops.pop()
    ops = [(n, op) for n, op in ops if op != "P"]
    return "".join(f"{n}{op}" for n, op in ops)

def build_clinical_twin_bams(tag, src_bam, region_str, anchor_pos):
    fa_txt = subprocess.check_output(
        ["samtools", "faidx", grch38_ref, f"{region_str.split(':')[0]}:{anchor_pos}-{anchor_pos + 149}"],
        text=True,
    )
    ref_seg = "".join(l.strip() for l in fa_txt.splitlines() if not l.startswith(">")).upper()
    if len(ref_seg) < 120:
        raise AssertionError(f"{tag}: short ref segment at {anchor_pos}")

    sam_raw = subprocess.check_output(["samtools", "view", "-h", src_bam, region_str], text=True)
    hdr_lines, aln_clean, aln_hclip = [], [], []
    for line in sam_raw.splitlines():
        if line.startswith("@"):
            hdr_lines.append(line)
            continue
        cols = line.split("\t")
        if len(cols) < 11 or cols[5] == "*":
            continue
        clean_cig = strip_outer_hp(cols[5])
        if not clean_cig:
            continue
        c_clean = list(cols)
        c_clean[5] = clean_cig
        c_hclip = list(cols)
        c_hclip[5] = f"15H{clean_cig}20H"
        pos_i = int(cols[3])
        aln_clean.append((pos_i, "\t".join(c_clean)))
        aln_hclip.append((pos_i, "\t".join(c_hclip)))

    # Inject 60 SEQ-matched clinical template reads with exact 90M3D15M3I vs 90M3D15M3I20H
    # (SEQ length = 90 + 15 + 3 = 108 bases; trailing 3I adjacent to 20H triggers
    # AstraZeneca-NGS/VarDictJava#341 & AstraZeneca-NGS/VarDictJava#103!)
    chrom = region_str.split(":")[0]
    ins_bases = "ACG" if ref_seg[108:111] != "ACG" else "TGC"
    twin_seq = ref_seg[0:90] + ref_seg[93:108] + ins_bases
    if len(twin_seq) != 108:
        raise AssertionError(f"Expected 108-bp SEQ, got {len(twin_seq)}")
    twin_qual = "I" * 108
    for idx in range(60):
        flag = 16 if (idx % 2 == 1) else 0
        rname = f"seqmatched_90M3D15M3I_{tag}_{idx}"
        aln_clean.append(
            (anchor_pos, f"{rname}\t{flag}\t{chrom}\t{anchor_pos}\t60\t90M3D15M3I\t*\t0\t0\t{twin_seq}\t{twin_qual}")
        )
        cig_hc = "90M3D15M3I20H" if idx < 30 else "12H90M3D15M3I20H"
        aln_hclip.append(
            (anchor_pos, f"{rname}\t{flag}\t{chrom}\t{anchor_pos}\t60\t{cig_hc}\t*\t0\t0\t{twin_seq}\t{twin_qual}")
        )

    aln_clean.sort(key=lambda x: x[0])
    aln_hclip.sort(key=lambda x: x[0])
    clean_bam = os.path.join(tmp, f"{tag}_clean.bam")
    hclip_bam = os.path.join(tmp, f"{tag}_hclip.bam")
    for out_bam, recs in ((clean_bam, aln_clean), (hclip_bam, aln_hclip)):
        sam_p = out_bam + ".sam"
        with open(sam_p, "w", encoding="utf-8") as fh:
            for h in hdr_lines:
                fh.write(h + "\n")
            for _, r in recs:
                fh.write(r + "\n")
        subprocess.run(["samtools", "view", "-b", "-o", out_bam, sam_p], check=True)
        subprocess.run(["samtools", "index", out_bam], check=True)
    return clean_bam, hclip_bam, len(aln_clean)

cases = [
    ("seqc2_hcc1395", os.path.join(g4_dir, "seqc2", "HCC1395.bam"), "chr17:74747665-74755000", 74750100, []),
    ("qiaseq_spe", os.path.join(w4_dir, "bam", "SampleA_QGN1.bam"), "chr1:11144542-11145012", 11144620, []),
]

for tag, src_bam, reg, anchor, extra in cases:
    clean_bam, hclip_bam, n_reads = build_clinical_twin_bams(tag, src_bam, reg, anchor)
    if n_reads < 150:
        raise AssertionError(f"{tag}: expected >= 150 real clinical reads + injected twins in {reg}, got {n_reads}")
    vj_clean = subprocess.check_output(
        [vdr, "-G", grch38_ref, "-b", clean_bam, "-R", reg, "-N", tag, "-f", "0.005", "--compat=vardictjava", "-v", *extra],
        text=True,
    )
    vj_hclip = subprocess.check_output(
        [vdr, "-G", grch38_ref, "-b", hclip_bam, "-R", reg, "-N", tag, "-f", "0.005", "--compat=vardictjava", "-v", *extra],
        text=True,
    )
    if vj_clean == vj_hclip:
        raise AssertionError(f"{tag}: expected default --compat=vardictjava to diverge on 90M3D15M3I20H")

    off_var_clean = subprocess.check_output(
        [vdr, "-G", grch38_ref, "-b", clean_bam, "-R", reg, "-N", tag, "-f", "0.005", "--compat=off", "--vcfall", "--emit-intermediate", "--fisher", *extra],
        text=True,
    )
    off_var_hclip = subprocess.check_output(
        [vdr, "-G", grch38_ref, "-b", hclip_bam, "-R", reg, "-N", tag, "-f", "0.005", "--compat=off", "--vcfall", "--emit-intermediate", "--fisher", *extra],
        text=True,
    )
    off_vcf_clean = subprocess.check_output(
        [vdr, "-G", grch38_ref, "-b", clean_bam, "-R", reg, "-N", tag, "-f", "0.005", "--compat=off", "--vcfall", "-v", "--fisher", *extra],
        text=True,
    )
    off_vcf_hclip = subprocess.check_output(
        [vdr, "-G", grch38_ref, "-b", hclip_bam, "-R", reg, "-N", tag, "-f", "0.005", "--compat=off", "--vcfall", "-v", "--fisher", *extra],
        text=True,
    )
    if off_var_clean != off_var_hclip:
        raise AssertionError(f"{tag}: --compat=off .var (including PMEAN) diverged between 90M3D15M3I and 90M3D15M3I20H")
    if off_vcf_clean != off_vcf_hclip:
        raise AssertionError(f"{tag}: --compat=off .vcf (including PMEAN) diverged between 90M3D15M3I and 90M3D15M3I20H")
    n_vars = sum(1 for l in off_vcf_clean.splitlines() if l and not l.startswith("#"))
    if n_vars < 2:
        raise AssertionError(f"{tag}: expected >= 2 called variants in twin slice, got {n_vars}")
    print(f"  PASS Gate 2 ({tag}): {n_reads} reads, {n_vars} variants — byte-identical --compat=off .var & .vcf + PMEAN")
PY

# =============================================================================
# Gate 3: 6-Cohort Accounting & Non-Negative QUAL Invariants
#         (AstraZeneca-NGS/VarDictJava#316, AstraZeneca-NGS/VarDict#168,
#          bcbio/bcbio-nextgen#1704)
# =============================================================================
echo "[3/7] Scanning 100% of --compat=off .var & .vcf rows across all 6 clinical cohorts..."
python3 - "$VARDICTRUST_BIN" "$GRCH38_REF" "$GIAB_REF" "$G4_DIR" "$WAVE4_DIR" "$WORK_TMP" <<'PY'
import os, subprocess, sys

vdr, grch38_ref, giab_ref, g4_dir, w4_dir, tmp = sys.argv[1:7]

seqc2_sub_bed = os.path.join(tmp, "seqc2_gate3.bed")
with open(seqc2_sub_bed, "w", encoding="utf-8") as fh:
    fh.write("chr11\t6400000\t6459439\tslice_1\nchr17\t74747665\t74952140\tslice_3\n")

hg002_chr20_bed = os.path.join(tmp, "hg002_chr20_gate3.bed")
with open(hg002_chr20_bed, "w", encoding="utf-8") as fh:
    fh.write("chr20\t10000000\t10060000\thg002_c20\n")

qgn_bed = os.path.join(w4_dir, "beds", "qiaseq_spe_grch38.amp8.bed")
ilm_bed = os.path.join(w4_dir, "beds", "ilm_twoprimer_grch38.amp8.bed")
tfs_bed = os.path.join(w4_dir, "beds", "tfs_iontorrent_grch38.amp8.bed")

cohorts = [
    # 1. GIAB HG002 (CMRG/Mosaic + chr20 WGS)
    ("GIAB_HG002_CMRG", "simple", giab_ref, os.path.join(w4_dir, "hg002_cmrg_mosaic", "HG002.cmrg_mosaic.slice.bam"),
     ["-c", "1", "-S", "2", "-E", "3", "-g", "4", "-f", "0.01", "-th", "4", os.path.join(w4_dir, "hg002_cmrg_mosaic", "hg002_cmrg_mosaic.bed")]),
    ("GIAB_HG002_chr20", "simple", giab_ref, os.path.join(g4_dir, "hg002", "HG002.chr20.bam"),
     ["-c", "1", "-S", "2", "-E", "3", "-g", "4", "-f", "0.01", "-th", "4", hg002_chr20_bed]),
    # 2. GIAB HG008 T/N
    ("GIAB_HG008_TN", "somatic", giab_ref, f"{g4_dir}/hg008/HG008-T.slice.bam|{g4_dir}/hg008/HG008-N.slice.bam",
     ["-c", "1", "-S", "2", "-E", "3", "-f", "0.01", "-th", "4", os.path.join(g4_dir, "hg008", "test5.bed")]),
    # 3. FDA SEQC2 HCC1395 T/N
    ("SEQC2_HCC1395_TN", "somatic", grch38_ref, f"{g4_dir}/seqc2/HCC1395.bam|{g4_dir}/seqc2/HCC1395BL.bam",
     ["-c", "1", "-S", "2", "-E", "3", "-g", "4", "-f", "0.01", "-th", "4", seqc2_sub_bed]),
    # 4. QIAseq SPE (single & 0.5% ctDNA SampleAC5 paired)
    ("QIAseq_SPE_Single", "amplicon", grch38_ref, os.path.join(w4_dir, "bam", "SampleA_QGN1.bam"),
     ["-a", "10:0.95", "-sp", "-f", "0.005", "-th", "4", qgn_bed]),
    ("QIAseq_ctDNA_Paired", "somatic", grch38_ref, f"{w4_dir}/bam/SampleAC5_QGN1.bam|{w4_dir}/bam/SampleB_QGN1.bam",
     ["-a", "10:0.95", "-sp", "-f", "0.0025", "-th", "4", qgn_bed]),
    # 5. Illumina PE Amplicon (single & paired)
    ("Illumina_PE_Single", "amplicon", grch38_ref, os.path.join(w4_dir, "bam", "SampleA_ILM1.bam"),
     ["-a", "10:0.95", "-f", "0.01", "-th", "4", ilm_bed]),
    ("Illumina_PE_Paired", "somatic", grch38_ref, f"{w4_dir}/bam/SampleA_ILM1.bam|{w4_dir}/bam/SampleB_ILM1.bam",
     ["-a", "10:0.95", "-f", "0.01", "-th", "4", ilm_bed]),
    # 6. Ion Torrent Amplicon (single & paired)
    ("IonTorrent_Single", "amplicon", grch38_ref, os.path.join(w4_dir, "bam", "SampleA_TFS1.bam"),
     ["-a", "10:0.95", "-f", "0.01", "-th", "4", tfs_bed]),
    ("IonTorrent_Paired", "somatic", grch38_ref, f"{w4_dir}/bam/SampleA_TFS1.bam|{w4_dir}/bam/SampleB_TFS1.bam",
     ["-a", "10:0.95", "-f", "0.01", "-th", "4", tfs_bed]),
]

total_var_rows = 0
total_vcf_rows = 0

for name, kind, ref_fa, bam_arg, extra in cohorts:
    var_txt = subprocess.check_output(
        [vdr, "-G", ref_fa, "-b", bam_arg, "--compat=off", "--vcfall", "--emit-intermediate", "--fisher", *extra],
        text=True,
    )
    vcf_txt = subprocess.check_output(
        [vdr, "-G", ref_fa, "-b", bam_arg, "--compat=off", "--vcfall", "-v", "--fisher", *extra],
        text=True,
    )
    v_lines = [l for l in var_txt.splitlines() if l.strip()]
    c_lines = [l for l in vcf_txt.splitlines() if l.strip() and not l.startswith("#")]
    if not v_lines or not c_lines:
        raise AssertionError(f"{name}: expected non-empty .var ({len(v_lines)}) and .vcf ({len(c_lines)})")

    # Audit 100% of .var rows
    for row in v_lines:
        cols = row.split("\t")
        if kind in ("simple", "amplicon"):
            exp_len = 38 if kind == "simple" else 40
            if len(cols) != exp_len:
                raise AssertionError(f"{name}: expected {exp_len} .var cols, got {len(cols)}")
            dp, vd = int(cols[7]), int(cols[8])
            rf, rr, vf, vr = int(cols[9]), int(cols[10]), int(cols[11]), int(cols[12])
            pmean, qmean, mq, nm = float(cols[16]), float(cols[18]), float(cols[22]), float(cols[29])
            hicnt = int(cols[30])
            if dp < rf + rr + vd:
                raise AssertionError(f"{name} .var: DP ({dp}) < RD ({rf}+{rr}) + VD ({vd}): {row}")
            if vf + vr != vd:
                raise AssertionError(f"{name} .var: Fwd+Rev ({vf}+{vr}) != VD ({vd}): {row}")
            if not (0 <= hicnt <= vd):
                raise AssertionError(f"{name} .var: HICNT ({hicnt}) outside [0, VD={vd}]: {row}")
            if min(pmean, qmean, mq, nm) < 0.0:
                raise AssertionError(f"{name} .var: negative mean metric ({pmean}, {qmean}, {mq}, {nm}): {row}")
        else:
            if len(cols) != 61:
                raise AssertionError(f"{name}: expected 61 somatic .var cols, got {len(cols)}")
            for blk_name, base in (("tumor", 7), ("normal", 27)):
                dp, vd = int(cols[base]), int(cols[base + 1])
                rf, rr, vf, vr = int(cols[base + 2]), int(cols[base + 3]), int(cols[base + 4]), int(cols[base + 5])
                pmean, qmean = float(cols[base + 9]), float(cols[base + 11])
                mq, nm = float(cols[base + 13]), float(cols[base + 17])
                if dp < rf + rr + vd:
                    raise AssertionError(f"{name} .var ({blk_name}): DP ({dp}) < RD ({rf}+{rr}) + VD ({vd}): {row}")
                if vf + vr != vd:
                    raise AssertionError(f"{name} .var ({blk_name}): Fwd+Rev ({vf}+{vr}) != VD ({vd}): {row}")
                if min(pmean, qmean, mq, nm) < 0.0:
                    raise AssertionError(f"{name} .var ({blk_name}): negative mean ({pmean}, {qmean}, {mq}, {nm}): {row}")
        total_var_rows += 1

    # Audit 100% of .vcf rows
    for row in c_lines:
        cols = row.split("\t")
        qual = float(cols[5])
        if qual < 0.0:
            raise AssertionError(f"{name} .vcf: negative QUAL ({qual}): {row}")
        info_kv = dict(x.split("=", 1) for x in cols[7].split(";") if "=" in x)
        if "HICNT" in info_kv:
            hicnt = int(info_kv["HICNT"])
            vd_info = int(info_kv["VD"])
            if not (0 <= hicnt <= vd_info):
                raise AssertionError(f"{name} .vcf: INFO HICNT ({hicnt}) outside [0, VD={vd_info}]: {row}")
        fmt_keys = cols[8].split(":")
        for smp in cols[9:]:
            vals = dict(zip(fmt_keys, smp.split(":")))
            dp = int(vals["DP"])
            vd = int(vals["VD"])
            ad = [int(x) for x in vals["AD"].split(",")]
            ald = [int(x) for x in vals["ALD"].split(",")]
            if dp < sum(ad):
                raise AssertionError(f"{name} .vcf: DP ({dp}) < AD[0]+AD[1] ({ad}): {row}")
            if sum(ald) != vd:
                raise AssertionError(f"{name} .vcf: ALD[0]+ALD[1] ({ald}) != VD ({vd}): {row}")
        total_vcf_rows += 1

print(f"  PASS Gate 3: Audited {total_var_rows} .var rows & {total_vcf_rows} .vcf rows across all 6 cohorts (0 violations)")
PY

# =============================================================================
# Gate 4: IUPAC FASTA Reduction & Symbolic <...> Bracket Safety
#         (AstraZeneca-NGS/VarDictJava#267, AstraZeneca-NGS/VarDictJava#287,
#          bcbio/bcbio-nextgen#3321)
# =============================================================================
echo "[4/7] Verifying IUPAC FASTA reduction & <DEL>/<DUP>/<INV> bracket safety via bcftools norm --check-ref e..."
python3 - "$VARDICTRUST_BIN" "$FA_HG19" "$VJ_BAMS" "$WORK_TMP" <<'PY'
import gzip, os, re, subprocess, sys

vdr, fa_hg19, vj_bams, tmp = sys.argv[1:5]

# Extract chr10 and chr20 from hg19.fa and inject IUPAC ambiguity bases at indel/SNV/SV anchor loci
iupac_fa = os.path.join(tmp, "hg19_iupac_chr10_chr20.fa")
raw_fa = subprocess.check_output(["samtools", "faidx", fa_hg19, "chr10", "chr20"], text=True)

seqs = {}
cur_name = None
buf = []
for line in raw_fa.splitlines():
    if line.startswith(">"):
        if cur_name:
            seqs[cur_name] = bytearray("".join(buf).encode("ascii"))
        cur_name = line[1:].split()[0]
        buf = []
    else:
        buf.append(line.strip())
if cur_name:
    seqs[cur_name] = bytearray("".join(buf).encode("ascii"))

# Inject IUPAC ambiguity codes (R, Y, M, K, S, W, B, D, H, V) across Colo829 DEL/DUP/INV windows
# including exact upstream anchors and interior positions where reads carry A/C/G/T bases
iupac_cycle = b"RYMKSWBDHV"
for cname, spans in (
    ("chr10", [(89700390, 89700650), (89713700, 89713960)]),
    ("chr20", [(46841560, 46841850), (47124940, 47125250), (47131950, 47132040)]),
):
    arr = seqs[cname]
    k = 0
    for s_pos, e_pos in spans:
        for p1 in range(s_pos, e_pos, 17):
            arr[p1 - 1] = iupac_cycle[k % len(iupac_cycle)]
            k += 1

with open(iupac_fa, "w", encoding="utf-8") as fh:
    for cname in ("chr10", "chr20"):
        fh.write(f">{cname}\n")
        data = seqs[cname]
        for off in range(0, len(data), 60):
            fh.write(data[off : off + 60].decode("ascii") + "\n")
subprocess.run(["samtools", "faidx", iupac_fa], check=True)

cases = [
    ("IUPAC_DEL", f"{vj_bams}/Colo829.chr10.DELDiscReverse.bam", ["-R", "chr10:89700380-89713976", "-f", "0.001", "-L", "500"]),
    ("IUPAC_DUP", f"{vj_bams}/Colo829.chr20.DUPDiscSV.bam", ["-R", "chr20:47124929-47131802", "-f", "0.001", "-L", "500"]),
    ("IUPAC_INV", f"{vj_bams}/Colo829.chr20.INVSV.bam", ["-R", "chr20:46841554-46845755", "-f", "0.001", "-L", "500"]),
    ("IUPAC_SHORT_DUP_MAT", f"{vj_bams}/Colo829.chr20.MarkingDUPSV.bam", ["-R", "chr20:47131954-47132045", "-f", "0.001"]),
]

valid_alt_re = re.compile(r"^(?:[ACGTN]+|<[A-Za-z0-9_-]+>)(?:,(?:[ACGTN]+|<[A-Za-z0-9_-]+>))*$")
corrupt_tokens = ("<INN>", "<NEL>", "<AUP>", "<INS_CORRUPT>")
seen_svs = set()
checked_rows = 0

for label, bam_p, extra in cases:
    out_gz = os.path.join(tmp, f"{label}.vcf.gz")
    subprocess.run(
        [vdr, "-G", iupac_fa, "-b", bam_p, "--compat=off", "--vcfall", "--output", out_gz, "--tbi", *extra],
        check=True,
    )
    # bcftools norm --check-ref e validates REF against IUPAC FASTA & rejects invalid VCF alleles
    subprocess.run(
        ["bcftools", "norm", "--check-ref", "e", "-f", iupac_fa, out_gz],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=True,
    )
    with gzip.open(out_gz, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            ref_a, alt_a = cols[3], cols[4]
            for bad_tok in corrupt_tokens:
                if bad_tok in alt_a:
                    raise AssertionError(f"{label}: bracket corruption {bad_tok} in ALT={alt_a}")
            if not valid_alt_re.match(alt_a):
                raise AssertionError(f"{label}: non-canonical/IUPAC base leaked into ALT={alt_a}: {line}")
            for single_alt in alt_a.split(","):
                if ref_a == single_alt:
                    raise AssertionError(f"{label}: non-pileup REF == ALT ({ref_a}) leaked: {line}")
                if single_alt.startswith("<"):
                    seen_svs.add(single_alt)
            checked_rows += 1

for req in ("<DEL>", "<DUP>", "<INV>"):
    if req not in seen_svs:
        raise AssertionError(f"Expected intact {req} across IUPAC FASTA runs, got {seen_svs}")
print(f"  PASS Gate 4: {checked_rows} rows passed bcftools norm --check-ref e (0 IUPAC ALT leaks, 0 REF==ALT, intact {sorted(seen_svs)})")
PY

# =============================================================================
# Gate 5: Bounded -L (--min-sv-len) Enforcement & MAX_EXPLICIT_SV_SEQ_BP Ceiling
#         (AstraZeneca-NGS/VarDictJava#328, AstraZeneca-NGS/VarDict#157)
# =============================================================================
echo "[5/7] Verifying bounded -L (--min-sv-len) materialization, discordant suppression, & -L 1000000 ceiling..."
python3 - "$VARDICTRUST_BIN" "$FA_HG19" "$FA_HS38DH" "$VJ_BAMS" "$WORK_TMP" <<'PY'
import gzip, os, re, subprocess, sys

vdr, fa_hg19, fa_hs38dh, vj_bams, tmp = sys.argv[1:6]
dna_re = re.compile(r"^[ACGTN]+$")

# 5a) Sub-1000 bp split-refined <DUP> and <INV> materialize explicit sequence alleles under default -L 1000
sub1k_split_cases = [
    ("Colo829_DUP843", fa_hg19, f"{vj_bams}/Colo829.chr20.DUPDiscSV.bam", "chr20:47124929-47131802"),
    ("Colo829_MarkingDUP", fa_hg19, f"{vj_bams}/Colo829.chr20.MarkingDUPSV.bam", "chr20:47131954-47132045"),
    ("NA12878_INVSimple", fa_hs38dh, f"{vj_bams}/NA12878.chrom20.INVSimpleSV.bam", "20:1961100-1961400"),
    ("NA12878_INVFull", fa_hs38dh, f"{vj_bams}/NA12878.chrom20.INVFullSV.bam", "20:5926940-5927179"),
]

materialized_total = 0
for name, ref_fa, bam_p, reg in sub1k_split_cases:
    out_gz = os.path.join(tmp, f"gate5_{name}_L1000.vcf.gz")
    subprocess.run(
        [vdr, "-G", ref_fa, "-b", bam_p, "-R", reg, "-f", "0.001", "--compat=off", "--vcfall", "--output", out_gz, "--tbi"],
        check=True,
    )
    subprocess.run(
        ["bcftools", "norm", "--check-ref", "e", "-f", ref_fa, out_gz],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=True,
    )
    explicit_long = 0
    with gzip.open(out_gz, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            alt = cols[4]
            if "<DUP" in alt or "<INV" in alt or "<DEL" in alt:
                raise AssertionError(f"{name}: sub-1000 bp SV remained symbolic under default -L 1000: {line}")
            if max(len(cols[3]), len(alt)) >= 20 and dna_re.match(cols[3]) and dna_re.match(alt):
                explicit_long += 1
    if explicit_long < 1:
        raise AssertionError(f"{name}: expected >= 1 explicit >=20 bp materialized SV allele under -L 1000")
    materialized_total += explicit_long

# 5b) Discordant-only clusters (SPLITREAD=0) shorter than -L are suppressed, and rescued when -L <= span
disc_cases = [
    ("NA21144_DUPDiscF_266bp", fa_hs38dh, f"{vj_bams}/NA21144.chrom11.DUPDiscF.bam", "11:77614410-77614671", "1000", "50"),
    ("NA21144_INV5Disc_11708bp", fa_hs38dh, f"{vj_bams}/NA21144.chrom11.INV5Disc.bam", "11:57970876-57971275", "15000", "1000"),
]
for name, ref_fa, bam_p, reg, suppress_l, rescue_l in disc_cases:
    vcf_suppressed = subprocess.check_output(
        [vdr, "-G", ref_fa, "-b", bam_p, "-R", reg, "-f", "0.001", "-L", suppress_l, "--compat=off", "--vcfall", "-v"],
        text=True,
    )
    if any(
        ("<DUP>" in l or "<INV>" in l)
        for l in vcf_suppressed.splitlines() if not l.startswith("#")
    ):
        raise AssertionError(f"{name}: discordant-only SV < -L {suppress_l} leaked:\n{vcf_suppressed}")
    vcf_rescued = subprocess.check_output(
        [vdr, "-G", ref_fa, "-b", bam_p, "-R", reg, "-f", "0.001", "-L", rescue_l, "--compat=off", "--vcfall", "-v"],
        text=True,
    )
    if not any(
        ("<DUP>" in l or "<INV>" in l)
        for l in vcf_rescued.splitlines() if not l.startswith("#")
    ):
        raise AssertionError(f"{name}: expected discordant-only SV to be emitted under -L {rescue_l}:\n{vcf_rescued}")

# 5c) -L 1000000 ceiling: > MAX_EXPLICIT_SV_SEQ_BP (32,768 bp) SVs stay cleanly bounded as symbolic <DEL>/<DUP>/<INV>
span3_bam = os.path.join(tmp, "span3_sv.bam")
ceil_gz = os.path.join(tmp, "gate5_L1000000.vcf.gz")
subprocess.run(
    [vdr, "-G", fa_hg19, "-b", span3_bam, "-R", "chr20:46841200-46895000",
     "-f", "0.001", "-L", "1000000", "--compat=off", "--vcfall", "--output", ceil_gz, "--tbi"],
    check=True,
)
subprocess.run(
    ["bcftools", "norm", "--check-ref", "e", "-f", fa_hg19, ceil_gz],
    stdout=subprocess.DEVNULL,
    stderr=subprocess.DEVNULL,
    check=True,
)
ceil_svs = set()
with gzip.open(ceil_gz, "rt", encoding="utf-8") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        cols = line.rstrip("\n").split("\t")
        if max(len(cols[3]), len(cols[4])) > 32769:
            raise AssertionError(f"-L 1000000 exceeded MAX_EXPLICIT_SV_SEQ_BP (32,768 bp): len={max(len(cols[3]), len(cols[4]))}")
        if cols[4] in ("<DEL>", "<DUP>", "<INV>"):
            ceil_svs.add(cols[4])
if ceil_svs != {"<DEL>", "<DUP>", "<INV>"}:
    raise AssertionError(f"Expected >32.7kb SVs to remain symbolic {{<DEL>,<DUP>,<INV>}} under -L 1000000, got {ceil_svs}")

print(f"  PASS Gate 5: {materialized_total} sub-1kb split SVs materialized, discordant <1kb suppressed, -L 1000000 bounded at 32,768 bp")
PY

# =============================================================================
# Gate 6: -p Multi-Deletion Homopolymer Coordinate Retention (Untiled == Tiled)
#         (AstraZeneca-NGS/VarDictJava#330)
# =============================================================================
echo "[6/7] Verifying -p multi-deletion homopolymer coordinate retention across untiled & split-tile runs..."
python3 - "$VARDICTRUST_BIN" "$FA_B37" "$GIAB_REF" "$VJ_BAMS" "$WAVE4_DIR" <<'PY'
import os, subprocess, sys

vdr, fa_b37, giab_ref, vj_bams, w4_dir = sys.argv[1:6]

pileup_cases = [
    (
        "N13532_TwoDels1Ins",
        fa_b37,
        os.path.join(vj_bams, "N13532.TwoDels1InsStringBoundError.bam"),
        "11",
        1279500,
        1281000,
        300,
        True,
    ),
    (
        "threeIndels_Homopolymer",
        fa_b37,
        os.path.join(vj_bams, "threeIndels_zero_deletion_tumor.bam"),
        "20",
        62361250,
        62361500,
        80,
        False,
    ),
    (
        "HG002_CMRG_Pileup",
        giab_ref,
        os.path.join(w4_dir, "hg002_cmrg_mosaic", "HG002.cmrg_mosaic.slice.bam"),
        "chr1",
        1019800,
        1021200,
        350,
        True,
    ),
]

for name, ref_fa, bam_p, chrom, s_pos, e_pos, tile_sz, full_read_span in pileup_cases:
    reg = f"{chrom}:{s_pos}-{e_pos}"
    # Determine 1-based positions covered by >= 1 non-duplicate MAPQ>=1 primary alignment
    depth_out = subprocess.check_output(["samtools", "depth", "-Q", "1", "-r", reg, bam_p], text=True)
    covered = {int(l.split("\t")[1]) for l in depth_out.splitlines() if l.strip() and int(l.split("\t")[2]) > 0}
    if len(covered) < 50:
        raise AssertionError(f"{name}: expected >= 50 covered bases in {reg}, got {len(covered)}")

    var_u = subprocess.check_output(
        [vdr, "-G", ref_fa, "-b", bam_p, "-R", reg, "-p", "--compat=off", "--vcfall", "--emit-intermediate"],
        text=True,
    )
    var_t = subprocess.check_output(
        [vdr, "-G", ref_fa, "-b", bam_p, "-R", reg, "-p", "--compat=off", "--vcfall", "--emit-intermediate",
         "--tile-size", str(tile_sz), "-th", "4"],
        text=True,
    )
    vcf_u = subprocess.check_output(
        [vdr, "-G", ref_fa, "-b", bam_p, "-R", reg, "-p", "--compat=off", "--vcfall", "-v"],
        text=True,
    )
    vcf_t = subprocess.check_output(
        [vdr, "-G", ref_fa, "-b", bam_p, "-R", reg, "-p", "--compat=off", "--vcfall", "-v",
         "--tile-size", str(tile_sz), "-th", "4"],
        text=True,
    )
    if full_read_span:
        if var_u != var_t:
            raise AssertionError(f"{name}: -p .var untiled vs --tile-size {tile_sz} mismatch")
        if vcf_u != vcf_t:
            raise AssertionError(f"{name}: -p .vcf untiled vs --tile-size {tile_sz} mismatch")

    for run_kind, vcf_body in (("untiled", vcf_u), (f"tiled-{tile_sz}", vcf_t)):
        positions = []
        dot_positions = set()
        for line in vcf_body.splitlines():
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split("\t")
            p = int(cols[1])
            if p <= 0 or p < s_pos or p > e_pos:
                raise AssertionError(f"{name} ({run_kind}): out-of-bounds or pos-0 sentinel coordinate {p} in {reg}")
            positions.append(p)
            if cols[4] == ".":
                if p in dot_positions:
                    raise AssertionError(f"{name} ({run_kind}): duplicate ALT=. reference row at {chrom}:{p}")
                dot_positions.add(p)

        if any(positions[i] > positions[i + 1] for i in range(len(positions) - 1)):
            raise AssertionError(f"{name} ({run_kind}): non-monotonic -p coordinates")
        emitted_set = set(positions)
        missing = sorted(covered - emitted_set)
        if missing:
            raise AssertionError(f"{name} ({run_kind}): -p dropped covered coordinates {missing[:10]} in {reg}")

print("  PASS Gate 6: 100% covered homopolymer coordinates retained monotonically (0 pos-0/dup rows, untiled == tiled)")
PY

# =============================================================================
# Gate 7: Strict VCF Schema (AD / END / GT / TYPE) & Contig-Scoped rtg vcfeval
#         Without --squash-ploidy Proving Subclonal Somatic TP Recovery
#         (clintval/vartovcf#99, clintval/vartovcf#101, clintval/vartovcf#105,
#          clintval/vartovcf#109, EBIvariation/vcf-validator#99)
# =============================================================================
echo "[7/7] Verifying strict VCF schema (AD/END/GT/TYPE) & contig-scoped rtg vcfeval without --squash-ploidy..."
python3 - "$VARDICTRUST_BIN" "$GRCH38_REF" "$FA_HG19" "$FA_B37" "$VJ_BAMS" "$G4_DIR" "$WAVE4_DIR" "$WORK_TMP" <<'PY'
import gzip, os, re, subprocess, sys

vdr, grch38_ref, fa_hg19, fa_b37, vj_bams, g4_dir, w4_dir, tmp = sys.argv[1:9]

def expected_vartype(ref_a, alt_a):
    if alt_a == ".":
        return ""
    if alt_a.startswith("<"):
        inner = alt_a[1:-1]
        return inner.split("-", 1)[0]
    if len(ref_a) == len(alt_a):
        return "SNV" if len(ref_a) == 1 else "Complex"
    if len(ref_a) == 1 and len(alt_a) > 1 and alt_a.startswith(ref_a):
        return "Insertion"
    if len(alt_a) == 1 and len(ref_a) > 1 and ref_a.startswith(alt_a):
        return "Deletion"
    # Shared prefix/suffix pure insertion or deletion vs Complex
    r, a = ref_a, alt_a
    while len(r) > 1 and len(a) > 1 and r[0] == a[0]:
        r, a = r[1:], a[1:]
    while len(r) > 1 and len(a) > 1 and r[-1] == a[-1]:
        r, a = r[:-1], a[:-1]
    if len(r) == 1 and len(a) > 1 and a.startswith(r):
        return "Insertion"
    if len(a) == 1 and len(r) > 1 and r.startswith(a):
        return "Deletion"
    return "Complex"

def audit_strict_vcf(label, vcf_text, paired):
    rows = 0
    sv_end_seen = 0
    for line in vcf_text.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        cols = line.split("\t")
        ref_a, alt_a, info, fmt_s = cols[3], cols[4], cols[7], cols[8]
        kv = dict(x.split("=", 1) for x in info.split(";") if "=" in x)
        is_sv = alt_a.startswith("<")
        if is_sv:
            if "END" not in kv:
                raise AssertionError(f"{label}: symbolic SV missing INFO/END: {line}")
            sv_end_seen += 1
        else:
            if "END" in kv:
                raise AssertionError(f"{label}: non-SV record leaked INFO/END: {line}")
        # Single-allele post-adj_complex TYPE check
        if "," not in alt_a and alt_a != "." and not is_sv:
            exp_t = expected_vartype(ref_a, alt_a)
            got_t = kv.get("TYPE", "")
            if got_t != exp_t:
                raise AssertionError(f"{label}: stale INFO/TYPE={got_t} (expected {exp_t}) for {ref_a}->{alt_a}: {line}")

        fmt = fmt_s.split(":")
        n_alt = 0 if alt_a == "." else len(alt_a.split(","))
        any_vd_pos = False
        any_alt_gt = False
        for smp in cols[9:]:
            vals = dict(zip(fmt, smp.split(":")))
            if vals["GT"] == "1/0":
                raise AssertionError(f"{label}: non-canonical 1/0 genotype leaked: {line}")
            if "1" in vals["GT"]:
                any_alt_gt = True
            if int(vals["VD"]) > 0:
                any_vd_pos = True
            ad_len = len(vals["AD"].split(","))
            exp_ad = 2 if paired else (1 + n_alt)
            if ad_len != exp_ad:
                raise AssertionError(f"{label}: AD length {ad_len} != {exp_ad} on ALT={alt_a}: {line}")
        if alt_a != "." and any_vd_pos and not any_alt_gt:
            raise AssertionError(f"{label}: unobserved ALT (all samples GT=0/0 with VD>0): {line}")
        rows += 1
    return rows, sv_end_seen

# 7a) Audit single-sample -R (omit_end_tag=false!), -p pileup (ALT=.), symbolic SV, and paired somatic VCFs
single_r_vcf = subprocess.check_output(
    [vdr, "-G", fa_hg19, "-b", os.path.join(vj_bams, "Colo829.chr10.DELDiscReverse.bam"),
     "-R", "chr10:89700380-89713976", "-f", "0.001", "--compat=off", "-v"],
    text=True,
)
n_sr, n_sv_end = audit_strict_vcf("single_R_sv", single_r_vcf, paired=False)
if n_sv_end < 1:
    raise AssertionError("Expected >= 1 symbolic SV with INFO/END in single_R_sv")

single_p_vcf = subprocess.check_output(
    [vdr, "-G", fa_b37, "-b", os.path.join(vj_bams, "N13532.TwoDels1InsStringBoundError.bam"),
     "-R", "11:1280138-1280341", "-p", "--compat=off", "-v"],
    text=True,
)
n_sp, _ = audit_strict_vcf("single_pileup", single_p_vcf, paired=False)

# 7b) Contig-scoped rtg format on samtools faidx-extracted target contigs + rtg vcfeval WITHOUT --squash-ploidy
target_contigs = ["chr7", "chr9", "chr11", "chr12", "chr17", "chr22"]
scoped_fa = os.path.join(tmp, "seqc2_target_contigs.fa")
with open(scoped_fa, "w", encoding="utf-8") as fh:
    subprocess.run(["samtools", "faidx", grch38_ref, *target_contigs], stdout=fh, check=True)
subprocess.run(["samtools", "faidx", scoped_fa], check=True)

scoped_sdf = os.path.join(tmp, "seqc2_target_contigs.sdf")
subprocess.run(["rtg", "format", "-o", scoped_sdf, scoped_fa], stdout=subprocess.DEVNULL, check=True)

# Intersect seqc2/slices.bed with seqc2/eval.bed so rtg vcfeval evaluates exact high-confidence intervals
slices_bed = os.path.join(g4_dir, "seqc2", "slices.bed")
hc_eval_bed = os.path.join(g4_dir, "seqc2", "eval.bed")
scoped_eval_bed = os.path.join(tmp, "seqc2_scoped_eval.bed")

slice_ivs = {}
with open(slices_bed, "r", encoding="utf-8") as fh:
    for line in fh:
        if line.strip() and not line.startswith("#"):
            p = line.split("\t")
            slice_ivs.setdefault(p[0], []).append((int(p[1]), int(p[2])))

with open(hc_eval_bed, "r", encoding="utf-8") as fh, open(scoped_eval_bed, "w", encoding="utf-8") as out:
    for line in fh:
        if not line.strip() or line.startswith("#"):
            continue
        p = line.split("\t")
        c = p[0]
        if c not in slice_ivs:
            continue
        s, e = int(p[1]), int(p[2])
        for ss, se in slice_ivs[c]:
            lo, hi = max(s, ss), min(e, se)
            if hi > lo:
                out.write(f"{c}\t{lo}\t{hi}\n")

# Build diploid somatic truth baseline VCF (GT=0/1 for TVAF < 0.80, GT=1/1 for TVAF >= 0.80) on target contigs
truth_src = os.path.join(g4_dir, "seqc2", "truth", "truth.vcf.gz")
truth_scoped_vcf = os.path.join(tmp, "seqc2_truth_diploid.vcf")
truth_scoped_gz = truth_scoped_vcf + ".gz"
raw_truth = subprocess.check_output(["bcftools", "view", "-R", slices_bed, truth_src], text=True)
tvaf_re = re.compile(r"(?:^|;)TVAF=([0-9.eE+-]+)")
contig_lens = {}
with open(scoped_fa + ".fai", "r", encoding="utf-8") as fh:
    for line in fh:
        p = line.split("\t")
        contig_lens[p[0]] = int(p[1])

subclonal_truth_n = 0
with open(truth_scoped_vcf, "w", encoding="utf-8") as out:
    out.write("##fileformat=VCFv4.2\n")
    for c in target_contigs:
        out.write(f"##contig=<ID={c},length={contig_lens[c]}>\n")
    out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\n")
    for line in raw_truth.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        p = line.split("\t")
        if p[0] not in contig_lens:
            continue
        m = tvaf_re.search(p[7])
        tvaf = float(m.group(1)) if m else 0.30
        if tvaf < 0.30:
            subclonal_truth_n += 1
        gt = "1/1" if tvaf >= 0.80 else "0/1"
        out.write(f"{p[0]}\t{p[1]}\t{p[2]}\t{p[3]}\t{p[4]}\t.\tPASS\t.\tGT\t{gt}\n")

subprocess.run(["bgzip", "-f", truth_scoped_vcf], check=True)
subprocess.run(["tabix", "-f", "-p", "vcf", truth_scoped_gz], check=True)

# Call SEQC2 HCC1395|HCC1395BL with discovery -f 0.01 and VCF AF threshold --vcf-af 0.30:
# Subclonal somatic calls (0.01 <= AF < 0.30) receive GT=0/0 under --compat=vardictjava
# (dropped as CALL=IGN by rtg vcfeval without --squash-ploidy!) vs GT=0/1 under --compat=off!
bam_pair = f"{g4_dir}/seqc2/HCC1395.bam|{g4_dir}/seqc2/HCC1395BL.bam"

def run_seqc2_vcfeval(mode_tag, compat_flag):
    raw_vcf = os.path.join(tmp, f"seqc2_{mode_tag}_raw.vcf")
    som_vcf = os.path.join(tmp, f"seqc2_{mode_tag}_somatic.vcf")
    som_gz = som_vcf + ".gz"
    subprocess.run(
        [vdr, "-G", scoped_fa, "-b", bam_pair, "-N", "TUMOR|NORMAL",
         "-c", "1", "-S", "2", "-E", "3", "-g", "4",
         "-f", "0.01", "--vcf-af", "0.30", "--fisher", compat_flag,
         "--output", raw_vcf, "-th", "4", slices_bed],
        check=True,
    )
    with open(raw_vcf, "r", encoding="utf-8") as fh:
        raw_txt = fh.read()
    if compat_flag == "--compat=off":
        audit_strict_vcf("seqc2_compat_off", raw_txt, paired=True)

    gt00_somatic = 0
    with open(som_vcf, "w", encoding="utf-8") as out:
        for line in raw_txt.splitlines():
            if line.startswith("#"):
                out.write(line + "\n")
                continue
            p = line.split("\t")
            if "STATUS=StrongSomatic" in p[7] or "STATUS=LikelySomatic" in p[7]:
                t_fields = p[9].split(":")
                n_fields = p[10].split(":")
                if t_fields[0] == "0/0" and n_fields[0] == "0/0" and (int(t_fields[2]) > 0 or int(n_fields[2]) > 0):
                    gt00_somatic += 1
                out.write(line + "\n")
    subprocess.run(["bgzip", "-f", som_vcf], check=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", som_gz], check=True)

    eval_dir = os.path.join(tmp, f"vcfeval_{mode_tag}")
    # Deliberately omit --squash-ploidy to verify diploid GT=0/1 agreement!
    subprocess.run(
        ["rtg", "vcfeval",
         "-b", truth_scoped_gz,
         "-c", som_gz,
         "-e", scoped_eval_bed,
         "-t", scoped_sdf,
         "--all-records",
         "--vcf-score-field=QUAL",
         "--sample=TUMOR,TUMOR",
         "-o", eval_dir],
        stdout=subprocess.DEVNULL,
        check=True,
    )
    tp_n = 0
    with gzip.open(os.path.join(eval_dir, "tp.vcf.gz"), "rt", encoding="utf-8") as fh:
        for line in fh:
            if not line.startswith("#"):
                tp_n += 1
    return tp_n, gt00_somatic

tp_vj, gt00_vj = run_seqc2_vcfeval("vj", "--compat=vardictjava")
tp_off, gt00_off = run_seqc2_vcfeval("off", "--compat=off")

if gt00_vj < 10:
    raise AssertionError(f"Expected >= 10 GT=0/0 subclonal somatic calls under --compat=vardictjava, got {gt00_vj}")
if gt00_off != 0:
    raise AssertionError(f"Expected 0 GT=0/0 somatic tumor calls under --compat=off, got {gt00_off}")
if tp_off <= tp_vj:
    raise AssertionError(
        f"Expected --compat=off rtg vcfeval (no --squash-ploidy) to recover subclonal somatic TPs over --compat=vardictjava: tp_off={tp_off} vs tp_vj={tp_vj}"
    )

print(
    f"  PASS Gate 7: Strict schema audited ({n_sr}+{n_sp}+SEQC2 rows); "
    f"rtg vcfeval (no --squash-ploidy, {subclonal_truth_n} subclonal truth sites): "
    f"TP {tp_vj} -> {tp_off} (+{tp_off - tp_vj} recovered; GT=0/0 somatic {gt00_vj} -> {gt00_off})"
)
PY

ELAPSED_S=$(( ( $(now_ms) - T0_MS ) / 1000 ))
echo "================================================================================"
echo " ALL 7 WAVE-5 PUBLIC CLINICAL & VALIDATOR GATES PASSED (${ELAPSED_S}s)"
echo "================================================================================"
