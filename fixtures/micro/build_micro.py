#!/usr/bin/env python3
"""Build the self-contained micro fixture (fixtures/micro/).

Deterministic: fixed RNG seed. Emits micro.ref.fa, micro.bed, micro.sam.
Convert SAM -> BAM with:
  samtools view -b micro.sam | samtools sort -o micro.bam - && samtools index micro.bam
"""
import random

random.seed(42)

CHR1_LEN = 300
CHR2_LEN = 260


def rand_seq(n: int) -> str:
    return "".join(random.choice("ACGT") for _ in range(n))


chr1 = rand_seq(CHR1_LEN)
chr2 = rand_seq(CHR2_LEN)

# --- variants -------------------------------------------------------------
# SNV at chr1 (0-based) pos 100 in ~half the reads covering it.
SNV_POS = 100
alt_base = next(b for b in "ACGT" if b != chr1[SNV_POS])
# 3 bp deletion of chr2 (0-based) [150, 153) in ~half the covering reads.
DEL_START, DEL_END = 150, 153

READ_LEN = 60

reads = []  # (name, chr, start0, seq, cigar, flag)


def add_reads(chr_seq: str, chrom: str, first: int, last: int):
    global read_id
    for start in range(first, last, 4):  # step 4 -> ~15x depth mid-window
        end = min(start + READ_LEN, len(chr_seq))
        if end - start < 30:
            continue
        seq = list(chr_seq[start:end])
        cig_last = end - start
        cigar = f"{cig_last}M"
        mut = False
        if chrom == "chr1" and start <= SNV_POS < end and (start // 4) % 2 == 0:
            seq[SNV_POS - start] = alt_base
            mut = True
        if chrom == "chr2" and start <= DEL_START and end >= DEL_END and (start // 4) % 2 == 0:
            seq = chr_seq[start:DEL_START] + chr_seq[DEL_END:end]
            cigar = f"{DEL_START - start}M{DEL_END - DEL_START}D{end - DEL_END}M"
            mut = True
        reads.append(
            (
                f"r{chrom}_{start}",
                chrom,
                start + 1,  # SAM is 1-based
                "".join(seq),
                cigar,
            )
        )


add_reads(chr1, "chr1", 40, 240)
add_reads(chr2, "chr2", 110, 200)

with open("micro.sam", "w") as sam:
    sam.write("@HD\tVN:1.6\tSO:coordinate\n")
    sam.write(f"@SQ\tSN:chr1\tLN:{CHR1_LEN}\n")
    sam.write(f"@SQ\tSN:chr2\tLN:{CHR2_LEN}\n")
    for name, chrom, pos, seq, cigar in sorted(reads, key=lambda r: (r[1], r[2])):
        qual = "I" * len(seq)
        sam.write(
            f"{name}\t0\t{chrom}\t{pos}\t60\t{cigar}\t*\t0\t0\t{seq}\t{qual}\n"
        )

with open("micro.ref.fa", "w") as fa:
    fa.write(">chr1\n")
    for i in range(0, CHR1_LEN, 60):
        fa.write(chr1[i : i + 60] + "\n")
    fa.write(">chr2\n")
    for i in range(0, CHR2_LEN, 60):
        fa.write(chr2[i : i + 60] + "\n")

with open("micro.bed", "w") as bed:
    bed.write("chr1\t40\t260\tmicro_chr1\n")
    bed.write("chr2\t110\t230\tmicro_chr2\n")

print(f"chr1[SNV_POS]={chr1[SNV_POS]} -> {alt_base}; del chr2[{DEL_START},{DEL_END})")
print(f"reads: {len(reads)}")
