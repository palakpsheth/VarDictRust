#!/usr/bin/env python3
"""Unit tests for harness/select_slice_regions.py."""
import contextlib
import io
import os
import sys
import tempfile
import unittest

HARNESS_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if HARNESS_DIR not in sys.path:
    sys.path.insert(0, HARNESS_DIR)

import select_slice_regions


def make_vcf(path, records):
    with open(path, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, ref, alt in records:
            fh.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t30\tPASS\t.\n")


def make_bed(path, intervals):
    with open(path, "w") as fh:
        for chrom, start, end in intervals:
            fh.write(f"{chrom}\t{start}\t{end}\n")


class TestSelectSliceRegions(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.d = self.tmpdir.name

    def tearDown(self):
        self.tmpdir.cleanup()

    def run_tool(self, args):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            select_slice_regions.main(args)
        out = buf.getvalue().strip()
        return out.splitlines() if out else []

    def test_min_count_and_padding(self):
        """6 chr1 variants within 1 kb plus 2 chr2 variants, with --min-count 5 -> exactly one chr1 region, padded."""
        # 6 chr1 variants at 0-based positions 1000, 1100, 1200, 1300, 1400, 1500 (VCF 1-based POS: 1001..1501)
        # 2 chr2 variants at 0-based positions 1000, 1200 (VCF 1-based POS: 1001, 1201)
        chr1_vars = [("chr1", p, "A", "T") for p in [1001, 1101, 1201, 1301, 1401, 1501]]
        chr2_vars = [("chr2", p, "A", "T") for p in [1001, 1201]]
        vcf_path = os.path.join(self.d, "truth.vcf")
        make_vcf(vcf_path, chr1_vars + chr2_vars)

        eval_bed = os.path.join(self.d, "eval.bed")
        make_bed(eval_bed, [("chr1", 0, 100000), ("chr2", 0, 100000)])

        args = [
            "--truth", vcf_path,
            "--eval", eval_bed,
            "--chain", "500",
            "--min-count", "5",
            "--pad", "1000",
            "--max-span", "5000",
            "--max-regions", "10",
        ]
        lines = self.run_tool(args)
        self.assertEqual(len(lines), 1)
        fields = lines[0].split("\t")
        self.assertEqual(fields[0], "chr1")
        # chr1: first pos 0-based = 1000, pad = 1000 -> 0; last pos 0-based = 1500, rl = 1 -> 1501 + 1000 = 2501
        self.assertEqual(fields[1], "0")
        self.assertEqual(fields[2], "2501")
        self.assertEqual(fields[3], "slice_1_6")

    def test_exclude_bed_covering_cluster(self):
        """An exclude BED covering the chr1 cluster -> SystemExit."""
        chr1_vars = [("chr1", p, "A", "T") for p in [1001, 1101, 1201, 1301, 1401, 1501]]
        chr2_vars = [("chr2", p, "A", "T") for p in [1001, 1201]]
        vcf_path = os.path.join(self.d, "truth.vcf")
        make_vcf(vcf_path, chr1_vars + chr2_vars)

        eval_bed = os.path.join(self.d, "eval.bed")
        make_bed(eval_bed, [("chr1", 0, 100000), ("chr2", 0, 100000)])

        exclude_bed = os.path.join(self.d, "exclude.bed")
        make_bed(exclude_bed, [("chr1", 500, 2000)])

        args = [
            "--truth", vcf_path,
            "--eval", eval_bed,
            "--exclude", exclude_bed,
            "--chain", "500",
            "--min-count", "5",
            "--pad", "1000",
            "--max-span", "5000",
            "--max-regions", "10",
        ]
        with self.assertRaises(SystemExit):
            self.run_tool(args)

    def test_eval_bed_drops_variants(self):
        """An eval BED that drops 2 of the 6 chr1 variants -> SystemExit."""
        chr1_vars = [("chr1", p, "A", "T") for p in [1001, 1101, 1201, 1301, 1401, 1501]]
        chr2_vars = [("chr2", p, "A", "T") for p in [1001, 1201]]
        vcf_path = os.path.join(self.d, "truth.vcf")
        make_vcf(vcf_path, chr1_vars + chr2_vars)

        # Eval BED covers only 1150-2000 on chr1, dropping 1000 and 1100 (leaving 4 variants, < min-count 5)
        eval_bed = os.path.join(self.d, "eval.bed")
        make_bed(eval_bed, [("chr1", 1150, 2000), ("chr2", 0, 100000)])

        args = [
            "--truth", vcf_path,
            "--eval", eval_bed,
            "--chain", "500",
            "--min-count", "5",
            "--pad", "1000",
            "--max-span", "5000",
            "--max-regions", "10",
        ]
        with self.assertRaises(SystemExit):
            self.run_tool(args)

    def test_equal_size_clusters_lexicographically_first(self):
        """Two equal-size clusters with --max-regions 1 -> the lexicographically first chromosome is chosen."""
        chr1_vars = [("chr1", p, "A", "T") for p in [1001, 1101, 1201, 1301, 1401]]
        chr2_vars = [("chr2", p, "A", "T") for p in [1001, 1101, 1201, 1301, 1401]]
        vcf_path = os.path.join(self.d, "truth.vcf")
        make_vcf(vcf_path, chr1_vars + chr2_vars)

        eval_bed = os.path.join(self.d, "eval.bed")
        make_bed(eval_bed, [("chr1", 0, 100000), ("chr2", 0, 100000)])

        args = [
            "--truth", vcf_path,
            "--eval", eval_bed,
            "--chain", "500",
            "--min-count", "5",
            "--pad", "1000",
            "--max-span", "5000",
            "--max-regions", "1",
        ]
        lines = self.run_tool(args)
        self.assertEqual(len(lines), 1)
        fields = lines[0].split("\t")
        self.assertEqual(fields[0], "chr1")
        self.assertNotIn("chr2", lines[0])


if __name__ == "__main__":
    unittest.main()
