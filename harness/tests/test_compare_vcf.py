import subprocess, sys, tempfile, os, unittest
HDR = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
def rec(pos, alt="T", qual="30", info="DP=10;AF=0.5"):
    return f"chr1\t{pos}\t.\tA\t{alt}\t{qual}\tPASS\t{info}\tGT:DP\t0/1:10\n"
def run(java, rust, *args):
    d = tempfile.mkdtemp(); j, r = os.path.join(d, "j.vcf"), os.path.join(d, "r.vcf")
    open(j, "w").write(HDR + "".join(java)); open(r, "w").write(HDR + "".join(rust))
    p = subprocess.run([sys.executable, "harness/compare_vcf.py", *args, j, r], capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr
def run_custom(j_hdr, r_hdr, java, rust, *args):
    d = tempfile.mkdtemp(); j, r = os.path.join(d, "j.vcf"), os.path.join(d, "r.vcf")
    with open(j, "w") as f: f.write(j_hdr + "".join(java))
    with open(r, "w") as f: f.write(r_hdr + "".join(rust))
    p = subprocess.run([sys.executable, "harness/compare_vcf.py", *args, j, r], capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr

class T(unittest.TestCase):
    def test_identical(self):
        rc, out = run([rec(1), rec(2)], [rec(2), rec(1)]); self.assertEqual(rc, 0); self.assertIn("MATCH (2 variants)", out)
    def test_duplicate_key_counts(self):
        rc, out = run([rec(1)], [rec(1), rec(1)]); self.assertEqual(rc, 1); self.assertIn("record counts differ", out)
    def test_min_variants_vacuous(self):
        self.assertEqual(run([], [], "--min-variants", "1")[0], 1)
    def test_vacuous_warns(self):
        rc, out = run([], []); self.assertEqual(rc, 0); self.assertIn("vacuous", out)
    def test_eps(self):
        # Under raw float comparison with eps=1e-6 (A5, D-52):
        self.assertEqual(run([rec(1, qual="1.0000004")], [rec(1, qual="1.0000000")])[0], 0)
        self.assertEqual(run([rec(1, qual="1.000001")], [rec(1, qual="1.0")])[0], 0)
        self.assertEqual(run([rec(1, qual="1.00001")], [rec(1, qual="1.0")])[0], 1)
    def test_set_difference(self):
        self.assertEqual(run([rec(1)], [rec(2)])[0], 1)
    def test_malformed_row(self):
        rc, out = run(["chr1\t1\t.\tA\n"], ["chr1\t1\t.\tA\n"]); self.assertEqual(rc, 1); self.assertIn("malformed", out)
    def test_header_mismatch(self):
        h1 = "##fileformat=VCFv4.2\n##INFO=<ID=FOO>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        h2 = "##fileformat=VCFv4.2\n##INFO=<ID=BAR>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        rc, out = run_custom(h1, h2, [rec(1)], [rec(1)])
        self.assertEqual(rc, 1)
        self.assertIn("HEADER line", out)
    def test_header_volatile_ignored(self):
        h1 = "##fileformat=VCFv4.2\n##fileDate=20261008\n##source=VarDict_v1.8.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        h2 = "##fileformat=VCFv4.2\n##fileDate=20261009\n##source=VarDictRust\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        rc, out = run_custom(h1, h2, [rec(1)], [rec(1)])
        self.assertEqual(rc, 0)
        self.assertIn("MATCH (1 variants)", out)
    def test_strict_bytes_reports(self):
        rc, out = run([rec(1, info="NM=1.0")], [rec(1, info="NM=1")], "--strict-bytes")
        self.assertEqual(rc, 1)
        self.assertIn("STRICT-BYTES: 1", out)
    def test_strict_bytes_vcf_formatting_diffs(self):
        for j_info, r_info in [
            ("DUPRATE=0", "DUPRATE=0.0"),
            ("MQ=2.", "MQ=2"),
            ("ODDRATIO=1.0E-4", "ODDRATIO=0.0001"),
        ]:
            self.assertEqual(run([rec(1, info=j_info)], [rec(1, info=r_info)])[0], 0)
            rc, out = run([rec(1, info=j_info)], [rec(1, info=r_info)], "--strict-bytes")
            self.assertEqual(rc, 1, f"expected strict-bytes failure for {j_info} vs {r_info}")
            self.assertIn("STRICT-BYTES: 1", out)
    def test_strict_bytes_var_formatting_diffs(self):
        base = ["S", "G", "chr1", "10", "10", "A", "T", "10", "5", "3", "2", "3", "2", "A/T", "0.5"] + ["0"] * 21
        for val_a, val_b in [("0", "0.0"), ("2.", "2"), ("1.0E-4", "0.0001")]:
            row_a = list(base); row_a[14] = val_a
            row_b = list(base); row_b[14] = val_b
            d = tempfile.mkdtemp()
            j, r = os.path.join(d, "j.var"), os.path.join(d, "r.var")
            with open(j, "w") as f: f.write("\t".join(row_a) + "\n")
            with open(r, "w") as f: f.write("\t".join(row_b) + "\n")
            p_nostrict = subprocess.run([sys.executable, "harness/compare_var.py", j, r], capture_output=True, text=True)
            self.assertEqual(p_nostrict.returncode, 0)
            p_strict = subprocess.run([sys.executable, "harness/compare_var.py", "--strict-bytes", j, r], capture_output=True, text=True)
            self.assertEqual(p_strict.returncode, 1, f"expected strict-bytes failure in compare_var.py for {val_a} vs {val_b}")

    def test_truncated_format_subfield(self):
        j_row = "chr1\t10\t.\tA\tT\t30\tPASS\tDP=100;AF=0.5\tGT:DP:VD\t0/1:100:50\n"
        r_row = "chr1\t10\t.\tA\tT\t30\tPASS\tDP=100;AF=0.5\tGT:DP:VD\t0/1:100\n"
        rc, out = run([j_row], [r_row])
        self.assertNotEqual(rc, 0)
        self.assertIn("MISMATCH", out)

    def test_extra_format_subfield(self):
        j_row = "chr1\t10\t.\tA\tT\t30\tPASS\tDP=100;AF=0.5\tGT:DP\t0/1:100\n"
        r_row = "chr1\t10\t.\tA\tT\t30\tPASS\tDP=100;AF=0.5\tGT:DP\t0/1:100:50\n"
        rc, out = run([j_row], [r_row])
        self.assertNotEqual(rc, 0)
        self.assertIn("MISMATCH", out)

    def test_dropped_sample_column(self):
        paired_hdr = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL\n"
        j_paired = "chr1\t10\t.\tA\tT\t30\tPASS\tDP=100\tGT:DP\t0/1:100\t0/0:90\n"
        r_dropped_normal = "chr1\t10\t.\tA\tT\t30\tPASS\tDP=100\tGT:DP\t0/1:100\n"
        rc1, out1 = run_custom(paired_hdr, paired_hdr, [j_paired], [r_dropped_normal])
        self.assertNotEqual(rc1, 0)
        self.assertIn("MISMATCH", out1)

        r_eight_cols = "chr1\t10\t.\tA\tT\t30\tPASS\tDP=100\n"
        rc2, out2 = run([rec(10)], [r_eight_cols])
        self.assertNotEqual(rc2, 0)
        self.assertIn("MISMATCH", out2)

    def test_malformed_nine_col_row(self):
        nine_col = "chr1\t10\t.\tA\tT\t30\tPASS\tDP=100\tGT:DP\n"
        rc, out = run([nine_col], [nine_col])
        self.assertNotEqual(rc, 0)
        self.assertIn("MISMATCH", out)

    def test_stray_header_after_data_rows(self):
        rc, out = run([rec(10), "#stray_header\n", rec(20)], [rec(10), rec(20)])
        self.assertNotEqual(rc, 0)
        self.assertIn("MISMATCH", out)

    def test_streaming_fifo_single_pass(self):
        import shutil, threading
        d = tempfile.mkdtemp()
        try:
            j_fifo = os.path.join(d, "j.fifo")
            r_fifo = os.path.join(d, "r.fifo")
            os.mkfifo(j_fifo)
            os.mkfifo(r_fifo)
            payload = HDR + rec(1) + rec(2)

            def writer(path):
                try:
                    with open(path, "w") as fh:
                        fh.write(payload)
                except BrokenPipeError:
                    pass

            tj = threading.Thread(target=writer, args=(j_fifo,), daemon=True)
            tr = threading.Thread(target=writer, args=(r_fifo,), daemon=True)
            tj.start()
            tr.start()
            p = subprocess.run(
                [sys.executable, "harness/compare_vcf.py", "--min-variants", "2", j_fifo, r_fifo],
                capture_output=True,
                text=True,
                timeout=3,
            )
            tj.join(timeout=1)
            tr.join(timeout=1)
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertIn("MATCH (2 variants)", p.stdout)
        finally:
            shutil.rmtree(d, ignore_errors=True)

if __name__ == "__main__": unittest.main()
