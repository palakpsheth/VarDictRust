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
if __name__ == "__main__": unittest.main()
