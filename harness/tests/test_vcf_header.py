import os
import re
import shutil
import subprocess
import tempfile
import unittest


class TestVcfHeaderParity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.keep_dir = tempfile.mkdtemp()
        env = dict(os.environ)
        env["KEEP_DIR"] = cls.keep_dir
        cmd = [
            "bash",
            "harness/compare.sh",
            "fixtures/micro/micro.ref.fa",
            "fixtures/micro/micro.bam",
            "fixtures/micro/micro.bed",
            "SAMPLE",
        ]
        cls.proc = subprocess.run(cmd, env=env, capture_output=True, text=True)

    @classmethod
    def tearDownClass(cls):
        if os.path.exists(cls.keep_dir):
            shutil.rmtree(cls.keep_dir)

    def test_headers_match_verbatim(self):
        java_vcf = os.path.join(self.keep_dir, "java.vcf")
        rust_vcf = os.path.join(self.keep_dir, "rust.vcf")
        self.assertTrue(os.path.exists(java_vcf), f"java.vcf missing: {self.proc.stderr}")
        self.assertTrue(os.path.exists(rust_vcf), f"rust.vcf missing: {self.proc.stderr}")

        with open(java_vcf) as f:
            j_lines = [line.rstrip("\r\n") for line in f if line.startswith("#")]
        with open(rust_vcf) as f:
            r_lines = [line.rstrip("\r\n") for line in f if line.startswith("#")]

        # Filter out volatile lines matching ^##(fileDate|source)=
        volatile_pattern = re.compile(r"^##(fileDate|source)=")
        j_filtered = [line for line in j_lines if not volatile_pattern.match(line)]
        r_filtered = [line for line in r_lines if not volatile_pattern.match(line)]

        # Check that required header sections are present
        self.assertTrue(any(line.startswith("##fileformat=") for line in r_filtered))
        self.assertTrue(any(line.startswith("##reference=") for line in r_filtered))
        self.assertTrue(any(line.startswith("##contig=") for line in r_filtered))
        self.assertTrue(any(line.startswith("##INFO=") for line in r_filtered))
        self.assertTrue(any(line.startswith("##FILTER=") for line in r_filtered))
        self.assertTrue(any(line.startswith("##FORMAT=") for line in r_filtered))
        self.assertTrue(any(line.startswith("#CHROM") for line in r_filtered))

        # Check PMEAN specifically against Java's exact string
        expected_pmean_desc = 'Description="The mean distance to the nearest 5 or 3 prime read end (whichever is closer) in all reads that support the variant call"'
        r_pmean_lines = [line for line in r_filtered if "ID=PMEAN" in line]
        self.assertTrue(len(r_pmean_lines) > 0, "PMEAN header line not found in rust.vcf")
        self.assertIn(
            expected_pmean_desc,
            r_pmean_lines[0],
            f"Rust PMEAN description does not match Java: {r_pmean_lines[0]}",
        )

        # Assert all header lines match verbatim between Java and Rust
        self.assertEqual(j_filtered, r_filtered)


if __name__ == "__main__":
    unittest.main()
