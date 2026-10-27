#!/usr/bin/env python3
"""Convert upstream VarDictJava testdata thin-FASTA CSVs to indexed FASTAs.

Upstream VarDictJava (AstraZeneca-NGS) only commits sparse reference sequence
chunks as CSVs (`VarDictJava/testdata/fastas/*.fa.csv`) to save repository space.
Their integration test harness parses these CSVs on the fly using Java classes.
VarDictRust and samtools require true FASTA files with companion `.fa.fai` indices.

This script reconstructs bitwise-identical FASTA files and companion `.fai` indices
directly in pure Python (zero external dependencies, runtime <10s).
"""

import fcntl
import os
import sys
import time

def convert_csv_to_fasta(csv_path: str, fa_path: str, fai_path: str):
    chrs = {}
    with open(csv_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(",")
            if len(parts) < 4:
                continue
            chrom, start, seq = parts[0], int(parts[1]), parts[3]
            if chrom not in chrs:
                chrs[chrom] = []
            chrs[chrom].append((start, seq))
    
    tmp_fa = f"{fa_path}.tmp.{os.getpid()}"
    tmp_fai = f"{fai_path}.tmp.{os.getpid()}"
    
    fai_lines = []
    with open(tmp_fa, "wb") as out_fa:
        for chrom, chunks in chrs.items():
            max_len = max(start + len(seq) - 1 for start, seq in chunks)
            seq_arr = bytearray(b"N" * max_len)
            for start, seq in chunks:
                start_idx = start - 1
                seq_bytes = seq.encode("ascii")
                seq_arr[start_idx:start_idx + len(seq_bytes)] = seq_bytes
            
            header = f">{chrom}\n".encode("ascii")
            offset = out_fa.tell() + len(header)
            out_fa.write(header)
            
            line_len = min(60, len(seq_arr))
            for i in range(0, len(seq_arr), 60):
                out_fa.write(seq_arr[i:i+60] + b"\n")
            
            fai_lines.append(f"{chrom}\t{max_len}\t{offset}\t{line_len}\t{line_len+1}\n")
        out_fa.flush()
        os.fsync(out_fa.fileno())
            
    with open(tmp_fai, "w", encoding="utf-8") as out_fai:
        out_fai.writelines(fai_lines)
        out_fai.flush()
        os.fsync(out_fai.fileno())

    # Atomic rename into final locations
    os.replace(tmp_fa, fa_path)
    os.replace(tmp_fai, fai_path)

def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    fasta_dir = os.path.join(root, "VarDictJava", "testdata", "fastas")
    if not os.path.isdir(fasta_dir):
        print(f"convert_fastas: directory not found: {fasta_dir}", file=sys.stderr)
        sys.exit(1)

    import tempfile
    lock_file = os.path.join(tempfile.gettempdir(), "vardictrust_convert_fastas.lock")
    with open(lock_file, "w") as lock_fd:
        try:
            fcntl.flock(lock_fd, fcntl.LOCK_EX)
        except OSError:
            pass  # If flock is not supported on target filesystem, continue

        t0 = time.time()
        converted = 0
        for fn in sorted(os.listdir(fasta_dir)):
            if not fn.endswith(".fa.csv"):
                continue
            csv_path = os.path.join(fasta_dir, fn)
            fa_name = fn[:-4]  # strip .csv -> .fa
            fa_path = os.path.join(fasta_dir, fa_name)
            fai_path = fa_path + ".fai"
            
            # Convert if missing or if csv is newer than either FASTA or FAI, or if FAI is older than FASTA
            if (
                not os.path.exists(fa_path)
                or not os.path.exists(fai_path)
                or os.path.getmtime(csv_path) > os.path.getmtime(fa_path)
                or os.path.getmtime(csv_path) > os.path.getmtime(fai_path)
                or os.path.getmtime(fai_path) < os.path.getmtime(fa_path)
            ):
                convert_csv_to_fasta(csv_path, fa_path, fai_path)
                converted += 1

        elapsed = round(time.time() - t0, 2)
        try:
            import subprocess
            excl = subprocess.check_output(
                ["git", "-C", fasta_dir, "rev-parse", "--git-path", "info/exclude"],
                text=True, stderr=subprocess.DEVNULL
            ).strip()
            if excl:
                excl_path = excl if os.path.isabs(excl) else os.path.join(fasta_dir, excl)
                os.makedirs(os.path.dirname(excl_path), exist_ok=True)
                existing = open(excl_path, "r", encoding="utf-8").read() if os.path.exists(excl_path) else ""
                with open(excl_path, "a", encoding="utf-8") as ef:
                    for pat in ("testdata/fastas/*.fa", "testdata/fastas/*.fa.fai"):
                        if pat not in existing:
                            ef.write(f"\n{pat}\n")
        except Exception:
            pass
        if converted > 0:
            print(f"convert_fastas: successfully materialized {converted} FASTAs in {elapsed}s")
        else:
            print(f"convert_fastas: all reference FASTAs up-to-date ({elapsed}s)")

if __name__ == "__main__":
    main()

