#!/usr/bin/env python3
"""Convert upstream VarDictJava testdata thin-FASTA CSVs to indexed FASTAs.

Upstream VarDictJava (AstraZeneca-NGS) only commits sparse reference sequence
chunks as CSVs (`VarDictJava/testdata/fastas/*.fa.csv`) to save repository space.
Their integration test harness parses these CSVs on the fly using Java classes.
VarDictRust and samtools require true FASTA files with companion `.fa.fai` indices.

This script reconstructs bitwise-identical FASTA files and companion `.fai` indices
directly in pure Python (zero external dependencies, runtime <10s).
"""

import argparse
import errno
import fcntl
import os
import stat
import subprocess
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

def _needs_conversion(csv_path: str, fa_path: str, fai_path: str) -> bool:
    try:
        csv_mtime = os.path.getmtime(csv_path)
        fa_mtime = os.path.getmtime(fa_path)
        fai_mtime = os.path.getmtime(fai_path)
    except OSError:
        return True
    return csv_mtime > fa_mtime or csv_mtime > fai_mtime or fai_mtime < fa_mtime


def _ensure_git_excludes(fasta_dir: str) -> None:
    try:
        excl = subprocess.check_output(
            ["git", "-C", fasta_dir, "rev-parse", "--git-path", "info/exclude"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
        if not excl:
            return
        excl_path = excl if os.path.isabs(excl) else os.path.join(fasta_dir, excl)
        os.makedirs(os.path.dirname(excl_path), exist_ok=True)
        existing = ""
        if os.path.exists(excl_path):
            with open(excl_path, "r", encoding="utf-8") as rf:
                existing = rf.read()
        missing = [
            pat
            for pat in (
                "testdata/fastas/*.fa",
                "testdata/fastas/*.fa.fai",
                "testdata/fastas/.convert_fastas.lock",
            )
            if pat not in existing
        ]
        if missing:
            with open(excl_path, "a", encoding="utf-8") as ef:
                for pat in missing:
                    ef.write(f"\n{pat}\n")
    except Exception:
        pass


def _open_lock_fd(fasta_dir: str) -> int:
    flags = os.O_CREAT | os.O_RDWR | getattr(os, "O_NOFOLLOW", 0)
    primary_lock = os.path.join(fasta_dir, ".convert_fastas.lock")
    try:
        fd = os.open(primary_lock, flags, 0o600)
        try:
            st = os.fstat(fd)
            if not stat.S_ISREG(st.st_mode):
                raise PermissionError(f"Lock file is not a regular file: {primary_lock}")
        except BaseException:
            os.close(fd)
            raise
        return fd
    except OSError as exc:
        if exc.errno not in (errno.EROFS, errno.EACCES, errno.EPERM):
            raise

    runtime_dir = os.environ.get("XDG_RUNTIME_DIR") or "/tmp"
    fallback_lock = os.path.join(
        runtime_dir, f"vardictrust-convert-fastas-{os.getuid()}.lock"
    )
    fd = os.open(fallback_lock, flags, 0o600)
    try:
        st = os.fstat(fd)
        if not stat.S_ISREG(st.st_mode) or st.st_uid != os.getuid():
            raise PermissionError(
                f"Unsafe fallback lock file ownership or mode: {fallback_lock}"
            )
    except BaseException:
        os.close(fd)
        raise
    return fd


def main():
    default_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    parser = argparse.ArgumentParser(description="Convert upstream VarDictJava thin-FASTA CSVs to indexed FASTAs.")
    parser.add_argument("--data-root", default=None, help="Root directory containing VarDictJava/testdata/fastas")
    args = parser.parse_args()
    if args.data_root is not None:
        root = os.path.abspath(args.data_root)
    else:
        root = default_root
        if not os.path.isdir(os.path.join(root, "VarDictJava", "testdata", "fastas")):
            try:
                common_dir = subprocess.check_output(
                    ["git", "-C", root, "rev-parse", "--path-format=absolute", "--git-common-dir"],
                    text=True,
                    stderr=subprocess.DEVNULL,
                ).strip()
                candidate = os.path.dirname(common_dir)
                if os.path.isdir(os.path.join(candidate, "VarDictJava", "testdata", "fastas")):
                    root = candidate
            except Exception:
                pass
    fasta_dir = os.path.join(root, "VarDictJava", "testdata", "fastas")
    if not os.path.isdir(fasta_dir):
        print(f"convert_fastas: directory not found: {fasta_dir}", file=sys.stderr)
        sys.exit(1)

    t0 = time.time()
    csv_entries = []
    for fn in sorted(os.listdir(fasta_dir)):
        if not fn.endswith(".fa.csv"):
            continue
        csv_path = os.path.join(fasta_dir, fn)
        fa_path = os.path.join(fasta_dir, fn[:-4])  # strip .csv -> .fa
        fai_path = fa_path + ".fai"
        csv_entries.append((csv_path, fa_path, fai_path))

    if not any(_needs_conversion(c, fa, fai) for c, fa, fai in csv_entries):
        _ensure_git_excludes(fasta_dir)
        elapsed = round(time.time() - t0, 2)
        print(f"convert_fastas: all reference FASTAs up-to-date ({elapsed}s)")
        return

    _ensure_git_excludes(fasta_dir)
    lock_fd = _open_lock_fd(fasta_dir)
    try:
        try:
            fcntl.flock(lock_fd, fcntl.LOCK_EX)
        except OSError:
            pass  # If flock is not supported on target filesystem, continue

        converted = 0
        for csv_path, fa_path, fai_path in csv_entries:
            if _needs_conversion(csv_path, fa_path, fai_path):
                convert_csv_to_fasta(csv_path, fa_path, fai_path)
                converted += 1

        elapsed = round(time.time() - t0, 2)
        _ensure_git_excludes(fasta_dir)
        if converted > 0:
            print(f"convert_fastas: successfully materialized {converted} FASTAs in {elapsed}s")
        else:
            print(f"convert_fastas: all reference FASTAs up-to-date ({elapsed}s)")
    finally:
        os.close(lock_fd)

if __name__ == "__main__":
    main()

