#!/usr/bin/env python3
"""Tier-1: compare two VarDict intermediate .var tables (TSV, no header).
Field-by-field, row order preserved: numeric within --eps (raw float comparison), else exact.
--cols A-B restricts to the 1-based inclusive column range A..B. Exit 0=MATCH, 1=MISMATCH."""
import sys, argparse, math

def eq(a, b, eps):
    if a == b: return True
    try:
        fa, fb = float(a), float(b)
        if math.isnan(fa) or math.isnan(fb): return False
        return abs(fa - fb) <= (eps + 1e-12)
    except (ValueError, TypeError):
        return False

def load(p):
    with open(p) as f:
        return [ln.rstrip("\n").split("\t") for ln in f if ln.strip()]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("java"); ap.add_argument("rust"); ap.add_argument("--eps", type=float, default=1e-6)
    ap.add_argument("--cols", default=None, help="1-based inclusive range, e.g. 8-13")
    a = ap.parse_args()
    J, R = load(a.java), load(a.rust)
    if len(J) != len(R):
        print(f"MISMATCH: row count java={len(J)} rust={len(R)}"); sys.exit(1)
    lo, hi = 0, None
    if a.cols:
        s, e = a.cols.split("-"); lo, hi = int(s) - 1, int(e)
    errs = []
    for i, (j, r) in enumerate(zip(J, R)):
        if len(j) != len(r):
            errs.append(f"row {i}: field count java={len(j)} rust={len(r)}")
        else:
            top = min(hi, len(j)) if hi else len(j)
            for c in range(lo, top):
                if not eq(j[c], r[c], a.eps):
                    errs.append(f"row {i} col {c+1}: {j[c]!r} != {r[c]!r}")
        if len(errs) >= 20: break
    if errs:
        print("MISMATCH:"); [print("  " + e) for e in errs[:20]]; sys.exit(1)
    print(f"MATCH ({len(J)} rows)"); sys.exit(0)

if __name__ == "__main__":
    main()
