#!/usr/bin/env python3
"""Tier-2: compare two VCFs for VarDictRust <-> VarDictJava equivalence.
Exact: CHROM/POS/REF/ALT, FILTER, INFO key-set, FORMAT key-set, per-sample non-numeric values (incl GT).
Tolerant: numeric values within --eps after 5-dp HALF_EVEN. Records aligned by key -> order ignored.
Exit 0=MATCH, 1=MISMATCH."""
import sys, argparse
from decimal import Decimal, ROUND_HALF_EVEN

def r5(s):
    try: return float(Decimal(s).quantize(Decimal("0.00001"), rounding=ROUND_HALF_EVEN))
    except Exception: return None

def eq(a, b, eps):
    if a == b: return True
    ra, rb = r5(a), r5(b)
    return ra is not None and rb is not None and abs(ra - rb) <= eps

def kv(field):
    d = {}
    if field in (".", ""): return d
    for p in field.split(";"):
        k, sep, v = p.partition("="); d[k] = v if sep else True
    return d

def load(path):
    out = []
    with open(path) as fh:
        for ln in fh:
            if ln.strip() and not ln.startswith("#"):
                out.append(ln.rstrip("\n").split("\t"))
    return out

def keyof(r): return (r[0], r[1], r[3], r[4])  # CHROM, POS, REF, ALT

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("java"); ap.add_argument("rust"); ap.add_argument("--eps", type=float, default=1e-5)
    a = ap.parse_args()
    J = {keyof(r): r for r in load(a.java)}
    R = {keyof(r): r for r in load(a.rust)}
    if set(J) != set(R):
        oj = sorted(map(str, set(J) - set(R)))[:10]; orr = sorted(map(str, set(R) - set(J)))[:10]
        print(f"MISMATCH: variant sets differ. java-only(<=10)={oj} rust-only(<=10)={orr}"); sys.exit(1)
    errs = []
    for k in sorted(J):
        j, r = J[k], R[k]; w = f"{k[0]}:{k[1]} {k[2]}>{k[3]}"
        if j[6] != r[6]: errs.append(f"{w} FILTER: {j[6]!r} != {r[6]!r}")
        if not eq(j[5], r[5], a.eps): errs.append(f"{w} QUAL: {j[5]!r} != {r[5]!r}")
        ij, ir = kv(j[7]), kv(r[7])
        if set(ij) != set(ir):
            errs.append(f"{w} INFO keys: only-java={sorted(set(ij)-set(ir))} only-rust={sorted(set(ir)-set(ij))}")
        else:
            for key in ij:
                if ij[key] is True or ir[key] is True:
                    if ij[key] != ir[key]: errs.append(f"{w} INFO {key}: {ij[key]!r} != {ir[key]!r}")
                elif not eq(str(ij[key]), str(ir[key]), a.eps):
                    errs.append(f"{w} INFO {key}: {ij[key]!r} != {ir[key]!r}")
        if len(j) > 9 and len(r) > 9:
            if j[8] != r[8]:
                errs.append(f"{w} FORMAT keys: {j[8]!r} != {r[8]!r}")
            else:
                fkeys = j[8].split(":")
                for si in range(9, min(len(j), len(r))):
                    for key, xj, xr in zip(fkeys, j[si].split(":"), r[si].split(":")):
                        if not eq(xj, xr, a.eps): errs.append(f"{w} sample{si-9} {key}: {xj!r} != {xr!r}")
        if len(errs) >= 20: break
    if errs:
        print("MISMATCH:"); [print("  " + e) for e in errs[:20]]; sys.exit(1)
    print(f"MATCH ({len(J)} variants)"); sys.exit(0)

if __name__ == "__main__":
    main()
