#!/usr/bin/env python3
"""Tier-2: compare two VCFs for VarDictRust <-> VarDictJava equivalence.
Exact: CHROM/POS/REF/ALT, FILTER, INFO key-set, FORMAT key-set, per-sample non-numeric values (incl GT).
Tolerant: numeric values within --eps (raw float comparison). Records aligned by key -> order ignored.
Exit 0=MATCH, 1=MISMATCH."""
import sys, argparse, math
from collections import defaultdict

def eq(a, b, eps):
    if a == b: return True
    try:
        fa, fb = float(a), float(b)
        if math.isnan(fa) or math.isnan(fb): return False
        return abs(fa - fb) <= (eps + 1e-12)
    except (ValueError, TypeError):
        return False

def kv(field):
    d = {}
    if field in (".", ""): return d
    for p in field.split(";"):
        k, sep, v = p.partition("="); d[k] = v if sep else True
    return d

def load(path):
    out = []
    with open(path) as fh:
        for n, ln in enumerate(fh, 1):
            if ln.strip() and not ln.startswith("#"):
                cols = ln.rstrip("\n").split("\t")
                if len(cols) < 8:
                    print(f"MISMATCH: malformed VCF row {path}:{n}: {ln.strip()[:80]!r}"); sys.exit(1)
                out.append(cols)
    return out

def keyof(r): return (r[0], r[1], r[3], r[4])  # CHROM, POS, REF, ALT

def group(rows):
    g = defaultdict(list)
    for r in rows: g[keyof(r)].append(r)
    for k in g: g[k].sort(key="\t".join)
    return g

def compare_pair(k, j, r, eps, errs):
    w = f"{k[0]}:{k[1]} {k[2]}>{k[3]}"
    if j[6] != r[6]: errs.append(f"{w} FILTER: {j[6]!r} != {r[6]!r}")
    if not eq(j[5], r[5], eps): errs.append(f"{w} QUAL: {j[5]!r} != {r[5]!r}")
    ij, ir = kv(j[7]), kv(r[7])
    if set(ij) != set(ir):
        errs.append(f"{w} INFO keys: only-java={sorted(set(ij)-set(ir))} only-rust={sorted(set(ir)-set(ij))}")
    else:
        for key in ij:
            if ij[key] is True or ir[key] is True:
                if ij[key] != ir[key]: errs.append(f"{w} INFO {key}: {ij[key]!r} != {ir[key]!r}")
            elif not eq(str(ij[key]), str(ir[key]), eps):
                errs.append(f"{w} INFO {key}: {ij[key]!r} != {ir[key]!r}")
    if len(j) > 9 and len(r) > 9:
        if j[8] != r[8]: errs.append(f"{w} FORMAT keys: {j[8]!r} != {r[8]!r}")
        else:
            fkeys = j[8].split(":")
            for si in range(9, min(len(j), len(r))):
                for key, xj, xr in zip(fkeys, j[si].split(":"), r[si].split(":")):
                    if not eq(xj, xr, eps): errs.append(f"{w} sample{si-9} {key}: {xj!r} != {xr!r}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("java"); ap.add_argument("rust"); ap.add_argument("--eps", type=float, default=1e-6)
    ap.add_argument("--min-variants", type=int, default=0)
    a = ap.parse_args()
    LJ, LR = load(a.java), load(a.rust)
    J, R = group(LJ), group(LR)
    if set(J) != set(R):
        oj = sorted(map(str, set(J) - set(R)))[:10]; orr = sorted(map(str, set(R) - set(J)))[:10]
        print(f"MISMATCH: variant sets differ. java-only(<=10)={oj} rust-only(<=10)={orr}"); sys.exit(1)
    bad = sorted(k for k in J if len(J[k]) != len(R[k]))
    if bad:
        print(f"MISMATCH: record counts differ for {len(bad)} key(s) (<=10): "
              f"{[(str(k), len(J[k]), len(R[k])) for k in bad[:10]]}"); sys.exit(1)
    if len(LJ) < a.min_variants:
        print(f"MISMATCH: only {len(LJ)} records compared (< --min-variants {a.min_variants})"); sys.exit(1)
    errs = []
    for k in sorted(J):
        for j, r in zip(J[k], R[k]): compare_pair(k, j, r, a.eps, errs)
        if len(errs) >= 20: break
    if errs:
        print("MISMATCH:"); [print("  " + e) for e in errs[:20]]; sys.exit(1)
    if not LJ: print("WARNING: vacuous comparison (0 records)", file=sys.stderr)
    print(f"MATCH ({len(LJ)} variants)"); sys.exit(0)

if __name__ == "__main__":
    main()
