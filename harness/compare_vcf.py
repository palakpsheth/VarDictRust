#!/usr/bin/env python3
"""Tier-2: compare two VCFs for VarDictRust <-> VarDictJava equivalence.
Exact: CHROM/POS/REF/ALT, FILTER, INFO key-set, FORMAT key-set, per-sample non-numeric values (incl GT).
Tolerant: numeric values within --eps (raw float comparison). Records aligned by key -> order ignored.
Exit 0=MATCH, 1=MISMATCH."""
import sys, argparse, math, re
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

def load_header(path):
    hdr = []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#"):
                hdr.append(ln.rstrip("\r\n"))
            elif ln.strip():
                break
    return hdr

def compare_headers(j_hdr, r_hdr, errs):
    vol_pat = re.compile(r"^##(fileDate|source)=")
    fj = [l.rstrip("\r\n") for l in j_hdr if not vol_pat.match(l.strip())]
    fr = [l.rstrip("\r\n") for l in r_hdr if not vol_pat.match(l.strip())]
    if fj != fr:
        for i in range(min(len(fj), len(fr))):
            if fj[i] != fr[i]:
                errs.append(f"HEADER line {i+1}: java={fj[i]!r} != rust={fr[i]!r}")
                if len(errs) >= 20:
                    return
        if len(fj) != len(fr):
            extra = "java" if len(fj) > len(fr) else "rust"
            first_extra = fj[min(len(fj), len(fr))] if len(fj) > len(fr) else fr[min(len(fj), len(fr))]
            errs.append(f"HEADER line count: java={len(fj)} != rust={len(fr)} (first extra {extra}: {first_extra!r})")

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
    for v in g.values():
        if len(v) > 1: v.sort(key="\t".join)
    return g

def compare_pair(k, j, r, eps, errs, strict=None):
    if j == r: return
    w = f"{k[0]}:{k[1]} {k[2]}>{k[3]}"
    if j[6] != r[6]: errs.append(f"{w} FILTER: {j[6]!r} != {r[6]!r}")
    if not eq(j[5], r[5], eps): errs.append(f"{w} QUAL: {j[5]!r} != {r[5]!r}")
    elif strict is not None and j[5] != r[5]: strict.append(f"{w} QUAL: {j[5]!r} vs {r[5]!r}")
    if j[7] != r[7]:
        ij, ir = kv(j[7]), kv(r[7])
        if set(ij) != set(ir):
            errs.append(f"{w} INFO keys: only-java={sorted(set(ij)-set(ir))} only-rust={sorted(set(ir)-set(ij))}")
        else:
            for key in ij:
                if ij[key] is True or ir[key] is True:
                    if ij[key] != ir[key]: errs.append(f"{w} INFO {key}: {ij[key]!r} != {ir[key]!r}")
                elif not eq(str(ij[key]), str(ir[key]), eps):
                    errs.append(f"{w} INFO {key}: {ij[key]!r} != {ir[key]!r}")
                elif strict is not None and str(ij[key]) != str(ir[key]):
                    strict.append(f"{w} INFO {key}: {ij[key]!r} vs {ir[key]!r}")
    if len(j) > 9 and len(r) > 9 and j[8:] != r[8:]:
        if j[8] != r[8]: errs.append(f"{w} FORMAT keys: {j[8]!r} != {r[8]!r}")
        else:
            fkeys = j[8].split(":")
            for si in range(9, min(len(j), len(r))):
                if j[si] == r[si]: continue
                for key, xj, xr in zip(fkeys, j[si].split(":"), r[si].split(":")):
                    if not eq(xj, xr, eps): errs.append(f"{w} sample{si-9} {key}: {xj!r} != {xr!r}")
                    elif strict is not None and xj != xr: strict.append(f"{w} sample{si-9} {key}: {xj!r} vs {xr!r}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("java"); ap.add_argument("rust"); ap.add_argument("--eps", type=float, default=1e-6)
    ap.add_argument("--min-variants", type=int, default=0)
    ap.add_argument("--strict-bytes", action="store_true")
    a = ap.parse_args()
    JH, RH = load_header(a.java), load_header(a.rust)
    LJ, LR = load(a.java), load(a.rust)
    J, R = group(LJ), group(LR)
    errs = []
    strict = [] if a.strict_bytes else None
    compare_headers(JH, RH, errs)
    if set(J) != set(R):
        oj = sorted(map(str, set(J) - set(R)))[:10]; orr = sorted(map(str, set(R) - set(J)))[:10]
        print(f"MISMATCH: variant sets differ. java-only(<=10)={oj} rust-only(<=10)={orr}"); sys.exit(1)
    bad = sorted(k for k in J if len(J[k]) != len(R[k]))
    if bad:
        print(f"MISMATCH: record counts differ for {len(bad)} key(s) (<=10): "
              f"{[(str(k), len(J[k]), len(R[k])) for k in bad[:10]]}"); sys.exit(1)
    if len(LJ) < a.min_variants:
        print(f"MISMATCH: only {len(LJ)} records compared (< --min-variants {a.min_variants})"); sys.exit(1)
    for k in sorted(J):
        for j, r in zip(J[k], R[k]): compare_pair(k, j, r, a.eps, errs, strict)
        if len(errs) >= 20: break
    if errs:
        print("MISMATCH:"); [print("  " + e) for e in errs[:20]]; sys.exit(1)
    if not LJ: print("WARNING: vacuous comparison (0 records)", file=sys.stderr)
    if strict is not None:
        print(f"STRICT-BYTES: {len(strict)} textual difference(s) within eps")
        for e in strict[:20]: print("  " + e)
    print(f"MATCH ({len(LJ)} variants)"); sys.exit(0)

if __name__ == "__main__":
    main()
