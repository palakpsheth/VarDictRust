#!/usr/bin/env python3
"""Pick dense truth-variant clusters as slice regions (deterministic).
Output: 4-col BED (0-based, half-open) `chrom start end slice_<i>_<count>`, merged, sorted."""
import argparse, bisect, gzip, sys

def merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]: out[-1] = (out[-1][0], max(out[-1][1], e))
        else: out.append((s, e))
    return out

def read_bed(path):
    d = {}
    with open(path) as fh:
        for ln in fh:
            if not ln.strip() or ln.startswith(("#", "track", "browser")): continue
            c, s, e = ln.rstrip("\n").split("\t")[:3]; d.setdefault(c, []).append((int(s), int(e)))
    return {c: merge(v) for c, v in d.items()}

def inside(iv, c, p0):
    lst = iv.get(c, []); i = bisect.bisect_right(lst, (p0, float("inf"))) - 1
    return i >= 0 and lst[i][0] <= p0 < lst[i][1]

def vcf_sites(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for ln in fh:
            if ln.startswith("#"): continue
            f = ln.split("\t", 7)
            if len(f) >= 7 and not (f[6] == "." or "PASS" in f[6].split(";")): continue
            yield f[0], int(f[1]) - 1, len(f[3])

def clusters(sites, chain, max_span):
    cur = []
    for c, p, rl in sites:
        if cur and (c != cur[0][0] or p - cur[-1][1] > chain or p + rl - cur[0][1] > max_span):
            yield cur; cur = []
        cur.append((c, p, rl))
    if cur: yield cur

def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", action="append", required=True)
    ap.add_argument("--eval", "--eval-bed", dest="eval", required=True)
    ap.add_argument("--exclude", "--exclude-bed", dest="exclude")
    ap.add_argument("--chain", "--max-gap", dest="chain", type=int, required=True)
    ap.add_argument("--min-count", type=int, required=True)
    ap.add_argument("--pad", type=int, default=1000)
    ap.add_argument("--max-span", type=int, required=True)
    ap.add_argument("--max-regions", type=int, required=True)
    ap.add_argument("--out")
    a = ap.parse_args(argv)
    ev = read_bed(a.eval); ex = read_bed(a.exclude) if a.exclude else {}
    sites = sorted({s for t in a.truth for s in vcf_sites(t) if inside(ev, s[0], s[1]) and not inside(ex, s[0], s[1])})
    cl = [k for k in clusters(sites, a.chain, a.max_span) if len(k) >= a.min_count]
    if not cl: sys.exit("select_slice_regions: no cluster meets --min-count; change parameters deliberately and record why")
    cl.sort(key=lambda k: (-len(k), k[0][0], k[0][1]))
    regs = {}
    for k in cl[:a.max_regions]:
        c = k[0][0]; regs.setdefault(c, []).append((max(0, k[0][1] - a.pad), max(p + rl for _, p, rl in k) + a.pad))
    n = 0
    out_fh = open(a.out, "w") if a.out else sys.stdout
    try:
        for c in sorted(regs):
            for s, e in merge(regs[c]):
                n += 1; cnt = sum(1 for cc, p, _ in sites if cc == c and s <= p < e)
                print(f"{c}\t{s}\t{e}\tslice_{n}_{cnt}", file=out_fh)
    finally:
        if a.out and out_fh is not sys.stdout:
            out_fh.close()

if __name__ == "__main__": main()
