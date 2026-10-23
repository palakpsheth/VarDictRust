#!/usr/bin/env python3
"""harness/profile_catalog.py: rank the hot spots in one or more profile cell dirs into a Markdown catalog.

Usage: profile_catalog.py <cell-dir>... -o catalog.md [--min-pct 2]

Reads whatever of these exists in each cell dir (graceful partial output; stdlib only):
  cpu-collapsed.txt + cpu-meta.json   pprof CPU stacks `frame;frame;... N`, root->leaf (profile-cpu build)
  dhat-heap.json                      dhat-rs 0.3.3 heap profile (profile-dhat build)
  heap.pb.gz                          gzipped pprof protobuf from jemalloc_pprof (profile-heap build)
and writes beside them: dhat-bytes-collapsed.txt / dhat-blocks-collapsed.txt / heap-collapsed.txt (full stacks,
root->leaf) and, if `inferno-flamegraph` is on PATH, dhat-bytes.svg / dhat-blocks.svg.

Attribution rule (the same for CPU, dhat and heap): a sample/allocation is charged to its INNERMOST frame whose
symbol contains `vardictrust_core::` (falling back to `vardictrust::`), never to the allocator/std leaf. dhat `fs`
lists are innermost-first (ftbl[0] = "[root]" never appears in them); CPU/heap collapsed stacks are root->leaf.
"""
import argparse
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
from collections import OrderedDict

HASH = re.compile(r"::h[0-9a-f]{16}$")
ADDR = re.compile(r"^0x[0-9a-f]+: ")
SRCLOC = re.compile(r" \([^()]*:\d+(?::\d+)?\)$")
NO_VD = "(no vardictrust frame)"
TOP_STACKS = 3
STACK_DEPTH = 8


def clean(sym):
    """dhat `0xADDR: symbol (file:line:col)` / `symbol::h<16 hex>` -> bare symbol (so one function is one row)."""
    return HASH.sub("", SRCLOC.sub("", ADDR.sub("", sym.strip())))


def vd_frame(frames):
    """Innermost vardictrust frame; `frames` is innermost-first."""
    for marker in ("vardictrust_core::", "vardictrust::"):
        for f in frames:
            if marker in f:
                return f
    return NO_VD


def md(sym):
    """frame as an inline code span: symbols hold `<..>`, `|` and `*` that Markdown would otherwise eat."""
    return "`%s`" % sym.replace("|", "\\|") if not sym.startswith("(") and not sym.startswith("other (") else sym


def pct(n, total):
    return "%.1f%%" % (100.0 * n / total) if total else "0.0%"


def table(rows, total, min_pct, count_col, pct_col, extra_cols=()):
    """rows: [(name, count, *extra)] -> markdown, sorted by count desc; rows under min_pct folded into 'other'."""
    rows = sorted(rows, key=lambda r: (-r[1], r[0]))
    keep = [r for r in rows if total and 100.0 * r[1] / total >= min_pct]
    rest = [r for r in rows if r not in keep]
    if rest:
        keep.append(tuple(["other (%d rows)" % len(rest)] + [sum(r[i] for r in rest) for i in range(1, len(rest[0]))]))
    heads = ["rank", "frame", count_col, pct_col] + list(extra_cols)
    out = ["| " + " | ".join(heads) + " |", "|---:|---|" + "---:|" * (len(heads) - 2)]
    for i, r in enumerate(keep, 1):
        cells = [str(i), md(r[0]), str(r[1]), pct(r[1], total)] + [str(x) for x in r[2:]]
        out.append("| " + " | ".join(cells) + " |")
    return out


def write_collapsed(path, stacks):
    """stacks: {stack-string: weight}; zero weights dropped, sorted."""
    lines = sorted("%s %d" % (k, v) for k, v in stacks.items() if v > 0)
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + ("\n" if lines else ""))


# ---- CPU ---------------------------------------------------------------------------------------
def cpu_section(d, min_pct):
    path = os.path.join(d, "cpu-collapsed.txt")
    if not os.path.exists(path):
        return None
    leaf, attr, total = {}, {}, 0
    with open(path) as fh:
        for ln in fh:
            ln = ln.rstrip("\n")
            if not ln:
                continue
            stack, _, n = ln.rpartition(" ")
            frames = [clean(f) for f in stack.split(";")]
            n = int(n)
            total += n
            leaf[frames[-1]] = leaf.get(frames[-1], 0) + n
            a = vd_frame(frames[::-1])
            attr[a] = attr.get(a, 0) + n
    out = ["### CPU (pprof, 499 Hz)", ""]
    meta = os.path.join(d, "cpu-meta.json")
    if os.path.exists(meta):
        with open(meta) as fh:
            m = json.load(fh)
        t, k = m["ticks_total"], m["samples_kept"]
        out += ["pprof kept %d of %d SIGPROF ticks (%s kept, %s dropped; the drop includes the address blocklist "
                "and pprof `try_write` contention drops)." % (k, t, pct(k, t), pct(t - k, t)), ""]
    out += ["Percentages are of the %d kept samples." % total, "", "**Self time by leaf frame**", ""]
    out += table(list(leaf.items()), total, min_pct, "samples", "% CPU")
    out += ["", "**Time under the innermost `vardictrust` frame (inclusive of non-vardict callees)**", ""]
    out += table(list(attr.items()), total, min_pct, "samples", "% CPU")
    return out


# ---- dhat --------------------------------------------------------------------------------------
def dhat_section(d, min_pct):
    path = os.path.join(d, "dhat-heap.json")
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        doc = json.load(fh)
    ftbl = doc.get("ftbl", [])
    per = OrderedDict()   # attributed frame -> [tb, tbk, mb, [(tb, frames)]]
    by_bytes, by_blocks = {}, {}
    for pp in doc.get("pps", []):
        frames = [clean(ftbl[i]) for i in pp.get("fs", []) if 0 <= i < len(ftbl)]   # innermost-first
        tb, tbk, mb = pp.get("tb", 0), pp.get("tbk", 0), pp.get("mb", 0)
        stack = ";".join(reversed(frames))
        by_bytes[stack] = by_bytes.get(stack, 0) + tb
        by_blocks[stack] = by_blocks.get(stack, 0) + tbk
        row = per.setdefault(vd_frame(frames), [0, 0, 0, []])
        row[0] += tb
        row[1] += tbk
        row[2] += mb
        row[3].append((tb, frames))
    write_collapsed(os.path.join(d, "dhat-bytes-collapsed.txt"), by_bytes)
    write_collapsed(os.path.join(d, "dhat-blocks-collapsed.txt"), by_blocks)
    tot_b = sum(r[0] for r in per.values())
    tot_k = sum(r[1] for r in per.values())
    out = ["### dhat (allocation sites)", "",
           "total_bytes %d, total_blocks %d. `max_bytes` is the sum of per-site peaks (an upper bound on any "
           "frame's true peak). Percentages are of total_bytes." % (tot_b, tot_k), ""]
    out += table([(k, v[0], v[1], v[2]) for k, v in per.items()], tot_b, min_pct, "total_bytes", "% bytes",
                 ["total_blocks", "max_bytes"])
    out += ["", "#### Top stacks (innermost frame first, by total_bytes)", ""]
    for name, v in sorted(per.items(), key=lambda kv: (-kv[1][0], kv[0])):
        if tot_b and 100.0 * v[0] / tot_b < min_pct:
            continue
        out.append("- `%s`" % name)
        for tb, frames in sorted(v[3], key=lambda x: -x[0])[:TOP_STACKS]:
            out.append("  - %d B: `%s`" % (tb, " ← ".join(frames[:STACK_DEPTH])))
    return out


# ---- heap (pprof profile.proto subset, hand-decoded) -----------------------------------------------
def _varint(b, i):
    r = s = 0
    while True:
        c = b[i]
        i += 1
        r |= (c & 0x7F) << s
        s += 7
        if c < 0x80:
            return r, i


def _fields(b):
    """yield (field number, wire type, value) of one protobuf message."""
    i = 0
    while i < len(b):
        k, i = _varint(b, i)
        w = k & 7
        if w == 0:
            v, i = _varint(b, i)
        elif w == 2:
            n, i = _varint(b, i)
            v, i = b[i:i + n], i + n
        elif w == 1:
            v, i = b[i:i + 8], i + 8
        elif w == 5:
            v, i = b[i:i + 4], i + 4
        else:
            raise ValueError("unsupported protobuf wire type %d" % w)
        yield k >> 3, w, v


def _ints(w, v):
    """repeated varint field, packed (wire 2) or not (wire 0)"""
    if w == 0:
        return [v]
    out, i = [], 0
    while i < len(v):
        x, i = _varint(v, i)
        out.append(x)
    return out


def decode_pprof(raw):
    """-> (value type names, samples [(loc ids leaf-first, values)], {loc id: (address, [function ids inner-first])},
    {function id: name}). Ports only the profile.proto subset the catalog needs."""
    strings, types, samples, locs, funcs = [], [], [], {}, {}
    for f, w, v in _fields(raw):
        if f == 6:
            strings.append(v.decode("utf-8", "replace"))
        elif f == 1:
            types.append(dict((ff, vv) for ff, _, vv in _fields(v)))
        elif f == 2:
            ids, vals = [], []
            for ff, ww, vv in _fields(v):
                if ff == 1:
                    ids += _ints(ww, vv)
                elif ff == 2:
                    vals += _ints(ww, vv)
            samples.append((ids, vals))
        elif f == 4:
            lid, addr, lines = 0, 0, []
            for ff, ww, vv in _fields(v):
                if ff == 1:
                    lid = vv
                elif ff == 3:
                    addr = vv
                elif ff == 4:
                    lines += [x for fff, _, x in _fields(vv) if fff == 1]
            locs[lid] = (addr, lines)
        elif f == 5:
            d = dict((ff, vv) for ff, _, vv in _fields(v))
            funcs[d.get(1, 0)] = d.get(2, 0)
    name = lambda i: strings[i] if 0 <= i < len(strings) else ""
    return [(name(t.get(1, 0)), name(t.get(2, 0))) for t in types], samples, locs, \
        {k: name(v) for k, v in funcs.items()}


def pick_value_index(types):
    """Weight = live heap bytes: 'inuse_space', else 'alloc_space', else the first 'bytes'-unit type, else the
    last type (the pprof convention puts the default view last)."""
    for want in ("inuse_space", "alloc_space"):
        for i, (t, _) in enumerate(types):
            if t == want:
                return i
    for i, (_, u) in enumerate(types):
        if u == "bytes":
            return i
    return len(types) - 1


def heap_section(d, min_pct):
    path = os.path.join(d, "heap.pb.gz")
    if not os.path.exists(path):
        return None
    with open(path, "rb") as fh:
        raw = fh.read()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    types, samples, locs, funcs = decode_pprof(raw)
    vi = pick_value_index(types)
    named = False
    stacks, per, total = {}, {}, 0
    for ids, vals in samples:
        w = vals[vi] if 0 <= vi < len(vals) else 0
        frames = []   # innermost-first: location order is leaf-first, each location's lines inner-first
        for lid in ids:
            addr, fids = locs.get(lid, (0, []))
            names = [funcs.get(x, "") for x in fids]
            if any(names):
                named = True
            frames += [clean(n) for n in names if n] or ["0x%x" % addr]
        key = ";".join(reversed(frames))
        stacks[key] = stacks.get(key, 0) + w
        a = vd_frame(frames)
        per[a] = per.get(a, 0) + w
        total += w
    write_collapsed(os.path.join(d, "heap-collapsed.txt"), stacks)
    out = ["### Heap (jemalloc_pprof, sample type `%s`)" % (types[vi][0] if types else "?"), "",
           "Weight is live in-use bytes at the end of the run (a snapshot taken at exit, not the peak). "
           "Total %d bytes over %d sampled stacks." % (total, len(samples)), ""]
    if not named:
        out += ["> no function names in this profile (addresses only); frames are raw addresses.", ""]
    out += table(list(per.items()), total, min_pct, "bytes", "% bytes")
    return out


# ---- svg hook -------------------------------------------------------------------------------------
_hinted = False


def render_dhat_svgs(d):
    global _hinted
    exe = shutil.which("inferno-flamegraph")
    if not exe:
        if not _hinted:
            sys.stderr.write("profile_catalog: inferno-flamegraph not on PATH; skipping dhat SVGs "
                             "(cargo install inferno --locked)\n")
            _hinted = True
        return
    for kind, count in (("bytes", "bytes"), ("blocks", "blocks")):
        src = os.path.join(d, "dhat-%s-collapsed.txt" % kind)
        dst = os.path.join(d, "dhat-%s.svg" % kind)
        if not os.path.exists(src):
            continue
        with open(src) as i, open(dst, "w") as o:
            r = subprocess.run([exe, "--title", "dhat %s (%s)" % (kind, os.path.basename(os.path.abspath(d))),
                                "--countname", count], stdin=i, stdout=o, stderr=subprocess.PIPE)
        if r.returncode != 0:
            sys.stderr.write("profile_catalog: inferno-flamegraph failed for %s: %s\n" %
                             (src, r.stderr.decode("utf-8", "replace").strip()))
            os.remove(dst)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("dirs", nargs="+")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--min-pct", type=float, default=2.0, help="fold rows below this %% into 'other' (default 2)")
    a = ap.parse_args()
    out_dir = os.path.dirname(os.path.abspath(a.out))
    doc = ["# Profile catalog", "",
           "Attribution: innermost `vardictrust_core::` frame (fallback `vardictrust::`); rows under %g%% are "
           "folded into `other`." % a.min_pct, ""]
    for d in a.dirs:
        rel = os.path.relpath(os.path.abspath(d), out_dir)
        doc += ["## " + (rel if rel != "." else os.path.basename(os.path.abspath(d))), ""]
        sections = [cpu_section(d, a.min_pct), dhat_section(d, a.min_pct), heap_section(d, a.min_pct)]
        if os.path.exists(os.path.join(d, "dhat-heap.json")):
            render_dhat_svgs(d)
        sections = [s for s in sections if s]
        if not sections:
            doc += ["_no profile data in this dir_", ""]
        for s in sections:
            doc += s + [""]
        svgs = sorted(f for f in os.listdir(d) if f.endswith(".svg"))
        if svgs:
            doc += ["### Flamegraphs", ""]
            doc += ["- [%s](%s)" % (f, os.path.normpath(os.path.join(rel, f))) for f in svgs] + [""]
    with open(a.out, "w") as fh:
        fh.write("\n".join(doc))


if __name__ == "__main__":
    main()
