#!/usr/bin/env python3
"""Pinned golden baseline store for java/cpp oracle outputs (AGENTS.md §14, issue #31).

Usage (always via the shim: bash harness/golden_store.sh <verb> ...):
  oracle LABEL REF BAM[|BAM] BED OUT [SAMPLE [FLAGS...]]   run_oracle.sh argv + LABEL
  cpp    LABEL REF BAM[|BAM] BED OUT [SAMPLE [FLAGS...]]   run_cpp.sh   argv + LABEL
  raw    --kind java|cpp --label L --out F [--err-out F] [--input P]... -- CMD...
  verify [--quiet]            re-hash every live artifact against its manifest
  stats  [--run ID] [--min-hit R] [--json]
  list                        one line per cell
  tombstone KEY REASON        retire a cell (append-only; the next run re-pins)

Hit = same oracle kind + pin + conda-env sha + wrapper sha + input sha256s + argv (+ AF_THR,
LC_ALL, LANG, CPP_CHUNK): the pinned output, its recorded stderr and its `time -v` file are
replayed and no engine runs. A `time -v` file is only replayed on the host that measured it
(AGENTS.md §14). Miss = run live, pin on exit 0 (never when the output is empty/implausible,
an input changed during the run, or the store cannot be written — the run result still stands).
Env: GOLDEN_STORE_DIR (default scratch/golden-store), GOLDEN_MODE=auto|regen|off,
GOLDEN_TIME_OUT=<file>, GOLDEN_RUN_ID, GOLDEN_REHASH=1, GOLDEN_ORACLE_CMD / GOLDEN_CPP_CMD
(test seams; folded into the key). Python stdlib only.
"""
import contextlib
import datetime
import gzip
import hashlib
import json
import os
import pathlib
import platform
import shlex
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import zlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
SCHEMA = 2
GZ_MIN_BYTES = 1_000_000
TIME_BIN = "/usr/bin/time"
LOCK_STALE_S = 600
TMP_GC_S = 86400
ENV_KEYS = {"java": ("AF_THR", "LC_ALL", "LANG"), "cpp": ("AF_THR", "CPP_CHUNK", "LC_ALL", "LANG")}
EVENTS_LIVE = ("MISS", "CORRUPT", "REGEN", "FAIL", "REFUSED", "UNSTABLE", "PINFAIL")


def store_dir():
    return pathlib.Path(os.environ.get("GOLDEN_STORE_DIR") or ROOT / "scratch" / "golden-store")


def mode():
    m = os.environ.get("GOLDEN_MODE", "auto")
    if m not in ("auto", "regen", "off"):
        sys.exit(f"golden_store: GOLDEN_MODE must be auto|regen|off, got {m!r}")
    return m


def now_utc():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def atomic_write(path, data):
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=".tmp-")
    with os.fdopen(fd, "wb") as f:
        f.write(data)
    os.replace(tmp, path)


# --- input identity -------------------------------------------------------------

def _hash(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_file(path):
    """sha256 of a file, cached by (resolved path, size, mtime_ns, ctime_ns); GOLDEN_REHASH=1 skips the cache."""
    p = pathlib.Path(path)
    if not p.is_file():
        sys.exit(f"golden_store: input not found: {path}")
    p = p.resolve()
    st = p.stat()
    cache = store_dir() / ".shacache.json"
    try:
        data = json.loads(cache.read_text())
    except (OSError, ValueError):
        data = {}
    k = f"{p}|{st.st_size}|{st.st_mtime_ns}|{st.st_ctime_ns}"
    if os.environ.get("GOLDEN_REHASH") != "1" and k in data:
        return data[k]
    data[k] = digest = _hash(p)
    if len(data) > 5000:  # temp BEDs add an entry per run: drop entries whose file is gone
        data = {kk: v for kk, v in data.items() if os.path.exists(kk.split("|")[0])}
    atomic_write(cache, json.dumps(data, sort_keys=True).encode())
    return digest


def oracle_pin(kind):
    name = {"java": "vardict-java", "cpp": "vardictcpp"}[kind]
    for line in (ROOT / "benchmarks" / "env" / "vardict-oracle.yml").read_text().splitlines():
        s = line.strip()
        if s.startswith("- " + name + "="):
            return s[2:]
    sys.exit(f"golden_store: {name} pin not found in benchmarks/env/vardict-oracle.yml")


def normalise_argv(args, inputs):
    """Replace every input path inside argv with @IN<i> so checkouts/worktrees share cells."""
    order = sorted(enumerate(map(str, inputs)), key=lambda t: -len(t[1]))
    out = []
    for a in args:
        for i, p in order:
            a = a.replace(p, f"@IN{i}")
        out.append(a)
    return out


def cell_key(kind, inputs, args, wrapper="-", literal=()):
    """Key = sha256 of everything that can change the oracle's output.

    `literal` inputs are hashed but their path string stays in argv: VarDict embeds the -G path
    verbatim in `##reference=` and compare_vcf.py does not exempt it, so a different spelling of
    the same FASTA must be a different cell. The thread count (-th) is deliberately kept: a time
    file measured at one -th must not be replayed into a run at another."""
    lit = {str(p) for p in literal}
    subst = [p for p in inputs if str(p) not in lit]
    ident = {
        "schema": SCHEMA,
        "kind": kind,
        "pin": oracle_pin(kind),
        "env_pin": sha256_file(ROOT / "benchmarks" / "env" / "vardict-oracle.yml"),
        "wrapper": wrapper,
        "inputs": [sha256_file(p) for p in inputs],
        "argv": normalise_argv(list(args), subst),
        "env": {k: os.environ.get(k, "") for k in ENV_KEYS[kind]},
    }
    blob = json.dumps(ident, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(blob).hexdigest()[:16], ident


def _sigs(inputs):
    out = []
    for p in inputs:
        st = os.stat(p)
        out.append((st.st_size, st.st_mtime_ns, st.st_ctime_ns))
    return out


# --- manifest + artifacts ---------------------------------------------------------

def cell_dir(key):
    return store_dir() / "cells" / key


@contextlib.contextmanager
def _cell_lock(key):
    """mkdir-based lock (atomic, stdlib, portable); a lock older than LOCK_STALE_S is a dead holder."""
    d = cell_dir(key)
    d.mkdir(parents=True, exist_ok=True)
    lock = d / ".lock"
    deadline = time.time() + LOCK_STALE_S
    while True:
        try:
            os.mkdir(lock)
            break
        except FileExistsError:
            try:
                if time.time() - lock.stat().st_mtime > LOCK_STALE_S:
                    shutil.rmtree(lock, ignore_errors=True)
                    continue
            except OSError:
                continue
            if time.time() > deadline:
                sys.exit(f"golden_store: timed out waiting for the lock on cell {key}")
            time.sleep(0.05)
    try:
        yield
    finally:
        shutil.rmtree(lock, ignore_errors=True)


def load_manifest(key):
    p = cell_dir(key) / "manifest.json"
    if not p.exists():
        return {"key": key, "entries": []}
    try:
        return json.loads(p.read_text())
    except (OSError, ValueError):
        sys.exit(f"golden_store: manifest of cell {key} is corrupt ({p}); inspect it, then tombstone or delete the cell")


def live_entry(key):
    """(index, entry) of the newest entry, or None when absent or tombstoned (append-only log)."""
    entries = load_manifest(key)["entries"]
    if not entries or entries[-1].get("tombstone"):
        return None
    return len(entries) - 1, entries[-1]


def host_info():
    cpu = "unknown"
    try:
        with open("/proc/cpuinfo") as f:
            for line in f:
                if line.startswith("model name"):
                    cpu = line.split(":", 1)[1].strip()
                    break
    except OSError:
        pass
    return {"cpu": cpu, "nproc": os.cpu_count(), "kernel": platform.release()}


def _host_matches(recorded):
    cur = host_info()
    return recorded.get("cpu") == cur["cpu"] and recorded.get("nproc") == cur["nproc"]


def _store_file(src, dest_dir, name):
    raw = pathlib.Path(src).stat().st_size
    sha = _hash(src)
    gz = raw >= GZ_MIN_BYTES
    fname = name + (".gz" if gz else "")
    with open(src, "rb") as fi, open(dest_dir / fname, "wb") as fo:
        if gz:
            with gzip.GzipFile(filename="", mode="wb", fileobj=fo, mtime=0) as z:
                shutil.copyfileobj(fi, z, 1 << 20)
        else:
            shutil.copyfileobj(fi, fo, 1 << 20)
    return {"file": fname, "sha256": sha, "bytes": raw, "gz": gz}


def _restore_file(entry_dir, rec, dest=None):
    """Stream an artifact (decompressing), optionally writing it to dest; True iff sha256 matches."""
    h = hashlib.sha256()
    opener = gzip.open if rec["gz"] else open
    fo = None
    try:
        if dest is not None:
            dest = pathlib.Path(dest)
            dest.parent.mkdir(parents=True, exist_ok=True)
            fo = open(dest, "wb")
        with opener(entry_dir / rec["file"], "rb") as fi:
            for chunk in iter(lambda: fi.read(1 << 20), b""):
                h.update(chunk)
                if fo:
                    fo.write(chunk)
    except (OSError, EOFError, zlib.error):
        return False
    finally:
        if fo:
            fo.close()
    return h.hexdigest() == rec["sha256"]


def put(key, ident, label, kind, files, time_file):
    with _cell_lock(key):
        d = cell_dir(key)
        m = load_manifest(key)
        n = len(m["entries"])
        tmp = d / f".tmp-{os.getpid()}"
        shutil.rmtree(tmp, ignore_errors=True)
        tmp.mkdir(parents=True)
        try:
            arts = {role: _store_file(path, tmp, role) for role, path in files.items()}
            if time_file and pathlib.Path(time_file).is_file() and os.path.getsize(time_file) > 0:
                arts["time"] = _store_file(time_file, tmp, "time")
            shutil.rmtree(d / str(n), ignore_errors=True)
            os.replace(tmp, d / str(n))
        except BaseException:
            shutil.rmtree(tmp, ignore_errors=True)
            raise
        m["entries"].append({
            "schema": SCHEMA, "label": label, "kind": kind,
            "pin": ident["pin"], "env_pin": ident["env_pin"], "wrapper": ident["wrapper"],
            "inputs": ident["inputs"], "argv": ident["argv"], "env": ident["env"],
            "artifacts": arts, "host": host_info(), "date_utc": now_utc(),
        })
        m["key"] = key
        atomic_write(d / "manifest.json", json.dumps(m, indent=1, sort_keys=True).encode())
        for t in d.glob(".tmp-*"):  # leftovers of killed runs
            try:
                if time.time() - t.stat().st_mtime > TMP_GC_S:
                    shutil.rmtree(t, ignore_errors=True)
            except OSError:
                pass


def log_event(event, kind, key, label):
    sd = store_dir()
    sd.mkdir(parents=True, exist_ok=True)
    fields = [now_utc(), os.environ.get("GOLDEN_RUN_ID", "-"), event, kind, key, label or "-"]
    with open(sd / "stats.log", "a") as f:
        f.write("\t".join(fields) + "\n")


# --- run-or-replay core -------------------------------------------------------------

def _tmp(prefix):
    f = tempfile.NamedTemporaryFile(delete=False, prefix=prefix)
    f.close()
    return f.name


def _stream(path, dest):
    with open(path, "rb") as f:
        shutil.copyfileobj(f, dest)


def _on_term(signum, frame):
    raise SystemExit(143)


def _kill_group(proc):
    for sig in (signal.SIGTERM, signal.SIGKILL):
        try:
            os.killpg(proc.pid, sig)
        except ProcessLookupError:
            return
        try:
            proc.wait(timeout=10)
            return
        except subprocess.TimeoutExpired:
            continue


def _run_live(cmd, stdout_path, err_path, time_path):
    full = ([TIME_BIN, "-v", "-o", time_path] if os.path.exists(TIME_BIN) else []) + list(cmd)
    so = open(stdout_path, "wb") if stdout_path else None
    try:
        with open(err_path, "wb") as err:
            proc = subprocess.Popen(full, stdout=so, stderr=err, start_new_session=True)
            try:
                return proc.wait()
            except BaseException:  # SIGTERM (from `timeout`) / Ctrl-C: take the whole command group down
                _kill_group(proc)
                raise
    finally:
        if so:
            so.close()


def _replay(key, got, out, err_out, time_out):
    """'hit' | 'corrupt' | 'host' (a time file from another host is never replayed)."""
    idx, e = got
    ed = cell_dir(key) / str(idx)
    arts = e["artifacts"]
    if time_out and "time" in arts and not _host_matches(e.get("host", {})):
        return "host"
    with tempfile.TemporaryDirectory() as td:
        td = pathlib.Path(td)
        if not _restore_file(ed, arts["out"], out):
            return "corrupt"
        if not _restore_file(ed, arts["err"], td / "err"):
            return "corrupt"
        if "time" in arts and not _restore_file(ed, arts["time"], td / "time"):
            return "corrupt"
        if err_out:
            shutil.copyfile(td / "err", err_out)
        else:
            _stream(td / "err", sys.stderr.buffer)
            sys.stderr.flush()
        if time_out:
            if "time" in arts:
                shutil.copyfile(td / "time", time_out)
            else:
                pathlib.Path(time_out).unlink(missing_ok=True)
    return "hit"


def _plausible(kind, out, err_path):
    if os.path.getsize(out) == 0 and os.path.getsize(err_path) == 0:
        return False, "empty output and empty stderr"
    if kind == "java" and str(out).endswith(".vcf"):
        with open(out, "rb") as f:
            if b"#CHROM" not in f.read(1 << 22):
                return False, "java VCF without a #CHROM header line"
    return True, ""


def _pin(m, kind, label, key, ident, out, ef, tf, rc, corrupt, sigs, inputs):
    if m == "off":
        return "OFF"
    if rc != 0:
        return "FAIL"
    if not pathlib.Path(out).is_file():
        print(f"golden_store: {label}: command exited 0 but {out} is missing; not pinned", file=sys.stderr)
        return "FAIL"
    ok, why = _plausible(kind, out, ef)
    if not ok:
        print(f"golden_store: {label}: not pinned ({why})", file=sys.stderr)
        return "REFUSED"
    if _sigs(inputs) != sigs:
        print(f"golden_store: {label}: an input changed during the run; not pinned", file=sys.stderr)
        return "UNSTABLE"
    try:
        put(key, ident, label, kind, {"out": out, "err": ef}, tf)
    except OSError as e:
        print(f"golden_store: {label}: pin failed ({e}); the result of this run is NOT stored", file=sys.stderr)
        return "PINFAIL"
    return "REGEN" if m == "regen" else ("CORRUPT" if corrupt else "MISS")


def run_cached(kind, label, inputs, key_args, cmd, out, err_out=None, capture_stdout=False,
               wrapper="-", literal=()):
    m = mode()
    time_out = os.environ.get("GOLDEN_TIME_OUT")
    if time_out and not os.path.exists(TIME_BIN):
        sys.exit(f"golden_store: GOLDEN_TIME_OUT needs {TIME_BIN}")
    key = ident = None
    corrupt = False
    if m != "off":
        key, ident = cell_key(kind, inputs, key_args, wrapper, literal)
    if m == "auto":
        got = live_entry(key)
        if got is not None:
            r = _replay(key, got, out, err_out, time_out)
            if r == "hit":
                log_event("HIT", kind, key, label)
                return 0
            corrupt = r == "corrupt"
    sigs = _sigs(inputs)
    ef, tf = _tmp("gs-err-"), _tmp("gs-time-")
    prev = signal.signal(signal.SIGTERM, _on_term)
    try:
        rc = _run_live(cmd, out if capture_stdout else None, ef, tf)
        if err_out:
            shutil.copyfile(ef, err_out)
        else:
            _stream(ef, sys.stderr.buffer)
            sys.stderr.flush()
        event = _pin(m, kind, label, key, ident, out, ef, tf, rc, corrupt, sigs, inputs)
        log_event(event, kind, key or "-", label)
    finally:
        signal.signal(signal.SIGTERM, prev)
        if time_out and os.path.exists(tf) and os.path.getsize(tf) > 0:
            shutil.copyfile(tf, time_out)  # also after a kill: the partial `time -v` is a diagnostic
        for p in (ef, tf):
            try:
                os.unlink(p)
            except OSError:
                pass
    return rc


def _split_dashdash(args):
    if "--" not in args:
        sys.exit("golden_store: missing '--' before the command")
    i = args.index("--")
    return args[:i], args[i + 1:]


def cmd_raw(args):
    opts, cmd = _split_dashdash(args)
    kind = label = out = err_out = None
    inputs = []
    it = iter(opts)
    for o in it:
        v = next(it, None)
        if v is None:
            sys.exit(f"golden_store: {o} needs a value")
        if o == "--kind":
            kind = v
        elif o == "--label":
            label = v
        elif o == "--out":
            out = v
        elif o == "--err-out":
            err_out = v
        elif o == "--input":
            inputs.append(v)
        else:
            sys.exit(f"golden_store: unknown option {o}")
    if kind not in ("java", "cpp") or not out or not cmd:
        sys.exit("golden_store: raw needs --kind java|cpp, --out, and a command after --")
    return run_cached(kind, label, inputs, cmd, cmd, out, err_out, capture_stdout=True)


def cmd_verify(args):
    quiet = "--quiet" in args
    bad, n = [], 0
    for mf in sorted((store_dir() / "cells").glob("*/manifest.json")):
        key = mf.parent.name
        try:
            got = live_entry(key)
        except SystemExit as e:
            bad.append(str(e))
            continue
        if got is None:
            continue
        idx, e = got
        n += 1
        for role, rec in e["artifacts"].items():
            if not _restore_file(mf.parent / str(idx), rec):
                bad.append(f"{key} ({e.get('label')}): artifact '{role}' fails its sha256")
    for line in bad:
        print(f"golden_store: CORRUPT {line}", file=sys.stderr)
    if not quiet or bad:
        print(f"golden-store verify: {n} cell(s) checked, {len(bad)} corrupt")
    return 1 if bad else 0


def cmd_stats(args):
    run = min_hit = None
    as_json = False
    it = iter(args)
    for o in it:
        if o == "--json":
            as_json = True
            continue
        v = next(it, None)
        if o == "--run" and v:
            run = v
        elif o == "--min-hit" and v:
            min_hit = float(v)
        else:
            sys.exit("golden_store: stats [--run ID] [--min-hit RATIO] [--json]")
    hit = live = 0
    events = []
    log = store_dir() / "stats.log"
    for line in (log.read_text().splitlines() if log.exists() else []):
        f = line.split("\t")
        if len(f) != 6 or (run and f[1] != run):
            continue
        if f[2] == "HIT":
            hit += 1
        elif f[2] in EVENTS_LIVE:
            live += 1
        else:
            continue
        events.append({"event": f[2], "kind": f[3], "key": f[4], "label": f[5]})
    total = hit + live
    ratio = hit / total if total else 0.0
    if as_json:
        print(json.dumps({"run": run, "hit": hit, "live": live, "events": events}, sort_keys=True))
    else:
        print(f"golden-store: {hit} hit / {live} live ({ratio:.1%} served from store)")
    return 1 if min_hit is not None and (total == 0 or ratio < min_hit) else 0


def cmd_list(args):
    for mf in sorted((store_dir() / "cells").glob("*/manifest.json")):
        m = load_manifest(mf.parent.name)
        last = m["entries"][-1] if m["entries"] else {}
        label = next((e.get("label") for e in reversed(m["entries"]) if e.get("label")), "-")
        kind = next((e.get("kind") for e in reversed(m["entries"]) if e.get("kind")), "-")
        status = "tombstoned" if last.get("tombstone") else "live"
        print(f"{m['key']}\t{kind}\t{label}\t{status}\t{len(m['entries'])}")
    return 0


def cmd_tombstone(args):
    if len(args) < 2:
        sys.exit("golden_store: tombstone KEY REASON")
    key, reason = args[0], " ".join(args[1:])
    if not (cell_dir(key) / "manifest.json").exists():
        print(f"golden_store: unknown cell key {key}", file=sys.stderr)
        return 1
    with _cell_lock(key):
        m = load_manifest(key)
        m["entries"].append({"tombstone": True, "reason": reason, "date_utc": now_utc()})
        atomic_write(cell_dir(key) / "manifest.json", json.dumps(m, indent=1, sort_keys=True).encode())
    return 0


def _cmd_wrapper(kind, env_name, script, args):
    if len(args) < 5:
        sys.exit(f"golden_store: {kind} needs LABEL REF BAM BED OUT [SAMPLE [FLAGS...]]")
    label, ref, bam, bed, out, *rest = args
    fai = ref + ".fai"
    inputs = [ref] + ([fai] if os.path.isfile(fai) else []) + bam.split("|") + [bed]
    override = os.environ.get(env_name)
    base = shlex.split(override) if override else ["bash", str(ROOT / "harness" / script)]
    wrapper = "override:" + override if override else sha256_file(ROOT / "harness" / script)
    return run_cached(kind, label, inputs, [ref] + rest, base + [ref, bam, bed, out] + rest, out,
                      wrapper=wrapper, literal=[ref])


def cmd_oracle(args):
    return _cmd_wrapper("java", "GOLDEN_ORACLE_CMD", "run_oracle.sh", args)


def cmd_cpp(args):
    return _cmd_wrapper("cpp", "GOLDEN_CPP_CMD", "run_cpp.sh", args)


VERBS = {"raw": cmd_raw, "oracle": cmd_oracle, "cpp": cmd_cpp, "verify": cmd_verify,
         "stats": cmd_stats, "list": cmd_list, "tombstone": cmd_tombstone}


def main(argv):
    if len(argv) < 2 or argv[1] not in VERBS:
        sys.exit(__doc__)
    return VERBS[argv[1]](argv[2:])


if __name__ == "__main__":
    sys.exit(main(sys.argv))
