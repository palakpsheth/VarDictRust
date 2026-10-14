# Scale-bar host-class recalibration proposal (T3D-01)

**Status:** SIGNED OFF (2026-10-13, see [§Maintainer sign-off](#maintainer-sign-off)) & IMPLEMENTED in T3D-02 (`harness/scale_bench.sh` + `benchmarks/evidence/scale.json` `PASS`).
*(Historical T3D-01 context: When drafted, the pre-registered bar in `harness/scale_bench.sh:1-10`
(`wall(th1)/wall(th8) >= 4.8 AND maxRSS(chr1 full) <= 1.5 x maxRSS(10 Mb)`)
was held frozen pending sign-off, and the initial FAIL stood until protocol §4 passed under the signed-off bar in Task 2 / T3D-02.)*

## 1. Derivation first (no observed numbers used)

The recalibrated bar is a pure function of pinned physical-core topology,
derived from the pre-registered bar alone:

1. The pre-registered speedup bar (4.8 at `-th 8`) was set on the Ryzen host
   with 8 full cores under the pin. Implied parallel efficiency: `4.8 / 8 = 0.6`.
2. Hyperthreads / E-cores do not add a full core's worth of throughput for this
   workload, so the bar must scale with **physical cores under the pin (P)**,
   not with logical CPUs or requested threads.
3. Proposed rule: **`bar(P) = 0.6 × min(requested_threads, P)`**,
   where `P` = distinct physical cores under the CPU pin (counting method §2).
   The RSS leg (`<= 1.5x`) is host-independent and unchanged.

Nothing in steps 1–3 uses any observed wall time, speedup, or RSS number
from `scratch/g4/scale/extra.json`. Observations appear only in §3 as a
post-hoc consistency check, never as an input to the rule.

## 2. Counting method (exact) + worked examples

Physical cores under the pin = distinct `(socket, core)` pairs among the
pinned logical CPUs:

```bash
taskset -pc $$            # shows the current pin, e.g. "pid $$'s current affinity list: 0-7"
lscpu -p=CPU,CORE,SOCKET  # machine-parseable topology; intersect CPU list with the pin
```

Count distinct `(CORE,SOCKET)` rows whose `CPU` is in the pin set. That count
is `P`. Use `lscpu` (not `/proc/cpuinfo` processor count, which reports
logical CPUs and overcounts under HT).

Worked examples (illustrative topology, not measurements):

- **8C/8T pin** (8 physical cores, no HT — the Ryzen-class case):
  pin CPUs `0-7` map to 8 distinct `(core, socket)` pairs → `P = 8`,
  `bar = 0.6 × min(8, 8) = 4.8`. Recovers the original bar exactly.
- **4C/8T pin** (4 physical cores with 2-way HT — the i7-13700H
  `taskset -c 0-7` case): pin CPUs `0-7` map to 4 distinct
  `(core, socket)` pairs → `P = 4`, `bar = 0.6 × min(8, 4) = 2.4`.

## 3. Why the 4.8x bar misses here (no product defect)

Observed (`scratch/g4/scale/extra.json`, pinned `taskset -c 0-7` on i7-13700H):
speedup 3.31x (th1 1984s / th8 ~600s, pre-rebase run) then 3.57x post-rebase
(th1 2771.98s / th8 776.36s); RSS leg passes (chr1-full vs 10 Mb ratio
1.365–1.45x, under the 1.5x bar). Ryzen original: 3436s / 709s = 4.85x (PASS).

Mechanism hypothesis (consistent with, not fitted to, the above):

- The `taskset -c 0-7` pin on the hybrid P/E i7-13700H captures ~4 physical
  cores with HT (8 logical CPUs), versus 8 full Ryzen cores. SMT thread pairs
  share execution resources, so scaling past P saturates early — the th4→th8
  step (881.26s → 776.36s) shows exactly this flattening.
- The i7 th1 (1984–2772s) is markedly faster than the Ryzen th1 (3436s)
  (higher single-thread boost), which compresses the `wall(th1)/wall(th8)`
  ratio even when absolute th8 throughput is competitive.
- Output identity across thread counts still holds (`cmp` cells in
  `scale_bench.sh`); the RSS leg passes. Both point away from a scaling
  defect in the engine and toward a host-class mismatch in the bar.

Post-hoc check only: 3.31x/3.57x ≥ 2.4 (the §1 rule's output for P=4) and
1.365x ≤ 1.5x. This does not set the bar — §1 does.

## 4. Re-measurement protocol (binding once signed off)

1. Same pins as the FAIL run (`taskset -c 0-7` wherever `scale_bench.sh`
   applies it); record `P` via the §2 counting method in `extra.json`.
2. Three full `harness/scale_bench.sh` runs, courtesy-wrapped
   (`harness/courtesy.sh`), same ref/BAM/BED/sample.
3. Take the **median** wall time per cell across the three runs; recompute
   `speedup_8 = median(th1)/median(th8)` and the RSS ratio from medians.
4. PASS iff median speedup ≥ `bar(P)` from §1 AND median RSS ratio ≤ 1.5
   AND all cross-thread `cmp` identity checks hold.
5. On PASS, Task 2 updates the script bar + `benchmarks/evidence/scale.json`
   via a live runner (never by hand, per integrity §0b).

Cited tokens for verification: `3.57`, `3.31`, `1.365`, `taskset`.

## 5. Addendum: RSS-leg variance watch (no bar change)

The RSS leg (`<= 1.5x`) is host-independent by construction (memory scaling,
not topology), so no host-class rule is proposed for it — the 1.5 bar stands
as written. What the signed campaign exposed instead is measurement noise:
the median RSS ratio passed at 1.448 while an individual run reached 1.489,
leaving ~1% headroom in which a noisy future run could flake the gate with
no engine change. Rule going forward (binding on the next scale campaign,
no re-measurement of past runs): retain all three runs' per-cell `wall_s`
and `maxrss_kb` tables (not just medians) alongside `extra.json`, and flag
in the campaign report any single run whose RSS ratio exceeds the bar even
when the median passes. A repeated flag across two campaigns — not a single
noisy run — is the trigger for reopening the 1.5 value itself.

## Maintainer sign-off

SIGN-OFF GRANTED by maintainer (recorded by controller), 2026-10-13: implement `bar(P) = 0.6 × min(8, P)` (P=4 → 2.4 on the `taskset -c 0-7` pin) with the RSS leg unchanged (≤ 1.5x) per §§1–2; protocol §4 binding on Task 2 (T3D-02).
