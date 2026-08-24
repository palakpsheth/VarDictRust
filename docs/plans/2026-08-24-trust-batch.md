# Trust Batch (2026-08-24) — post-audit adoption of ecosystem-survey findings

Derived from `docs/ecosystem-survey.md` gap analysis; approved amendments applied
(committed micro-fixtures before CI gate; honest-status vocabulary instead of second
flag doc; realistic fuzz targets; nightly scoping; triage policy).

## Tasks
- **G0** committed self-contained micro fixture (`fixtures/micro/` + `golden/micro.expected.var`)
  so CI can run a real case without untracked data. Gate: dual-engine MATCH vs live oracle.
- **G3** `harness/smoke.sh` (VARDICT_BIN byte-compare) + release.yml artifact gate.
  Gates: YAML valid; local SMOKE-PASS on release binary; corpus fail-set identity 91/1.
- **G1** `fuzz/` crate (workspace-excluded, own nightly pin): `fuzz_cigar`, `fuzz_regions`;
  seeds incl. crash-family inputs; `harness/fuzz.sh`; triage policy per CLAUDE.md §9.
  Gates: ≥300s/target clean; fmt/clippy/test green; corpus fail-set identity.
- **A** `docs/ecosystem-survey.md` (controller-direct).
- **B'** README Supported Flags table gains honest status column (controller-direct).
- **E** `docs/plans/giab-seqc2-validation.md`: deferred T167–T169 implementation plan (controller-direct).

## Status
G0 `3f845c8`, G3 `827de8f`, G1 `a182a69` — implementer complete; review pending.
