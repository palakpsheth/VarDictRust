# /snapshot-update

Update golden test snapshots after a justified, documented output change.

## Steps
1. Verify the output change is justified (not a regression)
2. Document the reason in the commit message
3. Run: `cargo insta test --accept` (for insta snapshots)
4. For golden/ files: copy new oracle output with provenance
5. Update `golden/PROVENANCE` with the change reason
6. Run `cargo test` and `harness/compare.sh` to confirm
7. Commit with provenance: `git add -A && git commit -m "snapshot-update: <reason>"`
