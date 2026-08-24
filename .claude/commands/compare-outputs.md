# /compare-outputs

Run differential testing against the VarDictJava oracle.

Usage: `/compare-outputs [case|--all]`

## Steps
1. Build release: `cargo build --release`
2. Run: `bash harness/compare.sh <ref> <bam> <bed>` for the specified case
3. Or for Tier-1: `python3 harness/compare_var.py golden/<case>.var <ours.var>`
4. Report MATCH or MISMATCH with the diff
