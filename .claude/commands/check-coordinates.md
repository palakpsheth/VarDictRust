# /check-coordinates

Audit coordinate usage for 0-based vs 1-based safety.

## Steps
1. Search for raw integer arithmetic on Pos0/Pos1 values
2. Verify all BAM positions use Pos0, all VCF/VarDict positions use Pos1
3. Check that conversions only go through `.to_one_based()` / `.to_zero_based()`
4. Flag any raw `.0` field access that could mix coordinate systems
5. Report findings with file:line citations
