#!/usr/bin/env Rscript
# Deterministic R fisher.test oracle generator for VarDictRust (Issues B, M & N).
# Executes on Linux x86_64 in the vardict-oracle conda environment:
#   "${CONDA_PREFIX:-$HOME/miniconda3}/envs/vardict-oracle/bin/Rscript" fixtures/jvm-oracles/RFisherOracle.R > crates/vardictrust-core/testdata/fisher/r_fisher_oracle.tsv

options(scipen = 999)

cat("# generator=fixtures/jvm-oracles/RFisherOracle.R\n")
cat(sprintf("# r.version=%s.%s\n", R.version$major, R.version$minor))
cat(sprintf("# platform=%s\n", R.version$platform))
cat(sprintf("# sizeof.longdouble=%d\n", .Machine$sizeof.longdouble))
cat("# columns: ref_fwd\tref_rev\talt_fwd\talt_rev\tp_two\tp_less\tp_greater\todd_ratio\n")

# Read the same (ref_fwd, ref_rev, alt_fwd, alt_rev) table coordinates from java_fisher_oracle.tsv
args <- commandArgs(trailingOnly = TRUE)
in_path <- if (length(args) >= 1) args[1] else "crates/vardictrust-core/testdata/fisher/java_fisher_oracle.tsv"
lines <- readLines(in_path)

for (line in lines) {
    if (startsWith(line, "#") || nchar(trimws(line)) == 0) {
        next
    }
    parts <- strsplit(line, "\t")[[1]]
    rf <- as.integer(parts[1])
    rr <- as.integer(parts[2])
    af <- as.integer(parts[3])
    ar <- as.integer(parts[4])

    # (145, 494, 1, 0) is a Java --fisher (cm3 saddle-point) pin in java_fisher_oracle.tsv.
    # In R mode, dhyper(145, 639, 1, 146) = 73/320 = 0.228125 (an exact 5-decimal midpoint
    # that is a repeating binary fraction), where glibc libm (0x3fcd333333333333) and
    # musl libm (0x3fcd333333333334) straddle the midpoint by 1 ULP.
    if (rf == 145L && rr == 494L && af == 1L && ar == 0L) {
        next
    }

    mat <- matrix(c(rf, rr, af, ar), nrow = 2)
    h_two <- fisher.test(mat, alternative = "two.sided")
    h_less <- fisher.test(mat, alternative = "less")
    h_greater <- fisher.test(mat, alternative = "greater")

    p_two <- round(h_two$p.value, 5)
    p_less <- round(h_less$p.value, 5)
    p_greater <- round(h_greater$p.value, 5)
    est <- unname(h_two$estimate)
    or_val <- if (is.infinite(est)) "Inf" else as.character(round(est, 5))

    cat(sprintf("%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n",
        rf, rr, af, ar,
        as.character(p_two),
        as.character(p_less),
        as.character(p_greater),
        or_val
    ))
}
