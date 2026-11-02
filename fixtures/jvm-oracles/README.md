# JVM & R Kernel Oracles (`fixtures/jvm-oracles/`)

Deterministic oracle generators for VarDictJava v1.8.3 (`0886261` / `009e017`) formatting, `java.util.HashMap` operation traces, and Fisher's exact test (`FisherExact.java` and R `stats::fisher.test`).

## Provenance & Regeneration

All committed oracle tables in `crates/vardictrust-core/testdata/` are generated on `x86_64-unknown-linux-gnu` using the pinned `vardict-oracle` conda environment:

```bash
ORACLE_ENV="${CONDA_PREFIX:-$HOME/miniconda3}/envs/vardict-oracle"
VARDICT_LIB="$ORACLE_ENV/share/vardict-java-1.8.3-0/lib"
CP="$VARDICT_LIB/VarDict-1.8.3.jar:$VARDICT_LIB/commons-math3-3.6.1.jar"

javac -cp "$CP" -d /tmp/jvm-oracles fixtures/jvm-oracles/*.java
java -cp "/tmp/jvm-oracles:$CP" JavaFmtOracle > crates/vardictrust-core/testdata/format/java_fmt_oracle.tsv
java -cp "/tmp/jvm-oracles:$CP" FisherOracle > crates/vardictrust-core/testdata/fisher/java_fisher_oracle.tsv
java -cp "/tmp/jvm-oracles" JavaMapOracle crates/vardictrust-core/testdata/javamap
"$ORACLE_ENV/bin/Rscript" fixtures/jvm-oracles/RFisherOracle.R > crates/vardictrust-core/testdata/fisher/r_fisher_oracle.tsv
```
