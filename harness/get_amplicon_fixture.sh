#!/usr/bin/env bash
# Setup multi-amplicon overlapping panel fixture exercising primer trimming and AMPBIAS.
set -euo pipefail

OUTDIR="scratch/g4/amplicon"
mkdir -p "$OUTDIR"

REF_FA="VarDictJava/testdata/fastas/hg19.fa"
SRC_BAM="VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/L861Q.bam"

[ -f "$REF_FA" ] || { echo "Missing reference $REF_FA" >&2; exit 1; }
[ -f "$SRC_BAM" ] || { echo "Missing BAM $SRC_BAM" >&2; exit 1; }

# Symlink reference and BAM into fixture dir
ln -sf "$(realpath "$REF_FA")" "$OUTDIR/ref.fa"
ln -sf "$(realpath "$REF_FA.fai")" "$OUTDIR/ref.fa.fai"
ln -sf "$(realpath "$SRC_BAM")" "$OUTDIR/amplicon.bam"
ln -sf "$(realpath "$SRC_BAM.bai")" "$OUTDIR/amplicon.bam.bai"

# Setup multi-amplicon overlapping panel (correct_8columns_amplicon.bed)
# Amplicon 1: 55259400-55259600, insert 55259420-55259580
# Amplicon 2: 55259500-55259650, insert 55259510-55259640
# Both amplicons overlap and exercise primer trimming (insert < full region).
# Variant at 55259519 is called in amp1 (GDAMP=1), covered by amp2 (TLAMP=2)
# but not called in amp2, triggering GDAMP < TLAMP - NCAMP -> AMPBIAS filter!
cat << 'EOF' > "$OUTDIR/correct_8columns_amplicon.bed"
chr7	55259400	55259600	amp1	100	+	55259420	55259580
chr7	55259500	55259650	amp2	100	+	55259510	55259640
EOF

echo "Amplicon fixture prepared at $OUTDIR"
