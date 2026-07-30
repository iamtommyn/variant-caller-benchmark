#!/bin/bash
BAM="NA12878.chrom20.ILLUMINA.bwa.CEU.low_coverage.20121211.bam"
REF="hs37_chr20.fa"
OUTDIR="results"
DEPTHS=(10 20 30)

mkdir -p $OUTDIR/bcftools $OUTDIR/stats

echo "=== Running bcftools mpileup/call ==="
bcftools mpileup -f $REF -r 20 \
    --min-MQ 20 --min-BQ 20 \
    --annotate FORMAT/DP,FORMAT/AD \
    $BAM | \
bcftools call --multiallelic-caller --variants-only \
    -o $OUTDIR/bcftools/raw.vcf.gz -O z
bcftools index $OUTDIR/bcftools/raw.vcf.gz
echo "bcftools mpileup done"

for D in "${DEPTHS[@]}"; do
    bcftools filter -i "INFO/DP>=${D} && QUAL>=20" \
        $OUTDIR/bcftools/raw.vcf.gz \
        -o $OUTDIR/bcftools/filtered_dp${D}.vcf.gz -O z
    bcftools index $OUTDIR/bcftools/filtered_dp${D}.vcf.gz
    COUNT=$(bcftools stats $OUTDIR/bcftools/filtered_dp${D}.vcf.gz | grep "number of records:" | cut -f4)
    echo "  bcftools DP>=${D}: ${COUNT} variants"
done

echo "=== Overlap analysis ==="
for D in "${DEPTHS[@]}"; do
    echo "--- Depth >= ${D} ---"
    mkdir -p $OUTDIR/stats/isec_dp${D}
    bcftools isec \
        $OUTDIR/gatk/pass_dp${D}.vcf.gz \
        $OUTDIR/freebayes/filtered_dp${D}.vcf.gz \
        $OUTDIR/bcftools/filtered_dp${D}.vcf.gz \
        -p $OUTDIR/stats/isec_dp${D} -n=3
done

echo "=== All done! Run: python3 analyze_results.py ==="
