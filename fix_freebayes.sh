#!/bin/bash
BAM="NA12878.chrom20.ILLUMINA.bwa.CEU.low_coverage.20121211.bam"
REF="hs37_chr20.fa"
OUTDIR="results"
DEPTHS=(10 20 30)

echo "=== Fixing FreeBayes filtering ==="
for D in "${DEPTHS[@]}"; do
    bcftools filter -i "INFO/DP>=${D} && QUAL>=20" \
        $OUTDIR/freebayes/raw.vcf.gz \
        -o $OUTDIR/freebayes/filtered_dp${D}.vcf.gz -O z
    bcftools index $OUTDIR/freebayes/filtered_dp${D}.vcf.gz
    COUNT=$(bcftools stats $OUTDIR/freebayes/filtered_dp${D}.vcf.gz | grep "number of records:" | cut -f4)
    echo "  FreeBayes DP>=${D}: ${COUNT} variants"
done

echo "=== Fixing bcftools filtering ==="
for D in "${DEPTHS[@]}"; do
    bcftools filter -i "INFO/DP>=${D} && QUAL>=20" \
        $OUTDIR/bcftools/raw.vcf.gz \
        -o $OUTDIR/bcftools/filtered_dp${D}.vcf.gz -O z
    bcftools index $OUTDIR/bcftools/filtered_dp${D}.vcf.gz
    COUNT=$(bcftools stats $OUTDIR/bcftools/filtered_dp${D}.vcf.gz | grep "number of records:" | cut -f4)
    echo "  bcftools DP>=${D}: ${COUNT} variants"
done

echo "=== Step 6: Overlap analysis ==="
for D in "${DEPTHS[@]}"; do
    echo "--- Depth >= ${D} ---"
    bcftools isec \
        $OUTDIR/gatk/pass_dp${D}.vcf.gz \
        $OUTDIR/freebayes/filtered_dp${D}.vcf.gz \
        $OUTDIR/bcftools/filtered_dp${D}.vcf.gz \
        -p $OUTDIR/stats/isec_dp${D} -n=3
done

echo "=== All done! Run: python3 analyze_results.py ==="
