#!/bin/bash
set -e

BAM="NA12878.chrom20.ILLUMINA.bwa.CEU.low_coverage.20121211.bam"
REF="hs37_chr20.fa"
REGION="20"
THREADS=4
OUTDIR="results"
DEPTHS=(10 20 30)

mkdir -p $OUTDIR/gatk $OUTDIR/freebayes $OUTDIR/bcftools $OUTDIR/stats

echo "=== Step 1: BAM stats ==="
samtools flagstat $BAM > $OUTDIR/stats/flagstat.txt
cat $OUTDIR/stats/flagstat.txt

echo "=== Step 2: Read depth ==="
samtools depth -r $REGION $BAM > $OUTDIR/stats/depth.txt
python3 summarize_depth.py $OUTDIR/stats/depth.txt

echo "=== Step 3: GATK HaplotypeCaller ==="
gatk HaplotypeCaller \
    -R $REF -I $BAM \
    -O $OUTDIR/gatk/raw.vcf.gz \
    -L $REGION \
    --native-pair-hmm-threads $THREADS \
    --verbosity ERROR
echo "GATK done"

for D in "${DEPTHS[@]}"; do
    gatk VariantFiltration \
        -R $REF \
        -V $OUTDIR/gatk/raw.vcf.gz \
        -O $OUTDIR/gatk/filtered_dp${D}.vcf.gz \
        --filter-expression "DP < ${D}" \
        --filter-name "LowDepth" \
        --verbosity ERROR
    bcftools view -f PASS $OUTDIR/gatk/filtered_dp${D}.vcf.gz \
        -o $OUTDIR/gatk/pass_dp${D}.vcf.gz -O z
    bcftools index $OUTDIR/gatk/pass_dp${D}.vcf.gz
    COUNT=$(bcftools stats $OUTDIR/gatk/pass_dp${D}.vcf.gz | grep "number of records:" | cut -f4)
    echo "  GATK DP>=${D}: ${COUNT} variants"
done

echo "=== Step 4: FreeBayes ==="
freebayes -f $REF -r $REGION \
    --min-mapping-quality 20 --min-base-quality 20 \
    $BAM > $OUTDIR/freebayes/raw.vcf
bgzip $OUTDIR/freebayes/raw.vcf
bcftools index $OUTDIR/freebayes/raw.vcf.gz
echo "FreeBayes done"

for D in "${DEPTHS[@]}"; do
    bcftools filter -i "DP>=${D} && QUAL>=20" \
        $OUTDIR/freebayes/raw.vcf.gz \
        -o $OUTDIR/freebayes/filtered_dp${D}.vcf.gz -O z
    bcftools index $OUTDIR/freebayes/filtered_dp${D}.vcf.gz
    COUNT=$(bcftools stats $OUTDIR/freebayes/filtered_dp${D}.vcf.gz | grep "number of records:" | cut -f4)
    echo "  FreeBayes DP>=${D}: ${COUNT} variants"
done

echo "=== Step 5: bcftools mpileup/call ==="
bcftools mpileup -f $REF -r $REGION \
    --min-MQ 20 --min-BQ 20 \
    --annotate FORMAT/DP,FORMAT/AD \
    $BAM | \
bcftools call --multiallelic-caller --variants-only \
    -o $OUTDIR/bcftools/raw.vcf.gz -O z
bcftools index $OUTDIR/bcftools/raw.vcf.gz
echo "bcftools done"

for D in "${DEPTHS[@]}"; do
    bcftools filter -i "DP>=${D} && QUAL>=20" \
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
