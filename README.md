# Variant Calling Pipeline Benchmark

Benchmarking three widely used variant callers — **GATK HaplotypeCaller**, **FreeBayes**, and **bcftools** — on low-coverage whole-genome sequencing data (NA12878, chromosome 20), evaluated across multiple read-depth filtering thresholds.

## Overview

Variant calling is a foundational step in genomic analysis, but different tools can produce substantially different results from the same input data. This project benchmarks the accuracy, sensitivity, agreement, and runtime of three common variant callers on a standard low-coverage sample (NA12878, chr20, ~4.83x mean depth), evaluated at minimum read depth thresholds of DP≥10, DP≥20, and DP≥30.

**Key findings:**
- At DP≥10, GATK detected the most variants (3,105), followed by bcftools (2,541) and FreeBayes (2,094), with 1,736 variants shared across all three tools.
- GATK showed the highest sensitivity but the lowest inter-tool agreement (55.9%); FreeBayes was the most conservative caller with the highest agreement rate (82.9%).
- Increasing the depth threshold sharply reduced variant counts across all tools, reflecting the limitations of the low-coverage dataset.
- bcftools was the fastest tool (~1 min); GATK was the slowest (~4.32 min).

## Repository Structure

```
.
├── variant_calling_pipeline.sh    # Main pipeline: alignment through variant calling
├── run_bcftools.sh                # bcftools calling step
├── fix_freebayes.sh               # Fix for FreeBayes depth-field filtering bug
├── analyze_results.py             # Computes summary statistics across callers
├── summarize_depth.py             # Summarizes per-base depth (from samtools depth output)
├── pipeline.log                   # Log from a full pipeline run
│
├── gatk/                          # GATK HaplotypeCaller output (VCFs)
├── freebayes/                     # FreeBayes output (VCFs)
├── bcftools/                      # bcftools output (VCFs)
│
├── isec_dp10/ isec_dp20/ isec_dp30/            # Intersection of all 3 callers, by depth threshold
├── isec_gatk_bc_dp10/ ... _dp30/               # GATK ∩ bcftools, by depth threshold
├── isec_gatk_fb_dp10/ ... _dp30/               # GATK ∩ FreeBayes, by depth threshold
├── isec_fb_bc_dp10/ ... _dp30/                 # FreeBayes ∩ bcftools, by depth threshold
│
├── fig1_variant_counts.png        # Variant counts by tool and depth threshold
├── fig2_snp_indel.png             # SNP vs. INDEL breakdown by tool
└── flagstat.txt                   # Alignment QC summary (samtools flagstat)
```

## Data

This project uses the **NA12878** sample (chromosome 20 subset), a well-characterized reference genome widely used for benchmarking in genomics.

- **Reference genome**: hg19/hs37, chr20 — not included in this repo due to size. Download from [Ensembl](https://ftp.ensembl.org/) or [UCSC Genome Browser](https://hgdownload.soe.ucsc.edu/).
- **Aligned reads (BAM)**: NA12878 chr20, Illumina low-coverage — available from [Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle) or the [1000 Genomes Project](https://www.internationalgenome.org/).

Raw reference and alignment files are excluded via `.gitignore` — see below.

## Pipeline

1. **Alignment QC** — `samtools flagstat` on the input BAM (`flagstat.txt`)
2. **Depth calculation** — `samtools depth`, summarized via `summarize_depth.py`
3. **Variant calling** — each of the three callers run independently:
   - GATK HaplotypeCaller
   - FreeBayes (with a fix for an ambiguous depth-field filtering bug — see `fix_freebayes.sh`)
   - bcftools
4. **Depth filtering** — each caller's output filtered at DP≥10, DP≥20, DP≥30
5. **Intersection analysis** — `bcftools isec` used to compare agreement between tools at each depth threshold
6. **Summary statistics & figures** — `analyze_results.py` generates variant counts, sensitivity, and agreement metrics

## Usage

```bash
# Run the full pipeline (requires reference genome + BAM in place — see Data section)
bash variant_calling_pipeline.sh

# Re-run bcftools calling step only
bash run_bcftools.sh

# Generate summary statistics and figures
python analyze_results.py
```

## Requirements

- `bwa`, `samtools`, `bcftools`
- GATK 4.x
- FreeBayes
- Python 3 (pandas, matplotlib)

## Notes

- Raw reference files (`.fa`, `.fai`, `.dict`), BAM/BAI files, and large intermediate outputs (e.g., raw per-base depth files) are excluded from version control — see `.gitignore`.
- This project accompanies a written analysis: *"Benchmarking Variant Calling Tools on Low-Coverage Whole-Genome Sequencing Data: A Comparison of GATK HaplotypeCaller, FreeBayes, and bcftools on Chromosome 20 of NA12878."*

## Author

Thomas Nguyen — Bioinformatics, Baylor University
