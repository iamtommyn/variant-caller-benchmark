import os, subprocess
import matplotlib.pyplot as plt
import pandas as pd

OUTDIR = "results"
DEPTHS = [10, 20, 30]
TOOLS  = {"GATK": "gatk/pass", "FreeBayes": "freebayes/filtered", "bcftools": "bcftools/filtered"}

def count_vcf(path):
    if not os.path.exists(path): return 0
    r = subprocess.run(["bcftools","stats",path], capture_output=True, text=True)
    for line in r.stdout.split('\n'):
        if line.startswith('SN') and 'number of records' in line:
            return int(line.strip().split('\t')[3])
    return 0

def get_stats(path):
    snps = indels = 0
    if not os.path.exists(path): return snps, indels
    r = subprocess.run(["bcftools","stats",path], capture_output=True, text=True)
    for line in r.stdout.split('\n'):
        if line.startswith('SN') and 'number of SNPs' in line:
            snps = int(line.strip().split('\t')[3])
        if line.startswith('SN') and 'number of indels' in line:
            indels = int(line.strip().split('\t')[3])
    return snps, indels

rows = []
for label, prefix in TOOLS.items():
    for d in DEPTHS:
        path = f"{OUTDIR}/{prefix}_dp{d}.vcf.gz"
        rows.append({'Tool': label, 'Depth': f'DP>={d}', 'Variants': count_vcf(path)})
df = pd.DataFrame(rows)
print(df.pivot(index='Depth', columns='Tool', values='Variants').to_string())

pivot = df.pivot(index='Depth', columns='Tool', values='Variants')
ax = pivot.plot(kind='bar', figsize=(9,5), width=0.7,
                color=['#2196F3','#4CAF50','#FF5722'])
ax.set_title("Variants Detected by Tool and Minimum Read Depth", fontsize=13)
ax.set_xlabel("Min Depth Threshold"); ax.set_ylabel("Number of Variants")
ax.tick_params(axis='x', rotation=0)
for c in ax.containers: ax.bar_label(c, fmt='%d', fontsize=8, padding=2)
plt.tight_layout()
plt.savefig(f"{OUTDIR}/fig1_variant_counts.png", dpi=150)
plt.close(); print("Saved fig1_variant_counts.png")

rows2 = []
for label, prefix in TOOLS.items():
    for d in DEPTHS:
        snps, indels = get_stats(f"{OUTDIR}/{prefix}_dp{d}.vcf.gz")
        rows2.append({'Tool': label, 'Depth': f'DP>={d}', 'SNPs': snps, 'INDELs': indels})
df2 = pd.DataFrame(rows2)
fig, axes = plt.subplots(1,2,figsize=(12,5))
for ax, col in zip(axes, ['SNPs','INDELs']):
    df2.pivot(index='Depth', columns='Tool', values=col).plot(kind='bar', ax=ax, width=0.7)
    ax.set_title(col); ax.tick_params(axis='x', rotation=0)
plt.suptitle("SNP vs INDEL Detection", fontsize=13)
plt.tight_layout()
plt.savefig(f"{OUTDIR}/fig2_snp_indel.png", dpi=150)
plt.close(); print("Saved fig2_snp_indel.png")

print("\nDone! Check results/ folder for your figures.")
