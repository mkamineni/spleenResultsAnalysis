import pandas as pd

dir = "../Data/fuma_results/grayemphasis_common_var/"
snp_annot_path = dir+"snps.txt"
sig_snp_path = dir+"leadSNPs.txt"
out_path = dir+"nearest_genes_lead_snps.csv"

#read in all SNP annotations
snp_annot = pd.read_csv(snp_annot_path, sep = "\t")
print(snp_annot.head(5))

#read in sig SNPs
sig_snps = pd.read_csv(sig_snp_path, sep = "\t")

snp_annot = snp_annot[["rsID", "GenomicLocus", "nearestGene"]]
sig_snps = sig_snps[["rsID", "GenomicLocus"]]
print(snp_annot.shape[0])
print(sig_snps.shape[0])

merged = pd.merge(snp_annot, sig_snps, on = ["rsID", "GenomicLocus"], how = "inner")
print(merged.shape[0])
print(len(merged["nearestGene"].unique()))
merged.to_csv(out_path, index=False)
