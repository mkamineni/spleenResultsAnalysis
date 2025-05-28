import pandas as pd
from argparse import ArgumentParser

if __name__=="__main__":
    parser = ArgumentParser()

    parser.add_argument('--trait', '-trait', 
    	help = "feature", 
    	default = None, 
    	type = str)

args = parser.parse_args()
trait = args.trait

dir = "../Data/fuma_results/"+trait+"/"
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
print("Unique nearest genes")
print(len(merged["nearestGene"].unique()))
merged.to_csv(out_path, index=False)
