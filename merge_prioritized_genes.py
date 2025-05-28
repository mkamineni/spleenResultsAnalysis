import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from argparse import ArgumentParser

if __name__=="__main__":
    parser = ArgumentParser()

    parser.add_argument('--trait', '-trait', 
    	help = "feature", 
    	default = None, 
    	type = str)

args = parser.parse_args()
trait = args.trait

all_genes_path = "../Data/fuma_results/"+trait+"/genes.txt"
nearest_path = "../Data/fuma_results/"+trait+"/nearest_genes_lead_snps.csv"
pops_path = "../Data/fuma_results/"+trait+"/genes_highest_pops.csv"
out_path =  "../Data/fuma_results/"+trait+"/all_prioritized_genes_lead.csv"
table_out_path = "../Data/fuma_results/"+trait+"/table_prioritized_genes_lead.csv"
filter_out_path = "../Data/fuma_results/"+trait+"/filtered_prioritized_genes_lead.csv"

def make_gene_to_prior_grid(df):
	filt_genes = set()
	gene_to_methods = dict()
	for index, row in df.iterrows():
		n, e, p = row['nearestGene'], row['eqtl_gene'], row['pops_gene']
		gene = None
		methods = set()
		if n==e:
			gene = n
			methods.add('Nearest')
			methods.add('eQTL')	
		if n==p:
			gene = n
			methods.add('Nearest')
			methods.add('PoPS')
		if p==e:
			gene = p
			methods.add('PoPS')
			methods.add('eQTL')
		if gene:
			filt_genes.add(gene)
			for method in methods:
				gene_to_methods.setdefault(gene, methods).add(method)
	print(len(filt_genes))
	sorted_keys = sorted(gene_to_methods, key=lambda key: len(gene_to_methods[key]), reverse = True)
	new_df = pd.DataFrame(columns = ['Nearest', 'PoPS', 'eQTL'], index = sorted_keys)
	for index, row in new_df.iterrows():
		methods = gene_to_methods[index]
		for elem in methods:
			row[elem] = 'X'
	new_df.to_csv(table_out_path)



nearest_genes = pd.read_csv(nearest_path, index_col=0)
nearest_genes = nearest_genes[["GenomicLocus", "nearestGene"]]
nearest_genes = nearest_genes.drop_duplicates()

pops_genes = pd.read_csv(pops_path)
pops_genes = pops_genes.rename(columns = {"symbol":"pops_gene"})
pops_genes = pops_genes[["GenomicLocus", "pops_gene"]]
pops_genes["GenomicLocus"] = pops_genes["GenomicLocus"].astype(str)
pops_genes = pops_genes.set_index('pops_gene').apply(lambda col:
    col.str.split(': ?', expand=True).stack())\
    .droplevel(1).reset_index().fillna('')
pops_genes["GenomicLocus"] = pops_genes["GenomicLocus"].astype(int)

all_genes = pd.read_csv(all_genes_path, sep="\t", index_col=0)
eqtl_genes = all_genes[all_genes["eqtlMapSNPs"]>0]
eqtl_genes = eqtl_genes[["GenomicLocus", "symbol"]]
eqtl_genes["GenomicLocus"] = eqtl_genes["GenomicLocus"].astype(str)

print("Unique eqtl Genes")
print(len(eqtl_genes.symbol.unique()))

eqtl_genes = eqtl_genes.rename(columns = {"symbol":"eqtl_gene"})
eqtl_genes = eqtl_genes.set_index('eqtl_gene').apply(lambda col:
    col.str.split(': ?', expand=True).stack())\
    .droplevel(1).reset_index().fillna('')
eqtl_genes["GenomicLocus"] = eqtl_genes["GenomicLocus"].astype(int)

merged = pd.merge(nearest_genes, pops_genes, on = "GenomicLocus", how="outer")
merged = pd.merge(merged, eqtl_genes, on = "GenomicLocus", how="outer")
print(merged.head(10))
make_gene_to_prior_grid(merged)

out_merged = merged.groupby(['GenomicLocus'], as_index=False)[['nearestGene', 'pops_gene', 'eqtl_gene']].agg(lambda x: ','.join(list(set(x.dropna()))))
out_merged.to_csv(out_path, index=False)

#new_merged = merged[["GenomicLocus", "pops_gene", "PoPS_Score"]].drop_duplicates().dropna()
#new_merged = merged[(merged["pops_gene"]==merged["eqtl_gene"]) & (merged["nearestGene"]==merged["eqtl_gene"])]
#new_merged.to_csv(filter_out_path, index=False)
