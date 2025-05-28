import pandas as pd

phenos = ["spleen_original_glszm_GrayLevelNonUniformity","spleen_original_shape_Sphericity","spleen_original_glcm_Correlation","spleen_original_firstorder_Energy","spleen_original_gldm_SmallDependenceHighGrayLevelEmphasis","spleen_original_glszm_LargeAreaLowGrayLevelEmphasis","spleen_original_glrlm_RunLengthNonUniformity", "spleen_original_glcm_Id","spleen_original_glcm_Idn"]
out_path = '../Data/fuma_results/top_all_phenos_genes_lead_atleast_one.csv'
phenos_len = len(phenos)+1

pheno_to_genes = {}
all_genes = set()

pops_genes = {}
eqtl_genes = {}
for pheno in phenos:
	df = pd.read_csv("../Data/fuma_results/"+pheno+"/table_prioritized_genes_lead.csv",  header='infer')
	print(df.head())
	for col in df.columns:
		print(col)
	df = df.rename(columns = {"Unnamed: 0": "Index"})
	#print(df['PoPS'])
	genes = set(df.Index.unique())	
	pheno_to_genes[pheno] = genes
	all_genes = all_genes.union(genes)
	for index, row in df.iterrows():
		print(row)
		gene = row["Index"]
		if row["PoPS"] == "X":
			pops_genes.setdefault(gene, 0)
			pops_genes[gene] +=1
		if row["eQTL"]=="X":
			eqtl_genes.setdefault(gene, 0)
			eqtl_genes[gene] +=1

sorted_pops = sorted(pops_genes.items(), reverse=True, key=lambda kv: kv[1])
sorted_eqtl = sorted(eqtl_genes.items(), reverse=True, key=lambda kv: kv[1])
print(sorted_pops)
print(sorted_eqtl)

final_df = pd.DataFrame(columns = phenos+['count'], index = list(all_genes))
print(final_df.index)
for index, row in final_df.iterrows():
	count = 0
	for pheno in phenos:
		if index in pheno_to_genes[pheno]:
			row[pheno] = 'X'
			count += 1
	row['count'] = count
	print(index)
#final_df = final_df[final_df["count"]>1]
final_df = final_df.sort_values(by = ['count']+phenos, ascending = [False]*phenos_len).drop('count', axis = 1)
print(final_df.head(5))
final_df.to_csv(out_path)


