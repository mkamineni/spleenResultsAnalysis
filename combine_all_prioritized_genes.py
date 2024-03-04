import pandas as pd

phenos = ["spleen_original_glszm_GrayLevelNonUniformity","spleen_original_shape_Sphericity","spleen_original_glcm_Correlation","spleen_original_firstorder_Energy","spleen_original_gldm_SmallDependenceHighGrayLevelEmphasis","spleen_original_glszm_LargeAreaLowGrayLevelEmphasis","spleen_original_glrlm_RunLengthNonUniformity", "spleen_original_glcm_Id","spleen_original_firstorder_90Percentile","spleen_original_shape_Maximum2DDiameterSlice","spleen_original_gldm_LargeDependenceHighGrayLevelEmphasis","spleen_original_glszm_SmallAreaLowGrayLevelEmphasis","spleen_original_shape_Flatness"]
out_path = '../Data/fuma_results/top_all_phenos_genes_lead_atleast_one.csv'

pheno_to_genes = {}
all_genes = set()
for pheno in phenos:
	df = pd.read_csv("../Data/fuma_results/"+pheno+"/table_prioritized_genes_lead.csv", index_col = 0)
	genes = set(df.index.unique())	
	pheno_to_genes[pheno] = genes
	all_genes = all_genes.union(genes)

final_df = pd.DataFrame(columns = phenos+['count'], index = list(all_genes))
for index, row in final_df.iterrows():
	count = 0
	for pheno in phenos:
		if index in pheno_to_genes[pheno]:
			row[pheno] = 'X'
			count += 1
	row['count'] = count
#final_df = final_df[final_df["count"]>1]
final_df = final_df.sort_values(by = ['count']+phenos, ascending = [False]*14).drop('count', axis = 1)
final_df.to_csv(out_path)


