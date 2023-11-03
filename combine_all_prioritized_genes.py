import pandas as pd

phenos = ["sphericity", "graynonuniform", "kurtosis", "graydependence", "grayemphasis", "grayinversediff"]
out_path = '../Data/fuma_results/top_all_phenos_genes_lead.csv'

pheno_to_genes = {}
all_genes = set()
for pheno in phenos:
	df = pd.read_csv("../Data/fuma_results/"+pheno+"_common_var/table_prioritized_genes_lead.csv", index_col = 0)
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
final_df = final_df[final_df["count"]>1]
final_df = final_df.sort_values(by = 'count', ascending = False).drop('count', axis = 1)
final_df.to_csv(out_path)


