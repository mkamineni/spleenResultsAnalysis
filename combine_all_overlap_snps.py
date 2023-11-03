import re
import pandas as pd

phenos = ["sphericity", "graynonuniform", "kurtosis", "graydependence", "grayemphasis", "grayinversediff"]
out_path = '../Data/overlap_CAD/all_overlap_snps.csv'

pheno_to_snps = {}
all_snps = set()
for pheno in phenos:
	with open("../Data/overlap_CAD/"+pheno+".csv", "r") as f:
		lines = f.readlines()
		print(lines)
		snps = [line.split()[2] for line in lines[1:-2]]
		print(snps)
		pheno_to_snps[pheno] = snps
		all_snps = all_snps.union(snps)

final_df = pd.DataFrame(columns = phenos+['count'], index = list(all_snps))
for index, row in final_df.iterrows():
        count = 0
        for pheno in phenos:
                if index in pheno_to_snps[pheno]:
                        row[pheno] = 'X'
                        count += 1
        row['count'] = count
final_df = final_df.sort_values(by = 'count', ascending = False).drop('count', axis = 1)
final_df.to_csv(out_path)

