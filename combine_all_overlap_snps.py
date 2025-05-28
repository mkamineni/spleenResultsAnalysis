import re
import pandas as pd

phenos = ["sphericity", "graynonuniform", "kurtosis", "graydependence", "grayemphasis", "grayinversediff"]
phenos = ["spleen_original_glszm_GrayLevelNonUniformity", "spleen_original_shape_Sphericity", "spleen_original_glcm_Correlation", "spleen_original_firstorder_Energy", "spleen_original_gldm_SmallDependenceHighGrayLevelEmphasis", "spleen_original_gldm_GrayLevelVariance", "spleen_original_glszm_LargeAreaLowGrayLevelEmphasis", "spleen_original_glrlm_RunLengthNonUniformity", "spleen_original_glcm_Id", "spleen_original_glcm_Idn"]

phenos = ["spleen_original_glszm_GrayLevelNonUniformity", "spleen_original_shape_Sphericity", "spleen_original_glcm_Correlation", "spleen_original_firstorder_Energy", "spleen_original_glszm_LargeAreaLowGrayLevelEmphasis", "spleen_original_glrlm_RunLengthNonUniformity", "spleen_original_glcm_Id", "spleen_original_glcm_Idn"]

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

