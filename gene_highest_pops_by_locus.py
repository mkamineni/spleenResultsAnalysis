import pandas as pd


trait = "spleen_original_glcm_Id"
pops_path = "../PoPs/out/r_"+trait+".preds"
fuma_path = "../Data/fuma_results/"+trait+"/genes.txt"
out_path = "../Data/fuma_results/"+trait+"/genes_highest_pops.csv"

fuma_res = pd.read_csv(fuma_path, sep = "\t")

pops_res = pd.read_csv(pops_path, sep = "\t")

fuma_res = fuma_res[["ensg", "symbol", "GenomicLocus"]] 

fuma_res = fuma_res.rename(columns = {"ensg":"ENSGID"})

pops_res = pops_res[["ENSGID", "PoPS_Score"]]

print("FUMA results %s" %str(fuma_res.shape))
print("PoPS results %s" %str(pops_res.shape))
merged = pd.merge(fuma_res, pops_res, on = "ENSGID", how = "inner")
print("merged %s" %str(merged.shape))

grouped = merged.groupby('GenomicLocus', as_index = False).apply(lambda x: x.sort_values('PoPS_Score', ascending = False))
print(grouped.head(20))
grouped = grouped.groupby('GenomicLocus', as_index = False).first().reset_index()

print("grouped %s" %str(grouped.shape))

print(grouped.head(20))
grouped.to_csv(out_path, index=False) 
