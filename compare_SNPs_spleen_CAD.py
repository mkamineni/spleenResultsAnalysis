import pandas as pd

#pval_cad = 2.52 * 10**(-5)
pval_cad = 5*10^(-8)
pval_spleen = 5*10^(-8)

spleen_gwas_path = "../Data/gwas_results/filt_small_spleen_chr_sphericity.regenie.gz"
cad_gwas_path = "../Data/gwas_results/cad_gwas.tsv"

spleen_gwas = pd.read_csv(spleen_gwas_path, compression = 'gzip', sep = "\t")
cad_gwas = pd.read_csv(cad_gwas_path, sep = "\t")

spleen_gwas_sig = set(spleen_gwas[spleen_gwas["Pvalue"]<pval_spleen].ID)
cad_gwas_sig = set(cad_gwas[cad_gwas["p_value"]<pval_cad].markername)

print(len(spleen_gwas_sig))
print(len(cad_gwas_sig.shape))

inter = spleen_gwas_sig.intersecton(cad_gwas_sig)
print(len(inter))
print(inter)





