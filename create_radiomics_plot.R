library(data.table)
library(ggplot2)
library(dplyr)
library(pheatmap)

file <- "../Data/radiomics_data/CADcohort_all_demo_rad_final_nov1.csv"
df <- read.csv(file, header=TRUE)
out <- "../Data/figures/rad_heatmap.png"
out_small <- "../Data/figures/rad_heatmap_select.png"

colnames(df) <- gsub('spleen_original_', '',colnames(df))
df = subset(df, select = -c(Prevalent_Coronary_Artery_Disease_INTERMEDIATE, ID, pce_goff, time_to_follow_up, Years_To_Coronary_Artery_Disease_INTERMEDIATE, train)) 
# “race_asian”, “race_black”, “race_mixed”,“race_other”, “race_white”, “sex_Female”, “sex_Male”))
df = df %>% relocate("age", .after="sex_Male")
pheatmap(cor(df),cluster_rows = FALSE, cluster_cols =FALSE, fontsize = 20, filename=out, width=30, height=30)

df = subset(df, select = c(shape_Sphericity, shape_MinorAxisLength, gldm_DependenceEntropy, glszm_GrayLevelNonUniformity, glcm_Correlation, race_asian, race_black, race_mixed, race_other, race_white, sex_Female, sex_Male, age))
pheatmap(cor(df),cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 10, filename=out_small)
