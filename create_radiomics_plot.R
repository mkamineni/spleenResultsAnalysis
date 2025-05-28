library(psych)
library(data.table)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(heatmaply)
library(seriation)
library(dendextend)
library(tidyverse)
library(Hmisc)
library(scales)
library(grid)
library(gridExtra)

file <- "../Data/radiomics_data/CADcohort_all_demo_rad_final_nov1.csv"
file2 <- "../Data/radiomics_data/CAD_pat_char.csv"
coefs <- read.csv(file, header=TRUE)
out <- "../Data/figures/rad_heatmap.pdf"
out_small <- "../Data/figures/rad_heatmap_select.pdf"

colnames(coefs) <- gsub("spleen_original_", "", colnames(coefs))
colnames(coefs) <- gsub("shape_", "", colnames(coefs))
colnames(coefs) <- gsub("firstorder_", "", colnames(coefs))
colnames(coefs) <- gsub("glcm_", "", colnames(coefs))
colnames(coefs) <- gsub("gldm_", "", colnames(coefs))
colnames(coefs) <- gsub("ngtdm_Strength", "Strength", colnames(coefs))
colnames(coefs) <- gsub("ngtdm_Complexity", "Complexity", colnames(coefs))
colnames(coefs) <- gsub("ngtdm_Busyness", "Busyness", colnames(coefs))
colnames(coefs) <- gsub("ngtdm_Coareness", "Coarseness", colnames(coefs))

colnames(coefs) <- gsub("glszm_ZoneVariance", "ZoneVariance", colnames(coefs))
colnames(coefs) <- gsub("glszm_ZonePercentage", "ZonePercentage", colnames(coefs))
colnames(coefs) <- gsub("glszm_ZoneEntropy", "ZoneEntropy", colnames(coefs))
colnames(coefs) <- gsub("glszm_SmallAreaLowGrayLevelEmphasis", "SmallAreaLowGrayLevelEmphasis", colnames(coefs))
colnames(coefs) <- gsub("glszm_SmallAreaHighGrayLevelEmphasis", "SmallAreaHighGrayLevelEmphasis", colnames(coefs))
colnames(coefs) <- gsub("glszm_SmallAreaEmphasis", "SmallAreaEmphasis", colnames(coefs))
colnames(coefs) <- gsub("glszm_SizeZoneNonUniformityNormalized", "SizeZoneNonUniformityNormalized", colnames(coefs))
colnames(coefs) <- gsub("glszm_LowGrayLevelZoneEmphasis", "LowGrayLevelZoneEmphasis", colnames(coefs))
colnames(coefs) <- gsub("glszm_LargeAreaLowGrayLevelEmphasis", "LargeAreaLowGrayLevelEmphasis", colnames(coefs))
colnames(coefs) <- gsub("glszm_LargeAreaHighGrayLevelEmphasis", "LargeAreaHighGrayLevelEmphasis", colnames(coefs))
colnames(coefs) <- gsub("glszm_LargeAreaEmphasis", "LargeAreaEmphasis", colnames(coefs))
colnames(coefs) <- gsub("glszm_GrayLevelZoneEmphasis", "GrayLevelZoneEmphasis", colnames(coefs))

colnames(coefs) <- gsub("glrlm_ShortRunLowGrayLevelEmphasis", "ShortRunLowGrayLevelEmphasis", colnames(coefs))
colnames(coefs) <- gsub("glrlm_ShortRunHighGrayLevelEmphasis", "ShortRunHighGrayLevelEmphasis", colnames(coefs))
colnames(coefs) <- gsub("glrlm_ShortRunEmphasis", "ShortRunEmphasis", colnames(coefs))
colnames(coefs) <- gsub("glrlm_RunVariance", "RunVariance", colnames(coefs))
colnames(coefs) <- gsub("glrlm_RunPercentage", "RunPercentage", colnames(coefs))
colnames(coefs) <- gsub("glrlm_RunLengthNonUniformity", "RunLengthNonUniformity", colnames(coefs))
colnames(coefs) <- gsub("glrlm_RunLengthNonUniformityNormalized", "RunLengthNonUniformityNormalized", colnames(coefs))
colnames(coefs) <- gsub("glrlm_RunEntropy", "RunEntropy", colnames(coefs))
colnames(coefs) <- gsub("glrlm_LowGrayLevelRunEmphasis", "LowGrayLevelRunEmphasis", colnames(coefs))
colnames(coefs) <- gsub("glrlm_LongRunLowGrayLevelEmphasis", "LongRunLowGrayLevelEmphasis", colnames(coefs))
colnames(coefs) <- gsub("glrlm_LongRunHighGrayLevelEmphasis", "LongRunHighGrayLevelEmphasis", colnames(coefs))
colnames(coefs) <- gsub("glrlm_LongRunEmphasis", "LongRunEmphasis", colnames(coefs))
colnames(coefs) <- gsub("glrlm_HighGrayLevelRunEmphasis", "HighGrayLevelRunEmphasis", colnames(coefs))

df<- coefs
splenic_feat_colors <- c(rep("Brown", 14), rep("Darkgreen", 18), rep("Red", 24), rep("Purple", 14), rep("Black", 16), rep("Darkturquoise", 16), rep("Blue", 5))

df = subset(df, select = -c(Prevalent_Coronary_Artery_Disease_INTERMEDIATE, ID, pce_goff, time_to_follow_up, Years_To_Coronary_Artery_Disease_INTERMEDIATE, train, race_asian, race_black, race_mixed,race_other, race_white, sex_Female, sex_Male, age))

o = rownames(cor(df))
od =  hclust(dist(cor(df)))$order
m2 = cor(df)[od, od]
m2[lower.tri(m2)]<- NA
m2 = m2[o, o]


pheatmap(m2,fillename=out, cluster_rows=hclust(dist(cor(df))), cluster_cols =hclust(dist(cor(df))), fontsize = 20, filename=out, width=30, height=30, na_col="white", border_color="white", device="pdf")
#t_df <- transpose(df)
#t_df$colors = splenic_feat_colors
#e$gtable$grobs[[4]]$gp <- gpar(col=splenic_feat_colors)
#cols=t_df[order(match(colnames(t_df), e$gtable$grobs[[5]]$label)), ]$colors
#e$gtable$grobs[[5]]$gp=gpar(col=splenic_feat_colors)
#png(out, width=30, height=30)
#grid::grid.newpage()
#grid::grid.draw(e$gtable)
#dev.off()


print("Mesh and Voxel Corr")
a = "Sphericity"
b = "Elongation"
cor_val = cor(df[a], df[b], method = "pearson", use="complete.obs")
print(cor_val)
corr_res = corr.test(df[a], df[b], method="pearson", alpha=0.05,ci=TRUE)
print(corr_res$p)


# Make a heat map for the important radiomics features

cors <- function(df) {
	M = Hmisc::rcorr(as.matrix(df)) 
	print(head(M, 3))
	Mdf = map(M, ~data.frame(.x)) 
	print(head(Mdf, 3))
	return(Mdf) 
}

formatted_cors <- function(df){
	new_df = cors(df)
	new_df = new_df %>% 
	map(~rownames_to_column(.x, var="measure1")) %>% 
	map(~pivot_longer(.x, -measure1, "measure2")) %>% 
	bind_rows(.id = "id") %>% 
	pivot_wider(names_from = id, values_from = value) %>%
	mutate(sig_p = ifelse(P < .05, T, F), p_if_sig = ifelse(P <.05, P, NA), r_if_sig = ifelse(P <.05, r, NA)) 
}

new_df <- read.csv(file2, header=TRUE)
colnames(new_df) <- gsub('spleen_original_', '',colnames(new_df))
print(colnames)

new_df = subset(new_df, select = c(glszm_GrayLevelNonUniformity, shape_Sphericity, glcm_Correlation, firstorder_Energy, gldm_SmallDependenceHighGrayLevelEmphasis, gldm_GrayLevelVariance, glszm_LargeAreaLowGrayLevelEmphasis, glrlm_RunLengthNonUniformity, glcm_Id, glcm_Idn, age, race_white, race_asian, race_mixed, race_other, ever_smoked,BMI,SBP,HDL.cholesterol,LDL.direct,Triglycerides,Prev_Diabetes_Type_2,Prev_Hypercholesterolemia,Prev_Hypertension,Total.Cholesterol))

colnames(new_df)= c("GLSZM Gray Level Non-Uniformity", "Sphericity", "GLCM Correlation", "Energy", "GLDM Small Dependence High Gray Level Emphasis", "GLDM Gray Level Variance", "GLSZM Large Area Low Gray Level Emphasis", "GLRLM Run Length Non-Uniformity", "GLCM Inverse Difference", "GLCM Inverse Difference Normalized", "Age", "Race White", "Race Asian", "Race Mixed", "Race Other", "Ever Smoked", "BMI", "SBP", "HDL", "LDL", "Triglycerides", "Prev Type 2 Diabetes", "Prev Hyper cholesterolemia", "Prev Hypertension", "Total Cholesterol")
splenic_feat_colors <- c('Darkturquoise', 'Brown', 'Red', 'Darkgreen', 'Purple', 'Red', 'Darkturquoise', 'Black', 'Red', 'Red')

#corr <- cor(new_df, use="complete.obs")
#corr <- corr[c("Gray Level Non-Uniformity", "Sphericity", "Correlation", "Energy", "Small Dependence High Gray Level Emphasis", "Gray Level Variance", "Large Area Low Gray Level Emphasis", "Run Length Non-Uniformity", "Inverse Difference", "90th Percentile", "Large Dependence High Gray Level Emphasis", "Small Area Low Gray Level Emphasis", "Flatness", "Max 2D Diameter Slice"),]
#pheatmap(corr,cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 10, filename=out_small)
wrap_20 <- wrap_format(25)
new_df <- formatted_cors(new_df)
print(head(new_df, 5))

new_df$measure1 <- factor(new_df$measure1, levels = unique(new_df$measure1))
new_df <- new_df[! new_df$measure2 %in% c("Age", "Race.White", "Race.Asian", "Race.Mixed", "Race.Other", "Ever.Smoked", "BMI", "SBP", "HDL", "LDL", "Triglycerides", "Prev.Type.2.Diabetes", "Prev.Hyper.cholesterolemia", "Prev.Hypertension", "Total.Cholesterol"), ]
print(unique(new_df$measure2))
new_df$measure2 <- factor(new_df$measure2, levels = unique(new_df$measure2))

new_df %>% 
 ggplot(aes(measure1, measure2, fill=r, label=round(r_if_sig,2))) +
 geom_tile() + 
 geom_text(size = 10) +
 #theme_classic() +
 labs(x = NULL, y = NULL, fill = "Pearson's\nCorrelation") + 
 scale_fill_gradient2(low="#888BC9",high="#3A68AE", limits=c(-1,1)) +
 scale_x_discrete(labels = wrap_20(unique(new_df$measure1)), expand=c(0,0)) + 
 scale_y_discrete(labels = wrap_20(unique(new_df$measure1)[1:14]), expand=c(0,0)) +
 theme(axis.text.y = element_text(size=30, color=splenic_feat_colors), axis.text.x = element_text(size=30, angle=90), legend.text = element_text(size=25), legend.title = element_text(size=30))


ggsave(out_small, height=15, width=35, device = "pdf")


