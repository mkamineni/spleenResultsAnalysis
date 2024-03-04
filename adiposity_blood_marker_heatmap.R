library(tidyverse)
library(data.table)
library(ggplot2)
library(dplyr)
library(reshape)
library(scales)
library(janitor)
library(cowplot)
wrap_10 <- wrap_format(15)
options(repr.plot.width=60, repr.plot.height=60)

trait = "blood"
file <- "../Data/radiomics_data/radspleen_wo_heme_cancer.csv"
file2 <- "../Data/radiomics_data/CAD_pat_char.csv"
enroll_time <- "../Data/radiomics_data/CAD_coh_time.csv"
mri_time <- "../Data/radiomics_data/CAD_coh_MRI_time.csv"
time_out <- "../Data/radiomics_data/CAD_coh_days_to_MRI.csv"


if (trait == "adiposity") {
  feat_cols = c("BMI", "vat", "asat", "WHRadjbmiInstance2")
  pretty_feat_cols = c("BMI", "VAT", "ASAT", "BMI-adjusted WHR")

  x_label = "Adiposity Marker"
  out <- "../Data/figures/adiposity_heatmap.png"
  data_out <- "../Data/figures/adiposity_feat_corr.csv"
  data_out_se <- "../Data/figures/adiposity_feat_corr_se.csv"

  width <- 6
} else if (trait == "blood") {
  feat_cols = c("C.reactive.protein", "Basophill.count", "Basophill.percentage", "Eosinophill.count", "Eosinophill.percentage", "Haematocrit.percentage", "Haemoglobin.concentration", "High.light.scatter.reticulocyte.count",
  "High.light.scatter.reticulocyte.percentage", "Immature.reticulocyte.fraction", "Lymphocyte.count", "Lymphocyte.percentage", "Mean.corpuscular.haemoglobin", "Mean.corpuscular.haemoglobin.concentration", "Mean.corpuscular.volume",
  "Mean.platelet..thrombocyte..volume", "Mean.reticulocyte.volume", "Mean.sphered.cell.volume", "Monocyte.count", "Monocyte.percentage", "Neutrophill.count", "Neutrophill.percentage", "Nucleated.red.blood.cell.count",
  "Nucleated.red.blood.cell.percentage", "Platelet.count", "Platelet.crit", "Platelet.distribution.width", "Red.blood.cell..erythrocyte..count", "Red.blood.cell..erythrocyte..distribution.width", "Reticulocyte.count",
  "Reticulocyte.percentage", "White.blood.cell..leukocyte..count")
  pretty_feat_cols = c("C Reactive Protein", "Basophill Count", "Basophill Percentage", "Eosinophill Count", "Eosinophill Percentage", "Haematocrit.percentage", "Haemoglobin Concentration", "High Light Scatter Reticulocyte Count",
  "High Light Scatter Reticulocyte Percentage", "Immature Reticulocyte Fraction", "Lymphocyte Count", "Lymphocyte Percentage", "Mean Corpuscular Haemoglobin", "Mean Corpuscular Haemoglobin Concentration", "Mean Corpuscular Volume",
  "Mean platelet Volume", "Mean Reticulocyte Volume", "Mean Sphered Cell Volume", "Monocyte Count", "Monocyte Percentage", "Neutrophill Count", "Neutrophill Percentage", "Nucleated RBC Count",
  "Nucleated RBC Percentage", "Platelet Count", "Plateletcrit", "Platelet Distribution Width", "RBC Count", "RBC Distribution Width", "Reticulocyte Count",
  "Reticulocyte Percentage", "WBC Count")

  x_label = "Hematological Marker"
  out <- "../Data/figures/blood_heatmap.png"
  data_out <- "../Data/figures/blood_feat_corr.csv"
  data_out_se <- "../Data/figures/blood_feat_corr_se.csv"

  width <- 30
}

cols <- append(feat_cols, "ID")
cols <- append(cols, "BMI")

# calculate the days to MRI covariate
enroll_time <- read.csv(enroll_time, header=TRUE)
mri_time <- read.csv(mri_time, header=TRUE)

enroll_time=enroll_time[c("f.eid", "Index_Date")]
colnames(enroll_time) = c("ID", "Enroll_Date")

mri_time=mri_time[c("f.eid", "Index_Date")]
colnames(mri_time)=c("ID", "MRI_Date")

times <- merge(enroll_time, mri_time, by="ID")
times$MRI_Date <- as.Date(times$MRI_Date, "%Y-%m-%d")
times$Enroll_Date <- as.Date(times$Enroll_Date, "%Y-%m-%d")

times["days_to_MRI"] <- as.numeric(difftime(times$MRI_Date, times$Enroll_Date, units="days"))
write.csv(times, time_out, row.names=FALSE, quote=FALSE)

# open main cohort
coh <- read.csv(file, header=TRUE)
pat_char <- read.csv(file2, header=TRUE)
times <- times[!duplicated(times["ID"]),]

pat_char=pat_char[cols]
pat_char=merge(pat_char, times, by="ID", all.x=TRUE)
print(nrow(coh))
print(nrow(pat_char))

# merge coh with patient characteristics
final_coh <- merge(coh, pat_char, by = "ID", all.x = TRUE)
print(nrow(final_coh))
final_coh<-final_coh[complete.cases(final_coh$BMI),]
print(nrow(final_coh))

coef_columns <- names(final_coh)[startsWith(names(final_coh), "spleen_")]

#coef_columns = append(c("Feature", "AIC"), vec)
coefs = data.frame(matrix(nrow=1, ncol = length(coef_columns))) 
colnames(coefs) <- coef_columns
sds = data.frame(matrix(nrow=1, ncol = length(coef_columns))) 
colnames(sds) <- coef_columns

# df with columns of splenic features and rows as biomarkers
for (marker in feat_cols) {
  coefs_for_marker = c()
  ses_for_marker = c()
  for (splenic_feat in coef_columns) {
    final_coh["y"] = scale(final_coh[splenic_feat])   #gives mean of 0 and sd of 1
    final_coh["marker"] = scale(final_coh[marker])   #gives mean of 0 and sd of 1
    m1 <- glm(y ~ marker+age+sex_Female+days_to_MRI+BMI, data = final_coh, family = gaussian())
    coef = summary(m1)$coefficients["marker", "Estimate"]
    std = summary(m1)$coefficients["marker", "Std. Error"]
    coefs_for_marker <- append(coefs_for_marker, list(coef/std))
    ses_for_marker <- append(ses_for_marker, list(std))

  }
    names(coefs_for_marker) <- coef_columns
    coefs <- rbind(coefs, coefs_for_marker)
    names(ses_for_marker) <- coef_columns
    sds <- rbind(sds, ses_for_marker)

}

coefs = coefs[-1,]
rownames(coefs) <- pretty_feat_cols
sds = sds[-1,]
rownames(sds) <- pretty_feat_cols

# NOTE that coef is currently z score so change before you write to CSV
#write.csv(coefs, data_out, row.names=TRUE)
#write.csv(sds, data_out_se, row.names=TRUE)


#colnames(coefs) <- gsub("^.*_", "", colnames(coefs))
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
splenic_feat_colors <- c(rep("Brown", 14), rep("Darkgreen", 18), rep("Red", 24), rep("Purple", 14), rep("Black", 16), rep("Darkturquoise", 16), rep("Blue", 5))

#aics = aics[1:]
#aics = aics[-1,]
data_melt <- melt(as.matrix(coefs))        
print(head(data_melt, 5))                                

ggp <- ggplot(data_melt, aes(X1, X2)) +                           # Create heatmap with ggplot2
  geom_tile(aes(fill = value, linewidth = 1)) +
  scale_fill_gradient2(low="#888BC9", mid = "white", high = "#3A68AE") +
  labs(x = x_label, y = "Splenic Feature", size = 30) +
  theme(axis.text.y = element_text(size=25, colour = splenic_feat_colors), axis.text.x = element_text(size=25, angle=90), legend.text = element_text(size=25), legend.title = element_text(size=30))


ggsave(out, width = width, height = 25, limitsize = FALSE)








