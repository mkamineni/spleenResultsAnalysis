library(tidyverse)
library(data.table)
library(ggplot2)
library(dplyr)
library(janitor)
library(cowplot)
options(repr.plot.width=50, repr.plot.height=15)

trait = "sex"
file <- "../Data/radiomics_data/radspleen_wo_heme_cancer.csv"
file2 <- "../Data/radiomics_data/CAD_pat_char.csv"
enroll_time <- "../Data/radiomics_data/CAD_coh_time.csv"
mri_time <- "../Data/radiomics_data/CAD_coh_MRI_time.csv"
time_out <- "../Data/radiomics_data/CAD_coh_days_to_MRI.csv"

if (trait=="BMI"){
  cols = c("BMI")
	out <- "../Data/figures/BMI_rad_forest_plot.png"
} else if (trait=="sex") {
  cols = c("BMI")
  out <- "../Data/figures/sex_rad_forest_plot.png"
} else if (trait == "adiposity") {
  cols = c("BMI", "vat", "asat", "WHRadjbmiInstance2")
  out <- "../Data/figures/adiposity_rad_forest_plot.png"
} else if (trait == "blood") {
  cols = c("C.reactive.protein", "Basophill.count", "Basophill.percentage", "Eosinophill.count", "Eosinophill.percentage", "Haematocrit.percentage", "Haemoglobin.concentration", "High.light.scatter.reticulocyte.count",
  "High.light.scatter.reticulocyte.percentage", "Immature.reticulocyte.fraction", "Lymphocyte.count", "Lymphocyte.percentage", "Mean.corpuscular.haemoglobin", "Mean.corpuscular.haemoglobin.concentration", "Mean.corpuscular.volume",
  "Mean.platelet..thrombocyte..volume", "Mean.reticulocyte.volume", "Mean.sphered.cell.volume", "Monocyte.count", "Monocyte.percentage", "Neutrophill.count", "Neutrophill.percentage", "Nucleated.red.blood.cell.count",
  "Nucleated.red.blood.cell.percentage", "Platelet.count", "Platelet.crit", "Platelet.distribution.width", "Red.blood.cell..erythrocyte..count", "Red.blood.cell..erythrocyte..distribution.width", "Reticulocyte.count",
  "Reticulocyte.percentage", "White.blood.cell..leukocyte..count")
  out <- "../Data/figures/blood_rad_forest_plot.png"
}

cols <- append(cols, "ID")

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
pat_char=left_join(pat_char, times, by="ID", all.x=FALSE, mult = "first")

# merge coh with patient characteristics
final_coh <- left_join(coh, pat_char, by = "ID", all.x=FALSE, mult = "first")
print(nrow(final_coh))
#final_coh<-final_coh[complete.cases(final_coh$BMI),]
print(nrow(final_coh))

print(median(final_coh$days_to_MRI))
print(IQR(final_coh$days_to_MRI))

vec <- names(final_coh)[startsWith(names(final_coh), "spleen_")]

columns = c("Feature", "Coefficient", "Low_95", "Upper_95")
coefs = data.frame(matrix(nrow=1, ncol = length(columns))) 
colnames(coefs) = columns

aic_columns = c("Feature", "AIC")
aics = data.frame(matrix(nrow=1, ncol = length(aic_columns))) 
colnames(aics) = aic_columns

for (splenic_feat in vec) {
    final_coh["y"] = scale(final_coh[splenic_feat])   #gives mean of 0 and sd of 1
    if (trait == "BMI" | trait == "sex") {
      m1 <- glm(y ~ BMI+age+sex_Female+days_to_MRI, data = final_coh, family = gaussian())
      if (trait == "BMI") {
        coef = summary(m1)$coefficients["BMI", "Estimate"]
        std = summary(m1)$coefficients["BMI", "Std. Error"]
      } else {
        coef = summary(m1)$coefficients["sex_Female", "Estimate"]
        std = summary(m1)$coefficients["sex_Female", "Std. Error"]
      }
      new_coef = c(splenic_feat, coef, coef-1.96*std, coef+1.96*std)
      names(new_coef) = columns
      coefs <- rbind(coefs, new_coef)
    } else {
      m1 <- glm(reformulate(termlabels = c(cols, c("age", "sex_Female", "days_to_MRI")), response='y'), data = final_coh, family = gaussian())
      new_aic = c(splenic_feat, with(summary(m1), 1 - deviance/null.deviance))
      names(new_aic) = aic_columns
      aics <- rbind(aics, new_aic)
    }

}

coefs$Feature <- gsub("spleen_original_", "", coefs$Feature)
coefs$Feature <- gsub("shape_", "", coefs$Feature)
coefs$Feature <- gsub("firstorder_", "", coefs$Feature)
coefs$Feature <- gsub("glcm_", "", coefs$Feature)
coefs$Feature <- gsub("gldm_", "", coefs$Feature)
coefs$Feature <- gsub("ngtdm_Strength", "Strength", coefs$Feature)
coefs$Feature <- gsub("ngtdm_Complexity", "Complexity", coefs$Feature)
coefs$Feature <- gsub("ngtdm_Busyness", "Busyness", coefs$Feature)
coefs$Feature <- gsub("ngtdm_Coareness", "Coarseness", coefs$Feature)

coefs$Feature <- gsub("glszm_ZoneVariance", "ZoneVariance", coefs$Feature)
coefs$Feature <- gsub("glszm_ZonePercentage", "ZonePercentage", coefs$Feature)
coefs$Feature <- gsub("glszm_ZoneEntropy", "ZoneEntropy", coefs$Feature)
coefs$Feature <- gsub("glszm_SmallAreaLowGrayLevelEmphasis", "SmallAreaLowGrayLevelEmphasis", coefs$Feature)
coefs$Feature <- gsub("glszm_SmallAreaHighGrayLevelEmphasis", "SmallAreaHighGrayLevelEmphasis", coefs$Feature)
coefs$Feature <- gsub("glszm_SmallAreaEmphasis", "SmallAreaEmphasis", coefs$Feature)
coefs$Feature <- gsub("glszm_SizeZoneNonUniformityNormalized", "SizeZoneNonUniformityNormalized", coefs$Feature)
coefs$Feature <- gsub("glszm_LowGrayLevelZoneEmphasis", "LowGrayLevelZoneEmphasis", coefs$Feature)
coefs$Feature <- gsub("glszm_LargeAreaLowGrayLevelEmphasis", "LargeAreaLowGrayLevelEmphasis", coefs$Feature)
coefs$Feature <- gsub("glszm_LargeAreaHighGrayLevelEmphasis", "LargeAreaHighGrayLevelEmphasis", coefs$Feature)
coefs$Feature <- gsub("glszm_LargeAreaEmphasis", "LargeAreaEmphasis", coefs$Feature)
coefs$Feature <- gsub("glszm_GrayLevelZoneEmphasis", "GrayLevelZoneEmphasis", coefs$Feature)

coefs$Feature <- gsub("glrlm_ShortRunLowGrayLevelEmphasis", "ShortRunLowGrayLevelEmphasis", coefs$Feature)
coefs$Feature <- gsub("glrlm_ShortRunHighGrayLevelEmphasis", "ShortRunHighGrayLevelEmphasis", coefs$Feature)
coefs$Feature <- gsub("glrlm_ShortRunEmphasis", "ShortRunEmphasis", coefs$Feature)
coefs$Feature <- gsub("glrlm_RunVariance", "RunVariance", coefs$Feature)
coefs$Feature <- gsub("glrlm_RunPercentage", "RunPercentage", coefs$Feature)
coefs$Feature <- gsub("glrlm_RunLengthNonUniformity", "RunLengthNonUniformity", coefs$Feature)
coefs$Feature <- gsub("glrlm_RunLengthNonUniformityNormalized", "RunLengthNonUniformityNormalized", coefs$Feature)
coefs$Feature <- gsub("glrlm_RunEntropy", "RunEntropy", coefs$Feature)
coefs$Feature <- gsub("glrlm_LowGrayLevelRunEmphasis", "LowGrayLevelRunEmphasis", coefs$Feature)
coefs$Feature <- gsub("glrlm_LongRunLowGrayLevelEmphasis", "LongRunLowGrayLevelEmphasis", coefs$Feature)
coefs$Feature <- gsub("glrlm_LongRunHighGrayLevelEmphasis", "LongRunHighGrayLevelEmphasis", coefs$Feature)
coefs$Feature <- gsub("glrlm_LongRunEmphasis", "LongRunEmphasis", coefs$Feature)
coefs$Feature <- gsub("glrlm_HighGrayLevelRunEmphasis", "HighGrayLevelRunEmphasis", coefs$Feature)
print(head(coefs, 10))

splenic_feat_colors <- c(rep("Brown", 14), rep("Darkgreen", 18), rep("Red", 24), rep("Purple", 14), rep("Black", 16), rep("Darkturquoise", 16), rep("Blue", 5))

if (trait == "BMI" | trait == "sex") {
  coefs = coefs[-1,]
  #print(head(coefs, 10))
  coefs$Color = splenic_feat_colors

  coefs$Coefficient = as.numeric(coefs$Coefficient)
  coefs =  coefs[order(coefs$Coefficient, decreasing = FALSE),]
  print(head(coefs, 130))

  coefs$Feature <- factor(coefs$Feature, levels = coefs$Feature)

  coef_plot <- ggplot(coefs, aes(x = as.numeric(Coefficient), y=Feature)) +
    geom_point(shape = 15, size = 2, color="#3A68AE") +  
    geom_errorbarh(aes(xmin = as.numeric(Low_95), xmax = as.numeric(Upper_95)), color="#3A68AE", height = 0.25) +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
    #scale_y_continuous(name = "", breaks=1:length(bmi_coefs$Feature), labels = bmi_coefs$Feature, trans = "reverse") +
    xlab("Beta (95% CI)") + 
    ylab("Splenic Radiomics Feature") +
    theme_bw() +
    theme(panel.border = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text.y = element_text(size = 5, colour = coefs$Color),
          axis.text.x.bottom = element_text(size = 8, colour = "black"),
          axis.title.x = element_text(size = 12, colour = "black"))
  ggsave(out)
} else {
  #aics = aics[1:]
  aics = aics[-1,]
  print(head(aics, 5))

  aics$AIC = as.numeric(aics$AIC)
  aics =  aics[order(aics$AIC, decreasing = FALSE),]
  print(head(aics, 10))

  print("before factor")
  aics$Feature <- factor(aics$Feature, levels = aics$Feature)
  print("after factor")

  aic_plot <- ggplot(aics, aes(x = as.numeric(AIC), y=Feature)) +
    geom_point(shape = 15, size = 2) +  
    xlab("McFadden's R2") + 
    ylab("Splenic Radiomics Feature") +
    ggtitle(trait) +
    theme_bw() +
    theme(panel.border = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text.y = element_text(size = 4, colour = "black"),
          axis.text.x.bottom = element_text(size = 8, colour = "black"),
          axis.title.x = element_text(size = 12, colour = "black"))
  ggsave(out)

}







