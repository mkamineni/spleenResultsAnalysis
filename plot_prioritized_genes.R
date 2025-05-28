library(data.table)
library(ggplot2)
library(dplyr)
library(scales)

file <- "../Data/fuma_results/top_all_phenos_genes_lead_atleast_one.csv"
df <- read.csv(file, header=TRUE)
print(colnames(df))
names(df) <- c('X', 'GLSZM Gray Level Non Uniformity', 'Sphericity', 
	'GLCM Correlation', 'Energy', 'GLDM Small Dependence High Gray Level Emphasis', 'GLSZM Large Area Low Gray Level Emphasis', 
	'GLRLM Run Length Non Uniformity', 'GLCM Inverse Difference', 'GLCM Inverse Difference Normalized')

#df <-df[order(nrow(df):1),]
df<-cbind(df,tag=1:nrow(df))
df_merge<-df[,c("X", "tag")]

print(head(df, 5))
df["Gray Level Variance"] <- "Y"
df["Small Dependence High Gray Level Emphasis"] <- "Y"
df["Small Area Low Gray Level Emphasis"] <- "Y"
df <- df[, c('X', 'tag','GLSZM Gray Level Non Uniformity', 'Sphericity', 
	'GLCM Correlation', 'Energy', 'GLDM Small Dependence High Gray Level Emphasis', 'GLSZM Large Area Low Gray Level Emphasis', 
	'GLRLM Run Length Non Uniformity', 'GLCM Inverse Difference', 'GLCM Inverse Difference Normalized')]
df$X <- NULL
splenic_colors <- c('Darkturquoise', 'Brown', 'Red', 'Darkgreen', 'Purple', 'Red', 'Darkturquoise', 'Black', 'Red', 'Red')

df <- melt(data.table(df), id.vars=c("tag"))
df$value <- as.numeric(df$value == "X")*1
df <- merge(df, df_merge, on="tag", how="right")
df<-df[order(-df$tag),]

df$tag <- NULL
print(head(df,10))

wrap_10 <- wrap_format(20)
options(repr.plot.width = 3, repr.plot.height =24)
df$X = factor(df$X, levels = unique(df$X))
plot = ggplot(df, aes(y=X, x=variable, fill=value))+geom_tile()+labs(y="Gene",x="Splenic Phenotype", size=30)+
theme(axis.text.x = element_text(angle=-90, size=30, color=splenic_colors), axis.text.y=element_text(angle=0, size=30, face="italic"), axis.title=element_text(size=20,face="bold"))+
scale_fill_gradient2(low="#888BC9", high="#3A68AE")+theme(legend.position="none")
ggsave("../Data/figures/prioritized_genes_grid_lead.pdf", height=50, width = 20, device="pdf", limitsize=FALSE)
