library(data.table)
library(ggplot2)
library(dplyr)
library(scales)

file <- "../Data/fuma_results/top_all_phenos_genes_lead_atleast_one.csv"
df <- read.csv(file, header=TRUE)
print(colnames(df))
names(df) <- c('X', 'Gray Level Non Uniformity', 'Sphericity', 
	'Correlation', 'Energy', 'Small Dependence High Gray Level Emphasis', 'Large Area Low Gray Level Emphasis', 
	'Run Length Non Uniformity', 'Inverse Difference', '90th Percentile', 'Max 2D Diameter Slice', 
	'Large Dependence High Gray Level Emphasis', 'Small Area Low Gray Level Emphasis', 
	'Flatness')
#df <-df[order(nrow(df):1),]
df<-cbind(df,tag=1:nrow(df))
df_merge<-df[,c("X", "tag")]

print(head(df, 5))
df["Gray Level Variance"] <- "Y"
df["Small Dependence High Gray Level Emphasis"] <- "Y"
df["Small Area Low Gray Level Emphasis"] <- "Y"
df <- df[, c('X', 'tag','Gray Level Non Uniformity', 'Sphericity', 
	'Correlation', 'Energy', 'Small Dependence High Gray Level Emphasis', 'Gray Level Variance', 'Large Area Low Gray Level Emphasis',
	'Run Length Non Uniformity', 'Inverse Difference', '90th Percentile', 'Max 2D Diameter Slice', 
	'Large Dependence High Gray Level Emphasis', 'Small Area Low Gray Level Emphasis', 
	'Flatness')]
df$X <- NULL
splenic_colors <- c('Darkturquoise', 'Brown', 'Red', 'Darkgreen', 'Purple', 'Red', 'Darkturquoise', 'Black', 'Red', 'Darkgreen', 'Brown', 'Purple', 'Darkturquoise', 'Brown')

df <- melt(data.table(df), id.vars=c("tag"))
df$value <- as.numeric(df$value == "X")*1
df <- merge(df, df_merge, on="tag", how="right")
df<-df[order(df$tag),]

df$tag <- NULL
print(head(df,10))

wrap_10 <- wrap_format(20)
options(repr.plot.width = 3, repr.plot.height =24)
df$X = factor(df$X, levels = unique(df$X))
plot = ggplot(df, aes(x=X, y=variable, fill=value))+geom_tile()+labs(x="Gene",y="Splenic Phenotype", size=24)+
theme(axis.text.y = element_text(angle=0, size=24, color=splenic_colors), axis.text.x=element_text(angle=-90, size=24, face="italic"), axis.title=element_text(size=20,face="bold"))+
scale_y_discrete(labels = wrap_10(df$variable))+scale_fill_gradient2(low="#888BC9", high="#3A68AE")+theme(legend.position="none")
ggsave("../Data/figures/prioritized_genes_grid_lead.png", height=20, width = 40, dpi=1000)
