library(data.table)
library(ggplot2)
library(dplyr)

file <- "../Data/fuma_results/top_all_phenos_genes_lead.csv"
df <- read.csv(file, header=TRUE)
names(df) <- c("X", "Sphericity", "Gray Non-Uniformity", "Kurtosis", "Gray Dependence", "Gray Emphasis", "Gray Inverse Diff")
#df <-df[order(nrow(df):1),]
df<-cbind(df,tag=1:nrow(df))
df_merge<-df[,c("X", "tag")]

print(head(df, 5))
df["Gray Emphasis"] <- "Y"
df$X <- NULL

df <- melt(data.table(df), id.vars=c("tag"))
df$value <- as.numeric(df$value == "X")*1
df <- merge(df, df_merge, on="tag", how="right")
df<-df[order(df$tag),]

df$tag <- NULL
print(head(df,10))

options(repr.plot.width = 3, repr.plot.height =24)
plot = ggplot(df, aes(x=variable, y=X, fill=value))+geom_tile()+labs(y="Gene",x="Splenic Phenotype", size=16)+theme(axis.text.x = element_text(angle=-45, size=16), axis.text.y=element_text(size=16, face="italic"), axis.title=element_text(size=20,face="bold"))+scale_x_discrete(position = "top")+scale_fill_gradient2()
ggsave("../Data/figures/prioritized_genes_grid_lead.png", height=14, dpi=1000)
