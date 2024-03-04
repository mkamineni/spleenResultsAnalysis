library(data.table)
library(ggplot2)
library(dplyr)
library(scales)
library(readxl)

file <- "../Data/figures/SupplementaryTables.xlsx"
df <- read_excel(file, sheet='S33')
names(df) <- c("SNP", "assoc", "Hypertension", "Diabetes", "Systolic Blood Pressure", "Diastolic Blood Pressure", "Smoker", "Total Cholesterol", "HDL", "LDL", "Triglyceride", "BMI/Weight")
print(head(df, 5))
df <- subset(df, select = -c(assoc))
df = df[-(1:3), , drop = FALSE]

df<-cbind(df, tag=1:nrow(df))
df_merge<-df[,c("SNP", "tag")]

print(head(df, 5))
df$SNP <- NULL

df <- melt(data.table(df), id.vars=c("tag"))
df$value <- as.numeric(df$value == "X")*1
df <- df %>% replace(is.na(.), 0)
df <- merge(df, df_merge, on="tag", how="right")
df<-df[order(df$tag),]

df$tag <- NULL
print(head(df,10))

wrap_10 <- wrap_format(20)
df$SNP = factor(df$SNP, levels = unique(df$SNP))
plot = ggplot(df, aes(x=SNP, y=variable, fill=value))+geom_tile()+labs(x="SNP",y="CAD Risk Factor", size=24)+
theme(axis.text.y = element_text(angle=0, size=24), axis.text.x=element_text(angle=-90, size=24, face="italic"), axis.title=element_text(size=20,face="bold"))+
scale_y_discrete(labels = wrap_10(df$variable))+scale_fill_gradient2(high="#3A68AE")+theme(legend.position="none")
ggsave("../Data/figures/overlap_SNP_figure.png", height=20, width = 40, dpi=1000)
