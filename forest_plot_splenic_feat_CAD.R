library(ggplot2)
library(scales)
library(cowplot)
#options(bitmapType = 'cairo', device = 'png')

wrap_10 <- wrap_format(15)
# create forest plot for logistic regression only
filename = "../Data/figures/logreg_forest.png"
data <- data.frame(study=c('GLSZM Gray Level Non-Uniformity', 'Sphericity', 
	'GLCM Correlation', 'Energy', 'GLDM Small Dependence High Gray Level Emphasis', 
	'GLDM Gray Level Variance', 'GLSZM Large Area Low Gray Level Emphasis', 
	'GLRLM Run Length Non Uniformity', 'GLCM Inverse Difference'), index=1:9, 
	result=c(1.58587305, 1.160977351, 1.301599008, 0.6400165257, 1.260875935, 
		0.7679382134, 1.133504941, 1.32717505, 0.8192561588), 
	error_lower = c(1.382777282, 1.091514919, 1.132132293, 0.5022430977, 
		1.107009105, 0.6512261324, 1.043332509, 1.077097269, 0.7016542747), 
	error_upper=c(1.818798561, 1.234860272, 1.496432871, 0.8155834396, 
		1.436129221, 0.905567314, 1.231470734, 1.635315271, 0.9565688944))

splenic_feat_colors <- c(rep("Brown", 14), rep("Darkgreen", 18), rep("Red", 24), rep("Purple", 14), rep("Black", 16), rep("Darkturquoise", 16), rep("Blue", 5))
splenic_feat_colors <- c("Darkturquoise", "Brown", "Red", "Darkgreen", "Purple", "Purple", "Darkturquoise", "Black", "Red")
ggplot(data=data, aes(y=index, x=result, xmin=error_lower, xmax=error_upper))+
theme_bw()+
theme(axis.text.y = element_text(size=12, color=splenic_feat_colors), axis.title=element_text(size=14,face="bold"), panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
geom_point(color="#3A68AE", size=3)+geom_errorbarh(height=.2, color="#3A68AE",lwd=1.4)+scale_y_continuous(breaks=1:9, labels=wrap_10(data$study))+labs(x='Odds Ratio (95% CI)', y='Splenic Feature')+geom_vline(xintercept=1, color='red', linetype='dashed', alpha=0.8)
ggsave(filename)


# create forest plot for cox only
filename2 = "../Data/figures/cox_forest.png"
data <- data.frame(study=c('90th Percentile', 'Max 2D Diameter Slice', 
	'GLDM Large Dependence High Gray Level Emphasis', 'GLSZM Small Area Low Gray Level Emphasis', 
	'Flatness'), index=1:5, 
	result=c(1.087198278, 1.064758486, 0.8639244121, 0.8606718274, 0.9243014647), 
	error_lower = c(1.007803853, 1.015806136, 0.7874275756, 0.7687594248, 0.8572748672), 
	error_upper=c(1.172847368, 1.116069882, 0.947852746, 0.9635732202, 0.9965685808))

splenic_feat_colors <- c("Darkgreen", "Brown", "Purple", "Darkturquoise", "Brown")

ggplot(data=data, aes(y=index, x=result, xmin=error_lower, xmax=error_upper))+geom_point(color="#3A68AE", size=3)+
geom_errorbarh(height=.2, color="#3A68AE",lwd=1.4)+scale_y_continuous(breaks=1:5, labels=wrap_10(data$study))+
theme_bw()+
theme(axis.text.y = element_text(size=12, color=splenic_feat_colors), axis.title=element_text(size=14,face="bold"), panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
labs(x='Hazard Ratio (95% CI)', y='Splenic Feature')+geom_vline(xintercept=1, color='red', linetype='dashed', alpha=0.8)
ggsave(filename2)

# create one plot together
#data <- data.frame(study=c('Gray Dependence','Gray Non-Uniformity','Gray Emphasis', 'Kurtosis'), index=1:4, result=c(0.8432835663, 1.105827585, 0.8748069898, 0.9102044813), error_lower = c(0.767599689, 1.025450512, 0.7861749568, 0.3742584041), error_upper=c(0.9264297307, 1.192504789, 0.9734312482, 2.213636858))
