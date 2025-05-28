library(ggplot2)
library(scales)
library(cowplot)
#options(bitmapType = 'cairo', device = 'png')

wrap_10 <- wrap_format(15)
# create forest plot for logistic regression only
filename = "../Data/figures/logreg_forest.pdf"
data <- data.frame(study=c('GLSZM Gray Level Non-Uniformity', 'Sphericity', 
	'GLCM Correlation', 'Energy', 'GLDM Small Dependence High Gray Level Emphasis', 
	'GLDM Gray Level Variance', 'GLSZM Large Area Low Gray Level Emphasis', 
	'GLRLM Run Length Non Uniformity', 'GLCM Inverse Difference'), index=1:9, 
	result=c(1.58593272, 1.1660401, 1.298730458, 0.6370196131, 1.2633772, 
		0.7667868262, 1.131855933, 1.322155496, 0.8171318018), 
	error_lower = c(1.382809046, 1.096603106, 1.129639766, 0.5001765358, 
		1.109241642, 0.6501341798, 1.041898178, 1.073549282, 0.7000382351), 
	error_upper=c(1.818893651, 1.23987385, 1.493131576, 0.8113015273, 
		1.438930788, 0.9043702901, 1.229580664, 1.628332471, 0.9538113034))

splenic_feat_colors <- c(rep("Brown", 14), rep("Darkgreen", 18), rep("Red", 24), rep("Purple", 14), rep("Black", 16), rep("Darkturquoise", 16), rep("Blue", 5))
splenic_feat_colors <- c("Darkturquoise", "Brown", "Red", "Darkgreen", "Purple", "Purple", "Darkturquoise", "Black", "Red")
ggplot(data=data, aes(y=index, x=result, xmin=error_lower, xmax=error_upper))+
theme_bw()+
theme(axis.text.y = element_text(size=12, color=splenic_feat_colors), axis.title=element_text(size=14,face="bold"), panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
geom_point(color="#3A68AE", size=3)+geom_errorbarh(height=.2, color="#3A68AE",lwd=1.4)+scale_y_continuous(breaks=1:9, labels=wrap_10(data$study))+labs(x='Odds Ratio (95% CI)', y='Splenic Feature')+geom_vline(xintercept=1, color='red', linetype='dashed', alpha=0.8)
ggsave(filename, device="pdf")


# create forest plot for cox only
filename2 = "../Data/figures/cox_forest.pdf"
data <- data.frame(study=c('GLRLM Run Length Non Uniformity', 
	'GLCM Inverse Difference Normalized'), index=1:2
, 
	result=c(1.17047307, 0.8981348993), 
	error_lower = c(1.091823856, 0.8482625502), 
	error_upper=c(1.254787757, 0.9509394198))

splenic_feat_colors <- c("Black", "Red")

ggplot(data=data, aes(y=index, x=result, xmin=error_lower, xmax=error_upper))+geom_point(color="#3A68AE", size=3)+
geom_errorbarh(height=.2, color="#3A68AE",lwd=1.4)+scale_y_continuous(breaks=1:2, labels=wrap_10(data$study))+
theme_bw()+
theme(axis.text.y = element_text(size=12, color=splenic_feat_colors), axis.title=element_text(size=14,face="bold"), panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
labs(x='Hazard Ratio (95% CI)', y='Splenic Feature')+geom_vline(xintercept=1, color='red', linetype='dashed', alpha=0.8)
ggsave(filename2, device="pdf")

# create one plot together
#data <- data.frame(study=c('Gray Dependence','Gray Non-Uniformity','Gray Emphasis', 'Kurtosis'), index=1:4, result=c(0.8432835663, 1.105827585, 0.8748069898, 0.9102044813), error_lower = c(0.767599689, 1.025450512, 0.7861749568, 0.3742584041), error_upper=c(0.9264297307, 1.192504789, 0.9734312482, 2.213636858))
