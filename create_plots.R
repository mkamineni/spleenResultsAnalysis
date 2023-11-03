library(ggplot2)
options(bitmapType = 'cairo', device = 'png')

# create forest plot for logistic regression only
filename = "../Data/figures/logreg_forest.png"
data <- data.frame(study=c('Gray Non-Uniformity','Gray Inverse Diff','Sphericity'), index=1:3, result=c(1.259904484, 0.8248733691, 1.110503303), error_lower = c(1.193233502, 0.7791491039, 1.048418274), error_upper=c(1.330300655, 0.8732809569, 1.176264872))

ggplot(data=data, aes(y=index, x=result, xmin=error_lower, xmax=error_upper))+theme(axis.text.y = element_text(size=12), axis.title=element_text(size=14,face="bold"))+geom_point(color="blue", size=3)+geom_errorbarh(height=.2, color="blue",lwd=1.4)+scale_y_continuous(breaks=c(1.0, 2.0, 3.0), labels=data$study)+labs(x='Odds Ratio (95% CI)', y='Splenic Feature')+geom_vline(xintercept=1, color='red', linetype='dashed', alpha=0.8)
ggsave(filename)


# create forest plot for cox only
filename2 = "../Data/figures/cox_forest.png"
data <- data.frame(study=c('Gray Dependence','Gray Non-Uniformity','Gray Emphasis', 'Kurtosis'), index=1:4, result=c(0.8432835663, 1.105827585, 0.8748069898, 0.9102044813), error_lower = c(0.767599689, 1.025450512, 0.7861749568, 0.8328028995), error_upper=c(0.9264297307, 1.192504789, 0.9734312482, 0.9947998478))

ggplot(data=data, aes(y=index, x=result, xmin=error_lower, xmax=error_upper))+geom_point(color="blue", size=3)+geom_errorbarh(height=.2, color="blue",lwd=1.4)+scale_y_continuous(breaks=c(1.0, 2.0, 3.0, 4.0), labels=data$study)+theme(axis.text.y = element_text(size=12), axis.title=element_text(size=14,face="bold"))+labs(x='Hazard Ratio (95% CI)', y='Splenic Feature')+geom_vline(xintercept=1, color='red', linetype='dashed', alpha=0.8)
ggsave(filename2)

# create one plot together
#data <- data.frame(study=c('Gray Dependence','Gray Non-Uniformity','Gray Emphasis', 'Kurtosis'), index=1:4, result=c(0.8432835663, 1.105827585, 0.8748069898, 0.9102044813), error_lower = c(0.767599689, 1.025450512, 0.7861749568, 0.3742584041), error_upper=c(0.9264297307, 1.192504789, 0.9734312482, 2.213636858))
