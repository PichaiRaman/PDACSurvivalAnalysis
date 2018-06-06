###############################################
#Code to get CCLE data & classify
#Pichai Raman
#4/26/2018
###############################################

#Call libraries
library(tidyr)
library(ggplot2)
source("../code/helper/pubTheme.R")
#Get data
data <- read.delim("../data/IHCData.txt");
data <- data[,c(1:6,8:14)];

#Tumor Differentiation vs Mean score
p <- ggplot(data, aes(Tumor.differentiation, Mean.Score))+geom_boxplot()+theme_Publication()+ theme(axis.text.x = element_text(angle = 90, hjust = 1));
p <- p+ylab("Mean Score")+xlab("Tumor Differentiation")+ggtitle("Tumor Differentiation vs Mean Score")
ggsave("../results/FigureMeanScore_vs_TumorDiff.eps", width=7, height=7);
t.test(Mean.Score ~ Tumor.differentiation, data)
#tumDifAnalysis <- aov(Mean.Score ~ Tumor.differentiation, data)
#summary(tumDifAnalysis)
#TukeyHSD(tumDifAnalysis)

#Vascular Invasion vs Mean score
p <- ggplot(data, aes(Vascual.invasion, Mean.Score))+geom_boxplot()+theme_Publication()+ theme(axis.text.x = element_text(angle = 90, hjust = 1));
p <- p+ylab("Mean Score")+xlab("Vascular Invasion")+ggtitle("Vascular Invasion vs Mean Score")
ggsave("../results/Figure3ALeft.eps", width=7, height=7);
VascInvAnalysis <- t.test(Mean.Score ~ Vascual.invasion, data)

#Perineural Invasion vs Mean score
p <- ggplot(data, aes(Perineural.invasion, Mean.Score))+geom_boxplot()+theme_Publication()+ theme(axis.text.x = element_text(angle = 90, hjust = 1));
p <- p+ylab("Mean Score")+xlab("Perineural Invasion")+ggtitle("Perineural Invasion vs Mean Score")
ggsave("../results/Figure3AMiddle.eps", width=7, height=7);
PerInvAnalysis <- t.test(Mean.Score~ Perineural.invasion, data)

#Stage vs Mean score
p <- ggplot(data, aes(AJCC.stage, Mean.Score))+geom_boxplot()+theme_Publication()+ theme(axis.text.x = element_text(angle = 90, hjust = 1));
p <- p+ylab("Mean Score")+xlab("AJCC Stage")+ggtitle("AJCC Stage vs Mean Score")
ggsave("../results/Figure3ARight.eps", width=7, height=7);
stageAnalysis <- aov(Mean.Score~ AJCC.stage, data)
summary(stageAnalysis)
TukeyHSD(stageAnalysis)

#Power analysis for Reviewer
d=.8 because large effect size

#Large effect size
pwr.t.test(d=.8, sig.level=0.05, power=.8)

#medium effect size
pwr.t.test(d=.5, sig.level=0.05, power=.8)






