####################################
#Comparison of our list of 707 genes
#to Copy number from TCGA
#Pichai Raman
#4/25/2018 
####################################

#Load libraries
library("ggplot2")
library("cowplot")
library("lumi")
library("limma");
source("helper/pubTheme.R");
source("helper/DemGGPlotSummary.R")


#Our list of genes
discList <- read.delim("../results/supplementalTable1.txt")
discListUp <-rownames(discList[discList[,"SRV_logFC"]>0,])
discListDown <-rownames(discList[discList[,"SRV_logFC"]<0,])

#Read in methlation & format & create M values for statistical analysis & remove duplicate patients
methlyationDat <- read.delim("../data/tcga/data_methylation_hm450.txt")
rownames(methlyationDat) <- methlyationDat[,1];
methlyationDat <- methlyationDat[-1:-2];
colnames(methlyationDat) <- gsub("\\.", "-", colnames(methlyationDat))
methlyationDat <- methlyationDat[,grep("-01", colnames(methlyationDat))]
colnames(methlyationDat) <- substr(colnames(methlyationDat), 1, 12)
methlyationDatKeep <- methlyationDat

methlyationDatMat <- as.matrix(methlyationDat)
methlyationDatMVal <- data.frame(beta2m(methlyationDatMat));
colnames(methlyationDatMVal) <- gsub("\\.", "-", colnames(methlyationDatMVal))

#Read in discovery analysis methylation data
sampComparison <- read.delim("../data/discAnalysisSamples.txt", stringsAsFactors=F);
rownames(sampComparison)<- sampComparison[,2];

#Now format sampComparison & methlationdata for limma analysis
intGenes <- intersect(c(discListUp,discListDown), rownames(methlyationDatMVal))
methlyationDatMVal <- methlyationDatMVal[intGenes,sampComparison[,2]]
methlyationDat <- methlyationDat[intGenes,sampComparison[,2]]

#Now get beta-value delta / choose .2 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4391864/
methlyationDat <- as.matrix(methlyationDat)
methlyationDat[is.na(methlyationDat)]<-.5
methlyationDat <- data.frame(methlyationDat)
colnames(methlyationDatMVal) <- gsub("\\.", "-", colnames(methlyationDatMVal))
betaDeltaDF <- data.frame(rowMeans(methlyationDat[,sampComparison[,1]=="Survival_poor"]), rowMeans(methlyationDat[,sampComparison[,1]=="Survival_good"]))
colnames(betaDeltaDF) <- c("meanBeta_SurvivalPoor", "meanBeta_SurvivalGood");
betaDeltaDF[,"BetaDelta"] <- betaDeltaDF[1]-betaDeltaDF[2]
betaDeltaDF <- betaDeltaDF[order(betaDeltaDF[,3]),]

#Run Limma & get P-values
sampComparison <- sampComparison[-2];
fTarget <- factor(sampComparison[,1]);
design <- model.matrix(~0+fTarget);
fit <- lmFit(methlyationDatMVal, design);
contrast.matrix <- makeContrasts(fTargetSurvival_poor-fTargetSurvival_good, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
output <- topTable(fit2, number=40000)
betaDeltaDF <- betaDeltaDF[rownames(output),]
mAnalysisOut <- data.frame(betaDeltaDF, output[,c("P.Value", "adj.P.Val")])
mAnalysisOut <- cbind(mAnalysisOut, discList[rownames(mAnalysisOut), "SRV_logFC"])
colnames(mAnalysisOut)[6]<- "DiscFC";
mAnalysisOut[,"GeneDiscAssoc"] <- ifelse(mAnalysisOut[,"DiscFC"]>0, "Poor", "Good");

#Look at concordance, looks like 560 are concordant and 116 aren't
sum(mAnalysisOut[,"BetaDelta"]*mAnalysisOut[,"DiscFC"]>0)



#Get's the correlation of a gene between methylation and GE
discAnalysis <- readRDS("../data/DiscoveryData.rds")
getCor <- function(x, retPlot=F)
{
tmpGE <- discAnalysis[[1]][x,]
colnames(methlyationDat) <- gsub("\\.", "-", colnames(methlyationDat)) 
tmpME <- methlyationDat[x,]
tmpIntNames <- intersect(names(tmpGE), names(tmpME));
tmpDFGEME <- data.frame(t(tmpGE[tmpIntNames]), t(tmpME[tmpIntNames]));
out <- cor.test(tmpDFGEME[,1], tmpDFGEME[,2], method="spearman")
myP <- c(out$p.value, out$estimate[[1]])

if(retPlot==T)
{
	tmpDFGEME <- na.omit(tmpDFGEME);
	tmpDFGEME <- cbind(tmpDFGEME, sampComparison[rownames(tmpDFGEME),])
	colnames(tmpDFGEME)<- c("GE", "ME", "Survival");
	tmpDFGEME[,"GE"]<- log2(tmpDFGEME[,"GE"]+1);
	tmpDFGEME[,"PointType"] <- "data";
	
	#Add mean point in per survival
	survGoodMean <- c(mean(tmpDFGEME[tmpDFGEME[,"Survival"]=="Survival_good","GE"]), mean(tmpDFGEME[tmpDFGEME[,"Survival"]=="Survival_good","ME"]), "Survival_good", "mean")
	survPoorMean <- c(mean(tmpDFGEME[tmpDFGEME[,"Survival"]=="Survival_poor","GE"]), mean(tmpDFGEME[tmpDFGEME[,"Survival"]=="Survival_poor","ME"]), "Survival_poor", "mean")
	tmpDFGEME <- rbind(tmpDFGEME, survGoodMean)
	tmpDFGEME <- rbind(tmpDFGEME, survPoorMean)
	tmpDFGEME[,1] <- as.numeric(tmpDFGEME[,1]);
	tmpDFGEME[,2] <- as.numeric(tmpDFGEME[,2]);

	p <- ggplot(tmpDFGEME)+geom_point(aes(GE, ME, color=Survival, size=PointType, shape=PointType))+geom_smooth(aes(GE, ME), method="lm")+theme_Publication();
	p <- p+ggtitle(x)+xlab("Gene Expression")+ylab("Gene Methlation")
	p <- p+guides(size=FALSE, shape=FALSE)

	if(x=="MET")
	{
	#plot(p)
	#ggsave(paste("../results/Figure2B.eps", sep=""), width=7, height=7);
	}
	if(x=="PTPRN")
	{
	#plot(p)
	#ggsave(paste("../results/Figure2C.eps", sep=""), width=7, height=7);
	}
	#plot(p)
	#ggsave(paste("../results/ExpVsMethPlots/", x, ".eps", sep=""), width=7, height=7);
}
return(myP)

}
corDat <- data.frame(t(sapply(rownames(mAnalysisOut), FUN=getCor)))
colnames(corDat) <- c("CorGEME_Pval", "CorGEME_Est")
mAnalysisOut <- cbind(mAnalysisOut, corDat);
mAnalysisOut[,"Hit"] <- mAnalysisOut[,"adj.P.Val"]<0.01&abs(mAnalysisOut[,"BetaDelta"])>.2&mAnalysisOut[,"BetaDelta"]*mAnalysisOut[,"DiscFC"]<0&mAnalysisOut[,"CorGEME_Pval"]<0.01

#Print out table and a graph(s)
p <- ggplot(mAnalysisOut, aes(BetaDelta, DiscFC, color=mAnalysisOut[,"BetaDelta"]*mAnalysisOut[,"DiscFC"]<0))+geom_point()+theme_Publication()
p <- p+geom_vline(xintercept = 0)+geom_hline(yintercept = 0)
p <- p+xlab("Beta Delta")+ylab("Discovery Log FC")+guides(color=guide_legend(title="Anti-correlated"))
p;
ggsave("../results/Figure2B.eps", plot=p, width=7, height=7);
mAnalysisOutHits <- mAnalysisOut[mAnalysisOut[,"Hit"]==T,]
sapply(rownames(mAnalysisOutHits), FUN=getCor, retPlot=T)
write.table(mAnalysisOutHits[1:7], "../results/SupplementalTable3.txt", sep="\t", row.names=T)












