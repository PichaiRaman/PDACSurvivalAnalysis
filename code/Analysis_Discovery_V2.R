###################################
#Code to determine survival of pancreatic cancer patients
#Pichai Raman
#4/25/2018
###################################


#Call Libraries
library("survival");
library("ggplot2");
library("reshape2");
library("stringr");
library("scales");
library("GEOquery");
library("Hmisc");
library("ggthemes");
library("VennDiagram");

#Source code
source("helper/voomLimma.R")
source("helper/DemGGPlotSummary.R")
source("helper/pubTheme.R")

#Create cutoffs
cutoffHigh <- 24
cutoffLow <- 12
options(warn=-1)


#################################
#Format TCGA Data
#################################

#Get Name
getGeneName <- function(x)
{
	x <- strsplit(x, "\\|")[[1]][1]
}

#Master clean function
cleanFormatDisc <- function(annot, exprs)
{
#Get appropriate columns for annot and recode
annot <- annot[,c("bcr_patient_barcode", "days_to_death", "vital_status", "days_to_last_followup", "days_to_birth", "pathologic_stage", "neoplasm_histologic_grade", "gender")];
colnames(annot) <- c("PATIENT_ID", "OS_MONTHS", "OS_STATUS", "FOLLOW_UP", "AGE", "STAGE", "GRADE", "GENDER");
annot[,"OS_STATUS"] <- ifelse(annot[,"OS_STATUS"]=="Alive", 0, 1);
annot[,"FOLLOW_UP"] <- as.numeric(annot[,"FOLLOW_UP"]);
annot[,"FOLLOW_UP"] <- round(annot[,"FOLLOW_UP"]/30)

annot[,"OS_MONTHS"] <- round(as.numeric(annot[,"OS_MONTHS"])/30);
annot[,"OS_MONTHS"] <- paste(annot[,"OS_MONTHS"] , annot[,"FOLLOW_UP"], sep="");
annot[,"OS_MONTHS"] <- as.numeric(as.character(gsub("NA", "",annot[,"OS_MONTHS"]))); 

annot[,"AGE"] <- round((-1)*annot[,"AGE"]/365);
rownames(annot) <- annot[,1];
annot <- annot[-1];

#Format expression
exprs <- unique(exprs);
exprs[,1] <- sapply(exprs[,1], FUN=getGeneName)

#Get max expression and sort and remove duplicates
exprs[,"maxExp"] <- apply(exprs[2:ncol(exprs)], FUN=max, MARGIN=1)
exprs <- exprs[order(-exprs[,"maxExp"]),];
exprs <- exprs[!duplicated(exprs[,1]),]
rownames(exprs) <- exprs[,1] 
exprs <- exprs[-1];

#print histogram of values and establish cutoff
dfForPlot <- data.frame(log10(exprs[,"maxExp"]+1));
colnames(dfForPlot) <- "myVal";
p <- ggplot(dfForPlot, aes(myVal))+geom_histogram(bins=500);
p <- p+geom_vline(xintercept = 2, color="red")
p <- p+theme_Publication()+ggtitle("Histogram of Maximum Value")
p <- p+xlab("Log10 Value")
ggsave("../results/SupplementalFigure_notused.eps", plot=p, width=7, height=7);


#Filter data & remove max value
exprs <- exprs[log10(exprs[,"maxExp"]+1)>3,]

#Rename Columns
colnames(exprs) <- gsub("\\.","-", colnames(exprs));
colnames(exprs) <- substring(colnames(exprs), 1, 12);

#Now get intersection of id's and conserve order
intSamps <- intersect(colnames(exprs), rownames(annot));
output <- list(exprs[,intSamps], annot[intSamps,]);
return(output);

}

annot_pr <- read.delim("../data/DiscoveryDataSet/TCGA_Panc_Annot.txt", stringsAsFactors=F);
annot_pr_old <- annot_pr;
exprs_pr <- read.delim("../data/DiscoveryDataSet/TCGA_Panc_Raw_Exp.txt", row.names=as.character(c(1:20330)))
pr <- cleanFormatDisc(annot_pr, exprs_pr);
saveRDS(pr, "../data/DiscoveryData.rds")

####################
#Write expression object for GSEA later
###################
write.table(data.frame(rownames(pr[[1]])), "../data/GeneUniverse.txt", sep="\t", row.names=F);

###################################################
#Examine people in the same age bracket
###################################################
annot_pr <- pr[[2]];
exprs_pr <- pr[[1]];

myZ <- function(x) { (x-mean(x))/sd(x)}
annot_pr[,"AGE_Z"] <- abs(myZ(annot_pr[,4]))
annot_pr <- annot_pr[annot_pr[,"AGE_Z"]<10,]
exprs_pr <- exprs_pr[,rownames(annot_pr)];

pr[[1]] <- exprs_pr;
pr[[2]] <- annot_pr;
#################################
#Create Groups based on survival time
#################################
annot_pr <- pr[[2]];

#Let's look at a plot of survival time and create logical cutoff's
ecFunc <- ecdf(annot_pr[,1])
annot_pr <- cbind(annot_pr, ecFunc(annot_pr[,1]));

p <- qplot(ecFunc(annot_pr[,1]), annot_pr[,1])
p <- p+geom_hline(yintercept=c(cutoffLow,cutoffHigh), color="red");
p <- p+labs(list(x="Percentile", y="Survival Time", title="Survival Time vs Survival Quantile"))
p <- p+theme_Publication()
p <- p+theme(plot.title = element_text(hjust = .5))
ggsave("../results/SupplementalFigure1.eps", plot=p, width=6, height=6);

lowSamps <- intersect(rownames(annot_pr[annot_pr[,1]<cutoffLow,]), rownames(annot_pr[annot_pr[,2]==1,]));
highSamps <- intersect(rownames(annot_pr[annot_pr[,1]>cutoffHigh,]), rownames(annot_pr[annot_pr[,2]==0,]));

discSampleDF <- rbind(setNames(data.frame("Survival_poor", lowSamps), c("Group", "Sample")), setNames(data.frame("Survival_good", highSamps), c("Group", "Sample")))
write.table(discSampleDF, "../data/discAnalysisSamples.txt", sep="\t", row.names=F);


print(paste("There are ", length(lowSamps), " samples in the low group", sep=""));
print(paste("There are ", length(highSamps), " samples in the high group", sep=""));

#################################
#Purity of groups 
#################################

purityData <- read.delim("../data/PancCellularityStudy/tumCellularity.txt");
rownames(purityData) <- substring(purityData[,1], 1, 12)

lowPurityData <- data.frame(table(purityData[lowSamps,"Purity.Class..high.or.low."]), "Survival-");
lowPurityData[,1] <- as.character(lowPurityData[,1])
lowPurityData <- rbind(lowPurityData, c("Excluded",1, "Survival-"))

highPurityData <- data.frame(table(purityData[highSamps,"Purity.Class..high.or.low."]), "Survival+");
highPurityData[,1] <- as.character(highPurityData[,1])
highPurityData <- rbind(highPurityData, c("Excluded",9, "Survival+"))

PurityDataBar <- rbind(lowPurityData, highPurityData)
PurityDataBar <- PurityDataBar[PurityDataBar[,"Freq"]>0,]
colnames(PurityDataBar) <- c("Purity", "Frequency", "Category");
PurityDataBar[,2] <- as.numeric(PurityDataBar[,2]);

p <- ggplot(PurityDataBar, aes(x = Category, y = Frequency, fill = Purity, label=Frequency));
p <- p+geom_bar(stat = "identity", position = "dodge")
p <- p+geom_text(size = 5, hjust=1.6, position = position_dodge(width=.9))+coord_flip();
ggsave("../results/SupplementalFigure_XX.eps", plot=p, width=6, height=6);


#################################
#Demographics on group 
#################################
annot_pr_demo <- annot_pr[,c(4:7)]
annot_pr_demo[,2] <- as.factor(annot_pr_demo[,2]);
annot_pr_demo[,3] <- as.factor(annot_pr_demo[,3]);
annot_pr_demo[,4] <- as.factor(annot_pr_demo[,4]);
groupsList <- list(highSamps, lowSamps)
demData <- annot_pr_demo
groupsLab <- c("Survival+", "Survival-");

#Update column names to be the right case
colnames(demData)<- capitalize(tolower(colnames(demData)))
p <- SummaryImage(demData, groupsList, groupsLab, "Population Cohort Characteritics")
ggsave("../results/SupplementalFigure2.eps", plot=p, width=7, height=7);


#Let's get stats on everything

#Age - p-value is .39
ks.test(annot_pr_demo[groupsList[[1]],"AGE"], annot_pr_demo[groupsList[[2]],"AGE"])

#Gender - P-vlaue is .34
genderMat <- rbind(as.numeric(table(annot_pr_demo[groupsList[[1]],"GENDER"])), as.numeric(table(annot_pr_demo[groupsList[[2]],"GENDER"])))
fisher.test(genderMat)

#Stage - P-vlaue is .01077
stageMat <- rbind(as.numeric(table(annot_pr_demo[groupsList[[1]],"STAGE"])), as.numeric(table(annot_pr_demo[groupsList[[2]],"STAGE"])))
fisher.test(stageMat)

#Grade - P-vlaue is 0.0006492
gradeMat <- rbind(as.numeric(table(annot_pr_demo[groupsList[[1]],"GRADE"])), as.numeric(table(annot_pr_demo[groupsList[[2]],"GRADE"])))
fisher.test(gradeMat)



##################################################################
#Run Voom Limma to get Large signature Figure 1A
##################################################################

SurvivalDEG <- voomLimma(exprs_pr, list(highSamps, lowSamps))[[2]];

#accessory for volcano plot
reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv, 
              log_breaks(base = base), 
              domain = c(1e-100, Inf))
}

#volcano plot, takes in limma analysis
plotVolcano <- function(result, title="Volcano Plot", otherCol="blue")
{
p <- ggplot(result, aes(x=logFC, y= adj.P.Val, color=Hit))+geom_point()+ scale_y_continuous(trans=reverselog_trans(10))+ggtitle(title)+scale_colour_manual(values =c("gray", otherCol));
p <- p+theme_Publication();
return(p);
}

#Create Data Frame for Volcano Plot
SurvivalDEG[,"absLogFC"] <- abs(SurvivalDEG[,"logFC"]);
SurvivalDEGHits <- unique(rownames(SurvivalDEG[SurvivalDEG[,"absLogFC"]>(.585)&SurvivalDEG[,"adj.P.Val"]<0.05,],));
SurvivalDEG[,"Hit"] <- rownames(SurvivalDEG)%in%SurvivalDEGHits

tmpString <- "Genes up-regulated in Survival- patients: "
print(paste(tmpString, nrow(SurvivalDEG[SurvivalDEG[,"Hit"]==T&SurvivalDEG[,"logFC"]>0,]), sep=""));
tmpString <- "Genes down-regulated in Survival- patients: "
print(paste(tmpString, nrow(SurvivalDEG[SurvivalDEG[,"Hit"]==T&SurvivalDEG[,"logFC"]<0,]), sep=""));

#png("../results/Figure2A.png", width=1440, height=1440,  res=324);
plotVolcano(SurvivalDEG, "Discovery Cohort Volcano Plot", "red");
#dev.off();
ggsave("../results/SupplementalFigure3A.eps", width=8, height=8);



##################################################################
#Run Voom Limma on just the high purity high and low samples
##################################################################

#5 samples with good survival and 8 with poor survival
highSampsPure <- intersect(highSamps, rownames(purityData[purityData[,"Purity.Class..high.or.low."]=="high",]));
lowSampsPure <- intersect(lowSamps, rownames(purityData[purityData[,"Purity.Class..high.or.low."]=="high",]));

SurvivalDEGHighPurity <- voomLimma(exprs_pr, list(highSampsPure, lowSampsPure))[[2]];


#Create Data Frame for Volcano Plot
SurvivalDEGHighPurity[,"absLogFC"] <- abs(SurvivalDEGHighPurity[,"logFC"]);
SurvivalDEGHitsPurity <- unique(rownames(SurvivalDEGHighPurity[SurvivalDEGHighPurity[,"absLogFC"]>(.585)&SurvivalDEGHighPurity[,"P.Value"]<0.05,],));
SurvivalDEGHighPurity[,"Hit"] <- rownames(SurvivalDEGHighPurity)%in%SurvivalDEGHitsPurity

tmpString <- "Genes up-regulated in Survival- High Purity Samples: "
print(paste(tmpString, nrow(SurvivalDEGHighPurity[SurvivalDEGHighPurity[,"Hit"]==T&SurvivalDEGHighPurity[,"logFC"]>0,]), sep=""));
tmpString <- "Genes down-regulated in Survival- High Purity Samples: "
print(paste(tmpString, nrow(SurvivalDEGHighPurity[SurvivalDEGHighPurity[,"Hit"]==T&SurvivalDEGHighPurity[,"logFC"]<0,]), sep=""));

cor.test(SurvivalDEGHighPurity[SurvivalDEGHits,"logFC"], SurvivalDEG[SurvivalDEGHits,"logFC"])
p <- qplot(SurvivalDEGHighPurity[SurvivalDEGHits,"logFC"], SurvivalDEG[SurvivalDEGHits,"logFC"])+geom_point()+geom_smooth(method="lm")+theme_bw();
p <- p+xlab("All Samples - LogFC")+ylab("High Purity Samples - LogFC")+geom_vline(xintercept=0)+geom_hline(yintercept=0)
ggsave("../results/SupplementaryFigure_XZ.eps", width=8, height=8);


##################################################################
#Tumor vs normal pancreatic cancer limma analysis Figure 1B
##################################################################
getGeneSymbol <- function(x)
{

xTmp <- strsplit(x, "\\//");
xTmp <- xTmp[[1]][2];
xTmp <- gsub(" ", "", xTmp);
return(xTmp);

}

cleanFormatValidationSetOne <- function(annot, exprs)
{
#Get appropriate columns for annot and recode
annot[,"TimeVar"] <- as.numeric(gsub("survival_month: ", "", annot[,"characteristics_ch1.1"]));
annot[,"EventVar"] <- as.numeric(gsub("cancer_death: ", "", annot[,"characteristics_ch1.2"]));
annot <- annot[,c("geo_accession", "title", "TimeVar", "EventVar")];
annot[,"group"] <- ifelse(grepl("nontumor", annot[,"title"]), "Normal", "Tumor");
colnames(annot)[1] <- c("Patient");

#Format expression
probeAnnot <- read.delim("../data/GPL6244-24073.txt", stringsAsFactors=F);
probeAnnot <- probeAnnot[,c("ID", "gene_assignment")];
probeAnnot[,2] <- sapply(probeAnnot[,2], FUN=getGeneSymbol);
probeAnnot <- na.omit(probeAnnot);
rownames(probeAnnot) <- probeAnnot[,1];
exprs <- exprs[rownames(probeAnnot),rownames(annot)]
exprs <- cbind(exprs, probeAnnot[2]);

#Filter out 
maxExprs <- apply(exprs[1:90], FUN=max, MARGIN=1)
exprs <- exprs[names(sort(maxExprs, T)),]
exprs <- exprs[!duplicated(exprs[,"gene_assignment"]),]
rownames(exprs) <- exprs[,"gene_assignment"];
exprs <- exprs[1:90];

#Now get intersection of id's and conserve order
intSamps <- intersect(colnames(exprs), rownames(annot));
output <- list(exprs[,intSamps], annot[intSamps,]);
return(output);

}

load("../data/ValidationDataSets/GSE28735.RData");
GSE28735_annot <- pData(dataGSE28735[[1]]);

#Demographic information on this cohort
demoTumNorm <- GSE28735_annot[,c("geo_accession", "source_name_ch1", "characteristics_ch1", "characteristics_ch1.1", "characteristics_ch1.2")]
summary(demoTumNorm)

GSE28735_exprs <- data.frame(exprs(dataGSE28735[[1]]));
pr2 <- cleanFormatValidationSetOne(GSE28735_annot, GSE28735_exprs);

GSE28735_exprs <- pr2[[1]];
GSE28735_annot <- pr2[[2]];

#############################################################
#Run Estimate and get purity of samples
#############################################################
library(estimate)
write.table(GSE28735_exprs, "../results/GSE28735_estimate.txt", sep="\t", row.names=T, quote=F)
filterCommonGenes(input.f="../results/GSE28735_estimate.txt", output.f="../results/GSE28735_10412genes.gct", id="GeneSymbol")
estimateScore("../results/GSE28735_10412genes.gct", "../results/GSE28735_10412_score.gct", platform="affymetrix")

GSE28735_purity <- read.delim("../results/GSE28735_10412_score.gct", sep="\t", skip=2)[2:92];
rownames(GSE28735_purity) <- GSE28735_purity[,1];
GSE28735_purity <- GSE28735_purity[-1];
GSE28735_purity <- data.frame(t(GSE28735_purity));
GSE28735_annot_purity <- cbind(GSE28735_annot, GSE28735_purity[,"TumorPurity"]);
colnames(GSE28735_annot_purity)[6] <- "Purity";
GSE28735_annot_purity <- GSE28735_annot_purity[GSE28735_annot_purity[,"group"]=="Tumor",]

puritySummary <- summary(GSE28735_annot_purity[,"Purity"])
p <- ggplot(GSE28735_annot_purity, aes(Purity, fill=2))+geom_density(alpha=.5)+theme_bw();
p <- p+geom_vline(xintercept=puritySummary[[3]], color="red", size=2, linetype=1)
p <- p+annotate("text", x = .7, y = 3.5, label = paste("Median: ", puritySummary[[3]], sep=""))
p <- p+ggtitle("ESTIMATE Purity Distribution - GSE28735", "center")
p <- p+scale_fill_continuous(guide=FALSE)
ggsave("../results/SupplementaryFigure_XAA.eps", width=8, height=8);

#For Reviewer - Confidence Intervals
a <- mean(GSE28735_annot_purity[,"Purity"])
s <- sd(GSE28735_annot_purity[,"Purity"])
n <- length(GSE28735_annot_purity[,"Purity"])
error <- qnorm(0.975)*s/sqrt(n)
print(a-error)
print(a+error)

purityABS <- na.omit(purityData[,"ABSOLUTE.Purity"])
a <- mean(purityABS)
s <- sd(purityABS)
n <- length(purityABS)
error <- qnorm(0.975)*s/sqrt(n)
print(a-error)
print(a+error)




#############################################################

targ <- GSE28735_annot["group"];
colnames(targ)[1] <- "targ"
exprs_pr_tmp <- GSE28735_exprs[,rownames(targ)]
targ <- targ[,1]
pair <- rep(c(1,2), 45)
design <- model.matrix(~0+targ+pair);
colnames(design) <- gsub("targ", "", colnames(design));
fit <- lmFit(exprs_pr_tmp, design);

cont <- "Tumor-Normal"
contrast.matrix <- eval(as.call(c(as.symbol("makeContrasts"),as.list
(cont),levels=list(design))))

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
tumNormDEGenes <- topTable(fit2, number=40000)
tumNormDEGenes[,"ID"] <- rownames(tumNormDEGenes);
tumNormDEGenes[,"absLogFC"] <- abs(tumNormDEGenes[,"logFC"]);
tumNormDEGenesHits <- unique(tumNormDEGenes[tumNormDEGenes[,"absLogFC"]>(.585)&tumNormDEGenes[,"adj.P.Val"]<0.05,"ID"]);
tumNormDEGenes[,"Hit"] <- rownames(tumNormDEGenes)%in%tumNormDEGenesHits;

tmpString <- "Genes up-regulated in Survival- patients: "
print(paste(nrow(tumNormDEGenes[tumNormDEGenes[,"Hit"]==T&tumNormDEGenes[,"logFC"]>0,]), sep=""));
tmpString <- "Genes down-regulated in Survival- patients: "
print(paste(nrow(tumNormDEGenes[tumNormDEGenes[,"Hit"]==T&tumNormDEGenes[,"logFC"]<0,]), sep=""));


#png("../results/Figure2B.png", width=1440, height=1440,  res=324);
plotVolcano(tumNormDEGenes, "Tumor vs normal Volcano Plot");
#dev.off();
ggsave("../results/SupplementalFigure3B.eps", width=8, height=8);

##################################################################
#Venn Diagram Figure 1C & Scatter Figure 1D
##################################################################

#Figure 1C Venn Diagram
cairo_ps("../results/Figure1_Venn_notused.eps", width=10, height=10);
myVenn <- draw.pairwise.venn(length(tumNormDEGenesHits), length(SurvivalDEGHits), length(intersect(SurvivalDEGHits, tumNormDEGenesHits)), c("Tumor\nvs\nNormal\nAnalysis", "Survival\nAnalysis"),  fill = c("blue", "red"), alpha=.8, cat.pos=c(-10, 0), label.col="black", cex = 2, cat.fontface = 1, cat.cex=2, cat.dist=0.05);
dev.off();

#Figure 1D Scatter plot of tum vs norm DE vs Survival DE
LargeSignature <- intersect(SurvivalDEGHits, tumNormDEGenesHits)
LargeSignatureDF <- cbind(tumNormDEGenes[LargeSignature, c("logFC", "adj.P.Val", "absLogFC")], SurvivalDEG[LargeSignature, c("logFC", "adj.P.Val", "absLogFC")])
colnames(LargeSignatureDF)[1:3] <- paste("TNC_", colnames(LargeSignatureDF)[1:3], sep="");
colnames(LargeSignatureDF)[4:6] <- paste("SRV_", colnames(LargeSignatureDF)[4:6], sep="");
LargeSignatureDF[,"Gene"] <- rownames(LargeSignatureDF);

#png("../results/Figure2D.png", width=1440, height=1440,  res=324);
p <- ggplot(LargeSignatureDF, aes(TNC_logFC, SRV_logFC, size=(-1)*log10(SRV_adj.P.Val)))+geom_point()+xlab("Tumor vs Normal Fold Change")+ylab("Survival- vs Survival+ Fold Change")+scale_colour_manual(values=c("gray", "purple"))
p <- p+ggtitle("Comparison of Survival Analysis & Tumor/Normal Analysis")+guides(size=guide_legend(title = "Fold Change"));
p <- p+theme_Publication();
p
#dev.off();
ggsave("../results/SupplementalFigure3C.eps", width=8, height=8);


#Write out supplemental table 1
LargeSignatureDF <- LargeSignatureDF[LargeSignatureDF[,"SRV_logFC"]*LargeSignatureDF[,"TNC_logFC"]>0,];
write.table(LargeSignatureDF, "../results/SupplementalTable1.txt", sep="\t", row.names=T);

#Get final number of genes up and down
print(paste("Final Signature up-regulated genes: ", sum(LargeSignatureDF[,"SRV_logFC"]>0), sep=""))
print(paste("Final Signature up-regulated genes: ", sum(LargeSignatureDF[,"SRV_logFC"]<0), sep=""))




