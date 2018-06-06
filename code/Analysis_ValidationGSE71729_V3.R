###################################
#Code to test whether gene ratio of GPC1 to any other membrane protein
#can determine survival of pancreatic cancer patients
#Pichai Raman
#4/1/2016
###################################

#Call Libraries
library("survival");
library("ggplot2");
library("reshape2");
library("limma");
library("AUC")
library("tidyverse");
source("helper/rocon.R")
source("helper/KaplanScan.R")
source("helper/pubTheme.R")


cutoffHigh <- 24
cutoffLow <- 12

#################################
#Load Validation GEO Dataset GSE71729
#################################
cleanFormatValidationSetOne <- function(annot, exprs)
{

#Get appropriate columns for annot and recode & filter
annot <- annot[grep("survival_months", annot[,"characteristics_ch2.2"]),];
annot[,"TimeVar"] <- as.numeric(gsub("survival_months: ", "", annot[,"characteristics_ch2.2"]));
annot[,"EventVar"] <- as.numeric(gsub("death_event_1death_0censor: ", "", annot[,"characteristics_ch2.3"]));
annot <- annot[,c("geo_accession", "TimeVar", "EventVar")];
annot <- na.omit(annot);
rownames(annot) <- annot[,1];

#Now get intersection of id's and conserve order
intSamps <- intersect(colnames(exprs), rownames(annot));
output <- list(exprs[,intSamps], annot[intSamps,]);
return(output);

}

load("../data/ValidationDataSets/GSE71729.RData");
pr <- cleanFormatValidationSetOne(GSE71729_annot, GSE71729_exprs);


#################################
#Normalize based on control genes
#################################
exprs_pr <- pr[[1]];
controlGenes <- c("TUBB", "ACTB", "UBC", "PPIA", "GUSB");
deltaCol <- colMeans(exprs_pr[controlGenes,])
deltaCol <- mean(deltaCol)-deltaCol;

normData <- function(x)
{
	x <- x+deltaCol;
}

exprs_pr <- data.frame(t(apply(exprs_pr, FUN=normData, MARGIN=1)));
pr[[1]] <- exprs_pr;

#################################
#Create Groups based on survival time
#################################
annot_pr <- pr[[2]];

#Let's look at a plot of survival time and create logical cutoff's
ecFunc <- ecdf(annot_pr[,2])
annot_pr <- cbind(annot_pr, ecFunc(annot_pr[,2]));

#pdf("../results/SurvivalDataValidation_GSE21501.pdf", width=3.5, height=5.5);
#p <- qplot(ecFunc(annot_pr[,2]), annot_pr[,2])+theme_bw()
#p <- p+geom_hline(yintercept=c(cutoffLow,cutoffHigh), color="red");
#p <- p+xlab("Percentile")+ylab("Survival Time");
#p;
#dev.off();

lowSamps <- intersect(rownames(annot_pr[annot_pr[,2]<cutoffLow,]), rownames(annot_pr[annot_pr[,3]==1,]));
highSamps <- intersect(rownames(annot_pr[annot_pr[,2]>cutoffHigh,]), rownames(annot_pr[annot_pr[,3]==0,]));

print(paste("There are ", length(lowSamps), " samples in the low group", sep=""));
print(paste("There are ", length(highSamps), " samples in the high group", sep=""));

exprs_pr <- pr[[1]]
############################################
#Classify Validation set 
#############################################

#Let's create normalized expression by columns first
expr_pr_norm <- exprs_pr


#Let's pull out our signature and test it
signatureList <- read.delim("../results/Table1.txt");
rownames(signatureList) <-signatureList[,1];
sigGenes <- unique(signatureList[,1]);

#Let's pull other signatures as well
virg13 <- as.character(read.delim("../data/GeneLists/thirteenGeneSig.txt")[,1]);
virg13Down <- c("MDM2", "TGFA");
virg13Up <- setdiff(virg13, c("MDM2", "TGFA"));

moffitt15 <- as.character(read.delim("../data/GeneLists/fifteenGeneSig.txt")[,1]);
barts36 <- as.character(read.delim("../data/GeneLists/thirtySixGeneSig.txt")[,1]);

bartsGeneDown <- c("ITGBL1","ARRB1","CADPS","EIF4E3","ICOSLG","NOSTRIN","B3GNT1","RPSAP58","CNNM3","QDPR","ZNF471","GTF2IRD2","GTF2IRD2B")
bartsGeneUp <- setdiff(barts36, bartsGeneDown)


#Z-score matrix
myZ <- function(x) {(x-mean(x))/sd(x)}
#Function to create signature based on a list
createSigScore <- function(mySigUp, sigName="NoName", mySigDown=NULL)
{
	#Genes in intersection
	mySig <- union(mySigUp, mySigDown)
	intGenes <- intersect(rownames(exprs_pr), mySig);
	exprs_pr_cmat <- expr_pr_norm[intGenes,];
	exprs_pr_cmat_z <- data.frame(t(apply(exprs_pr_cmat, FUN=myZ, MARGIN=1)));
	#Pull out up and down genes
	discHitsGenesUp <- intersect(mySigUp, intGenes)
	discHitsGenesDown <- intersect(mySigDown, intGenes);
	#Generate Scores
	#myScores <- colSums(exprs_pr_cmat_z[discHitsGenesUp,]*(-1*log10(signatureList[discHitsGenesUp,"P.Value"])))-colSums(exprs_pr_cmat_z[discHitsGenesDown,]*(-1*log10(signatureList[discHitsGenesDown,"P.Value"])))
	myScores <- colSums(exprs_pr_cmat_z[discHitsGenesUp,])-colSums(exprs_pr_cmat_z[discHitsGenesDown,])
	aucDF <- data.frame(myScores);
	aucDF[,"sample"] <- rownames(aucDF);
	classLabs <- data.frame(c(lowSamps, highSamps), c(rep(1, length(lowSamps)), rep(0, length(highSamps))))
	rownames(classLabs) <- classLabs[,1];
	aucDF <- aucDF[rownames(classLabs),]
	aucDF <- cbind(aucDF[1], classLabs[2]);
	colnames(aucDF) <- c("score", "class");
	aucDF[,"Signature"] <- sigName;
	return(list(aucDF, myScores));
}
myGeneSigRes <- createSigScore(sigGenes, sigName="5-Gene Signature")
virg13 <- createSigScore(virg13Up, sigName="Virginia-13", mySigDown=virg13Down)
moffitt15 <- createSigScore(moffitt15, sigName="Moffitt-15")
barts36 <- createSigScore(bartsGeneUp, sigName="Barts-36", mySigDown=bartsGeneDown)
rocDF <- rbind(myGeneSigRes[[1]], virg13[[1]], moffitt15[[1]], barts36[[1]])
roconMult(rocDF, myTitle="ROC Curves for GSE71729");
ggsave("../results/Figure4ARight.eps", width=7, height=7);


#Now let's create distribution of random signatures and put in our signature
myDist <- c();
for(i in 1:5000)
{
randomSig <- createSigScore(sample(rownames(exprs_pr), 5), sigName="randSig")
myDist <- c(rocon(randomSig[[1]])[[3]], myDist);
}

pvalvsRandom <- 1-sum(rocon(myGeneSigRes[[1]])[[3]]>myDist)/length(myDist)
p <- qplot(myDist, bins=200)+geom_vline(xintercept=rocon(myGeneSigRes[[1]])[[3]], color="red")+theme_bw()+xlab("AUC")+ggtitle("AUC Distribution of Random Signatures (GSE71729)")
p <- p+annotate("text", x = .7, y = 150, label = paste("P-value =", round(pvalvsRandom, 6)))
ggsave("../results/SupplementalFigure7Right.eps", width=7, height=7);



#Now let's do kaplan-meier plots
annot_pr <- cbind(annot_pr, myGeneSigRes[[2]]);
colnames(annot_pr)[[5]] <- "myScores";
kapmPlot("myScores", annot_pr, T, perc=0.25, tVar="TimeVar", eVar="EventVar")+ylab("Survival")+xlab("Time");
ggsave("../results/Figure4BRight.eps", width=7, height=7);

###################################################
#Let's see what pathways correlate with the signature
#Reviewer comment 12
###################################################

MYC <- as.character(read.delim("../data/MYC.txt", skip=2)[,1])
MYCActivity <- createSigScore(MYC, sigName="MYC")

EGF <- as.character(read.delim("../data/EGF.txt", skip=2)[,1])
EGFActivity <- createSigScore(EGF, sigName="EGF")

Hypoxia <- as.character(read.delim("../data/Hypoxia.txt", skip=2)[,1])
HypoxiaActivity <- createSigScore(Hypoxia, sigName="Hypoxia")

TGFB <- as.character(read.delim("../data/TGFB.txt", skip=2)[,1])
TGFBActivity <- createSigScore(TGFB, sigName="TGFB")

Autophagy <- as.character(read.delim("../data/Autophagy.txt", skip=2)[,1])
AutophagyActivity <- createSigScore(Autophagy, sigName="Autophagy")

Inflam <- as.character(read.delim("../data/Inflam.txt", skip=2)[,1])
InflamActivity <- createSigScore(Inflam, sigName="Inflammation")


sigScoreCor <- cbind(myGeneSigRes[[1]][1], MYCActivity[[1]][1], EGFActivity[[1]][1], HypoxiaActivity[[1]][1], TGFBActivity[[1]][1], AutophagyActivity[[1]][1], InflamActivity[[1]][1]);
colnames(sigScoreCor) <- c("SignatureScore", "MYC", "EGF", "Hypoxia", "TGFB", "Autophagy", "Inflammation")
sigScoreCor <- data.frame(sigScoreCor)
sigScoreCorTS <- sigScoreCor %>% gather(Pathway, Activity, 2:7)

p <- ggplot(sigScoreCorTS, aes(SignatureScore, Activity, color=Pathway))+geom_point()+geom_smooth(method="lm")+facet_grid(~Pathway)
p <- p+theme_Publication()+xlab("Signature Score")+theme(legend.position = "right", legend.direction = "vertical", legend.key.size= unit(0.6, "cm"))
ggsave("../results/FigureXXY_SigVsPathwayGSE71729.eps", width=14, height=7);

MYCSigLM <- summary(lm(MYC ~ SignatureScore, sigScoreCor))
HypoxiaSigLM <- summary(lm(Hypoxia ~ SignatureScore, sigScoreCor))
EGFSigLM <- summary(lm(EGF ~ SignatureScore, sigScoreCor))
TGFBSigLM <- summary(lm(TGFB ~ SignatureScore, sigScoreCor))
AutoSigLM <- summary(lm(Autophagy ~ SignatureScore, sigScoreCor))
InflamSigLM <- summary(lm(Inflammation ~ SignatureScore, sigScoreCor))

########################################################
#Let's also classify each sample by Bailey and do boxplot
########################################################


baileySig <- read.delim("../data/otherSigs/bailey.txt")

ADEX <- as.character(baileySig[baileySig[,2]=="ADEX",1])
ADEXActivity <- createSigScore(ADEX, sigName="ADEX")[[2]]

Squamous <- as.character(baileySig[baileySig[,2]=="Squamous",1])
SquamousActivity <- createSigScore(Squamous, sigName="Squamous")[[2]]

Pancreatic_Progenitor <- as.character(baileySig[baileySig[,2]=="Pancreatic_Progenitor",1])
Pancreatic_ProgenitorActivity <- createSigScore(Pancreatic_Progenitor, sigName="Pancreatic_Progenitor")[[2]]

Immune <- as.character(baileySig[baileySig[,2]=="Immune",1])
ImmuneActivity <- createSigScore(Immune, sigName="Immune")[[2]]

classDF <- data.frame(ADEXActivity, SquamousActivity, Pancreatic_ProgenitorActivity, ImmuneActivity, myGeneSigRes[[2]])
colnames(classDF)<- c("ADEX", "Squamous", "PP", "Immune", "Signature");


classifySample <- function(x)
{
	x <- as.numeric(x);
	maxClass <- max(x[1:4])
	myClass <- "ADEX";
	if(maxClass==x[2])
	{
		myClass <- "Squamous";
	}
	if(maxClass==x[3])
	{
		myClass <- "PP";
	}
	if(maxClass==x[4])
	{
		myClass <- "Immune";
	}
	return(myClass)

}
classDF[,"ADEX"] <- classDF[,"ADEX"]/length(ADEX)
classDF[,"Squamous"] <- classDF[,"Squamous"]/length(Squamous)
classDF[,"PP"] <- classDF[,"PP"]/length(Pancreatic_Progenitor)
classDF[,"Immune"] <- classDF[,"Immune"]/length(Immune)
classDF[,"BaileyClass"] <- apply(classDF, FUN=classifySample, MARGIN=1)

p <- ggplot(classDF, aes(BaileyClass, Signature, fill=BaileyClass))+geom_boxplot();
p <- p+xlab("Subtype")+ylab("Signature Score")+theme_Publication()+guides(fill=FALSE);
ggsave("../results/FigureXXYZ_ComparisonToBaileyGSE71729.eps", width=7, height=5);
summary(aov(Signature ~ BaileyClass, classDF))
TukeyHSD(aov(Signature ~ BaileyClass, classDF))


