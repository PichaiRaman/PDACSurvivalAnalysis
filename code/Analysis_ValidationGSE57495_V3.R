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
library("ggplot2");
library("limma");
library("tidyverse");
source("helper/rocon.R")
source("helper/KaplanScan.R")
source("helper/pubTheme.R")


cutoffHigh <- 24
cutoffLow <- 12


#################################
#Load Validation GEO Dataset GSE57495
#################################
cleanFormatValidationSetOne <- function(annot, exprs)
{
#Get appropriate columns for annot and recode
annot[,"characteristics_ch1.1"] <- as.character(annot[,"characteristics_ch1.1"]);
annot[,"characteristics_ch1"] <- as.character(annot[,"characteristics_ch1"]);
annot[,"TimeVar"] <- as.numeric(gsub("overall survival \\(\\month\\)\\: ", "", annot[,"characteristics_ch1"]));
annot[,"EventVar"] <- ifelse(gsub("vital\\.status\\: ", "", annot[,"characteristics_ch1.1"])=="DEAD", 1, 0);
annot <- annot[,c("geo_accession", "TimeVar", "EventVar")];
colnames(annot)[1] <- c("Patient");
annot <- na.omit(annot);

#Format expression
probeAnnot <- read.delim("../data/Affy2_custom.txt");
probeAnnot <- probeAnnot[,c("ID", "GeneSymbol")];
exprs <- exprs[probeAnnot[,1],]
exprs <- cbind(exprs, probeAnnot[2]);
exprs <- exprs[!exprs[,"GeneSymbol"]=="",]
maxExprs <- apply(exprs[1:63], FUN=max, MARGIN=1)
exprs <- exprs[names(sort(maxExprs, T)),]
exprs <- exprs[!duplicated(exprs[,"GeneSymbol"]),]
rownames(exprs) <- exprs[,"GeneSymbol"];
exprs <- exprs[1:63];

#Now get intersection of id's and conserve order
intSamps <- intersect(colnames(exprs), rownames(annot));
output <- list(exprs[,intSamps], annot[intSamps,]);
return(output);

}

load("../data/ValidationDataSets/GSE57495.RData");
pr <- cleanFormatValidationSetOne(GSE57495_annot, GSE57495_exprs);



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

#png("../results/SurvivalDataValidation_GSE57495.png", width=1440, height=1440,  res=324);
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

#############################################
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
roconMult(rocDF, myTitle="ROC Curves for GSE57495");
ggsave("../results/Figure4AMiddle.eps", width=7, height=7);





#Now let's create distribution of random signatures and put in our signature
myDist <- c();
for(i in 1:5000)
{
randomSig <- createSigScore(sample(rownames(exprs_pr), 5), sigName="randSig")
myDist <- c(rocon(randomSig[[1]])[[3]], myDist);
}

pvalvsRandom <- 1-sum(rocon(myGeneSigRes[[1]])[[3]]>myDist)/length(myDist)
p <- qplot(myDist, bins=200)+geom_vline(xintercept=rocon(myGeneSigRes[[1]])[[3]], color="red")+theme_bw()+xlab("AUC")+ggtitle("AUC Distribution of Random Signatures (GSE57495)")
p <- p+annotate("text", x = .7, y = 150, label = paste("P-value =", round(pvalvsRandom, 6)))
ggsave("../results/SupplementalFigure7.eps", width=7, height=7);


#Now let's do kaplan-meier plots
annot_pr <- cbind(annot_pr, myGeneSigRes[[2]]);
colnames(annot_pr)[[5]] <- "myScores";
kapmPlot("myScores", annot_pr, T, perc=0.25, tVar="TimeVar", eVar="EventVar")+ylab("Survival")+xlab("Time");
ggsave("../results/Figure4BMiddle.eps", width=7, height=7);



########################################################
#Let's also classify each sample by Bailey and do boxplot
########################################################


baileySig <- read.delim("../data/otherSigs/bailey.txt")

ADEX <- as.character(baileySig[baileySig[,2]=="ADEX",1])
ADEXActivity <- createSigScore(ADEX, sigName="ADEX")[2]

Squamous <- as.character(baileySig[baileySig[,2]=="Squamous",1])
SquamousActivity <- createSigScore(Squamous, sigName="Squamous")[2]

Pancreatic_Progenitor <- as.character(baileySig[baileySig[,2]=="Pancreatic_Progenitor",1])
Pancreatic_ProgenitorActivity <- createSigScore(Pancreatic_Progenitor, sigName="Pancreatic_Progenitor")[2]

Immune <- as.character(baileySig[baileySig[,2]=="Immune",1])
ImmuneActivity <- createSigScore(Immune, sigName="Immune")[2]

classDF <- data.frame(ADEXActivity, SquamousActivity, Pancreatic_ProgenitorActivity, ImmuneActivity, myGeneSigRes[2])
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
ggsave("../results/FigureXXYZ_ComparisonToBaileyGSE57495.eps", width=7, height=5);

summary(aov(Signature ~ BaileyClass, classDF))
TukeyHSD(aov(Signature ~ BaileyClass, classDF))































