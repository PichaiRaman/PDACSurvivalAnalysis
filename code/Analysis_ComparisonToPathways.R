##########################################
#Optimize signature, try to find weights
#of signature which give best AUC in 2 data sets
#Pichai Raman
#4/25/2018
##########################################

#Call libraries
library("GEOquery");
library("preprocessCore");
source("helper/rocon.R")
source("helper/KaplanScan.R")
source("helper/pubTheme.R")
library("tidyverse")


#Get GEO Data and save for


 
#################################
#Format ICGC Data
#################################
cleanFormatDisc <- function(sampleDat, donorDat, specDat, geneAnnot, exprs)
{

#Get Sample Level and Specimen Level Data
rownames(sampleDat) <- sampleDat[,1];
sampleDat <- sampleDat[colnames(exprs),] # 269 rows
sampleDat <- merge(sampleDat, specDat[,c("icgc_specimen_id","tumour_histological_type", "tumour_grade", "specimen_type", "specimen_donor_treatment_type")], by.x="icgc_specimen_id", by.y="icgc_specimen_id")


#Merge Sample and Donor data frames
annot <- merge(sampleDat, donorDat, by.x="icgc_donor_id", by.y="icgc_donor_id", all.x=T)
annot[,"OS_STATUS"] <- ifelse(annot[,"donor_vital_status"]=="alive", 0, 1);
annot[,"OS_MONTHS"] <- round(as.numeric(annot[,"donor_survival_time"])/30);
rownames(annot) <- annot[,3];
annot <- annot[colnames(exprs),] # 269 rows

#Filter to remove non PDAC
annot <- annot[annot[,"tumour_histological_type"]=="Pancreatic Ductal Adenocarcinoma",]


#Get max expression and sort and remove duplicates
exprs[,"maxExp"] <- apply(exprs[2:ncol(exprs)], FUN=max, MARGIN=1);
exprs <- exprs[order(-exprs[,"maxExp"]),];
rownames(geneAnnot) <- geneAnnot[,1]
geneAnnot <- geneAnnot[rownames(exprs),]

exprs <- cbind(exprs, geneAnnot[,"Symbol"])
colnames(exprs)[271] <- "Symbol";
exprs <- exprs[!duplicated(exprs[,"Symbol"]),]
rownames(exprs) <- exprs[,"Symbol"]; 
exprs <- exprs[1:269];
myProbes <- rownames(exprs);
myCols <- colnames(exprs);
exprs <- normalize.quantiles(as.matrix(exprs));
rownames(exprs) <- myProbes;
colnames(exprs) <- myCols;
exprs <- data.frame(exprs);

#Now get intersection of id's and conserve order
intSamps <- intersect(colnames(exprs), rownames(annot));
output <- list(exprs[,intSamps], annot[intSamps,]);
return(output);

}

#load expression and annotation data
load("../data/ValidationDataSets/ICGC_ArrayData.RData");
exprs_pr <- arrayData;
rm(arrayData);
sample_pr <- read.delim("../data/sample.PACA-AU.tsv");
donor_pr <- read.delim("../data/donor.PACA-AU.tsv");
spec_pr <- read.delim("../data/specimen.PACA-AU.tsv");
geneAnnot <- read.delim("../data/GPL10558-50081.txt");

#Matrix
icgcDat <- cleanFormatDisc(sample_pr, donor_pr, spec_pr, geneAnnot, exprs_pr);
exprs_icgc <- icgcDat[[1]];
annot_icgc <- icgcDat[[2]];

#Normalize by Control genes
controlGenes <- c("TUBB", "ACTB", "UBC", "PPIA", "GUSB");
deltaCol <- colMeans(exprs_icgc[controlGenes,])
deltaCol <- mean(deltaCol)-deltaCol;

normData <- function(x)
{
	x <- x+deltaCol;
}
exprs_icgc <- data.frame(t(apply(exprs_icgc, FUN=normData, MARGIN=1)));
icgcDat[[1]] <- exprs_icgc;


#Let's create cutoff's
lowSamps <- intersect(rownames(annot_icgc[annot_icgc[,"OS_MONTHS"]<cutoffLow,]), rownames(annot_icgc[annot_icgc[,"OS_STATUS"]==1,]));
highSamps <- intersect(rownames(annot_icgc[annot_icgc[,"OS_MONTHS"]>cutoffHigh,]), rownames(annot_icgc[annot_icgc[,"OS_STATUS"]==0,]));

#################################
#End Format ICGC Data
#################################

################################
#Optimization algorithm
################################

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
	intGenes <- intersect(rownames(exprs_icgc), mySig);
	exprs_pr_cmat <- exprs_icgc[intGenes,];
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
	return(list(aucDF, aucDF["score"], myScores));
}


####################
#Logistic Model First
####################

myGeneSigRes <- createSigScore(sigGenes, sigName="5-Gene Signature")
virg13 <- createSigScore(virg13Up, sigName="Virginia-13", mySigDown=virg13Down)
moffitt15 <- createSigScore(moffitt15, sigName="Moffitt-15")
barts36 <- createSigScore(bartsGeneUp, sigName="Barts-36", mySigDown=bartsGeneDown)
rocDF <- rbind(myGeneSigRes[[1]], virg13[[1]], moffitt15[[1]], barts36[[1]])
roconMult(rocDF);
ggsave("../results/Figure4ALeft.eps", width=7, height=7);






#Now let's create distribution of random signatures and put in our signature
myDist <- c();
for(i in 1:5000)
{
randomSig <- createSigScore(sample(rownames(exprs_icgc), 5), sigName="randSig")
myDist <- c(rocon(randomSig[[1]])[[3]], myDist);
}

pvalvsRandom <- 1-sum(rocon(myGeneSigRes[[1]])[[3]]>myDist)/length(myDist)
p <- qplot(myDist, bins=200)+geom_vline(xintercept=rocon(myGeneSigRes[[1]])[[3]], color="red")+theme_bw()+xlab("AUC")+ggtitle("AUC Distribution of Random Signatures (ICGC)")
p <- p+annotate("text", x = .7, y = 150, label = paste("P-value =", round(pvalvsRandom, 6)))
ggsave("../results/SupplementalFigure7Left.eps", width=7, height=7);


#Now let's do kaplan-meier plots
annot_icgc <- cbind(annot_icgc[c(lowSamps, highSamps),c("icgc_sample_id", "donor_age_at_enrollment", "OS_MONTHS", "OS_STATUS")], myGeneSigRes[[2]]);
colnames(annot_icgc)[5] <- "myScores";
kapmPlot("myScores", annot_icgc, T, perc=0.25, tVar="OS_MONTHS", eVar="OS_STATUS")+ylab("Survival")+xlab("Time");
ggsave("../results/Figure4BLeft.eps", width=7, height=7);



#############################################
#Updating bailey section based on reviewer comment
#On why so few samples : Because we didn't classify all samples previously.
#############################################

#Examine alignment of signature with bailey NMF subgroups
baileyData <- read.delim("../data/BaileyICGC_NMFClass.txt");
tmpSamp <- icgcDat[[2]][,c("icgc_sample_id", "submitted_donor_id.x")]
tmpSamp <- cbind(tmpSamp, myGeneSigRes[[3]])
tmpSamp <- na.omit(tmpSamp);
tmpSamp <- merge(tmpSamp, baileyData, by.x="submitted_donor_id.x", by.y="icgc_id");
colnames(tmpSamp)[3] <- "score"
p <- ggplot(tmpSamp, aes(membership.ordered, score, fill=membership.ordered))+geom_boxplot();
p <- p+xlab("Subtype")+ylab("Signature Score")+theme_Publication()+guides(fill=FALSE);
ggsave("../results/Figure3BLeft.eps", width=7, height=5);

#Let's also get p-values of difference
sqScore <- tmpSamp[tmpSamp[,"membership.ordered"]=="Squamous", "score"]
adScore <- tmpSamp[tmpSamp[,"membership.ordered"]=="ADEX", "score"]
imScore <- tmpSamp[tmpSamp[,"membership.ordered"]=="Immunogenic", "score"]
ppScore <- tmpSamp[tmpSamp[,"membership.ordered"]=="Pancreatic Progenitor", "score"]

sqVsAd <- t.test(sqScore, adScore) # 0.008447 (used to be 0.1712)
sqVsim <- t.test(sqScore, imScore) # 8.334e-08 (used to be 0.0003492)
sqVspp <- t.test(sqScore, ppScore) # 2.351e-06 (used to be 0.0002186)



################################################
#Respond to reviewer and do multivariate testing
#so we can see if signature has predictive power
#even with age, 
######################
multTestDF <- cbind(icgcDat[[2]][,c("tumour_histological_type", "tumour_grade", "donor_sex", "donor_age_at_enrollment", "OS_STATUS", "OS_MONTHS")], myGeneSigRes[[3]]);
colnames(multTestDF)[7] <- "SignatureScore"

#Clean data a bit
multTestDF <- multTestDF[multTestDF[,"tumour_grade"]!="",]
multTestDF <- multTestDF[multTestDF[,"tumour_grade"]!="X - Cannot be assessed",]
multTestDF[,"tumour_grade"] <- factor(multTestDF[,"tumour_grade"], levels=c("1 - Well differentiated","2 - Moderately differentiated","3 - Poorly differentiated","4 - Undifferentiated"))
multTestDF[,"donor_sex"]<- factor(multTestDF[,"donor_sex"], levels=c("male", "female"))

#Run test
coxSigAnalysis <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ SignatureScore+tumour_grade+donor_sex+donor_age_at_enrollment, data=multTestDF)
#cox.zph(coxSigAnalysis)



############################################################
#Reviewer asked to relate signature score to squamous biology
############################################################

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
ggsave("../results/Figure3BRight.eps", width=14, height=7);

MYCSigLM <- summary(lm(MYC ~ SignatureScore, sigScoreCor))
HypoxiaSigLM <- summary(lm(Hypoxia ~ SignatureScore, sigScoreCor))
EGFSigLM <- summary(lm(EGF ~ SignatureScore, sigScoreCor))
TGFBSigLM <- summary(lm(TGFB ~ SignatureScore, sigScoreCor))
AutoSigLM <- summary(lm(Autophagy ~ SignatureScore, sigScoreCor))
InflamSigLM <- summary(lm(Inflammation ~ SignatureScore, sigScoreCor))



