################################################
#Code to develop a small signature from the genes chosen
#
#
################################################


#Call Libraries
library("survival");
library("ggplot2");
library("reshape2");
library("stringr");
library("scales")
library("GEOquery")
library("preprocessCore");
library("limma")
library("pheatmap")
library("corrplot")
source("helper/voomLimma.R")
source("helper/pubTheme.R")
source("helper/DemGGPlotSummary.R")
source("helper/")
cutoffHigh <- 24
cutoffLow <- 12


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

pr <- cleanFormatDisc(sample_pr, donor_pr, spec_pr, geneAnnot, exprs_pr);

exprs_pr <- pr[[1]];
annot_pr <- pr[[2]];


controlGenes <- c("TUBB", "ACTB", "UBC", "PPIA", "GUSB");
deltaCol <- colMeans(exprs_pr[controlGenes,])
deltaCol <- mean(deltaCol)-deltaCol;

normData <- function(x)
{
	x <- x+deltaCol;
}

exprs_pr <- data.frame(t(apply(exprs_pr, FUN=normData, MARGIN=1)));
pr[[1]] <- exprs_pr;


#Filter to only genes from initial analysis
discGenes <- read.delim("../results/SupplementalTable1.txt")
discGenes <- unique(discGenes[,"Gene"]);
######################
exprs_pr <- exprs_pr[intersect(rownames(exprs_pr), discGenes),]; #697 genes


#Let's create cutoff's
lowSamps <- intersect(rownames(annot_pr[annot_pr[,"OS_MONTHS"]<cutoffLow,]), rownames(annot_pr[annot_pr[,"OS_STATUS"]==1,]));
highSamps <- intersect(rownames(annot_pr[annot_pr[,"OS_MONTHS"]>cutoffHigh,]), rownames(annot_pr[annot_pr[,"OS_STATUS"]==0,]));
print(paste("There are ", length(lowSamps), " samples in the low group", sep="")); #70 Samples
print(paste("There are ", length(highSamps), " samples in the high group", sep="")); #46 Samples


#Now let's print clinical summary information on the lowSamps & highSamps
annot_pr_demo <- annot_pr[c(highSamps,lowSamps) ,c("tumour_grade", "donor_sex", "donor_age_at_diagnosis", "donor_tumour_stage_at_diagnosis")]
colnames(annot_pr_demo) <- c("Grade", "Gender", "Age", "Stage");
annot_pr_demo[,1] <- as.factor(annot_pr_demo[,1]);
annot_pr_demo[,2] <- as.factor(annot_pr_demo[,2]);
annot_pr_demo[,3] <- as.numeric(annot_pr_demo[,3]);
annot_pr_demo[,4] <- as.factor(annot_pr_demo[,4]);

groupsList <- list(highSamps, lowSamps)
demData <- annot_pr_demo
groupsLab <- c("Survival+", "Survival-");
demData[demData[,"Stage"]=="","Stage"]<- "NA"
demData[demData[,"Grade"]=="","Grade"]<- "NA"

#Update column names to be the right case
p <- SummaryImage(demData, groupsList, groupsLab, "Population Cohort Characteritics")
ggsave("../results/SupplementalFigure4.eps", plot=p, width=5, height=10);

#Now code to come up with signature
runLimma <- function(lowSamps, highSamps)
{
	targ <- data.frame(c(lowSamps, highSamps), c(rep("Poor", length(lowSamps)), rep("Good", length(highSamps))))
	rownames(targ) <- targ[,1];
	targ <- targ[-1];
	colnames(targ)[1] <- "targ"
	exprs_pr_tmp <- exprs_pr[,rownames(targ)]
	targ <- targ[,1]
	design <- model.matrix(~0+targ);
	colnames(design) <- gsub("targ", "", colnames(design));
	fit <- lmFit(exprs_pr_tmp, design);
	cont <- "Poor-Good"
	contrast.matrix <- eval(as.call(c(as.symbol("makeContrasts"),as.list(cont),levels=list(design))))
	fit2 <- contrasts.fit(fit, contrast.matrix)
	fit2 <- eBayes(fit2)
	outputAll <- topTable(fit2, number=40000)
	outputAll[,"ID"] <- rownames(outputAll);
	return(outputAll)
}
discGenes <- rownames(exprs_pr);

#Set a seed
set.seed(6292)

lSamps <- sample(lowSamps, 15);
hSamps <- sample(highSamps, 15);
tmpOut <- runLimma(lSamps, hSamps);
myDFPVAL <- tmpOut[discGenes,"P.Value"]
numIter <- 9;
for(i in 1:numIter)
{
	lSamps <- sample(lowSamps, 15);
	hSamps <- sample(highSamps, 15);
	tmpOut <- runLimma(lSamps, hSamps);
	tmpOut <- tmpOut[discGenes,]
	myDFPVAL <- cbind(myDFPVAL, tmpOut[,"P.Value"])

}
rownames(myDFPVAL) <- discGenes;
colnames(myDFPVAL) <- paste("Iteration_", c(1:(numIter+1)), sep="");

#Let's find out the gene count per iteration
GeneCountAllIter <- rowSums(myDFPVAL<0.05);
dfForPlot <- data.frame(GeneCountAllIter)
colnames(dfForPlot) <- "myVal";
p <- ggplot(dfForPlot, aes(myVal))+geom_histogram();
p <- p+geom_vline(xintercept = 5, color="red")
p <- p+theme_Publication()+ggtitle("Number of Genes Signficant in N Iterations")
p <- p+xlab("Number of Iterations")+ylab("Number of Genes");
ggsave("../results/SupplementalFigure5A.eps", plot=p, width=7, height=7);
GeneCountAllIter <- GeneCountAllIter[GeneCountAllIter>5]
sigGenes <- names(GeneCountAllIter);


##########################
#Okay now let's do all pair regression
##########################
cairo_ps("../results/SupplementalFigure5B.eps", width=7, height=7)
corrplot.mixed(cor(t(exprs_pr[sigGenes,c(lowSamps, highSamps)])), title="Correlation Matrix of Filtered Genes", mar=c(0,0,1,0));
dev.off();

#Remove 3 genes that are highly correlated...removed MCM4, LOX, DKK1
sigGenes <- c("ADM", "ASPM", "KRT6A", "DCBLD2", "E2F7");

#Run on all data
sigGenes <- runLimma(lowSamps, highSamps)[sigGenes,];
sigGenes <- sigGenes[,c("logFC", "P.Value")]
sigGenes[,"Direction"] <- ifelse(sigGenes[,"logFC"]>0, "Up", "Down");
geneInfo <- read.delim("~/Documents/Data/AnnotationFiles/Entrez_gene/Homo_sapiens.gene_info", header=F);
sigGenes <- data.frame(rownames(sigGenes), sigGenes)
colnames(sigGenes)[1] <- "Gene";
sigGenes <- merge(sigGenes, geneInfo[,c("V3", "V8", "V9")], by.x="Gene", by.y="V3")
colnames(sigGenes)[5:6] <- c("Chromosomal Location", "Gene Description");
print(nrow(sigGenes))
write.table(sigGenes[-2], "../results/Table1.txt", sep="\t", row.names=F);











