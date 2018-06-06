######################################
#Code to compare gene list & Pathways
#to known pancreatic cancer gene signatures, lists
#and cancer gene census
#
#Pichai Raman
#4/25/2018
######################################

#plot
library("gplots");
library("ggplot2");
library("scales");
source("helper/pubTheme.R")

#Let's read in the original list
data <- read.delim("../results/SupplementalTable1.txt");


###########################################
#Let's create table 1 of known oncogenes
##########################################
#Comparison to the cancer gene census
cgcList <- read.csv("../data/GeneLists/cancer_gene_census.csv");
cgcListSkinny <- unique(cgcList[,c("Gene.Symbol", "Molecular.Genetics")]);
intersectedGenesPancCancer <- intersect(unique(cgcList[,1]), rownames(data))

#Hypergeometic test 
runHypGeom <- function(set, genes,n=20000)
{
#number of white balls
x <- length(intersect(genes, set));

#white balls
m <- length(genes);

#black balls
n2 <- n-m; 

#balls drawn from the urn 
k <- length(set);


out <- phyper(x-1, m, n2, k, lower.tail=F);
setSize <- k;
overLap <- x;
numGenes <- m;

myRet <- c(setSize, numGenes, overLap, out); 
return(myRet);

}
#39819 comes from HGNC : http://www.genenames.org/cgi-bin/statistics
runHypGeom(rownames(data), unique(cgcList[,1]), 25000);


###########################################
#Let's compare to PDAC List
##########################################
pancUpSig <- as.character(read.delim("../data/GruetzzmanPancGeneSetUp.txt", skip=2)[,1]);
print(paste("Comparison to PDAC sig Up", runHypGeom(rownames(data[data[,"TNC_logFC"]>0,]), pancUpSig, 25000), sep=": "));
pancDownSig <- as.character(read.delim("../data/GruetzzmanPancGeneSetDown.txt", skip=2)[,1]);
print(paste("Comparison to PDAC sig Down", runHypGeom(rownames(data[data[,"TNC_logFC"]<0,]), pancDownSig, 25000), sep=": "));



geneInfo <-read.delim("~/Documents/Data/AnnotationFiles/Entrez_gene/Homo_sapiens.gene_info", header=F);
tmpTable <- geneInfo[geneInfo[,"V3"]%in%intersectedGenesPancCancer,]
tmpTable <- tmpTable[,c("V3", "V8", "V9")]
colnames(tmpTable) <- c("Gene Symbol", "Location", "Description");
tmpTable <- tmpTable[order(tmpTable[,1]),];
tmpTable[,"Molecular.Genetics"] <- cgcListSkinny[cgcListSkinny[,1]%in%tmpTable[,1],2];
tmpTable <- tmpTable[tmpTable[,"Molecular.Genetics"]=="Dom",];
tmpTable <- tmpTable[1:3];
#write.table(tmpTable, "../results/Table1.txt", sep="\t", row.names=F);



###########################################
#Let's compare to known oncogenes
##########################################

fifteenGL <- gsub(" ", "", as.character(read.delim("../data/GeneLists/fifteenGeneSig.txt")[,1]));
thirteenGL <- gsub(" ", "", as.character(read.delim("../data/GeneLists/thirteenGeneSig.txt")[,1]));
thirtySixGL <- gsub(" ", "", as.character(read.delim("../data/GeneLists/thirtySixGeneSig.txt")[,1]));
fortyEightGeneGL <- gsub(" ", "", as.character(read.delim("../data/GeneLists/FourtyEightGeneSig.txt")[,1]));
ourSig <- rownames(data);
createListFV <- list(Moffit=fifteenGL, U.Virginia=thirteenGL, Barts=thirtySixGL, Indiana.U=fortyEightGeneGL, Drexel.U=ourSig)
cairo_ps("../results/SupplementalFigure_notused2.eps", width=10, height=10);
venn(createListFV);
dev.off();

#Hypergeometric test
runHypGeom(rownames(data), fifteenGL, 25000) #Moffitt P-value 2.02 * 10^-6
runHypGeom(rownames(data), thirtySixGL, 25000) #Barts P-value 5.90 * 10^-6

#Let's also do a bar chart with % of each signature our sig captures
fifPerc <- length(intersect(ourSig, fifteenGL))/length(fifteenGL)
thirPerc <- length(intersect(ourSig, thirteenGL))/length(thirteenGL)
thirSizPerc <- length(intersect(ourSig, thirtySixGL))/length(thirtySixGL)
fourEigPerc <- length(intersect(ourSig, fortyEightGeneGL))/length(fortyEightGeneGL)

myBar <- data.frame(c("Moffit", "Barts", "Indiana.U"), c(fifPerc, thirSizPerc, fourEigPerc));
colnames(myBar) <- c("Signature", "Percentage_Overlap")
p <- ggplot(myBar, aes(Signature, Percentage_Overlap))+geom_bar(stat="identity")+theme_Publication();
p <- p+ylab("Percentage Overlap")+ggtitle("Overlap Percentage of Public Signatures")+scale_y_continuous(labels=percent)
ggsave("../results/Figure2C.eps", width=7, height=5);


###########################
#Reviewer Comments
###########################

#mofAnalysis <- SurvivalDEG[fifteenGL,]
#mofAnalysis <- na.omit(mofAnalysis); # Only 13 genes now

#Number of genes that are hits
#table[mofAnalysis[,"Hit"]==T,] # 11 of these were hits
#mofAnalysis <- mofAnalysis[mofAnalysis[,"Hit"],]
#tumNormDEGenes[rownames(mofAnalysis),]









