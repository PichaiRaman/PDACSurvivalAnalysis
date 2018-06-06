############################################
#Code to run enrichment analysis on significant genes
#Pichai Raman
#4/25/2018
#############################################

#Call libraries
library("ggplot2")
library("org.Hs.eg.db")
library("clusterProfiler");

source("helper/EnrichmentHelper.R");
source("helper/pubTheme.R")


########################################################
#Run GSEA
########################################################
#Load all Data
data <- read.delim("../results/SupplementalTable1.txt")
data <- data[data[,"TNC_logFC"]>0,]
data <- data[data[,"SRV_logFC"]>0,]

poorSurvivalGenes <- rownames(data);

#Run Enrichment
poorSurvivalSets <- funcEnrichmentAll(poorSurvivalGenes, qval=.99)
poorSurvivalSets <- poorSurvivalSets[poorSurvivalSets[,"Collection"]=="ReactomeC2",]
poorSurvivalSets[,"ADJ_P_VALUE"] <- p.adjust(poorSurvivalSets[,"P_VAL"], method="BH")
cutoff <- 0.05

#filter sets
poorSurvivalSets <- poorSurvivalSets[poorSurvivalSets[,"ADJ_P_VALUE"]<cutoff,];
poorSurvivalSets[,"LogP"] <- (-1)*log10(poorSurvivalSets[,"P_VAL"])
poorSurvivalSets[,"SetName"]<- rownames(poorSurvivalSets)
poorSurvivalSets[,"SetName"]<- factor(poorSurvivalSets[,"SetName"], levels=poorSurvivalSets[,"SetName"]);
write.table(poorSurvivalSets, "../results/SupplementalTable2.txt", sep="\t", row.names=T)


#capitilize
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
      sep="", collapse=" ")
}

#Read in pathway results
data<- poorSurvivalSets[1:25,];
data[,"PathwayName"] <- gsub("_", " ", rownames(data));
data[,"PathwayName"] <- tolower(data[,"PathwayName"]);
data[,"PathwayName"] <- unlist(lapply(data[,"PathwayName"], FUN=simpleCap));
data[,"PathwayName"] <- unlist(lapply(data[,"PathwayName"], FUN=simpleCap));
data[14,"PathwayName"] <- "Reactome Transport Of Glucose And Other Sugars Bile Salts etc.."
data[,"PathwayName"] <- factor(data[,"PathwayName"], levels=data[,"PathwayName"] )

p <- ggplot(data, aes(PathwayName, LogP))+geom_bar(stat="identity")+theme_Publication()+theme(axis.text.x=element_text(angle=80,hjust=1));
p <- p+xlab("Reactome Pathway Name")+ylab("Pathway Score (-log10 P-value)")+ggtitle("Enriched Pathways from 602 Up-regulated Genes");
p
ggsave("../results/Figure2A.eps", width=7, height=10);












