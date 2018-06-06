
library("GSEABase")

#C1 Positional Sets
PositionalSetsC1 <- getGmt("~/Documents/Code/Gsea/GeneSetDatabases/Archive_051017/c1.all.v3.1.symbols.gmt", collectionType=BroadCollection(), geneIdType= SymbolIdentifier())
PositionalSetsC1 <- geneIds(PositionalSetsC1);

#C2 Canonical Pathway Sets
CanPathwaysC2 <- getGmt("~/Documents/Code/Gsea/GeneSetDatabases/Archive_051017/c2.cp.v3.1.symbols.gmt", collectionType=BroadCollection(), geneIdType= SymbolIdentifier())
CanPathwaysC2 <- geneIds(CanPathwaysC2);

#C2 Biocarta sets
BiocartaC2 <- getGmt("~/Documents/Code/Gsea/GeneSetDatabases/Archive_051017/c2.cp.biocarta.v3.1.symbols.gmt", collectionType=BroadCollection(), geneIdType= SymbolIdentifier())
BiocartaC2 <- geneIds(BiocartaC2);

#C2 KEGG Sets
KeggC2 <- getGmt("~/Documents/Code/Gsea/GeneSetDatabases/Archive_051017/c2.cp.kegg.v3.1.symbols.gmt", collectionType=BroadCollection(), geneIdType= SymbolIdentifier())
KeggC2 <- geneIds(KeggC2);

#C2 Reactome Sets
ReactomeC2 <- getGmt("~/Documents/Code/Gsea/GeneSetDatabases/Archive_051017/c2.cp.reactome.v3.1.symbols.gmt", collectionType=BroadCollection(), geneIdType= SymbolIdentifier())
ReactomeC2 <- geneIds(ReactomeC2);

#C3 MiRNA Sets
MirC3 <- getGmt("~/Documents/Code/Gsea/GeneSetDatabases/Archive_051017/c3.mir.v3.1.symbols.gmt", collectionType=BroadCollection(), geneIdType= SymbolIdentifier())
MirC3 <- geneIds(MirC3);

#C3 TF Sets
TfC3 <- getGmt("~/Documents/Code/Gsea/GeneSetDatabases/Archive_051017/c3.tft.v3.1.symbols.gmt", collectionType=BroadCollection(), geneIdType= SymbolIdentifier())
TfC3 <- geneIds(TfC3);

#C4 Cancer Gene Neighborhood Sets
CgnC4 <- getGmt("~/Documents/Code/Gsea/GeneSetDatabases/Archive_051017/c4.cgn.v3.1.symbols.gmt", collectionType=BroadCollection(), geneIdType= SymbolIdentifier())
CgnC4 <- geneIds(CgnC4);

#C4 Cancer Modules
CmC4 <- getGmt("~/Documents/Code/Gsea/GeneSetDatabases/Archive_051017/c4.cm.v3.1.symbols.gmt", collectionType=BroadCollection(), geneIdType= SymbolIdentifier())
CmC4 <- geneIds(CmC4);

#C5 Go Biological Process
GoBPC5 <- getGmt("~/Documents/Code/Gsea/GeneSetDatabases/Archive_051017/c5.bp.v3.1.symbols.gmt", collectionType=BroadCollection(), geneIdType= SymbolIdentifier())
GoBPC5 <- geneIds(GoBPC5);

#C5 Go Molecular Function
GoMFC5 <- getGmt("~/Documents/Code/Gsea/GeneSetDatabases/Archive_051017/c5.mf.v3.1.symbols.gmt", collectionType=BroadCollection(), geneIdType= SymbolIdentifier())
GoMFC5 <- geneIds(GoMFC5);

#C5 GO Cellular Component
GoCCC5 <- getGmt("~/Documents/Code/Gsea/GeneSetDatabases/Archive_051017/c5.cc.v3.1.symbols.gmt", collectionType=BroadCollection(), geneIdType= SymbolIdentifier())
GoCCC5 <- geneIds(GoCCC5);

#C6 Oncogenic Signatures
OncSigC6 <- getGmt("~/Documents/Code/Gsea/GeneSetDatabases/Archive_051017/c6.all.v3.1.symbols.gmt", collectionType=BroadCollection(), geneIdType= SymbolIdentifier())
OncSigC6 <- geneIds(OncSigC6);

print("Done Getting Sets for Functional Enrichment");


#Functions to run functional enrichment on any gene set 
#based on mSigDBsets, use 13173 as universe size as that's the
#smallest set but can be reset
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

#Accessory for functional enrichment
funcEnrichment <- function(genes, sets, qval=.05, numRet=5)
{

out <- lapply(sets, FUN = runHypGeom, genes = genes);
out <- data.frame(out);
out <- data.frame(t(out));
out$ADJ_P_VAL <- p.adjust(out[,4], method="BH");
colnames(out)[1:5] <- c("SET_SIZE", "NUM_GENES_INPUT", "OVERLAP", "P_VAL", "ADJ_P_VALUE");
out <- out[order(out[,5]),];
if(min(out[,5])>qval)
{
out <- out[1:numRet,];
}
if(min(out[,5])<qval)
{
out <- out[out[,5]<qval,];
}
return(out);

}



funcEnrichmentAll <- function(genes, qval=.05, numRet=5)
{

#C1 Positional Sets
#tmp1 <- funcEnrichment(genes, PositionalSetsC1, qval, numRet)
#tmp1[,"Collection"] <- "PositionalSetC1";

#C2 Canonical Pathway Sets
#tmp2 <- funcEnrichment(genes, CanPathwaysC2, qval)
#tmp2[,"Collection"] <- "CanPathwaysC2";

#C2 Biocarta sets
tmp3 <- funcEnrichment(genes, BiocartaC2, qval)
tmp3[,"Collection"] <- "BiocartaC2";

#C2 KEGG Sets
tmp4 <- funcEnrichment(genes, KeggC2, qval)
tmp4[,"Collection"] <- "KeggC2";

#C2 Reactome Sets
tmp5 <- funcEnrichment(genes, ReactomeC2, qval)
tmp5[,"Collection"] <- "ReactomeC2";

#C3 MiRNA Sets
#tmp6 <- funcEnrichment(genes, MirC3, qval)
#tmp6[,"Collection"] <- "MirC3";

#C3 TF Sets
tmp7 <- funcEnrichment(genes, TfC3, qval)
tmp7[,"Collection"] <- "TranscriptionFactor_C3";

#C4 Cancer Gene Neighborhood Sets
tmp8 <- funcEnrichment(genes, CgnC4, qval)
tmp8[,"Collection"] <- "CancerGeneNetC4";

#C4 Cancer Modules
tmp9 <- funcEnrichment(genes, CmC4, qval)
tmp9[,"Collection"] <- "CancerModulesC4";

#C5 Go Biological Process
tmp10 <- funcEnrichment(genes, GoBPC5, qval)
tmp10[,"Collection"] <- "GOBiologicalProcC5";

#C5 Go Molecular Function
tmp11 <- funcEnrichment(genes, GoMFC5, qval)
tmp11[,"Collection"] <- "GOMolecularFunctionC5";

#C5 GO Cellular Component
tmp12 <- funcEnrichment(genes, GoCCC5, qval)
tmp12[,"Collection"] <- "GOCellularCompC5";

#C6 Oncogenic Signatures
tmp13 <- funcEnrichment(genes, OncSigC6, qval)
tmp13[,"Collection"] <- "OncogenicSigsC6";

#finaloutput <- rbind(tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13);

finaloutput <- rbind(tmp3, tmp4, tmp5, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13);


return(finaloutput);


}



