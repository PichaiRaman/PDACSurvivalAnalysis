###########################################
#Code to run Voom and Limma on a htseq count matrix
#
#Function requires a matrix of counts and 2 or more groups
#returns back a list
#
#First object is gene, p-value, Log-FC
#Second object is the actual limma output
#
###########################################

#####################
#EXAMPLES
#
#2 Groups:
#load("../data/mouse.brain.htseq.rda");
#groupsName <- list(colnames(mouse.brain.htseq)[1:5], colnames(mouse.brain.htseq)[6:10]);
#groupsNum <- list(c(1:5), c(6:10));
#limmaOutput2<- voomLimma(mouse.brain.htseq, groupsNum);
#
#Multiple Groups:
#groupsName <- list(colnames(mouse.brain.htseq)[1:3], colnames(mouse.brain.htseq)[4:6], colnames(mouse.brain.htseq)[7:10]);
#limmaOutputMult<- voomLimma(mouse.brain.htseq, groupsName);
#
#####################

#Function to run voom and then limma, by default first group is control
#If more than one group is given the ANOVA p-value is outputter
#all other groups are always compared back to control or you get an explosion in number of comparisons
voomLimma <- function(mtrx=NULL, grps=NULL)
{
#Call libraries
require("edgeR");
require("limma");
dataMat <- mtrx;
groups <- grps;

#Handle if its a numeric vector
if(class(groups[[1]])=="integer")
{
tmpGroups <- colnames(dataMat)[groups[[1]]];
for(i in 2:length(groups))
{
tmpGroups <- list(tmpGroups, colnames(dataMat)[groups[[i]]]);
}
groups <- tmpGroups;
}

#choose only given columns
dataMat <- dataMat[,unlist(groups)];

#Generate targets file
targets <- data.frame(groups[[1]],paste("Group_", rep(1, length(groups[[1]])), sep=""))
colnames(targets) <- c("Sample", "Condition");
for(i in 2:length(groups))
{
tmpTargets <- data.frame(groups[[i]],paste("Group_", rep(i, length(groups[[i]])), sep=""))
colnames(tmpTargets) <- c("Sample", "Condition");
targets <- rbind(targets, tmpTargets);
}
rownames(targets) <- targets[,1];
targets <- targets[-1];

#Run Voom
y <- DGEList(counts=as.matrix(dataMat), genes=rownames(dataMat))
y <- calcNormFactors(y)
design <- model.matrix(~targets[,"Condition"])
v <- voom(y,design,plot=TRUE);
voomData <- v$E

#Run linear fit and Create the contrast
design <- model.matrix(~0+targets[,"Condition"])
colnames(design) <- levels(targets[,1]);
fit <- lmFit(voomData, design);
myConts <- paste(levels(targets[,1])[-1], "Group_1", sep="-");
contrast.matrix <- makeContrasts(contrasts=myConts, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

#Okay now combine all into one result data frame
myComps <- colnames(fit2$contrasts);

if(length(myComps)==1)
{
output <- topTable(fit2, number=60000)
genOutput <- output[,c("logFC", "P.Value")];
AllResult <- list(genOutput, output, voomData, design, myComps);
}

if(length(myComps)>1)
{
output <- topTable(fit2, number=60000)
genOutput <- output[,c(c(1:length(myComps)), grep("P.Value", colnames(output)))];
colnames(genOutput)[1:length(myComps)] <- myComps;
tmpOut <- topTable(fit2, number=60000, 1)[,c("logFC", "P.Value", "adj.P.Val")];
colnames(tmpOut) <- paste(myComps[1], colnames(tmpOut), sep="_");
IndComp <- tmpOut;

    for(i in 2:length(myComps))
    {
    tmpOut <- topTable(fit2, number=60000, i)[,c("logFC", "P.Value", "adj.P.Val")];
    tmpOut <- tmpOut[rownames(IndComp),];
    colnames(tmpOut) <- paste(myComps[i], colnames(tmpOut), sep="_");
    IndComp <- cbind(IndComp, tmpOut);
    }

AllResult <- list(genOutput, IndComp, voomData, design, myComps);
}
names(AllResult) <- c("effectPval", "fullResult", "voomMatrix", "design", "comparison")
return(AllResult);
}















