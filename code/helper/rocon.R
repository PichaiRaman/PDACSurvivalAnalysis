##################################################
#Author : Pichai Raman
#Date : 10/15/15
#Package : This package is for the creation
#and display of ROC plots that are customizable
#with GGplot. It is meant both for novice and expert
#GGPlot & R users and has a command line component.
##################################################

#Call libraries
require("AUC");
require("ggplot2");



############################################
#Main function to generate ROC Plots
#
#data : data frame, 1st column is metric, second column is label (0=False, 1=True)
#diagCol : Color of the diagnomal line
#roccol : Color of the ROC line
#isDecreasing : by default the lower the metric the better, i.e. if you have a score of 1 its better than a score of 2 (i.e. p-values), if using fold change set to T
#myTitle : title of the plot
############################################
rocon <- function(data, diagCol="black", roccol="red", isDecreasing=T, myTitle="ROC Curve")
{
    #Removes all NA
    data <- na.omit(data);
    data[,2] <- as.numeric(as.character(data[,2]));
    pn <- nrow(subset(data, data[, 2] == 1))
    fn <- nrow(data) - pn
    diag = data.frame(x = seq(0, 1, by = .01), y = seq(0, 1,
    by = .01))
    data <- data[order(data[,1], decreasing=isDecreasing),];
    
    x = 0
    y = 0
    for (i in 1:nrow(data))
    {
        tpr <- sum(data[1:i,2])/pn;
        fpr <- (length(c(1:i))-sum(data[1:i,2]))/fn
        x <- c(x, fpr)
        y <- c(y, tpr)
    }
    
    if(isDecreasing==F)
    {
        data[,1] <- 1/data[,1];
    }
    myAuc <- auc(roc(data[,1], factor(data[,2])))
    
    #Create and return object
    rocdata <- data.frame(FPR = x, TPR = y, method="")
    legLabs <- paste("AUC = ", round(myAuc, 3), sep="");
    p <- ggplot(data = rocdata, aes(x = FPR, y = TPR, color=method)) +geom_line(size=.5) + geom_line(data = diag, aes(x = x,y = y), color ="black")+theme_bw()+ labs(x = "False Positive Rate",y = "True Positive Rate", title = "ROC curve")+scale_color_discrete(labels=legLabs)+ theme(legend.title=element_blank())

    return(list(rocdata[1:2], p, myAuc));
    
}


############################################
#Function to generate multiple ROC Plots
#
#data : data frame, 1st column is metric, second column is label (0=False, 1=True), 3rd is the method
#diagCol : Color of the diagnomal line
#isDecreasing : by default the lower the metric the better, i.e. if you have a score of 1 its better than a score of 2 (i.e. p-values), if using fold change set to T
#myTitle : title of the plot
############################################

roconMult <- function(data, diagCol="black", colors=NULL, isDecreasing=T, myTitle="ROC Curve")
{
    
    dataList <- split(data, f=data[,3]);
    diag = data.frame(x = seq(0, 1, by = .01), y = seq(0, 1, by = .01))
    output <- lapply(dataList, FUN=rocon, isDecreasing=isDecreasing)
    dfAll <- data.frame();
    allAUC <- data.frame();
    legLabs <- c();
    for(i in 1:length(output))
    {
        dfAll <- rbind(dfAll, data.frame(output[[i]][1], method=as.character(names(output)[i])));
        allAUC <- rbind(allAUC, data.frame(AUCValue=as.numeric(as.character(output[[i]][3])), method=names(output)[i]));
        legLabs <- c(legLabs, paste(method=names(output)[i], " : AUC = ", round(as.numeric(as.character(output[[i]][3])), 3), sep=""));
     }
    
    #Create and return object
    
    p <- ggplot(data = dfAll, aes(x = FPR, y = TPR, color=method)) +geom_line(size=.5) + geom_line(data = diag, aes(x = x,y = y), color ="black")+theme_bw()+ labs(x = "False Positive Rate",y = "True Positive Rate", title = myTitle)+scale_color_discrete(labels=legLabs)+ theme(legend.title=element_blank())
    return(list(dfAll, p, allAUC));
    
}




