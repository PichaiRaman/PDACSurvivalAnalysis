###################################################
# DemGGPlotSummary Package
#
# Package will create an image out of demographic information in table
# A column of choice can be used to group the data acccordingly
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'
###################################################

#Call libraries
require("ggplot2")
require("Rmisc");
require("ggplot2")
require("gridExtra")
require("cowplot")

#tempRemove
#Function to create summary image
SummaryImage <- function(demData=NULL, groupsList=NULL, groupsLab=NULL, myTitle="")
{
  #Add Grouping column to dataframe
  demDataTmp <- demData;
  demDataTmp[,"Group"] <- "NOVALUE";
  for(i in 1:length(groupsLab))
  {
    demDataTmp[groupsList[[i]],"Group"] <- groupsLab[i]
  }
  demData <- demDataTmp[demDataTmp[,"Group"]!="NOVALUE",];

  #put plots into a list
  myList <- list();
  myList <- lapply(colnames(demData[1:(ncol(demData)-1)]), FUN=createPlot, demDataTmp=demData);
  output <- grid_arrange_plots(myPlots=myList, sharedLegend=T, myTitle=myTitle)
  #Get plots and put into vector
}

grid_arrange_plots <- function(myPlots=NULL, sharedLegend=F, myTitle="", myNrow=2) {
    if(sharedLegend==F)
    {
      out <- plot_grid(plotlist=myPlots, nrow=myNrow, labels = "AUTO", align = 'h')
    } 
    if(sharedLegend==T)
    {
      g <- ggplotGrob(myPlots[[1]] + theme(legend.position="right"))$grobs;
      legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]];
      myPlots <- lapply(myPlots, FUN=removeLegend)
      mypgGrid <- plot_grid(plotlist=myPlots, nrow=myNrow, labels = "AUTO", align = 'h')+ggtitle(myTitle)
      out <- plot_grid(mypgGrid, legend, ncol=1, rel_heights = c(6,.3))
    }
    return(out);
}

#Accessory to just remove legends
removeLegend <- function(x)
{
  x <- x+theme(legend.position="none");  
  return(x);
}

#Function to create different plots based on data types
createPlot <- function(x=NULL, demDataTmp=demData, myTheme=theme_Publication(), myColors=c("red", "blue"))
{
  type=ifelse(class(demDataTmp[,x])=="numeric", "box", "bar");

  if(type=="box")
  {
    p <- ggplot(demDataTmp, aes_string("Group", x, fill="Group"))+geom_boxplot()+myTheme+ scale_fill_brewer(palette = "Set1");
  }
  if(type=="bar")
  {
    p <- ggplot(demDataTmp, aes_string(x, fill="Group"))+geom_bar(position = "dodge")+myTheme+theme(axis.text.x=element_text(angle=75,hjust=1))+ scale_fill_brewer(palette = "Set1")+ylab("Count")
  }
  return(p);
}








