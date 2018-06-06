########################################
#Kaplan-meieir based on finding optimal cut-point
########################################

#Call libraries
library("survival");

#1 = alive
#2 = dead

kmScan <- function(x, perc, tVar, eVar)
{
    
    #First sort data frame
    timeVar <- x[,as.character(tVar)];
    eventVar <- x[,as.character(eVar)];
    bestI <- 1;
    pVal <- 1;
    myDF <- data.frame();
    
    minSamps <- round(perc*dim(x)[1])
    
    #Now we must loop
	for(i in minSamps:(dim(x)[1]-minSamps))
	{
        x[,"Gene"] <- 1;
        x[1:i,"Gene"] <- 0;
        x.Surv <- Surv(timeVar, eventVar);
        myP <- pchisq(survdiff(x.Surv~Gene, data=x, rho=0)$chisq, df=1, lower=F);
        
        if(myP<pVal)
        {
            pVal <- myP;
            bestI <- i;
        }
        
        #put all data in a frame to get adj p-value
        myDF <- rbind(myDF, c(i,myP));
	}
	
	#now p.adjust and return
	adjpval <- min(p.adjust(myDF[,2], method="bonferroni"));
    
	return(c(bestI,pVal, adjpval));
    
}


#This will yield a p-value and show a plot
kapmPlot <- function(metric, myData, createPlot=T, perc=.05,  tVar="time", eVar="event")
{

    #Get metadata
    tmpMeta <- myData
    tmpMeta[,metric] <- as.numeric(tmpMeta[,metric]);
    tmpMeta <- tmpMeta[order(tmpMeta[,metric]),]
    colnames(tmpMeta)[5] <- "Gene";
      
    #Run scan
    tmpMetaScan <- tmpMeta;
    tmpMetaScan <- tmpMetaScan[,as.character(c(tVar, eVar, "Gene"))];
    out <- kmScan(tmpMetaScan, perc, tVar, eVar);
    
    #Sort DF and set it
    tmpMetaScan[,"GeneBin"] <- 1;
    tmpMetaScan[1:out[[1]],"GeneBin"] <- 0;
    
    #time
    timeVar <- tmpMetaScan[,tVar];
    #event
    eventVar <- tmpMetaScan[,eVar];
    
    #createsurvival
    t.Surv <- Surv(timeVar, eventVar);
    t.survfit <- survfit(t.Surv~GeneBin, data=tmpMetaScan);
    
    #Change strata names
    myLowName <- paste("Low : n = ", t.survfit$n[[1]], sep="");
    myHighName <- paste("High : n = ", t.survfit$n[[2]], sep="");
    names(t.survfit$strata) <- c(myLowName, myHighName)
    t.survframe <- createSurvivalFrame(t.survfit)
    
    if(createPlot==T)
    {
        tmpTitle <- paste("KM Plot - Sig Score", "\nP-val(Adj) :", format(out[2], scientific=T, digits=3), "(", format(out[3], scientific=T, digits=3), ")");
        myReturn <- qplot_survival(t.survframe, f.CI=F, myTitle=tmpTitle)+theme_bw()+scale_colour_manual(values=c("red", "blue") );
    }
    
    if(createPlot==F)
    {
        
        myReturn <- c(genes, out[2], out[3]);
        
    }
    
    myReturn;
    
}

createSurvivalFrame <- function(f.survﬁt){
    # initialise frame variable
    f.frame <- NULL
    
    # check if more then one strata
    if(length(names(f.survﬁt$strata)) == 0){
        
        # create data.frame with data from survﬁt
        f.frame <- data.frame(time=f.survﬁt$time, n.risk=f.survﬁt$n.risk, n.event=f.survﬁt$n.event, n.censor = f.survﬁt
        $n.censor, surv=f.survﬁt$surv, upper=f.survﬁt$upper, lower=f.survﬁt$lower)
        
        # create ﬁrst two rows (start at 1)
        f.start <- data.frame(time=c(0, f.frame$time[1]), n.risk=c(f.survﬁt$n, f.survﬁt$n), n.event=c(0,0),
        n.censor=c(0,0), surv=c(1,1), upper=c(1,1), lower=c(1,1))
        
        # add ﬁrst row to dataset
        f.frame <- rbind(f.start, f.frame)
        
        # remove temporary data
        rm(f.start)
    }
    else {
        
        # create vector for strata identiﬁcation
        f.strata <- NULL
        for(f.i in 1:length(f.survﬁt$strata)){
            
            # add vector for one strata according to number of rows of strata
            f.strata <- c(f.strata, rep(names(f.survﬁt$strata)[f.i], f.survﬁt$strata[f.i]))
        }
        
        # create data.frame with data from survﬁt (create column for strata)
        f.frame <- data.frame(time=f.survﬁt$time, n.risk=f.survﬁt$n.risk, n.event=f.survﬁt$n.event, n.censor = f.survﬁt
        $n.censor, surv=f.survﬁt$surv, upper=f.survﬁt$upper, lower=f.survﬁt$lower, strata=factor(f.strata))
        
        # remove temporary data
        rm(f.strata)
        
        # create ﬁrst two rows (start at 1) for each strata
        for(f.i in 1:length(f.survﬁt$strata)){
            
            # take only subset for this strata from data
            f.subset <- subset(f.frame, strata==names(f.survﬁt$strata)[f.i])
            
            # create ﬁrst two rows (time: 0, time of ﬁrst event)
            f.start <- data.frame(time=c(0, f.subset$time[1]), n.risk=rep(f.survﬁt[f.i]$n, 2), n.event=c(0,0),
            n.censor=c(0,0), surv=c(1,1), upper=c(1,1), lower=c(1,1), strata=rep(names(f.survﬁt$strata)[f.i],
            2))
            
            # add ﬁrst two rows to dataset
            f.frame <- rbind(f.start, f.frame)
            
            # remove temporary data
            rm(f.start, f.subset)
        }
        
        # reorder data
        f.frame <- f.frame[order(f.frame$strata, f.frame$time), ]
        
        # rename row.names
        rownames(f.frame) <- NULL
    }
    
    # return frame
    return(f.frame)
}


qplot_survival <- function(f.frame, f.CI="default", f.shape=3, myTitle){
    # use different plotting commands dependig whether or not strata's are given
    if("strata" %in% names(f.frame) == FALSE){
        # confidence intervals are drawn if not specified otherwise
        if(f.CI=="default" | f.CI==TRUE ){
            # create plot with 4 layers (first 3 layers only events, last layer only censored)
            # hint: censoring data for multiple censoring events at timepoint are overplotted
            # (unlike in plot.survfit in survival package)
            ggplot(data=f.frame) + geom_step(aes(x=time, y=surv), direction="hv") + geom_step(aes(x=time,
            y=upper), directions="hv", linetype=2) + geom_step(aes(x=time,y=lower), direction="hv", linetype=2) +
            geom_point(data=subset(f.frame, n.censor==1), aes(x=time, y=surv), shape=f.shape)+scale_y_continuous(limits = c(0, 1))
        }
        else {
            # create plot without confidence intervalls
            ggplot(data=f.frame) + geom_step(aes(x=time, y=surv), direction="hv") +
            geom_point(data=subset(f.frame, n.censor==1), aes(x=time, y=surv), shape=f.shape)+scale_y_continuous(limits = c(0, 1))
            
        }
    }
    else {
        if(f.CI=="default" | f.CI==FALSE){
            # without CI
            ggplot(data=f.frame, aes(group=strata, colour=strata, shape=strata)) + geom_step(aes(x=time, y=surv),
            direction="hv") + geom_point(data=subset(f.frame, n.censor==1), aes(x=time, y=surv), shape=f.shape)+ggtitle(myTitle)+scale_y_continuous(limits = c(0, 1));
        }
        else {
            # with CI (hint: use alpha for CI)
            ggplot(data=f.frame, aes(colour=strata, group=strata)) + geom_step(aes(x=time, y=surv),
            direction="hv") + geom_step(aes(x=time, y=upper), directions="hv", linetype=2, alpha=0.5) +
            geom_step(aes(x=time,y=lower), direction="hv", linetype=2, alpha=0.5) +
            geom_point(data=subset(f.frame, n.censor==1), aes(x=time, y=surv), shape=f.shape)+scale_y_continuous(limits = c(0, 1))
        }
    }
}

