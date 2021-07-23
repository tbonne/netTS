
#' Plot function for graphTS outputs
#'
#' This function will plot the output from the graphTS function
#' @param data Output from the graphTS function.
#' @param plotCI Whether to plot the confidence intervals (Default: TRUE)
#' @import ggplot2
#' @export
#'
#'
graphTS.plot <- function(data, plotCI=FALSE){

  if(plotCI==TRUE){
    g<-ggplot(data, aes(x=windowstart, y=measure))+geom_line(color="blue")+geom_point(color="blue") + geom_ribbon(aes(ymin=CI.low,ymax=CI.high),  fill="red", alpha=0.1) + theme_classic()
  } else {
    g<-ggplot(data, aes(x=windowstart, y=measure))+geom_line(color="blue")+geom_point(color="blue")  + theme_classic()
  }

  return(g)


}



#' Plot function for nodeTS outputs
#'
#' This function will plot the output from the nodeTS function
#' @param data Output from the nodeTS function.
#' @param legend Whether to include the legend (Default: TRUE)
#' @importFrom ggplot2 ggplot
#' @importFrom reshape2 melt
#' @export
#'
#'
nodeTS.plot <- function(data, legend=TRUE, plotCI=FALSE){

  #plot just the observed?
  if(plotCI==FALSE){

    if(!is.data.frame(data))data=data$obs

    data.long <- reshape2::melt(data, id = c("windowend","windowstart", "nEvents"))
    names(data.long)[names(data.long)=="variable"] <- "ID"

    if(legend==TRUE){
      g <- ggplot(data=data.long, aes(x = as.Date(windowstart), y = value,group = ID, color = ID)) + geom_line()  + theme_classic() + xlab("windowstart")
    } else {
      g <- ggplot(data=data.long, aes(x = as.Date(windowstart), y = value,group = ID, color = ID)) + geom_line()  + theme_classic() + xlab("windowstart") + theme(legend.position="none")
    }
  } else {

    #get the observed
    data.long <- reshape2::melt(data$obs, id = c("windowend","windowstart", "nEvents"))
    names(data.long)[names(data.long)=="variable"] <- "ID"

    #get the lower CI
    data.lowCI <- reshape2::melt(data$lowCI, id = c("windowend","windowstart", "nEvents"))
    data.long$lowCI <- data.lowCI$value

    #get the upper CI
    data.highCI <- reshape2::melt(data$highCI, id = c("windowend","windowstart", "nEvents"))
    data.long$highCI <- data.highCI$value

    #create the ggplots
    if(legend==TRUE){
      g <- ggplot(data=data.long, aes(x = as.Date(windowstart), y = value,group = ID, ymin=lowCI,ymax=highCI)) + geom_line(aes(color = ID)) + geom_ribbon(aes(fill=ID),color=NA, alpha=0.1)  + theme_classic() + xlab("windowstart")
    } else {
      g <- ggplot(data=data.long, aes(x = as.Date(windowstart), y = value,group = ID, ymin=lowCI,ymax=highCI)) + geom_line(aes(color = ID)) + geom_ribbon(aes(fill=ID),color=NA, alpha=0.1)  + theme_classic() + xlab("windowstart") + theme(legend.position="none")
    }


  }



  return(g)


}


#' Plot function for dyadTS outputs
#'
#' This function will plot the output from the dyadTS function
#' @param data Output from the dyadTS function.
#' @param legend Whether to include the legend (Default: TRUE)
#' @importFrom ggplot2 ggplot
#' @importFrom reshape2 melt
#' @export
#'
dyadTS.plot <- function(data=output.net.dir, legend=FALSE, plotCI=FALSE){

  #plot just the observed?
  if(plotCI==FALSE){

    if(!is.data.frame(data))data=data$obs

    data.long <- reshape2::melt(data, id = c("windowend","windowstart", "nEvents"))
    names(data.long)[names(data.long)=="variable"] <- "ID"

    if(legend==TRUE){
      g <- ggplot(data=data.long, aes(x = as.Date(windowstart), y = value,group = ID, color = ID)) + geom_line()  + theme_classic() + xlab("windowstart")
    } else {
      g <- ggplot(data=data.long, aes(x = as.Date(windowstart), y = value,group = ID, color = ID)) + geom_line()  + theme_classic() + xlab("windowstart") + theme(legend.position="none")
    }
  } else {

    #get the observed
    data.long <- reshape2::melt(data$obs, id = c("windowend","windowstart", "nEvents"))
    data.long$joinID <- paste0(data.long$variable,"_",data.long$windowstart)


    #get the lower CI
    data.lowCI <- reshape2::melt(data$lowCI, id = c("windowend","windowstart", "nEvents"))
    data.lowCI$joinID <- paste0(data.lowCI$variable,"_",data.lowCI$windowstart)
    names(data.lowCI)[names(data.lowCI)=="value"] <- "lowCI"
    data.long <- dplyr::left_join(data.long,data.lowCI[c("lowCI","joinID")], by="joinID")

    #get the upper CI
    data.highCI <- reshape2::melt(data$highCI, id = c("windowend","windowstart", "nEvents"))
    data.highCI$joinID <- paste0(data.highCI$variable,"_",data.highCI$windowstart)
    names(data.highCI)[names(data.highCI)=="value"] <- "highCI"
    data.long <- dplyr::left_join(data.long,data.highCI[c("highCI","joinID")], by="joinID")

    #change the name of the variable to ID
    names(data.long)[names(data.long)=="variable"] <- "ID"

    #create the ggplots
    if(legend==TRUE){
      g <- ggplot(data=data.long, aes(x = as.Date(windowstart), y = value,group = ID, ymin=lowCI,ymax=highCI)) + geom_line(aes(color = ID)) + geom_ribbon(aes(fill=ID),color=NA, alpha=0.1)  + theme_classic() + xlab("windowstart")
    } else {
      g <- ggplot(data=data.long, aes(x = as.Date(windowstart), y = value,group = ID, ymin=lowCI,ymax=highCI)) + geom_line(aes(color = ID)) + geom_ribbon(aes(fill=ID),color=NA, alpha=0.1)  + theme_classic() + xlab("windowstart") + theme(legend.position="none")
    }


  }

  return(g)


}


#' Plot function for windowsize check outputs
#'
#' This function will plot the output from the window.check function
#' @param data Output from the window.check function.
#' @param legend Whether to include the legend (Default: TRUE)
#' @importFrom ggplot2 ggplot
#' @importFrom reshape2 melt
#' @export
#'
#'
check.windowsize.plot <- function(data, legend=TRUE){

  if(legend==TRUE){
    g <- ggplot(data=data, aes(x = as.Date(windowstart), y = mean,group = factor(fracData), color = factor(fracData))) + geom_line() + geom_ribbon(aes(ymin=CI.low,ymax=CI.high,fill=factor(fracData)),color=NA,  alpha=0.05)  + theme_classic() + xlab("windowstart") + ylab("Similarity to observed network")
  } else {
    g <- ggplot(data=data, aes(x = as.Date(windowstart), y = mean,group = factor(fracData), color = factor(fracData))) + geom_line() + geom_ribbon(aes(ymin=CI.low,ymax=CI.high,fill=factor(fracData)),color=NA,  alpha=0.05) + theme_classic() + xlab("windowstart") + ylab("Similarity to observed network") + theme(legend.position="none")
  }

  return(g)

}
