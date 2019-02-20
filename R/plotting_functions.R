
#' Plot function for graphTS outputs
#'
#' This function will plot the output from the graphTS function
#' @param data Output from the graphTS function.
#' @param plotCI Whether to plot the confidence intervals (Default: TRUE)
#' @importFrom ggplot2 ggplot
#' @export
#'
#'
graphTS_plot <- function(data, plotCI=TRUE){

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
nodeTS_plot <- function(data, legend=TRUE){

  data.long <- reshape2::melt(data, id = c("windowend","windowstart", "nEvents"))
  names(data.long)[names(data.long)=="variable"] <- "ID"

  if(legend==TRUE){
    g <- ggplot(data=data.long, aes(x = as.Date(windowstart), y = value,group = ID, color = ID)) + geom_line()  + theme_classic() + xlab("windowstart")
  } else {
    g <- ggplot(data=data.long, aes(x = as.Date(windowstart), y = value,group = ID, color = ID)) + geom_line()  + theme_classic() + xlab("windowstart") + theme(legend.position="none")
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
#'
dyadTS_plot <- function(data=output.net.dir, legend=FALSE){

  data.long <- reshape2::melt(data, id = c("windowend","windowstart", "nEvents"))
  names(data.long)[names(data.long)=="variable"] <- "ID"

  if(legend==TRUE){
    g <- ggplot(data=data.long, aes(x = as.Date(windowstart), y = value,group = ID, color = ID)) + geom_line()  + xlab("windowstart") + theme_classic()
  } else {
    g <- ggplot(data=data.long, aes(x = as.Date(windowstart), y = value,group = ID, color = ID)) + geom_line()  + xlab("windowstart") + theme_classic() + theme(legend.position="none")
  }

  return(g)


}
