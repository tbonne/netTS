
#' Plotting function for nodeTS dataframes
#'
#' This function will plot the output of the nodeWin function
#' @param df.ts output dataframe from the nodeWin function
#' @param nEvents Opional argument to plot the number of events
#' @param dates Optional argument to plot the date as opposed to the time since the first event
#' @import ggplot2
#' @importFrom reshape2 melt
#' @examples
#'
#' ts.out<-nodeTS(event.data=groomEvents[1:200,])
#' nodeTS.plot(ts.out)
#'
#' @export
nodeTS.plot <- function(df.ts, nEvents = FALSE, dates = FALSE){

  if(nEvents ==FALSE){

    if(dates==FALSE){

      df.melt<-melt(df.ts, id=c("windowEnd"))
      df.melt<-df.melt[complete.cases(df.melt),]
      df.melt<-filter(df.melt, variable != "windowStart" & variable != "windowEnd" & variable != "windowStartDate" & variable != "windowEndDate", variable != "nEvents")

      fig<-ggplot(df.melt, aes(x=windowEnd, y=value, group=variable, color=variable))+ geom_line()+
        geom_point(color="blue") +
        labs(x= "Time since start")

    } else {

      df.melt<-melt(df.ts, id=c("windowEndDate"))
      df.melt<-df.melt[complete.cases(df.melt),]
      df.melt<-filter(df.melt, variable != "windowStart" | variable != "windowEnd" | variable != "windowStartDate" | variable != "windowEndDate", variable != "nEvents")

      fig<-ggplot(df.melt, aes(x=windowEndDate, y=value, group=variable, color=variable))+ geom_line()+
        geom_point(color="blue") +
        labs(x= "Time since start")
    }

  } else{

    if(dates==FALSE){

      df.melt<-melt(df.ts, id=c("windowEnd"))
      df.melt<-df.melt[complete.cases(df.melt),]
      df.melt<-filter(df.melt, variable != "windowStart" | variable != "windowEnd" | variable != "windowStartDate" | variable != "windowEndDate")

      fig<-ggplot(df.melt, aes(x=windowEnd, y=nEvents, group=variable, color=variable))+ geom_line()+
        geom_point(color="blue") +
        labs(x= "Time since start")

    } else {

      df.melt<-melt(df.ts[,(ncol(df.ts)-1)], id=c("windowEndDate"))
      df.melt<-df.melt[complete.cases(df.melt),]
      df.melt<-filter(df.melt, variable != "windowStart" | variable != "windowEnd" | variable != "windowStartDate" | variable != "windowEndDate")

      fig<-ggplot(df.melt, aes(x=windowEndDate, y=nEvents, group=variable, color=variable))+ geom_line()+
        geom_point(color="blue") +
        labs(x= "Time since start")
    }
  }

  fig

  return(fig)
}


#' Estimate consine similarity between nodes in a graph
#'
#' This function will calculate the cosine similarity between nodes at two time periods.
#' @param graph1 first graph (igraph)
#' @param graph2 second graph (igraph)
#' @export
#' @import igraph
#' @importFrom dplyr full_join
#' @examples
#'
#' #two random graphs
#' library(igraph)
#' graph1 <-erdos.renyi.game(n=10,p=0.1)
#' graph2 <-erdos.renyi.game(n=10,p=0.1)
#' cosine_between_graphs_nodes(graph1,graph2)
#'
#' #moving windo with cosine similarity
#' ts.out<-nodeTS(event.data=groomEvents[1:100,], type='cosine')
#' nodeTS.plot(ts.out)
#'
cosine_between_graphs_nodes<- function(graph1, graph2){

  node.cosine <- vector()

  #create weighted edge list from first graph
  g1.edges<-as.data.frame(get.edgelist(graph1, names=TRUE))
  if(is.null(igraph::E(graph1)$weight)){
    g1.edges$weight <- rep(1,nrow(g1.edges))
  } else {
    g1.edges$weight<-igraph::E(graph1)$weight
  }
  g1.edges$joinC <- ifelse(as.character(g1.edges$V1) < as.character(g1.edges$V2), paste(g1.edges$V1, g1.edges$V2), paste(g1.edges$V2, g1.edges$V1))

  #create weighted edge list from second graph
  g2.edges<-as.data.frame(get.edgelist(graph2, names=TRUE))
  if(is.null(igraph::E(graph2)$weight)){
    g2.edges$weight <- rep(1,nrow(g2.edges))
  } else {
    g2.edges$weight<-igraph::E(graph2)$weight
  }
  g2.edges$joinC <- ifelse(as.character(g2.edges$V1) < as.character(g2.edges$V2), paste(g2.edges$V1, g2.edges$V2), paste(g2.edges$V2, g2.edges$V1))

  #join the weighted edge lists, setting NAs equal to 0
  comb<-full_join(g1.edges,g2.edges,by="joinC")
  comb$weight.x[is.na(comb$weight.x)]<-0
  comb$weight.y[is.na(comb$weight.y)]<-0

  names.unique<-unique(comb$V1.x)

  for(i in 1:(length(names.unique)-1)){
    temp.node.w <- dplyr::filter(comb, V1.x==names.unique[i])
    node.cosine[length(node.cosine)+1]<-lsa::cosine(temp.node.w$weight.x,temp.node.w$weight.y)
  }

  #node.cosine[is.nan(node.cosine)]<-0
  cos.df<-as.data.frame(t(node.cosine))
  names(cos.df)<-names.unique[-length(names.unique)]

  return(cos.df)

}

#' netWin function
#'
#' This function will take a dataframe with events between individuals/objects, and take node level measures using a moving window approach.
#' A time column is required.
#' @param event.data dataframe containing events between individuals/objects
#' @param windowSize size of the window in which to make network measures (should be the same scale as the time column)
#' @param windowShift the amount to shift the moving window for each measure
#' @param windowStart The time of the first window. This should corespond to the first events.
#' @param type is the type of measure to take in each window. Currently available are: betweennes, closeness, eigen, and cc (i.e., clustering coeficient)
#' @param directedNet Whether the events are directed or no: true or false.
#' @param threshold minimum number of events to calculate a network measure (otherwise NA is produced).
#' @param lag The lag at which to calculate cosine similarity in network structure.
#' @param startDate Optional argument to set the date of the first event.
#' @export
#' @import igraph
#' @importFrom plyr rbind.fill
#' @import lubridate

#' @examples
#'
#' ts.out<-nodeTS(event.data=groomEvents[1:200,])
#' nodeTS.plot(ts.out)
#'
nodeTS <- function (event.data,windowSize =30,windowShift= 1, type="cc",directedNet=T, threshold=30,windowStart=0,lag=1, startDate=NULL){

  #intialize
  windowEnd=windowStart+windowSize
  netValues <- data.frame()
  gplist <- rep(list(NA),lag)

  if(windowEnd>max(event.data$time))print("Error: the window size is set larger than the max time difference")

  #set global dataframe with proper names
  g.global <- create.a.network(event.data)
  netValues<- t(rep(1,vcount(g.global)))
  colnames(netValues)<-names(igraph::V(g.global))


  #for every window
  while (windowStart + windowSize<=max(event.data$time)) {

    #subset the data
    df.window<-create.window(event.data, windowStart, windowEnd)

    #if there is enough data in this window...
    if(nrow(df.window)>threshold){

      #if there is no previous network, and the measure requires one
      if(is.na(gplist[1]) & type=='cosine'){

        gplist[[length(gplist)+1]] <- create.a.network(df.window)
        gplist<-gplist[-1]
        df.measure <- as.data.frame(NA)
        measure.uncertainty<-c(NA,NA,NA)
        measure.random<-c(NA,NA,NA)

      } else {

        #create a network
        g <- create.a.network(df.window)

        #calculate measure
        if(type=='between')measure <- betweenness(g)
        if(type=='eigne')measure <- eigen_centrality(g)
        if(type=='close')measure <- closeness(g)
        if(type=='cc')measure <- transitivity(g, type=c('local'))
        if(type=='degree')measure <- degree(g)
        if(type=='strength')measure <- strength(g)
        if(type=='cosine')measure <- cosine_between_graphs_nodes(graph1=g,graph2=gplist[[1]])

        #create a dataframe with the measures
        if(type=='cosine'){
          df.measure <- measure
          df.measure$windowStart <- windowStart
          df.measure$windowEnd <- windowEnd
        } else {
          df.measure <- data.frame(t(measure))
          colnames(df.measure)<-names(igraph::V(g))
          df.measure$windowStart <- windowStart
          df.measure$windowEnd <- windowEnd
        }

        #estimate uncertainty of measure (bootstrap)
        #measure.uncertainty <- estimate.uncertainty.boot(df.window, nb=nBoot, type=type, directedNet = directedNet)

        #esitmate range of random (permutation)
        #measure.random <- estimate.random.range.perm(g, np=nPerm, type=type, directedNet = directedNet)

        #save last network
        gplist[[length(gplist)+1]] <- g
        gplist<-gplist[-1]

      }
    } else {
      df.measure <- as.data.frame(NA)
      #measure.uncertainty<-c(NA,NA,NA)
      #measure.random<-c(NA,NA,NA)
      gp=NA
    }

    #get window date range
    windowStartDate <- NA
    windowEndDate <- NA
    if(is.null(startDate)==FALSE){
      windowStartDate <- mdy(startDate) + days(windowStart)
      windowEndDate <- mdy(startDate) + days(windowEnd)
    }
    df.measure$windowStartDate <- windowStartDate
    df.measure$windowEndDate <- windowEndDate
    df.measure$nEvents <- nrow(df.window)

    #record values
    netValues <- rbind.fill(list(as.data.frame(netValues),df.measure ))

    #move window over
    windowEnd = windowEnd + windowShift
    windowStart = windowStart + windowShift

  }

  return (netValues[-1,])
}

