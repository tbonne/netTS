#' Create a network from an edge list
#'
#' This function will generate a network from an event dataframe
#' @param events dataframe containing all events.
#' @importFrom igraph graph_from_data_frame
#'
create.a.network<-function(events){

  elist<-create.an.edgeList(events)
  names(elist)<- c("from", "to", "weight")
  gg <- graph_from_data_frame(elist, directed = TRUE, vertices = NULL)

  return(gg)
}

#' Create a network from an edge list using the SRI.
#'
#' This function will generate an SRI network from an event dataframe
#' @param events dataframe containing all events.
#' @importFrom igraph graph_from_data_frame
#'
create.a.network.SRI<-function(events){

  elist<-create.an.edgeList(events)
  elist.sri <- vector()
  for(i in 1:nrow(elist)){

    Nab <- elist$n[i]

    Na <- elist %>% dplyr::filter(from==as.character(elist[i,]$from) ) %>% dplyr::filter(to!=as.character(elist[i,]$to) ) %>% summarise(x=sum(n))
    Nb <- elist %>% dplyr::filter(from==as.character(elist[i,]$to) ) %>% dplyr::filter(to!=as.character(elist[i,]$from) ) %>% summarise(x=sum(n))
    elist.sri[length(elist.sri)+1] <- Nab / (Nab + Na[[1]] + Nb[[1]])
  }

  elist$n <- elist.sri

  names(elist)<- c("from", "to", "weight")
  elist <- elist[complete.cases(elist),]
  gg <- graph_from_data_frame(elist, directed = TRUE, vertices = NULL)

  return(gg)
}

#' Create an edge list
#'
#' This function will generate an edge list from an events dataframe
#' @param events dataframe containing all events
#' @importFrom dplyr count
#'
create.an.edgeList<-function(events){

  edge.list.NN<-dplyr::count(events, from , to)

  return(edge.list.NN)
}

#' Create a window
#'
#' This function will generate a window from a larger event dataframe, based on the time of the events
#' @param df.total dataframe containing all events
#' @param start the starting time of the window
#' @param end the ending time of the window
#' @importFrom dplyr filter
#' @importFrom dplyr select
#'
create.window <- function(df.total, start, end){

  df.win <- dplyr::filter(df.total, time < end &  time >= start)
  window.sub <- dplyr::select(df.win, from, to )

  return (window.sub)
}

#' Create a window with a time column
#'
#' This function will generate a window from a larger event dataframe, based on the time of the events. The dataframe returned will have to, from, and time as columns.
#' @param df.total dataframe containing all events
#' @param start the starting time of the window
#' @param end the ending time of the window
#' @importFrom dplyr filter
#' @importFrom dplyr select
#'
create.window.time <- function(df.total, start, end){

  df.win <- dplyr::filter(df.total, time < end &  time >= start)
  window.sub <- dplyr::select(df.win, from, to, time )

  return (window.sub)
}


#' Bootstrap estimates
#'
#' This function will generate bootstrapped estimates for each network
#' @param dataSub dataframe containing events
#' @param nb number of bootstraps
#' @param type the measure to use for the bootstrap. Currently available are: between (i.e., mean betweenness centrality),
#'  close (i.e., mean closenes), eigen (i.e., eigen vector centrality), and cc (i.e., clustering coeficient). See the igraph package
#'  for details on each measure.
#' @param directedNet wheter the network is directed or not.
#' @param previousNet A second network to compare disimilariy (e.g., for cosine similariy).
#' @import igraph
#'
estimate.uncertainty.boot <- function(dataSub,nb, type, directedNet,previousNet=NULL){

  # Variables used during bootstrap
  bootstrap.number <- nb


  if(bootstrap.number>0){
    sample.size <- nrow(dataSub)

    #where to temporarily store bootstrapped samples from the original data
    boot.values <- vector('numeric')

    for (n in 1:bootstrap.number){

      #random sample with replacement from the observed interactions
      random.rows <- sample(sample.size, replace=TRUE)
      Boot.List <- dataSub[random.rows,]

      # Create the social network
      Boot.Network <- create.a.network(Boot.List)

      #calculate and store the network measure calculated from the bootstrapped sample
      if(type=='between')boot.values[length(boot.values)+1] <- mean(betweenness(Boot.Network,directed = directedNet))
      if(type=='eigen')boot.values[length(boot.values)+1] <- mean(eigen_centrality(Boot.Network,directed = directedNet))
      if(type=='close')boot.values[length(boot.values)+1] <- mean(closeness(Boot.Network))
      if(type=='cc')boot.values[length(boot.values)+1] <- transitivity(Boot.Network)
      if(type=='degree')boot.values[length(boot.values)+1] <- mean(degree(Boot.Network))
      if(type=='strength')boot.values[length(boot.values)+1] <- mean(strength(Boot.Network))
      if(type=='cosine')boot.values[length(boot.values)+1] <- cosine_between_graphs(Boot.Network,previousNet)
    }

    return (quantile(boot.values, probs = c(0.0275,0.5,0.975)))

  }else{
    return (c(NA,NA,NA))
  }

}

#' Permutation test
#'
#' This function will generate edge permutations for each network
#' @param graphW Network (igraph)
#' @param np number of permutations
#' @param type the measure to use for permutation. Currently available are: between (i.e., mean betweenness centrality),
#'  close (i.e., mean closenes), eigen (i.e., eigen vector centrality), and cc (i.e., clustering coeficient). See the igraph package
#'  for details on each measure.
#' @param directedNet wheter the network is directed or not.
#' @param previousNet A second network to compare disimilariy (e.g., for cosine similariy).
#' @import igraph
#'
estimate.random.range.perm <- function(graphW,np, type, directedNet,previousNet=NULL){

  #Parameters needed
  permutation.number <- np

  if(permutation.number>0){
    number.individuals <- vcount(graphW)
    number.edges <- ecount(graphW)
    weights <- E(graphW)$weight

    #store the permutation values
    permutation.values <- vector('numeric')

    for (n in 1:permutation.number){

      #create a random graph with the same number of nodes and edges as the observed graph
      Permute.Network <- sample_gnm(number.individuals, number.edges, directed = T,loops = F)
      E(Permute.Network)$weight <- sample(weights)

      # Calculate the metric
      if(type=="between")permutation.values[length(permutation.values)+1] <- mean(betweenness(Permute.Network,directed = directedNet))
      if(type=="eigen")permutation.values[length(permutation.values)+1] <- mean(eigen_centrality(Permute.Network,directed = directedNet))
      if(type=="close")permutation.values[length(permutation.values)+1] <- mean(closeness(Permute.Network))
      if(type=="cc")permutation.values[length(permutation.values)+1] <- transitivity(Permute.Network)
      if(type=="degree")permutation.values[length(permutation.values)+1] <- mean(degree(Permute.Network))
      if(type=="strength")permutation.values[length(permutation.values)+1] <- mean(strength(Permute.Network))
      if(type=="cosine")permutation.values[length(permutation.values)+1] <- cosine_between_graphs(Permute.Network,previousNet)
    }

    return (quantile(permutation.values, probs = c(0.0275,0.5,0.975)))

  }else{
    return (c(NA,NA,NA))
  }

}

#' graphTS function
#'
#' This function will take a dataframe with events between individuals/objects, and take network measures using a moving window approach.
#' @param event.data dataframe containing events between individuals/objects. This dataframe should have a 'to' and a 'from' column, as well as 'time' column representing the time between the start of observations and the events.
#' @param nBoot number of bootstrap estimates to take for each measure
#' @param nPerm number of permutation estimates to take for each measure
#' @param windowSize size of the window in which to make network measures (should be the same scale as the time column)
#' @param windowShift the amount to shift the moving window for each measure
#' @param windowStart The time of the first window. This should corespond to the first events.
#' @param type is the type of measure to take in each window. Currently available are: between (i.e., mean betweenness centrality),
#'  close (i.e., mean closenes), eigen (i.e., eigen vector centrality), and cc (i.e., clustering coeficient). See the igraph package
#'  for details on each measure.
#' @param lag The lag at which to calculate cosine similarity in network structure.
#' @param directedNet Whether the events are directed or no: true or false.
#' @param threshold minimum number of events to calculate a network measure (otherwise NA is produced)
#' @param startDate Optional argument used to define the start date of events. Use dd/mm/yyyy format (e.g., startDate = '24/07/2002').
#' @export
#' @import igraph
#' @importFrom lubridate dmy
#' @examples
#'
#' ts.out<-graphTS(event.data=groomEvents[1:200,])
#' graphTS.plot(ts.out)
#'
graphTS <- function (event.data,nBoot=100,nPerm=100,windowSize =30,windowShift= 1, type="cc",directedNet=T, threshold=30,windowStart=0, lag=1,startDate=NULL, method="none"){

  #intialize
  windowEnd=windowStart+windowSize
  netValues <- data.frame(t(rep(-1,12)))
  names(netValues)<-c(type,paste(type,".low95",sep=""),paste(type,".med50",sep=""),paste(type,".high95",sep=""),"perm.low95","perm.med50","perm.high95","nEvents","windowStart","windowEnd","windowStartDate","windowEndDate")
  gplist <- rep(list(NA),lag)

  if(windowEnd>max(event.data$time))print("Error: the window size is set larger than the max time difference")

  #for every window
  while (windowStart + windowSize<=max(event.data$time)) {

    #subset the data
    df.window<-create.window(event.data, windowStart, windowEnd)

    #if there is enough data in this window...
    if(nrow(df.window)>threshold ){

      #if there is no previous network, and the measure requires one
      if(is.na(gplist[1]) & type=='cosine'){

        if(method == "none"){
          gplist[[length(gplist)+1]] <- create.a.network(df.window)
        } else if (method == "SRI"){
          gplist[[length(gplist)+1]] <- create.a.network.SRI(df.window)
        }
        gplist<-gplist[-1]
        measure <- NA
        measure.uncertainty<-c(NA,NA,NA)
        measure.random<-c(NA,NA,NA)

      } else {

        #create a network
        if(method == "none"){
          g <- create.a.network(df.window)
        } else if (method == "SRI"){
          g <- create.a.network.SRI(df.window)
        }

        #calculate measure
        if(type=='between')measure <- mean(betweenness(g))
        if(type=='eigen')measure <- mean(eigen_centrality(g)$vector)
        if(type=='close')measure <- mean(closeness(g))
        if(type=='cc')measure <- transitivity(g)
        if(type=='degree')measure <- mean(degree(g))
        if(type=='strength')measure <- mean(strength(g))
        if(type=='cosine')measure <- cosine_between_graphs(graph1=g,graph2=gplist[[1]])
        if(type=='SRI')measure <- get_SRI(create.window.time(event.data, windowStart, windowEnd))

        #estimate uncertainty of measure (bootstrap)
        if(type=='cosine'){
          measure.uncertainty <- estimate.uncertainty.boot(df.window, nb=nBoot, type=type, directedNet = directedNet, previousNet = gplist[[1]])
        } else{
          measure.uncertainty <- estimate.uncertainty.boot(df.window, nb=nBoot, type=type, directedNet = directedNet)
        }

        #esitmate range of random (permutation)
        if(type=='cosine'){
          measure.random <- estimate.random.range.perm(g, np=nPerm, type=type, directedNet = directedNet, previousNet=gplist[[1]])
        } else {
          measure.random <- estimate.random.range.perm(g, np=nPerm, type=type, directedNet = directedNet)
        }

        #save last network
        gplist[[length(gplist)+1]] <- g
        gplist<-gplist[-1]

      }
    } else {
      measure <- NA
      measure.uncertainty<-c(NA,NA,NA)
      measure.random<-c(NA,NA,NA)
      gplist[[length(gplist)+1]] <- NA
      gplist<-gplist[-1]
    }

    #get window date range
    windowStartDate <- NA
    windowEndDate <- NA
    if(is.null(startDate)==FALSE){
      windowStartDate <- ymd(startDate) + days(windowStart)
      windowEndDate <- ymd(startDate) + days(windowEnd)
    }

    #record each measure as we go
    instance<-data.frame(measure, measure.uncertainty[1],measure.uncertainty[2], measure.uncertainty[3],measure.random[1],measure.random[2],measure.random[3],nrow(df.window),windowStart,windowEnd,windowStartDate, windowEndDate)
    names(instance) <- c(type,paste(type,".low95",sep=""),paste(type,".med50",sep=""),paste(type,".high95",sep=""),"perm.low95","perm.med50","perm.high95","nEvents","windowStart","windowEnd","windowStartDate","windowEndDate")
    netValues <- rbind(netValues,instance)

    #move window over
    windowEnd = windowEnd + windowShift
    windowStart = windowStart + windowShift
  }

  netValues<-netValues[-1,]
  return (netValues)
}

#' Plotting function for netTS dataframes
#'
#' This function will plot the output of the netTS function
#' @param df.ts output dataframe from the netTS function
#' @param nEvents Option to plot the number of events within each window.
#' @param dates Option to plot the figure with dates as opposed to days since start. Note this only applies if a start date was provided to the graphTS function.
#' @export
#' @import ggplot2
#' @examples
#'
#' ts.out<-graphTS(event.data=groomEvents[1:200,])
#' graphTS.plot(ts.out)
#'
graphTS.plot<-function(df.ts, nEvents=FALSE, dates=FALSE){

  if(nEvents==FALSE){

    if(dates==FALSE){

      fig<-ggplot(df.ts, aes(x=df.ts$windowEnd, y=df.ts[,1]))+ geom_line()+
        geom_ribbon(aes(ymin = df.ts[,2], ymax = df.ts[,4], fill="bootstrap"),alpha=0.2) +
        geom_ribbon(aes(ymin = df.ts[,5], ymax = df.ts[,7], fill="permutation"),alpha=0.2) +
        geom_point(color="blue") +
        labs(x= "Time since start", y=names(df.ts)[1])+
        scale_colour_manual(name="Shading", values=c(bootstrap="red", permutation="blue"))+theme_minimal()
    } else {
      fig<-ggplot(df.ts, aes(x=df.ts$windowEndDate, y=df.ts[,1]))+ geom_line()+
        geom_ribbon(aes(ymin = df.ts[,2], ymax = df.ts[,4], fill="bootstrap"),alpha=0.2) +
        geom_ribbon(aes(ymin = df.ts[,5], ymax = df.ts[,7], fill="permutation"),alpha=0.2) +
        geom_point(color="blue") +
        labs(x= "Time since start", y=names(df.ts)[1])+
        scale_colour_manual(name="Shading", values=c(bootstrap="red", permutation="blue"))+theme_minimal()
    }

  } else {

    if(dates==FALSE){
      fig<-ggplot(df.ts, aes(x=df.ts$windowEnd, y=df.ts$nEvents))+ geom_line()+
        geom_point(color="blue") +
        labs(x= "Time since start", y="Number of events")+theme_minimal()
    } else{
      fig<-ggplot(df.ts, aes(x=df.ts$windowEndDate, y=df.ts$nEvents))+ geom_line()+
        geom_point(color="blue") +
        labs(x= "Time since start", y="Number of events")+theme_minimal()
    }
  }

  fig
  return(fig)

}


#' Estimate consine similarity between graphs
#'
#' This function will calculate the cosine similarity between two graphs
#' @param graph1 first graph (igraph)
#' @param graph2 second graph (igraph)
#' @export
#' @import igraph
#' @importFrom dplyr full_join
#' @importFrom lsa cosine
#' @examples
#'
#' #two random graphs
#' library(igraph)
#' graph1 <-erdos.renyi.game(n=10,p=0.1)
#' graph2 <-erdos.renyi.game(n=10,p=0.1)
#' cosine_between_graphs(graph1,graph2)
#'
cosine_between_graphs <- function(graph1,graph2){

  #create weighted edge list from first graph
  g1.edges<-as.data.frame(get.edgelist(graph1, names=TRUE))
  if(is.null(igraph::E(graph1)$weight)){
    g1.edges$weight <- rep(1,nrow(g1.edges))
  } else {
    g1.edges$weight<-igraph::E(graph1)$weight
  }
  #create a join column (ordered)
  g1.edges$joinC <- ifelse(as.character(g1.edges$V1) < as.character(g1.edges$V2), paste(g1.edges$V1, g1.edges$V2), paste(g1.edges$V2, g1.edges$V1))

  #create weighted edge list from second graph
  g2.edges<-as.data.frame(get.edgelist(graph2, names=TRUE))
  if(is.null(igraph::E(graph2)$weight)){
    g2.edges$weight <- rep(1,nrow(g2.edges))
  } else {
    g2.edges$weight<-igraph::E(graph2)$weight
  }
  #create a join column (ordered)
  g2.edges$joinC <- ifelse(as.character(g2.edges$V1) < as.character(g2.edges$V2), paste(g2.edges$V1, g2.edges$V2), paste(g2.edges$V2, g2.edges$V1))

  #join the weighted edge lists, setting NAs equal to 0
  comb<-full_join(g1.edges,g2.edges,by="joinC")
  comb$weight.x[is.na(comb$weight.x)]<-0
  comb$weight.y[is.na(comb$weight.y)]<-0

  return(cosine(comb$weight.x,comb$weight.y))

}


#' Estimate SRI value for each dyad
#'
#' This function will calculate the SRI score for each dyad
#' @param graph1 first graph (igraph)
#' @param graph2 second graph (igraph)
#' @export
#' @import igraph
#' @importFrom dplyr full_join
#'
#' @examples
#'
#'
get_SRI <- function(df.window){

  #for each individual calculate SRI scores to others
  #@importFrom asnipe get_group_by_individual
  unique.ids<-unique(c(df.window$to,df.window$from))



}


