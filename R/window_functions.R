<<<<<<< HEAD
=======

>>>>>>> origin/master
#' Create a network from an edge list
#'
#' This function will generate a network from an event dataframe
#' @param events dataframe containing all events.
#'
create.a.network<-function(events){

  elist<-create.an.edgeList(events)
  names(elist)<- c("from", "to", "weight")
<<<<<<< HEAD
  gg <- igraph::graph_from_data_frame(elist, directed = TRUE, vertices = NULL)
=======
  gg <- graph_from_data_frame (elist, directed = TRUE, vertices = NULL)
>>>>>>> origin/master

  return(gg)
}

<<<<<<< HEAD
=======

>>>>>>> origin/master
#' Create an edge list
#'
#' This function will generate an edge list from an events dataframe
#' @param events dataframe containing all events
#'
create.an.edgeList<-function(events){

<<<<<<< HEAD
  edge.list.NN<-dplyr::count(events, from , to)
=======
  edge.list.NN<-dplyr::count(events, id1 , id2)
>>>>>>> origin/master

  return(edge.list.NN)
}

#' Create a window
#'
#' This function will generate a window from a larger event dataframe, based on the time of the events
#' @param df.total dataframe containing all events
#' @param start the starting time of the window
#' @param end the ending time of the window
#'
create.window <- function(df.total, start, end){

<<<<<<< HEAD
  df.win <- dplyr::filter(df.total, time < end &  time >= start)
  window.sub <- dplyr::select(df.win, from, to )
=======
  df.win <- filter(df.total, DaySinceStart < end &  DaySinceStart >= start)
  window.sub <- select(df.win, id , NN.Female )
  window.sub <- filter(window.sub, NN.Female!="XX" & NN.Female!="" & as.character(id)!= as.character(NN.Female))
>>>>>>> origin/master

  return (window.sub)
}


#' Bootstrap estimates
#'
#' This function will generate bootstrapped estimates for each network
#' @param dataSub dataframe containing events
#' @param nb number of bootstraps
#' @param type the measure to use for the bootstrap. Currently available are: betweennes, closeness, eigen, and cc (i.e., clustering coeficient)
#' @param directedNet wheter the network is directed or not.
#'
estimate.uncertainty.boot <- function(dataSub,nb, type, directedNet){

  # Variables used during bootstrap
  bootstrap.number <- nb
  sample.size <- nrow(dataSub)

  #where to temporarily store bootstrapped samples from the original data
  boot.values <- vector('numeric')

  for (n in 1:bootstrap.number){

    #random sample with replacement from the observed interactions
    random.rows <- sample(1:sample.size, sample.size, replace=T)
    Boot.List <- dataSub[random.rows,]

    # Create the social network
    Boot.Network <- create.a.network(Boot.List)

    #calculate and store the network measure calculated from the bootstrapped sample
<<<<<<< HEAD
    if(type=='betweennes')boot.values[length(boot.values)+1] <- mean(igraph::betweenness(Boot.Network,directed = directedNet))
    if(type=='eigen')boot.values[length(boot.values)+1] <- mean(igraph::eigen_centrality(Boot.Network,directed = directedNet))
    if(type=='closeness')boot.values[length(boot.values)+1] <- mean(igraph::closeness(Boot.Network,directed = directedNet))
    if(type=='cc')boot.values[length(boot.values)+1] <- igraph::transitivity(Boot.Network)
    if(type=='degree')boot.values[length(boot.values)+1] <- mean(igraph::degree(Boot.Network))
    if(type=='strength')boot.values[length(boot.values)+1] <- mean(igraph::strength(Boot.Network))
=======
    if(type=='betweennes')boot.values[length(boot.values)+1] <- mean(betweenness(Boot.Network,directed = directedNet))
    if(type=='eigen')boot.values[length(boot.values)+1] <- mean(eigen_centrality(Boot.Network,directed = directedNet))
    if(type=='closeness')boot.values[length(boot.values)+1] <- mean(closeness(Boot.Network,directed = directedNet))
    if(type=='cc')boot.values[length(boot.values)+1] <- transitivity(Boot.Network)
    if(type=='degree')boot.values[length(boot.values)+1] <- mean(degree(Boot.Network))
    if(type=='strength')boot.values[length(boot.values)+1] <- mean(strength(Boot.Network))
>>>>>>> origin/master
  }

  return (quantile(boot.values, probs = c(0.0275,0.5,0.975)))

}



#' Permutation test
#'
#' This function will generate edge permutations for each network
#' @param graphW Network (igraph)
#' @param np number of permutations
#' @param type the measure to use for permutation. Currently available are: betweennes, closeness, eigen, and cc (i.e., clustering coeficient)
#' @param directedNet wheter the network is directed or not.
#'
estimate.random.range.perm <- function(graphW,np, type, directedNet){

  #Parameters needed
  permutation.number <- np
<<<<<<< HEAD
  number.individuals <- igraph::vcount(graphW)
  number.edges <- igraph::ecount(graphW)
  weights <- igraph::E(graphW)$weight
=======
  number.individuals <- vcount(graphW)
  number.edges <- ecount(graphW)
  weights <- E(graphW)$weight
>>>>>>> origin/master

  #store the permutation values
  permutation.values <- vector('numeric')

  for (n in 1:permutation.number){

    #create a random graph with the same number of nodes and edges as the observed graph
<<<<<<< HEAD
    Permute.Network <- igraph::sample_gnm(number.individuals, number.edges, directed = T,loops = F)
    igraph::E(Permute.Network)$weight <- sample(weights)

    # Calculate the metric
    if(type=="betweennes")permutation.values[length(permutation.values)+1] <- mean(igraph::betweenness(Permute.Network,directed = directedNet))
    if(type=="eigen")permutation.values[length(permutation.values)+1] <- mean(igraph::eigen_centrality(Permute.Network,directed = directedNet))
    if(type=="closeness")permutation.values[length(permutation.values)+1] <- mean(igraph::closeness(Permute.Network,directed = directedNet))
    if(type=="cc")permutation.values[length(permutation.values)+1] <- igraph::transitivity(Permute.Network)
    if(type=="degree")permutation.values[length(permutation.values)+1] <- mean(igraph::degree(Permute.Network))
    if(type=="strength")permutation.values[length(permutation.values)+1] <- mean(igraph::strength(Permute.Network))
=======
    Permute.Network <- sample_gnm(number.individuals, number.edges, directed = T,loops = F)
    E(Permute.Network)$weight <- sample(weights)

    # Calculate the metric
    if(type=="betweennes")permutation.values[length(permutation.values)+1] <- mean(betweenness(Permute.Network,directed = directedNet))
    if(type=="eigen")permutation.values[length(permutation.values)+1] <- mean(eigen_centrality(Permute.Network,directed = directedNet))
    if(type=="closeness")permutation.values[length(permutation.values)+1] <- mean(closeness(Permute.Network,directed = directedNet))
    if(type=="cc")permutation.values[length(permutation.values)+1] <- transitivity(Permute.Network)
    if(type=="degree")permutation.values[length(permutation.values)+1] <- mean(degree(Permute.Network))
    if(type=="strength")permutation.values[length(permutation.values)+1] <- mean(strength(Permute.Network))
>>>>>>> origin/master
  }

  return (quantile(permutation.values, probs = c(0.0275,0.5,0.975)))

}

#' netWin function
#'
#' This function will take a dataframe with events between individuals/objects, and take network measures using a moving window approach.
#' A time column is required.
#' @param event.data dataframe containing events between individuals/objects
#' @param nBoot number of bootstrap estimates to take for each measure
#' @param nPerm number of permutation estimates to take for each measure
#' @param windowSize size of the window in which to make network measures (should be the same scale as the time column)
#' @param windowShift the amount to shift the moving window for each measure
#' @param type is the type of measure to take in each window. Currently available are: betweennes, closeness, eigen, and cc (i.e., clustering coeficient)
#' @param directedNet Whether the events are directed or no: true or false.
#' @param threshold minimum number of events to calculate a network measure (otherwise NA is produced)
<<<<<<< HEAD
#' @export
#' @examples
#'
#' #not run
#' library(netTS)
#' ts.out<-netWin(event.data=df)
=======
#' @examples
#'
#' #not run
#' ts.out<-netTS(event.data=df)
>>>>>>> origin/master
#' plot.netTS(ts.out)
#'
netWin <- function (event.data,nBoot=100,nPerm=100,windowSize =30,windowShift= 1, type="cc",directedNet=T, threshold=30){

  #intialize
  windowStart=0
  windowEnd=windowStart+windowSize
  netValues <- data.frame()

  if(windowEnd>max(event.data$time))print("Error: the window size is set larger than the max time difference")

  #for every window
  while (windowStart + windowSize<=max(event.data$time)) {

    #subset the data
    df.window<-create.window(event.data, windowStart, windowEnd)

    #if there is data in this window...
    if(nrow(df.window)>threshold){

      #create a network
      g <- create.a.network(df.window)

      #calculate measure
<<<<<<< HEAD
      if(type=='betweennes')measure <- mean(igraph::betweenness(g))
      if(type=='eigne')measure <- mean(igraph::eigen_centrality(g))
      if(type=='closeness')measure <- mean(igraph::closeness(g))
      if(type=='cc')measure <- igraph::transitivity(g)
      if(type=='degree')measure <- mean(igraph::degree(g))
      if(type=='strength')measure <- mean(igraph::strength(g))
=======
      if(type=='betweennes')measure <- mean(betweenness(g))
      if(type=='eigne')measure <- mean(eigen_centrality(g))
      if(type=='closeness')measure <- mean(closeness(g))
      if(type=='cc')measure <- transitivity(g)
      if(type=='degree')measure <- mean(degree(g))
      if(type=='strength')measure <- mean(strength(g))
>>>>>>> origin/master

      #estimate uncertainty of measure (bootstrap)
      measure.uncertainty <- estimate.uncertainty.boot(df.window, nb=nBoot, type=type, directedNet = directedNet)

      #esitmate range of random (permutation)
      measure.random <- estimate.random.range.perm(g, np=nPerm, type=type, directedNet = directedNet)

    } else {
      measure <- NA
      measure.uncertainty<-c(NA,NA,NA)
      measure.random<-c(NA,NA,NA)
    }

    #record each measure as we go
    netValues <- rbind(netValues,c(measure, measure.uncertainty,measure.random, windowStart,windowEnd))

    #move window over
    windowEnd = windowEnd + windowShift
    windowStart = windowStart + windowShift

  }

  names(netValues)<-c(type,paste(type,".low95",sep=""),paste(type,".med50",sep=""),paste(type,".high95",sep=""),"perm.low95","perm.med50","perm.high95","windowStart","windowEnd")
  return (netValues)
}

#' Plotting function for netTS dataframes
#'
#' This function will plot the output of the netTS function
#' @param df output dataframe from the netTS function
<<<<<<< HEAD
#' @export
#' @examples
#'
#' #not run
#' library(netTS)
#' ts.out<-netWin(event.data=df)
=======
#' @examples
#'
#' #not run
#' ts.out<-netTS(event.data=df)
>>>>>>> origin/master
#' plot.netTS(ts.out)
#'
plot.netTS<-function(df){

<<<<<<< HEAD
  fig<-ggplot2::ggplot(df, aes(x=df[,8], y=df[,1]))+ geom_line()+
    geom_ribbon(aes(ymin = df[,2], ymax = df[,4], fill="bootstrap"),alpha=0.2) +
    geom_ribbon(aes(ymin = df[,5], ymax = df[,7], fill="permutation"),alpha=0.2) +
    geom_point(color="blue") +
    labs(x= "Time since start", y=names(df)[1])+
    scale_colour_manual(name="Shading", values=c(bootstrap="red", permutation="blue"))+theme_minimal()
=======
  fig<-ggplot(df, aes(x=df[,8], y=df[,1]))+ geom_line()+
    geom_ribbon(aes(ymin = df[,2], ymax = df[,4], fill="bootstrap"),alpha=0.2) +
    geom_ribbon(aes(ymin = df[,5], ymax = df[,7], fill="permutation"),alpha=0.2) +
    geom_point(color="blue") +
    labs(x= "Days since start", y=names(df)[1])+
    scale_colour_manual(name="Shading", values=c(bootstrap="red", permutation="blue"))
>>>>>>> origin/master

  fig
  return(fig)

}
