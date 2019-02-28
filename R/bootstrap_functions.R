

#' Bootstrap convergence check
#'
#' This function will estimate the convergence of the chosen network measure using bootstrapped samples of the data.
#' @param data Dataframe with relational data in the first two rows, and a time stamp in the third row. Note: time stamps should be in ymd or ymd_hms format. The lubridate package can be very helpful in organizing times.
#' @param windowsize The size of each window in which to generate a network.
#' @param windowshift The amount of time to shift the window when generating networks.
#' @param directed Whether to consider the network as directed or not (TRUE/FALSE).
#' @param measureFun The measurment function to perform the bootstap on (should be at the node level).
#' @param corFun The method used to compare observed node/dyad values with bootstrapped values: 1-Cosine similarity, 2-pearsons correlation
#' @param boot.samples The number of bootstrapped samples to run (Default=100)
#' @param SRI Wether to use the simple ratio index (Default=FALSE)
#' @param probs The quantiles of the bootrap samples to return (Default=c(0.025,0.975)).
#' @importFrom stats cor.test quantile
#' @importFrom igraph set_graph_attr degree
#' @export
#'
#'
convergence.check.boot <- function(data, windowsize=days(30), windowshift=days(1), directed = FALSE, measureFun=degree, corFun = 1,boot.samples=100, SRI=FALSE, probs=c(0.025,0.975)){

  #intialize times
  windowstart <- min(data[,3])
  windowend=windowstart+windowsize
  if(windowend>max(data[,3]))print("warnning: the window size is set larger than the observed data.")

  #for every window generate a network
  conv.values <- data.frame(mean=-1,CI.low=-1,CI.high=-1,windowstart=as.Date("2001-12-30"),windowend=as.Date("2001-12-30"))
  while (windowstart + windowsize<=max(data[,3])) {

    #subset the data
    df.window<-create.window(data, windowstart, windowend)
    Observation.Events <- nrow(df.window)

    if(Observation.Events>0){

      #calculate the desired network measure
      g <- create.a.network(df.window, directed = directed, SRI=FALSE, effort = 1)
      g <- set_graph_attr(g, "nEvents", Observation.Events)
      g <- set_graph_attr(g, "windowstart", windowstart)
      g <- set_graph_attr(g, "windowend", windowend)

      #store measure
      obs.measures <- measureFun(g)

      if(length(obs.measures)<2)print("Warning: the convergence check with boot only works with network measures that return multiple values: e.g., one value per node/dyad.")

      #store correlation measures
      cor.measures <- vector()

      #Number of
      for(j in 1:boot.samples){

        #bootstrap sample from the data in this window
        df.sub<-df.window[sample(nrow(df.window),replace = T),]

        #create a network and add it to the list
        g.boot <- create.a.network(df.sub, directed = directed, SRI=SRI, effort = 1)
        g.boot <- set_graph_attr(g.boot, "nEvents", Observation.Events)
        g.boot <- set_graph_attr(g.boot, "windowstart", windowstart)
        g.boot <- set_graph_attr(g.boot, "windowend", windowend)

        #take correlation measure
        comb.by.names<-cbind(obs.measures,measureFun(g.boot)[names(obs.measures)])
        comb.by.names[is.na(comb.by.names[,2]),2]<-0
        if(corFun==2){
          cor.measures[length(cor.measures)+1] <- cor.test(comb.by.names[,1],comb.by.names[,2])$estimate
        } else if(corFun == 1){
          cor.measures[length(cor.measures)+1] <- lsa::cosine(comb.by.names[,1],comb.by.names[,2])
        }
      }

      #calculate convergence
      conv.values <- rbind(conv.values,data.frame(mean=mean(cor.measures),CI.low=quantile(cor.measures,probs = probs[1],na.rm = T),CI.high=quantile(cor.measures,probs = probs[2],na.rm = T),windowstart=windowstart, windowend=windowend))

    } else{
      #calculate convergence
      conv.values <- rbind(conv.values,data.frame(mean=NA,CI.low=NA,CI.high=NA,windowstart=windowstart, windowend=windowend))
    }



    #move the window
    windowend = windowend + windowshift
    windowstart = windowstart + windowshift

  }

  conv.values<-conv.values[-1,]

  return(conv.values)

}

#' Bootstrap convergence check for network level measures
#'
#' This function will estimate the convergence of the chosen network measure using bootstrapped samples of the data.
#' @param data Dataframe with relational data in the first two rows, and a time stamp in the third row. Note: time stamps should be in ymd or ymd_hms format. The lubridate package can be very helpful in organizing times.
#' @param windowsize The size of each window in which to generate a network.
#' @param windowshift The amount of time to shift the window when generating networks.
#' @param directed Whether to consider the network as directed or not (TRUE/FALSE).
#' @param measureFun The measurment function to perform the bootstap on (should be at the node level).
#' @param boot.samples The number of bootstrapped samples to run (Default=100)
#' @param SRI Wether to use the simple ratio index (Default=FALSE)
#' @param probs The quantiles of the bootrap samples to return (Default=c(0.025,0.975)).
#' @importFrom stats cor.test quantile
#' @importFrom igraph set_graph_attr degree
#' @export
#'
#'
convergence.check.boot.graph <- function(data, windowsize=days(30), windowshift=days(1), directed = FALSE, measureFun=eigen_mean,boot.samples=100, SRI=FALSE, probs=c(0.025,0.975)){

  #intialize times
  windowstart <- min(data[,3])
  windowend=windowstart+windowsize
  if(windowend>max(data[,3]))print("warnning: the window size is set larger than the observed data.")

  #for every window generate a network
  conv.values <- data.frame(mean=-1,sd=-1,CI.low=-1,CI.high=-1,windowstart=as.Date("2001-12-30"),windowend=as.Date("2001-12-30"))
  while (windowstart + windowsize<=max(data[,3])) {

    #subset the data
    df.window<-create.window(data, windowstart, windowend)
    Observation.Events <- nrow(df.window)

    if(Observation.Events>0){

      #calculate the desired network measure
      g <- create.a.network(df.window, directed = directed, SRI=SRI, effort = 1)
      g <- set_graph_attr(g, "nEvents", Observation.Events)
      g <- set_graph_attr(g, "windowstart", windowstart)
      g <- set_graph_attr(g, "windowend", windowend)

      #store measure
      obs.measures <- measureFun(g)

      if(length(obs.measures)>1)print("Warning: the convergence check with boot for graphs only works with network measures that return a single value.")

      #store network measure
      net.measures <- vector()

      #Number of
      for(j in 1:boot.samples){

        #bootstrap sample from the data in this window
        df.sub<-df.window[sample(nrow(df.window),replace = T),]

        #create a network and add it to the list
        g.boot <- create.a.network(df.sub, directed = directed, SRI=SRI, effort = 1)
        g.boot <- set_graph_attr(g.boot, "nEvents", Observation.Events)
        g.boot <- set_graph_attr(g.boot, "windowstart", windowstart)
        g.boot <- set_graph_attr(g.boot, "windowend", windowend)

        #take network measure
        net.measures[length(net.measures)+1] <- measureFun(g.boot)
      }

      #calculate convergence
      conv.values <- rbind(conv.values,data.frame(mean=mean(net.measures),sd=sd(net.measures),CI.low=quantile(net.measures,probs = probs[1],na.rm = T),CI.high=quantile(net.measures,probs = probs[2],na.rm = T),windowstart=windowstart, windowend=windowend))

    } else{
      #calculate convergence
      conv.values <- rbind(conv.values,data.frame(mean=NA,sd=NA,CI.low=NA,CI.high=NA,windowstart=windowstart, windowend=windowend))
    }



    #move the window
    windowend = windowend + windowshift
    windowstart = windowstart + windowshift

  }

  conv.values<-conv.values[-1,]

  return(conv.values)

}


#' Variance by window size check
#'
#' This function will estimate the variance in the time series based on the window size.
#' @param data Dataframe with relational data in the first two rows, and a time stamp in the third row. Note: time stamps should be in ymd or ymd_hms format. The lubridate package can be very helpful in organizing times.
#' @param windowsize_min The min size of windowsize to test.
#' @param windowsize_max The max size of windowsize to test.
#' @param by The resolution at which to test window sizes between the min and the max window sizes
#' @param windowshift The amount of time to shift the window when generating networks.
#' @param directed Whether to consider the network as directed or not (TRUE/FALSE).
#' @param measureFun The measurment function to perform the bootstap on (Default: density).
#' @param SRI Wether to use the simple ratio index (Default=FALSE)
#' @importFrom stats cor.test quantile
#' @importFrom igraph set_graph_attr degree edge_density
#' @export
#'
#'
convergence.check.var<-function(data, windowsize_min=days(10),windowsize_max=days(40),by=days(1), windowshift=days(1), directed = FALSE, measureFun=igraph::edge_density, SRI=FALSE){

  #setup dataframe to return variance values
  df.var <- data.frame(windowsize = -1, var=-1)

  #for each window size calculate the overall variance and add it to df.var
  windowsize_seq = windowsize_min
  while(windowsize_seq<=windowsize_max){

    #calculate time series
    graph.values<-graphTS(data, windowsize = windowsize_seq, windowshift= windowshift, measureFun=measureFun ,effortFun=NULL, permutationFun=perm.events,directed=directed, lagged=FALSE, lag=1, firstNet=FALSE, cores=1, nperm=0, probs=0.95, SRI=SRI)

    #record
    df.var<-rbind(df.var, data.frame(windowsize=as.numeric(as.duration(windowsize_seq),"days"), var=var(graph.values[,1], na.rm = T)) )

    #update window size tested
    windowsize_seq = windowsize_seq + by
  }

  df.var<-df.var[-1,]

  return(df.var)

}

#' Convergence check
#'
#' This function will estimate the convergence of the chosen network measure using random subsets of the data.
#' @param data Dataframe with relational data in the first two rows, and a time stamp in the third row. Note: time stamps should be in ymd or ymd_hms format. The lubridate package can be very helpful in organizing times.
#' @param windowsize The size of each window in which to generate a network.
#' @param windowshift The amount of time to shift the window when generating networks.
#' @param directed Whether to consider the network as directed or not (TRUE/FALSE).
#' @param measureFun The function used to take network measurements.
#' @param random.sample.size The maximum number of random samples to remove when recalculating the network metric.
#' @param SRI Wether to use the simple ratio index or not (Default=FALSE).
#' @importFrom stats lm complete.cases coef sigma
#' @importFrom igraph set_graph_attr
#' @export
#'
#'
convergence.check.value<-function(data, windowsize, windowshift, directed = FALSE, measureFun,random.sample.size=30, SRI=FALSE){

  #intialize times
  windowstart <- min(data[,3])
  windowend=windowstart+windowsize
  if(windowend>max(data[,3]))print("warnning: the window size is set larger than the observed data.")

  #for every window generate a network
  conv.values <- data.frame(slope=-1, sd=-1, windowstart=as.Date("2001-12-30"),windowend=as.Date("2001-12-30"))
  while (windowstart + windowsize<=max(data[,3])) {

    #subset the data
    df.window<-create.window(data, windowstart, windowend)
    Observation.Events <- nrow(df.window)

    #store network measures
    net.measures <- data.frame(value=-1,sample=-1)

    if(Observation.Events>1){

      #Number of
      for(j in seq(max(Observation.Events-random.sample.size,1),Observation.Events,by=1)){

        #subset window
        df.sub<-df.window[sample(nrow(df.window),j),]

        #create a network and add it to the list
        g <- create.a.network(df.sub, directed = directed, SRI)
        g <- set_graph_attr(g, "nEvents", Observation.Events)
        g <- set_graph_attr(g, "windowstart", windowstart)
        g <- set_graph_attr(g, "windowend", windowend)

        #take measure
        net.measures <- rbind(net.measures,data.frame(value=measureFun(g),sample=j))
      }

    } else {
      net.measures <- rbind(net.measures,data.frame(value=NA,sample=0))
    }

    net.measures<-net.measures[-1,]
    net.measures<-net.measures[complete.cases(net.measures),]
    #calculate convergence (right now just using the slope...)

    if(nrow(net.measures)>0){
      conv.values <- rbind(conv.values,data.frame( slope = coef(lm(value~sample, data = net.measures))["sample"], sd=sigma(lm(value~sample, data = net.measures)), windowstart=windowstart,windowend=windowend))
    } else {
      conv.values <- rbind(conv.values,data.frame( slope = NA, sd=NA, windowstart=windowstart,windowend=windowend))
    }

    #move the window
    windowend = windowend + windowshift
    windowstart = windowstart + windowshift


  }

  conv.values<-conv.values[-1,]

  return(conv.values)

}
