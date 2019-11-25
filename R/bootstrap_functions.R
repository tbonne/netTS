

#' Bootstrap convergence check with subsampling
#'
#' This function will estimate the convergence of the chosen network measure using bootstrapped samples of the data, while also checking for measurement sensitivity to data subsampling.
#' @param data Dataframe with relational data in the first two rows, and a time stamp in the third row. Note: time stamps should be in ymd or ymd_hms format. The lubridate package can be very helpful in organizing times.
#' @param windowsize The size of each window in which to generate a network.
#' @param windowshift The amount of time to shift the window when generating networks.
#' @param directed Whether to consider the network as directed or not (TRUE/FALSE).
#' @param measureFun The measurment function to be used with the bootstapped networks.
#' @param corFun The method used to compare observed values with bootstrapped values: 1-Cosine similarity, 2-pearsons correlation, 3-Euclidean distance
#' @param boot.samples The number of bootstrapped samples to run (Default=100)
#' @param SRI Wether to use the simple ratio index (Default=FALSE)
#' @param probs The quantiles of the bootrap samples to return (Default=c(0.025,0.975)).
#' @param subsamples A vector of values between 0-1 used to subsample the original dataframe.
#' @param plot Wether a plot of the results should be produced.
#' @importFrom stats cor.test quantile
#' @importFrom igraph set_graph_attr degree
#' @export
#'
#'
check.windowsize <- function(data, windowsize=days(30), windowshift=days(1), directed = FALSE, measureFun=degree, corFun = 1,boot.samples=100, SRI=FALSE, probs=c(0.025,0.975), subsamples=c(1,0.8,0.6), plot=TRUE){

  #dataframe to store the results
  df.results <- data.frame(mean=-1,CI.low=-1,CI.high=-1,windowstart=as.Date("2001-12-30"),windowend=as.Date("2001-12-30"), fracData = -1)

  for(i in subsamples){

    #subsample data
    data.sub = data[sample(1:nrow(data), size = round(i*nrow(data),digits = 0), replace = FALSE),]

    #run the bootstrap
    df.temp<-convergence.check.boot(data.sub, windowsize,windowshift,directed,measureFun,corFun,boot.samples,SRI,probs, fullData = data)

    #record the data
    df.temp$fracData = i
    df.results <- rbind(df.results, df.temp)
  }

  #remove the first row used to initialize the dataframe
  df.results<-df.results[-1,]

  #plot
  if(plot==T)window.check.plot(df.results)

  #return the results
  return(df.results)

}

#' Bootstrap convergence check
#'
#' This function will estimate the convergence of the chosen network measure using bootstrapped samples of the data.
#' @param data Dataframe with relational data in the first two rows, and a time stamp in the third row. Note: time stamps should be in ymd or ymd_hms format. The lubridate package can be very helpful in organizing times.
#' @param windowsize The size of each window in which to generate a network.
#' @param windowshift The amount of time to shift the window when generating networks.
#' @param directed Whether to consider the network as directed or not (TRUE/FALSE).
#' @param measureFun The measurment function to perform the bootstap on (should be at the node level).
#' @param effortFun This is a function that takes as input the data within a window of time and returns the total sampling effort.
#' @param effortData This is a dataframe containing the data used to calculate sampling effort.
#' @param corFun The method used to compare observed node/dyad values with bootstrapped values: 1-Cosine similarity, 2-pearsons correlation, 3-Euclidean distance
#' @param boot.samples The number of bootstrapped samples to run (Default=100)
#' @param SRI Wether to use the simple ratio index (Default=FALSE)
#' @param probs The quantiles of the bootrap samples to return (Default=c(0.025,0.975)).
#' @param fullData This is the full dataset, if a subset dataset is being used to compare bootstrap samples to the full dataset.
#' @importFrom stats cor.test quantile
#' @importFrom igraph set_graph_attr degree
#' @export
#'
#'
convergence.check.boot <- function(data, windowsize=days(30), windowshift=days(1), directed = FALSE, measureFun=degree, corFun = 1,boot.samples=100, SRI=FALSE, probs=c(0.025,0.975), effortFun=NULL, effortData=NULL,fullData=NULL){

  #intialize times
  windowstart <- min(data[,3])
  windowend=windowstart+windowsize
  if(windowend>max(data[,3]))print("warnning: the window size is set larger than the observed data.")

  #for every window generate a network
  conv.values <- data.frame(mean=-1,CI.low=-1,CI.high=-1,windowstart=as.Date("2001-12-30"),windowend=as.Date("2001-12-30"))

  while (windowstart + windowsize<=max(data[,3])) {

    if(is.null(fullData)){

      #subset the data
      df.window<-create.window(data, windowstart, windowend)
      Observation.Events <- nrow(df.window)

      #calculate sampling effort
      if(is.null(effortFun)==FALSE & is.null(effortData)==TRUE  ){ #there is an effort function and it requires no external data
        g = effortFun(df.window, directed = directed)
        g <- set_graph_attr(g, "nEvents", Observation.Events)
        g <- set_graph_attr(g, "windowstart", windowstart )
        g <- set_graph_attr(g, "windowend", windowend)

      }else if(is.null(effortFun)==FALSE & is.null(effortData)==FALSE ){ #there is an effort function and it requires some external data
        effortData.sub <- effortData[effortData[,1]>=windowstart & effortData[,1]<windowend,]
        g = effortFun(df.window, effortData.sub, directed = directed)
        g <- set_graph_attr(g, "nEvents", Observation.Events)
        g <- set_graph_attr(g, "windowstart", windowstart )
        g <- set_graph_attr(g, "windowend", windowend)

      } else { #there is no effort function

        #calculate the desired network measure
        g <- create.a.network(df.window, directed = directed, SRI=FALSE)
        g <- set_graph_attr(g, "nEvents", Observation.Events)
        g <- set_graph_attr(g, "windowstart", windowstart)
        g <- set_graph_attr(g, "windowend", windowend)

      }

      #store measure
      obs.measures <- measureFun(g)

    } else {

      #subset the data
      df.window.full<-create.window(fullData, windowstart, windowend)
      df.window<-create.window(data, windowstart, windowend)
      Observation.Events.full <- nrow(df.window.full)
      Observation.Events <- nrow(df.window)

      #calculate sampling effort
      if(is.null(effortFun)==FALSE & is.null(effortData)==TRUE  ){ #there is an effort function and it requires no external data
        g = effortFun(df.window.full, directed = directed)
        g <- set_graph_attr(g, "nEvents", Observation.Events)
        g <- set_graph_attr(g, "windowstart", windowstart )
        g <- set_graph_attr(g, "windowend", windowend)

      }else if(is.null(effortFun)==FALSE & is.null(effortData)==FALSE ){ #there is an effort function and it requires some external data
        effortData.sub <- effortData[effortData[,1]>=windowstart & effortData[,1]<windowend,]
        g = effortFun(df.window.full, effortData.sub, directed = directed)
        g <- set_graph_attr(g, "nEvents", Observation.Events)
        g <- set_graph_attr(g, "windowstart", windowstart )
        g <- set_graph_attr(g, "windowend", windowend)

      } else { #there is no effort function

        #calculate the desired network measure
        g.full <- create.a.network(df.window.full, directed = directed, SRI=FALSE)
        g.full <- set_graph_attr(g.full, "nEvents", Observation.Events)
        g.full <- set_graph_attr(g.full, "windowstart", windowstart)
        g.full <- set_graph_attr(g.full, "windowend", windowend)

      }

      #store measure
      obs.measures <- measureFun(g.full)

    }

    if(Observation.Events>0){

      if(length(obs.measures)<2){
        corFun=3
        if(corFun!=3)print("Warning: only one measure prduced by the measurement function. corFun set to euclidean distance")
      }

      #store correlation measures
      cor.measures <- vector()

      #Number of
      for(j in 1:boot.samples){

        #bootstrap sample from the data in this window
        df.sub<-df.window[sample(nrow(df.window),replace = T),]

        #calculate sampling effort
        if(is.null(effortFun)==FALSE & is.null(effortData)==TRUE  ){ #there is an effort function and it requires no external data
          g = effortFun(df.sub, directed = directed)
          g <- set_graph_attr(g, "nEvents", Observation.Events)
          g <- set_graph_attr(g, "windowstart", windowstart )
          g <- set_graph_attr(g, "windowend", windowend)

        }else if(is.null(effortFun)==FALSE & is.null(effortData)==FALSE ){ #there is an effort function and it requires some external data
          effortData.sub <- effortData[effortData[,1]>=windowstart & effortData[,1]<windowend,]
          g = effortFun(df.sub, effortData.sub, directed = directed)
          g <- set_graph_attr(g, "nEvents", Observation.Events)
          g <- set_graph_attr(g, "windowstart", windowstart )
          g <- set_graph_attr(g, "windowend", windowend)

        } else { #there is no effort function

          #create a network and add it to the list
          g.boot <- create.a.network(df.sub, directed = directed, SRI=SRI)
          g.boot <- set_graph_attr(g.boot, "nEvents", Observation.Events)
          g.boot <- set_graph_attr(g.boot, "windowstart", windowstart)
          g.boot <- set_graph_attr(g.boot, "windowend", windowend)

        }

        #take correlation measure

        if(corFun==2){
          comb.by.names<-cbind(obs.measures,measureFun(g.boot)[names(obs.measures)])
          comb.by.names[is.na(comb.by.names[,2]),2]<-0
          cor.measures[length(cor.measures)+1] <- cor.test(comb.by.names[,1],comb.by.names[,2])$estimate
        } else if(corFun == 1){
          comb.by.names<-cbind(obs.measures,measureFun(g.boot)[names(obs.measures)])
          comb.by.names[is.na(comb.by.names[,2]),2]<-0
          cor.measures[length(cor.measures)+1] <- lsa::cosine(comb.by.names[,1]-mean(comb.by.names[,1]),comb.by.names[,2]-mean(comb.by.names[,2]))
        }else if(corFun == 3){
          cor.measures[length(cor.measures)+1] <- sqrt(sum((measureFun(g.boot)-obs.measures) ^ 2))
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
      g <- create.a.network(df.window, directed = directed, SRI=SRI)
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
        g.boot <- create.a.network(df.sub, directed = directed, SRI=SRI)
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
#' @param max.subsample.size The maximum number of samples to remove when recalculating the network metric.
#' @param SRI Wether to use the simple ratio index or not (Default=FALSE).
#' @param n.boot Number of times to repeate the subsample procedure (i.e., bootstrap estimate of slope)
#' @param probs The probablility interval to calculate the confidence intervals.
#' @importFrom stats lm complete.cases coef sigma
#' @importFrom igraph set_graph_attr
#' @importFrom bootstrap bcanon
#' @export
#'
#'
convergence.check.value<-function(data, windowsize, windowshift, directed = FALSE, measureFun,max.subsample.size=30, SRI=FALSE, n.boot=100, probs=c(0.025,0.975)){

  #intialize times
  windowstart <- min(data[,3])
  windowend=windowstart+windowsize
  if(windowend>max(data[,3]))print("warnning: the window size is set larger than the observed data.")

  #for every window generate a network
  conv.values <- data.frame(slope=-1, CI.low=-1, CI.high=-1, windowstart=as.Date("2001-12-30"),windowend=as.Date("2001-12-30"))
  while (windowstart + windowsize<=max(data[,3])) {

    #subset the data
    df.window<-create.window(data, windowstart, windowend)
    Observation.Events <- nrow(df.window)

    #store the slope values
    estimated.slopes <- vector()

    #bootstrap samples to calculate mean and CI for slope
    for(k in 1:n.boot){

      #store network measures
      net.measures <- data.frame(value=-1,sample=-1)

      if(Observation.Events>1){

        #Number of
        for(j in seq(max(Observation.Events-max.subsample.size,1),Observation.Events,by=1)){

          #subset window
          df.sub<-df.window[sample(nrow(df.window),j, replace = FALSE),]

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

      #clean up the data
      net.measures<-net.measures[-1,]
      net.measures<-net.measures[complete.cases(net.measures),]

      #calculate slope of the line (effect of subsampling on the network measure)
      if(nrow(net.measures)>0){
        estimated.slopes[length(estimated.slopes)+1] <- coef(lm(value~sample, data = net.measures))["sample"]
      } else {
        estimated.slopes[length(estimated.slopes)+1] <- NA
      }

    }

    #calculate and store estiamtes from the bootstrapped sample
    if(abs(max(estimated.slopes, na.rm = T) - min(estimated.slopes, na.rm = T)) < .Machine$double.eps ^ 0.5){

      print("Note: no range in estimated network measure produced by sub-sampling")
      conv.values <- rbind(conv.values, data.frame(slope=mean(estimated.slopes), lci = mean(estimated.slopes),uci = mean(estimated.slopes),windowstart=windowstart,windowend=windowend ))

    } else{

      bca = bcanon(estimated.slopes,n.boot,mean,alpha=probs)
      conv.values <- rbind(conv.values, data.frame(slope=mean(estimated.slopes), CI.low = bca$conf[1,2],CI.high =  bca$conf[2,2],windowstart=windowstart,windowend=windowend ))

    }

    #move the window
    windowend = windowend + windowshift
    windowstart = windowstart + windowshift

  }

  #remove first row used for intialization
  conv.values<-conv.values[-1,]

  return(conv.values)

}
