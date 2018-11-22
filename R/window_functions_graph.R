

#' graphTS function
#'
#' This function will take a dataframe with events between individuals/objects, and take network measures using a moving window approach.
#' @param data A dataframe with relational data in the first two rows, with weights in the thrid row, and a time stamp in the fourth row. Note: time stamps should be in ymd or ymd_hms format. The lubridate package can be very helpful in organizing times.
#' @param windowsize The size of the moving window in which to take network measures. These should be provided as e.g., days(30), hours(5), ... etc.
#' @param windowshift The amount to shift the moving window for each measure. Again times should be provided as e.g., days(1), hours(1), ... etc.
#' @param measureFun This is a function that takes as an input a igraph network and returns a single value. There are functions within netTS (see details), and custom made functions can be used.
#' @param directed Whether the events are directed or no: true or false.
#' @param lagged Whether the network measure function used requires the comparison between two networks. e.g., comparing the current network to one lagged by 10 days. If TRUE the measureFun should take two graphs as input and return a single value. The order of inputs in the function is the lagged network followed by the current network.
#' @param lag If lagged is set to TRUE, this is the lag at which to compare networks.
#' @param firstNet If lagged is set to TRUE, this forces the comparisons between graphs to always be between the current and first graph.
#' @param cores This allows for multiple cores to be used while generating networks and calculating network measures.
#' @param nperm This allows for the estimation the network measure assuming random permutations. Currently the 95 percent quantiles are returned.
#' @param permutationFun This is a function that takes as input an events dataframe and a measurment function, and will return a quantile range of network values (see permutation vignette).
#' @param probs When nperm > 0 this will determine the probability of the permutation values returned from the permuations.
#' @param check.convergence If this is TRUE the function will calculate network measures for each window using random subsets of the data to measure the stability of the network measure. The value returned is the slope from the random subsets of decreasing size.
#' @param random.sample.size If check.convergence is TRUE this specifices the minimum size of the random subset used in calculating the convergence slope, i.e., minimum random subset size = actual sample size within a window - random.sample.size.
#' @param trim Whether nodes that are in windows beyond their first/last observation time are removed (i.e., only partially within a time window).
#' @param SRI Whether to convert edges to the simple ratio index: Nab / (Nab + Na + Nb). Default is set to FALSE.
#' @export
#' @importFrom lubridate days
#' @examples
#'
#' ts.out<-graphTS(data=groomEvents[1:200,])
#'
#'
graphTS <- function (data, windowsize = days(30), windowshift= days(1), measureFun=degree_mean,effortFun=NULL, permutationFun=perm.events,directed=FALSE, lagged=FALSE, lag=1, firstNet=FALSE, cores=1, nperm=0, probs=0.95, SRI=FALSE){

  #extract networks from the dataframe
  if(cores > 1){
    graphlist <- extract_networks_para(data, windowsize, windowshift, directed, cores = 2, SRI=SRI, effortFun = effortFun)
  } else {
    graphlist <- extract_networks(data, windowsize, windowshift, directed, SRI=SRI, effortFun = effortFun)
  }

  #extract measures from the network list
  if(lagged==FALSE){
    values <- extract_measure_network(graphlist, measureFun)
  } else {
    values <- extract_lagged_measure_network(graphlist, measureFun, lag, firstNet)
  }

  if(nperm>0){
    perm.values <- permutation.graph.values(data, windowsize, windowshift, directed, measureFun = measureFun, probs=probs, SRI=SRI)
    values <- cbind(data.frame(values),perm.values)
    values<-values[,c(1,5,6,2,3,4)]
  }

  return (values)

}


#' Convergence check
#'
#' This function will estimate the convergence of the chosen network measure using random subsets of the data.
#' @param data Dataframe with relational data in the first two rows, with weights in the thrid row, and a time stamp in the fourth row. Note: time stamps should be in ymd or ymd_hms format. The lubridate package can be very helpful in organizing times.
#' @param windowsize The size of each window in which to generate a network.
#' @param windowshift The amount of time to shift the window when generating networks.
#' @param directed Whether to consider the network as directed or not (TRUE/FALSE).
#' @importFrom stats coef
#' @importFrom stats lm
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

#' Bootstrap convergence check
#'
#' This function will estimate the convergence of the chosen network measure using bootstrapped samples of the data.
#' @param data Dataframe with relational data in the first two rows, with weights in the thrid row, and a time stamp in the fourth row. Note: time stamps should be in ymd or ymd_hms format. The lubridate package can be very helpful in organizing times.
#' @param windowsize The size of each window in which to generate a network.
#' @param windowshift The amount of time to shift the window when generating networks.
#' @param directed Whether to consider the network as directed or not (TRUE/FALSE).
#' @param boot.samples The number of bootstrapped samples to run (Default=100)
#' @importFrom stats cor.test
#' @importFrom igraph set_graph_attr
#' @export
#'
#'
convergence.check.boot<-function(data, windowsize, windowshift, directed = FALSE, measureFun=out_degree,boot.samples=100, SRI=FALSE, probs=c(0.025,0.975)){

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
      #cor.measures[length(cor.measures)+1] <- cor.test(comb.by.names[,1],comb.by.names[,2])$estimate
      cor.measures[length(cor.measures)+1] <- lsa::cosine(comb.by.names[,1],comb.by.names[,2])
    }

    #calculate convergence
    conv.values <- rbind(conv.values,data.frame(mean=mean(cor.measures),CI.low=quantile(cor.measures,probs = probs[1],na.rm = T),CI.high=quantile(cor.measures,probs = probs[2],na.rm = T),windowstart=windowstart, windowend=windowend))

    #move the window
    windowend = windowend + windowshift
    windowstart = windowstart + windowshift

  }

  conv.values<-conv.values[-1,]

  return(conv.values)

}


#' Extract networks from a moving window
#'
#' This function will create a time series of networks from a dataframe with relational events and a time stamp.
#' @param data Dataframe with relational data in the first two rows, with weights in the thrid row, and a time stamp in the fourth row. Note: time stamps should be in ymd or ymd_hms format. The lubridate package can be very helpful in organizing times.
#' @param windowsize The size of each window in which to generate a network.
#' @param windowshift The amount of time to shift the window when generating networks.
#' @param directed Whether to consider the network as directed or not (TRUE/FALSE).
#' @param trim Whether to remove nodes from the network if they are past the last observation time.
#' @importFrom igraph set_graph_attr
#' @export
#'
#'
extract_networks<-function(data, windowsize, windowshift, directed = FALSE,trim=FALSE, SRI=FALSE, effortFun=NULL){

  #intialize times
  windowstart <- min(data[,3])
  windowend=windowstart+windowsize
  if(windowend>max(data[,3]))print("warnning: the window size is set larger than the observed data.")

  #for every window generate a network
  netlist <- list()
  while (windowstart + windowsize<=max(data[,3])) {

    #subset the data
    df.window<-create.window(data, windowstart, windowend)
    if(trim==TRUE)df.window<-trim_graph(df.window,data,windowstart, windowend)
    Observation.Events <- nrow(df.window)

    #calculate effort
    if(is.null(effortFun)==FALSE){
      effort = effortFun(df.window)
    } else {
      effort = 1
    }

    #create a network and add it to the list
    g <- create.a.network(df.window, directed = directed, SRI, effort=effort)
    g <- set_graph_attr(g, "nEvents", Observation.Events)
    g <- set_graph_attr(g, "windowstart", windowstart )
    g <- set_graph_attr(g, "windowend", windowend)
    netlist[[length(netlist)+1]] <- g

    #move the window
    windowend = windowend + windowshift
    windowstart = windowstart + windowshift

  }

  print(paste0(length(netlist)," networks extracted"))
  return(netlist)

}

#' Extract networks from a moving window using multiple cores
#'
#' This function will create a time series of networks from a dataframe with relational events and a time stamp, using parallel processing.
#' @param data Dataframe with relational data in the first two rows, with weights in the thrid row, and a time stamp in the fourth row. Note: time stamps should be in ymd or ymd_hms format. The lubridate package can be very helpful in organizing times.
#' @param windowsize The size of each window in which to generate a network.
#' @param windowshift The amount of time to shift the window when generating networks.
#' @param directed Whether to consider the network as directed or not (TRUE/FALSE).
#' @param cores How many cores should be used.
#' @importFrom parallel makeCluster
#' @importFrom igraph set_graph_attr
#' @export
#'
#'
extract_networks_para<-function(data, windowsize, windowshift, directed = FALSE, cores=2,trim=FALSE, SRI, effortFun=NULL){

  #SRI not implimented yet
  if(SRI==TRUE)print("Warning SRI not yet available for parallel extraction of networks. Using SRI == FALSE.")

  #intialize times
  windowStart <- min(data[,3])
  windowEnd=windowStart+windowsize
  if(windowEnd>max(data[,3]))print("warnning: the window size is set larger than the observed data.")

  #generate a list of windows times
  window.ranges <- data.frame(start=windowStart, end=windowEnd)
  endDay=max(data[,3])
  while(windowEnd<=endDay){
    window.ranges <-  rbind(window.ranges,data.frame(start=windowStart, end=windowStart+windowsize))
    windowStart = windowStart + windowshift
    windowEnd = windowStart + windowsize
  }
  window.ranges<-window.ranges[-1,]

  #setup parallel backend
  cl <- parallel::makeCluster(cores)
  registerDoParallel(cl)

  #generate the networks
  final.net.list<-net.para(data, window.ranges, directed, trim = trim, effortFun=effortFun)

  #stop cluster
  parallel::stopCluster(cl)

  #report number of networks extracted
  print(paste0(length(final.net.list)," networks extracted"))

  return(final.net.list)
}

#' Extract networks in parallel using a dataframe of times
#'
#' This function will generate networks in parallel using a dataframe with time constraints.
#' @param data Dataframe with relational data in the first two rows, with weights in the thrid row, and a time stamp in the fourth row. Note: time stamps should be in ymd or ymd_hms format. The lubridate package can be very helpful in organizing times.
#' @param window.ranges The dataframe containing the start and end times of each window to create a network from.
#' @param directed Whether to consider the networks are directed or not.
#' @export
#'
#'
net.para<-function(data, window.ranges,directed=FALSE,trim, effortFun=NULL){

  #run the processes
  try(finalMatrix <- foreach(i=1:nrow(window.ranges), .export=c("effortFun","window.net","create.window", "create.a.network","window.net.para"), .packages = c("igraph", "dplyr") ) %dopar%

        net.window.para(data,windowstart = window.ranges[i,1], windowend = window.ranges[i,2], directed, trim=trim, effortFun=effortFun)

  )

  return(finalMatrix)
}


#' Extract one network within time constriants
#'
#' This function will generate one network from a dataframe with time constraints.
#' @param data Dataframe with relational data in the first two rows, with weights in the thrid row, and a time stamp in the fourth row. Note: time stamps should be in ymd or ymd_hms format. The lubridate package can be very helpful in organizing times.
#' @param windowstart The start of the window.
#' @param windowend The end of the window.
#' @param directed Whether to consider the network as directed or not.
#' @importFrom igraph set_graph_attr
#' @export
#'
#'
net.window<-function(data, windowstart, windowend,directed=FALSE, SRI){

  #subset the data
  df.window<-create.window(data, windowstart, windowend)
  Observation.Events <- nrow(df.window)

  #create a network and add it to the list
  g <- create.a.network(df.window, directed = directed, SRI)
  g <- set_graph_attr(g, "nEvents", Observation.Events)
  g <- set_graph_attr(g, "windowstart", windowstart)
  g <- set_graph_attr(g, "windowend", windowend)

  return(g)

}



#' Extract one network within time constriants
#'
#' This function will generate one network from a dataframe with time constraints.
#' @param data Dataframe with relational data in the first two rows, with weights in the thrid row, and a time stamp in the fourth row. Note: time stamps should be in ymd or ymd_hms format. The lubridate package can be very helpful in organizing times.
#' @param windowstart The start of the window.
#' @param windowend The end of the window.
#' @param directed Whether to consider the network as weighted. (default=FALSE)
#' @importFrom igraph set_graph_attr
#' @export
#'
#'
net.window.para<-function(data, windowstart, windowend,directed=FALSE, trim=FALSE, effortFun=NULL){

  #subset the data
  df.window <- data[data[[3]] >= windowstart & data[[3]] < windowend,]
  if(trim==TRUE)df.window<-trim_graph(df.window,data,windowstart, windowend)
  Observation.Events <- nrow(df.window)

  #calculate effort
  if(is.null(effortFun)==FALSE){
    effort = effortFun(df.window)
  } else {
    effort = 1
  }

  #create a network and add it to the list
  names(data)<-c("to","from","weight","date")
  elist<-data %>% dplyr::group_by(.dots=c("to","from")) %>% summarise(sum(weight)/effort)
  g <- graph_from_data_frame(elist, directed = directed, vertices = NULL)
  if(is.simple(g)==FALSE)g<-simplify(g, edge.attr.comb=list(weight="sum"))

  #add attributes
  g <- igraph::set_graph_attr(g, "nEvents", Observation.Events)
  g <- igraph::set_graph_attr(g, "windowstart", windowstart)
  g <- igraph::set_graph_attr(g, "windowend", windowend)

  return(g)

}


#' Extract network measures from a list of networks
#'
#' This function will estimate network measures from a list of networks.
#' @param netlist List of networks.
#' @param measureFun A function that takes a network as input and returns a single value.
#' @export
#' @importFrom igraph get.graph.attribute
#' @importFrom lubridate ymd
#'
#'
extract_measure_network<-function(netlist, measureFun, inOut=NULL){

  #store measures
  net.measure <- data.frame(measure=-1,nEvents=-1,windowstart=ymd("2000-01-01"), windowend=ymd("2000-01-01"))

  #extract measures
  if(exists('measureFun', mode='function')){

    for(i in 1:length(netlist)) {

      #get this network values
      obs.measure<-measureFun(netlist[[i]])
      windowstart.temp=igraph::get.graph.attribute(netlist[[i]], "windowstart" )
      windowend.temp=igraph::get.graph.attribute(netlist[[i]], "windowend" )

      if(length(obs.measure)>1)print("Warning: the measurement function is returing multiple values. In the case of graphTS the function should only return one value")

      #check if weighted average should be taken
      if(is.null(inOut)==FALSE){
        obs.measure<-weightedAvg(obs.measure,inOut,windowstart.temp,windowend.temp)
      }

      df.temp <- data.frame(measure=obs.measure,
                            nEvents=igraph::get.graph.attribute(netlist[[i]], "nEvents" ),
                            windowstart=windowstart.temp,
                            windowend=windowend.temp)
      net.measure <- rbind(net.measure,df.temp)
    }

  } else {
    print("Error: the measurment function was not found.")
  }

  net.measure<-net.measure[-1,]
  return(net.measure)

}


#' Calculate the weighted average of network measures
#'
#' This function will estimate the weighted average of node level measures based on their leaving and entering each network (e.g., due to births or deaths).
#' @param obs.values Vector of node level network measures, with column names for each individual measure.
#' @param inOut A dataframe with names in the first column, in a second column the dates these individuals first entered the network, and in a third column the dates these individuals left the network.
#' @param windowstart The starting date of the window.
#' @param windowend The ending date of the window.
#' @export
#' @importFrom igraph get.graph.attribute
#'
#'
weightedAvg<-function(obs.values, inOut, windowstart,windowend){

  IDs<-names(obs.values)
  window.interval <- lubridate::interval(windowstart,windowend)

  #get individual weights
  ID.weights<-vector()
  for (i in IDs){

    #get an individuals start and end time in the network
    ind.inOut<-inOut[inOut[,1]==IDs,]
    in.date <- ind.inOut[1,2]
    out.date <- ind.inOut[1,3]

    weight[i]<-1

    if ( (in.date %within% window.interval) & (out.date %within% window.interval) ){

      weight[i]<- as.numeric(as.duration(lubridate::interval(in.date,out.date)),"seconds") /lubridate::seconds(windowsize)

    } else if (in.date %within% window.interval){

      weight[i]<- as.numeric(as.duration(lubridate::interval(in.date,windowend)),"seconds")/lubridate::seconds(windowsize)

    } else if (out.date %within% window.interval){

      weight[i]<- as.numeric(as.duration(lubridate::interval(windowstart,out.date)),"seconds")/lubridate::seconds(windowsize)

    }

    if (weight[i]>1 | weight[i]< 0) print("something wrong with weights")

    ID.weights[length(ID.weights)] <- weight
  }

  #get new average network value
  total.weight <- sum(ID.weights)
  new.obs.values<- ID.weights*obs.values/total.weight

  return(new.obs.values)

}


#' Extract measures from a list of networks when the measure requires comparisons between networks
#'
#' This function will estimate network measures from a list of networks.
#' @param netlist List of networks.
#' @param measureFun A function that takes two networks as input and returns a single value. The first network is the lagged network, and the second is the current network.
#' @param lag At what lag should networks be compared? The number here will be based on the order of the network list generated. E.g., a list of networks generated using a window shift of 10 days, and a lag of 1, would compare networks 10days apart.
#' @param firstNet If TRUE the comparison between networks is always between the current and first network.
#' @export
#' @importFrom igraph get.graph.attribute
#'
#'
extract_lagged_measure_network<-function(netlist, measureFun, lag=1, firstNet=FALSE){

  #store measures
  net.measure <- data.frame(measure=-1,nEvents=-1,windowstart=ymd_hms("2000-01-01 12:00:00"), windowend=ymd_hms("2000-01-01 12:00:00"))

  if(exists('measureFun', mode='function')){



    for(i in 1:length(netlist)) {

      if(firstNet == FALSE){
        if(i-lag>1){

          df.temp <- data.frame(measure=measureFun(netlist[[i-lag]],netlist[[i]]),
                                nEvents=igraph::get.graph.attribute(netlist[[i]], "nEvents" ),
                                windowstart=igraph::get.graph.attribute(netlist[[i]], "windowstart" ),
                                windowend=igraph::get.graph.attribute(netlist[[i]], "windowend" ))

          net.measure <- rbind(net.measure,df.temp)

        } else {
          df.temp <- data.frame(measure=NA,
                                nEvents=igraph::get.graph.attribute(netlist[[i]], "nEvents" ),
                                windowstart=igraph::get.graph.attribute(netlist[[i]], "windowstart" ),
                                windowend=igraph::get.graph.attribute(netlist[[i]], "windowend" ))

          net.measure <- rbind(net.measure,df.temp)
        }

      } else {

        df.temp <- data.frame(measure=measureFun(netlist[[1]],netlist[[i]]),
                              nEvents=igraph::get.graph.attribute(netlist[[i]], "nEvents" ),
                              windowstart=igraph::get.graph.attribute(netlist[[i]], "windowstart" ),
                              windowend=igraph::get.graph.attribute(netlist[[i]], "windowend" ))

        net.measure <- rbind(net.measure,df.temp)

      }
    }


  } else {
    print("Error: the measurment function was not found.")
  }

  net.measure<-net.measure[-1,]
  return(net.measure)

}



#' Use permutation to extract uncertainty
#'
#' This function will estimate network measures given random permutations on the original data.
#' @param data Dataframe with relational data in the first two rows, with weights in the thrid row, and a time stamp in the fourth row. Note: time stamps should be in ymd or ymd_hms format. The lubridate package can be very helpful in organizing times.
#' @param windowsize The size of each window in which to generate a network.
#' @param windowshift The amount of time to shift the window when generating networks.
#' @param directed Whether to consider the network as directed or not (TRUE/FALSE).
#' @param measureFun This is a function that takes as an input a igraph network and returns a single value.
#' @param probs numeric vector of probabilities with values in [0,1].
#' @param nperm Number of permutations to perform before extracting network measures.
#' @export
#'
#'
permutation.graph.values<-function(data, windowsize, windowshift, directed = FALSE,measureFun, probs=0.95, nperm=1000, SRI=FALSE){

  #intialize times
  windowstart <- min(data[,3])
  windowend=windowstart+windowsize
  if(windowend>max(data[,3]))print("warnning: the window size is set larger than the observed data.")

  #monitor the progress
  pb <- txtProgressBar(min = as.numeric(windowstart + windowsize), max = as.numeric(max(data[,3])), style = 3)

  #for every window generate a network
  perm.values.high <- vector()
  perm.values.low <- vector()
  while (windowstart + windowsize<=max(data[,3])) {

    #subset the data
    df.window<-create.window(data, windowstart, windowend)
    Observation.Events <- nrow(df.window)

    #perform permutations
    perm.out<-perm.events(df.window, measureFun, directed, probs=probs,nperm= nperm, SRI=SRI)

    #record the high and low estimates
    perm.values.high[[length(perm.values.high)+1]] <- perm.out[2]
    perm.values.low[[length(perm.values.low)+1]] <- perm.out[1]

    #move the window
    windowend = windowend + windowshift
    windowstart = windowstart + windowshift

    #update progress bar
    setTxtProgressBar(pb,  as.numeric(windowend) )

  }

  perm.df<-data.frame(CI.low=perm.values.low,CI.high=perm.values.high)

  return(perm.df)

}





#' Trim nodes when taking network measures.
#'
#' This function removes node from the network when they are beyond their min and max observed times, then takes the network measure.
#' @param nodevalues Output from the nodeTS function.
#' @param data The events dataframe used in the nodeTS function.
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @export
#'
trim_graph<-function(df.window, data, windowstart, windowend){

  #ensure the names of the first four columns
  names(data)[1:3]<- c("from","to","date")

  #which names to keep
  names.kept<-unique( c(df.window[,1],df.window[,2]) )

  #Initialize trimed dataframes with important vars
  #df.trim<- data.frame(remove=rep(NA,nrow(nodevalues)))

  #loop through each ID and trim based on min and max date observed
  for(i in 1:length(names.kept)){

    #determine the min and max dates the focal was seen
    df.temp <- data %>% filter(from == names.kept[i] | to == names.kept[i])
    min.date<-min(df.temp$date)
    max.date<-max(df.temp$date)

    #check if individual was seen before/after this particular window
    if( (min.date<=windowstart & max.date>=windowend) == FALSE){

      #remove this individual from the window
      df.window<-df.window[df.window[,1]!=names.kept[i],]
      df.window<-df.window[df.window[,2]!=names.kept[i],]
    }

  }

  return(df.window)

}


