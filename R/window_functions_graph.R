

#' graphTS function
#'
#' This function will take a dataframe with events between individuals/objects, and take network measures using a moving window approach.
#' @param data A dataframe with relational data in the first two rows, with weights in the thrid row, and a time stamp in the fourth row. Note: time stamps should be in ymd or ymd_hms format. The lubridate package can be very helpful in organizing times.
#' @param windowsize The size of the moving window in which to take network measures. These should be provided as e.g., days(30), hours(5), ... etc.
#' @param windowshift The amount to shift the moving window for each measure. Again times should be provided as e.g., days(1), hours(1), ... etc.
#' @param measureFun This is a function that takes as an input a igraph network and returns a single value. There are functions within netTS (see details), and custom made functions can be used.
#' @param effortFun This is a function that takes as input the data within a window of time and returns the total sampling effort.
#' @param directed Whether the events are directed or no: true or false.
#' @param lagged Whether the network measure function used requires the comparison between two networks. e.g., comparing the current network to one lagged by 10 days. If TRUE the measureFun should take two graphs as input and return a single value. The order of inputs in the function is the lagged network followed by the current network.
#' @param lag If lagged is set to TRUE, this is the lag at which to compare networks.
#' @param firstNet If lagged is set to TRUE, this forces the comparisons between graphs to always be between the current and first graph.
#' @param cores This allows for multiple cores to be used while generating networks and calculating network measures.
#' @param nperm This allows for the estimation the network measure assuming random permutations. Currently the 95 percent quantiles are returned.
#' @param permutationFun This is a function that takes as input an events dataframe and a measurment function, and will return a quantile range of network values (see permutation vignette).
#' @param probs When nperm > 0 this will determine the probability of the permutation values returned from the permuations.
#' @param SRI Whether to convert edges to the simple ratio index: Nab / (Nab + Na + Nb). Default is set to FALSE.
#' @export
#' @importFrom lubridate days
#' @examples
#'
#' ts.out<-graphTS(data=groomEvents)
#'
#'
graphTS <- function (data, windowsize = days(30), windowshift= days(1), measureFun=degree_mean ,effortFun=NULL, permutationFun=perm.events,directed=FALSE, lagged=FALSE, lag=1, firstNet=FALSE, cores=1, nperm=0, probs=0.95, SRI=FALSE){

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
    perm.values <- permutation.graph.values(data, windowsize, windowshift, directed, measureFun = measureFun, probs=probs, SRI=SRI, graphlist = graphlist,permutationFun=permutationFun, nperm=nperm )
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

#' Bootstrap convergence check
#'
#' This function will estimate the convergence of the chosen network measure using bootstrapped samples of the data.
#' @param data Dataframe with relational data in the first two rows, with weights in the thrid row, and a time stamp in the fourth row. Note: time stamps should be in ymd or ymd_hms format. The lubridate package can be very helpful in organizing times.
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
convergence.check.boot <- function(data, windowsize, windowshift, directed = FALSE, measureFun=degree,boot.samples=100, SRI=FALSE, probs=c(0.025,0.975)){

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
      #cor.measures[length(cor.measures)+1] <- cor.test(comb.by.names[,1],comb.by.names[,2])$estimate
      cor.measures[length(cor.measures)+1] <- lsa::cosine(comb.by.names[,1],comb.by.names[,2])
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


#' Variance by window size check
#'
#' This function will estimate the variance in the time series based on the window size.
#' @param data Dataframe with relational data in the first two rows, with weights in the thrid row, and a time stamp in the fourth row. Note: time stamps should be in ymd or ymd_hms format. The lubridate package can be very helpful in organizing times.
#' @param windowsize_min The min size of windowsize to test.
#' @param windowsize_max The max size of windowsize to test.
#' @param by The resolution at which to test window sizes between the min and the max window sizes
#' @param windowshift The amount of time to shift the window when generating networks.
#' @param directed Whether to consider the network as directed or not (TRUE/FALSE).
#' @param measureFun The measurment function to perform the bootstap on (Default: density).
#' @param SRI Wether to use the simple ratio index (Default=FALSE)
#' @importFrom stats cor.test quantile
#' @importFrom igraph set_graph_attr degree density
#' @export
#'
#'
convergence.check.var<-function(data, windowsize_min=days(10),windowsize_max=days(40),by=days(1), windowshift=days(1), directed = FALSE, measureFun=igraph::edge_density, SRI=FALSE){

  #setup dataframe to return variance values
  df.var <- data.frame(windowsize = -1, var=-1)

  #for each window size calculate the overall variance and add it to df.var
  #pb <- txtProgressBar(0, length(x), style = 3)
  windowsize_seq = windowsize_min
  while(windowsize_seq<=windowsize_max){

    #calculate time series
    graph.values<-graphTS(data, windowsize = windowsize_seq, windowshift= windowshift, measureFun=measureFun ,effortFun=NULL, permutationFun=perm.events,directed=directed, lagged=FALSE, lag=1, firstNet=FALSE, cores=1, nperm=0, probs=0.95, SRI=SRI)

    #record
    df.var<-rbind(df.var, data.frame(windowsize=as.numeric(as.duration(windowsize_seq),"days"), var=var(graph.values[,1])) )

    #update window size tested
    windowsize_seq = windowsize_seq + by

    print(windowsize_seq)

  }

  df.var<-df.var[-1,]

  return(df.var)

}




#' Extract networks from a moving window
#'
#' This function will create a time series of networks from a dataframe with relational events and a time stamp.
#' @param data Dataframe with relational data in the first two rows, with weights in the thrid row, and a time stamp in the fourth row. Note: time stamps should be in ymd or ymd_hms format. The lubridate package can be very helpful in organizing times.
#' @param windowsize The size of each window in which to generate a network.
#' @param windowshift The amount of time to shift the window when generating networks.
#' @param directed Whether to consider the network as directed or not (TRUE/FALSE).
#' @param SRI Whether to use the simple ratio index (Default=FALSE).
#' @param effortFun This is a function that takes as input the data within a window of time and returns the total sampling effort.
#' @importFrom igraph set_graph_attr
#' @export
#'
#'
extract_networks<-function(data, windowsize, windowshift, directed = FALSE, SRI=FALSE, effortFun=NULL){

  #intialize times
  windowstart <- min(data[,3])
  windowend=windowstart+windowsize
  if(windowend>max(data[,3]))print("warnning: the window size is set larger than the observed data.")

  #for every window generate a network
  netlist <- list()
  while (windowstart + windowsize<=max(data[,3])) {

    #subset the data
    df.window<-create.window(data, windowstart, windowend)
    Observation.Events <- nrow(df.window)

    #calculate effort
    if(is.data.frame(effortFun)){
      effort = sum(effortFun[(effortFun[,1]>=windowstart & effortFun[,1]<windowend), ][,2])
    }else if(is.null(effortFun)==FALSE){
      effort = effortFun(df.window)
    } else {
      effort = 1
    }

    #create a network and add it to the list
    if(Observation.Events>0){
      g <- create.a.network(df.window, directed = directed, SRI, effort=effort)
      g <- set_graph_attr(g, "nEvents", Observation.Events)
      g <- set_graph_attr(g, "windowstart", windowstart )
      g <- set_graph_attr(g, "windowend", windowend)
      g <- set_graph_attr(g, "effort", effort)
      netlist[[length(netlist)+1]] <- g
    } else {
      g <- make_empty_graph(n=0, directed = directed)
      g <- set_graph_attr(g, "nEvents", Observation.Events)
      g <- set_graph_attr(g, "windowstart", windowstart )
      g <- set_graph_attr(g, "windowend", windowend)
      g <- set_graph_attr(g, "effort", effort)
      netlist[[length(netlist)+1]] <- g
    }


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
#' @param SRI Wether to use the simple ratio index (Default=FALSE).
#' @param effortFun This is a function that takes as input the data within a window of time and returns the total sampling effort.
#' @importFrom parallel makeCluster
#' @importFrom igraph set_graph_attr
#' @export
#'
#'
extract_networks_para<-function(data, windowsize, windowshift, directed = FALSE, cores=2, SRI, effortFun=NULL){

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
  final.net.list<-net.para(data, window.ranges, directed, effortFun=effortFun)

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
#' @param effortFun This is a function that takes as input the data within a window of time and returns the total sampling effort.
#' @importFrom foreach foreach %dopar%
#' @export
#'
#'
net.para<-function(data, window.ranges,directed=FALSE, effortFun=NULL){

  #run the processes
  try(finalMatrix <- foreach(i=1:nrow(window.ranges), .export=c("effortFun","window.net","create.window", "create.a.network","window.net.para"), .packages = c("igraph", "dplyr") ) %dopar%

        net.window.para(data,windowstart = window.ranges[i,1], windowend = window.ranges[i,2], directed, effortFun=effortFun)

  )

  return(finalMatrix)
}


#' Extract one network within time constriants
#'
#' This function will generate one network from a dataframe with time constraints.
#' @param data Dataframe with relational data in the first two rows, with weights in the thrid row, and a time stamp in the fourth row. Note: time stamps should be in ymd or ymd_hms format. The lubridate package can be very helpful in organizing times.
#' @param windowstart The start of the window.
#' @param windowend The end of the window.
#' @param directed Whether to consider the network as weighted. (default=FALSE)
#' @param effortFun This is a function that takes as input the data within a window of time and returns the total sampling effort.
#' @importFrom igraph set_graph_attr
#' @export
#'
#'
net.window.para<-function(data, windowstart, windowend,directed=FALSE, effortFun=NULL){

  #subset the data
  df.window <- data[data[[3]] >= windowstart & data[[3]] < windowend,]
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
  g <- igraph::set_graph_attr(g, "effort", effort)

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
extract_measure_network<-function(netlist, measureFun){

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


#' Calculate the weighted average correcting for nodes enter or leaving the network
#'
#' This function will estimate the weights of node based on their first and last observation and the time within a given window (e.g., due to births or deaths), and calculate the weighted average network value.
#' @param obs.values Vector of node level network measures, with column names for each individual measure.
#' @param inOut A dataframe with names in the first column, in a second column the dates these individuals first entered the network, and in a third column the dates these individuals left the network.
#' @param net The network used to calculate the weighted mean.
#' @export
#' @importFrom igraph get.graph.attribute
#' @importFrom lubridate %within% interval seconds as.duration
#'
#'
weighted_mean<-function(obs.values, inOut, net){

  #Individual nodes
  IDs<-names(obs.values)

  #start and end dates of this network
  windowstart = get.graph.attribute(net, "windowstart" )
  windowend   = get.graph.attribute(net, "windowend" )
  window.interval <- lubridate::interval(windowstart,windowend)
  windowsize = as.numeric(lubridate::seconds(window.interval))

  #get individual weights
  ID.weights<-vector()
  for (i in IDs){

    #get an individuals start and end time in the network
    ind.inOut<-inOut[inOut[,1]==i,]
    in.date <- ind.inOut[1,2]
    out.date <- ind.inOut[1,3]

    weight<-1

    if ( (in.date %within% window.interval) & (out.date %within% window.interval) ){

      weight<- as.numeric(as.duration(lubridate::interval(in.date,out.date)),"seconds") /windowsize

    } else if (in.date %within% window.interval){

      weight<- as.numeric(as.duration(lubridate::interval(in.date,windowend)),"seconds")/windowsize

    } else if (out.date %within% window.interval){

      weight<- as.numeric(as.duration(lubridate::interval(windowstart,out.date)),"seconds")/windowsize

    }

    if (weight>1 | weight< 0) print("something wrong with weights")

    ID.weights[length(ID.weights)+1] <- weight
  }

  #get new average network value
  total.weight <- sum(ID.weights)
  new.obs.values<- ID.weights*obs.values

  return(mean(new.obs.values))

}

#' First and last observation time of each node
#'
#' This calculates nodes min and max observed times.
#' @param data The events dataframe used in the nodeTS function.
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @export
#'
node_first_last<-function(data){

  #initialize inOut dataframe
  inOut <- data.frame(ID=data[1,1],min=data[1,3], max=data[1,3] )
  inOut<-inOut[-1,]

  #ensure the names of the first four columns
  names(data)[1:3]<- c("from","to","date")

  #which names to keep
  names.kept<-unique( c(as.character(data[,1]),as.character(data[,2]) ) )

  #loop through each ID and get min and max date observed
  for(i in 1:length(names.kept)){

    #determine the min and max dates the focal was seen
    df.temp <- data %>% dplyr::filter(from == names.kept[i] | to == names.kept[i])
    min.date<-min(df.temp$date)
    max.date<-max(df.temp$date)
    inOut <- rbind(inOut,data.frame(ID=names.kept[i],min=min.date,max=max.date))

  }

  return(inOut)
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
#' @importFrom lubridate ymd_hms
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
#' @param SRI Wether to use the simple ratio index (Default=FALSE).
#' @param graphlist A list of networks to perfrom permutations on.
#' @param permutationFun A function that will be used to perform permutations on the event data (i.e., before the network) or on the network itself. See vignette: Using_permutations.
#' @importFrom R.utils doCall
#' @export
#'
#'
permutation.graph.values<-function(data, windowsize, windowshift, directed = FALSE,measureFun, probs=0.95, nperm=1000, SRI=FALSE,graphlist=NULL, permutationFun=perm.events){

  print("perm")

  #vectors to record permutation results
  perm.values.high <- vector()
  perm.values.low <- vector()

  #monitor the progress
  #pb <- txtProgressBar(min = 1, max = length(graphlist), style = 3)

  #run the permutation on each network
  for(i in 1:length(graphlist)){

    #get the time bounds of each network
    windowstart <- graph_attr(graphlist[[i]],"windowstart")
    windowend   <- graph_attr(graphlist[[i]],"windowend")

    #get the data within this time range
    df.window<-create.window(data, windowstart, windowend)
    Observation.Events <- nrow(df.window)

    #perform permutations
    perm.out<-doCall(permutationFun,data=df.window, measureFun=measureFun, directed=directed, probs=probs,nperm= nperm, SRI=SRI, effort= graph_attr(graphlist[[i]],"effort") )

    #record the high and low estimates
    perm.values.high[[length(perm.values.high)+1]] <- perm.out[2]
    perm.values.low[[length(perm.values.low)+1]] <- perm.out[1]

    #update progress bar
    #setTxtProgressBar(pb,  as.numeric(windowend) )

  }

  #close progress bar
  #close(pb)

  #put the estimates into a data frame
  perm.df<-data.frame(CI.low=perm.values.low,CI.high=perm.values.high)

  return(perm.df)

}


