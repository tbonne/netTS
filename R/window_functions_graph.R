

#' graphTS function
#'
#' This function will take a dataframe with events between individuals/objects, and take network measures using a moving window approach.
#' @param data A dataframe with relational data in the first two columns, and a time stamp in the third column. An optional column with a weight can be added if there is a duration or magnitude for each interaction (column name for this should be set to 'weight'). Note: time stamps should be in ymd or ymd_hms format. The lubridate package can be very helpful in organizing times.
#' @param windowsize The size of the moving window in which to take network measures. These should be provided as e.g., days(30), hours(5), ... etc.
#' @param windowshift The amount to shift the moving window for each measure. Again times should be provided as e.g., days(1), hours(1), ... etc.
#' @param measureFun This is a function that takes as an input a igraph network and returns a single value. There are functions within netTS (see details), and custom made functions can be used.
#' @param effortFun This is a function that takes as input the data within a window of time and returns the total sampling effort.
#' @param effortData This is a dataframe containing the data used to calculate sampling effort.
#' @param directed Whether the events are directed or no: true or false.
#' @param lagged Whether the network measure function used requires the comparison between two networks. e.g., comparing the current network to one lagged by 10 days. If TRUE the measureFun should take two graphs as input and return a single value. The order of inputs in the function is the lagged network followed by the current network.
#' @param lag If lagged is set to TRUE, this is the lag at which to compare networks.
#' @param firstNet If lagged is set to TRUE, this forces the comparisons between graphs to always be between the current and first graph.
#' @param cores This allows for multiple cores to be used while generating networks and calculating network measures.
#' @param nperm This allows for the estimation the network measure assuming random permutations. Currently the 95 percent quantiles are returned.
#' @param permutationFun This is a function that takes as input an events dataframe and a measurment function, and will return a quantile range of network values (see permutation vignette).
#' @param probs When nperm > 0 this will determine the probability of the permutation values returned from the permuations.
#' @param SRI Whether to convert edges to the simple ratio index: Nab / (Nab + Na + Nb). Default is set to FALSE.
#' @param windowstart When the moving window should start. Default is the minimum observation time.
#' @param windowend When the moving window should stop. Default is the maximum observation time.
#' @export
#' @importFrom lubridate days
#' @examples
#'
#' ts.out<-graphTS(data=groomEvents)
#'
#'
graphTS <- function (data, windowsize = days(30), windowshift= days(1), measureFun=degree_mean ,effortFun=NULL,effortData=NULL, permutationFun=perm.events,directed=FALSE, lagged=FALSE, lag=1, firstNet=FALSE, cores=1, nperm=0, probs=0.95, SRI=FALSE, windowstart=NULL,windowend=NULL){

  #check for missing data
  if(sum(is.na(data)) > 0){
    print(paste0("Data contains NA, removing ",sum(!complete.cases(data)), " row(s)."))
    data<-data[complete.cases(data),]
  }

  #extract networks from the dataframe
  if(cores > 1){
    graphlist <- extract_networks_para(data, windowsize, windowshift, directed=directed, cores = cores, SRI=SRI, effortFun = effortFun, effortData=effortData, winstart=windowstart,winend=windowend)
  } else {
    graphlist <- extract_networks(data, windowsize, windowshift, directed=directed, SRI=SRI, effortFun = effortFun, effortData=effortData, winstart=windowstart,winend=windowend)
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


#' Extract networks from a moving window
#'
#' This function will create a time series of networks from a dataframe with relational events and a time stamp.
#' @param data A dataframe with relational data in the first two columns, and a time stamp in the third column. An optional column with a weight can be added if there is a duration or magnitude for each interaction (column name for this should be set to 'weight'). Note: time stamps should be in ymd or ymd_hms format. The lubridate package can be very helpful in organizing times.
#' @param windowsize The size of each window in which to generate a network.
#' @param windowshift The amount of time to shift the window when generating networks.
#' @param directed Whether to consider the network as directed or not (TRUE/FALSE).
#' @param SRI Whether to use the simple ratio index (Default=FALSE).
#' @param effortFun This is a function that takes as input the data within a window of time and returns the total sampling effort.
#' @param effortData This is a dataframe containing the data used to calculate sampling effort. The first column should contain timedate values.
#' @importFrom igraph set_graph_attr
#' @export
#'
#'
extract_networks<-function(data, windowsize, windowshift, directed = FALSE, SRI=FALSE, effortFun=NULL, effortData=NULL, winstart=NULL,winend=NULL){

  #intialize times
  if(is.null(winstart) ){
    windowstart <- min(data[,3])
  } else {
    windowstart <- winstart
  }

  if(is.null(winend) ){
    windowmax <- max(data[,3])
  } else {
    windowmax <- winend
  }

  windowend=windowstart+windowsize
  if(windowend>windowmax)print("warnning: the window size is set larger than the observed data.")

  #for every window generate a network
  netlist <- list()
  while (windowend<=windowmax) {

    #subset the data
    df.window<-create.window(data, windowstart, windowend)
    Observation.Events <- nrow(df.window)

    #calculate sampling effort
    if(is.null(effortFun)==FALSE & is.null(effortData)==TRUE  ){ #there is an effort function and it requires no external data
      g = effortFun(df.window, directed = directed)
      g <- set_graph_attr(g, "nEvents", Observation.Events)
      g <- set_graph_attr(g, "windowstart", windowstart )
      g <- set_graph_attr(g, "windowend", windowend)
      netlist[[length(netlist)+1]] <- g

    }else if(is.null(effortFun)==FALSE & is.null(effortData)==FALSE ){ #there is an effort function and it requires some external data
      effortData.sub <- effortData[effortData[,1]>=windowstart & effortData[,1]<windowend,]
      g = effortFun(df.window, effortData.sub, directed = directed)
      g <- set_graph_attr(g, "nEvents", Observation.Events)
      g <- set_graph_attr(g, "windowstart", windowstart )
      g <- set_graph_attr(g, "windowend", windowend)
      netlist[[length(netlist)+1]] <- g

    } else { #there is no effort function

      if(Observation.Events>0){
        g <- create.a.network(df.window, directed = directed, SRI)
        g <- set_graph_attr(g, "nEvents", Observation.Events)
        g <- set_graph_attr(g, "windowstart", windowstart )
        g <- set_graph_attr(g, "windowend", windowend)
        netlist[[length(netlist)+1]] <- g

      } else {
        g <- make_empty_graph(n=0, directed = directed)
        g <- set_graph_attr(g, "nEvents", Observation.Events)
        g <- set_graph_attr(g, "windowstart", windowstart )
        g <- set_graph_attr(g, "windowend", windowend)
        netlist[[length(netlist)+1]] <- g
      }
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
#' @param data A dataframe with relational data in the first two columns, and a time stamp in the third column. An optional column with a weight can be added if there is a duration or magnitude for each interaction (column name for this should be set to 'weight'). Note: time stamps should be in ymd or ymd_hms format. The lubridate package can be very helpful in organizing times.
#' @param windowsize The size of each window in which to generate a network.
#' @param windowshift The amount of time to shift the window when generating networks.
#' @param directed Whether to consider the network as directed or not (TRUE/FALSE).
#' @param cores How many cores should be used.
#' @param SRI Wether to use the simple ratio index (Default=FALSE).
#' @param effortFun This is a function that takes as input the data within a window of time and returns the total sampling effort.
#' @param effortData This is a dataframe containing the data used to calculate sampling effort. The first column should contain timedate values.
#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom igraph set_graph_attr
#' @export
#'
#'
extract_networks_para<-function(data, windowsize, windowshift, directed = FALSE, cores=2, SRI, effortFun=NULL, effortData=NULL, winstart=NULL,winend=NULL){

  #SRI not implimented yet
  if(SRI==TRUE)print("Warning SRI not yet available for parallel extraction of networks. Using SRI == FALSE.")

  #intialize times
  #intialize times
  if(is.null(winstart) ){
    windowStart <- min(data[,3])
  } else {
    windowStart <- winstart
  }

  if(is.null(winend) ){
    windowmax <- max(data[,3])
  } else {
    windowmax <- winend
  }


  windowEnd=windowStart+windowsize
  if(windowEnd>windowmax)print("warnning: the window size is set larger than the observed data.")

  #generate a list of windows times
  window.ranges <- data.frame(start=windowStart, end=windowEnd)
  while(windowEnd<=windowmax){
    window.ranges <-  rbind(window.ranges,data.frame(start=windowStart, end=windowStart+windowsize))
    windowStart = windowStart + windowshift
    windowEnd = windowStart + windowsize
  }
  window.ranges<-window.ranges[-1,]

  #setup parallel backend
  cl <- parallel::makeCluster(cores)
  registerDoParallel(cl)

  #generate the networks
  final.net.list<-net.para(data, window.ranges, directed, effortFun=effortFun, effortData=effortData)

  #stop cluster
  parallel::stopCluster(cl)

  #report number of networks extracted
  print(paste0(length(final.net.list)," networks extracted"))

  return(final.net.list)
}

#' Extract networks in parallel using a dataframe of times
#'
#' This function will generate networks in parallel using a dataframe with time constraints.
#' @param data A dataframe with relational data in the first two columns, and a time stamp in the third column. An optional column with a weight can be added if there is a duration or magnitude for each interaction (column name for this should be set to 'weight'). Note: time stamps should be in ymd or ymd_hms format. The lubridate package can be very helpful in organizing times.
#' @param window.ranges The dataframe containing the start and end times of each window to create a network from.
#' @param directed Whether to consider the networks are directed or not.
#' @param effortFun This is a function that takes as input the data within a window of time and returns the total sampling effort.
#' @importFrom foreach foreach %dopar%
#' @export
#'
#'
net.para<-function(data, window.ranges,directed=FALSE, effortFun=NULL, effortData=NULL){

  #finalMatrix <- NA

  #run the processes
  finalMatrix <- foreach(i=1:nrow(window.ranges), .export=c("create.window", "create.a.network","net.window.para"), .packages = c("igraph",  "dplyr") ) %dopar%

        net.window.para(data,windowstart = window.ranges[i,1], windowend = window.ranges[i,2], directed, effortFun=effortFun, effortData=effortData)


  return(finalMatrix)

}


#' Extract one network within time constriants
#'
#' This function will generate one network from a dataframe with time constraints.
#' @param data A dataframe with relational data in the first two columns, and a time stamp in the third column. An optional column with a weight can be added if there is a duration or magnitude for each interaction (column name for this should be set to 'weight'). Note: time stamps should be in ymd or ymd_hms format. The lubridate package can be very helpful in organizing times.
#' @param windowstart The start of the window.
#' @param windowend The end of the window.
#' @param directed Whether to consider the network as weighted. (default=FALSE)
#' @param effortFun This is a function that takes as input the data within a window of time and returns the total sampling effort.
#' @param effortData This is a dataframe containing the data used to calculate sampling effort. The first column should contain timedate values.
#' @importFrom igraph set_graph_attr
#' @export
#'
#'
net.window.para<-function(data, windowstart, windowend,directed=FALSE, effortFun=NULL, effortData=NULL){

  #subset the data
  df.window <- data[data[[3]] >= windowstart & data[[3]] < windowend,]
  Observation.Events <- nrow(df.window)

  #calculate sampling effort
  if(is.null(effortFun)==FALSE & is.null(effortData)==TRUE ){ #there is an effort function and it requires no external data
    g = effortFun(df.window, directed = directed)
    g <- igraph::set_graph_attr(g, "nEvents", Observation.Events)
    g <- igraph::set_graph_attr(g, "windowstart", windowstart)
    g <- igraph::set_graph_attr(g, "windowend", windowend)

  }else if(is.null(effortFun)==FALSE & is.null(effortData)==FALSE){ #there is an effort function and it requires some external data
    effortData.sub <- effortData[effortData[,1]>=windowstart & effortData[,1]<windowend,]
    g = effortFun(df.window, effortData.sub, directed = directed)
    g <- igraph::set_graph_attr(g, "nEvents", Observation.Events)
    g <- igraph::set_graph_attr(g, "windowstart", windowstart)
    g <- igraph::set_graph_attr(g, "windowend", windowend)

  } else { #there is no effort function

    #create a network and add it to the list
    names(data)[1:2]<-c("from","to")
    if(is.null(data$weight))data$weight=1
    #elist <- as.data.frame(as.data.table(data)[,.(sum(weight)/effort), by=list(from,to)])
    elist<-data %>% dplyr::group_by(.dots=c("to","from")) %>% summarise(sum(weight))
    names(elist)<-c("from","to","weight")
    g <- graph_from_data_frame(elist, directed = directed, vertices = NULL)
    if(is.simple(g)==FALSE)g<-simplify(g, edge.attr.comb=list(weight="sum"))

    #add attributes
    g <- igraph::set_graph_attr(g, "nEvents", Observation.Events)
    g <- igraph::set_graph_attr(g, "windowstart", windowstart)
    g <- igraph::set_graph_attr(g, "windowend", windowend)
  }

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
#' @param data A dataframe with relational data in the first two columns, and a time stamp in the third column. An optional column with a weight can be added if there is a duration or magnitude for each interaction (column name for this should be set to 'weight'). Note: time stamps should be in ymd or ymd_hms format. The lubridate package can be very helpful in organizing times.
#' @param windowsize The size of each window in which to generate a network.
#' @param windowshift The amount of time to shift the window when generating networks.
#' @param directed Whether to consider the network as directed or not (TRUE/FALSE).
#' @param measureFun This is a function that takes as an input a igraph network and returns a single value.
#' @param effortFun This is a function that takes as input the data within a window of time and returns the total sampling effort.
#' @param effortData This is a dataframe containing the data used to calculate sampling effort. The first column should contain timedate values.
#' @param probs numeric vector of probabilities with values in [0,1].
#' @param nperm Number of permutations to perform before extracting network measures.
#' @param SRI Wether to use the simple ratio index (Default=FALSE).
#' @param graphlist A list of networks to perfrom permutations on.
#' @param permutationFun A function that will be used to perform permutations on the event data (i.e., before the network) or on the network itself. See vignette: Using_permutations.
#' @importFrom R.utils doCall
#' @export
#'
#'
permutation.graph.values<-function(data, windowsize, windowshift, directed = FALSE,measureFun, effortFun=NULL, effortData=NULL, probs=0.95, nperm=1000, SRI=FALSE,graphlist=NULL, permutationFun=perm.events){

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
    perm.out<-doCall(permutationFun,data=df.window, measureFun=measureFun, directed=directed, windowstart=windowstart, windowend=windowend, probs=probs,nperm= nperm, SRI=SRI, effortFun=NULL, effortData=NULL)

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

