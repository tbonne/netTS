

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
#' @param cores This allows for multiple cores to be used while generating networks and calculating network measures.
#' @param nperm This allows for the estimation the network measure assuming random permutations. Currently the 95% quantiles are returned.
#' @param probs When nperm > 0 this will determine the probability of the permutation values returned from the permuations.
#' @export
#' @import lubridate
#' @examples
#'
#' ts.out<-graphTS(data=groomEvents[1:200,])
#'
graphTS <- function (data,windowsize =days(30), windowshift= days(1), measureFun=degree_mean,directed=FALSE, lagged=FALSE, lag=1, firstNet=FALSE, cores=1, nperm=0, probs=0.95){

  #extract networks from the dataframe
  if(cores > 1){
    graphlist <- extract_networks_para(data, windowsize, windowshift, directed, cores = 2)
  } else {
    graphlist <- extract_networks(data, windowsize, windowshift, directed)
  }

  #extract measures from the network list
  if(lagged==FALSE){
    values <- extract_measure_network(graphlist, measureFun)
  } else {
    values <- extract_lagged_measure_network(graphlist, measureFun, lag, firstNet)
  }

  if(nperm>0){
    perm.values <- permutation.graph.values(data, windowsize, windowshift, directed, measureFun = measureFun, probs=probs)
    values <- cbind(data.frame(values),perm.values)
    values<-values[,c(1,5,6,2,3,4)]
  }

  return (values)

}


#' Extract networks from a moving window
#'
#' This function will create a time series of networks from a dataframe with relational events and a time stamp.
#' @param data Dataframe with relational data in the first two rows, with weights in the thrid row, and a time stamp in the fourth row. Note: time stamps should be in ymd or ymd_hms format. The lubridate package can be very helpful in organizing times.
#' @param windowsize The size of each window in which to generate a network.
#' @param windowshift The amount of time to shift the window when generating networks.
#' @param directed Whether to consider the network as directed or not (TRUE/FALSE).
#' @importFrom igraph set_graph_attr
#' @export
#' @examples
#'
#'
#'
extract_networks<-function(data, windowsize, windowshift, directed = FALSE){

  #intialize times
  windowstart <- min(data[,4])
  windowend=windowstart+windowsize
  if(windowend>max(data[,4]))print("warnning: the window size is set larger than the observed data.")

  #for every window generate a network
  netlist <- list()
  while (windowstart + windowsize<=max(data[,4])) {

    #subset the data
    df.window<-create.window(data, windowstart, windowend)
    Observation.Events <- nrow(df.window)

    #create a network and add it to the list
    g <- create.a.network(df.window, directed = directed)
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
#' @import parallel
#' @import dplyr
#' @importFrom igraph set_graph_attr
#' @export
#' @examples
#'
#'
#'
extract_networks_para<-function(data, windowsize, windowshift, directed = FALSE, cores=2){

  #intialize times
  windowStart <- min(data[,4])
  windowEnd=windowStart+windowsize
  if(windowEnd>max(data[,4]))print("warnning: the window size is set larger than the observed data.")

  #generate a list of windows times
  window.ranges <- data.frame(start=windowStart, end=windowEnd)
  endDay=max(data[,4])
  while(windowEnd<=endDay){
    window.ranges <-  rbind(window.ranges,data.frame(start=windowStart, end=windowStart+windowsize))
    windowStart = windowStart + windowshift
    windowEnd = windowStart + windowsize
  }
  window.ranges<-window.ranges[-1,]

  #setup parallel backend
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)

  #generate the networks
  final.net.list<-net.para(data, window.ranges, directed)

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
net.para<-function(data, window.ranges,directed=FALSE){

  #run the processes
  try(finalMatrix <- foreach(i=1:nrow(window.ranges), .export=c("window.net","create.window", "create.a.network","window.net.para"), .packages = c("igraph", "dplyr") ) %dopar%

        net.window.para(data,windowstart = window.ranges[i,1], windowend = window.ranges[i,2], directed)

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
net.window<-function(data, windowstart, windowend,directed=FALSE){

  #subset the data
  df.window<-create.window(data, windowstart, windowend)
  Observation.Events <- nrow(df.window)

  #create a network and add it to the list
  g <- create.a.network(df.window, directed = directed)
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
net.window.para<-function(data, windowstart, windowend,directed=FALSE){

  #subset the data
  df.window <- data[data[[4]] >= windowstart & data[[4]] < windowend,]
  Observation.Events <- nrow(df.window)

  #create a network and add it to the list
  names(data)<-c("to","from","weight","date")
  elist<-data %>% dplyr::group_by(.dots=c("to","from")) %>% summarise(sum(weight))
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
#' @import igraph
#' @examples
#'
#'
#'
extract_measure_network<-function(netlist, measureFun){

  #store measures
  net.measure <- data.frame(measure=-1,nEvents=-1,windowstart=ymd("2000-01-01"), windowend=ymd("2000-01-01"))

  #extract measures
  if(exists('measureFun', mode='function')){

    for(i in 1:length(netlist)) {
      df.temp <- data.frame(measure=measureFun(netlist[[i]]),
                            nEvents=igraph::get.graph.attribute(netlist[[i]], "nEvents" ),
                            windowstart=igraph::get.graph.attribute(netlist[[i]], "windowstart" ),
                            windowend=igraph::get.graph.attribute(netlist[[i]], "windowend" ))
      net.measure <- rbind(net.measure,df.temp)
    }

  } else {
    print("Error: the measurment function was not found.")
  }

  net.measure<-net.measure[-1,]
  return(net.measure)

}




#' Extract measures from a list of networks when the measure requires comparisons between networks
#'
#' This function will estimate network measures from a list of networks.
#' @param netlist List of networks.
#' @param measureFun A function that takes two networks as input and returns a single value. The first network is the lagged network, and the second is the current network.
#' @param lag At what lag should networks be compared? The number here will be based on the order of the network list generated. E.g., a list of networks generated using a window shift of 10 days, and a lag of 1, would compare networks 10days apart.
#' @export
#' @import igraph
#' @examples
#'
#'
#'
extract_lagged_measure_network<-function(netlist, measureFun, lag=1, firstNet){

  #store measures
  net.measure <- data.frame(measure=-1,nEvents=-1,windowstart=ymd("2000-01-01"), windowend=ymd("2000-01-01"))

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
#' @param
#' @export
#' @examples
#'
#'
#'
permutation.graph.values<-function(data, windowsize, windowshift, directed = FALSE,measureFun, probs=0.95, nperm=1000){

  #intialize times
  windowstart <- min(data[,4])
  windowend=windowstart+windowsize
  if(windowend>max(data[,4]))print("warnning: the window size is set larger than the observed data.")

  #monitor the progress
  pb <- txtProgressBar(min = as.numeric(windowstart + windowsize), max = as.numeric(max(data[,4])), style = 3)

  #for every window generate a network
  perm.values.high <- vector()
  perm.values.low <- vector()
  while (windowstart + windowsize<=max(data[,4])) {

    #subset the data
    df.window<-create.window(data, windowstart, windowend)
    Observation.Events <- nrow(df.window)

    #perform permutations
    perm.out<-perm.interactions(df.window, measureFun, directed, probs=probs,nperm= nperm)

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




#' Perform permutation
#'
#' This function will permute a network by randomly switching individuals within events.
#' @param data Dataframe with relational data in the first two rows, with weights in the thrid row, and a time stamp in the fourth row. Note: time stamps should be in ymd or ymd_hms format. The lubridate package can be very helpful in organizing times.
#' @param measureFun This is a function that takes as an input a igraph network and returns a single value.
#' @param directed Whether to consider the network as directed or not (TRUE/FALSE).
#' @param nperm Number of permutations to perform before extracting network measures.
#' @param probs numeric vector of probabilities with values in [0,1].
#' @export
#' @examples
#'
#'
#'
perm.interactions <- function(data, measureFun, directed=FALSE, nperm=1000, probs=0.95){

  net.list <- list(create.a.network(data, directed))
  Perm.measure<-vector()

  for(i in 1:nperm){

    no.loops= FALSE

    #choose to or from grooming to permute
    if(0.5 > runif(1)){

      while(no.loops == FALSE){

        #choose two individuals to switch
        rows.to.switch <- sample(1:nrow(data),2,F)

        #record old order
        old.order <- data$to
        new.order <- data$to

        #update order
        new.order[rows.to.switch[1]] <- old.order[rows.to.switch[2]]
        new.order[rows.to.switch[2]] <- old.order[rows.to.switch[1]]

        #check to make sure there are no self loops
        if(sum(as.character(data$from)==as.character(new.order) )==0){

          data$to <- new.order
          NewData<- data
          no.loops=TRUE

        }
      }

    }  else {
      while(no.loops == FALSE){

        #choose two individuals to switch
        rows.to.switch <- sample(1:nrow(data),2,F)

        #record old order
        old.order <- data$from
        new.order <- data$from

        #update order
        new.order[rows.to.switch[1]] <- old.order[rows.to.switch[2]]
        new.order[rows.to.switch[2]] <- old.order[rows.to.switch[1]]

        #check to make sure there are no self loops
        if(sum(as.character(new.order)==as.character(data$to) )==0){

          data$from <- new.order
          NewData<- data
          no.loops=TRUE

        }
      }
    }

    #Create graph in order to get the measure
    Perm.network <- create.a.network(NewData, directed)

    # Get measure
    Perm.measure[length(Perm.measure)+1]<- measureFun(Perm.network)

  }

  probs.left<-1-probs
  return(quantile(Perm.measure, probs = c( (0+probs.left/2), (1-probs.left/2) ), na.rm=T))
}














