

#' nodeTS function
#'
#' This function will take a dataframe with events between individuals/objects, and take network measures using a moving window approach.
#' @param data A dataframe with relational data in the first two rows, with weights in the thrid row, and a time stamp in the fourth row. Note: time stamps should be in ymd_hms format. The lubridate package can be very helpful in organizing times. Note: names in the first two columns are case sensitive.
#' @param windowsize The size of the moving window in which to take network measures. These should be provided as e.g., days(30), hours(5), ... etc.
#' @param windowshift The amount to shift the moving window for each measure. Again times should be provided as e.g., days(1), hours(1), ... etc.
#' @param measureFun This is a function that takes as an input a igraph network and returns values for each node in the network. There are functions within netTS (see details), and custom made functions can be used.
#' @param directed Whether the events are directed or no: true or false.
#' @param lagged Whether the network measure function used requires the comparison between two networks. e.g., comparing the current network to one lagged by 10 days. If TRUE the measureFun should take two graphs as input and return a single value. The order of inputs in the function is the lagged network followed by the current network.
#' @param lag If lagged is set to TRUE, this is the lag at which to compare networks.
#' @param firstNet If lagged is set to TRUE, this compares the subsequent networks to the first network.
#' @param cores This allows for multiple cores to be used while generating networks and calculating network measures.
#' @export
#' @importFrom lubridate days
#' @examples
#'
#' ts.out<-nodeTS(data=groomEvents[1:200,])
#'
nodeTS <- function (data,windowsize =days(30), windowshift= days(1), measureFun=degree, effortFun=NULL,directed=FALSE, lagged=FALSE, lag=1, firstNet=FALSE, cores=1){

  #extract networks from the dataframe
  if(cores > 1){
    graphlist <- extract_networks_para(data, windowsize, windowshift, directed, cores = 2, effortFun = effortFun)
  } else {
    graphlist <- extract_networks(data, windowsize, windowshift, directed, effortFun = effortFun)
  }



  #extract measures from the network list
  all.unique.names <- unique(c(as.character(data[,1]),as.character(data[,2]) ))
  if(lagged==FALSE){
    values <- extract_measure_nodes(netlist=graphlist, measureFun, unique.names = all.unique.names)
  } else {
    values <- extract_lagged_measure_nodes(graphlist, measureFun, lag, unique.names = all.unique.names, firstNet=firstNet)
  }

  return (values)

}



#' Extract node level network measures from a list of networks
#'
#' This function will estimate node level network measures from a list of networks.
#' @param netlist List of networks.
#' @param measureFun A function that takes a network as input and returns a value for each node.
#' @param unique.names A list of all names/nodes in the networks.
#' @export
#' @importFrom igraph get.graph.attribute
#' @importFrom dplyr bind_rows
#'
#'
#'
extract_measure_nodes<-function(netlist, measureFun, unique.names){

  #store measures - set global dataframe with proper names
  netvalues <- data.frame(t(rep(-1,length(unique.names))))
  names(netvalues) <- unique.names
  net.measure <- data.frame(nEvents=-1,windowstart=ymd_hms("2000-01-01 12:00:00"), windowend=ymd_hms("2000-01-01 12:00:00"))
  netvalues<-cbind(netvalues,net.measure)

  #extract measures
  if(exists('measureFun', mode='function')){

    for(i in 1:length(netlist)) {
      df.temp.nodes <- data.frame(t(measureFun(netlist[[i]])))
      df.temp.graph <- data.frame(nEvents=igraph::get.graph.attribute(netlist[[i]], "nEvents" ),
                            windowstart=igraph::get.graph.attribute(netlist[[i]], "windowstart" ),
                            windowend=igraph::get.graph.attribute(netlist[[i]], "windowend" ))
      df.temp <- cbind(df.temp.nodes,df.temp.graph)
      netvalues <- dplyr::bind_rows(netvalues,df.temp)
    }

  } else {
    print("Error: the measurment function was not found.")
  }

  netvalues<-netvalues[-1,]
  return(netvalues)

}


#' Extract node level measures from a list of networks when the measure requires comparisons between networks.
#'
#' This function will estimate node level network measures from a list of networks.
#' @param netlist List of networks.
#' @param measureFun A function that takes two networks as input and returns a single value. The first network is the lagged network, and the second is the current network.
#' @param lag At what lag should networks be compared? The number here will be based on the order of the network list generated. E.g., a list of networks generated using a window shift of 10 days, and a lag of 1, would compare networks 10days apart.
#' @param unique.names A list of all names/nodes in the networks.
#' @export
#' @importFrom igraph get.graph.attribute
#' @importFrom dplyr bind_rows
#'
#'
#'
extract_lagged_measure_nodes<-function(netlist, measureFun, lag=1, unique.names, firstNet=FALSE){

  #store measures - set global dataframe with proper names
  netvalues <- data.frame(t(rep(-1,length(unique.names))))
  names(netvalues) <- unique.names
  net.measure <- data.frame(nEvents=-1,windowstart=ymd_hms("2000-01-01 12:00:00"), windowend=ymd_hms("2000-01-01 12:00:00"))
  netvalues<-cbind(netvalues,net.measure)

  if(exists('measureFun', mode='function')){

    if(firstNet == FALSE){
      for(i in 1:length(netlist)) {
        if(i-lag>1){

          df.temp.nodes <- data.frame((measureFun(netlist[[i-lag]],netlist[[i]])))
          df.temp.graph <- data.frame(nEvents=igraph::get.graph.attribute(netlist[[i]], "nEvents" ),
                                      windowstart=igraph::get.graph.attribute(netlist[[i]], "windowstart" ),
                                      windowend=igraph::get.graph.attribute(netlist[[i]], "windowend" ))
          df.temp <- cbind(df.temp.nodes,df.temp.graph)
          netvalues <- bind_rows(netvalues,df.temp)

        } else {
          df.temp.nodes <- data.frame(t(rep(NA,length(unique.names)) ))
          names(df.temp.nodes)<-unique.names
          df.temp.graph <- data.frame(nEvents=igraph::get.graph.attribute(netlist[[i]], "nEvents" ),
                                      windowstart=igraph::get.graph.attribute(netlist[[i]], "windowstart" ),
                                      windowend=igraph::get.graph.attribute(netlist[[i]], "windowend" ))
          df.temp <- cbind(df.temp.nodes,df.temp.graph)
          netvalues <- bind_rows(netvalues,df.temp)
        }

      }
    } else {

      for(i in 1:length(netlist)) {
        df.temp.nodes <- data.frame((measureFun(netlist[[1]],netlist[[i]])))
        df.temp.graph <- data.frame(nEvents=igraph::get.graph.attribute(netlist[[i]], "nEvents" ),
                                    windowstart=igraph::get.graph.attribute(netlist[[i]], "windowstart" ),
                                    windowend=igraph::get.graph.attribute(netlist[[i]], "windowend" ))
        df.temp <- cbind(df.temp.nodes,df.temp.graph)
        netvalues <- bind_rows(netvalues,df.temp)
      }

    }

  } else {
    print("Error: the measurment function was not found.")
  }

  netvalues<-netvalues[-1,]
  return(netvalues)

}


#' Trim nodes.
#'
#' This function removes node values beyond their min and max observed times.
#' @param nodevalues Output from the nodeTS function.
#' @param data The events dataframe used in the nodeTS function.
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @export
#'
trim_nodes<-function(nodevalues, data){

  #ensure the names of the first four columns
  names(data)[1:3]<- c("from","to","date")

  #which names to keep
  names.kept<-colnames(nodevalues)[1:(length(nodevalues)-3)]

  #Initialize trimed dataframes with important vars
  df.trim<- data.frame(remove=rep(NA,nrow(nodevalues)))

  #loop through each ID and trim based on min and max date observed
  for(i in 1:length(names.kept)){

    #determine the min and max dates the focal was seen
    df.temp <- data %>% filter(from == names.kept[i] | to == names.kept[i])
    min.date<-min(df.temp$date)
    max.date<-max(df.temp$date)

    #remove all window estimates outside those dates
    df.temp2 <- nodevalues %>% dplyr::select(names.kept[i], windowstart, windowend)
    df.temp2[,1] <- ifelse(df.temp2[,2]<min.date,NA,df.temp2[,1])
    df.temp2[,1] <- ifelse(df.temp2[,3]>max.date,NA,df.temp2[,1])

    #build new trimed dataframes
    df.temp3 <- data.frame((df.temp2[,1]))
    names(df.temp3) <- c(names.kept[i])
    df.trim <- cbind(df.trim, df.temp3)

  }

  df.trim<-df.trim[,-1]
  df.times <- nodevalues %>% dplyr::select("nEvents","windowstart","windowend")
  df.trim <- cbind(df.trim,df.times)


  return(df.trim)

}

