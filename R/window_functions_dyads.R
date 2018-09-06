


#' dyadTS function
#'
#' This function will take a dataframe with events between individuals/objects, and take network measures using a moving window approach.
#' @param data A dataframe with relational data in the first two rows, with weights in the thrid row, and a time stamp in the fourth row. Note: time stamps should be in ymd or ymd_hms format. The lubridate package can be very helpful in organizing times. Note: names in the first two columns are case sensitive.
#' @param windowsize The size of the moving window in which to take network measures. These should be provided as e.g., days(30), hours(5), ... etc.
#' @param windowshift The amount to shift the moving window for each measure. Again times should be provided as e.g., days(1), hours(1), ... etc.
#' @param measureFun This is a function that takes as an input a igraph network and returns values for each node in the network. There are functions within netTS (see details), and custom made functions can be used.
#' @param directed Whether the events are directed or no: true or false.
#' @param lagged Whether the network measure function used requires the comparison between two networks. e.g., comparing the current network to one lagged by 10 days. If TRUE the measureFun should take two graphs as input and return a single value. The order of inputs in the function is the lagged network followed by the current network.
#' @param lag If lagged is set to TRUE, this is the lag at which to compare networks.
#' @param cores This allows for multiple cores to be used while generating networks and calculating network measures.
#' @export
#' @importFrom lubridate days
#' @examples
#'
#' ts.out<-dyadTS(data=groomEvents[1:200,])
#'
dyadTS <- function (data, windowsize=days(30), windowshift=days(1), measureFun=dyad_weight, directed=FALSE, lagged=FALSE, lag=1, cores=1){

  #if undirected, organise the data to look at the data without order.
  if(directed==FALSE) data <- order_events(data)

  #extract networks from the dataframe
  if(cores > 1){
    graphlist <- extract_networks_para(data, windowsize, windowshift, directed, cores = 2)
  } else {
    graphlist <- extract_networks(data, windowsize, windowshift, directed)
  }

  #extract measures from the network list
  all.unique.names <- unique(paste0(as.character(data[,1]),"_",as.character(data[,2]) ))
  if(lagged==FALSE){
    values <- extract_measure_dyads(netlist=graphlist, measureFun, unique.names = all.unique.names)
  } else {
    values <- extract_lagged_measure_dyads(graphlist, measureFun, lag, unique.names = all.unique.names)
  }

  return (values)

}


#' Order events alphabetically
#'
#' This function will take an event dataframe and order the names alphabetically. Useful in undirected graphs to ensure all edges have the same name through time.
#' @param data Dataframe containing events between individuals/objects
#'
order_events <- function(data){

  data[,1] <- as.character(data[,1])
  data[,2] <- as.character(data[,2])

  #fill in and order the events
  for(i in 1:nrow(data)){

    #get names of individuals/nodes
    names.found <- data[i,1:2][1,]
    if(as.character(names.found[1,1])>as.character(names.found[1,2]) )names.found<-names.found[,c(2,1)]

    #re-order
    data[i,1]<-names.found[1,1]
    data[i,2]<-names.found[1,2]

  }

  return (data)
}


#' Extract dyad level network measures from a list of networks
#'
#' This function will estimate dyad level network measures from a list of networks.
#' @param netlist List of networks.
#' @param measureFun A function that takes a network as input and returns a value for each dyad.
#' @param unique.names A list of all dyads/edges in the networks.
#' @export
#' @importFrom igraph get.graph.attribute
#' @importFrom dplyr bind_rows
#'
#'
#'
extract_measure_dyads<-function(netlist, measureFun, unique.names){

  #store measures - set global dataframe with proper names
  netvalues <- data.frame(t(rep(-1,length(unique.names))))
  names(netvalues) <- unique.names
  net.measure <- data.frame(nEvents=-1,windowstart=ymd("2000-01-01"), windowend=ymd("2000-01-01"))
  netvalues<-cbind(netvalues,net.measure)

  #extract measures
  if(exists('measureFun', mode='function')){

    for(i in 1:length(netlist)) {
      df.temp.nodes <- data.frame((measureFun(netlist[[i]])))
      df.temp.graph <- data.frame(nEvents=igraph::get.graph.attribute(netlist[[i]], "nEvents" ),
                                  windowstart=igraph::get.graph.attribute(netlist[[i]], "windowstart" ),
                                  windowend=igraph::get.graph.attribute(netlist[[i]], "windowend" ))
      df.temp <- cbind(df.temp.nodes,df.temp.graph)
      netvalues <- bind_rows(netvalues,df.temp)
    }

  } else {
    print("Error: the measurment function was not found.")
  }

  netvalues<-netvalues[-1,]
  return(netvalues)

}


#' Extract dyad level measures from a list of networks when the measure requires comparisons between networks.
#'
#' This function will estimate dyad level network measures from a list of networks.
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
extract_lagged_measure_dyads<-function(netlist, measureFun, lag=1, unique.names){

  #store measures - set global dataframe with proper names
  netvalues <- data.frame(t(rep(-1,length(unique.names))))
  names(netvalues) <- unique.names
  net.measure <- data.frame(nEvents=-1,windowstart=ymd("2000-01-01"), windowend=ymd("2000-01-01"))
  netvalues<-cbind(netvalues,net.measure)

  if(exists('measureFun', mode='function')){

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
    print("Error: the measurment function was not found.")
  }

  netvalues<-netvalues[-1,]
  return(netvalues)

}


#' Trim dyads.
#'
#' This function removes dyad values beyond their min and max observed times.
#' @param dyadvalues Output from the dyadTS function.
#' @param data The events dataframe used in the dyadTS function.
#' @param directed Whether to treat the dyads as directed or not
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @export
#'
trim_dyads<-function(dyadvalues, data ,directed){

  #Ensure the names of the first four columns
  names(data)[1:4]<- c("from","to","weight","date")

  #which names to keep
  names.kept<-colnames(dyadvalues)[1:(length(dyadvalues)-3)]

  #Initialize trimed dataframes with important vars
  df.trim<- data.frame(remove=rep(NA,nrow(dyadvalues)))

  #loop through each ID and trim based on min and max date observed
  for(i in 1:length(names.kept)){

    #get both nodes from dyad
    node1 <- strsplit(names.kept[i], split="_")[[1]][1]
    node2 <- strsplit(names.kept[i], split="_")[[1]][2]

    #determine the min and max dates the focal was seen
    if(directed==FALSE){

      #get all occurances of the dyad
      df.temp <- data %>% filter( (from == node1 & to == node2) |  (from == node2 & to == node1) )
      if(nrow(df.temp)>0){
        min.date<-min(df.temp$date)
        max.date<-max(df.temp$date)
      } else {
        min.date<-NA
        max.date<-NA
      }


    } else {

      #get all occurances of the dyad
      df.temp <- data %>% filter( (from == node1 & to == node2))
      if(nrow(df.temp)>0){
        min.date<-min(df.temp$date)
        max.date<-max(df.temp$date)
      } else {
        min.date<-NA
        max.date<-NA
      }

    }

    if(is.na(min.date)==FALSE){

      #remove all window estimates outside those dates
      df.temp2 <- dyadvalues %>% dplyr::select(names.kept[i], windowstart, windowend)
      df.temp2[,1] <- ifelse(df.temp2[,2]<min.date,NA,df.temp2[,1])
      df.temp2[,1] <- ifelse(df.temp2[,3]>max.date,NA,df.temp2[,1])

      #build new trimed dataframes
      df.temp3 <- data.frame((df.temp2[,1]))
      names(df.temp3) <- c(names.kept[i])
      df.trim <- cbind(df.trim, df.temp3)

    } else {

      #build new trimed dataframes
      df.temp3 <- data.frame(rep(NA,nrow(df.trim) ))
      names(df.temp3) <- c(names.kept[i])
      df.trim <- cbind(df.trim, df.temp3)

    }


  }

  df.trim<-df.trim[,-1]
  df.times <- dyadvalues %>% dplyr::select("nEvents","windowstart","windowend")
  df.trim <- cbind(df.trim,df.times)


  return(df.trim)

}


