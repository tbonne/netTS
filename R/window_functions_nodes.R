

#' nodeTS function
#'
#' This function will take a dataframe with events between individuals/objects, and take network measures using a moving window approach.
#' @param data A dataframe with relational data in the first two columns, and a time stamp in the third columns Note: time stamps should be in ymd_hms format. The lubridate package can be very helpful in organizing times. Note: names in the first two columns are case sensitive.
#' @param windowsize The size of the moving window in which to take network measures. These should be provided as e.g., days(30), hours(5), ... etc.
#' @param windowshift The amount to shift the moving window for each measure. Again times should be provided as e.g., days(1), hours(1), ... etc.
#' @param measureFun This is a function that takes as an input a igraph network and returns values for each node in the network. There are functions within netTS (see details), and custom made functions can be used.
#' @param effortFun This is a function that takes as input the data within a window of time and returns the total sampling effort.
#' @param effortData This is a dataframe containing the data used to calculate sampling effort.
#' @param directed Whether the events are directed or no: true or false.
#' @param lagged Whether the network measure function used requires the comparison between two networks. e.g., comparing the current network to one lagged by 10 days. If TRUE the measureFun should take two graphs as input and return a single value. The order of inputs in the function is the lagged network followed by the current network.
#' @param lag If lagged is set to TRUE, this is the lag at which to compare networks.
#' @param firstNet If lagged is set to TRUE, this compares the subsequent networks to the first network.
#' @param cores This allows for multiple cores to be used while generating networks and calculating network measures.
#' @param SRI Wether to use the simple ratio index (Default=FALSE).
#' @param nperm This allows for the estimation the network measure assuming random permutations. Currently the 95 percent quantiles are returned.
#' @param permutationFun This is a function that takes as input an events dataframe and a measurment function, and will return a quantile range of network values (see permutation vignette).
#' @param probs When nperm > 0 this will determine the probability of the permutation values returned from the permuations.
#' @export
#' @importFrom lubridate days
#' @importFrom igraph degree
#' @examples
#'
#' ts.out<-nodeTS(data=groomEvents)
#'
nodeTS <- function (data,windowsize =days(30), windowshift= days(1), measureFun=degree, effortFun=NULL,effortData=NULL,directed=FALSE, lagged=FALSE, lag=1, firstNet=FALSE, cores=1, permutationFun=perm.events.multiple.outputs, nperm=0, probs=0.95, SRI=FALSE){

  #check for missing data
  if(sum(is.na(data)) > 0){
    print(paste0("Data contains NA, removing ",sum(!complete.cases(data)), " row(s)."))
    data<-data[complete.cases(data),]
  }

  #extract networks from the dataframe
  if(cores > 1){
    graphlist <- extract_networks_para(data, windowsize, windowshift, directed, cores = 2, effortFun = effortFun, effortData = effortData)
  } else {
    graphlist <- extract_networks(data, windowsize, windowshift, directed, effortFun = effortFun, effortData = effortData)
  }



  #extract measures from the network list
  all.unique.names <- unique(c(as.character(data[,1]),as.character(data[,2]) ))
  if(lagged==FALSE){
    values <- extract_measure_nodes(netlist=graphlist, measureFun, unique.names = all.unique.names)
  } else {
    values <- extract_lagged_measure_nodes(graphlist, measureFun, lag, unique.names = all.unique.names, firstNet=firstNet)
  }

  #run permutations and extract measures
  if(nperm>0){
    perm.values <- permutation.multi.values(data, windowsize, windowshift, directed, measureFun = measureFun, probs=probs, SRI=SRI, graphlist = graphlist,permutationFun=permutationFun, nperm=nperm, unique.names = all.unique.names )
    values <- list(obs=values,lowCI=cbind(perm.values$low,values[c("windowstart","windowend","nEvents") ]),highCI= cbind(perm.values$high,values[c("windowstart","windowend","nEvents") ]) )
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
  names(netvalues) <- as.character(unique.names)
  net.measure <- data.frame(nEvents=-1,windowstart=igraph::get.graph.attribute(netlist[[1]], "windowstart" ), windowend=igraph::get.graph.attribute(netlist[[1]], "windowend" ))
  netvalues<-cbind(netvalues,net.measure)

  #extract measures
  if(exists('measureFun', mode='function')){

    for(i in 1:length(netlist)) {
      df.temp.nodes <- as.data.frame(t(measureFun(netlist[[i]])))
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
#' @param firstNet Wether to make all network comparisons to the first network observed (Default=FALSE).
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
  net.measure <- data.frame(nEvents=-1,windowstart=igraph::get.graph.attribute(netlist[[1]], "windowstart" ), windowend=igraph::get.graph.attribute(netlist[[1]], "windowend" ))
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
#' @param enter_leave Optional: a data frame that specifies the entering and leaving dates of each node (ID, enterDate, leaveDate)
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @export
#'
trim_nodes<-function(nodevalues, data, enter_leave=NULL){

  #ensure the names of the first four columns
  names(data)[1:3]<- c("from","to","date")

  #which names to keep
  names.kept<-colnames(nodevalues)[1:(length(nodevalues)-3)]

  #Initialize trimed dataframes with important vars
  df.trim<- data.frame(remove=rep(NA,nrow(nodevalues)))

  #(Default) use first and last observations to set trim
  if(is.null(enter_leave)){

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

    #Use a user specified dates
  } else {

    #ensure column names
    names(enter_leave)[1:3]<- c("ID","enter","leave")

    #loop through each ID and trim based on min and max date specified
    for(i in 1:length(names.kept)){

      #determine the min and max dates the focal was seen
      df.temp <- enter_leave %>% filter(ID == names.kept[i])
      min.date<-df.temp[,2]
      max.date<-df.temp[,3]

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
  }

  return(df.trim)

}

