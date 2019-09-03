################# Some functions for measurment of sampling effort ###########################


############################################
#
#  Effort functions
#
############################################


#' Min/Max time per day
#'
#' This function will estimate the total time spent sampling, using the min and max observed values per day. It will export a network where the edges are corrected for sampling time.
#' @param df.window An events data frame with a column containing the date and time.
#' @importFrom dplyr group_by summarise
#' @importFrom lubridate minute hour date interval as.duration day
#' @importFrom stats time
#' @export
#'
#'
effort.time <- function(df.window, directed=FALSE){

  #get days and times
  names(df.window)[3] <- c("date")
  df.window$day <- time_length(lubridate::interval(lubridate::date(min(df.window$date)),lubridate::date(df.window$date)  ),"days")
  df.window$time <- time(df.window$date)

  #get min and max times per day
  hours.min <- df.window %>% dplyr::group_by(day) %>% dplyr::summarise(min=min(lubridate::hour(date) + (lubridate::minute(date) )/60 ))
  hours.max <- df.window %>% dplyr::group_by(day) %>% dplyr::summarise(max=max(lubridate::hour(date) + (lubridate::minute(date) )/60 ))

  #get the sum of the difference between max and min times within a day
  sample.time = sum(hours.max-hours.min)

  #create a network with the edge corrected values
  g <- create.a.network(df.window, directed)
  E(g)$weight <- E(g)$weight/sample.time
  g <- set_graph_attr(g, "effort", sample.time)

  return(g)
}


#' Unique scan IDs
#'
#' This function will correct edge weights in a network, using a count of the number of scans used to construct the network.
#' @param df.window The events data frame with interactions.
#' @df.scans A dataframe with the first column as a data column, and a second column with the number of scans durring that day.
#' @param directed Whether a directed network should be constructed (Default=FALSE).
#' @export
#'
#'
effort.scan <- function(df.window, df.scans, directed=FALSE){

  #calculate how many scans there are in the window
  total.samples <- sum(df.scans$scanID)

  #create a network with the edge corrected values
  g <- create.a.network(df.window, directed)
  E(g)$weight <- E(g)$weight/total.samples
  g <- set_graph_attr(g, "effort", total.samples)

  #return value
  return(g)
}


#' Focal sampling
#'
#' This function will estimate the sampling effort from focal data. It will calculate the total durration each dyad had the potential to be observed (i.e., Interactions between A and B could only be seen if A or B was the focal). This sampling correction should be applied post network construction.
#' @param df.window A data frame with focal duration per individual per unique focal, with a date time column.
#' @param effortData A data frame with the first column as a datatime stamp, a second column with the ID of the focal, and a thrid columbn with the total time of the focal.
#' @param directed Whether a directed network should be constructed (Default=FALSE).
#' @export
#'
#'
effort.focal <- function(df.window,effortData, directed=FALSE){

  #create a network
  g<-create.a.network(df.window,directed=directed)

  #for each dyad calculate the observation time and correct for sampling effort (i.e., focal time of the individuals that make up the dyad)
  list_of_edges <- E(g)
  cor.edge.weights = vector()
  tot.effort = 0
  for(e in 1:length(list_of_edges) ){
    ids = ends(g,list_of_edges[e])
    tot = sum((effortData[effortData[,2]==ids[1],])[,3]) + sum((effortData[effortData[,2]==ids[2],])[,3])
    tot.effort = tot.effort + tot
    cor.edge.weights[length(cor.edge.weights)+1] <- E(g)$weight[e]/tot
  }
  E(g)$weight = cor.edge.weights

  #add attribute
  g <- set_graph_attr(g, "effort", tot.effort)

  return(g)
}

