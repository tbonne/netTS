#' Create a network
#'
#' This function will generate a network from an events dataframe.
#' @param data Dataframe containing all events. The first two coloums should contain the 'to' and 'from' indentities involved in the interaction, while the third column should contain the 'weight' of the interaction.
#' @param directed Treat the network as directed or not (Default = FALSE)
#' @import igraph
#' @export
#'
create.a.network<-function(data, directed = FALSE){

  elist<-create.an.edgeList(data)
  gg <- graph_from_data_frame(elist, directed = directed, vertices = NULL)

  if(is.simple(gg)==FALSE)gg<-simplify(gg, edge.attr.comb=list(weight="sum"))

  return(gg)

}



#' Create an edge list
#'
#' This function will generate an edge list from an events dataframe
#' @param data Dataframe containing all events. The first two coloums should contain the 'to' and 'from' indentities involved in the interaction, while the third column should contain the 'weight' of the interaction.
#' @import dplyr
#' @export
#'
create.an.edgeList<-function(data){

  #create a network and add it to the list
  names(data)[1:3]<-c("to","from","weight")
  elist<-data %>% dplyr::group_by(.dots=c("to","from")) %>% summarise(sum(weight))
  names(elist)<-c("to","from","weight")

  return(elist)
}

#' Create a window
#'
#' This function will generate a window from a larger event dataframe, based on the time of the events
#' @param data dataframe containing all events
#' @param start the starting time of the window
#' @param end the ending time of the window
#' @export
#'
#'
create.window <- function(data, start, end){

  df.window <- data[data[[4]] >= start & data[[4]] < end,]

  return (df.window)
}







