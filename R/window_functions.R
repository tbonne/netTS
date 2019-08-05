#' Create a network
#'
#' This function will generate a network from an events dataframe.
#' @param data Dataframe containing all events. The first two coloums should contain the 'to' and 'from' indentities involved in the interaction, while the third column should contain the 'weight' of the interaction.
#' @param directed Treat the network as directed or not (Default=FALSE).
#' @param SRI Wether to use the simple ratio index or not (Default=FALSE).
#' @param effort The number of scans, or measure of sampling effort (e.g., hours sampling). The default of 1 assumes equal sampling effort througout.
#' @importFrom igraph graph_from_data_frame simplify
#' @export
#'
create.a.network<-function(data, directed = FALSE, SRI=FALSE, effort=1){

  if(SRI==FALSE){
    elist<-create.an.edgeList(data,effort)
    gg <- graph_from_data_frame(elist, directed = directed, vertices = NULL)

    if(is.simple(gg)==FALSE)gg<-simplify(gg, edge.attr.comb=list(weight="sum"))

  } else {

    gg<-create.a.network.SRI(data,directed=directed)

    if(is.simple(gg)==FALSE)gg<-simplify(gg, edge.attr.comb=list(weight="sum"))

  }

  return(gg)

}



#' Create an edge list
#'
#' This function will generate an edge list from an events dataframe
#' @param data Dataframe containing all events. The first two coloums should contain the 'to' and 'from' indentities involved in the interaction, while the third column should contain the 'weight' of the interaction.
#' @param effort The number of scans, or measure of sampling effort (e.g., hours sampling). The default of 1 assumes equal sampling effort througout.
#' @import data.table
#' @export
#'
create.an.edgeList<-function(data,effort=1){

  #create a network and add it to the list
  names(data)[1:2]<-c("from","to")
  if(is.null(data$weight))data$weight=1
  elist <- as.data.frame(as.data.table(data)[,.(sum(weight)/effort), by=list(from,to)])
  names(elist)<-c("from","to","weight")

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

  df.window <- data[data[[3]] >= start & data[[3]] < end,]

  return (df.window)
}


#' Create a network from an edge list using the SRI.
#'
#' This function will generate an SRI network from an event dataframe. Right now only directed graphs are considered.
#' @param events dataframe containing all events.
#' @param directed Wether the network is directed or not (Default=FALSE).
#' @importFrom stats complete.cases
#' @importFrom dplyr summarise filter
#' @importFrom igraph graph_from_data_frame
#' @export
#'
create.a.network.SRI <- function(events, directed=FALSE){

  if(directed==TRUE){

  elist<-as.data.frame(create.an.edgeList(events))
  elist.sri <-  data.frame(from="NA",to="NA",weight=-1, stringsAsFactors=FALSE)
  for(i in 1:nrow(elist)){

    #how much A interacts with B
    Nab <- elist %>% dplyr::filter(from==as.character(elist[i,][,1]) ) %>% dplyr::filter(to==as.character(elist[i,][,2]) ) %>% summarise(x=sum(weight))

    #how much A interacts with not B
    Na <- elist %>% dplyr::filter(from==as.character(elist[i,][,1]) ) %>% dplyr::filter(to!=as.character(elist[i,][,2]) ) %>% summarise(x=sum(weight))

    #how much B gets interacted with by not A
    Nb <- elist %>% dplyr::filter(to==as.character(elist[i,][,2]) ) %>% dplyr::filter(from!=as.character(elist[i,][,1]) ) %>% summarise(x=sum(weight))

    #if no events
    if(nrow(Na)==0)Na=0
    if(nrow(Nb)==0)Nb=0

    #add SRI to the list
    df.temp<-data.frame(from=elist[i,][,1],to=elist[i,][,2],weight=as.numeric(Nab / (Nab + Na + Nb)))
    elist.sri <- rbind(elist.sri,df.temp)

  }

  elist.sri <- elist.sri[-1,]

  names(elist.sri)<- c("from", "to", "weight")
  elist.sri <- elist.sri[complete.cases(elist.sri),]
  gg <- graph_from_data_frame(elist.sri, directed = TRUE, vertices = NULL)

  } else { #Undirected SRI

    elist<-as.data.frame(create.an.edgeList(events))
    elist.sri <- data.frame(from="-1",to="-1",weight=1 , stringsAsFactors=FALSE)

    elist<-order_events(as.data.frame(elist))

    for(i in 1:nrow(elist)){

      #how much A and B interact
      Nab <- elist %>% dplyr::filter(from==as.character(elist[i,][,1]) ) %>% dplyr::filter(to==as.character(elist[i,][,2]) ) %>% summarise(x=sum(weight))

      #how much A interacts with not B
      Na1 <- elist %>% dplyr::filter(from==as.character(elist[i,][,1]) ) %>% dplyr::filter(to!=as.character(elist[i,][,2]) ) %>% summarise(x=sum(weight))
      Na2 <- elist %>% dplyr::filter(to==as.character(elist[i,][,1]) ) %>% dplyr::filter(from!=as.character(elist[i,][,2]) ) %>% summarise(x=sum(weight))
      Na <- Na1+Na2

      #how much B interacts with not A
      Nb1 <- elist %>% dplyr::filter(to==as.character(elist[i,][,2]) ) %>% dplyr::filter(from!=as.character(elist[i,][,1]) ) %>% summarise(x=sum(weight))
      Nb2 <- elist %>% dplyr::filter(from==as.character(elist[i,][,2]) ) %>% dplyr::filter(to!=as.character(elist[i,][,1]) ) %>% summarise(x=sum(weight))
      Nb <- Nb1 + Nb2

      #if no events
      if(nrow(Na)==0)Na=0
      if(nrow(Nb)==0)Nb=0

      #add SRI to the list
      df.temp<-data.frame(from=elist[i,][,1],to=elist[i,][,2],weight=as.numeric(Nab / (Nab + Na + Nb)))
      elist.sri <- rbind(elist.sri,df.temp)

    }

    elist.sri <- elist.sri[-1,]

    names(elist.sri)<- c("from", "to", "weight")
    elist.sri <- elist.sri[complete.cases(elist.sri),]
    gg <- graph_from_data_frame(elist.sri, directed = FALSE, vertices = NULL)

  }

  return(gg)
}



