################# Some network functions for measurment extraction ###########################


############################################
#
#  Graph level measurement functions
#
############################################


#' Mean degree of a network
#'
#' This function will estimate mean degree of a network
#' @param net An igraph network.
#' @importFrom igraph degree
#' @export
#'
#'
degree_mean<-function(net){

  return(mean(degree(net)))

}


#' Mean eigenvector centrality of a network
#'
#' This function will estimate mean eigenvector centrality of a network
#' @param net An igraph network.
#' @importFrom igraph eigen_centrality
#' @export
#'
#'
eigen_mean<-function(net){

  return(mean(eigen_centrality(net)$vector))

}





#' Estimate consine similarity between graphs
#'
#' This function will calculate the cosine similarity between two graphs.
#' @param graph1 first graph (igraph)
#' @param graph2 second graph (igraph)
#' @param directed Whether the graph is directed or not (default = FALSE)
#' @param considerZeros Whether to treat edges present in one graph, but not the other, as zero (TRUE), or consider only those edges present in both graphs (FALSE).
#' @param center Whether to center the resulting vectors before measuring cosine similarity.
#' @export
#' @importFrom igraph E get.edgelist
#' @importFrom dplyr full_join
#' @importFrom lsa cosine
#' @examples
#'
#' #two random graphs
#' library(igraph)
#' graph1 <-erdos.renyi.game(n=10,p=0.1)
#' graph2 <-erdos.renyi.game(n=10,p=0.1)
#' cosine_between_graphs(graph1,graph2)
#'
cosine_between_graphs <- function(graph1,graph2, directed=FALSE, considerZeros=TRUE, center=FALSE){

  #create weighted edge list from first graph
  g1.edges<-as.data.frame(get.edgelist(graph1, names=TRUE))
  if(is.null(igraph::E(graph1)$weight)){
    g1.edges$weight <- rep(1,nrow(g1.edges))
  } else {
    g1.edges$weight<-igraph::E(graph1)$weight
  }
  #create a join column (ordered)
  if(directed==FALSE){
    g1.edges$joinC <- ifelse(as.character(g1.edges$V1) < as.character(g1.edges$V2), paste(g1.edges$V1, g1.edges$V2), paste(g1.edges$V2, g1.edges$V1))
  } else {
    g1.edges$joinC <- paste(g1.edges$V1, g1.edges$V2)
  }

  #create weighted edge list from second graph
  g2.edges<-as.data.frame(get.edgelist(graph2, names=TRUE))
  if(is.null(igraph::E(graph2)$weight)){
    g2.edges$weight <- rep(1,nrow(g2.edges))
  } else {
    g2.edges$weight<-igraph::E(graph2)$weight
  }
  #create a join column (ordered)
  if(directed==FALSE){
    g2.edges$joinC <- ifelse(as.character(g2.edges$V1) < as.character(g2.edges$V2), paste(g2.edges$V1, g2.edges$V2), paste(g2.edges$V2, g2.edges$V1))
  } else {
    g2.edges$joinC <- paste(g2.edges$V1, g2.edges$V2)
  }

  #join the weighted edge lists, setting NAs equal to 0
  comb<-full_join(g1.edges,g2.edges,by="joinC")

  if(considerZeros){
  comb$weight.x[is.na(comb$weight.x)]<-0
  comb$weight.y[is.na(comb$weight.y)]<-0
  } else {
    comb<-comb[complete.cases(comb),]
  }

  #center both vectors
  if(center){
  comb$weight.x = comb$weight.x - mean(comb$weight.x)
  comb$weight.y = comb$weight.y - mean(comb$weight.y)
  }

  return(cosine(comb$weight.x,comb$weight.y))

}


#' Estimate similarity between graphs using correlation
#'
#' This function will calculate pearson's correlation between two graphs.
#' @param graph1 first graph (igraph)
#' @param graph2 second graph (igraph)
#' @param directed Whether the graph is directed or not (default = FALSE)
#' @param considerZeros Whether to treat edges present in one graph, but not the other, as zero (TRUE), or consider only those edges present in both graphs (FALSE).
#' @export
#' @importFrom igraph E get.edgelist
#' @importFrom dplyr full_join
#' @examples
#'
#' #two random graphs
#' library(igraph)
#' graph1 <-erdos.renyi.game(n=10,p=0.1)
#' graph2 <-erdos.renyi.game(n=10,p=0.1)
#' cor_between_graphs(graph1,graph2)
#'
cor_between_graphs <- function(graph1,graph2, directed=FALSE, considerZeros=TRUE){

  #create weighted edge list from first graph
  g1.edges<-as.data.frame(get.edgelist(graph1, names=TRUE))
  if(is.null(igraph::E(graph1)$weight)){
    g1.edges$weight <- rep(1,nrow(g1.edges))
  } else {
    g1.edges$weight<-igraph::E(graph1)$weight
  }
  #create a join column (ordered)
  if(directed==FALSE){
    g1.edges$joinC <- ifelse(as.character(g1.edges$V1) < as.character(g1.edges$V2), paste(g1.edges$V1, g1.edges$V2), paste(g1.edges$V2, g1.edges$V1))
  } else {
    g1.edges$joinC <- paste(g1.edges$V1, g1.edges$V2)
  }

  #create weighted edge list from second graph
  g2.edges<-as.data.frame(get.edgelist(graph2, names=TRUE))
  if(is.null(igraph::E(graph2)$weight)){
    g2.edges$weight <- rep(1,nrow(g2.edges))
  } else {
    g2.edges$weight<-igraph::E(graph2)$weight
  }
  #create a join column (ordered)
  if(directed==FALSE){
    g2.edges$joinC <- ifelse(as.character(g2.edges$V1) < as.character(g2.edges$V2), paste(g2.edges$V1, g2.edges$V2), paste(g2.edges$V2, g2.edges$V1))
  } else {
    g2.edges$joinC <- paste(g2.edges$V1, g2.edges$V2)
  }

  #join the weighted edge lists, setting NAs equal to 0
  comb<-full_join(g1.edges,g2.edges,by="joinC")

  if(considerZeros){
    comb$weight.x[is.na(comb$weight.x)]<-0
    comb$weight.y[is.na(comb$weight.y)]<-0
  } else {
    comb<-comb[complete.cases(comb),]
  }

  return(cor.test(comb$weight.x,comb$weight.y))

}


#' Estimate similarity between graphs using Euclidean distance
#'
#' This function will calculate Euclidean distance between two graphs.
#' @param graph1 first graph (igraph)
#' @param graph2 second graph (igraph)
#' @param directed Whether the graph is directed or not (default = FALSE)
#' @param considerZeros Whether to treat edges present in one graph, but not the other, as zero (TRUE), or consider only those edges present in both graphs (FALSE).
#' @export
#' @importFrom igraph E get.edgelist
#' @importFrom dplyr full_join
#' @examples
#'
#' #two random graphs
#' library(igraph)
#' graph1 <-erdos.renyi.game(n=10,p=0.1)
#' graph2 <-erdos.renyi.game(n=10,p=0.1)
#' cor_between_graphs(graph1,graph2)
#'
dist_between_graphs <- function(graph1,graph2, directed=FALSE, considerZeros=TRUE){

  #create weighted edge list from first graph
  g1.edges<-as.data.frame(get.edgelist(graph1, names=TRUE))
  if(is.null(igraph::E(graph1)$weight)){
    g1.edges$weight <- rep(1,nrow(g1.edges))
  } else {
    g1.edges$weight<-igraph::E(graph1)$weight
  }
  #create a join column (ordered)
  if(directed==FALSE){
    g1.edges$joinC <- ifelse(as.character(g1.edges$V1) < as.character(g1.edges$V2), paste(g1.edges$V1, g1.edges$V2), paste(g1.edges$V2, g1.edges$V1))
  } else {
    g1.edges$joinC <- paste(g1.edges$V1, g1.edges$V2)
  }

  #create weighted edge list from second graph
  g2.edges<-as.data.frame(get.edgelist(graph2, names=TRUE))
  if(is.null(igraph::E(graph2)$weight)){
    g2.edges$weight <- rep(1,nrow(g2.edges))
  } else {
    g2.edges$weight<-igraph::E(graph2)$weight
  }
  #create a join column (ordered)
  if(directed==FALSE){
    g2.edges$joinC <- ifelse(as.character(g2.edges$V1) < as.character(g2.edges$V2), paste(g2.edges$V1, g2.edges$V2), paste(g2.edges$V2, g2.edges$V1))
  } else {
    g2.edges$joinC <- paste(g2.edges$V1, g2.edges$V2)
  }

  #join the weighted edge lists, setting NAs equal to 0
  comb<-full_join(g1.edges,g2.edges,by="joinC")

  if(considerZeros){
    comb$weight.x[is.na(comb$weight.x)]<-0
    comb$weight.y[is.na(comb$weight.y)]<-0
  } else {
    comb<-comb[complete.cases(comb),]
  }

  return(dist(rbind(comb$weight.x,comb$weight.y)))

}



