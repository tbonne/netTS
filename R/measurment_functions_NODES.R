
###################################################
#
# Functions for network measures around nodes
#
###################################################











#' Estimate consine similarity between nodes in a graph
#'
#' This function will calculate the cosine similarity between nodes at two time periods.
#' @param graph1 first graph (igraph)
#' @param graph2 second graph (igraph)
#' @param directed Boolean indicating if the graph is directed or not.
#' @param mode For directed graphs this allows for in, out, or total cosine similarity to be calculated for each node. This argument is ignored for undirected graphs.
#' @export
#' @import igraph
#' @importFrom dplyr full_join
#' @examples
#'
#' #two random graphs
#' library(igraph)
#' graph1 <-erdos.renyi.game(n=10,p=0.1)
#' graph2 <-erdos.renyi.game(n=10,p=0.1)
#' cosine_between_nodes(graph1,graph2)
#'
#'
cosine_between_nodes<- function(graph1, graph2, directed=FALSE, mode="out"){

  node.cosine <- vector()

  #create weighted edge list from first graph
  g1.edges<-as.data.frame(get.edgelist(graph1, names=TRUE))
  if(is.null(igraph::E(graph1)$weight)){
    g1.edges$weight <- rep(1,nrow(g1.edges))
  } else {
    g1.edges$weight<-igraph::E(graph1)$weight
  }
  g1.edges$joinC <- ifelse(as.character(g1.edges$V1) < as.character(g1.edges$V2), paste(g1.edges$V1, g1.edges$V2), paste(g1.edges$V2, g1.edges$V1))

  #create weighted edge list from second graph
  g2.edges<-as.data.frame(get.edgelist(graph2, names=TRUE))
  if(is.null(igraph::E(graph2)$weight)){
    g2.edges$weight <- rep(1,nrow(g2.edges))
  } else {
    g2.edges$weight<-igraph::E(graph2)$weight
  }
  g2.edges$joinC <- ifelse(as.character(g2.edges$V1) < as.character(g2.edges$V2), paste(g2.edges$V1, g2.edges$V2), paste(g2.edges$V2, g2.edges$V1))

  #join the weighted edge lists, setting NAs equal to 0
  comb<-full_join(g1.edges,g2.edges,by="joinC")
  comb$weight.x[is.na(comb$weight.x)]<-0
  comb$weight.y[is.na(comb$weight.y)]<-0

  names.unique<-unique(c(as.character(comb$V1.x),as.character(comb$V2.x) ))
  names.unique <- names.unique[is.na(names.unique)==FALSE]

  if(directed==FALSE){

    for(i in 1:(length(names.unique))){
      temp.node.w <- dplyr::filter(comb, V1.x==names.unique[i] | V2.x==names.unique[i])
      node.cosine[length(node.cosine)+1]<-lsa::cosine(temp.node.w$weight.x,temp.node.w$weight.y)
    }

  } else if (directed==TRUE){

    if(mode=="out"){
      for(i in 1:(length(names.unique))){
        temp.node.w <- dplyr::filter(comb, V1.x==names.unique[i])
        node.cosine[length(node.cosine)+1]<-lsa::cosine(temp.node.w$weight.x,temp.node.w$weight.y)
      }
    }

    if(mode=="in"){
      for(i in 1:(length(names.unique))){
        temp.node.w <- dplyr::filter(comb, V2.x==names.unique[i])
        node.cosine[length(node.cosine)+1]<-lsa::cosine(temp.node.w$weight.x,temp.node.w$weight.y)
      }
    }

    if(mode=="total"){
      for(i in 1:(length(names.unique))){
        temp.node.w <- dplyr::filter(comb, V1.x==names.unique[i] | V2.x==names.unique[i])
        node.cosine[length(node.cosine)+1]<-lsa::cosine(temp.node.w$weight.x,temp.node.w$weight.y)
      }
    }


  }

  #node.cosine[is.nan(node.cosine)]<-0
  cos.df<-as.data.frame(t(node.cosine))
  names(cos.df)<-names.unique

  return(cos.df)

}



#' Estimate skewness of the edge weight distribution for each node in a graph
#'
#' This function will calculate the skewness of edge weights from nodes in a graph.
#' @param graph1 An igraph object.
#' @param type The type of edges to pull from the graph: 'in', 'out', or 'all'. In non-directed graphs this is ignored.
#' @export
#' @import igraph
#' @import moments
#' @examples
#' library(igraph)
#'
#' #Random graph
#' graph1 <- random.graph.game(n=10, p.or.m = 0.3)
#' E(graph1)$weight <- c(1,2,3,4,5,6,7,7,7,7)
#'
#' #skewness
#' #edge.weight.skeweness(graph1)
#'
edge.weight.skewness <- function(graph1, type = "all"){

  if(type == "all"){
    edge.skewness <- vector()
    for(i in 1:vcount(graph1)){
      inc.edges <- incident(graph1,  V(graph1)[i], mode="all")
      edge.skewness[length(edge.skewness)+1]<-skewness(E(graph1)[inc.edges]$weight)
    }

  } else if (type == "out") {
    edge.skewness <- vector()
    for(i in 1:vcount(graph1)){
      inc.edges <- incident(graph1,  V(graph1)[i], mode="out")
      edge.skewness[length(edge.skewness)+1]<-skewness(E(graph1)[inc.edges]$weight)
    }

  } else if(type == "in"){
    edge.skewness <- vector()
    for(i in 1:vcount(graph1)){
      inc.edges <- incident(graph1,  V(graph1)[i], mode="in")
      edge.skewness[length(edge.skewness)+1]<-skewness(E(graph1)[inc.edges]$weight)
    }
  }

  return (edge.skewness)
}


