################# Some network functions for measurment extraction ###########################


############################################
#
#  Dyad level measurement functions
#
############################################



#' dyad_weight function
#'
#' This function will take a graph and take dyad level sum of weights.
#' @param g graph to extract dyad measures from
#' @export
#' @importFrom igraph E get.edgelist
#' @importFrom plyr rbind.fill
#'
dyad_weight <- function(g){

  weight.vector <- vector()

  for(i in 1:length(E(g))){

    weight.value=E(g)[i]$weight

    weight.vector[length(weight.vector)+1] <- weight.value
  }

  names(weight.vector)<-paste(get.edgelist(g)[,1],get.edgelist(g)[,2],sep="_")

  return (weight.vector)

}


#' dyad_mean function
#'
#' This function will take a graph and take dyad level mean weight.
#' @param g graph to extract dyad measures from
#' @export
#' @importFrom igraph E is.directed which_mutual ends get.edgelist
#' @importFrom plyr rbind.fill
#'
dyad_mean <- function(g){

  mean.vector <- vector()

  for(i in 1:length(E(g))){

    mean.value=E(g)[i]$weight

    #check to see if there is a reciprical edge
    if(is.directed(g) & which_mutual(g,E(g)[i])){
      mean.value = mean( E(g)[ends(g,E(g)[i])[,2] %--% ends(g,E(g)[i])[,1]]$weight)
    }

    mean.vector[length(mean.vector)+1] <- mean.value
  }

  names(mean.vector)<-paste(get.edgelist(g)[,1],get.edgelist(g)[,2],sep="_")

  return (mean.vector)

}

#' dyad_diff function
#'
#' This function will take a graph and take dyad level difference in weights.
#' @param g graph to extract dyad measures from
#' @export
#' @importFrom igraph E is.directed which_mutual ends get.edgelist
#' @importFrom plyr rbind.fill
#'
dyad_diff <- function(g){

  diff.vector <- vector()

  for(i in 1:length(E(g))){

    diff.value=0

    #check to see if there is a reciprical edge
    if(is.directed(g) & which_mutual(g,E(g)[i])){
      diff.values = ( E(g)[ends(g,E(g)[i])[,2] %--% ends(g,E(g)[i])[,1]]$weight)
      diff.value = abs(diff.values[1]-diff.values[2])
    }

    diff.vector[length(diff.vector)+1] <- diff.value
  }

  names(diff.vector)<-paste(get.edgelist(g)[,1],get.edgelist(g)[,2],sep="_")

  return (diff.vector)

}


#' dyad_change function
#'
#' This function will take two graphs and take the difference in dyad weight.
#' @param g1 First graph.
#' @param g2 Second graph.
#' @param directed Wether the networks are directed (Default=FALSE).
#' @export
#' @importFrom igraph E head_of tail_of as_ids is.directed which_mutual ends get.edgelist
#' @importFrom plyr rbind.fill
#'
dyad_change <- function(g1,g2, directed=FALSE){

  #get weights by dyad for g1
  g1.dyad.names <- vector()
  for(i in 1:length(E(g1))){
    g.head <- head_of(g1, E(g1)[i])
    g.tail <- tail_of(g1, E(g1)[i])
    dyad.name.temp <- paste0(as.character(as_ids(g.tail)),"_",as.character(as_ids(g.head)))

    #correct for direction if undirected (simply write the dyad in alphabetical order)
    if(directed==FALSE){
      if(g.head>g.tail){
        dyad.name.temp <- paste0(as.character(as_ids(g.tail)),"_",as.character(as_ids(g.head)))
      } else{
        dyad.name.temp <- paste0(as.character(as_ids(g.head)),"_",as.character(as_ids(g.tail)))
      }
    }

    #record dyad name
    g1.dyad.names[length(g1.dyad.names)+1] <- dyad.name.temp

  }

  g1.df <- data.frame(t(E(g1)$weight))
  names(g1.df) <- g1.dyad.names

  #get weights by dyad for g2
  g2.dyad.names <- vector()
  for(i in 1:length(E(g2))){
    g.head <- head_of(g2, E(g2)[i])
    g.tail <- tail_of(g2, E(g2)[i])
    dyad.name.temp <- paste0(as.character(as_ids(g.tail)),"_",as.character(as_ids(g.head)))

    #correct for direction if undirected
    if(directed==FALSE){
      if(g.head>g.tail){
        dyad.name.temp <- paste0(as.character(as_ids(g.tail)),"_",as.character(as_ids(g.head)))
      } else{
        dyad.name.temp <- paste0(as.character(as_ids(g.head)),"_",as.character(as_ids(g.tail)))
      }
    }

    #record dyad name
    g2.dyad.names[length(g2.dyad.names)+1] <- dyad.name.temp

  }

  g2.df <- data.frame(t(E(g2)$weight))
  names(g2.df) <- g2.dyad.names


  #compare the two
  df.total <- bind_rows(g1.df,g2.df)
  df.total[ is.na(df.total)] <- 0
  diff <- df.total[2,]-df.total[1,]
  diff<-as.vector(diff)

  return(diff)
}





