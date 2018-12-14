#' Perform permutation on the events dataframe
#'
#' This function will permute a network by randomly switching individuals between events.
#' @param data Dataframe with relational data in the first two rows, and a time stamp in the third row. Note: time stamps should be in ymd or ymd_hms format. The lubridate package can be very helpful in organizing times.
#' @param measureFun This is a function that takes as an input a igraph network and returns a single value.
#' @param directed Whether to consider the network as directed or not (TRUE/FALSE).
#' @param nperm Number of permutations to perform before extracting network measures.
#' @param probs numeric vector of probabilities with values in [0,1].
#' @export
#'
#'
perm.events <- function(data, measureFun, directed=FALSE, nperm=1000, probs=0.95, SRI, effort=1){

  Perm.measure<-vector()

  for(i in 1:nperm){

    no.loops= FALSE

    while(no.loops == FALSE){

      #choose two rows to permute individuals
      rows.to.switch <- sample(1:nrow(data),2,F)

      #record old order
      old.order <- data[rows.to.switch,1:2]
      new.order <- data[rows.to.switch,1:2]

      #choose who to permute
      r1<-sample(1:2,1)
      r2<-sample(1:2,1)
      perm.1<- old.order[1,r1]
      perm.2<- old.order[2,r2]

      #switch individuals
      new.order[1,r1] <- as.character(perm.2)
      new.order[2,r2] <- as.character(perm.1 )

      #check to make sure there are no self loops
      if(sum(as.character(new.order[,1])==as.character(new.order[,2]) )==0){

        #record switch in the dataframe
        data[rows.to.switch,1:2] <- new.order
        no.loops=TRUE

      }
    }

    #Create graph in order to get the measure
    Perm.network <- create.a.network(data, directed, SRI = SRI, effort=effort)

    # Get measure
    Perm.measure[length(Perm.measure)+1]<- measureFun(Perm.network)

  }

  probs.left<-1-probs
  return(quantile(Perm.measure, probs = c( (0+probs.left/2), (1-probs.left/2) ), na.rm=T))
}




#' Perform permutation on the events dataframe, maintaining in/out degree of nodes
#'
#' This function will permute a network by randomly switching individuals between events, keeping the oringinal order of the interactions (e.g., individuals in the from or to columns are permuted, but no permutations between from and to columns).
#' @param data Dataframe with relational data in the first two rows, with weights in the thrid row, and a time stamp in the fourth row. Note: time stamps should be in ymd or ymd_hms format. The lubridate package can be very helpful in organizing times.
#' @param measureFun This is a function that takes as an input a igraph network and returns a single value.
#' @param directed Whether to consider the network as directed or not (TRUE/FALSE).
#' @param nperm Number of permutations to perform before extracting network measures.
#' @param probs numeric vector of probabilities with values in [0,1].
#' @export
#'
#'
perm.events.directed <- function(data, measureFun, directed=FALSE, nperm=1000, probs=0.95, SRI, effort=1){

  net.list <- list(create.a.network(data, directed,SRI=SRI))
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
    Perm.network <- create.a.network(NewData, directed, SRI = SRI, effort=effort)

    # Get measure
    Perm.measure[length(Perm.measure)+1]<- measureFun(Perm.network)

  }

  probs.left<-1-probs
  return(quantile(Perm.measure, probs = c( (0+probs.left/2), (1-probs.left/2) ), na.rm=T))
}


#' Perform edge weight permutations
#'
#' This function will permute a network by randomly switching weights between edges.
#' @param data Dataframe with relational data in the first two rows, and a time stamp in the third row. Note: time stamps should be in ymd or ymd_hms format. The lubridate package can be very helpful in organizing times.
#' @param measureFun This is a function that takes as an input a igraph network and returns a single value.
#' @param directed Whether to consider the network as directed or not (TRUE/FALSE).
#' @param nperm Number of permutations to perform before extracting network measures.
#' @param probs numeric vector of probabilities with values in [0,1].
#' @importFrom igraph E
#' @export
#'
#'
perm.edge.weights <- function(data, measureFun, directed=FALSE, nperm=1000, probs=0.95, SRI=FALSE, effort=1){

  net.original <- create.a.network(data, directed,SRI=SRI, effort=effort)
  Perm.measure<-vector()

  for(i in 1:nperm){

    #permute edges
    net.per <- net.original
    E(net.per)$weight <- sample(E(net.original)$weight, replace = F)

    # Get measure
    Perm.measure[length(Perm.measure)+1]<- measureFun(net.per)

  }

  probs.left<-1-probs
  return(quantile(Perm.measure, probs = c( (0+probs.left/2), (1-probs.left/2) ), na.rm=T))

}


#' Perform edge permutations maintaining the degree distribution of original graph.
#'
#' This function will permute a network by randomly moving edges (see the keeping_degseq function in igraph).
#' @param data Dataframe with relational data in the first two rows, and a time stamp in the third row. Note: time stamps should be in ymd or ymd_hms format. The lubridate package can be very helpful in organizing times.
#' @param measureFun This is a function that takes as an input a igraph network and returns a single value.
#' @param directed Whether to consider the network as directed or not (TRUE/FALSE).
#' @param nperm Number of permutations to perform before extracting network measures.
#' @param probs numeric vector of probabilities with values in [0,1].
#' @importFrom igraph rewire
#' @export
#'
#'
perm.edge.degseq <- function(data, measureFun, directed=FALSE, nperm=1000, probs=0.95, SRI=FALSE, effort=1){

  net.original <- create.a.network(data, directed,SRI=SRI, effort=effort)
  Perm.measure<-vector()

  for(i in 1:nperm){

    #permute edges
    net.per <- rewire(net.original, with=keeping_degseq(niter = vcount(net.original) * 10) )

    # Get measure
    Perm.measure[length(Perm.measure)+1]<- measureFun(net.per)

  }

  probs.left<-1-probs
  return(quantile(Perm.measure, probs = c( (0+probs.left/2), (1-probs.left/2) ), na.rm=T))

}



