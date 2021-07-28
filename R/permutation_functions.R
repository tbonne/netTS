#' Perform permutation on the events dataframe
#'
#' This function will permute a network by randomly switching individuals between events.
#' @param data Dataframe with relational data in the first two rows, and a time stamp in the third row. Note: time stamps should be in ymd or ymd_hms format. The lubridate package can be very helpful in organizing times.
#' @param measureFun This is a function that takes as an input a igraph network and returns a single value.
#' @param directed Whether to consider the network as directed or not (TRUE/FALSE).
#' @param nperm Number of permutations to perform before extracting network measures.
#' @param probs numeric vector of probabilities with values in [0,1].
#' @param SRI Wether to use the simple ratio index (Default=FALSE).
#' @param effortFun This is a function that takes as input the data within a window of time and returns the total sampling effort.
#' @param effortData This is a dataframe containing the data used to calculate sampling effort.
#' @importFrom stats quantile
#' @export
#'
#'
perm.events <- function(data, measureFun, directed=FALSE, nperm=1000, probs=0.95, SRI=FALSE, effortFun=NULL, effortData=NULL, windowstart, windowend){

  #vector to store values
  Perm.measure<-vector()

  #force columns to be character
  data[,1]<- as.character(data[,1])
  data[,2]<- as.character(data[,2])

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
    if(is.null(effortFun)==FALSE & is.null(effortData)==TRUE  ){ #there is an effort function and it requires no external data
      Perm.network = effortFun(data, directed = directed)

    }else if(is.null(effortFun)==FALSE & is.null(effortData)==FALSE ){ #there is an effort function and it requires some external data
      effortData.sub <- effortData[effortData[,1]>=windowstart & effortData[,1]<windowend,]
      Perm.network = effortFun(data, effortData.sub, directed = directed)

    } else { #there is no effort function
      Perm.network <- create.a.network(data, directed = directed, SRI = SRI)
    }


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
#' @param SRI Wether to use the simple ratio index (Default=FALSE).
#' @param effortFun This is a function that takes as input the data within a window of time and returns the total sampling effort.
#' @param effortData This is a dataframe containing the data used to calculate sampling effort.
#' @importFrom stats quantile runif
#' @export
#'
#'
perm.events.directed <- function(data, measureFun, directed=FALSE, nperm=1000, probs=0.95, SRI=FALSE, effortFun=NULL, effortData=NULL, windowstart, windowend){

  #vector to store values
  Perm.measure<-vector()

  #force columns to be character
  data[,1]<- as.character(data[,1])
  data[,2]<- as.character(data[,2])

  for(i in 1:nperm){

    no.loops= FALSE

    #choose to or from grooming to permute
    if(0.5 > runif(1)){

      while(no.loops == FALSE){

        #choose two individuals to switch
        rows.to.switch <- sample(1:nrow(data),2,F)

        #record old order
        old.order <- data[,2]
        new.order <- data[,2]

        #update order
        new.order[rows.to.switch[1]] <- old.order[rows.to.switch[2]]
        new.order[rows.to.switch[2]] <- old.order[rows.to.switch[1]]

        #check to make sure there are no self loops
        if(sum(as.character(data[,1])==as.character(new.order) )==0){ #$from

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
        old.order <- data[,1]
        new.order <- data[,1]

        #update order
        new.order[rows.to.switch[1]] <- old.order[rows.to.switch[2]]
        new.order[rows.to.switch[2]] <- old.order[rows.to.switch[1]]

        #check to make sure there are no self loops
        if(sum(as.character(new.order)==as.character(data[,2]) )==0){ #$to

          data$from <- new.order
          NewData<- data
          no.loops=TRUE

        }
      }
    }

    #Create graph in order to get the measure
    if(is.null(effortFun)==FALSE & is.null(effortData)==TRUE  ){ #there is an effort function and it requires no external data
      Perm.network = effortFun(data, directed = directed)

    }else if(is.null(effortFun)==FALSE & is.null(effortData)==FALSE ){ #there is an effort function and it requires some external data
      effortData.sub <- effortData[effortData[,1]>=windowstart & effortData[,1]<windowend,]
      Perm.network = effortFun(data, effortData.sub, directed = directed)

    } else { #there is no effort function
      Perm.network <- create.a.network(data, directed = directed, SRI = SRI)
    }

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
#' @param SRI Wether to use the simple ratio index (Default=FALSE).
#' @param effortFun This is a function that takes as input the data within a window of time and returns the total sampling effort.
#' @param effortData This is a dataframe containing the data used to calculate sampling effort.
#' @importFrom igraph E
#' @importFrom stats quantile
#' @export
#'
#'
perm.edge.weights <- function(data, measureFun, directed=FALSE, nperm=1000, probs=0.95, SRI=FALSE, effortFun=NULL,effortData=NULL, windowstart, windowend){

  #Create graph in order to get the measure
  if(is.null(effortFun)==FALSE & is.null(effortData)==TRUE  ){ #there is an effort function and it requires no external data
    net.original = effortFun(data, directed = directed)

  }else if(is.null(effortFun)==FALSE & is.null(effortData)==FALSE ){ #there is an effort function and it requires some external data
    effortData.sub <- effortData[effortData[,1]>=windowstart & effortData[,1]<windowend,]
    net.original = effortFun(data, effortData.sub, directed = directed)

  } else { #there is no effort function
    net.original <- create.a.network(data, directed = directed, SRI = SRI)
  }

  Perm.measure<-vector()

  for(i in 1:nperm){

    #permute edges
    net.per <- net.original
    E(net.per)$weight <- sample(E(net.original)$weight, replace = F)

    # Get measure
    Perm.measure[length(Perm.measure)+1]<- measureFun(net.per)

  }

  probs.left<-1-probs
  return(quantile(Perm.measure, probs = c( (0+probs.left/2), (1-probs.left/2) ) , na.rm=T))

}


#' Perform edge permutations maintaining the degree distribution of original graph.
#'
#' This function will permute a network by randomly moving edges (see the keeping_degseq function in igraph).
#' @param data Dataframe with relational data in the first two rows, and a time stamp in the third row. Note: time stamps should be in ymd or ymd_hms format. The lubridate package can be very helpful in organizing times.
#' @param measureFun This is a function that takes as an input a igraph network and returns a single value.
#' @param directed Whether to consider the network as directed or not (TRUE/FALSE).
#' @param nperm Number of permutations to perform before extracting network measures.
#' @param probs numeric vector of probabilities with values in [0,1].
#' @param SRI Wether to use the simple ratio index (Default=FALSE).
#' @param effortFun This is a function that takes as input the data within a window of time and returns the total sampling effort.
#' @param effortData This is a dataframe containing the data used to calculate sampling effort.
#' @importFrom igraph rewire
#' @importFrom stats quantile
#' @export
#'
#'
perm.edge.degseq <- function(data, measureFun, directed=FALSE, nperm=1000, probs=0.95, SRI=FALSE, effortFun=NULL,effortData=NULL, windowstart, windowend){

  #Create graph in order to get the measure
  if(is.null(effortFun)==FALSE & is.null(effortData)==TRUE  ){ #there is an effort function and it requires no external data
    net.original = effortFun(data, directed = directed)

  }else if(is.null(effortFun)==FALSE & is.null(effortData)==FALSE ){ #there is an effort function and it requires some external data
    effortData.sub <- effortData[effortData[,1]>=windowstart & effortData[,1]<windowend,]
    net.original = effortFun(data, effortData.sub, directed = directed)

  } else { #there is no effort function
    net.original <- create.a.network(data, directed = directed, SRI = SRI)
  }

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


#' Perform permutation on the events dataframe
#'
#' This function will permute a network by randomly switching individuals between events.
#' @param data Dataframe with relational data in the first two rows, and a time stamp in the third row. Note: time stamps should be in ymd or ymd_hms format. The lubridate package can be very helpful in organizing times.
#' @param measureFun This is a function that takes as an input a igraph network and returns a single value.
#' @param directed Whether to consider the network as directed or not (TRUE/FALSE).
#' @param nperm Number of permutations to perform before extracting network measures.
#' @param probs numeric vector of probabilities with values in [0,1].
#' @param SRI Wether to use the simple ratio index (Default=FALSE).
#' @param effortFun This is a function that takes as input the data within a window of time and returns the total sampling effort.
#' @param effortData This is a dataframe containing the data used to calculate sampling effort.
#' @importFrom stats quantile
#' @export
#'
#'
perm.events.multiple.outputs <- function(data, measureFun, directed=FALSE, nperm=1000, probs=0.95, SRI=FALSE, effortFun=NULL, effortData=NULL, windowstart, windowend){

  #vector to store values
  Perm.measure<-data.frame()

  #force columns to be character
  data[,1]<- as.character(data[,1])
  data[,2]<- as.character(data[,2])

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
    if(is.null(effortFun)==FALSE & is.null(effortData)==TRUE  ){ #there is an effort function and it requires no external data
      Perm.network = effortFun(data, directed = directed)

    }else if(is.null(effortFun)==FALSE & is.null(effortData)==FALSE ){ #there is an effort function and it requires some external data
      effortData.sub <- effortData[effortData[,1]>=windowstart & effortData[,1]<windowend,]
      Perm.network = effortFun(data, effortData.sub, directed = directed)

    } else { #there is no effort function
      Perm.network <- create.a.network(data, directed = directed, SRI = SRI)
    }


    # Get measure
    Perm.measure <- bind_rows(Perm.measure, measureFun(Perm.network) )
    #i'm here right now ... need to figure out the best way to store these multiple measures...
    #measureFun = dyad_diff
  }


  probs.left<-1-probs

  #get the quantiles for each column
  df.quant<-sapply(Perm.measure, FUN=quantile, probs=c( (0+probs.left/2), (1-probs.left/2) ), na.rm=T)

  return(df.quant)
}



#' Use permutation to extract uncertainty
#'
#' This function will estimate network measures given random permutations on the original data.
#' @param data A dataframe with relational data in the first two columns, and a time stamp in the third column. An optional column with a weight can be added if there is a duration or magnitude for each interaction (column name for this should be set to 'weight'). Note: time stamps should be in ymd or ymd_hms format. The lubridate package can be very helpful in organizing times.
#' @param windowsize The size of each window in which to generate a network.
#' @param windowshift The amount of time to shift the window when generating networks.
#' @param directed Whether to consider the network as directed or not (TRUE/FALSE).
#' @param measureFun This is a function that takes as an input a igraph network and returns a single value.
#' @param effortFun This is a function that takes as input the data within a window of time and returns the total sampling effort.
#' @param effortData This is a dataframe containing the data used to calculate sampling effort. The first column should contain timedate values.
#' @param probs numeric vector of probabilities with values in [0,1].
#' @param nperm Number of permutations to perform before extracting network measures.
#' @param SRI Wether to use the simple ratio index (Default=FALSE).
#' @param graphlist A list of networks to perfrom permutations on.
#' @param permutationFun A function that will be used to perform permutations on the event data (i.e., before the network) or on the network itself. See vignette: Using_permutations.
#' @importFrom R.utils doCall
#' @export
#'
#'
permutation.multi.values<-function(data, windowsize, windowshift, directed = FALSE,measureFun, effortFun=NULL, effortData=NULL, probs=0.95, nperm=1000, SRI=FALSE,graphlist=NULL, permutationFun=perm.events.multiple.outputs, unique.names = all.unique.names ){

  print("perm")

  #dataframe to record permutation results
  perm.values.high <- data.frame(matrix(ncol = length(unique.names), nrow = 0))
  perm.values.low <- data.frame(matrix(ncol = length(unique.names), nrow = 0))
  colnames(perm.values.high) <- unique.names
  colnames(perm.values.low) <- unique.names

  #monitor the progress
  #pb <- txtProgressBar(min = 1, max = length(graphlist), style = 3)

  #run the permutation on each network
  for(i in 1:length(graphlist)){

    #get the time bounds of each network
    windowstart <- graph_attr(graphlist[[i]],"windowstart")
    windowend   <- graph_attr(graphlist[[i]],"windowend")

    #get the data within this time range
    df.window<-create.window(data, windowstart, windowend)
    Observation.Events <- nrow(df.window)

    if(Observation.Events >= 2 ){

      #perform permutations
      perm.out<-doCall(permutationFun,data=df.window, measureFun=measureFun, directed=directed, windowstart=windowstart, windowend=windowend, probs=probs,nperm= nperm, SRI=SRI, effortFun=NULL, effortData=NULL)

      #record the high and low estimates
      df.temp <- data.frame(matrix(ncol = length(unique.names), nrow = 0))
      colnames(df.temp) <- unique.names
      df.temp <- bind_rows(df.temp,as.data.frame(perm.out) )
      perm.values.low<- bind_rows(perm.values.low,df.temp[1,] )
      perm.values.high<- bind_rows(perm.values.high,df.temp[2,] )

    } else {

      perm.values.low.na <- data.frame(matrix(rep(NA,length(unique.names)),ncol = length(unique.names), nrow = 1))
      perm.values.high.na <- data.frame(matrix(rep(NA,length(unique.names)),ncol = length(unique.names), nrow = 1))
      colnames(perm.values.low.na) <- unique.names
      colnames(perm.values.high.na) <- unique.names

      perm.values.low<- bind_rows(perm.values.low,perm.values.low.na )
      perm.values.high<- bind_rows(perm.values.high,perm.values.high.na )

    }

    #update progress bar
    #setTxtProgressBar(pb,  as.numeric(windowend) )

  }

  #close progress bar
  #close(pb)

  #put the dataframes into a list
  perm.list<-list(low=perm.values.low,high=perm.values.high)

  return(perm.list)

}

