#' Simulate some observation data
#'
#' This function will simulate events data.
#' @param nodes Number of individuals within a simulated group.
#' @param sampling.periods The number of times the simulated group is observed.
#' @param sampling.periods.per.day The number of sampling perids per day.
#' @param true.net (Optional) A true underlying network describing the probability of each individual interacting.
#' @importFrom igraph make_full_graph
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom stats runif
#' @importFrom corpcor make.positive.definite
#' @importFrom MASS mvrnorm
#' @export
#'
sim.events.data <- function(nodes, sampling.periods, sampling.periods.per.day=1, true.net=NULL, ind.probs=NULL, cor.mat=NULL, ind.sd=NULL){

  #setup dataframe to capture simulated data
  day=lubridate::ymd_hms("2002/07/24 00:00:00")
  df.sim <- data.frame(from = 1, to=2, date=day, sampleID=1 )

  #convert correlation matrix to covariance matrix
  if( !is.null(cor.mat) & !is.null(ind.sd) ){
    change_vec <- rep(ind.sd,nrow(cor.mat))
    b <- change_vec %*% t(change_vec)
    a_covariance <- b * cor.mat
    a_covariance <- make.positive.definite(a_covariance)
    print("correlation matrix used for updating node behaivour")
  }

  #Set probability each individual will show the behaviour
  if(is.null(ind.probs)){
    ind.behav <- runif(nodes) #random
    print("Uniform random draws used for updating node behaivour")
  } else{
    ind.behav <- ind.probs #user set
    if(is.null(cor.mat) | is.null(ind.sd)){
      print("Fixed node behaivour")
    }
  }

  #monitor the progress
  pb <- txtProgressBar(min = 1, max = sampling.periods, style = 3)

  for(i in 1:sampling.periods){

    #update probability each individual will show the behaviour
    if(is.null(ind.probs)){
      ind.behav <- runif(nodes)
    } else if (!is.null(cor.mat) & !is.null(ind.sd)){
      ind.behav <- ind.behav +  mvrnorm(1,rep(0,length(ind.behav)),a_covariance)
    }

    #Each individual perform the behaviour
    for(j in 1:length(ind.behav)){
      if(runif(1)<ind.behav[j]){

        #who with?
        if(is.null(true.net)){
          g.random <- igraph::make_full_graph(nodes)
          E(g.random)$weight <- runif(ecount(g.random))
          m.random <- igraph::get.adjacency(g.random, attr="weight", sparse = T)
          chosen.node<-sample(1:nodes,1,prob = m.random[,j])
          if(j!=chosen.node)df.sim<-rbind(df.sim,data.frame(from = j, to=chosen.node, date=day,sampleID=i ))

        } else {
          m.random <- igraph::get.adjacency(true.net, attr="weight", sparse = T)
          chosen.node<-sample(1:nodes,1,prob = m.random[,j])
          df.sim<-rbind(df.sim,data.frame(from = j, to=chosen.node, date=day,sampleID=i ))
        }
      }
    }

    #add day
    if( i %% sampling.periods.per.day == 0){
      day = day+lubridate::days(1)
    }

    #update progress bar
    setTxtProgressBar(pb,  i )

  }

  close(pb)

  return(df.sim)
}

