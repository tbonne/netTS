#' Simulate some observation data
#'
#' This function will simulate events data.
#' @param nodes Number of individuals within a simulated group.
#' @param sampling.periods The number of times the simulated group is observed.
#' @param sampling.periods.per.day The number of sampling perids per day.
#' @param true.net (Optional) A true underlying network describing the probability of each individual interacting.
#' @param ind.probs (Optional) A vector specifying the probability of observing the individual performing the behaivour. Should be the same length as the number of nodes.
#' @param ind.sd (Optional) A value for the standard deviation around observed probability of observing a behaviour. Used with the cor.mat option to determine the rate of change in probability of behaviour.
#' @param cor.mat (Optional) A correlation matrix (size nodes x nodes) describing the dependence between individual changes in behaviour.
#' @importFrom igraph make_full_graph
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom stats runif
#' @importFrom corpcor make.positive.definite is.positive.definite
#' @importFrom MASS mvrnorm
#' @export
#'
sim.events.data <- function(nodes, sampling.periods, sampling.periods.per.day=1, true.net=NULL, ind.probs=NULL, e.sd=NULL, cor.mat=NULL, covariates=NULL, covariates.beta=NULL){

  ####create time series of behaivour probabilities

  #mean values
  A <- matrix(runif(nodes), nrow=nodes, ncol=1)
  if(!is.null(ind.probs) ) A <- matrix(ind.probs, nrow=nodes, ncol=1)

  #correlation matrix
  B <- diag(nodes)
  if(!is.null(cor.mat))B<-cor.mat

  #error/change rate
  ind.sd=1
  if(!is.null(ind.sd))ind.sd=e.sd

  #covariate estiamtes
  C <- matrix(0, nrow=nodes, ncol=1)
  if(!is.null(covariates.beta)) C <- matrix(covariates.beta,nrow=nodes)

  #setup for sim
  X <- matrix(NA, nrow=sampling.periods, ncol=nodes, dimnames=list(NULL, as.character(1:nodes)))
  U <- matrix(0, nrow=sampling.periods, ncol=1)
  if(!is.null(covariates)) U <- covariates
  E <- matrix(rnorm(sampling.periods*nodes) * ind.sd, nrow=sampling.periods, ncol=nodes)

  #initial values
  X[1,] <- runif(nodes,0,1)

  #simulate the probability time series
  for(i in 2:sampling.periods){
    p_scale<-A + B%*%matrix(X[i-1,],ncol=1) + C%*%matrix(U[i,],ncol=1) + E[i,]
    X[i,] <- exp(p_scale)/(1+exp(p_scale))
  }

  #setup dataframe to capture simulated data
  day=lubridate::ymd_hms("2002/07/24 00:00:00")
  df.sim <- data.frame(from = 1, to=2, date=day, sampleID=1 )

  #monitor the progress
  pb <- txtProgressBar(min = 1, max = sampling.periods, style = 3)

  for(i in 1:sampling.periods){

    #update probability each individual will show the behaviour
    ind.behav <- X[i,]

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

