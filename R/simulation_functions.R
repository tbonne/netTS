#' Simulate some observation data
#'
#' This function simulates interaction data based on a group scaning methodology. For a given number of sampling periods, the simulation will generate interactions based on a user defined network, and individual interaction probabilities.
#' @param nodes Number of individuals within a simulated group.
#' @param sampling.periods The number of times the simulated group is observed.
#' @param sampling.periods.per.day The number of sampling periods per day.
#' @param true.net (Optional) A network describing who can interact with whom.
#' @param A (Optional) A vector specifying the average probability interaction for each individual. Should be the same length as the number of nodes. Default: runif(node,0,1)
#' @param B (Optional) A matrix (node x node) specifying the inter-dependence of individual interaction probabilities. Default: matrix of zeros
#' @param C (Optional) A matrix (node X covariates) specifying the effect of covariate values on each individual's interaction probability. Default: matrix of zeros
#' @param D (Optional) A matrix (node x node) specifying the inter-dependence of noise around individual variations probabilities. Default: matrix of zeros
#' @param E (Optional) A value for the standard deviation (noise) around individual's mean interaction probability. Default: 0.1
#' @param covariates A matrix (sampling period x covariates) with covariate values for each sampling period.
#' @importFrom igraph make_full_graph
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom stats runif
#' @importFrom corpcor make.positive.definite is.positive.definite
#' @importFrom MASS mvrnorm
#' @export
#'
sim.events.data <- function(nodes, sampling.periods, sampling.periods.per.day=1, true.net=NULL, A=NULL, B=NULL, C=NULL, D=NULL, E=NULL, covariates=NULL){

  #initialize individual probability of behaviour
  if(is.null(A))A<-matrix(runif(nodes,0,1), nrow=nodes, ncol=1)
  A.logit <- log(A/(1-A))

  #correlation matrix: means
  if(is.null(B))B<-diag(nodes)*0

  #correlation matrix: error
  if(is.null(D))D<-diag(nodes)*0.0

  #error/change rate
  if(is.null(E))E=0.1

  #covariate effects on individual probabilities
  if(is.null(C)) C <- matrix(0, nrow=nodes, ncol=1)

  #covariate values
  if(is.null(covariates)) covariates <- matrix(0, nrow=sampling.periods, ncol=1)

  #probabilties
  X <- matrix(NA, nrow=sampling.periods, ncol=nodes, dimnames=list(NULL, as.character(1:nodes)))

  #setup dataframe to capture simulated data
  day=lubridate::ymd_hms("2002/07/24 00:00:00")
  df.sim <- data.frame(from = 1, to=2, date=day, sampleID=1 )

  #monitor the progress
  pb <- txtProgressBar(min = 1, max = sampling.periods, style = 3)

  #individuals start at mean probability
  ind.behav <- A

  #Start simulations
  for(i in 1:sampling.periods){

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

    ###update probability each individual will show the behaviour
    #convert probability to logit scale
    ind.behav<-ifelse(ind.behav<0.001,0.001,ifelse(ind.behav>0.999,0.999,ind.behav))
    p.logit<-log(ind.behav/(1-ind.behav))
    #adjust probability (AR model)
    p.logit.mean <- A.logit + B%*%matrix(p.logit,ncol=1) + C%*%matrix(covariates[i,],ncol=1) + D%*%matrix(p.logit-A.logit,ncol=1) + matrix(mvrnorm(n = 1, rep(0,nodes), Sigma=diag(nodes)),ncol=1)*E

    #convert back to probability
    p.mean <- exp(p.logit.mean)/(exp(p.logit.mean)+1)
    ind.behav <- p.mean
    X[i,] <- p.mean

    #update progress bar
    setTxtProgressBar(pb,  i )

  }

  close(pb)

  return(list(beha=df.sim,probs=X))
}


