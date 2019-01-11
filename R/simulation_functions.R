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
#' @export
#'
sim.events.data <- function(nodes, sampling.periods, sampling.periods.per.day=1, true.net=NULL){

  day=lubridate::ymd_hms("2002/07/24 00:00:00")
  df.sim <- data.frame(from = 1, to=2, date=day, sampleID=1 )

  #monitor the progress
  pb <- txtProgressBar(min = 1, max = sampling.periods, style = 3)

  for(i in 1:sampling.periods){

    #probability each individual will show the behaviour
    ind.behav <- runif(nodes)

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
