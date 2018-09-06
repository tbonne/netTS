context("dyadTS tests")
library(netTS)
library(lubridate)
library(igraph)

#setup one network
subdata <- groomEvents[groomEvents$date<(min(groomEvents$date)+days(30)),]
fixed.net.dir <- create.a.network(subdata, directed = TRUE)
fixed.net.undir <- create.a.network(subdata, directed = FALSE)


###############directed
output.net.dir  <- dyadTS(groomEvents,windowsize = days(30), windowshift = days(30), measureFun = dyad_weight, directed=TRUE)

#get an edge from the known graph
dyad.test<-ends(fixed.net.dir, E(fixed.net.dir)[1])
dyad.test<-paste0(dyad.test[1],"_",dyad.test[2])
dyad.weight<-edge_attr(fixed.net.dir, name = "weight")[1]

#get same dyad from the dyad.ts output
dyad.out.weight<-output.net.dir[,colnames(output.net.dir)==dyad.test][1]


#################undirected
output.net.undir<- dyadTS(groomEvents,windowsize = days(30), windowshift = days(30), measureFun = dyad_weight, directed=FALSE)

#get an edge from the known graph
dyad.test.undir<-ends(fixed.net.undir, E(fixed.net.undir)[1])
dyad.test.undir<-paste0(dyad.test.undir[1],"_",dyad.test.undir[2])
dyad.weight.undir<-edge_attr(fixed.net.undir, name = "weight")[1]

#get same dyad from the dyad.ts output
dyad.out.weight.undir<-output.net.undir[,colnames(output.net.undir)==dyad.test.undir][1]


#####################################################

test_that("subsetted network is equal", {
  expect_equal(dyad.out.weight, dyad.weight )
  expect_equal(dyad.out.weight.undir, dyad.weight.undir )
})
