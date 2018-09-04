context("graphTS tests")
library(netTS)
library(lubridate)
library(igraph)

#setup one network
subdata <- groomEvents[groomEvents$date<(min(groomEvents$date)+days(30)),]
fixed.net.dir <- create.a.network(subdata, directed = TRUE)
fixed.net.undir <- create.a.network(subdata, directed = FALSE)

#test if the correct network is extracted
dir.nets <- extract_networks(groomEvents,windowsize = days(30),windowshift = days(30),directed = TRUE)
undir.nets <- extract_networks(groomEvents,windowsize = days(30),windowshift = days(30),directed = FALSE)

#test if graphTS code gets the same measures
output.net.dir<-graphTS(groomEvents,windowsize = days(30), windowshift = days(30), measureFun = degree_mean, directed=TRUE)
output.net.undir<-graphTS(groomEvents,windowsize = days(30), windowshift = days(30), measureFun = degree_mean, directed=FALSE)



test_that("subsetted network is equal", {
  expect_equal(vcount(dir.nets[[1]]), vcount(fixed.net.dir) )
  expect_equal(vcount(undir.nets[[1]]), vcount(fixed.net.undir) )
  expect_equal(output.net.dir[1,1], mean(degree(fixed.net.dir)) )
  expect_equal(output.net.undir[1,1], mean(degree(fixed.net.undir)) )
})
