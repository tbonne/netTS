context("nodeTS tests")
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
output.net.dir<-nodeTS(groomEvents,windowsize = days(30), windowshift = days(30), measureFun = degree, directed=TRUE)

output.net.undir<-nodeTS(groomEvents,windowsize = days(30), windowshift = days(30), measureFun = degree, directed=FALSE)

out_degree <- function(x){
  degree(x, mode="out")
}
output.net.dir.out<-nodeTS(groomEvents,windowsize = days(30), windowshift = days(30), measureFun = out_degree, directed=TRUE)


test_that("subsetted network is equal", {
  expect_equal(vcount(dir.nets[[1]]), vcount(fixed.net.dir) )
  expect_equal(vcount(undir.nets[[1]]), vcount(fixed.net.undir) )
  expect_equal(output.net.dir[1,"Malc"], as.vector(degree(fixed.net.dir)["Malc"]) )
  expect_equal(output.net.undir[1,"Malc"], as.vector(degree(fixed.net.undir)["Malc"]) )
  expect_equal(output.net.dir.out[1,"Laur"], as.vector(degree(fixed.net.dir, mode="out")["Laur"]) )
})
