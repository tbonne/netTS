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

#extract measures test
test.m <-extract_measure_nodes(dir.nets, strength,unique(c(groomEvents$ID,groomEvents$PartnerID)) )
test.m2 <- extract_lagged_measure_nodes(dir.nets, cosine_between_nodes,lag=1,unique.names=unique(c(groomEvents$ID,as.character(groomEvents$PartnerID) )))

out_degree <- function(x){
  degree(x, mode="out")
}
output.net.dir.out<-nodeTS(groomEvents,windowsize = days(30), windowshift = days(30), measureFun = out_degree, directed=TRUE)

df_enter_leave <- data.frame(ID = unique(groomEvents$ID), startDate = min(groomEvents$date), endDate = min(groomEvents$date)  )
df_enter_leave$startDate <- ymd_hms(df_enter_leave$startDate)
df_enter_leave$endDate <- ymd_hms(df_enter_leave$endDate)
trimed_df <- trim_nodes(output.net.dir, groomEvents, enter_leave = df_enter_leave)

#node measure tests
out.edgeS<-edge.weight.skewness(fixed.net.undir)
out.cosNo<-cosine_between_nodes(fixed.net.dir,fixed.net.dir)

test_that("subsetted network is equal", {
  expect_equal(vcount(dir.nets[[1]]), vcount(fixed.net.dir) )
  expect_equal(vcount(undir.nets[[1]]), vcount(fixed.net.undir) )
  expect_equal(output.net.dir[1,"Malc"], as.vector(degree(fixed.net.dir)["Malc"]) )
  expect_equal(output.net.undir[1,"Malc"], as.vector(degree(fixed.net.undir)["Malc"]) )
  expect_equal(output.net.dir.out[1,"Laur"], as.vector(degree(fixed.net.dir, mode="out")["Laur"]) )
  expect_equal(test.m[1,1], 17)
  expect_equal(test.m2[3,1], 0.5728558)
  expect_equal(trimed_df[1,1], NA)
  expect_equal(out.cosNo[1,1], 1)
  expect_equal(out.edgeS[2],NaN)
})
