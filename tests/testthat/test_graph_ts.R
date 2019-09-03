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

mean_out_degree <- function(x){
  mean(degree(x, mode="out"))
}
output.net.dir.out<-graphTS(groomEvents,windowsize = days(30), windowshift = days(30), measureFun = mean_out_degree, directed=TRUE)


#second test set:
data<- data.frame ( from=c("laur", "mori", "daff", "lore", "elto"),
                    to = c("mori", "daff", "elto", "daff", "daff"),
                    date= c("2015-01-01", "2015-01-10", "2015-01-15","2015-01-20","2015-01-30"))
data$from<- as.character(data$from)
data$to<- as.character(data$to)
data$date<- ymd(data$date)
graph<- create.a.network(data)
values<- graphTS(data, windowsize = days(10), windowshift= days(5), measureFun=degree_mean)
netlist<-extract_networks(data, windowsize = days(10), windowshift= days(5), directed=FALSE)
net.measure<-extract_measure_network(netlist, measureFun=degree_mean)


#third test set (parallel code):
para.netlist <- extract_networks_para(data, windowsize= days(10), windowshift= days(5), directed = FALSE, cores=2, SRI=FALSE, effortFun=NULL)

## check net.para function
# need window.range first
#generate a list of windows times
# windowStart=as_date("2019-01-01")
# windowShift = days(round(365/4))
# windowEnd=windowStart+windowShift
# window.ranges <- data.frame(start=windowStart, end=windowEnd)
# endDay=as_date("2019-12-01")
# while(windowEnd<=endDay){
#   window.ranges <-  rbind(window.ranges,data.frame(start=windowStart, end=windowStart+windowShift))
#   windowStart = windowStart + windowShift
#   windowEnd = windowStart + windowShift
# }


test_that("subsetted network is equal", {
  expect_equal(vcount(dir.nets[[1]]), vcount(fixed.net.dir) )
  expect_equal(vcount(undir.nets[[1]]), vcount(fixed.net.undir) )
  expect_equal(output.net.dir[1,1], mean(degree(fixed.net.dir)) )
  expect_equal(output.net.undir[1,1], mean(degree(fixed.net.undir)) )
  expect_equal(output.net.dir.out[1,1], mean(degree(fixed.net.dir, mode="out")) )
  expect_equal(values[1,1], 4/3)
  expect_equal(length(netlist), 4)
  expect_equal(length(para.netlist), 4)
  expect_equal(net.measure[1,1], 4/3)#check measure
  expect_equal(net.measure[1,2], 2) #check nb of extracted EVENTS
})
