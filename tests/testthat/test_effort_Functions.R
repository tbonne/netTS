context("Effort function tests")
library(testthat)


## Measurment function effort
df.effort <-data.frame(from=c("daff","damo","lore","elto","damo","gizm", "daff","lore"),
                       to= c("gizm", "ella", "ella", "ella","daff","daff","ella", "elto"),
                       date=c("2017-07-09 08:00:00", "2017-07-09 08:08:00","2017-07-09 08:10:00","2017-07-09 16:00:00","2017-07-09  15:10:00","2017-07-09 16:11:00","2017-07-09  18:00:00","2017-07-09  09:00:00"),
                       sampleID = c(1,1,1,2,2,2,3,4))
df.effort$date<-ymd_hms(df.effort$date)

##focal data
df.focal <- data.frame(date=c("2017-07-09 08:00:00", "2017-07-09 08:08:00","2017-07-09 08:10:00","2017-07-09 16:00:00","2017-07-09  15:10:00","2017-07-09 16:11:00","2017-07-09  18:00:00","2017-07-09  09:00:00"),
                       ID = c("daff","damo","lore","elto","damo","gizm", "daff","lore"),
                       duration = c(10,10,10,10,10,10,10,10))


### controling for entering and leaving

## Check weighted mean: test it via the extract measure netwrok function
#here elto shows up later in network

firstLast<- data.frame( ID = c( "laur", "mori", "daff", "lore", "elto"),
                        In= c("2015-01-01", "2015-01-01", "2015-01-01,","2015-01-01","2015-01-01"),
                        Out= c("2015-01-06", "2015-01-30", "2015-01-30", "2015-01-30", "2015-01-30"))

firstLast$ID<- as.character(firstLast$ID)
firstLast$In<- lubridate::ymd(firstLast$In)
firstLast$Out<- lubridate::ymd(firstLast$Out)
netlist<-extract_networks(firstLast, windowsize = days(10), windowshift= days(5), directed=FALSE)
obs.values<- degree(netlist[[1]])





test_that("effort.scan is good I guess", {
  expect_equal(as.numeric(graph_attr(effort.scan(df.window=df.effort),"effort")), 4)
  expect_equal(as.numeric(graph_attr(effort.time(df.window=df.effort),"effort")), 10)
  expect_equal(as.numeric(E(effort.focal(df.window=df.effort,effortData=df.focal))$weight[1]), 1/40)
  skip('skip')
  expect_equal(weighted_mean(netlist, inOut = firstLast, netlist[[1]]), 3.5/3)
  expect_equal(node_first_last(data)[1,2],as.Date("2015-01-01"))
})

