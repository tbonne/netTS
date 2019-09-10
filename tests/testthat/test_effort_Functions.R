context("Effort function tests")
library(testthat)


## Measurment function effort
df.effort <-data.frame(from=c("daff","damo","lore","elto","damo","gizm", "daff","lore"),
                       to= c("gizm", "ella", "ella", "ella","daff","daff","ella", "elto"),
                       date=c("2017-07-09 08:00:00", "2017-07-09 08:08:00","2017-07-09 08:10:00","2017-07-09 16:00:00","2017-07-09  15:10:00","2017-07-09 16:11:00","2017-07-09  18:00:00","2017-07-09  09:00:00"),
                       scanID = c(1,1,1,2,2,2,3,4))
df.effort$date<-ymd_hms(df.effort$date)

df.effort.scans <- data.frame(date=c("2017-07-09 08:00:00", "2017-07-09 08:08:00","2017-07-09 08:10:00","2017-07-09 16:00:00","2017-07-09  15:10:00","2017-07-09 16:11:00","2017-07-09  18:00:00","2017-07-09  09:00:00"),
                              scanID= c(1,1,1,2,2,2,3,4))

##focal data
df.focal <- data.frame(date=c("2017-07-09 08:00:00", "2017-07-09 08:08:00","2017-07-09 08:10:00","2017-07-09 16:00:00","2017-07-09  15:10:00","2017-07-09 16:11:00","2017-07-09  18:00:00","2017-07-09  09:00:00"),
                       ID = c("daff","damo","lore","elto","damo","gizm", "daff","lore"),
                       duration = c(10,10,10,10,10,10,10,10))


### controling for entering and leaving

#calculate when individuals were first and last seen in a network
df.inOut<-node_first_last(df.effort)

#use these to get weighted measures
net.list <- extract_networks(df.effort, windowsize = hours(1), windowshift = hours(1))
measures.deg <- degree(net.list[[2]])
mean.corrected.w <- weighted_mean(measures.deg, df.inOut, net.list[[2]])
graph_attr(net.list[[2]], "windowstart")
graph_attr(net.list[[2]], "windowend")

test_that("effort.scan is good I guess", {
  expect_equal(as.numeric(graph_attr(effort.scan(df.window=df.effort, df.scans = df.effort.scans),"effort")), 16)
  expect_equal(as.numeric(graph_attr(effort.time(df.window=df.effort),"effort")), 10)
  expect_equal(as.numeric(E(effort.focal(df.window=df.effort,effortData=df.focal))$weight[1]), 1/40)
  expect_equal(df.inOut[1,3],ymd_hms("2017-07-09 18:00:00"))
  expect_equal(mean.corrected.w, (0*1+1*1)/(2)  )
  #expect_equal(eff_time[1,1], 10/60  )
})


