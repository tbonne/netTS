context("Measurment tests")
library(testthat)

## Test: measurment functions


#generate data
df1 <-data.frame(from=c("daff","damo","lore","elto","damo","gizm", "daff","lore"),
                 to= c("gizm", "ella", "ella", "ella","daff","daff","ella", "elto"),
                 weight= c(2,4,6,8,10,1,1,1),
                 date=c("2017-07-09", "2017-07-09","2017-07-09","2017-07-09","2017-07-09","2017-07-09","2017-07-09","2017-07-09"))

df <-data.frame(from=c("wall","mori","panc","nige","schm","bone", "raje","ring"),
                to= c("daff", "ella", "lore", "gizm","damo","daff","daff", "elto"),
                weight= c(2,4,6,8,10,1,1,1),
                date=c("2017-07-09", "2017-07-09","2017-07-09","2017-07-09","2017-07-09","2017-07-09","2017-07-09","2017-07-09"))

#create networks
net1 <-create.a.network(df1, directed = FALSE)
net1.d <-create.a.network(df1, directed = TRUE)
net3<- create.a.network(df)
net3.d <- create.a.network(df, directed = TRUE)

#create identical network
net2<-net1
net2.d <- net1.d

# Test Cosine
cosine_between_graphs(net2,net3, directed=FALSE,center=F,considerZeros = F)
cosine_between_graphs(net1.d,net2.d, directed=TRUE)
cosine_between_graphs(net1.d,net2.d, directed=FALSE)


# Test Correlation
cor_between_graphs(net2,net3, directed=FALSE)
cor_between_graphs(net1.d,net2.d, directed=TRUE)
cor_between_graphs(net1,net2, directed=FALSE)

# Test distance
dist_between_graphs(net2,net3, directed=FALSE)
dist_between_graphs(net1.d,net2.d, directed=TRUE)
dist_between_graphs(net1,net2, directed=FALSE)

#cosine_between_graphs(net1,net3) # when completely different

## Test: measurment_functions_DYADS.R

## for dyad_weight and dyad_sum they have the same description on github.
# However if network is directed and you have reciprical edges, dyad_weight will distinguish them
# with dyad_sum, the egde direction isn't taken into account so that the total sum of weight (for the reciprical dyads) is given for each reciprical dyad
# you can play around with these 2 netowrks and you'll understand what I mean haha.

## definition: take the sum of weight at dyad level

## def:not sure..at directed level, the weight sum of the dyad is duplicated when reciprical dyads

#maybe delete dyad_sum?

## here take the mean weight at dyad level:

## ADD: give better description and say that if net non directed then = dyad_weight


## dyad difference in weight within the SAME network (so only look @ reciprical dyads)


## different in weight BETWEEN networks
# here take same net1 but multiply weights by 2
df1.bis <-data.frame(from=c("daff","damo","lore","elto","damo","gizm", "daff","lore"),
                     to= c("gizm", "ella", "ella", "ella","daff","daff","ella", "elto"),
                     weight= c(4,8,12,16,20,2,2,2),
                     date=c("2017-07-09", "2017-07-09","2017-07-09","2017-07-09","2017-07-09","2017-07-09","2017-07-09","2017-07-09"))

net1.bis<- create.a.network(df1.bis)
net1.bis.d<- create.a.network(df1.bis,directed = TRUE)

### The user HAS to specify directed. otherwise won't check

#dyad_weight(net1.bis.d)


##check extract lagged measure network
#create first df that will give similar networks
same.net.df<- data.frame ( from=c("laur", "mori", "laur", "mori", "laur", "mori"),
                           to = c("mori", "daff", "mori", "daff","mori", "daff"),
                           date= c("2015-01-01", "2015-01-10", "2015-01-15","2015-01-20","2015-01-28","2015-01-30"))
same.net.df$from<- as.character(same.net.df$from)
same.net.df$to<- as.character(same.net.df$to)
same.net.df$date<- ymd(same.net.df$date)
lagged.netlist<-extract_networks(same.net.df, windowsize = days(10), windowshift= days(9), directed=FALSE)
diff.net.df<- data.frame ( from=c("laur", "laur", "laur", "laur", "mori", "laur"),
                           to = c("mori", "mori", "mori", "mori","daff", "mori"),
                           date= c("2015-01-01", "2015-01-11", "2015-01-18","2015-01-26","2015-02-03","2015-02-12"))
diff.net.df$from<- as.character(diff.net.df$from)
diff.net.df$to<- as.character(diff.net.df$to)
diff.net.df$date<- ymd(diff.net.df$date)
lagged.diff.netlist<-extract_networks(diff.net.df, windowsize = days(8), windowshift= days(8), directed=FALSE)




# Tests
test_that("Degree mean is good I guess", {
  expect_equal(as.numeric(degree_mean(net1)), 14/6)
  expect_equal(as.numeric(degree_mean(net1.d)),16/6)
  expect_equal(as.numeric(eigen_mean(net1)), mean(eigen_centrality(net1)$vector))
  expect_equal(as.numeric(eigen_mean(net1.d)),mean(eigen_centrality(net1.d)$vector))
  expect_equal(as.numeric(cosine_between_graphs(net1,net2)), 1)
  expect_equal(as.numeric(cosine_between_graphs(net2,net3)), 0)
  expect_equal(as.numeric(cosine_between_graphs(net2.d,net3.d)), 0)
  expect_equal(as.numeric(dyad_weight(net1)[1]), 10)
  expect_equal(as.numeric(dyad_weight(net1.d)[1]),2)
  expect_equal(as.numeric(dyad_change(net1,net2)[1]), 0)
  expect_equal(as.numeric(dyad_change(net1,net1.bis)[1]), 10)
  expect_equal(as.numeric(dyad_change(net1.d,net2.d)[1]), 0)
  expect_equal(as.numeric(dyad_change(net1.d,net1.bis.d, directed = TRUE)[1]), 2)
  expect_equal(as.numeric(dyad_diff(net1)[1]), 0) # will get 0s if non directed (to add)
  expect_equal(as.numeric(dyad_diff(net1.d)[1]),1)
  expect_equal(as.numeric(dyad_mean(net1)[2]), 3)
  expect_equal(as.numeric(dyad_mean(net1.d)[1]),1.5)
  expect_equal(extract_lagged_measure_network(lagged.netlist, lag=1,measureFun = cosine_between_graphs, firstNet = TRUE)[2,1],1) # check that 2nd network same as 1st one
  expect_equal(extract_lagged_measure_network(lagged.diff.netlist, lag=1,measureFun = cosine_between_graphs, firstNet = FALSE)[5,1],0) ## check that the 4th andf 5th network are different
  expect_equal(cosine_between_nodes (net1,net2, directed=FALSE)[1,1], 1) ## when net1 = net 2
  expect_equal(as.numeric(cosine_between_nodes(net2,net3)[1,1] ), 0) #when net2 completely different from  net3
  expect_equal(as.numeric(cosine_between_nodes(net1.d,net2.d,directed=TRUE)[1,1]), 1) ##with directed network
  expect_equal(as.numeric(cosine_between_nodes(net2.d,net3.d)[1,1] ), 0)
  expect_equal(as.numeric(cosine_between_graphs(net1.d,net2.d, directed=TRUE)), 1)
  expect_equal(as.numeric(cor_between_graphs(net1,net2)$estimate), 1)
  expect_equal(as.numeric(round(cor_between_graphs(net2,net3)$estimate,digits = 6)), -0.476419)
  expect_equal(as.numeric(cor_between_graphs(net2.d,net1.d, directed = T)$estimate), 1)
  expect_equal(as.numeric(cor_between_graphs(net2,net1, directed = F)$estimate), 1)
  expect_equal(as.numeric(dist_between_graphs(net1,net2)), 0)
  expect_equal(as.numeric(round(dist_between_graphs(net2,net3),digits = 4)), 21.2132)
  expect_equal(as.numeric(dist_between_graphs(net2.d,net1.d, directed = T)), 0)
  expect_equal(as.numeric(dist_between_graphs(net2,net1, directed = F)), 0)
  expect_equal(as.numeric(edge.weight.skewness(net1))[2], 0)
})
