context("input to igraph tests")
library(netTS)
library(lubridate)
library(igraph)
library(dplyr)

#make up a simple dataframe of observations
test.df <- data.frame(from = c("amy", "greg", "sarah", "olga"),
                      to = c("greg","amy", "amy", "greg"),
                      weight= c(1,2,3,4),
                      date = c("2014-10-23","2014-10-23","2014-10-23","2014-11-23" )
                      )

#create graph
out.graph.undirected <- create.a.network(test.df, directed=FALSE)
out.graph.directed <- create.a.network(test.df, directed=TRUE)

#test that it worked properly
test_that("Out network is equal to expectations from custom dataframe", {
  expect_equal(vcount(out.graph.undirected), 4)
  expect_equal(vcount(out.graph.directed), 4)
  expect_equal(ecount(out.graph.undirected), 3)
  expect_equal(ecount(out.graph.directed), 4)
  expect_equal(sum(E(out.graph.undirected)$weight - c(3,3,4)),0 )
  expect_equal(sum(E(out.graph.directed)$weight - c(1,2,4,3)) , 0)
})
