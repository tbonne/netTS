context("Cosine dataframe structure")
library(netTS)

graph1 <- igraph::erdos.renyi.game(15,0.3,directed = T)
graph2 <- igraph::erdos.renyi.game(15,0.3,directed = T)
graph3 <- igraph::erdos.renyi.game(15,0.3,directed = F)
graph4 <- igraph::erdos.renyi.game(15,0.3,directed = F)


test_that("Number of nodes equals number of columns", {
  expect_equal(ncol(netTS::cosine_between_graphs_nodes(graph1,graph2)), 15)
  expect_equal(ncol(netTS::cosine_between_graphs_nodes(graph3,graph4)), 15)
})
