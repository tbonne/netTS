context("Cosine measure tests")
library(netTS)

test.df.1 <- data.frame(from = c("amy", "greg", "sarah", "olga"),
                      to =     c("greg","amy", "amy", "greg"),
                      weight= c(1,2,3,4),
                      date = c("2014-10-23","2014-10-23","2014-10-23","2014-11-23" )
)

test.df.2 <- data.frame(from = c("amy", "greg", "sarah", "olga"),
                        to = c("greg","amy", "amy", "greg"),
                        weight= c(1,2,3,4),
                        date = c("2014-10-23","2014-10-23","2014-10-23","2014-11-23" )
)

test.df.3 <- data.frame(from = c("amy", "greg","olga", "olga"),
                        to = c("greg","amy", "arron","greg"),
                        weight= c(1,2,3,4),
                        date = c("2014-10-23","2014-10-23","2014-10-23","2014-11-23" )
)

#create networks
graph.1 <- create.a.network(test.df.1, directed=FALSE)
graph.2 <- create.a.network(test.df.2, directed=FALSE)
graph.3 <- create.a.network(test.df.3, directed=FALSE)

graph1 <- graph.1
graph2 <- graph.2

#create networks - directed
graph.1.d <- create.a.network(test.df.1, directed=TRUE)
graph.2.d <- create.a.network(test.df.2, directed=TRUE)
graph.3.d <- create.a.network(test.df.3, directed=TRUE)

graph1 <- graph.1.d
graph2 <- graph.3.d

test_that("Number of nodes equals number of columns", {
  expect_equal(as.numeric(cosine_between_nodes(graph.1,graph.2)[1]), 1)
  expect_equal(as.numeric(round(cosine_between_nodes(graph.2,graph.3)[1],7)), 0.7071068)
  expect_equal(as.numeric(cosine_between_nodes(graph.1.d,graph.2.d, directed = T)[1]), 1)
  expect_equal(as.numeric(round(cosine_between_nodes(graph.2,graph.3, directed = T)[1],7)), 0.7071068)
})
