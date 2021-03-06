% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/measurment_functions_GRAPH.R
\name{dist_between_graphs}
\alias{dist_between_graphs}
\title{Estimate similarity between graphs using Euclidean distance}
\usage{
dist_between_graphs(graph1, graph2, directed = FALSE, considerZeros = TRUE)
}
\arguments{
\item{graph1}{first graph (igraph)}

\item{graph2}{second graph (igraph)}

\item{directed}{Whether the graph is directed or not (default = FALSE)}

\item{considerZeros}{Whether to treat edges present in one graph, but not the other, as zero (TRUE), or consider only those edges present in both graphs (FALSE).}
}
\description{
This function will calculate Euclidean distance between two graphs.
}
\examples{

#two random graphs
library(igraph)
graph1 <-erdos.renyi.game(n=10,p=0.1)
graph2 <-erdos.renyi.game(n=10,p=0.1)
cor_between_graphs(graph1,graph2)

}
