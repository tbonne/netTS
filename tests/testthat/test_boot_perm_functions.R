context("Bootstrap and permutation function tests")
library(testthat)

## /permutation_functions.R

# test perm.events function
set.seed(123)
df1<-groomEvents[1:1000,]


## test that the swapping of individuals work by creating a df of only 2 individuals, using a non directed approach
df.perm<- data.frame(from = c("daff","daff","daff","daff","daff","daff","daff","daff","daff","daff","daff","daff"),
                     to = c("elto","elto","elto","elto","elto","elto","elto","elto","elto","elto","elto","elto"),
                     date = c("2019-01-01","2019-01-01","2019-01-01","2019-01-02","2019-01-02","2019-01-02","2019-01-03","2019-01-03","2019-01-03","2019-01-04","2019-01-04","2019-01-04") )
df.perm$from<- as.character(df.perm$from)
df.perm$to<-as.character(df.perm$to)
df.perm$date<-ymd(df.perm$date)
data2<-df.perm

df.perm3<- data.frame(from = c("daff","damo","saff","lori","band","rock"),
                     to = c("elto","elto","elto","elto","elto","elto"),
                     weight = c(2,2,2,2,2,2))

#calculate some permutations
perm.noDiff<-perm.events(data2, nperm =100,measureFun = eigen_mean)
perm.noDiff.directed <- perm.events.directed(data2, nperm =100,measureFun = eigen_mean)
perm.noDiff.edge <-perm.edge.weights(df.perm3, nperm =100,measureFun = eigen_mean)
perm.noDiff.degseq <- perm.edge.degseq(df.perm3, nperm =100,measureFun = eigen_mean)

graphlist.temp<-extract_networks(data2, windowsize=days(1),windowshift = days(1))
perm.ci<-permutation.graph.values(data2, windowsize=days(1),windowshift = days(1), measureFun = eigen_mean, graphlist = graphlist.temp, nperm = 10)

#bootstrap test
mean_str <- function(x){
  mean(strength(x))
}
mean_deg <- function(x){
  mean(degree(x))
}

out.w1 <-check.windowsize(data=df1, windowsize=days(30),windowshift = days(30), measureFun=NULL,corFun = 1,boot.samples=10)
out.w2 <-check.windowsize(data=df1, windowsize=days(30),windowshift = days(30), measureFun=NULL,corFun = 2,boot.samples=10)
out.w3 <-check.windowsize(data=df1, windowsize=days(30),windowshift = days(30), measureFun=NULL,corFun = 3,boot.samples=10)

out.b1 <-convergence.check.boot(data=df1, windowsize=days(30),windowshift = days(30), measureFun = strength ,corFun = 1,boot.samples=10)
out.b2 <-convergence.check.boot(data=df1, windowsize=days(30),windowshift = days(30), measureFun = strength ,corFun = 2,boot.samples=10)
out.b3 <-convergence.check.boot(data=df1, windowsize=days(30),windowshift = days(30), measureFun = mean_str ,corFun = 3,boot.samples=10)

out.bg1 <- convergence.check.boot.graph(data=df1, windowsize=days(30),windowshift = days(30), measureFun = eigen_mean,boot.samples=10)

out.v1<-convergence.check.var(df1, windowsize_min=days(1),windowsize_max=days(10),by=days(1), windowshift=days(30), directed = FALSE, measureFun=igraph::edge_density, SRI=FALSE)

out.t1<-check.timescale(df1, windowsize_min=days(1),windowsize_max=days(10),by=days(1), windowshift=days(30), directed = FALSE, measureFun=igraph::edge_density, summaryFun = var, SRI=FALSE)

out.reg <- convergence.check.value(df1, windowsize=days(120),windowshift = days(120),max.subsample.size = 30, SRI = FALSE, n.boot = 10, probs = c(0.025, 0.975))

test_that("perm.events works good", {
  expect_equal( as.numeric(perm.noDiff[1]), as.numeric(perm.noDiff[2]) )
  expect_equal( as.numeric(perm.noDiff.directed[1]), as.numeric(perm.noDiff.directed[2]))
  expect_equal( as.numeric(perm.noDiff.edge[1]), as.numeric(perm.noDiff.edge[2]))
  expect_equal( as.numeric(perm.noDiff.degseq[1]), as.numeric(perm.noDiff.degseq[2]))
  expect_equal( as.numeric(perm.ci[1,1]), as.numeric(perm.ci[1,2]))
  expect_equal( round(out.w1[1,1],3), 0.771)
  expect_equal( round(out.w2[1,1],3), 0.785)
  expect_equal( round(out.w3[1,1],3), 12.967)
  expect_equal( round(out.b1[1,1],3), 0.882)
  expect_equal( round(out.b2[1,1],3), 0.927)
  expect_equal( round(out.b3[1,1],3), 0.041)
  expect_equal( round(out.bg1[1,1],3), 0.161)
  expect_equal( round(out.v1[1,2],3), 0.185)
  expect_equal( round(out.t1[6,2],3), 0.012)
  expect_equal( round(out.reg[1,1],3), 0.067)

})   ## test that permutation is occuring and that the original df and sampled df differ

