context("Bootstrap and permutation function tests")
library(testthat)

## /permutation_functions.R

# test perm.events function ... more specifically test the sample part (gist of the function)
set.seed(123)
df1<-groomEvents
df.original<-df1[,1:2]
df.original[,1]<- as.character(df.original[,1])
df.original[,2]<- as.character(df.original[,2])
data<- df1[,1:2]

## test that the swapping of individuals work by creating a df of only 2 individuals, using a non directed approach
df.perm<- data.frame(from = c("daff","daff","daff","daff","daff","daff"),
                     to = c("elto","elto","elto","elto","elto","elto"),
                     date = c("2019-01-01","2019-01-01","2019-01-01","2019-01-02","2019-01-02","2019-01-02") )
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

test_that("perm.events works good", {
  expect_equal( as.numeric(perm.noDiff[1]), as.numeric(perm.noDiff[2]) )
  expect_equal( as.numeric(perm.noDiff.directed[1]), as.numeric(perm.noDiff.directed[2]))
  expect_equal( as.numeric(perm.noDiff.edge[1]), as.numeric(perm.noDiff.edge[2]))
  expect_equal( as.numeric(perm.noDiff.degseq[1]), as.numeric(perm.noDiff.degseq[2]))
  expect_equal( as.numeric(perm.ci[1]), as.numeric(perm.ci[2]))
})   ## test that permutation is occuring and that the original df and sampled df differ

