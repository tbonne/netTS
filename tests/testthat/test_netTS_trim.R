context("Trim tests")
library(netTS)

# node.ts <- nodeTS(groomEvents[1:200,])
# dyad.ts.unDir <- dyadTS(groomEvents[1:200,], directed = FALSE)
# dyad.ts.Dir <- dyadTS(groomEvents[1:200,], directed = TRUE)
#
# node.ts.trim <- trim_nodes(node.ts, groomEvents[1:200,])
# dyad.ts.trim.unDir <- trim_dyads(dyad.ts.unDir, groomEvents[1:200,], directed = FALSE)
# dyad.ts.trim.Dir <- trim_dyads(dyad.ts.Dir, groomEvents[1:200,], directed = TRUE)
#
# #test that it worked properly
# test_that("Trimmed network is equal to expectations", {
#   expect_equal(ncol(node.ts.trim),ncol(node.ts))
#   expect_equal(ncol(dyad.ts.trim.unDir),ncol(dyad.ts.unDir))
#   expect_equal(ncol(dyad.ts.trim.Dir),ncol(dyad.ts.Dir))
#   expect_equal(is.na(dyad.ts.trim.Dir[1,1]),TRUE)
#   expect_equal(dyad.ts.trim.unDir[6,2],14)
#   expect_equal(is.na(dyad.ts.trim.unDir[5,2]),TRUE)
# })
