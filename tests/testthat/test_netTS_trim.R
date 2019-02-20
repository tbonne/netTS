context("Trim tests")
library(netTS)
library(lubridate)

#set random
set.seed(123)

#simulate some data
events1 <-sim.events.data(nodes=3,sampling.periods = 100,sampling.periods.per.day = 1)
events1 <- events1$beha
events2 <-sim.events.data(nodes=2,sampling.periods = 100,sampling.periods.per.day = 1)
events2 <- events2$beha
events2$sampleID<-events2$sampleID+100
events2$date<-events2$date+days(100)
eventsAll <- rbind(events1,events2)

str(eventsAll)
eventsAll$from <- as.character(eventsAll$from)
eventsAll$to <- as.character(eventsAll$to)

#extract measures
node.ts <- nodeTS(eventsAll)
dyad.ts.unDir <- dyadTS(eventsAll, directed = FALSE)
dyad.ts.Dir <- dyadTS(eventsAll, directed = TRUE)

#trim (default)
node.ts.trim <- trim_nodes(node.ts, eventsAll)
dyad.ts.trim.unDir <- trim_dyads(dyad.ts.unDir, eventsAll, directed = FALSE)
dyad.ts.trim.Dir <- trim_dyads(dyad.ts.Dir, eventsAll, directed = TRUE)

#trim (by user specification)

enter_leave <- data.frame(ID=c("1","2","3"), enter=c("2002-07-24","2002-07-24","2002-07-24"),leave=c("2003-02-08","2003-02-08","2002-10-30"))
enter_leave$enter <- ymd(as.character(enter_leave$enter))
enter_leave$leave <- ymd(as.character(enter_leave$leave))

enter_leave_dyad <- data.frame(ID=c("1_2","2_3","1_3"), enter=c("2002-07-24","2002-07-24","2002-07-24"),leave=c("2003-02-08","2002-10-30","2002-10-30"))
enter_leave_dyad$enter <- ymd(as.character(enter_leave_dyad$enter))
enter_leave_dyad$leave <- ymd(as.character(enter_leave_dyad$leave))

enter_leave_dyad_dir <- data.frame(ID=c("1_2","3_2","2_3","3_1","2_1","1_3"), enter=c("2002-07-24","2002-07-24","2002-07-24","2002-07-24","2002-07-24","2002-07-24"),leave=c("2003-02-08","2002-10-30","2002-10-30","2002-10-30","2003-02-08","2002-10-30"))
enter_leave_dyad_dir$enter <- ymd(as.character(enter_leave_dyad_dir$enter))
enter_leave_dyad_dir$leave <- ymd(as.character(enter_leave_dyad_dir$leave))

node.ts.trim.user <- trim_nodes(node.ts, eventsAll, enter_leave=enter_leave)
dyad.ts.trim.unDir.user <- trim_dyads(dyad.ts.unDir, eventsAll, directed = FALSE,enter_leave=enter_leave_dyad)
dyad.ts.trim.Dir.user <- trim_dyads(dyad.ts.Dir, eventsAll, directed = TRUE, enter_leave=enter_leave_dyad_dir)

#true value
true.leave = ymd("2002-07-24")+days(100)


#test that it worked properly
test_that("Trimmed network is equal to expectations", {
   expect_equal(is.na(node.ts.trim[70,2]),FALSE )
   expect_equal(is.na(node.ts.trim[71,2]),TRUE )
   expect_equal(is.na(node.ts.trim.user[69,2]),FALSE )
   expect_equal(is.na(node.ts.trim.user[70,2]),TRUE )
   expect_equal(is.na(dyad.ts.trim.unDir[70,2]),FALSE )
   expect_equal(is.na(dyad.ts.trim.unDir[71,2]),TRUE )
   expect_equal(is.na(dyad.ts.trim.unDir.user[70,2]),FALSE )
   expect_equal(is.na(dyad.ts.trim.unDir.user[71,2]),TRUE )
   expect_equal(is.na(dyad.ts.trim.Dir[70,2]),FALSE )
   expect_equal(is.na(dyad.ts.trim.Dir[71,2]),TRUE )
   expect_equal(is.na(dyad.ts.trim.Dir.user[70,2]),FALSE )
   expect_equal(is.na(dyad.ts.trim.Dir.user[71,2]),TRUE )
 })
