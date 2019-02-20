################# Some functions for measurment of sampling effort ###########################


############################################
#
#  Effort functions
#
############################################


#' Min/Max time per day
#'
#' This function will estimate the total time spent sampling, using the min and max observed values per day
#' @param df.window An events data frame with a column containing the date and time.
#' @importFrom dplyr group_by summarise
#' @importFrom lubridate minute hour date interval as.duration day
#' @importFrom stats time
#' @export
#'
#'
effort.time <- function(df.window){
  df.window$day <- as.numeric(as.duration(lubridate::interval(lubridate::date(min(df.window$date)),lubridate::date(df.window$date)  )),"days")
  df.window$time <- time(df.window$date)
  hours.min <- df.window %>% dplyr::group_by(day) %>% dplyr::summarise(min=min(lubridate::hour(date) + (lubridate::minute(date) )/60 ))
  hours.max <- df.window %>% dplyr::group_by(day) %>% dplyr::summarise(max=max(lubridate::hour(date) + (lubridate::minute(date) )/60 ))
  return(sum(hours.max-hours.min))
}

#' Unique sampeling IDs
#'
#' This function will estimate the total sampling, using a count of the time a sample was performed
#' @param df.window An events data frame with a column (sampleID) containing a unique id for every sampling event (i.e., group scan in a gambit of the group sampling procedure).
#' @export
#'
#'
effort.scan <- function(df.window){
  if(is.null(df.window$sampleID) )print("Please rename the column containing the sample IDs: sampleID")
  total.samples <- length(unique(df.window$sampleID))
  return(total.samples)
}
