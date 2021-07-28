### In an attempt to speed up the analysis all modifications will be done using data.table set function
### This breaks with the standard assumptions in R that functions have no side effects
### I may change this in the future if it becomes too confusing/annoying

#' Convert glucose concentrations from mg/dL to mmol/L
#'
#' @param x x glucose monitoring data.table
#'
#' @return glucose monitoring data.table (also updated as a side effect)
#' @export
#'
#' @examples
#' \dontrun{
#' to-do
#' }
convert_mgdl_2_molL <- function(x){
    mgDl_2_mmolL <- 18.018018
    data.table::set(x, j = "Glucose", value = x[["Glucose"]]/mgDl_2_mmolL)
    x[]
}

#' Update Date column to include time zone information
#'
#' @param x glucose monitoring data.table
#' @param tz timezone to use
#'
#' @return glucose monitoring data.table (also updated as a side effect)
#' @export
#'
#' @import data.table
#'
#' @examples
#' \dontrun{
#' to-do
#' }
### Excel does not handle time zones
fix_date <- function(x, tz) {
  Date <- ElapsedTime <- NULL # Silence build notes
  x[, Date:=lubridate::force_tz(Date[1], tz) + ElapsedTime - ElapsedTime[1]]
  x[]
}


#' Update data to be sampled every minute
#'
#' @param x glucose monitoring data.table
#' @param sample_rate sampling rate in seconds
#'
#' @return glucose monitoring data.table (NOT updated as a side effect)
#' @export
#' 
#' @import data.table
#'
#' @examples
#' \dontrun{
#' to-do
#' }
# Sometime there will be datapoints taken at odd intervals, discard those
# Sometimes datapoints will be missing, add rows to contain those
uniformize_sample_rate <- function(x, sample_rate = 60)
{
  # Only very few timepoints should be affected by dropping timepoints where a full minute is missing
  before <- nrow(x)
  time_correct <- c(x[["ElapsedTime"]][1], diff(x[["ElapsedTime"]])) == sample_rate
  idx_bad <- which(!time_correct)
  
  x <- x[(time_correct), ]
  excluded_message <- paste("Datapoints excluded due to uneven sampling rate:", before - nrow(x), "\nThese dropouts happen at rows:", paste(idx_bad, collapse = "; "))

  # Create all datapoints if every timepoint was measured
  y <- data.table::data.table(Date = seq(x[["Date"]][1], x[["Date"]][nrow(x)], by = sample_rate))
  y <- data.table::merge.data.table(y, x, by = "Date", all = TRUE)

  Date <- ElapsedTime <- Sample_ID <- NULL # Silence build notes
  # Some timepoints are missing data
  y[is.na(ElapsedTime), ElapsedTime:=lubridate::as.duration(Date - y[["Date"]][1] + y[["ElapsedTime"]][1])]
  y[, Sample_ID:=x[["Sample_ID"]][1]]
  
  if(!all(round(diff(y[["ElapsedTime"]]), digits = 4) == sample_rate)) stop("Failed to make uniform sampling rate")
  
  y[]
}


#' Add a number of derived time related data
#'
#' @param x glucose monitoring data.table
#' @param light_on Period object with time the light is turned on
#' @param light_off Period object with time the light is turned off 
#' @param dst string either "sync" or "independent". If "sync" light is turned on and off 
#' when the computer clock are the the timepoints specified in light_on and light_off. If
#' "independent" light is running on a 24 hour cycle independent of computer time.
#' These options only matter if the experiment covers a daylight savings time change.
#'
#' @return glucose monitoring data.table (also updated as a side effect)
#' @export
#' 
#' @import data.table
#'
#' @examples
#' \dontrun{
#' to-do
#' }
add_derived_date_info <- function(x, light_on, light_off, dst)
{
  Date <- ElapsedTime <- .N <- Light_on <- PhaseId <- Day <- Week <- ZT_exact <- ZT <- NULL # Silence build notes
  first_timepoint <- x[1, Date - ElapsedTime]
  last_time_point <- x[.N, Date]
  experimental_days <- lubridate::day(lubridate::as.period(last_time_point - first_timepoint)) + 1
  
  
  first_light_on <- lubridate::ymd_hms(paste0(lubridate::year(first_timepoint), "-",
                                              lubridate::month(first_timepoint), "-",
                                              lubridate::day(first_timepoint),  " ",
                                              lubridate::hour(light_on), ":",
                                              lubridate::minute(light_on), ":",
                                              "0"), tz = lubridate::tz(x$Date))
  
  # Adding 24 hours will, for some reason, add 23 or 25 when crossing DST
  # We need to work in UTC for lubridate to behave properly when adding 24 hours
  first_light_on <- lubridate::with_tz(first_light_on, "UTC")
  
  first_light_off <- first_light_on + (light_off - light_on)
  # The plus 0.1 second is so that measurements taken at the time the
  # light goes on is counted as dark.
  first_light_on <- first_light_on + lubridate::seconds(0.1)
  if (dst == "independent"){
    # If we add hours while in UTC and *then* convert back light on/off
    # cycles every 24 hours, independent of computer time.
    light_on_all <- first_light_on + lubridate::hours((0:experimental_days)*24)
    light_off_all <- first_light_off + lubridate::hours((0:experimental_days)*24)
    
    light_on_all <- lubridate::with_tz(light_on_all, lubridate::tz(x$Date))
    light_off_all <- lubridate::with_tz(light_off_all, lubridate::tz(x$Date))
  } else if (dst == "sync"){
    # If we go back to whatever timezone we were in and add days, timezone is 
    # respected.
    first_light_on <- lubridate::with_tz(first_light_on, lubridate::tz(x$Date))
    first_light_off <- lubridate::with_tz(first_light_off, lubridate::tz(x$Date))
    
    light_on_all <- first_light_on + lubridate::days(0:experimental_days)
    light_off_all <- first_light_off + lubridate::days(0:experimental_days)
  } else {
    stop("Daylight savings time setting must be 'independent' or 'sync'")
  }
  light_periods <- as.list(lubridate::interval(light_on_all, light_off_all))
  
  x[, Light_on:=FALSE]
  x[lubridate::`%within%`(Date, light_periods), Light_on:=TRUE]
  
  x[, PhaseId:=cumsum(c(1, abs(diff(Light_on))))]
  x[, Day:=c(0, diff(Light_on))]
  x[Day == -1, Day:=0]
  x[, Day:=cumsum(Day)]
  x[, Week:=Day %/% 7]
  x[, Week:=Week + 1]
  x[, Day:=Day + 1]
  x[, ZT_exact:=lubridate::time_length(Date - light_on_all[[.GRP]], unit = "hours") %% 24, by = "Day"]
  x[, ZT:=floor(ZT_exact)]
  x[]
}

#' Exclude timepoints
#'
#' @param x glucose monitoring data.table
#' @param exclusions data.frame-like object with POSIXct columns Start and End denoting exclusion intervals
#'
#' @return x glucose monitoring data.table (also updated as a side effect)
#' @export
#'
#' @import data.table
#'
#' @examples
#' \dontrun{
#' to-do
#' }
exclude_timepoints <- function(x, exclusions, invert)
{
  included <- Date <- NULL # Silence build notes
  if (is.null(x[["included"]])) x[, included:=TRUE]

  if (nrow(exclusions) > 0){
    exclusion_durations <- as.list(lubridate::interval(exclusions$Start, exclusions$End + 1)) #The plus one adds one second and ensures that the ending minute is also excluded.
    
    if (invert) {
      x[!lubridate::`%within%`(Date, exclusion_durations), included:=FALSE]
    } else {
      x[lubridate::`%within%`(Date, exclusion_durations), included:=FALSE]
    }
    
  }
  x[]
}


#' Linear interpolation of glucose for small gaps
#'
#' @param x glucose monitoring data.table
#' @param max_gap maximum number of missing datapoints to impute
#' @param column character, column to impute
#'
#' @return x glucose monitoring data.table (NOT updated as a side effect)
#' @export
#'
#' @import data.table
#' 
#' @examples
#' \dontrun{
#' to-do
#' }
linear_imputation <- function(x, column, max_gap)
{
  # Sometimes there are no e.g. activity data.
  if(all(is.na(x[[column]]))) return(x)
  
  # Which runs have a missing value
  value_na_rle <- rle(is.na(x[[column]]))

  # If the first or the last element(s) are missing we cannot interpolate them,
  # so they are discarded instead
  n_runs <- length(value_na_rle$values)

  bad_measurements <- integer(0)
  if (value_na_rle$values[1]) {
    bad_measurements <- seq(from = 1, to = value_na_rle$lengths[1])
  }

  if (value_na_rle$values[n_runs]) {
    start_of_last_run <- cumsum(value_na_rle$lengths)[n_runs - 1] + 1
    bad_measurements <- c(bad_measurements,
                         seq(from = start_of_last_run, to = nrow(x)))
  }

  if (length(bad_measurements) > 0) {
    x <- x[-bad_measurements]
  }

  # Which runs have a missing value and is shorter than the max gap?
  value_na_rle <- rle(is.na(x[[column]]))
  to_interpolate <- value_na_rle$values & value_na_rle$lengths <= max_gap

  # Which indexes should be interpolated?
  ends <- cumsum(value_na_rle$lengths)
  starts <- c(0, ends[-length(ends)]) + 1
  idxs <- Map(seq, starts[to_interpolate], ends[to_interpolate])

  # Simple linear interpolation
  for (i in idxs){
    before_val <- x[[column]][min(i) - 1]
    after_val <- x[[column]][max(i) + 1]
    interpolated_values <- seq(from = before_val,
                              to = after_val,
                              length.out = length(i) + 2)
    interpolated_values <- interpolated_values[-c(1, length(interpolated_values))]
    data.table::set(x, i, column, interpolated_values)
    data.table::set(x, i, "included", TRUE)
  }
  x[]
}

#' Mark timepoints without glucose data as excluded
#'
#' @param x glucose monitoring data.table
#'
#' @return glucose monitoring data.table
#' @export
#' @import data.table
#'
#' @examples
#' \dontrun{
#' to-do
#' }
exclude_missing <- function(x) {
  Glucose <- included <- NULL
  if (is.null(x[["included"]])) x[, included:=TRUE]
  x[is.na(Glucose), included:=FALSE]
  x[]
}

#' Find baseline glucose expression
#'
#' @param x glucose monitoring data.table
#' @param baseline_window number of timepoints to include in baseline calculations
#' @param max_missing_baseline maximum number of excluded datapoints before baseline is set to NA
#'
#' @return glucose monitoring data.table
#' @export
#' @import data.table
#'
#' @examples
#' \dontrun{
#' to-do
#' }
find_baseline <- function(x, baseline_window, max_missing_baseline)
{
  if (is.null(x[["included"]])) x[, included:=TRUE]
  included <- baseline <- Glucose <- NULL
  changePoints <- c(0, 0, diff(sign(diff(x[["Glucose"]]))))
  x[, baseline:=NA_real_]
  x[(included) & changePoints != 0, baseline:=Glucose]
  # data.table::frollapply is twice as fast as zoo::rollapply
  # but lacks the 'partial' option. These extra steps replicate the 'partial' functionality
  
  n_missing <- data.table::frollsum(c(rep(NA, baseline_window/2), !x$included, rep(NA, baseline_window/2)),
                       n = baseline_window,
                       na.rm = TRUE,
                       fill = NA,
                       align = "center")
  
  median_baseline <- data.table::frollapply(c(rep(NA, baseline_window/2), x$baseline, rep(NA, baseline_window/2)),
                                           n = baseline_window,
                                           FUN = stats::median,
                                           na.rm = TRUE,
                                           fill = NA,
                                           align = "center")
  
  if (baseline_window %% 2 == 1) {
    idx_included <- (ceiling(baseline_window/2)):(length(median_baseline) - floor(baseline_window/2))
  } else {
    idx_included <- (baseline_window/2 + 1):(length(median_baseline) - baseline_window/2)
  }
  
  median_baseline <- median_baseline[idx_included]
  n_missing <- n_missing[idx_included]
  x[, baseline:=median_baseline]
  
  x[n_missing > max_missing_baseline, baseline:=NA]

  x[]
}

#' Find baseline glucose expression, ignoring high glucose values
#'
#' @param x glucose monitoring data.table
#' @param baseline_window number of timepoints to include in baseline calculations
#' @param max_missing_baseline maximum number of excluded datapoints before baseline is set to NA
#'
#' @return glucose monitoring data.table
#' @export
#' @import data.table
#'
#' @examples
#' \dontrun{
#' to-do
#' }
find_baseline_ignore_high <- function(x, baseline_window, max_missing_baseline)
{
  if (is.null(x[["included"]])) x[, included:=TRUE]
  included <- baseline <- Glucose <- NULL
  changePoints <- c(0, 0, diff(sign(diff(x[["Glucose"]]))))
  x[, baseline:=NA_real_]
  
  # High turning points are likely caused by spikes in glucose, there are no corresponding dips,
  # so if we exclude these spikes we get a more stable estimate
  
  # data.table::frollapply is twice as fast as zoo::rollapply
  # but lacks the 'partial' option. These extra steps replicate the 'partial' functionality
  glucose <- x$Glucose
  glucose[!x$included] <- NA_real_
  roll_mean_glucose <- data.table::frollmean(c(rep(NA, baseline_window/2), glucose, rep(NA, baseline_window/2)),
                                    n = baseline_window,
                                    na.rm = TRUE,
                                    fill = NA,
                                    align = "center")
  
  if (baseline_window %% 2 == 1) {
    idx_no_edges <- (ceiling(baseline_window/2)):(length(roll_mean_glucose) - floor(baseline_window/2))
  } else {
    idx_no_edges <- (baseline_window/2 + 1):(length(roll_mean_glucose) - baseline_window/2)
  }
  
  x[(included) & changePoints != 0 & Glucose < roll_mean_glucose[idx_no_edges], baseline:=Glucose]
  n_missing <- data.table::frollsum(c(rep(NA, baseline_window/2), !x$included, rep(NA, baseline_window/2)),
                                    n = baseline_window,
                                    na.rm = TRUE,
                                    fill = NA,
                                    align = "center")
  
  median_baseline <- data.table::frollmean(c(rep(NA, baseline_window/2), x$baseline, rep(NA, baseline_window/2)),
                                            n = baseline_window,
                                            na.rm = TRUE,
                                            fill = NA,
                                            align = "center")
  
  if (baseline_window %% 2 == 1) {
    idx_included <- (ceiling(baseline_window/2)):(length(median_baseline) - floor(baseline_window/2))
  } else {
    idx_included <- (baseline_window/2 + 1):(length(median_baseline) - baseline_window/2)
  }
  
  median_baseline <- median_baseline[idx_included]
  n_missing <- n_missing[idx_included]
  x[, baseline:=median_baseline]
  
  x[n_missing > max_missing_baseline, baseline:=NA]
  
  x[]
}

#' Smooth baseline
#'
#' @param x glucose monitoring data.table
#' @param baseline_window number of timepoints to include in baseline calculations
#' 
#' @return glucose monitoring data.table
#' @export
#' @import data.table
#'
#' @examples
#' \dontrun{
#' to-do
#' }
smooth_baseline <- function(x, baseline_window)
{
  if(is.null(x[["baseline"]])) stop("baseline must be calculated before smoothing")
  baseline <- NULL
  excluded_background <- is.na(x$baseline)
  mean_baseline <- data.table::frollmean(c(rep(NA, baseline_window/2), x$baseline, rep(NA, baseline_window/2)),
                              n = baseline_window,
                              na.rm = TRUE,
                              fill = NA,
                              align = "center")
  
  if (baseline_window %% 2 == 1) {
    idx_included <- (ceiling(baseline_window/2)):(length(mean_baseline) - floor(baseline_window/2))
  } else {
    idx_included <- (baseline_window/2 + 1):(length(mean_baseline) - baseline_window/2)
  }
  
  mean_baseline <- mean_baseline[idx_included]
  
  x[, baseline:=mean_baseline]
  x[excluded_background, baseline:=NA]
  x
}

#' Tag excursions
#'
#' @param x glucose monitoring data.table
#' @param excursion_low increase below baseline before excursion is tagged
#' @param excursion_high increase above baseline before excursion is tagged 
#'
#' @return glucose monitoring data.table
#' @import data.table
#' @export
#'
#' @examples
#' \dontrun{
#' to-do
#' }
tag_excursions <- function(x, excursion_low, excursion_high) {
  if(is.null(x[["baseline"]])) stop("baseline must be calculated before tagging excursions")
  
  excursion <- Glucose <- baseline <- deviation_from_baseline <- NULL
  
  x[, excursion:=Glucose - baseline]
  x[, deviation_from_baseline:=Glucose - baseline]
  x[excursion > excursion_low & excursion < excursion_high, excursion:=NA]
  x[]
}

 
#' Find peaks and nadirs
#'
#' @param x glucose monitoring data.table
#' @param max_min_window integer, window to scan for local peaks/nadirs. Must be odd.
#'
#' @return glucose monitoring data.table
#' @import data.table
#' @export
#'
#' @examples
#' \dontrun{
#' to-do
#' }
find_peaks_and_nadirs <- function(x, max_min_window)
{
  if(is.null(x[["excursion"]])) stop("excursions must be calculated before tagging peaks and nadirs")
  if(max_min_window %% 2 != 1) stop("max_min_window must be odd")
  
  peak <- included <- nadir <- excursion <- NULL
  
  roll_fun <- function(which_fun){
    which_fun_na_fix <- function(x) {
      out <- which_fun(x)
      if (length(out) == 0) {
        NA
      } else {
        out
      }
    }
    out <- data.table::frollapply(x = c(rep(NA, max_min_window/2), x$excursion, rep(NA, max_min_window/2)),
                                        n = max_min_window, 
                                        FUN = which_fun_na_fix, 
                                        fill = NA,
                                        align = "center"
    )
    # returns 8.236138e-312 instead of 0 on NA
    out <- floor(out)
    
    #max_min_window is always odd
    idx_included <- (ceiling(max_min_window/2)):(length(out) - floor(max_min_window/2))
    out[idx_included]
  }

  # roll_fun returns the position in the defined window, centered on the
  # observation which has the highest/lowest value
  # when the observation with the highest/lowest value is the centerpoint a
  # local minimum/maximum is found
  # in cases with ties, which.min / which.max returns the first
  midpoint <- ceiling(max_min_window/2)

  peaks <- roll_fun(which.max) == midpoint

  nadirs <- roll_fun(which.min) == midpoint

  # A peak is a peak if:
  #     1. it is included in the data
  #     2. it has a positive excursion
  #     3. it is the highest value in the window
  x[, peak:=FALSE]
  x[(included) & sign(excursion) == 1 & peaks, peak:=TRUE]

  # Similar for nadirs
  x[, nadir:=FALSE]
  x[(included) & sign(excursion) == -1 & nadirs, nadir:=TRUE]

  # Sometimes two peaks in the same window has exactly the same value
  # In those cases keep only the first

  scan_fn <- function(window){
    any(abs(window[-midpoint] - window[midpoint]) < 10^-4)
  }
  
  duplicate_peak <- data.table::frollapply(x = c(rep(NA, midpoint), x$excursion),
                                          n = midpoint, 
                                          FUN = scan_fn, 
                                          fill = NA,
                                          align = "right"
  )
  idx_included <- (midpoint + 1):length(duplicate_peak)
  duplicate_peak <- duplicate_peak[idx_included]
  duplicate_peak <- as.logical(duplicate_peak)
  
  x[peak & duplicate_peak, peak:=FALSE]
  x[nadir & duplicate_peak, nadir:=FALSE]
  x[]
}

#' Tag events
#'
#' @param x glucose monitoring data.table
#' @param event character (regex) to denoting event(s) to detect
#' @param before int, minutes before event to tag
#' @param after int, minutes after event to tag
#'
#' @return glucose monitoring data.table
#' @export
#'
#' @examples
#' \dontrun{
#' to-do
#' }
tag_events <- function(x, event, before, after) {
  x[, Event_ID:=stringr::str_detect(Event, event)]
  x[is.na(Event_ID), Event_ID:=FALSE]
  x[, Event_ID:=cumsum(Event_ID)]
  
  x[, time_to_event:=(.N):1, by = "Event_ID"]
  x[, time_since_event:=1:(.N), by = "Event_ID"]
  
  x[time_to_event > before & time_since_event > after, Event_ID:=NA]
  x[time_since_event < after, Event_ID:=Event_ID-1]
  x[Event_ID == -1, Event_ID:=NA]
  x[, Event_ID:=as.numeric(as.factor(Event_ID))]
  x[, is_event:=!is.na(Event_ID)]
  
  x[, ZT_event_exact:=lubridate::time_length(lubridate::duration(time_since_event, units = "minute"), unit = "hours") %% 24, by = "Event_ID"]
  x[, ZT_event:=floor(ZT_event_exact)]
  
  x[, time_to_event:=NULL]
  x[, time_since_event:=NULL]
  
  x[]
}

#' Run the standard pipeline for a single sample
#'
#' @param sample_id sample to be processed
#' @param cge cgm data object
#'
#' @return glucose monitoring data.table
#' @import data.table
#' @export
#'
#' @examples
#' \dontrun{
#' to-do
#' }
run_standard_preprocess_pipeline <- function(sample_id, cge) {
  # Silence no visible binding warnings
  Date <- NULL
  
  x <- cge$data[[sample_id]]
  
  if (get_option(cge, "mgdl_2_mmolL") == "y") {
    x <- convert_mgdl_2_molL(x)
  }
  
  x <- fix_date(x, get_option(cge, "time_zone"))
  x <- uniformize_sample_rate(x)
  x <- add_derived_date_info(x = x, 
                             light_on = get_option(cge, "light_on"), 
                             light_off = get_option(cge, "light_off"), 
                             dst = get_option(cge, "DST"))
  
  exclusions <- get_exclusions(cge, sample_id)
  x <- exclude_timepoints(x, exclusions, invert = get_option(cge, "invert_exclusions") == "y")
  
  x <- linear_imputation(x, column = "Glucose", max_gap = get_option(cge, "max_gap"))
  #x <- linear_imputation(x, column = "Temperature", max_gap = get_option(cge, "max_gap")) # Maybe reactivate at some point
  
  x <- find_baseline_ignore_high(x = x, 
                     baseline_window = get_option(cge, "baseline_window"), 
                     max_missing_baseline = get_option(cge, "max_missing_baseline")
  )
  
  x <- smooth_baseline(x = x, 
                         baseline_window = get_option(cge, "baseline_window"))
  
  x <- tag_excursions(x = x, 
                      excursion_low = get_option(cge, "excursion_low"), 
                      excursion_high = get_option(cge, "excursion_high"))
  
  x <- find_peaks_and_nadirs(x, max_min_window = get_option(cge, "max_min_window"))
  
  if (!is.na(get_option(cge, "event_letter"))){
    x <- tag_events(x = x, 
                   event = get_option(cge, "event_letter"), 
                   before = get_option(cge, "pre_event_window"), 
                   after = get_option(cge, "post_event_window"))
  }
  
  # UTC used throughout processing, reset prior to returning object
  x[, Date:=lubridate::with_tz(Date, get_option(cge, "time_zone"))]
  
  x[]
}