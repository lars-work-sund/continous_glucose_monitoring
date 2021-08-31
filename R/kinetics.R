#' Create empty kinetics table
#'
#' @return data.table
#'
#' @examples
#' \dontrun{
#' to-do
#' }
empty_kinetics <- function() {
  data.table::data.table( 
    excursion_ID         = integer(0),
    excursion            = numeric(0),
    excursion_duration   = integer(0), 
    single_peak          = logical(0),
    area_under_excursion = numeric(0),
    mean_uptake          = numeric(0),
    mean_clearance       = numeric(0),
    max_uptake           = numeric(0),
    max_clearance        = numeric(0),
    mean_excursion       = numeric(0),
    Sample_ID            = character(0),
    ElapsedTime          = lubridate::duration(),
    Light_on             = logical(0),
    phase_ID             = numeric(0),
    Day                  = numeric(0),
    Week                 = numeric(0),
    ZT                   = numeric(0),
    Event_ID             = integer(0),
    is_event             = logical(0),
    ZT_event             = numeric(0),
    peak_type            = character(0),
    n_peaks              = integer(0),
    peak_number          = numeric(0),
    Group                = character(0)
  )
}

#' Single peak kinetics
#'
#' @param x glucose monitoring data.table
#' @param excursion_high increase above baseline before excursion is tagged
#'
#' @return kinetics data.table
#' @import data.table
#' @export
#'
#' @examples
#' \dontrun{
#' to-do
#' }
single_peak_kinetics <- function(x, excursion_high){
  # Silence no visible binding warnings
  single_peak <- . <- excursion <- peak <- cumulative_uptake_time <- 
    cumulative_clearence_time <- uptake_slope <- clearance_slope <- Sample_ID <- 
    ElapsedTime <- Light_on <- phase_ID <- Day <- Week <- ZT <- Group <- NULL
  
  if (any(x$peak & x$single_peak)) {
    x[(single_peak), 
      .(excursion            = excursion[peak],
        excursion_duration   = .N, 
        single_peak          = TRUE,
        area_under_excursion = sum(excursion) - (excursion[[1]] + excursion[[.N]])/2, #To get area under curve first and last point are weighted by half 
        mean_uptake          = (excursion[peak] - excursion_high)/cumulative_uptake_time[peak],
        mean_clearance       = (excursion[peak] - excursion_high)/cumulative_clearence_time[peak],
        max_uptake           = max(uptake_slope),
        max_clearance        = min(clearance_slope),
        mean_excursion       = mean(excursion),
        Sample_ID            = Sample_ID[1],
        ElapsedTime          = ElapsedTime[peak],
        Light_on             = Light_on[peak],
        phase_ID             = phase_ID[peak],
        Day                  = Day[peak],
        Week                 = Week[peak],
        ZT                   = ZT[peak],
        Event_ID             = Event_ID[peak],
        is_event             = is_event[peak],
        ZT_event             = ZT_event[peak],
        peak_type            = "Single",
        n_peaks              = 1,
        peak_number          = 1,
        Group                = Group[peak]
      ), 
      by = c("excursion_ID")]
  } else {
    empty_kinetics()
  }
}

#' Single peak kinetics
#'
#' @param x glucose monitoring data.table
#' @param excursion_high increase above baseline before excursion is tagged
#'
#' @return kinetics data.table
#' @import data.table
#' @export
#'
#' @examples
#' \dontrun{
#' to-do
#' }
multi_peak_kinetics <- function(x, excursion_high){
  # Silence no visible binding warnings
  single_peak <- . <- excursion <- peak <- cumulative_uptake_time <- uptake_slope <- 
    clearance_slope <- time_to_next_peak <- Sample_ID <- ElapsedTime <- Light_on <- 
    phase_ID <- Day <- Week <- ZT <- Group <- excursion_ID <- peak_type <- 
    time_since_last_peak <- peak_number <- cumulative_clearence_time <- NULL
  
  x <- x[(!single_peak)]
  
  first_peak <- function(x) min(which(x))
  last_peak <- function(x) max(which(x))
  
  if(nrow(x) > 0){
    first_peak_dt <- x[, 
                     .(excursion            = excursion[first_peak(peak)],
                       excursion_duration   = .N, 
                       single_peak          = FALSE,
                       area_under_excursion = sum(excursion) - (excursion[[1]] + excursion[[.N]])/2, #To get area under curve first and last point are weighted by half ,
                       mean_uptake          = (excursion[first_peak(peak)] - excursion_high)/cumulative_uptake_time[first_peak(peak)],
                       mean_clearance       = NA,
                       max_uptake           = max(uptake_slope[1:(cumulative_uptake_time[first_peak(peak)])]),
                       max_clearance        = min(clearance_slope[(cumulative_uptake_time[first_peak(peak)]):(cumulative_uptake_time[first_peak(peak)] + time_to_next_peak[first_peak(peak)])]),
                       mean_excursion       = mean(excursion),
                       Sample_ID            = Sample_ID[1],
                       ElapsedTime          = ElapsedTime[first_peak(peak)],
                       Light_on             = Light_on[first_peak(peak)],
                       phase_ID             = phase_ID[first_peak(peak)],
                       Day                  = Day[first_peak(peak)],
                       Week                 = Week[first_peak(peak)],
                       ZT                   = ZT[first_peak(peak)],
                       Event_ID             = Event_ID[first_peak(peak)],
                       is_event             = is_event[first_peak(peak)],
                       ZT_event             = ZT_event[first_peak(peak)],
                       peak_type            = "First",
                       n_peaks              = sum(peak),
                       peak_number          = 1,
                       Group                = Group[first_peak(peak)]
                     ), 
                     by = c("excursion_ID")]
  } else {
    first_peak_dt <- empty_kinetics()
  }
  
  par_seq <- function(from, to){
    Map(seq, from = from, to = to)
  }
  
  max_uptake_fn <- function(slope, peak, cumulative_uptake_time, time_since_last_peak){
    idx_peaks <- which(peak)
    idx_peaks <- idx_peaks[-c(1, length(idx_peaks))]
    
    # This bit can be calculated directly from the peaks, but this is much more readable
    ends <- cumulative_uptake_time[idx_peaks]
    starts <- ends - time_since_last_peak[idx_peaks]
    
    idxs <- par_seq(starts, ends)
    helper <- function(idx){
      max(slope[idx])
    }
    unlist(lapply(idxs, helper))
  }
  
  max_clerance_fn <- function(slope, peak, cumulative_uptake_time, time_to_next_peak){
    idx_peaks <- which(peak)
    idx_peaks <- idx_peaks[-c(1, length(idx_peaks))]
    
    # This bit can be calculated directly from the peaks, but this is much more readable
    starts <- cumulative_uptake_time[idx_peaks]
    ends <- starts + time_to_next_peak[idx_peaks]
    
    idxs <- par_seq(starts, ends)
    helper <- function(idx){
      min(slope[idx])
    }
    unlist(lapply(idxs, helper))
  }
  
  get_value_at_peak_fn <- function(value, peak){
    idx_peaks <- which(peak)
    idx_peaks <- idx_peaks[-c(1, length(idx_peaks))]
    value[idx_peaks]
  }
  
  excursion_with_internal_peaks <- stats::na.omit(x[, excursion_ID[peak_type == "Internal"]])
  if (length(excursion_with_internal_peaks) > 0){
    x_nested <- x[excursion_ID %in% excursion_with_internal_peaks]
    internal_peak_dt <- x_nested[, 
                                 .(excursion            = get_value_at_peak_fn(excursion, peak),
                                   excursion_duration   = .N,
                                   single_peak          = FALSE,
                                   area_under_excursion = sum(excursion) - (excursion[[1]] + excursion[[.N]])/2, #To get area under curve first and last point are weighted by half ,
                                   mean_uptake          = NA,
                                   mean_clearance       = NA,
                                   max_uptake           = max_uptake_fn(uptake_slope, peak, cumulative_uptake_time, time_since_last_peak),
                                   max_clearance        = max_clerance_fn(clearance_slope, peak, cumulative_uptake_time, time_to_next_peak),
                                   mean_excursion       = mean(excursion),
                                   Sample_ID            = Sample_ID[1],
                                   ElapsedTime          = get_value_at_peak_fn(ElapsedTime, peak),
                                   Light_on             = get_value_at_peak_fn(Light_on, peak),
                                   phase_ID             = get_value_at_peak_fn(phase_ID, peak),
                                   Day                  = get_value_at_peak_fn(Day, peak),
                                   Week                 = get_value_at_peak_fn(Week, peak),
                                   ZT                   = get_value_at_peak_fn(ZT, peak),
                                   Event_ID             = get_value_at_peak_fn(Event_ID, peak),
                                   is_event             = get_value_at_peak_fn(is_event, peak),
                                   ZT_event             = get_value_at_peak_fn(ZT_event, peak),
                                   peak_type            = "Internal",
                                   n_peaks              = sum(peak),
                                   peak_number          = get_value_at_peak_fn(peak_number, peak) + 1,
                                   Group                = get_value_at_peak_fn(Group, peak)
                                 ), 
                                 by = c("excursion_ID")]
  } else {
    internal_peak_dt <- empty_kinetics()
  }
  
  if (nrow(x) > 0){
    last_peak_dt <- x[, 
                    .(excursion            = excursion[last_peak(peak)],
                      excursion_duration   = .N, 
                      single_peak          = FALSE,
                      area_under_excursion = sum(excursion) - (excursion[[1]] + excursion[[.N]])/2, #To get area under curve first and last point are weighted by half ,
                      mean_uptake          = NA,
                      mean_clearance       = (excursion[last_peak(peak)] - excursion_high)/cumulative_clearence_time[last_peak(peak)],
                      max_uptake           = max(uptake_slope[(cumulative_uptake_time[last_peak(peak)] - time_since_last_peak[last_peak(peak)]):cumulative_uptake_time[last_peak(peak)]]),
                      max_clearance        = min(clearance_slope[cumulative_uptake_time[last_peak(peak)]:(.N)]),
                      mean_excursion       = mean(excursion),
                      Sample_ID            = Sample_ID[1],
                      ElapsedTime          = ElapsedTime[last_peak(peak)],
                      Light_on             = Light_on[last_peak(peak)],
                      phase_ID             = phase_ID[last_peak(peak)],
                      Day                  = Day[last_peak(peak)],
                      Week                 = Week[last_peak(peak)],
                      ZT                   = ZT[last_peak(peak)],
                      Event_ID             = Event_ID[last_peak(peak)],
                      is_event             = is_event[last_peak(peak)],
                      ZT_event             = ZT_event[last_peak(peak)],
                      peak_type            = "Last",
                      n_peaks              = sum(peak),
                      peak_number          = sum(peak),
                      Group                = Group[last_peak(peak)]
                    ), 
                    by = c("excursion_ID")]
    
  } else {
    last_peak_dt <- empty_kinetics()
  }
  list(First = first_peak_dt,
       Internal = internal_peak_dt,
       Last = last_peak_dt)
}

#' Add various peak timers
#'
#' @param x glucose monitoring data.table
#' @param min_peak_duration minimum duration for peak to be counted
#'
#' @return glucose monitoring data.table
#' @import data.table
#' @export
#'
#' @examples
#' \dontrun{
#' to-do
#' }
add_peak_timers <- function(x, min_peak_duration){
  if(is.null(x[["excursion"]])) stop("excursions must be calculated before adding peak timers")
  if(is.null(x[["peak"]])) stop("peak must be tagged before adding peak timers")
  if(all(is.na(x[["excursion"]]))) stop("Sample ", x$Sample_ID[[1]], " has no excursions. Exclude sample or change excursion threshold.")
  # Uptake and clearence timers
  excursion_tmp_grp <- excursion <- cumulative_uptake_time <- cumulative_clearence_time <- 
    tmp_peak <- peak <- time_to_next_peak <- time_since_last_peak <- included <- time_since_last_excluded <- NULL
  
  x[, excursion_tmp_grp:=data.table::rleid(is.na(excursion))]
  x[!is.na(excursion), cumulative_uptake_time:=1:(.N), by = "excursion_tmp_grp"]
  x[!is.na(excursion), cumulative_clearence_time:=(.N):1, by = "excursion_tmp_grp"]
  
  # Only reset peak counter if excursion lasts longer than minPeakDuration
  x[, tmp_peak:=data.table::fifelse(rep(.N >= min_peak_duration, .N), peak, FALSE), by = "excursion_tmp_grp"]
  
  x[, excursion_tmp_grp:=cumsum(tmp_peak)]
  x[, time_to_next_peak:=(.N):1, by = "excursion_tmp_grp"]
  x[, excursion_tmp_grp:=c(0, excursion_tmp_grp[-.N])]
  x[, time_since_last_peak:=1:(.N), by = "excursion_tmp_grp"]
  
  x[, excursion_tmp_grp:=cumsum(!included)]
  x[, time_since_last_excluded:=1:(.N), by = "excursion_tmp_grp"]
  
  x[, excursion_tmp_grp:=NULL]
  x[, tmp_peak:=NULL]
  x[]
}

#' Calculate slopes
#'
#' @param x glucose monitoring data.table
#' @param n_points data points used for slope calculation
#'
#' @return glucose monitoring data.table
#' @import data.table
#' @export
#'
#' @examples
#' \dontrun{
#' to-do
#' }
calc_slopes <- function(x, n_points){
  uptake_slope <- clearance_slope <- NULL
  # Use standard equations to calculate slope
  x_sum <- sum(seq_len(n_points))
  x_2_sum <- sum(seq_len(n_points)^2)
  y_sum <- data.table::frollsum(x$Glucose, n_points, align = "right")
  xy_sum <- data.table::frollapply(x$Glucose, n_points, function(x){sum(x*seq_len(n_points))}, align = "right")
  x[, uptake_slope:=((n_points*xy_sum) - (x_sum*y_sum))/((n_points*x_2_sum) - (x_sum^2))]
  x[, clearance_slope:=c(uptake_slope[-c(seq_len(n_points - 1))], rep(NA, n_points - 1))]
  x[]
}

#' Select excursions to include in kinetics calculations
#'
#' @param x glucose monitoring data.table
#' @param min_peak_duration integer, minimum excursion duration before peak is used for kinetics calculations
#'
#' @return glucose monitoring data.table
#' @import data.table
#' @export
#'
#' @examples
#' \dontrun{
#' to-do
#' }
select_excursions <- function(x, min_peak_duration){
  if(is.null(x[["excursion"]])) stop("excursions must be calculated before selecting excursions")
  if(is.null(x[["time_since_last_excluded"]]) | is.null(x[["time_since_last_peak"]])) stop("peak times must be added before selecting excursions")
  
  grp <- excursion <- included <- time_since_last_excluded <- time_since_last_peak <- 
    time_since_last_peak <- peak <- V1 <- excursion_ID <- NULL
  
  x[, grp:=rleid(is.na(excursion))]
  # We are currently only interested in positive excursions
  x <- x[, included_in_kinetics:=excursion > 0]
  # Remove peaks that contains excluded values
  # Remove first peak after exclusion (could be due to handling etc.)
  excludedExcursions <- x[(included_in_kinetics), any(!included) |
                            any(time_since_last_excluded <= time_since_last_peak) |
                            .N < min_peak_duration |
                            any(is.na(excursion)) |
                            sum(peak) == 0, # small excursions close to a larger excursion may have no peaks
                          by = "grp"][(V1), grp]
  x <- x[(grp %in% excludedExcursions), included_in_kinetics:=FALSE]
  x[, excursion_ID:=NA_integer_]
  
  x[, excursion_ID:=as.integer(rleid(grp))]
  x[, grp:=NULL]
  x[]
}

#' Tag excursion with multiple peaks
#'
#' @param x glucose monitoring data.table
#'
#' @return glucose monitoring data.table
#' @import data.table
#' @export
#'
#' @examples
#' \dontrun{
#' to-do
#' }
tag_multi_peaks <- function(x){
  single_peak <- peak <- peak_type <- NULL
  
  x[, single_peak:=rep(sum(peak) == 1, .N), by = "excursion_ID"]
  x[, peak_type:=NA_character_]
  if (nrow(x[(!single_peak) & (peak)]) > 0){
    x[(!single_peak) & (peak), c("peak_type", "peak_number"):=list(c("First", rep("Internal", sum(peak) - 2), "Last"),
                                                                            c(NA, seq_len(sum(peak) - 2), NA)), by = "excursion_ID"]
  } else {
    x[, c("peak_type", "peak_number"):=list(NA, NA)]
  }
  x[(single_peak) & (peak), c("peak_type", "peak_number"):=list("Single", NA)]
  x[]
}

#' Tag excursion with multiple nadirs
#'
#' @param x glucose monitoring data.table
#'
#' @return glucose monitoring data.table
#' @import data.table
#' @export
#'
#' @examples
#' \dontrun{
#' to-do
#' }
tag_multi_nadirs <- function(x){
  nadir <- nadir_type <- excursion <- NULL
  x[(nadir), nadir_type:=fifelse(.N > 1, "Local", "Global"), by = "excursion_ID"]
  x[(nadir) & nadir_type == "Local", nadir_type:=fifelse(excursion > min(excursion), "Local", "Global"), by = "excursion_ID"]
  x[]
}

#' Secondary filtering of peaks
#'
#' @param x glucose monitoring data.table
#' @param peak_ratio minimum ratio between peak value and minimum obtained since last peak
#'
#' @return glucose monitoring data.table
#' @import data.table
#' @export
#'
#' @examples
#' \dontrun{
#' to-do
#' }
secondary_filtering <- function(x, peak_ratio) {
  all_nadir <- nadir <- all_peak <- peak <- excursion <- peak_type <- NULL
  x[, all_nadir:=nadir]
  
  #x[, nadir:=nadir_type == "Global"]
  x[nadir_type == "Local", nadir:=FALSE]
  
  x[, all_peak:=peak]
  
  helper_fn <- function(excursion, peak) {
    if (sum(peak) < 3) {
      return(peak)
    }
    
    
    intervals <- split(excursion, cumsum(peak))

    # special case where the first time point is a peak (simply add fake interval to be removed below)
    if (peak[[1]]) {
      intervals[["0"]] <- c(0)
    }
    
    # first interval contains values up to first peak, which is always included.
    # last interval contains values after last peak, which is never interesting.
    # second to last interval contains values after last peak, which is always included.
    intervals <- intervals[-c(1, length(intervals) - 1, length(intervals))]
    
    minimums <- unlist(lapply(intervals, min))
    
    idx <- which(peak)
    retain_peak <- (excursion[idx[-c(1, length(idx))]]/minimums) >= peak_ratio
    
    peak[idx[-c(1, length(idx))]] <- retain_peak
    peak
  }
  
  iteration <- 0
  converged <- FALSE
  old_peaks <- x$peak
  
  while (!converged & iteration < 500) {
    x[excursion > 0, peak:=helper_fn(excursion, peak), by = "excursion_ID"]
    
    if (all(old_peaks == x$peak, na.rm = TRUE)) {
      converged <- TRUE
    }
    old_peaks <- x$peak
    iteration <- iteration + 1
  }
  
  x[!peak & all_peak, peak_type:=NA]
  
  x[, peak_number:=NA]
  
  if (nrow(x[(!single_peak) & (peak)]) > 0){
    x[(!single_peak) & (peak), peak_number:=c(NA, seq_len(sum(peak) - 2), NA), by = "excursion_ID"]
  } else {
    x[, c("peak_type", "peak_number"):=list(NA, NA)]
  }
  x[(single_peak) & (peak), c("peak_type", "peak_number"):=list("Single", NA)]
  
  x[]
}

#' get kinetics calculations of a single sample
#'
#' @param x glucose monitoring data.table
#' @param excursion_high increase above baseline before excursion is tagged
#' @param min_peak_duration mimumim peak duration
#' @param datapoints_for_slope datapoints used for slope calculation
#'
#' @return data.table
#' @export
#'
#' @examples
#' \dontrun{
#' to-do
#' }
get_sample_kinetics <- function(x, excursion_high, min_peak_duration, datapoints_for_slope, peak_ratio) {
  # Silence no visible binding warnings
  uptakeSpring <- max_uptake <- excursion <- clearanceSpring <- max_clearance <- NULL
  
  # Sometimes there are no peaks
  x <- tryCatch(add_peak_timers(x, min_peak_duration = min_peak_duration),
           error = function(e) {NULL})
  
  if (is.null(x)) {
    out <- empty_kinetics()
    out$uptakeSpring <- numeric(0)
    out$clearanceSpring <- numeric(0)
    return(out)
  }
  
  x <- calc_slopes(x, n_points = datapoints_for_slope)
  
  x <- select_excursions(x, min_peak_duration = min_peak_duration)
  
  x <- tag_multi_peaks(x)
  
  x <- tag_multi_nadirs(x)
  
  x <- secondary_filtering(x, peak_ratio)
  
  kinetics_single   <- single_peak_kinetics(x[(included_in_kinetics)], excursion_high)
  kinetics_multiple <- multi_peak_kinetics(x[(included_in_kinetics)], excursion_high)
  
  kinetics_multiple$Single <- kinetics_single
  kinetics_all <- do.call(what = "rbind", kinetics_multiple)
  data.table::setkeyv(kinetics_all, "ElapsedTime")
  
  kinetics_all[, uptakeSpring:=max_uptake/excursion]
  kinetics_all[, clearanceSpring:=abs(max_clearance)/excursion]
  
  kinetics_all[]
}

#' Get spring constants
#'
#' @param kinetics_all data.table with kinetics data
#'
#' @return list with spring constants
#' @export
#'
#' @examples
#' \dontrun{
#' to-do
#' }
get_spring_constants <- function(kinetics_all) {
  
  data.table::set(kinetics_all, j = "peak_type", value = factor(kinetics_all$peak_type, levels = c("Single", "First", "Internal", "Last")))
  uptakes <- kinetics_all[, c(
    nObservations = .N,
    as.list(stats::coef(stats::lm(max_uptake ~ excursion, data = .SD))),
    rsquared = summary(stats::lm(max_uptake ~ excursion, data = .SD))$r.squared), 
    .SDcols = c("max_uptake", "excursion"), keyby = c("peak_type", "Sample_ID", "Group")] 
  
  uptakes_both <- kinetics_all[peak_type %in% c("Single", "First"), c(
    peak_type = "Single+First",
    nObservations = .N,
    as.list(stats::coef(stats::lm(max_uptake ~ excursion, data = .SD))),
    rsquared = summary(stats::lm(max_uptake ~ excursion, data = .SD))$r.squared), 
    .SDcols = c("max_uptake", "excursion"), keyby = c("Sample_ID", "Group")] 
  setcolorder(uptakes_both, colnames(uptakes))
  
  uptakes <- data.table::rbindlist(list(uptakes, uptakes_both))
  
  
  clearance <- kinetics_all[, c(
    nObservations = .N,
    as.list(stats::coef(stats::lm(-max_clearance ~ excursion, data = .SD))),
    rsquared = summary(stats::lm(-max_clearance ~ excursion, data = .SD))$r.squared), 
    .SDcols = c("max_clearance", "excursion"), keyby = c("peak_type","Sample_ID", "Group")]
  
  clearance_both <- kinetics_all[peak_type %in% c("Single", "Last"), c(
    peak_type = "Single+Last",
    nObservations = .N,
    as.list(stats::coef(stats::lm(-max_clearance ~ excursion, data = .SD))),
    rsquared = summary(stats::lm(-max_clearance ~ excursion, data = .SD))$r.squared), 
    .SDcols = c("max_clearance", "excursion"), keyby = c("Sample_ID", "Group")]
  
  setcolorder(clearance_both, colnames(clearance))
  
  clearance <- data.table::rbindlist(list(clearance, clearance_both))
  
  list(uptakes = uptakes, clearance = clearance)
}