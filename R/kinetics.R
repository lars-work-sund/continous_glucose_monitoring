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
    excursionId        = integer(0),
    excursion          = numeric(0),
    excursionDuration  = integer(0), 
    singlePeak         = logical(0),
    areaUnderExcursion = numeric(0),
    meanUptake         = numeric(0),
    meanClearance      = numeric(0),
    maxUptake          = numeric(0),
    maxClearance       = numeric(0),
    Sample_ID          = character(0),
    ElapsedTime        = lubridate::duration(),
    Light_on           = logical(0),
    PhaseId            = numeric(0),
    Day                = numeric(0),
    Week               = numeric(0),
    ZT                 = numeric(0),
    Event_ID           = integer(0),
    is_event           = logical(0),
    ZT_event           = numeric(0),
    nestedPeakType     = character(0),
    nPeaks             = integer(0),
    peakNumber         = numeric(0),
    Group              = character(0)
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
  singlePeak <- . <- excursion <- peak <- cumulativeUptakeTime <- 
    cumulativeClearenceTime <- uptakeSlope <- clearanceSlope <- Sample_ID <- 
    ElapsedTime <- Light_on <- PhaseId <- Day <- Week <- ZT <- Group <- NULL
  
  if (any(x$peak & x$singlePeak)) {
    x[(singlePeak), 
      .(excursion          = excursion[peak],
        excursionDuration  = .N, 
        singlePeak         = TRUE,
        areaUnderExcursion = sum(excursion) - (excursion[[1]] + excursion[[.N]])/2, #To get area under curve first and last point are weighted by half 
        meanUptake         = (excursion[peak] - excursion_high)/cumulativeUptakeTime[peak],
        meanClearance      = (excursion[peak] - excursion_high)/cumulativeClearenceTime[peak],
        maxUptake          = max(uptakeSlope),
        maxClearance       = min(clearanceSlope),
        Sample_ID          = Sample_ID[1],
        ElapsedTime        = ElapsedTime[peak],
        Light_on           = Light_on[peak],
        PhaseId            = PhaseId[peak],
        Day                = Day[peak],
        Week               = Week[peak],
        ZT                 = ZT[peak],
        Event_ID           = Event_ID[peak],
        is_event           = is_event[peak],
        ZT_event           = ZT_event[peak],
        nestedPeakType     = "Single",
        nPeaks             = 1,
        peakNumber         = 1,
        Group              = Group[peak]
      ), 
      by = c("excursionId")]
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
  singlePeak <- . <- excursion <- peak <- cumulativeUptakeTime <- uptakeSlope <- 
    clearanceSlope <- timeToNextPeak <- Sample_ID <- ElapsedTime <- Light_on <- 
    PhaseId <- Day <- Week <- ZT <- Group <- excursionId <- nestedPeakType <- 
    timeSinceLastPeak <- nestedPeakNumber <- cumulativeClearenceTime <- NULL
  
  x <- x[(!singlePeak)]
  
  first_peak <- function(x) min(which(x))
  last_peak <- function(x) max(which(x))
  
  if(nrow(x) > 0){
    first_peak_dt <- x[, 
                     .(excursion          = excursion[first_peak(peak)],
                       excursionDuration  = .N, 
                       singlePeak         = FALSE,
                       areaUnderExcursion = NA,
                       meanUptake         = (excursion[first_peak(peak)] - excursion_high)/cumulativeUptakeTime[first_peak(peak)],
                       meanClearance      = NA,
                       maxUptake          = max(uptakeSlope[1:(cumulativeUptakeTime[first_peak(peak)])]),
                       maxClearance       = min(clearanceSlope[(cumulativeUptakeTime[first_peak(peak)]):(cumulativeUptakeTime[first_peak(peak)] + timeToNextPeak[first_peak(peak)])]),
                       Sample_ID          = Sample_ID[1],
                       ElapsedTime        = ElapsedTime[first_peak(peak)],
                       Light_on           = Light_on[first_peak(peak)],
                       PhaseId            = PhaseId[first_peak(peak)],
                       Day                = Day[first_peak(peak)],
                       Week               = Week[first_peak(peak)],
                       ZT                 = ZT[first_peak(peak)],
                       Event_ID           = Event_ID[first_peak(peak)],
                       is_event           = is_event[first_peak(peak)],
                       ZT_event           = ZT_event[first_peak(peak)],
                       nestedPeakType     = "First",
                       nPeaks             = sum(peak),
                       peakNumber         = 1,
                       Group              = Group[first_peak(peak)]
                     ), 
                     by = c("excursionId")]
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
  
  excursion_with_internal_peaks <- stats::na.omit(x[, excursionId[nestedPeakType == "Internal"]])
  if (length(excursion_with_internal_peaks) > 0){
    x_nested <- x[excursionId %in% excursion_with_internal_peaks]
    internal_peak_dt <- x_nested[, 
                                 .(excursion          = get_value_at_peak_fn(excursion, peak),
                                   excursionDuration  = .N,
                                   singlePeak         = FALSE,
                                   areaUnderExcursion = NA,
                                   meanUptake         = NA,
                                   meanClearance      = NA,
                                   maxUptake          = max_uptake_fn(uptakeSlope, peak, cumulativeUptakeTime, timeSinceLastPeak),
                                   maxClearance       = max_clerance_fn(clearanceSlope, peak, cumulativeUptakeTime, timeToNextPeak),
                                   Sample_ID          = Sample_ID[1],
                                   ElapsedTime        = get_value_at_peak_fn(ElapsedTime, peak),
                                   Light_on           = get_value_at_peak_fn(Light_on, peak),
                                   PhaseId            = get_value_at_peak_fn(PhaseId, peak),
                                   Day                = get_value_at_peak_fn(Day, peak),
                                   Week               = get_value_at_peak_fn(Week, peak),
                                   ZT                 = get_value_at_peak_fn(ZT, peak),
                                   Event_ID           = get_value_at_peak_fn(Event_ID, peak),
                                   is_event           = get_value_at_peak_fn(is_event, peak),
                                   ZT_event           = get_value_at_peak_fn(ZT_event, peak),
                                   nestedPeakType     = "Internal",
                                   nPeaks             = sum(peak),
                                   peakNumber         = get_value_at_peak_fn(nestedPeakNumber, peak) + 1,
                                   Group              = get_value_at_peak_fn(Group, peak)
                                 ), 
                                 by = c("excursionId")]
  } else {
    internal_peak_dt <- empty_kinetics()
  }
  
  if (nrow(x) > 0){
    last_peak_dt <- x[, 
                    .(excursion          = excursion[last_peak(peak)],
                      excursionDuration  = .N, 
                      singlePeak         = FALSE,
                      areaUnderExcursion = NA,
                      meanUptake         = NA,
                      meanClearance      = (excursion[last_peak(peak)] - excursion_high)/cumulativeClearenceTime[last_peak(peak)],
                      maxUptake          = max(uptakeSlope[(cumulativeUptakeTime[last_peak(peak)] - timeSinceLastPeak[last_peak(peak)]):cumulativeUptakeTime[last_peak(peak)]]),
                      maxClearance       = min(clearanceSlope[cumulativeUptakeTime[last_peak(peak)]:(.N)]),
                      Sample_ID          = Sample_ID[1],
                      ElapsedTime        = ElapsedTime[last_peak(peak)],
                      Light_on           = Light_on[last_peak(peak)],
                      PhaseId            = PhaseId[last_peak(peak)],
                      Day                = Day[last_peak(peak)],
                      Week               = Week[last_peak(peak)],
                      ZT                 = ZT[last_peak(peak)],
                      Event_ID           = Event_ID[last_peak(peak)],
                      is_event           = is_event[last_peak(peak)],
                      ZT_event           = ZT_event[last_peak(peak)],
                      nestedPeakType     = "Last",
                      nPeaks             = sum(peak),
                      peakNumber         = sum(peak),
                      Group              = Group[last_peak(peak)]
                    ), 
                    by = c("excursionId")]
    
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
  excursion_tmp_grp <- excursion <- cumulativeUptakeTime <- cumulativeClearenceTime <- 
    tmp_peak <- peak <- timeToNextPeak <- timeSinceLastPeak <- included <- timeSinceLastExcluded <- NULL
  
  x[, excursion_tmp_grp:=data.table::rleid(is.na(excursion))]
  x[!is.na(excursion), cumulativeUptakeTime:=1:(.N), by = "excursion_tmp_grp"]
  x[!is.na(excursion), cumulativeClearenceTime:=(.N):1, by = "excursion_tmp_grp"]
  
  # Only reset peak counter if excursion lasts longer than minPeakDuration
  x[, tmp_peak:=data.table::fifelse(rep(.N >= min_peak_duration, .N), peak, FALSE), by = "excursion_tmp_grp"]
  
  x[, excursion_tmp_grp:=cumsum(tmp_peak)]
  x[, timeToNextPeak:=(.N):1, by = "excursion_tmp_grp"]
  x[, excursion_tmp_grp:=c(0, excursion_tmp_grp[-.N])]
  x[, timeSinceLastPeak:=1:(.N), by = "excursion_tmp_grp"]
  
  x[, excursion_tmp_grp:=cumsum(!included)]
  x[, timeSinceLastExcluded:=1:(.N), by = "excursion_tmp_grp"]
  
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
  uptakeSlope <- clearanceSlope <- NULL
  # Use standard equations to calculate slope
  x_sum <- sum(seq_len(n_points))
  x_2_sum <- sum(seq_len(n_points)^2)
  y_sum <- data.table::frollsum(x$Glucose, n_points, align = "right")
  xy_sum <- data.table::frollapply(x$Glucose, n_points, function(x){sum(x*seq_len(n_points))}, align = "right")
  x[, uptakeSlope:=((n_points*xy_sum) - (x_sum*y_sum))/((n_points*x_2_sum) - (x_sum^2))]
  x[, clearanceSlope:=c(uptakeSlope[-c(seq_len(n_points - 1))], rep(NA, n_points - 1))]
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
  if(is.null(x[["timeSinceLastExcluded"]]) | is.null(x[["timeSinceLastPeak"]])) stop("peak times must be added before selecting excursions")
  
  grp <- excursion <- included <- timeSinceLastExcluded <- timeSinceLastPeak <- 
    timeSinceLastPeak <- peak <- V1 <- excursionId <- NULL
  
  x[, grp:=rleid(is.na(excursion))]
  # We are currently only interested in positive excursions
  x <- x[excursion > 0]
  # Remove peaks that contains excluded values
  # Remove first peak after exclusion (could be due to handling etc.)
  excludedExcursions <- x[, any(!included) |
                            any(timeSinceLastExcluded <= timeSinceLastPeak) |
                            .N < min_peak_duration |
                            any(is.na(excursion)) |
                            sum(peak) == 0, # small excursions close to a larger excursion may have no peaks
                          by = "grp"][(V1), grp]
  x <- x[!(grp %in% excludedExcursions)]
  x[, excursionId:=as.integer(rleid(grp))]
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
  singlePeak <- peak <- nestedPeakType <- NULL
  
  x[, singlePeak:=rep(sum(peak) == 1, .N), by = "excursionId"]
  x[, nestedPeakType:=NA_character_]
  if (nrow(x[(!singlePeak) & (peak)]) > 0){
    x[(!singlePeak) & (peak), c("nestedPeakType", "nestedPeakNumber"):=list(c("First", rep("Internal", sum(peak) - 2), "Last"),
                                                                            c(NA, seq_len(sum(peak) - 2), NA)), by = "excursionId"]
  } else {
    x[, c("nestedPeakType", "nestedPeakNumber"):=list(NA, NA)]
  }
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
get_sample_kinetics <- function(x, excursion_high, min_peak_duration, datapoints_for_slope) {
  # Silence no visible binding warnings
  uptakeSpring <- maxUptake <- excursion <- clearanceSpring <- maxClearance <- NULL
  
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
  
  kinetics_single   <- single_peak_kinetics(x, excursion_high)
  kinetics_multiple <- multi_peak_kinetics(x, excursion_high)
  
  kinetics_multiple$Single <- kinetics_single
  kinetics_all <- do.call(what = "rbind", kinetics_multiple)
  data.table::setkeyv(kinetics_all, "ElapsedTime")
  
  kinetics_all[, uptakeSpring:=maxUptake/excursion]
  kinetics_all[, clearanceSpring:=abs(maxClearance)/excursion]
  
  kinetics_all[]
}
