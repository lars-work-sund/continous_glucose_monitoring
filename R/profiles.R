#' Make breaks
#'
#' @param x glucose monitoring data.table
#' @param low mimimum value
#' @param high maximum value
#' @param step interval between bins
#' @param divide_by number to divide by. If data is sampled every minute and 
#' this is set to 60, then results are returned as per hour.
#'
#' @return list of breaks
make_breaks <- function(x, low, high, step, divide_by){
  x[x < low] <- low
  x[x > high] <- high
  out <- cut(x, seq(low, high, by = step), include.lowest = TRUE, dig.lab = 4)
  out <- table(out) / divide_by
  as.list(out)
}

#' Calculate profiles
#'
#' @param x glucose monitoring data.table
#' @param by_row what should be in the rows (often Light_on)
#' @param by_col what should be in the rows (often Day or Week)
#' @param stat unevaluated expression that will be performed within x
#' @param low lowest value of profile range
#' @param high highest value of profile range 
#' @param min_frac_summaries minimum fraction of data needed before profiling is performed
#' @param subset_expression optional expression used to subset the data (often (peak | nadir))
#' @param update_names logical, should column names have the name of the variable prepended?
#' @param step numeric, interval between bins
#' @param as_percent logical, report results as percent instead of minutes per hour
#'
#' @return data.table
#' @export
#' 
#' @examples
#' \dontrun{
#' to-do
#' }
make_profile <- function(x, by_row, by_col, stat, low, high, step, min_frac_summaries, subset_expression, as_percent, update_names = TRUE, excursion_duration = FALSE) {
  # Silence no visible binding warnings
  tmp_included <- included <- . <- value <- Interval <- NULL
  
  group_by <- c(by_row, by_col)
  
  
  x[, tmp_included:=(included) & !is.na(eval(stat))]
  if (as_percent) {
    normalizing_constant <- 100
  } else {
    normalizing_constant <- 60
  }
  
  if (excursion_duration) {
    if (is.null(x[["excursion_ID"]])) # No peaks
    {
      x[, tmp_included:=(included) & !is.na(eval(stat))]
      return(data.table())
    }
    
    profile <- x[subset_expression, .(n_minutes = .N, 
                 excursion = suppressWarnings(as.numeric(max(filter_max_missing(eval(stat), tmp_included, min_frac_summaries), na.rm = TRUE)))), 
             by = c("excursion_ID", group_by)]
    profile[is.infinite(excursion), excursion:=NA]
    #profile <- profile[!is.na(excursion)]
    profile[excursion < low, excursion:=low]
    profile[excursion > high, excursion:=high]
    profile[, Interval:=cut(excursion, seq(low, high, by = step), include.lowest = TRUE, dig.lab = 4)]
    profile <- profile[, .(value = mean(n_minutes)), by = c(group_by, "Interval")]
    
    
    #tmp <- data.table::dcast(tmp, stats::as.formula(paste(by_row, "+ Interval ~", by_col)), value.var = "n_minutes", fun.aggregate = mean, drop = FALSE)
  } else {
    profile <- x[,
                 make_breaks(filter_max_missing(eval(stat), tmp_included, min_frac_summaries)[filter_max_missing(eval(subset_expression) & tmp_included, tmp_included, min_frac_summaries)], 
                             low = low, 
                             high = high, 
                             step = step, 
                             divide_by = sum(tmp_included)/normalizing_constant),
                 by = group_by]
    profile <- data.table::melt(profile, id.vars = group_by, variable.name = "Interval")
  }
  
  
  values_mean <- profile[, .("__tmp_col__" = "mean", value = mean(value, na.rm = TRUE)), by = c("Interval", by_row)]
  setnames(values_mean, "__tmp_col__", by_col)
  values_sem <- profile[, .("__tmp_col__" = "SEM", value = stats::sd(value, na.rm = TRUE)/length(stats::na.omit(value))), by = c("Interval", by_row)]
  setnames(values_sem, "__tmp_col__", by_col)
  
  long_form <- rbind(profile, values_mean, values_sem)
  
  n_obs <- x[(tmp_included), .N/60, by = group_by]
  n_obs[, Interval:="Hours included"]
  n_obs <- data.table::dcast(n_obs, stats::as.formula(paste(by_row, "+ Interval ~", by_col)), value.var = "V1", drop = FALSE)
  n_obs$mean <- NA_real_
  n_obs$SEM <- NA_real_
  setcolorder(n_obs, c(by_row, "Interval", "mean", "SEM"))
  
  out <- data.table::dcast(long_form, stats::as.formula(paste(by_row, "+ Interval ~", by_col)), drop = FALSE)
  out <- out[, colnames(n_obs), with = FALSE]
  
  #simple sanity check
  #all_ok <- profile[, sum(value), by = group_by]
  #all_ok[V1 == 0, V1:=NaN] #Is situations where some data is included but less than min_frac_summaries. If none are included NaN is already returned
  
  #if(!all(stats::na.omit(all_ok[, V1 - 60] < 0.0001))) {stop("Something went wrong when calculating profiles. There are not 60 minutes for every hour of measurement.")}
  
  x[, tmp_included:=NULL]
  out <- rbind(out, n_obs)
  
  if (update_names){
    old <- unique(profile[[by_col]])
    data.table::setnames(out, as.character(old), paste0(by_col, old), skip_absent = TRUE)
  }
  out
}

#' Get individual and group-wise profiles
#'
#' @param cge cgm_experiment object
#' @param stat unevaluated expression that will be performed within x
#' @param low lowest value of profile range
#' @param high highest value of profile range 
#' @param subset_expression optional expression used to subset the data (often (peak | nadir))
#' @param step numeric, interval between bins
#' @param as_percent logical, report results as percent instead of minutes per hour
#'
#' @return data.frame
#' @export
#'
#' @examples
#' \dontrun{
#' to-do
#' }
get_profiles <- function(cge, stat, low, high, step, as_percent, subset_expression = TRUE, excursion_duration = FALSE) {
  get_sample_kinetics(cge)
  
  stat <- substitute(stat)
  subset_expression <- substitute(subset_expression)
  
  out_part1 <- lapply(cge$data, make_profile, by_row = "Light_on", 
                      by_col = "Day", stat = stat, low = low, high = high, step = step,
                      min_frac_summaries = get_option(cge, "min_frac_summaries"), 
                      as_percent = as_percent, subset_expression = subset_expression,
                      excursion_duration = excursion_duration)
  
  all_data <- data.table::rbindlist(cge$data, fill = TRUE)
  data.table::setkey(all_data, "Group")
  
  groups <- unique(all_data$Group)
  out_part2 <- lapply(groups, function(x) {make_profile(all_data[x], "Light_on", 
                                                        "Week", stat = stat, low, step = step,
                                                        high, get_option(cge, "min_frac_summaries"), 
                                                        as_percent = as_percent,
                                                        subset_expression = subset_expression,
                                                        excursion_duration = excursion_duration)})
  names(out_part2) <- groups
  
  c(out_part1, out_part2)
}
