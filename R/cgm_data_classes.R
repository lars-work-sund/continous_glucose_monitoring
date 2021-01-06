#' Glucose data constructor
#'
#' @param x list of data.tables with (potentially a subset of) the columns 
#' needed for the full analysis
#'
#' @return a cgm_glucose_data object
#' @export
#'
#' @examples
new_glucose_data <- function(x = data.table::data.table()) {
  vctrs::new_list_of(x, ptype = data.table::data.table(), 
                     class = "cgm_glucose_data")
}

#' Glucose data validator
#'
#' @param x a cgm_glucose_data object
#'
#' @return a cgm_glucose_data object
#' @export
#'
#' @examples
validate_glucose_data <- function(x) {
  # It is a valid data structure if the column names and types match
  
  is_valid_format <- function(x, sample) {
    prototype <- list(
      "ElapsedTime" = "Duration", 
      "Date" = c("POSIXct", "POSIXt"),
      "Event" = "character", 
      "Glucose" = "numeric", 
      "Temperature" = "numeric", 
      "Activity" = "numeric", 
      "Sample_ID" = "character", 
      "Light_on" = "logical", 
      "Phase_ID" = "numeric", 
      "Day" = "numeric", 
      "Week" = "numeric", 
      "ZT" = "numeric", 
      "ZT_exact" = "numeric", 
      "Included" = "logical", 
      "Group" = "character", 
      "Baseline" = "numeric", 
      "Excursion" = "numeric", 
      "Peak" = "logical", 
      "Nadir" = "logical"
    )
    
    # the prototype defines the expected columns and their types
    proper_names <- all(names(x) %in% names(prototype))
    if(!proper_names) {
      stop(sprintf("%s: %s not recognized column name. Must be one of %s",
                   sample, 
                   paste(setdiff(names(x), names(prototype)), collapse = ", "), 
                   paste(setdiff(names(prototype), names(x)), collapse = ", ")), 
           call. = FALSE)
    }
    
    proper_class <- unlist(Map(inherits, x, prototype[colnames(x)]))
    
    if (!all(proper_class)) {
      stop(sprintf("\n%s: %s does not inherit from %s", 
                   sample, 
                   colnames(x)[!proper_class], 
                   prototype[colnames(x)[!proper_class]]),
           call. = FALSE)
    }
  }
  for (i in names(x)) {
    is_valid_format(x = x[[i]], sample = i)
  }
  unnamed <- is.null(names(x))
  if (unnamed) stop("Samples must be named")
  
  duplicate_names <- anyDuplicated(names(x))
  if (duplicate_names) stop("Samples must be uniquely named")
  
  x
}

#' Glucose data helper 
#'
#' @param x list of data.tables with (potentially a subset of) the columns 
#' needed for the full analysis
#'
#' @return a cgm_glucose_data object 
#' @export 
#'
#' @examples
glucose_data <- function(x) {
  validate_glucose_data(new_glucose_data(x))
}



#' Configuration constructor
#'
#' @param settings data.frame with experimental settings
#' @param groupings data.frame with sample groupings
#' @param exclusions list of data.frames with start and end of excluded timepoints
#'
#' @return configuration class list
#' @export
#'
#' @examples
new_config <- function(settings, groupings, exclusions) {
  out <- list(
    settings = settings,
    groupings = groupings,
    exclusions = exclusions
  )
  class(out) <- "cgm_config"
  out
}

#' Settings validator
#'
#' @param x settings part of a configuration object
#'
#' @return
#'
#' @examples
validate_settings <- function(x) {
  # Light on / off
  proper_time_format <- function(x, parameter) {
    if(!inherits(x, "Period")) stop(paste(parameter, "does not inherit from Period"))
    
    h <- lubridate::hour(x)
    m <- lubridate::minute(x)
    if(h < 0 | h > 23) stop(paste(parameter, "hour must be between 0 and 23"))
    if(m < 0 | m > 59) stop(paste(parameter, "minute must be between 0 and 59"))
  }
  proper_time_format(x$settings[["light_on"]], "light_on")
  proper_time_format(x$settings[["light_off"]], "light_off")
  
  # DST
  if(!(x$settings$DST %in% c("independent", "sync"))) stop("DST must be either independent or sync")
  
  # mgdl_2_mmolL
  if(!(x$settings$mgdl_2_mmolL %in% c("y", "n"))) stop("mgdl_2_mmolL must be either y or n")
  
  # max_gap
  if(x$settings$max_gap < 0) stop("max_gap must be 0 or greater")

  # baseline_window
  if(x$settings$baseline_window < 0) stop("baseline_window must be 0 or greater")
  
  # excursion_low & excursion_high
  if(x$settings$excursion_low > x$settings$excursion_high) stop("excursion_low must be less than excursion_high")
  
  # max_min_window
  if(x$settings$max_min_window < 3 | (x$settings$max_min_window %% 2 != 1)) stop("max_min_window must be an odd integer greater than 1")
  
  # datapoints_for_slope
  if(x$settings$datapoints_for_slope < 0) stop("datapoints_for_slope must be 0 or greater")
  
  # datapoints_for_slope
  if(x$settings$min_frac_summaries < 0 | x$settings$min_frac_summaries > 1) stop("min_frac_summaries must be between 0 and 1")
  
  x
}

#' Groupings validator
#'
#' @param x groupings part of settings
#'
#' @return
#'
#' @examples
validate_grouping <- function(x) {
  
  expected_colnames <- c("SampleID", "Group", "Include (Y/N)")
  
  proper_names <- all(colnames(x$groupings) %in% expected_colnames)
  if(!proper_names) {
    stop(sprintf("%s not recognized column name. Must be one of %s\n",
                 paste(setdiff(colnames(x$groupings), expected_colnames), collapse = ", "), 
                 paste(setdiff(expected_colnames, colnames(x$groupings)), collapse = ", ")))
  }
  
  if(class(x$groupings[["SampleID"]]) != "character") stop("SampleID must be character")
  if(!class(x$groupings[["Group"]]) %in% c("character", "factor")) stop("Group must be character or factor")
  if(!all(x$groupings[["Include (Y/N)"]] %in% c("y", "n"))) stop("Include (Y/N) must be y or n")
  
  x
}

#' Exclusions validator
#'
#' @param x exclusions part of settings
#'
#' @return
#'
#' @examples
validate_exclusions <- function(x) {
  valid_exclusion <- function(y, sample) {
    
    expected_colnames <- c("Start", "End")
    proper_names <- all(colnames(y) %in% expected_colnames)
    if(!proper_names) {
      stop(sprintf("%s not recognized column name. Must be one of %s\n",
                   paste(setdiff(colnames(y), expected_colnames), collapse = ", "), 
                   paste(setdiff(expected_colnames, colnames(y)), collapse = ", ")))
    }
    
    start_ok <- lubridate::is.instant(y$Start)
    if(!start_ok) stop(paste0(sample,": Start is not recognized as a timepoint"))
    
    end_ok <- lubridate::is.instant(y$End)
    if(!end_ok) stop(paste0(sample, ": End is not recognized as a timepoint"))
    
    proper_order <- all(y$Start < y$End)
    if(!proper_order) stop("Some exclusions stop before they start")
    
    invisible(y)
  }
  
  if(is.null(names(x$exclusions))) stop("Exclusions must be named")
  Map(valid_exclusion, x$exclusions, names(x$exclusions))
  
  x
}

#' Configuration validator
#'
#' @param x configuration object
#'
#' @return
#' @export
#'
#' @examples
validate_config <- function(x) {
  validate_settings(x)
  validate_grouping(x)
  validate_exclusions(x)
  x
}

#' Configuration helper
#'
#' @param settings data.frame with experimental settings
#' @param groupings data.frame with sample groupings
#' @param exclusions list of data.frames with start and end of excluded timepoints
#'
#' @return
#' @export
#'
#' @examples
config <- function(settings, groupings, exclusions) {
  tryCatch(expr = {
    exclusions <- lapply(exclusions, dplyr::rename, 
         "Start" = "ExclusionStart (DD-MM-YYYY  HH:MM:SS)",
         "End" = "ExclusionEnd (DD-MM-YYYY  HH:MM:SS)"
         )
  }, error = function(e){invisible(NULL)}
  )
  if (!is.null(groupings$`Include (Y/N)`)) {
    groupings$`Include (Y/N)` <- tolower((groupings$`Include (Y/N)`))
  }
  if (!is.null(groupings$Group) & !(is.factor(groupings$Group) | is.character(groupings$Group))) {
    groupings$Group <- as.character(groupings$Group)
  }
  if (!is.null(settings$mgdl_2_mmolL)) {
    settings$mgdl_2_mmolL <- tolower(settings$mgdl_2_mmolL)
  }
  
  validate_config(new_config(settings = settings, 
                             groupings = groupings, 
                             exclusions = exclusions))
}

#' Continous glucose monitoring experiment constructor
#'
#' @param glucose_data gcm_data object
#' @param config gcm_config object
#'
#' @return cgm_experiment object
#' @export
#'
#' @examples
new_cgm_experiment <- function(glucose_data, config) {
  out <- list(data = glucose_data,
              config = config)
  class(out) <- "cgm_experiment"
  out
}

#' Continous glucose monitoring experiment validator
#'
#' @param x cgm_experiment object
#'
#' @return
#' @export
#'
#' @examples
validate_cgm_experiment <- function(x) {
  extra_samples <- setdiff(names(x$data), names(x$config$exclusions))
  if (length(extra_samples) != 0) stop(paste("Samples not present in config detected: ", paste(extra_samples, collapse = ", ")))
  
  extra_config <- setdiff(names(x$config$exclusions), c(names(x$data), "all"))
  if (length(extra_config) != 0) stop(paste("Samples present in config but not in data: ", paste(extra_config, collapse = ", ")))
  
  if (!("all" %in% names(x$config$exclusions))) stop("'all' is missing from config")
  
  proper_specific <- function(data, exclusion, sample) {
    exclusion_interval <- lubridate::interval(exclusion$Start, exclusion$End)
    all_measured <- lubridate::interval(min(data$Date), max(data$Date))
    specific_overlaps <- all(lubridate::int_overlaps(exclusion_interval, all_measured))
    
    if (!specific_overlaps) stop(paste0(sample, ": Not all exclusions overlap data"))
  }
  Map(proper_specific, x$data, x$config$exclusions[names(x$data)], names(x$data))
  
  proper_all <- function(data) {
    exclusion_interval <- lubridate::interval(x$config$exclusions$all$Start, x$config$exclusions$all$End)
    all_measured <- lubridate::interval(min(data$Date), max(data$Date))
    lubridate::int_overlaps(exclusion_interval, all_measured)
  }
  all_overlap <- lapply(x$data, proper_all)
  all_overlap <- all(Reduce(`&`, all_overlap))
  if (!all_overlap) stop("Some exclusions in 'all' overlaps no data")
  
  
  grouped_samples <- x$config$groupings$SampleID
  present_samples <- names(x$data)
  
  all_present <- all(grouped_samples %in% present_samples)
  all_grouped <- all(present_samples %in% grouped_samples)
  
  if(!all_grouped) {
    stop(sprintf("%s Not present in sample grouping\n",
                 paste(setdiff(present_samples, grouped_samples), collapse = ", ")))
  }
  if(!all_present) {
    stop(sprintf("%s present in sample grouping but not in data\n",
                 paste(setdiff(grouped_samples, present_samples), collapse = ", ")))
  }
  
  all_samples_grouped <- all(x$config$groupings$SampleID %in% names(x$data))
  if (!all_samples_grouped) stop("Some samples are missing from the grouping")
  
  x
}

#' Continous glucose monitoring experiment helper
#'
#' @param glucose_data gcm_data object
#' @param config gcm_config object
#'
#' @return
#' @export
#'
#' @examples
cgm_experiment <- function(glucose_data, config) {
  validate_cgm_experiment(new_cgm_experiment(glucose_data, config))
}

#' Print Values
#'
#' @param x cgm_experiment object
#'
#' @return
#' @export
#'
#' @examples
print.cgm_experiment <- function(x) {
  nSamples <- length(x$data)
  nObs <- range(unlist(lapply(x$data, nrow)))
  nExclusions <- sum(unlist(lapply(x$config$exclusions, nrow)))
  cat(sprintf("Continous Glucose Monitoring experiment\nSamples: %d\nTimepoints: %d to %d\nExclusions: %d", nSamples, nObs[1], nObs[2], nExclusions))
  invisible(x)
}

