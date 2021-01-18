#' Helper function for creating an analysis configuration file
#'
#' @param path Path where config file should be located
#' @param sample_names Vector of sample names
#' @param events events as reported in the data file
#'
#' @return NULL
#' @export
#'
#' @examples
#' \dontrun{
#' to-do
#' }
create_config <- function(path, sample_names, events) {
  sheets <- c("README", "sample_grouping", "settings", "events", "all", sample_names)
  
  out <- vector(mode = "list", length = length(sheets))
  names(out) <- sheets
  
  out[["events"]] <- events
  
  out[["sample_grouping"]] <- data.table::data.table(SampleID = sample_names, 
                                   Group = rep("Please Specify", length(sample_names)), 
                                   "Include (Y/N)" = rep("Please Specify", length(sample_names)))
  template <- data.table::data.table("ExclusionStart (DD-MM-YYYY  HH:MM:SS)" = character(0), 
                         "ExclusionEnd (DD-MM-YYYY  HH:MM:SS)" = character(0), 
                         "Notes" = character(0))
  for (i in sample_names){
    out[[i]] <- template
  }
  out[["all"]] <- template
  out[["Time Zones"]] <- data.table::data.table(OlsonNames())
  out[["settings"]] <- data.table::data.table(Parameter = c("light_on", "light_off", "time_zone", "DST", "mgdl_2_mmolL", "max_gap", "baseline_window", "max_missing_baseline", "excursion_low", "excursion_high", "max_min_window", "min_peak_duration", "datapoints_for_slope", "min_frac_summaries"), 
                                  value = c("", "", "Please Specify", "independent", "Please specify", "5", "1440", "600", "-1", "1", "17", "3", "4", "0.5"), 
                                  notes = c("HH:MM", "HH:MM", "Change to sync if light cycle follows clock-time after DST", "Script assumes mmol/L, if values are mg/dl set to Y", "Maximum gap to interpolate", "Number of datapoints used for baseline calculations", "Maximum number of missing datapoints where baseline calculation is performed", "Excursion (in mmol/L) needed to flag nadir", "Excursion (in mmol/L) needed to flag peak", "Interval to search for local peaks and nadirs (must be odd). Scale is the same as measurement interval, 17 works for minutes, change as needed.", "Minimum duration before peak is used for kinetics calculations", "Number of datapoints used when calculating uptake and clearance slopes. Remember to rescale if data is not in one minute interval", "Minimum fraction of observations that needs to be included in order to calculate summary statistics"))
  openxlsx::write.xlsx(out, file = path)
  invisible(NULL)
}



#' Read settings
#'
#' @param file path to config file
#'
#' @return list of settings
#'
#' @examples
#' \dontrun{
#' to-do
#' }
read_settings <- function(file){
  excel_formatted <- readxl::read_xlsx(sheet = "settings", 
                                path = file, 
                                range = readxl::cell_cols(c("A", "B")), 
                                col_types = c("text", "text"))
  
  settings <- as.list(excel_formatted$value)
  names(settings) <- excel_formatted$Parameter
  
  settings[["light_on"]] <- lubridate::hms(
    format(openxlsx::convertToDateTime(
      as.numeric(settings[["light_on"]])), "%H:%M:%S"))
  settings[["light_off"]] <- lubridate::hms(
    format(openxlsx::convertToDateTime(
      as.numeric(settings[["light_off"]])), "%H:%M:%S"))
  settings[["DST"]] <- tolower(settings[["DST"]])
  settings[["mgdl_2_mmolL"]] <- tolower(settings[["mgdl_2_mmolL"]])
  settings[["max_gap"]] <- as.integer(settings[["max_gap"]])
  settings[["baseline_window"]] <- as.integer(settings[["baseline_window"]])
  settings[["max_missing_baseline"]] <- as.integer(settings[["max_missing_baseline"]])
  settings[["excursion_low"]] <- as.integer(settings[["excursion_low"]])
  settings[["excursion_high"]] <- as.integer(settings[["excursion_high"]])
  settings[["max_min_window"]] <- as.integer(settings[["max_min_window"]])
  settings[["min_peak_duration"]] <- as.integer(settings[["min_peak_duration"]])
  settings[["datapoints_for_slope"]] <- as.integer(settings[["datapoints_for_slope"]])
  settings[["min_frac_summaries"]] <- as.numeric(settings[["min_frac_summaries"]])
  settings
}


#' Read exclusions
#'
#' @param sheet name of sheet to read
#' @param file path to config file
#'
#' @return list of exclusions
#'
#' @examples
#' \dontrun{
#' to-do
#' }
read_exclusion <- function(sheet, file) {
  exclusions <- readxl::read_xlsx(sheet,
                                  path = file, 
                                  range = readxl::cell_cols(c("A", "B")),
                                  col_types = c("date", "date")
  )
}

#' Read sample groupings
#'
#' @param file path to config file
#'
#' @return list of groupings
#'
#' @examples
#' \dontrun{
#' to-do
#' }
read_groupings <- function(file) {
  readxl::read_xlsx(sheet = "sample_grouping", path = file)
}

#' Read configurations from config file
#'
#' @param file path to configurations file
#' @param samples names of samples
#'
#' @return cgm_config object
#' @export
#'
#' @examples
#' \dontrun{
#' to-do
#' }
read_config <- function(file, samples) {
  settings <- read_settings(file)
  group <- read_groupings(file)
  
  samples <- c("all", samples)
  exclusions <- lapply(samples, read_exclusion, file = file)
  names(exclusions) <- samples
  
  config(settings = settings, groupings = group, exclusions = exclusions)
}

#' Read continuous glucose monitoring data
#'
#' @param file file with continuous glucose monitoring data
#' @param pattern regex to recognize sheets containing data
#' @param parallel logical, should data be loaded in parallel?
#'
#' @importFrom foreach %dopar%
#' 
#' @return cgm_glucose_data object
#' @export
#'
#' @examples
#' \dontrun{
#' to-do
#' }
read_data <- function(file, pattern = "Parameters", parallel = FALSE) {
  sheets <- stringr::str_subset(readxl::excel_sheets(file), pattern)
  
  if (parallel) {
    cores <- parallel::detectCores()
    cl <- parallel::makeCluster(cores[1] - 1)
    doParallel::registerDoParallel(cl)
    i <- NULL # To silence warning that i is undefined
    glc_data <- foreach::foreach(i = sheets, .packages = "continousGlucoseMonitoring") %dopar% 
      data_reader(i, file)
    
    parallel::stopCluster(cl)
  } else {
    glc_data <- lapply(sheets, data_reader, file)
  }
  names(glc_data) <- sheets
  glucose_data(glc_data)
}

#' data reading helper
#'
#' @param sheet_name name of sheet to read
#' @param file path to data file
#'
#' @return data.table
#'
#' @examples
#' \dontrun{
#' to-do
#' }
data_reader <- function(sheet_name, file)
{
  raw_data <- readxl::read_xlsx(file, 
                       sheet = sheet_name, 
                       range = readxl::cell_cols(c("A", "M")),
                       col_types = c("date", 
                                     "skip", 
                                     "skip", 
                                     "text", 
                                     "text", 
                                     "skip", 
                                     "skip", 
                                     "numeric",
                                     "numeric",
                                     "numeric",
                                     "numeric",
                                     "numeric",
                                     "numeric"
                                     )
  )
  data.table::setDT(raw_data)
  
  data.table::setnames(raw_data, c("Gavg(mmol/L):Glucose", "Gavg(mg/dL):Glucose", "T_Mean(Celsius):Temp", "T_Mean(Celsius):Temperat", "A_TA(Counts):Activity", "A_TA(Counts.s):Activity"), 
                       c("Glucose", "Glucose", "Temperature", "Temperature", "Activity", "Activity"), skip_absent = TRUE)
  
  unwanted_columns <- setdiff(colnames(raw_data), c("Date", "ElapsedTime", "Event", "Glucose", "Temperature", "Activity"))
  data.table::set(raw_data, j = unwanted_columns, value = NULL)
  data.table::set(raw_data, j = "Glucose", value = as.numeric(raw_data$Glucose))
  data.table::set(raw_data, j = "Sample_ID", value = sheet_name)
  data.table::set(raw_data, j = "ElapsedTime", value = lubridate::as.duration(lubridate::hms(raw_data$ElapsedTime)))
  
  raw_data
}

#' Helper script for preparing for data analysis
#' 
#' @description Helper script to prepare data generated using this data generation platform.
#'
#' @param data_file path to data file
#' @param configuration_file path to configuration file, if missing one will be created
#' @param pattern Regex used to recognize data sheets
#' @param parallel logical, should data be loaded in parallel?
#'
#' @return cgm_experiment object
#' @export
#'
#' @examples 
#' \dontrun{
#' to-do
#' }
prepare_experiment <- function(data_file, configuration_file, pattern = "Parameters", parallel = TRUE) {
  if (!file.exists(data_file)) stop(paste(data_file, "not found"))
  samples <- stringr::str_subset(readxl::excel_sheets(data_file), pattern)
  if (!file.exists(configuration_file)) {
    
    events <- readxl::read_xlsx(data_file, sheet = "Events")

    create_config(path = configuration_file, sample_names = samples, events = events)
    stop("A settings file was not found, so one was generated, please fill out with subject inclusion/exclusion, grouping information, light cycle and excluded timepoints and re-run script.")
  }
  
  configuration <- read_config(configuration_file, samples)
  glucose_data <- read_data(data_file, pattern = pattern, parallel = parallel)
  
  cgm_experiment(glucose_data, configuration)
}
