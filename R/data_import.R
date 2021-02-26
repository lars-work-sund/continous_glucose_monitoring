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
  sheets <- c("README", "sample_grouping", "settings", "Time Zones", "events", "all", sample_names)
  
  out <- vector(mode = "list", length = length(sheets))
  names(out) <- sheets
  
  out[["events"]] <- events
  
  out[["sample_grouping"]] <- data.table::data.table(SampleID = sample_names, 
                                                     Group = rep("Please Specify", length(sample_names)), 
                                                     "Include (Y/N)" = rep("y", length(sample_names)),
                                                     Alias = rep("", length(sample_names)))
  template <- data.table::data.table("ExclusionStart (DD-MM-YYYY  HH:MM:SS)" = character(0), 
                                     "ExclusionEnd (DD-MM-YYYY  HH:MM:SS)" = character(0), 
                                     "Notes" = character(0))
  out[["all"]] <- template
  out[["Time Zones"]] <- data.table::data.table(OlsonNames())
  
  out[["settings"]] <- data.table::fread(system.file("extdata", "config_file.csv", package = "continousGlucoseMonitoring", mustWork = TRUE))
  
  
  for (i in sample_names){
    out[[i]] <- template
  }
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
  settings[["time_zone"]] <- settings[["time_zone"]] # for consistency
  settings[["activity_col"]] <- as.character(settings[["activity_col"]])
  settings[["temperature_col"]] <- as.character(settings[["temperature_col"]])
  settings[["glucose_col"]] <- as.character(settings[["glucose_col"]])
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
  settings[["profile_glucose_bins"]] <- as.numeric(stringr::str_split(as.character(settings[["profile_glucose_bins"]]), ";")[[1]])
  settings[["profile_peak_iso_bins"]] <- as.numeric(stringr::str_split(as.character(settings[["profile_peak_iso_bins"]]), ";")[[1]])
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

#' data reading helper
#'
#' @param sheet_name name of sheet to read
#' @param file path to data file
#' @param glucose_col character, column name for glucose data (typically H)
#' @param temperature_col character, column name for temperature data 
#' @param activity_col character, column name for activity data 
#'
#' @return data.table
#'
#' @examples
#' \dontrun{
#' to-do
#' }
#' 
data_reader <- function(sheet_name, file, glucose_col, temperature_col, activity_col)
{
  # Temperature and activity columns are optional
  if (missing(temperature_col)) temperature_col <- ""
  if (missing(activity_col)) activity_col <- ""
  
  # Cells to be read have to be specified as a contiguous block of columns, fortunately
  # helper functions can figure out exactly which columns need to be loaded.
  col_range <- readxl::cell_cols(c("A", glucose_col, temperature_col, activity_col))
  
  glucose_col <- cellranger::letter_to_num(glucose_col)
  temperature_col <- cellranger::letter_to_num(temperature_col) 
  activity_col <- cellranger::letter_to_num(activity_col)
  
  col_types <- rep("skip", max(glucose_col, temperature_col, activity_col, na.rm = TRUE))
  col_types[1] <- "date"
  col_types[4] <- "text"
  col_types[5] <- "text"
  col_types[glucose_col] <- "numeric"
  col_types[temperature_col] <- "numeric"
  col_types[activity_col] <- "numeric"
  
  
  raw_data <- readxl::read_xlsx(file, 
                                sheet = sheet_name, 
                                range = col_range,
                                col_types = col_types
  )
  data.table::setDT(raw_data)
  
  # We need to standardize the names for glucose, activity and temperature
  idx <- order(c(glucose_col, temperature_col, activity_col), na.last = NA)
  new_names <- c("Date", "ElapsedTime", "Event", c("Glucose", "Temperature", "Activity")[idx])
  
  data.table::setnames(raw_data, new_names, skip_absent = TRUE)
  
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
#'
#' @return cgm_experiment object
#' @importFrom foreach %dopar%
#' @export
#'
#' @examples 
#' \dontrun{
#' to-do
#' }
prepare_experiment <- function(data_file, configuration_file, pattern = "Parameters") {
  if (!file.exists(data_file)) stop(paste(data_file, "not found"))
  samples <- stringr::str_subset(readxl::excel_sheets(data_file), pattern)
  if (!file.exists(configuration_file)) {
    
    events <- readxl::read_xlsx(data_file, sheet = "Events")
    
    create_config(path = configuration_file, sample_names = samples, events = events)
    warning("A settings file was not found, so one was generated.\nPlease fill it out with subject inclusion/exclusion, grouping information, light cycle and excluded timepoints and re-run script.", call. = FALSE)
    return(NULL)
  }
  
  configuration <- read_config(configuration_file, samples)
  
  samples_to_load <- configuration$groupings$SampleID[tolower(configuration$groupings$`Include (Y/N)`) == "y"]
  
  glc_data <- foreach::foreach(i = samples_to_load, .packages = "continousGlucoseMonitoring") %dopar% 
    data_reader(sheet_name = i, file = data_file, 
                glucose_col = configuration$settings$glucose_col, 
                temperature_col = configuration$settings$temperature_col,
                activity_col = configuration$settings$activity_col)
  
  names(glc_data) <- samples_to_load
  glucose_data <- glucose_data(glc_data)
  
  configuration$exclusions <- configuration$exclusions[c("all", samples_to_load)]
  configuration$groupings <- configuration$groupings[match(samples_to_load, configuration$groupings$SampleID), ]
  
  cgm_experiment(glucose_data, configuration)
}
