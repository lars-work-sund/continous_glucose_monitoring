#' Main wrapper for analyzing a continuous glucose monitoring experiment
#'
#' @param data_file path to data file
#' @param configuration_file path to configuration file, if missing one will be created
#' @param out_folder path to folder where results will be saved, if missing one will be created
#' @param pattern Regex used to recognize data sheets
#' @param parallel logical, should data be loaded in parallel?
#' @param reload logical, should previously preprocessed data be reloaded
#'
#' @return cge object
#' @importFrom foreach %dopar%
#' @importFrom ggplot2 %+%
#' @import data.table
#' @export
#'
#' @examples
#' \dontrun{
#' to-do
#' }
analyse_experiment <- function(data_file, configuration_file, out_folder, pattern = "Parameters", parallel = TRUE, reload = TRUE) {
  # Silence no visible binding warnings
  Activity <- Temperature <- Glucose <- baseline <- peak <- nadir <- NULL
  
  ###############################
  ##### Loading and preprocessing
  ###############################
  time_start <- proc.time()
  # If preprocessing has already been performed, then consider reloading
  if (parallel) {
    cores <- parallel::detectCores()
    cl <- parallel::makeCluster(cores[1] - 1)
    doParallel::registerDoParallel(cl)
  
    on.exit({
    parallel::stopCluster(cl)
    foreach::registerDoSEQ()
    })
  }
  
  if(file.exists(file.path(out_folder, "preprocessed_data.RDS")) & reload) {
    stored_data <- file.path(out_folder, "preprocessed_data.RDS")
    cge <- readRDS(stored_data)
    cge$data <- lapply(cge$data, data.table::setDT)
    message(paste("Preprocessed data loaded from:", stored_data))
    return(cge)
  } else {
    message("Loading glucose data")
    cge <- prepare_experiment(data_file = data_file, 
                              configuration_file = configuration_file, 
                              pattern = pattern
                              )
    if (is.null(cge)) return(NULL) # Config file did not exist
    cge$data <- lapply(cge$data, data.table::setDT) #I should figure out why I need this
    
    # Update names with an alias
    idx_alias <- !is.na(cge$config$groupings$Alias)
    if (any(idx_alias)){
      old_names <- names(cge$data)
      cge$config$groupings$SampleID[idx_alias] <- cge$config$groupings$Alias[idx_alias]
      cge <- update_names(cge, cge$config$groupings$SampleID)
      cge$config$groupings$Alias[idx_alias] <- old_names[idx_alias]
      
      Map(data.table::set, cge$data, j = "Sample_ID", value = names(cge$data))
    }
    
    # Run preprocessing
    message("Running pre-processing")
    glc_preprocessed <- foreach::foreach(i = names(cge$data), .packages = "continousGlucoseMonitoring") %dopar% {
      run_standard_preprocess_pipeline(sample_id = i, cge = cge)
    }
    
    names(glc_preprocessed) <- names(cge$data)
    glc_preprocessed <- lapply(glc_preprocessed, data.table::setDT) #I should figure out why I need this
    
    # Add group information
    Map(data.table::set, glc_preprocessed, j = "Group", value = cge$config$groupings$Group[match(names(cge$data), cge$config$groupings$SampleID)])
    
    # To avoid cge$data getting changes when we modify Date in glc_preprocessed
    cge$data <- lapply(glc_preprocessed, copy)
    
    dir.create(file.path(out_folder), showWarnings = FALSE)
    message("Writing pre-processed files to disk")
    
    # Dates get converted to UTC when saving to excel. We format them before that and save as character
    for (i in glc_preprocessed) {
      set(i, j = "Date", value = format(i$Date, format = "%F %H:%M"))
      #set(i, j = "Date", value = format(i$Date, format = "%FT%H:%M:%S%z")) # ISO
    }
    
    writexl::write_xlsx(glc_preprocessed, file.path(out_folder, "preprocessed_samples.xlsx"))
    
  }
  
  
  ###########################
  ##### Kinetics calculations
  ###########################
  message("Performing kinetics calculations")
  
  # Faster not to use parallel processing
  kinetics <- lapply(cge$data, get_sample_kinetics, 
                        excursion_high = get_option(cge, "excursion_high"), 
                        min_peak_duration = get_option(cge, "min_peak_duration"),
                        datapoints_for_slope = get_option(cge, "datapoints_for_slope"))
  
  kinetics_all <- do.call(what = "rbind", kinetics)
  cge$kinetics <- kinetics_all
  
  data.table::set(kinetics_all, j = "nestedPeakType", value = factor(kinetics_all$nestedPeakType, levels = c("Single", "First", "Internal", "Last")))
  uptakes <- kinetics_all[, c(as.list(stats::coef(stats::lm(maxUptake ~ excursion, data = .SD))),
                                rsquared = summary(stats::lm(maxUptake ~ excursion, data = .SD))$r.squared), 
                            .SDcols = c("maxUptake", "excursion"), keyby = c("nestedPeakType", "Sample_ID")] 
  
  clearance <- kinetics_all[, c(as.list(stats::coef(stats::lm(-maxClearance ~ excursion, data = .SD))),
                                   rsquared = summary(stats::lm(-maxClearance ~ excursion, data = .SD))$r.squared), 
                               .SDcols = c("maxClearance", "excursion"), keyby = c("nestedPeakType","Sample_ID")]
  
  Sample_ID <- nestedPeakType <- NULL
  pbase <- ggplot2::ggplot(kinetics_all, ggplot2::aes_string(x = "excursion")) +
    ggplot2::geom_point(ggplot2::aes_string(colour = "nestedPeakType"), size = 1) +
    ggplot2::geom_smooth(data = kinetics_all[nestedPeakType == "Single"], method='lm',formula=y~x) +
    ggplot2::facet_wrap(~Sample_ID)
  
  p1 <- pbase %+% ggplot2::aes_string(y = "maxUptake")
  p2 <- pbase %+% ggplot2::aes_string(y = "maxClearance")
  
  message("Writing kinetics results to disk")
  writexl::write_xlsx(kinetics, file.path(out_folder, "kinetics.xlsx"))
  writexl::write_xlsx(list(uptakes = uptakes, clearance = clearance), file.path(out_folder, "spring_constants.xlsx"))
  
  #
  suppressWarnings({
    ggplot2::ggsave(p1, filename = file.path(out_folder, "excursion_vs_max_uptake.pdf"), width = 15, height = 15, units = "cm", dpi = 320)
    ggplot2::ggsave(p2, filename = file.path(out_folder, "excursion_vs_max_clearance.pdf"), width = 15, height = 15, units = "cm", dpi = 320)
  })
  
  ###############
  ##### Summaries
  ###############
  message("Summarizing data")
  all_stats <- list(
    glucose = glucose_statistics(cge),
    peak = peak_statistics(cge),
    isoglycemic = isoglycemic_statistics(cge),
    kinetics = kinetics_statistics(cge)
  )
  if (length(get_option(cge, "activity_col")) > 0) {
    all_stats[["activity"]] <- other_statistics(cge, Activity)
  }
  if (length(get_option(cge, "temperature_col")) > 0) {
    all_stats[["temperature"]] <- other_statistics(cge, Temperature)
  }
  
  cge$stats <- all_stats
  
  message("Writing summaries")
  writexl::write_xlsx(all_stats$glucose, file.path(out_folder, "glucose_statistics.xlsx"))
  writexl::write_xlsx(all_stats$peak, file.path(out_folder, "peak_statistics.xlsx"))
  writexl::write_xlsx(all_stats$isoglycemic, file.path(out_folder, "isoglycemic_statistics.xlsx"))
  writexl::write_xlsx(all_stats$kinetics, file.path(out_folder, "kinetics_statistics.xlsx"))
  
  if (length(get_option(cge, "activity_col")) > 0) {
    writexl::write_xlsx(all_stats$activity, file.path(out_folder, "activity_statistics.xlsx"))
  }
  if (length(get_option(cge, "temperature_col"))> 0) {
    writexl::write_xlsx(all_stats$temperature, file.path(out_folder, "temperature_statistics.xlsx"))
  }

  message("Preparing profiles")
  glucose_profile <- get_profiles(cge, stat = Glucose, 
                                  low = get_option(cge, "profile_glucose_bins")[2], 
                                  high = get_option(cge, "profile_glucose_bins")[3],
                                  step = get_option(cge, "profile_glucose_bins")[1],
                                  as_percent = TRUE
                                  )
  isoglycemic_profile <- get_profiles(cge, stat = Glucose - baseline, 
                                      low = get_option(cge, "profile_peak_iso_bins")[2], 
                                      high = get_option(cge, "profile_peak_iso_bins")[3],
                                      step = get_option(cge, "profile_peak_iso_bins")[1],
                                      as_percent = TRUE
                                      )
  peak_frequency_profile <- get_profiles(cge, stat = (Glucose - baseline), subset_expression = (peak | nadir), 
                                         low = get_option(cge, "profile_peak_iso_bins")[2], 
                                         high = get_option(cge, "profile_peak_iso_bins")[3],
                                         step = get_option(cge, "profile_peak_iso_bins")[1],
                                         as_percent = FALSE
                                         )

  cge$profiles <- list(
    glucose = glucose_profile,
    peak = peak_frequency_profile,
    isoglycemic = isoglycemic_profile
  )
  
  message("Writing profiles")
  writexl::write_xlsx(glucose_profile, file.path(out_folder, "Absolute BG Profile.xlsx")) # Previoulsy Time in Absolute BG Ranges.xlsx
  writexl::write_xlsx(isoglycemic_profile, file.path(out_folder, "Isoglycemic BG Profile.xlsx")) # Previously Isoglycemic Profile.xlsx
  writexl::write_xlsx(peak_frequency_profile, file.path(out_folder, "Peak Frequency.xlsx")) # Previously Excursion Frequence.xlsx
  
  
  
  message("Saving cge object")
  saveRDS(cge, file = file.path(out_folder, "preprocessed_data.RDS"))
  
  run_time <- (proc.time() - time_start)["elapsed"]
  message(paste("Dataset analysed in", round(lubridate::duration(run_time))))
  
  cge
}
