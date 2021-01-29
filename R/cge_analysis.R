#' Main wrapper for analyzing a continuous glucose monitoring experiment
#'
#' @param data_file path to data file
#' @param configuration_file path to configuration file, if missing one will be created
#' @param out_folder path to folder where results will be saved, if missing one will be created
#' @param pattern Regex used to recognize data sheets
#' @param parallel logical, should data be loaded in parallel?
#' @param reload logical, should previously preprocessed data be reloaded
#'
#' @return
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
    message(paste("Preprocessed data loaded from:", stored_data))
  } else {
    message("Loading glucose data")
    cge <- prepare_experiment(data_file = data_file, 
                              configuration_file = configuration_file, 
                              pattern = pattern
                              )
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
    cge <- preprocess_samples(cge)
    cge$data <- lapply(cge$data, data.table::setDT) #I should figure out why I need this
    
    # Add group information
    Map(data.table::set, cge$data, j = "Group", value = cge$config$groupings$Group[match(names(cge$data), cge$config$groupings$SampleID)])
    
    dir.create(file.path(out_folder), showWarnings = FALSE)
    message("Writing pre-processed files to disk")
    writexl::write_xlsx(cge$data, file.path(out_folder, "preprocessed_samples.xlsx"))
    saveRDS(cge, file = file.path(out_folder, "preprocessed_data.RDS"))
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
  message("Summarizing glucose data")
  glucose_stats <- glucose_statistics(cge)
  peak_stats <- peak_statistics(cge)
  isoglycemic_stats <- isoglycemic_statistics(cge)

  message("Writing summaries")
  writexl::write_xlsx(glucose_stats, file.path(out_folder, "glucose_statistics.xlsx"))
  writexl::write_xlsx(peak_stats, file.path(out_folder, "peak_statistics.xlsx"))
  writexl::write_xlsx(isoglycemic_stats, file.path(out_folder, "isoglycemic_statistics.xlsx"))

  message("Preparing profiles")
  glucose_profile <- get_profiles(cge, stat = Glucose, 1, 16)
  isoglycemic_profile <- get_profiles(cge, stat = Glucose - baseline, -4.75, 9)
  peak_frequency_profile <- get_profiles(cge, stat = (Glucose - baseline), subset_expression = (peak | nadir), -4.75, 9)

  message("Writing profiles")
  writexl::write_xlsx(glucose_profile, file.path(out_folder, "Time in Absolute BG Ranges.xlsx"))
  writexl::write_xlsx(isoglycemic_profile, file.path(out_folder, "Isoglycemic Profile.xlsx"))
  writexl::write_xlsx(peak_frequency_profile, file.path(out_folder, "Excursion Frequency.xlsx"))
  
  cge$kinetics <- kinetics_all
  cge$stats <- list(
    glucose = glucose_stats,
    peak = peak_stats,
    isoglycemic = isoglucose_stats
  )
  cge$profiles <- list(
    glucose = glucose_profile,
    peak = peak_frequency_profile,
    isoglycemic = isoglycemic_profile
  )
  
  run_time <- (proc.time() - time_start)["elapsed"]
  message(paste("Dataset analysed in", round(lubridate::duration(run_time))))
  
  cge
}
