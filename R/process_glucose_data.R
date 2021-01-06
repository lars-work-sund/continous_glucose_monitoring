convertUnits <- function(x, convertUnits){
  if (convertUnits){
    mgDl_2_mmolL <- 18.018018
    x[, Glucose:=Glucose/mgDl_2_mmolL]
  }
  x
}

# Sometime there will be datapoints taken at odd intervals, discard those
# Sometimes datapoints will be missing, add rows to contain those
addMissingMinutes <- function(x)
{
  toTime <- function(z){
    duration(hour = as.integer(z[, 1]), minute = as.integer(z[, 2]), second = as.integer(z[, 3]))
  }
  x[, experimentTime:=toTime(str_split_fixed(ElapsedTime, ":|,", n = 4)[, -4])]
  
  # Only very few timepoints should be affected by dropping timepoints where a full minute is missing
  before <- nrow(x)
  ind <- x[, which(c(60, diff(experimentTime)) != 60)]
  x <- x[c(60, diff(experimentTime)) == 60, ]
  message(paste("Datapoints excluded due to uneven sampling rate:", before - nrow(x), "\nThese dropouts happen at rows:", paste(ind, collapse = "; ")))
  
  # This convoluted way of finding the end date ensures that daylight savings time is handled
  # I suspect seq should be able to handle seq(expStartTime, expEndTime, by = 60), but that currently doesn't work
  expStartTime <- x[1, as.integer(experimentTime)]
  expEndTime <- x[.N, as.integer(experimentTime)]
  
  y <- x[, .(experimentTime = duration(seq(expStartTime, expEndTime, by = 60)))]
  y <- merge(y, x, by = "experimentTime", all = TRUE)
  
  # Some timepoints are missing data
  startDate <- x[1, Date]
  y[is.na(Date), Date:=startDate + experimentTime]
  
  y[, MouseID:=x$MouseID[[1]]]
  y[, ElapsedTime:=NULL]
  y
}

# light containes the hour light is turned on, and the hour light is turned off
addDerivedDateInfo <- function(x, light, dst = "independent")
{
  # if dst = "independent", light is turned on and off independently of daylight savings time (e.g. switches every 12 hours)
  # if dst = "sync", light follows the time of day (e.g. turns on at 6, no matter how many hours of dark)
  
  # if dst = "independent", light should contain on/off times on the first experimental day
  
  # if the experiment does not cover a daylight savings time shift, both options return the same result. The option sync may be slightly faster
  lightHour <- floor(light)
  lightMinute <- (light - lightHour) * 60
  
  if (dst == "independent"){
    minutesOfLight <- (light[2] - light[1]) * 60
    minutesOfDark <- (24 - (light[2] - light[1])) * 60
    
    firstTimepoint <- x[1, Date - experimentTime]
    firstLightOn <- paste0(year(firstTimepoint), "-",
                           month(firstTimepoint), "-",
                           day(firstTimepoint),  " ",
                           floor(light[1]), ":",
                           (light[1] %% 1) * 60, ":",
                           "0"
    ) %>% ymd_hms(tz = tz(x$Date)) 
    
    # The plus 0.1 second is so that measurements taken at the time the 
    # light goes on is counted as dark.
    firstLightOn <- firstLightOn + seconds(0.1)
    
    firstLightOff <- firstLightOn + minutes(minutesOfLight)
    firstLight <- interval(firstLightOn, firstLightOff)
    
    lastTimePoint <- x[.N, Date]
    
    expDays <- interval(firstTimepoint, lastTimePoint) %>%
      int_length %>%
      divide_by(60 * 60 * 24)
    
    cumulativeHours <- hours(cumsum(c(0, rep(24, ceiling(expDays) + 1))))
    lightOnPeriods <- as.list(int_shift(firstLight, cumulativeHours))
    
    x[, Darkness:=1]
    x[Date %within% lightOnPeriods, Darkness:=0]
    # firstInHours <- hour(firstTimepoint) + minute(firstTimepoint)/60
    # startingPhase <- as.integer(!(firstInHours > lightOn[1] & firstInHours < lightOn[2]))
    # 
    # if (startingPhase == 0){ # Starts in darkness
    #   nDarkLeft <- 
    # }
    # 
    # firstDay <- x[1, c(year(Date), month(Date), day(Date))]
    # 
    # x[year(Date) == firstDay[1] & 
    #     month(Date) == firstDay[2] & 
    #     day(Date) == firstDay[3],
    #   Darkness:=1
    #   ]
    # x[year(Date) == firstDay[1] & 
    #     month(Date) == firstDay[2] & 
    #     day(Date) == firstDay[3] &
    #     (((hour(Date) > lightHour[1]) & (hour(Date) < (lightHour[2]))) | # In the whole hours of light
    #     ((hour(Date) == lightHour[1] & minute(Date) > lightMinute[1]) | # In the minutes after light is on
    #        (hour(Date) == (lightHour[2]) & minute(Date) <= lightMinute[2]))), # In the minutes before light is off
    #   Darkness:=0
    #   ]
    # 
    # startingPhase <- x[1, Darkness]
    # x[Darkness != startingPhase, Darkness:=NA]
    # 
    # minutesOfLight <- (light[2] - light[1]) * 60
    # minutesOfDark <- (24 - (light[2] - light[1])) * 60
    # if (startingPhase == 1){
    #   lightDark <- rep(c(0, 1), c(minutesOfLight, minutesOfDark))
    # } else {
    #   lightDark <- rep(c(1, 0), c(minutesOfDark, minutesOfLight))
    # }
    # lightDarkAll <- rep(lightDark, length.out = x[, sum(is.na(Darkness))])
    # 
    # x[is.na(Darkness), Darkness:=lightDarkAll]
    
  } else if (dst == "sync"){
    x[, Darkness:=1]
    x[((hour(Date) > lightHour[1]) & (hour(Date) < (lightHour[2]))) | # In the whole hours of light
        ((hour(Date) == lightHour[1] & minute(Date) > lightMinute[1]) | # In the minutes after light is on
           (hour(Date) == (lightHour[2]) & minute(Date) <= lightMinute[2])), # In the minutes before light is off
      Darkness:=0]
  } else {
    stop("Daylight savings time setting must be 'independent' or 'sync'")
  }
  zt <- (hour(x$Date) + minute(x$Date)/60 - light[1]) %% 24
  x[, PhaseId:=cumsum(c(1, abs(diff(Darkness))))]
  x[, Day:=c(0, diff(Darkness))]
  x[Day == 1, Day:=0]
  x[, Day:=cumsum(abs(Day))]
  x[, Week:=Day %/% 7]
  x[, ZT:=floor(zt)]
  x[, ZT_exact:=zt]
  x
}

excludeTimepoints <- function(x, settingsFile = NULL)
{
  x[, included:=TRUE]
  if (!is.null(settingsFile)){
    sheets <- readxl::excel_sheets(settingsFile)
    exclusionReader <- function(x){
      read_xlsx(x,
                path = settingsFile, 
                range = cell_cols(c("A", "B")),
                col_types = c("date", "date")
      ) %>%
        data.table
    }
    if(x$MouseID[[1]] %in% sheets){
      relevantSheets <- c("all", x$MouseID[[1]])
    } else {
      relevantSheets <- "all"
    }
    exclusions <- lapply(relevantSheets, exclusionReader) %>%
      do.call(what = "rbind")
    
    if (nrow(exclusions) > 0){
      for (i in seq_len(nrow(exclusions))){
        x[Date %between% exclusions[i, c(`ExclusionStart (DD-MM-YYYY  HH:MM:SS)`, `ExclusionEnd (DD-MM-YYYY  HH:MM:SS)`)], included:=FALSE]
      }
    }
  }
  x
}

addGroups <- function(x, settingsFile){
  groupInfo <- read_xlsx("SampleGrouping", 
                         path = settingsFile
  ) %>% 
    data.table(key = "MouseID")
  groupNames <- colnames(groupInfo)
  groupNames <- groupNames[!(groupNames %in% c("MouseID", "Include (Y/N)"))]
  
  groupInfo <- groupInfo[x$MouseID, ..groupNames]
  x[, (groupNames):=groupInfo]
  x
}

interpolateGlucose <- function(x, maxGap)
{
  # Which runs have a missing value
  glucoseRle <- rle(x[, is.na(Glucose)])
  
  # If the first or the last element(s) are missing we cannot interpolate them,
  # so they are discarded instead
  nRuns <- length(glucoseRle$values)
  badMeasurements <- integer(0)
  
  if (glucoseRle$values[1]){
    badMeasurements <- c(badMeasurements, seq(from = 1, to = glucoseRle$lengths[1]))
  }
  if (glucoseRle$values[nRuns]){
    startOfLastRun <- cumsum(glucoseRle$lengths)[nRuns - 1] + 1
    badMeasurements <- c(badMeasurements, 
                         seq(from = startOfLastRun, to = nrow(x)))
  }
  
  if (length(badMeasurements) > 0)
  {
    x <- x[-badMeasurements]
  }
  
  # Which runs have a missing value and is shorter than the max gap?
  glucoseRle <- rle(x[, is.na(Glucose)])
  toInterpolate <- glucoseRle$values & glucoseRle$lengths <= maxGap
  
  # Which indexes should be interpolated?
  ends <- cumsum(glucoseRle$lengths)
  starts <- c(0, ends[-length(ends)]) + 1
  inds <- Map(seq, starts[toInterpolate], ends[toInterpolate])
  
  # Simple linear interpolation
  for (i in inds){
    beforeVal <- x[min(i) - 1, Glucose]
    afterVal <- x[max(i) + 1, Glucose]
    interpolatedValues <- seq(from = beforeVal, 
                              to = afterVal, 
                              length.out = length(i) + 2)
    interpolatedValues <- interpolatedValues[-c(1, length(interpolatedValues))]
    x[i, Glucose:=interpolatedValues]
  }
  x
}

findBaseline <- function(x, baselineWindow, excursionLow, excursionHigh)
{
  changePoints <- c(0, 0, diff(sign(diff(x$Glucose))))
  x[, baseline:=NA_real_]
  x[(included) & changePoints != 0, baseline:=Glucose]
  
  medianBaseline <- rollapply(x$baseline,
                              width = baselineWindow,
                              FUN = median,
                              na.rm = TRUE,
                              partial = TRUE, # Are we sure we want this?
                              fill = NA,
                              align = "center")
  x[, baseline:=medianBaseline]
  x[!(included), baseline:=NA]
  
  x[, excursion:=Glucose - baseline]
  x[excursion > excursionLow & excursion < excursionHigh, excursion:=NA]
  
  x
}

smoothBackground <- function(x, baselineWindow)
{
  meanBackground <- rollapply(x$baseline,
                              width = baselineWindow,
                              FUN = mean,
                              na.rm = TRUE,
                              partial = TRUE, # Are we sure we want this?
                              fill = NA,
                              align = "center")
  x[, baseline:=meanBackground]
  x[!(included), baseline:=NA]
  x
}

findExcursions <- function(x, excursionLow, excursionHigh)
{
  x[, excursion:=Glucose - baseline]
  x[excursion > excursionLow & excursion < excursionHigh, excursion:=NA]
  x
}

findPeaksAndNadirs <- function(x, maxMinWindow)
{
  # in cases with ties, which.min / which.max returns the first
  
  myMinMax <- function(x, whichFun){
    # if there are only NAs, which.max/min returns integer(0)
    # this confuses rollapply, instead we want a NA
    out <- whichFun(x)
    if (length(out) == 0) return(NA)
    out
  }
  
  rollFun <- function(whichFun)
  {
    zoo::rollapply(x$excursion, 
                   width = maxMinWindow, 
                   FUN = myMinMax,
                   whichFun = whichFun,
                   partial = TRUE,
                   fill = NA,
                   align = "center")
  }
  
  # rollFun returns the position in the defined window, centered on the 
  # observation which has the highest/lowest value
  # when the observation with the highes/lowest value is the centerpoint a 
  # local minimum/maximum is found
  midpoint <- ceiling(maxMinWindow/2)
  
  peaks <- rollFun(which.max) == midpoint
  
  nadirs <- rollFun(which.min) == midpoint
  
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
  
  scanFn <- function(window){
    any(abs(window[-midpoint] - window[midpoint]) < 10^-4)
  }
  duplicatePeak <- zoo::rollapply(x$excursion, 
                                  width = midpoint, 
                                  FUN = scanFn,
                                  fill = NA,
                                  partial = TRUE,
                                  align = "right")
  
  x[peak & duplicatePeak, peak:=FALSE]
  x[nadir & duplicatePeak, nadir:=FALSE]
  x
}

################################################
######## Kinetics
# imputeBloodGlucose <- function(x){
#   ind <- x[, which(is.na(Glucose))]
#   ind <- ind[c(TRUE, diff(ind) != 1)]
#   for (i in ind){
#     x[i, Glucose:=x[c(i-1, i+1), mean(Glucose)]]
#   }
#   x
# }

addPeakTimers <- function(x, minPeakDuration){
  # Uptake and clearence timers
  x[, grp:=rleid(is.na(excursion))]
  x[!is.na(excursion), cumulativeUptakeTime:=1:(.N), by = "grp"]
  x[!is.na(excursion), cumulativeClearenceTime:=(.N):1, by = "grp"]
  
  # Only reset peak counter if excursion lasts longer than minPeakDuration
  x[, tmpPeak:=fifelse(rep(.N >= minPeakDuration, .N), peak, FALSE), by = "grp"]
  
  # # For convenience exclude too short peaks from further analysis (should probably go into seperate function)
  # x[!(tmpPeak), peak:=FALSE]
  
  x[, grp:=cumsum(tmpPeak)]
  x[, timeToNextPeak:=(.N):1, by = "grp"]
  x[, grp:=c(0, grp[-.N])]
  x[, timeSinceLastPeak:=1:(.N), by = "grp"]
  
  x[, grp:=cumsum(!included)]
  x[, timeSinceLastExcluded:=1:(.N), by = "grp"]
  
  x[, grp:=NULL]
  x[, tmpPeak:=NULL]
  x
}

calcSlopes <- function(x, nPoints){
  # Division by 60 rescales so slopes are in change/minute rather than change/second
  #splineFit <- x[!is.na(Glucose) & !is.na(background), 
  # splineFit <- x[!is.na(Glucose) & !is.na(baseline), 
  #                smooth.spline(x = as.integer(experimentTime)/60, 
  #                              y = Glucose - baseline, 
  #                              all.knots = TRUE)]
  # x[, slope:=predict(splineFit, as.integer(experimentTime)/60, deriv = 1)$y]
  # x[, acceleration:=predict(splineFit, as.integer(experimentTime)/60, deriv = 2)$y]
  calcSlope <- function(val, align){
    helper <- function(val, idx) {tryCatch(coef(lm(val ~ idx))[2], error = function(e){NA_real_})}
    frollapply(val, 
               n = nPoints, 
               FUN = helper, 
               idx = seq_len(nPoints),
               align = align)
  }
  slope <- calcSlope(x$Glucose, "right")
  
  x[, uptakeSlope:=slope]
  x[, clearanceSlope:=c(slope[-c(seq_len(nPoints - 1))], rep(NA, nPoints - 1))]
  x
}

selectExcursions <- function(x, minPeakDuration){
  x[, grp:=rleid(is.na(excursion))]
  # We are currently only interested in positive excursions
  x <- x[excursion > 0]
  # Remove peaks that contains excluded values
  # Remove first peak after exclusion (could be due to handling etc.)
  excludedExcursions <- x[, any(!included) | 
                            any(timeSinceLastExcluded <= timeSinceLastPeak) |
                            .N < minPeakDuration |
                            any(is.na(excursion)) |
                            sum(peak) == 0, # small excursions close to a larger excursion may have no peaks 
                          by = "grp"][(V1), grp]
  x <- x[!(grp %in% excludedExcursions)]
  x[, excursionId:=as.integer(rleid(grp))]
  x[, grp:=NULL]
  x
}

tagMultiPeaks <- function(x){
  x[, singlePeak:=rep(sum(peak) == 1, .N), by = "excursionId"]
  x[, nestedPeakType:=NA_character_]
  if (nrow(x[(!singlePeak) & (peak)]) > 0){
    x[(!singlePeak) & (peak), c("nestedPeakType", "nestedPeakNumber"):=list(c("First", rep("Internal", sum(peak) - 2), "Last"),
                                                                            c(NA, seq_len(sum(peak) - 2), NA)), by = "excursionId"]
  } else {
    x[, c("nestedPeakType", "nestedPeakNumber"):=list(NA, NA)]
  }
  x
}
