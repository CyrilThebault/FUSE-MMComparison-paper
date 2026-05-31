# This R script was adapted from Knoben (2024) 
#
# Knoben, W.: CH-Earth/multi-model-mosaic-paper: Peer review release, Zenodo [code], 
# https://doi.org/10.5281/zenodo.13515769, 2024. a, b


#! ---------------------------------------------------------------------------------------
#!
#! Description       :
#!
#! Authors           : Cyril Thebault <cyril.thebault@ucalgary.ca>
#!
#! Creation date     : 2024-09-14 13:37:55
#! Modification date : 2025-09-16 09:17:23
#!
#! Comments          : Script modified to work within Cyril's framework
#!
#! ---------------------------------------------------------------------------------------


# Load required libraries
library(dplyr)
library(tidyr)
library(lpSolve)

##########################
#   FUNCTIONS
##########################

# Helper functions

loadRData <- function(file_name) {
  load(file_name)
  get(ls()[ls() != "file_name"])
}

# Modified bootjack function from gumboot package (Clark et al., 2021) to use KGEcomp score
#
# Clark, M. P., Vogel, R. M., Lamontagne, J. R., Mizukami, N., Knoben, W. J. M., Tang, G., 
# Gharari, S., Freer, J. E., Whitfield, P. H., Shook, K. R., & Papalexiou, S. M. (2021). 
# The Abuse of Popular Performance Metrics in Hydrologic Modeling. Water Resources Research, 57(9), 
# e2020WR029001. https://doi.org/10.1029/2020WR029001

bootjack_v2 <- function (flows, GOF_stat = c("NSE", "KGE", "KGEcomp"), nSample = 1000, 
                         waterYearMonth = 10, startYear = NULL, endYear = NULL, minDays = 100, 
                         minYears = 10, returnSamples = FALSE, seed = NULL, bootYearFile = NULL) {
  if ("KGE" %in% GOF_stat) {
    KGE_is_present <- TRUE
  } else {
    KGE_is_present <- FALSE
  }
  
  if ("NSE" %in% GOF_stat) {
    NSE_is_present <- TRUE
  } else {
    NSE_is_present <- FALSE
  }
  
  if ("KGEcomp" %in% GOF_stat) {
    KGEcomp_is_present <- TRUE
  } else {
    KGEcomp_is_present <- FALSE
  }
  
  options(dplyr.summarise.inform = F)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  write_bootyears <- FALSE
  read_bootyears <- FALSE
  if (!is.null(bootYearFile)) {
    if (!file.exists(bootYearFile)){ 
      write_bootyears <- TRUE
    } else {
      read_bootyears <- TRUE
    }
  }
  flows$year <- as.numeric(format(flows$date, format = "%Y"))
  flows$month <- as.numeric(format(flows$date, format = "%m"))
  flows$day <- as.numeric(format(flows$date, format = "%d"))
  iyUnique <- unique(flows$year)
  nTrials <- length(endYear)
  flows$iyWater <- ifelse(flows$month >= waterYearMonth, flows$year + 
                            1, flows$year)
  iyWater <- ifelse(flows$month >= waterYearMonth, flows$year + 
                      1, flows$year)
  nYears <- length(unique(iyWater))
  
  
  yearsJack <- array(-9999, c(nYears))
  yearsBoot <- array(-9999, c(nYears, nSample))
  zeroVal <- -1e-10
  ixValid <- which((flows$obs > zeroVal) & (flows$sim > zeroVal))
  good_flows <- flows[ixValid, ]
  valid_days <- good_flows %>% group_by(iyWater) %>% summarise(good_days = n_distinct(date))
  valid_years <- valid_days[valid_days$good_days > minDays, ]
  if (!is.null(startYear)) {
    valid_years <- valid_years[valid_years$iyWater >= startYear, ]
  }
  if (!is.null(endYear)) {
    valid_years <- valid_years[valid_years$iyWater <= endYear, ]
  }
  nyValid <- nrow(valid_years)
  if (nyValid < minYears) {
    errorStats <- data.frame(GOF_stat = "", seJack = NA_real_, 
                             seBoot = NA_real_, p05 = NA_real_, p50 = NA_real_, 
                             p95 = NA_real_, score = NA_real_, biasJack = NA_real_, 
                             biasBoot = NA_real_, seJab = NA_real_)
    return(errorStats)
  }
  izUnique <- valid_years$iyWater
  samplingStrategy <- c("jack", "boot")
  n_sampling_strategies <- length(samplingStrategy)
  mSample <- c(nyValid + 1, nSample + 1)
  
  statsJack <- lapply(1:3, function(x) data.frame(matrix(NA, nrow = mSample[1], ncol = 6)))
  names(statsJack) <- GOF_stat
  colnames_list <- c("meanSim", "meanObs", "varSim", "varObs", "rProd", "score")
  statsJack <- lapply(statsJack, function(df) { colnames(df) <- colnames_list; df })
  
  statsBoot <- lapply(1:3, function(x) data.frame(matrix(NA, nrow = mSample[2], ncol = 6)))
  names(statsBoot) <- GOF_stat
  colnames_list <- c("meanSim", "meanObs", "varSim", "varObs", "rProd", "score")
  statsBoot <- lapply(statsBoot, function(df) { colnames(df) <- colnames_list; df })
  
  if("KGEcomp" %in% GOF_stat){
    
    statsJack$KGEcomp$kge <- NA
    statsJack$KGEcomp$meanSiminv <- NA
    statsJack$KGEcomp$meanObsinv <- NA
    statsJack$KGEcomp$varSiminv <- NA
    statsJack$KGEcomp$varObsinv <- NA
    statsJack$KGEcomp$rProdinv <- NA
    statsJack$KGEcomp$kgeinv <- NA
    
    statsJack$KGEcomp <- statsJack$KGEcomp[, c(
      "meanSim", "meanObs", "varSim", "varObs", "rProd",
      "kge", "meanSiminv", "meanObsinv", "varSiminv", "varObsinv", "rProdinv", "kgeinv",
      "score"
    )]
    
    
    statsBoot$KGEcomp$kge <- NA
    statsBoot$KGEcomp$meanSiminv <- NA
    statsBoot$KGEcomp$meanObsinv <- NA
    statsBoot$KGEcomp$varSiminv <- NA
    statsBoot$KGEcomp$varObsinv <- NA
    statsBoot$KGEcomp$rProdinv <- NA
    statsBoot$KGEcomp$kgeinv <- NA
    
    statsBoot$KGEcomp <- statsBoot$KGEcomp[, c(
      "meanSim", "meanObs", "varSim", "varObs", "rProd",
      "kge", "meanSiminv", "meanObsinv", "varSiminv", "varObsinv", "rProdinv", "kgeinv",
      "score"
    )]
  }
  
  
  
  if (read_bootyears) {
    bootYears <- read.csv(bootYearFile, header = FALSE)
  }
  for (iStrategy in 1:n_sampling_strategies) {
    for (iSample in 1:mSample[iStrategy]) {
      if (iSample == 1) {
        jxValid <- ixValid
      } else {
        if (samplingStrategy[iStrategy] == "jack") {
          yearsJack[iSample - 1] = izUnique[iSample - 
                                              1]
          iyIndex <- which(iyWater[ixValid] != izUnique[iSample - 
                                                          1])
          jxValid <- ixValid[iyIndex]
        }
        else {
          if (!read_bootyears) {
            uRand <- runif(nyValid)
            ixYear <- floor(uRand * nyValid) + 1
            iyYear <- izUnique[ixYear]
          }
          else {
            iyYear <- bootYears[, iSample - 1]
          }
          yearsBoot[1:nyValid, iSample - 1] <- iyYear
          iyIndex <- which(iyWater[ixValid] == iyYear[1])
          for (iYear in 2:nyValid) iyIndex <- c(iyIndex, 
                                                which(iyWater[ixValid] == iyYear[iYear]))
          jxValid <- ixValid[iyIndex]
        }
      }
      
      
      if (NSE_is_present) {
        
        qSimValid <- flows$sim[jxValid]
        qObsValid <- flows$obs[jxValid]
        wyValid = iyWater[jxValid]
        meanSim <- mean(qSimValid, na.rm = TRUE)
        meanObs <- mean(qObsValid, na.rm = TRUE)
        varSim <- var(qSimValid, na.rm = TRUE)
        varObs <- var(qObsValid, na.rm = TRUE)
        rProd <- cor(qSimValid, qObsValid)
        
        yBeta = (meanObs - meanSim)/sqrt(varObs)
        alpha = sqrt(varSim)/sqrt(varObs)
        nse = 2 * alpha * rProd - yBeta^2 - alpha^2
        if (samplingStrategy[iStrategy] == "jack") {
          statsJack$NSE[iSample,] = c(meanSim, meanObs,
                                      varSim, varObs,
                                      rProd,
                                      nse)
        }
        if (samplingStrategy[iStrategy] == "boot") {
          statsBoot$NSE[iSample,] = c(meanSim, meanObs,
                                      varSim, varObs,
                                      rProd,
                                      nse)
        }
      }
      
      if (KGE_is_present) {
        
        qSimValid <- flows$sim[jxValid]
        qObsValid <- flows$obs[jxValid]
        wyValid = iyWater[jxValid]
        meanSim <- mean(qSimValid, na.rm = TRUE)
        meanObs <- mean(qObsValid, na.rm = TRUE)
        varSim <- var(qSimValid, na.rm = TRUE)
        varObs <- var(qObsValid, na.rm = TRUE)
        rProd <- cor(qSimValid, qObsValid)
        
        xBeta <- meanSim/meanObs
        alpha <- sqrt(varSim)/sqrt(varObs)
        kge <- 1 - sqrt((xBeta - 1)^2 + (alpha - 1)^2 + 
                          (rProd - 1)^2)
        
        if (samplingStrategy[iStrategy] == "jack") {
          statsJack$KGE[iSample,] = c(meanSim, meanObs,
                                      varSim, varObs,
                                      rProd,
                                      kge)
        }
        if (samplingStrategy[iStrategy] == "boot") {
          statsBoot$KGE[iSample,] = c(meanSim, meanObs,
                                      varSim, varObs,
                                      rProd,
                                      kge)
        }
      }
      
      if (KGEcomp_is_present) {
        
        epsilon <- mean(flows$obs, na.rm = TRUE )/100
        transfo <- -1
        qSimValidinv <- (epsilon+flows$sim[jxValid])^transfo
        qObsValidinv <- (epsilon+flows$obs[jxValid])^transfo
        wyValid = iyWater[jxValid]
        meanSiminv <- mean(qSimValidinv, na.rm = TRUE)
        meanObsinv <- mean(qObsValidinv, na.rm = TRUE)
        varSiminv <- var(qSimValidinv, na.rm = TRUE)
        varObsinv <- var(qObsValidinv, na.rm = TRUE)
        rProdinv <- cor(qSimValidinv, qObsValidinv)
        
        xBetainv <- meanSiminv/meanObsinv
        alphainv <- sqrt(varSiminv)/sqrt(varObsinv)
        kgeinv <- 1 - sqrt((xBetainv - 1)^2 + (alphainv - 1)^2 + 
                             (rProdinv - 1)^2)
        
        
        qSimValid <- flows$sim[jxValid]
        qObsValid <- flows$obs[jxValid]
        wyValid = iyWater[jxValid]
        meanSim <- mean(qSimValid, na.rm = TRUE)
        meanObs <- mean(qObsValid, na.rm = TRUE)
        varSim <- var(qSimValid, na.rm = TRUE)
        varObs <- var(qObsValid, na.rm = TRUE)
        rProd <- cor(qSimValid, qObsValid)
        
        xBeta <- meanSim/meanObs
        alpha <- sqrt(varSim)/sqrt(varObs)
        kge <- 1 - sqrt((xBeta - 1)^2 + (alpha - 1)^2 + 
                          (rProd - 1)^2)
        
        
        kgecomp = (kge+kgeinv)/2
        
        
        
        if (samplingStrategy[iStrategy] == "jack") {
          statsJack$KGEcomp[iSample,] = c(meanSim, meanObs,
                                          varSim, varObs,
                                          rProd,
                                          kge,
                                          meanSiminv, meanObsinv,
                                          varSiminv, varObsinv,
                                          rProdinv,
                                          kgeinv,
                                          kgecomp)
        }
        if (samplingStrategy[iStrategy] == "boot") {
          statsBoot$KGEcomp[iSample,] = c(meanSim, meanObs,
                                          varSim, varObs,
                                          rProd,
                                          kge,
                                          meanSiminv, meanObsinv,
                                          varSiminv, varObsinv,
                                          rProdinv,
                                          kgeinv,
                                          kgecomp)
        }
        
      }
      
      
    }
  }
  if (write_bootyears & (samplingStrategy[iStrategy] == "boot")) {
    write.table(yearsBoot, file = bootYearFile, row.names = FALSE, 
                col.names = FALSE, sep = ",")
  }
  if (returnSamples) {
    return_vals <- list(statsBoot = statsBoot, statsJack = statsJack)
    return(return_vals)
  }
  errorStats <- data.frame(GOF_stat = GOF_stat, seJack = NA_real_, 
                           seBoot = NA_real_, p05 = NA_real_, p50 = NA_real_, p95 = NA_real_, 
                           score = NA_real_, biasJack = NA_real_, biasBoot = NA_real_, 
                           seJab = NA_real_)
  # if (any(sapply(statsJack, function(df) any(df <= -9998, na.rm = TRUE)))) {
  #   return(errorStats)
  # }
  numstats <- length(GOF_stat)
  colnames <- names(statsJack)
  if (KGE_is_present) {
    kge_col <- "KGE"
  }
  if (NSE_is_present) {
    nse_col <- "NSE"
  }
  if (KGEcomp_is_present) {
    kgecomp_col <- "KGEcomp"
  }
  
  for (iPlot in 1:numstats) {
    if (GOF_stat[iPlot] == "NSE") {
      ixPos <- nse_col
    }
    if (GOF_stat[iPlot] == "KGE") {
      ixPos <- kge_col
    }
    if (GOF_stat[iPlot] == "KGEcomp") {
      ixPos <- kgecomp_col
    }
    xJack <- statsJack[[ixPos]][, "score"]
    xBoot <- statsBoot[[ixPos]][, "score"]
    iSort <- order(xJack)
    score <- xJack[1]
    zJack <- xJack[2:(nYears + 1)]
    zBoot <- xBoot[2:(nSample + 1)]
    ixJack <- which(zJack > -9998 & (!is.na(zJack)))
    nJack <- length(ixJack)
    jackMean <- mean(zJack, na.rm = TRUE)
    jackScore <- (nJack * score) - (nJack - 1) * jackMean
    sumSqErr <- (nJack - 1) * sum((jackMean - zJack[ixJack])^2)
    seJack <- sqrt(sumSqErr/nJack)
    ySample <- zBoot[order(zBoot)]
    seBoot <- sd(zBoot)
    p05 <- ySample[floor(0.05 * nSample) + 1]
    p50 <- ySample[floor(0.5 * nSample) + 1]
    p95 <- ySample[floor(0.95 * nSample) + 1]
    biasJack <- (nJack - 1) * (jackMean - score)
    biasBoot <- mean(zBoot) - score
    jabData <- vector("numeric", nYears)
    for (iYear in 2:(nYears + 1)) {
      matchYear <- vector("integer", nSample)
      for (iSample in 1:nSample) {
        ixMatch <- which(yearsBoot[, iSample] == iyUnique[iYear])
        nMatch <- length(ixMatch)
        matchYear[iSample] <- nMatch
      }
      ixMissing <- which(matchYear == 0)
      nMissing <- length(ixMissing)
      xSample <- zBoot[ixMissing]
      ySample <- xSample[order(xSample)]
      p05jack_R <- quantile(xSample, 0.05, type = 3)
      p95jack_R <- quantile(xSample, 0.95, type = 3)
      p05jack <- ySample[floor(0.05 * nMissing) + 1]
      p95jack <- ySample[floor(0.95 * nMissing) + 1]
      jabData[iYear - 1] <- p95jack - p05jack
    }
    jabMean <- mean(jabData)
    sumSqErr <- (nYears - 1) * sum((jabMean - jabData)^2)
    seJab <- sqrt(sumSqErr/nYears)
    errorStats[iPlot, ] <- c(GOF_stat[iPlot], seJack, seBoot, 
                             p05, p50, p95, score, biasJack, biasBoot, seJab)
  }
  errorStats$seJack <- as.numeric(errorStats$seJack)
  errorStats$seBoot <- as.numeric(errorStats$seBoot)
  errorStats$p05 <- as.numeric(errorStats$p05)
  errorStats$p50 <- as.numeric(errorStats$p50)
  errorStats$p95 <- as.numeric(errorStats$p95)
  errorStats$score <- as.numeric(errorStats$score)
  errorStats$biasJack <- as.numeric(errorStats$biasJack)
  errorStats$biasBoot <- as.numeric(errorStats$biasBoot)
  errorStats$seJab <- as.numeric(errorStats$seJab)
  return(errorStats)
}

# function to minimize the number of models used across the catchments based on sampling uncertainty
min_set_cover <- function(constraint_matrix) {
  
  obj <- rep(1, ncol(constraint_matrix))
  
  constraint_matrix[is.na(constraint_matrix)] <- FALSE
  
  constraint_matrix <- as.data.frame(lapply(constraint_matrix, as.numeric))
  
  # Check if any elements are not covered by any subset
  uncovered_elements <- as.character(df_eval_gumboot$gauge_id[rowSums(constraint_matrix) == 0])
  
  # Remove rows corresponding to uncovered elements
  if (length(uncovered_elements) > 0) {
    constraint_matrix <- constraint_matrix[rowSums(constraint_matrix) > 0, , drop = FALSE]
  }
  
  # Set up the linear programming problem
  direction <- rep(">=", nrow(constraint_matrix))
  rhs <- rep(1, nrow(constraint_matrix))
  
  # Solve the problem
  result <- lp("min", obj, constraint_matrix, direction, rhs, all.bin = TRUE)
  
  # Check if a solution was found
  if (result$status != 0) {
    stop("No feasible solution found")
  }
  
  covered_elements = as.character(df_eval_gumboot$gauge_id[!as.character(df_eval_gumboot$gauge_id) %in% uncovered_elements])
  # Return the results
  list(
    num_selected = sum(result$solution),
    selected_indices = which(result$solution == 1),
    covered_elements = covered_elements,
    uncovered_elements = uncovered_elements,
    coverage_percentage = round(length(covered_elements) / nrow(df_eval_gumboot) * 100, 2)
  )
}

##########################
#   VARIABLES
##########################

FUSE_path <-  "/work/comphyd_lab/users/cyril.thebault/Postdoc_Ucal/02_DATA/FUSE"
CAMELS_path <- "/home/cyril.thebault/Postdoc_Ucal/02_DATA/CAMELS"

dataset = "CAMELS"
spatialisation = "Lumped"
inputdata = "daymet"
CritCal = "KGECOMP"
TransfoCal="1" #needed for the tree even if not used with KGECOMP; should match the value set during FUSE execution

SimStart = "1989-01-01"
SimEnd = "2009-12-31"

CalStart = "1989-01-01"
CalEnd = "1998-12-31"

WU = 2

##########################
#   MAIN
##########################

print("---- Get metrics for FUSE")

# Load gauge metadata
df <- read.csv(file.path(CAMELS_path, "camels_topo.txt"), sep = ";") %>%
  select(gauge_id, gauge_lat, gauge_lon, area_gages2, area_geospa_fabric)


Eval_FUSE_long_Comp = loadRData(file.path(FUSE_path, dataset, spatialisation, inputdata, CritCal, TransfoCal, "Eval_FUSE.Rdata"))%>%
  mutate(`Cal : KGEcomp` = (`Cal : KGE :  1` + `Cal : KGE : -1`) / 2,
         `Eval : KGEcomp` = (`Eval : KGE :  1` + `Eval : KGE : -1`) / 2)

Eval_FUSE_Comp <- Eval_FUSE_long_Comp %>%
  mutate(
    gauge_id = as.integer(sub("USA_", "", Codes))
  ) %>%
  select(
    gauge_id,
    ModelDecisions,
    cal_kge = `Cal : KGEcomp`,
    val_kge = `Eval : KGEcomp`
  ) %>%
  pivot_longer(
    cols = c(cal_kge, val_kge),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    var = paste0("Comp_", ModelDecisions, "_", metric)
  ) %>%
  select(gauge_id, var, value) %>%
  pivot_wider(
    names_from = var,
    values_from = value
  ) %>%
  select(
    gauge_id,
    all_of(sort(names(.)[names(.) != "gauge_id"]))
  )

list_names <- names(Eval_FUSE_Comp) %>%
  grep("^Comp_", ., value = TRUE) %>%
  sub("_(cal|val)_kge$", "", .) %>%
  unique()

df_ini <- df %>%
  inner_join(Eval_FUSE_Comp, by = "gauge_id") 

write.csv(df_ini, file.path(FUSE_path, dataset, spatialisation, inputdata, CritCal, TransfoCal, "brief_analysis_updated_modeling_results.csv"), row.names = FALSE)

print("---- Find the best model")
# Create df_eval
df_eval <- df_ini

# Find best model for each basin
df_eval$best_cal_column <- apply(df_eval %>% select(ends_with("cal_kge")), 1, function(x) {
  # Check if x is not empty and contains valid values
  if (length(x) == 0 || all(is.na(x))) {
    return(NA)
  }
  
  # Find the name of the maximum value
  max_index <- which.max(x)
  if (length(max_index) == 0 || is.na(max_index)) {
    return(NA)
  }
  
  max_name <- names(x)[max_index]
  if (is.na(max_name)) {
    return(NA)
  } else {
    return(max_name)
  }
})

df_eval$best_cal_model <- sapply(strsplit(df_eval$best_cal_column, "_"), function(x) {
  if(NA %in% x){
    return(NA)
  } else{
    return(paste(x[1:2], collapse = "_"))
  }
  
})

df_eval$best_cal_score <- apply(df_eval %>% select(ends_with("cal_kge")), 1, function(x){
  if(all(is.na(x))){
    return(NA)
  } else{
    return(max(x, na.rm = TRUE))
  }
})


write.csv(df_eval, file.path(FUSE_path, dataset, spatialisation, inputdata, CritCal, TransfoCal, "brief_analysis_updated_modeling_results_eval.csv"), row.names = FALSE)

print("---- Sampling uncertainty (gumboot)")
# Create df_gumboot
df_gumboot = df_eval

# Prepare for gumboot analysis
gumboot_cols <- c("GOF_stat", "seJack", "seBoot", "p05", "p50", "p95", "score", "biasJack", "biasBoot", "seJab")
for (col in gumboot_cols) {
  df_gumboot[,col] <- NA
}

qobs <- loadRData(file.path(FUSE_path, "ObservedStreamflow.Rdata"))$Qobs_mm
qsim <- loadRData(file.path(FUSE_path, dataset, spatialisation, inputdata, CritCal, TransfoCal, "SimArray.Rdata"))

# Gumboot analysis
for (i in 1:nrow(df_gumboot)) {
  basin <- df_gumboot[i,"gauge_id"]
  print(paste("Running", basin))
  
  qobs_sub <- tibble(
    date = as.Date(rownames(qobs)),
    qobs_mm_d = qobs[[paste0("USA_", sprintf("%08d", basin))]]
  ) %>%
    filter(date >= as.Date(CalStart) & date <= as.Date(CalEnd))
  
  best_cal_model <- sub("Comp_", "",df_gumboot$best_cal_model[i])

  if(is.na(best_cal_model)){
    print(paste("No best model found, skipping basin", basin))
    next
  }
  
  qsim_sub <- tibble(
    date = as.Date(dimnames(qsim)$Dates),
    !!as.character(basin) := qsim[,best_cal_model,paste0("USA_", sprintf("%08d", basin))]
  ) %>%
    filter(date >= as.Date(CalStart) & date <= as.Date(CalEnd))
  
  
  
  foot <- inner_join(qobs_sub, qsim_sub, by = "date") %>%
    select(date, obs = qobs_mm_d, sim = as.character(basin)) %>%
    na.omit()
  
  clean_foot <- na.omit(foot)
  
  # Run gumboot, catching errors if needed
  result <- try({
    result <- bootjack_v2(clean_foot, GOF_stat = "KGEcomp", seed = 1)
  })
  
  # Check if we encountered an error, and return all NaNs if so
  if (inherits(result, "try-error")) {
    print("An error occurred, returning all NA")
    next
  } else if (all(is.na(result))) {
    print("Not enough data, minimum 10 years required")
    next
  }
  
  df_gumboot[i,colnames(result)] = result
}

# Save results with gumboot analysis
write.csv(df_gumboot, file.path(FUSE_path, dataset, spatialisation, inputdata, CritCal, TransfoCal, "brief_analysis_modeling_results_with_gumboot.csv"), row.names = FALSE)

print("---- Analysis of models within uncertainty bounds")

df_eval_gumboot = df_gumboot

# Calculate range
df_eval_gumboot$range_5_95 <- df_eval_gumboot$p95 - df_eval_gumboot$p05

# Analysis of models within uncertainty bounds
columns_cal <- grep("_cal_kge$", names(df_eval_gumboot), value = TRUE)

for (column_cal in columns_cal) {
  model <- paste(strsplit(column_cal, "_")[[1]][1:2], collapse = "_")
  flag_column <- paste0(model, "_above_p5")
  df_eval_gumboot[[flag_column]] <- df_eval_gumboot[[column_cal]] >= df_eval_gumboot$p05
}

columns_flags = grep("_above_p5$", names(df_eval_gumboot), value = TRUE)
df_eval_gumboot$similar_model_count =  rowSums(df_eval_gumboot[,columns_flags], na.rm = TRUE)

constraint_matrix <- df_eval_gumboot[,columns_flags]

# Function to solve the set cover problem
result <- min_set_cover(constraint_matrix)
num_selected <- result$num_selected
selected_indices <- result$selected_indices
selected_models <- list_names[selected_indices]
covered_basins <- colSums(df_eval_gumboot[,columns_flags[selected_indices]], na.rm = TRUE)

# Print results
cat("Number of subsets selected:", num_selected, "\n")
cat("Indices of selected subsets:", selected_indices, "\n")
cat("Selected models:", selected_models, "\n")
cat("Basins covered:", covered_basins, "\n")


df_eval_gumboot$mod_need = NA

tmp_df <- df_eval_gumboot[! df_eval_gumboot$similar_model_count == 0, paste0(selected_models, "_above_p5")]
basins_left <- nrow(df)

while (basins_left > 0) {
  # Perform a greedy selection for the basins we have left
  model_needed <- names(which.max(colSums(tmp_df, na.rm = TRUE)))
  
  # Subset the dataframe to keep only those basins for which we have no model yet
  
  ind_tmp <- tmp_df[[model_needed]]
  ind_tmp[is.na(ind_tmp)] <- FALSE
  tmp_df <- tmp_df[!ind_tmp, ]
  
  ind_gumboot = is.na(df_eval_gumboot$mod_need) & df_eval_gumboot[[model_needed]]
  ind_gumboot[is.na(ind_gumboot)] <- FALSE
  df_eval_gumboot$mod_need[ind_gumboot] <- sapply(strsplit(model_needed, "_"), function(x) paste(x[1:2], collapse = "_"))
  
  # Update the count of basins we still need - this will terminate the loop
  basins_left <- nrow(tmp_df)
}

print(sort(table(df_eval_gumboot$mod_need), decreasing = TRUE))

# Save results with gumboot analysis
write.csv(df_eval_gumboot, file.path(FUSE_path, dataset, spatialisation, inputdata, CritCal, TransfoCal, "brief_analysis_modeling_results_with_gumboot_eval.csv"), row.names = FALSE)

