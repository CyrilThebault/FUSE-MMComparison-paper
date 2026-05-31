# libraries

library(abind)

# functions 

loadRData <- function(file_name) {
  e <- new.env(parent = emptyenv())
  load(file_name, envir = e)
  objs <- ls(envir = e)
  if (length(objs) == 0) return(NULL)
  e[[objs[1]]]
}

# Scalar KGE (kept for window-to-window KGE on vectors)
KGE_vec <- function(sim, obs) {
  ok <- !is.na(obs) & !is.na(sim)
  qsim <- sim[ok]
  qobs <- obs[ok]
  if (length(qobs) < 2) return(list(KGE = NA_real_, alpha = NA_real_, beta = NA_real_, ere = NA_real_))
  
  sigma_Obs <- sd(qobs)
  sigma_Sim <- sd(qsim)
  if (!is.finite(sigma_Obs) || sigma_Obs == 0) return(list(KGE = NA_real_, alpha = NA_real_, beta = NA_real_, ere = NA_real_))
  
  alpha <- sigma_Sim / sigma_Obs
  ere   <- suppressWarnings(cor(qobs, qsim, method = "pearson"))
  beta  <- sum(qsim) / sum(qobs)
  
  list(
    KGE   = 1 - sqrt((ere - 1)^2 + (beta - 1)^2 + (alpha - 1)^2),
    alpha = alpha,
    beta  = beta,
    ere   = ere
  )
}

MAE_vec <- function(sim, obs) {
  ok <- !is.na(obs) & !is.na(sim)
  if (!any(ok)) return(NA_real_)
  mean(abs(obs[ok] - sim[ok]))
}

# Vectorized KGE for all columns in a matrix vs one obs vector (pairwise complete)
KGE_matrix <- function(sim_mat, obs_vec) {
  tau <- length(obs_vec)
  m <- ncol(sim_mat)
  
  obs_mat <- matrix(obs_vec, nrow = tau, ncol = m)
  
  ok <- !is.na(sim_mat) & !is.na(obs_mat)
  n <- colSums(ok)
  bad <- (n < 2L)
  
  sim0 <- sim_mat; sim0[!ok] <- 0
  obs0 <- obs_mat; obs0[!ok] <- 0
  
  sx  <- colSums(sim0)
  sy  <- colSums(obs0)
  sxx <- colSums(sim0 * sim0)
  syy <- colSums(obs0 * obs0)
  sxy <- colSums(sim0 * obs0)
  
  mx <- sx / n
  my <- sy / n
  
  # sample variance/cov with (n-1) denominator, matching sd()/cor()
  vx  <- (sxx - n * mx * mx) / (n - 1)
  vy  <- (syy - n * my * my) / (n - 1)
  cov <- (sxy - n * mx * my) / (n - 1)
  
  sd_x <- sqrt(vx)
  sd_y <- sqrt(vy)
  
  ere   <- cov / (sd_x * sd_y)
  alpha <- sd_x / sd_y
  beta  <- sx / sy
  
  kge <- 1 - sqrt((ere - 1)^2 + (beta - 1)^2 + (alpha - 1)^2)
  
  kge[bad | !is.finite(kge)] <- NA_real_
  kge
}

KGEinv_matrix <- function(sim_mat, obs_vec, epsilon) {
  KGE_matrix((epsilon + sim_mat)^-1, (epsilon + obs_vec)^-1)
}

# control file

FUSE_path = "/project/6079554/thebault/Postdoc_Ucal/02_DATA/FUSE"

dataset = "CAMELS"
spatialisation = "Lumped"
inputdata = "daymet"
CritCal = "KGECOMP"

outputfolder = "Comp_WA_q_KGEcomp"

# Simulations
Qsim = loadRData(file.path(FUSE_path,dataset,spatialisation,inputdata,CritCal,"1", "SimArray.Rdata"))


# Obervations
Qobs = loadRData(file = file.path(FUSE_path,dataset,spatialisation,inputdata,CritCal,"1","ObservedStreamflow.Rdata"))$Qobs

CalStart = "1989-01-01"
CalEnd = "1998-12-31"

EvalStart = "1999-01-01"
EvalEnd = "2009-12-31"

time_windows = seq(1,30,3)
k_models = seq(1,20,2)
k_neighbours = seq(1,20,2)

DatesR = as.Date(dimnames(Qsim)[[1]])
models = dimnames(Qsim)[[2]]
catchments = dimnames(Qsim)[[3]]


mydf = expand.grid(catchments, time_windows, k_models, k_neighbours)
colnames(mydf) = c("Code", "Tau", "Km", "Kn")
mydf = mydf[order(mydf$Code),]

if(file.exists(file.path(FUSE_path,dataset,spatialisation,inputdata,CritCal, outputfolder, "indexfail.Rdata"))){
  
  indexfail = loadRData(file.path(FUSE_path,dataset,spatialisation,inputdata,CritCal, outputfolder, "indexfail.Rdata"))
  
  if(is.null(indexfail)){
    stop("All simulations completed")
  }
  
  mydf = mydf[indexfail,]
  
}

############

args <- commandArgs(trailingOnly = TRUE)
id <- as.numeric(args[1])

catchment = mydf$Code[id]
tau = mydf$Tau[id]
km = mydf$Km[id]
knn = mydf$Kn[id]


Ind_Cal = which(DatesR >= CalStart & DatesR <= CalEnd)
Ind_Eval = which(DatesR >= EvalStart & DatesR <= EvalEnd)

# Subset simulated
simulated = Qsim[,models,catchment]

# Subset observed
qmm = Qobs[,catchment]

#epsilon (used for KGEinv calculation)
epsilon = mean(qmm[Ind_Cal], na.rm = TRUE)/100

#median Cal period
medCal= median(qmm[Ind_Cal], na.rm = TRUE)

# Create a matrix to store the dynamic weights for each time step and each simulation
dynamic_weights <- matrix(NA, nrow = nrow(simulated), ncol = ncol(simulated))
rownames(dynamic_weights) <- dimnames(simulated)[["Dates"]]
colnames(dynamic_weights) <- dimnames(simulated)[["Decisions"]]

# Create a vector to store the weighted average for each time step
dynamic_weighted_avg <- rep(NA,nrow(simulated))
names(dynamic_weighted_avg) <- dimnames(simulated)[["Dates"]]

# Create a matrix to store the metric used for each time step for neighbors and models selection
selection_metric <- matrix(NA, nrow = nrow(simulated), ncol = 2)
rownames(selection_metric) <- dimnames(simulated)[["Dates"]]
colnames(selection_metric) <- c("Neighbors", "Models")

# Main loop
for (t in c(Ind_Cal, Ind_Eval)) {
  
  window_indices <- (t - tau):(t - 1)
  
  observed_window <- qmm[window_indices]
  dynamic_window  <- dynamic_weighted_avg[window_indices]
  
  best_window <- dynamic_window
  ind_na <- is.na(best_window)
  
  # Warm-up: fill missing best_window from obs, else from median(sim)
  if (any(ind_na)) {
    best_window[ind_na] <- observed_window[ind_na]
    ind_na <- is.na(best_window)
    if (any(ind_na)) {
      if (sum(ind_na) == 1) {
        best_window[ind_na] <- median(simulated[window_indices[ind_na], ], na.rm = TRUE)
      } else {
        best_window[ind_na] <- apply(simulated[window_indices[ind_na], ], 1, median, na.rm = TRUE)
      }
    }
  }
  
  fv_best_window <- best_window
  medTarget <- median(best_window, na.rm = TRUE)
  
  # Build Ind_Search quickly (exclude a contiguous range)
  lo <- min(window_indices) - 30 - tau
  hi <- max(window_indices) + 30 + tau
  
  in_excl <- (Ind_Cal >= lo) & (Ind_Cal <= hi)
  Ind_Search <- Ind_Cal[!in_excl]
  nS <- length(Ind_Search)
  
  # ---------- PASS 1: window-to-window metrics only ----------
  kge_windows    <- rep(NA_real_, nS)
  kgeinv_windows <- rep(NA_real_, nS)
  kgecomp_windows <- rep(NA_real_, nS)
  mae_windows    <- rep(NA_real_, nS)
  
  for (inc in seq_len(nS)) {
    i <- Ind_Search[inc]
    obs_cal <- qmm[(i - tau):(i - 1)]
    
    # window scores
    kge_windows[inc]    <- suppressWarnings(KGE_vec(fv_best_window, obs_cal)$KGE)
    kgeinv_windows[inc] <- suppressWarnings(KGE_vec((epsilon + fv_best_window)^-1, (epsilon + obs_cal)^-1)$KGE)
    kgecomp_windows[inc] <- (kge_windows[inc] + kgeinv_windows[inc])/2
    mae_windows[inc]    <- MAE_vec(fv_best_window, obs_cal)
  }
  
  # Choose neighbor metric
  metric_neighbors <- "MAE"
  kgeOK_windows <- is.finite(kge_windows)
  if (sum(kgeOK_windows) >= 0.8 * length(kge_windows)) {
    best_indices_windows <- order(kgecomp_windows, decreasing = TRUE)[1:knn]
    metric_neighbors <- "KGEcomp"
  } else {
    best_indices_windows <- order(mae_windows, decreasing = FALSE)[1:knn]
    metric_neighbors <- "MAE"
  }
  
  # ---------- PASS 2: model metrics only for selected windows ----------
  knn_eff <- length(best_indices_windows)
  sub_kge_models    <- matrix(NA_real_, nrow = knn_eff, ncol = ncol(simulated))
  sub_kgeinv_models <- matrix(NA_real_, nrow = knn_eff, ncol = ncol(simulated))
  sub_kgecomp_models <- matrix(NA_real_, nrow = knn_eff, ncol = ncol(simulated))
  sub_mae_models    <- matrix(NA_real_, nrow = knn_eff, ncol = ncol(simulated))
  
  for (j in seq_len(knn_eff)) {
    inc <- best_indices_windows[j]
    i <- Ind_Search[inc]
    
    sim_cal <- simulated[(i - tau):(i - 1), , drop = FALSE]
    obs_cal <- qmm[(i - tau):(i - 1)]
    
    # model scores for this neighbor window (vectorized across columns)
    sub_kge_models[j, ]    <- KGE_matrix(sim_cal, obs_cal)
    sub_kgeinv_models[j, ] <- KGEinv_matrix(sim_cal, obs_cal, epsilon)
    sub_kgecomp_models[j, ] <-(sub_kge_models[j, ]+sub_kgeinv_models[j, ])/2
    
    # vectorized MAE across models (fast)
    sub_mae_models[j, ] <- colMeans(abs(sim_cal - obs_cal), na.rm = TRUE)
  }
  
  
  # Average across selected windows
  if (knn_eff > 1) {
    sub_kgecomp_models_mean <- colMeans(sub_kgecomp_models, na.rm = TRUE)
    sub_mae_models_mean     <- colMeans(sub_mae_models, na.rm = TRUE)
  } else {
    sub_kgecomp_models_mean <- as.numeric(sub_kgecomp_models)
    sub_mae_models_mean     <- as.numeric(sub_mae_models)
  }
  
  # Choose model metric (KGE/KGEinv if available, else MAE)
  metric_models <- "MAE"
  kgeOK_models <- is.finite(sub_kgecomp_models_mean)
  if (sum(kgeOK_models) >= 0.8 * length(sub_kgecomp_models_mean)) {
    best_indices_mod <- order(sub_kgecomp_models_mean, decreasing = TRUE)[1:km]
    metric_models <- "KGEcomp"
  } else {
    best_indices_mod <- order(sub_mae_models_mean, decreasing = FALSE)[1:km]
    metric_models <- "MAE"
  }
  
  # Weights and weighted average
  mod_weight <- numeric(ncol(simulated))
  mod_weight[best_indices_mod] <- 1 / km
  
  dynamic_weights[t, ] <- mod_weight
  dynamic_weighted_avg[t] <- sum(simulated[t, ] * mod_weight, na.rm = TRUE)
  
  selection_metric[t, "Neighbors"] <- metric_neighbors
  selection_metric[t, "Models"]    <- metric_models
}

simulated_WA <- list(
  dates            = DatesR,
  weighted_avg     = dynamic_weighted_avg,
  weights          = dynamic_weights,
  selection_metric = selection_metric,
  searchPeriod     = c(CalStart, CalEnd),
  runPeriod        = c(CalStart, EvalEnd)
)

#####################
# Save results
#####################

dir_out = file.path(FUSE_path,dataset,spatialisation,inputdata,CritCal, outputfolder, catchment, "output")
if(! dir.exists(dir_out)){dir.create(dir_out, recursive = TRUE)}

save(simulated_WA, file = paste0(dir_out,"/", catchment,"_",tau,"_",km,"_",knn,".Rdata" ))

print(paste0("Saved file: ",dir_out,"/", catchment,"_",tau,"_",km,"_",knn,".Rdata" ))
