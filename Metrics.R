#! ---------------------------------------------------------------------------------------
#!
#! Description       :
#!
#! Authors           : Cyril Thebault <cyril.thebault@ucalgary.ca>
#!
#! Creation date     : 2025-04-29 16:54:28
#! Modification date :
#!
#! Comments          :
#!
#! ---------------------------------------------------------------------------------------

#! ----------------------------- function


loadRData <- function(file_name) {
  load(file_name)
  get(ls()[ls() != "file_name"])
}


KGE <- function(sim, obs) {
  
  suis_NA <- is.na(obs) | is.na(sim)
  
  qsim = sim[!suis_NA]
  qobs = obs[!suis_NA]
  
  sigma_Obs <- sd(qobs)
  sigma_Sim <- sd(qsim)
  alpha <- sigma_Sim/sigma_Obs
  ere <- cor(qobs, 
             qsim, method = "pearson")
  beta <- sum(qsim)/sum(qobs)
  return(list(KGE = 1-sqrt((ere-1)^2+(beta-1)^2+(alpha-1)^2),
              alpha = alpha,
              beta = beta,
              ere = ere))
} 

KGEp <- function(sim, obs) {
  
  suis_NA <- is.na(obs) | is.na(sim)
  
  qsim = sim[!suis_NA]
  qobs = obs[!suis_NA]
  
  sigma_Obs <- sd(qobs)
  sigma_Sim <- sd(qsim)
  alphap <- (sigma_Sim/mean(qsim))/(sigma_Obs/mean(qobs))
  ere <- cor(qobs, 
             qsim, method = "pearson")
  beta <- sum(qsim)/sum(qobs)
  return(list(KGEp = 1-sqrt((ere-1)^2+(beta-1)^2+(alphap-1)^2),
              alphap = alphap,
              beta = beta,
              ere = ere))
} 


NSE <- function(sim, obs) {
  
  suis_NA <- is.na(obs) | is.na(sim)
  
  qsim = sim[!suis_NA]
  qobs = obs[!suis_NA]
  
  numerator = sum((qobs - qsim)^2)
  denominator = sum((qobs - mean(qobs))^2)
  
  nse = 1-(numerator/denominator)
  return(nse)
} 


MAE <- function(sim, obs) {
  
  suis_NA <- is.na(obs) | is.na(sim)
  
  qsim = sim[!suis_NA]
  qobs = obs[!suis_NA]
  
  ae = abs(qobs - qsim)
  
  mae = mean(ae)
  
  return(mae)
} 

RMSE <- function(sim, obs) {
  
  suis_NA <- is.na(obs) | is.na(sim)
  
  qsim = sim[!suis_NA]
  qobs = obs[!suis_NA]
  
  se = (qobs - qsim)^2
  mse = mean(se)
  rmse = sqrt(mse)
  
  return(rmse)
} 

