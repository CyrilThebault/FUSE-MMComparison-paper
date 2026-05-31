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


#! ----------------------------- workbook directory

loadRData <- function(file_name) {
  load(file_name)
  get(ls()[ls() != "file_name"])
}

source("/home/cyril.thebault/Postdoc_Ucal/03_CODES/Metrics.R")

FUSE_path = "/work/comphyd_lab/users/cyril.thebault/Postdoc_Ucal/02_DATA/FUSE"

dataset = "CAMELS"
spatialisation = "Lumped"
inputdata = "daymet"
CritCal = "KGECOMP"

# Simulations
Qsim = loadRData(file.path(FUSE_path,dataset,spatialisation,inputdata,CritCal,"1", "SimArray.Rdata"))

# Obervations
Qobs = loadRData(file = file.path(FUSE_path,dataset,spatialisation,inputdata,CritCal,"1","ObservedStreamflow.Rdata"))$Qobs


CalStart = "1989-01-01"
CalEnd = "1998-12-31"

EvalStart = "1999-01-01"
EvalEnd = "2009-12-31"


DatesR = as.Date(dimnames(Qsim)[[1]])
models = dimnames(Qsim)[[2]]
catchments = dimnames(Qsim)[[3]]

MetricsEval = c('KGE', 'NSE')
TransfoEval = c(1,-1)


# Combinations
comb2 <- combn(models, 2, simplify = FALSE)
comb3 <- combn(models, 3, simplify = FALSE)
# comb4 <- combn(models, 4, simplify = FALSE)
allcombinations <- c(comb2, comb3)
allcombinations_names = sapply(allcombinations, paste, collapse = "_")

############

args <- commandArgs(trailingOnly = TRUE)
id <- as.numeric(args[1])

catchment = catchments[id]
  
Ind_Cal = which(DatesR >= CalStart & DatesR <= CalEnd)
Ind_Eval = which(DatesR >= EvalStart & DatesR <= EvalEnd)

#prepare evaluation
TableEval <- array(
  NA, 
  dim = c(length(MetricsEval), length(TransfoEval) + 1,2, length(allcombinations_names)),
  dimnames = list(
    metric = MetricsEval,
    transfo = c(as.character(TransfoEval), "Comp"),
    period = c("Cal", "Eval"),
    names = allcombinations_names
  )
)


# Loop over each combination of decisions
for (comb_idx in seq_along(allcombinations)) {
  
  selected_decisions <- allcombinations[[comb_idx]]
  selected_decisions_name <- allcombinations_names[[comb_idx]]
  
  # cat(sprintf("Processing combination %d/%d: %s\n", comb_idx, length(allcombinations_names), selected_decisions_name))
  
  # Subset simulated
  simulated = rowMeans(Qsim[,selected_decisions,catchment])
  
  # Subset observed
  observed = Qobs[,catchment]
  
  
  for (metric in MetricsEval) {
    
    myFun <- get(metric)
    
    for (transfo in TransfoEval) {
      
      if (as.numeric(transfo) < 0) {
        epsilon <- mean(observed[Ind_Cal], na.rm = TRUE) / 100
      } else {
        epsilon <- 0
      }
      
      sim_cal_transf <- (epsilon + simulated[Ind_Cal])^transfo
      obs_cal_transf <- (epsilon + observed[Ind_Cal])^transfo
      
      sim_eval_transf <- (epsilon + simulated[Ind_Eval])^transfo
      obs_eval_transf <- (epsilon + observed[Ind_Eval])^transfo
      
      TableEval[metric, as.character(transfo), "Cal", selected_decisions_name] <- ifelse(metric == "KGE",
                                                                           myFun(sim_cal_transf, obs_cal_transf)$KGE,
                                                                           myFun(sim_cal_transf, obs_cal_transf))
      
      TableEval[metric, as.character(transfo), "Eval", selected_decisions_name] <- ifelse(metric == "KGE",
                                                                            myFun(sim_eval_transf, obs_eval_transf)$KGE,
                                                                            myFun(sim_eval_transf, obs_eval_transf))
    } # loop transfo
    
    TableEval[metric, "Comp", "Cal", selected_decisions_name] <- (TableEval[metric, "1", "Cal", selected_decisions_name]+TableEval[metric, "-1", "Cal", selected_decisions_name])/2
    
    TableEval[metric, "Comp", "Eval", selected_decisions_name] <- (TableEval[metric, "1", "Eval", selected_decisions_name]+TableEval[metric, "-1", "Eval", selected_decisions_name])/2
    
  } # loop metric
} # loop combinations


best_name <- dimnames(TableEval)$names[
  which.max(TableEval["KGE", "Comp", "Cal", ])
]

best_value <- max(TableEval["KGE", "Comp", "Cal", ])

# cat(sprintf("Best name: %s, KGE (Comp, Cal): %f\n", best_name, best_value))

best_decisions = unlist(strsplit(best_name, "_"))

bestsim = rowMeans(Qsim[,best_decisions,catchment])

#####################
# Save results
#####################

dir_out_sim = file.path(FUSE_path,dataset,spatialisation,inputdata,CritCal, "SA", catchment, "output")
if(! dir.exists(dir_out_sim)){dir.create(dir_out_sim, recursive = TRUE)}

save(bestsim, file = paste0(dir_out_sim,"/", catchment,".Rdata" ))

print(paste0("Saved file: ",dir_out_sim,"/", catchment,".Rdata" ))


dir_out_eval = file.path(FUSE_path,dataset,spatialisation,inputdata,CritCal, "SA", catchment, "evaluation")
if(! dir.exists(dir_out_eval)){dir.create(dir_out_eval, recursive = TRUE)}

save(TableEval, file = paste0(dir_out_eval,"/", catchment,"_eval.Rdata" ))

print(paste0("Saved file: ",dir_out_eval,"/", catchment,"_eval.Rdata" ))
