
#! ---------------------------------------------------------------------------------------
#!
#! Description       :
#!
#! Authors           : Cyril Thebault <cyril.thebault@ucalgary.ca>
#!
#! Creation date     : 2024-04-29 16:54:28
#! Modification date :
#!
#! Comments          :
#!
#! ---------------------------------------------------------------------------------------

#! ----------------------------- package loading

library(abind)
source("Metrics.R")

#! ---------------------------- bassins

##-----------------------------------------
##-------------- VARIABLES ----------------
##-----------------------------------------


dir_FUSE = "/work/comphyd_lab/users/cyril.thebault/Postdoc_Ucal/02_DATA/FUSE"

dir_BV = "/home/cyril.thebault/Postdoc_Ucal/02_DATA/BDD"

catchments = unlist(read.table(file.path(dir_BV, "liste_BV_CAMELS_559.txt")), use.names = FALSE)

dataset = "CAMELS"
spatialisation = "Lumped"
inputdata = "daymet"
CritCal = "KGECOMP"
TransfoCal="SA"


list_eval <- list()

for(catchment in catchments){
  
  print(catchment)
  
  ##-----------------------------------------
  ##------------- MERGE EVAL ----------------
  ##-----------------------------------------

  dir_catchment= file.path(dir_FUSE, dataset, spatialisation, inputdata,CritCal, TransfoCal, catchment)
  dir_eval= paste0(dir_catchment,"/evaluation")
  myfile = paste0(catchment, "_eval.Rdata")
  if(!file.exists(file.path(dir_eval,myfile))){next}
  TableEval = loadRData(file.path(dir_eval, myfile))
  
  # Append to list
  list_eval[[length(list_eval) + 1]] <- TableEval

}

TableEvalAll <- abind::abind(list_eval, along = 5)

# Add dimension names
dimnames(TableEvalAll) <- list(
  metric  = dimnames(list_eval[[1]])[[1]],
  transfo = dimnames(list_eval[[1]])[[2]],
  period  = dimnames(list_eval[[1]])[[3]],
  model   = dimnames(list_eval[[1]])[[4]],
  code    = catchments
)

save(TableEvalAll, file = paste0(dir_catchment, "/../Eval_SA.Rdata"))
