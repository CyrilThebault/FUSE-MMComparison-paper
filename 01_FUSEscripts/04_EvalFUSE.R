
#! ---------------------------------------------------------------------------------------
#!
#! Description       :
#!
#! Authors           : Cyril Thebault <cyril.thebault@inrae.fr>
#!
#! Creation date     : 2024-10-26 16:54:28
#! Modification date :
#!
#! Comments          :
#!
#! ---------------------------------------------------------------------------------------


#! ----------------------------- path definition

#! -------------- sources

#! source directory
path_src <- "/home/cyril.thebault/Postdoc_Ucal/02_DATA"
# path_src <- "../02_DATA"

#! -------------- results

#! result directory
path_res <- "/home/cyril.thebault/Postdoc_Ucal/02_DATA"
# path_res <- "../02_DATA"


#! ----------------------------- workbook directory

# setwd("D:/cyril.thebault/03_CODES")
# setwd("/nfs/home/thebaultc/03_CODES")

#! ----------------------------- package loading

library(ncdf4)

source(file = "Metrics.R")

#! ---------------------------- bassins

##-----------------------------------------
##-------------- VARIABLES ----------------
##-----------------------------------------

dir_FUSE = "/work/comphyd_lab/users/cyril.thebault/Postdoc_Ucal/02_DATA/FUSE"

dir_BV = "/home/cyril.thebault/Postdoc_Ucal/02_DATA/BDD"

catchments = unlist(read.table(file.path(dir_BV, "liste_BV_CAMELS_559.txt")), use.names = FALSE)

models = read.table(paste0(dir_FUSE, "/list_decision_78.txt"), sep = ";", header = TRUE)[,1]

MetricsEval = c('KGE', 'NSE')
TransfoEval = c(1,-1)

WU = 2

CalStart = "1989-01-01"
CalEnd = "1998-12-31"

SimStart = "1989-01-01"
SimEnd = "2009-12-31"

EvalStart = "1999-01-01"
EvalEnd = "2009-12-31"

dataset = "CAMELS"
spatialisation = "Lumped"
inputdata = "daymet"
CritCal = "KGECOMP"
TransfoCal="1" #needed for the tree even if not used with KGECOMP; should match the value set during FUSE execution

##-----------------------------------------
##----------------- HPC -------------------
##-----------------------------------------

args <- commandArgs(trailingOnly = TRUE)
id <- as.numeric(args[1])
# id = 1

catchment = catchments[id]

print("---------------")
print(catchment)
print("---------------")

##-----------------------------------------
##----------------- DATA ------------------
##-----------------------------------------

# Dates
DatesR <- seq(from = as.Date(gsub(format(as.Date(min(CalStart, SimStart, EvalStart)), "%Y"),
                                  as.numeric(format(as.Date(min(CalStart, SimStart, EvalStart)), "%Y"))-WU,
                                  min(CalStart, SimStart, EvalStart))),
              to = as.Date(max(CalEnd, SimEnd, EvalEnd)), 
              by = 'day')


# Observed streamflow
myfileobs = paste0(dir_FUSE,"/",dataset, "/",spatialisation, "/",inputdata, "/",CritCal, "/",TransfoCal,"/",catchment,"/input/",catchment,"_input.nc")
if(!file.exists(myfileobs)){stop("The observation file does not exist.")}
obs = nc_open(myfileobs)

Qobs = ncvar_get(obs, "q_obs")

nc_close(obs)

Qobs[is.nan(Qobs)] = NA
Qobs[Qobs< 0] = NA

Qobs = data.frame(Qobs)
colnames(Qobs) = catchment
rownames(Qobs) = DatesR

for(mod in models){
  ##-----------------------------------------
  ##-------------- SIMULATION ---------------
  ##-----------------------------------------
  
  myfilesim = paste0(dir_FUSE,"/",dataset, "/",spatialisation, "/",inputdata, "/",CritCal, "/",TransfoCal,"/",catchment,"/output/",catchment,"_",mod,"_runs_best.nc")
  if(!file.exists(myfilesim)){
    print(paste0("The simulation file for model decision ", mod, " does not exist."))
    next
  }
  sim = nc_open(myfilesim)
  
  Qsim = ncvar_get(sim, "q_routed")
  
  nc_close(sim)
  
  if(length(Qsim)==0){
    print(paste0("The simulation file for model decision ", mod, " is empty."))
    next
  }
  
  Qsim[is.nan(Qsim)] = NA
  Qsim[Qsim< 0] = NA
  
  Qsim = data.frame(Qsim)
  colnames(Qsim) = catchment
  rownames(Qsim) = DatesR
  
  ##-----------------------------------------
  ##-------------- EVALUATION ---------------
  ##-----------------------------------------
  IndPeriod_WarmUp <- seq(
    which(as.POSIXct(DatesR) == as.POSIXct(as.Date(gsub(format(as.Date(SimStart), "%Y"),
                                                        as.numeric(format(as.Date(SimStart), "%Y")) - WU,
                                                        SimStart)))), # Set aside warm-up period
    which(as.POSIXct(DatesR) == as.POSIXct(as.Date(SimStart)-1)) # Until the end of the time series
  )
  
  
  IndPeriod_Cal <- seq(
    which(as.POSIXct(DatesR) == as.POSIXct(as.Date(CalStart))), # Set aside warm-up period
    which(as.POSIXct(DatesR) == as.POSIXct(as.Date(CalEnd))) # Until the end of the time series
  )
  
  IndPeriod_Eval <- seq(
    which(as.POSIXct(DatesR) == as.POSIXct(as.Date(EvalStart))), # Set aside warm-up period
    which(as.POSIXct(DatesR) == as.POSIXct(as.Date(EvalEnd))) # Until the end of the time series
  )
  
  
  mygrid = expand.grid(MetricsEval, TransfoEval)
  colnames(mygrid) = c("Metrics", "Transformations")
  combinations <- apply(mygrid, 1, function(row) paste(row, collapse = " : "))
  TableEval = data.frame(matrix(NA, nrow = length(combinations), ncol = 2))
  colnames(TableEval) = c("Cal", "Eval")
  rownames(TableEval) = combinations
  
  for(i in 1:nrow(mygrid)){
    
    metric = as.character(mygrid$Metrics[i])
    myFun = get(metric)
    
    transfo = as.numeric(mygrid$Transformations[i])
    
    if(transfo < 0 ){
      epsilon = mean(Qobs[IndPeriod_Cal, catchment], na.rm = TRUE)/100
    } else{
      epsilon = 0
    }
    
    TableEval[i, 'Cal'] = ifelse(metric == "KGE",
                                 myFun(sim = (epsilon + Qsim[IndPeriod_Cal, catchment])^transfo,
                                       obs = (epsilon + Qobs[IndPeriod_Cal, catchment])^transfo)$KGE,
                                 myFun(sim = (epsilon + Qsim[IndPeriod_Cal,catchment])^transfo,
                                       obs = (epsilon + Qobs[IndPeriod_Cal, catchment])^transfo)
    )
    
    TableEval[i, 'Eval'] = ifelse(metric == "KGE",
                                  myFun(sim = (epsilon + Qsim[IndPeriod_Eval, catchment])^transfo,
                                        obs = (epsilon + Qobs[IndPeriod_Eval, catchment])^transfo)$KGE,
                                  myFun(sim = (epsilon + Qsim[IndPeriod_Eval, catchment])^transfo,
                                        obs = (epsilon + Qobs[IndPeriod_Eval, catchment])^transfo)
    )
    
  }
  
  
  dir_eval= paste0(dir_FUSE, "/",dataset, "/",spatialisation, "/",inputdata, "/",CritCal, "/",TransfoCal,"/",catchment, "/evaluation")
  if (!dir.exists(dir_eval)) {
    dir.create(dir_eval, recursive = TRUE)
  }
  save(TableEval, file = paste0(dir_eval, "/",catchment,"_",mod,"_eval.Rdata"))
  
  print(mod)
}