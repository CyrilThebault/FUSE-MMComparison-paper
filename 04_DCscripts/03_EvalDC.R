
#! ---------------------------------------------------------------------------------------
#!
#! Description       :
#!
#! Authors           : Cyril Thebault <cyril.thebault@ucalgary.ca>
#!
#! Creation date     : 2025-01-08 16:54:28
#! Modification date :
#!
#! Comments          :
#!
#! ---------------------------------------------------------------------------------------

#! ----------------------------- package loading

source("Metrics.R")


#! ---------------------------- bassins

##-----------------------------------------
##-------------- VARIABLES ----------------
##-----------------------------------------

dir_FUSE = "/project/6079554/thebault/Postdoc_Ucal/02_DATA/FUSE"

dir_BV = "/home/thebault/Postdoc_Ucal/02_DATA/BDD"

catchments = unlist(read.table(file.path(dir_BV, "liste_BV_CAMELS_559.txt")), use.names = FALSE)

# models = read.table(paste0(dir_FUSE, "/list_decision_78.txt"), sep = ";", header = TRUE)[,1]

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
TransfoCal="Comp_WA_q_KGEcomp"

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


# Observed streamflow
Qobs_all = loadRData(file = file.path(dir_FUSE,dataset,spatialisation,inputdata,CritCal,"1","ObservedStreamflow.Rdata"))$Qobs

Qobs = Qobs_all[,catchment, drop=FALSE]
Qobs$Date = DatesR


pathWA = paste0(dir_FUSE,"/",dataset, "/",spatialisation, "/",inputdata, "/",CritCal, "/",TransfoCal,"/",catchment,"/output")
filesWA = list.files(pathWA)

for(myfilesim in filesWA){
  ##-----------------------------------------
  ##-------------- SIMULATION ---------------
  ##-----------------------------------------
  
  sim = loadRData(file.path(pathWA, myfilesim))
  
  
  Qsim_ini = sim$weighted_avg
  
  Qsim_ini[is.nan(Qsim_ini)] = NA
  Qsim_ini[Qsim_ini< 0] = NA
  
  Qsim_ini = data.frame(Qsim_ini)
  colnames(Qsim_ini) = catchment
  Qsim_ini$Date = sim$dates
  
  Qsim = merge(Qobs["Date"], Qsim_ini, by = "Date", all.x = TRUE)
  colnames(Qsim) = c("Date", catchment)
  
  
  ##-----------------------------------------
  ##-------------- EVALUATION ---------------
  ##-----------------------------------------
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
  myfileeval <- sub("\\.Rdata$", "_eval.Rdata", myfilesim)
  save(TableEval, file = file.path(dir_eval, myfileeval))
  
  print(myfileeval)
}