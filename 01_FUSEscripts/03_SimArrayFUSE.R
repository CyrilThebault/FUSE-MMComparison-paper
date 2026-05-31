#! ---------------------------------------------------------------------------------------
#!
#! Description       :
#!
#! Authors           : Cyril Thebault <cyril.thebault@ucalgary.ca>
#!
#! Creation date     : 2024-10-21 17:35:01
#! Modification date :
#!
#! Comments          :
#!
#! ---------------------------------------------------------------------------------------

library(ncdf4)

#! ----------------------------- path

FUSE_path <-  "/work/comphyd_lab/users/cyril.thebault/Postdoc_Ucal/02_DATA/FUSE"
BDD_path <- "/home/cyril.thebault/Postdoc_Ucal/02_DATA/BDD"

#! ----------------------------  function

loadRData <- function(file_name) {
  load(file_name)
  get(ls()[ls() != "file_name"])
}

#! ----------------------------  main

catchments <-  unlist(read.table(file.path(BDD_path, "liste_BV_CAMELS_559.txt")))

decisions <- read.table(file.path(FUSE_path, "list_decision_78.txt"), header = TRUE, sep = ";")

decisionsID = as.character(decisions$ID)

dataset = "CAMELS"
spatialisation = "Lumped"
inputdata = "daymet"
CritCal = "KGECOMP"
TransfoCal="1" #needed for the tree even if not used with KGECOMP; should match the value set during FUSE execution

for(catchment in catchments){
  
  print("############")
  print(catchment)
  print("############")
  
  
  for(decision in decisionsID){
    
    print(decision)
    
    SIM_path = file.path(FUSE_path, dataset, spatialisation, inputdata, CritCal, TransfoCal, catchment, "output")
    myfile = paste0(catchment, "_",decision,"_runs_best.nc")
    
    if(!file.exists(file.path(SIM_path, myfile))){next}
      
    Qsim = nc_open(filename = file.path(SIM_path, myfile))
    
    if(!exists("Qsim_array")){
      
      # get time
      time_var <- ncvar_get(Qsim, "time")
      time_units <- ncatt_get(Qsim, "time", "units")$value
      time_split <- unlist(strsplit(time_units, " "))
      time_step <- time_split[1]
      
      if (time_step == "days") {
        ts = 86400  # 86400 seconds in a day
      } else if (time_step == "hours") {
        ts = 3600   # 3600 seconds in an hour
      } else if (time_step == "minutes") {
        ts = 60     # 60 seconds in a minute
      } else if (time_step == "seconds") {
        ts = 1
      } else {
        stop("Unsupported time step: ", time_step)
      }
      
      origin_date <- as.POSIXct(time_split[3], format = "%Y-%m-%d", tz = "UTC")
      
      DatesR = as.Date(origin_date + (time_var * ts))
      
      # Create an empty array with dimensions (dates, decisions, catchments)
      Qsim_array <- array(NA, dim = c(length(DatesR), length(decisionsID), length(catchments)),
                              dimnames = list(Dates = as.character(DatesR),
                                              Decisions = decisionsID,
                                              Catchments = catchments))
      
    }
    
    values = ncvar_get(Qsim, "q_routed")
    
    if(length(values) != length(DatesR)){next}
    
    Qsim_array[,decision,catchment] = values

    
    nc_close(Qsim)
  }
  
}


save(Qsim_array, file = file.path(FUSE_path, dataset, spatialisation, inputdata, CritCal, TransfoCal, "SimArray.Rdata"))