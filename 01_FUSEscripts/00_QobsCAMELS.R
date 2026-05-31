#! ---------------------------------------------------------------------------------------
#!
#! Description       :
#!
#! Authors           : Cyril Thebault <cyril.thebault@ucalgary.ca>
#!
#! Creation date     : 2024-07-15 14:54:28
#! Modification date :
#!
#! Comments          :
#!
#! ---------------------------------------------------------------------------------------

library(ncdf4)
library(sf)
library(reticulate)

#! ----------------------------- path

CAMELS_path = "/Users/cyrilthebault/Postdoc_Ucal/02_DATA/CAMELS" 

#! ----------------------------  function

loadRData <- function(file_name) {
  load(file_name)
  get(ls()[ls() != "file_name"])
}

#! ----------------------------  main

dir_BV = paste0("/Users/cyrilthebault/Postdoc_Ucal/02_DATA/BDD")
catchments = unlist(read.table(file.path(dir_BV, "liste_BV_CAMELS_559.txt")))

DatesR <- seq(from = as.Date("1987-01-01"),
              to = as.Date("2009-12-31"), 
              by = 'day')


# shp_catchments = read_sf(file.path(CAMELS_path, "basin_set_full_res", "HCDN_nhru_final_671.shp"))


Qobs = data.frame(matrix(nrow = length(DatesR), ncol = length(catchments)))
colnames(Qobs) = catchments
rownames(Qobs) = as.character(DatesR)

area <- setNames(rep(NA, length(catchments)), catchments)


for(catchment in catchments){
  
  #####################
  # Observed streamflow
  #####################
  
  # Get the ID number of the catchment as in CAMELS
  IDnumber = substr(catchment, 5,nchar(catchment))
  
  # Get catchment area 
  
  BM_forcing_path = file.path(CAMELS_path,"basin_timeseries_v1p2_metForcing_obsFlow","basin_dataset_public_v1p2","basin_mean_forcing", "daymet")
  BM_forcing_file = NULL
  for(folder in list.files(BM_forcing_path)){
    BM_forcing_file = list.files(file.path(BM_forcing_path, folder), pattern = IDnumber)
    if(length(BM_forcing_file) != 0){break}
  }
  
  catchment_area = as.numeric(readLines(paste0(BM_forcing_path,"/",folder,"/", BM_forcing_file))[3])*10^-6
  
  # Path of the observations
  observed_path = file.path(CAMELS_path,"basin_timeseries_v1p2_metForcing_obsFlow","basin_dataset_public_v1p2","usgs_streamflow")
  
  # Select the file corresponding to the catchment
  observed_file = NULL
  for(folder in list.files(observed_path)){
    observed_file = list.files(file.path(observed_path, folder), pattern = IDnumber)
    if(length(observed_file) != 0){break}
  }
  
  # Read the file
  observed = read.table(paste0(observed_path,"/",folder,"/", observed_file), header = FALSE, quote = "")
  colnames(observed) = c("Id", "Year", "Month", "Day", "Qcfs", "Flag")
  
  # Create column with Date
  observed$Date = as.Date(paste(observed$Year, observed$Month, observed$Day, sep = '-'))
  
  # Convert streamflow
  observed$qm3s =  observed$Qcfs*0.028317
  observed$qmm = observed$qm3s*86400/(catchment_area*10^3) 
  
  qmm = observed$qmm[match(DatesR, observed$Date)]
  
  Qobs[,catchment] = qmm
  
  area[catchment] = catchment_area

  print(catchment)
}

list_obs = list(Qobs_mm = Qobs, area_km2 = area)

save(list_obs, file = "/Users/cyrilthebault/Postdoc_Ucal/02_DATA/FUSE/ObservedStreamflow.Rdata")
