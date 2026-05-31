#! ---------------------------------------------------------------------------------------
#!
#! Description       :
#!
#! Authors           : Cyril Thebault <cyril.thebault@ucalgary.ca>
#!
#! Creation date     : 2024-02-26 16:54:28
#! Modification date :
#!
#! Comments          :
#!
#! ---------------------------------------------------------------------------------------

#! ----------------------------- package loading

library(ncdf4)
library(CFtime)
library(convertr)
library(sf)
library(stars)
library(dplyr)
library(airGR)

#! ----------------------------  File manager

FUSE_path = "/work/comphyd_lab/users/cyril.thebault/Postdoc_Ucal/02_DATA/FUSE"   # (01) Path of the workflow
CAMELS_path = "/home/cyril.thebault/Postdoc_Ucal/02_DATA/CAMELS/basin_timeseries_v1p2_metForcing_obsFlow/basin_dataset_public_v1p2"         # (02) Path of the CAMELS-spat database\
dir_BV = "/home/cyril.thebault/Postdoc_Ucal/02_DATA/BDD"

models = unlist(read.table(file.path(FUSE_path, "list_decision_78.txt"), sep = ";", header = TRUE)$ID)

catchments = unlist(read.table(file.path(dir_BV, "liste_BV_CAMELS_559.txt")))

SimStart = "1989-01-01"
SimEnd = "2009-12-31"

CalStart = "1989-01-01"
CalEnd = "1998-12-31"

WU = 2

shp_catchments = read_sf(paste0(CAMELS_path,"/../../basin_set_full_res/HCDN_nhru_final_671.shp"))

CritCal = "KGECOMP"
TransfoCal = "1"
dataset = "CAMELS"
spatialisation = "Lumped"
inputdata = "daymet"

#! ----------------------------  catchments

args <- commandArgs(trailingOnly = TRUE)
id <- as.numeric(args[1])

catchment = catchments[id]
  
#! ----------------------------  create folder
mydemo = "fuse_template"
DemoFolder = file.path(FUSE_path, mydemo)
CatchmentFolder = file.path(FUSE_path, dataset, spatialisation, inputdata, CritCal, TransfoCal, catchment)

#create folder tree
if(!dir.exists(CatchmentFolder)){
  dir.create(path = file.path(CatchmentFolder, "input"), recursive = TRUE)
  dir.create(path = file.path(CatchmentFolder, "output"), recursive = TRUE)
  dir.create(path = file.path(CatchmentFolder, "settings/fuse_zDecisions"), recursive = TRUE)
}

if(length(list.files(file.path(CatchmentFolder, "settings")))==1){
  #copy settings files
  settings_files = list.files(file.path(DemoFolder,"settings"))
  
  invisible(file.copy(from = paste0(DemoFolder,"/settings/", settings_files), 
                      to = file.path(CatchmentFolder,"settings"),
                      recursive = FALSE)
  )
  
}

DecisionFiles = paste0("fuse_zDecisions_",models,".txt")
if(!any(DecisionFiles %in% list.files(file.path(CatchmentFolder, "settings/fuse_zDecisions")))){
  
  invisible(file.copy(from = paste0(DemoFolder,"/settings/fuse_zDecisions/", DecisionFiles), 
                      to = paste0(CatchmentFolder,"/settings/fuse_zDecisions"), recursive = TRUE)
  )
  
}

#! ----------------------------  get input data

if((!paste0(catchment,"_elev_bands.nc") %in% list.files(file.path(CatchmentFolder, "input")))|
   (!paste0(catchment,"_input.nc") %in% list.files(file.path(CatchmentFolder, "input")))){
  
  IDnumber = substr(catchment, 5,nchar(catchment))
  
  
  shp_catchment = shp_catchments[shp_catchments$hru_id == as.numeric(IDnumber),]
  lat = shp_catchment$lat_cen
  lon = shp_catchment$lon_cen
  
  if(!paste0(catchment,"_input.nc") %in% list.files(file.path(CatchmentFolder, "input"))){
    
    ####### Dates
    DatesR <- seq(from = as.Date(gsub(format(as.Date(min(SimStart, CalStart)), "%Y"),
                                      as.numeric(format(as.Date(min(SimStart, CalStart)), "%Y"))-WU,
                                      min(SimStart, CalStart))),
                  to = as.Date(max(SimEnd, CalEnd)), 
                  by = 'day')
    
    ###### Precipitation, temperature & potential evapotranspiration
    
    # Basin mean
    
    BM_forcing_path = paste0(CAMELS_path,"/basin_mean_forcing/", inputdata)
    BM_forcing_file = NULL
    for(folder in list.files(BM_forcing_path)){
      BM_forcing_file = list.files(file.path(BM_forcing_path, folder), pattern = IDnumber)
      if(length(BM_forcing_file) != 0){break}
    }
    
    BM_forcing = read.table(paste0(BM_forcing_path,"/",folder,"/", BM_forcing_file), header = FALSE, quote = "", skip = 4)
    colnames(BM_forcing) = c('Year', 'Month', 'Day', 'Hour', 'Dayl[s]', 'P[mm/d]', 'Srad[W/m2]', 'SWE[mm]', 'Tmax[C]', 'Tmin[C]', 'Vp[Pa]')
    BM_forcing$Date = as.Date(paste(BM_forcing$Year, BM_forcing$Month, BM_forcing$Day, sep = '-'))
    
    BM_forcing[,'Tmean[C]'] = rowMeans(BM_forcing[, c('Tmax[C]', 'Tmin[C]')])
    BM_forcing[,'P[mm/d]'][BM_forcing[,'P[mm/d]'] < 0] = 0
    BM_forcing[,'PET_Oudin[mm/d]'] <- PE_Oudin(JD = as.POSIXlt(BM_forcing[,"Date"])$yday + 1,
                                               Temp = BM_forcing[,'Tmean[C]'],
                                               Lat = lat, LatUnit = "deg")
    
    BM_subforcing = BM_forcing[BM_forcing$Date %in% DatesR,]
    
    catchment_area = as.numeric(readLines(paste0(BM_forcing_path,"/",folder,"/", BM_forcing_file))[3])*10^-6
    
    ###### Observed streamflow
    observed_path = paste0(CAMELS_path,"/usgs_streamflow")
    observed_file = paste0(IDnumber, "_streamflow_qc.txt")
    observed = read.table(paste0(observed_path,"/",folder,"/", observed_file), header = FALSE, quote = "")
    
    colnames(observed) = c("Id", "Year", "Month", "Day", "Qcfs", "Flag")
    
    observed$Date = as.Date(paste(observed$Year, observed$Month, observed$Day, sep = '-'))
    
    observed$qm3s =  observed$Qcfs*0.028317
    
    catchment_area = as.numeric(readLines(paste0(BM_forcing_path,"/",folder,"/", BM_forcing_file))[3])*10^-6
    
    observed$qmm = observed$qm3s*86400/(catchment_area*10^3) 
    
    qmm = observed$qmm[match(DatesR, observed$Date)]
    
    qmm[is.nan(qmm)] = NA
    qmm[qmm< 0] = NA
    
    BM_subforcing[,'Qobs[mm/d]'] = qmm
    
    
    #! ----------------------------  create input file
    
    inputname <- paste0(CatchmentFolder, "/input/",catchment,"_input.nc")
    
    
    # define dimensions
    latdim <- ncdim_def(name = "latitude", units = "degreesN", vals = lat, longname = "latitude", create_dimvar=TRUE)
    londim <- ncdim_def(name = "longitude", units = "degreesE", vals = lon, longname = "longitude", create_dimvar=TRUE)
    timedim <-  ncdim_def(name = "time", units = paste0("days since ", min(DatesR)), vals = (1:nrow(BM_subforcing))-1, longname = "time",create_dimvar=TRUE, unlim = TRUE)
    
    # define variables
    PET_def <- ncvar_def(name = "pet", 
                         units = "mm/day",
                         dim = list(londim, latdim, timedim),
                         longname = "Potential evaportanspiration estimated using Oudin et al., 2005, JoH",
                         prec = "double"
    )
    
    P_def <- ncvar_def(name = "pr", 
                       units = "mm/day",
                       dim = list(londim, latdim, timedim),
                       longname = "Mean daily precipitation",
                       prec = "double"
    )
    
    Qobs_def <- ncvar_def(name = "q_obs", 
                          units = "mm/day",
                          dim = list(londim, latdim, timedim),
                          longname = "Mean observed daily discharge",
                          prec = "double"
    )
    
    T_def <- ncvar_def(name = "temp", 
                       units = "degC",
                       dim = list(londim, latdim, timedim),
                       longname = "Mean daily temperature",
                       prec = "double"
    )
    
    
    ncinput <- nc_create(inputname,list(PET_def, P_def, Qobs_def, T_def), force_v4 = TRUE)
    
    # put variables
    ncvar_put(ncinput,PET_def, BM_subforcing[,'PET_Oudin[mm/d]'])
    ncvar_put(ncinput,P_def, BM_subforcing[,'P[mm/d]'])
    ncvar_put(ncinput,Qobs_def,BM_subforcing[,'Qobs[mm/d]'])
    ncvar_put(ncinput,T_def,BM_subforcing[,'Tmean[C]'])
    
    
    # add global attributes
    ncatt_put(ncinput,0,"author", "Cyril Thébault")
    ncatt_put(ncinput,0,"date",paste0("Created ", format(Sys.time(), format = "%Y/%m/%d %H:%M:%S")))
    ncatt_put(ncinput,0,"institution",'University of Calgary')
    
    nc_close(ncinput)
  }
  
  #! ----------------------------  get elevation band data
  
  if(!paste0(catchment,"_elev_bands.nc") %in% list.files(file.path(CatchmentFolder, "input"))){
    
    # Elevation bands
    
    EB_forcing_path = paste0(CAMELS_path,"/elev_bands_forcing/", inputdata)
    EB_forcing_file_list = list.files(file.path(EB_forcing_path, folder), pattern = IDnumber)
    EB_forcing_file_list = EB_forcing_file_list[grepl(".list", EB_forcing_file_list)]
    
    EB_forcing_file_df = read.table(paste0(EB_forcing_path,"/",folder,"/", EB_forcing_file_list), header = FALSE, skip = 1)
    colnames(EB_forcing_file_df) = c("FileName","Area[m2]")
    
    subcatchment_area = EB_forcing_file_df[,"Area[m2]"]*10^-6
    
    
    #! ----------------------------  create elevation band file
    
    elevname <- paste0(CatchmentFolder, "/input/",catchment,"_elev_bands.nc")
    
    
    # define dimensions
    latdim <- ncdim_def(name = "latitude", units = "degreesN", vals = lat, longname = "latitude", create_dimvar=TRUE)
    londim <- ncdim_def(name = "longitude", units = "degreesE", vals = lon, longname = "longitude", create_dimvar=TRUE)
    elevdim <-  ncdim_def(name = "elevation_band", units = "-", vals = 1:nrow(EB_forcing_file_df), longname = "elevation_band",create_dimvar=TRUE)
    
    
    # define variables
    AreaFrac_def <- ncvar_def(name = "area_frac", 
                              units = "-",
                              dim = list(londim, latdim, elevdim),
                              longname = "Fraction of the catchment covered by each elevation band",
                              prec = "double"
    )
    
    MeanElev_def <- ncvar_def(name = "mean_elev", 
                              units = "m asl",
                              dim = list(londim, latdim, elevdim),
                              longname = "Mid-point elevation of each elevation band",
                              prec = "double"
    )
    
    PrecFrac_def <- ncvar_def(name = "prec_frac", 
                              units = "-",
                              dim = list(londim, latdim, elevdim),
                              longname = "Fraction of catchment precipitation that falls on each elevation band - same as area_frac",
                              prec = "double"
    )
    
    
    ncelev <- nc_create(elevname,list(AreaFrac_def, MeanElev_def, PrecFrac_def), force_v4 = TRUE)
    
    # put variables
    ncvar_put(ncelev,AreaFrac_def, subcatchment_area/sum(subcatchment_area))
    ncvar_put(ncelev,MeanElev_def, as.numeric(substr(EB_forcing_file_df[,"FileName"], 20,22))*100+50)
    ncvar_put(ncelev,PrecFrac_def,subcatchment_area/sum(subcatchment_area))
    
    # add global attributes
    ncatt_put(ncelev,0,"author", "Cyril Thébault")
    ncatt_put(ncelev,0,"date",paste0("Created ", format(Sys.time(), format = "%Y/%m/%d %H:%M:%S")))
    ncatt_put(ncelev,0,"institution",'University of Calgary')
    
    nc_close(ncelev)
  }
}

#! ----------------------------  fm_catch file

for(mod in models){
  
  FM = readLines( file.path(DemoFolder, "fm_catch.txt"))
  
  FM[which(grepl("! SETNGS_PATH", FM))] = gsub( "/my/path/to/fuse/catchment/settings/" , file.path("/dev/shm",catchment,"settings/") , FM[which(grepl("! SETNGS_PATH", FM))] )
  FM[which(grepl("! INPUT_PATH", FM))] = gsub( "/my/path/to/fuse/catchment/input/" , file.path("/dev/shm",catchment,"input/") , FM[which(grepl("! INPUT_PATH", FM))] )
  FM[which(grepl("! OUTPUT_PATH", FM))] = gsub( "/my/path/to/fuse/catchment/output/" , file.path("/dev/shm",catchment,"output/") , FM[which(grepl("! OUTPUT_PATH", FM))] )
  
  FM[which(grepl("! Q_ONLY", FM))] = gsub( "FALSE" , "TRUE" , FM[which(grepl("! Q_ONLY", FM))] )
  
  FM[which(grepl("! M_DECISIONS", FM))] = gsub( "902" , mod , FM[which(grepl("! M_DECISIONS", FM))] )
  FM[which(grepl("! FMODEL_ID", FM))] = gsub( "902" , mod , FM[which(grepl("! FMODEL_ID", FM))] )
  
  FM[which(grepl("! date_start_sim", FM))] = gsub( "2000-10-01" , 
                                                   as.Date(gsub(format(as.Date(min(SimStart, CalStart)), "%Y"),
                                                                as.numeric(format(as.Date(min(SimStart, CalStart)), "%Y"))-WU,
                                                                min(SimStart, CalStart))), 
                                                   FM[which(grepl("! date_start_sim", FM))] )
  FM[which(grepl("! date_end_sim", FM))] = gsub( "2005-09-30" , SimEnd , FM[which(grepl("! date_end_sim", FM))] )
  FM[which(grepl("! date_start_eval", FM))] = gsub( "2001-10-01" , CalStart , FM[which(grepl("! date_start_eval", FM))] )
  FM[which(grepl("! date_end_eval", FM))] = gsub( "2005-09-30" , CalEnd , FM[which(grepl("! date_end_eval", FM))] )
  
  FM[which(grepl("! METRIC", FM))] = gsub( "RMSE" , CritCal , FM[which(grepl("! METRIC", FM))] )
  FM[which(grepl("! TRANSFO", FM))] = gsub( "1" , TransfoCal , FM[which(grepl("! TRANSFO", FM))] )
  
  FM[which(grepl("! MAXN", FM))] = gsub( "20" , 10000 , FM[which(grepl("! MAXN", FM))] )
  
  writeLines( FM , paste0(CatchmentFolder,"/",catchment,"_",mod,".txt") )
  
}

print(paste0(catchment, " DONE"))
