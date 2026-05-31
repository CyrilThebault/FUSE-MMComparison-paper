
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


#! ----------------------------- path definition

#! -------------- sources

#! source directory
path_src <- "/home/cyril.thebault/Postdoc_Ucal/02_DATA"

#! -------------- results

#! result directory
path_res <- "/home/cyril.thebault/Postdoc_Ucal/02_DATA"


#! ----------------------------- workbook directory


#! ----------------------------- package loading

source("Metrics.R")

#! ---------------------------- bassins

##-----------------------------------------
##-------------- VARIABLES ----------------
##-----------------------------------------


dir_FUSE = "/work/comphyd_lab/users/cyril.thebault/Postdoc_Ucal/02_DATA/FUSE"

dir_BV = "/home/cyril.thebault/Postdoc_Ucal/02_DATA/BDD"

catchments = unlist(read.table(file.path(dir_BV, "liste_BV_CAMELS_559.txt")), use.names = FALSE)
models = read.table(paste0(dir_FUSE, "/list_decision_78.txt"), sep = ";", header = TRUE)[,1]

dataset = "CAMELS"
spatialisation = "Lumped"
inputdata = "daymet"
CritCal = "KGECOMP"
TransfoCal="1" #needed for the tree even if not used with KGECOMP; should match the value set during FUSE execution

CritEval = c("KGE", "NSE")
TransfoEval = c(1,-1)

##-----------------------------------------
##----------------- TREE ------------------
##-----------------------------------------

mydf = data.frame(expand.grid(models, catchments))
colnames(mydf) = c('ModelDecisions', 'Codes')
mydf = mydf[,c('Codes', 'ModelDecisions')]

combinations <- expand.grid(c("Cal", "Eval"), CritEval, TransfoEval)
column_names <- apply(combinations, 1, function(row) paste(row, collapse = " : "))
mydf[,column_names] <- NA

for(id in 1:nrow(mydf)){
  
  mod = as.character(mydf$ModelDecisions[id])
  catchment = as.character(mydf$Codes[id])

  ##-----------------------------------------
  ##------------- MERGE EVAL ----------------
  ##-----------------------------------------
  
  dir_merge = dir_FUSE
  
  dir_catchment= file.path(dir_merge, dataset, spatialisation, inputdata,CritCal, TransfoCal, catchment)
  
  dir_eval= paste0(dir_catchment,"/evaluation")
  
  myfile = paste0(catchment,"_", mod, "_eval.Rdata")
  
  if(!file.exists(file.path(dir_eval,myfile))){next}
  
  TableEval = loadRData(file.path(dir_eval, myfile))
  
  # Create an empty list to store the new row names and values
  new_data <- list()
  
  # Loop through each column and each row to combine names
  for (col in names(TableEval)) {
    for (row in rownames(TableEval)) {
      # Combine column name and row name
      new_row_name <- paste(col, row, sep = " : ")
      # Extract the value
      new_data[[new_row_name]] <- TableEval[row, col]
    }
  }
  
  # Convert the list to a single-column data frame
  single_column_df <- data.frame(Value = unlist(new_data))
  
  
  mydf[id, column_names] = single_column_df[column_names,]
  
  print(paste0(catchment, " --- ", mod))
}

save(mydf, file = paste0(dir_catchment, "/../Eval_FUSE.Rdata"))
