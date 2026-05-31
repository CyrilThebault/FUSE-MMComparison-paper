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

setwd("/Users/cyrilthebault/Postdoc_Ucal/02_DATA/FUSE-MMComparison-paper")

#! ----------------------------- package loading

library(sf)
library(dplyr)
library(ggplot2)
library(ggspatial)
library(ggpubr)
library(grid)
library(dplyr)
library(tidyr)
library(viridis)
library(patchwork)
library(tibble)


source(file = "Metrics.R")

#! ---------------------------- bassins

##-----------------------------------------
##-------------- VARIABLES ----------------
##-----------------------------------------

######## Code catchments

catchments = as.character(unname(unlist(read.table(file.path("00_DATA", "liste_BV_CAMELS_559.txt")))))

######## Shapefiles

## catchments outlets
csv_outlets_DB <- readLines(file.path("00_DATA", "gauge_information.txt"))
csv_outlets_DB = strsplit(csv_outlets_DB, "\t")
csv_outlets_DB[[1]] = c("HUC_02", "GAGE_ID", "GAGE_NAME", "LAT", "LON", "DRAINAGE AREA (KM^2)")
csv_outlets_DB = data.frame(matrix(unlist(csv_outlets_DB), nrow=length(csv_outlets_DB), byrow=TRUE),stringsAsFactors=FALSE)
colnames(csv_outlets_DB) = csv_outlets_DB[1,]
csv_outlets_DB = csv_outlets_DB[! csv_outlets_DB$HUC_02 == "HUC_02",]
csv_outlets_DB$Station_lon = csv_outlets_DB$LON
csv_outlets_DB$Station_lat = csv_outlets_DB$LAT
csv_outlets_DB$Station_id = csv_outlets_DB$GAGE_ID
csv_outlets_DB$Country = "USA"

shp_outlets_DB <- sf::st_as_sf(csv_outlets_DB, coords = c('Station_lon','Station_lat'), crs = 4326)

modif = nchar(shp_outlets_DB$Station_id) == 7
shp_outlets_DB$Station_id[shp_outlets_DB$Country == 'USA' & modif] = paste0("0", shp_outlets_DB$Station_id[shp_outlets_DB$Country == 'USA' & modif])

shp_outlets_DB$Codes = paste(shp_outlets_DB$Country, shp_outlets_DB$Station_id, sep='_')

## rivers

## north-america
NorthAm = file.path("Shp","boundaries_p_2021_v3.shp")
NorthAm = read_sf(NorthAm)
NorthAm = st_transform(NorthAm, crs = 4326)
NorthAm = NorthAm[NorthAm$COUNTRY == "USA",]

drop_these = c('US-AK', 'US-HI', 'US-PR', 'US-VI')
NorthAm = NorthAm[! NorthAm$STATEABB %in% drop_these,]

######## Subset

shp_outlets = subset(shp_outlets_DB, Codes %in% catchments)

####### /!\ Remove catchments that fails during multi-model approaches
# 3 catchments fails in the benchmark (unable to calculate KGEcomp metric)
# 1 catchment  fail  in the mosaic based on performance (unable to calculate KGEcomp metric)
# 7 catchments fails in the mosaic based on performance equivalence (unable to do the bootstrap)
# 9 catchments fails in the spatially and temporally static combination (unable to calculate KGEcomp metric on one of the models)
# 1 catchment  fail  in the spatially variable and temporally static combination (unable to calculate KGEcomp metric on one of the models)
# 1 catchment  fail  in the dynamic combination (unable to calculate KGEcomp metric)

catchments_fail = c("USA_02427250", "USA_04056500", "USA_05062500", "USA_05412500", 
                    "USA_06354000", "USA_06360500", "USA_06441500", "USA_06447000",
                    "USA_06452000", "USA_06468250", "USA_07263295", "USA_07362587",
                    "USA_09484600", "USA_11151300", "USA_12141300")

shp_outlets = subset(shp_outlets, !Codes %in% catchments_fail)

###### Add performance to the subset

dir_FUSE = file.path("00_DATA")

# FUSE Comp
Eval_FUSE_long_Comp = loadRData(file.path(dir_FUSE, "KGEcomp", "Eval_FUSE.Rdata"))%>%
  mutate(`Cal : KGEcomp` = (`Cal : KGE :  1` + `Cal : KGE : -1`) / 2,
         `Eval : KGEcomp` = (`Eval : KGE :  1` + `Eval : KGE : -1`) / 2)

Eval_FUSE_Comp <- data.frame(Eval_FUSE_long_Comp %>%
                               dplyr::select(Codes, ModelDecisions, `Cal : KGEcomp`) %>%  
                               pivot_wider(names_from = ModelDecisions, values_from = `Cal : KGEcomp`)) %>%
  column_to_rownames(var = "Codes")

best_Comp = names(which.max(apply(Eval_FUSE_Comp, 2, median, na.rm = TRUE)))
best_decision_Comp = as.numeric(gsub("X", "", best_Comp))

Eval_FUSE_best_Comp_Cal <- Eval_FUSE_long_Comp %>%
  filter(ModelDecisions == best_decision_Comp) %>%
  dplyr::select(Codes,`Cal : KGEcomp`) %>%
  rename(`KGEcomp` = `Cal : KGEcomp`)

Eval_FUSE_best_Comp_Eval <- Eval_FUSE_long_Comp %>%
  filter(ModelDecisions == best_decision_Comp) %>%
  dplyr::select(Codes,`Eval : KGEcomp`) %>%
  rename(`KGEcomp` = `Eval : KGEcomp`)




shp_outlets_joined <-shp_outlets %>%
  left_join(
    Eval_FUSE_best_Comp_Cal %>%
      select(Codes, Cal = KGEcomp),
    by = "Codes"
  ) %>%
  left_join(
    Eval_FUSE_best_Comp_Eval %>%
      select(Codes, Eval = KGEcomp),
    by = "Codes"
  ) %>%
  drop_na(Cal, Eval)

###### Plot

gg1 <- ggplot() +
  geom_sf(data = NorthAm, fill = "grey90", color = "grey40") +
  geom_sf(data = shp_outlets_joined, aes(fill = Cal), size = 3, shape = 21, color = "black") +
  scale_fill_gradientn(
    limits = c(-0.42, 1),
    colors = c(rep("purple",2), viridis(80)),
    breaks = seq(-0.4,1,0.1),
    labels = c("< -0.4", seq(-3,10,1)/10),
    oob = scales::squish,
    guide = guide_colorbar(title = expression(KGE[comp]), barwidth = 1.5, barheight = 10)
  ) +
  labs(title = "(a) Calibration period")+
  annotation_scale(location = "bl", style = "bar") +
  annotation_north_arrow(location = "br", style = north_arrow_fancy_orienteering, height = unit(1.1, "cm"), width = unit(1.1, "cm"))+
  theme_bw()

gg2 <- ggplot() +
  geom_sf(data = NorthAm, fill = "grey90", color = "grey40") +
  geom_sf(data = shp_outlets_joined, aes(fill = Eval), size = 3, shape = 21, color = "black") +
  scale_fill_gradientn(
    limits = c(-0.42, 1),
    colors = c(rep("purple",2), viridis(80)),
    breaks = seq(-0.4,1,0.1),
    labels = c("< -0.4", seq(-3,10,1)/10),
    oob = scales::squish,
    guide = guide_colorbar(expression(KGE[comp]), barwidth = 1.5, barheight = 10)
  ) +
  labs(title = "(b) Evaluation period")+
  annotation_scale(location = "bl", style = "bar") +
  annotation_north_arrow(location = "br", style = north_arrow_fancy_orienteering, height = unit(1.1, "cm"), width = unit(1.1, "cm"))+
  theme_bw()



# Define the layout
layout <- "
A
B
"

# Combine the plots
combined_plot <- gg1 + gg2 +
  plot_layout(design = layout)


ggsave(plot = combined_plot, 
       filename = "99_Figures/f04.png",
       width = 6, height = 6, dpi = 300)

