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
library(mapview)
library(ggplot2)
library(ggspatial)
library(ggpubr)
library(grid)
library(dplyr)
library(tidyr)
library(viridis)
library(patchwork)


source(file = "Metrics.R")

#! ---------------------------- bassins

##-----------------------------------------
##-------------- VARIABLES ----------------
##-----------------------------------------

######## Code catchments

catchments = as.character(unname(unlist(read.table(file.path("00_DATA", "liste_BV_CAMELS_559.txt")))))

######## Shapefiles

## catchments outlets
csv_outlets_DB <- readLines(file.path("00_DATA","gauge_information.txt"))
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

###### Get the evaluation data

dir_FUSE = file.path('00_DATA')

# FUSE comp
Eval_FUSE_long_Comp = loadRData(file.path(dir_FUSE, "KGEcomp", "Eval_FUSE.Rdata"))%>%
  mutate(`Cal : KGEcomp` = (`Cal : KGE :  1` + `Cal : KGE : -1`) / 2,
         `Eval : KGEcomp` = (`Eval : KGE :  1` + `Eval : KGE : -1`) / 2)

decision_medians_Comp <- Eval_FUSE_long_Comp %>%
  group_by(ModelDecisions) %>%
  summarise(median_Cal_KGEcomp = median(`Cal : KGEcomp`, na.rm = TRUE))

best_decision_Comp <- decision_medians_Comp %>%
  filter(median_Cal_KGEcomp == max(median_Cal_KGEcomp)) %>%
  pull(ModelDecisions)

Eval_FUSE_best_Comp <- Eval_FUSE_long_Comp %>%
  filter(ModelDecisions == best_decision_Comp) %>%
  dplyr::select(Codes,`Eval : KGEcomp`) %>%
  rename(`KGEcomp` = `Eval : KGEcomp`) %>%
  mutate(Type = "Benchmark (top-performing single model for everywhere)")


# FUSE Mosa perf
df_eval_gumboot = read.table(file.path(dir_FUSE, "KGEcomp","brief_analysis_modeling_results_with_gumboot_eval.csv"),
                             header = TRUE, sep = ",")
df_eval_gumboot$Codes = ifelse(nchar(df_eval_gumboot$gauge_id) == 7, paste0("USA_0", df_eval_gumboot$gauge_id), paste0("USA_", df_eval_gumboot$gauge_id))
df_eval_gumboot$best_cal_model =gsub(x = df_eval_gumboot$best_cal_model, pattern = "Comp_", replacement = "")

df_lookup <- data.frame(Codes = df_eval_gumboot$Codes, best_cal_model = df_eval_gumboot$best_cal_model)

Eval_FUSE_best_Mosa_P <- Eval_FUSE_long_Comp %>%
  left_join(df_lookup, by = "Codes") %>%
  group_by(Codes) %>%
  filter(if_else(is.na(best_cal_model), row_number() == 1, ModelDecisions == best_cal_model)) %>% 
  ungroup()

Eval_FUSE_best_Mosa_P[is.na(Eval_FUSE_best_Mosa_P$best_cal_model),3:12] = NA

Eval_FUSE_best_Mosa_P <- Eval_FUSE_best_Mosa_P %>%
  dplyr::select(Codes, `Eval : KGEcomp`) %>%
  rename(`KGEcomp` = `Eval : KGEcomp`) %>%
  mutate(Type = "Mosaic based on performance")

# FUSE Mosa perf-equivalence
df_eval_gumboot = read.table(file.path(dir_FUSE, "KGEcomp","brief_analysis_modeling_results_with_gumboot_eval.csv"),
                             header = TRUE, sep = ",")
df_eval_gumboot$Codes = ifelse(nchar(df_eval_gumboot$gauge_id) == 7, paste0("USA_0", df_eval_gumboot$gauge_id), paste0("USA_", df_eval_gumboot$gauge_id))
df_eval_gumboot$mod_need = gsub(x = df_eval_gumboot$mod_need, pattern = "Comp_", replacement = "")

df_lookup <- data.frame(Codes = df_eval_gumboot$Codes, mod_need = df_eval_gumboot$mod_need)

Eval_FUSE_best_Mosa_PE <- Eval_FUSE_long_Comp %>%
  left_join(df_lookup, by = "Codes") %>%
  group_by(Codes) %>%
  filter(if_else(is.na(mod_need), row_number() == 1, ModelDecisions == mod_need)) %>% 
  ungroup()

Eval_FUSE_best_Mosa_PE[is.na(Eval_FUSE_best_Mosa_PE$mod_need),3:12] = NA

Eval_FUSE_best_Mosa_PE <- Eval_FUSE_best_Mosa_PE %>%
  dplyr::select(Codes,`Eval : KGEcomp`) %>%
  rename(`KGEcomp` = `Eval : KGEcomp`) %>%
  mutate(Type = "Mosaic based on performance-equivalence")


# FUSE SA (static in time and space)

Eval_FUSE_long_SA = loadRData(file.path(dir_FUSE,"KGEcomp", "Eval_SA.Rdata"))
Eval_FUSE_long_SA_Cal <- Eval_FUSE_long_SA["KGE", "Comp", "Cal", , ]
Eval_FUSE_long_SA_Eval <- Eval_FUSE_long_SA["KGE", "Comp", "Eval", , ]

decision_medians_SA_TS <- apply(Eval_FUSE_long_SA_Cal, 1, median, na.rm = TRUE)
best_decision_SA_TS <- names(decision_medians_SA_TS)[which.max(decision_medians_SA_TS)]

Eval_FUSE_best_SA_TS <- tibble(
  Codes = dimnames(Eval_FUSE_long_SA)$code,
  KGEcomp = as.numeric(Eval_FUSE_long_SA_Eval[best_decision_SA_TS, ] ),
  Type = "Spatially and temporally static combination"
)


# FUSE SA (static in time only)

# Eval_FUSE_long_SA = loadRData(file.path(dir_FUSE,"KGEcomp", "Eval_SA.Rdata"))
# Eval_FUSE_long_SA_Cal <- Eval_FUSE_long_SA["KGE", "Comp", "Cal", , ]
# Eval_FUSE_long_SA_Eval <- Eval_FUSE_long_SA["KGE", "Comp", "Eval", , ]

best_decision_indices_SA_T <- vapply(
  seq_len(ncol(Eval_FUSE_long_SA_Cal)),
  function(i) {
    col <- Eval_FUSE_long_SA_Cal[, i]
    if (all(is.na(col))) NA_integer_ else which.max(col)
  },
  FUN.VALUE = integer(1)
)
best_decision_SA_T <- dimnames(Eval_FUSE_long_SA)$model[best_decision_indices_SA_T]

Eval_FUSE_best_SA_T <- tibble(
  Codes = dimnames(Eval_FUSE_long_SA)$code,
  KGEcomp = as.numeric(Eval_FUSE_long_SA_Eval[cbind(best_decision_indices_SA_T, seq_along(best_decision_indices_SA_T))]),
  Type = "Spatially variable and temporally static combination"
)

# FUSE WA

Eval_FUSE_long_WA = loadRData(file.path(dir_FUSE,"KGEcomp", "Eval_WA.Rdata"))%>%
  mutate(`Cal : KGEcomp` = (`Cal : KGE :  1` + `Cal : KGE : -1`) / 2,
         `Eval : KGEcomp` = (`Eval : KGE :  1` + `Eval : KGE : -1`) / 2)

best_decision_WA <- Eval_FUSE_long_WA %>%
  filter(!startsWith(as.character(ModelDecisions), "1_")) %>%
  group_by(Codes) %>%
  summarise(
    best_model = {
      vals <- `Cal : KGEcomp`
      if (all(is.na(vals))) NA_character_
      else as.character(ModelDecisions[which.max(vals)])
    },
    .groups = "drop"
  )%>%
  pull(best_model)


Eval_FUSE_best_WA <- tibble(
  Codes = unique(as.character(Eval_FUSE_long_WA$Codes)),
  ModelDecisions = best_decision_WA
) %>%
  left_join(
    Eval_FUSE_long_WA %>%
      select(Codes, ModelDecisions, `Eval : KGEcomp`),
    by = c("Codes", "ModelDecisions")
  ) %>%
  transmute(
    Codes = factor(Codes, levels = unique(as.character(Eval_FUSE_long_WA$Codes))),
    KGEcomp = `Eval : KGEcomp`,
    Type = "Dynamic combination"
  )



# Add the columns
shp_outlets_joined <-shp_outlets %>%
  left_join(
    Eval_FUSE_best_Comp %>%
      select(Codes, Comp = KGEcomp),
    by = "Codes"
  ) %>%
  left_join(
    Eval_FUSE_best_Mosa_P %>%
      select(Codes, Mosa_P = KGEcomp),
    by = "Codes"
  ) %>%
  left_join(
    Eval_FUSE_best_Mosa_PE %>%
      select(Codes, Mosa_PE = KGEcomp),
    by = "Codes"
  ) %>%
  left_join(
    Eval_FUSE_best_SA_TS %>%
      select(Codes, SA_TS = KGEcomp),
    by = "Codes"
  ) %>%
  left_join(
    Eval_FUSE_best_SA_T %>%
      select(Codes, SA_T = KGEcomp),
    by = "Codes"
  ) %>%
  left_join(
    Eval_FUSE_best_WA %>%
      select(Codes, WA = KGEcomp),
    by = "Codes"
  ) %>%
  drop_na(Comp, Mosa_P, Mosa_PE, SA_TS, SA_T, WA)


# Add the differences

shp_outlets_joined <- shp_outlets_joined %>%
  mutate(
    diff_Mosa_P  = Mosa_P - Comp,
    diff_Mosa_PE = Mosa_PE - Comp,
    diff_SA_TS   = SA_TS - Comp,
    diff_SA_T    = SA_T - Comp,
    diff_WA      = WA - Comp
  )

###### Plot

# Mosa_P
gg1_1 <- ggplot() +
  geom_sf(data = NorthAm, fill = "grey90", color = "grey40") +
  geom_sf(data = shp_outlets_joined, aes(fill = Mosa_P), size = 3, shape = 21, color = "black") +
  scale_fill_gradientn(
    limits = c(-0.42, 1),
    colors = c(rep("purple",2), viridis(80)),
    breaks = seq(-0.4,1,0.1),
    labels = c("< -0.4", seq(-3,10,1)/10),
    oob = scales::squish,
    guide = guide_colorbar(title = expression(KGE[comp]), barwidth = 1.5, barheight = 10)
  ) +
  labs(title = "(a) Mosaic based on performance")+
  annotation_scale(location = "bl", style = "bar") +
  annotation_north_arrow(location = "br", style = north_arrow_fancy_orienteering, height = unit(1.1, "cm"), width = unit(1.1, "cm"))+
  theme_bw()

gg1_2 <- ggplot() +
  geom_sf(data = NorthAm, fill = "grey90", color = "grey40") +
  geom_sf(
    data = shp_outlets_joined,
    aes(fill = diff_Mosa_P),
    size = 3, shape = 21, color = "black"
  ) +
  scale_fill_gradient2(
    low = "red",      # For negative values
    mid = "white",    # For zero
    high = "blue",    # For positive values
    midpoint = 0,
    limits = c(-0.5, 0.5),
    oob = scales::squish,  # Values beyond limits will be squished to end colors
    name = expression(ΔKGE[comp]),     # Legend title (you can customize it)
    breaks = c(-0.5, 0, 0.5),
    labels = c("< -0.5", "0", "> 0.5")
  ) +
  labs(title = "Difference with the benchmark") +
  annotation_scale(location = "bl", style = "bar") +
  annotation_north_arrow(
    location = "br",
    style = north_arrow_fancy_orienteering,
    height = unit(1.1, "cm"),
    width = unit(1.1, "cm")
  ) +
  theme_bw()

# Mosa_PE
gg2_1 <- ggplot() +
  geom_sf(data = NorthAm, fill = "grey90", color = "grey40") +
  geom_sf(data = shp_outlets_joined, aes(fill = Mosa_PE), size = 3, shape = 21, color = "black") +
  scale_fill_gradientn(
    limits = c(-0.42, 1),
    colors = c(rep("purple",2), viridis(80)),
    breaks = seq(-0.4,1,0.1),
    labels = c("< -0.4", seq(-3,10,1)/10),
    oob = scales::squish,
    guide = guide_colorbar(title = expression(KGE[comp]), barwidth = 1.5, barheight = 10)
  ) +
  labs(title = "(b) Mosaic based on performance-equivalence")+
  annotation_scale(location = "bl", style = "bar") +
  annotation_north_arrow(location = "br", style = north_arrow_fancy_orienteering, height = unit(1.1, "cm"), width = unit(1.1, "cm"))+
  theme_bw()

gg2_2 <- ggplot() +
  geom_sf(data = NorthAm, fill = "grey90", color = "grey40") +
  geom_sf(
    data = shp_outlets_joined,
    aes(fill = diff_Mosa_PE),
    size = 3, shape = 21, color = "black"
  ) +
  scale_fill_gradient2(
    low = "red",      # For negative values
    mid = "white",    # For zero
    high = "blue",    # For positive values
    midpoint = 0,
    limits = c(-0.5, 0.5),
    oob = scales::squish,  # Values beyond limits will be squished to end colors
    name = expression(ΔKGE[comp]),     # Legend title (you can customize it)
    breaks = c(-0.5, 0, 0.5),
    labels = c("< -0.5", "0", "> 0.5")
  ) +
  labs(title = "Difference with the benchmark") +
  annotation_scale(location = "bl", style = "bar") +
  annotation_north_arrow(
    location = "br",
    style = north_arrow_fancy_orienteering,
    height = unit(1.1, "cm"),
    width = unit(1.1, "cm")
  ) +
  theme_bw()

# SA_TS
gg3_1 <- ggplot() +
  geom_sf(data = NorthAm, fill = "grey90", color = "grey40") +
  geom_sf(data = shp_outlets_joined, aes(fill = SA_TS), size = 3, shape = 21, color = "black") +
  scale_fill_gradientn(
    limits = c(-0.42, 1),
    colors = c(rep("purple",2), viridis(80)),
    breaks = seq(-0.4,1,0.1),
    labels = c("< -0.4", seq(-3,10,1)/10),
    oob = scales::squish,
    guide = guide_colorbar(title = expression(KGE[comp]), barwidth = 1.5, barheight = 10)
  ) +
  labs(title = "(c) Spatially and temporally static combination")+
  annotation_scale(location = "bl", style = "bar") +
  annotation_north_arrow(location = "br", style = north_arrow_fancy_orienteering, height = unit(1.1, "cm"), width = unit(1.1, "cm"))+
  theme_bw()

gg3_2 <- ggplot() +
  geom_sf(data = NorthAm, fill = "grey90", color = "grey40") +
  geom_sf(
    data = shp_outlets_joined,
    aes(fill = diff_SA_TS),
    size = 3, shape = 21, color = "black"
  ) +
  scale_fill_gradient2(
    low = "red",      # For negative values
    mid = "white",    # For zero
    high = "blue",    # For positive values
    midpoint = 0,
    limits = c(-0.5, 0.5),
    oob = scales::squish,  # Values beyond limits will be squished to end colors
    name = expression(ΔKGE[comp]),     # Legend title (you can customize it)
    breaks = c(-0.5, 0, 0.5),
    labels = c("< -0.5", "0", "> 0.5")
  ) +
  labs(title = "Difference with the benchmark") +
  annotation_scale(location = "bl", style = "bar") +
  annotation_north_arrow(
    location = "br",
    style = north_arrow_fancy_orienteering,
    height = unit(1.1, "cm"),
    width = unit(1.1, "cm")
  ) +
  theme_bw()

# SA_T
gg4_1 <- ggplot() +
  geom_sf(data = NorthAm, fill = "grey90", color = "grey40") +
  geom_sf(data = shp_outlets_joined, aes(fill = SA_T), size = 3, shape = 21, color = "black") +
  scale_fill_gradientn(
    limits = c(-0.42, 1),
    colors = c(rep("purple",2), viridis(80)),
    breaks = seq(-0.4,1,0.1),
    labels = c("< -0.4", seq(-3,10,1)/10),
    oob = scales::squish,
    guide = guide_colorbar(title = expression(KGE[comp]), barwidth = 1.5, barheight = 10)
  ) +
  labs(title = "(d) Spatially variable and temporally static combination")+
  annotation_scale(location = "bl", style = "bar") +
  annotation_north_arrow(location = "br", style = north_arrow_fancy_orienteering, height = unit(1.1, "cm"), width = unit(1.1, "cm"))+
  theme_bw()

gg4_2 <- ggplot() +
  geom_sf(data = NorthAm, fill = "grey90", color = "grey40") +
  geom_sf(
    data = shp_outlets_joined,
    aes(fill = diff_SA_T),
    size = 3, shape = 21, color = "black"
  ) +
  scale_fill_gradient2(
    low = "red",      # For negative values
    mid = "white",    # For zero
    high = "blue",    # For positive values
    midpoint = 0,
    limits = c(-0.5, 0.5),
    oob = scales::squish,  # Values beyond limits will be squished to end colors
    name = expression(ΔKGE[comp]),     # Legend title (you can customize it)
    breaks = c(-0.5, 0, 0.5),
    labels = c("< -0.5", "0", "> 0.5")
  ) +
  labs(title = "Difference with the benchmark") +
  annotation_scale(location = "bl", style = "bar") +
  annotation_north_arrow(
    location = "br",
    style = north_arrow_fancy_orienteering,
    height = unit(1.1, "cm"),
    width = unit(1.1, "cm")
  ) +
  theme_bw()

# WA
gg5_1 <- ggplot() +
  geom_sf(data = NorthAm, fill = "grey90", color = "grey40") +
  geom_sf(data = shp_outlets_joined, aes(fill = WA), size = 3, shape = 21, color = "black") +
  scale_fill_gradientn(
    limits = c(-0.42, 1),
    colors = c(rep("purple",2), viridis(80)),
    breaks = seq(-0.4,1,0.1),
    labels = c("< -0.4", seq(-3,10,1)/10),
    oob = scales::squish,
    guide = guide_colorbar(title = expression(KGE[comp]), barwidth = 1.5, barheight = 10)
  ) +
  labs(title = "(e) Dynamic combination")+
  annotation_scale(location = "bl", style = "bar") +
  annotation_north_arrow(location = "br", style = north_arrow_fancy_orienteering, height = unit(1.1, "cm"), width = unit(1.1, "cm"))+
  theme_bw()

gg5_2 <- ggplot() +
  geom_sf(data = NorthAm, fill = "grey90", color = "grey40") +
  geom_sf(
    data = shp_outlets_joined,
    aes(fill = diff_WA),
    size = 3, shape = 21, color = "black"
  ) +
  scale_fill_gradient2(
    low = "red",      # For negative values
    mid = "white",    # For zero
    high = "blue",    # For positive values
    midpoint = 0,
    limits = c(-0.5, 0.5),
    oob = scales::squish,  # Values beyond limits will be squished to end colors
    name = expression(ΔKGE[comp]),     # Legend title (you can customize it)
    breaks = c(-0.5, 0, 0.5),
    labels = c("< -0.5", "0", "> 0.5")
  ) +
  labs(title = "Difference with the benchmark") +
  annotation_scale(location = "bl", style = "bar") +
  annotation_north_arrow(
    location = "br",
    style = north_arrow_fancy_orienteering,
    height = unit(1.1, "cm"),
    width = unit(1.1, "cm")
  ) +
  theme_bw()

# Define the layout
layout <- "
AB
CD
EF
GH
IJ
"

# Combine the plots
combined_plot <- gg1_1 + gg1_2 + 
  gg2_1 + gg2_2 + 
  gg3_1 + gg3_2 + 
  gg4_1 + gg4_2 + 
  gg5_1 + gg5_2 + 
  plot_layout(design = layout)


ggsave(plot = combined_plot, filename = "99_Figures/f06.png",
       width = 12, height = 15, dpi = 300)




df <- sf::st_drop_geometry(shp_outlets_joined)

diff_cols <- grep("^diff_", names(df), value = TRUE)

percent_summary <- function(x, tol = 0) {
  n <- sum(!is.na(x))
  gt <- sum(x >  tol, na.rm = TRUE)
  eq <- sum(abs(x) <= tol, na.rm = TRUE)   # use tol if you want a near-zero band
  lt <- sum(x < -tol, na.rm = TRUE)
  c(
    n          = n,
    greater_n  = gt, greater_pct = 100 * gt / n,
    equal_n    = eq, equal_pct   = 100 * eq / n,
    less_n     = lt, less_pct    = 100 * lt / n
  )
}


diff_stats <- t(sapply(df[diff_cols], percent_summary))
diff_stats <- as.data.frame(diff_stats)
diff_stats <- round(diff_stats, 2)
print(diff_stats)



