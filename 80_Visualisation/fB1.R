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

library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(purrr)
library(ggnewscale)
library(ggpubr)
library(sf)

source(file = "Metrics.R")

##-----------------------------------------
##---------------- MAIN ------------------
##-----------------------------------------

dir_FUSE = file.path("00_DATA")

# FUSE comp
Eval_FUSE_long_Comp = loadRData(file.path(dir_FUSE,"KGEcomp", "Eval_FUSE.Rdata"))%>%
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
  mutate(Type = "Benchmark")


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
Eval_FUSE_long_WA = loadRData(file.path(dir_FUSE, "KGEcomp", "Eval_WA.Rdata"))%>%
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



combined_df <-rbind(Eval_FUSE_best_Comp,
                    Eval_FUSE_best_Mosa_P, Eval_FUSE_best_Mosa_PE,
                    Eval_FUSE_best_SA_TS, Eval_FUSE_best_SA_T, 
                    Eval_FUSE_best_WA
)

combined_df$Failed <- is.na(combined_df$KGEcomp)




####
border_path <- "Shp/boundaries_p_2021_v3.shp"

# Load and process CONUS shapefile
gdf_borders <- st_read(border_path)
conus <- gdf_borders %>%
  filter(COUNTRY == "USA") %>%
  filter(!STATEABB %in% c("US-AK", "US-HI", "US-PR", "US-VI"))
conus <- st_transform(conus, 4326)

df_eval_gumboot = read.table(file.path(dir_FUSE, "KGEcomp","brief_analysis_modeling_results_with_gumboot_eval.csv"),
                             header = TRUE, sep = ",")

# Create spatial data
gdf_sub <- st_as_sf(df_eval_gumboot, coords = c("gauge_lon", "gauge_lat"), crs = 4326)
gdf_sub <- st_transform(gdf_sub, st_crs(conus))
gdf_sub$Codes = ifelse(nchar(gdf_sub$gauge_id)==7, paste0("USA_0", gdf_sub$gauge_id), paste0("USA_", gdf_sub$gauge_id))


color_dict <- c(
  "Benchmark" = "darkorange",
  "Mosaic based on performance" = "lightpink1",
  "Mosaic based on performance-equivalence" = "lightpink4",
  "Spatially and temporally static combination" = "olivedrab2",
  "Spatially variable and temporally static combination" = "olivedrab4",
  "Dynamic combination" = "firebrick4"
)


gdf_fail <- gdf_sub %>%
  left_join(combined_df, by = "Codes")

gdf_fail$Type <- factor(gdf_fail$Type, levels = names(color_dict))



gg_fail <- ggplot() +
  geom_sf(data = conus, fill = "grey60", color = "white") +
  geom_sf(
    data = gdf_fail %>% filter(Failed),
    aes(fill = Type), color = "black",
    size = 3, shape = 21
  ) +
  scale_fill_manual(values = color_dict) +
  facet_wrap(~ Type) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 9)
  )

gg_fail

ggsave(filename = "99_Figures/fB1.png",
       plot = gg_fail, 
       width = 10, height = 4, dpi = 300, units = "in",
       bg = "white")
