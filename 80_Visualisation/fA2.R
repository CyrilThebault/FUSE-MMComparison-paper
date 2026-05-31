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

library(tidyverse)
library(sf)
library(ggplot2)
library(ggpubr)
library(viridis)
library(RColorBrewer)

generate_random_colors <- function(num_colors, seed = 1997) {
  set.seed(seed)
  colors <- sample(colors(), num_colors)
  return(colors)
}

make_alternating_labels <- function(labels) {
  if (!is.character(labels)) labels <- as.character(labels)
  ifelse(seq_along(labels) %% 2 == 0,
         paste0("\n ", labels),   # even indices: shifted slightly down
         paste0(labels,"\n"))    # odd indices: shifted slightly up
}

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



border_path <- file.path("Shp","boundaries_p_2021_v3.shp")

# Load and process CONUS shapefile
gdf_borders <- st_read(border_path)
conus <- gdf_borders %>%
  filter(COUNTRY == "USA") %>%
  filter(!STATEABB %in% c("US-AK", "US-HI", "US-PR", "US-VI"))
conus <- st_transform(conus, 4326)

df_eval_gumboot = read.table("00_DATA/KGEcomp/brief_analysis_modeling_results_with_gumboot_eval.csv",
               header = TRUE, sep = ",")

df_eval_gumboot$Codes = ifelse(nchar(df_eval_gumboot$gauge_id) == 7, paste0("USA_0", df_eval_gumboot$gauge_id), paste0("USA_", df_eval_gumboot$gauge_id))

df_eval_gumboot = df_eval_gumboot[!df_eval_gumboot$Codes %in% catchments_fail,]

# Create spatial data
gdf_sub <- st_as_sf(df_eval_gumboot, coords = c("gauge_lon", "gauge_lat"), crs = 4326)
gdf_sub <- st_transform(gdf_sub, st_crs(conus))

# Remove basins without uncertainty info
gdf_gumboot <- gdf_sub %>% filter(score > -990) 

list_names = gsub("_above_p5","",grep("_above_p5$", names(df_eval_gumboot), value = TRUE))

# Create a color dictionary
color_dict <- setNames(generate_random_colors(length(unique(list_names))),
                       unique(list_names))

df <- gdf_gumboot %>%
  st_drop_geometry() %>%
  count(best_cal_model, name = "n") %>%
  mutate(num = as.numeric(sub("Comp_", "", best_cal_model))) %>%
  arrange(desc(n), num)

models_best <- df$best_cal_model
models_best_num <- df$n
models_best_lab <- df$num

####### MOSAIC PERF

# Loop over the models we identified and create incremental plots
ix =  length(models_best) 
model <- models_best[ix]

# Define bar plotting stats
bar_x <- models_best # always the same but kept here for clarity
# bar_y <- cumsum(models_best_num[1:(ix)])
# bar_y <- c(bar_y, rep(-1, length(bar_x) - length(bar_y)))
bar_y <- as.numeric(models_best_num)
  
# Initialize plot

gdf_plot = gdf_gumboot[gdf_gumboot$mod_need %in%  head(models_best,ix),]

# Make the actual figure

# axis_labels <- ifelse(
#   models_best_lab %in% models_needed_lab,
#   paste0("<span style='color:red;'>", models_best_lab, "</span>"),
#   models_best_lab
# )

p1 <- ggplot() +
  geom_sf(data = conus, fill = "grey60", color = "white", size = 0.5) +
  geom_sf(data = gdf_plot, aes(fill = factor(best_cal_model, levels = models_best)), size = 3, shape = 21, color = "black") +
  scale_fill_manual(values = color_dict,
                    guide = 'none') +
  labs(x = "", y = "", title = "a) Mosaic based on performance")+
  theme_bw()


p2 <- ggplot(data.frame(x = bar_x, y = bar_y), aes(x = factor(x, levels = models_best), y = y, fill = factor(x, levels = models_best))) +
  geom_bar(stat = "identity",  col = "black") +
  scale_fill_manual(values = color_dict,
                    guide = 'none') +
  #scale_x_discrete(labels = models_best_lab) +
  # scale_x_discrete(labels = axis_labels) +
  scale_x_discrete(labels = make_alternating_labels(models_best_lab)) +
  labs(title = paste0("b) Number of models selected: ", length(bar_y)), x = "Model structure ID", y = "Number of catchments") +
  theme_bw()
  #theme(axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1, size = 8)) 
  # ylim(0, nrow(gdf_gumboot))


plot1 <- ggarrange(
  p1, p2,
  ncol = 2, nrow = 1,
  widths = c(1, 1.4)  # Set width ratio
)

# Save the plot
ggsave(filename = "99_Figures/fA2.png",
       plot = plot1, 
       width = 14, height = 3.25, dpi = 300, units = "in",
       bg = "white")

