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
library(forcats)
library(ggplot2)
library(tibble)
library(purrr)
library(ggnewscale)
library(ggpubr)
library(scales)
library(viridis)
library(stringr)



source(file = "Metrics.R")

generate_random_colors <- function(num_colors, seed = 1997) {
  set.seed(seed)
  colors <- sample(colors(), num_colors)
  return(colors)
}

collapse_models <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA_character_)
  x <- sort(as.integer(x))        
  paste(x, collapse = "_")       
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


##-----------------------------------------
##---------------- PARAM ------------------
##-----------------------------------------

dir_FUSE = file.path("00_DATA")

# FUSE WA

Eval_FUSE_long_WA = loadRData(file.path(dir_FUSE,"KGEcomp", "Eval_WA.Rdata"))%>%
  mutate(`Cal : KGEcomp` = (`Cal : KGE :  1` + `Cal : KGE : -1`) / 2,
         `Eval : KGEcomp` = (`Eval : KGE :  1` + `Eval : KGE : -1`) / 2)%>%
  filter(!Codes %in% catchments_fail)

df_param_WA <- Eval_FUSE_long_WA %>%
  filter(!startsWith(as.character(ModelDecisions), "1_")) %>%
  group_by(Codes) %>%
  summarise(
    best_cal_param = {
      vals <- `Cal : KGEcomp`
      if (all(is.na(vals))) NA_character_
      else as.character(ModelDecisions[which.max(vals)])
    },
    .groups = "drop"
  )



# 1) Split parameters and pivot longer
param_long <- df_param_WA %>%
  separate(best_cal_param, into = c("τ [day]", "m [-]", "k [-]"), sep = "_", convert = TRUE) %>%
  pivot_longer(c(`τ [day]`, `m [-]`, `k [-]`), names_to = "parameter", values_to = "value") %>%
  filter(!is.na(value))

# 2) Compute counts per parameter
param_counts <- param_long %>%
  count(parameter, value)

# 3) Define the full range of possible parameter values
# (You can adjust these based on the known model parameter grids)
full_ranges <- list(
  `τ [day]` = seq(4, 30, 3),
  `m [-]`  = seq(1, 20, 2),
  `k [-]`   = seq(1, 20, 2)
)

# 4) Create a complete data frame including zeros for missing combinations
param_complete <- bind_rows(
  lapply(names(full_ranges), function(p) {
    tibble(parameter = p, value = full_ranges[[p]])
  })
) %>%
  left_join(param_counts, by = c("parameter", "value")) %>%
  mutate(n = ifelse(is.na(n), 0, n))  # Replace missing counts with 0

# 5) Make 'value' a factor with numeric ordering
param_complete <- param_complete %>%
  mutate(value = factor(value, levels = sort(unique(value))))

# 6) Plot
p_params <- ggplot(param_complete, aes(x = value, y = n)) +
  geom_col(col = "black", linewidth = 0.3, fill = "skyblue") +
  facet_wrap(~ parameter, ncol = 3, scales = "free_x") +
  labs(
    x = "Parameter value",
    y = "Number of catchments"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )


p_params


ggsave(plot = p_params, filename = "99_Figures/fA8.png",
       width = 10, height = 5, dpi = 300)




##-----------------------------------------
##------------- MOD INDIV -----------------
##-----------------------------------------

dir_FUSE = file.path("00_DATA")

# FUSE WA

SelecMod_df = loadRData(file.path(dir_FUSE,"KGEcomp", "selected_models_WA.Rdata"))


##----------------- Overall -------------------

# Flatten all model values into a single vector
all_models <- unlist(
  unlist(SelecMod_df, recursive = TRUE, use.names = FALSE),
  recursive = TRUE, use.names = FALSE
)

all_models <- all_models[!is.na(all_models) & all_models != ""]

# Count occurrences and calculate usage percentage
model_usage <- data.frame(model = all_models) %>%
  count(model, name = "n") %>%
  mutate(pct = n / sum(n)) %>%
  arrange(desc(pct)) %>%
  mutate(model = factor(model, levels = model))

color_dict <- setNames(generate_random_colors(nrow(model_usage)),
                       unique(sort(as.character(model_usage$model))))

# Plot all models (sorted descending)
gg <- ggplot(model_usage, aes(x = model, y = pct, fill = model)) +
  geom_col(col = "black") +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_fill_manual(values = color_dict,
                    guide = 'none') +
  labs(
    x = "Model structure ID (individual model within the combinations)",
    y = "Usage Percentage"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

gg

ggsave(plot = gg, filename = "99_Figures/fA9.png",
       width = 14, height = 3.25, dpi = 300)

##----------------- Per catchment analysis -------------------


#Percentage of model use per catchment
catchment_usage <- SelecMod_df %>%
  pivot_longer(everything(), names_to = "catchment", values_to = "models") %>%
  unnest(models) %>%
  filter(models != "" & !is.na(models)) %>%
  count(catchment, models, name = "n") %>%
  group_by(catchment) %>%
  mutate(pct = n / sum(n)) %>%
  ungroup()%>%
  mutate(models = factor(models, levels = sort(as.numeric(unique(models)))))

# Heatmap
hm <- ggplot(catchment_usage, aes(x = models, y = catchment, fill = pct)) +
  geom_tile() +
  scale_fill_gradientn(
    limits = c(0, 10)/100,
    colors = c(viridis(80), rep("red",2)),
    breaks = seq(0,10,1)/100,
    labels = c(paste0(seq(0,9,1), "%"),"> 10%"),
    oob = scales::squish,
    guide = guide_colorbar(title = "Usage Percentage", barwidth = 1.5, barheight = 10)
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    panel.grid = element_blank()) +
  labs(x = "Model Structure ID", y = "Catchment ID")


hm

ggsave(plot = hm, filename = "99_Figures/fA10.png",
       width = 11, height = 7, dpi = 300)


# ##-----------------------------------------
# ##------------- MOD COMBI -----------------
# ##-----------------------------------------
# 
# dir_FUSE = "00_DATA"
# 
# # FUSE WA
# 
# SelecMod_df = loadRData(file.path(dir_FUSE,"KGEcomp", "selected_models_WA.Rdata"))
# SelecMod_combo <- SelecMod_df %>%
#   mutate(across(everything(), ~ map_chr(.x, collapse_models)))
# 
# 
# ##----------------- Overall -------------------
# 
# # Flatten all model values into a single vector
# all_models <- unlist(SelecMod_combo)
# 
# all_models <- all_models[!is.na(all_models) & all_models != ""]
# 
# # Count occurrences and calculate usage percentage
# model_usage <- data.frame(model = all_models) %>%
#   count(model, name = "n") %>%
#   mutate(pct = n / sum(n)) %>%
#   arrange(desc(pct)) %>%
#   mutate(model = factor(model, levels = model))
# 
# 
# model_usage <- model_usage %>%
#   mutate(
#     n_models = str_count(as.character(model), "_") + 1
#   )
# 
# # Check the distribution
# table(model_usage$n_models)
# 
# # Plot top models faceted by number of models in combination
# top_models <- model_usage %>%
#   arrange(desc(pct)) %>%
#   group_by(n_models) %>%
#   slice_head(n = 20) %>%  # Top 20 per group
#   ungroup()
# 
# gg <- ggplot(top_models, aes(x = reorder(model, -pct), y = pct)) +
#   geom_linerange(aes(ymin = 0.0001/100, ymax = pct),
#                  color = "skyblue", linewidth = 2) +
#   scale_y_log10(
#     labels = percent_format(accuracy = 0.001),
#     breaks = c(0.0001, 0.001, 0.01, 0.1, 1)/100,
#     limits = c(0.0001/100, NA),
#     expand = expansion(mult = c(0, 0.05))
#   ) +
#   facet_wrap(~ n_models, scales = "free_x",
#              labeller = labeller(n_models = function(x) paste(x, "model(s)"))) +
#   labs(
#     x = "Model structure ID",
#     y = "Usage Percentage (log scale)",
#     title = "Top Models by Combination Size"
#   ) +
#   theme_bw(base_size = 12) +
#   theme(
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank()
#   )
# 
# gg
# 
# 
# ##----------------- Per catchment analysis -------------------
# 
# #Percentage of model use per catchment
# catchment_usage <- SelecMod_combo %>%
#   pivot_longer(everything(), names_to = "catchment", values_to = "models") %>%
#   unnest(models) %>%
#   filter(models != "" & !is.na(models)) %>%
#   count(catchment, models, name = "n") %>%
#   group_by(catchment) %>%
#   mutate(pct = n / sum(n)) %>%
#   ungroup()%>%
#   arrange(desc(pct))
# 
# catchment_usage <- catchment_usage %>%
#   mutate(
#     n_models = str_count(as.character(models), "_") + 1
#   )
# 
# # Heatmap
# hm <- ggplot(catchment_usage, aes(x = models, y = catchment, fill = pct)) +
#   geom_tile() +
#   scale_fill_gradientn(
#     limits = c(0, 10)/100,
#     colors = c(viridis(80), rep("red",2)),
#     breaks = seq(0,10,1)/100,
#     labels = c(paste0(seq(0,9,1), "%"),"> 10%"),
#     oob = scales::squish,
#     guide = guide_colorbar(title = "Usage Percentage", barwidth = 1.5, barheight = 10)
#   ) +
#   facet_wrap(~ n_models, scales = "free",
#              labeller = labeller(n_models = function(x) paste(x, "model(s)"))) +
#   theme_bw() +
#   theme(
#     axis.text = element_blank(),
#     axis.ticks = element_blank(),
#     panel.grid = element_blank()) +
#   labs(x = "Model Structure ID", y = "Catchment ID")
# 
# 
# hm

