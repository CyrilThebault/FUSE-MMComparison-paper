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
library(ggpubr)
library(tibble)
library(purrr)

source(file = "Metrics.R")

##-----------------------------------------
##-------------- PARENTS MOD --------------
##-----------------------------------------

# TOPMODEL: 64;multiplc_e;tension2_1;unlimpow_2;tmdl_param;perc_f2sat;sequential;intflwnone;rout_gamma;temp_index
# VIC: 96;multiplc_e;onestate_1;fixedsiz_2;arno_x_vic;perc_w2sat;sequential;intflwnone;rout_gamma;temp_index
# PRMS: 34;multiplc_e;tension2_1;unlimfrc_2;prms_varnt;perc_f2sat;sequential;intflwnone;rout_gamma;temp_index
# SACRAMENTO: 170;multiplc_e;tension1_1;tens2pll_2;prms_varnt;perc_lower;sequential;intflwnone;rout_gamma;temp_index
parents <- c("PRMS" = "X34", 
             "TOPMODEL" = "X64", 
             "VIC" = "X96", 
             "SACRAMENTO" = "X170")

##-----------------------------------------
##---------------- BEST MOD ---------------
##-----------------------------------------

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

best_models <- c(
  "Composite ensemble: KGEcomp"   = best_Comp
)


##-----------------------------------------
##------------ Calibration ----------------
##-----------------------------------------

# # FUSE Comp
# dir_FUSE_Comp = file.path(path_res, 'FUSE')
# Eval_FUSE_long_Comp = loadRData(file.path(dir_FUSE_Comp,"01_Paper", "Comp", "Eval_FUSE.Rdata"))%>%
#   mutate(`Cal : KGEcomp` = (`Cal : KGE :  1` + `Cal : KGE : -1`) / 2,
#          `Eval : KGEcomp` = (`Eval : KGE :  1` + `Eval : KGE : -1`) / 2)
# 
# Eval_FUSE_Comp <- data.frame(Eval_FUSE_long_Comp %>%
#                                dplyr::select(Codes, ModelDecisions, `Cal : KGEcomp`) %>%  # dplyr::select only the relevant columns
#                                pivot_wider(names_from = ModelDecisions, values_from = `Cal : KGEcomp`)) %>%
#   column_to_rownames(var = "Codes")



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

Eval_FUSE_Comp   <- Eval_FUSE_Comp[!rownames(Eval_FUSE_Comp) %in% catchments_fail, ]



# Ensure the rownames are preserved as a column for merging
Eval_FUSE_Comp$rownames <- rownames(Eval_FUSE_Comp)


data_frames <- list(Eval_FUSE_Comp)


# Merge the data frames by rownames
combined_df <- reduce(data_frames, full_join, by = "rownames")

colnames(combined_df)[colnames(combined_df) == "rownames"] = "Catchment"

# Pivot the data frame into a long format
long_df <- combined_df %>%
  pivot_longer(
    cols = -Catchment,            # Pivot all columns except 'Catchment'
    names_to = "Model",           # New column for the model names
    values_to = "KGE_Cal"        # New column for the model values
  )

long_df$KGE_Cal[is.na(long_df$KGE_Cal)] = -9999

# Create a CDF dataframe for each ModelDecision
cdf_combined <- long_df %>%
  group_by(Model) %>%
  arrange(as.numeric(KGE_Cal)) %>%
  mutate(CDF = seq(0, 1, length.out = n()),
         Model = as.character(Model)) %>%
  ungroup()



# Modify the color_dict and add a new column for legend grouping

cdf_combined <- cdf_combined %>%
  mutate(IsBest = ifelse(Model %in% best_models, TRUE, FALSE))%>%
  mutate(IsParent = ifelse(Model %in% parents, TRUE, FALSE))


color_dict <- c(
  "FUSE ensemble" = "darkorange",
  "Top-performing single model" = "darkorange",
  "Parents models" = "gray50"
)

targets = c(0.125, 0.25, 0.375, 0.5, 0.625 ,0.75, 0.875)

find_closest_cdf <- function(cdf_values, targets) {
  sapply(targets, function(target) {
    which.min(abs(cdf_values - target))  # Index of the closest CDF value
  })
}

CDF_values = unique(cdf_combined$CDF)

index = find_closest_cdf(CDF_values, targets)

approx_target = CDF_values[index]

cdf_boxplot_data <- cdf_combined %>%
  filter(CDF %in% approx_target)

# Create a small offset for the boxplots to avoid overlap (shift by LegendGroup)
cdf_boxplot_data <- cdf_boxplot_data %>%
  mutate(CDF_offset = case_when(
    TRUE ~ CDF
  ))

cdf_boxplot_data$LegendGroup = TRUE

# Create the plot
gg1 <- ggplot(cdf_combined, aes(x = KGE_Cal, y = CDF, group = Model)) +
  # FUSE ensemble CDF
  geom_line(data = subset(cdf_combined, !IsBest), 
            color = color_dict["FUSE ensemble"], linewidth = 0.4, alpha = 0.2) +
  
  # Parents models CDF
  geom_line(data = subset(cdf_combined, IsParent), color = color_dict["Parents models"], linewidth = 0.5, alpha = 1) +
  
  # Best models CDF
  geom_line(data = subset(cdf_combined, IsBest), color = "black", linewidth = 1.8, alpha = 1) +
  geom_line(data = subset(cdf_combined, IsBest), 
            color = color_dict["Top-performing single model"], linewidth = 1.5, alpha = 1) +
  
  # Add vertical reference line
  geom_vline(xintercept = 1, col = "deepskyblue4", size = 0.8, linetype = "dashed") +
  
  # Add boxplots at CDF = 0.25, 0.5, and 0.75 for each LegendGroup
  geom_boxplot(data = cdf_boxplot_data, aes(x = KGE_Cal, y = CDF_offset, group = interaction(CDF, LegendGroup)),
               fill = color_dict["FUSE ensemble"],
               width = 0.05, alpha = 0.5, outlier.shape = NA) +  # `width` controls horizontal width
  
  # Labels and theme
  labs(x = expression(KGE[comp]), y = "Cumulative distribution function (CDF)", title = "(a) Calibration period") +
  theme_bw() +
  coord_cartesian(xlim = c(-0.41, 1.02), ylim = c(-0.02, 1.02), expand = FALSE) +
  
  scale_x_continuous(
    breaks = seq(-0.4, 1.0, by = 0.2),   # major tick marks every 0.2
    minor_breaks = seq(-0.4, 1.0, by = 0.1)  # faint minor grid every 0.1
  )+
  
  # Remove legend for cleaner visualization
  theme(legend.position = "none")   

gg1



medians <- apply(Eval_FUSE_Comp[ , setdiff(names(Eval_FUSE_Comp), "rownames")],
                 2, median, na.rm = TRUE)
min_median <- min(medians)
max_median <- max(medians)
min_col <- names(medians)[which.min(medians)]
max_col <- names(medians)[which.max(medians)]
cat("Minimum median:", min_median, " (column:", min_col, ")\n")
cat("Maximum median:", max_median, " (column:", max_col, ")\n")

##-----------------------------------------
##------------ Evaluation -----------------
##-----------------------------------------

# FUSE Comp
Eval_FUSE_long_Comp = loadRData(file.path(dir_FUSE, "KGEcomp", "Eval_FUSE.Rdata"))%>%
  mutate(`Cal : KGEcomp` = (`Cal : KGE :  1` + `Cal : KGE : -1`) / 2,
         `Eval : KGEcomp` = (`Eval : KGE :  1` + `Eval : KGE : -1`) / 2)

Eval_FUSE_Comp <- data.frame(Eval_FUSE_long_Comp %>%
                               dplyr::select(Codes, ModelDecisions, `Eval : KGEcomp`) %>%  # dplyr::select only the relevant columns
                               pivot_wider(names_from = ModelDecisions, values_from = `Eval : KGEcomp`)) %>%
  column_to_rownames(var = "Codes")



# Ensure the rownames are preserved as a column for merging
Eval_FUSE_Comp$rownames <- rownames(Eval_FUSE_Comp)


data_frames <- list(Eval_FUSE_Comp)


# Merge the data frames by rownames
combined_df <- reduce(data_frames, full_join, by = "rownames")

colnames(combined_df)[colnames(combined_df) == "rownames"] = "Catchment"

# Pivot the data frame into a long format
long_df <- combined_df %>%
  pivot_longer(
    cols = -Catchment,            # Pivot all columns except 'Catchment'
    names_to = "Model",           # New column for the model names
    values_to = "KGE_Eval"        # New column for the model values
  )

long_df$KGE_Eval[is.na(long_df$KGE_Eval)] = -9999

# Create a CDF dataframe for each ModelDecision
cdf_combined <- long_df %>%
  group_by(Model) %>%
  arrange(as.numeric(KGE_Eval)) %>%
  mutate(CDF = seq(0, 1, length.out = n()),
         Model = as.character(Model)) %>%
  ungroup()



cdf_combined <- cdf_combined %>%
  mutate(IsBest = ifelse(Model %in% best_models, TRUE, FALSE))%>%
  mutate(IsParent = ifelse(Model %in% parents, TRUE, FALSE))


color_dict <- c(
  "FUSE ensemble" = "darkorange",
  "Top-performing single model" = "darkorange",
  "Parents models" = "gray50"
)

targets = c(0.125, 0.25, 0.375, 0.5, 0.625 ,0.75, 0.875)

find_closest_cdf <- function(cdf_values, targets) {
  sapply(targets, function(target) {
    which.min(abs(cdf_values - target))  # Index of the closest CDF value
  })
}

CDF_values = unique(cdf_combined$CDF)

index = find_closest_cdf(CDF_values, targets)

approx_target = CDF_values[index]

cdf_boxplot_data <- cdf_combined %>%
  filter(CDF %in% approx_target)

# Create a small offset for the boxplots to avoid overlap (shift by LegendGroup)
cdf_boxplot_data <- cdf_boxplot_data %>%
  mutate(CDF_offset = case_when(
    TRUE ~ CDF
  ))

cdf_boxplot_data$LegendGroup = TRUE

# Create the plot
gg2 <- ggplot(cdf_combined, aes(x = KGE_Eval, y = CDF, group = Model)) +
  # FUSE ensemble CDF
  geom_line(data = subset(cdf_combined, !IsBest), 
            color = color_dict["FUSE ensemble"], linewidth = 0.4, alpha = 0.2) +
  
  # Parents models CDF
  geom_line(data = subset(cdf_combined, IsParent), color = color_dict["Parents models"], linewidth = 0.5, alpha = 1) +
  
  # Best models CDF
  geom_line(data = subset(cdf_combined, IsBest), color = "black", linewidth = 1.8, alpha = 1) +
  geom_line(data = subset(cdf_combined, IsBest), 
            color = color_dict["Top-performing single model"], linewidth = 1.5, alpha = 1) +
  
  # Add vertical reference line
  geom_vline(xintercept = 1, col = "deepskyblue4", size = 0.8, linetype = "dashed") +
  
  # Add boxplots at CDF = 0.25, 0.5, and 0.75 for each LegendGroup
  geom_boxplot(data = cdf_boxplot_data, aes(x = KGE_Eval, y = CDF_offset, group = interaction(CDF, LegendGroup)),
                                            fill = color_dict["FUSE ensemble"],
               width = 0.05, alpha = 0.5, outlier.shape = NA) +  # `width` controls horizontal width
  
  # Labels and theme
  labs(x = expression(KGE[comp]), y = "Cumulative distribution function (CDF)", title = "(b) Evaluation period") +
  theme_bw() +
  coord_cartesian(xlim = c(-0.41, 1.02), ylim = c(-0.02, 1.02), expand = FALSE) +
  
  scale_x_continuous(
    breaks = seq(-0.4, 1.0, by = 0.2),   # major tick marks every 0.2
    minor_breaks = seq(-0.4, 1.0, by = 0.1)  # faint minor grid every 0.1
  )+
  
  # Remove legend for cleaner visualization
  theme(legend.position = "none")   

gg2


medians <- apply(Eval_FUSE_Comp[ , setdiff(names(Eval_FUSE_Comp), "rownames")],
                 2, median, na.rm = TRUE)
min_median <- min(medians)
max_median <- max(medians)
min_col <- names(medians)[which.min(medians)]
max_col <- names(medians)[which.max(medians)]
cat("Minimum median:", min_median, " (column:", min_col, ")\n")
cat("Maximum median:", max_median, " (column:", max_col, ")\n")

##-----------------------------------------
##---------------- Plot ------------------
##-----------------------------------------

combined_plot <- ggarrange(
  gg1, gg2,
  ncol = 2, nrow = 1,
  widths = c(1, 1)
)

combined_plot
# Save the plot
ggsave(filename = "99_Figures/f03.png",
       plot = combined_plot, 
       width = 10, height = 5, dpi = 600, units = "in",
       bg = "white")


