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

source(file = "Metrics.R")

##-----------------------------------------
##---------------- MAIN ------------------
##-----------------------------------------

dir_FUSE = file.path("00_DATA")

# FUSE Comp
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
df_eval_gumboot = read.table(file.path(dir_FUSE, "KGEcomp", "brief_analysis_modeling_results_with_gumboot_eval.csv"),
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
df_eval_gumboot = read.table(file.path(dir_FUSE, "KGEcomp", "brief_analysis_modeling_results_with_gumboot_eval.csv"),
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
# 
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

Eval_FUSE_best_Comp   <- Eval_FUSE_best_Comp[!Eval_FUSE_best_Comp$Codes %in% catchments_fail, ]
Eval_FUSE_best_Mosa_P <- Eval_FUSE_best_Mosa_P[!Eval_FUSE_best_Mosa_P$Codes %in% catchments_fail, ]
Eval_FUSE_best_Mosa_PE <- Eval_FUSE_best_Mosa_PE[!Eval_FUSE_best_Mosa_PE$Codes %in% catchments_fail, ]
Eval_FUSE_best_SA_TS  <- Eval_FUSE_best_SA_TS[!Eval_FUSE_best_SA_TS$Codes %in% catchments_fail, ]
Eval_FUSE_best_SA_T   <- Eval_FUSE_best_SA_T[!Eval_FUSE_best_SA_T$Codes %in% catchments_fail, ]
Eval_FUSE_best_WA     <- Eval_FUSE_best_WA[!Eval_FUSE_best_WA$Codes %in% catchments_fail, ]


# Only to get statistics

Eval_FUSE_best_All <- Reduce(function(x, y) merge(x, y, by = "Codes"),
                             list(
                               Eval_FUSE_best_Comp[, c("Codes", "KGEcomp")],
                               Eval_FUSE_best_Mosa_P[, c("Codes", "KGEcomp")],
                               Eval_FUSE_best_Mosa_PE[, c("Codes", "KGEcomp")],
                               Eval_FUSE_best_SA_TS[, c("Codes", "KGEcomp")],
                               Eval_FUSE_best_SA_T[, c("Codes", "KGEcomp")],
                               Eval_FUSE_best_WA[, c("Codes", "KGEcomp")]
                             ))
colnames(Eval_FUSE_best_All) <- c("Codes", "Comp", "Mosa_P", "Mosa_PE", "SA_TS", "SA_T", "WA")

summary_stats <- data.frame(
  Min  = apply(Eval_FUSE_best_All[, -1], 2, min,  na.rm = TRUE),
  Q10  = apply(Eval_FUSE_best_All[, -1], 2, quantile, probs = 0.10, na.rm = TRUE),
  Q25  = apply(Eval_FUSE_best_All[, -1], 2, quantile, probs = 0.25, na.rm = TRUE),
  Mean = apply(Eval_FUSE_best_All[, -1], 2, mean, na.rm = TRUE),
  Q50  = apply(Eval_FUSE_best_All[, -1], 2, quantile, probs = 0.50, na.rm = TRUE),
  Q75  = apply(Eval_FUSE_best_All[, -1], 2, quantile, probs = 0.75, na.rm = TRUE),
  Q90  = apply(Eval_FUSE_best_All[, -1], 2, quantile, probs = 0.90, na.rm = TRUE),
  Max  = apply(Eval_FUSE_best_All[, -1], 2, max,  na.rm = TRUE)
)
print(summary_stats)


# Create data frame for the plots

Eval_FUSE_best_Comp = Eval_FUSE_best_Comp[,c("Type",  "KGEcomp")]
Eval_FUSE_best_Mosa_P = Eval_FUSE_best_Mosa_P[,c("Type",  "KGEcomp")]
Eval_FUSE_best_Mosa_PE = Eval_FUSE_best_Mosa_PE[,c("Type",  "KGEcomp")]
Eval_FUSE_best_SA_TS = Eval_FUSE_best_SA_TS[,c("Type",  "KGEcomp")]
Eval_FUSE_best_SA_T = Eval_FUSE_best_SA_T[,c("Type",  "KGEcomp")]
Eval_FUSE_best_WA = Eval_FUSE_best_WA[,c("Type",  "KGEcomp")]


# Ensure the rownames are preserved as a column for merging
combined_df <-rbind(Eval_FUSE_best_Comp,
                    Eval_FUSE_best_Mosa_P, Eval_FUSE_best_Mosa_PE,
                    Eval_FUSE_best_SA_TS, Eval_FUSE_best_SA_T, 
                    Eval_FUSE_best_WA
)



##########################
#--------Boxplot---------#
##########################

median_benchmark = median(Eval_FUSE_best_Comp$KGEcomp)

# Define color palette
color_dict <- c(
  "Benchmark (top-performing single model for everywhere)" = "darkorange",
  "Mosaic based on performance" = "lightpink1",
  "Mosaic based on performance-equivalence" = "lightpink4",
  "Spatially and temporally static combination" = "olivedrab2",
  "Spatially variable and temporally static combination" = "olivedrab4",
  "Dynamic combination" = "firebrick4"
)


# Convert Type to factor with ordered levels
combined_df$Type <- factor(combined_df$Type, 
                           levels = c("Dynamic combination", "Spatially variable and temporally static combination", "Spatially and temporally static combination",
                                          "Mosaic based on performance-equivalence", "Mosaic based on performance", 
                                          "Benchmark (top-performing single model for everywhere)"))

# Create gg
gg <- ggplot(combined_df, aes(x = `KGEcomp`, y = Type, fill = Type)) +
  stat_summary(geom = "boxplot", 
               fun.data = function(x) {
                 data.frame(
                   y = median(x),
                   ymin = quantile(x, 0.10),
                   lower = quantile(x, 0.25),
                   middle = median(x),
                   upper = quantile(x, 0.75),
                   ymax = quantile(x, 0.90)
                 )
               },
               position = position_dodge2()) +
  scale_fill_manual(values = color_dict) +
  # coord_cartesian(ylim = c(-0.41,1))+
  geom_vline(xintercept = 1, col = "deepskyblue4", size = 0.8, linetype = "dashed") +
  geom_vline(xintercept = median_benchmark, col = "black", size = 0.5, linetype = "dashed") +
  labs(y = NULL, x = expression(KGE[comp]), title = "") +
  guides(fill = guide_legend(title = "Modelling approach"))+
  theme_bw() +
  theme(legend.position = "none")
  # theme(axis.text.x = element_blank())
  
# Arrange vertically
gg


ggsave(plot = gg, filename = "99_Figures/f05.png",
       width = 10, height = 5, dpi = 300)
