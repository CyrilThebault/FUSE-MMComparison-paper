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

source(file = "Metrics.R")

generate_random_colors <- function(num_colors, seed = 1997) {
  set.seed(seed)
  colors <- sample(colors(), num_colors)
  return(colors)
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
##---------------- MAIN ------------------
##-----------------------------------------

dir_FUSE = file.path("00_DATA")

# FUSE SA (static in time and space)

Eval_FUSE_long_SA = loadRData(file.path(dir_FUSE,"KGEcomp", "Eval_SA.Rdata"))
Eval_FUSE_long_SA_Cal <- Eval_FUSE_long_SA["KGE", "Comp", "Cal", , !dimnames(Eval_FUSE_long_SA)$code %in% catchments_fail]
Eval_FUSE_long_SA_Eval <- Eval_FUSE_long_SA["KGE", "Comp", "Eval", , !dimnames(Eval_FUSE_long_SA)$code %in% catchments_fail]

decision_medians_SA_TS <- apply(Eval_FUSE_long_SA_Cal, 1, median, na.rm = TRUE)
best_decision_SA_TS <- names(decision_medians_SA_TS)[which.max(decision_medians_SA_TS)]

df_lookup_SA_TS = data.frame(Codes = dimnames(Eval_FUSE_long_SA)[["code"]], best_cal_model = best_decision_SA_TS)


# Reorder ModelDecisions by median Cal:KGEcomp
data_plot <- as.data.frame(Eval_FUSE_long_SA_Cal) %>%
  mutate(model = rownames(Eval_FUSE_long_SA_Cal)) %>%
  pivot_longer(
    cols = -model,
    names_to = "code",
    values_to = "performance"
  )%>%
  mutate(model = fct_reorder(model, performance, .fun = median, .desc = TRUE))%>%
  rename(
    ModelDecisions = model,
    Codes = code,
    `Cal : KGEcomp` = performance
  )


# Compute quantiles per model
summary_stats <- data_plot %>%
  group_by(ModelDecisions) %>%
  summarise(
    q10 = quantile(`Cal : KGEcomp`, 0.10, na.rm = TRUE),
    q25 = quantile(`Cal : KGEcomp`, 0.25, na.rm = TRUE),
    q50 = quantile(`Cal : KGEcomp`, 0.50, na.rm = TRUE),
    q75 = quantile(`Cal : KGEcomp`, 0.75, na.rm = TRUE),
    q90 = quantile(`Cal : KGEcomp`, 0.90, na.rm = TRUE)
  ) %>%
  mutate(ModelDecisions = fct_reorder(ModelDecisions, q50, .desc = TRUE))


# Plot quantile lines
gg <- ggplot(summary_stats, aes(x = ModelDecisions)) +
  geom_linerange(aes(ymin = q10, ymax = q90), color = "grey40", size = 0.3) +
  geom_linerange(aes(ymin = q25, ymax = q75), color = "steelblue", size = 0.6) +
  geom_point(aes(y = q50), color = "black", size = 0.5) +
  geom_line(aes(y = q50, group = 1), color = "black", size = 0.3) +
  geom_hline(yintercept = 1, color = "deepskyblue4", size = 0.8, linetype = "dashed") +
  scale_x_discrete(
    breaks = function(x) x[seq(1, length(x), by = 10000)]
  ) +
  
  geom_segment(
    x = summary_stats$ModelDecisions[2000],
    xend = summary_stats$ModelDecisions[1],
    y = 0.91,         
    yend = summary_stats$q50[1],
    arrow = arrow(length = unit(0.15, "cm")),
    color = "red",
    size = 0.7
  ) +
  
  # Label next to arrow
  geom_text(
    x = summary_stats$ModelDecisions[2000], 
    y = 0.92, label = gsub("_", " & ",summary_stats$ModelDecisions[1]),
    color = "red",
    size = 3,
    hjust = 0
  ) +
  
  
  theme_bw() +
  labs(
    x = "Model structure combination ID (two or three models)",
    y = expression(KGE[comp]*" (calibration period)")
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

ggsave(plot = gg, filename = "99_Figures/fA5.png",
       width = 10, height = 5, dpi = 300)

