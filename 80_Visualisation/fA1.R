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
library(forcats)

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

# FUSE comp

Eval_FUSE_long_Comp = loadRData(file.path(dir_FUSE,"KGEcomp", "Eval_FUSE.Rdata"))%>%
  mutate(`Cal : KGEcomp` = (`Cal : KGE :  1` + `Cal : KGE : -1`) / 2,
         `Eval : KGEcomp` = (`Eval : KGE :  1` + `Eval : KGE : -1`) / 2)%>%
  filter(!Codes %in% catchments_fail)

decision_medians_Comp <- Eval_FUSE_long_Comp %>%
  group_by(ModelDecisions) %>%
  summarise(median_Cal_KGEcomp = median(`Cal : KGEcomp`, na.rm = TRUE))

best_decision_Comp <- decision_medians_Comp %>%
  filter(median_Cal_KGEcomp == max(median_Cal_KGEcomp)) %>%
  pull(ModelDecisions)

df_lookup_Comp = data.frame(Codes = unique(Eval_FUSE_long_Comp$Codes), best_cal_model = best_decision_Comp)


# Reorder ModelDecisions by median Cal:KGEcomp
data_plot <- Eval_FUSE_long_Comp %>%
  mutate(ModelDecisions = factor(ModelDecisions)) %>%
  mutate(ModelDecisions = fct_reorder(ModelDecisions, `Cal : KGEcomp`, .fun = median, .desc = TRUE))

color_dict <- setNames(generate_random_colors(length(unique(data_plot$ModelDecisions))),
                       unique(sort(as.character(data_plot$ModelDecisions))))

gg <- ggplot(data_plot, aes(x = ModelDecisions, y = `Cal : KGEcomp`, fill = ModelDecisions)) +
  stat_summary(
    geom = "boxplot",
    fun.data = function(x) {
      data.frame(
        ymin   = quantile(x, 0.10, na.rm = TRUE),
        lower  = quantile(x, 0.25, na.rm = TRUE),
        middle = median(x, na.rm = TRUE),
        upper  = quantile(x, 0.75, na.rm = TRUE),
        ymax   = quantile(x, 0.90, na.rm = TRUE)
      )
    },
    position = position_dodge2(width = 0.9),
    width = 0.7,
    fill = "skyblue"
  ) +
  geom_hline(yintercept = 1, col = "deepskyblue4", size = 0.8, linetype = "dashed") +
  # scale_fill_manual(values = color_dict,
  #                   guide = 'none') +
  theme_bw() +
  labs(
    x = "Model structure ID",
    y = expression(KGE[comp]* " (calibration period)")
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )


ggsave(plot = gg, filename = "99_Figures/fA1.png",
       width = 10, height = 5, dpi = 300)

