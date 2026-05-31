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

dir_FUSE =  file.path("00_DATA")

# FUSE Mosa perf
Eval_FUSE_long_Mosa_P = loadRData(file.path(dir_FUSE,"KGEcomp", "Eval_FUSE.Rdata"))%>%
  mutate(`Cal : KGEcomp` = (`Cal : KGE :  1` + `Cal : KGE : -1`) / 2,
         `Eval : KGEcomp` = (`Eval : KGE :  1` + `Eval : KGE : -1`) / 2)

best_rows <- Eval_FUSE_long_Mosa_P %>%
  group_by(Codes) %>%
  slice_max(order_by = `Cal : KGEcomp`, n = 1, with_ties = FALSE) %>%
  ungroup()

best_decision_Mosa_P <- best_rows %>%
  pull(ModelDecisions)

Eval_FUSE_best_Mosa_P <- best_rows %>%
  transmute(
    Codes = as.character(Codes),
    KGEcomp = `Cal : KGEcomp`,
    Type = "Mosaic based on performance"
  )


# FUSE SA (static in time only)

Eval_FUSE_long_SA = loadRData(file.path(dir_FUSE,"KGEcomp", "Eval_SA.Rdata"))
Eval_FUSE_long_SA_Cal <- Eval_FUSE_long_SA["KGE", "Comp", "Cal", , ]
Eval_FUSE_long_SA_Eval <- Eval_FUSE_long_SA["KGE", "Comp", "Eval", , ]

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
  KGEcomp = as.numeric(Eval_FUSE_long_SA_Cal[cbind(best_decision_indices_SA_T, seq_along(best_decision_indices_SA_T))]),
  Type = "Spatially variable and temporally static combination"
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

Eval_FUSE_best_Mosa_P <- Eval_FUSE_best_Mosa_P[!Eval_FUSE_best_Mosa_P$Codes %in% catchments_fail, ]
Eval_FUSE_best_SA_T   <- Eval_FUSE_best_SA_T[!Eval_FUSE_best_SA_T$Codes %in% catchments_fail, ]

#######
df_diff <- Eval_FUSE_best_Mosa_P %>%
  rename(KGEcomp_Mosa = KGEcomp) %>%
  inner_join(
    Eval_FUSE_best_SA_T %>% rename(KGEcomp_SA_T = KGEcomp),
    by = "Codes"
  ) %>%
  mutate(diff = KGEcomp_Mosa - KGEcomp_SA_T)

pct_left  <- mean(df_diff$diff < 0) * 100
pct_left
pct_right <- mean(df_diff$diff > 0) * 100
pct_right

cut_pos <- sum(df_diff$diff < 0)

gg <- ggplot(df_diff, aes(x = reorder(Codes, diff), y = diff, group = 1)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = cut_pos, linetype = "dashed", color = "red", linewidth = 1.2) +
  annotate("text", 
           x = cut_pos / 2, 
           y = 0.09, 
           label = sprintf("Combination of 2 or 3 models is better\n(%.1f%%)", pct_left),
           hjust = 0.5) +
  annotate("text", 
           x = cut_pos + (nrow(df_diff) - cut_pos) / 2, 
           y = 0.09, 
           label = sprintf("Single model\nis better\n(%.1f%%)", pct_right),
           hjust = 0.5) +
  labs(
    title = expression(KGE[comp]~difference~
                       "(Mosaic based on performance - Spatially variable and temporally static combination)"),
    x = "Catchments",
    y = expression(Delta~KGE[comp]~"(Calibration period)")
  ) + 
  coord_cartesian(ylim = c(-0.1, 0.1)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),  
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

gg

ggsave(plot = gg, filename = "99_Figures/fC1.png",
       width = 10, height = 5, dpi = 300)
