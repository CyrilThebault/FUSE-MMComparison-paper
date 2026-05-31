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

best_decision_indices_SA_T <- vapply(
  seq_len(ncol(Eval_FUSE_long_SA_Cal)),
  function(i) {
    col <- Eval_FUSE_long_SA_Cal[, i]
    if (all(is.na(col))) NA_integer_ else which.max(col)
  },
  FUN.VALUE = integer(1)
)
best_decision_SA_T <- dimnames(Eval_FUSE_long_SA_Cal)$model[best_decision_indices_SA_T]

df_lookup_SA_T = data.frame(Codes = dimnames(Eval_FUSE_long_SA_Cal)[["code"]], best_cal_model = best_decision_SA_T)


# ---- 1) Frequency (in %) of each combination (e.g., "72_126", "48_50_72") ----
comb_freq <- df_lookup_SA_T %>%
  filter(!is.na(best_cal_model) & best_cal_model != "") %>%
  count(best_cal_model, name = "n") %>%
  mutate(pct = n / sum(n)) %>%
  arrange(desc(pct)) %>%
  mutate(best_cal_model = fct_reorder(best_cal_model, pct, .desc = TRUE))

p_comb <- ggplot(comb_freq, aes(x = best_cal_model, y = n)) +
  geom_col(fill="skyblue", col = "black", linewidth = 0.3) +
  labs(
    x = "Model structure combination ID (two or three models)",
    y = "Number of catchments"
  ) +
  scale_x_discrete(
    breaks = function(x) x[seq(1, length(x), by = 50)]
  ) +
  theme_bw() +
  theme(
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

p_comb


ggsave(filename = "99_Figures/fA6.png",
       plot = p_comb, 
       width = 14, height = 3.25, dpi = 300, units = "in",
       bg = "white")

# ---- 2) Frequency (in %) of each model ID within combinations ----
# Splits "72_126" into "72" and "126" (and handles combos with 3+ parts too)
model_freq <- df_lookup_SA_T %>%
  filter(!is.na(best_cal_model) & best_cal_model != "") %>%
  separate_rows(best_cal_model, sep = "_", convert = FALSE) %>%
  rename(model_id = best_cal_model) %>%
  count(model_id, name = "n") %>%
  mutate(
    model_num = as.numeric(as.character(model_id)),
    pct = n / sum(n)
  ) %>%
  arrange(desc(n), model_num) %>%
  mutate(model_id = factor(model_id, levels = model_id))


indiv <- unique(unlist(strsplit(attr(Eval_FUSE_long_SA_Cal, "dimnames")$model,"_")))

color_dict <- setNames(generate_random_colors(length(indiv)),
                       unique(sort(as.character(indiv))))

p_model <- ggplot(model_freq, aes(x = model_id, y = n, fill = model_id)) +
  geom_col(col = "black", linewidth = 0.3) +
  labs(
    x = "Model structure ID (individual model within the combinations)",
    y = "Number of catchments"
  ) +
  scale_fill_manual(values = color_dict,
                    guide = 'none') +
  theme_bw() +
  theme(
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

p_model

ggsave(filename = "99_Figures/fA7.png",
       plot = p_model, 
       width = 14, height = 3.25, dpi = 300, units = "in",
       bg = "white")







