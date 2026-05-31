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
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(abind)
library(reshape2)
library(viridis)
library(patchwork)

source(file = "Metrics.R")


dir_FUSE = file.path('00_DATA')


##-----------------------------------------
##--------------- Gumboot -----------------
##-----------------------------------------

# Sampling uncertainty over the evaluation period for all multi-model approaches
su_array = loadRData(file.path(dir_FUSE, "KGEcomp", "SU_MM.Rdata"))

su_numeric <- su_array[, , , setdiff(dimnames(su_array)[[4]], "GOF_stat"), drop=FALSE]
storage.mode(su_numeric) <- "numeric"

# Compute range
range_array <- su_numeric[,,, "p95", drop = FALSE] -
  su_numeric[,,, "p05", drop = FALSE]

dimnames(range_array)[[4]] <- "range"

# Attach using abind
su <- abind(su_numeric, range_array, along = 4)

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

su <- su[!(dimnames(su)[[1]] %in% catchments_fail), , , ,drop=FALSE]

##-----------------------------------------
##------------- Evaluation ----------------
##-----------------------------------------

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

Eval_FUSE_best_Comp <- Eval_FUSE_best_Comp[!Eval_FUSE_best_Comp$Codes %in% catchments_fail,]

##################################
##        ScaterPlot            ##
##################################

range_data <- su[,c("Comp","MosaP", "MosaPE", "SATS","SAT","WA"),c("KGEcomp"),"range"]

# Convert array to data frame (you already did this)
df <- as.data.frame(range_data)
methods <- setdiff(names(df), "Comp")

df$Codes = rownames(df)

method_labels <- c(
  WA      = "Sampling uncertainty of the\ndynamic combination",
  SAT     = "Sampling uncertainty of the\nspatially variable and temporally static combination",
  SATS    = "Sampling uncertainty of the\nspatially and temporally static combination",
  MosaPE  = "Sampling uncertainty of the\nmosaic based on performance-equivalence",
  MosaP   = "Sampling uncertainty of the\nmosaic based on performance",
  Comp    = "Sampling uncertainty of the\nbenchmark"
)

df <- df %>%
  left_join(
    Eval_FUSE_best_Comp %>% select(Codes, KGE = KGEcomp),
    by = "Codes"
  )




#### Calculate percentages

comp_col="Comp"

df_compare <- df
for (m in methods) {
  df_compare[[paste0(m, "_vs_Comp")]] <- ifelse(
    df[[m]] > df[[comp_col]], ">",
    ifelse(df[[m]] < df[[comp_col]], "<", "=")
  )
}

head(df_compare[ , c(comp_col, methods, paste0(methods, "_vs_Comp"))])

levels3 <- c(">", "=", "<")

count_mat <- sapply(methods, function(m) {
  tab <- table(factor(df_compare[[paste0(m, "_vs_Comp")]], levels = levels3))
  as.integer(tab)
})
rownames(count_mat) <- levels3

# Percentages
pct_mat <- round(100 * sweep(count_mat, 2, colSums(count_mat), "/"), 1)

counts_tbl <- data.frame(Case = rownames(count_mat), count_mat, row.names = NULL)
colnames(counts_tbl)[-1] <- methods

pct_tbl <- data.frame(Case = rownames(pct_mat), pct_mat, row.names = NULL)
colnames(pct_tbl)[-1] <- methods

counts_tbl   
pct_tbl      


### PLOT
plot_tags <- paste0("(", letters[seq_along(methods)], ")")

plots <- lapply(seq_along(methods), function(i) {
  
  m <- methods[i]
  
  print(m)
  
  pct_less <- pct_tbl[pct_tbl$Case == "<", m]
  pct_greater <- pct_tbl[pct_tbl$Case == ">", m]
  
  ggplot(df, aes_string(x = m, y = "Comp", color = "KGE")) +
    scale_colour_gradientn(
      limits = c(-0.42, 1),
      colors = c(rep("purple",2), viridis(80)),
      breaks = seq(-0.4,1,0.1),
      labels = c("< -0.4", seq(-3,10,1)/10),
      oob = scales::squish,
      guide = guide_colorbar(title = expression("Benchmark" ~ KGE[comp]), barwidth = 1.5, barheight = 10)
    ) +
    geom_point(alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    annotate("text", x = 0.05, y = 1.32, label = paste0(round(pct_less), "%"), hjust = 0, size = 5) +
    annotate("text", x = 1.35, y = 0.05, label = paste0(round(pct_greater), "%"), hjust = 1, size = 5) +
    coord_cartesian(xlim = c(0, 1.4), ylim = c(0, 1.4)) +
    labs(
      title = plot_tags[i],
      x = method_labels[m],
      y = method_labels["Comp"]
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(size = 12),
      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 9)
    )
})

# Combine all 5 plots in 2 rows and add the legend in the empty spot

combined_plot <-
  wrap_plots(
    c(plots, list(guide_area())),  
    nrow = 2
  ) +
  plot_layout(guides = "collect") 

# Display it
combined_plot


ggsave(combined_plot, filename = "99_Figures/f08.png",
       width = 10, height = 7, dpi = 300)





