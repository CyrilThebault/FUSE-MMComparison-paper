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

source(file = "Metrics.R")


dir_FUSE = file.path("00_DATA")


##-----------------------------------------
##-------------- VARIABLES ----------------
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

##################################
##           Boxplot            ##
##################################

range_data <- su[,c("Comp","MosaP", "MosaPE", "SATS","SAT","WA"),c("KGEcomp"),"range"]

median_benchmark = median(range_data[,"Comp"])

df_range <- melt(range_data, varnames = c("catchment", "model"), value.name = "range")

df_range$model <- factor(df_range$model,
                         levels = c("WA","SAT", "SATS", "MosaPE","MosaP","Comp"),
                         labels = c("Dynamic combination", "Spatially variable and temporally static combination", "Spatially and temporally static combination",
                                    "Mosaic based on performance-equivalence", "Mosaic based on performance", 
                                    "Benchmark (top-performing single model for everywhere)"))


# Define color palette
color_dict <- c(
  "Benchmark (top-performing single model for everywhere)" = "darkorange",
  "Mosaic based on performance" = "lightpink1",
  "Mosaic based on performance-equivalence" = "lightpink4",
  "Spatially and temporally static combination" = "olivedrab2",
  "Spatially variable and temporally static combination" = "olivedrab4",
  "Dynamic combination" = "firebrick4"
)


# Create plot 
gg1 <- ggplot(df_range, aes(x = range, y = model, fill = model)) +
  stat_summary(
    geom = "boxplot",
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
    position = position_dodge2()
  ) +
  geom_vline(xintercept = 0, col = "deepskyblue4", size = 0.8, linetype = "dashed") +
  geom_vline(xintercept = median_benchmark, col = "black", size = 0.5, linetype = "dashed") +
  scale_fill_manual(values = color_dict) +
  coord_cartesian(xlim = c(0, 0.45)) +
  labs(
    y = NULL,
    x = "Sampling uncertainty (p95 - p05 from gumboot)"
  ) +
  guides(fill = "none") + 
  theme_bw()


gg1

ggsave(gg1, filename = "99_Figures/f07.png",
       width = 10, height = 5, dpi = 300)
