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
library(ggspatial)
library(ggpubr)
library(grid)
library(dplyr)
library(tidyr)
library(abind)
library(gridExtra)
library(patchwork)
library(cowplot)

source(file = "Metrics.R")


dir_FUSE = file.path('00_DATA')

##-----------------------------------------
##-------------- VARIABLES ----------------
##-----------------------------------------

su_array = loadRData(file.path(dir_FUSE, "KGEcomp", "SU_MM.Rdata"))

su_numeric <- su_array[, , , setdiff(dimnames(su_array)[[4]], "GOF_stat"), drop=FALSE]
storage.mode(su_numeric) <- "numeric"

# Compute range
range_array <- su_numeric[,,, "p95", drop = FALSE] -
  su_numeric[,,, "p05", drop = FALSE]

dimnames(range_array)[[4]] <- "range"

# Attach using abind
su <- abind(su_numeric, range_array, along = 4)

# Extract relevant dimension names
catchments <- dimnames(su)[[1]]
decisions <- dimnames(su)[[2]]
metrics <- dimnames(su)[[3]]
variables <- dimnames(su)[[4]]


# Create an array to store the classification results
result_array <- array(NA, dim = c(length(catchments), length(metrics), length(decisions), length(decisions)),
                      dimnames = list(catchments, metrics, decisions, decisions))

# Classification function
classify_performance <- function(p05_M, p95_M, score_M, p05_W, p95_W, score_W) {
  if (is.na(p05_M) | is.na(p95_M) | is.na(p05_W) | is.na(p95_W)) {
    return(NA)
  }
  
  highest_score = max(c(score_M, score_W))
  lowest_score = min(c(score_M, score_W))
  
  p05 = ifelse(score_M == highest_score, p05_M, p05_W)
  p95 = ifelse(score_M == highest_score, p95_M, p95_W)
  
  if (score_W == highest_score & p05 > lowest_score) {
    return("Better")
  } else if (score_M == highest_score & p05 > lowest_score) {
    return("Worse")
  } else {
    return("Equivalent")
  }
}

# Loop through catchments, metrics, and other decisions
for (i in seq_along(catchments)) {
  for (j in seq_along(metrics)) {
    for (k in seq_along(decisions)) {
      for (l in seq_along(decisions)){
        
        # Extract required values
        p05_W <- su[i, k, j, which(variables == "p05")]
        p95_W <- su[i, k, j, which(variables == "p95")]
        score_W <- su[i, k, j, which(variables == "p50")]
        
        p05_M <- su[i, l, j, which(variables == "p05")]
        p95_M <- su[i, l, j, which(variables == "p95")]
        score_M <- su[i, l, j, which(variables == "p50")]
        
        # Apply classification
        result_array[i, j, k, l] <- classify_performance(p05_M, p95_M, score_M, p05_W, p95_W, score_W)
      }
    }
  }
}

# Identify rows (catchments) that contain any NA across the full 4D array
na_catchments <- dimnames(result_array)[[1]][apply(result_array, 1, function(x) any(is.na(x)))]

# Drop those rows
result_array <- result_array[!(dimnames(result_array)[[1]] %in% na_catchments), , ,]


df_long <- as.data.frame.table(
  apply(result_array, c(2, 3), function(x) {
    round(prop.table(table(factor(x, levels = c("Better", "Equivalent", "Worse")))) * 100)
  })
)

names(df_long) <- c("category", "reference", "compared", "percent")

df_long <- df_long[, c("reference", "compared", "category", "percent")]

df_summary_named <- df_long %>%
  pivot_wider(
    names_from = category,
    values_from = percent
  )

labels <- c(
  Comp   = "Benchmark\n(top-performing single\nmodel for everywhere)",
  MosaP  = "Mosaic based on\nperformance",
  MosaPE = "Mosaic based on\nperformance-equivalence",
  SATS   = "Spatially and\ntemporally static\ncombination",
  SAT    = "Spatially variable\nand temporally\nstatic combination",
  WA     = "Dynamic combination"
)

df_summary_named$reference <- factor(
  df_summary_named$reference,
  levels = names(labels),
  labels = labels
)

df_summary_named$compared <- factor(
  df_summary_named$compared,
  levels = names(labels),
  labels = labels
)


approaches = unique(df_summary_named$reference)

plot_tags <- paste0("(", letters[seq_along(approaches)], ")")

plots <- lapply(seq_along(approaches), function(i) {
  
  ref = approaches[i]
  
  df_tmp <- df_summary_named %>%
    filter(reference == ref) %>%
    mutate(
      compared = factor(compared, levels = rev(approaches))
    ) %>%
    pivot_longer(
      cols = c(Better, Equivalent, Worse),
      names_to = "category",
      values_to = "percent"
    ) %>%
    group_by(reference, compared) %>%
    mutate(
      percent_plot = percent / sum(percent) * 100
    ) %>%
    arrange(category) %>%
    mutate(
      cum = cumsum(percent_plot),
      mid = cum - percent_plot / 2, 
      small = percent < 9
    ) %>%
    ungroup() %>%
    mutate(
      category = factor(category, levels = c("Worse", "Equivalent", "Better"))
    )
  
  ggplot(df_tmp, aes(x = compared, y = percent_plot, fill = category)) +
    geom_bar(stat = "identity", color = "black") +
    coord_flip() +
    
    scale_fill_manual(values = c(
      Better = "#4C6EDB",
      Equivalent = "white",
      Worse = "red"
    )) +
    
    # Labels normaux
    geom_text(
      data = subset(df_tmp, !small),
      aes(y = mid, label = ifelse(percent == 0, "", paste0(percent, "%"))),
      size = 3
    ) +
    
    # Labels déportés
    geom_text(
      data = subset(df_tmp, small & percent > 0),
      aes(
        y = ifelse(category == "Better", mid + 10, mid - 10),  
        label = paste0(percent, "%")
      ),
      size = 3
    ) +
    
    # Flèches horizontales
    geom_segment(
      data = subset(df_tmp, small & percent > 0),
      aes(
        x = compared,
        xend = compared,
        y = ifelse(category == "Better", mid + 5, mid - 5),
        yend = mid
      ),
      arrow = arrow(length = unit(0.2, "cm")),
      inherit.aes = FALSE
    ) +
    
    labs(
      title = paste(plot_tags[i],gsub("\n", " ", ref)),
      x = "",
      y = "Percentage"
    ) +
    
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 10),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    )
})

# afficher
plots[[1]]


legend_plot <- ggplot() +
  
  # Texte gauche
  annotate("text",x = 0, y = 1, label = "The approach in the title is...", hjust = 0, size = 4) +
  
  # Points + labels
  geom_point(aes(x = 5.8, y = 1), shape = 21, size = 6,
             fill = "#4C6EDB", color = "black") +
  annotate("text", x = 6.5, y = 1, label = "Better", hjust = 0, size = 4)+
  
  geom_point(aes(x = 9.2, y = 1), shape = 21, size = 6,
             fill = "white", color = "black") +
  annotate("text", x = 9.9, y = 1, label = "Equivalent", hjust = 0, size = 4)+
  
  geom_point(aes(x = 13.2, y = 1), shape = 21, size = 6,
             fill = "red", color = "black") +
  annotate("text", x = 13.9, y = 1, label = "Worse", hjust = 0, size = 4)+
  
  # Texte droite
  annotate(
    "text",
    x = 16,
    y = 1,
    label = "...than the one in the y-axis",
    hjust = 0,
    size = 4
  ) +
  
  xlim(0, 20) +
  ylim(0.5, 1.5) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  )

legend_plot

gg_main <- wrap_plots(plots[1:6], nrow = 3, ncol = 2)

gg <- plot_grid(
  gg_main,
  legend_plot,
  ncol = 1,
  rel_heights = c(1, 0.08)
)

gg


ggsave(plot = gg, filename = "99_Figures/f09.png",
       width = 9.5, height = 10, dpi = 300)
