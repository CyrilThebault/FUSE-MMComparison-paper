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
library(cowplot)

source(file = "Metrics.R")

dir_FUSE = file.path('00_DATA')



##-----------------------------------------
##-------------- VARIABLES ----------------
##-----------------------------------------

su_array = loadRData(file.path(dir_FUSE, "KGEcomp","SU_MM.Rdata"))

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

######## Shapefiles

## catchments outlets
csv_outlets_DB <- readLines(file.path(dir_FUSE,"gauge_information.txt"))
csv_outlets_DB = strsplit(csv_outlets_DB, "\t")
csv_outlets_DB[[1]] = c("HUC_02", "GAGE_ID", "GAGE_NAME", "LAT", "LON", "DRAINAGE AREA (KM^2)")
csv_outlets_DB = data.frame(matrix(unlist(csv_outlets_DB), nrow=length(csv_outlets_DB), byrow=TRUE),stringsAsFactors=FALSE)
colnames(csv_outlets_DB) = csv_outlets_DB[1,]
csv_outlets_DB = csv_outlets_DB[! csv_outlets_DB$HUC_02 == "HUC_02",]
csv_outlets_DB$Station_lon = csv_outlets_DB$LON
csv_outlets_DB$Station_lat = csv_outlets_DB$LAT
csv_outlets_DB$Station_id = csv_outlets_DB$GAGE_ID
csv_outlets_DB$Country = "USA"

shp_outlets_DB <- sf::st_as_sf(csv_outlets_DB, coords = c('Station_lon','Station_lat'), crs = 4326)

modif = nchar(shp_outlets_DB$Station_id) == 7
shp_outlets_DB$Station_id[shp_outlets_DB$Country == 'USA' & modif] = paste0("0", shp_outlets_DB$Station_id[shp_outlets_DB$Country == 'USA' & modif])

shp_outlets_DB$Codes = paste(shp_outlets_DB$Country, shp_outlets_DB$Station_id, sep='_')

shp_outlets = subset(shp_outlets_DB, Codes %in% dimnames(result_array)[[1]])

## rivers

## north-america
NorthAm = file.path("Shp","boundaries_p_2021_v3.shp")
NorthAm = read_sf(NorthAm)
NorthAm = st_transform(NorthAm, crs = 4326)
NorthAm = NorthAm[NorthAm$COUNTRY == "USA",]

drop_these = c('US-AK', 'US-HI', 'US-PR', 'US-VI')
NorthAm = NorthAm[! NorthAm$STATEABB %in% drop_these,]


##################################
##            MAP 1             ##
##################################

# Define fill colors
fill_colors <- c("Better" = "blue", "Equivalent" = "white", "Worse" = "red")


labels <- c(
  Comp   = "Benchmark\n(top-performing single\nmodel for everywhere)",
  MosaP  = "Mosaic based on\nperformance",
  MosaPE = "Mosaic based on\nperformance-\nequivalence",
  SATS   = "Spatially and\ntemporally static\ncombination",
  SAT    = "Spatially variable\nand temporally\nstatic combination",
  WA     = "Dynamic combination"
)

txt_size <- 4.08
txt_face <- "bold"
txt_lineheight <- 1.0

num_size <- txt_size
num_face <- txt_face

bbox <- st_bbox(NorthAm)

plot_list <- list()

for (i in seq_along(decisions)) {
  l <- decisions[i]
  
  for (j in seq_along(decisions)) {
    k <- decisions[j]
    
    if (l == k) {
      
      cell_plot <- ggplot() + theme_void()
      
    } else {
      
      new_col_name <- paste0(l, "_vs_", k)
      
      shp_outlets[[new_col_name]] <- factor(
        result_array[, l, k],
        levels = c("Worse", "Equivalent", "Better")
      )
      
      n_total <- sum(!is.na(shp_outlets[[new_col_name]]))
      
      pct_worse <- round(
        100 * sum(shp_outlets[[new_col_name]] == "Worse", na.rm = TRUE) / n_total
      )
      
      pct_equiv <- round(
        100 * sum(shp_outlets[[new_col_name]] == "Equivalent", na.rm = TRUE) / n_total
      )
      
      pct_better <- round(
        100 * sum(shp_outlets[[new_col_name]] == "Better", na.rm = TRUE) / n_total
      )
      
      label_plot <- ggplot() +
        annotate("text", x = 0.70, y = 0., label = pct_better,
                 color = "#4C6EDB", size = num_size, fontface = num_face, hjust = 1) +
        annotate("text", x = 0.72, y = 0., label = "/",
                 color = "black", size = num_size, fontface = num_face, hjust = 0.5) +
        annotate("text", x = 0.78, y = 0., label = pct_equiv,
                 color = "black", size = num_size, fontface = num_face, hjust = 0.5) +
        annotate("text", x = 0.83, y = 0., label = "/",
                 color = "black", size = num_size, fontface = num_face, hjust = 0.5) +
        annotate("text", x = 0.85, y = 0., label = pct_worse,
                 color = "red", size = num_size, fontface = num_face, hjust = 0) +
        coord_cartesian(xlim = c(0, 1), ylim = c(-0.3, 1), clip = "off") +
        theme_void() +
        theme(
          plot.margin = margin(0, 0, -6, 0)
        )
      
      map_plot <- ggplot() +
        geom_sf(
          data = NorthAm,
          fill = "grey90",
          color = "grey40",
          linewidth = 0.25
        ) +
        geom_sf(
          data = shp_outlets,
          aes(fill = .data[[new_col_name]]),
          size = 1.6,
          shape = 21,
          color = "black",
          linewidth = 0.25
        ) +
        scale_fill_manual(
          values = fill_colors,
          limits = c("Worse", "Equivalent", "Better"),
          na.value = "grey80",
          drop = FALSE
        ) +
        coord_sf(
          xlim = c(bbox["xmin"], bbox["xmax"]),
          ylim = c(bbox["ymin"], bbox["ymax"]),
          expand = TRUE
        ) +
        theme_void() +
        theme(
          legend.position = "none",
          plot.margin = margin(2, 2, 2, 2)
        )
      
      cell_plot <- plot_grid(
        label_plot,
        map_plot,
        ncol = 1,
        rel_heights = c(0.15, 1)
      )
    }
    
    plot_list[[paste0(l, "_vs_", k)]] <- cell_plot
  }
}

col_titles <- lapply(decisions, function(x) {
  ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = labels[[x]],
      fontface = txt_face,
      size = txt_size,
      lineheight = txt_lineheight
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
})

row_titles <- lapply(decisions, function(x) {
  ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = labels[[x]],
      angle = 90,
      fontface = txt_face,
      size = txt_size,
      lineheight = txt_lineheight
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
})

top_row <- plot_grid(
  ggplot() + theme_void(),
  plotlist = col_titles,
  nrow = 1,
  rel_widths = c(0.22, rep(1, length(decisions)))
)

map_rows <- list()

for (i in seq_along(decisions)) {
  l <- decisions[i]
  
  row_maps <- lapply(decisions, function(k) {
    plot_list[[paste0(l, "_vs_", k)]]
  })
  
  map_rows[[i]] <- plot_grid(
    row_titles[[i]],
    plotlist = row_maps,
    nrow = 1,
    rel_widths = c(0.22, rep(1, length(decisions)))
  )
}

main_grid <- plot_grid(
  top_row,
  plotlist = map_rows,
  ncol = 1,
  rel_heights = c(0.55, rep(1, length(decisions)))
)

legend_plot <- ggplot() +
  annotate(
    "text",
    x = 0.5,
    y = 1,
    label = "The approach in the row is...",
    hjust = 0,
    size = txt_size+0.25*txt_size,
    fontface = txt_face
  ) +
  geom_point(
    aes(x = 4.7, y = 1),
    shape = 21,
    size = txt_size+0.5*txt_size,
    fill = "#4C6EDB",
    color = "black"
  ) +
  annotate(
    "text",
    x = 5.25,
    y = 1,
    label = "Better",
    hjust = 0,
    size = txt_size+0.25*txt_size,
    fontface = txt_face
  ) +
  geom_point(
    aes(x = 8.0, y = 1),
    shape = 21,
    size = txt_size+0.5*txt_size,
    fill = "white",
    color = "black"
  ) +
  annotate(
    "text",
    x = 8.55,
    y = 1,
    label = "Equivalent",
    hjust = 0,
    size = txt_size+0.25*txt_size,
    fontface = txt_face
  ) +
  geom_point(
    aes(x = 12.0, y = 1),
    shape = 21,
    size = txt_size+0.5*txt_size,
    fill = "red",
    color = "black"
  ) +
  annotate(
    "text",
    x = 12.55,
    y = 1,
    label = "Worse",
    hjust = 0,
    size = txt_size+0.25*txt_size,
    fontface = txt_face
  ) +
  annotate(
    "text",
    x = 14.5,
    y = 1,
    label = "...than the one in the colum",
    hjust = 0,
    size = txt_size+0.25*txt_size,
    fontface = txt_face
  ) +
  xlim(0, 18) +
  ylim(0.5, 1.5) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(0, 0, 0, 0)
  )

gg1 <- plot_grid(
  main_grid,
  legend_plot,
  ncol = 1,
  rel_heights = c(1, 0.08)
)

ggsave(
  plot = gg1,
  filename = "99_Figures/fD1.png",
  width = 18,
  height = 12,
  dpi = 300,
  bg = "white"
)