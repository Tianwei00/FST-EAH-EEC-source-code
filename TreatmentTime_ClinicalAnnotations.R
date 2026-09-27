suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(grid)
})

font_family <- "Arial"

dat <- read.csv(
  "clinical_data.csv",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_columns <- c(
  "ID",
  "Treatment outcome",
  "PR",
  "ER",
  "Age",
  "Histology",
  "Molecular classification",
  "Treatment time"
)

stopifnot(all(required_columns %in% colnames(dat)))

dat[["Molecular classification"]] <- factor(
  dat[["Molecular classification"]],
  levels = c("POLEmut", "MMRd", "NSMP", "p53abn")
)

dat[["Treatment outcome"]] <- factor(
  dat[["Treatment outcome"]],
  levels = c("CR", "No-CR")
)

dat$Histology <- factor(
  dat$Histology,
  levels = c("EAH", "EEC")
)

dat$Age <- factor(
  dat$Age,
  levels = c("<30", "≥30")
)

dat$ER <- factor(
  dat$ER,
  levels = c("High", "Low", "Deficient")
)

dat$PR <- factor(
  dat$PR,
  levels = c("High", "Low", "Deficient")
)

dat[["Treatment time"]] <- as.numeric(dat[["Treatment time"]])

dat <- dat[
  order(
    dat[["Molecular classification"]],
    dat[["Treatment time"]],
    dat$ID
  ),
]

rownames(dat) <- dat$ID

molecular_col <- c(
  "POLEmut" = "#C1BD38",
  "MMRd"    = "#DE8E69",
  "NSMP"    = "#F9A825",
  "p53abn"  = "#71702F"
)

histology_col <- c(
  "EAH" = "#99B898",
  "EEC" = "#C0392B"
)

age_col <- c(
  "<30" = "#D7ACA1",
  "≥30" = "#C2697F"
)

er_col <- c(
  "High"      = "#2B8CBE",
  "Low"       = "#9ECAE1",
  "Deficient" = "#CCCCCC"
)

pr_col <- er_col

outcome_col <- c(
  "CR"    = "#08519C",
  "No-CR" = "#CB181D"
)

annotation_colours <- list(
  "Molecular classification" = molecular_col,
  "Histology" = histology_col,
  "Age" = age_col,
  "ER" = er_col,
  "PR" = pr_col,
  "Treatment outcome" = outcome_col
)

y_max <- ceiling(
  max(dat[["Treatment time"]], na.rm = TRUE) * 1.05 / 5
) * 5

y_breaks <- pretty(
  c(0, y_max),
  n = 5
)

y_breaks <- y_breaks[
  y_breaks >= 0 & y_breaks <= y_max
]

point_colours <- outcome_col[
  as.character(dat[["Treatment outcome"]])
]

treatment_time_plot <- AnnotationFunction(
  which = "column",
  height = unit(45, "mm"),
  
  fun = function(index) {
    
    pushViewport(
      viewport(
        xscale = c(0.5, length(index) + 0.5),
        yscale = c(0, y_max),
        clip = "off"
      )
    )
   
    grid.rect(
      gp = gpar(
        fill = "#F0F0F0",
        col = NA
      )
    )
    
    grid.points(
      x = seq_along(index),
      y = dat[["Treatment time"]][index],
      default.units = "native",
      pch = 16,
      size = unit(1.05, "mm"),
      gp = gpar(
        col = point_colours[index],
        fill = point_colours[index]
      )
    )
    
    grid.yaxis(
      at = y_breaks,
      label = y_breaks,
      gp = gpar(
        fontfamily = font_family,
        fontsize = 9,
        col = "black",
        lwd = 0.8
      )
    )
    
    grid.text(
      "Treatment time (months)",
      x = unit(-13, "mm"),
      y = unit(0.5, "npc"),
      rot = 90,
      gp = gpar(
        fontfamily = font_family,
        fontface = "bold",
        fontsize = 10
      )
    )
    
    popViewport()
  }
)

legend_parameters <- lapply(
  names(annotation_colours),
  function(x) {
    list(
      title_gp = gpar(
        fontfamily = font_family,
        fontface = "bold",
        fontsize = 9
      ),
      labels_gp = gpar(
        fontfamily = font_family,
        fontsize = 9
      ),
      grid_width = unit(4, "mm"),
      grid_height = unit(4, "mm")
    )
  }
)

names(legend_parameters) <- names(annotation_colours)

top_annotation <- HeatmapAnnotation(
  "Treatment time" = treatment_time_plot,
  
  "Molecular classification" =
    dat[["Molecular classification"]],
  
  "Histology" = dat$Histology,
  "Age" = dat$Age,
  "ER" = dat$ER,
  "PR" = dat$PR,
  
  "Treatment outcome" =
    dat[["Treatment outcome"]],
  
  col = annotation_colours,
  
  annotation_height = unit.c(
    unit(45, "mm"),
    rep(unit(4.2, "mm"), 6)
  ),
  
  gap = unit(
    c(2.5, rep(0.6, 5)),
    "mm"
  ),
  
  border = FALSE,
  na_col = "#D9D9D9",

  annotation_name_side = "left",
  annotation_name_rot = 0,

  show_annotation_name = c(
    FALSE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE
  ),
  
  annotation_name_gp = gpar(
    fontfamily = font_family,
    fontface = "bold",
    fontsize = 9
  ),
  
  annotation_legend_param = legend_parameters
)

make_figure <- function() {
  
  dummy_matrix <- matrix(
    0,
    nrow = 1,
    ncol = nrow(dat)
  )
  
  colnames(dummy_matrix) <- dat$ID
  
  ht <- Heatmap(
    dummy_matrix,
    top_annotation = top_annotation,
    
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    
    show_row_names = FALSE,
    show_column_names = FALSE,
    show_heatmap_legend = FALSE,
    
    rect_gp = gpar(type = "none"),
    height = unit(0.1, "mm"),
    border = FALSE
  )
  
  draw(
    ht,
    annotation_legend_side = "right",
    heatmap_legend_side = "right",
    merge_legends = FALSE,
    
    padding = unit(
      c(5, 8, 5, 20),
      "mm"
    )
  )
}

make_figure()