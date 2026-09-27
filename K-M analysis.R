library(survival)
library(survminer)
library(ggplot2)
library(ggpubr)
library(svglite)

theme_set(
  theme_bw(
    base_family = "Arial",
    base_size = 13
  )
)
update_geom_defaults("text", list(family = "Arial"))
update_geom_defaults("label", list(family = "Arial"))

dat <- read.csv(
  "KM_data.csv",
  check.names = FALSE
)

dat$histology <- factor(
  dat$histology,
  levels = c(1, 2),
  labels = c("EAH", "EEC")
)

dat_surv <- dat[
  complete.cases(
    dat[, c("treatment time", "treatment outcomes", "histology")]
  ),
]

fit <- survfit(
  Surv(treatment time, treatment outcomes) ~ histology,
  data = dat_surv
)

logrank_test <- survdiff(
  Surv(treatment time, treatment outcomes) ~ histology,
  data = dat_surv
)

df <- length(logrank_test$n) - 1

p_value <- pchisq(
  logrank_test$chisq,
  df = df,
  lower.tail = FALSE
)

pval_text <- if (p_value < 0.001) {
  "Log-rank test p < 0.001"
} else {
  sprintf("Log-rank test p = %.3f", p_value)
}

p <- ggsurvplot(
  fit,
  data = dat_surv,
  
  fun = "event",
  
  pval = pval_text,
  pval.size = 4.5,
  pval.coord = c(20, 0.10),
  
  risk.table = TRUE,
  risk.table.col = "strata",
  risk.table.fontsize = 6,
  risk.table.height = 0.30,
  
  xlab = "Fertility-sparing treatment time (months)",
  ylab = "Cumulative CR rate",
  ylim = c(0, 1.0),
  
  legend = "top",
  legend.title = "Histology",
  legend.labs = c("EAH", "EEC"),
  
  break.x.by = 5,
  break.y.by = 0.25,
  
  palette = c("#E7B800", "#2E9FDF"),
  
  ggtheme = theme_bw(
    base_family = "Arial",
    base_size = 13
  ) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 13),
      
      legend.text = element_text(size = 13),
      legend.title = element_text(size = 13),
      
      legend.direction = "horizontal",
      legend.justification = "center",
      
      plot.margin = grid::unit(
        c(1, 1, 1, 1),
        "lines"
      )
    )
)

p$table <- p$table +
  theme(
    text = element_text(family = "Arial"),
    legend.position = "none"
  )

combined_plot <- ggarrange(
  p$plot,
  p$table,
  ncol = 1,
  heights = c(2.5, 1),
  align = "v"
)

combined_plot
