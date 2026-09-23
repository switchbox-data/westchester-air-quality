# Executive-summary COMPOSITE charts for the report, in the Switchbox house style.
#
# Builds the diesel-first health-benefit curves: avoided cases per 100,000 per
# year (region-level mean) vs. share of fleet electrified, one line per health
# outcome. Produces a PM2.5 composite (4 outcomes; primary) and an NO2 companion
# (3 outcomes).
#
# Reads ONLY cached results (Results/Mitigation_EV_diesel_health_benefit.RData) —
# no recomputation, no Census API. Run inside the project container with Code,
# Data, Results and Figures bind-mounted.
#
# Data shape (see Code/04_Mitigation_strategies.R, sec. (2)/(3)):
#   HB_EV_diesel_rate_tot_stat[[k]] is the region-level avoided RATE per 100,000,
#     a matrix with one row per fleet-replacement fraction and columns
#     Mean, SD, P2.5, P25, P50, P75, P97.5.
#   k indexes the 7 pollutant-outcome pairs in fixed order:
#     1 PM2.5 Mortality, 2 PM2.5 Stroke, 3 PM2.5 Lung cancer, 4 PM2.5 Asthma,
#     5 NO2 Mortality, 6 NO2 Lung cancer, 7 NO2 Asthma.
#   The replacement-fraction grid is EV_diesel$replacement =
#     c(0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5); full diesel replacement lands
#     at ~0.10, which is the meaningful kink the report narrative hinges on.

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

source("Code/_theme_switchbox.R")
theme_set(theme_sb())

load("Results/Mitigation_EV_diesel_health_benefit.RData")

# Fixed pollutant-outcome order matching HB_EV_diesel_* list index k.
frac_grid <- c(0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5)
full_diesel_frac <- 0.1  # full diesel replacement / narrative kink

outcome_meta <- data.frame(
  k         = 1:7,
  pollutant = c(rep("PM2.5", 4), rep("NO2", 3)),
  outcome   = c("All-cause mortality", "Stroke", "Lung cancer", "Asthma",
                "All-cause mortality", "Lung cancer", "Asthma"),
  stringsAsFactors = FALSE
)

# Assemble a long data frame of region-level mean avoided rate (+ 95% UI) per
# fraction for the requested pollutant.
build_curve_df <- function(pollutant) {
  meta <- outcome_meta[outcome_meta$pollutant == pollutant, ]
  do.call(rbind, lapply(seq_len(nrow(meta)), function(i) {
    k   <- meta$k[i]
    mat <- HB_EV_diesel_rate_tot_stat[[k]]
    data.frame(
      frac    = frac_grid,
      rate    = mat[, "Mean"],
      lower   = mat[, "P2.5"],
      upper   = mat[, "P97.5"],
      outcome = meta$outcome[i],
      stringsAsFactors = FALSE
    )
  }))
}

# Build the composite chart. Legend ordered by magnitude at full diesel
# replacement (largest avoided rate first).
build_chart <- function(df, title) {
  order_at <- df[df$frac == full_diesel_frac, ]
  lvls <- order_at$outcome[order(order_at$rate, decreasing = TRUE)]
  df$outcome <- factor(df$outcome, levels = lvls)

  ggplot(df, aes(x = frac, y = rate, color = outcome)) +
    geom_vline(xintercept = full_diesel_frac,
               linetype = "dashed", color = "grey60", linewidth = 0.4) +
    annotate("text", x = full_diesel_frac, y = Inf,
             label = "Full diesel replacement (~10%)",
             hjust = -0.03, vjust = 1.6, size = 3, color = "grey40") +
    geom_line(aes(color = outcome), linewidth = 1) +
    scale_color_sb(name = "Health outcome") +
    scale_x_continuous(
      labels = scales::percent_format(accuracy = 1),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      x = "Share of fleet electrified (diesel-first)",
      y = "Avoided cases per 100,000 per year",
      title = title
    )
}

dir.create("Figures", showWarnings = FALSE)

# ---- PM2.5 composite (primary, 4 outcomes) ---------------------------------
df_pm <- build_curve_df("PM2.5")
p_pm  <- build_chart(df_pm, "Avoided health burden under diesel-first electrification (PM2.5)")
save_sb_fig(p_pm, "Figures/ExecSummary_diesel_avoided_cases_PM2.5.png")

# ---- NO2 companion (3 outcomes) --------------------------------------------
df_no2 <- build_curve_df("NO2")
p_no2  <- build_chart(df_no2, "Avoided health burden under diesel-first electrification (NO2)")
save_sb_fig(p_no2, "Figures/ExecSummary_diesel_avoided_cases_NO2.png")

# ---- Cross-check: values at full diesel replacement (~10%) -----------------
chk <- rbind(df_pm, df_no2)
chk <- chk[chk$frac == full_diesel_frac, c("outcome", "rate")]
cat("\nAvoided rate per 100,000 at full diesel replacement (frac = 0.10):\n")
print(chk, row.names = FALSE)
