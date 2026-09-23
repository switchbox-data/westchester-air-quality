# Executive-summary COMPOSITE charts for the report, in the Switchbox house style.
#
# Report Figure 6 replacement. Builds two independent electrification-scenario
# curves on one absolute x-axis: vehicle-miles electrified (in billions),
# converting each scenario's own replacement fraction into miles via that
# fuel's total fleet VMT (gasoline vs. diesel). One line per health outcome,
# solid = diesel scenario, dashed = gasoline scenario. Produces a PM2.5
# composite (4 outcomes; primary) and an NO2 companion (3 outcomes).
#
# Legends: the "Health outcome" color legend stays on the right; the "Fleet
# replaced" linetype legend sits at the bottom of the plot, horizontal, with
# its key lines drawn in neutral grey45 (per-guide legend positioning,
# ggplot2 >= 3.5.0) so it reads as a distinct, secondary legend rather than
# competing with the outcome colors.
#
# Each chart carries two grey endpoint labels, both computed from the plotted
# data (not hard-coded): "full diesel fleet" at the diesel line's endpoint
# (grows rightward, hjust = 0) and "full gasoline fleet" at the gasoline
# line's endpoint (grows leftward, hjust = 1, nudged up slightly to clear the
# line). The x scale gets a small right-side expansion so the diesel label
# is not clipped at the panel edge.
#
# Below the main panel (and its own vehicle-miles x-axis), each chart adds two
# thin "ruler" strips built with patchwork so a reader can convert absolute
# vehicle-miles to share-of-fleet at a glance: a SHORT diesel-fleet ruler that
# ends exactly at the diesel lines' endpoint (0/50/100% ticks), and a
# near-full-width gasoline-fleet ruler (0/25/50/75/100% ticks). All three rows
# (main panel, diesel ruler, gasoline ruler) share an identical x scale
# (limits, expansion) so the rulers' 100% ticks land on the same horizontal
# pixel position as each scenario's line endpoint in the panel above.
#
# The grey "Fleet replaced" linetype key is NOT drawn as a ggplot guide --
# patchwork's guides = "collect" merges per-guide positions ("right" for
# color, "bottom" for linetype) into a single shared guide area, which pulled
# the linetype key up next to the color legend instead of leaving it at the
# bottom. Instead it's hand-drawn (two short grey45 line segments + labels) as
# its own composition row, pinned below the gasoline ruler, so its position is
# deterministic. The color legend stays attached to the main plot alone (its
# own `position = "right"` guide), right of the main panel only.
#
# Reads ONLY cached results:
#   Results/Mitigation_EV_dieselonly_health_benefit.RData
#   Results/Mitigation_EV_gasonly_health_benefit.RData
# No recomputation, no Census API. Run inside the project container with Code,
# Data, Results and Figures bind-mounted, e.g.:
#   docker run --rm -v "$PWD/Code":/project/Code -v "$PWD/Results":/project/Results \
#     -v "$PWD/Figures":/project/Figures -w /project westchester-air-quality:latest \
#     Rscript Code/report_figures_composite_gasdiesel.R
#
# Data shape (see Code/04_Mitigation_strategies.R and its gas/diesel-only
# variants):
#   HB_EV_diesel_rate_tot_stat[[k]] / HB_EV_gas_rate_tot_stat[[k]] is the
#     region-level avoided RATE per 100,000, a matrix with one row per
#     fleet-replacement fraction and columns Mean, SD, P2.5, P25, P50, P75,
#     P97.5.
#   k indexes the 7 pollutant-outcome pairs in fixed order:
#     1 PM2.5 All-cause mortality, 2 PM2.5 Stroke, 3 PM2.5 Lung cancer,
#     4 PM2.5 Asthma, 5 NO2 All-cause mortality, 6 NO2 Lung cancer, 7 NO2 Asthma.
#   The replacement-fraction grid is EV_diesel$replacement == EV_gas$replacement
#     == c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
#     each fraction being of THAT fuel's own fleet (not the total fleet). The
#     two scenarios are independent — a diesel fleet ~10% of vehicle-miles vs.
#     a gasoline fleet ~89% — so converting each to absolute vehicle-miles
#     electrified (via each fuel's own total VMT) puts them on one comparable
#     axis instead of a shared percent axis.

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(patchwork)
})

source("Code/_theme_switchbox.R")
theme_set(theme_sb())

# grey45 does not collide with any brand color in sb_colors (sky/midnight/
# carrot/saffron/pistachio/black/white/midnight_text/pistachio_text) -- all of
# those are either saturated brand hues or pure black/white, none of them a
# neutral mid-grey.
legend_key_grey <- "grey45"

load("Results/Mitigation_EV_dieselonly_health_benefit.RData")
load("Results/Mitigation_EV_gasonly_health_benefit.RData")

# Total fleet VMT for each fuel; used to convert each scenario's replacement
# fraction into absolute vehicle-miles electrified. Units are THOUSANDS of
# miles: Ningrui's county source-type table (in millions of VMT) sums to
# exactly these values / 1e3, so dividing by 1e6 below yields billions.
VMT_GAS_TOTAL    <- 6310050.048
VMT_DIESEL_TOTAL <- 683544.208

# Fixed pollutant-outcome order matching HB_EV_*_rate_tot_stat list index k.
outcome_meta <- data.frame(
  k         = 1:7,
  pollutant = c(rep("PM2.5", 4), rep("NO2", 3)),
  outcome   = c("All-cause mortality", "Stroke", "Lung cancer", "Asthma",
                "All-cause mortality", "Lung cancer", "Asthma"),
  stringsAsFactors = FALSE
)

# Assemble a long data frame of region-level mean avoided rate (+ 95% UI) per
# fraction, for one pollutant and one scenario (diesel or gasoline).
build_curve_df <- function(pollutant, HB_list, EV_df, scenario_label) {
  meta <- outcome_meta[outcome_meta$pollutant == pollutant, ]
  do.call(rbind, lapply(seq_len(nrow(meta)), function(i) {
    k   <- meta$k[i]
    mat <- HB_list[[k]]
    data.frame(
      frac     = EV_df$replacement,
      rate     = mat[, "Mean"],
      lower    = mat[, "P2.5"],
      upper    = mat[, "P97.5"],
      outcome  = meta$outcome[i],
      scenario = scenario_label,
      stringsAsFactors = FALSE
    )
  }))
}

# Combine both scenarios (diesel + gasoline) into one long df for a pollutant.
build_pollutant_df <- function(pollutant) {
  rbind(
    build_curve_df(pollutant, HB_EV_diesel_rate_tot_stat, EV_diesel, "Diesel"),
    build_curve_df(pollutant, HB_EV_gas_rate_tot_stat, EV_gas, "Gasoline")
  )
}

# Build the composite chart on an absolute vehicle-miles-electrified axis.
# Color legend ordered by magnitude of the diesel scenario at full (100%)
# replacement of that fuel's fleet (largest first). Linetype legend moved to
# the bottom in neutral grey; both fleet endpoints get a data-driven grey
# label.
build_chart <- function(df, title) {
  order_at <- df[df$frac == 1 & df$scenario == "Diesel", ]
  lvls <- order_at$outcome[order(order_at$rate, decreasing = TRUE)]
  df$outcome  <- factor(df$outcome, levels = lvls)
  df$scenario <- factor(df$scenario, levels = c("Diesel", "Gasoline"))

  vmt_total_bn <- ifelse(df$scenario == "Diesel",
                          VMT_DIESEL_TOTAL / 1e6,
                          VMT_GAS_TOTAL / 1e6)
  df$vmt_bn <- df$frac * vmt_total_bn

  diesel_end_x <- VMT_DIESEL_TOTAL / 1e6
  diesel_end_y <- max(df$rate[df$scenario == "Diesel" & df$frac == 1])

  gas_end_x <- VMT_GAS_TOTAL / 1e6
  gas_end_y <- max(df$rate[df$scenario == "Gasoline" & df$frac == 1])

  ggplot(df, aes(x = vmt_bn, y = rate, color = outcome, linetype = scenario)) +
    geom_line(linewidth = 1) +
    scale_color_sb(name = "Health outcome",
                    guide = guide_legend(order = 1, position = "right")) +
    # The linetype scale still distinguishes Diesel (solid) vs. Gasoline
    # (dashed) lines, but draws NO ggplot guide here: patchwork's guides =
    # "collect" cannot be trusted to keep a per-guide bottom position separate
    # from the color legend's right position (it merges both into one guide
    # area). Instead, build_legend_strip() hand-draws this same "Fleet
    # replaced" key as its own composition row, pinned below the gasoline
    # ruler -- see build_composite().
    scale_linetype_manual(
      name   = "Fleet replaced",
      values = c(Diesel = "solid", Gasoline = "22"),
      guide  = "none"
    ) +
    scale_x_continuous(limits = c(0, gas_end_x),
                        expand = expansion(mult = c(0, 0.02))) +
    annotate("text", x = diesel_end_x, y = diesel_end_y,
             label = "full diesel fleet",
             hjust = 0, vjust = -0.4, size = 3, color = "grey40") +
    annotate("text", x = gas_end_x, y = gas_end_y,
             label = "full gasoline fleet",
             hjust = 1, vjust = -0.4, size = 3, color = "grey40") +
    labs(
      x = "Vehicle-miles electrified (billions)",
      y = "Avoided cases per 100,000 per year",
      title = title
    ) -> p

  list(plot = p, diesel_end_x = diesel_end_x, gas_end_x = gas_end_x)
}

# Build one thin "ruler" strip: a spine from x = 0 to fleet_max (data units,
# billions of VMT), with perpendicular ticks + percent labels at ticks_frac *
# fleet_max, and a small fleet-name label. panel_max sets the x scale's upper
# limit and MUST match the main panel's (and the other ruler's) so the two
# rulers and the main panel line up pixel-for-pixel on x. label_side controls
# whether the fleet-name label grows away from ("right") or back along
# ("left") the spine from its right end, so it never gets cramped or clipped.
build_ruler <- function(fleet_max, panel_max, ticks_frac, fleet_label,
                         label_side = c("right", "left")) {
  label_side <- match.arg(label_side)
  tick_df <- data.frame(
    x     = ticks_frac * fleet_max,
    label = paste0(round(ticks_frac * 100), "%")
  )
  # Dodge the endpoint tick labels off the ends of the ruler (0% grows left,
  # 100% grows right, interior ticks stay centered) so they don't collide on
  # the short diesel ruler, where three centered labels would overlap.
  n <- nrow(tick_df)
  tick_df$hjust <- 0.5
  tick_df$hjust[1] <- 1
  tick_df$hjust[n] <- 0

  label_hjust <- if (label_side == "right") 0 else 1
  label_x     <- if (label_side == "right") fleet_max + panel_max * 0.01 else fleet_max

  ggplot(tick_df) +
    annotate("segment", x = 0, xend = fleet_max, y = 0, yend = 0,
             linewidth = 0.5, color = "black") +
    geom_segment(aes(x = x, xend = x, y = -0.3, yend = 0.3),
                 linewidth = 0.5, color = "black") +
    geom_text(aes(x = x, y = -1.1, label = label, hjust = hjust),
              family = "IBM-Plex-Sans", size = 3.5, color = "black") +
    annotate("text", x = label_x, y = 1.1, label = fleet_label,
             hjust = label_hjust, vjust = 0, family = "IBM-Plex-Sans",
             size = 3.5, fontface = "bold", color = "grey30") +
    scale_x_continuous(limits = c(0, panel_max),
                        expand = expansion(mult = c(0, 0.02))) +
    scale_y_continuous(limits = c(-2, 2)) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(text = element_text(family = "IBM-Plex-Sans", size = 12))
}

# Hand-drawn "Fleet replaced" legend key: two short line segments (solid /
# dashed) in neutral grey45 with labels, centered under the panel width. Built
# by hand -- not a ggplot guide -- so it renders as its own composition row,
# below both rulers, regardless of how patchwork's guide collection would
# otherwise merge it with the color legend (see note in build_chart()).
build_legend_strip <- function() {
  key_df <- data.frame(
    scenario = factor(c("Diesel", "Gasoline"), levels = c("Diesel", "Gasoline")),
    x0       = c(0.36, 0.57),
    x1       = c(0.45, 0.66)
  )
  ggplot(key_df) +
    geom_segment(aes(x = x0, xend = x1, y = 0, yend = 0, linetype = scenario),
                 linewidth = 0.7, color = legend_key_grey) +
    geom_text(aes(x = x1 + 0.015, y = 0, label = scenario),
              hjust = 0, vjust = 0.5, family = "IBM-Plex-Sans", size = 3.8,
              color = "grey30") +
    annotate("text", x = 0.34, y = 0, label = "Fleet replaced:",
             hjust = 1, vjust = 0.5, family = "IBM-Plex-Sans", size = 3.8,
             fontface = "bold", color = "grey30") +
    scale_linetype_manual(values = c(Diesel = "solid", Gasoline = "22"), guide = "none") +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(-1, 1)) +
    coord_cartesian(clip = "off") +
    theme_void()
}

# Stack the main panel over the two fleet-share rulers and the hand-drawn
# "Fleet replaced" key (patchwork). The main panel, diesel ruler, and gasoline
# ruler share an identical x scale so ticks align with line endpoints above.
# No guide collection is used: the color legend stays attached to the main
# plot alone (its own `position = "right"` guide, right of row 1 only), and
# the linetype key is plain geoms in the bottom row -- so nothing can merge
# the two into one guide area.
build_composite <- function(p_main, diesel_end_x, gas_end_x) {
  p_diesel_ruler <- build_ruler(
    fleet_max = diesel_end_x, panel_max = gas_end_x,
    ticks_frac = c(0, 0.5, 1), fleet_label = "Diesel fleet",
    label_side = "right"
  )
  p_gas_ruler <- build_ruler(
    fleet_max = gas_end_x, panel_max = gas_end_x,
    ticks_frac = c(0, 0.25, 0.5, 0.75, 1), fleet_label = "Gasoline fleet",
    label_side = "left"
  )
  p_legend <- build_legend_strip()
  (p_main / p_diesel_ruler / p_gas_ruler / p_legend) +
    plot_layout(heights = c(1, 0.10, 0.10, 0.12))
}

dir.create("Figures", showWarnings = FALSE)

# ---- PM2.5 composite (primary, 4 outcomes) ---------------------------------
df_pm    <- build_pollutant_df("PM2.5")
ch_pm    <- build_chart(df_pm, "Avoided health burden per vehicle-mile electrified: diesel vs. gasoline (PM2.5)")
comp_pm  <- build_composite(ch_pm$plot, ch_pm$diesel_end_x, ch_pm$gas_end_x)
save_sb_fig(comp_pm, "Figures/ExecSummary_gasdiesel_avoided_cases_PM2.5.png", width = 9, height = 5.6)
save_sb_fig(comp_pm, "Figures/ExecSummary_gasdiesel_avoided_cases_PM2.5.svg", width = 9, height = 5.6)

# ---- NO2 companion (3 outcomes) --------------------------------------------
df_no2   <- build_pollutant_df("NO2")
ch_no2   <- build_chart(df_no2, "Avoided health burden per vehicle-mile electrified: diesel vs. gasoline (NO2)")
comp_no2 <- build_composite(ch_no2$plot, ch_no2$diesel_end_x, ch_no2$gas_end_x)
save_sb_fig(comp_no2, "Figures/ExecSummary_gasdiesel_avoided_cases_NO2.png", width = 9, height = 5.6)
save_sb_fig(comp_no2, "Figures/ExecSummary_gasdiesel_avoided_cases_NO2.svg", width = 9, height = 5.6)

# ---- Cross-check: values at 100% replacement of each fuel's own fleet ------
chk <- rbind(df_pm, df_no2)
chk <- chk[chk$frac == 1, c("outcome", "scenario", "rate")]
chk <- chk[order(chk$outcome, chk$scenario), ]
cat("\nAvoided rate per 100,000 at 100% replacement of that fuel's fleet:\n")
print(chk, row.names = FALSE)
