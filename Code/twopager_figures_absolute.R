# Absolute-number exhibits for the advocacy 2-pager.
#
# Two figures, both in ABSOLUTE cases per year (not per 100,000), with PM2.5 and
# NO2 combined per outcome so the exhibits match the 2-pager's stat callouts
# (~230 attributable deaths/yr, ~280 child asthma cases/yr; ~80 deaths and ~90
# asthma cases avoided under full diesel electrification):
#   Figures/TwoPager_health_number_map.png        deaths | child asthma, by tract
#   Figures/TwoPager_gasdiesel_avoided_cases_absolute.png
#                                                 avoided cases vs. vehicle-miles
#
# Combination rule: some cases are attributable to both pollutants, so combined
# counts assume independent effects, PAF_comb = 1 - (1-PAF_pm)(1-PAF_no2),
# evaluated per Monte Carlo draw (the same rule as the 2-pager footnotes).
# Stroke has no NO2 concentration-response function, so its combined value is
# the PM2.5 value. For the avoided-cases chart the burden draws (Code/02) and
# mitigation draws (Code/05) come from separate simulation runs, so pairing
# them draw-by-draw is an independence approximation; only the plotted means
# are used here.
#
# Reads ONLY cached results + Data. Run inside the project container:
#   docker run --rm -v "$PWD/Code":/project/Code -v "$PWD/Data":/project/Data \
#     -v "$PWD/Results":/project/Results -v "$PWD/Figures":/project/Figures \
#     -w /project westchester-air-quality:latest Rscript Code/twopager_figures_absolute.R

source("Code/_theme_switchbox.R")

suppressPackageStartupMessages({
  library(sf)
  library(ggplot2)
  library(ggnewscale)
  library(scales)
  library(patchwork)
})

# Baseline incidence rates per 100,000 of the at-risk age group, from
# Code/02_Health_risk.R baseline_data (CDC Wonder / GBD 2023).
RATE_MORT   <- 941.22   # all-cause mortality, 18+
RATE_ASTHMA <- 1233.65  # asthma incidence, 0-20
RATE_LC     <- 135.70   # lung cancer incidence, 50+

# attr_* are attributable-count draws (matrix, units x draws); N is the
# baseline case count (length = rows). Matrix/vector ops recycle down columns,
# so element [i, j] pairs with N[i].
combine_attr <- function(attr_pm, attr_no2, N) {
  N * (1 - (1 - attr_pm / N) * (1 - attr_no2 / N))
}
combine_avoided <- function(attr_pm, attr_no2, av_pm, av_no2, N) {
  N * ((1 - (attr_pm - av_pm) / N) * (1 - (attr_no2 - av_no2) / N) -
       (1 - attr_pm / N) * (1 - attr_no2 / N))
}


## 1. Tract map: absolute attributable deaths and child asthma cases ----------

WEST_tract <- st_read("Data/Westchester_tracts_selected.shp", quiet = TRUE)
WEST_tract <- WEST_tract[order(WEST_tract$GEOID), ]

load("Data/Westchester_tract_population_by_age.RData")  # WEST_tract_popu
pop <- WEST_tract_popu[order(WEST_tract_popu$GEOID), ]
stopifnot(identical(WEST_tract$GEOID, pop$GEOID))

# risk_number rows align positionally with the GEOID-sorted shapefile, same as
# risk_rate_stat in Code/report_figures_maps.R. List index k: 1 PM2.5 mortality,
# 2 PM2.5 stroke, 3 PM2.5 lung cancer, 4 PM2.5 asthma, 5 NO2 mortality,
# 6 NO2 lung cancer, 7 NO2 asthma.
load("Results/Health_risk_number_tract.RData")  # risk_number (draws), risk_number_stat

N_mort_t   <- RATE_MORT   * pop$Popu_over18 / 1e5
N_asthma_t <- RATE_ASTHMA * pop$Popu_0_20   / 1e5

mort_comb   <- combine_attr(as.matrix(risk_number[[1]]), as.matrix(risk_number[[5]]), N_mort_t)
asthma_comb <- combine_attr(as.matrix(risk_number[[4]]), as.matrix(risk_number[[7]]), N_asthma_t)

cat(sprintf("Map totals (mean/yr): deaths %.1f, child asthma %.1f\n",
            sum(rowMeans(mort_comb)), sum(rowMeans(asthma_comb))))

load("Results/Exposure_PM2.5_tract.RData")  # PM_tract_annual_s_w_h, for the bbox
shared_bm <- sb_basemap(sf::st_as_sfc(sf::st_bbox(PM_tract_annual_s_w_h)), zoom = 13)
rm(PM_tract_annual_s_w_h)

# Same builder as Code/report_figures_maps.R sb_choropleth.
sb_choropleth <- function(tr, name, panel_title = NULL) {
  p <- ggplot()
  if (!is.null(shared_bm)) {
    p <- p +
      geom_raster(data = shared_bm, aes(x = x, y = y, fill = hex)) +
      scale_fill_identity() +
      ggnewscale::new_scale_fill()
  }
  p <- p +
    geom_sf(data = tr, aes(fill = Mean),
            color = "black", linewidth = 0.4, alpha = 0.9) +
    scale_fill_sb_map(name = name) +
    theme_void(base_family = "IBM-Plex-Sans") +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "right",
      legend.title = element_text(size = 12, hjust = 0.5),
      legend.text = element_text(size = 12)
    )
  if (!is.null(panel_title)) {
    p <- p +
      ggtitle(panel_title) +
      theme(plot.title = element_text(
        family = "IBM-Plex-Sans", size = 13,
        colour = sb_colors[["midnight"]], hjust = 0.5,
        margin = margin(b = 4)
      ))
  }
  p
}

panel_w <- 5
panel_h <- 4

p_mort <- sb_choropleth(cbind(WEST_tract, Mean = rowMeans(mort_comb)),
                        name = "Deaths\nper year",
                        panel_title = "Premature deaths")
p_asthma <- sb_choropleth(cbind(WEST_tract, Mean = rowMeans(asthma_comb)),
                          name = "Cases\nper year",
                          panel_title = "New childhood asthma cases")

save_sb_fig(p_mort + p_asthma, "Figures/TwoPager_health_number_map.png",
            width = 2 * panel_w, height = panel_h)


## 2. Avoided-cases chart: absolute counts, PM2.5 + NO2 combined --------------

load("Results/Health_risk_number_and_rate_region.RData")     # risk_number_tot (7 x 1000)
load("Results/Mitigation_EV_dieselonly_health_benefit.RData") # HB_EV_diesel_number_tot, EV_diesel
load("Results/Mitigation_EV_gasonly_health_benefit.RData")    # HB_EV_gas_number_tot, EV_gas

N_MORT   <- RATE_MORT   * sum(pop$Popu_over18) / 1e5
N_ASTHMA <- RATE_ASTHMA * sum(pop$Popu_0_20)   / 1e5
N_LC     <- RATE_LC     * sum(pop$Popu_over50) / 1e5

# Mean combined avoided cases per replacement fraction (rows of the 14-step
# grid). k_pm/k_no2 index the pollutant-outcome pairs; k_no2 = NA -> PM2.5 only.
combined_curve <- function(HB_list, k_pm, k_no2, N) {
  if (is.na(k_no2)) return(rowMeans(HB_list[[k_pm]]))
  attr_pm  <- matrix(risk_number_tot[k_pm, ],  nrow = 14, ncol = 1000, byrow = TRUE)
  attr_no2 <- matrix(risk_number_tot[k_no2, ], nrow = 14, ncol = 1000, byrow = TRUE)
  rowMeans(combine_avoided(attr_pm, attr_no2, HB_list[[k_pm]], HB_list[[k_no2]], N))
}

outcome_meta <- data.frame(
  outcome = c("All-cause mortality", "Asthma", "Lung cancer", "Stroke"),
  k_pm    = c(1, 4, 3, 2),
  k_no2   = c(5, 7, 6, NA),
  N       = c(N_MORT, N_ASTHMA, N_LC, NA),  # N unused for PM2.5-only outcomes
  stringsAsFactors = FALSE
)

build_scenario_df <- function(HB_list, EV_df, scenario_label) {
  do.call(rbind, lapply(seq_len(nrow(outcome_meta)), function(i) {
    m <- outcome_meta[i, ]
    data.frame(
      frac     = EV_df$replacement,
      cases    = combined_curve(HB_list, m$k_pm, m$k_no2, m$N),
      outcome  = m$outcome,
      scenario = scenario_label,
      stringsAsFactors = FALSE
    )
  }))
}

df <- rbind(
  build_scenario_df(HB_EV_diesel_number_tot, EV_diesel, "Diesel"),
  build_scenario_df(HB_EV_gas_number_tot, EV_gas, "Gasoline")
)

# Chart + ruler machinery adapted from Code/report_figures_composite_gasdiesel.R
# (same layout; y is now absolute cases per year).
theme_set(theme_sb())
legend_key_grey <- "grey45"

VMT_GAS_TOTAL    <- 6310050.048  # thousands of miles (see gasdiesel script)
VMT_DIESEL_TOTAL <- 683544.208

build_chart <- function(df, title) {
  order_at <- df[df$frac == 1 & df$scenario == "Diesel", ]
  lvls <- order_at$outcome[order(order_at$cases, decreasing = TRUE)]
  df$outcome  <- factor(df$outcome, levels = lvls)
  df$scenario <- factor(df$scenario, levels = c("Diesel", "Gasoline"))

  vmt_total_bn <- ifelse(df$scenario == "Diesel",
                          VMT_DIESEL_TOTAL / 1e6,
                          VMT_GAS_TOTAL / 1e6)
  df$vmt_bn <- df$frac * vmt_total_bn

  diesel_end_x <- VMT_DIESEL_TOTAL / 1e6
  diesel_end_y <- max(df$cases[df$scenario == "Diesel" & df$frac == 1])
  gas_end_x <- VMT_GAS_TOTAL / 1e6
  gas_end_y <- max(df$cases[df$scenario == "Gasoline" & df$frac == 1])

  p <- ggplot(df, aes(x = vmt_bn, y = cases, color = outcome, linetype = scenario)) +
    geom_line(linewidth = 1) +
    scale_color_sb(name = "Health outcome",
                    guide = guide_legend(order = 1, position = "right")) +
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
      y = "Avoided cases per year",
      title = title
    )
  list(plot = p, diesel_end_x = diesel_end_x, gas_end_x = gas_end_x)
}

build_ruler <- function(fleet_max, panel_max, ticks_frac, fleet_label,
                         label_side = c("right", "left")) {
  label_side <- match.arg(label_side)
  tick_df <- data.frame(
    x     = ticks_frac * fleet_max,
    label = paste0(round(ticks_frac * 100), "%")
  )
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
  (p_main / p_diesel_ruler / p_gas_ruler / build_legend_strip()) +
    plot_layout(heights = c(1, 0.10, 0.10, 0.12))
}

ch <- build_chart(df, "Avoided cases per year: electrifying diesel vs. gasoline (PM2.5 + NO2)")
comp <- build_composite(ch$plot, ch$diesel_end_x, ch$gas_end_x)
save_sb_fig(comp, "Figures/TwoPager_gasdiesel_avoided_cases_absolute.png",
            width = 9, height = 5.6)

# Cross-check against the 2-pager callouts (~80 deaths, ~90 asthma at full diesel).
chk <- df[df$frac == 1, ]
cat("\nCombined avoided cases per year at 100% replacement of that fuel's fleet:\n")
print(chk[order(chk$scenario, -chk$cases), c("scenario", "outcome", "cases")],
      row.names = FALSE)
