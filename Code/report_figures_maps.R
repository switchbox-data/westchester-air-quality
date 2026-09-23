# Re-plot the report's tract-level choropleth maps from cached Results in the
# Switchbox house style (no recomputation, no Census API, no NMF refit).

source("Code/_theme_switchbox.R")

suppressPackageStartupMessages({
  library(sf)
  library(ggplot2)
  library(ggnewscale)
  library(scales)
  library(patchwork)
})

# Per-panel canvas (inches), held fixed across every combined figure so all
# panels are identical in size: 1x2 combos = 2*panel_w x 1*panel_h, the 2x2 =
# 2*panel_w x 2*panel_h (exactly twice the height of each 1x2 combo).
panel_w <- 5
panel_h <- 4


## Shared basemap: every tract map covers the same study-area bbox, so fetch the
## greyscale tiles ONCE and reuse. Holding a single copy (rather than one per
## figure) keeps memory bounded now that combined figures retain several panels.
## zoom 13 is sharp enough for the 3000px-wide individual exports while small
## enough that the panels we keep around fit comfortably in container memory.
load("Results/Exposure_PM2.5_tract.RData")  # PM_tract_annual_s_w_h
shared_bm <- sb_basemap(sf::st_as_sfc(sf::st_bbox(PM_tract_annual_s_w_h)), zoom = 13)
rm(PM_tract_annual_s_w_h)


## Shared map builder: greyscale basemap (or white fallback) + theme ----------
## `panel_title` adds a short heading for use inside combined multi-panel figs.
sb_choropleth <- function(tr, name, legend_labels = waiver(), panel_title = NULL) {
  bm <- shared_bm

  p <- ggplot()
  if (!is.null(bm)) {
    p <- p +
      geom_raster(data = bm, aes(x = x, y = y, fill = hex)) +
      scale_fill_identity() +
      ggnewscale::new_scale_fill()
  }
  p <- p +
    geom_sf(data = tr, aes(fill = Mean),
            color = "black", linewidth = 0.4, alpha = 0.9) +
    scale_fill_sb_map(name = name, labels = legend_labels) +
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


## 1. Annual average exposure concentration (PM2.5, NO2) --------------------
## Source objects are already sf with the value column merged onto geometry.
exposure_panels <- list()
for (pollut_std in c("PM2.5", "NO2")) {
  load(paste0("Results/Exposure_", pollut_std, "_tract.RData"))  # PM_tract_annual_s_w_h
  tr <- PM_tract_annual_s_w_h
  tr$Mean <- tr$Wt_SWH
  legend_title <- ifelse(pollut_std == "PM2.5", "μg/m³", "ppb")
  p <- sb_choropleth(tr, name = legend_title)
  save_sb_fig(p, paste0("Figures/Exposure_", pollut_std, "_tract.png"),
              width = 10, height = 8)
  exposure_panels[[pollut_std]] <- sb_choropleth(
    tr, name = legend_title, panel_title = pollut_std
  )
  rm(PM_tract_annual_s_w_h, tr)
}


## 2. Traffic source contribution, percentage (PM2.5, NO2) ------------------
## Source objects are already sf; Fac1_perc is the traffic fraction.
traffic_panels <- list()
for (pollut_std in c("PM2.5", "NO2")) {
  load(paste0("Results/Exposure_source_", pollut_std, "_tract.RData"))  # PM_tract_annual_s_w_h
  tr <- PM_tract_annual_s_w_h
  tr$Mean <- tr$Fac1_perc
  p <- sb_choropleth(tr, name = "Percentage",
                     legend_labels = scales::percent_format(accuracy = 1))
  save_sb_fig(p, paste0("Figures/Exposure_traffic_", pollut_std, "_tract_perc.png"),
              width = 10, height = 8)
  traffic_panels[[pollut_std]] <- sb_choropleth(
    tr, name = "Percentage",
    legend_labels = scales::percent_format(accuracy = 1),
    panel_title = pollut_std
  )
  rm(PM_tract_annual_s_w_h, tr)
}


## 3. Attributable health-risk rate per 100,000 (7 pollutant/outcome pairs) --
## risk_rate_stat rows align positionally with the GEOID-sorted shapefile,
## matching the original cbind(WEST_tract, Mean = ...) join in 02_Health_risk.R.
WEST_tract <- st_read("Data/Westchester_tracts_selected.shp", quiet = TRUE)
WEST_tract <- WEST_tract[order(WEST_tract$GEOID), ]

load("Results/Health_risk_rate_tract.RData")  # risk_rate_stat (list of 7)

pollutant_list <- c(rep("PM2.5", 4), rep("NO2", 3))
outcome_list   <- c("Mortality", "Stroke", "Lung_cancer", "Asthma",
                    "Mortality", "Lung_cancer", "Asthma")

# Panel titles for the 2x2 health composite (mortality + asthma only).
health_panels <- list()
health_panel_titles <- list(
  "PM2.5_Mortality" = "PM2.5 — premature mortality",
  "NO2_Mortality"   = "NO2 — premature mortality",
  "PM2.5_Asthma"    = "PM2.5 — childhood asthma",
  "NO2_Asthma"      = "NO2 — childhood asthma"
)

for (k in seq_along(pollutant_list)) {
  tr <- cbind(WEST_tract, Mean = risk_rate_stat[[k]][, "Mean"])
  p <- sb_choropleth(tr, name = "Per 100,000")
  save_sb_fig(p, paste0("Figures/Health_risk_rate_",
                        pollutant_list[k], "_", outcome_list[k], "_tract.png"),
              width = 10, height = 8)
  key <- paste0(pollutant_list[k], "_", outcome_list[k])
  if (!is.null(health_panel_titles[[key]])) {
    health_panels[[key]] <- sb_choropleth(
      tr, name = "Per 100,000", panel_title = health_panel_titles[[key]]
    )
  }
}


## 4. Combined multi-panel report figures (patchwork) -----------------------
## Every panel uses the same per-panel canvas, so panel sizes are identical
## across all three combined figures.

# Combined exposure (Fig 3+4): PM2.5 | NO2, 1 row x 2 cols.
combined_exposure <- exposure_panels[["PM2.5"]] + exposure_panels[["NO2"]]
save_sb_fig(combined_exposure, "Figures/Combined_exposure_tract.png",
            width = 2 * panel_w, height = panel_h)

# Combined traffic (Fig 5+6): traffic%->PM2.5 | traffic%->NO2, 1 row x 2 cols.
combined_traffic <- traffic_panels[["PM2.5"]] + traffic_panels[["NO2"]]
save_sb_fig(combined_traffic, "Figures/Combined_traffic_tract_perc.png",
            width = 2 * panel_w, height = panel_h)

# Combined health 2x2 (Fig 7a/7b/8a/8b):
#   row 1 = premature mortality (PM2.5 | NO2)
#   row 2 = childhood asthma    (PM2.5 | NO2)
combined_health <-
  (health_panels[["PM2.5_Mortality"]] | health_panels[["NO2_Mortality"]]) /
  (health_panels[["PM2.5_Asthma"]]    | health_panels[["NO2_Asthma"]])
save_sb_fig(combined_health,
            "Figures/Combined_health_rate_mortality_asthma_2x2.png",
            width = 2 * panel_w, height = 2 * panel_h)
