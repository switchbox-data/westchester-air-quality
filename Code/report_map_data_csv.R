# Export one CSV per tract-level choropleth map shown in the report, using the
# SAME cached Results objects and joins as Code/report_figures_maps.R, so each
# CSV's values are exactly what the corresponding map draws. No recomputation.
#
# Output: Figures/map_data/<map>.csv, each with columns GEOID, <value>.
# For handoff to a designer rendering the maps from the tract shapefile
# (Data/Westchester_tracts_selected.shp), joined on GEOID.

suppressPackageStartupMessages({
  library(sf)
})

out_dir <- "Figures/map_data"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# GEOID as a plain character string (drop geometry, keep the join key + value).
write_map_csv <- function(geoid, value, value_name, file) {
  df <- data.frame(GEOID = as.character(geoid), stringsAsFactors = FALSE)
  df[[value_name]] <- value
  path <- file.path(out_dir, file)
  write.csv(df, path, row.names = FALSE)
  cat(sprintf("wrote %-52s %d tracts\n", file, nrow(df)))
}


## 1. Annual average exposure concentration (Fig: Combined_exposure_tract) ----
##    value = Wt_SWH (PM2.5 in ug/m3, NO2 in ppb)
for (pollut_std in c("PM2.5", "NO2")) {
  load(paste0("Results/Exposure_", pollut_std, "_tract.RData"))  # PM_tract_annual_s_w_h
  tr <- PM_tract_annual_s_w_h
  value_name <- if (pollut_std == "PM2.5") "exposure_ug_m3" else "exposure_ppb"
  write_map_csv(tr$tract, tr$Wt_SWH, value_name,
                paste0("exposure_", pollut_std, "_tract.csv"))
  rm(PM_tract_annual_s_w_h, tr)
}


## 2. Traffic source contribution, percentage (Fig: Combined_traffic_tract) --
##    value = Fac1_perc (fraction 0-1; the map's legend shows it as a percent)
for (pollut_std in c("PM2.5", "NO2")) {
  load(paste0("Results/Exposure_source_", pollut_std, "_tract.RData"))  # PM_tract_annual_s_w_h
  tr <- PM_tract_annual_s_w_h
  write_map_csv(tr$tract, tr$Fac1_perc, "traffic_fraction",
                paste0("traffic_perc_", pollut_std, "_tract.csv"))
  rm(PM_tract_annual_s_w_h, tr)
}


## 3. Attributable health-risk rate per 100,000 (Fig: Combined_health 2x2) ---
##    risk_rate_stat rows align positionally with the GEOID-sorted shapefile,
##    matching cbind(WEST_tract, Mean = ...) in report_figures_maps.R.
##    The 2x2 report figure shows exactly these 4 pollutant/outcome panels.
WEST_tract <- st_read("Data/Westchester_tracts_selected.shp", quiet = TRUE)
WEST_tract <- WEST_tract[order(WEST_tract$GEOID), ]

load("Results/Health_risk_rate_tract.RData")  # risk_rate_stat (list of 7)

# Positional index into the fixed 7-pair sequence used throughout the pipeline:
#   1 PM2.5 Mortality, 2 PM2.5 Stroke, 3 PM2.5 Lung_cancer, 4 PM2.5 Asthma,
#   5 NO2 Mortality,   6 NO2 Lung_cancer, 7 NO2 Asthma
panels_in_report <- list(
  list(k = 1, file = "health_rate_PM2.5_Mortality_tract.csv"),
  list(k = 5, file = "health_rate_NO2_Mortality_tract.csv"),
  list(k = 4, file = "health_rate_PM2.5_Asthma_tract.csv"),
  list(k = 7, file = "health_rate_NO2_Asthma_tract.csv")
)

for (p in panels_in_report) {
  write_map_csv(WEST_tract$GEOID, risk_rate_stat[[p$k]][, "Mean"],
                "rate_per_100k", p$file)
}

cat("\nDone. CSVs in", out_dir, "\n")
