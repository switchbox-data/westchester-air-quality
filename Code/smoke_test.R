# Minimal smoke test: verifies the environment can load every package
# referenced by the project scripts and can read the input data files.
# Skips the full pipeline (1000-iteration bootstrap, tidycensus network calls,
# NMF factorization) — those are integration-level workloads.

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(reshape2)
  library(lubridate)
  library(clipr)
  library(maptiles)
  library(RColorBrewer)
  library(ggplot2)
  library(ggnewscale)
  library(splines)
  library(DescTools)
  library(tidyr)
  library(tigris)
  library(tidycensus)
  library(NMF)
})

cat("== package load: OK ==\n")
cat("R version:", R.version.string, "\n\n")

# Shapefile reads (small, fast)
seg <- st_read("Data/Mean_concentration_road_segments.shp", quiet = TRUE)
colnames(seg)[1:10] <- c("segment_id","blackcarbon","c2h6","ch4","co",
                        "co2","no","no2","o3","pm_2.5")
cat("road segments:", nrow(seg), "rows\n")
cat("  PM2.5 mean:", round(mean(seg$pm_2.5, na.rm = TRUE), 3), "\n")
cat("  NO2   mean:", round(mean(seg$no2,   na.rm = TRUE), 3), "\n\n")

tracts <- st_read("Data/Westchester_tracts_selected.shp", quiet = TRUE)
cat("census tracts:", nrow(tracts), "rows\n\n")

# Mobile-monitoring RData (the big one — ~164 MB on disk)
load("Data/Mobile_monitoring_1second_data_Aclima.RData")
stopifnot(exists("dat"))
cat("mobile monitoring rows:", nrow(dat), "\n")
cat("modalities:", paste(unique(as.character(dat$modality)), collapse = ", "), "\n\n")

# Tiny slice of the actual exposure-assessment workflow: per-block-group
# hourly median for PM2.5 only, with no bootstrap.
pm <- dat %>%
  filter(modality == "pm_2.5", quality_flag == 0) %>%
  mutate(timestamp = with_tz(ymd_hms(timestamp, tz = "UTC"), "America/New_York"))

pm_sf <- st_as_sf(pm, coords = c("lon","lat"), crs = 4326, remove = FALSE) |>
  st_transform(crs = 32115)

tracts_proj <- st_transform(tracts, 32115)
nearest <- st_nearest_feature(pm_sf, tracts_proj)
pm_sf$tract <- tracts_proj$GEOID[nearest]

per_tract <- st_drop_geometry(pm_sf) |>
  mutate(hour = floor_date(timestamp, "hour")) |>
  group_by(tract, hour) |>
  summarise(median = pmax(median(value, na.rm = TRUE), 0),
            .groups = "drop") |>
  group_by(tract) |>
  summarise(annual_pm25 = mean(median, na.rm = TRUE), .groups = "drop") |>
  arrange(desc(annual_pm25))

cat("== PM2.5 tract-level annual avg (top 5) ==\n")
print(head(per_tract, 5))

cat("\n== smoke test PASSED ==\n")

# Persist a small artifact so a host bind-mount can verify Results/ writes work.
saveRDS(per_tract, file = "Results/smoke_test_pm25_by_tract.rds")
cat("wrote Results/smoke_test_pm25_by_tract.rds\n")
