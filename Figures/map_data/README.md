# Tract-level map data (designer handoff)

One CSV per choropleth map shown in the report. Each row is one census tract.
Join to the tract geometries in `Data/Westchester_tracts_selected.shp` on the
11-digit 2020 Census tract **GEOID** (the shapefile's `GEOID` column). 70 tracts
each. Values are exactly what the corresponding published map draws — generated
from the same cached results by `Code/report_map_data_csv.R`.

| CSV | Value column | Meaning | Report figure |
|---|---|---|---|
| `exposure_PM2.5_tract.csv` | `exposure_ug_m3` | Annual avg PM₂.₅ (µg/m³) | Combined exposure — PM₂.₅ panel |
| `exposure_NO2_tract.csv` | `exposure_ppb` | Annual avg NO₂ (ppb) | Combined exposure — NO₂ panel |
| `traffic_perc_PM2.5_tract.csv` | `traffic_fraction` | Traffic source share of PM₂.₅, **fraction 0–1** (map legend shows as %) | Combined traffic — PM₂.₅ panel |
| `traffic_perc_NO2_tract.csv` | `traffic_fraction` | Traffic source share of NO₂, **fraction 0–1** (map legend shows as %) | Combined traffic — NO₂ panel |
| `health_rate_PM2.5_Mortality_tract.csv` | `rate_per_100k` | Attributable premature-mortality rate per 100,000 | Health 2×2 — top-left |
| `health_rate_NO2_Mortality_tract.csv` | `rate_per_100k` | Attributable premature-mortality rate per 100,000 | Health 2×2 — top-right |
| `health_rate_PM2.5_Asthma_tract.csv` | `rate_per_100k` | Attributable childhood-asthma rate per 100,000 | Health 2×2 — bottom-left |
| `health_rate_NO2_Asthma_tract.csv` | `rate_per_100k` | Attributable childhood-asthma rate per 100,000 | Health 2×2 — bottom-right |

Notes:
- **Traffic columns are fractions (0–1)**, not percentages — multiply by 100 for
  the legend the report uses (e.g. `0.7574` → `76%`).
- The report's Combined-health figure shows only these 4 outcome panels. The
  pipeline also produces stroke and lung-cancer rate maps (PM₂.₅ stroke, lung
  cancer; NO₂ lung cancer); rerun `Code/report_map_data_csv.R` with those
  indices added if those are wanted too.
- House-style ramp used in the report is white → sky → midnight (sequential).
