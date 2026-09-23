# Westchester Air Pollution Health Impact Assessment

This repository has the code for the Westchester Air Pollution Health Impact Assessment. The study uses one year of street-level air monitoring data from Aclima, Inc. to estimate:

1. the annual exposure to PM2.5 and NO₂ in each census tract of lower Westchester County, New York;
2. the deaths and illnesses that this exposure causes;
3. the part of the exposure that comes from traffic;
4. the health benefits of traffic measures, for example vehicle electrification and a low-emission zone.

Switchbox, the Energy Justice Law and Policy Center (EJLPC) and the University of Washington (UW) prepared the study. Dr. Magali Blanco and Dr. Ningrui Liu of the UW Department of Environmental and Occupational Health Sciences led the technical analysis. The full method will be in a peer-reviewed paper from the UW team.

## Quickstart

You need Docker, `curl` and a free [Census API key](https://api.census.gov/data/key_signup.html).

```sh
git clone https://github.com/switchbox-data/westchester-air-quality.git
cd westchester-air-quality
make data                         # download the input data (157 MB)
export CENSUS_API_KEY=<your-key>
make                              # build the Docker image, then run stages 01–04
```

The full run is slow, because the bootstrap and NMF steps take a long time. Stage 01 needs about 6 GB of memory. To run one stage only, use `make 01`, `make 02`, `make 03` or `make 04`. The `Makefile` has more options.

You do not need to run the pipeline to see the results. `Results/` and `Figures/` have the outputs of the last run, except for three files that are too large for GitHub (see `.gitignore`).

## Data

### Where the data is

The Aclima data is **public domain**. It is in the public S3 bucket `s3://switchbox-ny-air-monitoring` (AWS region `us-west-2`). You do not need an AWS account to get it.

`scripts/fetch_data.sh` downloads the files and checks each file against the checksums in `scripts/data_manifest.sha256`.

```sh
scripts/fetch_data.sh                  # Westchester input only, to Data/ (same as `make data`)
scripts/fetch_data.sh segments tvoc    # also the statewide folders, to Data/aclima/
scripts/fetch_data.sh all              # everything, about 3 GB
```

To see all the files in the bucket:

```sh
aws s3 ls --no-sign-request --recursive s3://switchbox-ny-air-monitoring/
```

### What the data is

Aclima collected the data for the New York State Department of Environmental Conservation (NYSDEC) in 10 areas of New York State, between July 2022 and August 2023. Aclima vehicles drove the public roads in each area many times and measured nine pollutants once per second.

| Bucket folder | What it holds |
|---|---|
| `westchester/` | The 1-second data for the Westchester study area, as an R data frame. **This is the only file that the pipeline needs.** |
| `raw_1s/` | The 1-second data for all 10 areas. |
| `segments/` | One-year summary values for each road segment (about 100 m long), for all 10 areas. |
| `tvoc/` | Peaks of volatile organic compounds and air toxics, and a biogenic-methane indicator, for all 10 areas. |

The bucket's `README.md` gives the columns of each file, the measurement period of each area and the changes from the original Aclima delivery.

### Files that are in this repository

These small inputs are committed to `Data/`:

| File | What it is |
|---|---|
| `Mean_concentration_road_segments.*` | The Aclima 100-m road segments for the study area (8,391 segments). UW removed some segments at the edge of the area. |
| `Westchester_tracts_selected.*` | The census tracts of the study area. |
| `Westchester_tract_population_by_age.RData` | The population of each tract by age group. The health-rate calculations use these age groups. |

The pipeline also gets block groups and American Community Survey data from the Census Bureau at run time. This is why it needs `CENSUS_API_KEY`.

## Pipeline

Each script reads from `Data/` and `Results/`, and writes to `Results/` and `Figures/`.

| Script | What it does |
|---|---|
| `Code/01_Exposure_assessment.R` | Makes tract-level annual exposure estimates for PM2.5 and NO₂ from the 1-second data (point → block group → tract), with bootstrap uncertainty. |
| `Code/02_Health_risk.R` | Estimates the deaths and illnesses that the exposure causes, and the benefits of general concentration reductions. |
| `Code/03_Source_apportionment.R` | Finds the part of the exposure that comes from traffic, with non-negative matrix factorization (NMF). |
| `Code/04_Mitigation_strategies_NEW.R` | Estimates the health benefits of vehicle electrification and a low-emission zone. `04_Mitigation_strategies_old.R` is an earlier version, kept for reference. |
| `Code/05_Mitigation_gas_diesel_compute.R` | Estimates the benefits of the electrification of only gasoline or only diesel vehicles. It is not part of `make`. Run it by hand (see the comment at the top of the file). |
| `Code/report_*.R`, `Code/twopager_figures_absolute.R` | Make the report figures again from the saved `Results/`. They do not calculate again and do not call the Census API. |
| `Code/smoke_test.R` | Checks that the Docker image can load all packages and read the input files. |
| `Code/_helpers.R`, `Code/_theme_switchbox.R` | Shared functions, and the Switchbox chart style. |

### Order of health outcomes

In the result tables for the seven health outcomes, the order is:

1. PM2.5 all-cause mortality
2. PM2.5 stroke
3. PM2.5 lung cancer
4. PM2.5 asthma
5. NO₂ all-cause mortality
6. NO₂ lung cancer
7. NO₂ asthma

### RStudio

If you prefer RStudio to Docker, open `New_York_MM_project.Rproj`, run `make data` (or `scripts/fetch_data.sh`), and run the scripts in `Code/` in order. The `Dockerfile` lists the R packages that you must install.

## Citation

Switchbox, Energy Justice Law and Policy Center, and the University of Washington. 2026. *Westchester Air Pollution Health Impact Assessment: A blueprint for turning mobile monitoring data into community-based mitigation planning.* Switchbox.

## License

- Code: MIT (see `LICENSE`).
- Aclima data: public domain.
