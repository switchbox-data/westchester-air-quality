
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
.census_key <- Sys.getenv("CENSUS_API_KEY")
if (!nzchar(.census_key)) stop("CENSUS_API_KEY env var is not set")
census_api_key(.census_key)



###### This code shows the final version: Temporally weighted average estimates (point -> block group -> census tract)
###### Comparisons with previous versions (GAM; point -> road segment -> census tract) are removed






## Read 1-second point-level calibrated data from mobile monitoring
# dat <- read.csv("Data/mvy_clipped_f.csv")
# dat <- dat[!duplicated(dat),]
load("Data/Mobile_monitoring_1second_data_Aclima.RData")


## Read 100-m road segments from Aclima (removing some edge roads by Ningrui)
seg_conc_mean <- st_read("Data/Mean_concentration_road_segments.shp") # 8391 road segments
colnames(seg_conc_mean)[1:10] <- c("segment_id","blackcarbon","c2h6","ch4","co","co2","no","no2","o3","pm_2.5")



## Read shapefiles of census tracts and block groups in our study area
WEST_tract <- st_read("Data/Westchester_tracts_selected.shp")
WEST_tract <- WEST_tract[order(WEST_tract$GEOID),]

WEST_block <- block_groups(state = "36", county = "119", cb = TRUE, year = 2023)
WEST_block <- WEST_block %>% filter(substr(GEOID, 1, 11) %in% WEST_tract$GEOID)
WEST_block <- st_transform(WEST_block, crs = 32115)
summary(WEST_block$ALAND)/1000/1000
quantile(WEST_block$ALAND/1000/1000, probs = c(0.025,0.05,0.1,0.9,0.95,0.975))




## For loop for PM2.5 and NO2 for exposure assessment.
## Override via `POLLUTANTS=pm_2.5 Rscript Code/01_...` to run one pollutant
## at a time — keeps peak memory under ~6 GB instead of ~12 GB.
for (pollut_name in strsplit(Sys.getenv("POLLUTANTS", "pm_2.5,no2"), ",")[[1]]) {
  
  pollut_std <- ifelse(pollut_name == "pm_2.5","PM2.5","NO2")
  
  ## Filter data about the focused pollutant
  PM_data <- dat %>% filter((modality == pollut_name) & (quality_flag == 0))
  PM_data$timestamp <- ymd_hms(PM_data$timestamp, tz = "UTC")
  PM_data$timestamp <- with_tz(PM_data$timestamp, tzone = "America/New_York")
  
  PM_data <- st_as_sf(
    PM_data,
    coords = c("lon", "lat"),
    crs = 4326,      # WGS84 lon/lat
    remove = FALSE  # keep lon/lat columns
  )
  PM_data <- st_transform(PM_data, crs = 32115)
  
  
  ## Link points to the nearest block group
  block_id <- st_nearest_feature(PM_data, WEST_block)
  PM_data$block <- WEST_block$GEOID[block_id]
  PM_data$tract <- substr(PM_data$block,1,11)
  
  
  ## Calculate visit-level (hour-level) median for each block group
  PM_block_hour_median <- st_drop_geometry(PM_data) %>%
    mutate(hour = floor_date(timestamp, "hour")) %>%
    group_by(block, hour) %>%
    summarise(
      median = pmax(median(value, na.rm = T),0),
      N = n(),
      .groups = "drop"
    )
  
  
  
  
  
  ###### 1. Version without bootstrap uncertainty #####
  ## This version is for exposure mapping
  
  ## Calculate block-group-level stratum-level average exposure
  ## 20 strata (season * weekday/weekend * time periods in a day) * 195 block groups
  PM_block_s_w_h <- PM_block_hour_median %>%
    group_by(block) %>%
    mutate(median_wins = Winsorize(median)) %>%
    ungroup() %>%
    mutate(month = month(hour), weekday = wday(hour), hour_id = hour(hour)) %>%
    mutate(season = ifelse(month %in% c(3,4,5),"Spring",
                           ifelse(month %in% c(6,7,8),"Summer",
                                  ifelse(month %in% c(9,10,11),"Fall","Winter"))),
           TOW = ifelse(weekday %in% c(1,7),"Weekend","Weekday"),
           
           TOD = ifelse(TOW == "Weekday",
                        ifelse(((hour_id >= 6) & (hour_id < 10)) | ((hour_id >= 16) & (hour_id < 20)),"Rush",
                               ifelse((hour_id >= 10) & (hour_id < 16),"Nonrush","Night")),
                        ifelse((hour_id >= 6) & (hour_id < 20),"Nonrush","Night"))) %>%
    group_by(block, season, TOW, TOD) %>%
    summarise(
      Mean = mean(median_wins),
      N_pass = n(),
      .groups = "drop"
    ) %>%
    complete(
      block, season, TOW, TOD,
      fill = list(Mean = NA)
    ) %>%
    filter(!((TOW == "Weekend") & (TOD == "Rush"))) %>%
    mutate(wt = 1/4 * ifelse(TOW == "Weekday", 5/7, 2/7) * 
             ifelse(TOW == "Weekday",
                    ifelse(TOD == "Rush", 8/24, ifelse(TOD == "Nonrush", 6/24, 10/24)),
                    ifelse(TOD == "Nonrush", 14/24, 10/24)))
  
  
  ## Calculate block-group-level annual average exposure
  PM_block_annual_s_w_h <- PM_block_s_w_h %>%
    group_by(block) %>%
    summarise(
      N_pass = sum(N_pass, na.rm = T),
      N_season = n_distinct(season[!is.na(Mean)]),
      N_strata = sum(!is.na(Mean)),
      Wt_SWH = sum(Mean * wt, na.rm = T) / sum(wt[!is.na(Mean)]),
      .groups = "drop"
    )
  
  summary(PM_block_annual_s_w_h[,-1])

  
  
  ## Calculate tract-level annual average exposure
  popu_block <- get_acs(
    state = "NY",
    county = "Westchester",
    geography = "block group",
    variables = "B01003_001",
    geometry = FALSE,
    year = 2023) %>%
    filter(GEOID %in% PM_block_annual_s_w_h$block) %>%
    select(GEOID, estimate)
  
  if (sum(popu_block$GEOID == PM_block_annual_s_w_h$block) != nrow(PM_block_annual_s_w_h)) stop("Order of block group GEOIDs is wrong!") # 195
  PM_block_annual_s_w_h <- cbind(PM_block_annual_s_w_h, popu = popu_block$estimate)
  
  PM_tract_annual_s_w_h <- PM_block_annual_s_w_h %>%
    mutate(tract = substr(block,1,11)) %>%
    group_by(tract) %>%
    summarise(
      Wt_SWH = sum(Wt_SWH * popu) / sum(popu),
      .groups = "drop"
    )
  
  if (sum(PM_tract_annual_s_w_h$tract == WEST_tract$GEOID) != nrow(PM_tract_annual_s_w_h)) stop("Order of tract GEOIDs is wrong!") # 70
  PM_tract_annual_s_w_h <- st_as_sf(merge(PM_tract_annual_s_w_h, WEST_tract[,"GEOID"], by.x = "tract", by.y = "GEOID"))
  
  save(PM_tract_annual_s_w_h, 
       file = paste0("Results/Exposure_",pollut_std,"_tract.RData"))
  
  
  
  ## Plot tract-level annual average exposure map
  bb  <- st_bbox(WEST_tract)
  bb_sf <- st_as_sfc(bb)
  bm <- get_tiles(bb_sf, provider = "OpenStreetMap", zoom = 14)
  bm_df <- as.data.frame(bm, xy = TRUE)
  bm_df$hex <- rgb(bm_df[[3]], bm_df[[4]], bm_df[[5]],
                   maxColorValue = 255)
  

  ggplot() +
    geom_raster(data = bm_df, aes(x = x, y = y, fill = hex)) +
    scale_fill_identity() +
    new_scale_fill() +
    geom_sf(data = PM_tract_annual_s_w_h,
            aes(fill = Wt_SWH), alpha = 0.8,
            color = "black",
            linewidth = 0.4) +
    scale_fill_distiller(
      name = ifelse(pollut_name == "pm_2.5", "μg/m³", "ppb"),
      palette = "RdYlBu",
      direction = -1,
    ) +
    guides(
      fill = guide_colorbar(
        barheight = unit(6, "cm"),
        barwidth  = unit(0.6, "cm"),
        frame.colour = "black",
        frame.linewidth = 0.4,
        ticks.colour = "black"
      )
    ) +
    theme_minimal() +
    theme_void() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      legend.position = "right",
      legend.box.spacing = unit(0.05, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14)
    )

  
  ggsave(paste0("Figures/Exposure_",pollut_std,"_tract.png"),
         width = 10, height = 8, dpi = 300)
  
  
  
  
  
  ###### 2. Temporal weighted average for each block group (with uncertainty version) #####
  ## This version is prepared for health risk assessment which needs uncertainty analysis
  
  set.seed(123)
  
  PM_tract_boot_tot <- list()
  
  for (b in 1:1000) {
    if (b %% 10 == 0) print(paste0("b = ",b,"; Time = ", Sys.time()))
    
    ## Calculate block-group-level stratum-level average exposure
    PM_block_s_w_h <- PM_block_hour_median %>%
      mutate(month = month(hour), weekday = wday(hour), hour_id = hour(hour)) %>%
      mutate(season = ifelse(month %in% c(3,4,5),"Spring",
                             ifelse(month %in% c(6,7,8),"Summer",
                                    ifelse(month %in% c(9,10,11),"Fall","Winter"))),
             TOW = ifelse(weekday %in% c(1,7),"Weekend","Weekday"),
             
             TOD = ifelse(TOW == "Weekday",
                          ifelse(((hour_id >= 6) & (hour_id < 10)) | ((hour_id >= 16) & (hour_id < 20)),"Rush",
                                 ifelse((hour_id >= 10) & (hour_id < 16),"Nonrush","Night")),
                          ifelse((hour_id >= 6) & (hour_id < 20),"Nonrush","Night"))) %>%
      
      group_by(block, season, TOW, TOD) %>%
      slice_sample(prop = 1, replace = TRUE) %>%  ## Bootstrap resampling stratified by temporal strata to keep temporal structure
      ungroup() %>%
      
      group_by(block) %>%
      mutate(median_wins = Winsorize(median)) %>%
      ungroup() %>%
      
      group_by(block, season, TOW, TOD) %>%
      summarise(
        Mean = mean(median_wins),
        N_pass = n(),
        .groups = "drop"
      ) %>%
      complete(
        block, season, TOW, TOD,
        fill = list(Mean = NA)
      ) %>%
      filter(!((TOW == "Weekend") & (TOD == "Rush"))) %>%
      mutate(wt = 1/4 * ifelse(TOW == "Weekday", 5/7, 2/7) * 
               ifelse(TOW == "Weekday",
                      ifelse(TOD == "Rush", 8/24, ifelse(TOD == "Nonrush", 6/24, 10/24)),
                      ifelse(TOD == "Nonrush", 14/24, 10/24)))
    
    
    ## Calculate block-group-level annual average exposure
    PM_block_annual_s_w_h <- PM_block_s_w_h %>%
      group_by(block) %>%
      summarise(
        N_pass = sum(N_pass, na.rm = T),
        N_season = n_distinct(season[!is.na(Mean)]),
        N_strata = sum(!is.na(Mean)),
        Wt_SWH = sum(Mean * wt, na.rm = T) / sum(wt[!is.na(Mean)]),
        .groups = "drop"
      )
    
    
    PM_block_annual_s_w_h <- cbind(PM_block_annual_s_w_h, popu = popu_block$estimate)
    
    
    ## Calculate tract-level annual average exposure
    PM_tract_boot_block <- PM_block_annual_s_w_h %>%
      mutate(tract = substr(block,1,11)) %>%
      group_by(tract) %>%
      reframe(
        boot = Wt_SWH[sample.int(
          n = dplyr::n(),
          size = 1000,
          prob = popu / sum(popu),
          replace = TRUE
        )],  ## Monte Carlo simulation, ready for health risk assessment
      ) %>%
      mutate(draw_id = row_number(), .by = tract)
    PM_tract_boot_block <- dcast(PM_tract_boot_block, tract ~ draw_id, value.var = "boot")
    PM_tract_boot_block$run_id <- b
    
    PM_tract_boot_tot[[b]] <- PM_tract_boot_block
  }
  
  PM_tract_boot <- do.call(rbind, PM_tract_boot_tot)
  PM_tract_boot <- PM_tract_boot %>% arrange(tract, run_id)
  
  
  save(PM_tract_boot, 
       file = paste0("Results/Exposure_boot_",pollut_std,"_tract.RData"))
  
}






