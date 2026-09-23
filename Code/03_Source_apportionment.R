
library(sf)
library(dplyr)
library(ggplot2)
library(ggnewscale)
library(lubridate)
library(reshape2)
library(DescTools)
library(tidyr)
library(maptiles)
library(RColorBrewer)


library(NMF)

library(tigris)
library(tidycensus)
.census_key <- Sys.getenv("CENSUS_API_KEY")
if (!nzchar(.census_key)) stop("CENSUS_API_KEY env var is not set")
census_api_key(.census_key)



########## 1. Input data preparation ##########

## Read Aclima's mobile monitoring 1-second point-level data
load("Data/Mobile_monitoring_1second_data_Aclima.RData")
dat$timestamp <- ymd_hms(dat$timestamp, tz = "UTC")
dat$timestamp <- with_tz(dat$timestamp, tzone = "America/New_York")

dat_wide <- dat %>% filter(quality_flag == 0)
dat_wide <- dcast(data = dat_wide,
                  timestamp + lat + lon ~ modality,
                  value.var = "value")

dat_NMF <- st_as_sf(
  dat_wide,
  coords = c("lon", "lat"),
  crs = 4326,      # WGS84 lon/lat
  remove = FALSE  # keep lon/lat columns
)


## Read Aclima's 100-m road segment
seg_conc_mean <- st_read("Data/Mean_concentration_road_segments.shp")
colnames(seg_conc_mean)[1] <- "segment_id"



## Calculate the median for each visit of each road segment
seg_id <- st_nearest_feature(dat_NMF, seg_conc_mean)
dat_NMF$segment_id <- seg_conc_mean$segment_id[seg_id]


dat_seg_hour_median <- st_drop_geometry(dat_NMF) %>%
  mutate(hour = floor_date(timestamp, "hour")) %>%
  group_by(segment_id, hour) %>%
  summarise(
    PM2.5 = pmax(median(pm_2.5, na.rm = T),0),
    BC = pmax(median(blackcarbon, na.rm = T),0),
    NO2 = pmax(median(no2, na.rm = T),0),
    NO = pmax(median(no, na.rm = T),0),
    CO = pmax(median(co, na.rm = T),0),
    CO2 = pmax(median(co2, na.rm = T),0),
    Ozone = pmax(median(o3, na.rm = T),0),
    CH4 = pmax(median(ch4, na.rm = T),0),
    C2H6 = pmax(median(c2h6, na.rm = T),0),
    N = n(),
    .groups = "drop"
  )

## Free the 27.9M-row raw frames before NMF — they're no longer referenced.
rm(dat, dat_wide, dat_NMF); gc()


## Background subtraction for CO2, CH4, and C2H6
dat_seg_hour_median$date <- date(dat_seg_hour_median$hour)

dat_seg_hour_bg_sub <- dat_seg_hour_median %>%
  group_by(date) %>%
  mutate(CO2 = CO2 - min(CO2, na.rm = T),
         CH4 = CH4 - min(CH4, na.rm = T),
         C2H6 = C2H6 - min(C2H6, na.rm = T)) %>%
  ungroup() %>%
  mutate(PM2.5 = Winsorize(PM2.5, val = quantile(PM2.5, probs = c(0.05,0.95), na.rm = T)),
         BC = Winsorize(BC, val = quantile(BC, probs = c(0.05,0.95), na.rm = T)),
         NO2 = Winsorize(NO2, val = quantile(NO2, probs = c(0.05,0.95), na.rm = T)),
         NO = Winsorize(NO, val = quantile(NO, probs = c(0.05,0.95), na.rm = T)),
         CO = Winsorize(CO, val = quantile(CO, probs = c(0.05,0.95), na.rm = T)),
         CO2 = Winsorize(CO2, val = quantile(CO2, probs = c(0,0.95), na.rm = T)),
         Ozone = Winsorize(Ozone, val = quantile(Ozone, probs = c(0.05,0.95), na.rm = T)),
         CH4 = Winsorize(CH4, val = quantile(CH4, probs = c(0,0.95), na.rm = T)),
         C2H6 = Winsorize(C2H6, val = quantile(C2H6, probs = c(0,0.95), na.rm = T)))
# Winsorize overall (only 95th percentile for CO2, CH4, and C2H6), 312,541 rows


dat_seg_hour_bg_sub$ratio <- apply(dat_seg_hour_bg_sub,1,function(x) {9 - sum(is.na(x))})
sum(dat_seg_hour_bg_sub$ratio >= 9) # 183,140 rows (58.6%)
sum(dat_seg_hour_bg_sub$ratio >= 8) # 223,210 rows (71.4%)
sum(dat_seg_hour_bg_sub$ratio >= 7) # 289,409 rows (recommended, 92.6%)
sum(dat_seg_hour_bg_sub$ratio >= 6) # 303,811 rows (97.2%)


dat_seg_hour_NMF <- dat_seg_hour_bg_sub %>% filter(ratio >= 7) # 289,409 rows

dat_seg_pass <- dat_seg_hour_NMF %>%
  group_by(segment_id) %>%
  summarise(
    N_pass = n(),
    .groups = "drop"
  )



## Weight matrix and modified conc matrix
# MDL: PM2.5, BC, NO2, NO, CO, CO2, Ozone, CH4, C2H6
MDL <- c(0.9, 0.6, 3.1, 4.0, 0.026, 3.3, 1.8, 34, 2.4)

conc <- dat_seg_hour_NMF[,3:11]
wt <- conc


for (k in 1:9) {conc[which(conc[,k] <= MDL[k]),k] <- 0.5 * MDL[k]}
conc[is.na(conc)] <- 0
conc <- as.matrix(conc)

for (k in 1:9) {
  wt[which(wt[,k] <= MDL[k]),k] <- 1 / (5/6 * MDL[k])^2
  wt[which(wt[,k] > MDL[k]),k] <- 1 / ((0.1 * wt[which(wt[,k] > MDL[k]),k])^2 + (0.5 * MDL[k])^2)
}
wt[is.na(wt)] <- 0
wt <- as.matrix(wt)

save(dat_seg_hour_NMF, conc, wt, file = "Results/NMF_input_matrices.RData")





########## 2. Weighted NMF analysis ##########

########## (2-1) Try different number of factors: 2-6 factors ##########

## Factor-count sweep. The downstream analysis only requires the 3-factor
## fit, so the default skips the 2,4,5,6 fits (each takes ~22 min sequential).
## Override to do the full sweep with `NMF_FACTORS=2,3,4,5,6`.
.nmf_factors <- as.integer(strsplit(Sys.getenv("NMF_FACTORS", "3"), ",")[[1]])
.full_sweep  <- all(2:6 %in% .nmf_factors)

## Fit NMF models
for (k in .nmf_factors) {
  print(paste0("Number of factors = ",k,"; Time = ",Sys.time()))
  
  fit <- nmf(
    x = conc,
    rank = k,
    method = "ls-nmf",
    weight = wt,
    seed = 123,
    nrun = 20,
    .options = list(parallel = as.integer(Sys.getenv("NMF_PARALLEL", "20")))
  )

  save(fit, file = paste0("Results/NMF_model_factor_",k,".RData"))
}



## Calculate IM and IS (model-selection diagnostic across the 2..6 sweep).
## Only runs when the full sweep was actually fit; otherwise we skip
## straight to the final 3-factor analysis below.
if (.full_sweep) {
IM <- c()
IS <- c()
for (k in 2:6) {
  load(paste0("Results/NMF_model_factor_",k,".RData"))

  G_contrib <- basis(fit)   # Factor contribution, n x k
  F_profile <- coef(fit)    # Factor profile, k x p

  G_mean <- apply(G_contrib,2,mean)
  G_contrib <- sweep(G_contrib, 2, G_mean, "/")  # Keep the column average of G matrix as 1, same as PMF
  F_profile <- F_profile * G_mean

  res_scaled <- sqrt(wt) * (conc - fitted(fit))   # Scaled residual


  # Calculate IM and IS
  IM <- c(IM, max(apply(res_scaled,2,mean)))
  IS <- c(IS, max(apply(res_scaled,2,sd)))
}




IM_IS <- data.frame(Factor = 2:6, IM = IM, IS = IS) %>%
  pivot_longer(cols = c(IM, IS),
               names_to = "Metric",
               values_to = "Value")

ggplot(IM_IS, aes(x = Factor, y = Value, color = Metric)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    x = "Number of factors",
    y = "IM or IS",
    color = NULL
  ) +
  
  scale_color_manual(values = c(
    "IM" = "#B22222",
    "IS" = "#7A9CC6"
  )) +
  
  scale_x_continuous(breaks = seq(2, 6, by = 1)) +
  
  theme_classic(base_size = 16) +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.text = element_text(size = 16),
    axis.title.x = element_text(size = 18, margin = margin(t = 12)),
    axis.title.y = element_text(size = 18, margin = margin(r = 10)),
    axis.text = element_text(size = 16, color = "black"),
    axis.ticks = element_line(linewidth = 1, color = "black"),
    axis.ticks.length = unit(0.18, "cm"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major = element_line(color = "#D0D0D0", linewidth = 0.5, linetype = "dashed"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(12, 18, 12, 12)
  )

ggsave("Figures/NMF_model_IM_IS.png",
       width = 8, height = 6, dpi = 600)
}  # end of if (.full_sweep)





## Final model: Three factor version
load("Results/NMF_model_factor_3.RData")
G_contrib <- basis(fit)   # Factor contribution, n x k
F_profile <- coef(fit)    # Factor profile, k x p

G_mean <- apply(G_contrib,2,mean)
G_contrib <- sweep(G_contrib, 2, G_mean, "/")  # Keep the column average of G matrix as 1, same as PMF
F_profile <- F_profile * G_mean







########## (2-2) Source contribution for interpretation ##########
G_contrib <- data.frame(dat_seg_hour_NMF[,1:2], G_contrib)
colnames(G_contrib)[-(1:2)] <- paste0("Fac",1:3)


##### (1) Weighted annual average factor contribution
G_seg_s_w_h <- G_contrib %>%
  mutate(month = month(hour), weekday = wday(hour), hour_id = hour(hour)) %>%
  mutate(season = ifelse(month %in% c(3,4,5),"Spring",
                         ifelse(month %in% c(6,7,8),"Summer",
                                ifelse(month %in% c(9,10,11),"Fall","Winter"))),
         TOW = ifelse(weekday %in% c(1,7),"Weekend","Weekday"),
         
         TOD = ifelse(TOW == "Weekday",
                      ifelse(((hour_id >= 6) & (hour_id < 10)) | ((hour_id >= 16) & (hour_id < 20)),"Rush",
                             ifelse((hour_id >= 10) & (hour_id < 16),"Nonrush","Night")),
                      ifelse((hour_id >= 6) & (hour_id < 20),"Nonrush","Night"))) %>%
  group_by(segment_id, season, TOW, TOD) %>%
  summarise(
    Fac1 = mean(Fac1),
    Fac2 = mean(Fac2),
    Fac3 = mean(Fac3),
    N_pass = n(),
    .groups = "drop"
  ) %>%
  complete(
    segment_id, season, TOW, TOD,
    fill = list(Fac1 = NA, Fac2 = NA, Fac3 = NA, Fac4 = NA)
  ) %>%
  filter(!((TOW == "Weekend") & (TOD == "Rush"))) %>%
  mutate(wt = 1/4 * ifelse(TOW == "Weekday", 5/7, 2/7) * 
           ifelse(TOW == "Weekday",
                  ifelse(TOD == "Rush", 8/24, ifelse(TOD == "Nonrush", 6/24, 10/24)),
                  ifelse(TOD == "Nonrush", 14/24, 10/24)))



## Annual average contribution

G_seg_annual_s_w_h <- G_seg_s_w_h %>%
  group_by(segment_id) %>%
  summarise(
    N_season = n_distinct(season[!is.na(Fac1)]),
    N_strata = sum(!is.na(Fac1)),
    Fac1 = sum(Fac1 * wt, na.rm = T) / sum(wt[!is.na(Fac1)]),
    Fac2 = sum(Fac2 * wt, na.rm = T) / sum(wt[!is.na(Fac2)]),
    Fac3 = sum(Fac3 * wt, na.rm = T) / sum(wt[!is.na(Fac3)]),
    .groups = "drop"
  ) %>%
  left_join(
    dat_seg_pass %>% select(segment_id, N_pass),
    by = c("segment_id" = "segment_id")
  ) %>%
  filter(N_pass >= 20, N_season >= 3, N_strata >= 12)

G_seg_annual_s_w_h <- merge(G_seg_annual_s_w_h, seg_conc_mean[,"segment_id"],
                            by = "segment_id", all.x = T)
G_seg_annual_s_w_h <- st_as_sf(G_seg_annual_s_w_h)




for (k in 1:3) {
  ggplot() +
    geom_sf(data = G_seg_annual_s_w_h,
            aes(color = .data[[paste0("Fac",k)]]), alpha = 1,
            linewidth = 0.4) +
    scale_color_distiller(
      limits = c(0.5, 1.5),
      palette = "RdYlBu",
      direction = -1,
      oob = scales::squish,
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
    ) + 
    labs(color = "Contribution",
         title = paste0("Weighted annual average contribution of Factor ",k))
  
  ggsave(paste0("Figures/NMF_model_annual_contribution_Fac",k,".png"),
         width = 10, height = 8, dpi = 600)
  
}







##### (2) Seasonal variation
G_seg_seasonal <- G_seg_s_w_h %>%
  group_by(segment_id, season) %>%
  summarise(
    N_strata = sum(!is.na(Fac1)),
    Fac1 = sum(Fac1 * wt, na.rm = T) / sum(wt[!is.na(Fac1)]),
    Fac2 = sum(Fac2 * wt, na.rm = T) / sum(wt[!is.na(Fac2)]),
    Fac3 = sum(Fac3 * wt, na.rm = T) / sum(wt[!is.na(Fac3)]),
    N_pass = sum(N_pass, na.rm = T),
    .groups = "drop"
  ) %>%
  filter(N_pass >= 5, N_strata >= 3)

G_seg_seasonal <- merge(G_seg_seasonal, seg_conc_mean[,"segment_id"],
                        by = "segment_id", all.x = T)
G_seg_seasonal <- st_as_sf(G_seg_seasonal)



for (k in 1:3) {
  for (season_label in c("Spring","Summer","Fall","Winter")) {
    ggplot() +
      geom_sf(data = G_seg_seasonal %>% filter(season == season_label),
              aes(color = .data[[paste0("Fac",k)]]), alpha = 1,
              linewidth = 0.4) +
      scale_color_distiller(
        limits = c(0.5, 1.5),
        palette = "RdYlBu",
        direction = -1,
        oob = scales::squish,
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
      ) + 
      labs(color = "Contribution",
           title = paste0("Average contribution of Factor ",k," in ",season_label))
    
    
    ggsave(paste0("Figures/NMF_model_seasonal_contribution_Fac",k,"_",season_label,".png"),
           width = 10, height = 8, dpi = 600)
  }
}







##### (3) Rush vs non-rush hours variation
G_seg_rush <- G_seg_s_w_h %>%
  group_by(segment_id, TOD) %>%
  summarise(
    N_strata = sum(!is.na(Fac1)),
    Fac1 = sum(Fac1 * wt, na.rm = T) / sum(wt[!is.na(Fac1)]),
    Fac2 = sum(Fac2 * wt, na.rm = T) / sum(wt[!is.na(Fac2)]),
    Fac3 = sum(Fac3 * wt, na.rm = T) / sum(wt[!is.na(Fac3)]),
    # Fac4 = sum(Fac4 * wt, na.rm = T) / sum(wt[!is.na(Fac4)]),
    N_pass = sum(N_pass, na.rm = T),
    .groups = "drop"
  ) %>%
  filter(N_pass >= 5, N_strata >= 3)

G_seg_rush <- merge(G_seg_rush, seg_conc_mean[,"segment_id"],
                    by = "segment_id", all.x = T)
G_seg_rush <- st_as_sf(G_seg_rush)






for (k in 1:3) {
  for (rush_label in c("Rush","Nonrush","Night")) {
    ggplot() +
      geom_sf(data = G_seg_rush %>% filter(TOD == rush_label),
              aes(color = .data[[paste0("Fac",k)]]), alpha = 1,
              linewidth = 0.4) +
      scale_color_distiller(
        limits = c(0.5, 1.5),
        palette = "RdYlBu",
        direction = -1,
        oob = scales::squish,
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
      ) + 
      labs(color = "Contribution",
           title = paste0("Average contribution of Factor ",k," in ",rush_label))
    
    
    ggsave(paste0("Figures/NMF_model_diurnal_contribution_Fac",k,"_",rush_label,".png"),
           width = 10, height = 8, dpi = 600)
  }
}







########## 3. Source-specific exposure ##########
##### Follow Visit -> Block group -> Tract approach, similar to exposure assessment


## Read source profile and contribution matrices
# load("Results/NMF_model_factor_3.RData")
# load("Results/NMF_input_matrices.RData")

G_contrib <- basis(fit)   # Factor contribution, n x k
F_profile <- coef(fit)    # Factor profile, k x p

G_mean <- apply(G_contrib,2,mean)
G_contrib <- sweep(G_contrib, 2, G_mean, "/")  # Keep the column average of G matrix as 1, same as PMF
F_profile <- F_profile * G_mean

G_contrib <- data.frame(dat_seg_hour_NMF[,1:2], G_contrib)
colnames(G_contrib)[-(1:2)] <- paste0("Fac",1:3)



## Read shapefiles of tracts and blocks
WEST_tract <- st_read("Data/Westchester_tracts_selected.shp")
WEST_tract <- WEST_tract[order(WEST_tract$GEOID),]

WEST_block <- block_groups(state = "36", county = "119", cb = TRUE, year = 2023)
WEST_block <- WEST_block %>% filter(substr(GEOID, 1, 11) %in% WEST_tract$GEOID)
WEST_block <- st_transform(WEST_block, crs = 32115)




## Main code for annual traffic-related exposure

for (pollut_name in c("PM2.5", "NO2")) {
  
  ## Road-segment-level visit-level source-specific exposure
  PM_data <- as.matrix(G_contrib[,-(1:2)]) * matrix(rep(as.numeric(F_profile[,pollut_name]), nrow(G_contrib)),
                                                    ncol = 3, byrow = T)
  PM_data <- data.frame(G_contrib[,1:2], PM_data)
  
  
  ## Load tract-level total exposure estimates
  load(paste0("Results/Exposure_",pollut_name,"_tract.RData"))
  PM_tract_tot <- PM_tract_annual_s_w_h
  
  
  ## Link 100-m road segments to the nearest block group
  PM_data <- merge(PM_data, seg_conc_mean[,"segment_id"],
                   by = "segment_id", all.x = T) %>% st_as_sf()
  PM_data <- st_transform(PM_data, crs = 32115)
  
  
  block_id <- st_nearest_feature(PM_data, WEST_block)
  PM_data$block <- WEST_block$GEOID[block_id]
  PM_data$tract <- substr(PM_data$block,1,11)
  
  
  ## Calculate block-group-level weighted averages
  PM_block_s_w_h <- st_drop_geometry(PM_data) %>%
    group_by(block) %>%
    mutate(Fac1 = Winsorize(Fac1),
           Fac2 = Winsorize(Fac2),
           Fac3 = Winsorize(Fac3)) %>%
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
      Fac1 = mean(Fac1),
      Fac2 = mean(Fac2),
      Fac3 = mean(Fac3),
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
  
  
  PM_block_annual_s_w_h <- PM_block_s_w_h %>%
    group_by(block) %>%
    summarise(
      N_pass = sum(N_pass, na.rm = T),
      N_season = n_distinct(season[!is.na(Fac1)]),
      N_strata = sum(!is.na(Fac1)),
      Fac1 = sum(Fac1 * wt, na.rm = T) / sum(wt[!is.na(Fac1)]),
      Fac2 = sum(Fac2 * wt, na.rm = T) / sum(wt[!is.na(Fac2)]),
      Fac3 = sum(Fac3 * wt, na.rm = T) / sum(wt[!is.na(Fac3)]),
      # Fac4 = sum(Fac4 * wt, na.rm = T) / sum(wt[!is.na(Fac4)]),
      .groups = "drop"
    )
  
  
  
  ## Calculate tract-level weighted averages
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
      Fac1 = sum(Fac1 * popu) / sum(popu),
      Fac2 = sum(Fac2 * popu) / sum(popu),
      Fac3 = sum(Fac3 * popu) / sum(popu),
      .groups = "drop"
    )
  
  
  ## Calculate tract-level source contribution percentages
  PM_tract_annual_s_w_h <- merge(PM_tract_annual_s_w_h, st_drop_geometry(PM_tract_tot),
                                 by = "tract", all.x = T)
  colnames(PM_tract_annual_s_w_h)[5] <- "Conc_tot"
  
  PM_tract_annual_s_w_h <- PM_tract_annual_s_w_h %>%
    mutate(Fac1_perc = Fac1 / pmax(Fac1 + Fac2 + Fac3, Conc_tot),
           Fac2_perc = Fac2 / pmax(Fac1 + Fac2 + Fac3, Conc_tot),
           Fac3_perc = Fac3 / pmax(Fac1 + Fac2 + Fac3, Conc_tot))
  
  
  
  if (sum(PM_tract_annual_s_w_h$tract == WEST_tract$GEOID) != nrow(PM_tract_annual_s_w_h)) stop("Order of tract GEOIDs is wrong!") # 70
  PM_tract_annual_s_w_h <- st_as_sf(merge(PM_tract_annual_s_w_h, WEST_tract[,"GEOID"], by.x = "tract", by.y = "GEOID"))
  
  
  save(PM_tract_annual_s_w_h, file = paste0("Results/Exposure_source_",pollut_name,"_tract.RData"))
  

  
  
  ## Plot traffic-related PM2.5 exposure
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
            aes(fill = Fac1), alpha = 0.8,
            color = "black",
            linewidth = 0.4) +
    scale_fill_distiller(
      name = ifelse(pollut_name == "PM2.5","μg/m³","ppb"),
      palette = "RdYlBu",
      direction = -1,
      oob = scales::squish,
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
  
  
  ggsave(paste0("Figures/Exposure_traffic_",pollut_name,"_tract_conc.png"),
         width = 10, height = 8, dpi = 300)
  
  
  
  ## Plot traffic-related PM2.5 percentage
  ggplot() +
    geom_raster(data = bm_df, aes(x = x, y = y, fill = hex)) +
    scale_fill_identity() +
    new_scale_fill() +
    geom_sf(data = PM_tract_annual_s_w_h,
            aes(fill = Fac1_perc), alpha = 0.8,
            color = "black",
            linewidth = 0.4) +
    scale_fill_distiller(
      name = "Percentage",
      palette = "RdYlBu",
      direction = -1,
      labels = scales::percent_format(accuracy = 1)
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
  
  
  ggsave(paste0("Figures/Exposure_traffic_",pollut_name,"_tract_perc.png"),
         width = 10, height = 8, dpi = 300)
  

}









