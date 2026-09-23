

library(sf)
library(dplyr)
library(lubridate)
library(maptiles)
library(RColorBrewer)
library(ggplot2)
library(ggnewscale)
library(tidyr)
library(tigris)
library(tidycensus)
.census_key <- Sys.getenv("CENSUS_API_KEY")
if (!nzchar(.census_key)) stop("CENSUS_API_KEY env var is not set")
census_api_key(.census_key)




########## Notes ##########

### Section 1-3: Health risk assessment based on current tract-level exposure
### Section 4: Health benefit evaluation for hypothetical concentration reduction scenarios





########## 1. PAF estimation ##########

## Load tract-level annual average exposure with bootstrap uncertainties from Code 01
load("Results/Exposure_boot_NO2_tract.RData")
NO2_tract_boot <- PM_tract_boot # NO2 boot data, names as NO2_tract_boot

load("Results/Exposure_boot_PM2.5_tract.RData") # PM2.5 boot data, names as PM_tract_boot


## Exposure-response relationship
ER_data <- data.frame(pollutant = c(rep("PM2.5",6), rep("NO2",3)),
                      outcome = c("All-cause mortality","Stroke incidence","Lung cancer incidence","Asthma incidence",
                                  "Alzheimer hospitalization","Dementia hospitalization",
                                  "All-cause mortality","Lung cancer incidence","Asthma incidence"),
                      beta = c(0.01133,0.00343,0.037844,0.043672,0.139762,0.076961, log(1.045), log(1.02), log(1.082)),
                      SE = c(0.0016,0.001265,0.013121,0.000885,0.017753,0.018905, log(1.052/1.037)/1.96/2, log(1.04/1)/1.96/2, log(1.16/1.01)/1.96/2),
                      unit = c(rep(1,6), 8.1, 10.4, 5))


## Baseline health data
baseline_data <- data.frame(outcome = c("Mortality", "Stroke", "LC", "Asthma", "Parkinson", "Dementia"),
                            rate = c(941.22, 503.88, 135.70, 1233.65, 80.45, 229.49))

## Pollutant-outcome pairs
pair <- c(101,202,303,404,701,803,904)

# The first digit represents the row number in the ER_data, i.e., pollutant-outcome pair ID
# The last digit represents thr row number in the baseline_data, i.e., outcome ID
# Therefore, each element in pair means:
# 101 PM2.5 with all-cause mortality
# 202 PM2.5 with stroke incidence
# 303 PM2.5 with lung cancer incidence
# 404 PM2.5 with asthma incidence
# 701 NO2 with all-cause mortality
# 803 NO2 with lung cancer incidence
# 904 NO2 with asthma incidence
# This is also the order for the following list elements.


## Risk assessment function
risk_cal <- function(risk_outcome) {
  ER_id <- risk_outcome %/% 100
  baseline_id <- risk_outcome %% 100
  
  pollutant_id <- ER_data$pollutant[ER_id]
  if (pollutant_id == "PM2.5") {exposure <- PM_tract_boot} else {exposure <- NO2_tract_boot}
  beta <- ER_data$beta[ER_id]
  SE <- ER_data$SE[ER_id]
  unit <- ER_data$unit[ER_id]
  
  set.seed(12345)
  RR0 <- exp(beta + SE * rnorm(1000))^(1/unit)
  # if (pollutant_id == "PM2.5") {TMREL <- runif(1000, min = 2.4, max = 5.9)} else (TMREL <- rep(0, 1000))
  TMREL <- rep(0, 1000) # Now we assume TMREL=0 for both PM2.5 and NO2
  
  PAF <- lapply(unique(exposure$tract),
                function(i) {
                  tract_exp <- exposure %>% filter(tract == i) %>% select(-c(tract, run_id))
                  RR <- RR0^tract_exp
                  RR_avg <- apply(RR, 1, mean)
                  PAF <- pmax(1 - RR0^TMREL/RR_avg, 0)
                  return(PAF)
                })
  PAF <- as.data.frame(do.call(rbind, PAF))
  
  return(PAF)
}



## PAF for each risk-outcome pair: Monte Carlo results
PAF_pair <- list()
for (k in 1:length(pair)) {
  print(paste0("k = ",k,"; Time = ", Sys.time()))
  PAF_pair[[k]] <- risk_cal(pair[k])
}



## PAF for each risk-outcome pair: Summary statistics
summary_stat <- function(v) {
  stat <- c(mean(v),
            sd(v),
            quantile(v, probs = 0.025),
            quantile(v, probs = 0.25),
            quantile(v, probs = 0.50),
            quantile(v, probs = 0.75),
            quantile(v, probs = 0.975))
  return(stat)
}


PAF_stat <- list()
for (k in 1:length(PAF_pair)) {
  PAF_stat[[k]] <- t(apply(PAF_pair[[k]],1,summary_stat))
  colnames(PAF_stat[[k]]) <- c("Mean","SD","P2.5","P25","P50","P75","P97.5")
}



save(PAF_pair, PAF_stat, file = "Results/PAF_results_for_health_risk.RData")









########## 2. Attributable health risks ##########

#### (1) Attributable risk: Tract-level rate per 100,000
risk_rate <- list()
risk_rate_stat <- list()

for (k in 1:length(pair)) {
  # Monte Carlo results
  risk_rate[[k]] <- PAF_pair[[k]] * baseline_data$rate[pair[k] %% 100]  # Risk rate per 100,000
  
  # Summary statistics
  risk_rate_stat[[k]] <- t(apply(risk_rate[[k]],1,summary_stat))
  colnames(risk_rate_stat[[k]]) <- c("Mean","SD","P2.5","P25","P50","P75","P97.5")
}

save(risk_rate, risk_rate_stat, file = "Results/Health_risk_rate_tract.RData")




#### (2) Attributable risk: Tract-level number

## Extract the age-group-specific population data from ACS
popu_count <- get_acs(
  state = "NY",
  geography = "tract",
  variables = paste0("B01001_",sprintf("%03d", 1:49)),
  geometry = FALSE,
  year = 2023)

WEST_tract <- st_read("Data/Westchester_tracts_selected.shp")
WEST_tract <- WEST_tract[order(WEST_tract$GEOID),]
popu_count <- popu_count[popu_count$GEOID %in% WEST_tract$GEOID,]

popu_0_20 <- popu_count[popu_count$variable %in% paste0("B01001_",sprintf("%03d", c(3:8,27:32))),]
popu_0_20 <- aggregate(estimate ~ GEOID,data = popu_0_20,sum)
colnames(popu_0_20)[2] <- "Popu_0_20"

popu_over18 <- popu_count[popu_count$variable %in% paste0("B01001_",sprintf("%03d", c(7:25,31:49))),]
popu_over18 <- aggregate(estimate ~ GEOID,data = popu_over18,sum)
colnames(popu_over18)[2] <- "Popu_over18"

popu_over50 <- popu_count[popu_count$variable %in% paste0("B01001_",sprintf("%03d", c(16:25,40:49))),]
popu_over50 <- aggregate(estimate ~ GEOID,data = popu_over50,sum)
colnames(popu_over50)[2] <- "Popu_over50"

popu_over65 <- popu_count[popu_count$variable %in% paste0("B01001_",sprintf("%03d", c(20:25,44:49))),]
popu_over65 <- aggregate(estimate ~ GEOID,data = popu_over65,sum)
colnames(popu_over65)[2] <- "Popu_over65"

popu_over70 <- popu_count[popu_count$variable %in% paste0("B01001_",sprintf("%03d", c(22:25,46:49))),]
popu_over70 <- aggregate(estimate ~ GEOID,data = popu_over70,sum)
colnames(popu_over70)[2] <- "Popu_over70"

popu_stat <- cbind(popu_0_20,popu_over18[,-1,drop = F],popu_over50[,-1,drop = F],
                   popu_over65[,-1,drop = F],popu_over70[,-1,drop = F])


## Link population data to census tract GEOID
WEST_tract_popu <- merge(WEST_tract, popu_stat, by = "GEOID", all.x = T)
WEST_tract_popu <- st_drop_geometry(WEST_tract_popu) %>%
  select(GEOID, Popu_0_20, Popu_over18, Popu_over50, Popu_over65)
save(WEST_tract_popu, file = "Data/Westchester_tract_population_by_age.RData")


## Calculate attributable risk number
popu_group <- c("Popu_over18","Popu_over65","Popu_over50","Popu_0_20",
                "Popu_over18","Popu_over50","Popu_0_20") # Age group for each pollutant-outcome pair

risk_number <- list()
risk_number_stat <- list()

for (k in 1:length(pair)) {
  # Monte Carlo results
  risk_number[[k]] <- risk_rate[[k]] * WEST_tract_popu[[popu_group[k]]]/100000  # Risk number
  
  # Summary statistics
  risk_number_stat[[k]] <- t(apply(risk_number[[k]],1,summary_stat))
  colnames(risk_number_stat[[k]]) <- c("Mean","SD","P2.5","P25","P50","P75","P97.5")
}


save(risk_number, risk_number_stat, file = "Results/Health_risk_number_tract.RData")






#### (3) Attributable risk: Total number and rate for the whole study region
risk_number_tot <- c()
risk_rate_tot <- c()

for (k in 1:length(pair)) {
  risk_number_tot <- rbind(risk_number_tot, apply(risk_number[[k]],2,sum))
  risk_rate_tot <- rbind(risk_rate_tot, 
                         apply(risk_number[[k]],2,sum) / (sum(WEST_tract_popu[[popu_group[k]]])/100000))
}

risk_number_tot_stat <- t(apply(risk_number_tot, 1, summary_stat))
risk_rate_tot_stat <- t(apply(risk_rate_tot, 1, summary_stat))
colnames(risk_number_tot_stat) <- colnames(risk_rate_tot_stat) <- c("Mean","SD","P2.5","P25","P50","P75","P97.5")

save(risk_number_tot, risk_number_tot_stat, risk_rate_tot, risk_rate_tot_stat,
     file = "Results/Health_risk_number_and_rate_region.RData")








########## 3. Visualization ##########
bb  <- st_bbox(WEST_tract)
bb_sf <- st_as_sfc(bb)
bm <- get_tiles(bb_sf, provider = "OpenStreetMap", zoom = 14)
bm_df <- as.data.frame(bm, xy = TRUE)
bm_df$hex <- rgb(bm_df[[3]], bm_df[[4]], bm_df[[5]],
                 maxColorValue = 255)


pollutant_list <- c(rep("PM2.5",4), rep("NO2",3))
outcome_list <- c("Mortality","Stroke","Lung_cancer","Asthma",
                  "Mortality","Lung_cancer","Asthma")



## Plot tract-level health risk rate per 100,000
for (k in 1:length(pair)) {
  WEST_tract_risk <- cbind(WEST_tract, Mean = risk_rate_stat[[k]][,"Mean"])

  ggplot() +
    geom_raster(data = bm_df, aes(x = x, y = y, fill = hex)) +
    scale_fill_identity() +
    new_scale_fill() +
    geom_sf(data = WEST_tract_risk,
            aes(fill = Mean), alpha = 0.8,
            color = "black",
            linewidth = 0.4) +
    scale_fill_distiller(
      name = "Per 100,000",
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
  
  
  ggsave(paste0("Figures/Health_risk_rate_",pollutant_list[k],"_",outcome_list[k],"_tract.png"),
         width = 10, height = 8, dpi = 300)
}



## Plot tract-level health risk rate per 100,000
for (k in 1:length(pair)) {
  WEST_tract_PAF <- cbind(WEST_tract, Mean = PAF_stat[[k]][,"Mean"])
  
  ggplot() +
    geom_raster(data = bm_df, aes(x = x, y = y, fill = hex)) +
    scale_fill_identity() +
    new_scale_fill() +
    geom_sf(data = WEST_tract_PAF,
            aes(fill = Mean), alpha = 0.8,
            color = "black",
            linewidth = 0.4) +
    scale_fill_distiller(
      name = "PAF",
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
  
  
  ggsave(paste0("Figures/Health_risk_PAF_",pollutant_list[k],"_",outcome_list[k],"_tract.png"),
         width = 10, height = 8, dpi = 300)
  
}











########## 4. Health benefit ##########

##### (1) Attributable health benefit: Tract-level

### PIF function — defined in Code/_helpers.R so it can be reused by 04.
source("Code/_helpers.R")



### Health benefits at tract level when reducing concentrations by percentages
PIF_result <- PIF_stat <- HB_rate <- HB_number <- HB_rate_stat <- HB_number_stat <- list()

for (k in 1:length(pair)) {
  PIF_result[[k]] <- PIF_stat[[k]] <- HB_rate[[k]] <- HB_number[[k]] <- HB_rate_stat[[k]] <- HB_number_stat[[k]] <- list()

  for (delta in 1:10) {
    print(paste0("k = ",k,"; delta = ",delta,"; Time = ", Sys.time()))
    
    ## PIF
    PIF_result[[k]][[delta]] <- PIF_cal(pair[k], (delta)/10, "Percent")
    PIF_stat[[k]][[delta]] <- t(apply(PIF_result[[k]][[delta]],1,summary_stat))
    
    ## Health benefit rate per 100,000
    HB_rate[[k]][[delta]] <- PIF_result[[k]][[delta]] * baseline_data$rate[pair[k] %% 100]
    HB_rate_stat[[k]][[delta]] <- t(apply(HB_rate[[k]][[delta]],1,summary_stat))
    
    ## Health benefit number
    HB_number[[k]][[delta]] <- HB_rate[[k]][[delta]] * WEST_tract_popu[[popu_group[k]]]/100000
    HB_number_stat[[k]][[delta]] <- t(apply(HB_number[[k]][[delta]],1,summary_stat))
    
    
    ## Update column names
    colnames(PIF_stat[[k]][[delta]]) <- colnames(HB_rate_stat[[k]][[delta]]) <- 
      colnames(HB_number_stat[[k]][[delta]]) <- c("Mean","SD","P2.5","P25","P50","P75","P97.5")
    
  }
}


save(PIF_result, PIF_stat, HB_rate, HB_number, HB_rate_stat, HB_number_stat, 
     file = "Results/Health_benefit_conc_reduction_scenario_tract.RData")






##### (2) Attributable health benefit: For the whole study region

HB_number_tot <- HB_rate_tot <- list()
HB_number_tot_stat <- HB_rate_tot_stat <- HB_perc_tot_stat <- list()

for (k in 1:length(pair)) {
  
  HB_number_tot[[k]] <- HB_rate_tot[[k]] <- rep(0,1000)

  for (s in 1:length(HB_number[[k]])) {
    
    ## Region-level total benefit number
    sum_tmp <- apply(HB_number[[k]][[s]],2,sum)
    HB_number_tot[[k]] <- rbind(HB_number_tot[[k]],
                                sum_tmp)
    
    ## Region-level total benefit rate per 100,000
    HB_rate_tot[[k]] <- rbind(HB_rate_tot[[k]],
                              sum_tmp / (sum(WEST_tract_popu[[popu_group[k]]])/100000))
    
  }
  
  HB_number_tot_stat[[k]] <- data.frame(seq(0,1,by = 0.1),
                                        t(apply(HB_number_tot[[k]],1,summary_stat)))
  HB_rate_tot_stat[[k]] <- data.frame(seq(0,1,by = 0.1),
                                      t(apply(HB_rate_tot[[k]],1,summary_stat)))
  
  ## Region-level total percentage of benefit in total baseline risk
  HB_perc_tot_stat[[k]] <- data.frame(seq(0,1,by = 0.1),
                                      t(apply(HB_rate_tot[[k]],1,summary_stat)) / baseline_data$rate[pair[k] %% 100])
  
  colnames(HB_number_tot_stat[[k]]) <- colnames(HB_rate_tot_stat[[k]]) <- 
    colnames(HB_perc_tot_stat[[k]]) <- c("Reduction","Mean","SD","P2.5","P25","P50","P75","P97.5")
  
}


save(HB_number_tot, HB_rate_tot, HB_number_tot_stat, HB_rate_tot_stat, HB_perc_tot_stat,
     file = "Results/Health_benefit_conc_reduction_scenario_region.RData")







##### (3) Plot curves of attributable health benefit rate for the whole study region

outcome_full_name <- c("all-cause mortality","stroke incidence","lung cancer incidence","asthma incidence",
                       "all-cause mortality","lung cancer incidence","asthma incidence")


for (k in 1:length(pair)) {
  print(k)
  
  ggplot(HB_rate_tot_stat[[k]], aes(x = Reduction)) +
    
    # Uncertainty band
    geom_ribbon(aes(ymin = `P2.5`, ymax = `P97.5`),
                fill = "#B22222",
                alpha = 0.2,
                color = NA) +
    
    # Mean line
    geom_line(aes(y = Mean), linewidth = 1.2, color = "#B22222") +
    
    # Optional points (can remove if too busy)
    geom_point(aes(y = Mean), size = 2, color = "#B22222") +
    
    labs(
      x = ifelse(pollutant_list[k] == "PM2.5",
                 expression("Reduction percentage of PM"[2.5]),
                 expression("Reduction percentage of NO"[2])),
      y = paste0("Avoidable ",outcome_full_name[k]," rate\n(per 100,000)"),
      color = NULL,
      fill = NULL
    ) +
    
    scale_x_continuous(
      breaks = seq(0, 1, by = 0.2),
      labels = scales::percent_format(accuracy = 1),
      expand = expansion(mult = c(0, 0))
    ) +
    
    scale_y_continuous(
      expand = expansion(mult = c(0, 0))
    ) +
    
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
      plot.margin = margin(12, 24, 12, 12)
    )
  

  ggsave(paste0("Figures/Health_benefit_conc_reduc_",pollutant_list[k],"_",outcome_list[k],"_region.png"),
         width = 8, height = 6, dpi = 600)
}





##### (4) Plot prioritized tracts by avoidable numbers

ref_perc <- 0.2  ## The reduction percentage you want to emphasize to rank tracts that benefit most


for (k in 1:length(pair)) {
  print(k)
  
  df <- do.call(cbind, lapply(1:length(HB_number_stat[[k]]), 
                              function(x) HB_number_stat[[k]][[x]][,"Mean"]))
  df <- cbind(Reduction = seq(0, 1, by = 0.1), t(cbind(0, df)))
  colnames(df)[-1] <- paste0("Tract",1:(ncol(df)-1))
  
  top5_tract <- order(df[which(df[,1] == ref_perc),-1], decreasing = TRUE)[1:5]
  
  df_long <- as.data.frame(df) %>%
    pivot_longer(cols = -Reduction,
                 names_to = "tract",
                 values_to = "value") %>%
    mutate(highlight = ifelse(tract %in% paste0("Tract",top5_tract), "Top 5", "Other tracts"))
  
  
  ## Plot mean curves of attributable risk number for 70 tracts in one figure
  ggplot() +
    geom_line(
      data = df_long %>% filter(highlight == "Other tracts"),
      aes(x = Reduction, y = value, group = tract),
      color = "grey70",
      linewidth = 0.6
    ) +
    geom_line(
      data = df_long %>% filter(highlight == "Top 5"),
      aes(x = Reduction, y = value, group = tract),
      color = "#D94801",
      linewidth = 0.6
    ) +
    labs(
      x = ifelse(pollutant_list[k] == "PM2.5",
                 expression("Reduction percentage of PM"[2.5]),
                 expression("Reduction percentage of NO"[2])),
      y = paste0("Avoidable ",outcome_list[k]," number"),
      color = NULL
    ) +
    scale_x_continuous(breaks = seq(0, 1, 0.2),
                       labels = function(x) paste0(x * 100, "%"),
                       expand = expansion(mult = c(0, 0))) +
    
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
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
      plot.margin = margin(12, 24, 12, 12)
    )
  
  ggsave(paste0("Figures/Health_benefit_conc_reduc_",
                pollutant_list[k],"_",outcome_list[k],"_tract.png"),
         width = 8, height = 6, dpi = 600)
  
  
  ## Plot top 5 tracts that benefit most from PM2.5/NO2 reduction by 10%
  WEST_tract_top5 <- WEST_tract[top5_tract,]
  WEST_tract_other <- WEST_tract %>%
    filter(!GEOID %in% WEST_tract_top5$GEOID)
  
  ggplot() +
    geom_raster(data = bm_df, aes(x = x, y = y, fill = hex)) +
    scale_fill_identity() +
    
    ggnewscale::new_scale_fill() +
    
    geom_sf(
      data = WEST_tract_other,
      fill = "grey65",
      alpha = 0.6,
      color = "grey35",
      linewidth = 0.4
    ) +
    
    geom_sf(
      data = WEST_tract_top5,
      fill = "yellow",
      color = "grey20",
      linewidth = 0.7,
      alpha = 0.9,
      show.legend = FALSE
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
    coord_sf(datum = NA) +
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
  
  
  ggsave(paste0("Figures/Health_benefit_conc_reduc_",
                pollutant_list[k],"_",outcome_list[k],"_tract_top5_map.png"),
         width = 8, height = 6, dpi = 600)
  
  
}














