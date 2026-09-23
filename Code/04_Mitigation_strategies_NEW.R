
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(ggnewscale)
library(maptiles)


########## 1. Vehicle electrification ##########

### MOVES model output
MOVES_data <- data.frame(
  Fuel = c("Gasoline", "Diesel", "CNG", "Ethanol", "Electric"),
  VMT = c(6310050.048, 683544.208, 8013.151, 15283.118, 105657.032),
  NO2 = c(140800.539, 304993.8, 1471.393, 182.665, 0),
  PM2.5_exhaust = c(34062.615, 54305.411, 107.988, 63.238, 0),
  PM2.5_brakewear = c(17303.548, 18793.611, 431.523, 41.768, 85.502),
  PM2.5_tirewear = c(8066.363, 2150.267, 19.4, 19.354, 136.724)
)


MOVES_data <- MOVES_data %>%
  mutate(PM2.5 = PM2.5_exhaust + PM2.5_brakewear + PM2.5_tirewear,
         EF_NO2 = NO2 / VMT,
         EF_PM2.5 = PM2.5 / VMT,
         VMT_perc = VMT / sum(VMT))





########## (1-1) Emission reduction ##########
EV_unif <- data.frame(replacement = c(0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5))
EV_unif <- EV_unif %>%
  mutate(NO2 = replacement * sum(MOVES_data$VMT[1:4] * (MOVES_data$EF_NO2[1:4] - MOVES_data$EF_NO2[5])),
         NO2_reduction = NO2 / sum(MOVES_data$NO2),
         PM2.5 = replacement * sum(MOVES_data$VMT[1:4] * (MOVES_data$EF_PM2.5[1:4] - MOVES_data$EF_PM2.5[5])),
         PM2.5_reduction = PM2.5 / sum(MOVES_data$PM2.5))

EV_unif_long <- EV_unif %>%
  select(-NO2, -PM2.5) %>%
  rename(
    NO2 = NO2_reduction,
    PM2.5 = PM2.5_reduction
  ) %>%
  pivot_longer(cols = c(NO2, PM2.5),
               names_to = "Pollutant",
               values_to = "Reduction") %>%
  mutate(Scenario = "Uniform")



EV_diesel <- data.frame(replacement = c(0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5))
EV_diesel <- EV_diesel %>%
  mutate(NO2 = pmin(replacement/MOVES_data$VMT_perc[2],1) * MOVES_data$VMT[2] * (MOVES_data$EF_NO2[2] - MOVES_data$EF_NO2[5]) + 
           pmax((replacement - MOVES_data$VMT_perc[2])/(1 - MOVES_data$VMT_perc[2]),0) * sum(MOVES_data$VMT[c(1,3,4)] * (MOVES_data$EF_NO2[c(1,3,4)] - MOVES_data$EF_NO2[5])),
         NO2_reduction = NO2 / sum(MOVES_data$NO2),
         PM2.5 = pmin(replacement/MOVES_data$VMT_perc[2],1) * MOVES_data$VMT[2] * (MOVES_data$EF_PM2.5[2] - MOVES_data$EF_PM2.5[5]) + 
           pmax((replacement - MOVES_data$VMT_perc[2])/(1 - MOVES_data$VMT_perc[2]),0) * sum(MOVES_data$VMT[c(1,3,4)] * (MOVES_data$EF_PM2.5[c(1,3,4)] - MOVES_data$EF_PM2.5[5])),
         PM2.5_reduction = PM2.5 / sum(MOVES_data$PM2.5))


EV_diesel_long <- EV_diesel %>%
  select(-NO2, -PM2.5) %>%
  rename(
    NO2 = NO2_reduction,
    PM2.5 = PM2.5_reduction
  ) %>%
  pivot_longer(cols = c(NO2, PM2.5),
               names_to = "Pollutant",
               values_to = "Reduction") %>%
  mutate(Scenario = "Diesel")


EV_emission <- rbind(EV_unif_long, EV_diesel_long)


ggplot(EV_emission, aes(x = replacement, y = Reduction, color = Pollutant, linetype = Scenario)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    x = "EV replacement",
    y = "Emission reduction",
    color = NULL,
    linetype = NULL
  ) +
  
  # Custom colors and linetypes
  scale_color_manual(values = c(
    "PM2.5" = "#B22222",
    "NO2" = "#7A9CC6"
  )) +
  
  scale_linetype_manual(values = c(
    "Uniform" = "solid",
    "Diesel"  = "22"
  )) + 
  
  # Axis formatting
  scale_x_continuous(breaks = seq(0, 0.5, by = 0.1),
                     labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0))) +
  
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0))) +
  
  # Theme
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

ggsave(paste0("Figures/Mitigation_EV_emission_reduction.png"),
       width = 8, height = 6, dpi = 600)









########## (1-2) Concentration reduction ##########

EV_unif_conc_result <- c()
EV_unif_conc_list <- list()

for (pollut_name in c("PM2.5", "NO2")) {
  load(paste0("Results/Exposure_source_",pollut_name,"_tract.RData"))
  
  EV_unif_conc_reduc <- PM_tract_annual_s_w_h[,c("tract","Fac1_perc")]
  
  for (k in 1:nrow(EV_unif)) {
    EV_unif_conc_reduc[[paste0("Case",k)]] <- EV_unif_conc_reduc$Fac1_perc * EV_unif[[paste0(pollut_name,"_reduction")]][k]
  }
  
  EV_unif_conc_list[[pollut_name]] <- EV_unif_conc_reduc
  
  EV_unif_conc_result <- rbind(EV_unif_conc_result,
                               data.frame(replacement = EV_unif$replacement,
                                          Pollutant = pollut_name,
                                          Mean = apply(st_drop_geometry(EV_unif_conc_reduc[,-(1:2)]), 2, mean),
                                          P2.5 = apply(st_drop_geometry(EV_unif_conc_reduc[,-(1:2)]), 2, function(x) quantile(x, probs = 0.025)),
                                          P97.5 = apply(st_drop_geometry(EV_unif_conc_reduc[,-(1:2)]), 2, function(x) quantile(x, probs = 0.975))))
}


EV_unif_conc_result$Scenario <- "Uniform"





EV_diesel_conc_result <- c()
EV_diesel_conc_list <- list()


for (pollut_name in c("PM2.5", "NO2")) {
  load(paste0("Results/Exposure_source_",pollut_name,"_tract.RData"))
  
  EV_diesel_conc_reduc <- PM_tract_annual_s_w_h[,c("tract","Fac1_perc")]
  
  for (k in 1:nrow(EV_diesel)) {
    EV_diesel_conc_reduc[[paste0("Case",k)]] <- EV_diesel_conc_reduc$Fac1_perc * EV_diesel[[paste0(pollut_name,"_reduction")]][k]
  }
  
  EV_diesel_conc_list[[pollut_name]] <- EV_diesel_conc_reduc
  
  EV_diesel_conc_result <- rbind(EV_diesel_conc_result,
                                 data.frame(replacement = EV_diesel$replacement,
                                            Pollutant = pollut_name,
                                            Mean = apply(st_drop_geometry(EV_diesel_conc_reduc[,-(1:2)]), 2, mean),
                                            P2.5 = apply(st_drop_geometry(EV_diesel_conc_reduc[,-(1:2)]), 2, function(x) quantile(x, probs = 0.025)),
                                            P97.5 = apply(st_drop_geometry(EV_diesel_conc_reduc[,-(1:2)]), 2, function(x) quantile(x, probs = 0.975))))
}

EV_diesel_conc_result$Scenario <- "Diesel"

EV_concentration <- rbind(EV_unif_conc_result, EV_diesel_conc_result)



ggplot(EV_concentration, 
       aes(x = replacement, color = Pollutant, fill = Pollutant, linetype = Scenario)) +
  
  # Uncertainty band
  geom_ribbon(aes(ymin = P2.5, ymax = P97.5),
              alpha = 0.2,
              color = NA) +
  
  # Mean line
  geom_line(aes(y = Mean), linewidth = 1.2) +
  
  # Optional points (can remove if too busy)
  geom_point(aes(y = Mean), size = 2) +
  
  labs(
    x = "EV replacement",
    y = "Concentration reduction",
    color = NULL,
    fill = NULL,
    linetype = NULL
  ) +
  
  scale_color_manual(values = c(
    "PM2.5" = "#B22222",
    "NO2" = "#7A9CC6"
  )) +
  
  scale_fill_manual(values = c(
    "PM2.5" = "#B22222",
    "NO2" = "#7A9CC6"
  )) +
  
  scale_linetype_manual(values = c(
    "Uniform" = "solid",
    "Diesel"  = "22"
  )) + 
  
  scale_x_continuous(
    breaks = seq(0, 0.5, by = 0.1),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0))
  ) +
  
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
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
    plot.margin = margin(12, 18, 12, 12)
  )



ggsave(paste0("Figures/Mitigation_EV_concentration_reduction.png"),
       width = 8, height = 6, dpi = 600)




########## (1-3) Health benefit ##########

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


## Population data
WEST_tract <- st_read("Data/Westchester_tracts_selected.shp")
WEST_tract <- WEST_tract[order(WEST_tract$GEOID),]
load("Data/Westchester_tract_population_by_age.RData")
popu_group <- c("Popu_over18","Popu_over65","Popu_over50","Popu_0_20",
                "Popu_over18","Popu_over50","Popu_0_20") # Age group for each pollutant-outcome pair



## PIF function
PIF_mitigation <- function(risk_outcome, C_perc) {
  ER_id <- risk_outcome %/% 100
  baseline_id <- risk_outcome %% 100
  
  pollutant_id <- ER_data$pollutant[ER_id]
  if (pollutant_id == "PM2.5") {exposure <- PM_tract_boot} else {exposure <- NO2_tract_boot}
  beta <- ER_data$beta[ER_id]
  SE <- ER_data$SE[ER_id]
  unit <- ER_data$unit[ER_id]
  
  set.seed(12345)
  RR0 <- exp(beta + SE * rnorm(1000))^(1/unit)
  colnames(C_perc) <- c("tract","perc")
  
  PIF <- lapply(unique(exposure$tract),
                function(i) {
                  tract_exp <- exposure %>% filter(tract == i) %>% select(-c(tract, run_id))
                  RR <- RR0^tract_exp
                  RR_limit <- RR0^as.matrix(tract_exp * (1 - C_perc$perc[C_perc$tract == i]))
                  PIF <- apply(RR - RR_limit, 1, mean) / apply(RR, 1, mean)
                  return(PIF)
                })
  PIF <- as.data.frame(do.call(rbind, PIF))
  
  return(PIF)
}



## Summary statistic function
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





##### (1) Uniform replacement
HB_EV_unif_rate <- HB_EV_unif_rate_stat <- list() # Tract-level rate
HB_EV_unif_number <- HB_EV_unif_number_stat <- list() # Tract-level number

HB_EV_unif_number_tot <- HB_EV_unif_number_tot_stat <- list() # Regional-level number
HB_EV_unif_rate_tot <- HB_EV_unif_rate_tot_stat <- list() # Regional-level rate
HB_EV_unif_PIF_stat <- list() # Regional-level percentage


for (k in 1:length(pair)) {
  
  HB_EV_unif_rate[[k]] <- HB_EV_unif_rate_stat[[k]] <- list()
  HB_EV_unif_number[[k]] <- HB_EV_unif_number_stat[[k]] <- list()
  
  HB_EV_unif_number_tot[[k]] <- HB_EV_unif_number_tot_stat[[k]] <- numeric(0)
  HB_EV_unif_rate_tot[[k]] <- HB_EV_unif_rate_tot_stat[[k]] <- numeric(0)
  
  
  for (s in 1:nrow(EV_unif)) {
    print(paste0("k = ",k,"; s = ",s,"; Time = ", Sys.time()))
    
    ## Tract level
    HB_EV_unif_rate[[k]][[s]] <- PIF_mitigation(pair[k], st_drop_geometry(EV_unif_conc_list[[ER_data$pollutant[pair[k] %/% 100]]][,c("tract",paste0("Case",s))])) * 
      baseline_data$rate[pair[k] %% 100]
    
    HB_EV_unif_number[[k]][[s]] <- HB_EV_unif_rate[[k]][[s]] * WEST_tract_popu[[popu_group[k]]]/100000
    

    
    
    ## Region level
    HB_EV_unif_number_tot[[k]] <- rbind(HB_EV_unif_number_tot[[k]],
                                        apply(HB_EV_unif_number[[k]][[s]], 2, sum))
    HB_EV_unif_rate_tot[[k]] <- rbind(HB_EV_unif_rate_tot[[k]],
                                      HB_EV_unif_number_tot[[k]][s,] / (sum(WEST_tract_popu[[popu_group[k]]])/100000))
    
    
    ## Summary statistics
    HB_EV_unif_rate_stat[[k]][[s]] <- t(apply(HB_EV_unif_rate[[k]][[s]], 1, summary_stat))
    HB_EV_unif_number_stat[[k]][[s]] <- t(apply(HB_EV_unif_number[[k]][[s]], 1, summary_stat))
    colnames(HB_EV_unif_rate_stat[[k]][[s]]) <- colnames(HB_EV_unif_number_stat[[k]][[s]]) <- 
      c("Mean","SD","P2.5","P25","P50","P75","P97.5")
  }
  
  HB_EV_unif_rate_tot_stat[[k]] <- t(apply(HB_EV_unif_rate_tot[[k]], 1, summary_stat))
  HB_EV_unif_number_tot_stat[[k]] <- t(apply(HB_EV_unif_number_tot[[k]], 1, summary_stat))
  HB_EV_unif_PIF_stat[[k]] <- HB_EV_unif_rate_tot_stat[[k]] / baseline_data$rate[pair[k] %% 100]
  
  colnames(HB_EV_unif_rate_tot_stat[[k]]) <- colnames(HB_EV_unif_number_tot_stat[[k]]) <- 
   colnames(HB_EV_unif_PIF_stat[[k]]) <- c("Mean","SD","P2.5","P25","P50","P75","P97.5")
  
  
}



save(HB_EV_unif_rate, HB_EV_unif_rate_stat,
     HB_EV_unif_number, HB_EV_unif_number_stat,
     HB_EV_unif_number_tot, HB_EV_unif_number_tot_stat,
     HB_EV_unif_rate_tot, HB_EV_unif_rate_tot_stat,
     HB_EV_unif_PIF_stat,
     file = "Results/Mitigation_EV_unif_health_benefit.RData")








##### (2) Diesel replacement first
HB_EV_diesel_rate <- HB_EV_diesel_rate_stat <- list() # Tract-level rate
HB_EV_diesel_number <- HB_EV_diesel_number_stat <- list() # Tract-level number

HB_EV_diesel_number_tot <- HB_EV_diesel_number_tot_stat <- list() # Regional-level number
HB_EV_diesel_rate_tot <- HB_EV_diesel_rate_tot_stat <- list() # Regional-level rate
HB_EV_diesel_PIF_stat <- list() # Regional-level percentage


for (k in 1:length(pair)) {
  
  HB_EV_diesel_rate[[k]] <- HB_EV_diesel_rate_stat[[k]] <- list()
  HB_EV_diesel_number[[k]] <- HB_EV_diesel_number_stat[[k]] <- list()
  
  HB_EV_diesel_number_tot[[k]] <- HB_EV_diesel_number_tot_stat[[k]] <- numeric(0)
  HB_EV_diesel_rate_tot[[k]] <- HB_EV_diesel_rate_tot_stat[[k]] <- numeric(0)
  
  
  for (s in 1:nrow(EV_diesel)) {
    print(paste0("k = ",k,"; s = ",s,"; Time = ", Sys.time()))
    
    ## Tract level
    HB_EV_diesel_rate[[k]][[s]] <- PIF_mitigation(pair[k], st_drop_geometry(EV_diesel_conc_list[[ER_data$pollutant[pair[k] %/% 100]]][,c("tract",paste0("Case",s))])) * 
      baseline_data$rate[pair[k] %% 100]
    
    HB_EV_diesel_number[[k]][[s]] <- HB_EV_diesel_rate[[k]][[s]] * WEST_tract_popu[[popu_group[k]]]/100000
    
    
    
    
    ## Region level
    HB_EV_diesel_number_tot[[k]] <- rbind(HB_EV_diesel_number_tot[[k]],
                                        apply(HB_EV_diesel_number[[k]][[s]], 2, sum))
    HB_EV_diesel_rate_tot[[k]] <- rbind(HB_EV_diesel_rate_tot[[k]],
                                      HB_EV_diesel_number_tot[[k]][s,] / (sum(WEST_tract_popu[[popu_group[k]]])/100000))
    
    
    ## Summary statistics
    HB_EV_diesel_rate_stat[[k]][[s]] <- t(apply(HB_EV_diesel_rate[[k]][[s]], 1, summary_stat))
    HB_EV_diesel_number_stat[[k]][[s]] <- t(apply(HB_EV_diesel_number[[k]][[s]], 1, summary_stat))
    colnames(HB_EV_diesel_rate_stat[[k]][[s]]) <- colnames(HB_EV_diesel_number_stat[[k]][[s]]) <- 
      c("Mean","SD","P2.5","P25","P50","P75","P97.5")
  }
  
  HB_EV_diesel_rate_tot_stat[[k]] <- t(apply(HB_EV_diesel_rate_tot[[k]], 1, summary_stat))
  HB_EV_diesel_number_tot_stat[[k]] <- t(apply(HB_EV_diesel_number_tot[[k]], 1, summary_stat))
  HB_EV_diesel_PIF_stat[[k]] <- HB_EV_diesel_rate_tot_stat[[k]] / baseline_data$rate[pair[k] %% 100]
  
  colnames(HB_EV_diesel_rate_tot_stat[[k]]) <- colnames(HB_EV_diesel_number_tot_stat[[k]]) <- 
    colnames(HB_EV_diesel_PIF_stat[[k]]) <- c("Mean","SD","P2.5","P25","P50","P75","P97.5")
  
  
}



save(HB_EV_diesel_rate, HB_EV_diesel_rate_stat,
     HB_EV_diesel_number, HB_EV_diesel_number_stat,
     HB_EV_diesel_number_tot, HB_EV_diesel_number_tot_stat,
     HB_EV_diesel_rate_tot, HB_EV_diesel_rate_tot_stat,
     HB_EV_diesel_PIF_stat,
     file = "Results/Mitigation_EV_diesel_health_benefit.RData")










##### (3) Plot curves of health benefits for the whole study region
pollutant_list <- c(rep("PM2.5",4), rep("NO2",3))
outcome_list <- c("all-cause mortality","stroke incidence","lung cancer incidence","asthma incidence",
                  "all-cause mortality","lung cancer incidence","asthma incidence")


for (k in 1:length(pair)) {
  print(k)
  
  ## Plot avoidable risk rate
  EV_health_benefit <- rbind(data.frame(HB_EV_unif_rate_tot_stat[[k]], Scenario = "Uniform"),
                             data.frame(HB_EV_diesel_rate_tot_stat[[k]], Scenario = "Diesel"))
  EV_health_benefit <- cbind(replacement = rep(EV_unif$replacement,2), EV_health_benefit)
  
  ggplot(EV_health_benefit,
         aes(x = replacement, color = Scenario, fill = Scenario)) +
    
    # Uncertainty band
    geom_ribbon(aes(ymin = `P2.5`, ymax = `P97.5`),
                alpha = 0.2,
                color = NA) +
    
    # Mean line
    geom_line(aes(y = Mean), linewidth = 1.2) +
    
    # Optional points (can remove if too busy)
    geom_point(aes(y = Mean), size = 2) +
    
    labs(
      x = "EV replacement",
      y = paste0("Avoidable ",outcome_list[k]," rate\n(per 100,000)"),
      color = NULL,
      fill = NULL
    ) +
    
    scale_color_manual(values = c(
      "Uniform" = "#00B050",
      "Diesel" = "#FF8C00"
    )) +
    
    scale_fill_manual(values = c(
      "Uniform" = "#00B050",
      "Diesel" = "#FF8C00"
    )) +
    
    scale_x_continuous(
      breaks = seq(0, 0.5, by = 0.1),
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
      plot.margin = margin(12, 18, 12, 12)
    )
  
  
  ggsave(paste0("Figures/Mitigation_EV_HB_rate_",pollutant_list[k],"_",outcome_list[k],"_region.png"),
         width = 8, height = 6, dpi = 600)
  
  
  
  
  ## Plot avoidable percentage for baseline risk
  EV_health_percent <- rbind(data.frame(HB_EV_unif_PIF_stat[[k]], Scenario = "Uniform"),
                             data.frame(HB_EV_diesel_PIF_stat[[k]], Scenario = "Diesel"))
  EV_health_percent <- cbind(replacement = rep(EV_unif$replacement,2), EV_health_percent)
  
  
  ggplot(EV_health_percent, 
         aes(x = replacement, color = Scenario, fill = Scenario)) +
    
    # Uncertainty band
    geom_ribbon(aes(ymin = `P2.5`, ymax = `P97.5`),
                alpha = 0.2,
                color = NA) +
    
    # Mean line
    geom_line(aes(y = Mean), linewidth = 1.2) +
    
    # Optional points (can remove if too busy)
    geom_point(aes(y = Mean), size = 2) +
    
    labs(
      x = "EV replacement",
      y = paste0("Avoidable ",outcome_list[k]),
      color = NULL,
      fill = NULL
    ) +
    
    scale_color_manual(values = c(
      "Uniform" = "#00B050",
      "Diesel" = "#FF8C00"
    )) +
    
    scale_fill_manual(values = c(
      "Uniform" = "#00B050",
      "Diesel" = "#FF8C00"
    )) +
    
    scale_x_continuous(
      breaks = seq(0, 0.5, by = 0.1),
      labels = scales::percent_format(accuracy = 1),
      expand = expansion(mult = c(0, 0))
    ) +
    
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 0.1),
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
      plot.margin = margin(12, 18, 12, 12)
    )
  
  
  ggsave(paste0("Figures/Mitigation_EV_HB_perc_",pollutant_list[k],"_",outcome_list[k],"_region.png"),
         width = 8, height = 6, dpi = 600)
  
  
}






##### (4) Plot top 5 tracts that benefit most
bb  <- st_bbox(WEST_tract)
bb_sf <- st_as_sfc(bb)
bm <- get_tiles(bb_sf, provider = "OpenStreetMap", zoom = 14)
bm_df <- as.data.frame(bm, xy = TRUE)
bm_df$hex <- rgb(bm_df[[3]], bm_df[[4]], bm_df[[5]],
                 maxColorValue = 255)



ref_perc <- 0.2  ## The EV replacement percentage you want to emphasize to rank tracts that benefit most



## Uniform replacement
for (k in 1:length(pair)) {
  print(k)
  
  df <- do.call(cbind, lapply(1:length(HB_EV_unif_number_stat[[k]]), 
                              function(x) HB_EV_unif_number_stat[[k]][[x]][,1]))
  df <- t(df)
  colnames(df) <- paste0("Tract",1:70)
  
  top5_tract <- order(df[which(EV_unif$replacement == ref_perc),], decreasing = TRUE)[1:5]
  
  df_long <- as.data.frame(df) %>%
    mutate(reduction = EV_unif$replacement) %>%
    pivot_longer(cols = -reduction,
                 names_to = "tract",
                 values_to = "value") %>%
    mutate(highlight = ifelse(tract %in% paste0("Tract",top5_tract), "Top 5", "Other tracts"))
  
  
  
  ggplot() +
    geom_line(
      data = df_long %>% filter(highlight == "Other tracts"),
      aes(x = reduction, y = value, group = tract),
      color = "grey70",
      linewidth = 0.6
    ) +
    geom_line(
      data = df_long %>% filter(highlight == "Top 5"),
      aes(x = reduction, y = value, group = tract),
      color = "#D94801",
      linewidth = 0.6
    ) +
    labs(
      x = "EV replacement",
      y = paste0("Avoidable ",outcome_list[k]," number"),
      color = NULL
    ) +
    scale_x_continuous(breaks = seq(0, 0.5, 0.1),
                       labels = function(x) paste0(x * 100, "%"),
                       expand = expansion(mult = c(0, 0))) +
    
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.justification = "center",
      legend.text = element_text(size = 12),
      axis.title.x = element_text(size = 14, margin = margin(t = 12)),
      axis.title.y = element_text(size = 14, margin = margin(r = 10)),
      axis.text = element_text(size = 12, color = "black"),
      axis.ticks = element_line(linewidth = 1, color = "black"),
      axis.ticks.length = unit(0.18, "cm"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.grid.major = element_line(color = "#D0D0D0", linewidth = 0.5, linetype = "dashed"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(12, 18, 12, 12)
    )
  
  ggsave(paste0("Figures/Mitigation_EV_HB_unif_number_",
                pollutant_list[k],"_",outcome_list[k],"_tract.png"),
         width = 8, height = 6, dpi = 600)
  
  
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
  
  
  ggsave(paste0("Figures/Mitigation_EV_HB_unif_number_",
                pollutant_list[k],"_",outcome_list[k],"_tract_top5_map.png"),
         width = 8, height = 6, dpi = 600)
  
  
}



## Diesel-first replacement
for (k in 1:length(pair)) {
  print(k)
  
  df <- do.call(cbind, lapply(1:length(HB_EV_diesel_number_stat[[k]]), 
                              function(x) HB_EV_diesel_number_stat[[k]][[x]][,1]))
  df <- t(df)
  colnames(df) <- paste0("Tract",1:70)
  
  top5_tract <- order(df[which(EV_diesel$replacement == ref_perc),], decreasing = TRUE)[1:5]
  
  df_long <- as.data.frame(df) %>%
    mutate(reduction = EV_diesel$replacement) %>%
    pivot_longer(cols = -reduction,
                 names_to = "tract",
                 values_to = "value") %>%
    mutate(highlight = ifelse(tract %in% paste0("Tract",top5_tract), "Top 5", "Other tracts"))
  
  
  
  ggplot() +
    geom_line(
      data = df_long %>% filter(highlight == "Other tracts"),
      aes(x = reduction, y = value, group = tract),
      color = "grey70",
      linewidth = 0.6
    ) +
    geom_line(
      data = df_long %>% filter(highlight == "Top 5"),
      aes(x = reduction, y = value, group = tract),
      color = "#D94801",
      linewidth = 0.6
    ) +
    labs(
      x = "EV replacement",
      y = paste0("Avoidable ",outcome_list[k]," number"),
      color = NULL
    ) +
    scale_x_continuous(breaks = seq(0, 0.5, 0.1),
                       labels = function(x) paste0(x * 100, "%"),
                       expand = expansion(mult = c(0, 0))) +
    
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.justification = "center",
      legend.text = element_text(size = 12),
      axis.title.x = element_text(size = 14, margin = margin(t = 12)),
      axis.title.y = element_text(size = 14, margin = margin(r = 10)),
      axis.text = element_text(size = 12, color = "black"),
      axis.ticks = element_line(linewidth = 1, color = "black"),
      axis.ticks.length = unit(0.18, "cm"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.grid.major = element_line(color = "#D0D0D0", linewidth = 0.5, linetype = "dashed"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(12, 18, 12, 12)
    )
  
  ggsave(paste0("Figures/Mitigation_EV_HB_diesel_number_",
                pollutant_list[k],"_",outcome_list[k],"_tract.png"),
         width = 8, height = 6, dpi = 600)
  
  
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
  
  
  ggsave(paste0("Figures/Mitigation_EV_HB_diesel_number_",
                pollutant_list[k],"_",outcome_list[k],"_tract_top5_map.png"),
         width = 8, height = 6, dpi = 600)
  
  
}








########## 1A. Vehicle electrification NEW: gasoline or diesel ##########

### MOVES model output
MOVES_data <- data.frame(
  Fuel = c("Gasoline", "Diesel", "CNG", "Ethanol", "Electric"),
  VMT = c(6310050.048, 683544.208, 8013.151, 15283.118, 105657.032),
  NO2 = c(140800.539, 304993.8, 1471.393, 182.665, 0),
  PM2.5_exhaust = c(34062.615, 54305.411, 107.988, 63.238, 0),
  PM2.5_brakewear = c(17303.548, 18793.611, 431.523, 41.768, 85.502),
  PM2.5_tirewear = c(8066.363, 2150.267, 19.4, 19.354, 136.724)
)


MOVES_data <- MOVES_data %>%
  mutate(PM2.5 = PM2.5_exhaust + PM2.5_brakewear + PM2.5_tirewear,
         EF_NO2 = NO2 / VMT,
         EF_PM2.5 = PM2.5 / VMT,
         VMT_perc = VMT / sum(VMT))





########## (1-1) Emission reduction ##########
EV_gas <- data.frame(replacement = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
EV_gas <- EV_gas %>%
  mutate(NO2 = replacement * MOVES_data$VMT[1] * (MOVES_data$EF_NO2[1] - MOVES_data$EF_NO2[5]),
         NO2_reduction = NO2 / sum(MOVES_data$NO2),
         PM2.5 = replacement * MOVES_data$VMT[1] * (MOVES_data$EF_PM2.5[1] - MOVES_data$EF_PM2.5[5]),
         PM2.5_reduction = PM2.5 / sum(MOVES_data$PM2.5))

EV_gas_long <- EV_gas %>%
  select(-NO2, -PM2.5) %>%
  rename(
    NO2 = NO2_reduction,
    PM2.5 = PM2.5_reduction
  ) %>%
  pivot_longer(cols = c(NO2, PM2.5),
               names_to = "Pollutant",
               values_to = "Reduction") %>%
  mutate(Scenario = "Gasoline")




ggplot(EV_gas_long, aes(x = replacement, y = Reduction, color = Pollutant)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    x = "Fraction of gasoline vehicles replaced",
    y = "Emission reduction",
    color = NULL,
    linetype = NULL
  ) +
  
  # Custom colors and linetypes
  scale_color_manual(values = c(
    "PM2.5" = "#B22222",
    "NO2" = "#7A9CC6"
  )) +
  
  # Axis formatting
  scale_x_continuous(breaks = seq(0, 1, by = 0.1),
                     labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0))) +
  
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0))) +
  
  # Theme
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

ggsave(paste0("Figures/Mitigation_EV_emission_reduction_gas.png"),
       width = 8, height = 6, dpi = 600)












EV_diesel <- data.frame(replacement = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
EV_diesel <- EV_diesel %>%
  mutate(NO2 = replacement * MOVES_data$VMT[2] * (MOVES_data$EF_NO2[2] - MOVES_data$EF_NO2[5]),
         NO2_reduction = NO2 / sum(MOVES_data$NO2),
         PM2.5 = replacement * MOVES_data$VMT[2] * (MOVES_data$EF_PM2.5[2] - MOVES_data$EF_PM2.5[5]),
         PM2.5_reduction = PM2.5 / sum(MOVES_data$PM2.5))


EV_diesel_long <- EV_diesel %>%
  select(-NO2, -PM2.5) %>%
  rename(
    NO2 = NO2_reduction,
    PM2.5 = PM2.5_reduction
  ) %>%
  pivot_longer(cols = c(NO2, PM2.5),
               names_to = "Pollutant",
               values_to = "Reduction") %>%
  mutate(Scenario = "Diesel")



ggplot(EV_diesel_long, aes(x = replacement, y = Reduction, color = Pollutant)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    x = "Fraction of diesel vehicles replaced",
    y = "Emission reduction",
    color = NULL,
    linetype = NULL
  ) +
  
  # Custom colors and linetypes
  scale_color_manual(values = c(
    "PM2.5" = "#B22222",
    "NO2" = "#7A9CC6"
  )) +
  
  # Axis formatting
  scale_x_continuous(breaks = seq(0, 1, by = 0.1),
                     labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0))) +
  
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0))) +
  
  # Theme
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

ggsave(paste0("Figures/Mitigation_EV_emission_reduction_diesel.png"),
       width = 8, height = 6, dpi = 600)









########## (1-2) Concentration reduction ##########

EV_gas_conc_result <- c()
EV_gas_conc_list <- list()

for (pollut_name in c("PM2.5", "NO2")) {
  load(paste0("Results/Exposure_source_",pollut_name,"_tract.RData"))
  
  EV_gas_conc_reduc <- PM_tract_annual_s_w_h[,c("tract","Fac1_perc")]
  
  for (k in 1:nrow(EV_gas)) {
    EV_gas_conc_reduc[[paste0("Case",k)]] <- EV_gas_conc_reduc$Fac1_perc * EV_gas[[paste0(pollut_name,"_reduction")]][k]
  }
  
  EV_gas_conc_list[[pollut_name]] <- EV_gas_conc_reduc
  
  EV_gas_conc_result <- rbind(EV_gas_conc_result,
                               data.frame(replacement = EV_gas$replacement,
                                          Pollutant = pollut_name,
                                          Mean = apply(st_drop_geometry(EV_gas_conc_reduc[,-(1:2)]), 2, mean),
                                          P2.5 = apply(st_drop_geometry(EV_gas_conc_reduc[,-(1:2)]), 2, function(x) quantile(x, probs = 0.025)),
                                          P97.5 = apply(st_drop_geometry(EV_gas_conc_reduc[,-(1:2)]), 2, function(x) quantile(x, probs = 0.975))))
}


EV_gas_conc_result$Scenario <- "Gasoline"





ggplot(EV_gas_conc_result, 
       aes(x = replacement, color = Pollutant, fill = Pollutant)) +
  
  # Uncertainty band
  geom_ribbon(aes(ymin = P2.5, ymax = P97.5),
              alpha = 0.2,
              color = NA) +
  
  # Mean line
  geom_line(aes(y = Mean), linewidth = 1.2) +
  
  # Optional points (can remove if too busy)
  geom_point(aes(y = Mean), size = 2) +
  
  labs(
    x = "Fraction of gasoline vehicles replaced",
    y = "Concentration reduction",
    color = NULL,
    fill = NULL,
    linetype = NULL
  ) +
  
  scale_color_manual(values = c(
    "PM2.5" = "#B22222",
    "NO2" = "#7A9CC6"
  )) +
  
  scale_fill_manual(values = c(
    "PM2.5" = "#B22222",
    "NO2" = "#7A9CC6"
  )) +
  
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.1),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0))
  ) +
  
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
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
    plot.margin = margin(12, 18, 12, 12)
  )



ggsave(paste0("Figures/Mitigation_EV_concentration_reduction_gasoline.png"),
       width = 8, height = 6, dpi = 600)














EV_diesel_conc_result <- c()
EV_diesel_conc_list <- list()


for (pollut_name in c("PM2.5", "NO2")) {
  load(paste0("Results/Exposure_source_",pollut_name,"_tract.RData"))
  
  EV_diesel_conc_reduc <- PM_tract_annual_s_w_h[,c("tract","Fac1_perc")]
  
  for (k in 1:nrow(EV_diesel)) {
    EV_diesel_conc_reduc[[paste0("Case",k)]] <- EV_diesel_conc_reduc$Fac1_perc * EV_diesel[[paste0(pollut_name,"_reduction")]][k]
  }
  
  EV_diesel_conc_list[[pollut_name]] <- EV_diesel_conc_reduc
  
  EV_diesel_conc_result <- rbind(EV_diesel_conc_result,
                                 data.frame(replacement = EV_diesel$replacement,
                                            Pollutant = pollut_name,
                                            Mean = apply(st_drop_geometry(EV_diesel_conc_reduc[,-(1:2)]), 2, mean),
                                            P2.5 = apply(st_drop_geometry(EV_diesel_conc_reduc[,-(1:2)]), 2, function(x) quantile(x, probs = 0.025)),
                                            P97.5 = apply(st_drop_geometry(EV_diesel_conc_reduc[,-(1:2)]), 2, function(x) quantile(x, probs = 0.975))))
}

EV_diesel_conc_result$Scenario <- "Diesel"




ggplot(EV_diesel_conc_result, 
       aes(x = replacement, color = Pollutant, fill = Pollutant)) +
  
  # Uncertainty band
  geom_ribbon(aes(ymin = P2.5, ymax = P97.5),
              alpha = 0.2,
              color = NA) +
  
  # Mean line
  geom_line(aes(y = Mean), linewidth = 1.2) +
  
  # Optional points (can remove if too busy)
  geom_point(aes(y = Mean), size = 2) +
  
  labs(
    x = "Fraction of diesel vehicles replaced",
    y = "Concentration reduction",
    color = NULL,
    fill = NULL,
    linetype = NULL
  ) +
  
  scale_color_manual(values = c(
    "PM2.5" = "#B22222",
    "NO2" = "#7A9CC6"
  )) +
  
  scale_fill_manual(values = c(
    "PM2.5" = "#B22222",
    "NO2" = "#7A9CC6"
  )) +
  
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.1),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0))
  ) +
  
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
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
    plot.margin = margin(12, 18, 12, 12)
  )



ggsave(paste0("Figures/Mitigation_EV_concentration_reduction_diesel.png"),
       width = 8, height = 6, dpi = 600)





########## (1-3) Health benefit ##########

## Load tract-level annual average exposure with bootstrap uncertainties from Code 01
load("Results/Exposure_boot_NO2_tract.RData") # NO2 boot data, names as NO2_tract_boot

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


## Population data
WEST_tract <- st_read("Data/Westchester_tracts_selected.shp")
WEST_tract <- WEST_tract[order(WEST_tract$GEOID),]
load("Data/Westchester_tract_population_by_age.RData")
popu_group <- c("Popu_over18","Popu_over65","Popu_over50","Popu_0_20",
                "Popu_over18","Popu_over50","Popu_0_20") # Age group for each pollutant-outcome pair



## PIF function
PIF_mitigation <- function(risk_outcome, C_perc) {
  ER_id <- risk_outcome %/% 100
  baseline_id <- risk_outcome %% 100
  
  pollutant_id <- ER_data$pollutant[ER_id]
  if (pollutant_id == "PM2.5") {exposure <- PM_tract_boot} else {exposure <- NO2_tract_boot}
  beta <- ER_data$beta[ER_id]
  SE <- ER_data$SE[ER_id]
  unit <- ER_data$unit[ER_id]
  
  set.seed(12345)
  RR0 <- exp(beta + SE * rnorm(1000))^(1/unit)
  colnames(C_perc) <- c("tract","perc")
  
  PIF <- lapply(unique(exposure$tract),
                function(i) {
                  tract_exp <- exposure %>% filter(tract == i) %>% select(-c(tract, run_id))
                  RR <- RR0^tract_exp
                  RR_limit <- RR0^as.matrix(tract_exp * (1 - C_perc$perc[C_perc$tract == i]))
                  PIF <- apply(RR - RR_limit, 1, mean) / apply(RR, 1, mean)
                  return(PIF)
                })
  PIF <- as.data.frame(do.call(rbind, PIF))
  
  return(PIF)
}



## Summary statistic function
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





##### (1) Gasoline vehicle replacement
HB_EV_gas_rate <- HB_EV_gas_rate_stat <- list() # Tract-level rate
HB_EV_gas_number <- HB_EV_gas_number_stat <- list() # Tract-level number

HB_EV_gas_number_tot <- HB_EV_gas_number_tot_stat <- list() # Regional-level number
HB_EV_gas_rate_tot <- HB_EV_gas_rate_tot_stat <- list() # Regional-level rate
HB_EV_gas_PIF_stat <- list() # Regional-level percentage


for (k in 1:length(pair)) {
  
  HB_EV_gas_rate[[k]] <- HB_EV_gas_rate_stat[[k]] <- list()
  HB_EV_gas_number[[k]] <- HB_EV_gas_number_stat[[k]] <- list()
  
  HB_EV_gas_number_tot[[k]] <- HB_EV_gas_number_tot_stat[[k]] <- numeric(0)
  HB_EV_gas_rate_tot[[k]] <- HB_EV_gas_rate_tot_stat[[k]] <- numeric(0)
  
  
  for (s in 1:nrow(EV_gas)) {
    print(paste0("k = ",k,"; s = ",s,"; Time = ", Sys.time()))
    
    ## Tract level
    HB_EV_gas_rate[[k]][[s]] <- PIF_mitigation(pair[k], st_drop_geometry(EV_gas_conc_list[[ER_data$pollutant[pair[k] %/% 100]]][,c("tract",paste0("Case",s))])) * 
      baseline_data$rate[pair[k] %% 100]
    
    HB_EV_gas_number[[k]][[s]] <- HB_EV_gas_rate[[k]][[s]] * WEST_tract_popu[[popu_group[k]]]/100000
    
    
    
    
    ## Region level
    HB_EV_gas_number_tot[[k]] <- rbind(HB_EV_gas_number_tot[[k]],
                                        apply(HB_EV_gas_number[[k]][[s]], 2, sum))
    HB_EV_gas_rate_tot[[k]] <- rbind(HB_EV_gas_rate_tot[[k]],
                                      HB_EV_gas_number_tot[[k]][s,] / (sum(WEST_tract_popu[[popu_group[k]]])/100000))
    
    
    ## Summary statistics
    HB_EV_gas_rate_stat[[k]][[s]] <- t(apply(HB_EV_gas_rate[[k]][[s]], 1, summary_stat))
    HB_EV_gas_number_stat[[k]][[s]] <- t(apply(HB_EV_gas_number[[k]][[s]], 1, summary_stat))
    colnames(HB_EV_gas_rate_stat[[k]][[s]]) <- colnames(HB_EV_gas_number_stat[[k]][[s]]) <- 
      c("Mean","SD","P2.5","P25","P50","P75","P97.5")
  }
  
  HB_EV_gas_rate_tot_stat[[k]] <- t(apply(HB_EV_gas_rate_tot[[k]], 1, summary_stat))
  HB_EV_gas_number_tot_stat[[k]] <- t(apply(HB_EV_gas_number_tot[[k]], 1, summary_stat))
  HB_EV_gas_PIF_stat[[k]] <- HB_EV_gas_rate_tot_stat[[k]] / baseline_data$rate[pair[k] %% 100]
  
  colnames(HB_EV_gas_rate_tot_stat[[k]]) <- colnames(HB_EV_gas_number_tot_stat[[k]]) <- 
    colnames(HB_EV_gas_PIF_stat[[k]]) <- c("Mean","SD","P2.5","P25","P50","P75","P97.5")
  
  
}



save(HB_EV_gas_rate, HB_EV_gas_rate_stat,
     HB_EV_gas_number, HB_EV_gas_number_stat,
     HB_EV_gas_number_tot, HB_EV_gas_number_tot_stat,
     HB_EV_gas_rate_tot, HB_EV_gas_rate_tot_stat,
     HB_EV_gas_PIF_stat,
     file = "Results/Mitigation_EV_gas_health_benefit.RData")








##### (2) Diesel replacement first
HB_EV_diesel_rate <- HB_EV_diesel_rate_stat <- list() # Tract-level rate
HB_EV_diesel_number <- HB_EV_diesel_number_stat <- list() # Tract-level number

HB_EV_diesel_number_tot <- HB_EV_diesel_number_tot_stat <- list() # Regional-level number
HB_EV_diesel_rate_tot <- HB_EV_diesel_rate_tot_stat <- list() # Regional-level rate
HB_EV_diesel_PIF_stat <- list() # Regional-level percentage


for (k in 1:length(pair)) {
  
  HB_EV_diesel_rate[[k]] <- HB_EV_diesel_rate_stat[[k]] <- list()
  HB_EV_diesel_number[[k]] <- HB_EV_diesel_number_stat[[k]] <- list()
  
  HB_EV_diesel_number_tot[[k]] <- HB_EV_diesel_number_tot_stat[[k]] <- numeric(0)
  HB_EV_diesel_rate_tot[[k]] <- HB_EV_diesel_rate_tot_stat[[k]] <- numeric(0)
  
  
  for (s in 1:nrow(EV_diesel)) {
    print(paste0("k = ",k,"; s = ",s,"; Time = ", Sys.time()))
    
    ## Tract level
    HB_EV_diesel_rate[[k]][[s]] <- PIF_mitigation(pair[k], st_drop_geometry(EV_diesel_conc_list[[ER_data$pollutant[pair[k] %/% 100]]][,c("tract",paste0("Case",s))])) * 
      baseline_data$rate[pair[k] %% 100]
    
    HB_EV_diesel_number[[k]][[s]] <- HB_EV_diesel_rate[[k]][[s]] * WEST_tract_popu[[popu_group[k]]]/100000
    
    
    
    
    ## Region level
    HB_EV_diesel_number_tot[[k]] <- rbind(HB_EV_diesel_number_tot[[k]],
                                          apply(HB_EV_diesel_number[[k]][[s]], 2, sum))
    HB_EV_diesel_rate_tot[[k]] <- rbind(HB_EV_diesel_rate_tot[[k]],
                                        HB_EV_diesel_number_tot[[k]][s,] / (sum(WEST_tract_popu[[popu_group[k]]])/100000))
    
    
    ## Summary statistics
    HB_EV_diesel_rate_stat[[k]][[s]] <- t(apply(HB_EV_diesel_rate[[k]][[s]], 1, summary_stat))
    HB_EV_diesel_number_stat[[k]][[s]] <- t(apply(HB_EV_diesel_number[[k]][[s]], 1, summary_stat))
    colnames(HB_EV_diesel_rate_stat[[k]][[s]]) <- colnames(HB_EV_diesel_number_stat[[k]][[s]]) <- 
      c("Mean","SD","P2.5","P25","P50","P75","P97.5")
  }
  
  HB_EV_diesel_rate_tot_stat[[k]] <- t(apply(HB_EV_diesel_rate_tot[[k]], 1, summary_stat))
  HB_EV_diesel_number_tot_stat[[k]] <- t(apply(HB_EV_diesel_number_tot[[k]], 1, summary_stat))
  HB_EV_diesel_PIF_stat[[k]] <- HB_EV_diesel_rate_tot_stat[[k]] / baseline_data$rate[pair[k] %% 100]
  
  colnames(HB_EV_diesel_rate_tot_stat[[k]]) <- colnames(HB_EV_diesel_number_tot_stat[[k]]) <- 
    colnames(HB_EV_diesel_PIF_stat[[k]]) <- c("Mean","SD","P2.5","P25","P50","P75","P97.5")
  
  
}



save(HB_EV_diesel_rate, HB_EV_diesel_rate_stat,
     HB_EV_diesel_number, HB_EV_diesel_number_stat,
     HB_EV_diesel_number_tot, HB_EV_diesel_number_tot_stat,
     HB_EV_diesel_rate_tot, HB_EV_diesel_rate_tot_stat,
     HB_EV_diesel_PIF_stat,
     file = "Results/Mitigation_EV_diesel_health_benefit.RData")










##### (3) Plot curves of health benefits for the whole study region
pollutant_list <- c(rep("PM2.5",4), rep("NO2",3))
outcome_list <- c("all-cause mortality","stroke incidence","lung cancer incidence","asthma incidence",
                  "all-cause mortality","lung cancer incidence","asthma incidence")


for (k in 1:length(pair)) {
  print(k)
  
  ## Plot avoidable risk rate
  EV_health_benefit <- rbind(data.frame(HB_EV_gas_rate_tot_stat[[k]], Scenario = "Gasoline"),
                             data.frame(HB_EV_diesel_rate_tot_stat[[k]], Scenario = "Diesel"))
  EV_health_benefit <- cbind(replacement = rep(EV_gas$replacement,2), EV_health_benefit)
  
  ggplot(EV_health_benefit,
         aes(x = replacement, color = Scenario, fill = Scenario)) +
    
    # Uncertainty band
    geom_ribbon(aes(ymin = `P2.5`, ymax = `P97.5`),
                alpha = 0.2,
                color = NA) +
    
    # Mean line
    geom_line(aes(y = Mean), linewidth = 1.2) +
    
    # Optional points (can remove if too busy)
    geom_point(aes(y = Mean), size = 2) +
    
    labs(
      x = "Fraction of vehicles replaced",
      y = paste0("Avoidable ",outcome_list[k]," rate\n(per 100,000)"),
      color = NULL,
      fill = NULL
    ) +
    
    scale_color_manual(values = c(
      "Gasoline" = "#00B050",
      "Diesel" = "#FF8C00"
    )) +
    
    scale_fill_manual(values = c(
      "Gasoline" = "#00B050",
      "Diesel" = "#FF8C00"
    )) +
    
    scale_x_continuous(
      breaks = seq(0, 1, by = 0.1),
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
  
  
  ggsave(paste0("Figures/Mitigation_EV_HB_rate_",pollutant_list[k],"_",outcome_list[k],"_region.png"),
         width = 8, height = 6, dpi = 600)
  
  
  
  
  ## Plot avoidable percentage for baseline risk
  EV_health_percent <- rbind(data.frame(HB_EV_gas_PIF_stat[[k]], Scenario = "Gasoline"),
                             data.frame(HB_EV_diesel_PIF_stat[[k]], Scenario = "Diesel"))
  EV_health_percent <- cbind(replacement = rep(EV_gas$replacement,2), EV_health_percent)
  
  
  ggplot(EV_health_percent, 
         aes(x = replacement, color = Scenario, fill = Scenario)) +
    
    # Uncertainty band
    geom_ribbon(aes(ymin = `P2.5`, ymax = `P97.5`),
                alpha = 0.2,
                color = NA) +
    
    # Mean line
    geom_line(aes(y = Mean), linewidth = 1.2) +
    
    # Optional points (can remove if too busy)
    geom_point(aes(y = Mean), size = 2) +
    
    labs(
      x = "EV replacement",
      y = paste0("Avoidable ",outcome_list[k]),
      color = NULL,
      fill = NULL
    ) +
    
    scale_color_manual(values = c(
      "Gasoline" = "#00B050",
      "Diesel" = "#FF8C00"
    )) +
    
    scale_fill_manual(values = c(
      "Gasoline" = "#00B050",
      "Diesel" = "#FF8C00"
    )) +
    
    scale_x_continuous(
      breaks = seq(0, 1, by = 0.1),
      labels = scales::percent_format(accuracy = 1),
      expand = expansion(mult = c(0, 0))
    ) +
    
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 0.1),
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
      plot.margin = margin(12, 18, 12, 12)
    )
  
  
  ggsave(paste0("Figures/Mitigation_EV_HB_perc_",pollutant_list[k],"_",outcome_list[k],"_region.png"),
         width = 8, height = 6, dpi = 600)
  
  
}







##### (4) Plot top 5 tracts that benefit most
bb  <- st_bbox(WEST_tract)
bb_sf <- st_as_sfc(bb)
bm <- get_tiles(bb_sf, provider = "OpenStreetMap", zoom = 14)
bm_df <- as.data.frame(bm, xy = TRUE)
bm_df$hex <- rgb(bm_df[[3]], bm_df[[4]], bm_df[[5]],
                 maxColorValue = 255)



ref_perc <- 0.2  ## The EV replacement percentage you want to emphasize to rank tracts that benefit most



## Gasoline vehicle replacement
for (k in 1:length(pair)) {
  print(k)
  
  df <- do.call(cbind, lapply(1:length(HB_EV_gas_number_stat[[k]]), 
                              function(x) HB_EV_gas_number_stat[[k]][[x]][,1]))
  df <- t(df)
  colnames(df) <- paste0("Tract",1:70)
  
  top5_tract <- order(df[which(EV_gas$replacement == ref_perc),], decreasing = TRUE)[1:5]
  
  df_long <- as.data.frame(df) %>%
    mutate(reduction = EV_gas$replacement) %>%
    pivot_longer(cols = -reduction,
                 names_to = "tract",
                 values_to = "value") %>%
    mutate(highlight = ifelse(tract %in% paste0("Tract",top5_tract), "Top 5", "Other tracts"))
  
  
  
  ggplot() +
    geom_line(
      data = df_long %>% filter(highlight == "Other tracts"),
      aes(x = reduction, y = value, group = tract),
      color = "grey70",
      linewidth = 0.6
    ) +
    geom_line(
      data = df_long %>% filter(highlight == "Top 5"),
      aes(x = reduction, y = value, group = tract),
      color = "#D94801",
      linewidth = 0.6
    ) +
    labs(
      x = "EV replacement",
      y = paste0("Avoidable ",outcome_list[k]," number"),
      color = NULL
    ) +
    scale_x_continuous(breaks = seq(0, 1, 0.1),
                       labels = function(x) paste0(x * 100, "%"),
                       expand = expansion(mult = c(0, 0))) +
    
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.justification = "center",
      legend.text = element_text(size = 12),
      axis.title.x = element_text(size = 14, margin = margin(t = 12)),
      axis.title.y = element_text(size = 14, margin = margin(r = 10)),
      axis.text = element_text(size = 12, color = "black"),
      axis.ticks = element_line(linewidth = 1, color = "black"),
      axis.ticks.length = unit(0.18, "cm"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.grid.major = element_line(color = "#D0D0D0", linewidth = 0.5, linetype = "dashed"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(12, 18, 12, 12)
    )
  
  ggsave(paste0("Figures/Mitigation_EV_HB_gas_number_",
                pollutant_list[k],"_",outcome_list[k],"_tract.png"),
         width = 8, height = 6, dpi = 600)
  
  
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
  
  
  ggsave(paste0("Figures/Mitigation_EV_HB_gas_number_",
                pollutant_list[k],"_",outcome_list[k],"_tract_top5_map.png"),
         width = 8, height = 6, dpi = 600)
  
  
}



## Diesel vehicle replacement
for (k in 1:length(pair)) {
  print(k)
  
  df <- do.call(cbind, lapply(1:length(HB_EV_diesel_number_stat[[k]]), 
                              function(x) HB_EV_diesel_number_stat[[k]][[x]][,1]))
  df <- t(df)
  colnames(df) <- paste0("Tract",1:70)
  
  top5_tract <- order(df[which(EV_diesel$replacement == ref_perc),], decreasing = TRUE)[1:5]
  
  df_long <- as.data.frame(df) %>%
    mutate(reduction = EV_diesel$replacement) %>%
    pivot_longer(cols = -reduction,
                 names_to = "tract",
                 values_to = "value") %>%
    mutate(highlight = ifelse(tract %in% paste0("Tract",top5_tract), "Top 5", "Other tracts"))
  
  
  
  ggplot() +
    geom_line(
      data = df_long %>% filter(highlight == "Other tracts"),
      aes(x = reduction, y = value, group = tract),
      color = "grey70",
      linewidth = 0.6
    ) +
    geom_line(
      data = df_long %>% filter(highlight == "Top 5"),
      aes(x = reduction, y = value, group = tract),
      color = "#D94801",
      linewidth = 0.6
    ) +
    labs(
      x = "EV replacement",
      y = paste0("Avoidable ",outcome_list[k]," number"),
      color = NULL
    ) +
    scale_x_continuous(breaks = seq(0, 1, 0.1),
                       labels = function(x) paste0(x * 100, "%"),
                       expand = expansion(mult = c(0, 0))) +
    
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.justification = "center",
      legend.text = element_text(size = 12),
      axis.title.x = element_text(size = 14, margin = margin(t = 12)),
      axis.title.y = element_text(size = 14, margin = margin(r = 10)),
      axis.text = element_text(size = 12, color = "black"),
      axis.ticks = element_line(linewidth = 1, color = "black"),
      axis.ticks.length = unit(0.18, "cm"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.grid.major = element_line(color = "#D0D0D0", linewidth = 0.5, linetype = "dashed"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(12, 18, 12, 12)
    )
  
  ggsave(paste0("Figures/Mitigation_EV_HB_diesel_number_",
                pollutant_list[k],"_",outcome_list[k],"_tract.png"),
         width = 8, height = 6, dpi = 600)
  
  
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
  
  
  ggsave(paste0("Figures/Mitigation_EV_HB_diesel_number_",
                pollutant_list[k],"_",outcome_list[k],"_tract_top5_map.png"),
         width = 8, height = 6, dpi = 600)
  
  
}


















########## 2. Low emission zone ##########

########## (2-1) Emission reduction ##########
LEZ <- data.frame(perc = seq(0, 1, by = 0.2))
LEZ <- LEZ %>%
  mutate(NO2 = perc * MOVES_data$NO2[2],
         NO2_reduction = NO2 / sum(MOVES_data$NO2),
         PM2.5 = perc * MOVES_data$PM2.5[2],
         PM2.5_reduction = PM2.5 / sum(MOVES_data$PM2.5))


LEZ_long <- LEZ %>%
  select(-NO2, -PM2.5) %>%
  rename(
    NO2 = NO2_reduction,
    PM2.5 = PM2.5_reduction
  ) %>%
  pivot_longer(cols = c(NO2, PM2.5),
               names_to = "Pollutant",
               values_to = "Reduction")


ggplot(LEZ_long, aes(x = perc, y = Reduction, color = Pollutant)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    x = "Diesel vehicle reduction",
    y = "Emission reduction",
    color = NULL
  ) +
  
  # Custom colors
  scale_color_manual(values = c(
    "PM2.5" = "#B22222",
    "NO2" = "#7A9CC6"
  )) +
  
  # Axis formatting
  scale_x_continuous(breaks = seq(0, 1, by = 0.2),
                     labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0))) +
  
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0))) +
  
  # Theme
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

ggsave(paste0("Figures/Mitigation_LEZ_emission_reduction.png"),
       width = 8, height = 6, dpi = 600)




########## (2-2) Concentration reduction ##########

LEZ_conc_result <- c()
LEZ_conc_list <- list()


for (pollut_name in c("PM2.5", "NO2")) {
  load(paste0("Results/Exposure_source_",pollut_name,"_tract.RData"))
  
  LEZ_conc_reduc <- PM_tract_annual_s_w_h[,c("tract","Fac1_perc")]
  
  for (k in 1:nrow(LEZ)) {
    LEZ_conc_reduc[[paste0("Case",k)]] <- LEZ_conc_reduc$Fac1_perc * LEZ[[paste0(pollut_name,"_reduction")]][k]
  }
  
  LEZ_conc_list[[pollut_name]] <- LEZ_conc_reduc
  
  LEZ_conc_result <- rbind(LEZ_conc_result,
                           data.frame(perc = LEZ$perc,
                                      Pollutant = pollut_name,
                                      Mean = apply(st_drop_geometry(LEZ_conc_reduc[,-(1:2)]), 2, mean),
                                      P2.5 = apply(st_drop_geometry(LEZ_conc_reduc[,-(1:2)]), 2, function(x) quantile(x, probs = 0.025)),
                                      P97.5 = apply(st_drop_geometry(LEZ_conc_reduc[,-(1:2)]), 2, function(x) quantile(x, probs = 0.975))))
}




ggplot(LEZ_conc_result, 
       aes(x = perc, color = Pollutant, fill = Pollutant)) +
  
  # Uncertainty band
  geom_ribbon(aes(ymin = P2.5, ymax = P97.5),
              alpha = 0.2,
              color = NA) +
  
  # Mean line
  geom_line(aes(y = Mean), linewidth = 1.2) +
  
  # Optional points (can remove if too busy)
  geom_point(aes(y = Mean), size = 2) +
  
  labs(
    x = "Diesel vehicle reduction",
    y = "Concentration reduction",
    color = NULL,
    fill = NULL,
    linetype = NULL
  ) +
  
  scale_color_manual(values = c(
    "PM2.5" = "#B22222",
    "NO2" = "#7A9CC6"
  )) +
  
  scale_fill_manual(values = c(
    "PM2.5" = "#B22222",
    "NO2" = "#7A9CC6"
  )) +
  
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.2),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0))
  ) +
  
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
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
    plot.margin = margin(12, 18, 12, 12)
  )



ggsave(paste0("Figures/Mitigation_LEZ_concentration_reduction.png"),
       width = 8, height = 6, dpi = 600)




########## (2-3) Health benefit ##########

##### (1) Calculate health benefits
HB_LEZ_rate <- HB_LEZ_rate_stat <- list() # Tract-level rate
HB_LEZ_number <- HB_LEZ_number_stat <- list() # Tract-level number

HB_LEZ_number_tot <- HB_LEZ_number_tot_stat <- list() # Regional-level number
HB_LEZ_rate_tot <- HB_LEZ_rate_tot_stat <- list() # Regional-level rate
HB_LEZ_PIF_stat <- list() # Regional-level percentage


for (k in 1:length(pair)) {
  
  HB_LEZ_rate[[k]] <- HB_LEZ_rate_stat[[k]] <- list()
  HB_LEZ_number[[k]] <- HB_LEZ_number_stat[[k]] <- list()
  
  HB_LEZ_number_tot[[k]] <- HB_LEZ_number_tot_stat[[k]] <- numeric(0)
  HB_LEZ_rate_tot[[k]] <- HB_LEZ_rate_tot_stat[[k]] <- numeric(0)
  
  
  for (s in 1:nrow(LEZ)) {
    print(paste0("k = ",k,"; s = ",s,"; Time = ", Sys.time()))
    
    ## Tract level
    HB_LEZ_rate[[k]][[s]] <- PIF_mitigation(pair[k], st_drop_geometry(LEZ_conc_list[[ER_data$pollutant[pair[k] %/% 100]]][,c("tract",paste0("Case",s))])) * 
      baseline_data$rate[pair[k] %% 100]
    
    HB_LEZ_number[[k]][[s]] <- HB_LEZ_rate[[k]][[s]] * WEST_tract_popu[[popu_group[k]]]/100000
    
    
    
    
    ## Region level
    HB_LEZ_number_tot[[k]] <- rbind(HB_LEZ_number_tot[[k]],
                                        apply(HB_LEZ_number[[k]][[s]], 2, sum))
    HB_LEZ_rate_tot[[k]] <- rbind(HB_LEZ_rate_tot[[k]],
                                      HB_LEZ_number_tot[[k]][s,] / (sum(WEST_tract_popu[[popu_group[k]]])/100000))
    
    
    ## Summary statistics
    HB_LEZ_rate_stat[[k]][[s]] <- t(apply(HB_LEZ_rate[[k]][[s]], 1, summary_stat))
    HB_LEZ_number_stat[[k]][[s]] <- t(apply(HB_LEZ_number[[k]][[s]], 1, summary_stat))
    colnames(HB_LEZ_rate_stat[[k]][[s]]) <- colnames(HB_LEZ_number_stat[[k]][[s]]) <- 
      c("Mean","SD","P2.5","P25","P50","P75","P97.5")
  }
  
  HB_LEZ_rate_tot_stat[[k]] <- t(apply(HB_LEZ_rate_tot[[k]], 1, summary_stat))
  HB_LEZ_number_tot_stat[[k]] <- t(apply(HB_LEZ_number_tot[[k]], 1, summary_stat))
  HB_LEZ_PIF_stat[[k]] <- HB_LEZ_rate_tot_stat[[k]] / baseline_data$rate[pair[k] %% 100]
  
  colnames(HB_LEZ_rate_tot_stat[[k]]) <- colnames(HB_LEZ_number_tot_stat[[k]]) <- 
    colnames(HB_LEZ_PIF_stat[[k]]) <- c("Mean","SD","P2.5","P25","P50","P75","P97.5")
  
  
}



save(HB_LEZ_rate, HB_LEZ_rate_stat,
     HB_LEZ_number, HB_LEZ_number_stat,
     HB_LEZ_number_tot, HB_LEZ_number_tot_stat,
     HB_LEZ_rate_tot, HB_LEZ_rate_tot_stat,
     HB_LEZ_PIF_stat,
     file = "Results/Mitigation_LEZ_health_benefit.RData")









##### (2) Plot curves for the whole study region
pollutant_list <- c(rep("PM2.5",4), rep("NO2",3))
outcome_list <- c("all-cause mortality","stroke incidence","lung cancer incidence","asthma incidence",
                  "all-cause mortality","lung cancer incidence","asthma incidence")

for (k in 1:length(pair)) {
  print(k)
  
  ## Plot avoidable risk rate
  df <- HB_LEZ_rate_tot_stat[[k]]
  df <- data.frame(perc = LEZ$perc, df)
  
  ggplot(df, aes(x = perc)) +
    
    # Uncertainty band
    geom_ribbon(aes(ymin = `P2.5`, ymax = `P97.5`),
                fill = "#00B050",
                alpha = 0.2,
                color = NA) +
    
    # Mean line
    geom_line(aes(y = Mean), linewidth = 1.2, color = "#00B050") +
    
    # Optional points (can remove if too busy)
    geom_point(aes(y = Mean), size = 2, color = "#00B050") +
    
    labs(
      x = "Diesel vehicle reduction",
      y = paste0("Avoidable ",outcome_list[k]," rate\n(per 100,000)"),
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
      plot.margin = margin(12, 18, 12, 12)
    )
  
  
  ggsave(paste0("Figures/Mitigation_LEZ_HB_rate_",pollutant_list[k],"_",outcome_list[k],"_region.png"),
         width = 8, height = 6, dpi = 600)
  
  
  
  
  ## Plot avoidable percentage for baseline risk
  df <- HB_LEZ_PIF_stat[[k]]
  df <- data.frame(perc = LEZ$perc, df)
  
  ggplot(df, aes(x = perc)) +
    
    # Uncertainty band
    geom_ribbon(aes(ymin = `P2.5`, ymax = `P97.5`),
                fill = "#FF8C00",
                alpha = 0.2,
                color = NA) +
    
    # Mean line
    geom_line(aes(y = Mean), linewidth = 1.2, color = "#FF8C00") +
    
    # Optional points (can remove if too busy)
    geom_point(aes(y = Mean), size = 2, color = "#FF8C00") +
    
    labs(
      x = "Diesel vehicle reduction",
      y = paste0("Avoidable ",outcome_list[k]),
      color = NULL,
      fill = NULL
    ) +
    
    scale_x_continuous(
      breaks = seq(0, 1, by = 0.2),
      labels = scales::percent_format(accuracy = 1),
      expand = expansion(mult = c(0, 0))
    ) +
    
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 0.1),
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
      plot.margin = margin(12, 18, 12, 12)
    )
  
  
  ggsave(paste0("Figures/Mitigation_LEZ_HB_perc_",pollutant_list[k],"_",outcome_list[k],"_region.png"),
         width = 8, height = 6, dpi = 600)
  
}






##### (3) Plot top 5 tracts

ref_perc <- 0.2

for (k in 1:length(pair)) {
  print(k)
  
  df <- do.call(cbind, lapply(1:length(HB_LEZ_number_stat[[k]]), 
                              function(x) HB_LEZ_number_stat[[k]][[x]][,1]))
  df <- t(df)
  colnames(df) <- paste0("Tract",1:70)
  
  top5_tract <- order(df[which(LEZ$perc == ref_perc),], decreasing = TRUE)[1:5]
  
  df_long <- as.data.frame(df) %>%
    mutate(reduction = LEZ$perc) %>%
    pivot_longer(cols = -reduction,
                 names_to = "tract",
                 values_to = "value") %>%
    mutate(highlight = ifelse(tract %in% paste0("Tract",top5_tract), "Top 5", "Other tracts"))
  
  
  
  ggplot() +
    geom_line(
      data = df_long %>% filter(highlight == "Other tracts"),
      aes(x = reduction, y = value, group = tract),
      color = "grey70",
      linewidth = 0.6
    ) +
    geom_line(
      data = df_long %>% filter(highlight == "Top 5"),
      aes(x = reduction, y = value, group = tract),
      color = "#D94801",
      linewidth = 0.6
    ) +
    labs(
      x = "Diesel vehicle reduction",
      y = paste0("Avoidable ",outcome_list[k]," number"),
      color = NULL
    ) +
    scale_x_continuous(breaks = seq(0, 1, 0.2),
                       labels = function(x) paste0(x * 100, "%"),
                       expand = expansion(mult = c(0, 0))) +
    
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.justification = "center",
      legend.text = element_text(size = 12),
      axis.title.x = element_text(size = 14, margin = margin(t = 12)),
      axis.title.y = element_text(size = 14, margin = margin(r = 10)),
      axis.text = element_text(size = 12, color = "black"),
      axis.ticks = element_line(linewidth = 1, color = "black"),
      axis.ticks.length = unit(0.18, "cm"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.grid.major = element_line(color = "#D0D0D0", linewidth = 0.5, linetype = "dashed"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(12, 18, 12, 12)
    )
  
  ggsave(paste0("Figures/Mitigation_LEZ_HB_number_",
                pollutant_list[k],"_",outcome_list[k],"_tract.png"),
         width = 8, height = 6, dpi = 600)
  
  
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
  
  
  ggsave(paste0("Figures/Mitigation_LEZ_HB_number_",
                pollutant_list[k],"_",outcome_list[k],"_tract_top5_map.png"),
         width = 8, height = 6, dpi = 600)
  
  
}







########## 3. Congestion pricing ##########

PIF_cal <- function(risk_outcome, C_delta, type) {
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
  TMREL <- rep(0, 1000)
  
  PIF <- lapply(unique(exposure$tract),
                function(i) {
                  tract_exp <- exposure %>% filter(tract == i) %>% select(-c(tract, run_id))
                  RR <- RR0^tract_exp
                  if (type == "Absolute") {
                    RR_limit <- RR0^pmax(as.matrix(tract_exp - C_delta), TMREL)
                  } else if (type == "Percent") {
                    RR_limit <- RR0^pmax(as.matrix(tract_exp * (1 - C_delta)), TMREL)
                  }
                  
                  PIF <- apply(RR - RR_limit, 1, mean) / apply(RR, 1, mean)
                  return(PIF)
                })
  PIF <- as.data.frame(do.call(rbind, PIF))
  
  return(PIF)
}


## Health benefit rate by tracts and draws
HB_CP <- list()
for (k in 1:length(pair)) {
  print(paste0("k = ",k,"; Time = ", Sys.time()))
  HB_CP[[k]] <- PIF_cal(pair[k], 0.22, "Percent") * baseline_data$rate[pair[k] %% 100]
}


## Health benefit rate for the whole study region
HB_CP_tot <- c()
for (k in 1:length(HB_CP)) {
  dat <- HB_CP[[k]]
  dat[dat < 0] <- 0
  dat <- cbind(dat, popu = WEST_tract_popu[,popu_group[k]])
  dat[,1:1000] <- dat[,1:1000] * dat$popu/100000
  dat_summary <- apply(dat[,1:1000],2,sum) / (sum(dat$popu)/100000)
  HB_CP_tot <- rbind(HB_CP_tot,
                     c(Mean = mean(dat_summary),
                       LCI = quantile(dat_summary, probs = 0.025),
                       UCI = quantile(dat_summary, probs = 0.975)))
}

HB_CP_perc <- HB_CP_tot / baseline_data$rate[pair %% 100]

save(HB_CP, HB_CP_tot, HB_CP_perc,
     file = "Results/Mitigation_CP_health_benefit.RData")


