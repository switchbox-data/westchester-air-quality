# COMPUTE-ONLY extraction of section "1A. Vehicle electrification NEW: gasoline
# or diesel" from Code/04_Mitigation_strategies_NEW.R (lines 928-1534).
#
# Computes region- and tract-level health benefits for two SEPARATE
# electrification scenarios:
#   (1) Gasoline vehicle replacement only (fraction of the gasoline fleet
#       replaced by EVs; diesel fleet untouched).
#   (2) Diesel vehicle replacement only (fraction of the diesel fleet
#       replaced by EVs; gasoline fleet untouched).
# Both use the replacement grid c(0, 0.05, 0.1, ..., 0.9, 1). This is NOT the
# same as the legacy "diesel-first" scenario (Code/04_Mitigation_strategies_NEW.R
# section 1, saved to Results/Mitigation_EV_diesel_health_benefit.RData), which
# replaces the diesel fleet first and then spills over into gasoline/CNG/ethanol
# replacement once diesel is exhausted. To avoid clobbering those legacy
# results, this script's outputs are saved under distinct filenames:
#   Results/Mitigation_EV_gasonly_health_benefit.RData
#   Results/Mitigation_EV_dieselonly_health_benefit.RData
#
# Reads ONLY cached Results/ and Data/ inputs (Exposure_source_*_tract.RData,
# Exposure_boot_*_tract.RData, Westchester_tracts_selected.shp,
# Westchester_tract_population_by_age.RData) — no recomputation upstream, no
# Census API. NO plotting is done here (see Code/04_Mitigation_strategies_NEW.R
# for the ggplot figures); this script only produces the saved .RData objects
# for downstream figure scripts to consume.
#
# Run inside the project container with Code, Data, Results bind-mounted:
#
#   docker run --rm \
#     -v "$PWD/Code":/project/Code -v "$PWD/Data":/project/Data \
#     -v "$PWD/Results":/project/Results \
#     -w /project westchester-air-quality:latest Rscript Code/05_Mitigation_gas_diesel_compute.R
#
# Takes on the order of an hour to run (bootstrap PIF over 1000 draws x 7
# pollutant-outcome pairs x 14 replacement fractions x 2 scenarios).

library(dplyr)
library(tidyr)
library(sf)


########## 1A. Vehicle electrification: gasoline or diesel ##########

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




EV_diesel <- data.frame(replacement = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
EV_diesel <- EV_diesel %>%
  mutate(NO2 = replacement * MOVES_data$VMT[2] * (MOVES_data$EF_NO2[2] - MOVES_data$EF_NO2[5]),
         NO2_reduction = NO2 / sum(MOVES_data$NO2),
         PM2.5 = replacement * MOVES_data$VMT[2] * (MOVES_data$EF_PM2.5[2] - MOVES_data$EF_PM2.5[5]),
         PM2.5_reduction = PM2.5 / sum(MOVES_data$PM2.5))




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




########## (1-3) Health benefit ##########

## Load tract-level annual average exposure with bootstrap uncertainties from Code 01
## NOTE: Exposure_boot_NO2_tract.RData loads an object named PM_tract_boot (not
## NO2_tract_boot, despite what a comment in the source section claimed) -- it
## must be renamed BEFORE loading the PM2.5 boot file, which also creates an
## object named PM_tract_boot and would otherwise silently overwrite the NO2
## data. This mirrors the (correct) pattern in section 1 of
## Code/04_Mitigation_strategies_NEW.R (lines 273-276).
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
     EV_gas, EV_gas_conc_list,
     file = "Results/Mitigation_EV_gasonly_health_benefit.RData")




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
     EV_diesel, EV_diesel_conc_list,
     file = "Results/Mitigation_EV_dieselonly_health_benefit.RData")
