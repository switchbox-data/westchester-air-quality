# Shared helpers used by multiple analysis scripts.
# Source from each script that needs them: source("Code/_helpers.R")
#
# Assumes the caller has already loaded dplyr and set up:
#   ER_data, PM_tract_boot, NO2_tract_boot

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
