# *****************************************************************************
#
# Script: population_natural_history.R
#
# Purpose: Minimal population model (Fernando, 22 Oct 2026), natural history
#          only. Cohort results per 1,000 people alive at each age, times the
#          FONASA inflation factor at that age (beneficiaries / 1,000), give
#          the number of people in Chile at that age in one year.
#
# Author: Jorge Roa
#
# Email: jorgeroa@stanford.edu
#
# ISS: #28   Version: #24
#
# Date Created: 23 October 2026
#
# *****************************************************************************
#
# Notes:
#   Run from the simcrc_chile root, after
#   analysis/data/wrangle_fonasa_population.R:
#     Rscript analysis/cea/population_natural_history.R
#   Reads the v0.15.0 Base no-screening run and FONASA 2025 beneficiaries by
#   single age. Writes outputs/cea/population_demand/natural_history_by_age.csv.
#
# *****************************************************************************

suppressMessages(library(data.table))

dt_fonasa <- fread("data-raw/fonasa_population_by_single_age_2025.csv")[
  , .(age, inflation_factor)]
dt_nh <- fread(paste0("outputs/cea/2026Chile/Base/RawModelOutput_SimCRC_R/",
                      "2026Chile_NoScreening.csv"), skip = "age,n_alive")

dt_nh <- dt_nh[, .(age,
                   crc_cases_per_1000  = (n_CRCstageI_sym + n_CRCstageII_sym +
                     n_CRCstageIII_sym + n_CRCstageIV_sym) / n_alive * 1000,
                   crc_deaths_per_1000 = n_CRC_dth / n_alive * 1000)]
dt_nh <- merge(dt_nh, dt_fonasa, by = "age")
dt_nh[, `:=`(crc_cases_chile  = crc_cases_per_1000 * inflation_factor,
             crc_deaths_chile = crc_deaths_per_1000 * inflation_factor)]

dir.create("outputs/cea/population_demand", showWarnings = FALSE)
fwrite(dt_nh, "outputs/cea/population_demand/natural_history_by_age.csv")

print(dt_nh[age %% 10 == 0])
cat("Chile, one year, ages 40-100: CRC cases", round(sum(dt_nh$crc_cases_chile)),
    "| CRC deaths", round(sum(dt_nh$crc_deaths_chile)), "\n")
