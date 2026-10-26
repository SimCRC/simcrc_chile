# *****************************************************************************
#
# Script: pop_weights_chile.R
#
# Purpose: Single-age FONASA population (ages 40-100) and single-age weights
#          to rescale the uspstf_summary outputs to the Chilean population.
#
# Author: Jorge Roa and Claudia Seguin
#
# Email: jorgeroa@stanford.edu/cseguin@stanford.edu
#
# Date Created: 26 Oct 2026
#
# *****************************************************************************

# Notes:
# 1. Prepare the pop projections
# convert from 5-year age groups to single-year age using a linear interpolation for now
# 2. Create single-age weights for the Chile population based on Fonasa data
#    (applied to the uspstf_summary files in a separate script)

library(data.table)

# *****************************************************************************
#### 1. FONASA population by single year of age (40-100) ####
# *****************************************************************************
# Each 5-year band's average count per year (n / 5) is placed at the band
# midpoint (42, 47, ..., 97); single ages are read off the straight line
# between neighbouring midpoints:
#   n(a) = y_k + (a - m_k) * (y_k+1 - y_k) / (m_k+1 - m_k),  m_k <= a <= m_k+1
# Within each band the values are then rescaled so they sum to the band total:
#   n_adj(a) = n(a) * N_band / sum(n(a) for a in band)
# Age 100: the open "Más de 99 años" band is treated as 100-104 (knot at 102)
# and age 100 is interpolated between the 97 and 102 knots, not rescaled.
# Extending the 90-94 -> 95-99 trend instead gives a negative count.
# "S.I." (unknown age, ~0.08%) is dropped.

age_min <- 40
age_max <- 100

dt_fonasa_groups <- fread("data-raw/fonasa_population_by_age_group_2025.csv",
                          encoding = "UTF-8")

# Band bounds; "Más de 99 años" -> 100-104, drop "S.I."
dt_fonasa_groups <- dt_fonasa_groups[grupo_edad != "S.I."]
dt_fonasa_groups[grepl("^\\d+ a \\d+", grupo_edad), `:=`(
  age_lower = as.integer(sub("^(\\d+) a .*", "\\1", grupo_edad)),
  age_upper = as.integer(sub("^\\d+ a (\\d+) .*", "\\1", grupo_edad))
)]
dt_fonasa_groups[grepl("^Más de", grupo_edad),
                 `:=`(age_lower = 100L, age_upper = 104L)]
stopifnot(!anyNA(dt_fonasa_groups$age_lower),
          !anyNA(dt_fonasa_groups$age_upper))

# Knots: average count per single year at each band midpoint
dt_fonasa_groups[, `:=`(
  age_mid    = (age_lower + age_upper) / 2,
  n_per_year = n_beneficiarios / (age_upper - age_lower + 1)
)]

# Linear interpolation between knots
dt_fonasa <- data.table(age = age_min:age_max)
dt_fonasa[, n_linear := approx(x = dt_fonasa_groups$age_mid,
                               y = dt_fonasa_groups$n_per_year,
                               xout = age)$y]
dt_fonasa[dt_fonasa_groups, on = .(age >= age_lower, age <= age_upper),
          `:=`(grupo_edad = i.grupo_edad, n_group = i.n_beneficiarios)]

# Rescale within each complete band so single ages sum to the band total
dt_fonasa[, n_beneficiarios := n_linear * n_group / sum(n_linear),
          by = grupo_edad]
dt_fonasa[age == 100, n_beneficiarios := n_linear]
stopifnot(!anyNA(dt_fonasa$n_beneficiarios),
          all(dt_fonasa$n_beneficiarios > 0))

dt_fonasa <- dt_fonasa[, .(grupo_edad, age, n_beneficiarios)]

# *****************************************************************************
#### 2. Single-age weights: Chilean population / simulated n_alive ####
# *****************************************************************************
# Goal: rescale the uspstf_summary outputs (RawModelOutput_SimCRC_R/*.csv, one
# row per age 40-100) from the simulated cohort to the actual Chilean (FONASA)
# population.
#
#   weight(a) = n_beneficiarios_FONASA(a) / n_alive_simulated(a)
#
# - The denominator is n_alive from the No Screening uspstf_summary, so every
#   strategy gets the same weight at a given age. This keeps differences
#   between strategies (e.g., CRC deaths averted) intact. Dividing by each
#   strategy's own n_alive would force all strategies to the same n_alive and
#   erase those differences.
# - The weights apply to ALL uspstf_summary outcome columns except age.
#   After weighting, n_alive in the No Screening file equals the FONASA
#   population at each age.
# - This script ONLY creates and saves the weights. A separate script, run
#   before process_uspstf_output(), applies them to every uspstf_summary file.
#   CRC cases, CRC deaths and the other outcomes are summarized only AFTER the
#   uspstf_summary files are updated to reflect the actual Chilean population.

path_uspstf_noscreening <- paste0(
  "outputs/cea/2026Chile/Base/RawModelOutput_SimCRC_R/",
  "2026Chile_NoScreening.csv"
)
path_weights_out <- "outputs/cea/2026Chile/Base/pop_weights_chile_single_age.csv"

dt_nh <- fread(path_uspstf_noscreening, skip = "age,n_alive")

dt_pop_weights <- merge(dt_fonasa, dt_nh[, .(age, n_alive)], by = "age")
stopifnot(nrow(dt_pop_weights) == age_max - age_min + 1)
dt_pop_weights[, weight := n_beneficiarios / n_alive]

fwrite(dt_pop_weights[, .(age, n_beneficiarios, n_alive, weight)],
       path_weights_out)