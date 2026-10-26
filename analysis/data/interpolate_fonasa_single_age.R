# *****************************************************************************
#
# Script: interpolate_fonasa_single_age.R
#
# Purpose: Interim single-year-of-age FONASA population (ages 40-100) for the
#          P_a term of the colonoscopy demand model, interpolated linearly
#          from the 5-year bands in fonasa_population_by_age_group_2025.csv.
#          Replace once Nicolas sends a FONASA breakdown by single year of age.
#
# Author: Jorge Roa and Claudia Seguin
#
# Email: jorgeroa@stanford.edu/cseguin@stanford.edu
#
# Date Created: 26 Oct 2026
#
# *****************************************************************************
#
# Notes:
#   Method: each band's average count per year (n / 5) is placed at the band
#   midpoint (42, 47, ..., 97) and single ages are read off the straight line
#   between neighbouring midpoints. Within each band the interpolated values
#   are then rescaled so they sum exactly to the band's FONASA count.
#
#   Age 100: the open "Más de 99 años" band is treated as spanning 100-104
#   (knot at 102, value n / 5), so age 100 is interpolated between the 97 and
#   102 knots. It is NOT rescaled, because only part of the 100+ band is kept.
#   Extrapolating the 90-94 -> 95-99 trend instead gives a negative count.
#
#   "S.I." (unknown age, ~0.08% of beneficiaries) is dropped.
#
# *****************************************************************************

library(data.table)

path_fonasa_groups <- "data-raw/fonasa_population_by_age_group_2025.csv"
path_fonasa_out    <- "data-raw/fonasa_population_by_single_age_40_100_2025.csv"

age_min <- 40
age_max <- 100

dt_fonasa <- fread(path_fonasa_groups, encoding = "UTF-8")

# Parse band bounds; "Más de 99 años" -> 100-104 (see Notes), drop "S.I."
dt_fonasa <- dt_fonasa[grupo_edad != "S.I."]
dt_fonasa[grepl("^\\d+ a \\d+", grupo_edad), `:=`(
  age_lower = as.integer(sub("^(\\d+) a .*", "\\1", grupo_edad)),
  age_upper = as.integer(sub("^\\d+ a (\\d+) .*", "\\1", grupo_edad))
)]
dt_fonasa[grepl("^Más de", grupo_edad), `:=`(age_lower = 100L, age_upper = 104L)]
stopifnot(!anyNA(dt_fonasa$age_lower), !anyNA(dt_fonasa$age_upper))

# Knots: average count per single year at each band midpoint
dt_fonasa[, `:=`(
  age_mid    = (age_lower + age_upper) / 2,
  n_per_year = n_beneficiarios / (age_upper - age_lower + 1)
)]

# Linear interpolation between knots for ages 40-100
dt_single_age <- data.table(age = age_min:age_max)
dt_single_age[, n_linear := approx(
  x = dt_fonasa$age_mid, y = dt_fonasa$n_per_year, xout = age
)$y]
dt_single_age[dt_fonasa, on = .(age >= age_lower, age <= age_upper),
              `:=`(grupo_edad = i.grupo_edad, n_group = i.n_beneficiarios)]

# Rescale within each complete band so single ages sum to the band total
dt_single_age[, n_beneficiarios := n_linear * n_group / sum(n_linear),
              by = grupo_edad]
dt_single_age[age == 100, n_beneficiarios := n_linear]

stopifnot(
  !anyNA(dt_single_age$n_beneficiarios),
  all(dt_single_age$n_beneficiarios > 0)
)

dt_single_age <- dt_single_age[, .(grupo_edad, age, n_beneficiarios)]

fwrite(dt_single_age, path_fonasa_out)
