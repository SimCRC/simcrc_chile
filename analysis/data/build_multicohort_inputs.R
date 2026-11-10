# *****************************************************************************
#
# Script: build_multicohort_inputs.R
#
# Purpose: Inputs for a multi-cohort SimCRC run in Chile, following Claudia
#          Seguin's Mexico City method. (1) Birth-year grid: birth year =
#          program year - age, for ages 45 to 99. (2) One life table per
#          birth cohort and sex, read along the diagonal of the INE period
#          tables: cohort b at age a gets the death rate of year b + a. Each
#          table starts at the cohort's age in the first program year, or 45.
#
# Author: Jorge Roa
#
# Email: jorgeroa@stanford.edu
#
# ISS: #28   Version: #24
#
# Date Created: 21 October 2026
#
# *****************************************************************************
#
# Notes:
#   Run from the simcrc_chile root:
#     Rscript analysis/data/build_multicohort_inputs.R
#   Program years 2000-2019 are PLACEHOLDERS until Fernando and Nicolás pick
#   the program. Rates: national (País) m(x) of the INE tables 1992-2050;
#   years after 2050 reuse 2050, and age 100 gets rate 1000 so everyone dies
#   by 100. Writes to data-raw/multicohort/, with the INE national rates
#   with l(x), T(x) and e(x), in ine_national_rates.csv for the life
#   expectancy check.
#
# *****************************************************************************

# 1. Packages ----------------------------------------------------------------

suppressMessages({
  library(data.table)
  library(readxl)
})

# 2. Program settings (placeholders) -----------------------------------------

first_program_year <- 2000
last_program_year  <- 2019
first_age          <- 45
last_age           <- 100
dir_out            <- "data-raw/multicohort"
dir.create(dir_out, showWarnings = FALSE)

# 3. Birth-year grid ----------------------------------------------------------

dt_grid <- CJ(age = first_age:(last_age - 1),
              year = first_program_year:last_program_year)
dt_grid[, birth_year := year - age]
fwrite(dcast(dt_grid, age ~ year, value.var = "birth_year"),
       file.path(dir_out, "birth_year_grid.csv"))

# 4. INE national death rates -------------------------------------------------

dt_ine <- as.data.table(suppressMessages(read_excel(
  "data-raw/life table/tablas-de-mortalidad-de-chile-1992-2050_Ch.xlsx",
  sheet = "BD Tablas de Mortalidad", skip = 1)))[, c(2:6, 8, 11, 13)]
setnames(dt_ine, c("region", "year", "sex", "age", "mortality_rate",
                   "lx", "Tx", "ex"))
dt_ine <- dt_ine[region == "País", .(
  year = as.integer(year),
  sex  = fifelse(sex == "Mujer", "female", "male"),
  age  = as.integer(sub("-.*", "", age)),
  mortality_rate = as.numeric(mortality_rate),
  lx = as.numeric(lx),
  Tx = as.numeric(Tx),
  ex_ine = as.numeric(ex)
)]

# 5. Cohort life tables -------------------------------------------------------

dt_lt <- CJ(sex = c("female", "male"),
            birth_year = min(dt_grid$birth_year):max(dt_grid$birth_year),
            age = first_age:last_age)

dt_lt <- dt_lt[age >= pmax(first_age, first_program_year - birth_year)]

dt_lt[, calendar_year := birth_year + age]
dt_lt[, year := pmin(calendar_year, max(dt_ine$year))]
dt_lt[dt_ine, mortality_rate := i.mortality_rate, on = .(sex, age, year)]
dt_lt[age == last_age, mortality_rate := 1000]
dt_lt[, p_cum_HD := 1 - exp(-cumsum(mortality_rate)), by = .(sex, birth_year)]
dt_lt[, p_HD := diff(c(0, p_cum_HD)), by = .(sex, birth_year)]
dt_lt[, year := NULL]



stopifnot(!anyNA(dt_lt$mortality_rate))

# 6. Write outputs -------------------------------------------------------------

fwrite(dt_ine, file.path(dir_out, "ine_national_rates.csv"))
fwrite(dt_lt, file.path(dir_out, "cohort_life_tables.csv"))
fwrite(dcast(dt_lt[sex == "female"], age ~ birth_year,
             value.var = "mortality_rate"),
       file.path(dir_out, "cohort_mortality_rates_female_wide.csv"))

cat("Birth cohorts", min(dt_lt$birth_year), "-", max(dt_lt$birth_year),
    "|", nrow(dt_lt), "life-table rows\n")
