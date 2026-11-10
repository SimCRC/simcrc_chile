# *****************************************************************************
#
# Script: multicohort_population.R
#
# Purpose: Show the multi-cohort life tables in use: one female population
#          per birth cohort from simcrc::get_dt_population(), stacked, with
#          women alive at 45+ by calendar year and survival from other
#          causes for five cohorts. No natural history or calibration needed.
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
#   Run from the simcrc_chile root, after build_multicohort_inputs.R:
#     Rscript analysis/data/multicohort_population.R
#   100,000 women per cohort, unweighted until the population projections are
#   added. Writes outputs/figures/multicohort_survival_by_cohort.png.
#
# *****************************************************************************

# 1. Packages ----------------------------------------------------------------

suppressMessages({
  library(simcrc)
  library(data.table)
  library(ggplot2)
})

# 2. Cohort life tables -------------------------------------------------------

dt_lt <- fread("data-raw/multicohort/cohort_life_tables.csv")[sex == "female"]
v_program_years <- as.integer(names(
  fread("data-raw/multicohort/birth_year_grid.csv"))[-1])

# 3. Population, one per birth cohort -----------------------------------------

set.seed(2026)
dt_pop <- dt_lt[, get_dt_population(year = birth_year, byear = birth_year,
                                    p_female = 1, n_pop = 1e5,
                                    dt_life_table_F = .SD),
                by = birth_year, .SDcols = c("age", "p_HD", "p_cum_HD")]

# 4. Women alive at 45+ by calendar year -------------------------------------

dt_alive <- data.table(year = v_program_years)
dt_alive[, n_alive_45plus := sapply(year, function(y) dt_pop[
  byear + 45 <= y & byear + floor(age_death_oc) >= y, .N])]

# 5. Survival figure ----------------------------------------------------------

v_show <- c(1930, 1940, 1950, 1960, 1970)
dt_surv <- dt_pop[byear %in% v_show,
                  .(age = seq(dt_lt[birth_year == .BY$byear, min(age)], 100)),
                  by = byear]
dt_surv[, alive := dt_pop[byear == .BY$byear, mean(age_death_oc >= age)],
        by = .(byear, age)]

fig_surv <- ggplot(dt_surv, aes(age, alive, colour = factor(byear))) +
  geom_line(linewidth = 0.9) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "Survival from other causes, by birth cohort (INE Chile, women)",
       subtitle = paste0("Each cohort starts at its age in ",
                         v_program_years[1], " (or 45). ",
                         "Years after 2050 use the 2050 rates."),
       x = "age", y = "proportion alive", colour = "birth year") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank())
ggsave("outputs/figures/multicohort_survival_by_cohort.png", fig_surv,
       width = 9, height = 5, dpi = 300, bg = "white")

print(dt_alive)

# 6. Check simulated life expectancy against INE ------------------------------

dt_ine <- fread("data-raw/multicohort/ine_national_rates.csv")
dt_ine[, ex_ine_100 := (Tx - Tx[age == 100]) / lx, by = .(year, sex)]
dt_chk <- CJ(year = seq(1995, 2050, 5), sex = c("female", "male"),
             age = c(45, 65))
dt_chk[, ex_simcrc := {
  d <- dt_ine[year == .BY$year & sex == .BY$sex & age >= .BY$age]
  d[age == 100, mortality_rate := 1000]
  d[, p_cum_HD := 1 - exp(-cumsum(mortality_rate))]
  d[, p_HD := diff(c(0, p_cum_HD))]
  p <- get_dt_population(year = .BY$year, byear = .BY$year,
                         p_female = as.numeric(.BY$sex == "female"),
                         n_pop = 1e5,
                         dt_life_table_F = if (.BY$sex == "female") d,
                         dt_life_table_M = if (.BY$sex == "male") d)
  mean(p$age_death_oc) - .BY$age
}, by = .(year, sex, age)]
dt_chk[dt_ine, ex_ine_100 := i.ex_ine_100, on = .(year, sex, age)]
cat("Simulated minus INE years left to 100: from",
    round(min(dt_chk$ex_simcrc - dt_chk$ex_ine_100), 3), "to",
    round(max(dt_chk$ex_simcrc - dt_chk$ex_ine_100), 3), "\n")

fig_ex <- ggplot(dt_ine[age %in% c(45, 65)],
                 aes(year, ex_ine_100, colour = sex)) +
  geom_line(linewidth = 0.9) +
  geom_point(data = dt_chk, aes(year, ex_simcrc), size = 2.2) +
  facet_wrap(~ age, scales = "free_y",
             labeller = labeller(
               age = function(a) paste("years left at", a, "(to 100)"))) +
  labs(title = "Life expectancy in Chile: INE published vs SimCRC simulated",
       subtitle = paste0("Lines: INE years left up to age 100, ",
                         "(T(x) - T(100)) / l(x). Points: mean years left in ",
                         "100,000 deaths sampled from our tables."),
       x = "calendar year", y = "years left up to age 100", colour = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank())
ggsave("outputs/figures/multicohort_life_expectancy_check.png", fig_ex,
       width = 10, height = 4.8, dpi = 300, bg = "white")
