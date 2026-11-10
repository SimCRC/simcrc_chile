# *****************************************************************************
#
# Script: population_colonoscopy_demand.R
#
# Purpose: National colonoscopy demand in one year, inflating the single-cohort
#          SimCRC results to the FONASA population (Fernando, 22 Oct 2026).
#          (1) Population table: at each age a, weight w_a = FONASA
#          beneficiaries at a / people alive at a in the simulated cohort
#          (no screening). (2) Demand = sum over ages of the cohort's
#          colonoscopies at a times w_a, split into screening, follow-up
#          (after a positive stool test) and surveillance.
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
#     Rscript analysis/cea/population_colonoscopy_demand.R
#   Reads the v0.15.0 Base run (outputs/cea/2026Chile/Base) and FONASA 2025
#   beneficiaries by single age (5-year bands split evenly). Perfect
#   adherence, as in the CEA; every age is at its place in the strategy's
#   schedule, so this is a ballpark, not a projection. Writes
#   outputs/cea/population_demand/.
#
# *****************************************************************************

# 1. Packages ----------------------------------------------------------------

suppressMessages(library(data.table))

dir_raw <- "outputs/cea/2026Chile/Base/RawModelOutput_SimCRC_R"
dir_out <- "outputs/cea/population_demand"
dir.create(dir_out, showWarnings = FALSE)

# 2. Population table: FONASA people per simulated person, by age ----------

dt_fonasa <- fread("data-raw/fonasa_population_by_single_age_2025.csv")[
  , .(age, n_fonasa = n_beneficiarios)]
dt_cohort <- fread(file.path(dir_raw, "2026Chile_NoScreening.csv"),
                   skip = "age,n_alive")[, .(age, n_cohort = n_alive)]
dt_population <- merge(dt_fonasa, dt_cohort, by = "age")
dt_population[, weight := n_fonasa / n_cohort]

# 3. Cohort colonoscopies by strategy and age ------------------------------

dt_col <- rbindlist(lapply(list.files(dir_raw, full.names = TRUE), function(f) {
  dt <- fread(f, skip = "age,n_alive")
  dt[, .(strategy = sub("^2026Chile_(COL|FIT)(..)(..)q(.*)[.]csv$",
                        "\\1_\\2_\\3_\\4", basename(f)),
         age,
         screening    = n_pos_scrcol + n_neg_scrcol,
         followup     = n_pos_fucol + n_neg_fucol,
         surveillance = n_pos_survcol + n_neg_survcol)]
}))
dt_col[strategy == "2026Chile_NoScreening.csv", strategy := "NoScreening"]

# 4. Inflate to the FONASA population ---------------------------------------

v_type <- c("screening", "followup", "surveillance")
dt_age <- merge(dt_col, dt_population[, .(age, n_fonasa, weight)], by = "age")
dt_age[, (v_type) := lapply(.SD, `*`, weight), .SDcols = v_type]
dt_age[, total := screening + followup + surveillance]

dt_strategy <- dt_age[, lapply(.SD, sum), by = strategy,
                      .SDcols = c(v_type, "total")]
dt_icer <- fread("outputs/cea/ce_results/df_icer_all_strategies.csv")
dt_strategy[dt_icer, efficient := i.Status == "ND", on = .(strategy = Strategy)]
setorder(dt_strategy, total)

# 5. Write outputs -------------------------------------------------------------

fwrite(dt_population, file.path(dir_out, "population_weights_by_age.csv"))
fwrite(dt_age[, .(strategy, age, n_fonasa, screening, followup, surveillance,
                  total)],
       file.path(dir_out, "colonoscopy_demand_by_age.csv"))
fwrite(dt_strategy, file.path(dir_out, "colonoscopy_demand_by_strategy.csv"))

cat("Strategies:", nrow(dt_strategy), "| FONASA ages", min(dt_age$age), "-",
    max(dt_age$age), "\n\nEfficient strategies, colonoscopies per year:\n")
print(dt_strategy[efficient == TRUE, lapply(.SD, round), by = strategy,
                  .SDcols = c("screening", "followup", "surveillance",
                              "total")])
