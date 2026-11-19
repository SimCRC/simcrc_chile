# *****************************************************************************
#
# Script: population_colonoscopy_demand.R
#
# Purpose: National colonoscopy demand in one year for the FONASA population
#          (Fernando, 22 Oct 2026). Sums the colonoscopies in the
#          FONASA-weighted uspstf_summary files (Claudia's weights, section
#          3.5 of cea_analysis_chile_parallel.R), split into screening,
#          follow-up (after a positive stool test) and surveillance.
#          Symptomatic diagnoses are reported but left out of the total.
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
#   Run from the simcrc_chile root, after cea_analysis_chile_parallel.R:
#     Rscript analysis/cea/population_colonoscopy_demand.R
#   Program colonoscopies get the CEA's 5% repeat rate (col_infl_rate = 1.05
#   in R/process_uspstf_output.R). Perfect adherence, as in the CEA. Writes
#   outputs/cea/2026Chile/Base_ChilePop/colonoscopies_by_age_chilepop.csv and
#   outputs/cea/population_demand/colonoscopy_demand_by_strategy.csv.
#
# *****************************************************************************

# 1. Packages ----------------------------------------------------------------

suppressMessages(library(data.table))

dir_pop <- "outputs/cea/2026Chile/Base_ChilePop"
dir_out <- "outputs/cea/population_demand"
dir.create(dir_out, showWarnings = FALSE)
col_infl_rate <- 1.05

# 2. Colonoscopies by strategy and age, FONASA population -----------------

dt_age <- rbindlist(lapply(
  list.files(file.path(dir_pop, "RawModelOutput_SimCRC_R"), full.names = TRUE),
  function(f) {
    dt <- fread(f, skip = "age,n_alive")
    dt[, .(strategy     = sub("^2026Chile_(.*)[.]csv$", "\\1", basename(f)),
           age,
           screening    = (n_pos_scrcol + n_neg_scrcol) * col_infl_rate,
           followup     = (n_pos_fucol + n_neg_fucol) * col_infl_rate,
           surveillance = (n_pos_survcol + n_neg_survcol) * col_infl_rate,
           symptomatic  = n_CRCstageI_sym + n_CRCstageII_sym +
             n_CRCstageIII_sym + n_CRCstageIV_sym)]
  }))
dt_age[, total := screening + followup + surveillance]

# 3. Totals by strategy -------------------------------------------------------

v_type <- c("screening", "followup", "surveillance", "total", "symptomatic")
dt_strategy <- dt_age[, lapply(.SD, sum), by = strategy, .SDcols = v_type]
dt_icer <- fread("outputs/cea/ce_results/df_icer_all_strategies.csv")
dt_icer[, strategy := gsub("_", "", Strategy)]
dt_icer[, strategy := sub("^(COL|FIT)(\\d{4})(\\d+)$", "\\1\\2q\\3", strategy)]
dt_strategy[dt_icer, efficient := i.Status == "ND", on = "strategy"]
setorder(dt_strategy, total)

# 4. Write outputs -------------------------------------------------------------

fwrite(dt_age, file.path(dir_pop, "colonoscopies_by_age_chilepop.csv"))
fwrite(dt_strategy, file.path(dir_out, "colonoscopy_demand_by_strategy.csv"))

cat("Strategies:", nrow(dt_strategy), "| ages", min(dt_age$age), "-",
    max(dt_age$age), "\n\nEfficient strategies, colonoscopies per year:\n")
print(dt_strategy[efficient == TRUE, lapply(.SD, round), by = strategy,
                  .SDcols = v_type])
