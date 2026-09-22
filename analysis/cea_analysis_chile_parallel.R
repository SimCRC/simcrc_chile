# *****************************************************************************
#
# Script: cea_analysis_chile_parallel.R
#
# Purpose: Run Cost-Effectiveness Analysis (CEA) for CRC screening strategies
#          (parallelized version)
#
# Author: Jorge Roa, Claudia Seguin and Carlos Pineda
#
# Email: jorgeroa@stanford.edu/cseguin@stanford.edu/cpinedaa@uw.edu
#
# Date Created: 01 August 2025
#
# *****************************************************************************


remove(list = ls())
gc()

# *****************************************************************************
#### Load packages and functions ####
# *****************************************************************************

library(data.table)
library(simcrc)
library(dplyr)
library(remotes)
library(CRCmortality)
library(dplyr)
library(ggplot2)
library(readxl)
library(dampack)
library(openxlsx)
library(tibble)
library(stringr)
library(parallel)
library(doParallel)
library(foreach)


# *****************************************************************************
#### 1.0 Natural history ####
# *****************************************************************************

#install prerelease simcrc

# remotes::install_github("simcrc/simcrc", ref = "v0.13.002", force = TRUE)
# ^ disabled: tag "v0.13.002" does not exist (real tags are v0.13.0.02/.03/.04), so this
#   line always errored; had it resolved it would downgrade simcrc on every CEA run.

# identify_strategies_to_remove() lives in this repo, not in simcrc; source it so the
# script runs in a fresh session instead of relying on an interactive global env.
source("R/identify_strategies_to_remove.R")
source("R/crc_allocation_time.R")
source("R/uspstf_summary.R")
source("R/process_uspstf_output.R")





# Get the model version
simcrc_model_version <- paste0("SimCRC v",as.character(packageVersion("simcrc")))

# Load calibrated parameters
# ---- Calibration source ------------------------------------------------------
# Set calibration_run_target to the exact run folder this CEA belongs to. Pinning
# it is deliberate: with auto-select (below) a CEA started before the intended
# calibration has written its l_params_calibrated_sets_*.RData silently falls back
# to the previous run's posteriors -- no error, wrong results. Pinned, that case
# stops the script instead.
# Set to NULL to fall back to "most recent run that has a best-parameter-set file"
# (folders are named vX.Y.Z.YYYYMMDD.HHMM, so a descending name sort is newest).
calibration_run_target <- "v0.15.0.20260921.1018"

calibration_base <- "outputs/BayCANN_versions/Chile/Adenoma/F"
cal_sets_files <- list.files(calibration_base,
                             pattern = "^l_params_calibrated_sets_.*\\.RData$",
                             recursive = TRUE, full.names = TRUE)
if (length(cal_sets_files) == 0) {
  stop("No l_params_calibrated_sets_*.RData found under '", calibration_base,
       "'. Run 12_best_param_set.R for a calibration before the CEA.")
}
cal_runs <- basename(dirname(cal_sets_files))

if (is.null(calibration_run_target)) {
  sel_idx <- order(cal_runs, decreasing = TRUE)[1]
  message("CEA calibration source: auto-selected most recent run.")
} else {
  sel_idx <- match(calibration_run_target, cal_runs)
  if (is.na(sel_idx)) {
    stop("Pinned calibration run '", calibration_run_target,
         "' has no l_params_calibrated_sets_*.RData yet under '", calibration_base,
         "'.\n  The calibration has not finished 12_best_param_set.R -- wait for it, or",
         " set calibration_run_target <- NULL to use the most recent available run.",
         "\n  Runs currently available: ", paste(sort(cal_runs, decreasing = TRUE),
                                                 collapse = ", "))
  }
}

calibration_run    <- cal_runs[sel_idx]
calibration_folder <- dirname(cal_sets_files[sel_idx])
message("CEA using calibration run: ", calibration_run)
load(cal_sets_files[sel_idx])
l_params_Min_AbsolutErr <- l_params_calibrated_sets$Min_AbsolutErr
l_params_all <- load_params_init(fromFile = TRUE, filename = l_params_Min_AbsolutErr)

# Update the model start age and the survival by race defaults. (CHANGE THE DEFAULTS IN SIMCRC)
l_params_all$min_age_lesion_onset <- 15
l_params_all$mort_by_race <- FALSE
# l_params_all$year_surv_improv <- 2003    # We haven't adjusted this for Chile

# Define the simulation population size
n_pop <- 1e6 # Run at least 1 mil for publications, 10mil if possible for stable results

# Define the cohort age
cohort_age <- 40

# Sample the age of death from life table
df_lt_chile <- read.csv("data-raw/df_lifetable_2017_CH.csv")
colnames(df_lt_chile)[colnames(df_lt_chile) == "age"] <- "Age"
colnames(df_lt_chile)[colnames(df_lt_chile) == "mortality_rate"] <- "mortality.rates"



dt_pop <- simcrc::get_dt_population(year = 1980,
                                    byear = 1980,
                                    p_female = 1,
                                    p_white = 0.8,
                                    n_pop = n_pop,
                                    dt_life_table_F = df_lt_chile,
                                    dt_life_table_M = NULL)

# Run SimCRC natural history
l_out_simcrc <- simcr_nathist_ssp_DES(l_params_all = l_params_all,
                                      dt_pop = dt_pop,
                                      SSP_pathway = FALSE,
                                      optimize_memory = TRUE)

dt_crc_pop <- l_out_simcrc$dt_crc_pop


# *****************************************************************************
#### 2.0 Screening strategies ####
# *****************************************************************************

# Load the screening strategies for this project
df_screening_strategies <- readr::read_csv("data-raw/df_CH_2026_strategies.csv", show_col_types = FALSE)

# It's okay to not run this, only prevents duplicated runs
# Add function to remove additional screening strategies which are effectively the same
# (e.g. COL4575q10 & COL4580q10, screening happens at 45,55,65,75, so we will remove one of them)
l_strategies_to_remove <- identify_strategies_to_remove(modality = "COL",
                                                        start_ages = c(45,50,55),
                                                        stop_ages  = c(70,75,80,85),
                                                        intervals  = c(5,10,15))
df_screening_strategies <- df_screening_strategies %>% filter(!strategy %in% l_strategies_to_remove)

# Check that there are no duplicate ids to prevent overwriting
if(any(duplicated(df_screening_strategies$id))){
  stop("There are duplicate ids in the screening strategies. Please check the input file and remove duplicates.")
}

# Pipeline smoke subset. Set to NULL to run the full Chile set (68 strategies:
# NoScreening + 31 COL + 36 FIT). Named rather than positional so the subset always
# spans both modalities -- section 5.0 needs >= 2 rows and section 6.0 needs >= 2
# NoScreening/FIT rows, so a COL-only subset silently skips the FIT frontier.
# ALWAYS clear RawModelOutput_SimCRC_R before a subset run: ProcessUSPSTFOutput globs
# that folder, so leftover CSVs from a previous run get blended in silently.
smoke_subset <- NULL

if (!is.null(smoke_subset)) {
  df_screening_strategies <- df_screening_strategies %>%
    filter(strategy %in% smoke_subset)
  stopifnot(nrow(df_screening_strategies) == length(smoke_subset))
}


n_ids <- nrow(df_screening_strategies)

output_template_year <- 2028

# *****************************************************************************
#### 3.0 Run screening and surveillance (parallel) ####
# *****************************************************************************

n_cores <- min(8L, parallel::detectCores())

cat(sprintf("Running %d strategies across %d cores\n", n_ids, n_cores))

# Use absolute path so parallel workers can find the file regardless of their working directory
log_file <- normalizePath(
  paste0("output/parallel_run_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"),
  mustWork = FALSE
)
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
cat(sprintf("Progress log -> %s\n", log_file))

log_msg <- function(msg) {
  line <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", msg, "\n")
  cat(line)
  cat(line, file = log_file, append = TRUE)
}

log_msg(sprintf("Starting: %d strategies on %d cores", n_ids, n_cores))

cl <- makeForkCluster(n_cores)   # default: worker startup banners stay hidden (no outfile="")
registerDoParallel(cl)

t_parallel_start <- proc.time()
set.seed(3)

foreach(
  i = 1:n_ids,
  .packages     = c("data.table", "simcrc", "dplyr"),
  .export       = c("uspstf_summary", "crc_allocation_time")
) %dopar% {

  data.table::setDTthreads(1L)

  run_id <- df_screening_strategies[i, ]

  # Run info
  id              <- run_id$id
  strategy_name   <- run_id$strategy
  screening_years <- seq(run_id$age_to_begin_screening,
                         run_id$age_to_end_screening,
                         run_id$frequency_screening)
  screening_modality <- run_id$modality
  follow_up          <- run_id$follow_up

  # Record start to the log file only (disk, for an optional `tail -f`); console output
  # is printed master-side by the completion callback to avoid worker interleaving.
  cat(sprintf("[%s] START strategy %d/%d: %s\n",
              format(Sys.time(), "%H:%M:%S"), i, n_ids, strategy_name),
      file = log_file, append = TRUE)

  start_time <- Sys.time()

  capture.output({

    # Compact path (simcrc >= 0.15.0): screening_surveillance_counts() counts
    # every round before releasing its per-round columns, freeing 23 columns in
    # screening and 24 in surveillance, and returns the USPSTF summary directly.
    l_screening_args <- list(
      p_coverage               = run_id$p_coverage,
      p_adherence              = run_id$p_adherence,
      age_to_begin_screening   = run_id$age_to_begin_screening,
      age_to_end_screening     = run_id$age_to_end_screening,
      frequency_screening      = run_id$frequency_screening,
      sens_small_adenoma       = run_id$sens_small_adenoma,
      sens_medium_adenoma      = run_id$sens_medium_adenoma,
      sens_large_adenoma       = run_id$sens_large_adenoma,
      sens_crc                 = run_id$sens_crc,
      sens_by_ad               = run_id$sens_by_ad,
      spec                     = run_id$spec,
      screening_reach          = run_id$screening_reach,
      p_reach_cecum            = run_id$p_reach_cecum,
      p_reach_ascending        = run_id$p_reach_ascending,
      p_reach_transverse       = run_id$p_reach_transverse,
      p_reach_descending       = run_id$p_reach_descending,
      p_reach_sigmoid          = run_id$p_reach_sigmoid,
      p_reach_rectum           = run_id$p_reach_rectum,
      p_death_scr              = run_id$p_death_scr,
      surveillance             = run_id$surveillance,
      confirmation             = run_id$follow_up,
      p_adherence_confirmation = run_id$p_adherence_confirmation,
      sens_small_adenoma_col   = run_id$sens_small_adenoma_col,
      sens_medium_adenoma_col  = run_id$sens_medium_adenoma_col,
      sens_large_adenoma_col   = run_id$sens_large_adenoma_col,
      sens_crc_col             = run_id$sens_crc_col,
      sens_by_ad_col           = run_id$sens_by_ad_col,
      spec_col                 = run_id$spec_col,
      confirmation_reach       = run_id$confirmation_reach,
      p_reach_cecum_conf       = run_id$p_reach_cecum_conf,
      p_reach_ascending_conf   = run_id$p_reach_ascending_conf,
      p_reach_transverse_conf  = run_id$p_reach_transverse_conf,
      p_reach_descending_conf  = run_id$p_reach_descending_conf,
      p_reach_sigmoid_conf     = run_id$p_reach_sigmoid_conf,
      p_reach_rectum_conf      = run_id$p_reach_rectum_conf,
      p_death_conf             = run_id$p_death_conf
    )

    l_surveillance_args <- list(
      p_adherence             = run_id$p_adherence_surv,
      sens_small_adenoma      = run_id$sens_small_adenoma_surv,
      sens_medium_adenoma     = run_id$sens_medium_adenoma_surv,
      sens_large_adenoma      = run_id$sens_large_adenoma_surv,
      sens_crc                = run_id$sens_crc_surv,
      sens_by_ad              = run_id$sens_by_ad_surv,
      spec                    = run_id$spec_surv,
      surveillance_reach      = run_id$surveillance_reach,
      p_reach_cecum_surv      = run_id$p_reach_cecum_surv,
      p_reach_ascending_surv  = run_id$p_reach_ascending_surv,
      p_reach_transverse_surv = run_id$p_reach_transverse_surv,
      p_reach_descending_surv = run_id$p_reach_descending_surv,
      p_reach_sigmoid_surv    = run_id$p_reach_sigmoid_surv,
      p_reach_rectum_surv     = run_id$p_reach_rectum_surv,
      p_death_surv            = run_id$p_death_surv
    )

    l_out <- screening_surveillance_counts(
      dt_crc_pop           = dt_crc_pop,
      l_screening_args     = l_screening_args,
      l_surveillance_args  = l_surveillance_args,
      copy_input           = TRUE,
      screening_modality   = screening_modality,
      output_template_year = output_template_year,
      min_age              = cohort_age,
      max_age              = 100
    )

    dt_export_final <- l_out$summary

  })

  end_time <- Sys.time()
  execution_time <- difftime(end_time, start_time, units = "mins")

  # Create a header string
  time <- format(Sys.time(), "%Y-%m-%d_%H:%M")
  # We never updated the output template string for 2020 since we don't need it moving forward
  # But the fields themselves are correct.
  if(output_template_year == 2020){
    if(screening_modality == "COL" | screening_modality == "SIG" | screening_modality == "CTC"){
      header_string <- paste0("date,", time,"\n",
                              "model,", simcrc_model_version,"\n",
                              "risk_scenario,", "BASE","\n",
                              "population,", "TOTAL","\n",
                              "adherence,", "PERFECT","\n",
                              "stooltest_type,", NA,"\n",
                              "stooltest_startage,", NA,"\n",
                              "stooltest_stopage,", NA,"\n",
                              "stooltest_interval,", NA,"\n",
                              "structuralexam_type,", run_id$modality,"\n",
                              "structuralexam_startage,", run_id$age_to_begin_screening,"\n",
                              "structuralexam_stopage,", run_id$age_to_end_screening,"\n",
                              "structuralexam_interval,", run_id$frequency_screening,"\n")
    }else if(screening_modality == "FIT" | screening_modality == "sDNA-FIT"){
      header_string <- paste0("date,", time,"\n",
                              "model,", simcrc_model_version,"\n",
                              "risk_scenario,", "BASE","\n",
                              "population,", "TOTAL","\n",
                              "adherence,", "PERFECT","\n",
                              "stooltest_type,", run_id$modality,"\n",
                              "stooltest_startage,", run_id$age_to_begin_screening,"\n",
                              "stooltest_stopage,", run_id$age_to_end_screening,"\n",
                              "stooltest_interval,", run_id$frequency_screening,"\n",
                              "structuralexam_type,", NA,"\n",
                              "structuralexam_startage,", NA,"\n",
                              "structuralexam_stopage,", NA,"\n",
                              "structuralexam_interval,", NA,"\n")
    }else{
      header_string <- paste0("date,", time,"\n",
                              "model,", simcrc_model_version,"\n",
                              "risk_scenario,", "BASE","\n",
                              "population,", "TOTAL","\n",
                              "adherence,", "PERFECT","\n",
                              "stooltest_type,", NA,"\n",
                              "stooltest_startage,", NA,"\n",
                              "stooltest_stopage,", NA,"\n",
                              "stooltest_interval,", NA,"\n",
                              "structuralexam_type,", NA,"\n",
                              "structuralexam_startage,", NA,"\n",
                              "structuralexam_stopage,", NA,"\n",
                              "structuralexam_interval,", NA,"\n")
    }
  }else if(output_template_year == 2028){
    if(screening_modality == "COL" | screening_modality == "SIG" | screening_modality == "CTC"){
      header_string <- paste0("date,", time,"\n",
                              "model,", simcrc_model_version,"\n",
                              "risk_scenario,", run_id$scenario,"\n",
                              "population,", "TOTAL","\n",
                              "adherence_fucol,", run_id$p_adherence_confirmation, "\n",
                              "screening_startage,", run_id$age_to_begin_screening,"\n",
                              "screening_stopage,", run_id$age_to_end_screening,"\n",
                              "stoolorbloodtest_type,", NA,"\n",
                              "stoolorbloodtest_interval,", NA,"\n",
                              "structuralexam_type,", run_id$modality,"\n",
                              "structuralexam_interval,", run_id$frequency_screening,"\n")
    }else if(screening_modality == "FIT" | screening_modality == "sDNA-FIT"| screening_modality == "mtsDNA"|
             screening_modality == "SRNA"| screening_modality == "BLOOD"){
      header_string <- paste0("date,", time,"\n",
                              "model,", simcrc_model_version,"\n",
                              "risk_scenario,", run_id$scenario,"\n",
                              "population,", "TOTAL","\n",
                              "adherence_fucol,", run_id$p_adherence_confirmation,"\n",
                              "screening_startage,", run_id$age_to_begin_screening,"\n",
                              "screening_stopage,", run_id$age_to_end_screening,"\n",
                              "stoolorbloodtest_type,", run_id$modality,"\n",
                              "stoolorbloodtest_interval,", run_id$frequency_screening,"\n",
                              "structuralexam_type,", NA,"\n",
                              "structuralexam_interval,", NA,"\n")
    }else{
      header_string <- paste0("date,", time,"\n",
                              "model,", simcrc_model_version,"\n",
                              "risk_scenario,", run_id$scenario,"\n",
                              "population,", "TOTAL","\n",
                              "adherence_fucol,", run_id$p_adherence_confirmation,"\n",
                              "screening_startage,", NA,"\n",
                              "screening_stopage,", NA,"\n",
                              "stoolorbloodtest_type,", NA,"\n",
                              "stoolorbloodtest_interval,", NA,"\n",
                              "structuralexam_type,", NA,"\n",
                              "structuralexam_interval,", NA,"\n")
    }
  }else{
    stop("Invalid output_template_year. Please choose either 2020 or 2028.")
  }

  # File name (absolute path so workers resolve it correctly)
  uspstf_folder <- normalizePath(
    paste0("output/", run_id$project, "/", run_id$scenario, "/RawModelOutput_SimCRC_R"),
    mustWork = FALSE
  )
  if (!dir.exists(uspstf_folder)) {
    dir.create(uspstf_folder, recursive = TRUE)
  }
  uspstf_filename <- file.path(uspstf_folder, paste0(run_id$project, "_", strategy_name, ".csv"))

  # Write the header to the file
  writeLines(header_string, uspstf_filename)

  # Append the data to the same file
  write.table(dt_export_final,
              file = uspstf_filename,
              col.names = TRUE, row.names = FALSE, quote = FALSE, append = TRUE, sep = ",")

  # Record completion to the log file only (disk); the console completion line is printed
  # master-side by the progress callback (in finish order, no interleaving).
  cat(sprintf("[%s] DONE  strategy %d/%d: %s (%.2f min)\n",
              format(Sys.time(), "%H:%M:%S"), i, n_ids, strategy_name, execution_time),
      file = log_file, append = TRUE)

  NULL
}

cat("\n")  # spacer after the per-strategy completion lines
stopCluster(cl)

t_parallel_total <- (proc.time() - t_parallel_start)[["elapsed"]]
log_msg(sprintf("All %d strategies completed in %.1f minutes (%.0f seconds) on %d cores.",
                n_ids, t_parallel_total / 60, t_parallel_total, n_cores))


# *****************************************************************************
###  4.0 Process the model output ---------------------------------------------
# *****************************************************************************

project <- unique(df_screening_strategies$project)
scenarios <- unique(df_screening_strategies$scenario)

df_uspstf_output <- ProcessUSPSTFOutput(analysis_folder = paste0("output/", project, "/", scenarios[1]),
                                        input_folder = "data-raw/cea_inputs",
                                        input_prefix = project,
                                        first_age_of_interest = cohort_age,
                                        col_infl_rate = 1.05,
                                        discount_rate = 0.03,
                                        col_spec_adj = 0.86,
                                        crc_care_costs_file = "crc_care_costs.csv",
                                        screen_costs_file = "screen_costs_v3.csv",
                                        crc_care_disutility_file = "crc_care_utility_loss.csv",
                                        screen_disutility_file = "screen_utility_loss_WithStoolTestValues.csv",
                                        general_health_utility_weights_file = "GeneralHealthUtilityWeightsByAge.csv",
                                        selected_outcomes_for_model_data_file = "model_data_outcomes_boolean_addDscQALY.csv",
                                        model_run_data_tag = "_Base",
                                        folder_for_output = "ce_results",
                                        include_SimCRC = F, ## Are you including SimCRC model results?
                                        include_SimCRC_R = T, ## Are you including SimCRC-R model results?
                                        include_MISCAN = F, ## Are you including MISCAN model results?
                                        include_CRCSPIN = F,
                                        output_template_year = output_template_year) ## either 2020 or 2028 depending on the output you are generating) ## Are you including CRC-SPIN model results?)


# *****************************************************************************
###  5.0 Perform the CEA for all strategies -----------------------------------
# *****************************************************************************

# Update the strategy label for plotting
df_uspstf_output$Strategy <- gsub("2026Chile_", "", df_uspstf_output$Strategy)

# Pick your cost variable (discounted costs) from the uspstf_output
v_crc_costs <- df_uspstf_output$DiscountedTotalCostsper1000

# Pick your benefit variable (discounted QALYG) from the uspstf_output
v_crc_qalys <- df_uspstf_output$DiscountedQALYGainedper1000

# Calculate icers using the vector of costs and vector of qalys
icer_all_stategies <- calculate_icers(cost = v_crc_costs,
                                      effect = v_crc_qalys,
                                      strategies = df_uspstf_output$Strategy)

write.csv(icer_all_stategies, file = "ce_results/df_icer_all_strategies.csv")

# Plot the efficient frontier
plot_ce <- dampack:::plot.icers(icer_all_stategies,
                                label = c("frontier"))

#replace y axis label to Discounted TotalCosts per 1000

plot_ce <- plot_ce + ylab("Discounted Total Costs per 1000") + xlab("Discounted QALYs Gained per 1000")

plot_ce

ggsave(filename = "ce_results/plot_ce_all_strategies.png", width = 6.5, height = 4, units = "in", dpi = 300)


# *****************************************************************************
###  6.0 Perform the CEA for only No Screening and FIT strategies -------------
# *****************************************************************************

# Add a column for the modality
df_uspstf_output <- df_uspstf_output %>% mutate(Modality = case_when(
  str_detect(Strategy, "NoScreening") ~ "NoScreening",
  str_detect(Strategy, "FIT") ~ "FIT",
  str_detect(Strategy, "COL") ~ "COL"))

# Filter only "NoScreening" and "FIT"
df_uspstf_output_FIT <- df_uspstf_output %>% filter(Modality == "NoScreening" | Modality == "FIT")

# Needs >= 2 strategies to compute an ICER frontier; skip cleanly on subsets without FIT.
if (nrow(df_uspstf_output_FIT) < 2) {
  message("Section 6.0 (NoScreening + FIT CEA) skipped: needs >= 2 matching strategies, found ",
          nrow(df_uspstf_output_FIT), " in the current subset.")
} else {

# Pick your cost variable (discounted costs) from the uspstf_output
v_crc_costs <- df_uspstf_output_FIT$DiscountedTotalCostsper1000

# Pick your benefit variable (discounted QALYG) from the uspstf_output
v_crc_qalys <- df_uspstf_output_FIT$DiscountedQALYGainedper1000

# Calculate icers using the vector of costs and vector of qalys
icer_FIT <- calculate_icers(cost = v_crc_costs,
                            effect = v_crc_qalys,
                            strategies = df_uspstf_output_FIT$Strategy)

write.csv(icer_FIT, file = "ce_results/df_icer_FIT.csv")

# Plot the efficient frontier
plot_ce_FIT <- dampack:::plot.icers(icer_FIT,
                                    label = c("frontier"))

plot_ce_FIT <- plot_ce_FIT + ylab("Discounted Total Costs per 1000") + xlab("Discounted QALYs Gained per 1000")




ggsave(filename = "ce_results/plot_ce_FIT.png", width = 6.5, height = 4, units = "in", dpi = 300)
}
