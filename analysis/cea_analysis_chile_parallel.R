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
library(ggpop)
library(ggrepel)
library(openxlsx)
library(tibble)
library(stringr)
library(parallel)
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
# CEA_CAL_RUN overrides the pinned run, so the same script can produce the
# no-recalibration control Fernando asked for (v0.14.0.2 parameters under the
# installed v0.15.0). Empty string restores the auto-select path below.
.cal_run_env <- Sys.getenv("CEA_CAL_RUN", "v0.15.0.20260921.1018")
calibration_run_target <- if (nzchar(.cal_run_env)) .cal_run_env else NULL

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
# Must match the calibration that produced these parameters: the September run
# fitted at 15, the August v0.14.0.2 run at 10. Mismatching them is not a
# control, it is a third model.
l_params_all$min_age_lesion_onset <-
  as.numeric(Sys.getenv("CEA_MIN_AGE_ONSET", "15"))
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
                                      SSP_pathway = FALSE)

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
.smoke_env <- Sys.getenv("CEA_SMOKE", "")
smoke_subset <- if (nzchar(.smoke_env)) strsplit(.smoke_env, ",")[[1]] else NULL

if (!is.null(smoke_subset)) {
  df_screening_strategies <- df_screening_strategies %>%
    filter(strategy %in% smoke_subset)
  stopifnot(nrow(df_screening_strategies) == length(smoke_subset))
}


# A control run writes everywhere the production run does -- raw CSVs, the
# xlsx, the ICER tables and the two plots -- so both are redirected together.
.scenario_tag <- Sys.getenv("CEA_SCENARIO", "")
if (nzchar(.scenario_tag)) df_screening_strategies$scenario <- .scenario_tag
ce_out <- Sys.getenv("CEA_CE_RESULTS", "ce_results")
dir.create(ce_out, recursive = TRUE, showWarnings = FALSE)
message("CEA scenario: ", unique(df_screening_strategies$scenario),
        " | ce_results -> ", ce_out,
        " | min_age_lesion_onset = ", l_params_all$min_age_lesion_onset)

n_ids <- nrow(df_screening_strategies)

output_template_year <- 2028

# *****************************************************************************
#### 3.0 Run screening and surveillance (parallel) ####
# *****************************************************************************

# CEA_CORES caps the worker count. 8 is fine for the September parameters, but
# min_age_lesion_onset = 10 carries far more lesions per person and 8 workers on
# the annual-FIT strategies exhausted 64 GB and got the run killed. Drop to 4
# for that configuration.
n_cores <- min(as.integer(Sys.getenv("CEA_CORES", "8")), parallel::detectCores())

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

# Backend is chosen by where the script runs, because the two constraints conflict.
# Fork shares the population copy-on-write and is much faster, but forking from
# RStudio is unreliable; and a live progress bar needs doSNOW's progress callback,
# which only exists on a PSOCK cluster. So: RStudio gets PSOCK plus the bar,
# a terminal gets fork plus the per-strategy lines in log_file.
force_backend <- Sys.getenv("CEA_BACKEND", "")   # "fork" or "psock" to override
in_rstudio    <- Sys.getenv("RSTUDIO") == "1" || .Platform$GUI == "RStudio"
use_fork      <- if (nzchar(force_backend)) {
  force_backend == "fork"
} else {
  !in_rstudio && .Platform$OS.type != "windows"
}

if (use_fork) {
  library(doParallel)
  cl <- makeForkCluster(n_cores)
  registerDoParallel(cl)
  log_msg("Backend: fork + doParallel (fastest). Progress -> tail -f the log file.")
} else {
  library(doSNOW)
  cl <- makeCluster(n_cores)
  registerDoSNOW(cl)
  log_msg("Backend: PSOCK + doSNOW (RStudio-safe, slower). Progress bar below.")
}

t_parallel_start <- proc.time()

fe_args <- list(
  i             = 1:n_ids,
  .packages     = c("data.table", "simcrc", "dplyr"),
  .export       = c("uspstf_summary", "crc_allocation_time")
)

if (!use_fork) {
  pb <- utils::txtProgressBar(min = 0, max = n_ids, style = 3)
  fe_args$.options.snow <- list(progress = function(n) {
    utils::setTxtProgressBar(pb, n)
    utils::flush.console()
  })
}

do.call(foreach::foreach, fe_args) %dopar% {

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

    set.seed(3)
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
if (!use_fork) close(pb)
stopCluster(cl)

t_parallel_total <- (proc.time() - t_parallel_start)[["elapsed"]]
log_msg(sprintf("All %d strategies completed in %.1f minutes (%.0f seconds) on %d cores.",
                n_ids, t_parallel_total / 60, t_parallel_total, n_cores))


# *****************************************************************************
###  4.0 Process the model output ---------------------------------------------
# *****************************************************************************
source("R/identify_strategies_to_remove.R")
source("R/crc_allocation_time.R")
source("R/uspstf_summary.R")
source("R/process_uspstf_output.R")
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

write.csv(icer_all_stategies, file = file.path(ce_out, "df_icer_all_strategies.csv"))

# ---- Prepare the plotting frame ---------------------------------------------
wtp_threshold <- 16e6   # willingness to pay per QALY, CLP

v_icons  <- c(Colonoscopy = "user-doctor", FIT = "vial",     `No screening` = "users-slash")
v_colors <- c(Colonoscopy = "#2a78d6",     FIT = "#eb6834",  `No screening` = "#1baf7a")
v_sel    <- "#a01b1b"   # the selected strategy, and the frontier line

# Label plaques take an 18% tint of their modality hue, so the readout is tied to
# its icon by colour. Kept this light on purpose: dark text over each tint stays
# above 12:1, where the hue at full strength would not.
v_fills  <- vapply(v_colors,
                   function(x) grDevices::colorRampPalette(c("#ffffff", x))(100)[18],
                   character(1))

df_ce <- data.frame(
  Strategy = as.character(icer_all_stategies$Strategy),
  cost     = as.numeric(icer_all_stategies$Cost) / 1e6,
  effect   = as.numeric(icer_all_stategies$Effect),
  status   = as.character(icer_all_stategies$Status),
  icer     = as.numeric(icer_all_stategies$ICER),
  stringsAsFactors = FALSE)

df_ce$modality <- factor(
  ifelse(grepl("NoScreening", df_ce$Strategy), "No screening",
         ifelse(grepl("^FIT", df_ce$Strategy), "FIT", "Colonoscopy")),
  levels = names(v_icons))
df_ce$icon_name <- unname(v_icons[as.character(df_ce$modality)])
df_ce$efficient <- df_ce$status == "ND"
df_ce$alpha_eff <- ifelse(df_ce$efficient, 1, 0.40)

# COL_45_70_15 -> "COL 45-70 q15"
df_ce$label <- vapply(df_ce$Strategy, function(s) {
  if (grepl("NoScreening", s)) return("No screening")
  p <- strsplit(s, "_")[[1]]
  if (length(p) == 4L) sprintf("%s %s-%s q%s", p[1], p[2], p[3], p[4]) else s
}, character(1), USE.NAMES = FALSE)

df_frontier <- df_ce[df_ce$efficient, ]
df_frontier <- df_frontier[order(df_frontier$effect), ]

# Optimal at WTP = last frontier point still reached by a step at or below wtp.
# Deliberately not a ray from no screening: that reads as AVERAGE
# cost-effectiveness, and every Chile strategy clears 16M on that basis while
# incremental ICERs run past 400M.
v_icer <- df_frontier$icer
v_icer[is.na(v_icer)] <- 0
i_opt <- max(which(cumsum(v_icer > wtp_threshold) == 0))
df_frontier$is_opt <- seq_len(nrow(df_frontier)) == i_opt
df_frontier$label[i_opt] <- sprintf(
  "%s\n%sM  -  %.1f QALYs\nICER %.1fM/QALY",
  df_frontier$label[i_opt],
  formatC(df_frontier$cost[i_opt], format = "f", digits = 0, big.mark = ","),
  df_frontier$effect[i_opt], v_icer[i_opt] / 1e6)

library(ggpop)
library(ggrepel)
# ---- Plot the efficient frontier --------------------------------------------
# Icon marks carry modality by shape as well as hue, so the figure survives
# greyscale and colour-vision deficiency. icon is mapped in aes() rather than
# fixed per layer because ggplot draws every show.legend layer's glyph into
# every key, which would overlay a stethoscope and a vial in both.
plot_ce <-
  ggplot(df_ce, aes(x = effect, y = cost)) +
  geom_line(data = df_frontier, colour = "#a01b1b", linewidth = 0.7) +
  # The layer that owns the legend must contain EVERY modality. With a level
  # present only in this layer and not in the other icon layer, ggpop draws the
  # previous key's glyph for it -- No screening came out as a green vial. So one
  # layer holds all the points and owns the legend, and the frontier is redrawn
  # bolder on top. alpha must be numeric here; ggpop rejects a logical column.
  geom_icon_point(aes(colour = modality, icon = icon_name, alpha = alpha_eff),
                  size = 1.3, show.legend = TRUE, legend_icons = TRUE,dpi = 500) +
  # the selected strategy is dropped here and redrawn in red below, so its
  # modality-coloured icon does not show through from underneath
  geom_icon_point(data = df_frontier[-i_opt, ],
                  aes(colour = modality, icon = icon_name),
                  size = 2.3, show.legend = FALSE,dpi = 500) +
  # the selected strategy is redrawn in the frontier red so it is found at a
  # glance; colour is fixed rather than mapped, so it cannot disturb the legend
  geom_icon_point(data = df_frontier[i_opt, ],
                  aes(icon = icon_name), colour = v_sel,
                  size = 2.3, show.legend = FALSE,dpi = 500) +
  scale_alpha_identity(guide = "none") +
  scale_fill_manual(values = v_fills, guide = "none") +
  geom_label(
    data = df_frontier,
    aes(label = label, fill = modality,
        fontface = ifelse(is_opt, "bold", "plain")),
    hjust = 0,
    # the selected readout is pushed further down-right into open space so it
    # clears the labels stacked along the steep part of the frontier
    nudge_x = ifelse(df_frontier$is_opt, 0.055, 0.014) *
              diff(range(df_ce$effect)),
    nudge_y = ifelse(df_frontier$is_opt, -0.075, 0) *
              diff(range(df_ce$cost)),
    size = 2.4, colour = "#26261f", lineheight = 1.05,
    label.size = 0.18, label.r = unit(0.1, "lines"),
    label.padding = unit(0.16, "lines")) +
  scale_colour_manual(values = v_colors, name = NULL) +
  scale_legend_icon(size = 7) +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 10),
                     minor_breaks = NULL,
                     expand = expansion(mult = c(0.02, 0.20))) +
  scale_y_continuous(labels = scales::label_number(accuracy = 1, big.mark = ",",
                                                   suffix = "M"),
                     breaks = scales::breaks_pretty(n = 10),
                     minor_breaks = NULL) +
  labs(title    = "Cost-effectiveness of CRC screening strategies, Chile",
       subtitle = paste0(simcrc_model_version, " - ", n_ids, " strategies, ",
                         format(n_pop, big.mark = ",", scientific = FALSE),
                         " cohort, 3% discounting"),
       x        = "Discounted QALYs gained per 1,000",
       y        = "Discounted total costs per 1,000 (CLP millions)",
       caption  = paste0(
         "Icons on the frontier are solid; dominated strategies are faded.",
         "\nBold label = optimal at a willingness to pay of 16.0M per QALY (",
         df_ce$label[match(df_frontier$Strategy[i_opt], df_ce$Strategy)],
         "). Chosen on the incremental ICER, not the ratio to no screening.")) +
  theme_pop(base_size = 10) +
  theme(axis.title       = element_text(colour = "#57574f", size = 9.5),
        axis.text        = element_text(colour = "#77776e", size = 8.5),
        axis.ticks       = element_line(colour = "#c9c9c2", linewidth = 0.3),
        axis.line        = element_line(colour = "#c9c9c2", linewidth = 0.3),
        panel.grid.major = element_line(colour = "#ececE8", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        plot.caption     = element_text(colour = "#77776e", hjust = 0, size = 7),
        plot.caption.position = "plot")

plot_ce

ggsave(filename = file.path(ce_out, "plot_ce_all_strategies.png"), plot = plot_ce,
       width = 8.5, height = 5.4, units = "in", dpi = 300, bg = "white")


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

write.csv(icer_FIT, file = file.path(ce_out, "df_icer_FIT.csv"))

# ---- Prepare the plotting frame ---------------------------------------------
# (wtp_threshold, v_icons and v_colors are set in section 5.0 above)


df_ce_FIT <- data.frame(
  Strategy = as.character(icer_FIT$Strategy),
  cost     = as.numeric(icer_FIT$Cost) / 1e6,
  effect   = as.numeric(icer_FIT$Effect),
  status   = as.character(icer_FIT$Status),
  icer     = as.numeric(icer_FIT$ICER),
  stringsAsFactors = FALSE)

df_ce_FIT$modality <- factor(
  ifelse(grepl("NoScreening", df_ce_FIT$Strategy), "No screening",
         ifelse(grepl("^FIT", df_ce_FIT$Strategy), "FIT", "Colonoscopy")),
  levels = names(v_icons))
df_ce_FIT$icon_name <- unname(v_icons[as.character(df_ce_FIT$modality)])
df_ce_FIT$efficient <- df_ce_FIT$status == "ND"
df_ce_FIT$alpha_eff <- ifelse(df_ce_FIT$efficient, 1, 0.40)

# COL_45_70_15 -> "COL 45-70 q15"
df_ce_FIT$label <- vapply(df_ce_FIT$Strategy, function(s) {
  if (grepl("NoScreening", s)) return("No screening")
  p <- strsplit(s, "_")[[1]]
  if (length(p) == 4L) sprintf("%s %s-%s q%s", p[1], p[2], p[3], p[4]) else s
}, character(1), USE.NAMES = FALSE)

df_frontier_FIT <- df_ce_FIT[df_ce_FIT$efficient, ]
df_frontier_FIT <- df_frontier_FIT[order(df_frontier_FIT$effect), ]

# Optimal at WTP = last frontier point still reached by a step at or below wtp.
# Deliberately not a ray from no screening: that reads as AVERAGE
# cost-effectiveness, and every Chile strategy clears 16M on that basis while
# incremental ICERs run past 400M.
v_icer_FIT <- df_frontier_FIT$icer
v_icer_FIT[is.na(v_icer_FIT)] <- 0
i_opt_FIT <- max(which(cumsum(v_icer_FIT > wtp_threshold) == 0))
df_frontier_FIT$is_opt <- seq_len(nrow(df_frontier_FIT)) == i_opt_FIT
df_frontier_FIT$label[i_opt_FIT] <- sprintf(
  "%s\n%sM  -  %.1f QALYs\nICER %.1fM/QALY",
  df_frontier_FIT$label[i_opt_FIT],
  formatC(df_frontier_FIT$cost[i_opt_FIT], format = "f", digits = 0, big.mark = ","),
  df_frontier_FIT$effect[i_opt_FIT], v_icer_FIT[i_opt_FIT] / 1e6)

# ---- Plot the efficient frontier --------------------------------------------
# Icon marks carry modality by shape as well as hue, so the figure survives
# greyscale and colour-vision deficiency. icon is mapped in aes() rather than
# fixed per layer because ggplot draws every show.legend layer's glyph into
# every key, which would overlay a stethoscope and a vial in both.
plot_ce_FIT <-
  ggplot(df_ce_FIT, aes(x = effect, y = cost)) +
  geom_line(data = df_frontier_FIT, colour = "#a01b1b", linewidth = 0.7) +
  # The layer that owns the legend must contain EVERY modality. With a level
  # present only in this layer and not in the other icon layer, ggpop draws the
  # previous key's glyph for it -- No screening came out as a green vial. So one
  # layer holds all the points and owns the legend, and the frontier is redrawn
  # bolder on top. alpha must be numeric here; ggpop rejects a logical column.
  geom_icon_point(aes(colour = modality, icon = icon_name, alpha = alpha_eff),
                  size = 1.3, show.legend = TRUE, legend_icons = TRUE, dpi = 500) +
  # the selected strategy is dropped here and redrawn in red below, so its
  # modality-coloured icon does not show through from underneath
  geom_icon_point(data = df_frontier_FIT[-i_opt_FIT, ],
                  aes(colour = modality, icon = icon_name),
                  size = 2.3, show.legend = FALSE, dpi = 500) +
  # the selected strategy is redrawn in the frontier red so it is found at a
  # glance; colour is fixed rather than mapped, so it cannot disturb the legend
  geom_icon_point(data = df_frontier_FIT[i_opt_FIT, ],
                  aes(icon = icon_name), colour = v_sel,
                  size = 2.3, show.legend = FALSE, dpi = 500) +
  scale_alpha_identity(guide = "none") +
  scale_fill_manual(values = v_fills, guide = "none") +
  geom_label(
    data = df_frontier_FIT,
    aes(label = label, fill = modality,
        fontface = ifelse(is_opt, "bold", "plain")),
    hjust = 0, nudge_x = 0.014 * diff(range(df_ce_FIT$effect)),
    size = 2.4, colour = "#26261f", lineheight = 1.05,
    label.size = 0.18, label.r = unit(0.1, "lines"),
    label.padding = unit(0.16, "lines")) +
  scale_colour_manual(values = v_colors, name = NULL) +
  scale_legend_icon(size = 7) +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 10),
                     minor_breaks = NULL,
                     expand = expansion(mult = c(0.02, 0.20))) +
  scale_y_continuous(labels = scales::label_number(accuracy = 1, big.mark = ",",
                                                   suffix = "M"),
                     breaks = scales::breaks_pretty(n = 10),
                     minor_breaks = NULL) +
  labs(title    = "Cost-effectiveness of FIT strategies, Chile",
       subtitle = paste0(simcrc_model_version, " - no screening + ",
                         nrow(df_uspstf_output_FIT) - 1, " FIT strategies"),
                         
       x        = "Discounted QALYs gained per 1,000",
       y        = "Discounted total costs per 1,000 (CLP millions)",
       caption  = paste0(
         "Icons on the frontier are solid; dominated strategies are faded.",
         "\nBold label = optimal at a willingness to pay of 16.0M per QALY (",
         df_ce_FIT$label[match(df_frontier_FIT$Strategy[i_opt_FIT], df_ce_FIT$Strategy)],
         "). Chosen on the incremental ICER, not the ratio to no screening.")) +
  theme_pop(base_size = 10) +
  theme(axis.title       = element_text(colour = "#57574f", size = 9.5),
        axis.text        = element_text(colour = "#77776e", size = 8.5),
        axis.ticks       = element_line(colour = "#c9c9c2", linewidth = 0.3),
        axis.line        = element_line(colour = "#c9c9c2", linewidth = 0.3),
        panel.grid.major = element_line(colour = "#ececE8", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        plot.caption     = element_text(colour = "#77776e", hjust = 0, size = 7),
        plot.caption.position = "plot")


plot_ce_FIT

ggsave(filename = file.path(ce_out, "plot_ce_FIT.png"), plot = plot_ce_FIT,
       width = 8.5, height = 5.4, units = "in", dpi = 300, bg = "white")
}

