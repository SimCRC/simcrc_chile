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
library(ceacrc)
library(parallel)
library(doSNOW)
library(foreach)


# *****************************************************************************
#### 1.0 Natural history ####
# *****************************************************************************

# Get the model version
simcrc_model_version <- paste0("SimCRC v",as.character(packageVersion("simcrc")))

# Load calibrated parameters
load("outputs/BayCANN_versions/Chile/Adenoma/F/v0.13.0/v0.13.0.20260406.1214/l_params_calibrated_sets_SimCRC_v0.13.0.20260406.1214_Adenoma_F.RData")
l_params_Min_MSE <- l_params_calibrated_sets$Min_MSE
l_params_all <- load_params_init(fromFile = TRUE, filename = l_params_Min_MSE)

# Update the model start age and the survival by race defaults. (CHANGE THE DEFAULTS IN SIMCRC)
l_params_all$min_age_lesion_onset <- 10
l_params_all$mort_by_race <- FALSE
# l_params_all$year_surv_improv <- 2003    # We haven't adjusted this for Chile

# Define the simulation population size
n_pop <- 1e6 # Run at least 1 mil for publications, 10mil if possible for stable results

# Define the cohort age
cohort_age <- 40

# Sample the age of death from life table
df_lt_chile <- read.csv("~/Documents/GitHub/simcrc_chile/data-raw/df_lifetable_2017_CH.csv")
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

n_ids <- nrow(df_screening_strategies)


# *****************************************************************************
#### 3.0 Run screening and surveillance (parallel) ####
# *****************************************************************************

n_cores <- ceiling(parallel::detectCores()/ 2)

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

cl <- makeCluster(n_cores)
registerDoSNOW(cl)

pb   <- txtProgressBar(min = 0, max = n_ids, style = 3, file = stderr())
opts <- list(progress = function(n) {
  setTxtProgressBar(pb, n)
  flush.console()
})

t_parallel_start <- proc.time()
set.seed(3)

foreach(
  i = 1:n_ids,
  .packages     = c("data.table", "simcrc", "ceacrc", "dplyr"),
  .options.snow = opts
) %dopar% {

  run_id <- df_screening_strategies[i, ]

  # Run info
  id              <- run_id$id
  strategy_name   <- run_id$strategy
  screening_years <- seq(run_id$age_to_begin_screening,
                         run_id$age_to_end_screening,
                         run_id$frequency_screening)
  screening_modality <- run_id$modality
  follow_up          <- run_id$follow_up

  start_time <- Sys.time()

  capture.output({

    # Run the screening module
    results_screening <- screening_detection(dt_crc_pop = dt_crc_pop,
                                             p_coverage      = run_id$p_coverage,
                                             p_adherence     = run_id$p_adherence,
                                             age_to_begin_screening = run_id$age_to_begin_screening,
                                             age_to_end_screening = run_id$age_to_end_screening,
                                             frequency_screening = run_id$frequency_screening,
                                             sens_small_adenoma = run_id$sens_small_adenoma,
                                             sens_medium_adenoma = run_id$sens_medium_adenoma,
                                             sens_large_adenoma = run_id$sens_large_adenoma,
                                             sens_crc = run_id$sens_crc,
                                             sens_by_ad = run_id$sens_by_ad,
                                             spec = run_id$spec,
                                             screening_reach = run_id$screening_reach,
                                             p_reach_cecum = run_id$p_reach_cecum,
                                             p_reach_ascending = run_id$p_reach_ascending,
                                             p_reach_transverse = run_id$p_reach_transverse,
                                             p_reach_descending = run_id$p_reach_descending,
                                             p_reach_sigmoid = run_id$p_reach_sigmoid,
                                             p_reach_rectum = run_id$p_reach_rectum,
                                             p_death_scr = run_id$p_death_scr,
                                             surveillance = run_id$surveillance,
                                             confirmation = run_id$follow_up,
                                             p_adherence_confirmation = run_id$p_adherence_confirmation,
                                             sens_small_adenoma_col = run_id$sens_small_adenoma_col,
                                             sens_medium_adenoma_col = run_id$sens_medium_adenoma_col,
                                             sens_large_adenoma_col = run_id$sens_large_adenoma_col,
                                             sens_crc_col = run_id$sens_crc_col,
                                             sens_by_ad_col = run_id$sens_by_ad_col,
                                             spec_col = run_id$spec_col,
                                             confirmation_reach = run_id$confirmation_reach,
                                             p_reach_cecum_conf = run_id$p_reach_cecum_conf,
                                             p_reach_ascending_conf = run_id$p_reach_ascending_conf,
                                             p_reach_transverse_conf = run_id$p_reach_transverse_conf,
                                             p_reach_descending_conf = run_id$p_reach_descending_conf,
                                             p_reach_sigmoid_conf = run_id$p_reach_sigmoid_conf,
                                             p_reach_rectum_conf = run_id$p_reach_rectum_conf,
                                             p_death_conf = run_id$p_death_conf,
                                             optimize_memory = FALSE)

    dt_pop_scr <- results_screening$dt_pop_screening

    # Run the surveillance module
    results_surveillance <- surveillance_detection(dt_pop_scr = dt_pop_scr,
                                                   p_adherence     = run_id$p_adherence_surv,
                                                   sens_small_adenoma = run_id$sens_small_adenoma_surv,
                                                   sens_medium_adenoma = run_id$sens_medium_adenoma_surv,
                                                   sens_large_adenoma = run_id$sens_large_adenoma_surv,
                                                   sens_crc = run_id$sens_crc_surv,
                                                   sens_by_ad = run_id$sens_by_ad_surv,
                                                   spec = run_id$spec_surv,
                                                   surveillance_reach = run_id$surveillance_reach,
                                                   p_reach_cecum_surv = run_id$p_reach_cecum_surv,
                                                   p_reach_ascending_surv = run_id$p_reach_ascending_surv,
                                                   p_reach_transverse_surv = run_id$p_reach_transverse_surv,
                                                   p_reach_descending_surv = run_id$p_reach_descending_surv,
                                                   p_reach_sigmoid_surv = run_id$p_reach_sigmoid_surv,
                                                   p_reach_rectum_surv = run_id$p_reach_rectum_surv,
                                                   p_death_surv = run_id$p_death_surv,
                                                   optimize_memory = TRUE)

    dt_pop_surv <- results_surveillance

    # Summarize the lesion-level file by age (in the standard uspstf format)
    dt_export_final <- uspstf_summary(datatable = dt_pop_surv,
                                      screening_modality = screening_modality,
                                      follow_up = follow_up,
                                      screening_years = screening_years,
                                      min_age = cohort_age,
                                      max_age = 100)

  })

  end_time <- Sys.time()
  execution_time <- difftime(end_time, start_time, units = "mins")

  # Create a header string
  time <- format(Sys.time(), "%Y-%m-%d_%H:%M")
  if (screening_modality == "COL" | screening_modality == "SIG" | screening_modality == "CTC") {
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
  } else if (screening_modality == "FIT" | screening_modality == "sDNA-FIT") {
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
  } else {
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

  # File name (absolute path so workers resolve it correctly)
  uspstf_folder <- normalizePath(
    paste0("output/", run_id$project, "/", run_id$scenario, "/RawModelOutput_SimCRC"),
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

  # Write progress to log file (cat() is unreliable inside parallel workers)
  msg <- sprintf("Strategy %s completed in %.2f minutes", strategy_name, execution_time)
  cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", msg, "\n"),
      file = log_file, append = TRUE)

  NULL
}

close(pb)
cat("\n")  # newline after progress bar
stopCluster(cl)

t_parallel_total <- (proc.time() - t_parallel_start)[["elapsed"]]
log_msg(sprintf("All %d strategies completed in %.1f minutes (%.0f seconds) on %d cores.",
                n_ids, t_parallel_total / 60, t_parallel_total, n_cores))


# *****************************************************************************
###  4.0 Process the model output ---------------------------------------------
# *****************************************************************************

project <- unique(df_screening_strategies$project)
scenarios <- unique(df_screening_strategies$scenario)

df_uspstf_output <- ceacrc::ProcessUSPSTFOutput(analysis_folder = paste0("output/", project, "/", scenarios[1]),
                                        input_folder = "cea_inputs",
                                        first_age_of_interest = cohort_age,
                                        col_infl_rate = 1.05,
                                        discount_rate = 0.03,
                                        col_spec_adj = 0.86,
                                        crc_care_costs_file = "crc_care_costs.csv",
                                        screen_costs_file = "screen_costs_v2.csv",
                                        crc_care_disutility_file = "crc_care_utility_loss.csv",
                                        screen_disutility_file = "screen_utility_loss_WithStoolTestValues.csv",
                                        general_health_utility_weights_file = "GeneralHealthUtilityWeightsByAge.csv",
                                        selected_outcomes_for_model_data_file = "model_data_outcomes_boolean_addDscQALY.csv",
                                        model_run_data_tag = "_Base",
                                        folder_for_output = "ce_results")


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

ggsave(filename = "ce_results/plot_ce_all_strategies.svg", width = 6.5, height = 4, units = "in")


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




ggsave(filename = "ce_results/plot_ce_FIT.svg", width = 6.5, height = 4, units = "in")


# *****************************************************************************
###  7.0 Single-strategy runtime test (FIT4585q1) -----------------------------
# *****************************************************************************

run_id <- df_screening_strategies %>% filter(strategy == "FIT4585q1")

strategy_name      <- run_id$strategy
screening_years    <- seq(run_id$age_to_begin_screening,
                          run_id$age_to_end_screening,
                          run_id$frequency_screening)
screening_modality <- run_id$modality
follow_up          <- run_id$follow_up

t_single_start <- Sys.time()
cat(sprintf("[%s] Starting single-strategy test: %s\n",
            format(t_single_start, "%H:%M:%S"), strategy_name))

results_screening <- screening_detection(dt_crc_pop = dt_crc_pop,
                                         p_coverage      = run_id$p_coverage,
                                         p_adherence     = run_id$p_adherence,
                                         age_to_begin_screening = run_id$age_to_begin_screening,
                                         age_to_end_screening = run_id$age_to_end_screening,
                                         frequency_screening = run_id$frequency_screening,
                                         sens_small_adenoma = run_id$sens_small_adenoma,
                                         sens_medium_adenoma = run_id$sens_medium_adenoma,
                                         sens_large_adenoma = run_id$sens_large_adenoma,
                                         sens_crc = run_id$sens_crc,
                                         sens_by_ad = run_id$sens_by_ad,
                                         spec = run_id$spec,
                                         screening_reach = run_id$screening_reach,
                                         p_reach_cecum = run_id$p_reach_cecum,
                                         p_reach_ascending = run_id$p_reach_ascending,
                                         p_reach_transverse = run_id$p_reach_transverse,
                                         p_reach_descending = run_id$p_reach_descending,
                                         p_reach_sigmoid = run_id$p_reach_sigmoid,
                                         p_reach_rectum = run_id$p_reach_rectum,
                                         p_death_scr = run_id$p_death_scr,
                                         surveillance = run_id$surveillance,
                                         confirmation = run_id$follow_up,
                                         p_adherence_confirmation = run_id$p_adherence_confirmation,
                                         sens_small_adenoma_col = run_id$sens_small_adenoma_col,
                                         sens_medium_adenoma_col = run_id$sens_medium_adenoma_col,
                                         sens_large_adenoma_col = run_id$sens_large_adenoma_col,
                                         sens_crc_col = run_id$sens_crc_col,
                                         sens_by_ad_col = run_id$sens_by_ad_col,
                                         spec_col = run_id$spec_col,
                                         confirmation_reach = run_id$confirmation_reach,
                                         p_reach_cecum_conf = run_id$p_reach_cecum_conf,
                                         p_reach_ascending_conf = run_id$p_reach_ascending_conf,
                                         p_reach_transverse_conf = run_id$p_reach_transverse_conf,
                                         p_reach_descending_conf = run_id$p_reach_descending_conf,
                                         p_reach_sigmoid_conf = run_id$p_reach_sigmoid_conf,
                                         p_reach_rectum_conf = run_id$p_reach_rectum_conf,
                                         p_death_conf = run_id$p_death_conf,
                                         optimize_memory = FALSE)

dt_pop_scr <- results_screening$dt_pop_screening

results_surveillance <- surveillance_detection(dt_pop_scr = dt_pop_scr,
                                               p_adherence     = run_id$p_adherence_surv,
                                               sens_small_adenoma = run_id$sens_small_adenoma_surv,
                                               sens_medium_adenoma = run_id$sens_medium_adenoma_surv,
                                               sens_large_adenoma = run_id$sens_large_adenoma_surv,
                                               sens_crc = run_id$sens_crc_surv,
                                               sens_by_ad = run_id$sens_by_ad_surv,
                                               spec = run_id$spec_surv,
                                               surveillance_reach = run_id$surveillance_reach,
                                               p_reach_cecum_surv = run_id$p_reach_cecum_surv,
                                               p_reach_ascending_surv = run_id$p_reach_ascending_surv,
                                               p_reach_transverse_surv = run_id$p_reach_transverse_surv,
                                               p_reach_descending_surv = run_id$p_reach_descending_surv,
                                               p_reach_sigmoid_surv = run_id$p_reach_sigmoid_surv,
                                               p_reach_rectum_surv = run_id$p_reach_rectum_surv,
                                               p_death_surv = run_id$p_death_surv,
                                               optimize_memory = TRUE)

t_single_end <- Sys.time()
t_single_elapsed <- difftime(t_single_end, t_single_start, units = "mins")
cat(sprintf("[%s] %s completed in %.2f minutes\n",
            format(t_single_end, "%H:%M:%S"), strategy_name, t_single_elapsed))


# *****************************************************************************
###  8.0 Resource monitor — how much RAM and CPU does one strategy use? -------
# Run FIT4585q1 (the heaviest) sequentially and sample RAM + CPU every second.
# From this we derive how many parallel workers your machine can safely support.
# *****************************************************************************

run_id_bench <- df_screening_strategies %>% filter(strategy == "FIT4585q1")

# Snapshot function: RAM used by this R session (MB) + logical CPU count
snapshot_resources <- function() {
  ram_mb <- sum(gc()[, 2]) * 8 / 1024   # gc() returns cells; each cell = 8 bytes on 64-bit
  list(ram_mb = ram_mb, timestamp = Sys.time())
}

# Sample in background every 2 seconds while strategy runs
monitor_file <- tempfile(fileext = ".csv")
monitor_pid  <- sys::exec_background(
  "Rscript",
  args = c("--vanilla", "-e", sprintf(
    'while(TRUE) { cat(format(Sys.time(), \"%%H:%%M:%%S\"), \",\", mem.maxVSize(), \"\\n\", sep=\"\", file=\"%s\", append=TRUE); Sys.sleep(2) }',
    monitor_file
  ))
)

cat("Resource monitor started. Running FIT4585q1 sequentially...\n")
t_res_start <- Sys.time()

invisible(capture.output({
  res_scr <- screening_detection(dt_crc_pop = dt_crc_pop,
                                 p_coverage      = run_id_bench$p_coverage,
                                 p_adherence     = run_id_bench$p_adherence,
                                 age_to_begin_screening = run_id_bench$age_to_begin_screening,
                                 age_to_end_screening = run_id_bench$age_to_end_screening,
                                 frequency_screening = run_id_bench$frequency_screening,
                                 sens_small_adenoma = run_id_bench$sens_small_adenoma,
                                 sens_medium_adenoma = run_id_bench$sens_medium_adenoma,
                                 sens_large_adenoma = run_id_bench$sens_large_adenoma,
                                 sens_crc = run_id_bench$sens_crc,
                                 sens_by_ad = run_id_bench$sens_by_ad,
                                 spec = run_id_bench$spec,
                                 screening_reach = run_id_bench$screening_reach,
                                 p_reach_cecum = run_id_bench$p_reach_cecum,
                                 p_reach_ascending = run_id_bench$p_reach_ascending,
                                 p_reach_transverse = run_id_bench$p_reach_transverse,
                                 p_reach_descending = run_id_bench$p_reach_descending,
                                 p_reach_sigmoid = run_id_bench$p_reach_sigmoid,
                                 p_reach_rectum = run_id_bench$p_reach_rectum,
                                 p_death_scr = run_id_bench$p_death_scr,
                                 surveillance = run_id_bench$surveillance,
                                 confirmation = run_id_bench$follow_up,
                                 p_adherence_confirmation = run_id_bench$p_adherence_confirmation,
                                 sens_small_adenoma_col = run_id_bench$sens_small_adenoma_col,
                                 sens_medium_adenoma_col = run_id_bench$sens_medium_adenoma_col,
                                 sens_large_adenoma_col = run_id_bench$sens_large_adenoma_col,
                                 sens_crc_col = run_id_bench$sens_crc_col,
                                 sens_by_ad_col = run_id_bench$sens_by_ad_col,
                                 spec_col = run_id_bench$spec_col,
                                 confirmation_reach = run_id_bench$confirmation_reach,
                                 p_reach_cecum_conf = run_id_bench$p_reach_cecum_conf,
                                 p_reach_ascending_conf = run_id_bench$p_reach_ascending_conf,
                                 p_reach_transverse_conf = run_id_bench$p_reach_transverse_conf,
                                 p_reach_descending_conf = run_id_bench$p_reach_descending_conf,
                                 p_reach_sigmoid_conf = run_id_bench$p_reach_sigmoid_conf,
                                 p_reach_rectum_conf = run_id_bench$p_reach_rectum_conf,
                                 p_death_conf = run_id_bench$p_death_conf,
                                 optimize_memory = FALSE)
  dt_pop_scr_bench <- res_scr$dt_pop_screening
  peak_ram_mb <- sum(gc()[, 6]) * 8 / 1024  # peak RAM during screening (max used column)

  surveillance_detection(dt_pop_scr = dt_pop_scr_bench,
                         p_adherence     = run_id_bench$p_adherence_surv,
                         sens_small_adenoma = run_id_bench$sens_small_adenoma_surv,
                         sens_medium_adenoma = run_id_bench$sens_medium_adenoma_surv,
                         sens_large_adenoma = run_id_bench$sens_large_adenoma_surv,
                         sens_crc = run_id_bench$sens_crc_surv,
                         sens_by_ad = run_id_bench$sens_by_ad_surv,
                         spec = run_id_bench$spec_surv,
                         surveillance_reach = run_id_bench$surveillance_reach,
                         p_reach_cecum_surv = run_id_bench$p_reach_cecum_surv,
                         p_reach_ascending_surv = run_id_bench$p_reach_ascending_surv,
                         p_reach_transverse_surv = run_id_bench$p_reach_transverse_surv,
                         p_reach_descending_surv = run_id_bench$p_reach_descending_surv,
                         p_reach_sigmoid_surv = run_id_bench$p_reach_sigmoid_surv,
                         p_reach_rectum_surv = run_id_bench$p_reach_rectum_surv,
                         p_death_surv = run_id_bench$p_death_surv,
                         optimize_memory = TRUE)
}))

tools::pskill(monitor_pid)  # stop the background monitor

t_res_elapsed <- difftime(Sys.time(), t_res_start, units = "mins")

# Peak RAM: read actual resident memory of this R process from the OS
pid          <- Sys.getpid()
rss_kb       <- as.numeric(system(paste("ps -o rss= -p", pid), intern = TRUE))
peak_ram_mb  <- rss_kb / 1024
total_ram_mb <- as.numeric(system("sysctl -n hw.memsize", intern = TRUE)) / 1024^2
safe_ram_mb  <- total_ram_mb * 0.80        # leave 20% headroom
phys_cores   <- detectCores(logical = FALSE)
logi_cores   <- detectCores(logical = TRUE)

recommended_cores <- max(1, min(phys_cores, floor(safe_ram_mb / peak_ram_mb)))

cat(sprintf("\n===== Resource Report =====\n"))
cat(sprintf("  Strategy:              FIT4585q1 (heaviest)\n"))
cat(sprintf("  Sequential runtime:    %.2f minutes\n", t_res_elapsed))
cat(sprintf("  Peak RAM (1 worker):   %.0f MB (%.1f GB)\n", peak_ram_mb, peak_ram_mb / 1024))
cat(sprintf("  Total system RAM:      %.0f MB (%.1f GB)\n", total_ram_mb, total_ram_mb / 1024))
cat(sprintf("  Safe RAM budget (80%%): %.0f MB (%.1f GB)\n", safe_ram_mb, safe_ram_mb / 1024))
cat(sprintf("  Physical cores:        %d\n", phys_cores))
cat(sprintf("  Logical cores:         %d\n", logi_cores))
cat(sprintf("  ----------------------------------\n"))
cat(sprintf("  Recommended n_cores:   %d\n", recommended_cores))
cat(sprintf("  (= floor(%.0f MB / %.0f MB per worker), capped at %d physical cores)\n",
            safe_ram_mb, peak_ram_mb, phys_cores))
