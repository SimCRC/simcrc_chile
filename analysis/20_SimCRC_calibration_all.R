#==============================================================================#
#
#  SimCRC Calibration wrapper 
#  
#  Author: Carlos Pineda-Antunez
#  Date: Dec 2025
#  Description: This script runs the full calibration process for SimCRC
#               using BayCANN approach. It includes LHS design, coverage
#               analysis, ANN training, Bayesian calibration using Stan,
#               posterior validation, and selection of best parameter sets.
#               Finally, it generates a calibration summary report.
#
#=============================================================================#

###### 1.Libraries and functions  =============================================

Sys.setenv(RETICULATE_PYTHON = "~/.virtualenvs/r-tensorflow/bin/python")
library(reticulate)
use_virtualenv("r-tensorflow", required = TRUE)
library(keras3)  #Install previously tensorflow
library(rstan)
library(reshape2)
library(tidyverse)
library(ggplot2)
library(doParallel)
library(dplyr)
library(data.table)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

library(devtools)
library(caTools)
library(plyr)
library(simcrc)
library(extraDistr) # To get Gompertz distribution
library(parallel)
library(doParallel)
library(microbenchmark)
library(beepr)

library(doBy)
library(readr)
library(stringr)  # to change graph labels

data.table::getDTthreads()

###### 2. Load functions =================================================

source("analysis/baycann_functions.R")
source("analysis/01_calibration_setup.R")
source("R/01_model_input_functions.R")
source("R/02_calibration_functions.R")
source("R/02_calibration_functions_ssp.R")
source("R/04_BayCANN_setup.R")  #Add specific ANN settings for each type of model
source("R/05_calibration_targets.R")
source("R/06_validation_functions.R")


###### 3. Set version of BayCANN ==============================================
Country           <- "Chile"  # USA or Chile
Machine            <- "Local" # Local or Argonne Sherlock
Model_name         <- "SimCRC"
main_version       <- paste0("v",packageVersion("simcrc"))
Date_version       <- format(Sys.time(), "%Y%m%d.%H%M")


models_to_calibrate <- c("model_adenoma_Chile"
                        #"model_female_adenoma"
                         #"model_male_adenoma",
                         #"model_female_both"
                        # "model_male_both"
                         )

for(models in models_to_calibrate) {

  current_model <- "model_adenoma_Chile" # For testing purposes
  #current_model <- models
  
  # Check what is the current model to calibrate
  if(current_model == "model_female_adenoma") {
    calibration_setup <- l_model_female_adenoma
  }  
  if(current_model == "model_male_adenoma") {
    calibration_setup <- l_model_male_adenoma
  } 
  if(current_model == "model_female_both") {
    calibration_setup <- l_model_female_both
  } 
  if(current_model == "model_male_both") {
    calibration_setup <- l_model_male_both
  } 
  if (current_model == "model_adenoma_Chile") {
    calibration_setup <- l_model_adenoma_Chile
  }
  
  #Get parameters
  
  Population_type    <- calibration_setup$population_type  
  Lesion_type        <- calibration_setup$lesion_type
  Calibration_number <- Date_version # Date version works to match multiple calibrations
  cal_version        <- paste0(main_version,".",Calibration_number)
  BayCANN_version    <- paste0(Model_name,"_",cal_version, "_",Lesion_type,"_",Population_type)
  folder_version     <- paste0("outputs/BayCANN_versions/",Country, "/",Lesion_type,"/",Population_type,"/",main_version)
  folder             <- paste0(folder_version, "/", cal_version) #Folder for specific version
  
  ## Set all the required paths for this version of BayCANN and save them in a list
  paths_calibration <- list(
    Model_name              = Model_name,
    Machine                 = Machine,
    Country                 = Country,
    main_version            = main_version,
    Lesion_type             = Lesion_type,
    Population_type         = Population_type,
    Date_version            = Date_version,
    Calibration_number      = Calibration_number,
    cal_version             = cal_version,
    BayCANN_version         = BayCANN_version,
    folder                  = folder,
    path_description        = paste0(folder,"/Calibration_description_",BayCANN_version,".txt"),
    path_priors             = paste0(folder,"/l_params_priors_SimCRC_",BayCANN_version,".RData"),
    path_lhs                = paste0(folder,"/l_LHS_SimCRC_",BayCANN_version,".RData"),
    path_coverage_num       = paste0(folder,"/fig_coverage_num_",BayCANN_version,".png"),
    path_coverage_cat       = paste0(folder,"/fig_coverage_cat_",BayCANN_version,".png"),
    path_keras_model        = paste0(folder,"/o_keras_model_",BayCANN_version,".keras"),
    path_ann_perform        = paste0(folder,"/fig_ANN_prediction_performance_",BayCANN_version,".png"),
    path_targets            = calibration_setup$targest_file,
    path_stan_model         = paste0(folder,"/o_stan_model_",BayCANN_version,".rds"),
    path_posteriors         = paste0(folder,"/dt_calibrated_posteriors_",BayCANN_version,".csv"),
    path_post_chains        = paste0(folder,"/fig_posterior_distribution_chains_",BayCANN_version,".png"),
    path_post_pairs         = paste0(folder,"/fig_posterior_distribution_pairwise_corr_",BayCANN_version,".png"),
    path_prior_post         = paste0(folder,"/fig_prior-vs-posterior_",BayCANN_version,".png"),
    path_val_ANN            = paste0(folder,"/prediction_ANN_posteriors_",BayCANN_version,".csv"),
    path_baycann_params     = paste0(folder,"/l_BayCANN_parameters_",BayCANN_version,".RData"),
    path_posterior_outputs  = paste0(folder,"/df_posterior_outputs_SimCRC_",BayCANN_version,".rda"),
    path_posterior_val_num  = paste0(folder,"/fig_internall_validation_",BayCANN_version,"_num.png"),
    path_posterior_val_cat  = paste0(folder,"/fig_internall_validation_",BayCANN_version,"_cat.png"),
    path_best_params_sets   = paste0(folder,"/l_params_calibrated_sets_", BayCANN_version, ".RData"),
    path_set_param_val_num  = paste0(folder,"/fig_parameter_set_validation_num_",BayCANN_version,".png"),
    path_set_param_val_cat  = paste0(folder,"/fig_parameter_set_validation_cat_",BayCANN_version,".png"),
    path_summary_template   = "outputs/BayCANN_versions/Both/Calibration_summary.qmd",
    path_calibration_summary = paste0(folder,"/Calibration_summary_",BayCANN_version,".qmd")
  )
  
  ###### 4. Create folder for current calibration version =======================
  
  
  # --- 1) Ensure natural history version folder exists --------------------------
  if (!dir.exists(folder_version)) {
    ok <- dir.create(folder_version, recursive = TRUE, showWarnings = FALSE)
    if (!ok) stop(sprintf("Failed to create folder_version '%s'.", folder_version))
    cat(sprintf("The folder version '%s' was created.\n", folder_version))
  } else {
    cat(sprintf("The folder version '%s' already exists.\n", folder_version))
  }
  
  # --- 2) Ensure calibration folder exists --------------------------------------
  if (!dir.exists(folder)) {
    ok <- dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    if (!ok) stop(sprintf("Failed to create folder '%s'.", folder))
    cat(sprintf("New calibration folder '%s' was created.\n", folder))
  } else {
    cat(sprintf("The folder '%s' already exists. Verify calibration number.\n", folder))
  }
  
  # --- 3) Ensure parent folder for description exists (extra safety) ------------
  desc_dir <- dirname(paths_calibration$path_description)
  if (!dir.exists(desc_dir)) {
    ok <- dir.create(desc_dir, recursive = TRUE, showWarnings = FALSE)
    if (!ok) stop(sprintf("Failed to create description directory '%s'.", desc_dir))
  }
  
  # --- 4) Write characteristics of this version of LHS/Calibration --------------
  version_particularity <- paste0(
    "First Chilean SimCRC model for Chilean population with priors from the USA corrected in April from Nicolas suggestions\n",
    "We will use the priors as our priors from the US and modify parameters regarding alpha_lesion_adenoma, hazard rates and probabilites from preclinical to clinical detection.\n",
    "Targets used are from Chilean data as described in the calibration setup.\n"
  )
  
  text_log <- c(
    paste0("Model name: ", Model_name),
    "",
    paste0("Model version: ", main_version),
    "",
    paste0("Date of calibration: ", Date_version),
    "",
    paste0("BayCANN version: ", BayCANN_version),
    "",
    paste0("Folder: ", folder),
    "",
    paste0("Models calibrated in this version: ", models_to_calibrate),
    "",
    version_particularity
  )
  
  # writeLines can take a file path directly; no need to open/close a connection
  writeLines(text_log, con = paths_calibration$path_description)
  
  cat(sprintf("Wrote calibration description to '%s'.\n", paths_calibration$path_description))
  
  ###### 5. Get base population =======================
  
  SSP_cal <- calibration_setup$SSP_pathway #TRUE if we want to include SSP pathway
  
  # Get base population
  n_pop <- calibration_setup$n_pop
   
  dt_pop <- get_dt_population(year = calibration_setup$cohort_year,
                              byear = calibration_setup$birth_year,
                              p_female = calibration_setup$p_female,
                              p_white = calibration_setup$p_white,
                              n_pop = calibration_setup$n_pop,
                              dt_life_table_F = calibration_setup$life_table_fem,
                              dt_life_table_M = calibration_setup$life_table_male)
  
  ###### 6. Generate LHS =========================================================
  
  # General parameters for LHS
  
  # Load all initial parameters and set the calibrated parameters from previous versions
  l_params_init <- calibration_setup$l_params_init
  
  #Set minimum age of lesion onset according to the calibrated model
  l_params_init$min_age_lesion_onset <- calibration_setup$min_age_lesion_onset
  
  # Do not use race to estimate mortality
  l_params_init$mort_by_race <- FALSE
  
  # Number of simulations
  n_sim <- 20000
  
  #Do Parallel
  parallel <- TRUE
  
  #Save LHS
  saveLHS <- TRUE
  
  #Number of cores for LHS
  #check number of cores
  if (parallel) {
    no_cores <- ceiling(parallel::detectCores()/ 2)
  } else {
    no_cores <- 1
  }
  
  # Set seed
  set.seed(20220906)



  calibration_setup <- l_model_adenoma_Chile
  
  source("analysis/02_LHS_design_all.R")  #Run LHS design
  
  #Save results
  char_data_LHS_filename <- paths_calibration$path_lhs
  if (saveLHS){
    data_LHS <- list(df_simcrc_outputs, df_params_lhs)
    save(data_LHS, file =  char_data_LHS_filename) 
    print(paste("LHS saved in: ",char_data_LHS_filename))
  }
  
  
  ###### 7. Coverage analysis ====================================================
  
  save_graph <- TRUE
  
  targets_file <- calibration_setup$targest_file
  
  source("analysis/03_Coverage_analysis.R")
  
  
  ###### 8. BayCANN calibration ==================================================

  # Save session before TensorFlow step (restart R & load before continuing)
  save.image(file = paste0(folder, "/session_pre_baycann_", BayCANN_version, ".RData"))

  load(paste0("outputs/BayCANN_versions/Chile/Adenoma/F/v0.13.0/", "v0.13.0.20260401.1325", "/session_pre_baycann_", "SimCRC_v0.13.0.20260401.1325_Adenoma_F", ".RData"))


  source("analysis/06_BayCANN_calibration_all_k3.R")
  
  ###### 9. Posterior validation =================================================
  
  #Validation specifications
  
  # Set seed
  set.seed(20220906)
  
  #Do Parallel
  parallel <- TRUE
  
  #Save validation output
  saveOutput <- TRUE
  
  #Save in calibration folder
  save_in_BayCANN <- TRUE
  
  #check number of cores
  if (parallel) {
    no_cores <- ceiling(parallel::detectCores())/ 2
  } else {
    no_cores <- 1
  }
  
  source("analysis/07_posterior_validations.R")
  
  
  # Get validation graphs
  
  chains_to_include <- c(1, 2,3,4)  #Specify the chains to include
  
  source("analysis/07_1_posterior_validations_graphs.R")
  
  ###### 9.1 Dwell, sojourn time 
  
  
  # Select data
  df_dwell_sojourn_CH <- df_model_outputs %>% 
    dplyr::select(chain, Adenoma_dwell, Sojourn_time, Total_dwell)
  
  
  load("data-raw/df_posterior_outputs_SimCRC_SimCRC_v0.12.0.1_Ad_F.rda")
  
  df_dwell_sojourn_US <- df_simcrc_outputs %>% 
    dplyr::select(chain, Adenoma_dwell, Sojourn_time, Total_dwell)
  
  # Function to summarize dwell/sojourn data
  summarize_dwell <- function(df, country_name) {
    df %>%
      dplyr::summarise(
        dplyr::across(c(Adenoma_dwell, Sojourn_time, Total_dwell),
                      list(
                        mean = ~mean(.x, na.rm = TRUE),
                        median = ~median(.x, na.rm = TRUE),
                        lower_95 = ~quantile(.x, 0.025, na.rm = TRUE),
                        upper_95 = ~quantile(.x, 0.975, na.rm = TRUE),
                        lower_50 = ~quantile(.x, 0.25, na.rm = TRUE),
                        upper_50 = ~quantile(.x, 0.75, na.rm = TRUE)
                      ),
                      .names = "{.col}||{.fn}")
      ) %>% 
      tidyr::pivot_longer(everything(),
                          names_to = c("metric", "stat"),
                          names_sep = "\\|\\|",
                          values_to = "value") %>% 
      tidyr::pivot_wider(names_from = stat, values_from = value) %>%
      mutate(
        country = country_name,
        metric_label = case_when(
          metric == "Adenoma_dwell" ~ "Adenoma Dwell",
          metric == "Sojourn_time" ~ "Sojourn Time",
          metric == "Total_dwell" ~ "Total Dwell"
        )
      )
  }
  
  # Summarize both datasets
  df_dwell_sojourn_CH_sum <- summarize_dwell(df_dwell_sojourn_CH, "Chile (SimCRC Chile v0.12.1 All)")
  df_dwell_sojourn_US_sum <- summarize_dwell(df_dwell_sojourn_US, "United States (SimCRC v0.12.0 Female)")
  
  # Combine
  df_dwell_sojourn_combined <- bind_rows(df_dwell_sojourn_CH_sum, 
                                         df_dwell_sojourn_US_sum)
  
  print(df_dwell_sojourn_combined)
  
  # Plot comparison with facet wrap by metric
  ggplot(df_dwell_sojourn_combined, 
         aes(x = country, y = mean, color = country)) +
    geom_point(size = 4, position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = lower_95, ymax = upper_95), 
                  width = 0.2, linewidth = 1, alpha = 0.6,
                  position = position_dodge(width = 0.5)) +
    # Add mean values with white background ON the point
    geom_label(aes(label = sprintf("%.1f", mean)),
               position = position_dodge(width = 0.5),
               vjust = 0.5, hjust = 0.5, size = 3, fontface = "bold",
               label.padding = unit(0.15, "lines"),
               label.size = 0,
               fill = "white",
               show.legend = FALSE) +
    # Add 95% CI labels (upper) with white background
    geom_label(aes(y = upper_95, label = sprintf("%.1f", upper_95)),
               position = position_dodge(width = 0.5),
               vjust = -0.5, size = 2.8,
               label.padding = unit(0.1, "lines"),
               label.size = 0,
               fill = "white",
               show.legend = FALSE) +
    # Add 95% CI labels (lower) with white background
    geom_label(aes(y = lower_95, label = sprintf("%.1f", lower_95)),
               position = position_dodge(width = 0.5),
               vjust = 1.5, size = 2.8,
               label.padding = unit(0.1, "lines"),
               label.size = 0,
               fill = "white",
               show.legend = FALSE) +
    facet_wrap(~metric_label, ncol = 3, scales = "free_y") +
    labs(
      title = "Dwell and Sojourn Time Estimates: Chile vs United States",
      subtitle = "Mean with 95% CrI | Values shown in years",
      x = "",
      y = "Time (years)",
      color = "Country"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "bottom",
      strip.text = element_text(face = "bold", size = 11),
      strip.background = element_rect(fill = "grey90", color = NA)
    ) +
    scale_color_manual(values = c("Chile (SimCRC Chile v0.12.1 All)" = "#E41A1C", 
                                  "United States (SimCRC v0.12.0 Female)" = "#377EB8")) +
    scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
    ggview::canvas(width = 10, height = 5)
  
  ggsave(filename = paste0(folder,"/fig_dwell_sojourn_comparison_",BayCANN_version,".png"),
         width = 10, height = 5, dpi = 300)
  
  ggsave(filename = paste0("outputs/BayCANN_versions/Chile/Adenoma/F/v0.12.1/v0.12.1.20260114.1804/fig_dwell_sojourn_comparison_",BayCANN_version,".png"),
         width = 10, height = 5, dpi = 300)
  ###### 10. Get the calibrated set of parameters ================================
  
  
  calibrated_params <- read.csv("outputs/BayCANN_versions/Chile/Adenoma/F/v0.12.1/v0.12.1.20260114.1804/dt_calibrated_posteriors_SimCRC_v0.12.1.20260114.1804_Adenoma_F.csv")
  Baycann_version <- BayCANN_version
  
  source("analysis/12_best_param_set.R")
  
   # Save "paths_calibration" list as .RData file in the calibration folder
  save(paths_calibration, file = paste0(folder,"/paths_calibration.RData"))
  
  # Save "calibration_setup" list as .RData file in the calibration folder
  save(calibration_setup, file = paste0(folder,"/calibration_setup.RData"))

} #End of loop for models to calibrate

Sys.time()
# 11. Combined outcomes -------------------------------------------------------


simcrc_targets <- read.csv(paths_calibration$path_targets)
simcrc_targets <- as.data.table(simcrc_targets)

source("analysis/16_combined_outcomes_SSP.R")

# 12. Calibration summary (Report on Quarto) ----------------------------------

# Load necessary library
library(fs)

# Define the source and destination paths
source_file <- paths_calibration$path_summary_template
destination_folder <- folder
destination_file <- file.path(destination_folder, basename(source_file))

# Copy the file
file_copy(source_file, destination_file)

#Change name of the file
file_move(destination_file, paths_calibration$path_calibration_summary)

# Open the file
# This will open the file in the default application associated with .qmd files
browseURL(paths_calibration$path_calibration_summary)

