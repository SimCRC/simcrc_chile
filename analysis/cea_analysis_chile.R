# *****************************************************************************
#
# Script: cea_analysis_chile.R
#
# Purpose: Run Cost-Effectiveness Analysis (CEA) for CRC screening strategies
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

library(ceacrc)


# *****************************************************************************
#### 1.0 Natural history ####
# *****************************************************************************

# Get the model version
simcrc_model_version <- paste0("SimCRC v",as.character(packageVersion("simcrc")))

# Load calibrated parameters 
load("data-raw/l_params_calibrated_sets_SimCRC_v0.12.1.20260114.1804_Adenoma_F.RData")
l_params_Min_AbsolutErr <- l_params_calibrated_sets$Min_AbsolutErr
l_params_all <- load_params_init(fromFile = TRUE, filename = l_params_Min_AbsolutErr)

# Update the model start age and the survival by race defaults. (CHANGE THE DEFAULTS IN SIMCRC)
l_params_all$min_age_lesion_onset <- 10
l_params_all$mort_by_race <- FALSE

# Define the simulation population size
n_pop <- 1e6#1e7

# Define the cohort age
cohort_age <- 40

# Sample the age of death from life table
df_lt_chile <- readRDS("data-raw/df_life_table_2017_CH.rds")
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
df_screening_strategies <- readr::read_csv("data-raw/df_CH_2026_strategies.csv",show_col_types = FALSE)

# Add function to remove duplicated strategies (e.g. COL_45_75_10 & COL_45_80_10)


# Check that there are no duplicate ids
if(any(duplicated(df_screening_strategies$id))){
  stop("There are duplicate ids in the screening strategies. Please check the input file and remove duplicates.")
}

n_ids <- nrow(df_screening_strategies)


# *****************************************************************************
#### 3.0 Run screening and surveillance ####
# *****************************************************************************
n_ids <- 2 # Memory limit reached when running the FIT strategy
pb <- txtProgressBar(min = 0, max = n_ids, style = 3)
set.seed(3)
for (i in c(1,4)) { #1:10) { #1:n_ids) {
  #i=4
  
  run_id <- df_screening_strategies[i,]
  
  # Run info
  id <- run_id$id
  strategy_name <- run_id$strategy
  screening_years <- seq(run_id$age_to_begin_screening,
                         run_id$age_to_end_screening,
                         run_id$frequency_screening)
  screening_modality <- run_id$modality
  follow_up <- run_id$follow_up
  
  start_time <- Sys.time() # Start time for this iteration
  
  cat("\nEvaluating", id, "\n")
  
  capture.output({
    
    # Run the screening module
    system.time(
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
                                               p_death_scr = run_id$p_death_scr, #probability of death from screening
                                               surveillance = run_id$surveillance, 
                                               # Variables below are for follow-up colonoscopy 
                                               confirmation = run_id$follow_up, #Set to TRUE to FIT tests
                                               p_adherence_confirmation = run_id$p_adherence_confirmation, # Updated lists from Fernando's tutorial.
                                               sens_small_adenoma_col = run_id$sens_small_adenoma_col, #sensitivity from follow-up colonoscopy
                                               sens_medium_adenoma_col = run_id$sens_medium_adenoma_col, #sensitivity from follow-up colonoscopy
                                               sens_large_adenoma_col = run_id$sens_large_adenoma_col, #sensitivity from follow-up colonoscopy
                                               sens_crc_col = run_id$sens_crc_col, #sensitivity crc from follow-up colonoscopy
                                               sens_by_ad_col = run_id$sens_by_ad_col,
                                               spec_col = run_id$spec_col, #spec from follow-up colonoscopy (Please do not change this assumption. Specificity is applied later in the code)
                                               confirmation_reach = run_id$confirmation_reach,
                                               p_reach_cecum_conf = run_id$p_reach_cecum_conf,                      
                                               p_reach_ascending_conf = run_id$p_reach_ascending_conf,   
                                               p_reach_transverse_conf = run_id$p_reach_transverse_conf,   
                                               p_reach_descending_conf = run_id$p_reach_descending_conf,         
                                               p_reach_sigmoid_conf = run_id$p_reach_sigmoid_conf,                
                                               p_reach_rectum_conf = run_id$p_reach_rectum_conf,
                                               p_death_conf = run_id$p_death_conf,
                                               optimize_memory = FALSE))
    
    dt_pop_scr <- results_screening$dt_pop_screening
    
    # Run the surveillance module
    system.time(
      results_surveillance <- surveillance_detection(dt_pop_scr = dt_pop_scr, 
                                                     p_adherence     = run_id$p_adherence_surv,
                                                     # Sensitivities of colonoscopy
                                                     sens_small_adenoma = run_id$sens_small_adenoma_surv,
                                                     sens_medium_adenoma = run_id$sens_medium_adenoma_surv,
                                                     sens_large_adenoma = run_id$sens_large_adenoma_surv,
                                                     sens_crc = run_id$sens_crc_surv,
                                                     sens_by_ad = run_id$sens_by_ad_surv, 
                                                     spec = run_id$spec_surv, #spec from follow-up colonoscopy (Please do not change this assumption. Specificity is applied later in the code)
                                                     surveillance_reach = run_id$surveillance_reach,
                                                     p_reach_cecum_surv = run_id$p_reach_cecum_surv,  
                                                     p_reach_ascending_surv = run_id$p_reach_ascending_surv,  
                                                     p_reach_transverse_surv = run_id$p_reach_transverse_surv,  
                                                     p_reach_descending_surv = run_id$p_reach_descending_surv,  
                                                     p_reach_sigmoid_surv = run_id$p_reach_sigmoid_surv,  
                                                     p_reach_rectum_surv = run_id$p_reach_rectum_surv,
                                                     p_death_surv = run_id$p_death_surv,
                                                     optimize_memory = TRUE))
    
    dt_pop_surv <- results_surveillance
    
    # Summarize the lesion-level file by age (in the standard uspstf format)
    dt_export_final <- uspstf_summary(datatable = dt_pop_surv,
                                      screening_modality = screening_modality,
                                      follow_up = follow_up,
                                      screening_years = screening_years,
                                      min_age = cohort_age, 
                                      max_age = 100)
    
  })
  
  end_time <- Sys.time() # End time for this iteration
  execution_time <- difftime(end_time, start_time, units = "mins") # Time difference

  # Create a header string
  time <- format(Sys.time(), "%Y-%m-%d_%H:%M")
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

  # File name
  uspstf_folder <- paste0("output/",run_id$project,"/",run_id$scenario,"/RawModelOutput_SimCRC")
  if(!dir.exists(uspstf_folder)) {
    dir.create(uspstf_folder, recursive = TRUE)
  }
  uspstf_filename <- paste0(uspstf_folder,"/",run_id$project ,"_",strategy_name, ".csv")
  
  # Write the header to the file
  writeLines(header_string, uspstf_filename)
  
  # Append the data to the same file
  write.table(dt_export_final,
              file = uspstf_filename,
              col.names = T, row.names = F, quote = F, append = T, sep = ',')
  
  # Display execution time
  cat(sprintf("Strategy %s completed in %.2f minutes\n", strategy_name, execution_time))
  
  setTxtProgressBar(pb, i)
  
}

close(pb)

# *****************************************************************************
###  4.0 Process the model output ---------------------------------------------
# *****************************************************************************

project <- unique(df_screening_strategies$project)
scenarios <- unique(df_screening_strategies$scenario)

df_uspstf_output <- ceacrc::ProcessUSPSTFOutput(analysis_folder = paste0("2020SDA/R1"), ## This should be the folder where subfolder "RawModelOutput..." lives. A
                                        input_folder = "cea_inputs",
                                        first_age_of_interest = cohort_age, # What is the 1st age of interest for analyses? 
                                        col_infl_rate = 1.05,  # CISNET models do not account for incomplete colonoscopies. Instead we assume 5% need to be repeated (due to poor prep, to try again at reaching cecum, etc)
                                        discount_rate = 0.03,
                                        col_spec_adj = 0.86, # CISNET models do not currently simulate non-adenomatous polyps. We account for their detection and removal in post-processing using the lack of specificity.
                                        crc_care_costs_file = "crc_care_costs.csv", # care costs input  (within input_folder)
                                        screen_costs_file = "screen_costs.csv", # screen costs input (within input_folder)
                                        crc_care_disutility_file = "crc_care_utility_loss.csv", # care utilities (within input_folder)
                                        screen_disutility_file = "screen_utility_loss_WithStoolTestValues.csv", # screen utilities (within input_folder)
                                        general_health_utility_weights_file = "GeneralHealthUtilityWeightsByAge.csv", # age utilities (within input_folder)
                                        selected_outcomes_for_model_data_file = "model_data_outcomes_boolean_addDscQALY.csv", # what outcomes do you want (within input_folder)            
                                        model_run_data_tag = "_basecosts", # You can add a tag to the model_run_data file!
                                        folder_for_output = "output")


# *****************************************************************************
###  5.0 Perform the CEA ------------------------------------------------
# *****************************************************************************

# Pick your cost variable (discounted costs) from the uspstf_output
v_crc_costs <- df_uspstf_output$DiscountedTotalCostsper1000

# Pick your benefit variable (discounted QALYG) from the uspstf_output
v_crc_qalys <- df_uspstf_output$DiscountedQALYGainedper1000

# Calculate icers using the vector of costs and vector of qalys
icer_hiv_mx_ar <- calculate_icers(cost = v_crc_costs,
                                  effect = v_crc_qalys,
                                  strategies = df_uspstf_output$Strategy)

writexl::write_xlsx(icer_hiv_mx_ar, path = "R_output/ICERs.xlsx")

# Plot the efficient frontier two ways
dampack:::plot.icers(icer_hiv_mx_ar, label = c("all"))


