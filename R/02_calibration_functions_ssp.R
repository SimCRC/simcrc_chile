#' Generate model outputs for calibration from a parameter set
#'
#' \code{calibration_out_ssp} computes model outputs to be used for calibration 
#' routines.
#'
#' @param v_params_calib Vector of parameters that need to be calibrated.
#' @param l_params_all List with all parameters of the decision model.
#' @param dt_pop Datatable dataset with the cohort information
#' @return 
#' 
#' A datatable with the outputs of the model
#' 
#' @export
calibration_out_ssp <- function(v_params_calib, l_params_all, dt_pop, id_draw = i, dt_long = FALSE){ # User defined
  # Substitute values of calibrated parameters in base-case with 
  # calibrated values
  l_params_all <- update_param_list(l_params_all   = l_params_all, 
                                    params_updated = v_params_calib)
  

  # Run simcrc nathist WITH SSP Pathway
  l_out_simcrc <- simcrc::simcr_nathist_ssp_DES(l_params_all = l_params_all, 
                                                dt_pop       = dt_pop, 
                                                SSP_pathway = TRUE)
  
  
  # Get main outputs of simcrc_nathist
  dt_crc_pop <- l_out_simcrc$dt_crc_pop
  
  # Get living population data by age, sex, smoking and year
  dt_living_pop <- get_living_pop(dt_crc_pop= dt_crc_pop,
                                  min_age = 0, 
                                  max_age = 100)
  # Calculate adenoma prevalence
  dt_adenoma_prev <- simcrc::get_adenoma_prevalence_ssp(
    dt_crc_pop = dt_crc_pop,
    dt_living_pop = dt_living_pop,
    v_lesion_location = c("P", "D", "R"),
    v_age = c(25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95),
    lesion_pathway = "Adenoma",
    age_groups = TRUE
  )
  
  # Calculate SSP prevalence considering misclassification rates
  dt_SSP_prev <- simcrc::get_adenoma_prevalence_ssp(
    dt_crc_pop = dt_crc_pop,
    dt_living_pop = dt_living_pop,
    v_lesion_location = c("P", "D", "R"),
    v_age = c(25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95),
    lesion_pathway = "SSP",
    age_groups = TRUE,
    misrate = TRUE,
    p_small  = 0.69,
    p_medium = 0.81,
    p_large  = 0.91
  )
  
  # Calculate CRC incidence overall
  
  dt_crc_inc_overall <- get_incidence_crc(
    dt_crc_pop = dt_crc_pop,
    dt_living_pop = dt_living_pop,
    ad_location = c("P", "D", "R"),
    v_age = c(20, 40, 50, 60, 70, 80),
    age_groups = TRUE)
  
  
  # Calculate CRC incidence in P
  
  dt_crc_inc_P <- get_incidence_crc(
    dt_crc_pop = dt_crc_pop,
    dt_living_pop = dt_living_pop,
    ad_location = c("P"),
    v_age = c(20, 40, 50, 60, 70, 80),
    age_groups = TRUE)
  
  # Calculate Stage distribution in P
  
  dt_stage_dist_P  <- get_stage_distribution(dt_crc_inc = dt_crc_inc_P)
  
  # Calculate CRC incidence in D
  
  dt_crc_inc_D <- get_incidence_crc(
    dt_crc_pop = dt_crc_pop,
    dt_living_pop = dt_living_pop,
    ad_location = c("D"),
    v_age = c(20, 40, 50, 60, 70, 80),
    age_groups = TRUE)
  
  # Calculate Stage distribution in D
  
  dt_stage_dist_D  <- get_stage_distribution(dt_crc_inc = dt_crc_inc_D) 
  
  # Calculate CRC incidence in R
  
  dt_crc_inc_R <- get_incidence_crc(
    dt_crc_pop = dt_crc_pop,
    dt_living_pop = dt_living_pop,
    ad_location = c("R"),
    v_age = c(20, 40, 50, 60, 70, 80),
    age_groups = TRUE)
  
  # Calculate Stage distribution in R
  
  dt_stage_dist_R  <- get_stage_distribution(dt_crc_inc = dt_crc_inc_R) 
  
  # Calculate preclinical prevalence for adenomas
  
  dt_prev_preclinical <- simcrc::get_prevalence_preclinical_ssp (dt_crc_pop= dt_crc_pop, 
                                                     dt_living_pop = dt_living_pop,
                                                     v_lesion_location = c("P", "D", "R"),
                                                     v_age = c(40,50,60,70,80,90),
                                                     lesion_pathway = "All",
                                                     age_groups  = TRUE)
  
  # Calculate preclinical prevalence for SSP
  
  dt_prev_preclinical_SSP <- simcrc::get_prevalence_preclinical_ssp (dt_crc_pop= dt_crc_pop, 
                                                     dt_living_pop = dt_living_pop,
                                                     v_lesion_location = c("P", "D", "R"),
                                                     v_age = c(40,50,60,70,80,90),
                                                     lesion_pathway = "SSP",
                                                     age_groups  = TRUE)
  
  # Calculate size distribution 
  
  dt_prop_size_P <- simcrc::get_size_distribution_ssp (dt_crc_pop = dt_crc_pop, 
                                           dt_living_pop = dt_living_pop,
                                           ad_location = c("P"),
                                           v_age = c(50:75),
                                           lesion_pathway = "Adenoma",
                                           age_groups  = FALSE)
  
  dt_prop_size_P_SSP <- simcrc::get_size_distribution_ssp (dt_crc_pop = dt_crc_pop, 
                                                       dt_living_pop = dt_living_pop,
                                                       ad_location = c("P"),
                                                       v_age = c(50:75),
                                                       lesion_pathway = "SSP",
                                                       age_groups  = FALSE)
  
  dt_prop_size_D <- simcrc::get_size_distribution_ssp (dt_crc_pop = dt_crc_pop, 
                                           dt_living_pop = dt_living_pop,
                                           ad_location = c("D"),
                                           v_age = c(50:75),
                                           lesion_pathway = "Adenoma",
                                           age_groups  = FALSE)
  
  dt_prop_size_D_SSP <- simcrc::get_size_distribution_ssp (dt_crc_pop = dt_crc_pop, 
                                                       dt_living_pop = dt_living_pop,
                                                       ad_location = c("D"),
                                                       v_age = c(50:75),
                                                       lesion_pathway = "SSP",
                                                       age_groups  = FALSE)
  
  dt_prop_size_R <- simcrc::get_size_distribution_ssp (dt_crc_pop = dt_crc_pop, 
                                           dt_living_pop = dt_living_pop,
                                           ad_location = c("R"),
                                           v_age = c(50:75),
                                           lesion_pathway = "Adenoma",
                                           age_groups  = FALSE)
  
  dt_prop_size_R_SSP <- simcrc::get_size_distribution_ssp (dt_crc_pop = dt_crc_pop, 
                                                       dt_living_pop = dt_living_pop,
                                                       ad_location = c("R"),
                                                       v_age = c(50:75),
                                                       lesion_pathway = "SSP",
                                                       age_groups  = FALSE)
  
  # proportion of CRC by lesion type
  
  prop_SSP <- prop.table(table(dt_crc_pop[crc_dx_i==1, lesion_type]))[2]
  
  # Calculates adenoma dwell time
  
  Adenoma_dwell <- mean(dt_crc_pop[lesion_type=="Adenoma",]$t_Adonset_Preclin, na.rm = TRUE)
  
  Sojourn_time  <- mean(dt_crc_pop[lesion_type=="Adenoma",]$t_Preclin_Sxdet, na.rm = TRUE)
  
  Total_dwell   <- mean(dt_crc_pop[lesion_type=="Adenoma",]$t_Adonset_Sxdet, na.rm = TRUE)
  
  
  dt_output <- data.table::data.table(id_draw                    = id_draw,
                                      Prev0AdOrPreClin_age27     = dt_adenoma_prev[age_group == '[25,30)' , prev_0_adenoma],
                                      Prev0AdOrPreClin_age32     = dt_adenoma_prev[age_group == '[30,35)' , prev_0_adenoma],
                                      Prev0AdOrPreClin_age37     = dt_adenoma_prev[age_group == '[35,40)' , prev_0_adenoma],
                                      Prev0AdOrPreClin_age42     = dt_adenoma_prev[age_group == '[40,45)' , prev_0_adenoma],
                                      Prev0AdOrPreClin_age47     = dt_adenoma_prev[age_group == '[45,50)' , prev_0_adenoma],
                                      Prev0AdOrPreClin_age52     = dt_adenoma_prev[age_group == '[50,55)' , prev_0_adenoma],
                                      Prev0AdOrPreClin_age57     = dt_adenoma_prev[age_group == '[55,60)' , prev_0_adenoma],
                                      Prev0AdOrPreClin_age62     = dt_adenoma_prev[age_group == '[60,65)' , prev_0_adenoma],
                                      Prev0AdOrPreClin_age67     = dt_adenoma_prev[age_group == '[65,70)' , prev_0_adenoma],
                                      Prev0AdOrPreClin_age72     = dt_adenoma_prev[age_group == '[70,75)' , prev_0_adenoma],
                                      Prev0AdOrPreClin_age77     = dt_adenoma_prev[age_group == '[75,80)' , prev_0_adenoma],
                                      Prev0AdOrPreClin_age82     = dt_adenoma_prev[age_group == '[80,85)' , prev_0_adenoma],
                                      Prev0AdOrPreClin_age87     = dt_adenoma_prev[age_group == '[85,90)' , prev_0_adenoma],
                                      Prev0AdOrPreClin_age92     = dt_adenoma_prev[age_group == '[90,95)' , prev_0_adenoma],
                                      Prev0AdOrPreClin_age97     = dt_adenoma_prev[age_group == '[95,Inf)' , prev_0_adenoma],
                                      Prev1AdOrPreClin_age27     = dt_adenoma_prev[age_group == '[25,30)' , prev_1_adenoma],
                                      Prev1AdOrPreClin_age32     = dt_adenoma_prev[age_group == '[30,35)' , prev_1_adenoma],
                                      Prev1AdOrPreClin_age37     = dt_adenoma_prev[age_group == '[35,40)' , prev_1_adenoma],
                                      Prev1AdOrPreClin_age42     = dt_adenoma_prev[age_group == '[40,45)' , prev_1_adenoma],
                                      Prev1AdOrPreClin_age47     = dt_adenoma_prev[age_group == '[45,50)' , prev_1_adenoma],
                                      Prev1AdOrPreClin_age52     = dt_adenoma_prev[age_group == '[50,55)' , prev_1_adenoma],
                                      Prev1AdOrPreClin_age57     = dt_adenoma_prev[age_group == '[55,60)' , prev_1_adenoma],
                                      Prev1AdOrPreClin_age62     = dt_adenoma_prev[age_group == '[60,65)' , prev_1_adenoma],
                                      Prev1AdOrPreClin_age67     = dt_adenoma_prev[age_group == '[65,70)' , prev_1_adenoma],
                                      Prev1AdOrPreClin_age72     = dt_adenoma_prev[age_group == '[70,75)' , prev_1_adenoma],
                                      Prev1AdOrPreClin_age77     = dt_adenoma_prev[age_group == '[75,80)' , prev_1_adenoma],
                                      Prev1AdOrPreClin_age82     = dt_adenoma_prev[age_group == '[80,85)' , prev_1_adenoma],
                                      Prev1AdOrPreClin_age87     = dt_adenoma_prev[age_group == '[85,90)' , prev_1_adenoma],
                                      Prev1AdOrPreClin_age92     = dt_adenoma_prev[age_group == '[90,95)' , prev_1_adenoma],
                                      Prev1AdOrPreClin_age97     = dt_adenoma_prev[age_group == '[95,Inf)' , prev_1_adenoma],
                                      Prev2AdOrPreClin_age27     = dt_adenoma_prev[age_group == '[25,30)' , prev_2_adenoma],
                                      Prev2AdOrPreClin_age32     = dt_adenoma_prev[age_group == '[30,35)' , prev_2_adenoma],
                                      Prev2AdOrPreClin_age37     = dt_adenoma_prev[age_group == '[35,40)' , prev_2_adenoma],
                                      Prev2AdOrPreClin_age42     = dt_adenoma_prev[age_group == '[40,45)' , prev_2_adenoma],
                                      Prev2AdOrPreClin_age47     = dt_adenoma_prev[age_group == '[45,50)' , prev_2_adenoma],
                                      Prev2AdOrPreClin_age52     = dt_adenoma_prev[age_group == '[50,55)' , prev_2_adenoma],
                                      Prev2AdOrPreClin_age57     = dt_adenoma_prev[age_group == '[55,60)' , prev_2_adenoma],
                                      Prev2AdOrPreClin_age62     = dt_adenoma_prev[age_group == '[60,65)' , prev_2_adenoma],
                                      Prev2AdOrPreClin_age67     = dt_adenoma_prev[age_group == '[65,70)' , prev_2_adenoma],
                                      Prev2AdOrPreClin_age72     = dt_adenoma_prev[age_group == '[70,75)' , prev_2_adenoma],
                                      Prev2AdOrPreClin_age77     = dt_adenoma_prev[age_group == '[75,80)' , prev_2_adenoma],
                                      Prev2AdOrPreClin_age82     = dt_adenoma_prev[age_group == '[80,85)' , prev_2_adenoma],
                                      Prev2AdOrPreClin_age87     = dt_adenoma_prev[age_group == '[85,90)' , prev_2_adenoma],
                                      Prev2AdOrPreClin_age92     = dt_adenoma_prev[age_group == '[90,95)' , prev_2_adenoma],
                                      Prev2AdOrPreClin_age97     = dt_adenoma_prev[age_group == '[95,Inf)' , prev_2_adenoma],
                                      Prev3AdOrPreClin_age27     = dt_adenoma_prev[age_group == '[25,30)' , prev_3_adenoma],
                                      Prev3AdOrPreClin_age32     = dt_adenoma_prev[age_group == '[30,35)' , prev_3_adenoma],
                                      Prev3AdOrPreClin_age37     = dt_adenoma_prev[age_group == '[35,40)' , prev_3_adenoma],
                                      Prev3AdOrPreClin_age42     = dt_adenoma_prev[age_group == '[40,45)' , prev_3_adenoma],
                                      Prev3AdOrPreClin_age47     = dt_adenoma_prev[age_group == '[45,50)' , prev_3_adenoma],
                                      Prev3AdOrPreClin_age52     = dt_adenoma_prev[age_group == '[50,55)' , prev_3_adenoma],
                                      Prev3AdOrPreClin_age57     = dt_adenoma_prev[age_group == '[55,60)' , prev_3_adenoma],
                                      Prev3AdOrPreClin_age62     = dt_adenoma_prev[age_group == '[60,65)' , prev_3_adenoma],
                                      Prev3AdOrPreClin_age67     = dt_adenoma_prev[age_group == '[65,70)' , prev_3_adenoma],
                                      Prev3AdOrPreClin_age72     = dt_adenoma_prev[age_group == '[70,75)' , prev_3_adenoma],
                                      Prev3AdOrPreClin_age77     = dt_adenoma_prev[age_group == '[75,80)' , prev_3_adenoma],
                                      Prev3AdOrPreClin_age82     = dt_adenoma_prev[age_group == '[80,85)' , prev_3_adenoma],
                                      Prev3AdOrPreClin_age87     = dt_adenoma_prev[age_group == '[85,90)' , prev_3_adenoma],
                                      Prev3AdOrPreClin_age92     = dt_adenoma_prev[age_group == '[90,95)' , prev_3_adenoma],
                                      Prev3AdOrPreClin_age97     = dt_adenoma_prev[age_group == '[95,Inf)' , prev_3_adenoma],
                                      CRCincPer100K_P_ages20_39  = dt_crc_inc_P[age_group == '[20,40)' , inc_overall], 
                                      CRCincPer100K_P_ages40_49  = dt_crc_inc_P[age_group == '[40,50)' , inc_overall],
                                      CRCincPer100K_P_ages50_59  = dt_crc_inc_P[age_group == '[50,60)' , inc_overall],
                                      CRCincPer100K_P_ages60_69  = dt_crc_inc_P[age_group == '[60,70)' , inc_overall],
                                      CRCincPer100K_P_ages70_79  = dt_crc_inc_P[age_group == '[70,80)' , inc_overall],
                                      CRCincPer100K_P_ages80_99  = dt_crc_inc_P[age_group == '[80,Inf)' , inc_overall],
                                      CRCincPer100K_D_ages20_39  = dt_crc_inc_D[age_group == '[20,40)' , inc_overall],
                                      CRCincPer100K_D_ages40_49  = dt_crc_inc_D[age_group == '[40,50)' , inc_overall],
                                      CRCincPer100K_D_ages50_59  = dt_crc_inc_D[age_group == '[50,60)' , inc_overall],
                                      CRCincPer100K_D_ages60_69  = dt_crc_inc_D[age_group == '[60,70)' , inc_overall],
                                      CRCincPer100K_D_ages70_79  = dt_crc_inc_D[age_group == '[70,80)' , inc_overall],
                                      CRCincPer100K_D_ages80_99  = dt_crc_inc_D[age_group == '[80,Inf)' , inc_overall],
                                      CRCincPer100K_R_ages20_39  = dt_crc_inc_R[age_group == '[20,40)' , inc_overall],
                                      CRCincPer100K_R_ages40_49  = dt_crc_inc_R[age_group == '[40,50)' , inc_overall],
                                      CRCincPer100K_R_ages50_59  = dt_crc_inc_R[age_group == '[50,60)' , inc_overall],
                                      CRCincPer100K_R_ages60_69  = dt_crc_inc_R[age_group == '[60,70)' , inc_overall],
                                      CRCincPer100K_R_ages70_79  = dt_crc_inc_R[age_group == '[70,80)' , inc_overall],
                                      CRCincPer100K_R_ages80_99  = dt_crc_inc_R[age_group == '[80,Inf)' , inc_overall],
                                      CRCincPer100K_ages20_39   = dt_crc_inc_overall[age_group == '[20,40)' , inc_overall],
                                      CRCincPer100K_ages40_49   = dt_crc_inc_overall[age_group == '[40,50)' , inc_overall],
                                      CRCincPer100K_ages50_59   = dt_crc_inc_overall[age_group == '[50,60)' , inc_overall],
                                      CRCincPer100K_ages60_69   = dt_crc_inc_overall[age_group == '[60,70)' , inc_overall],
                                      CRCincPer100K_ages70_79   = dt_crc_inc_overall[age_group == '[70,80)' , inc_overall],
                                      CRCincPer100K_ages80_99   = dt_crc_inc_overall[age_group == '[80,Inf)' , inc_overall],
                                      StageDistribution_P_S1    = dt_stage_dist_P[ , stage_dist_S1],
                                      StageDistribution_P_S2    = dt_stage_dist_P[ , stage_dist_S2],
                                      StageDistribution_P_S3    = dt_stage_dist_P[ , stage_dist_S3],
                                      StageDistribution_P_S4    = dt_stage_dist_P[ , stage_dist_S4],
                                      StageDistribution_D_S1    = dt_stage_dist_D[ , stage_dist_S1],
                                      StageDistribution_D_S2    = dt_stage_dist_D[ , stage_dist_S2],
                                      StageDistribution_D_S3    = dt_stage_dist_D[ , stage_dist_S3],
                                      StageDistribution_D_S4    = dt_stage_dist_D[ , stage_dist_S4],
                                      StageDistribution_R_S1    = dt_stage_dist_R[ , stage_dist_S1],
                                      StageDistribution_R_S2    = dt_stage_dist_R[ , stage_dist_S2],
                                      StageDistribution_R_S3    = dt_stage_dist_R[ , stage_dist_S3],
                                      StageDistribution_R_S4    = dt_stage_dist_R[ , stage_dist_S4],
                                      PrevPreclinical_ages40_49 = dt_prev_preclinical[age_group == '[40,50)' , prev_precl],
                                      PrevPreclinical_ages50_59 = dt_prev_preclinical[age_group == '[50,60)' , prev_precl],
                                      PrevPreclinical_ages60_69 = dt_prev_preclinical[age_group == '[60,70)' , prev_precl],
                                      PrevPreclinical_ages70_79 = dt_prev_preclinical[age_group == '[70,80)' , prev_precl],
                                      PrevPreclinical_ages80_89 = dt_prev_preclinical[age_group == '[80,90)' , prev_precl],
                                      SizeLRGivenAdInP_wtd      = dt_prop_size_P[ , prop_S],
                                      SizeMRGivenAdInP_wtd      = dt_prop_size_P[ , prop_M],
                                      SizeHRGivenAdInP_wtd      = dt_prop_size_P[ , prop_L],
                                      SizeLRGivenAdInD_wtd      = dt_prop_size_D[ , prop_S],
                                      SizeMRGivenAdInD_wtd      = dt_prop_size_D[ , prop_M],
                                      SizeHRGivenAdInD_wtd      = dt_prop_size_D[ , prop_L],
                                      SizeLRGivenAdInR_wtd      = dt_prop_size_R[ , prop_S],
                                      SizeMRGivenAdInR_wtd      = dt_prop_size_R[ , prop_M],
                                      SizeHRGivenAdInR_wtd      = dt_prop_size_R[ , prop_L],
                                      prop_SSP                  = prop_SSP,
                                      Prev0SSPOrPreClin_age27     = dt_SSP_prev[age_group == '[25,30)'  , prev_0_adenoma],
                                      Prev0SSPOrPreClin_age32     = dt_SSP_prev[age_group == '[30,35)'  , prev_0_adenoma],
                                      Prev0SSPOrPreClin_age37     = dt_SSP_prev[age_group == '[35,40)'  , prev_0_adenoma],
                                      Prev0SSPOrPreClin_age42     = dt_SSP_prev[age_group == '[40,45)'  , prev_0_adenoma],
                                      Prev0SSPOrPreClin_age47     = dt_SSP_prev[age_group == '[45,50)'  , prev_0_adenoma],
                                      Prev0SSPOrPreClin_age52     = dt_SSP_prev[age_group == '[50,55)'  , prev_0_adenoma],
                                      Prev0SSPOrPreClin_age57     = dt_SSP_prev[age_group == '[55,60)'  , prev_0_adenoma],
                                      Prev0SSPOrPreClin_age62     = dt_SSP_prev[age_group == '[60,65)'  , prev_0_adenoma],
                                      Prev0SSPOrPreClin_age67     = dt_SSP_prev[age_group == '[65,70)'  , prev_0_adenoma],
                                      Prev0SSPOrPreClin_age72     = dt_SSP_prev[age_group == '[70,75)'  , prev_0_adenoma],
                                      Prev0SSPOrPreClin_age77     = dt_SSP_prev[age_group == '[75,80)'  , prev_0_adenoma],
                                      Prev0SSPOrPreClin_age82     = dt_SSP_prev[age_group == '[80,85)'  , prev_0_adenoma],
                                      Prev0SSPOrPreClin_age87     = dt_SSP_prev[age_group == '[85,90)'  , prev_0_adenoma],
                                      Prev0SSPOrPreClin_age92     = dt_SSP_prev[age_group == '[90,95)'  , prev_0_adenoma],
                                      Prev0SSPOrPreClin_age97     = dt_SSP_prev[age_group == '[95,Inf)' , prev_0_adenoma],
                                      Prev1SSPOrPreClin_age27     = dt_SSP_prev[age_group == '[25,30)'  , prev_1_adenoma],
                                      Prev1SSPOrPreClin_age32     = dt_SSP_prev[age_group == '[30,35)'  , prev_1_adenoma],
                                      Prev1SSPOrPreClin_age37     = dt_SSP_prev[age_group == '[35,40)'  , prev_1_adenoma],
                                      Prev1SSPOrPreClin_age42     = dt_SSP_prev[age_group == '[40,45)'  , prev_1_adenoma],
                                      Prev1SSPOrPreClin_age47     = dt_SSP_prev[age_group == '[45,50)'  , prev_1_adenoma],
                                      Prev1SSPOrPreClin_age52     = dt_SSP_prev[age_group == '[50,55)'  , prev_1_adenoma],
                                      Prev1SSPOrPreClin_age57     = dt_SSP_prev[age_group == '[55,60)'  , prev_1_adenoma],
                                      Prev1SSPOrPreClin_age62     = dt_SSP_prev[age_group == '[60,65)'  , prev_1_adenoma],
                                      Prev1SSPOrPreClin_age67     = dt_SSP_prev[age_group == '[65,70)'  , prev_1_adenoma],
                                      Prev1SSPOrPreClin_age72     = dt_SSP_prev[age_group == '[70,75)'  , prev_1_adenoma],
                                      Prev1SSPOrPreClin_age77     = dt_SSP_prev[age_group == '[75,80)'  , prev_1_adenoma],
                                      Prev1SSPOrPreClin_age82     = dt_SSP_prev[age_group == '[80,85)'  , prev_1_adenoma],
                                      Prev1SSPOrPreClin_age87     = dt_SSP_prev[age_group == '[85,90)'  , prev_1_adenoma],
                                      Prev1SSPOrPreClin_age92     = dt_SSP_prev[age_group == '[90,95)'  , prev_1_adenoma],
                                      Prev1SSPOrPreClin_age97     = dt_SSP_prev[age_group == '[95,Inf)' , prev_1_adenoma],
                                      Prev2SSPOrPreClin_age27     = dt_SSP_prev[age_group == '[25,30)'  , prev_2_adenoma],
                                      Prev2SSPOrPreClin_age32     = dt_SSP_prev[age_group == '[30,35)'  , prev_2_adenoma],
                                      Prev2SSPOrPreClin_age37     = dt_SSP_prev[age_group == '[35,40)'  , prev_2_adenoma],
                                      Prev2SSPOrPreClin_age42     = dt_SSP_prev[age_group == '[40,45)'  , prev_2_adenoma],
                                      Prev2SSPOrPreClin_age47     = dt_SSP_prev[age_group == '[45,50)'  , prev_2_adenoma],
                                      Prev2SSPOrPreClin_age52     = dt_SSP_prev[age_group == '[50,55)'  , prev_2_adenoma],
                                      Prev2SSPOrPreClin_age57     = dt_SSP_prev[age_group == '[55,60)'  , prev_2_adenoma],
                                      Prev2SSPOrPreClin_age62     = dt_SSP_prev[age_group == '[60,65)'  , prev_2_adenoma],
                                      Prev2SSPOrPreClin_age67     = dt_SSP_prev[age_group == '[65,70)'  , prev_2_adenoma],
                                      Prev2SSPOrPreClin_age72     = dt_SSP_prev[age_group == '[70,75)'  , prev_2_adenoma],
                                      Prev2SSPOrPreClin_age77     = dt_SSP_prev[age_group == '[75,80)'  , prev_2_adenoma],
                                      Prev2SSPOrPreClin_age82     = dt_SSP_prev[age_group == '[80,85)'  , prev_2_adenoma],
                                      Prev2SSPOrPreClin_age87     = dt_SSP_prev[age_group == '[85,90)'  , prev_2_adenoma],
                                      Prev2SSPOrPreClin_age92     = dt_SSP_prev[age_group == '[90,95)'  , prev_2_adenoma],
                                      Prev2SSPOrPreClin_age97     = dt_SSP_prev[age_group == '[95,Inf)' , prev_2_adenoma],
                                      Prev3SSPOrPreClin_age27     = dt_SSP_prev[age_group == '[25,30)'  , prev_3_adenoma],
                                      Prev3SSPOrPreClin_age32     = dt_SSP_prev[age_group == '[30,35)'  , prev_3_adenoma],
                                      Prev3SSPOrPreClin_age37     = dt_SSP_prev[age_group == '[35,40)'  , prev_3_adenoma],
                                      Prev3SSPOrPreClin_age42     = dt_SSP_prev[age_group == '[40,45)'  , prev_3_adenoma],
                                      Prev3SSPOrPreClin_age47     = dt_SSP_prev[age_group == '[45,50)'  , prev_3_adenoma],
                                      Prev3SSPOrPreClin_age52     = dt_SSP_prev[age_group == '[50,55)'  , prev_3_adenoma],
                                      Prev3SSPOrPreClin_age57     = dt_SSP_prev[age_group == '[55,60)'  , prev_3_adenoma],
                                      Prev3SSPOrPreClin_age62     = dt_SSP_prev[age_group == '[60,65)'  , prev_3_adenoma],
                                      Prev3SSPOrPreClin_age67     = dt_SSP_prev[age_group == '[65,70)'  , prev_3_adenoma],
                                      Prev3SSPOrPreClin_age72     = dt_SSP_prev[age_group == '[70,75)'  , prev_3_adenoma],
                                      Prev3SSPOrPreClin_age77     = dt_SSP_prev[age_group == '[75,80)'  , prev_3_adenoma],
                                      Prev3SSPOrPreClin_age82     = dt_SSP_prev[age_group == '[80,85)'  , prev_3_adenoma],
                                      Prev3SSPOrPreClin_age87     = dt_SSP_prev[age_group == '[85,90)'  , prev_3_adenoma],
                                      Prev3SSPOrPreClin_age92     = dt_SSP_prev[age_group == '[90,95)'  , prev_3_adenoma],
                                      Prev3SSPOrPreClin_age97     = dt_SSP_prev[age_group == '[95,Inf)' , prev_3_adenoma],
                                      PrevPreclinical_SSP_ages40_49 = dt_prev_preclinical_SSP[age_group == '[40,50)' , prev_precl],
                                      PrevPreclinical_SSP_ages50_59 = dt_prev_preclinical_SSP[age_group == '[50,60)' , prev_precl],
                                      PrevPreclinical_SSP_ages60_69 = dt_prev_preclinical_SSP[age_group == '[60,70)' , prev_precl],
                                      PrevPreclinical_SSP_ages70_79 = dt_prev_preclinical_SSP[age_group == '[70,80)' , prev_precl],
                                      PrevPreclinical_SSP_ages80_89 = dt_prev_preclinical_SSP[age_group == '[80,90)' , prev_precl],
                                      SizeLRGivenSSPInP_wtd      = dt_prop_size_P_SSP[ , prop_S],
                                      SizeMRGivenSSPInP_wtd      = dt_prop_size_P_SSP[ , prop_M],
                                      SizeHRGivenSSPInP_wtd      = dt_prop_size_P_SSP[ , prop_L],
                                      SizeLRGivenSSPInD_wtd      = dt_prop_size_D_SSP[ , prop_S],
                                      SizeMRGivenSSPInD_wtd      = dt_prop_size_D_SSP[ , prop_M],
                                      SizeHRGivenSSPInD_wtd      = dt_prop_size_D_SSP[ , prop_L],
                                      SizeLRGivenSSPInR_wtd      = dt_prop_size_R_SSP[ , prop_S],
                                      SizeMRGivenSSPInR_wtd      = dt_prop_size_R_SSP[ , prop_M],
                                      SizeHRGivenSSPInR_wtd      = dt_prop_size_R_SSP[ , prop_L],
                                      Adenoma_dwell             = Adenoma_dwell,
                                      Sojourn_time              = Sojourn_time, 
                                      Total_dwell               = Total_dwell
                                      
  )
  
  
  if (dt_long) { 
    
    dt <- t(as.matrix(dt_output, ))
    dt_output <- data.table( output_name = row.names(dt), output_value = dt[,1])
  }
  
  return( dt_output)
  
}


#' Sample from uniform distributions of calibrated parameters from a Latin-
#' hypercube sampling (LHS) design.
#'
#' \code{sample_uniform_lhs_ssp} generates a sample of parameter sets from uniform
#' distributions following a Latin hypercube sampling (LHS) inlcuding parameter for SSP.
#' @param n_samp Number of samples.
#' @param v_param_names Vector with parameter names.
#' @param v_lb Vector with lower bounds for each parameter.
#' @param v_ub Vector with upper bounds for each parameter.
#' @return 
#' A matrix with as many columns as calibrated parameters and \code{n_samp} rows
#' . Each row corresponds to a parameter set sampled from uniform distributions
#' using a uniform LHS design.
#' @export
sample_uniform_lhs_ssp <- function(n_samp,
                               v_param_names = c(
                                 "IndividuakRiskMultVariance", 
                                 "alpha_lesion_adenoma",
                                 "beta_age",
                                 "AdNaturalHistoryPropensity_Gaussian_Variance",
                                 "AdGrowth_Exp_Rate_1_to_6_P",
                                 "AdGrowth_Exp_Rate_1_to_6_D",
                                 "AdGrowth_Exp_Rate_1_to_6_R",
                                 "AdGrowth_Exp_Rate_6_to_10_P",
                                 "AdGrowth_Exp_Rate_6_to_10_D",
                                 "AdGrowth_Exp_Rate_6_to_10_R", 
                                 "PreclinCancerProg_Exp_Rate_S1S2_P",
                                 "PreclinCancerProg_Exp_Rate_S1S2_D",
                                 "PreclinCancerProg_Exp_Rate_S1S2_R",
                                 "PreclinCancerProg_Exp_Rate_S2S3_P",
                                 "PreclinCancerProg_Exp_Rate_S2S3_D",
                                 "PreclinCancerProg_Exp_Rate_S2S3_R",
                                 "PreclinCancerProg_Exp_Rate_S3S4_P",
                                 "PreclinCancerProg_Exp_Rate_S3S4_D",
                                 "PreclinCancerProg_Exp_Rate_S3S4_R",
                                 "CancerOnset_Gompertz_Shape_P",
                                 "CancerOnset_Gompertz_Shape_D",
                                 "CancerOnset_Gompertz_Shape_R",
                                 "CancerOnset_Gompertz_Rate_P",
                                 "CancerOnset_Gompertz_Rate_D",
                                 "CancerOnset_Gompertz_Rate_R",
                                 "pSxDetS1_P",
                                 "pSxDetS1_D",
                                 "pSxDetS1_R",
                                 "hr_SxDetS2S1_P",
                                 "hr_SxDetS2S1_D",
                                 "hr_SxDetS2S1_R",
                                 "hr_SxDetS3S2_P",
                                 "hr_SxDetS3S2_D",
                                 "hr_SxDetS3S2_R",
                                 "hr_SxDetS4S3_P",
                                 "hr_SxDetS4S3_D",
                                 "hr_SxDetS4S3_R",           
                                 "alpha_lesion_ssp", #SSP 
                                 "beta_age_ssp",  # SSP
                                 "hr_SSPGrowth_1_to_6_P",  ##SSP
                                 "hr_SSPGrowth_6_to_10_P",  ##SSP
                                 "hr_SSPGrowth_1_to_6_D",  ##SSP
                                 "hr_SSPGrowth_6_to_10_D",  ##SSP
                                 "hr_SSPGrowth_1_to_6_R",  ##SSP
                                 "hr_SSPGrowth_6_to_10_R",  ##SSP
                                 "hr_SSPCancerOnset_P"  ##SSP
                               ),
                               # Lower bounds. Uniform distributions are defined.
                               #how to define the lower bounds
                               #Outputs of the model cover
                               v_lb = c(
                                 IndividuakRiskMultVariance = 0.01, 
                                 alpha_lesion_adenoma = -8, 
                                 beta_age =  0.02, 
                                 AdNaturalHistoryPropensity_Gaussian_Variance = 0.30,
                                 AdGrowth_Exp_Rate_1_to_6_P = 0.001, 
                                 AdGrowth_Exp_Rate_1_to_6_D = 0.001,
                                 AdGrowth_Exp_Rate_1_to_6_R = 0.001, 
                                 AdGrowth_Exp_Rate_6_to_10_P = 0.001,
                                 AdGrowth_Exp_Rate_6_to_10_D = 0.001,
                                 AdGrowth_Exp_Rate_6_to_10_R = 0.001, 
                                 PreclinCancerProg_Exp_Rate_S1S2_P = 0.20, #0.228536835,
                                 PreclinCancerProg_Exp_Rate_S1S2_D = 0.20, #0.270014625,
                                 PreclinCancerProg_Exp_Rate_S1S2_R = 0.30, #0.260442195,
                                 PreclinCancerProg_Exp_Rate_S2S3_P = 0.20, #0.0000001,
                                 PreclinCancerProg_Exp_Rate_S2S3_D = 0.20, #0.0000001,
                                 PreclinCancerProg_Exp_Rate_S2S3_R = 0.20, #0.0000001,
                                 PreclinCancerProg_Exp_Rate_S3S4_P = 0.30, #0.0000001,
                                 PreclinCancerProg_Exp_Rate_S3S4_D = 0.30, #0.0000001,
                                 PreclinCancerProg_Exp_Rate_S3S4_R = 0.30, #0.0000001,
                                 CancerOnset_Gompertz_Shape_P = 0.00001,
                                 CancerOnset_Gompertz_Shape_D = 0.0001,
                                 CancerOnset_Gompertz_Shape_R = 0.0001,
                                 CancerOnset_Gompertz_Rate_P = 0.001,
                                 CancerOnset_Gompertz_Rate_D = 0.001,
                                 CancerOnset_Gompertz_Rate_R = 0.001,
                                 pSxDetS1_P = 0.001,
                                 pSxDetS1_D = 0.001,
                                 pSxDetS1_R = 0.001,
                                 hr_SxDetS2S1_P = 1,
                                 hr_SxDetS2S1_D = 1,
                                 hr_SxDetS2S1_R = 1,
                                 hr_SxDetS3S2_P = 1,
                                 hr_SxDetS3S2_D = 1,
                                 hr_SxDetS3S2_R = 1,
                                 hr_SxDetS4S3_P = 1,
                                 hr_SxDetS4S3_D = 1,
                                 hr_SxDetS4S3_R = 1,
                                 alpha_lesion_ssp = -10, #SSP 
                                 beta_age_ssp  = -0.02,  # SSP 
                                 hr_SSPGrowth_1_to_6_P = 0.01,  ##SSP
                                 hr_SSPGrowth_6_to_10_P = 0.01,  ##SSP
                                 hr_SSPGrowth_1_to_6_D = 0.01,  ##SSP
                                 hr_SSPGrowth_6_to_10_D = 0.01,  ##SSP
                                 hr_SSPGrowth_1_to_6_R = 0.01,  ##SSP
                                 hr_SSPGrowth_6_to_10_R = 0.01,  ##SSP
                                 hr_SSPCancerOnset_P = 0.8  ##SSP
                               ),
                               v_ub = c(
                                 IndividuakRiskMultVariance = 0.8, 
                                 alpha_lesion_adenoma = -5, 
                                 beta_age =  0.06, 
                                 AdNaturalHistoryPropensity_Gaussian_Variance = 0.5,
                                 AdGrowth_Exp_Rate_1_to_6_P = 0.1, 
                                 AdGrowth_Exp_Rate_1_to_6_D = 0.2,
                                 AdGrowth_Exp_Rate_1_to_6_R = 0.3, 
                                 AdGrowth_Exp_Rate_6_to_10_P = 0.06,
                                 AdGrowth_Exp_Rate_6_to_10_D = 0.2,
                                 AdGrowth_Exp_Rate_6_to_10_R = 0.45,
                                 PreclinCancerProg_Exp_Rate_S1S2_P = 0.35, 
                                 PreclinCancerProg_Exp_Rate_S1S2_D = 0.35, 
                                 PreclinCancerProg_Exp_Rate_S1S2_R = 0.40, 
                                 PreclinCancerProg_Exp_Rate_S2S3_P = 0.30, 
                                 PreclinCancerProg_Exp_Rate_S2S3_D = 0.40, 
                                 PreclinCancerProg_Exp_Rate_S2S3_R = 0.50, 
                                 PreclinCancerProg_Exp_Rate_S3S4_P = 0.52, 
                                 PreclinCancerProg_Exp_Rate_S3S4_D = 0.7, 
                                 PreclinCancerProg_Exp_Rate_S3S4_R = 0.9, 
                                 CancerOnset_Gompertz_Shape_P = 0.0001,   #change for approach with exponential
                                 CancerOnset_Gompertz_Shape_D = 0.02,  #change for approach with exponential
                                 CancerOnset_Gompertz_Shape_R = 0.025,  #change for approach with exponential
                                 CancerOnset_Gompertz_Rate_P = 0.1, #change for approach with exponential
                                 CancerOnset_Gompertz_Rate_D = 0.04, #0.1 change for approach with exponential
                                 CancerOnset_Gompertz_Rate_R = 0.2, #change for approach with exponential
                                 pSxDetS1_P = 0.192,
                                 pSxDetS1_D = 0.192,
                                 pSxDetS1_R = 0.192,
                                 hr_SxDetS2S1_P = 24,
                                 hr_SxDetS2S1_D = 24,
                                 hr_SxDetS2S1_R = 24,
                                 hr_SxDetS3S2_P = 7,
                                 hr_SxDetS3S2_D = 7,
                                 hr_SxDetS3S2_R = 7,
                                 hr_SxDetS4S3_P = 7,
                                 hr_SxDetS4S3_D = 7,
                                 hr_SxDetS4S3_R = 7,
                                 alpha_lesion_ssp = -5, #SSP 
                                 beta_age_ssp  = 0.06,  # SSP 
                                 hr_SSPGrowth_1_to_6_P = 2.5,  ##SSP
                                 hr_SSPGrowth_6_to_10_P = 2.5,  ##SSP
                                 hr_SSPGrowth_1_to_6_D = 2.5,  ##SSP
                                 hr_SSPGrowth_6_to_10_D = 2.5,  ##SSP
                                 hr_SSPGrowth_1_to_6_R = 2.5,  ##SSP
                                 hr_SSPGrowth_6_to_10_R = 2.5,  ##SSP
                                 hr_SSPCancerOnset_P = 2.5  ##SSP
                               )){ 
  # Obtain number of parameters
  n_param <- length(v_param_names)
  
  ### Draw LHS from Uniform[0,1] distributions
  m_lhs_unit   <- lhs::randomLHS(n = n_samp, k = n_param)
  m_param_samp <- matrix(nrow = n_samp, ncol = n_param)
  colnames(m_param_samp) <- v_param_names
  for (i in 1:n_param){
    m_param_samp[, i] <- qunif(m_lhs_unit[, i],
                               min = v_lb[i],
                               max = v_ub[i])
  }
  
  l_results <- list(lower_bounds=v_lb , upper_bounds = v_ub, m_param_samp = m_param_samp)
  
  return(l_results)
  
}



