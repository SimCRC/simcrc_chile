#' Generate model outputs for calibration from a parameter set
#'
#' \code{calibration_out} computes model outputs to be used for calibration 
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
calibration_out <- function(v_params_calib, l_params_all, dt_pop, id_draw = i, dt_long = FALSE){ # User defined
  # Substitute values of calibrated parameters in base-case with 
  # calibrated values
  l_params_all <- update_param_list(l_params_all   = l_params_all, 
                                    params_updated = v_params_calib)
  
  # Run simcrc nathist
  # l_out_simcrc <- simcrc_nathist_DES(l_params_all = l_params_all,
  #                              dt_pop       = dt_pop)
  
  # Run simcrc nathist WITH SSP Pathway
  l_out_simcrc <- simcrc::simcr_nathist_ssp_DES(l_params_all = l_params_all,
                                                dt_pop       = dt_pop,
                                                SSP_pathway = FALSE)
  
  
  # Get main outputs of simcrc_nathist
  dt_crc_pop <- l_out_simcrc$dt_crc_pop
  dt_crc_adenomas_i <- l_out_simcrc$dt_crc_adenomas
  
  # Get living population data by age, sex, smoking and year
  dt_living_pop <- get_living_pop(dt_crc_pop= dt_crc_pop,
                                  min_age = 0, 
                                  max_age = 100)
  
  
  # Calculate adenoma prevalence
  dt_adenoma_prev <- simcrc::get_adenoma_prevalence_ssp(
    dt_crc_pop = dt_crc_pop,
    dt_living_pop = dt_living_pop,
    v_lesion_location = c("P", "D", "R"),
    v_age = c(25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99),
    lesion_pathway = "Adenoma",
    age_groups = TRUE
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
  
  # Calculate preclinical prevalence
  
  dt_prev_preclinical <- simcrc::get_prevalence_preclinical_ssp (dt_crc_pop= dt_crc_pop, 
                                                                 dt_living_pop = dt_living_pop,
                                                                 v_lesion_location = c("P", "D", "R"),
                                                                 v_age = c(40,50,60,70,80,90),
                                                                 lesion_pathway = "Adenoma",
                                                                 age_groups  = TRUE)
  
  # Calculate size distribution 
  
  dt_prop_size_P <- simcrc::get_size_distribution_ssp (dt_crc_pop = dt_crc_pop, 
                                                       dt_living_pop = dt_living_pop,
                                                       ad_location = c("P"),
                                                       v_age = c(50:75),
                                                       lesion_pathway = "Adenoma",
                                                       age_groups  = FALSE)
  
  dt_prop_size_D <- simcrc::get_size_distribution_ssp (dt_crc_pop = dt_crc_pop, 
                                                       dt_living_pop = dt_living_pop,
                                                       ad_location = c("D"),
                                                       v_age = c(50:75),
                                                       lesion_pathway = "Adenoma",
                                                       age_groups  = FALSE)
  
  dt_prop_size_R <- simcrc::get_size_distribution_ssp (dt_crc_pop = dt_crc_pop, 
                                                       dt_living_pop = dt_living_pop,
                                                       ad_location = c("R"),
                                                       v_age = c(50:75),
                                                       lesion_pathway = "Adenoma",
                                                       age_groups  = FALSE)
  
  
  # Calculates adenoma dwell time
  
  Adenoma_dwell <- mean(dt_crc_pop$t_Adonset_Preclin, na.rm = TRUE)
  
  Sojourn_time  <- mean(dt_crc_pop$t_Preclin_Sxdet, na.rm = TRUE)
  
  Total_dwell   <- mean(dt_crc_pop$t_Adonset_Sxdet, na.rm = TRUE)
  
  
  
  dt_output <- data.table::data.table(id_draw = id_draw,
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
                                      Prev0AdOrPreClin_age97     = dt_adenoma_prev[age_group == '[95,99)' , prev_0_adenoma],
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
                                      Prev1AdOrPreClin_age97     = dt_adenoma_prev[age_group == '[95,99)' , prev_1_adenoma],
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
                                      Prev2AdOrPreClin_age97     = dt_adenoma_prev[age_group == '[95,99)' , prev_2_adenoma],
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
                                      Prev3AdOrPreClin_age97     = dt_adenoma_prev[age_group == '[95,99)' , prev_3_adenoma],
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


#' Generate dweel time from a set of input parameters
#'
#' \code{calibration_out_dwell} computes dwell time of adenomas 
#'
#' @param v_params_calib Vector of parameters that need to be calibrated.
#' @param l_params_all List with all parameters of the decision model.
#' @param dt_pop Datatable dataset with the cohort information.
#' @param dt_long To select the structure of datatable, default is wide.
#' @return 
#' 
#' A datatable with the dwell times
#' 
#' @export
calibration_out_dwell <- function(v_params_calib, l_params_all, dt_pop, dt_long = FALSE){ # User defined
  # Substitute values of calibrated parameters in base-case with 
  # calibrated values
  l_params_all <- update_param_list(l_params_all   = l_params_all, 
                                    params_updated = v_params_calib)
  
  # Run simcrc nathist
  l_out_simcrc <- simcrc_nathist_DES(l_params_all = l_params_all, 
                                     dt_pop       = dt_pop)
  
  
  # Get main outputs of simcrc_nathist
  dt_crc_pop <- l_out_simcrc$dt_crc_pop
  dt_crc_adenomas_i <- l_out_simcrc$dt_crc_adenomas
  
  dt_dwell <- simcrc::get_adenoma_dwell_times(dt_crc_pop = dt_crc_pop)
  
  if (dt_long) { 
    
    dt <- t(as.matrix(dt_dwell, ))
    dt_dwell <- data.table( output_name = row.names(dt), output_value = dt[,1])
  }
  
  return(dt_dwell)
  
}


#' Sample from uniform distributions of calibrated parameters from a Latin-
#' hypercube sampling (LHS) design.
#'
#' \code{sample_uniform_lhs} generates a sample of parameter sets from uniform
#' distributions following a Latin hypercube sampling (LHS) design.
#' @param n_samp Number of samples.
#' @param v_param_names Vector with parameter names.
#' @param v_lb Vector with lower bounds for each parameter.
#' @param v_ub Vector with upper bounds for each parameter.
#' @return 
#' A matrix with as many columns as calibrated parameters and \code{n_samp} rows
#' . Each row corresponds to a parameter set sampled from uniform distributions
#' using a uniform LHS design.
#' @export

sample_uniform_lhs <- function(n_samp, l_params_priors
                               ){ 
  # Obtain number of parameters
  n_param <- length(l_params_priors$names)
  
  ### Draw LHS from Uniform[0,1] distributions
  m_lhs_unit   <- lhs::randomLHS(n = n_samp, k = n_param)
  m_param_samp <- matrix(nrow = n_samp, ncol = n_param)
  colnames(m_param_samp) <- l_params_priors$names
  
  #Check that upper and lower bounds are provided for all parameters and are numeric
  if (length(l_params_priors$lb) != n_param | length(l_params_priors$ub) != n_param){
    stop("Error in parameter bounds: number of lower and upper bounds must match number of parameters")
  }
  if (!is.numeric(l_params_priors$lb) | !is.numeric(l_params_priors$ub)){
    stop("Error in parameter bounds: lower and upper bounds must be numeric")
  }
  
  for (i in 1:n_param){
    m_param_samp[, i] <- qunif(m_lhs_unit[, i],
                               min = l_params_priors$lb[i],
                               max = l_params_priors$ub[i])
  
    #Confirm that upper bound is greater than lower bound
    if (l_params_priors$ub[i] <= l_params_priors$lb[i]){
      stop(paste0("Error in sampling parameter ", l_params_priors$names[i], ": upper bound not greater than lower bound"))
    }
    
  }
  
  #Test to avoid wrong bounds
  if (any(m_param_samp[, i] < l_params_priors$lb[i]) | any(m_param_samp[, i] > l_params_priors$ub[i])){
    stop(paste0("Error in sampling parameter ", l_params_priors$names[i], ": sampled values outside bounds"))
  }
  
  l_results <- list(lower_bounds=l_params_priors$lb , upper_bounds = l_params_priors$ub, m_param_samp = m_param_samp)
  
  return(l_results)
}

#test <- sample_uniform_lhs(100000,l_params_priors_Adenoma_female)


#' Get operating system
#' 
#' @return 
#' A string with the operating system.
#' @export
get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "MacOSX"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}
