# *****************************************************************************
#
# Script: uspstf_summary.R
#
# Purpose: Template that summarizes key outputs from the SimCRC model
#
# Author: Claudia Seguin and Jorge Roa
#
# Email: cseguin@stanford.edu/jorgeroa@stanford.edu
#
# Date Created: 01-Aug-2025
# *****************************************************************************

## USPSTF summary ------------------------------------------------------
#'
#'
#' \code{uspstf_summary} USPSTF summary template
#'
#' @param datatable Data table at the individual level with age of death, age of diagnosis, stage at diagnosis and cause of death.
#' @param screening_modality Screening modality simulated.
#' @param follow_up Boolean if the screening modality has a follow-up colonoscopy or not.
#' @param screening_years Years in which screening happens.
#' @param min_age The cohort age at which to start summarizing outcomes.
#' @param max_age The maximum age we care about.
#' @param output_template_year the year of the template we want to output. Options are 2020 and 2028. The 2028 template has more columns than the 2020 template, related to the new tests and SSLs
#' @return 
#' Matrix with the allocation of years by phases \code{m_allocation_time}.
#' @export
#' 
uspstf_summary <- function(datatable = dt_pop_surv, 
                           screening_modality = "COL",
                           follow_up = FALSE,
                           screening_years = c(50),
                           min_age = 40, 
                           max_age = 100,
                           output_template_year = 2028 # options are 2020 and 2028
                           ) {
  
  # Error in case the you pass an unexpected variable
  if(!(screening_modality %in% c("NoScreening", "COL", "FIT", "SIG", "HSgFOBT", "mtsDNA","SRNA", "BLOOD", "CTC")) ){
    stop("The screening_modality variable must be one of the following: 'NoScreening', 'COL', 'FIT', 'SIG', 'HSgFOBT', 'mtsDNA', 'SRNA', 'BLOOD', 'CTC'. The function uspstf_summary.R cannot yet handle other screening modalities.")
  }
  if(!(follow_up %in% c(TRUE, FALSE)) ){
    stop("The follow_up variable must be one of the following: TRUE or FALSE.")
  }
  if(!is.numeric(screening_years)) {
    stop("The screening_years variable must be a numeric vector.")
  }
  
  # *****************************************************************************
  #### 0.0 Create skeleton for the template called dt_export ####
  # *****************************************************************************
  
  dt_export <-  data.table(age = seq(min_age, max_age))
  
  # *****************************************************************************
  #### 1.0 Subset the correct population for dt_export ####
  # *****************************************************************************

  dt_pop_surv <- copy(datatable)  
  
  # Extra variables needed for the crc_allocation_time function
  dt_pop_surv$age_dx <- dt_pop_surv$age_CRC_dx
  dt_pop_surv$stage_at_dx <- dt_pop_surv$S_stage_at_dx
  
  # Keep only one row per individual since all summarized counts are individual-level
  dt_pop_surv_id <- dt_pop_surv[, .SD[1], by = "id"]
  
  # From all outcomes in the USPSTF template, remove people who die or get diagnosed cancer before age 40  
  dt_pop_surv_id <- dt_pop_surv_id[(age_CRC_dx >= min_age | is.na(age_CRC_dx)) & (age_death >= min_age | is.na(age_death)), ]

  
  # *****************************************************************************
  ##### 1.1 Create indicators and integers for certain variables ####
  # *****************************************************************************

  dt_pop_surv_id[, age_death_int := floor(age_death)]
  dt_pop_surv_id[, age_ind_crc:= pmin(age_CRC_dx, age_death, na.rm = TRUE)] 
  dt_pop_surv_id[, age_ind_crc_int := floor(age_ind_crc)] 
  
  # Check that every simulate person has an age of death
  if(any(is.na(dt_pop_surv_id$age_death))){
    stop("Every simulated individual must have an age of death.")
  }
  
  # *****************************************************************************
  #### 2.0  Count Variables for USPSTF ####
  # *****************************************************************************
  
  # *****************************************************************************
  ##### 2.1 Deaths by cause ####
  # *****************************************************************************
 
  # Number of deaths from crc
  dt_deaths_crc <- dt_pop_surv_id[death_cause == "crc" , .(n_CRC_dth = .N), by = age_death_int]
  dt_export[dt_deaths_crc, on = .(age = age_death_int), n_CRC_dth := i.n_CRC_dth]
  
  # Number of deaths from other causes 
  dt_deaths_oc <- dt_pop_surv_id[death_cause == "oc" , .(n_oc_dth = .N), by = age_death_int]
  dt_export[dt_deaths_oc, on = .(age = age_death_int), n_oc_dth := i.n_oc_dth]
  
  # Number of deaths from colonoscopies
  dt_deaths_col <- dt_pop_surv_id[death_cause %in% c("scr","surv") , .(n_compl_dth = .N), by = age_death_int]
  dt_export[dt_deaths_col, on = .(age = age_death_int), n_compl_dth := i.n_compl_dth]
  
  # *****************************************************************************
  ##### 2.2 People alive ####
  # *****************************************************************************
  
  # To calculate the number of people alive
  # We want to keep the double counting of the number of people alive and the number of people who died within a given age (aka the shift) to match the original definitions in the template that is used by all cisnet models.
  # For example, we assume if a person dies at age 40.5, we want them to count as both alive at the beginning of age 45, and dead during age 45.
  # If we remove the shift, there will be errors in the process_uspstf function and the implementation of screening will be different across models.
  
  ## Total number of deaths
  dt_deaths <- dt_pop_surv_id[, .(n_deaths = .N), by = age_death_int]
  dt_export[dt_deaths, on = .(age = age_death_int), `:=`(n_deaths = fifelse(is.na(i.n_deaths), 0, i.n_deaths))]
  dt_export[, n_cumsum_death := cumsum(n_deaths)]
  
  ## n_alive
  dt_export[, `:=` (n_alive = nrow(dt_pop_surv_id) - shift(n_cumsum_death, n = 1, fill = 0))]
  
  # *****************************************************************************
  ##### 2.3 CRC cases ####
  # *****************************************************************************

  # Use the simcrc summary function to calculate the number of CRC cases
  dt_CRC_by_stage_and_detection <- calc_counts_Dx(dt_crc_pop = dt_pop_surv, 
                                                  min_age = min_age, 
                                                  max_age = max_age,
                                                  ad_location = c("P", "D", "R"),
                                                  dt_long = FALSE,
                                                  extended = TRUE )
  
  # Only keep the relevant columns
  dt_crc <- dt_CRC_by_stage_and_detection[, .(age, 
                                              S1_scr,  S2_scr,  S3_scr,  S4_scr,
                                              S1_symp, S2_symp, S3_symp, S4_symp,
                                              S1_surv, S2_surv, S3_surv, S4_surv)]
  
  # Warning if you are simulating both sexes at the same time
  if(nrow(dt_crc[age == min_age,]) > 1){
    # If there are multiple sexes simulated, this should work to collapse
    dt_crc <- dt_crc[, lapply(.SD, sum), by = age, .SDcols = c("S1_scr",  "S2_scr",  "S3_scr",  "S4_scr",
                                                                                                              "S1_symp", "S2_symp", "S3_symp", "S4_symp",
                                                                                                              "S1_surv", "S2_surv", "S3_surv", "S4_surv")]
    warning("This function (calc_counts_Dx) is not collapsing over sex and will causes issues with merges to dt_export later in the code." )
  }
  
  
  # Rename S1_scr to n_CRCstageI_scr, S2_scr to n_CRCstageII_scr, etc. to match naming conventions in the template
  setnames(dt_crc, old = c("S1_scr", "S2_scr", "S3_scr", "S4_scr",
                           "S1_symp", "S2_symp", "S3_symp", "S4_symp",
                           "S1_surv", "S2_surv", "S3_surv", "S4_surv"),
           new = c("n_CRCstageI_scr", "n_CRCstageII_scr", "n_CRCstageIII_scr", "n_CRCstageIV_scr",
                   "n_CRCstageI_sym", "n_CRCstageII_sym", "n_CRCstageIII_sym", "n_CRCstageIV_sym",
                   "n_CRCstageI_surv", "n_CRCstageII_surv", "n_CRCstageIII_surv", "n_CRCstageIV_surv"))
  
  dt_export <- merge(dt_export, dt_crc, by = "age", all = TRUE)
  
  # *****************************************************************************
  ##### 2.4 Alive without CRC ####
  # *****************************************************************************
  
  dt_CRC <- dt_pop_surv_id[, .(n_CRC = .N), by = age_ind_crc_int]
  
  dt_export[dt_CRC, on = .(age = age_ind_crc_int), n_CRC := i.n_CRC]
  
  setorder(dt_export, age)
  
  dt_export[, n_cumsum_CRC := cumsum(n_CRC)]
  
  dt_export[, `:=` (n_alive_nocrc = nrow(dt_pop_surv_id) - shift(n_cumsum_CRC, n = 1, fill = 0))]
  
  # *****************************************************************************
  ##### 2.5 Count the number of screening and follow-up tests and findings  ####
  # *****************************************************************************
  dt_screening <- data.table()
  
  system.time(
    for (n_screening in 1:length(screening_years)) {
      # screening_years <- c(50,60,70,80)
      screen <- as.character(n_screening)
    
      # NEW: dynamically pull all ages that occur for this screen (including the deferred age, when applicable)
      ages_this_screen <- unique(
        dt_pop_surv_id[,
          get(paste0("age_at_", screen))
        ]
      )
      
      # inner loop over possible ages for this screen
      for (age in ages_this_screen) {
        
        # optional: restrict to people screened at exactly this age
        dt_age <- dt_pop_surv_id[
          get(paste0("age_at_", screen)) == age
        ]
        
        if (follow_up == FALSE) {
          
          n_pos <- dt_age[get(paste0("is_screened_at_", screen)) == 1 &
                            get(paste0("is_detected_at_", screen)) == 1, .N]
          
          # Most advanced lesion
          n_pos_small  <- dt_age[get(paste0("is_screened_at_", screen)) == 1 &
                                   get(paste0("is_detected_at_", screen)) == 1 &
                                   get(paste0("largest_size_dx_at_", screen)) == "small", .N]
          
          n_pos_medium <- dt_age[get(paste0("is_screened_at_", screen)) == 1 &
                                   get(paste0("is_detected_at_", screen)) == 1 &
                                   get(paste0("largest_size_dx_at_", screen)) == "medium", .N]
          
          n_pos_large  <- dt_age[get(paste0("is_screened_at_", screen)) == 1 &
                                   get(paste0("is_detected_at_", screen)) == 1 &
                                   get(paste0("largest_size_dx_at_", screen)) == "large", .N]
          
          n_pos_S1 <- dt_age[get(paste0("is_screened_at_", screen)) == 1 &
                               get(paste0("is_detected_at_", screen)) == 1 &
                               get(paste0("largest_size_dx_at_", screen)) == "S1", .N]
          n_pos_S2 <- dt_age[get(paste0("is_screened_at_", screen)) == 1 &
                               get(paste0("is_detected_at_", screen)) == 1 &
                               get(paste0("largest_size_dx_at_", screen)) == "S2", .N]
          n_pos_S3 <- dt_age[get(paste0("is_screened_at_", screen)) == 1 &
                               get(paste0("is_detected_at_", screen)) == 1 &
                               get(paste0("largest_size_dx_at_", screen)) == "S3", .N]
          n_pos_S4 <- dt_age[get(paste0("is_screened_at_", screen)) == 1 &
                               get(paste0("is_detected_at_", screen)) == 1 &
                               get(paste0("largest_size_dx_at_", screen)) == "S4", .N]
          
          if (n_pos != (n_pos_small+n_pos_medium+n_pos_large+n_pos_S1+n_pos_S2+n_pos_S3+n_pos_S4)) {
            warning(paste0("Check counts at age ", age,
                           ": sum of adenoma+CRC does not match total positives."))
          }
          
          n_neg <- dt_age[get(paste0("is_screened_at_", screen)) == 1 &
                            (get(paste0("is_detected_at_", screen)) == 0 |
                               is.na(get(paste0("is_detected_at_", screen)))), .N]
          
          n_pos_fucol <- 0
          n_neg_fucol <- 0
          
          # Most advanced adenoma
          n_pos_small_ad  <- dt_age[get(paste0("is_screened_at_", screen)) == 1 &
                                   get(paste0("is_detected_at_", screen)) == 1 &
                                   get(paste0("largest_ad_size_dx_at_", screen)) == "small", .N]
          
          n_pos_medium_ad <- dt_age[get(paste0("is_screened_at_", screen)) == 1 &
                                   get(paste0("is_detected_at_", screen)) == 1 &
                                   get(paste0("largest_ad_size_dx_at_", screen)) == "medium", .N]
          
          n_pos_large_ad  <- dt_age[get(paste0("is_screened_at_", screen)) == 1 &
                                   get(paste0("is_detected_at_", screen)) == 1 &
                                   get(paste0("largest_ad_size_dx_at_", screen)) == "large", .N]
          
          # Most advanced SSL
          n_pos_small_ssp  <- dt_age[get(paste0("is_screened_at_", screen)) == 1 &
                                      get(paste0("is_detected_at_", screen)) == 1 &
                                      get(paste0("largest_ssl_size_dx_at_", screen)) == "small", .N]
          
          n_pos_medium_ssp <- dt_age[get(paste0("is_screened_at_", screen)) == 1 &
                                      get(paste0("is_detected_at_", screen)) == 1 &
                                      get(paste0("largest_ssl_size_dx_at_", screen)) == "medium", .N]
          
          n_pos_large_ssp  <- dt_age[get(paste0("is_screened_at_", screen)) == 1 &
                                      get(paste0("is_detected_at_", screen)) == 1 &
                                      get(paste0("largest_ssl_size_dx_at_", screen)) == "large", .N]
          
          
        } else if (follow_up == TRUE) {
          
          # Total number of positive stool-based tests
          n_pos <- dt_age[get(paste0("is_screened_at_", screen)) == 1 & #Individuals get a stool-based test
                                    (get(paste0("is_pre_detected_at_", screen)) == 1 #Individuals have a positive stool-based test
                                    ), .N]

          # Among all of the stool-based test, calculate the findings at follow-up colonoscopy
          # Total number of positive follow-up colonoscopies
          n_pos_fucol <- dt_age[get(paste0("is_sent_to_fu_at_", screen)) == 1 & # Individuals who are sent to follow-up colonoscopy
                                          (get(paste0("is_detected_at_", screen)) == 1), .N] # Individuals with a positive follow-up colonoscopy

            # The number of test positive follow-up colonoscopies, by most advanced adenoma detected at this follow-up colonoscopy (excluding patients with cancer detected at this screen)
            n_pos_small <- dt_age[get(paste0("is_sent_to_fu_at_", screen)) == 1 & # Individuals who are sent to follow-up colonoscopy
                                                  (get(paste0("is_detected_at_", screen)) == 1) &
                                            get(paste0("largest_size_dx_at_", screen)) == "small", .N]

            n_pos_medium <- dt_age[get(paste0("is_sent_to_fu_at_", screen)) == 1 & # Individuals who are sent to follow-up colonoscopy
                                                   (get(paste0("is_detected_at_", screen)) == 1)&
                                             get(paste0("largest_size_dx_at_", screen)) == "medium", .N]

            n_pos_large <- dt_age[get(paste0("is_sent_to_fu_at_", screen)) == 1 & # Individuals who are sent to follow-up colonoscopy
                                                  (get(paste0("is_detected_at_", screen)) == 1)&
                                            get(paste0("largest_size_dx_at_", screen)) == "large", .N]

            n_pos_S1 <- dt_age[get(paste0("is_sent_to_fu_at_", screen)) == 1 & # Individuals who are sent to follow-up colonoscopy
                                            (get(paste0("is_detected_at_", screen)) == 1) &
                                            get(paste0("largest_size_dx_at_", screen)) == "S1", .N]

            n_pos_S2 <- dt_age[get(paste0("is_sent_to_fu_at_", screen)) == 1 & # Individuals who are sent to follow-up colonoscopy
                                             (get(paste0("is_detected_at_", screen)) == 1)&
                                             get(paste0("largest_size_dx_at_", screen)) == "S2", .N]

            n_pos_S3 <- dt_age[get(paste0("is_sent_to_fu_at_", screen)) == 1 & # Individuals who are sent to follow-up colonoscopy
                                            (get(paste0("is_detected_at_", screen)) == 1)&
                                            get(paste0("largest_size_dx_at_", screen)) == "S3", .N]

            n_pos_S4 <- dt_age[get(paste0("is_sent_to_fu_at_", screen)) == 1 & # Individuals who are sent to follow-up colonoscopy
                                         (get(paste0("is_detected_at_", screen)) == 1)&
                                         get(paste0("largest_size_dx_at_", screen)) == "S4", .N]

            if(n_pos_fucol != (n_pos_small+n_pos_medium+n_pos_large+n_pos_S1+n_pos_S2+n_pos_S3+n_pos_S4) ){
              warning(paste0("Check the calculation of n_pos at age ", age, ". The sum of positive tests by adenoma size and cancer stage does not equal the total number of positive tests."))
            }
          # Total number of negative follow-up colonoscopies
          n_neg_fucol <- n_pos - n_pos_fucol

          # Total number of negative stool-based tests
          n_neg <- dt_age[get(paste0("is_screened_at_", screen)) == 1 & #Individuals get a stool-based test
                                    (get(paste0("is_pre_detected_at_", screen)) == 0 | #Individuals have a negative stool-based test (0 or NA)
                                       is.na(get(paste0("is_pre_detected_at_", screen)))), .N]

          # Most advanced adenoma
          n_pos_small_ad  <- dt_age[get(paste0("is_sent_to_fu_at_", screen)) == 1 &
                                      get(paste0("is_detected_at_", screen)) == 1 &
                                      get(paste0("largest_ad_size_dx_at_", screen)) == "small", .N]
          
          n_pos_medium_ad <- dt_age[get(paste0("is_sent_to_fu_at_", screen)) == 1 &
                                      get(paste0("is_detected_at_", screen)) == 1 &
                                      get(paste0("largest_ad_size_dx_at_", screen)) == "medium", .N]
          
          n_pos_large_ad  <- dt_age[get(paste0("is_sent_to_fu_at_", screen)) == 1 &
                                      get(paste0("is_detected_at_", screen)) == 1 &
                                      get(paste0("largest_ad_size_dx_at_", screen)) == "large", .N]
          
          # Most advanced SSL
          n_pos_small_ssp  <- dt_age[get(paste0("is_sent_to_fu_at_", screen)) == 1 &
                                       get(paste0("is_detected_at_", screen)) == 1 &
                                       get(paste0("largest_ssl_size_dx_at_", screen)) == "small", .N]
          
          n_pos_medium_ssp <- dt_age[get(paste0("is_sent_to_fu_at_", screen)) == 1 &
                                       get(paste0("is_detected_at_", screen)) == 1 &
                                       get(paste0("largest_ssl_size_dx_at_", screen)) == "medium", .N]
          
          n_pos_large_ssp  <- dt_age[get(paste0("is_sent_to_fu_at_", screen)) == 1 &
                                       get(paste0("is_detected_at_", screen)) == 1 &
                                       get(paste0("largest_ssl_size_dx_at_", screen)) == "large", .N]
          
        }else{
          warning("The strategy dataframe does not specify if this modality has a follow-up colonoscopy (aka confirmation = TRUE)")
        }
        
            # Assign the counts
            if (screening_modality == "COL"){

              n_pos_scrcol <- n_pos
              n_neg_scrcol <- n_neg
              n_pos_fit <- 0
              n_neg_fit <- 0
              n_pos_sensa <- 0
              n_neg_sensa <- 0
              n_pos_sdna <- 0
              n_neg_sdna <- 0
              n_pos_srna <- 0
              n_neg_srna <- 0
              n_pos_blood <- 0
              n_neg_blood <- 0
              n_pos_sig <- 0
              n_neg_sig <- 0
              n_pos_ctc <- 0
              n_neg_ctc <- 0
              

            }else if (screening_modality == "FIT"){

              n_pos_scrcol <- 0
              n_neg_scrcol <- 0
              n_pos_fit <- n_pos
              n_neg_fit <- n_neg
              n_pos_sensa <- 0
              n_neg_sensa <- 0
              n_pos_sdna <- 0
              n_neg_sdna <- 0
              n_pos_srna <- 0
              n_neg_srna <- 0
              n_pos_blood <- 0
              n_neg_blood <- 0
              n_pos_sig <- 0
              n_neg_sig <- 0
              n_pos_ctc <- 0
              n_neg_ctc <- 0

            }else if (screening_modality == "HSgFOBT"){

              n_pos_scrcol <- 0
              n_neg_scrcol <- 0
              n_pos_fit <- 0
              n_neg_fit <- 0
              n_pos_sensa <- n_pos
              n_neg_sensa <- n_neg
              n_pos_sdna <- 0
              n_neg_sdna <- 0
              n_pos_srna <- 0
              n_neg_srna <- 0
              n_pos_blood <- 0
              n_neg_blood <- 0
              n_pos_sig <- 0
              n_neg_sig <- 0
              n_pos_ctc <- 0
              n_neg_ctc <- 0

            }else if (screening_modality == "mtsDNA"){

              n_pos_scrcol <- 0
              n_neg_scrcol <- 0
              n_pos_fit <- 0
              n_neg_fit <- 0
              n_pos_sensa <- 0
              n_neg_sensa <- 0
              n_pos_sdna <- n_pos
              n_neg_sdna <- n_neg
              n_pos_srna <- 0
              n_neg_srna <- 0
              n_pos_blood <- 0
              n_neg_blood <- 0
              n_pos_sig <- 0
              n_neg_sig <- 0
              n_pos_ctc <- 0
              n_neg_ctc <- 0

            }else if (screening_modality == "SRNA"){
              
              n_pos_scrcol <- 0
              n_neg_scrcol <- 0
              n_pos_fit <- 0
              n_neg_fit <- 0
              n_pos_sensa <- 0
              n_neg_sensa <- 0
              n_pos_sdna <- 0
              n_neg_sdna <- 0
              n_pos_srna <- n_pos
              n_neg_srna <- n_neg
              n_pos_blood <- 0
              n_neg_blood <- 0
              n_pos_sig <- 0
              n_neg_sig <- 0
              n_pos_ctc <- 0
              n_neg_ctc <- 0
              
            }else if (screening_modality == "BLOOD"){
              
              n_pos_scrcol <- 0
              n_neg_scrcol <- 0
              n_pos_fit <- 0
              n_neg_fit <- 0
              n_pos_sensa <- 0
              n_neg_sensa <- 0
              n_pos_sdna <- 0
              n_neg_sdna <- 0
              n_pos_srna <- 0
              n_neg_srna <- 0
              n_pos_blood <- n_pos
              n_neg_blood <- n_neg
              n_pos_sig <- 0
              n_neg_sig <- 0
              n_pos_ctc <- 0
              n_neg_ctc <- 0
              
            }else if (screening_modality == "SIG"){

              n_pos_scrcol <- 0
              n_neg_scrcol <- 0
              n_pos_fit <- 0
              n_neg_fit <- 0
              n_pos_sensa <- 0
              n_neg_sensa <- 0
              n_pos_sdna <- 0
              n_neg_sdna <- 0
              n_pos_srna <- 0
              n_neg_srna <- 0
              n_pos_blood <- 0
              n_neg_blood <- 0
              n_pos_sig <- n_pos
              n_neg_sig <- n_neg
              n_pos_ctc <- 0
              n_neg_ctc <- 0

            }else if (screening_modality == "CTC"){

              n_pos_scrcol <- 0
              n_neg_scrcol <- 0
              n_pos_fit <- 0
              n_neg_fit <- 0
              n_pos_sensa <- 0
              n_neg_sensa <- 0
              n_pos_sdna <- 0
              n_neg_sdna <- 0
              n_pos_srna <- 0
              n_neg_srna <- 0
              n_pos_blood <- 0
              n_neg_blood <- 0
              n_pos_sig <- 0
              n_neg_sig <- 0
              n_pos_ctc <- n_pos
              n_neg_ctc <- n_neg

            }else{

              # Print a warning if the user is trying a modality that we have not adjusted for
              if(screening_modality != "NoScreening"){
                warning(paste0("The modality '",screening_modality,"' has not been implemented yet. "))
              }

              n_pos_scrcol <- 0
              n_neg_scrcol <- 0
              n_pos_fit <- 0
              n_neg_fit <- 0
              n_pos_sensa <- 0
              n_neg_sensa <- 0
              n_pos_sdna <- 0
              n_neg_sdna <- 0
              n_pos_srna <- 0
              n_neg_srna <- 0
              n_pos_blood <- 0
              n_neg_blood <- 0
              n_pos_sig <- 0
              n_neg_sig <- 0
              n_pos_ctc <- 0
              n_neg_ctc <- 0

              n_pos_fucol <- 0
              n_neg_fucol <- 0
              n_pos_small <- 0
              n_pos_medium <- 0
              n_pos_large <- 0

            }
        
        dt_screening <- rbind(
          dt_screening,
          data.table(
            age = age,
            n_pos_scrcol = n_pos_scrcol,
            n_neg_scrcol = n_neg_scrcol,
            n_pos_fit = n_pos_fit,
            n_neg_fit = n_neg_fit,
            n_pos_sensa = n_pos_sensa,
            n_neg_sensa = n_neg_sensa,
            n_pos_sdna = n_pos_sdna,
            n_neg_sdna = n_neg_sdna,
            n_pos_srna = n_pos_srna,
            n_neg_srna = n_neg_srna,
            n_pos_blood = n_pos_blood,
            n_neg_blood = n_neg_blood,
            n_pos_sig = n_pos_sig,
            n_neg_sig = n_neg_sig,
            n_pos_ctc = n_pos_ctc,
            n_neg_ctc = n_neg_ctc,
            n_pos_fucol = n_pos_fucol,
            n_neg_fucol = n_neg_fucol,
            n_adn1to5_scr = n_pos_small_ad,
            n_adn6to9_scr = n_pos_medium_ad,
            n_adn10plus_scr = n_pos_large_ad,
            n_ssp1to5_scr = n_pos_small_ssp,
            n_ssp6to9_scr = n_pos_medium_ssp,
            n_ssp10plus_scr = n_pos_large_ssp,
            n_lesion1to5_scr = n_pos_small,
            n_lesion6to9_scr = n_pos_medium,
            n_lesion10plus_scr = n_pos_large,
            n_ovrdg_adn = 0,
            n_ovrdg_crc = 0,
            n_3ormorenonlargeadn_scr = 0
          ),
          fill = TRUE
        )
        
      } # end inner age loop
      
    } # end screening loop
  )
  
  # If there are rows where the same age appears multiple times, sum them together
  dt_screening <- dt_screening[, lapply(.SD, sum), by = age]
  
  dt_export <- merge(dt_export, dt_screening, by = "age", all = TRUE)

  # *****************************************************************************
  ##### 2.6 Count the number of surveillance colonoscopy tests and findings ####
  # *****************************************************************************
  
  # Remove people who do not have any surveillance colonoscopies to speed up the loop
  dt_surv_results_subset <- dt_pop_surv_id[n_surv > 0]
  
  # If at least one individual has at least one surveillance colonoscopy, then calculate the number of surveillance colonoscopy,
  # Otherwise assume the counts are zero
  if(nrow(dt_surv_results_subset) != 0){
    
    # Initialize empty summary table for storing the age-specific counts
    dt_surv_counts <- data.table(age = integer(), n_pos_survcol = integer(), n_neg_survcol = integer(),
                                 n_adn1to5_surv = integer(), n_adn6to9_surv = integer(), n_adn10plus_surv = integer(),
                                 n_ssp1to5_surv = integer(), n_ssp6to9_surv = integer(), n_ssp10plus_surv = integer(),
                                 n_lesion1to5_surv = integer(), n_lesion6to9_surv = integer(), n_lesion10plus_surv = integer())
    
    # Determine the maximum number of rounds of surveillance
    max_surv_rounds <- max(dt_surv_results_subset$n_surv)
    
    # Loop over the rounds of surveillance
    for (round in 1:max_surv_rounds) {
      #round = 1
      adherent_col <- paste0("is_surveilled_at_", round)
      detection_col <- paste0("is_detected_at_surv_", round)
      age_col <- paste0("age_surv_", round)
      largest_size_col <- paste0("largest_size_dx_at_surv_", round)
      largest_ad_size_col <- paste0("largest_ad_size_dx_at_surv_", round)
      largest_ssl_size_col <- paste0("largest_ssl_size_dx_at_surv_", round)
      
      # Count positives and negatives by age for this round
      dt_temp <- dt_surv_results_subset[(get(adherent_col)) == 1, # A surveillance colonoscopy happened this round
                                        .(n_pos_survcol    = sum(get(detection_col) == 1, na.rm = TRUE), # If something was found at this surveillance colonoscopy, it was a positive test
                                          n_neg_survcol    = sum(get(detection_col) == 0, na.rm = TRUE),
                                          n_adn1to5_surv   = sum(get(detection_col) == 1 & get(largest_ad_size_col) == "small", na.rm = TRUE),
                                          n_adn6to9_surv   = sum(get(detection_col) == 1 & get(largest_ad_size_col) == "medium", na.rm = TRUE),
                                          n_adn10plus_surv = sum(get(detection_col) == 1 & get(largest_ad_size_col) == "large", na.rm = TRUE),
                                          n_ssp1to5_surv   = sum(get(detection_col) == 1 & get(largest_ssl_size_col) == "small", na.rm = TRUE),
                                          n_ssp6to9_surv   = sum(get(detection_col) == 1 & get(largest_ssl_size_col) == "medium", na.rm = TRUE),
                                          n_ssp10plus_surv = sum(get(detection_col) == 1 & get(largest_ssl_size_col) == "large", na.rm = TRUE),
                                          n_lesion1to5_surv   = sum(get(detection_col) == 1 & get(largest_size_col) == "small", na.rm = TRUE),
                                          n_lesion6to9_surv   = sum(get(detection_col) == 1 & get(largest_size_col) == "medium", na.rm = TRUE),
                                          n_lesion10plus_surv = sum(get(detection_col) == 1 & get(largest_size_col) == "large", na.rm = TRUE)
                                          ),by = .(age = get(age_col))
                                        ]
      
      # Merge with the cumulative summary of the surveillance counts
      dt_surv_counts <- merge(dt_surv_counts, dt_temp, by = "age", all = TRUE, suffixes = c("", "_new"))
      
      # For each row in dt_surv_counts, add together the values in the columns "n_pos_survcol" and "n_pos_survcol_new" (ignoring missing values), 
      # and store this total back into the "n_pos_survcol" column. Repeat with all other outcome columns
      dt_surv_counts[, n_pos_survcol := rowSums(.SD, na.rm = TRUE), .SDcols = c("n_pos_survcol", "n_pos_survcol_new")]
      dt_surv_counts[, n_neg_survcol := rowSums(.SD, na.rm = TRUE), .SDcols = c("n_neg_survcol", "n_neg_survcol_new")]
      dt_surv_counts[, n_adn1to5_surv := rowSums(.SD, na.rm = TRUE), .SDcols = c("n_adn1to5_surv", "n_adn1to5_surv_new")]
      dt_surv_counts[, n_adn6to9_surv := rowSums(.SD, na.rm = TRUE), .SDcols = c("n_adn6to9_surv", "n_adn6to9_surv_new")]
      dt_surv_counts[, n_adn10plus_surv := rowSums(.SD, na.rm = TRUE), .SDcols = c("n_adn10plus_surv", "n_adn10plus_surv_new")]
      dt_surv_counts[, n_ssp1to5_surv := rowSums(.SD, na.rm = TRUE), .SDcols = c("n_ssp1to5_surv", "n_ssp1to5_surv_new")]
      dt_surv_counts[, n_ssp6to9_surv := rowSums(.SD, na.rm = TRUE), .SDcols = c("n_ssp6to9_surv", "n_ssp6to9_surv_new")]
      dt_surv_counts[, n_ssp10plus_surv := rowSums(.SD, na.rm = TRUE), .SDcols = c("n_ssp10plus_surv", "n_ssp10plus_surv_new")]
      dt_surv_counts[, n_lesion1to5_surv := rowSums(.SD, na.rm = TRUE), .SDcols = c("n_lesion1to5_surv", "n_lesion1to5_surv_new")]
      dt_surv_counts[, n_lesion6to9_surv := rowSums(.SD, na.rm = TRUE), .SDcols = c("n_lesion6to9_surv", "n_lesion6to9_surv_new")]
      dt_surv_counts[, n_lesion10plus_surv := rowSums(.SD, na.rm = TRUE), .SDcols = c("n_lesion10plus_surv", "n_lesion10plus_surv_new")]
      
      # Drop temp columns now that you have cumulated them
      dt_surv_counts[, c("n_pos_survcol_new", "n_neg_survcol_new","n_adn1to5_surv_new","n_adn6to9_surv_new","n_adn10plus_surv_new","n_ssp1to5_surv_new","n_ssp6to9_surv_new","n_ssp10plus_surv_new" ,"n_lesion1to5_surv_new","n_lesion6to9_surv_new","n_lesion10plus_surv_new") := NULL]
    }
  } else{
    dt_surv_counts <- data.table(age = min_age:max_age)
    dt_surv_counts[, n_pos_survcol := 0]
    dt_surv_counts[, n_neg_survcol := 0]
    dt_surv_counts[, n_adn1to5_surv := 0]
    dt_surv_counts[, n_adn6to9_surv := 0]
    dt_surv_counts[, n_adn10plus_surv := 0]
    dt_surv_counts[, n_ssp1to5_surv := 0]
    dt_surv_counts[, n_ssp6to9_surv := 0]
    dt_surv_counts[, n_ssp10plus_surv := 0]
    dt_surv_counts[, n_lesion1to5_surv := 0]
    dt_surv_counts[, n_lesion6to9_surv := 0]
    dt_surv_counts[, n_lesion10plus_surv := 0]
  }
  
  dt_export <- merge(dt_export, dt_surv_counts, by = "age", all = TRUE)
  
  dt_export[, `:=` (n_3ormorenonlargeadn_surv = 0 )] # not implemented yet
  
  
  # *****************************************************************************
  ##### 2.7 Life-years ####
  # *****************************************************************************
  
  ###### ly -------------------------------------------------------------
  
  ## calculation using aggregate counts using dt_export 
  dt_extra_ly <- dt_pop_surv_id[,.(extra_ly = sum(age_death - age_death_int)),
                                by = .(age_death_int)]
  
  dt_export <- merge(dt_export, dt_extra_ly, by.x = "age", by.y = "age_death_int", all = TRUE)
  
  dt_export[, ly := n_alive - n_deaths + extra_ly]
  
  ## calculation using individual level data using dt_pop_surv_id
  dt_pop_surv_id[, ly_i := age_death - min_age]
  
  ## sanity check
  n_ly_individual_data <- round(sum(as.numeric(dt_pop_surv_id$ly_i, na.rm = TRUE)),0)
  n_ly_aggregate_data <- round(sum(as.numeric(dt_export$ly, na.rm = TRUE)),0)
  
  if(n_ly_individual_data != n_ly_aggregate_data){
    warning("Check the calculation of ly's on the aggregate counts.")
  }
  
  
  ###### ly_nocrc -------------------------------------------------------------
  
  ## calculation using aggregate counts using dt_export 
  dt_extra_ly_nocrc <- dt_pop_surv_id[,.(extra_ly_nocrc = sum(age_ind_crc - age_ind_crc_int)),
                                      by = .(age_ind_crc_int)]
  
  dt_export <- merge(dt_export, dt_extra_ly_nocrc, by.x = "age", by.y = "age_ind_crc_int", all = TRUE)
  
  dt_export[, ly_nocrc := n_alive_nocrc - n_CRC + extra_ly_nocrc]
  
  dt_export[is.na(dt_export)] <- 0

  ## calculation using individual level data using dt_pop_surv_id
  dt_pop_surv_id[, ly_nocrc_i := age_ind_crc - min_age]
  
  ## sanity check
  n_ly_nocrc_individual_data <- round(sum(as.numeric(dt_pop_surv_id$ly_nocrc_i, na.rm = TRUE)),0)
  n_ly_nocrc_aggregate_data <- round(sum(as.numeric(dt_export$ly_nocrc, na.rm = TRUE)),0)
  
  if(n_ly_nocrc_individual_data != n_ly_nocrc_aggregate_data){
    warning("Check the calculation of ly_nocrc's on the aggregate counts.")
  }
  
  
  ###### ly_crc -------------------------------------------------------------
  
  ## calculation using individual level data using dt_pop_surv_id
  dt_pop_surv_id[!is.na(age_CRC_dx), ly_crc_i := age_death - age_CRC_dx]
  
  dt_pop_surv_id[is.na(ly_crc_i), ly_crc_i := 0]
  
  dt_export[is.na(dt_export)] <- 0
  
  # Necessary adjustment so that the crc_allocation_time function is using the correct stage of detection
  dt_pop_surv_id$stage_at_dx <- replace(dt_pop_surv_id$stage_at_dx, dt_pop_surv_id$stage_at_dx == 1, "I")
  dt_pop_surv_id$stage_at_dx <- replace(dt_pop_surv_id$stage_at_dx, dt_pop_surv_id$stage_at_dx == 2, "II")
  dt_pop_surv_id$stage_at_dx <- replace(dt_pop_surv_id$stage_at_dx, dt_pop_surv_id$stage_at_dx == 3, "III")
  dt_pop_surv_id$stage_at_dx <- replace(dt_pop_surv_id$stage_at_dx, dt_pop_surv_id$stage_at_dx == 4, "IV")
  
  dt_NS <- crc_allocation_time(stages = c("I", "II", "III", "IV"),
                               death_cause = c("crc", "oc"),
                               datatable = dt_pop_surv_id[!is.na(age_CRC_dx),],
                               min_age = min_age, 
                               max_age = max_age)
  
  dt_export <- merge(dt_export, dt_NS, by = "age", all = TRUE)
  
  dt_export[is.na(dt_export)] <- 0
  
  ####  3.0 Final export---------------------------------------------------------
  if(output_template_year == 2020){
    final_names <-  c("age", "n_alive", "n_alive_nocrc", 
                      "n_pos_fit", "n_pos_sensa", "n_pos_sdna", "n_pos_sig", "n_pos_ctc", "n_pos_scrcol",
                      "n_neg_fit", "n_neg_sensa", "n_neg_sdna", "n_neg_sig", "n_neg_ctc", "n_neg_scrcol",
                      "n_pos_fucol", "n_neg_fucol", "n_pos_survcol","n_neg_survcol",
                      "n_adn1to5_scr", "n_adn6to9_scr","n_adn10plus_scr",
                      "n_CRCstageI_scr", "n_CRCstageII_scr", "n_CRCstageIII_scr","n_CRCstageIV_scr",
                      "n_adn1to5_surv", "n_adn6to9_surv", "n_adn10plus_surv",
                      "n_CRCstageI_surv", "n_CRCstageII_surv", "n_CRCstageIII_surv","n_CRCstageIV_surv",
                      "n_CRCstageI_sym", "n_CRCstageII_sym","n_CRCstageIII_sym", "n_CRCstageIV_sym",
                      "n_CRC_dth","n_compl_dth","n_oc_dth",
                      "n_ovrdg_crc","n_ovrdg_adn", # did not implement yet
                      "ly","ly_nocrc", 
                      "ly_crcI_initial", "ly_crcII_initial", "ly_crcIII_initial","ly_crcIV_initial",
                      "ly_crcI_contin", "ly_crcII_contin", "ly_crcIII_contin","ly_crcIV_contin",
                      "ly_crcI_termcrc", "ly_crcII_termcrc", "ly_crcIII_termcrc","ly_crcIV_termcrc",
                      "ly_crcI_termoc", "ly_crcII_termoc", "ly_crcIII_termoc","ly_crcIV_termoc",
                      "n_3ormorenonlargeadn_scr","n_3ormorenonlargeadn_surv") # did not implement yet
    
    
  } else if(output_template_year == 2028){
    final_names <-  c("age", "n_alive", "n_alive_nocrc", 
                      "n_pos_fit", "n_pos_sdna", "n_pos_srna", "n_pos_blood", "n_pos_sig", "n_pos_ctc", "n_pos_scrcol",
                      "n_neg_fit", "n_neg_sdna", "n_neg_srna", "n_neg_blood", "n_neg_sig", "n_neg_ctc", "n_neg_scrcol",
                      "n_pos_fucol", "n_neg_fucol", "n_pos_survcol","n_neg_survcol",
                      "n_adn1to5_scr", "n_adn6to9_scr","n_adn10plus_scr",
                      "n_ssp1to5_scr", "n_ssp6to9_scr","n_ssp10plus_scr",
                      "n_lesion1to5_scr", "n_lesion6to9_scr","n_lesion10plus_scr",
                      "n_CRCstageI_scr", "n_CRCstageII_scr", "n_CRCstageIII_scr","n_CRCstageIV_scr",
                      "n_adn1to5_surv", "n_adn6to9_surv", "n_adn10plus_surv",
                      "n_ssp1to5_surv", "n_ssp6to9_surv", "n_ssp10plus_surv",
                      "n_lesion1to5_surv", "n_lesion6to9_surv", "n_lesion10plus_surv",
                      "n_CRCstageI_surv", "n_CRCstageII_surv", "n_CRCstageIII_surv","n_CRCstageIV_surv",
                      "n_CRCstageI_sym", "n_CRCstageII_sym","n_CRCstageIII_sym", "n_CRCstageIV_sym",
                      "n_CRC_dth","n_compl_dth","n_oc_dth",
                      "ly","ly_nocrc", 
                      "ly_crcI_initial", "ly_crcII_initial", "ly_crcIII_initial","ly_crcIV_initial",
                      "ly_crcI_contin", "ly_crcII_contin", "ly_crcIII_contin","ly_crcIV_contin",
                      "ly_crcI_termcrc", "ly_crcII_termcrc", "ly_crcIII_termcrc","ly_crcIV_termcrc",
                      "ly_crcI_termoc", "ly_crcII_termoc", "ly_crcIII_termoc","ly_crcIV_termoc",
                      "n_3ormorenonlargeadn_scr","n_3ormorenonlargeadn_surv") # did not implement yet
    
    
  } else{
    stop("The output_template_year variable must be one of the following: 2020 or 2028.")
  }
  
  dt_export_final <- dt_export[, ..final_names]
  
  dt_export_final[is.na(dt_export_final)] <- 0

  return(dt_export_final)
}
