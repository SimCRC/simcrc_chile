#' ProcessUSPSTFOutput
#' 
#' \code{ProcessUSPSTFOutput} is a function that reads in the USPSTF model output templates, applies disutilities and costs, and combines all summarized data into one table.
#' 
#' @param analysis_folder This should be the folder where subfolder "RawModelOutput..." lives.
#' @param input_folder Set the input_folder (refers to inputs for this post-processing code). If the input_folder is within the analysis folder you do not need to specify the entire path. However, if not you can specify a complete path.
#' @param folder_for_output Where do you want the R Output to go? This can be a new or existing directory.
#' @param include_SimCRC Are you including SimCRC model results?
#' @param include_MISCAN Are you including MISCAN model results?
#' @param include_CRCSPIN Are you including CRC-SPIN model results?
#' @param input_prefix What is the prefix for the input files? (e.g., USPSTF2020)
#' @param col_infl_rate CISNET models do not account for incomplete colonoscopies. Instead we assume 5% need to be repeated (due to poor prep, to try again at reaching cecum, etc)
#' @param discount_rate Discount rate for costs and utilities
#' @param col_spec_adj CISNET models do not currently simulate non-adenomatous polyps. We account for their detection and removal in post-processing using the lack of specificity.
#' @param crc_care_costs_file Care costs input (within input_folder)
#' @param screen_costs_file Screen costs input (within input_folder)
#' @param crc_care_disutility_file Care utilities (within input_folder)
#' @param screen_disutility_file Screen utilities (within input_folder)
#' @param general_health_utility_weights_file Age utilities (within input_folder)
#' @param selected_outcomes_for_model_data_file What outcomes do you want (within input_folder) 
#' @param model_run_data_tag You can add a tag to the model_run_data file!
#' @param ProcessUSPSTFOutput_script_location = "" Where is the script located?
#' 
#' @import data.table
#' @import dplyr
#' @import magrittr
#' @import stringr
#' @import tibble
#' @import writexl
#' 
#' @return A data.frame containing summarized data.
#' @examples
#' ProcessUSPSTFOutput(analysis_folder = "R1-Basecase", input_folder = "R_Processing_Inputs", folder_for_output = "output")
#' @export
ProcessUSPSTFOutput <- function(analysis_folder = c("R1-Basecase", "R2-ElevatedRisk", "R3-ElevatedRisk", "R-model"), # Options for analysis folder
                                output_template_year = 2020, ## either 2020 or 2028 depending on the output you are generating
                                psa_enabled = FALSE, # Not able to run this yet!
                                psa_params = df_psa_cea, # Not set up yet!
                                input_folder = "R_Processing_Inputs", ## Set the input_folder (refers to inputs for this post-processing code)
                                ## If the input_folder is within the analysis folder you do not need to specify the entire path. However, if not you 
                                ## can specify a complete path. 
                                folder_for_output = "R_Output",## Where do you want the R Output to go? This can be a new or existing directory. 
                                include_SimCRC = T, ## Are you including SimCRC-C++ model results? 
                                include_SimCRC_R = T, ## Are you including SimCRC-R model results? 
                                include_MISCAN = F, ## Are you including MISCAN model results?
                                include_CRCSPIN = F , ## Are you including CRC-SPIN model results?
                                first_age_of_interest = 40, # What is the 1st age of interest for analyses? 
                                input_prefix = "2020SDA", # Prefix to the files in need of processing.
                                col_infl_rate = 1.05,  # CISNET models do not account for incomplete colonoscopies. Instead we assume 5% need to be repeated (due to poor prep, to try again at reaching cecum, etc)
                                discount_rate = 0.03,
                                col_spec_adj = 0.86, # CISNET models do not currently simulate non-adenomatous polyps. We account for their detection and removal in post-processing using the lack of specificity.
                                crc_care_costs_file = "crc_care_costs.csv", # care costs input  (within input_folder)
                                screen_costs_file = "screen_costs.csv", # screen costs input (within input_folder)
                                crc_care_disutility_file = "crc_care_utility_loss.csv", # care utilities (within input_folder)
                                general_health_utility_weights_file = "GeneralHealthUtilityWeightsByAge.csv", # age utilities (within input_folder)
                                screen_disutility_file = "screen_utility_loss_WithStoolTestValues.csv", # screen utilities (within input_folder)
                                selected_outcomes_for_model_data_file = "model_data_outcomes_boolean.csv", # what outcomes do you want (within input_folder)            
                                model_run_data_tag = "_R1", # You can add a tag to the model_run_data file!
                                ProcessUSPSTFOutput_script_location = ""){
  
  
  #analysis_folder <- match.arg(analysis_folder)
  
  # Check if folder_for_output is missing or empty
  if (missing(folder_for_output) || folder_for_output == "") {
    stop("Error: You must specify a folder_for_output.")
  }
  
  # Check if required files are present in the input folder
  
  source("R/complications.R")
  
  # If the output template year is not one of the pre-specified years then throw an error.
  if(!(output_template_year %in% c(2020, 2028))){
    stop("output_template_year must be either 2020 or 2028")
  }
  
  ## Do you have all the necessary packages? If not, install. 
  list.of.packages <- c("data.table", "dplyr", "stringr", "tibble")
  missing.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(missing.packages)){install.packages(missing.packages)}
  # Load the necessary package libraries
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tibble)
  #library(xlsx)
  
  ## Set the working directory to the analysis folder
  #setwd(analysis_folder)
  time <- gsub(':','',Sys.time())
  
  ## Create the directory for the R output if it does not exist. If the folder already exists it will not write the folder. 
  folder_for_output <- paste0(analysis_folder,"/",folder_for_output)
  if(!(dir.exists(folder_for_output))){
    dir.create(folder_for_output, showWarnings = FALSE)
  }
  
  ## Read in costs and disutilites for quality-of-life (these files should be located in input_folder)
  crc_care_costs <- read.csv(file.path(input_folder, crc_care_costs_file), row.names = 1)
  screen_costs <- read.csv(file.path(input_folder, screen_costs_file), row.names = 1)
  crc_care_utility_loss <- read.csv(file.path(input_folder, crc_care_disutility_file), row.names = 1)
  screen_utility_loss <- read.csv(file.path(input_folder, screen_disutility_file), row.names = 1)
  
  # If PSA is enabled, overwrite the tables you read in with the values from the psa_params table. 
  # This will allow you to do PSA by just changing the values in the costs and utility tables at the beginning of the code.
  if(psa_enabled){
    crc_care_costs["Initial Phase (<65)","Stage_I"] <- psa_params$CRC_stage_I_inital_cost_under65
    crc_care_costs["Initial Phase (65+)","Stage_I"] <- psa_params$CRC_stage_I_inital_cost
    crc_care_costs["Initial Phase (<65)","Stage_II"] <- psa_params$CRC_stage_II_inital_cost_under65
    crc_care_costs["Initial Phase (65+)","Stage_II"] <- psa_params$CRC_stage_II_inital_cost
    crc_care_costs["Initial Phase (<65)","Stage_III"] <- psa_params$CRC_stage_III_inital_cost_under65
    crc_care_costs["Initial Phase (65+)","Stage_III"] <- psa_params$CRC_stage_III_inital_cost
    crc_care_costs["Initial Phase (<65)","Stage_IV"] <- psa_params$CRC_stage_IV_inital_cost_under65
    crc_care_costs["Initial Phase (65+)","Stage_IV"] <- psa_params$CRC_stage_IV_inital_cost
    
    crc_care_costs["Continuing Phase (<65)","Stage_I"] <- psa_params$CRC_stage_I_cont_cost_under65
    crc_care_costs["Continuing Phase (65+)","Stage_I"] <- psa_params$CRC_stage_I_cont_cost
    crc_care_costs["Continuing Phase (<65)","Stage_II"] <- psa_params$CRC_stage_II_cont_cost_under65
    crc_care_costs["Continuing Phase (65+)","Stage_II"] <- psa_params$CRC_stage_II_cont_cost
    crc_care_costs["Continuing Phase (<65)","Stage_III"] <- psa_params$CRC_stage_III_cont_cost_under65
    crc_care_costs["Continuing Phase (65+)","Stage_III"] <- psa_params$CRC_stage_III_cont_cost
    crc_care_costs["Continuing Phase (<65)","Stage_IV"] <- psa_params$CRC_stage_IV_cont_cost_under65
    crc_care_costs["Continuing Phase (65+)","Stage_IV"] <- psa_params$CRC_stage_IV_cont_cost
    
    crc_care_costs["Terminal Phase (Death CRC) (<65)","Stage_I"] <- psa_params$CRC_stage_I_termcrc_cost_under65
    crc_care_costs["Terminal Phase (Death CRC) (65+)","Stage_I"] <- psa_params$CRC_stage_I_termcrc_cost
    crc_care_costs["Terminal Phase (Death CRC) (<65)","Stage_II"] <- psa_params$CRC_stage_II_termcrc_cost_under65
    crc_care_costs["Terminal Phase (Death CRC) (65+)","Stage_II"] <- psa_params$CRC_stage_II_termcrc_cost
    crc_care_costs["Terminal Phase (Death CRC) (<65)","Stage_III"] <- psa_params$CRC_stage_III_termcrc_cost_under65
    crc_care_costs["Terminal Phase (Death CRC) (65+)","Stage_III"] <- psa_params$CRC_stage_III_termcrc_cost
    crc_care_costs["Terminal Phase (Death CRC) (<65)","Stage_IV"] <- psa_params$CRC_stage_IV_termcrc_cost_under65
    crc_care_costs["Terminal Phase (Death CRC) (65+)","Stage_IV"] <- psa_params$CRC_stage_IV_termcrc_cost
    
    crc_care_costs["Terminal Phase (Death Other) (<65)","Stage_I"] <- psa_params$CRC_stage_I_termoc_cost_under65
    crc_care_costs["Terminal Phase (Death Other) (65+)","Stage_I"] <- psa_params$CRC_stage_I_termoc_cost
    crc_care_costs["Terminal Phase (Death Other) (<65)","Stage_II"] <- psa_params$CRC_stage_II_termoc_cost_under65
    crc_care_costs["Terminal Phase (Death Other) (65+)","Stage_II"] <- psa_params$CRC_stage_II_termoc_cost
    crc_care_costs["Terminal Phase (Death Other) (<65)","Stage_III"] <- psa_params$CRC_stage_III_termoc_cost_under65
    crc_care_costs["Terminal Phase (Death Other) (65+)","Stage_III"] <- psa_params$CRC_stage_III_termoc_cost
    crc_care_costs["Terminal Phase (Death Other) (<65)","Stage_IV"] <- psa_params$CRC_stage_IV_termoc_cost_under65
    crc_care_costs["Terminal Phase (Death Other) (65+)","Stage_IV"] <- psa_params$CRC_stage_IV_termoc_cost
    
    screen_costs["FIT","Cost..LessThan65."] <- psa_params$FIT_cost_under65
    screen_costs["FIT","Cost..65Plus."] <- psa_params$FIT_cost
    
    screen_costs["SDNA","Cost..LessThan65."] <- psa_params$SDNA_cost_under65
    screen_costs["SDNA","Cost..65Plus."] <- psa_params$SDNA_cost
    
    screen_costs["SIG","Cost..LessThan65."] <- psa_params$SIG_cost_under65
    screen_costs["SIG","Cost..65Plus."] <- psa_params$SIG_cost
    
    screen_costs["CTC","Cost..LessThan65."] <- psa_params$CTC_cost_under65
    screen_costs["CTC","Cost..65Plus."] <- psa_params$CTC_cost
    
    screen_costs["Screening colonoscopy without polypectomy","Cost..LessThan65."]<- psa_params$scr_COL_wo_polyp_cost_under65
    screen_costs["Screening colonoscopy without polypectomy","Cost..65Plus."]<- psa_params$scr_COL_wo_polyp_cost
    
    screen_costs["Diagnostic colonoscopy without polypectomy","Cost..LessThan65."]<- psa_params$fu_COL_wo_polyp_cost_under65
    screen_costs["Diagnostic colonoscopy without polypectomy","Cost..65Plus."]<- psa_params$fu_COL_wo_polyp_cost
    
    screen_costs["Surveillance colonoscopy without polypectomy","Cost..LessThan65."]<- psa_params$surv_COL_wo_polyp_cost_under65
    screen_costs["Surveillance colonoscopy without polypectomy","Cost..65Plus."]<- psa_params$surv_COL_wo_polyp_cost
    
    screen_costs["Colonoscopy with polypectomy","Cost..LessThan65."]<- psa_params$COL_w_polyp_cost_under65
    screen_costs["Colonoscopy with polypectomy","Cost..65Plus."]<- psa_params$COL_w_polyp_cost
    
    screen_costs["Colonoscopy for diagnosis of a cancer by symptoms","Cost..LessThan65."]<- psa_params$symp_COL_w_polyp_cost_under65
    screen_costs["Colonoscopy for diagnosis of a cancer by symptoms","Cost..65Plus."]<- psa_params$symp_COL_w_polyp_cost
    
    screen_costs["Fatal perforation","Cost..LessThan65."] <- psa_params$fatal_perf_cost_under65
    screen_costs["Fatal perforation","Cost..65Plus."] <- psa_params$fatal_perf_cost
    
    screen_costs["Serious GI complication","Cost..LessThan65."] <- psa_params$serious_gi_cost_under65
    screen_costs["Serious GI complication","Cost..65Plus."] <- psa_params$serious_gi_cost
    
    screen_costs["Other GI complication","Cost..LessThan65."] <- psa_params$nonserious_gi_cost_under65
    screen_costs["Other GI complication","Cost..65Plus."] <- psa_params$nonserious_gi_cost
    
    screen_costs["Cardiovascular complication","Cost..LessThan65."] <- psa_params$cardio_cost_under65
    screen_costs["Cardiovascular complication","Cost..65Plus."] <- psa_params$cardio_cost
    
    
    crc_care_utility_loss["Initial Phase","Stage_I"] <- psa_params$CRC_stage_I_inital_disutil
    crc_care_utility_loss["Initial Phase","Stage_II"] <- psa_params$CRC_stage_II_inital_disutil
    crc_care_utility_loss["Initial Phase","Stage_III"] <- psa_params$CRC_stage_III_inital_disutil
    crc_care_utility_loss["Initial Phase","Stage_IV"] <- psa_params$CRC_stage_IV_inital_disutil
    
    crc_care_utility_loss["Continuing Phase","Stage_I"] <- psa_params$CRC_stage_I_cont_disutil
    crc_care_utility_loss["Continuing Phase","Stage_II"] <- psa_params$CRC_stage_II_cont_disutil
    crc_care_utility_loss["Continuing Phase","Stage_III"] <- psa_params$CRC_stage_III_cont_disutil
    crc_care_utility_loss["Continuing Phase","Stage_IV"] <- psa_params$CRC_stage_IV_cont_disutil
    
    crc_care_utility_loss["Terminal Phase (Death CRC)","Stage_I"] <- psa_params$CRC_stage_I_termcrc_disutil
    crc_care_utility_loss["Terminal Phase (Death CRC)","Stage_II"] <- psa_params$CRC_stage_II_termcrc_disutil
    crc_care_utility_loss["Terminal Phase (Death CRC)","Stage_III"] <- psa_params$CRC_stage_III_termcrc_disutil
    crc_care_utility_loss["Terminal Phase (Death CRC)","Stage_IV"] <- psa_params$CRC_stage_IV_termcrc_disutil
    
    crc_care_utility_loss["Terminal Phase (Death Other)","Stage_I"] <- psa_params$CRC_stage_I_termoc_disutil
    crc_care_utility_loss["Terminal Phase (Death Other)","Stage_II"] <- psa_params$CRC_stage_II_termoc_disutil
    crc_care_utility_loss["Terminal Phase (Death Other)","Stage_III"] <- psa_params$CRC_stage_III_termoc_disutil
    crc_care_utility_loss["Terminal Phase (Death Other)","Stage_IV"] <- psa_params$CRC_stage_IV_termoc_disutil
    
    screen_utility_loss["FIT positive",] <- psa_params$positive_FIT_disutil
    screen_utility_loss["FIT negative",] <- psa_params$negative_FIT_disutil
    
    screen_utility_loss["SDNA positive",] <- psa_params$positive_SDNA_disutil
    screen_utility_loss["SDNA negative",] <- psa_params$negative_SDNA_disutil
    
    screen_utility_loss["SIG positive",] <- psa_params$positive_SIG_disutil
    screen_utility_loss["SIG negative",] <- psa_params$negative_SIG_disutil
    
    screen_utility_loss["CTC positive",] <- psa_params$positive_CTC_disutil
    screen_utility_loss["CTC negative",] <- psa_params$negative_CTC_disutil
    
    screen_utility_loss["Any COL positive",] <- psa_params$positive_COL_disutil
    screen_utility_loss["Any COL negative",] <- psa_params$negative_COL_disutil
    
    screen_utility_loss["Serious GI complication (COL)",] <- psa_params$serious_gi_disutil
    screen_utility_loss["Other GI complication (COL)",] <- psa_params$nonserious_gi_disutil
    screen_utility_loss["Cardiovascular complication (COL)",] <- psa_params$cardio_disutil
    
    if(output_template_year == 2020){
      
      screen_costs["Hemoccult SENSA","Cost..LessThan65."] <- psa_params$SENSA_cost_under65
      screen_costs["Hemoccult SENSA","Cost..65Plus."] <- psa_params$SENSA_cost
      
      screen_utility_loss["SEN positive",] <- psa_params$positive_SENSA_disutil
      screen_utility_loss["SEN negative",] <- psa_params$negative_SENSA_disutil
      
    } else if(output_template_year == 2028){
      
      screen_costs["SRNA","Cost..LessThan65."] <- psa_params$SRNA_cost_under65
      screen_costs["SRNA","Cost..65Plus."] <- psa_params$SRNA_cost
      
      screen_costs["BLOOD","Cost..LessThan65."] <- psa_params$BLOOD_cost_under65
      screen_costs["BLOOD","Cost..65Plus."] <- psa_params$BLOOD_cost
      
      screen_utility_loss["SRNA positive",] <- psa_params$positive_SRNA_disutil
      screen_utility_loss["SRNA negative",] <- psa_params$negative_SRNA_disutil
      
      screen_utility_loss["BLOOD positive",] <- psa_params$positive_BLOOD_disutil
      screen_utility_loss["BLOOD negative",] <- psa_params$negative_BLOOD_disutil
      
    }
  }
  
  ## Read in the health utility (by age) weights but only keep the lines only for the ages we are interested in.
  general_health_utility_weights <- read.csv(file.path(input_folder, general_health_utility_weights_file)) %>% 
    filter(Age >= first_age_of_interest)  
  ## Create an empty vector and then fill with the names of the CISNET models that we want to include in the analyses. 
  model_list <- vector()  
  if(include_CRCSPIN){model_list <- c(model_list, "CRCSPIN")}
  if(include_MISCAN) {model_list <- c(model_list, "MISCAN")}
  if(include_SimCRC) {model_list <- c(model_list, "SimCRC")}
  if(include_SimCRC_R) {model_list <- c(model_list, "SimCRC_R")}
  
  if(length(model_list) < 1){
    stop("You must specify at least one CISNET Model")
  }
  ## Create an empty object named frontier. This will become the frontier_data table that contains one line of outcomes for each tested screening strategy from each model 
  model_data <- NULL
  #full_QALY_table <- NULL
  
  for(model_x in model_list){
    
    mismatch_population_counter <- 0
    ## The output from each model should be located in a folder "RawModelOutput_MODELNAME". Within that folder we then 
    ## search for all files with structure "...input_prefix....csv" in all subfolders of "RawModelOutput_MODELNAME"
    model_output_files <- list.files(
      path = file.path(analysis_folder, paste0("RawModelOutput_", model_x,"")),
      #pattern = "\\.csv$",  # matches files ending with .csv
      pattern = sprintf(".*%s.*\\.csv", input_prefix),
      full.names = TRUE,
      recursive = TRUE
    )
    
    
    
    ## Find where in the list of files the NoScreening file is located.  
    no_screening_location = grep(pattern = ".*NoScr.*", x = model_output_files)
    
    ## Make sure there is only one no screening file. If there is more than one or it is missing then stop! 
    if(length(no_screening_location) > 1){
      stop(paste("There is more than one NoScreening file (contains 'NoScr' in its name) for model",model_x))
    }else if(length(no_screening_location) < 1){
      stop(paste("There is no NoScreening file (contains 'NoScr' in its name) for model",model_x)) 
    }
    
    ## Move the no screening file to the top of the file list so that it can be processed first. 
    model_output_files = c(model_output_files[no_screening_location],model_output_files[-no_screening_location])
    
    ## For each of the files (tested screening strategies), process the output. 
    for(i in 1:length(model_output_files)){
      ## The first file that is processed should be the no screening file. If this is not the case then STOP! We just set the NS file to be the 1st file so this should never happen. 
      if (i == 1){
        if(!grepl("NoScr",model_output_files[i])){
          ## This isn't the NoScreening file! Throw!
          stop("First file must be the NoScreening file (and contain 'NoScr' in its name")
        }
      }  
      
      ## Keep track of where we are in the processing by printing the current strategy to the screen
      print(paste("Working on",model_output_files[i], sep = " "))
      cat(paste("Working on",model_output_files[i],"\n", sep = " "), file = paste0(folder_for_output, "/", time, "_ProcessUSPSTF_Log.txt"), append = TRUE)
      
      ## Read in the run information for the screening strategy
      if(output_template_year == 2020){
        run_info <- as.data.frame(
          t(read.csv(model_output_files[i],nrows = 14,header = FALSE,row.names = 1,blank.lines.skip = FALSE)), stringsAsFactors = FALSE
        )
      }else if(output_template_year == 2028){
        run_info <- as.data.frame(
          t(read.csv(model_output_files[i],nrows = 12,header = FALSE,row.names = 1,blank.lines.skip = FALSE)), stringsAsFactors = FALSE
        )
      }
      
      if(nrow(run_info)> 1){
        run_info <- run_info[1,]
      }
      
      # We want to change some of the run types. This will help with consistancy across models and will also set them to be exactly how we want them
      # for future output
      if(output_template_year == 2020){
        
        if(!is.na(run_info[,"stooltest_type"])){
          if(run_info[,"stooltest_type"] == "FDNA" |run_info[,"stooltest_type"] == "SDNA" | run_info[,"stooltest_type"] == " SDNA"| run_info[,"stooltest_type"] == "mtsDNA"){
            run_info[,"stooltest_type"] <- "sDNA-FIT"
          }else if(run_info[,"stooltest_type"] == "SEN" |run_info[,"stooltest_type"] == "HSFOBT" | run_info[,"stooltest_type"] == "FOBTHS"|
                   run_info[,"stooltest_type"] == " SEN"){
            run_info[,"stooltest_type"] <- "HSgFOBT"
          }
        }
        
        ## Assign a label to the strategy!
        ## If both stool and structural exams are NAs then this is a no screening run. 
        if(is.na(run_info[,"stooltest_type"]) & is.na(run_info[,"structuralexam_type"])
        ){
          strategy_label <- "NoScreening"
          
          ## If stool test is NA but a structural exam is present then the label should be structuralstrat_startage_stopage_interval
        }else if(is.na(run_info[,"stooltest_type"]) & !is.na(run_info[,"structuralexam_type"])){
          strategy_label <- paste(run_info[,paste0("structuralexam_", c("type", "startage", "stopage", "interval"))], collapse = "_")
          
          ## If stool test is present and structural is NA then the label should be stoolstrat_startage_stopage_interval
        }else if(!is.na(run_info[,"stooltest_type"]) & is.na(run_info[,"structuralexam_type"])){
          strategy_label <- paste(run_info[,paste0("stooltest_", c("type", "startage", "stopage", "interval"))], collapse = "_")
          
          ## Otherwise is must be some combined strategy. Make the label for both stool and structural strategies and then bind them together. Bind them in order of which has 
          ## an earlier start age. 
        }else {
          structural_strategy_label<- paste(run_info[,paste0("structuralexam_", c("type", "startage", "stopage", "interval"))], collapse = "_")
          stool_strategy_label <- paste(run_info[,paste0("stooltest_", c("type", "startage", "stopage", "interval"))], collapse = "_")
          if(as.numeric(run_info[,"structuralexam_startage"]) <= as.numeric(run_info[,"stooltest_startage"])){
            strategy_label <- paste(structural_strategy_label, stool_strategy_label, sep = "_")
          }else{strategy_label <- paste(stool_strategy_label, structural_strategy_label, sep = "_")}
        }
        
        # Skip the number of rows. CRC spin added a seed to the output so need in order to get the correct remaining output we need this. 
        if("seed" %in% colnames(run_info)){
          rows_to_skip <- 15
        }else if (colnames(run_info)[7] == "V7"){ # SimCRC has a bug which adds a line to the SDNA runs in the run info. This just allows us to ignore that row. 
          rows_to_skip <- 15
        }else{
          rows_to_skip <- 14
        }
        
        
      } else if (output_template_year == 2028){
        
        # For consistent naming across models
        if(!is.na(run_info[,"stoolorbloodtest_type"])){
          if(run_info[,"stoolorbloodtest_type"] == "FDNA" |run_info[,"stoolorbloodtest_type"] == "SDNA" | run_info[,"stoolorbloodtest_type"] == " SDNA"| run_info[,"stoolorbloodtest_type"] == "mtsDNA"| run_info[,"stoolorbloodtest_type"] == "dna" | run_info[,"stoolorbloodtest_type"] == "sDNA+FIT"){
            run_info[,"stoolorbloodtest_type"] <- "SDNA"
          }else if(run_info[,"stoolorbloodtest_type"] == "SEN" |run_info[,"stoolorbloodtest_type"] == "HSFOBT" | run_info[,"stoolorbloodtest_type"] == "FOBTHS"|
                   run_info[,"stoolorbloodtest_type"] == " SEN"){
            run_info[,"stoolorbloodtest_type"] <- "SEN"
          }else if(run_info[,"stoolorbloodtest_type"] == "SRNA" | run_info[,"stoolorbloodtest_type"] == "rna"| run_info[,"stoolorbloodtest_type"] == "sRNA+FIT"){
            run_info[,"stoolorbloodtest_type"] <- "SRNA"
          }else if(run_info[,"stoolorbloodtest_type"] == "FIT" | run_info[,"stoolorbloodtest_type"] == "fit"){
            run_info[,"stoolorbloodtest_type"] <- "FIT"
          }else if(run_info[,"stoolorbloodtest_type"] == "BLOOD"){
            run_info[,"stoolorbloodtest_type"] <- "BLOOD"
          }
        }
        
        if(!is.na(run_info[,"structuralexam_type"])){
          if(run_info[,"structuralexam_type"] == "sig" |run_info[,"structuralexam_type"] == "SIG" ){
            run_info[,"structuralexam_type"] <- "SIG"
          }
        }
        
        ## Assign a label to the strategy!
        ## If both stool and structural exams are NAs then this is a no screening run. 
        if(is.na(run_info[,"stoolorbloodtest_type"]) & is.na(run_info[,"structuralexam_type"])
        ){
          strategy_label <- "NoScreening"
          
          ## If stool test is NA but a structural exam is present then the label should be structuralstrat_startage_stopage_interval
        }else if(is.na(run_info[,"stoolorbloodtest_type"]) & !is.na(run_info[,"structuralexam_type"])){
          strategy_label <- do.call(paste, c(run_info[, c("structuralexam_type", "screening_startage", 
                                                          "screening_stopage", "structuralexam_interval")], 
                                             sep = "_"))
          ## If stool test is present and structural is NA then the label should be stoolstrat_startage_stopage_interval
        }else if(!is.na(run_info[,"stoolorbloodtest_type"]) & is.na(run_info[,"structuralexam_type"])){
          strategy_label <- do.call(paste, c(run_info[, c("stoolorbloodtest_type", "screening_startage", 
                                                          "screening_stopage", "stoolorbloodtest_interval")], 
                                             sep = "_"))
          ## Otherwise is must be some combined strategy. Make the label for both stool and structural strategies and then bind them together. Bind them in order of which has 
          ## an earlier start age. 
        }else {
          structural_strategy_label <- do.call(paste, c(run_info[, c("structuralexam_type", "screening_startage", 
                                                                     "screening_stopage", "structuralexam_interval")], 
                                                        sep = "_"))
          stool_strategy_label <- do.call(paste, c(run_info[, c("stoolorbloodtest_type", "screening_startage", 
                                                                "screening_stopage", "stoolorbloodtest_interval")], 
                                                   sep = "_"))
          strategy_label <- paste(structural_strategy_label, stool_strategy_label, sep = "_")
        }
        
        # Skip the number of rows. CRC spin added a seed to the output so need in order to get the correct remaining output we need this. 
        if("seed" %in% colnames(run_info)){
          rows_to_skip <- 15-2
        }else{
          rows_to_skip <- 14-2
        }
        
      }
      
      
      ## Now get the rest of the data from the model output (csv file). Skip the lines pertaining to the run info. Filter so the dataframe only contains the ages of interest. 
      #strategy_undsc <- as.data.frame(fread(model_output_files[i], skip = rows_to_skip)) 
      strategy_undsc <- as.data.frame(fread(model_output_files[i], skip = rows_to_skip)) 
      
      # If any values are NA, then fill with 0s (fix for CRCSPIN)
      strategy_undsc[is.na(strategy_undsc)] <- 0
      
      # If the column name is 'Age' then change it to 'age' (fix for MISCAN)
      if("Age" %in% colnames(strategy_undsc)){
        colnames(strategy_undsc)[colnames(strategy_undsc) == 'Age'] <- 'age'
      }
      
      
      filename <- basename(model_output_files[i])
      #strategy_label <- tools::file_path_sans_ext(filename)
      strategy_label <- ifelse(
        stringr::str_detect(strategy_label, "NoScr"),
        "NoScreening",
        strategy_label
      )
      
      
      colnames(strategy_undsc)[colnames(strategy_undsc) == 'n_3plusnonlargeadn_scr'] <- 'n_3ormorenonlargeadn_scr'
      colnames(strategy_undsc)[colnames(strategy_undsc) == 'n_3plusnonlargeadn_surv'] <- 'n_3ormorenonlargeadn_surv'
      
      dim(strategy_undsc)
      ## ERROR CHECK ## 
      ## Are the correct ages presented in the file? 
      if(nrow(filter(strategy_undsc, strategy_undsc$age > 120)) > 0){
        print(paste("WARNING: The ages presented in the",model_output_files[i],"are incorrect"))
        cat(paste("WARNING: The ages presented in the",model_output_files[i],"are incorrect", "\n"),file = paste0(folder_for_output, "/",time, "_ProcessUSPSTF_Log.txt"), append = TRUE )
      }
      
      strategy_undsc <- strategy_undsc %>% filter(age >= first_age_of_interest)
      strategy_undsc <- strategy_undsc %>% filter(age <= 100)
      ## Determine the number of people alive at the 1st age of interest who are free of colorectal cancer
      n_alive_initial_no_CRC <- strategy_undsc[1,"n_alive_nocrc"] 
      
      ## Determine the number of colonoscopies that would have had a polypectomy. That is all positive colonoscopies (cancer & adenomas)
      ## and some true negative colonoscopies that would have positives for polyps not deemed precancerous. 
      n_cols_with_polypectomy <- rowSums(strategy_undsc[,c("n_pos_scrcol", "n_pos_fucol", "n_pos_survcol", "n_CRCstageI_sym", "n_CRCstageII_sym", "n_CRCstageIII_sym", "n_CRCstageIV_sym")]) + 
        (1-col_spec_adj) * (rowSums(strategy_undsc[,c("n_neg_scrcol", "n_neg_fucol", "n_neg_survcol")])) #This accounts for polypectomy for polyps deemed not to be precancerous
      ## Calculate the number of complications (functions 'cardio', 'serious_GI', and 'other_GI' are at the top of this script)
      strategy_undsc[,"cardio_compl"] = n_cols_with_polypectomy*cardio(strategy_undsc["age"])
      strategy_undsc[,"serious_GI_compl"] = n_cols_with_polypectomy*serious_GI(strategy_undsc["age"])
      strategy_undsc[,"other_GI_compl"] = n_cols_with_polypectomy*other_GI(strategy_undsc["age"])   
      
      ## INFLATE THE NUMBER OF COLONOSCOPIES PERFORMED. 
      ## Loop through screening, follow up, and surveillance colonoscopies and 
      ## calculate the total number of colonoscopies (per type) by adding the positive cols and the negative cols and inflating according to the col_infl_rate. 
      for(col_type in c("scr", "fu", "surv")){
        strategy_undsc[,paste0("n_", col_type, "col_infl")] <- rowSums(strategy_undsc[,c(paste0("n_pos_", col_type, "col"), paste0("n_neg_", col_type, "col"))]) * col_infl_rate
        
        ## determine how many cols (of each type) would have had polypectomy. This includes all positive colonoscopies as well as some negative colonoscopies due to lack of specificity.
        strategy_undsc[,paste0("n_", col_type, "col_with_polypectomy")] <- strategy_undsc[,paste0("n_pos_", col_type, "col")] + (1-col_spec_adj) * strategy_undsc[,paste0("n_neg_", col_type, "col")]
        
        ## determine how many cols would NOT have have polypectomy. (TOTAL - w/ polypectomy)
        strategy_undsc[,paste0("n_", col_type, "col_without_polypectomy")] <- strategy_undsc[,paste0("n_", col_type, "col_infl")] - strategy_undsc[,paste0("n_", col_type, "col_with_polypectomy")]
      } 
      
      ## loop through cancer detection methods and calculate the number of cancers detected by colonoscopy, separate calculations for
      ## clinically (symptom) detected, screening detected, or surveillance detected and add as a new variable to the dataframe. 
      for (dx in c("scr", "surv", "sym")){
        strategy_undsc[,paste0("CRC_", dx ,"_dx")] <- rowSums(strategy_undsc[,c(paste0("n_CRCstage", c("I", "II", "III", "IV"), "_", dx))])
      }
      
      ## Adjust all the numerical output to account the dis-utility of age. Does not apply the utility to the age, or the number of persons alive!
      strategy_undsc_QALY <- cbind(strategy_undsc[,c("age", "n_alive", "n_alive_nocrc")], strategy_undsc[,!names(strategy_undsc) %in% c("age", "n_alive", "n_alive_nocrc")] * general_health_utility_weights[,"GeneralAgeUtilityWeight"])
      
      ## APPLY DISCOUNT
      ## Apply the discount to all the columns except Year, Age, n_alive, and n_alive_nocrc
      strategy_dsc <- cbind(strategy_undsc[,c("age", "n_alive", "n_alive_nocrc")],strategy_undsc[,!names(strategy_undsc) %in% c("age", "n_alive", "n_alive_nocrc") ]/(1+discount_rate)^(strategy_undsc[,"age"]-strategy_undsc[1,"age"]))
      
      strategy_dsc_QALY <- cbind(strategy_undsc_QALY[,c("age", "n_alive", "n_alive_nocrc")],strategy_undsc_QALY[,!names(strategy_undsc_QALY) %in% c("age", "n_alive", "n_alive_nocrc")]/(1+discount_rate)^(strategy_undsc_QALY[,"age"]-strategy_undsc_QALY[1,"age"]))
      
      ## Sum each column (both age discounted and not) to get the total number of each events for the population
      totals_undsc <- colSums(strategy_undsc)
      totals_undsc_QALY <- colSums(strategy_undsc_QALY)
      totals_dsc <- colSums(strategy_dsc)
      totals_dsc_QALY <- colSums(strategy_dsc_QALY)
      
      ## calculate the total number of cancer cases by adding screen detected, surveillance detected, and clinically detected. Additionally
      ## check that the population size is the same among runs. 
      total_crc_cases_undsc <- sum(totals_undsc[paste0("CRC_",c("scr", "surv", "sym"), "_dx")])
      if(i == 1){
        total_crc_cases_no_screening_undsc <- total_crc_cases_undsc
        population_no_screening <- strategy_undsc[1,"n_alive"]
      }else {
        population <- strategy_undsc[1,"n_alive"]
        if(population != population_no_screening){
          mismatch_population_counter <- mismatch_population_counter + 1
        }
      }
      ## calculate the total number of colonscopies that are preformed. 
      n_col_tests_total_infl_undsc <- sum(totals_undsc[c(paste0("n_", c("scr", "fu", "surv"), "col_infl"), "CRC_sym_dx")])
      
      ## For costs we may want to apply different costs for individuals aged less than 65 vs those aged 65 and older. Therefore we create tables with only the relevant ages. 
      totals_under65_undsc <- colSums(strategy_undsc %>% filter(age < 65))
      totals_under65_undsc_QALY <- colSums(strategy_undsc_QALY %>% filter(age < 65))
      totals_under65_dsc <- colSums(strategy_dsc %>% filter(age < 65))
      totals_under65_dsc_QALY <- colSums(strategy_dsc_QALY %>% filter(age < 65))
      
      totals_65plus_undsc <- colSums(strategy_undsc %>% filter(age >= 65))
      totals_65plus_undsc_QALY <- colSums(strategy_undsc_QALY %>% filter(age >= 65))
      totals_65plus_dsc <- colSums(strategy_dsc %>% filter(age >= 65))
      totals_65plus_dsc_QALY <- colSums(strategy_dsc_QALY %>% filter(age >= 65))
      
      ## Now we want to do some calculations (both undiscounted and discounted). 
      for (discounting in c("undsc", "dsc")){
        ## Get the relevant tables (undiscounted or discounted)  
        totals_table <- get(paste0("totals_", discounting)) 
        QALY_table <- get(paste0("totals_", discounting, "_QALY"))
        
        ## Determine the number of screening deaths, CRC deaths, and crc related deaths. 
        assign(paste0("screen_deaths_", discounting), totals_table["n_compl_dth"]) 
        assign(paste0("crc_deaths_", discounting), totals_table["n_CRC_dth"])
        assign(paste0("crc_and_screen_deaths_", discounting), get(paste0("crc_deaths_", discounting)) + get(paste0("screen_deaths_", discounting)))
        assign(paste0("crc_cases_", discounting), totals_table[""])
        ## Determine number of each type of col
        for (col_type in c("scr", "fu", "surv")){
          assign(paste0("n_",col_type , "col_infl_", discounting), totals_table[paste0("n_",col_type , "col_infl")])
        }
        
        ## Determine the total LY lived.
        assign(paste0("LY_",discounting), totals_table["ly"])
        
        
        
        ## Get the number of complications among the age group.
        assign(paste0("n_complications_", discounting), totals_table["n_compl_dth"] + totals_table["cardio_compl"] + totals_table["serious_GI_compl"] + totals_table["other_GI_compl"])
        
        ## Determine the number of stool tests. This is the sum of all the stool tests. 
        if(output_template_year == 2020){
          assign(paste0("n_stool_tests_" ,discounting), sum(totals_table[paste0(c("n_neg_", "n_pos_"), rep(c("fit", "sdna", "sensa"), 2))]))
          ## Determine the total number of each non-colonoscopy screening type.
          for (test in c("sig", "ctc", "fit", "sdna", "sensa")){
            assign(paste0("n_", test, "_tests_" ,discounting), sum(totals_table[paste0(c("n_pos_", "n_neg_"), test)]))
          }
        }else if (output_template_year == 2028){
          stool_cols <- unlist(lapply(
            c("fit", "sdna", "srna", "blood"),
            function(x) c(paste0("n_pos_", x), paste0("n_neg_", x))
          ))
          
          assign(paste0("n_stoolorblood_tests_", discounting),sum(totals_table[stool_cols]))
          ## Determine the total number of each non-colonoscopy screening type.
          for (test in c("sig", "ctc", "fit", "sdna", "srna", "blood")){
            assign(paste0("n_", test, "_tests_" ,discounting), sum(totals_table[paste0(c("n_pos_", "n_neg_"), test)]))
          }
        }
        
        
        
        ## Determine the number of symptom detected CRCs. 
        assign(paste0("n_CRC_sym_dx_", discounting), totals_table["CRC_sym_dx"])
        
        # if (discounting == "undsc"){
        #   full_QALY_table <- rbind(full_QALY_table, cbind(strategy_label, cbind(model_x, t(as.data.frame(QALY_table)))))
        # }
        
        ## Calculate QALYs-- accounting for different disutilities.
        if(output_template_year == 2020){
          novel_stool_test_qaly_loss <-
            screen_utility_loss["SEN positive",] * QALY_table["n_pos_sensa"] +
            screen_utility_loss["SEN negative",] * QALY_table["n_neg_sensa"]
          
        } else {
          novel_stool_test_qaly_loss <-
            screen_utility_loss["SRNA positive",] * QALY_table["n_pos_srna"] +
            screen_utility_loss["SRNA negative",] * QALY_table["n_neg_srna"] +
            screen_utility_loss["BLOOD positive",] * QALY_table["n_pos_blood"] +
            screen_utility_loss["BLOOD negative",] * QALY_table["n_neg_blood"]
        }
        
        assign(paste0("QALY_", discounting), QALY_table["ly"] - 
                 screen_utility_loss["SIG positive",]*QALY_table["n_pos_sig"] -
                 screen_utility_loss["SIG negative",]*QALY_table["n_neg_sig"] -
                 screen_utility_loss["FIT positive",]*QALY_table["n_pos_fit"] -
                 screen_utility_loss["FIT negative",]*QALY_table["n_neg_fit"] -
                 novel_stool_test_qaly_loss -
                 screen_utility_loss["SDNA positive",]*QALY_table["n_pos_sdna"] - 
                 screen_utility_loss["SDNA negative",]*QALY_table["n_neg_sdna"] - 
                 screen_utility_loss["CTC positive",]*QALY_table["n_pos_ctc"] - 
                 screen_utility_loss["CTC negative",]*QALY_table["n_neg_ctc"] -
                 screen_utility_loss["Any COL positive",]*QALY_table["n_scrcol_with_polypectomy"] -
                 screen_utility_loss["Any COL positive",]*QALY_table[ "n_fucol_with_polypectomy"] -
                 screen_utility_loss["Any COL positive",]*QALY_table["n_survcol_with_polypectomy"] -
                 screen_utility_loss["Any COL positive",]*QALY_table["CRC_sym_dx"] -
                 screen_utility_loss["Any COL negative",]*QALY_table["n_scrcol_without_polypectomy" ] -
                 screen_utility_loss["Any COL negative",]*QALY_table["n_fucol_without_polypectomy"] -
                 screen_utility_loss["Any COL negative",]*QALY_table["n_survcol_without_polypectomy" ] -
                 screen_utility_loss["Cardiovascular complication (COL)",]*QALY_table["cardio_compl"] -
                 screen_utility_loss["Serious GI complication (COL)",]*QALY_table["serious_GI_compl"] - 
                 screen_utility_loss["Other GI complication (COL)",]*QALY_table["other_GI_compl"]  -
                 (crc_care_utility_loss["Initial Phase","Stage_I"]*QALY_table["ly_crcI_initial"] + 
                    crc_care_utility_loss["Initial Phase","Stage_II"]*QALY_table["ly_crcII_initial"] +
                    crc_care_utility_loss["Initial Phase","Stage_III"]*QALY_table["ly_crcIII_initial"] +
                    crc_care_utility_loss["Initial Phase","Stage_IV"]*QALY_table["ly_crcIV_initial"] +
                    crc_care_utility_loss["Continuing Phase","Stage_I"]*QALY_table["ly_crcI_contin"] +
                    crc_care_utility_loss["Continuing Phase","Stage_II"]*QALY_table["ly_crcII_contin"] +
                    crc_care_utility_loss["Continuing Phase","Stage_III"]*QALY_table["ly_crcIII_contin"] +
                    crc_care_utility_loss["Continuing Phase","Stage_IV"]*QALY_table["ly_crcIV_contin"] +
                    crc_care_utility_loss["Terminal Phase (Death CRC)","Stage_I"]*QALY_table["ly_crcI_termcrc"] +
                    crc_care_utility_loss["Terminal Phase (Death CRC)","Stage_II"]*QALY_table["ly_crcII_termcrc"] +
                    crc_care_utility_loss["Terminal Phase (Death CRC)","Stage_III"]*QALY_table["ly_crcIII_termcrc"] +
                    crc_care_utility_loss["Terminal Phase (Death CRC)","Stage_IV"]*QALY_table["ly_crcIV_termcrc"] +
                    crc_care_utility_loss["Terminal Phase (Death Other)","Stage_I"]*QALY_table["ly_crcI_termoc"] + 
                    crc_care_utility_loss["Terminal Phase (Death Other)","Stage_II"]*QALY_table["ly_crcII_termoc"] + 
                    crc_care_utility_loss["Terminal Phase (Death Other)","Stage_III"]*QALY_table["ly_crcIII_termoc"] + 
                    crc_care_utility_loss["Terminal Phase (Death Other)","Stage_IV"]*QALY_table["ly_crcIV_termoc"] )) 
        if( strategy_label == "NoScreening" ){
          ## We compare other screening strategies to no screening so we want to save some values from the no screening
          ## strategy (colorectal cancer deaths, LY).
          assign(paste0("crc_deaths_no_screening_", discounting), get(paste0("crc_and_screen_deaths_", discounting)))
          assign(paste0("LY_no_screening_", discounting), totals_table["ly"])
          assign(paste0("QALY_no_screening_", discounting), get(paste0("QALY_", discounting)))
        } 
        ## COMPARE TO NO SCREENING ##
        ## Determine the LYs gained over no screening
        assign(paste0("LY_gained_over_NS_", discounting), get(paste0("LY_", discounting))-get(paste0("LY_no_screening_", discounting)))
        
        ## Determine the number of QALYs gained over no screening
        assign(paste0("QALY_gained_over_NS_", discounting), get(paste0("QALY_", discounting))- get(paste0("QALY_no_screening_", discounting)))
        
        ## For costs, we separate costs for those younger and those older than 65 years. If the row names (for crc care costs) or 
        ## column names (for screen costs) contain the pattern '(65+)' or (65Plus), respectively, we assume those values for the older population. 
        crc_care_costs_under65 <- crc_care_costs %>% rownames_to_column('phase') %>% slice(-grep("(65[+])",row.names(crc_care_costs))) %>% column_to_rownames('phase')
        crc_care_costs_65plus <- as.data.frame(crc_care_costs %>% rownames_to_column('phase') %>% slice(grep("(65[+])",row.names(crc_care_costs)))%>% column_to_rownames('phase'))
        screen_costs_under65 <- screen_costs %>% select(-grep("(65Plus)",colnames(screen_costs)))
        screen_costs_65plus <- screen_costs %>% select(grep("(65Plus)",colnames(screen_costs)))
        
        ## For calculating costs we separate the calculations for those younger than 65 and those 65 and older.
        age_groups <- c("under65", "65plus")
        if(first_age_of_interest >= 65){age_groups <- "65plus"}
        for(group in age_groups){
          
          ## Get the tables that are relevant for the age of interest 
          group_table <- get(paste0("totals_", group, "_", discounting))
          screen_costs_age_group <- get(paste0("screen_costs_", group))
          crc_care_costs_age_group <- get(paste0("crc_care_costs_", group))
          
          ## Calculate the costs by multiplying the cost times the number of events
          assign(paste0("fit_costs_", group, "_" ,discounting), screen_costs_age_group["FIT",] * sum(group_table[paste0(c("n_neg_", "n_pos_"),"fit")]))
          assign(paste0("sdna_costs_", group, "_", discounting), screen_costs_age_group["SDNA",] * sum(group_table[paste0(c("n_neg_", "n_pos_"),"sdna")]))
          assign(paste0("sig_costs_", group, "_", discounting), screen_costs_age_group["SIG",] * sum(group_table[paste0(c("n_neg_", "n_pos_"),"sig")]))
          assign(paste0("ctc_costs_", group, "_", discounting), screen_costs_age_group["CTC",] * sum(group_table[paste0(c("n_neg_", "n_pos_"),"ctc")]))
          assign(paste0("col_scr_costs_",group, "_" , discounting), screen_costs_age_group["Screening colonoscopy without polypectomy",]*group_table["n_scrcol_without_polypectomy"] + 
                   screen_costs_age_group["Colonoscopy with polypectomy",]*group_table["n_scrcol_with_polypectomy"] )
          assign(paste0("col_fu_costs_",group, "_" , discounting), screen_costs_age_group["Diagnostic colonoscopy without polypectomy",]*group_table["n_fucol_without_polypectomy"] + 
                   screen_costs_age_group["Colonoscopy with polypectomy",]*group_table["n_fucol_with_polypectomy"] )
          assign(paste0("col_surv_costs_", group, "_" ,discounting), screen_costs_age_group["Surveillance colonoscopy without polypectomy",]*group_table["n_survcol_without_polypectomy"] + 
                   screen_costs_age_group["Colonoscopy with polypectomy",]*group_table["n_survcol_with_polypectomy"] )
          assign(paste0("complication_costs_",group, "_" , discounting), screen_costs_age_group["Fatal perforation",] *group_table["n_compl_dth"] + 
                   screen_costs_age_group["Cardiovascular complication",] * group_table["cardio_compl"] +
                   screen_costs_age_group["Serious GI complication",] *group_table["serious_GI_compl"]  +
                   screen_costs_age_group["Other GI complication",] *group_table["other_GI_compl"] )
          
          if(output_template_year == 2020){
            assign(paste0("sensa_costs_", group, "_", discounting), screen_costs_age_group["Hemoccult SENSA",] * sum(group_table[paste0(c("n_neg_", "n_pos_"),"sensa")]))
            assign(paste0("screening_costs_",group, "_" , discounting), get(paste0("fit_costs_", group, "_" , discounting)) + get(paste0("sensa_costs_", group, "_" , discounting)) + 
                     get(paste0("sdna_costs_", group, "_" , discounting)) + get(paste0("sig_costs_", group, "_" , discounting)) + 
                     get(paste0("ctc_costs_", group, "_" , discounting)) + get(paste0("col_scr_costs_", group, "_" , discounting)))
            
          }else{
            assign(paste0("srna_costs_", group, "_", discounting), screen_costs_age_group["SRNA",] * sum(group_table[paste0(c("n_neg_", "n_pos_"),"srna")]))
            assign(paste0("blood_costs_", group, "_", discounting), screen_costs_age_group["BLOOD",] * sum(group_table[paste0(c("n_neg_", "n_pos_"),"blood")]))
            assign(paste0("screening_costs_",group, "_" , discounting), get(paste0("fit_costs_", group, "_" , discounting)) + get(paste0("srna_costs_", group, "_" , discounting)) + 
                     get(paste0("blood_costs_", group, "_" , discounting)) + 
                     get(paste0("sdna_costs_", group, "_" , discounting)) + get(paste0("sig_costs_", group, "_" , discounting)) + 
                     get(paste0("ctc_costs_", group, "_" , discounting)) + get(paste0("col_scr_costs_", group, "_" , discounting)))
            
          }
          
          ## For each cancer stage determine the cost. Costs are different by both stage and phase. We use pattern matching for getting the correct row (i.e. it ignores the '(65+)' for
          ## that age group)
          for(stage in c("I", "II", "III", "IV")){
            
            assign(paste0("crc_care_costs_stage_", stage, "_",group, "_" , discounting), crc_care_costs_age_group[grep("Initial Phase",row.names(crc_care_costs_age_group)),paste0("Stage_", stage)]*group_table[paste0("ly_crc", stage, "_initial")] +
                     crc_care_costs_age_group[grep("Continuing Phase",row.names(crc_care_costs_age_group)),paste0("Stage_", stage)]*group_table[paste0("ly_crc", stage, "_contin")] + 
                     crc_care_costs_age_group[grep("Death CRC",row.names(crc_care_costs_age_group)),paste0("Stage_", stage)]*group_table[paste0("ly_crc", stage, "_termcrc")] +
                     crc_care_costs_age_group[grep("Death Other",row.names(crc_care_costs_age_group)),paste0("Stage_", stage)]*group_table[paste0("ly_crc", stage, "_termoc")] +  
                     screen_costs_age_group[grep("Colonoscopy for diagnosis of a cancer by symptoms",row.names(screen_costs_age_group)),]*group_table[paste0("n_CRCstage", stage, "_sym")])
            
          } ## For crc care costs loop
          
          ## Get the costs of crc care and total costs for among the age group. 
          assign(paste0("crc_care_costs_all_stages_",group, "_" , discounting), sum(unlist(mget(paste0("crc_care_costs_stage_", c("I", "II", "III", "IV"), "_", group, "_" , discounting)))))
          
          assign(paste0("total_costs_", group, "_", discounting), get(paste0("screening_costs_",group, "_" , discounting)) + get(paste0("col_fu_costs_",group, "_" , discounting)) + get(paste0("col_surv_costs_",group, "_" , discounting)) +
                   get(paste0("complication_costs_",group, "_" , discounting)) + get(paste0("crc_care_costs_all_stages_",group, "_" , discounting)))
          
        }## For Age group-- cost loop
        
        ## Get the costs of crc care and total costs for the screening strategy (i.e. add the age groups together). mget allows for 'getting' of multiple
        ## objects and puts them in a list. We then need to unlist and sum them. This works better the adding two individual get functions for each age group
        ## because in the case that the 1st age of interest is greater than or equal to 65 we do not create values for less than 65 so trying to add them 
        ## together would result in an error. 
        assign(paste0("crc_care_costs_all_stages_", discounting), sum(unlist(mget(paste0("crc_care_costs_all_stages_", age_groups, "_", discounting)))))
        assign(paste0("total_costs_", discounting), sum(unlist(mget(paste0("total_costs_", age_groups, "_", discounting)))))
        
        ## Isolated Costs 
        assign(paste0("fit_costs_", discounting), sum(unlist(mget(paste0("fit_costs_", age_groups, "_", discounting)))))
        assign(paste0("sdna_costs_", discounting), sum(unlist(mget(paste0("sdna_costs_", age_groups, "_", discounting)))))
        assign(paste0("col_scr_costs_", discounting), sum(unlist(mget(paste0("col_scr_costs_", age_groups, "_", discounting)))))
        assign(paste0("col_fu_costs_", discounting), sum(unlist(mget(paste0("col_fu_costs_", age_groups, "_", discounting)))))
        assign(paste0("col_surv_costs_", discounting), sum(unlist(mget(paste0("col_surv_costs_", age_groups, "_", discounting)))))
        assign(paste0("sig_costs_", discounting), sum(unlist(mget(paste0("sig_costs_", age_groups, "_", discounting)))))
        assign(paste0("ctc_costs_", discounting), sum(unlist(mget(paste0("ctc_costs_", age_groups, "_", discounting)))))
        
        if(output_template_year == 2020){
          assign(paste0("sensa_costs_", discounting), sum(unlist(mget(paste0("sensa_costs_", age_groups, "_", discounting)))))
        }else{
          assign(paste0("srna_costs_", discounting), sum(unlist(mget(paste0("srna_costs_", age_groups, "_", discounting)))))
          assign(paste0("blood_costs_", discounting), sum(unlist(mget(paste0("blood_costs_", age_groups, "_", discounting)))))
        }
        
        ## Cost of complications
        assign(paste0("complication_costs_", discounting), sum(unlist(mget(paste0("complication_costs_", age_groups,"_", discounting)))))
        
      }## FOR DISCOUNTING LOOP!
      ## Determine the number of colorectal cancer cases avoided (subtract current strategy from no screening)
      crc_deaths_avoided_undsc <- crc_deaths_no_screening_undsc - crc_and_screen_deaths_undsc
      
      ## DETERMINE THE NUMBER OF CASES AVOIDED!
      crc_cases_avoided_undsc <- total_crc_cases_no_screening_undsc - total_crc_cases_undsc
      
      ## ERROR CHECKS ##
      
      ## Are there any screens BEFORE the age to begin screening
      if(output_template_year == 2020){
        if(strategy_label == "NoScreening"){
          if(sum(n_scrcol_infl_undsc, n_stool_tests_undsc, n_sig_tests_undsc, n_ctc_tests_undsc) > 0){
            print("There are screening tests for the 'No Screening' run")
            cat(paste("WARNING: There are screening tests for the 'No Screening' run\n"), file = paste0(folder_for_output, "/",time, "_ProcessUSPSTF_Log.txt"), append = TRUE)
          }
        }else{
          start_age <- min(as.numeric(run_info[c("stooltest_startage", "structuralexam_startage")]), na.rm = TRUE)
          totals_underStartAge_undsc <- colSums(strategy_undsc %>% filter(age < start_age))
          if(sum(totals_underStartAge_undsc[c("n_pos_fit", "n_neg_fit", 	"n_pos_sensa", "n_neg_sensa", "n_pos_sdna",  "n_neg_sdna",	
                                              "n_pos_sig",	"n_neg_sig",  "n_pos_ctc", "n_neg_ctc", "n_pos_scrcol", "n_neg_scrcol")]) > 0){
            print(paste("There are screens before screening is supposed to initiate for file", model_output_files[i]))
            cat(paste("WARNING: There are screens before screening is supposed to initiate for file", model_output_files[i], "\n"), file = paste0(folder_for_output, "/",time, "_ProcessUSPSTF_Log.txt"), append = TRUE )
          }
        }
      }else if(output_template_year == 2028){
        if(strategy_label == "NoScreening"){
          if(sum(n_scrcol_infl_undsc, n_stoolorblood_tests_undsc, n_sig_tests_undsc, n_ctc_tests_undsc) > 0){
            print("There are screening tests for the 'No Screening' run")
            cat(paste("WARNING: There are screening tests for the 'No Screening' run\n"), file = paste0(folder_for_output, "/",time, "_ProcessUSPSTF_Log.txt"), append = TRUE)
          }
        }else{
          start_age <- as.numeric(run_info[c("screening_startage")])
          totals_underStartAge_undsc <- colSums(strategy_undsc %>% filter(age < start_age))
          if(sum(totals_underStartAge_undsc[c("n_pos_fit", "n_neg_fit",  "n_pos_sdna",  "n_neg_sdna",	 "n_pos_srna",  "n_neg_srna", "n_pos_blood",  "n_neg_blood",
                                              "n_pos_sig",	"n_neg_sig",  "n_pos_ctc", "n_neg_ctc", "n_pos_scrcol", "n_neg_scrcol")]) > 0){
            print(paste("There are screens before screening is supposed to initiate for file", model_output_files[i]))
            cat(paste("WARNING: There are screens before screening is supposed to initiate for file", model_output_files[i], "\n"), file = paste0(folder_for_output, "/",time, "_ProcessUSPSTF_Log.txt"), append = TRUE )
          }
        }
      }
      
      
      ## The total number of stool tests should not exceed the number of tests for the strategy. (i.e. in a FIT strategy the number of FITs should be the same as the number 
      ## of stool tests.) 
      if(output_template_year == 2020){
        stool_test_lower <- str_to_lower(run_info[,"stooltest_type"])
        if(!is.na(stool_test_lower)){
          if(stool_test_lower == "sen"| stool_test_lower == "fobths" |stool_test_lower == "hsgfobt" | stool_test_lower == " sen"){
            stool_test_lower <- "sensa"
          }
          if(stool_test_lower == "fit-dna" | stool_test_lower == " sdna"| stool_test_lower == "sdna-fit"){
            stool_test_lower <- "sdna"
          }
        }
        if(!is.na(run_info[,"stooltest_type"])){
          if(get(paste0("n_",stool_test_lower,"_tests_", discounting)) != get(paste0("n_stool_tests_" ,discounting))){
            print(paste("The number of stool tests does not equal the number of ", run_info[,"stooltest_type"], " tests for the ", strategy_label, "strategy"))
            cat(paste("WARNING: The number of stool tests does not equal the number of ", run_info[,"stooltest_type"], " tests for the ", strategy_label, "strategy", "\n"), file = paste0(folder_for_output, "/",time, "_ProcessUSPSTF_Log.txt"), append = TRUE)
          }
        }
      }else if (output_template_year == 2028){
        stool_test_lower <- str_to_lower(run_info[,"stoolorbloodtest_type"])
        if(!is.na(stool_test_lower)){
          if(stool_test_lower == "fit-dna" | stool_test_lower == " sdna"| stool_test_lower == "sdna-fit"){
            stool_test_lower <- "sdna"
          }
        }
        if(!is.na(run_info[,"stoolorbloodtest_type"])){
          if(get(paste0("n_",stool_test_lower,"_tests_", discounting)) != get(paste0("n_stoolorblood_tests_" ,discounting))){
            print(paste("The number of stool tests does not equal the number of ", run_info[,"stoolorbloodtest_type"], " tests for the ", strategy_label, "strategy"))
            cat(paste("WARNING: The number of stool tests does not equal the number of ", run_info[,"stoolorbloodtest_type"], " tests for the ", strategy_label, "strategy", "\n"), file = paste0(folder_for_output, "/",time, "_ProcessUSPSTF_Log.txt"), append = TRUE)
          }
        }
      }
      
      ## Are all screening tests other than the test/test of interest 0? 
      if(output_template_year == 2020){
        screening_test_outcomes <- c("n_scrcol_infl_undsc", "n_sig_tests_undsc", "n_ctc_tests_undsc", "n_fit_tests_undsc", "n_sdna_tests_undsc", "n_sensa_tests_undsc")
      }else{
        screening_test_outcomes <- c("n_scrcol_infl_undsc", "n_sig_tests_undsc", "n_ctc_tests_undsc", "n_fit_tests_undsc", "n_sdna_tests_undsc", "n_srna_tests_undsc", "n_blood_tests_undsc")
      }
      
      if(strategy_label != "NoScreening"){
        tests_of_interest <- c(str_to_lower(run_info[,"structuralexam_type"]), stool_test_lower)
        tests_of_interest <- paste(tests_of_interest[!is.na(tests_of_interest)], collapse = "|")
        screening_test_outcomes <- screening_test_outcomes[!(grepl(tests_of_interest, screening_test_outcomes))]
        if(sum(unlist(mget(screening_test_outcomes)))> 0){
          print(paste("There are screening tests of the wrong modality occurring according to the run info for file",  model_output_files[i]))
          cat(paste("WARNING: There are screening tests of the wrong modality occurring according to the run info for file",  model_output_files[i], "\n"),file = paste0(folder_for_output, "/",time, "_ProcessUSPSTF_Log.txt"), append = TRUE )
        }
      }
      
      ## Do all positive (non-col) screening strategies have a diagnostic follow up exam. Check only if adherence is perfect. 
      if(output_template_year == 2020){
        if( run_info[,"adherence"] == "PERFECT" & round(sum(totals_undsc[paste0("n_pos_", c("fit", "sdna", "sensa", "sig", "ctc"))])- sum(totals_undsc[c("n_pos_fucol", "n_neg_fucol")])) > 0 ){
          print(paste("The number of positive screening tests is not equal to the number of diagnostic tests for file", model_output_files[i]))
          cat(paste("WARNING: The number of positive screening tests is not equal to the number of diagnostic tests for file", model_output_files[i], "\n"), file = paste0(folder_for_output, "/",time, "_ProcessUSPSTF_Log.txt"), append = TRUE)
        }
      }else{
        if( round(sum(totals_undsc[paste0("n_pos_", c("fit", "sdna", "srna","blood", "sig", "ctc"))])- sum(totals_undsc[c("n_pos_fucol", "n_neg_fucol")])) > 0 ){
          print(paste("The number of positive screening tests is not equal to the number of diagnostic tests for file", model_output_files[i]))
          cat(paste("WARNING: The number of positive screening tests is not equal to the number of diagnostic tests for file", model_output_files[i], "\n"), file = paste0(folder_for_output, "/",time, "_ProcessUSPSTF_Log.txt"), append = TRUE)
        }
      }
      
      ## If this is a COL only strategy is the number of FU cols 0? 
      if(output_template_year == 2020){
        if(!(is.na(run_info[,"structuralexam_type"])) & run_info[,"structuralexam_type"] == "COL" & is.na(run_info[,"stooltest_type"])){
          if(n_fucol_infl_undsc > 0){
            print(paste("There are 'follow-up' colonoscopies for a colonoscopy only strategy, file", model_output_files[i]))
            cat(paste("WARNING: There are 'follow-up' colonoscopies for a colonoscopy only strategy, file",model_output_files[i],"\n" ), file = paste0(folder_for_output, "/",time, "_ProcessUSPSTF_Log.txt"), append = TRUE)
          }
        }
      }else{
        if(!(is.na(run_info[,"structuralexam_type"])) & run_info[,"structuralexam_type"] == "COL" & is.na(run_info[,"stoolorbloodtest_type"])){
          if(n_fucol_infl_undsc > 0){
            print(paste("There are 'follow-up' colonoscopies for a colonoscopy only strategy, file", model_output_files[i]))
            cat(paste("WARNING: There are 'follow-up' colonoscopies for a colonoscopy only strategy, file",model_output_files[i],"\n" ), file = paste0(folder_for_output, "/",time, "_ProcessUSPSTF_Log.txt"), append = TRUE)
          }
        }
      }
      
      
      ## Is the number of surv cols after 85 go way down? 
      before_85 <- colSums(strategy_undsc %>% filter(age >= 79, age <=85))
      after_85 <- colSums(strategy_undsc %>% filter(age >= 86, age <=92))
      if(sum(before_85[c("n_pos_scrcol", "n_pos_survcol", "n_pos_fucol")]) < sum(after_85[c("n_pos_survcol", "n_neg_survcol")])){
        print(paste("There are too many surveillance colonoscopies for people ages 85 and older in file", model_output_files[i]))
        cat(paste("WARNING: There are too many surveillance colonoscopies for people ages 85 and older in file", model_output_files[i], "\n"), file = paste0(folder_for_output, "/",time, "_ProcessUSPSTF_Log.txt"), append = TRUE)
      }
      
      if(output_template_year == 2020){
        ## Does the number of adenomas on screening and screen detected cancers add up to the number of positive screening colonoscopies and positive diagnostic (FU) colonoscopies. 
        if(!are.equal(sum(totals_undsc[c("n_adn1to5_scr", "n_adn6to9_scr", "n_adn10plus_scr", "n_CRCstageI_scr", "n_CRCstageII_scr", "n_CRCstageIII_scr", "n_CRCstageIV_scr")]), sum(totals_undsc[c("n_pos_scrcol", "n_pos_fucol")]))){
          print(paste("The number of positive col screens and positive diagnostic tests is not equal to the number adenomas and screen detected cancers for file", model_output_files[i]))
          cat(paste("WARNING: The number of positive col screens and positive diagnostic tests is not equal to the number adenomas and screen detected cancers for file",model_output_files[i],"\n", sep = " "), file = paste0(folder_for_output, "/",time, "_ProcessUSPSTF_Log.txt"), append = TRUE)
        }
        ## Does the number of adenomas and surveillance detected cancers add up to the number of positive surveillance colonoscopies. 
        if(!are.equal(sum(totals_undsc[c("n_adn1to5_surv", "n_adn6to9_surv", "n_adn10plus_surv", "n_CRCstageI_surv", "n_CRCstageII_surv", "n_CRCstageIII_surv", "n_CRCstageIV_surv")]), sum(totals_undsc[c("n_pos_survcol")]))){
          print(paste("The number of positive surveillance cols is not equal to the number of adenomas and detected cancers found on surveillance for file", model_output_files[i]))
          cat(paste("WARNING: The number of positive surveillance cols is not equal to the number of adenomas and detected cancers found on surveillance for file", model_output_files[i],"\n"), file = paste0(folder_for_output, "/",time, "_ProcessUSPSTF_Log.txt"), append = TRUE)
        }
      }else{
        ## Does the number of lesions on screening and screen detected cancers add up to the number of positive screening colonoscopies and positive diagnostic (FU) colonoscopies. 
        if(!are.equal(sum(totals_undsc[c("n_lesion1to5_scr", "n_lesion6to9_scr", "n_lesion10plus_scr", "n_CRCstageI_scr", "n_CRCstageII_scr", "n_CRCstageIII_scr", "n_CRCstageIV_scr")]), sum(totals_undsc[c("n_pos_scrcol", "n_pos_fucol")]))){
          print(paste("The number of positive col screens and positive diagnostic tests is not equal to the number lesions and screen detected cancers for file", model_output_files[i]))
          cat(paste("WARNING: The number of positive col screens and positive diagnostic tests is not equal to the number lesion and screen detected cancers for file",model_output_files[i],"\n", sep = " "), file = paste0(folder_for_output, "/",time, "_ProcessUSPSTF_Log.txt"), append = TRUE)
        }
        ## Does the number of lesions and surveillance detected cancers add up to the number of positive surveillance colonoscopies. 
        if(!are.equal(sum(totals_undsc[c("n_lesion1to5_surv", "n_lesion6to9_surv", "n_lesion10plus_surv", "n_CRCstageI_surv", "n_CRCstageII_surv", "n_CRCstageIII_surv", "n_CRCstageIV_surv")]), sum(totals_undsc[c("n_pos_survcol")]))){
          print(paste("The number of positive surveillance cols is not equal to the number of lesions and detected cancers found on surveillance for file", model_output_files[i]))
          cat(paste("WARNING: The number of positive surveillance cols is not equal to the number of lesions and detected cancers found on surveillance for file", model_output_files[i],"\n"), file = paste0(folder_for_output, "/",time, "_ProcessUSPSTF_Log.txt"), append = TRUE)
        }
      }
      
      ## Does the LY with No CRC and ly in each of the stages and phases of cancer add up to the total LY. 
      if(!are.equal(round(as.numeric(sum(totals_undsc[c("ly_nocrc", "ly_crcI_initial" , "ly_crcII_initial", "ly_crcIII_initial", "ly_crcIV_initial", "ly_crcI_contin", "ly_crcII_contin", "ly_crcIII_contin", "ly_crcIV_contin",
                                                        "ly_crcI_termcrc", "ly_crcII_termcrc", "ly_crcIII_termcrc", "ly_crcIV_termcrc", "ly_crcI_termoc", "ly_crcII_termoc", "ly_crcIII_termoc", "ly_crcIV_termoc")])) - as.numeric(totals_undsc["ly"]),0),0)){
        print(paste("Life-years with CRC and life-years without CRC does not equal total life -years for file", model_output_files[i]))
        cat(paste("WARNING: Life-years with CRC and life-years without CRC does not equal total life -years for file", model_output_files[i],"\n"), file = paste0(folder_for_output, "/",time, "_ProcessUSPSTF_Log.txt"), append = TRUE)
      }
      ## Does the number alive at the start age equal the total number of people who died by final age. 
      if(!are.equal(strategy_undsc[1,"n_alive"], sum(crc_and_screen_deaths_undsc, totals_undsc["n_oc_dth"]))){
        print(paste("The number of deaths by age 100 is not equal to the number alive in at the 1st age of interest for file",model_output_files[i]))
        cat(paste("WARNING: The number of deaths by age 100 is not equal to the number alive in at the 1st age of interest for file",model_output_files[i],"\n"), file = paste0(folder_for_output, "/",time, "_ProcessUSPSTF_Log.txt"), append = TRUE)
      }
      
      ## Is the number alive without CRC equal to the total number of people alive at the 1st age of interest
      if(!are.equal(strategy_undsc[1,"n_alive"], strategy_undsc[1,"n_alive_nocrc"])){
        print(paste("The number of people alive in at the 1st age of interest is not the same as the number of people with no CRC for file",model_output_files[i]))
        cat(paste("WARNING: The number of people alive in at the 1st age of interest is not the same as the number of people with no CRC for file",model_output_files[i],"\n"), file = paste0(folder_for_output, "/",time, "_ProcessUSPSTF_Log.txt"), append = TRUE)
      }
      
      ## Is the average time in the initial and terminal phases in [0, 1) year? For each stage.
      ## Both bounds matter. A negative mean is as wrong as one over a year, and this
      ## check used to test only the upper one: in the broken screening arm the
      ## terminal-CRC column totalled -11,159 against a correct +3,904 and nothing
      ## fired, because every negative is < 1. na.rm covers a stage with no cases,
      ## where the mean is NaN and the old form errored instead of passing.
      v_mean_phaseI <- totals_undsc[c("ly_crcI_initial", "ly_crcI_termcrc",
                                    "ly_crcI_termoc")] /
        sum(totals_undsc[c("n_CRCstageI_scr", "n_CRCstageI_surv", "n_CRCstageI_sym")])
      if(any(v_mean_phaseI < 0 | v_mean_phaseI >= 1, na.rm = TRUE)){
        print(paste("Time in initial or terminal phase is incorrect for stage I"))
        cat(paste("WARNING: Time in initial or terminal phase is incorrect for stage I","\n"), file = paste0(folder_for_output, "/",time, "_ProcessUSPSTF_Log.txt"), append = TRUE)
      }
      v_mean_phaseII <- totals_undsc[c("ly_crcII_initial", "ly_crcII_termcrc",
                                    "ly_crcII_termoc")] /
        sum(totals_undsc[c("n_CRCstageII_scr", "n_CRCstageII_surv", "n_CRCstageII_sym")])
      if(any(v_mean_phaseII < 0 | v_mean_phaseII >= 1, na.rm = TRUE)){
        print(paste("Time in initial or terminal phase is incorrect for stage II"))
        cat(paste("WARNING: Time in initial or terminal phase is incorrect for stage II","\n"), file = paste0(folder_for_output, "/",time, "_ProcessUSPSTF_Log.txt"), append = TRUE)
      }
      v_mean_phaseIII <- totals_undsc[c("ly_crcIII_initial", "ly_crcIII_termcrc",
                                    "ly_crcIII_termoc")] /
        sum(totals_undsc[c("n_CRCstageIII_scr", "n_CRCstageIII_surv", "n_CRCstageIII_sym")])
      if(any(v_mean_phaseIII < 0 | v_mean_phaseIII >= 1, na.rm = TRUE)){
        print(paste("Time in initial or terminal phase is incorrect for stage III"))
        cat(paste("WARNING: Time in initial or terminal phase is incorrect for stage III","\n"), file = paste0(folder_for_output, "/",time, "_ProcessUSPSTF_Log.txt"), append = TRUE)
      }
      v_mean_phaseIV <- totals_undsc[c("ly_crcIV_initial", "ly_crcIV_termcrc",
                                    "ly_crcIV_termoc")] /
        sum(totals_undsc[c("n_CRCstageIV_scr", "n_CRCstageIV_surv", "n_CRCstageIV_sym")])
      if(any(v_mean_phaseIV < 0 | v_mean_phaseIV >= 1, na.rm = TRUE)){
        print(paste("Time in initial or terminal phase is incorrect for stage IV"))
        cat(paste("WARNING: Time in initial or terminal phase is incorrect for stage IV","\n"), file = paste0(folder_for_output, "/",time, "_ProcessUSPSTF_Log.txt"), append = TRUE)
      }
      
      
      ## Create the dataframe.
      if(output_template_year == 2020){
        model_data.df <-data.frame("Model" = model_x,
                                   "RiskScenario" = run_info[,"risk_scenario"],
                                   "StoolTestType"= run_info[,"stooltest_type"],
                                   "StoolStartAge" = run_info[,"stooltest_startage"], 
                                   "StoolStopAge" =run_info[,"stooltest_stopage"] ,
                                   "StoolInterval" = run_info[,"stooltest_interval"],
                                   "StructuralExamType" = run_info[,"structuralexam_type"],
                                   "StructuralStartAge" = run_info[,"structuralexam_startage"],
                                   "StructuralStopAge" = run_info[,"structuralexam_stopage"],
                                   "StructuralInterval" =run_info[,"structuralexam_interval"],
                                   "Adherence" = run_info[,"adherence"], 
                                   "Strategy" = strategy_label,
                                   "Colonoscopiesper1000" = n_col_tests_total_infl_undsc / n_alive_initial_no_CRC * 1000,
                                   #"ColonoscopiesperPerson" = n_col_tests_total_infl_undsc / n_alive_initial_no_CRC,
                                   "ScreeningColsper1000" = n_scrcol_infl_undsc / n_alive_initial_no_CRC * 1000, 
                                   "DiagnosticFUper1000" = n_fucol_infl_undsc / n_alive_initial_no_CRC * 1000, 
                                   "SurveillanceColsper1000" = n_survcol_infl_undsc / n_alive_initial_no_CRC * 1000, 
                                   "ColsForSymptDetectionper1000" = n_CRC_sym_dx_undsc / n_alive_initial_no_CRC * 1000 , # Number of cols for symptom detection is equal to the number of clinically detected CRCs 
                                   "NonScreeningCols" = (n_fucol_infl_undsc + n_survcol_infl_undsc + n_CRC_sym_dx_undsc) / n_alive_initial_no_CRC * 1000 , 
                                   "StoolTestsper1000"=  n_stool_tests_undsc/ n_alive_initial_no_CRC * 1000,
                                   "FITper1000"=  n_fit_tests_undsc/ n_alive_initial_no_CRC * 1000,
                                   "FITDNAper1000"=  n_sdna_tests_undsc/ n_alive_initial_no_CRC * 1000,
                                   "HSgFOBTper1000"=  n_sensa_tests_undsc / n_alive_initial_no_CRC * 1000,
                                   "SIGper1000"= n_sig_tests_undsc / n_alive_initial_no_CRC * 1000,
                                   "CTCper1000"= n_ctc_tests_undsc / n_alive_initial_no_CRC * 1000,
                                   "Complicationsper1000" = n_complications_undsc / n_alive_initial_no_CRC * 1000, 
                                   "CRCCasesper1000" = total_crc_cases_undsc / n_alive_initial_no_CRC * 1000 , 
                                   'DeathsFromColper1000' = screen_deaths_undsc/ n_alive_initial_no_CRC * 1000 , 
                                   "CRCandComplDeathsper1000"= crc_and_screen_deaths_undsc / n_alive_initial_no_CRC * 1000,
                                   "LYper1000" = LY_undsc / n_alive_initial_no_CRC * 1000,
                                   "LYGainedper1000"=  LY_gained_over_NS_undsc / n_alive_initial_no_CRC * 1000,  
                                   "DiscountedLYGainedper1000" =  LY_gained_over_NS_dsc/ n_alive_initial_no_CRC * 1000,
                                   "QALYper1000" = QALY_undsc /n_alive_initial_no_CRC * 1000,
                                   "QALYGainedper1000"=  QALY_gained_over_NS_undsc/ n_alive_initial_no_CRC * 1000,
                                   "DiscountedQALYGainedper1000"=  QALY_gained_over_NS_dsc / n_alive_initial_no_CRC * 1000,
                                   "LDGperPerson" = (LY_gained_over_NS_undsc*365.24) / n_alive_initial_no_CRC , 
                                   "CRCCasesAvertedper1000" = crc_cases_avoided_undsc / n_alive_initial_no_CRC * 1000, 
                                   "CRCDeathsAvertedper1000"= crc_deaths_avoided_undsc / n_alive_initial_no_CRC * 1000,
                                   "MortalityReduction" = crc_deaths_avoided_undsc / crc_deaths_no_screening_undsc,
                                   "FIT_costsper1000" = fit_costs_undsc/ n_alive_initial_no_CRC* 1000,
                                   "FITDNA_costsper1000" = sdna_costs_undsc/ n_alive_initial_no_CRC* 1000,
                                   "HSgFOBT_costsper1000" = sensa_costs_undsc/ n_alive_initial_no_CRC* 1000,
                                   "SIG_costsper1000" = sig_costs_undsc/ n_alive_initial_no_CRC* 1000,
                                   "CTC_costsper1000" = ctc_costs_undsc/ n_alive_initial_no_CRC* 1000,
                                   "COL_costsper1000" = (col_scr_costs_undsc + col_fu_costs_undsc + col_surv_costs_undsc)/ n_alive_initial_no_CRC* 1000,
                                   "Complication_costsper1000" = complication_costs_undsc / n_alive_initial_no_CRC* 1000,
                                   "CostsCancerCareper1000" = crc_care_costs_all_stages_undsc / n_alive_initial_no_CRC* 1000,
                                   "TotalCostsper1000"= total_costs_undsc / n_alive_initial_no_CRC* 1000,
                                   "DiscountedTotalCostsper1000"= total_costs_dsc / n_alive_initial_no_CRC* 1000,
                                   row.names = NULL
        )
      }else if (output_template_year == 2028){
        
        model_data.df <-data.frame("Model" = model_x,
                                   "RiskScenario" = run_info[,"risk_scenario"],
                                   "StoolTestType"= run_info[,"stoolorbloodtest_type"],
                                   "StoolStartAge" = run_info[,"screening_startage"], 
                                   "StoolStopAge" =run_info[,"screening_stopage"] ,
                                   "StoolInterval" = run_info[,"stoolorbloodtest_interval"],
                                   "StructuralExamType" = run_info[,"structuralexam_type"],
                                   "StructuralStartAge" = run_info[,"screening_startage"],
                                   "StructuralStopAge" = run_info[,"screening_stopage"],
                                   "StructuralInterval" =run_info[,"structuralexam_interval"],
                                   "AdherenceFollowUpCOL" = run_info[,"adherence_fucol"], 
                                   "Strategy" = strategy_label,
                                   "Colonoscopiesper1000" = n_col_tests_total_infl_undsc / n_alive_initial_no_CRC * 1000,
                                   #"ColonoscopiesperPerson" = n_col_tests_total_infl_undsc / n_alive_initial_no_CRC,
                                   "ScreeningColsper1000" = n_scrcol_infl_undsc / n_alive_initial_no_CRC * 1000, 
                                   "DiagnosticFUper1000" = n_fucol_infl_undsc / n_alive_initial_no_CRC * 1000, 
                                   "SurveillanceColsper1000" = n_survcol_infl_undsc / n_alive_initial_no_CRC * 1000, 
                                   "ColsForSymptDetectionper1000" = n_CRC_sym_dx_undsc / n_alive_initial_no_CRC * 1000 , # Number of cols for symptom detection is equal to the number of clinically detected CRCs 
                                   "NonScreeningCols" = (n_fucol_infl_undsc + n_survcol_infl_undsc + n_CRC_sym_dx_undsc) / n_alive_initial_no_CRC * 1000 , 
                                   "StoolorBloodTestsper1000"=  n_stoolorblood_tests_undsc/ n_alive_initial_no_CRC * 1000,
                                   "FITper1000"=  n_fit_tests_undsc/ n_alive_initial_no_CRC * 1000,
                                   "FITDNAper1000"=  n_sdna_tests_undsc/ n_alive_initial_no_CRC * 1000,
                                   "FITRNAper1000"=  n_srna_tests_undsc / n_alive_initial_no_CRC * 1000,
                                   "BLOODper1000"=  n_blood_tests_undsc / n_alive_initial_no_CRC * 1000,
                                   "SIGper1000"= n_sig_tests_undsc / n_alive_initial_no_CRC * 1000,
                                   "CTCper1000"= n_ctc_tests_undsc / n_alive_initial_no_CRC * 1000,
                                   "Complicationsper1000" = n_complications_undsc / n_alive_initial_no_CRC * 1000, 
                                   "CRCCasesper1000" = total_crc_cases_undsc / n_alive_initial_no_CRC * 1000 , 
                                   'DeathsFromColper1000' = screen_deaths_undsc/ n_alive_initial_no_CRC * 1000 , 
                                   "CRCandComplDeathsper1000"= crc_and_screen_deaths_undsc / n_alive_initial_no_CRC * 1000,
                                   "LYper1000" = LY_undsc / n_alive_initial_no_CRC * 1000,
                                   "LYGainedper1000"=  LY_gained_over_NS_undsc / n_alive_initial_no_CRC * 1000,  
                                   "DiscountedLYGainedper1000" =  LY_gained_over_NS_dsc/ n_alive_initial_no_CRC * 1000,
                                   "QALYper1000" = QALY_undsc /n_alive_initial_no_CRC * 1000,
                                   "QALYGainedper1000"=  QALY_gained_over_NS_undsc/ n_alive_initial_no_CRC * 1000,
                                   "DiscountedQALYGainedper1000"=  QALY_gained_over_NS_dsc / n_alive_initial_no_CRC * 1000,
                                   "LDGperPerson" = (LY_gained_over_NS_undsc*365.24) / n_alive_initial_no_CRC , 
                                   "CRCCasesAvertedper1000" = crc_cases_avoided_undsc / n_alive_initial_no_CRC * 1000, 
                                   "CRCDeathsAvertedper1000"= crc_deaths_avoided_undsc / n_alive_initial_no_CRC * 1000,
                                   "MortalityReduction" = crc_deaths_avoided_undsc / crc_deaths_no_screening_undsc,
                                   "FIT_costsper1000" = fit_costs_undsc/ n_alive_initial_no_CRC* 1000,
                                   "FITDNA_costsper1000" = sdna_costs_undsc/ n_alive_initial_no_CRC* 1000,
                                   "FITRNA_costsper1000" = srna_costs_undsc/ n_alive_initial_no_CRC* 1000,
                                   "BLOOD_costsper1000" = blood_costs_undsc/ n_alive_initial_no_CRC* 1000,
                                   "SIG_costsper1000" = sig_costs_undsc/ n_alive_initial_no_CRC* 1000,
                                   "CTC_costsper1000" = ctc_costs_undsc/ n_alive_initial_no_CRC* 1000,
                                   "COL_costsper1000" = (col_scr_costs_undsc + col_fu_costs_undsc + col_surv_costs_undsc)/ n_alive_initial_no_CRC* 1000,
                                   "Complication_costsper1000" = complication_costs_undsc / n_alive_initial_no_CRC* 1000,
                                   "CostsCancerCareper1000" = crc_care_costs_all_stages_undsc / n_alive_initial_no_CRC* 1000,
                                   "TotalCostsper1000"= total_costs_undsc / n_alive_initial_no_CRC* 1000,
                                   "DiscountedTotalCostsper1000"= total_costs_dsc / n_alive_initial_no_CRC* 1000,
                                   row.names = NULL
        )
      }
      ## Bind the temporary model_data data frame to the  other strategies 
      model_data <- rbind(model_data,  model_data.df)
    } #Loop with strategy 
    if(mismatch_population_counter >= 1){
      print(paste("Populations are not the same for model:", model_x, "variables such as LYGainedper1000, DiscountedLYGainedper1000, QALYGainedper1000 \n, 
                  DiscountedQALYGainedper1000, CRCCasesAvertedper1000, CRCDeathsAvertedper1000 and MortalityReduction are not reliable."))
      cat(paste("WARNING: Populations are not the same for model:", model_x, "variables such as LYGainedper1000, DiscountedLYGainedper1000, QALYGainedper1000,\n", 
                "DiscountedQALYGainedper1000, CRCCasesAvertedper1000, CRCDeathsAvertedper1000 and MortalityReduction are not reliable.\n"), file = paste0(folder_for_output, "/",time, "_ProcessUSPSTF_Log.txt"), append = TRUE) 
    }
  } #Loop with model_x
  
  ## Read the csv with the list of outcomes that you want and then select the appropriate columns. 
  desired_model_data_outcomes <- as.data.frame(read.csv(paste0(input_folder,"/", selected_outcomes_for_model_data_file), colClasses = c("character", "logical"))) %>% 
    filter(Boolean == TRUE)
  
  # model_data <- model_data[,(names(model_data) %in% c("Model", "RiskScenario", "StoolTestType","StoolStartAge", "StoolStopAge" , "StoolInterval",
  #                                                     "StructuralExamType", "StructuralStartAge", "StructuralStopAge", "StructuralInterval","Adherence", "Strategy", 
  #                                                     unlist(desired_model_data_outcomes["Outcome"])))]
  # model_data <- model_data[,c("Model", "RiskScenario", "StoolTestType","StoolStartAge", "StoolStopAge" , "StoolInterval",
  #                             "StructuralExamType", "StructuralStartAge", "StructuralStopAge", "StructuralInterval","Adherence", "Strategy", 
  #                             unlist(desired_model_data_outcomes["Outcome"]))]
  ## Get the arguments of the function. This will go on a second sheet in the excel file that is written. This was added to help
  ## keep track of files and when looking back we will have a better sense as to what the results in the file actually pertain to. 
  arguments <- data.frame("date" = Sys.time(),
                          "analysis_folder" = analysis_folder, 
                          "input_folder" = input_folder,  
                          "folder_for_output" = folder_for_output,
                          "include_SimCRC" = include_SimCRC,
                          "include_SimCRC_R" = include_SimCRC_R, 
                          "include_MISCAN" = include_MISCAN, 
                          "include_CRCSPIN" = include_CRCSPIN, 
                          "first_age_of_interest" = first_age_of_interest,
                          "input_prefix" = input_prefix, 
                          "col_infl_rate" = col_infl_rate,  
                          "discount_rate" = discount_rate,
                          "col_spec_adj" = col_spec_adj, 
                          "crc_care_costs_file" = crc_care_costs_file,
                          "screen_costs_file" = screen_costs_file,
                          "crc_care_disutility_file" = crc_care_disutility_file,
                          "screen_disutility_file" = screen_disutility_file,
                          "general_health_utility_weights_file" = general_health_utility_weights_file, 
                          "model_run_data_tag" = model_run_data_tag, 
                          "ProcessUSPSTFOutput_script_location" = ProcessUSPSTFOutput_script_location)
  arguments <- as.data.frame(t(arguments))
  # Add row names as a separate column to 'arguments'
  arguments$row_names <- row.names(arguments)
  # Set the new row names column as the first column
  file_name <- paste0(folder_for_output,"/",time, "_model_run_data",model_run_data_tag,".xlsx")
  
  # Create a list of data frames for each sheet in the Excel file
  sheets <- list(
    model_run_data = model_data,  # Data for 'model_run_data' sheet
    info = arguments # Data for 'info' sheet, with row names included
  )
  writexl::write_xlsx(sheets, path = file_name)
  
  ## Write the file! 
  #write.xlsx(model_data, file_name, sheetName = "model_run_data", row.names = FALSE)
  #write.xlsx(arguments, file_name, sheetName = "info", row.names = TRUE, col.names = FALSE, append =  TRUE)
  # write.xlsx(full_QALY_table, paste0(folder_for_output,"/", "QALYSums.xlsx"))
  ## Tell us that the model_data data file has been written by printing to the screen. 
  print("Done writing 'model_run_data'")
  cat("Done writing 'model_run_data'", file = paste0(folder_for_output, "/",time, "_ProcessUSPSTF_Log.txt"), append = TRUE)
  return(model_data)
}
