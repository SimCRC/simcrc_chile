
true_target_simcrc <- read.csv(param_BayCANN$targets_file)
true_target_simcrc$sd <- (true_target_simcrc$stopping_upper_bounds - true_target_simcrc$stopping_lower_bounds)/(2*1.96)


# 4. Select parameter set based on the max lp value ---------------------------

type_param_set <- "Max_lp"

max_lp <- max(calibrated_params$lp)

#get parameters set which maximize the lp
l_params_opt <- calibrated_params[which.max(calibrated_params$lp),]
# drop the variable lp
l_params_opt <- l_params_opt[,1:dim(l_params_opt)[2] - 1]

l_params_calibrated <- as.list(l_params_opt)

l_params_calibrated_Max_lp <- l_params_calibrated

#path_calibrated_set <- paste0(folder,"/l_params_calibrated_",type_param_set,"_",BayCANN_version,".rda")

#save(l_params_calibrated, file = path_calibrated_set)

if(!SSP_cal){ 
  #Adenoma
  output_simCRC_max_lp <- calibration_out(l_params_all = l_params_init,v_params_calib = l_params_calibrated, dt_pop, dt_long = TRUE)
} else {
  #SSP
  output_simCRC_max_lp <- calibration_out_ssp(l_params_all = l_params_init,v_params_calib = l_params_calibrated, dt_pop, dt_long = TRUE)
}
colnames(output_simCRC_max_lp) <- c("target_names","output_model")
output_simCRC_max_lp$type_param_set <- type_param_set
output_simCRC_max_lp <- inner_join(output_simCRC_max_lp , true_target_simcrc ,by=c("target_names"))  #Adding group


# 5. Select parameter set based on the log-likelihood maximization --------

type_param_set <- "Max_log_likelihood"

#Log-likelihood function
log_likelihood_single <- function(x, target_mean, target_sd) {
  # Check if input parameters are in the correct format
  if (!is.numeric(x) || !is.numeric(target_mean) || !is.numeric(target_sd)) {
    stop("Observed value, target mean, and target standard deviation should be numeric")
  }
  
  # Ensure that the standard deviation is positive
  if (target_sd <= 0) {
    stop("Standard deviation must be positive")
  }
  
  # Calculate the log-likelihood for the observation
  log_likelihood <- -0.5 * log(2 * pi * target_sd^2) -
    0.5 * ((x - target_mean)^2 / target_sd^2)
  
  return(as.numeric(log_likelihood))
}


v_targets_names <- param_BayCANN$outputs_names

true_target_simcrc <- read.csv(param_BayCANN$targets_file)
#true_target_simcrc <- read.csv("data-raw/20220909_simcrc_targets_ssp.csv")   #Reading manually

true_target_simcrc$sd <- (true_target_simcrc$stopping_upper_bounds - true_target_simcrc$stopping_lower_bounds)/(2*1.96)

#Model outputs included in the target file
true_targets_names <- colnames(df_simcrc_outputs)[colnames(df_simcrc_outputs) %in% true_target_simcrc$target_names]

n_targets <- length(true_targets_names)  #-1 to discount de id row
n_outputs <- dim(df_simcrc_outputs)[1]

m_ll_outputs <- matrix(data = NA, nrow = n_outputs , ncol = n_targets )
cont=1
for (target in true_targets_names) {
  y_true_mean = true_target_simcrc[true_target_simcrc$target_names==target,]$targets
  y_true_sd   = true_target_simcrc[true_target_simcrc$target_names==target,]$sd
  y_pred      = df_simcrc_outputs[,target]
  
  m_ll_outputs[,cont] = log_likelihood_single(x = y_pred, target_mean = y_true_mean, target_sd = y_true_sd)
  cont = cont + 1
}

df_ll_outputs <- as.data.frame(m_ll_outputs)
df_ll_outputs$sum <- rowSums(df_ll_outputs)

index_max_ll <- which.max(df_ll_outputs$sum)

l_params_opt <- calibrated_params[index_max_ll,]

l_params_calibrated <- as.list(l_params_opt)

l_params_calibrated_Max_log_likelihood <- l_params_calibrated

# path_calibrated_set <- paste0(folder,"/l_params_calibrated_",type_param_set,"_",BayCANN_version,".rda")
# 
# save(l_params_calibrated, file = path_calibrated_set)

if(!SSP_cal){ 
  
  #Adenoma
  output_simCRC_Max_log_likelihood <- calibration_out(l_params_all = l_params_init,v_params_calib = l_params_calibrated, dt_pop, dt_long = TRUE)
} else {
  
  #SSP
  output_simCRC_Max_log_likelihood <- calibration_out_ssp(l_params_all = l_params_init,v_params_calib = l_params_calibrated, dt_pop, dt_long = TRUE)
  
}
colnames(output_simCRC_Max_log_likelihood) <- c("target_names","output_model")
output_simCRC_Max_log_likelihood$type_param_set <- type_param_set
output_simCRC_Max_log_likelihood <- inner_join(output_simCRC_Max_log_likelihood , true_target_simcrc ,by=c("target_names"))  #Adding group

# 6. Select parameter set based on the Minimum Absolute Error --------

type_param_set <- "Min_AbsolutErr"

v_targets_names <- param_BayCANN$outputs_names

true_target_simcrc <- read.csv(param_BayCANN$targets_file)
#true_target_simcrc <- read.csv("data-raw/20220909_simcrc_targets_ssp.csv")   #Reading manually

true_target_simcrc$sd <- (true_target_simcrc$stopping_upper_bounds - true_target_simcrc$stopping_lower_bounds)/(2*1.96)

#Model outputs included in the target file
true_targets_names <- colnames(df_simcrc_outputs)[colnames(df_simcrc_outputs) %in% true_target_simcrc$target_names]

n_targets <- length(true_targets_names)  #-1 to discount de id row
n_outputs <- dim(df_simcrc_outputs)[1]


# Function to calculate the Standardized Mean Squared Error (SMSE)
smse <- function(y_true_mean, y_true_sd, y_pred) {

  # Calculate the Standardized Mean Squared Error (SMSE)
  sse_value <- abs(y_true_mean - y_pred) / y_true_sd
  
  return(sse_value)
}


m_smse_outputs <- matrix(data = NA, nrow = n_outputs , ncol = n_targets )
cont=1
for (target in true_targets_names) {
  
  y_true_mean = true_target_simcrc[true_target_simcrc$target_names==target,]$targets
  y_true_sd   = true_target_simcrc[true_target_simcrc$target_names==target,]$sd
  y_pred      = df_simcrc_outputs[,target]
  
  m_smse_outputs[,cont] = smse(y_true_mean = y_true_mean, y_true_sd = y_true_sd, y_pred = y_pred)
  cont = cont + 1
}

df_smse_outputs <- as.data.frame(m_smse_outputs)
df_smse_outputs$sum <- rowSums(df_smse_outputs)

index_min_ae <- which.min(df_smse_outputs$sum)

l_params_opt <- calibrated_params[index_min_ae,1:dim(calibrated_params)[2] - 1]

l_params_calibrated <- as.list(l_params_opt)

l_params_calibrated_Min_AbsolutErr <- l_params_calibrated

# BayCANN_version <- param_BayCANN$BayCANN_version
# 
# path_calibrated_set <- paste0(folder,"/l_params_calibrated_",type_param_set,"_",BayCANN_version,".rda")
# 
# save(l_params_calibrated, file = path_calibrated_set)

if(!SSP_cal){ 
  
  #Adenoma
  output_simCRC_Min_AbsolutErr <- calibration_out(l_params_all = l_params_init,v_params_calib = l_params_calibrated, dt_pop, dt_long = TRUE)
} else {
  
  #SSP
  output_simCRC_Min_AbsolutErr <- calibration_out_ssp(l_params_all = l_params_init,v_params_calib = l_params_calibrated, dt_pop, dt_long = TRUE)
  
}
colnames(output_simCRC_Min_AbsolutErr) <- c("target_names","output_model")
output_simCRC_Min_AbsolutErr$type_param_set <- type_param_set
output_simCRC_Min_AbsolutErr <- inner_join(output_simCRC_Min_AbsolutErr , true_target_simcrc ,by=c("target_names"))  #Adding group

# 7. Select parameter set based on the minimum of mean squared error --------

type_param_set <- "Min_MSE"

v_targets_names <- param_BayCANN$outputs_names

true_target_simcrc <- read.csv(param_BayCANN$targets_file)
#true_target_simcrc <- read.csv("data-raw/20220909_simcrc_targets_ssp.csv")   #Reading manually

true_target_simcrc$sd <- (true_target_simcrc$stopping_upper_bounds - true_target_simcrc$stopping_lower_bounds)/(2*1.96)

#Model outputs included in the target file
true_targets_names <- colnames(df_simcrc_outputs)[colnames(df_simcrc_outputs) %in% true_target_simcrc$target_names]

n_targets <- length(true_targets_names)  #-1 to discount de id row
n_outputs <- dim(df_simcrc_outputs)[1]

# Function to calculate the  Squared Error
get_mse <- function(y_true_mean, y_true_sd, y_pred) {
  
  # Calculate the Standardized Mean Squared Error (SMSE)
  mse_value <- (y_true_mean - y_pred)^2
  
  return(mse_value)
}

m_smse_outputs <- matrix(data = NA, nrow = n_outputs , ncol = n_targets )
cont=1
for (target in true_targets_names) {
  
  y_true_mean = true_target_simcrc[true_target_simcrc$target_names==target,]$targets
  y_true_sd   = true_target_simcrc[true_target_simcrc$target_names==target,]$sd
  y_pred      = df_simcrc_outputs[,target]
  
  m_smse_outputs[,cont] = smse(y_true_mean = y_true_mean, y_true_sd = y_true_sd, y_pred = y_pred)
  cont = cont + 1
}

df_smse_outputs <- as.data.frame(m_smse_outputs)
df_smse_outputs$mean <- rowMeans(df_smse_outputs)

index_min_mse <- which.min(df_smse_outputs$mean)

l_params_opt <- calibrated_params[index_min_mse,]

l_params_calibrated <- as.list(l_params_opt)

l_params_calibrated_Min_MSE <- l_params_calibrated

# BayCANN_version <- param_BayCANN$BayCANN_version
# 
# file_calibrated_params <- paste0(folder,"/l_params_calibrated_",type_param_set,"_",BayCANN_version,".rda")
# 
# save(l_params_calibrated, file = file_calibrated_params)

if(!SSP_cal){ 
  
  #Adenoma
  output_simCRC_Min_MSE <- calibration_out(l_params_all = l_params_init,v_params_calib = l_params_calibrated, dt_pop, dt_long = TRUE)
} else {
  
  #SSP
  output_simCRC_Min_MSE <- calibration_out_ssp(l_params_all = l_params_init,v_params_calib = l_params_calibrated, dt_pop, dt_long = TRUE)
  
}
colnames(output_simCRC_Min_MSE) <- c("target_names","output_model")
output_simCRC_Min_MSE$type_param_set <- type_param_set
output_simCRC_Min_MSE <- inner_join(output_simCRC_Min_MSE , true_target_simcrc ,by=c("target_names"))  #Adding group


# 8. select parameters based on the parameter means -------------------------------------------------

type_param_set <- "Post_mean"

l_params_opt <- colMeans(calibrated_params[,1:dim(calibrated_params)[2]-1], na.rm = TRUE)
#l_params_opt <- colMeans(df_posterior_param[1:1500,1:dim(df_posterior_param)[2]-1], na.rm = TRUE)
l_params_calibrated <- as.list(l_params_opt) 

l_params_calibrated_Post_mean <- l_params_calibrated

# file_calibrated_params <- paste0(folder,"/l_params_calibrated_",type_param_set,"_",BayCANN_version,".rda")
# 
# save(l_params_calibrated, file = file_calibrated_params)

if(!SSP_cal){ 
  #Adenoma
  output_simCRC_Post_mean <- calibration_out(l_params_all = l_params_init,v_params_calib = l_params_calibrated, dt_pop, dt_long = TRUE)
} else {
  #SSP
  output_simCRC_Post_mean <- calibration_out_ssp(l_params_all = l_params_init,v_params_calib = l_params_calibrated, dt_pop, dt_long = TRUE)
}
colnames(output_simCRC_Post_mean) <- c("target_names","output_model")
output_simCRC_Post_mean$type_param_set <- type_param_set
output_simCRC_Post_mean <- inner_join(output_simCRC_Post_mean , true_target_simcrc ,by=c("target_names"))  #Adding group


# 8. select parameters based on the parameter medians -------------------------------------------------

type_param_set <- "Post_median"



l_params_opt <- matrixStats::colMedians(as.matrix(calibrated_params[,1:dim(calibrated_params)[2]-1]), na.rm = TRUE)
l_params_calibrated <- as.list(l_params_opt) 

l_params_calibrated_Post_median <- l_params_calibrated

# file_calibrated_params <- paste0(folder,"/l_params_calibrated_",type_param_set,"_",BayCANN_version,".rda")
# 
# save(l_params_calibrated, file = file_calibrated_params)

if(!SSP_cal){ 
  
  #Adenoma
  output_simCRC_Post_median <- calibration_out(l_params_all = l_params_init,v_params_calib = l_params_calibrated, dt_pop, dt_long = TRUE)
} else {
  
  #SSP
  output_simCRC_Post_median <- calibration_out_ssp(l_params_all = l_params_init,v_params_calib = l_params_calibrated, dt_pop, dt_long = TRUE)
  
}
colnames(output_simCRC_Post_median) <- c("target_names","output_model")
output_simCRC_Post_median$type_param_set <- type_param_set
output_simCRC_Post_median <- inner_join(output_simCRC_Post_median , true_target_simcrc ,by=c("target_names"))  #Adding group


#####
v_type_param_set <- c("Max_lp", "Max_log_likelihood", "Min_AbsolutErr", "Min_MSE", "Post_mean", "Post_median")


l_params_calibrated_sets <- list()


for (type_param_set in v_type_param_set) {
  param_set_name <- paste0("l_params_calibrated_", type_param_set)
  l_params_calibrated_sets[[type_param_set]] <- get(param_set_name)
}


save(l_params_calibrated_sets, file = paths_calibration$path_best_params_sets)

# Append all outputs from different parameters set ------------------------

# Create a data frame to hold all selected sets

df_outputs_selected_sets <- rbind(     output_simCRC_max_lp, 
                                       output_simCRC_Max_log_likelihood, 
                                       output_simCRC_Min_AbsolutErr, 
                                       output_simCRC_Min_MSE,
                                       output_simCRC_Post_mean) 
                                       #output_simCRC_Post_median)


#----------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------#
# 10. Check the parameter set selected against the targets -------------------------------------------------


df_outputs_selected_sets$lesion_type <- factor(df_outputs_selected_sets$lesion_type, c("Adenoma", "SSP", "All"))

df_outputs_selected_sets[ df_outputs_selected_sets$target_groups=="Prev0Ad" , ]$target_groups  <- "Prevalence 0 lesions"
df_outputs_selected_sets[ df_outputs_selected_sets$target_groups=="Prev1Ad" , ]$target_groups  <- "Prevalence 1 lesion"
df_outputs_selected_sets[ df_outputs_selected_sets$target_groups=="Prev2Ad" , ]$target_groups  <- "Prevalence 2 lesions"
df_outputs_selected_sets[ df_outputs_selected_sets$target_groups=="Prev3PlusAd" , ]$target_groups  <- "Prevalance 3+ lesions"


#Set the order we want to have on the facets
df_outputs_selected_sets$target_groups<- factor(df_outputs_selected_sets$target_groups, c("CRCInc_Overall",
                                                                                                     "CRCInc_P", 
                                                                                                     "CRCInc_D",
                                                                                                     "CRCInc_R",
                                                                                                     "Prevalence 0 lesions",
                                                                                                     "Prevalence 1 lesion",
                                                                                                     "Prevalence 2 lesions",
                                                                                                     "Prevalance 3+ lesions",
                                                                                                     "StageDist_P",
                                                                                                     "StageDist_D",
                                                                                                     "StageDist_R",
                                                                                                     "PrevPreclin",
                                                                                                     "Size_P",
                                                                                                     "Size_D",
                                                                                                     "Size_R",
                                                                                                     "SSP"
                                                                                                     ))
color_values = c("Adenoma" = "#56B4E9", "SSP" = "#E69F00", "All" = "#009E73")

color_sets <- c("Max_lp" = "#D55E00", 
                 "Max_log_likelihood" = "#0072B2", 
                 "Min_AbsolutErr" = "#009E73", 
                 "Min_MSE" = "#CC79A7",
                 "Post_mean" = "#F0E442",
                 "Post_median" = "#999999")

#identification of categorical variables

cat_groups <- c("Size_P", "Size_D", "Size_R", "StageDist_D","StageDist_P","StageDist_R", "SSP")
df_outputs_selected_sets$categorical <- ifelse(df_outputs_selected_sets$target_groups %in% cat_groups,1,0)


plot_val_num <- ggplot(data = df_outputs_selected_sets[df_outputs_selected_sets$categorical==0,], 
                aes(x    = age, 
                    y    = targets, 
                    ymin = stopping_lower_bounds, 
                    ymax = stopping_upper_bounds
                    ))+ 
  geom_errorbar(width=1.2, size=0.5, color="red") +
  geom_ribbon(data = df_outputs_selected_sets[df_outputs_selected_sets$categorical==0,],
              aes(x    = age,
                  y    = output_model,
                  ymin = output_model,
                  ymax = output_model,
                  color=type_param_set),
              alpha = 0.3) +
  facet_wrap(~ target_groups + lesion_type,scales="free", ncol = 4) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(), legend.position="none") +
  #scale_color_manual(values = color_values)+
  # scale_color_manual(
  #   name = "Set selection", # Legend title
  #   values = c("#D55E00",
  #              "#0072B2",
  #              "#009E73",
  #              "#CC79A7",
  #              "#F0E442",
  #              "#999999") # Colors manually assigned
  #   #labels = c("Max_lp", "Max_log_likelihood","Min_AbsolutErr", "Min_MSE","Post_mean","Post_median") # Labels for legend
  # )+
  #scale_fill_manual(values =color_values)+
  #scale_y_continuous(breaks= 0:.2, labels=c(0,1))+
  scale_y_continuous(breaks = number_ticks(5))+
  #scale_x_continuous(breaks= 1:lim, labels=x_label) +
  #add legend for the ribbon colors
  theme_bw(base_size = 8) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 8, angle = 90),
        axis.title = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0)) +
  labs(title = "SimCRC-R calibrated parameter set validation ", 
       x     = "")

plot_val_num

ggsave(plot_val_num,
       filename = paths_calibration$path_set_param_val_num,
       width = 10, height = 6)


# First, create numeric positions for your groups
df_outputs_selected_sets$x_pos <- as.numeric(factor(df_outputs_selected_sets$stage))
shift <- 0.3

plot_val_cat <- ggplot(data = df_outputs_selected_sets[df_outputs_selected_sets$categorical==1,], 
                       aes(x    = x_pos - shift, 
                           y    = targets, 
                           ymin = stopping_lower_bounds, 
                           ymax = stopping_upper_bounds,
                           shape = "target"))+ 
  geom_point()+
  geom_errorbar(width=.2, size=1, alpha=0.5, color="red") +
  theme(legend.position="none") +
  geom_errorbar(data = df_outputs_selected_sets[df_outputs_selected_sets$categorical==1,],
                aes(x    = x_pos + shift,
                    y    = output_model,
                    ymin = output_model,
                    ymax = output_model, color= type_param_set), width=.2, size=0.9, alpha = 1) +
  geom_point(aes(x    = x_pos + shift,
                 y    = output_model , shape="model"))+
  facet_wrap(~ target_groups + lesion_type,scales="free", ncol = 3) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(), legend.position="none") +
  # scale_fill_manual(values = color_values)+
  # scale_color_manual(values = color_values)+
  scale_shape_manual(values = c(8 , 16)) +
  #scale_y_continuous(breaks= 0:.2, labels=c(0,1))+
  scale_y_continuous(breaks = number_ticks(5))+
  # Adjust x-axis labels back to original group labels
  scale_x_continuous(breaks = df_outputs_selected_sets$x_pos, labels = df_outputs_selected_sets$stage) +
  theme_bw(base_size = 8) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 8, angle = 0),
        axis.title = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0)) +
  labs(title = "SimCRC-R calibrated parameter set validation ", 
       x     = "")

plot_val_cat

ggsave(plot_val_cat,
       filename = paths_calibration$path_set_param_val_cat,
       width = 10, height = 6)



