


true_target_simcrc <- read.csv(param_BayCANN$targets_file)

true_target_simcrc$sd <- (true_target_simcrc$stopping_upper_bounds - true_target_simcrc$stopping_lower_bounds)/(2*1.96)


calibrated_params <- read.csv(param_BayCANN$path_posterior)

# --- Guard against metadata-leak contamination in the posteriors file ----------
# Pre-fix calibration runs could leak non-parameter columns into the posteriors CSV
# (a spurious duplicate 'lp__.1' alongside the real log-posterior 'lp__'). Drop the
# duplicate so `$lp` resolves unambiguously, and move 'lp__'/'chain' to the end so the
# downstream `[, 1:ncol-1]` parameter extraction never drops a real parameter.
# No-op on clean files, so it is safe to keep for future runs.
if ("lp__.1" %in% colnames(calibrated_params)) {
  print("Dropped spurious 'lp__.1' column from posteriors (metadata-leak artifact)")
  calibrated_params$lp__.1 <- NULL
}
v_tail_cols  <- intersect(c("lp__", "chain"), colnames(calibrated_params))
v_other_cols <- setdiff(colnames(calibrated_params), v_tail_cols)
calibrated_params <- calibrated_params[, c(v_other_cols, v_tail_cols), drop = FALSE]

# Select chains to include (all four chains)
selected_chains <- c(1, 2, 3, 4)
calibrated_params <- calibrated_params[calibrated_params$chain %in% selected_chains, ]
# Keep df_simcrc_outputs aligned with calibrated_params: the best-set selectors below
# (which.max / which.min) return ROW POSITIONS into df_simcrc_outputs, then index
# calibrated_params. If only calibrated_params is chain-filtered, those positions point
# at dropped/nonexistent rows and yield an all-NA parameter set. Both come from the same
# posteriors file in the same row order, so filtering both by the same mask keeps them aligned.
if (exists("df_simcrc_outputs") && "chain" %in% colnames(df_simcrc_outputs)) {
  df_simcrc_outputs <- df_simcrc_outputs[df_simcrc_outputs$chain %in% selected_chains, ]
}

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
  output_simCRC_max_lp <- calibration_out(l_params_all = l_params_init, v_params_calib = l_params_calibrated, id_draw = l_params_calibrated$id_draw, dt_pop, dt_long = TRUE)
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
  output_simCRC_Max_log_likelihood <- calibration_out(l_params_all = l_params_init,v_params_calib = l_params_calibrated,  id_draw = l_params_calibrated$id_draw, dt_pop, dt_long = TRUE)
} else {
  
  #SSP
  output_simCRC_Max_log_likelihood <- calibration_out_ssp(l_params_all = l_params_init,v_params_calib = l_params_calibrated, dt_pop, dt_long = TRUE)
  
}
colnames(output_simCRC_Max_log_likelihood) <- c("target_names","output_model")
output_simCRC_Max_log_likelihood$type_param_set <- type_param_set
output_simCRC_Max_log_likelihood <- inner_join(output_simCRC_Max_log_likelihood , true_target_simcrc ,by=c("target_names"))  #Adding group

# 6. Select parameter set based on the Minimum Absolute Error --------

type_param_set <- "Min_AbsolutErr"

v_targets_names <- true_target_simcrc$target_names


true_target_simcrc <- read.csv(param_BayCANN$targets_file)

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


m_smse_outputs <- matrix(data = NA, nrow = n_outputs , ncol = length(true_targets_names))
for (j in seq_along(true_targets_names)) {
  target      <- true_targets_names[j]
  y_true_mean <- true_target_simcrc[true_target_simcrc$target_names == target, ]$targets
  y_true_sd   <- true_target_simcrc[true_target_simcrc$target_names == target, ]$sd
  y_pred      <- df_simcrc_outputs[[target]]

  m_smse_outputs[, j] <- smse(y_true_mean = y_true_mean, y_true_sd = y_true_sd, y_pred = y_pred)
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
  output_simCRC_Min_AbsolutErr <- calibration_out(l_params_all = l_params_init,v_params_calib = l_params_calibrated, id_draw = l_params_calibrated$id_draw, dt_pop, dt_long = TRUE)
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

m_smse_outputs <- matrix(data = NA, nrow = n_outputs , ncol = length(true_targets_names))
for (j in seq_along(true_targets_names)) {
  target      <- true_targets_names[j]
  y_true_mean <- true_target_simcrc[true_target_simcrc$target_names == target, ]$targets
  y_true_sd   <- true_target_simcrc[true_target_simcrc$target_names == target, ]$sd
  y_pred      <- df_simcrc_outputs[[target]]

  m_smse_outputs[, j] <- smse(y_true_mean = y_true_mean, y_true_sd = y_true_sd, y_pred = y_pred)
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
  output_simCRC_Min_MSE <- calibration_out(l_params_all = l_params_init,v_params_calib = l_params_calibrated, id_draw = l_params_calibrated$id_draw, dt_pop, dt_long = TRUE)
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
  output_simCRC_Post_mean <- calibration_out(l_params_all = l_params_init,v_params_calib = l_params_calibrated, id_draw = l_params_calibrated$id_draw, dt_pop, dt_long = TRUE)
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
  output_simCRC_Post_median <- calibration_out(l_params_all = l_params_init,v_params_calib = l_params_calibrated, id_draw = l_params_calibrated$id_draw, dt_pop, dt_long = TRUE)
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

l_params_calibrated_sets$Min_AbsolutErr
save(l_params_calibrated_sets, file = paths_calibration$path_best_params_sets)

# Append all outputs from different parameters set ------------------------

# Create a data frame to hold all selected sets

df_outputs_selected_sets <- rbind(     output_simCRC_max_lp, 
                                       output_simCRC_Max_log_likelihood, 
                                       output_simCRC_Min_AbsolutErr, 
                                       output_simCRC_Min_MSE,
                                       output_simCRC_Post_mean,
                                       output_simCRC_Post_median)


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

color_sets <- c("Max_lp" = "#E41A1C", 
                 "Max_log_likelihood" = "#377EB8", 
                 "Min_AbsolutErr" = "#4DAF4A", 
                 "Min_MSE" = "#984EA3",
                 "Post_mean" = "#FF7F00",
                 "Post_median" = "#E69F00")

#identification of categorical variables

cat_groups <- c("Size_P", "Size_D", "Size_R", "StageDist_D","StageDist_P","StageDist_R", "SSP")
df_outputs_selected_sets$categorical <- ifelse(df_outputs_selected_sets$target_groups %in% cat_groups,1,0)


plot_val_num <- ggplot(data = df_outputs_selected_sets[df_outputs_selected_sets$categorical == 0, ], 
                       aes(x    = age, 
                           y    = targets, 
                           ymin = stopping_lower_bounds, 
                           ymax = stopping_upper_bounds)) + 
  geom_errorbar(width = 1.2, size = 0.5, color = "red") +
  geom_line(data = df_outputs_selected_sets[df_outputs_selected_sets$categorical == 0, ],
            aes(x     = age,
                y     = output_model,
                color = type_param_set),
            linewidth = 0.5) +
  facet_wrap(~ target_groups + lesion_type, scales = "free", ncol = 3) +
  scale_color_manual(
    name   = "Set selection",
    values = color_sets
  ) +
  scale_y_continuous(breaks = number_ticks(5)) +
  theme_bw(base_size = 8) +
  theme(
    plot.title       = element_text(size = 12, face = "bold"),
    axis.text.x      = element_text(size = 8, angle = 90),
    axis.title       = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.border     = element_rect(colour = "black", fill = NA),
    strip.background = element_blank(),
    strip.text       = element_text(hjust = 0, size=10)  ) +
  labs(title = "SimCRC-R calibrated parameter set validation", 
       x     = "")

plot_val_num

ggsave(plot_val_num,
       filename = paste0(folder, "/fig_internal_validation_", BayCANN_version, "_num_all_sets.png"),
       width = 10, height = 6)



# First, create numeric positions for your groups
df_outputs_selected_sets$x_pos <- as.numeric(factor(df_outputs_selected_sets$stage))
shift <- 0.2

# Recode stage labels, handling NAs safely
df_outputs_selected_sets$stage <- ifelse(
  is.na(df_outputs_selected_sets$stage), NA,
  ifelse(df_outputs_selected_sets$stage == "LR", "Low risk (LR)",
         ifelse(df_outputs_selected_sets$stage == "MR", "Medium risk (MR)",
                ifelse(df_outputs_selected_sets$stage == "HR", "High risk (HR)",
                       df_outputs_selected_sets$stage)))  # keeps "1","2","3","4" as-is
)

# Define ordered factor levels (NAs will fall outside and be dropped naturally)
df_outputs_selected_sets$stage <- factor(
  df_outputs_selected_sets$stage,
  levels = c("1", "2", "3", "4",
             "Low risk (LR)", "Medium risk (MR)", "High risk (HR)")
)

# Recompute x_pos after releveling (NAs in stage → NA in x_pos, safe to plot)
df_outputs_selected_sets$x_pos <- as.numeric(df_outputs_selected_sets$stage)

plot_val_cat <- ggplot(
  data = df_outputs_selected_sets[df_outputs_selected_sets$categorical == 1, ],
  aes(x    = x_pos - shift,
      y    = targets,
      ymin = stopping_lower_bounds,
      ymax = stopping_upper_bounds)) +
  # Target point + error bar (red)
  geom_point(aes(shape = "Target"), color = "red", size = 1.5) +
  geom_errorbar(width = 0.2, linewidth = 0.5, alpha = 0.5, color = "red") +
  # Model crossbar colored by parameter set
  geom_errorbar(
    data = df_outputs_selected_sets[df_outputs_selected_sets$categorical == 1, ],
    aes(x    = x_pos + shift,
        y    = output_model,
        ymin = output_model,
        ymax = output_model,
        color = type_param_set),
    width = 0.2, linewidth = 0.9, alpha = 1) +
  # Model point colored by parameter set
  geom_point(
    data = df_outputs_selected_sets[df_outputs_selected_sets$categorical == 1, ],
    aes(x     = x_pos + shift,
        y     = output_model,
        color = type_param_set,        # ← inherits parameter set color
        shape = "Model"),
    size = 1.5) +
  facet_wrap(~ target_groups + lesion_type, scales = "free", ncol = 3) +
  scale_color_manual(
    name   = "Set selection",
    values = color_sets) +
  scale_shape_manual(
    name   = "",
    values = c("Target" = 8, "Model" = 16),
    guide  = guide_legend(
      override.aes = list(
        color = c("red", "black"),     # ← Target = red star, Model = single black dot
        size  = c(2, 2)
      )
    )) +
  scale_y_continuous(breaks = number_ticks(5)) +
  scale_x_continuous(
    breaks = df_outputs_selected_sets$x_pos,
    labels = df_outputs_selected_sets$stage) +
  theme_bw(base_size = 8) +
  theme(
    plot.title       = element_text(size = 12, face = "bold"),
    axis.text.x      = element_text(size = 8, angle = 0),
    axis.title       = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.border     = element_rect(colour = "black", fill = NA),
    strip.background = element_blank(),
    strip.text       = element_text(hjust = 0, size = 10)) +
  labs(title = "SimCRC-R calibrated parameter set validation",
       x     = "")

plot_val_cat

  ggsave(plot_val_cat,
       filename = paste0(folder, "/fig_internal_validation_", BayCANN_version, "_cat_all_sets.png"),
       width = 10, height = 6)


#Now lets do it only for Min_AbsolutErr set, which is the one with the best fit to the targets
  
#plot_val_num 
  
ggplot(data = df_outputs_selected_sets[df_outputs_selected_sets$categorical == 0 & df_outputs_selected_sets$type_param_set == "Min_AbsolutErr", ], 
       aes(x    = age, 
           y    = targets, 
           ymin = stopping_lower_bounds, 
           ymax = stopping_upper_bounds)) + 
  geom_errorbar(width = 1.2, size = 0.5, color = "red") +
  geom_line(data = df_outputs_selected_sets[df_outputs_selected_sets$categorical == 0 & df_outputs_selected_sets$type_param_set == "Min_AbsolutErr", ],
            aes(x     = age,
                y     = output_model),
            color = color_sets["Min_AbsolutErr"], linewidth = 0.5) +
  facet_wrap(~ target_groups + lesion_type, scales = "free", ncol = 3) +
  scale_y_continuous(breaks = number_ticks(5)) +
  theme_bw(base_size = 8) +
  theme(
    plot.title       = element_text(size = 12, face = "bold"),
    axis.text.x      = element_text(size = 8, angle = 90),
    axis.title       = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.border     = element_rect(colour = "black", fill = NA),
    strip.background = element_blank(),
    strip.text       = element_text(hjust = 0, size=10)  ) +
  labs(title = "SimCRC-R calibrated parameter set validation (Set: Min_AbsolutErr)", 
       x     = "")


ggsave(
  filename = paste0(folder, "/fig_internal_validation_", BayCANN_version, "_num_Min_AbsolutErr.png"),
  width = 10, height = 6
)


#plot_val_cat

ggplot(
  data = df_outputs_selected_sets[df_outputs_selected_sets$categorical == 1 & df_outputs_selected_sets$type_param_set == "Min_AbsolutErr", ],
  aes(x    = x_pos - shift,
      y    = targets,
      ymin = stopping_lower_bounds,
      ymax = stopping_upper_bounds)) +
  # Target point + error bar (red)
  geom_point(aes(shape = "Target"), color = "red", size = 1.5) +
  geom_errorbar(width = 0.2, linewidth = 0.5, alpha = 0.5, color = "red") +
  # Model crossbar colored by parameter set
  geom_errorbar(
    data = df_outputs_selected_sets[df_outputs_selected_sets$categorical == 1 & df_outputs_selected_sets$type_param_set == "Min_AbsolutErr", ],
    aes(x    = x_pos + shift,
        y    = output_model,
        ymin = output_model,
        ymax = output_model),
    width = 0.2, linewidth = 0.9, alpha = 1, color = color_sets["Min_AbsolutErr"]) +
  # Model point colored by parameter set
  geom_point(
    data = df_outputs_selected_sets[df_outputs_selected_sets$categorical == 1 & df_outputs_selected_sets$type_param_set == "Min_AbsolutErr", ],
    aes(x     = x_pos + shift,
        y     = output_model,
        shape = "Model"),
    size = 1.5, color = color_sets["Min_AbsolutErr"]) +
  facet_wrap(~ target_groups + lesion_type, scales = "free", ncol = 3) +
  scale_shape_manual(
    name   = "",
    values = c("Target" = 8, "Model" = 16),
    guide  = guide_legend(
      override.aes = list(
        color = c("red", color_sets["Min_AbsolutErr"]),     # ← Target = red star, Model colored by parameter set
        size  = c(2, 2)
      )
    )) +
  scale_y_continuous(breaks = number_ticks(5)) +
  scale_x_continuous(
    breaks = df_outputs_selected_sets$x_pos,
    labels = df_outputs_selected_sets$stage) +
  theme_bw(base_size = 8) +
  theme(
    plot.title       = element_text(size = 12, face = "bold"),
    axis.text.x      = element_text(size = 8, angle = 0),
    axis.title       = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.border     = element_rect(colour = "black", fill = NA),
    strip.background = element_blank(),
    strip.text       = element_text(hjust = 0, size = 10)) +
  labs(title = "SimCRC-R calibrated parameter set validation (Set: Min_AbsolutErr)",
       x     = "")


ggsave(
  filename = paste0(folder, "/fig_internal_validation_", BayCANN_version, "_cat_Min_AbsolutErr.png"),
  width = 10, height = 6
)




# Compare the symptom detection probabilities for distal colon cancer
df_D <- NULL

# Parameters that fit the US-F stage distribution
pSxDetS1_D <- 0.0485
hr_SxDetS2S1_D <- 5.5
hr_SxDetS3S2_D <- 3.1
hr_SxDetS4S3_D <- 4.7

# Parameters that fit the Chile stage distribution
pSxDetS1_D <- l_params_calibrated_sets$Min_AbsolutErr$pSxDetS1_D
hr_SxDetS2S1_D <- l_params_calibrated_sets$Min_AbsolutErr$hr_SxDetS2S1_D
hr_SxDetS3S2_D <- l_params_calibrated_sets$Min_AbsolutErr$hr_SxDetS3S2_D
hr_SxDetS4S3_D <- l_params_calibrated_sets$Min_AbsolutErr$hr_SxDetS4S3_D

# Convert parameters to rates
rSxDetS1_D <- -log(1 - pSxDetS1_D)
rSxDetS2_D <- hr_SxDetS2S1_D*rSxDetS1_D
rSxDetS3_D <- hr_SxDetS3S2_D*rSxDetS2_D
rSxDetS4_D <- hr_SxDetS4S3_D*rSxDetS3_D

# Convert rates to probs
pSxDetS2_D <- 1 - exp(-rSxDetS2_D)
pSxDetS3_D <- 1 - exp(-rSxDetS3_D)
pSxDetS4_D <- 1 - exp(-rSxDetS4_D)

# Store probs in a table that you can keep appending to
df_D <- rbind(df_D, data.frame(
  id = "Chile",
  pSxDetS1_D = pSxDetS1_D,
  pSxDetS2_D = pSxDetS2_D,
  pSxDetS3_D = pSxDetS3_D,
  pSxDetS4_D = pSxDetS4_D
))

# df_D
# id pSxDetS1_D pSxDetS2_D pSxDetS3_D pSxDetS4_D
# US-F     0.0485  0.2392388  0.5715798  0.9813882
# Chile     0.0190  0.1172280  0.5608626  0.5608626



## Repeat for proximal colon (P)
df_P <- NULL

# Parameters that fit the US-F stage distribution
pSxDetS1_P <- 0.014
hr_SxDetS2S1_P <- 8.7
hr_SxDetS3S2_P <- 5.4
hr_SxDetS4S3_P <- 4.4

l_params_calibrated_sets$Min_AbsolutErr$pSxDetS1_P

# Parameters that fit the Chile stage distribution
pSxDetS1_P <- l_params_calibrated_sets$Min_AbsolutErr$pSxDetS1_P
hr_SxDetS2S1_P <- l_params_calibrated_sets$Min_AbsolutErr$hr_SxDetS2S1_P
hr_SxDetS3S2_P <- l_params_calibrated_sets$Min_AbsolutErr$hr_SxDetS3S2_P
hr_SxDetS4S3_P <- l_params_calibrated_sets$Min_AbsolutErr$hr_SxDetS4S3_P

# Convert parameters to rates
rSxDetS1_P <- -log(1 - pSxDetS1_P)
rSxDetS2_P <- hr_SxDetS2S1_P*rSxDetS1_P
rSxDetS3_P <- hr_SxDetS3S2_P*rSxDetS2_P
rSxDetS4_P <- hr_SxDetS4S3_P*rSxDetS3_P

# Convert rates to probs
pSxDetS2_P <- 1 - exp(-rSxDetS2_P)
pSxDetS3_P <- 1 - exp(-rSxDetS3_P)
pSxDetS4_P <- 1 - exp(-rSxDetS4_P)

# Store probs in a table that you can keep appending to
df_P <- rbind(df_P, data.frame(
  id = "Chile",
  pSxDetS1_P = pSxDetS1_P,
  pSxDetS2_P = pSxDetS2_P,
  pSxDetS3_P = pSxDetS3_P,
  pSxDetS4_P = pSxDetS4_P
))

# df_P
# id pSxDetS1_P pSxDetS2_P pSxDetS3_P pSxDetS4_P
# US-F      0.014 0.11543620  0.4843708  0.9457644
# Chile      0.008900666 0.05163168  0.4535654  0.4535654


## Repeat for rectal (R)
df_R <- NULL

# Parameters that fit the US-F stage distribution
pSxDetS1_R <- 0.07
hr_SxDetS2S1_R <- 2.04
hr_SxDetS3S2_R <- 6.5
hr_SxDetS4S3_R <- 1.8

# Parameters that fit the Chile stage distribution
pSxDetS1_R <- l_params_calibrated_sets$Min_AbsolutErr$pSxDetS1_R
hr_SxDetS2S1_R <- l_params_calibrated_sets$Min_AbsolutErr$hr_SxDetS2S1_R
hr_SxDetS3S2_R <- l_params_calibrated_sets$Min_AbsolutErr$hr_SxDetS3S2_R
hr_SxDetS4S3_R <- l_params_calibrated_sets$Min_AbsolutErr$hr_SxDetS4S3_R

# Convert parameters to rates
rSxDetS1_R <- -log(1 - pSxDetS1_R)
rSxDetS2_R <- hr_SxDetS2S1_R*rSxDetS1_R
rSxDetS3_R <- hr_SxDetS3S2_R*rSxDetS2_R
rSxDetS4_R <- hr_SxDetS4S3_R*rSxDetS3_R

# Convert rates to probs
pSxDetS2_R <- 1 - exp(-rSxDetS2_R)
pSxDetS3_R <- 1 - exp(-rSxDetS3_R)
pSxDetS4_R <- 1 - exp(-rSxDetS4_R)

# Store probs in a table that you can keep appending to
df_R <- rbind(df_R, data.frame(
  id = "Chile",
  pSxDetS1_R = pSxDetS1_R,
  pSxDetS2_R = pSxDetS2_R,
  pSxDetS3_R = pSxDetS3_R,
  pSxDetS4_R = pSxDetS4_R
))

