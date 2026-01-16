###########################  BayCANN   #########################################
#
#  Objective: Script to perform Bayesian calibration using ANN 
########################### <<<<<>>>>> ##############################################


#### 2. General parameters ========================================================

###### 2.1 parameters for Data preparation 
scale_type = 1  ## 1: for scale from -1 to 1; 2: for standardization ; 3: for scale from 0 to 1
set.seed(1234)

###### 2.1 parameters for ANN 
verbose          <- l_ANN_settings$verbose
n_batch_size     <- l_ANN_settings$n_batch_size
n_epochs         <- l_ANN_settings$n_epochs
patience         <- l_ANN_settings$patience
n_hidden_nodes   <- l_ANN_settings$n_hidden_nodes
n_hidden_layers  <- l_ANN_settings$n_hidden_layers
activation_fun   <- l_ANN_settings$activation_fun
init_W           <- l_ANN_settings$init_W


###### 2.2 parameters for Bayesian calibration
n_iter   <- l_ANN_settings$n_iter
n_thin   <- l_ANN_settings$n_thin
n_chains <- l_ANN_settings$n_chains

#### 3. Pre-processing actions  ===========================================

Normalize_inputs    <- l_ANN_settings$Normalize_inputs
Normalize_outputs   <- l_ANN_settings$Normalize_outputs
Scale_inputs        <- l_ANN_settings$Scale_inputs
Scale_outputs       <- l_ANN_settings$Scale_outputs
Remove_outliers     <- l_ANN_settings$Remove_outliers
Standardize_targets <- l_ANN_settings$Standardize_targets
Saved_data          <- l_ANN_settings$Saved_data
Selected_targets    <- l_ANN_settings$Selected_targets

#### 4. Load the input and output data for the simulations =======================


params_file <- outputs_file <- paths_calibration$path_lhs #LHS file with parameters and outputs

selected_params  <- NULL

load(outputs_file)

data_sim_param   <- data_LHS[[2]]  #LHS (prior parameter sampling)
data_sim_target  <- data_LHS[[1]]  # LHS-corresponding outputs

data_sim_param  <- as.data.frame(data_sim_param)
data_sim_target <- as.data.frame(data_sim_target)

# Load calibration targets
data_true_targets  <- read.csv(targets_file)

#keep only parameters with corresponding targets
id_draw_with_targets <- unique(data_sim_target$id_draw)
data_sim_param <- data_sim_param[(data_sim_param$id_draw %in% id_draw_with_targets),]

#Get and save priors
l_priors <- sample_uniform_lhs(n_samp = 1)

l_prior_lb <- l_priors$lower_bounds
l_prior_ub <- l_priors$upper_bounds
l_priors <- list(l_prior_lb = l_prior_lb , l_prior_ub=l_prior_ub)

path_baycann_priors <- paths_calibration$path_priors
save(l_priors, file = path_baycann_priors)

#Check if selected targets that are not included in the LHS outputs
if(Selected_targets){
  v_missing_targets <- setdiff(selected_targets, colnames(data_sim_target))
  if(length(v_missing_targets)>0){
    stop(paste("The following selected targets are not present in the LHS outputs:", paste(v_missing_targets, collapse = ", ")))
  }
}

#Drop variables in data_sim_target that are not included in data_true_targets
v_true_target_names <- data_true_targets$target_names
v_cols_to_drop <- setdiff(colnames(data_sim_target), v_true_target_names)
data_sim_target <- data_sim_target[, !(colnames(data_sim_target) %in% v_cols_to_drop)]

#Print removed targets
if(length(v_cols_to_drop)>0){
  print(paste("The following targets were removed from the LHS outputs as they are not present in the true targets:", paste(v_cols_to_drop, collapse = ", ")))
}

#Keep only selected targets
if(Selected_targets){
  data_sim_target <- subset(data_sim_target, select = selected_targets)
  print("Only selected targets kept")
}
#Keep only parameters to calibrate 
if("id_draw" %in% colnames(data_sim_param)) {
  data_sim_param <- subset(data_sim_param, select = -c(id_draw))
  print("id_draw column removed from parameters")
}


# check for NAs
v_NAs <- apply(data_sim_target, 1, function(x) any(is.na(x)))
table(v_NAs)

# Remove NAs from parameters and outputs 
data_sim_param  <- data_sim_param[!v_NAs,]
data_sim_target <- data_sim_target[!v_NAs,]

#Remove prevalence higher than 1
v_prevalence <- data_sim_target$Prev0AdOrPreClin_age27<=1
data_sim_param  <- data_sim_param[v_prevalence,]
data_sim_target <- data_sim_target[v_prevalence,]

#Confirm that there are no remaining NAs
v_NAs <- apply(data_sim_target, 1, function(x) any(is.na(x)))
table(v_NAs) #Check if there are NAs


###### 4.1 Get list of params and targets ####

###### 4.1 Input Parameters ####

data_sim_param <- as.matrix(data_sim_param)

###### 4.2 Model Outputs ####

data_sim_target <- as.matrix(data_sim_target)

#### 5. Removing outliers from model outputs ####

if (Remove_outliers) {
  vec_out <- outlier_vector(data_sim_target)
  data_sim_param  <- data_sim_param[!vec_out,]
  data_sim_target <- data_sim_target[!vec_out,]
}


#### 6. Normalize distributions ===================================================

#(parameters)
if (Normalize_inputs) {
  param_norm        <- best_normal_dataset(data_sim_param)
  data_sim_param    <- param_norm[["data_normal"]]
  inv_param_transf  <- param_norm[["inverse_dist"]]
}

#(model outputs)
if (Normalize_outputs) {
  target_norm          <- best_normal_dataset(data_sim_target)
  data_sim_target      <- target_norm[["data_normal"]]
  inv_target_transform <- target_norm[["inverse_dist"]]
}


#### 7. Train/test partition ======================================================

library(caTools)
set.seed(1223)

train.rows<-sample.split(data_sim_param[,1],SplitRatio=0.8)

data_sim_param_train  <- data_sim_param[train.rows,]
data_sim_param_test   <- data_sim_param[!train.rows,]

data_sim_target_train <- data_sim_target[train.rows,]
data_sim_target_test  <- data_sim_target[!train.rows,]

prepared_data <- prepare_data(xtrain = data_sim_param_train,
                              ytrain = data_sim_target_train,
                              xtest  = data_sim_param_test,
                              ytest  = data_sim_target_test,
                              scale  = scale_type)

list2env(prepared_data, envir = .GlobalEnv)

#### 8. Unscale to keep the original scale of output model ========================

if (!Scale_inputs) {
  xtrain_scaled <- unscale_data(xtrain_scaled, vec.mins = xmins, vec.maxs = xmaxs, vec.means=xmeans, vec.sds = xsds, type = scale_type)
  xtest_scaled <- unscale_data(xtest_scaled, vec.mins = xmins, vec.maxs = xmaxs, vec.means=xmeans, vec.sds = xsds, type = scale_type)
}

if (!Scale_outputs) {
  ytrain_scaled <- unscale_data(ytrain_scaled, vec.mins = ymins, vec.maxs = ymaxs, vec.means=ymeans, vec.sds = ysds, type = scale_type)
  ytest_scaled <- unscale_data(ytest_scaled, vec.mins = ymins, vec.maxs = ymaxs, vec.means=ymeans, vec.sds = ysds, type = scale_type)
}

####  9. Load the targets and their se ============================================

#keep targets of interest
v_targets_names  <- colnames(data_sim_target) 
data_true_targets  <- data_true_targets[(data_true_targets$target_names%in%v_targets_names),]   

if(Selected_targets) {
  data_true_targets  <- data_true_targets[(data_true_targets$target_names%in%selected_targets),]   #Selection of 56 targets
}

true_targets_mean  <- data_true_targets$targets
true_targets_upper <- data_true_targets$stopping_upper_bounds
true_targets_lower <- data_true_targets$stopping_lower_bounds
true_targets_se   <- (true_targets_upper - true_targets_lower)/(2*1.96)
#true_targets_se   <- data_true_targets$se   #when targets already had the se data


#standardize with respect to the true targets

if(Standardize_targets) {
  
  if(Scale_outputs) {
    true_targets_mean <- 2 * (true_targets_mean - ymins) / (ymaxs - ymins) - 1   ## range from -1 to 1
    true_targets_se <- 2 * (true_targets_se) / (ymaxs - ymins)
  }  
  
  for (i in 1:length(true_targets_mean)) {
    ytrain_scaled[,i] <- (ytrain_scaled[,i] - true_targets_mean[i]) / true_targets_se[i] 
    ytest_scaled[,i]  <- (ytest_scaled[,i] - true_targets_mean[i]) / true_targets_se[i] 
  }
}

#### 10. Scale the targets and their SE  ####

if (scale_type==1) {
  y_targets <- 2 * (true_targets_mean - ymins) / (ymaxs - ymins) - 1   ## range from -1 to 1
  y_targets_se <- 2 * (true_targets_se) / (ymaxs - ymins)
}

if (scale_type==2) {
  y_targets <- (true_targets_mean - ymeans)/ysds   ## Standardization
  y_targets_se <-(true_targets_se)/ysds
}

if (scale_type==3) {
  y_targets <- (true_targets_mean - ymins) / (ymaxs - ymins)   ## range from 0 to 1
  y_targets_se <-(true_targets_se) / (ymaxs - ymins)
}

y_targets <- t(as.matrix(y_targets))
y_targets_se <- t(as.matrix(y_targets_se))   

# converting "y_targets" to 0 and "y_targets_se" to 1
if( Standardize_targets) {
  y_targets <- y_targets - y_targets
  y_targets_se <- y_targets_se / y_targets_se
}

#### 11. Keras Section BayCANN SimCRC ==============================================

# File name of keras model

path_keras_model <- paste0(folder,"/model_keras_SimCRC_", cal_version, "_", Lesion_type, "_",Population_type,".h5")    ##File path for the compiled model

#Initializers
#init_W=initializer_random_uniform(minval = -0.7, maxval = 0.7,seed = 2312)   ###initialization of weights with uniform distribution
#init_W=initializer_random_normal(mean = 0, stddev = 0.1, seed = 2312)  ###initialization of weights with normal distribution

model <- keras_model_sequential()
mdl_string <- paste("model %>% layer_dense(units = n_hidden_nodes, kernel_initializer=init_W, activation = activation_fun, input_shape = n_inputs) %>%",
                    paste(rep(x = "layer_dense(units = n_hidden_nodes, activation = activation_fun) %>%",
                              (n_hidden_layers)), collapse = " "),
                    "layer_dense(units = n_outputs)")

inputs <- layer_input(shape = c(n_inputs))
x <- layer_dense(inputs, units = n_hidden_nodes, kernel_initializer = init_W, activation = activation_fun)

for (i in seq_len(n_hidden_layers)) {
  x <- layer_dense(x, units = n_hidden_nodes, activation = activation_fun)
}

outputs <- layer_dense(x, units = n_outputs)
model <- keras_model(inputs = inputs, outputs = outputs)



summary(model)

#model %>% compile(
#  loss = 'mean_squared_error',
#  optimizer = 'adam'  ,
#  metrics = list('mae',"accuracy")
#)

model$compile(
  loss = 'mean_squared_error',
  optimizer = 'adam',
  metrics = list('mae', 'accuracy')
)


keras.time <- proc.time()


#history <- model %>% fit(
#  xtrain_scaled, ytrain_scaled,
#  epochs = n_epochs,
#  batch_size = n_batch_size,
#  validation_data = list(xtest_scaled, ytest_scaled),
#  verbose = verbose,
#  callback_early_stopping(
#    monitor = "val_loss",
#    patience = patience,
#    verbose = 0,
#    restore_best_weights = TRUE
#  )
#)

n_batch_size <- as.integer(n_batch_size)
n_epochs     <- as.integer(n_epochs)
patience     <- as.integer(patience)


history <- model$fit(
  x = xtrain_scaled,
  y = ytrain_scaled,
  epochs = n_epochs,
  batch_size = n_batch_size,
  validation_data = list(xtest_scaled, ytest_scaled),
  verbose = verbose,
  callbacks = list(
    callback_early_stopping(
      monitor = "val_loss",
      patience = patience,
      verbose = 0,
      restore_best_weights = TRUE
    )
  )
)


t_training <- proc.time() - keras.time #keras ann fitting time

t_training <- t_training[3]/60

acc_err<-model$evaluate(xtest_scaled,ytest_scaled) # Model performance evaluation
metric_loss      <- acc_err[1]
metric_mae       <- acc_err[2]
metric_accuracy  <- acc_err[3]

save_model_hdf5(model,path_keras_model)  #Save the model
model <- load_model_hdf5(path_keras_model)


###### 11.4 History Graph ####

#plot(history$history)   #Plot loss function and accuracy function

# Extract history into dataframe
df_history <- as.data.frame(history$history)
df_history$epoch <- seq_len(nrow(df_history))

library(tidyr)
df_long <- pivot_longer(df_history, 
                        cols = -epoch, 
                        names_to = "variable", 
                        values_to = "value")


# Plot like old keras::plot(history)
library(ggplot2)

ggplot(df_long, aes(x = epoch, y = value, color = variable)) +
  geom_line(size = 1.2) +
  labs(title = "Training History", x = "Epoch", y = "Value") +
  theme_minimal()


###### 11.5 Prediction Graph  ####

pred <- model$predict(xtest_scaled)
ytest_scaled_pred <- data.frame(pred)
colnames(ytest_scaled_pred) <- y_names
head(ytest_scaled_pred)    #

ann_valid <- rbind(data.frame(sim = 1:n_test, ytest_scaled, type = "model"),
                   data.frame(sim = 1:n_test, ytest_scaled_pred, type = "pred"))
ann_valid_transpose <- ann_valid %>%
  pivot_longer(cols = -c(sim, type)) %>%
  pivot_wider(id_cols = c(sim, name), names_from = type, values_from = value)

##Partition of validation data for Graph (3 parts)
n_partition <-1
n_part_bach <-floor(n_outputs/n_partition)

ann_valid_transpose <- arrange(ann_valid_transpose,desc(name))

ann_valid_transpose1 <- ann_valid_transpose[(1):(n_part_bach*n_test),]
#ann_valid_transpose2 <- ann_valid_transpose[(n_part_bach*n_test+1):(2*n_part_bach*n_test),]
#ann_valid_transpose3 <- ann_valid_transpose[(2*n_part_bach*n_test+1):(3*n_part_bach*n_test),]
#ann_valid_transpose4 <- ann_valid_transpose[(3*n_part_bach*n_test+1):dim(ann_valid_transpose)[1],]

#part 1
pperformance <- ggplot(data = ann_valid_transpose1, aes(x = model, y = pred)) +
  geom_point(alpha = 0.5, color = "tomato") +
  facet_wrap(~name,  ncol = 10) +
  xlab("Model outputs (scaled)") +
  ylab("ANN predictions (scaled)") +
  #coord_equal() +
  theme_bw()

# #part 2
# ggplot(data = ann_valid_transpose2, aes(x = model, y = pred)) +
#   geom_point(alpha = 0.5, color = "tomato") +
#   facet_wrap(~name, scales="free", ncol = 7) +
#   xlab("Model outputs (scaled)") +
#   ylab("ANN predictions (scaled)") +
#   #coord_equal() +
#   theme_bw()
# 
# #part 3
# ggplot(data = ann_valid_transpose3, aes(x = model, y = pred)) +
#   geom_point(alpha = 0.5, color = "tomato") +
#   facet_wrap(~name, scales="free", ncol = 7) +
#   xlab("Model outputs (scaled)") +
#   ylab("ANN predictions (scaled)") +
#   #coord_equal() +
#   theme_bw()
# 
# #part 4
# ggplot(data = ann_valid_transpose4, aes(x = model, y = pred)) +
#   geom_point(alpha = 0.5, color = "tomato") +
#   facet_wrap(~name, scales="free", ncol = 7) +
#   xlab("Model outputs (scaled)") +
#   ylab("ANN predictions (scaled)") +
#   #coord_equal() +
#   theme_bw()

ggsave(filename = paths_calibration$path_ann_perform,
       pperformance,
       width = 12, height = 12)

#### 12. Stan section ==============================================================

path_posterior <- paths_calibration$path_posteriors

weights <- get_weights(model) #get ANN weights

n_hidden_layers <- length(weights)/2-2    #Removing bias layers and input and output layers  (Carlos P)
n_hidden_nodes  <- ncol(weights[[1]])     #Get number of hidden nodes from the firs layers (Carlos P)

# pass the weights and biases to Stan for Bayesian calibration
n_layers <- length(weights)
weight_first <- weights[[1]]
beta_first <- 1 %*% weights[[2]]
weight_last <- weights[[n_layers-1]]
beta_last <- 1 %*% weights[[n_layers]]
weight_middle <- array(0, c(n_hidden_layers, n_hidden_nodes, n_hidden_nodes))
beta_middle <- array(0, c(n_hidden_layers, 1, n_hidden_nodes))
for (l in 1:n_hidden_layers){
  weight_middle[l,,] <- weights[[l*2+1]]
  beta_middle[l,,] <- weights[[l*2+2]]
}

###Información para inferencia en Stan
stan.dat=list(
  num_hidden_nodes = n_hidden_nodes,
  num_hidden_layers= n_hidden_layers,
  num_inputs=n_inputs,
  num_outputs=n_outputs,
  num_targets=1,
  y_targets = y_targets,
  y_targets_se = y_targets_se,
  beta_first = beta_first,
  beta_middle = beta_middle,
  beta_last = beta_last,
  weight_first = weight_first,
  weight_middle = weight_middle,
  weight_last = weight_last)

# Select the stan file based on data transformation

if (Normalize_inputs) {
  file_perceptron <- "analysis/post_multi_perceptron_normal.stan"
} else {
  file_perceptron <- "analysis/post_multi_perceptron.stan"  
}

# Run stan file
stan.time <- proc.time()
m <- stan(file = file_perceptron,
          data = stan.dat,
          iter = n_iter,
          chains = n_chains,
          thin = n_thin,
          pars = c("Xq"),
          warmup = floor(n_iter/2),   ## (cp)
          seed = 12345) #for reproducibility. R's set.seed() will not work for stan
t_calibration <- proc.time() - stan.time # stan sampling time
t_calibration <- t_calibration[3] / 60
summary(m)

path_stan_model <- paths_calibration$path_stan_model

param_names    <- colnames(data_sim_param)

names(m)[1:n_inputs] <- param_names

saveRDS(m, path_stan_model)

m <- readRDS(path_stan_model) 

###### 12.1 Stan Diagnose  ----

stan_trace(m,pars=param_names,inc_warmup = FALSE)

stan_plot(m,pars=param_names, point_est = "mean", show_density = TRUE, fill_color = "maroon", ncol=2)

stan_hist(m,pars=param_names, inc_warmup = FALSE)

standensity <- stan_dens(m,pars=param_names, inc_warmup = FALSE, separate_chains=TRUE)

ggsave(filename = paths_calibration$path_post_chains,
       standensity,
       width = 24, height = 16)

stan_dens(m,pars=param_names, inc_warmup = FALSE, separate_chains=FALSE)

stan_ac(m,pars=param_names[1:15], inc_warmup = FALSE, separate_chains=TRUE)
stan_ac(m,pars=param_names[16:33], inc_warmup = FALSE, separate_chains=TRUE)

stan_rhat(m,pars=param_names)          # Rhat statistic 
stan_par(m,par=param_names[1])         # Mean metrop. acceptances, sample step size
stan_ess(m,pars=param_names)           # Effective sample size / Sample size
stan_mcse(m,pars=param_names)          # Monte Carlo SE / Posterior SD
stan_diag(m,)

###### 12.2 Stan extraction  ----

params <- rstan::extract(m, permuted=TRUE, inc_warmup = FALSE)
lp <- params$lp__
Xq <- params$Xq
Xq_df = as.data.frame(Xq)


# Scale the posteriors
if (Scale_inputs) {
  Xq_unscaled <- unscale_data(Xq_df, vec.mins = xmins, vec.maxs = xmaxs, vec.means = xmeans, vec.sds = xsds, type = scale_type)
} else {
  Xq_unscaled <- Xq_df
}

##### Inverse normalization  
if (Normalize_inputs) {
  Xq_unscaled <- normalize_predict(df=Xq_unscaled, list_transf = inv_param_transf,inverse = TRUE)
  data_sim_param_train <- normalize_predict(df=data_sim_param_train, list_transf = inv_param_transf,inverse = TRUE)
}

Xq_lp <- cbind(Xq_unscaled, lp)
colnames(Xq_lp) <- m@sim[["fnames_oi"]]

# Add id_draw
Xq_lp <- data.frame(id_draw = 1:nrow(Xq_lp), Xq_lp)

# Add chain id
chain_window = m@sim[["warmup2"]][1]
chain_id <- rep(1:n_chains, each = chain_window)
Xq_lp$chain <- chain_id

# Save the unscaled posterior samples
write.csv(Xq_lp,
          file = path_posterior,
          row.names = FALSE)

cal_mdl_1 <- path_posterior
### Load posterior
Xq_lp <- read.csv(file = cal_mdl_1)
n_col <- ncol(Xq_lp)
lp <- Xq_lp$lp__
#Remove id_draw, lp__ and chain columns
if("id_draw" %in% colnames(Xq_lp)) Xq_lp$id_draw <- NULL
if("chain" %in% colnames(Xq_lp)) Xq_lp$chain <- NULL

map_baycann <- Xq_unscaled[which.max(lp), ]     ### MAP for first BayCANN model

df_post_ann <- read.csv(file = cal_mdl_1)
if("lp__" %in% colnames(df_post_ann)) df_post_ann$lp__ <- NULL
if("id_draw" %in% colnames(df_post_ann)) df_post_ann$id_draw <- NULL
if("chain" %in% colnames(df_post_ann)) df_post_ann$chain <- NULL
colnames(df_post_ann) <- x_names


#gg_calib_post_pair_corr
# 
# ggsave(filename = paths_calibration$path_post_pairs,
#        gg_calib_post_pair_corr,
#        width = 24, height = 20)


#### Prior and prior graph
n_samp <- 2000
df_samp_prior <- reshape2::melt(cbind(Distribution = "Prior",
                                      as.data.frame(data_sim_param_train[1:n_samp, ])),
                                variable.name = "Parameter")

df_samp_post_ann   <- reshape2::melt(cbind(Distribution = "Posterior SimCRC SSP",
                                           as.data.frame(df_post_ann[1:n_samp, ])),
                                     variable.name = "Parameter")


#remove lp
#df_post_adenoma <- df_post_adenoma[, -ncol(df_post_adenoma)]
#
#df_samp_post_adenoma <- reshape2::melt(cbind(Distribution = "Posterior SimCRC Adenoma only",
#                                             as.data.frame(df_post_adenoma[1:n_samp, ])),
#                                       variable.name = "Parameter")


df_samp_prior_post <- rbind(df_samp_prior,
                            df_samp_post_ann)




#df_samp_prior_post$Distribution <- ordered(df_samp_prior_post$Distribution,
#                                           levels = levels)



df_samp_prior_post$Parameter <- factor(df_samp_prior_post$Parameter,
                                       levels = levels(df_samp_prior_post$Parameter),
                                       ordered = TRUE)

df_maps_n_true_params <- data.frame(Type = ordered(rep(c( "BayCANN MAP"), each = n_inputs),
                                                   levels = c("BayCANN MAP")),
                                    value = c(t(map_baycann)))
df_maps_n_true_params


### Plot priors and ANN and IMIS posteriors


df_maps_n_true_params$Parameter<-as.factor(x_names)


library(dampack)

gg_ann_vs_imis <- ggplot(df_samp_prior_post,
                         aes(x = value, y = ..density.., fill = Distribution)) +
  facet_wrap(~Parameter, scales = "free",
             ncol = 5,
             labeller = label_parsed) +
  # geom_vline(data = data.frame(Parameter = as.character(v_names_params_greek),
  #                             value = x_true_data$x, row.names = v_names_params_greek),
  #          aes(xintercept = value)) +
  #geom_vline(data = data.frame(Parameter = as.character(v_names_params_greek),
  #            value = c(t(map_baycann)), row.names = v_names_params_greek),
  #aes(xintercept = value), color = "tomato") +
  geom_vline(data = df_maps_n_true_params,
             aes(xintercept = value, label="MAP"), color = "tomato") +
  #scale_x_continuous(breaks = (5)) +
  scale_x_continuous(breaks = number_ticks(5)) +
  scale_color_manual("", values = c("black", "navy blue", "tomato","green")) +
  geom_density(alpha=0.5) +
  theme_bw(base_size = 16) +
  guides(fill = guide_legend(title = "", order = 1),
         linetype = guide_legend(title = "", order = 2),
         color = guide_legend(title = "", order = 2)) +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.margin=margin(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.background = element_rect(fill = "white",
                                        color = "white"),
        strip.text = element_text(hjust = 0))
gg_ann_vs_imis
ggsave(gg_ann_vs_imis,
       filename = paths_calibration$path_prior_post,
       width = 16, height = 12)

# 13. Internal validation  -------------------------------------------------

###### 13.1 Load Stan file -----
m <- readRDS(path_stan_model) 

#Extraction of parameters 
params <- rstan::extract(m)
Xq <- params$Xq
Xq_df = as.data.frame(Xq)
colnames(Xq_df) <- x_names
Xq_df <- as.matrix(Xq_df)

model <- load_model_hdf5(path_keras_model)

pred_posterior <- model$predict(Xq_df)

pred_posterior <- data.frame(pred_posterior)


colnames(pred_posterior) <- y_names

if (Scale_outputs) {
  pred_posterior_unsc <- unscale_data(pred_posterior, vec.mins = ymins, vec.maxs = ymaxs, vec.means = ymeans, vec.sds = ysds, type = scale_type)
}

path_validation_ANN <- paths_calibration$path_val_ANN
write.csv(pred_posterior_unsc,
          file = path_validation_ANN,
          row.names = FALSE)

# Emulator validation graph
true_target <- read.csv(targets_file)

#Load the model outputs
df_model_outputs <- pred_posterior_unsc


# Summarize the dataframe
summary_df <- data.frame(
  target_names = colnames(df_model_outputs),
  model_mean = apply(df_model_outputs, 2, mean),
  model_UB_95 = apply(df_model_outputs, 2, function(x) quantile(x, 0.975)),
  model_LB_95 = apply(df_model_outputs, 2, function(x) quantile(x, 0.025)),
  model_UB_50 = apply(df_model_outputs, 2, function(x) quantile(x, 0.75)),
  model_LB_50 = apply(df_model_outputs, 2, function(x) quantile(x, 0.25))
)

# Call the function to generate the plots
plots <- graph_internal_validation(model_outputs = summary_df, 
                                   calibration_targets = true_target,
                                   index = "target_names",
                                   categorial_groups = c("Size_P", "Size_D", "Size_R", "StageDist_D","StageDist_P","StageDist_R", "SSP"),
                                   subtype= "lesion_type",
                                   subtype_levels = c("Adenoma", "SSP", "All"),
                                   subtype_color = c("Adenoma" = "#56B4E9", "SSP" = "#E69F00", "All" = "#009E73"),
                                   target_color  = "gray",
                                   model_mean = "model_mean",
                                   model_LB_95 = "model_LB_95",
                                   model_UB_95 = "model_UB_95",
                                   model_LB_50 = "model_LB_50",
                                   model_UB_50 = "model_UB_50",
                                   target_mean = "targets",
                                   target_LB_95 = "stopping_lower_bounds",
                                   target_UB_P5 = "stopping_upper_bounds",
                                   target_groups = "target_groups",
                                   numeric = "age",
                                   categoric = "stage",
                                   plot_title = "SimCRC-R Posterior validation")


#Plot numeric targets
plot1 <- plots$plot_val_num
plot1



#plot categorial targets
plot2 <- plots$plot_val_cat
plot2 



# 14. Save BayCANN parameters  -------------------------------------------------
date_save <- Date_version
##Save BayCANN parameters
path_baycann_params <- paths_calibration$path_baycann_params
param_BayCANN <- list(date_save = date_save ,
                      Model_name = Model_name,               
                      BayCANN_version = BayCANN_version, 
                      scale_type = scale_type,
                      parameters_names = x_names,
                      outputs_names = y_names,
                      verbose = verbose,
                      Selected_targets = Selected_targets,
                      n_batch_size = n_batch_size,
                      n_chains = n_chains,
                      n_epochs = n_epochs,
                      patience = patience,
                      n_hidden_nodes = n_hidden_nodes,
                      n_hidden_layers = n_hidden_layers,
                      activation_fun = activation_fun,
                      n_iter = n_iter,
                      n_thin = n_thin,
                      Normalize_inputs = Normalize_inputs,
                      Normalize_outputs = Normalize_outputs,
                      Scale_inputs = Scale_inputs,
                      Scale_outputs = Scale_outputs,
                      Remove_outliers = Remove_outliers, 
                      Standardize_targets = Standardize_targets,
                      Saved_data = Saved_data,
                      params_file = params_file,
                      outputs_file = outputs_file,
                      targets_file = targets_file,
                      path_keras_model = path_keras_model,
                      t_training = t_training,
                      metric_loss = metric_loss,
                      metric_mae = metric_mae,
                      metric_accuracy = metric_accuracy, 
                      path_posterior = path_posterior, 
                      file_perceptron = file_perceptron,
                      t_calibration = t_calibration,
                      path_stan_model = path_stan_model,
                      path_validation_ANN = path_validation_ANN,
                      path_baycann_params = path_baycann_params
)

save(param_BayCANN, file = path_baycann_params)

if(!is.null(param_BayCANN)) {
  print("BayCANN calibration process finished successfully.\n BayCANN parameters saved successfully.")
} else {
  print("Error: BayCANN parameters not saved.")
}


#**************END OF FILE