
# Posterior validation graphs  --------------------------------------------

chains_str <- paste(chains_to_include, collapse = "-")

data_outputs_simcrc  <- df_simcrc_outputs[df_simcrc_outputs$chain %in% chains_to_include,]  #filter only chains selected

#data_outputs_simcrc$index <- "Model R"

#remove observations with missing values

v_NAs <- apply(data_outputs_simcrc, 1, function(x) any(is.na(x)))
table(v_NAs)

data_outputs_simcrc  <- data_outputs_simcrc[!v_NAs,]

# -------------------------------------------------------------------------
# Example of use ----------------------------------------------------------
# -------------------------------------------------------------------------

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)


#Load the targets
true_target <- read.csv(targets_file)

#Load the model outputs
df_model_outputs <- data_outputs_simcrc

#drop the id_draw column
df_model_outputs <- df_model_outputs[, !colnames(df_model_outputs) %in% c("id_draw")]

df_model_outputs <- as.data.frame(df_model_outputs)

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
                                   plot_title = paste0("SimCRC-R Posterior validation ", "(Chains: ", chains_str, ")"))


#Plot numeric targets
plot1 <- plots$plot_val_num
plot1

ggsave(plot1,
               filename = paste0(folder, "/fig_internal_validation_", BayCANN_version, "_num.png"),
               width = 10, height = 8)

#plot categorial targets
plot2 <- plots$plot_val_cat
plot2 

ggsave(plot2,
       filename = paste0(folder, "/fig_internal_validation_", BayCANN_version, "_cat.png"),
       width = 10, height = 8)




