
# Model outputs validation function
# This function generates plots for coverage analysis and internal validation of model outputs against true targets.
 

# Function to create a sequence of numbers for axis ticks
number_ticks <- function(n) {
  function(limits) pretty(seq(limits[1], limits[2], length.out = n))
}


#' Function to generate a graph that compare the model outputs against the true targets
#' @param model_outputs Data frame containing the model outputs
#' @param calibration_targets Data frame containing the true targets for comparison
#' @param index Column name to join the two data frames
#' @param categorial_groups Vector of categorical groups to be used in the analysis
#' @param subtype Column name for the subtype of lesions, e.g. "lesion_type"
#' @param subtype_levels Vector of levels for the subtype factor
#' @param subtype_color Named vector of colors for each subtype
#' @param target_color Color for the target points
#' @param model_mean Column name for the model mean values
#' @param model_LB_95 Column name for the lower bound of the 95% confidence interval
#' @param model_UB_95 Column name for the upper bound of the 95% confidence interval
#' @param model_LB_50 Column name for the lower bound of the 50% confidence interval
#' @param model_UB_50 Column name for the upper bound of the 50% confidence interval
#' @param target_mean Column name for the target mean values
#' @param target_LB_95 Column name for the lower bound of the 95% confidence interval for targets
#' @param target_UB_P5 Column name for the upper bound of the 95% confidence interval for targets
#' @param target_groups Column name for the target groups
#' @param numeric Column name for the numeric variable (e.g., age)
#' @param categoric Column name for the categorical variable (e.g., stage)
#' @param plot_title Title for the plot
#' @return A list containing two ggplot objects: one for numeric targets and one for categorical targets

graph_internal_validation <- function(model_outputs, 
                                      calibration_targets,
                                      index = "target_names",
                                      categorial_groups = c("Size_P", "Size_D", "Size_R", "StageDist_D","StageDist_P","StageDist_R", "SSP"),
                                      subtype= "lesion_type",
                                      subtype_levels = c("Adenoma", "SSP", "All"),
                                      subtype_color = c("Adenoma" = "#56B4E9", "SSP" = "#E69F00", "All" = "#009E73"),
                                      target_color  = "#D56E40",
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
                                      plot_title = "SimCRC-R coverage analysis") {
  
  # Ensure the input data frames are not empty
  if (nrow(model_outputs) == 0 || nrow(calibration_targets) == 0) {
    stop("Input data frames cannot be empty.")
  }
  # Check if the specified columns exist in the data frames
  required_model_columns <- c(index, model_mean, model_LB_95, model_UB_95, model_LB_50, model_UB_50)
  if (!all(required_model_columns %in% colnames(model_outputs))) {
    stop(paste("The following columns are missing in model outputs dataframe:", 
               paste(setdiff(required_model_columns, colnames(model_outputs)), collapse = ", ")))
  }
  
  required_target_columns <- c(index, target_mean, target_LB_95, target_UB_P5, target_groups, numeric, categoric, subtype)
  if (!all(required_target_columns %in% colnames(calibration_targets))) {
    stop(paste("The following columns are missing in calibration targets dataframe:", 
               paste(setdiff(required_target_columns, colnames(calibration_targets)), collapse = ", ")))
  }
  
  # Rename columns for clarity
  #rename variables in targets names
  colnames(calibration_targets)[colnames(calibration_targets) == target_mean] <- "target_mean"
  colnames(calibration_targets)[colnames(calibration_targets) == target_LB_95] <- "target_LB_95"
  colnames(calibration_targets)[colnames(calibration_targets) == target_UB_P5] <- "target_UB_P5"
  colnames(calibration_targets)[colnames(calibration_targets) == target_groups] <- "target_groups"
  colnames(calibration_targets)[colnames(calibration_targets) == numeric] <- "numeric"
  colnames(calibration_targets)[colnames(calibration_targets) == categoric] <- "categoric"
  colnames(calibration_targets)[colnames(calibration_targets) == subtype] <- "subtype"
  
  # Rename columns in model outputs
  colnames(model_outputs)[colnames(model_outputs) == model_mean] <- "model_mean"
  colnames(model_outputs)[colnames(model_outputs) == model_LB_95] <- "model_LB_95"
  colnames(model_outputs)[colnames(model_outputs) == model_UB_95] <- "model_UB_95"
  colnames(model_outputs)[colnames(model_outputs) == model_LB_50] <- "model_LB_50"
  colnames(model_outputs)[colnames(model_outputs) == model_UB_50] <- "model_UB_50"

  # Ensure the lesion_type is a factor with specified levels
  calibration_targets$subtype <- factor(calibration_targets$subtype, levels = subtype_levels)
  
  # Join the model outputs with the true targets
  model_targets <- inner_join(model_outputs, calibration_targets, by = index)
  
  # Check if the join was successful
  if (nrow(model_targets) == 0) {
    stop("No matching rows found between model outputs and calibration targets.")
  }
  
  # Keep only the necessary columns
  model_targets <- model_targets %>%
    dplyr::select(index, model_mean, model_LB_95, model_UB_95, model_LB_50, model_UB_50,
           target_mean, target_LB_95, target_UB_P5, target_groups, numeric, categoric, subtype)
  
  model_targets$categorical <- ifelse(model_targets$target_groups %in% categorial_groups,1,0) 
  
  #Create plot with numeric outputs
  
  plot_val_num <- ggplot(data = model_targets[model_targets$categorical==0,], 
                         aes(x    = numeric, 
                             y    = target_mean, 
                             ymin = target_LB_95, 
                             ymax = target_UB_P5,
                             color= subtype))+
    geom_errorbar(width=1.5, size=0.9, alpha = 1, color=target_color) +
    geom_point(aes(x    = numeric,
                   y    = target_mean , color=target_color, shape="target"))+
    theme(legend.position="none") +
    geom_ribbon(data = model_targets[model_targets$categorical==0,],
                aes(x    = numeric,
                    y    = model_mean,
                    ymin = model_LB_95,
                    ymax = model_UB_95,
                    fill = subtype),
                alpha = 0.3) +
    geom_ribbon(data = model_targets[model_targets$categorical==0,],
                aes(x    = numeric,
                    y    = model_mean,
                    ymin = model_LB_50,
                    ymax = model_UB_50,
                    fill = subtype),
                alpha = 0.4) +
    facet_wrap(~ target_groups + subtype,scales="free", ncol = 4) +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(), legend.position="none") +
    scale_color_manual(values = subtype_color)+
    scale_fill_manual(values =subtype_color)+
    scale_shape_manual(values = c("target" = 16 , "model" = 8)) +
    #scale_y_continuous(breaks= 0:.2, labels=c(0,1))+
    scale_y_continuous(breaks = number_ticks(5))+
    #scale_x_continuous(breaks= 1:lim, labels=x_label) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(size = 16, face = "bold"),
          axis.text.x = element_text(size = 12, angle = 90),
          axis.title = element_text(size = 14),
          panel.grid.major = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0)) +
    labs(title = plot_title, 
         x     = "")
  
  
  #Create plot with categorical outputs
  # First, create numeric positions for your groups
  model_targets$x_pos <- as.numeric(factor(model_targets$categoric))
  shift <- 0.2
  
  # Create the plot for categorical outputs
  plot_val_cat <- ggplot(data = model_targets[model_targets$categorical==1,], 
                         aes(x    = x_pos - shift, 
                             y    = target_mean, 
                             ymin = target_LB_95, 
                             ymax = target_UB_P5,
                             color=subtype))+ 
    geom_errorbar(width=.2, size=1, alpha=0.5, color=target_color) +
    geom_point(aes(x    = x_pos - shift,
                   y    = target_mean , shape="target", color=target_color))+
    theme(legend.position="none") +
    geom_errorbar(data = model_targets[model_targets$categorical==1,],
                  aes(x    = x_pos + shift,
                      y    = model_mean,
                      ymin = model_LB_95,
                      ymax = model_UB_95, color= subtype), width=.2, size=0.9, alpha = 1) +
    geom_point(aes(x    = x_pos + shift,
                   y    = model_mean , shape="model"))+
    geom_rect(aes(xmin = x_pos + shift/2, xmax = x_pos + shift*3/2, 
                     ymin = model_LB_50, ymax = model_UB_50, fill = subtype), 
               alpha = 0.3) +
    
    facet_wrap(~ target_groups + subtype,scales="free", ncol = 3) +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(), legend.position="none") +
    scale_fill_manual(values = subtype_color)+
    scale_color_manual(values = subtype_color)+
    scale_shape_manual(values = c("target" = 16 , "model" = 8)) +
    #scale_y_continuous(breaks= 0:.2, labels=c(0,1))+
    scale_y_continuous(breaks = number_ticks(5))+
    # Adjust x-axis labels back to original group labels
    scale_x_continuous(breaks = model_targets$x_pos, labels = model_targets$categoric) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(size = 16, face = "bold"),
          axis.text.x = element_text(size = 12, angle = 0),
          axis.title = element_text(size = 14),
          panel.grid.major = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0)) +
    labs(title   = plot_title, 
         x       ="")
  
  return(list(plot_val_num = plot_val_num, plot_val_cat = plot_val_cat))
  
}

