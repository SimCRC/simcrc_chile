

true_target <- read.csv(paths_calibration$path_targets)


#Select non-missing outputs
v_NAs <- apply(df_simcrc_outputs, 1, function(x) any(is.na(x)))
table(v_NAs)
df_model_outputs <- df_simcrc_outputs[!v_NAs, ]

# Summarize the dataframe
summary_df <- data.frame(
  target_names = colnames(df_model_outputs),
  model_mean = apply(df_model_outputs, 2, mean),
  model_UB_95 = apply(df_model_outputs, 2, function(x) quantile(x, 0.975)),
  model_LB_95 = apply(df_model_outputs, 2, function(x) quantile(x, 0.025)),
  model_UB_50 = apply(df_model_outputs, 2, function(x) quantile(x, 0.75)),
  model_LB_50 = apply(df_model_outputs, 2, function(x) quantile(x, 0.25))
)
model_targets <- inner_join(summary_df , true_target ,by=c("target_names"))  #Adding group


model_targets$lesion_type <- factor(model_targets$lesion_type, c("Adenoma", "SSP", "All"))

# Improve target group names for plotting
model_targets[ model_targets$target_groups=="Prev0Ad" , ]$target_groups  <- "Prevalence 0 lesions"
model_targets[ model_targets$target_groups=="Prev1Ad" , ]$target_groups  <- "Prevalence 1 lesion"
model_targets[ model_targets$target_groups=="Prev2Ad" , ]$target_groups  <- "Prevalence 2 lesions"
model_targets[ model_targets$target_groups=="Prev3PlusAd" , ]$target_groups  <- "Prevalance 3+ lesions"


#Set the order we want to have on the facets
model_targets$target_groups<- factor(model_targets$target_groups, c("CRCInc_Overall",
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
#Uncomment in case of two models
#model_targets$n_target <- rep(1:110, 2)

#identification of categorical variables

cat_groups <- c("Size_P", "Size_D", "Size_R", "StageDist_D","StageDist_P","StageDist_R", "SSP")
model_targets$categorical <- ifelse(model_targets$target_groups %in% cat_groups,1,0)


plot_val_num <- ggplot(data = model_targets[model_targets$categorical==0,], 
                       aes(x    = age, 
                           y    = targets, 
                           ymin = stopping_lower_bounds, 
                           ymax = stopping_upper_bounds,
                           color=lesion_type))+ 
  geom_errorbar(width=.4, size=0.9) +
  theme(legend.position="none") +
  geom_ribbon(data = model_targets[model_targets$categorical==0,],
              aes(x    = age,
                  y    = model_mean,
                  ymin = model_LB_95,
                  ymax = model_UB_95,
                  fill = lesion_type),
              alpha = 0.3) +
  geom_ribbon(data = model_targets[model_targets$categorical==0,],
              aes(x    = age,
                  y    = model_mean,
                  ymin = model_LB_50,
                  ymax = model_UB_50,
                  fill = lesion_type),
              alpha = 0.5) +
  facet_wrap(~ target_groups + lesion_type,scales="free", ncol = 4) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(), legend.position="none") +
  scale_color_manual(values = color_values)+
  scale_fill_manual(values =color_values)+
  #scale_y_continuous(breaks= 0:.2, labels=c(0,1))+
  scale_y_continuous(breaks = number_ticks(5))+
  #scale_x_continuous(breaks= 1:lim, labels=x_label) +
  theme_bw(base_size = 8) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 8, angle = 90),
        axis.title = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0)) +
  labs(title = "SimCRC-R coverage analysis", 
       x     = "")

plot_val_num

ggsave(plot_val_num,
       filename = paths_calibration$path_coverage_num,
       width = 10, height = 6)


# First, create numeric positions for your groups
model_targets$x_pos <- as.numeric(factor(model_targets$stage))
shift <- 0.3

plot_val_cat <- ggplot(data = model_targets[model_targets$categorical==1,], 
                       aes(x    = x_pos, 
                           y    = targets, 
                           ymin = stopping_lower_bounds, 
                           ymax = stopping_upper_bounds,
                           color=lesion_type, shape = "target"))+ 
  geom_point()+
  geom_errorbar(width=.2, size=1, alpha=0.5) +
  theme(legend.position="none") +
  geom_errorbar(data = model_targets[model_targets$categorical==1,],
                aes(x    = x_pos + shift,
                    y    = model_mean,
                    ymin = model_LB_95,
                    ymax = model_UB_95, color= lesion_type), width=.2, size=0.9, alpha = 1) +
  geom_point(aes(x    = x_pos + shift,
                 y    = model_mean , shape="model"))+
  facet_wrap(~ target_groups + lesion_type,scales="free", ncol = 3) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(), legend.position="none") +
  scale_fill_manual(values = color_values)+
  scale_color_manual(values = color_values)+
  scale_shape_manual(values = c(8 , 16)) +
  #scale_y_continuous(breaks= 0:.2, labels=c(0,1))+
  scale_y_continuous(breaks = number_ticks(5))+
  # Adjust x-axis labels back to original group labels
  scale_x_continuous(breaks = model_targets$x_pos, labels = model_targets$stage) +
  theme_bw(base_size = 8) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 8, angle = 0),
        axis.title = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0)) +
  labs(title = "SimCRC-R LHS coverage analysis ", 
       x     = "")

plot_val_cat

ggsave(plot_val_cat,
       filename = paths_calibration$path_coverage_cat,
       width = 10, height = 6)


