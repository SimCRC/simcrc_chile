# *****************************************************************************
#
# Script: 01_wrangle_targets.R
#
# Purpose: Wrangle targets data to compare countries and plot
#
# Author: Jorge Roa
#
# Date Created: 13 Jan 2026
#
# *****************************************************************************
#
# Notes:
#   
#
# *****************************************************************************



# 01 Load libraries ----------------------------------------------------------


library(tidyverse)
library(readxl)

# 02 Load and wrangle data ---------------------------------------------------

df_CH_targets <- read_xlsx("data-raw/true_target_simcrcRvCH.xlsx")

df_simcrc_targets <- read_csv("data-raw/20220909_simcrc_targets.csv")

df_targets_US_CH <- bind_rows(df_simcrc_targets %>% mutate(country = "USA"), 
                              df_CH_targets%>% mutate(country = "Chile"))

# Add names

df_targets_US_CH <- df_targets_US_CH %>%
  mutate(names = case_when(
    target_groups == "CRCInc_Overall" ~ "CRC Incidence Overall",
    target_groups == "CRCInc_P" ~ "CRC Incidence Proximal",
    target_groups == "CRCInc_D" ~ "CRC Incidence Distal",
    target_groups == "CRCInc_R" ~ "CRC Incidence Rectal",
    target_groups == "Prev0Ad" ~ "Prevalence 0 Adenomas",
    target_groups == "Prev1Ad" ~ "Prevalence 1 Adenoma",
    target_groups == "Prev2Ad" ~ "Prevalence 2 Adenomas",
    target_groups == "Prev3PlusAd" ~ "Prevalence 3+ Adenomas",
    target_groups == "PrevPreclin" ~ "Prevalence Preclinical CRC",
    target_groups == "Size_D" ~ "Adenoma Size Distal",
    target_groups == "Size_P" ~ "Adenoma Size Proximal",
    target_groups == "Size_R" ~ "Adenoma Size Rectal",
    target_groups == "StageDist_D" ~ "CRC Stage Distal",
    target_groups == "StageDist_P" ~ "CRC Stage Proximal",
    target_groups == "StageDist_R" ~ "CRC Stage Rectal",
    TRUE ~ target_groups
  ))

ggplot(df_targets_US_CH %>% filter(is.na(stage)), aes(x = age, y = targets, group = country, color = country)) +
  geom_ribbon(aes(ymin = stopping_lower_bounds, ymax = stopping_upper_bounds, fill = country), alpha = 0.2, color = NA) +
  geom_point() +
  geom_line() +
  facet_wrap(names~target_groups, scales = "free_y") +
  labs(title = "Numerical targets from SimCRC and Chile data",
       x = "\nAge\n",
       y = "\nCRC cases per 100k\n") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 15),
        legend.title = element_blank(),
        axis.text = element_text(size = 15),
        strip.text = element_text(size = 14)
  )

ggsave("figs/SimCRC_CH_NUM_targets.png", width = 15, height = 9, dpi = 300)


#Plot Categorical targets

cat_data <- df_targets_US_CH %>% filter(!is.na(stage))

unique(cat_data$stage)

ggplot(df_targets_US_CH %>% filter(!is.na(stage)), aes(x = stage, y = targets, group = country, color = country, linetype = stage)) +
  geom_point() +
  #add error bars
  geom_errorbar(aes(ymin = stopping_lower_bounds, ymax = stopping_upper_bounds), width = 0.2) +
  facet_wrap(names~target_groups, scales = "free_y") +
  labs(title = "Categorical targets from SimCRC and Chile data",
       x = "\nAge\n",
       y = "\nProportion (%)\n") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 15),
        legend.title = element_blank(),
        axis.text = element_text(size = 15),
        strip.text = element_text(size = 14)
  )

library(dplyr)
library(ggplot2)

cat_data <- df_targets_US_CH %>%
  filter(!is.na(stage)) %>%
  mutate(
    stage = as.character(stage),
    stage_type = if_else(stage %in% c("LR","MR","HR"), "Risk group", "Cancer stage"),
    stage_x = factor(
      stage,
      levels = c("LR","MR","HR","1","2","3","4")
    )
  )

ggplot(
  cat_data,
  aes(x = stage_x, y = targets, group = country, color = country)
) +
  geom_point(position = position_dodge(width = 0.35), size = 2) +
  geom_errorbar(
    aes(ymin = stopping_lower_bounds, ymax = stopping_upper_bounds),
    width = 0.2,
    position = position_dodge(width = 0.35)
  ) +
  facet_grid(stage_type ~ names + target_groups, scales = "free_y", space = "free_x") +
  labs(
    title = "Categorical targets from SimCRC and Chile data",
    x = "\nCategory\n",
    y = "\nProportion (%)\n"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 15),
    legend.title = element_blank(),
    axis.text = element_text(size = 15),
    strip.text = element_text(size = 14)
  )

library(dplyr)
library(ggplot2)
library(patchwork)

cat_data <- df_targets_US_CH %>%
  filter(!is.na(stage)) %>%
  mutate(
    stage = as.character(stage),
    stage = factor(stage, levels = c("LR","MR","HR","1","2","3","4")),
    names = as.character(names)
  )

target_names <- sort(unique(cat_data$names))

plots <- vector("list", length(target_names))

for (i in seq_along(target_names)) {
  
  nm <- target_names[[i]]
  d_i <- cat_data %>% filter(names == nm)
  
  plots[[i]] <-
    ggplot(d_i, aes(x = stage, y = targets, group = country, color = country)) +
    geom_point(position = position_dodge(width = 0.35), size = 2) +
    geom_errorbar(
      aes(ymin = stopping_lower_bounds, ymax = stopping_upper_bounds),
      width = 0.2,
      position = position_dodge(width = 0.35)
    ) +
    labs(
      title =  nm,
      subtitle = unique(d_i$target_groups),
      x = "\nStage\n",
      y = "\nProportion (%)\n",
      color = NULL
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",   # keep legends so patchwork can collect
      legend.text = element_text(size = 12),
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15),
      plot.title  = element_text(size = 12)
    )
}


  wrap_plots(plots, ncol = 3) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Categorical targets from SimCRC and Chile data") &
  theme(legend.position = "bottom")

  ggsave("figs/SimCRC_CH_CAT_targets.png", width = 15, height = 9, dpi = 300)





