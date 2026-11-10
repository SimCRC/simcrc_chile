# Chile vs USA cheat sheet, adapted from the original cheat_sheet_chile.R.
# Keeps the original simulation, outcome definitions, country life tables,
# colours, axes and layout. Uses the latest v0.15.0 Min_AbsolutErr sets.
# Run from the simcrc_chile root:
#   Rscript analysis/figures/11_cheat_sheet_Chile_vs_USA.R
# Inputs and outputs are in this repository; no sibling checkout is required.

path <- "outputs/comparison/USA_20261026_Chile_20260921/cheat_sheet/"
dir.create(path, recursive = TRUE, showWarnings = FALSE)
gc()

# *****************************************************************************
#### Load packages ####
# *****************************************************************************

# Load Libraries -------------------------------------------------------

library(data.table)
library(ggplot2)
library(ggtext)
library(tidyr)
library(dplyr)
library(simcrc)
library(scales)
library(patchwork)

setDTthreads(threads = 6)


# Generate SimCRC-R data --------------------------------------------------


path <- "outputs/comparison/USA_20261026_Chile_20260921/cheat_sheet/"
n_pop <- 1e6

# 1. Create population -------------------------------------------------
## data.table ----
set.seed(20220906)


#unit test idea.time to 6 filtrado. desmembrado dwell time. 
#where crc_dx_i == 1 , means they have dwell time. dwell_time_small = Age6mm - AgeLesionOnset. dwell_time_medium = Age10mm - Age6mm. dwell_time_large = AgePreclin - Age10mm.



# Define the simulation population size
n_pop <- 1e6 # Run at least 1 mil for publications, 10mil if possible for stable results

# Define the cohort age
cohort_age <- 40

# Define model-specific inputs
df_lt_chile <- read.csv("data-raw/df_lifetable_2017_CH.csv")
colnames(df_lt_chile)[colnames(df_lt_chile) == "age"] <- "Age"
colnames(df_lt_chile)[colnames(df_lt_chile) == "mortality_rate"] <- "mortality.rates"

model_specs <- list(
  "R v0.15.0-C MAE (Chile)" = list(
    param_file = "outputs/BayCANN_versions/Chile/Adenoma/F/v0.15.0/v0.15.0.20260921.1018/l_params_calibrated_sets_SimCRC_v0.15.0.20260921.1018_Adenoma_F.RData",
    dt_life_table_F = df_lt_chile,
    metric = "Min_AbsolutErr"
  ),
  "R v0.15.0-C MAE (USA)" = list(
    param_file = "data-raw/l_params_calibrated_sets_SimCRC_v0.15.0.20261026.130537_Adenoma_F.RData",
    dt_life_table_F = simcrc::l_life_tables$dt_LFT1977_females_0,
    metric = "Min_AbsolutErr"
  )
)

# Initialize results list
simcrc_v_nthist <- list()

# Loop through each version
for (v in names(model_specs)) {
  
  l_params_filename <- model_specs[[v]]$param_file

  dt_pop <- simcrc::get_dt_population(
    year = 1980,
    byear = 1980,
    p_female = 1,
    p_white = 0.8,
    n_pop = n_pop,
    dt_life_table_F = model_specs[[v]]$dt_life_table_F,
    dt_life_table_M = NULL
  )
  
  #if R v0.9.0.5-C, then modify this scale l_params[["AdOnset_Weibull_Scale"]] to 1.0e-08
  
  # Load model-specific parameters from each file
  param_env <- new.env(parent = emptyenv())
  loaded_objects <- load(l_params_filename, envir = param_env)

  metric <- model_specs[[v]]$metric

  if ("l_params_calibrated_sets" %in% loaded_objects) {
    l_params <- load_params_init(
      fromFile = TRUE,
      filename = param_env$l_params_calibrated_sets[[metric]]
    )
  } else if ("l_params_Min_AbsolutErr" %in% loaded_objects) {
    l_params <- load_params_init(
      fromFile = TRUE,
      filename = param_env$l_params_Min_AbsolutErr
    )
  } else {
    stop(
      paste0(
        "No supported parameter object found in file: ",
        l_params_filename,
        ". Expected `l_params_calibrated_sets` or `l_params_Min_AbsolutErr`."
      )
    )
  }

  # Update defaults consistently across versions
  # Apply Chile-specific override only to Chile model
  if (v == "R v0.15.0-C MAE (Chile)") {
    l_params$min_age_lesion_onset <- 10
  } else {
    l_params$min_age_lesion_onset <- 10
  }


  l_params$mort_by_race <- FALSE
  # l_params$year_surv_improv <- 2003    # We haven't adjusted this for Chile
  #drop $AdOnsetPropensity_Gamma_Alpha, AdNaturalHistoryPropensity_PreclinCancerProg_TampingTerm, AdNaturalHistoryPropensity_CancerOnset_TampingTerm
  
  
  # Run simulation
  sim_time <- system.time({
    result <- simcr_nathist_ssp_DES(
      l_params_all = l_params,
      dt_pop = dt_pop,
      SSP_pathway = FALSE
    )
  })
  
  # Store result
  simcrc_v_nthist[[v]] <- list(
    version = v,
    result = result,
    time = sim_time
  )
}

names(simcrc_v_nthist)



lesion_vars <- list(
  age_to_onset    = "age_lesion_onset",
  age_to_6mm      = "age_to_6mm",
  age_to_10mm     = "age_to_10mm",
  age_preclin     = "age_crc",
  age_clin        = "age_dx_symp",
  age_death_crc   = "age_death_crc"
)

l_df_age_r_long <- list()



#filter this names(simcrc_v_nthist) to version v0.7.0.1-C and v0.7.0-NC

# Loop through each version







for (version in names(simcrc_v_nthist)) {
  
  dt_crc_pop <- simcrc_v_nthist[[version]]$result$dt_crc_pop
  
  # Fix id
  if (!"id" %in% names(dt_crc_pop)) {
    dt_crc_pop[, id := seq_len(.N)]
  } else {
    dt_crc_pop[, id := zoo::na.locf(id, na.rm = FALSE)]
  }
  
  # Fix id_lesion
  if (!"id_lesion" %in% names(dt_crc_pop) || all(is.na(dt_crc_pop$id_lesion))) {
    dt_crc_pop[, id_lesion := 1L]
  } else {
    dt_crc_pop[, id_lesion := zoo::na.locf(id_lesion, na.rm = FALSE)]
  }
  
  # Create n_lesions: total lesions per person (including 0)
  dt_crc_pop[, n_lesions := sum(develop_lesion, na.rm = TRUE), by = id]
  dt_n_ade <- unique(dt_crc_pop[, .(id, n_lesions)])
  
  # Create output dataset with zeros included
  l_df_age_r_long[[paste0("R_", version, "_n_adenomas")]] <- data.table(
    age_raw = dt_n_ade$n_lesions,
    age = dt_n_ade$n_lesions,
    total = 1,
    type = paste0("R_", version),
    stage = "n_adenomas"
  )
  
  # dwell_time: age_dx_symp - age_lesion_onset
  dwell_values <- dt_crc_pop[
    crc_dx_i == 1 & is_crc_dx == 1 & !is.na(age_dx_symp) & !is.na(age_lesion_onset),
    age_dx_symp - age_lesion_onset
  ]
  valid_dwell_idx <- which(dwell_values >= 0 & !is.na(dwell_values))
  
  l_df_age_r_long[[paste0(version, "_dwell_time")]] <- data.table(
    age_raw = dwell_values[valid_dwell_idx],
    age = floor(dwell_values[valid_dwell_idx]),
    total = 1,
    type = paste0("R_", version),
    stage = "dwell_time"
  )
  
  # sojourn_time: age_dx_symp - age_crc
  sojourn_values <- dt_crc_pop[
    crc_dx_i == 1 & is_crc_dx == 1 & !is.na(age_dx_symp) & !is.na(age_crc),
    age_dx_symp - age_crc
  ]
  valid_soj_idx <- which(sojourn_values >= 0 & !is.na(sojourn_values))
  
  l_df_age_r_long[[paste0(version, "_sojourn_time")]] <- data.table(
    age_raw = sojourn_values[valid_soj_idx],
    age = floor(sojourn_values[valid_soj_idx]),
    total = 1,
    type = paste0("R_", version),
    stage = "sojourn_time"
  )
  
  # regular lesion stages
  age_onset_vec <- dt_crc_pop$age_lesion_onset
  age_6mm_vec   <- dt_crc_pop$age_to_6mm
  
  for (stage in names(lesion_vars)) {
    var <- lesion_vars[[stage]]
    
    values <- switch(stage,
                     "age_clin" = dt_crc_pop[
                       crc_dx_i == 1 & is_crc_dx == 1 & !is.na(get(var)), 
                       get(var)
                     ],
                     "age_death_crc" = dt_crc_pop[
                       crc_death_i == 1 & death_cause == "crc" & !is.na(get(var)),
                       get(var)
                     ],
                     dt_crc_pop[
                       get(var) < age_death & !is.na(get(var)),
                       get(var)
                     ]
    )
    
    indices <- switch(stage,
                      "age_clin" = which(dt_crc_pop$crc_dx_i == 1 & dt_crc_pop$is_crc_dx == 1 & !is.na(dt_crc_pop[[var]])),
                      "age_death_crc" = which(dt_crc_pop$crc_death_i == 1 & dt_crc_pop$death_cause == "crc" & !is.na(dt_crc_pop[[var]])),
                      which(dt_crc_pop[[var]] < dt_crc_pop$age_death & !is.na(dt_crc_pop[[var]]))
    )
    
    # Add the original stage entry
    l_df_age_r_long[[paste0(version, "_", stage)]] <- data.table(
      age_raw = values,
      age = floor(values),
      total = 1,
      type = paste0("R_", version),
      stage = stage
    )
    
    # Derived time stages
    if (stage == "age_to_6mm") {
      time_to_6mm_values <- values - age_onset_vec[indices]
      l_df_age_r_long[[paste0(version, "_time_to_6mm")]] <- data.table(
        age_raw = time_to_6mm_values,
        age = floor(time_to_6mm_values),
        total = 1,
        type = paste0("R_", version),
        stage = "time_to_6mm"
      )
    }
    
    if (stage == "age_to_10mm") {
      time_to_10mm_values <- values - age_6mm_vec[indices]
      l_df_age_r_long[[paste0(version, "_time_to_10mm")]] <- data.table(
        age_raw = time_to_10mm_values,
        age = floor(time_to_10mm_values),
        total = 1,
        type = paste0("R_", version),
        stage = "time_to_10mm"
      )
    }
    if (stage == "age_preclin") {
      age_preclin_vec <- values  # same as dt_crc_pop[[var]] at indices
      
      #Time from 6mm to Preclinical
      time_6mm_to_preclin <- age_preclin_vec - age_6mm_vec[indices]
      valid_6mm_idx <- which(!is.na(time_6mm_to_preclin) & time_6mm_to_preclin >= 0)
      
      l_df_age_r_long[[paste0(version, "_time_from_6mm_to_preclin")]] <- data.table(
        age_raw = time_6mm_to_preclin[valid_6mm_idx],
        age = floor(time_6mm_to_preclin[valid_6mm_idx]),
        total = 1,
        type = paste0("R_", version),
        stage = "time_from_6mm_to_preclin"
      )
      
      #Time from 10mm to Preclinical
      time_10mm_to_preclin <- age_preclin_vec - dt_crc_pop$age_to_10mm[indices]
      valid_10mm_idx <- which(!is.na(time_10mm_to_preclin) & time_10mm_to_preclin >= 0)
      
      l_df_age_r_long[[paste0(version, "_time_from_10mm_to_preclin")]] <- data.table(
        age_raw = time_10mm_to_preclin[valid_10mm_idx],
        age = floor(time_10mm_to_preclin[valid_10mm_idx]),
        total = 1,
        type = paste0("R_", version),
        stage = "time_from_10mm_to_preclin"
      )
    }
  }
  
  valid_idx <- with(dt_crc_pop,
                    crc_dx_i == 1)
  
  dt_valid <- dt_crc_pop[valid_idx]
  
  #Small
  dw_sm <- dt_valid[age_to_6mm >= age_lesion_onset]
  l_df_age_r_long[[paste0(version, "_dw_t_sm")]] <- data.table(
    age_raw = dw_sm$age_to_6mm - dw_sm$age_lesion_onset,
    age = floor(dw_sm$age_to_6mm - dw_sm$age_lesion_onset),
    total = 1,
    type = paste0("R_", version),
    stage = "dw_t_sm"
  )
  
  #Medium
  dw_md <- dt_valid[age_to_10mm >= age_to_6mm]
  l_df_age_r_long[[paste0(version, "_dw_t_md")]] <- data.table(
    age_raw = dw_md$age_to_10mm - dw_md$age_to_6mm,
    age = floor(dw_md$age_to_10mm - dw_md$age_to_6mm),
    total = 1,
    type = paste0("R_", version),
    stage = "dw_t_md"
  )
  
  #Large 
  dw_lg <- dt_valid[age_crc >= age_to_10mm]
  l_df_age_r_long[[paste0(version, "_dw_t_lg")]] <- data.table(
    age_raw = dw_lg$age_crc - dw_lg$age_to_10mm,
    age = floor(dw_lg$age_crc - dw_lg$age_to_10mm),
    total = 1,
    type = paste0("R_", version),
    stage = "dw_t_lg"
  )
  
}


df_r_long <- rbindlist(l_df_age_r_long, fill = TRUE)



df_plot_data <- df_r_long[
  , .(n_ade = sum(total)), by = .(type, stage, age)
]

#drop NA for age

df_plot_data <- df_plot_data[!is.na(age)]
#rename R versions to SimCRC v0.7.0 and so on

unique(df_plot_data$type)

df_plot_data[, 
  type := gsub("R_R v0.15.0-C MAE (Chile)", "R v0.15.0-C MAE (Chile)", type, fixed = TRUE)
]
df_plot_data[, 
  type := gsub("R_R v0.15.0-C MAE (USA)", "R v0.15.0-C MAE (USA)", type, fixed = TRUE)
]
df_plot_data[, type := as.character(type)]

df_plot_data[, type := gsub("\\s+", " ", type)]


unique(df_plot_data$type)

df_plot_data[, type := factor(
  type,
  levels = c("R v0.15.0-C MAE (Chile)", "R v0.15.0-C MAE (USA)"))]

#add total of adenomas to the data frame. add to stage n_adenomas_total

df_plot_data <- df_plot_data %>%
  mutate(
    stage = recode(stage,
                   "n_adenomas" = "n_adenomas_total")
  )
saveRDS(df_plot_data, file.path(path, "natural_history_plot_data.rds"))

# Define your palette
palette <- c(
  "R v0.15.0-C MAE (Chile)" = "#F94144",
  "R v0.15.0-C MAE (USA)" = "#F8961E"
)


# Load necessary package
library(data.table)

#ombine df_r_long








# List of all stages to plot individually
all_stages <- unique(df_plot_data$stage)



custom_stage_order <- c(
  "age_to_onset",      # Age lesion onset
  "time_to_6mm",       # Time to 6mm
  "age_to_6mm",        # Age to 6mm
  "time_to_10mm",      # Time to 10mm
  "age_to_10mm",       # Age to 10mm
  "time_from_6mm_to_preclin", # Time from 6mm to preclinical
  "time_from_10mm_to_preclin", # Time from 10mm to preclinical
  "age_preclin",       # Age to preclinical
  "age_clin",          # Age to clinical
  "age_death_crc",     # Age to death
  "dwell_time",        # Dwell time
  "sojourn_time",       # Sojourn
  "dw_t_sm",         # Dwell time small
  "dw_t_md",        # Dwell time medium
  "dw_t_lg"       # Dwell time large
)

stage_title_map <- c(
  age_to_onset                 = "Age to Onset",
  time_to_6mm                  = "Time to 6mm",
  age_to_6mm                   = "Age to 6mm",
  time_to_10mm                 = "Ade(s) time from\n6 mm to 10 mm",
  age_to_10mm                  = "Age to 10mm",
  time_from_6mm_to_preclin     = "Time from 6mm to Preclinical",
  time_from_10mm_to_preclin    = "Time from 10mm to Preclinical",
  age_preclin                  = "Age Preclinical",
  age_clin                     = "Age Clinical",
  age_death_crc                = "Age at CRC Death",
  dwell_time                   = "Dwell Time",
  sojourn_time                 = "Sojourn Time",
  dw_t_sm                      = "Dwell Time Small",
  dw_t_md                      = "Dwell Time Medium",
  dw_t_lg                      = "Dwell Time Large"
)


available_stages <- unique(df_plot_data$stage)
all_stages <- custom_stage_order[custom_stage_order %in% available_stages]

# Create legend labels using <span style="color:...">text</span>
legend_labels <- sapply(names(palette), function(x) {
  paste0("<span style='color:", palette[[x]], "'>", x, "</span>")
})
names(legend_labels) <- names(palette)

library(ggtext)

plot_list <- lapply(all_stages, function(stg) {
  
  #stg = "dw_t_md"
  
  df_stage <- df_plot_data %>%
    filter(stage == stg)
  stage_label <- stage_title_map[[stg]]  # Lookup readable label
  
  # Dynamically set x-axis label with white background
  x_axis_label <- if (grepl("age", stg)) {
    "Age"
  } else if (grepl("time", stg)) {
    "Time"
  } else {
    "Time"
  }
  
  
  
  # Compute label placement
  min_x <- min(df_stage$age, na.rm = TRUE)
  max_x <- max(df_stage$age, na.rm = TRUE)
  max_y <- max(df_stage$n_ade, na.rm = TRUE) / 1000
  
  label_x <- if (grepl("age", stg)) min_x else max_x
  hjust_val <- if (grepl("age", stg)) -0.1 else 1.1
  
  ggplot(df_stage, aes(x = age, y = n_ade / 1000, color = type)) +
    #geom_ribbon(
    #  data = df_ribbon %>%
    #    pivot_wider(names_from = type, values_from = n_ade) %>%
    #    arrange(age),
    #  aes(x = age, ymin = `Post. LB CI (2.5%)` / 1000, ymax = `Post. UB CI (97.5%)` / 1000),
    #  fill = "#F9C74F",
    #  alpha = 0.3,
    #  inherit.aes = FALSE
    #) +
    geom_line(aes(group = type), size = 2) +
    annotate(
      "text",
      x = label_x, y = max_y,
      label = stage_label,
      hjust = hjust_val, vjust = 1.5,
      fontface = "bold",
      size = 7,
      family = "Arial",
      color = "#455A64"  
    ) +
    labs(x = x_axis_label) +
    scale_x_continuous(breaks = seq(0, 100, 10)) +
    scale_y_continuous(
      breaks = scales::breaks_pretty(n = 5),
      labels = scales::label_comma()
    ) +
    theme_bw(base_size = 24) +
    scale_color_manual(values = palette, labels = legend_labels) +
    theme(
      axis.title.x = ggtext::element_textbox(
        size = 20,
        face = "bold",
        family = "Arial",
        padding = margin(2, 6, 2, 6),  # Slight padding for box shape
        margin = margin(t = -17, l=450),    # Move upward
        fill = "white",               # Box background
        box.color = NA,               # No visible border (set to "red" if debugging)
        halign = 0.5,                 # Center the box (not full-width)
        width = NULL                  # 🔥 Shrinks box to fit text
      ),
      axis.title.y = element_blank(),
      legend.position = "none",
      plot.title = element_blank(),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 18),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 10, r = 10, b = 40, l = 10)
    )
  
  
})





for (i in 13:15) {
  plot_list[[i]] <- plot_list[[i]] +
    scale_y_continuous(limits = c(0, 60), expand = expansion(mult = c(0, 0.05))) +
    theme(
      axis.title.x = ggtext::element_textbox(
        size = 20,
        face = "bold",
        family = "Arial",
        padding = margin(2, 6, 2, 6),
        margin = margin(t = -80, l = 450),
        fill = "white",
        box.color = NA,
        halign = 0.5,
        width = NULL
      ),
      plot.margin = margin(0, 0, 0, 0)
    )
}

# Extract all unique versions and stages
all_versions <- unique(gsub("^R_", "", df_r_long$type))
all_stages <- unique(df_r_long$stage)

# Function to extract stats per stage and version
extract_stats <- function(stage_name, version) {
  df_r <- df_r_long[stage == stage_name & type == paste0("R_", version)]
  
  values_r <- df_r[!is.na(age_raw), age_raw]
  
  mean_r <- mean(values_r)
  n_r <- length(values_r)
  
  list(
    mean_r = round(mean_r, 2),
    n_r = n_r
  )
}



# Initialize result table
df_stats <- data.table(stage = all_stages)

# Add R version results
for (v in all_versions) {
  stats_v <- lapply(all_stages, function(s) extract_stats(s, v))
  df_stats[[paste0("Mean (R_", v, ")")]] <- sapply(stats_v, function(x) x$mean_r)
  df_stats[[paste0("N (R_", v, ")")]] <- sapply(stats_v, function(x) x$n_r)
}


# Set stage as rownames if desired
setDF(df_stats)

#drop rownames

rownames(df_stats) <- NULL

df_stats_raw <- copy(df_stats)

#reorder stage

df_stats_raw$stage <- factor(df_stats_raw$stage,
                             levels = c(
                               "age_to_onset",
                               "time_to_6mm",
                               "age_to_6mm",
                               "time_to_10mm",
                               "age_to_10mm",
                               "time_from_6mm_to_preclin",
                               "time_from_10mm_to_preclin",
                               "age_preclin",
                               "age_clin",
                               "age_death_crc",
                               "dwell_time",
                               "sojourn_time"
                             )
)

df_stats_raw <- df_stats_raw[order(df_stats_raw$stage), ]


df_stats_raw$stage <- dplyr::recode(df_stats_raw$stage,
                                    "age_to_onset"   = "Age to Lesion Onset",
                                    "time_to_6mm"    = "Time to 6mm",
                                    "age_to_6mm"     = "Age to 6mm",
                                    "time_to_10mm"   = "Time to 10mm",
                                    "age_to_10mm"    = "Age to 10mm",
                                    "time_from_6mm_to_preclin" = "Time from 6mm to Preclinical",
                                    "time_from_10mm_to_preclin" = "Time from 10mm to Preclinical",
                                    "age_preclin"    = "Age to Preclinical CRC",
                                    "age_clin"       = "Age to Clinical CRC",
                                    "age_death_crc"  = "Age to CRC Death",
                                    "dwell_time"     = "Dwell Time",
                                    "sojourn_time"   = "Sojourn Time"
)




# Create two rows: one for mean, one for N
df_means <- df_stats_raw %>%
  transmute(
    Stage = stage,
    `R v0.15.0-C MAE (Chile)` = round(`Mean (R_R v0.15.0-C MAE (Chile))`, 1),
    `R v0.15.0-C MAE (USA)` = round(`Mean (R_R v0.15.0-C MAE (USA))`, 1))

# Step 1: Pivot to long format
df_means_long <- df_means %>%
  pivot_longer(
    cols = -Stage,
    names_to = "Model",
    values_to = "Mean"
  ) %>%
  mutate(Mean = as.numeric(Mean))

df_means_long$Model <- factor(df_means_long$Model, levels = names(palette))

#rename Time from 6mm to Preclinical to Time from 6mm\nto Preclinical

df_means_long$Stage <- dplyr::recode(df_means_long$Stage,
                                     "Time from 6mm to Preclinical" = "Time from 6mm\nto Preclinical",
                                     "Time from 10mm to Preclinical" = "Time from 10mm\nto Preclinical",
                                     "Time to 10mm" = "Ade(s) time from\n6mm to 10mm",
)

# Combine the two plots

unique(df_means_long$Model)

plot_mean <- ggplot(
  df_means_long %>%
    drop_na() %>%
    filter(Model %in% c("R v0.15.0-C MAE (Chile)", "R v0.15.0-C MAE (USA)")),
  aes(x = Model, y = 1, label = sprintf("%.1f", Mean), color = Model)
) +
  geom_text(
    size = 8,
    fontface = "bold"
  ) +
  facet_wrap(~Stage, ncol = 1) +
  scale_color_manual(values = palette) +
  labs(
    title = "Mean",
    x = "",
    y = ""
  ) +
  scale_y_continuous(limits = c(0.9, 1.1)) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20,
                              margin = margin(b = 8)),
    strip.text = element_text(face = "bold", size = 17),
    plot.margin = margin(t = 0, r = -100, b = 0, l = 0),
    legend.position = "none"
  )


# Step 1: Create clean df_ns
df_ns_clean <- df_stats_raw %>%
  #drop Time to 6mm, Time to 10mm
  filter(!grepl("Time to", stage)) %>%
  filter(!grepl("Time from", stage)) %>%
  #filter Dwell Time and Sojourn Time
  filter(!grepl("Dwell Time|Sojourn Time", stage)) %>%
  transmute(
    Stage = stage,
    `R v0.15.0-C MAE (Chile)` = as.numeric(gsub("[(),]", "", `N (R_R v0.15.0-C MAE (Chile))`)),
    `R v0.15.0-C MAE (USA)` = as.numeric(gsub("[(),]", "", `N (R_R v0.15.0-C MAE (USA))`)))
      
# Step 2: Pivot to long format
df_ns_long <- df_ns_clean %>%
  pivot_longer(
    cols = -Stage,
    names_to = "Model",
    values_to = "N"
  )


df_ns_long <- df_ns_long %>% drop_na()

df_ns_long$Model <- factor(df_ns_long$Model, levels = c(
  
  "R v0.15.0-C MAE (Chile)",
  "R v0.15.0-C MAE (USA)"))


# Generate colored labels for Models
model_colors <- palette
model_labels_colored <- sapply(levels(df_ns_long$Model), function(m) {
  paste0("<span style='color:", model_colors[[m]], ";'><b>", m, "</b></span>")
})
names(model_labels_colored) <- levels(df_ns_long$Model)




#put at the beginning of the data frame



df_ns_long$Stage <- dplyr::case_when(
  df_ns_long$Stage == "Age to Lesion Onset"     ~ "Total Number\nAdenomas",
  df_ns_long$Stage == "Age to 6mm"              ~ "Ade(s) reaching\n6mm",
  df_ns_long$Stage == "Age to 10mm"             ~ "Ade(s) reaching\n10mm",
  df_ns_long$Stage == "Age to Preclinical CRC"  ~ "Ade(s) reaching\nPreclinical CRC",
  df_ns_long$Stage == "Age to Clinical CRC"     ~ "Ade(s) reaching\nClinical CRC",
  df_ns_long$Stage == "Age to CRC Death"        ~ "Ade(s) reaching\nCRC Death",
  TRUE                                          ~ df_ns_long$Stage
)


unique(df_ns_long$Stage)

df_ns_long$Stage <- factor(df_ns_long$Stage, levels = c(
  "Total Number\nAdenomas",
  "Ade(s) reaching\n6mm",
  "Ade(s) reaching\n10mm",
  "Ade(s) reaching\nPreclinical CRC",
  "Ade(s) reaching\nClinical CRC",
  "Ade(s) reaching\nCRC Death"
))

df_plot_n <- df_ns_long %>%
  drop_na() %>%
  filter(Model %in% c("R v0.15.0-C MAE (Chile)", "R v0.15.0-C MAE (USA)")) %>%
  mutate(
    is_adenomas = Stage == "Total Number\nAdenomas",
    hjust_label = ifelse(is_adenomas, 1.1, NA_real_)  # only for label inside line
  )

# Order the two model labels consistently.

df_plot_n$Model <- factor(df_plot_n$Model, levels = c(
  
  "R v0.15.0-C MAE (Chile)",
  "R v0.15.0-C MAE (USA)"
))

# Now the plot:
plot_n <- ggplot(df_plot_n, aes(y = Model, x = N, color = Model)) +
  geom_segment(aes(x = 0, xend = N, y = Model, yend = Model), size = 1.5) +
  geom_point(size = 5) +
  
  # 🔹 Move label "inside" only for Number of Adenomas
  geom_label(
    data = df_plot_n %>% filter(is_adenomas),
    aes(label = scales::comma(N), hjust = hjust_label),
    size = 7,
    fontface = "bold",
    fill = "white",
    label.size = NA
  ) +
  
  # 🔹 Keep normal geom_text for all other stages
  geom_text(
    data = df_plot_n %>% filter(!is_adenomas),
    aes(label = scales::comma(N)),
    hjust = -0.2,
    size = 7,
    fontface = "bold"
  ) +
  
  facet_wrap(~Stage, ncol = 1, scales = "free_y") +
  scale_color_manual(values = palette) +
  scale_y_discrete(labels = model_labels_colored) +
  labs(
    title = "Ade(s) progress\nby model",
    x = "",
    y = ""
  ) +
  scale_x_continuous(
    breaks = scales::breaks_pretty(n = 3),
    labels = label_number(scale_cut = cut_short_scale()),
    expand = expansion(mult = c(0, 0.2))
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.spacing = unit(1, "lines"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20,
                              margin = margin(b = 8)),
    strip.text = element_text(face = "bold", size = 16),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_blank(),
    legend.position = "none"
  ) +
  coord_cartesian(clip = "off")



df_dummy_v_020_070 <- data.frame(
  type = factor(rep(names(palette), each = 2),
                levels = c("R v0.15.0-C MAE (Chile)", "R v0.15.0-C MAE (USA)")),
  x = rep(c(1, 2), times = length(palette)),
  y = 1
)

# Legend-only plot
legend_only_plot_v_020_070 <- ggplot(df_dummy_v_020_070 %>% 
                                       filter(type %in% c("R v0.15.0-C MAE (Chile)", "R v0.15.0-C MAE (USA)")
                                       ),
                                     aes(x = x, y = y, color = type, group = type)) +
  geom_line(linewidth = 0) +  # line actually visible now
  scale_color_manual(values = palette, labels = legend_labels) +
  guides(color = guide_legend(
    title = "", nrow = 1, byrow = TRUE,
    override.aes = list(linewidth = 4)
  )) +
  theme_void() +
  theme(
    plot.margin = margin(1, 1, 1, 1, "pt"),  # ultra-thin margins
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = c(0.5, 0.5),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.box.background = element_blank(),
    legend.box = "horizontal",
    legend.key.size = unit(2, "line"),
    legend.text = element_markdown(size = 25, face = "bold", family = "Arial"),
    legend.title = element_blank()
  )


# 2. Final ggplot

df_plot <- df_plot_data %>%
  filter(stage == "n_adenomas_total") %>%
  mutate(
    vjust_pos = ifelse(
      type %in% c("R v0.15.0-C MAE (Chile)", "R v0.15.0-C MAE (USA)") & age >= 3,  # 🔥 combined condition
      .5,  # move label above
      0.5    # otherwise beside
    ),
    hjust_pos = ifelse(
      type %in% c("R v0.15.0-C MAE (Chile)", "R v0.15.0-C MAE (USA)") & age >= 3,  # 🔥 combined condition
      1.2,  # centered
      1.2   # otherwise right aligned
    )
  )

#number of adenomas by model. group it

df_plot %>%
  group_by(type) %>%
  summarise(n_ade = sum(n_ade)) %>%
  ungroup() %>%
  mutate(
    type = factor(type, levels = c("R v0.15.0-C MAE (Chile)", "R v0.15.0-C MAE (USA)"))
  )

# Now the plot
plot_n_adenomas <- ggplot(df_plot %>% drop_na() %>% filter(type %in% c("R v0.15.0-C MAE (Chile)", "R v0.15.0-C MAE (USA)")
), aes(x = type, y = n_ade / 1000, color = type)) +
  geom_segment(aes(xend = type, y = 0, yend = n_ade / 1000), size = 1.5) +
  geom_point(size = 5) +
  geom_label(
    aes(
      label = scales::comma(n_ade),
      vjust = vjust_pos,    # 🔥 dynamic vjust
      hjust = hjust_pos     # 🔥 dynamic hjust
    ),
    angle = 90,
    size = 7,
    fontface = "bold",
    fill = "white",
    label.size = NA
  ) +
  facet_wrap(~age, ncol = 2, scales = "free_y") +
  scale_color_manual(values = palette) +
  scale_y_continuous(
    labels = scales::label_comma(),
    expand = expansion(mult = c(.20, .25))
  ) +
  labs(
    title = "Number of Adenomas",
    subtitle = "",
    x = "",
    y = "Lesions (per thousands)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title.y = element_text(size = 20),
    plot.title = element_text(vjust=-1, hjust = 0.5, face = "bold", size = 20),
    strip.text = element_text(face = "bold", size = 20),
    plot.caption = element_text(size = 10, hjust = 0),
    axis.text.x = ggtext::element_markdown(angle = 90, hjust = 1, size = 15),
    axis.text.y = element_text(size = 15),
    legend.position = "none"
  )

# Combine the plots
main_plot_mean_n <- wrap_plots(
  plot_mean,
  plot_n,
  ncol = 2
) + plot_layout(heights = c(1, 0.3), 
                widths = c(.7, 0.5))



main_plot_mean_n <- wrap_plots(
  plot_mean,
  plot_n,
  #ggleft_priors, 
  #ggright_priors,
  ncol = 2
) + plot_layout(heights = c(.8,.7), 
                widths = c(.5, .5))& theme(plot.margin = margin(0, 0, 0, 0))

# The original dwell panels already calculate these durations within the
# diagnosing lesion. Their means use diagnosed cases, unlike the general
# growth panels, which include all lesions reaching each size threshold.
# Small: age_to_6mm - age_lesion_onset
# Medium: age_to_10mm - age_to_6mm
# Large: age_crc - age_to_10mm
# Sojourn: age_dx_symp - age_crc
phase_rows <- match(c("dw_t_sm", "dw_t_md", "dw_t_lg", "sojourn_time"), df_stats$stage)
df_phase_means <- data.frame(
  Phase = c("Small adenoma (<6 mm)", "Medium adenoma (6-10 mm)",
            "Large adenoma before cancer", "Cancer before diagnosis (sojourn)"),
  Chile = df_stats[["Mean (R_R v0.15.0-C MAE (Chile))"]][phase_rows],
  USA = df_stats[["Mean (R_R v0.15.0-C MAE (USA))"]][phase_rows]
)
phase_counts <- c(
  Chile = df_stats[["N (R_R v0.15.0-C MAE (Chile))"]][phase_rows[1]],
  USA = df_stats[["N (R_R v0.15.0-C MAE (USA))"]][phase_rows[1]]
)
write.csv(df_phase_means, file.path(path, "phase_duration_means.csv"), row.names = FALSE)

# Use the existing Mean column's labels, country positions and typography.
df_phase_long <- df_phase_means %>%
  pivot_longer(c(Chile, USA), names_to = "Country", values_to = "Mean") %>%
  mutate(
    Stage = factor(Phase, levels = df_phase_means$Phase,
                   labels = c("Small adenoma (<6 mm)", "Medium adenoma (6-10 mm)",
                              "Large adenoma\nbefore cancer",
                              "Cancer before diagnosis\n(sojourn)")),
    Model = factor(paste0("R v0.15.0-C MAE (", Country, ")"), levels = names(palette))
  )
# Append phase rows to the same plot so both countries stay exactly aligned.
mean_rows <- plot_mean$data %>%
  mutate(Stage = as.character(Stage), Display = sprintf("%.1f", Mean))
phase_mean_rows <- df_phase_long %>%
  mutate(Stage = as.character(Stage), Display = sprintf("%.2f", Mean))
# Give the heading its own row; keep ordinary facet labels as plain text.
phase_heading <- data.frame(
  Stage = "phase_heading", Model = factor(names(palette)[1], levels = names(palette)),
  Mean = NA_real_, Display = ""
)
mean_order <- unique(c(mean_rows$Stage, "phase_heading", phase_mean_rows$Stage))
df_all_means <- bind_rows(mean_rows, phase_heading, phase_mean_rows) %>%
  mutate(Stage = factor(Stage, levels = mean_order))
phase_heading$Stage <- factor(phase_heading$Stage, levels = mean_order)
mean_column <- plot_mean %+% df_all_means +
  aes(label = Display) +
  facet_wrap(~Stage, ncol = 1,
             labeller = labeller(Stage = function(labels) {
               ifelse(labels == "phase_heading", "", labels)
             })) +
  geom_text(data = phase_heading, aes(x = 1.5, y = 1),
            label = "Phases dwell time", size = 9, fontface = "bold",
            inherit.aes = FALSE) +
  labs(caption = sprintf("Phase cases: Chile %s | USA %s",
                         scales::comma(phase_counts[["Chile"]]),
                         scales::comma(phase_counts[["USA"]]))) +
  theme(plot.caption = element_text(size = 12, hjust = .5),
        strip.text = element_text(size = 17, face = "bold"),
        plot.margin = margin(0, 0, 0, 0))

# Keep the event-count panel at its original height beside the extended means.
count_column <- (plot_n / plot_spacer()) + plot_layout(heights = c(.8, .7))
main_plot_mean_n <- wrap_plots(mean_column, count_column, ncol = 2) +
  plot_layout(widths = c(.5, .5)) & theme(plot.margin = margin(0, 0, 0, 0))


# Combine all plots into a grid



main_plot_grid <- wrap_plots(
  plot_list,
  ncol = 3
) +
  plot_layout(
    guides = "collect",
    widths = c(1, 1, 1),
    heights = rep(3, 5)  # 🔥 Increase row height (was 1, now 2)
  ) & 
  theme(
    plot.margin = margin(3, 3, 3, 3),
    panel.spacing = unit(0.2, "lines")
  )



df_w <- wrap_plots(plot_n_adenomas,
                   main_plot_grid, main_plot_mean_n,
                   ncol = 3
) + plot_layout(widths = c(.2,1.7,.7),
                heights = c(30, -1.1,-1)) & 
  theme(
    plot.margin = margin(3, 3, 3, 3),
    panel.spacing = unit(0.2, "lines")
  )



sheet <- (df_w / legend_only_plot_v_020_070) + plot_layout(heights = c(200, 1, 1.5)) +
  plot_annotation(
    title = "SimCRC Chile and USA v0.15.0",
    subtitle = paste0("Population: 1 million people","\n","Generated on: ", format(Sys.Date(), "%d-%b-%Y")),
    theme = theme(
      plot.title = element_text(
        size = 26,
        face = "bold",
        hjust = 0,       # left
        vjust = 1
      ),
      plot.subtitle = element_text(
        size = 20,
        face = "bold",
        hjust = 0,
        vjust = 1
      ),
      plot.caption = element_text(
        size = 22,
        face = "plain",
        hjust = 0,       # left-align caption
        vjust = 0,       # stay at bottom
        margin = margin(t = 10)  # tiny spacing above caption
      )
    )
  )

ggsave(paste0(path, "cheat_sheet_Chile_USA_v0.15.0.png"),
       plot = sheet, device = ragg::agg_png,
       width = 36, height = 24, dpi = 200)
write.csv(df_stats, file.path(path, "natural_history_summary.csv"), row.names = FALSE)
message("Saved Chile-USA cheat sheet to ", path)
