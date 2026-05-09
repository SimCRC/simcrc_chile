
# Compare priors: Chile vs USA  -----------------------------------------------
# Plots the prior range (lower bound – upper bound) for each calibrated
# parameter, side by side for Chile and USA.

library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Load priors ---------------------------------------------------------------

# Chile — latest v0.13.0 calibration
load("outputs/BayCANN_versions/Chile/Adenoma/F/v0.13.0/v0.13.0.20260406.1214/l_params_priors_SimCRC_SimCRC_v0.13.0.20260406.1214_Adenoma_F.RData")
priors_chile <- l_priors

# USA — v0.13.0
load("outputs/BayCANN_versions/USA/Adenoma/F/v0.13.0/v0.13.0.20260313.2334/l_params_priors_SimCRC_SimCRC_v0.13.0.20260313.2334_Adenoma_F.RData")
priors_usa <- l_priors

rm(l_priors)

# 2. Build tidy data frame -----------------------------------------------------

build_prior_df <- function(priors, country) {
  data.frame(
    parameter = priors$names,
    lb        = as.numeric(priors$lb),
    ub        = as.numeric(priors$ub),
    country   = country,
    stringsAsFactors = FALSE
  )
}

# Helper: compute derived symptom-detection rates from prior bounds
# All transformations are monotonically increasing on positive values,
# so lb→lb and ub→ub.
build_derived_rates <- function(priors, country) {
  # Look up lb/ub by name
  get_bound <- function(name, bound = "lb") {
    idx <- match(name, priors$names)
    if (is.na(idx)) return(NA_real_)
    as.numeric(if (bound == "lb") priors$lb[idx] else priors$ub[idx])
  }

  locations <- c("P", "D", "R")
  rows <- list()

  for (loc in locations) {
    # Stage 1: probability → rate
    p_lb <- get_bound(paste0("pSxDetS1_", loc), "lb")
    p_ub <- get_bound(paste0("pSxDetS1_", loc), "ub")
    r1_lb <- -log(1 - p_lb)
    r1_ub <- -log(1 - p_ub)
    # Stage 2: hr × rate_S1 → probability
    hr2_lb <- get_bound(paste0("hr_SxDetS2S1_", loc), "lb")
    hr2_ub <- get_bound(paste0("hr_SxDetS2S1_", loc), "ub")
    r2_lb <- hr2_lb * r1_lb
    r2_ub <- hr2_ub * r1_ub
    rows[[paste0("pSxDetS2_", loc)]] <- c(1 - exp(-r2_lb), 1 - exp(-r2_ub))

    # Stage 3: hr × rate_S2 → probability
    hr3_lb <- get_bound(paste0("hr_SxDetS3S2_", loc), "lb")
    hr3_ub <- get_bound(paste0("hr_SxDetS3S2_", loc), "ub")
    r3_lb <- hr3_lb * r2_lb
    r3_ub <- hr3_ub * r2_ub
    rows[[paste0("pSxDetS3_", loc)]] <- c(1 - exp(-r3_lb), 1 - exp(-r3_ub))

    # Stage 4: hr × rate_S3 → probability
    hr4_lb <- get_bound(paste0("hr_SxDetS4S3_", loc), "lb")
    hr4_ub <- get_bound(paste0("hr_SxDetS4S3_", loc), "ub")
    r4_lb <- hr4_lb * r3_lb
    r4_ub <- hr4_ub * r3_ub
    rows[[paste0("pSxDetS4_", loc)]] <- c(1 - exp(-r4_lb), 1 - exp(-r4_ub))
  }

  data.frame(
    parameter = names(rows),
    lb        = sapply(rows, `[`, 1),
    ub        = sapply(rows, `[`, 2),
    country   = country,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

df_priors <- rbind(
  build_prior_df(priors_chile, "Chile"),
  build_prior_df(priors_usa,   "USA"),
  build_derived_rates(priors_chile, "Chile"),
  build_derived_rates(priors_usa,   "USA")
)

# Midpoint for sorting / reference
df_priors$mid <- (df_priors$lb + df_priors$ub) / 2

# 3. Assign parameter groups for faceting --------------------------------------

df_priors$group <- case_when(
  grepl("IndividuakRiskMult|alpha_lesion|beta_age|Gaussian_Variance", df_priors$parameter) ~ "Risk & Onset",
  grepl("AdGrowth",           df_priors$parameter) ~ "Adenoma Growth",
  grepl("PreclinCancerProg",  df_priors$parameter) ~ "Cancer Progression",
  grepl("CancerOnset",        df_priors$parameter) ~ "Cancer Onset (Gompertz)",
  grepl("pSxDetS1",            df_priors$parameter) ~ "Sx Detection Prob (Stage 1)",
  grepl("pSxDetS2",            df_priors$parameter) ~ "Sx Detection Prob (Stage 2)",
  grepl("pSxDetS3",            df_priors$parameter) ~ "Sx Detection Prob (Stage 3)",
  grepl("pSxDetS4",            df_priors$parameter) ~ "Sx Detection Prob (Stage 4)",
  grepl("hr_SxDet",            df_priors$parameter) ~ "Symptom Detection HR",
  TRUE                                              ~ "Other"
)

# Set facet group order
group_order <- c(
  "Risk & Onset", "Adenoma Growth",
  "Cancer Onset (Gompertz)", "Cancer Progression",
  "Sx Detection Prob (Stage 1)", "Symptom Detection HR",
  "Sx Detection Prob (Stage 2)", "Sx Detection Prob (Stage 3)",
  "Sx Detection Prob (Stage 4)",
  "Other"
)
df_priors$group <- factor(df_priors$group, levels = group_order)

# Define logical parameter order (grouped, then within-group order)
# Includes both prior names and posterior CSV names (which differ)
param_order <- c(
  # Risk & Onset
  "IndividuakRiskMult_Mean", "IndividuakRiskMult_SD",
  "alpha_lesion_adenoma",
  "beta_age_adenoma",
  "AdNaturalHistoryPropensity_Gaussian_Variance",
  # Adenoma Growth
  "AdGrowth_A", "AdGrowth_B", "AdGrowth_C",
  # Cancer Onset (Gompertz) — Rate then Shape, P/D/R within each
  "CancerOnset_Gompertz_Rate_P", "CancerOnset_Gompertz_Rate_D", "CancerOnset_Gompertz_Rate_R",
  "CancerOnset_Gompertz_Shape_P", "CancerOnset_Gompertz_Shape_D", "CancerOnset_Gompertz_Shape_R",
  # Cancer Progression — stage by stage, P/D/R within each
  "PreclinCancerProg_S1S2_P", "PreclinCancerProg_S1S2_D", "PreclinCancerProg_S1S2_R",
  "PreclinCancerProg_S2S3_P", "PreclinCancerProg_S2S3_D", "PreclinCancerProg_S2S3_R",
  "PreclinCancerProg_S3S4_P", "PreclinCancerProg_S3S4_D", "PreclinCancerProg_S3S4_R",
  # Symptom Detection (probabilities)
  "pSxDetS1_P", "pSxDetS1_D", "pSxDetS1_R",
  # Symptom Detection HR — stage by stage, P/D/R within each
  "hr_SxDetS2S1_P", "hr_SxDetS2S1_D", "hr_SxDetS2S1_R",
  "hr_SxDetS3S2_P", "hr_SxDetS3S2_D", "hr_SxDetS3S2_R",
  "hr_SxDetS4S3_P", "hr_SxDetS4S3_D", "hr_SxDetS4S3_R",
  # Symptom Detection Probabilities (derived) — stage by stage, P/D/R within each
  "pSxDetS2_P", "pSxDetS2_D", "pSxDetS2_R",
  "pSxDetS3_P", "pSxDetS3_D", "pSxDetS3_R",
  "pSxDetS4_P", "pSxDetS4_D", "pSxDetS4_R"
)

# Keep only parameters that exist in the data, preserve order
param_order_present <- param_order[param_order %in% unique(df_priors$parameter)]
# Add any remaining parameters not in the explicit list
param_order_present <- c(param_order_present,
                         setdiff(unique(df_priors$parameter), param_order_present))

df_priors$parameter <- factor(df_priors$parameter,
                              levels = rev(param_order_present))

# 4. Detect narrow ranges (labels would overlap on top → push to sides) --------

# For each facet group, compute the full x-axis span and flag narrow ranges
df_priors <- df_priors %>%
  mutate(range = ub - lb) %>%
  group_by(group) %>%
  mutate(axis_span = max(ub) - min(lb),
         narrow    = range < 0.15 * axis_span) %>%
  ungroup()

# For narrow ranges: lb label goes left (hjust=1), ub label goes right (hjust=0)
# For wide ranges: both labels go on top (vjust=-1)
df_priors$lb_hjust <- ifelse(df_priors$narrow, 1.1,  0.5)
df_priors$lb_vjust <- ifelse(df_priors$narrow, 0.5, -1.0)
df_priors$ub_hjust <- ifelse(df_priors$narrow, -0.1, 0.5)
df_priors$ub_vjust <- ifelse(df_priors$narrow, 0.5, -1.0)

# 5. Plot — horizontal segments per parameter ----------------------------------

p <- ggplot(df_priors, aes(y = parameter, color = country)) +
  geom_linerange(aes(xmin = lb, xmax = ub),
                 linewidth = 2, alpha = 0.6,
                 position = position_dodge(width = 0.6)) +
  geom_point(aes(x = lb), size = 2, shape = "|",
             position = position_dodge(width = 0.6)) +
  geom_point(aes(x = ub), size = 2, shape = "|",
             position = position_dodge(width = 0.6)) +
  geom_text(aes(x = lb, label = sprintf("%.4g", lb),
                hjust = lb_hjust, vjust = lb_vjust),
            size = 2.5,
            position = position_dodge(width = 0.6), show.legend = FALSE) +
  geom_text(aes(x = ub, label = sprintf("%.4g", ub),
                hjust = ub_hjust, vjust = ub_vjust),
            size = 2.5,
            position = position_dodge(width = 0.6), show.legend = FALSE) +
  facet_wrap(~ group, scales = "free", ncol = 4) +
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  scale_color_manual(values = c("Chile" = "#E41A1C", "USA" = "#377EB8")) +
  theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(size = 16, face = "bold"),
    axis.text.y   = element_text(size = 9),
    strip.background = element_blank(),
    strip.text    = element_text(face = "bold", size = 11),
    legend.position = "top",
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.minor   = element_blank()
  ) +
  labs(
    title = "Prior Ranges: Chile vs USA (Adenoma F, v0.13.0)",
    x     = "Parameter value",
    y     = "",
    color = "Country"
  )

p

# 6. Save ---------------------------------------------------------------------

output_folder <- "outputs/BayCANN_versions/Chile/Adenoma/F/v0.13.0/v0.13.0.20260406.1214"

ggsave(p,
       filename = file.path(output_folder, "fig_priors_Chile_vs_USA_v0.13.0_Adenoma_F.png"),
       width = 20, height = 12, dpi = 300)


# 7. Load posteriors -----------------------------------------------------------

post_chile <- read.csv(file.path(output_folder,
                                 "dt_calibrated_posteriors_SimCRC_v0.13.0.20260406.1214_Adenoma_F.csv"))
post_usa   <- read.csv("outputs/BayCANN_versions/USA/Adenoma/F/v0.13.0/v0.13.0.20260313.2334/dt_calibrated_posteriors_SimCRC_v0.13.0.20260313.2334_Adenoma_F.csv")

# Drop non-parameter columns
drop_cols <- c("id_draw", "lp__", "chain")
post_chile <- post_chile[, !colnames(post_chile) %in% drop_cols]
post_usa   <- post_usa[,   !colnames(post_usa)   %in% drop_cols]

# 8. Compute derived symptom-detection rates from posterior draws --------------

# Chile — compute rates then convert S2-S4 to probabilities
rSxDetS1_P_chile <- -log(1 - post_chile$pSxDetS1_P)
rSxDetS1_D_chile <- -log(1 - post_chile$pSxDetS1_D)
rSxDetS1_R_chile <- -log(1 - post_chile$pSxDetS1_R)

rSxDetS2_P_chile <- post_chile$hr_SxDetS2S1_P * rSxDetS1_P_chile
rSxDetS2_D_chile <- post_chile$hr_SxDetS2S1_D * rSxDetS1_D_chile
rSxDetS2_R_chile <- post_chile$hr_SxDetS2S1_R * rSxDetS1_R_chile
post_chile$pSxDetS2_P <- 1 - exp(-rSxDetS2_P_chile)
post_chile$pSxDetS2_D <- 1 - exp(-rSxDetS2_D_chile)
post_chile$pSxDetS2_R <- 1 - exp(-rSxDetS2_R_chile)

rSxDetS3_P_chile <- post_chile$hr_SxDetS3S2_P * rSxDetS2_P_chile
rSxDetS3_D_chile <- post_chile$hr_SxDetS3S2_D * rSxDetS2_D_chile
rSxDetS3_R_chile <- post_chile$hr_SxDetS3S2_R * rSxDetS2_R_chile
post_chile$pSxDetS3_P <- 1 - exp(-rSxDetS3_P_chile)
post_chile$pSxDetS3_D <- 1 - exp(-rSxDetS3_D_chile)
post_chile$pSxDetS3_R <- 1 - exp(-rSxDetS3_R_chile)

rSxDetS4_P_chile <- post_chile$hr_SxDetS4S3_P * rSxDetS3_P_chile
rSxDetS4_D_chile <- post_chile$hr_SxDetS4S3_D * rSxDetS3_D_chile
rSxDetS4_R_chile <- post_chile$hr_SxDetS4S3_R * rSxDetS3_R_chile
post_chile$pSxDetS4_P <- 1 - exp(-rSxDetS4_P_chile)
post_chile$pSxDetS4_D <- 1 - exp(-rSxDetS4_D_chile)
post_chile$pSxDetS4_R <- 1 - exp(-rSxDetS4_R_chile)

# USA — compute rates then convert S2-S4 to probabilities
rSxDetS1_P_usa <- -log(1 - post_usa$pSxDetS1_P)
rSxDetS1_D_usa <- -log(1 - post_usa$pSxDetS1_D)
rSxDetS1_R_usa <- -log(1 - post_usa$pSxDetS1_R)

rSxDetS2_P_usa <- post_usa$hr_SxDetS2S1_P * rSxDetS1_P_usa
rSxDetS2_D_usa <- post_usa$hr_SxDetS2S1_D * rSxDetS1_D_usa
rSxDetS2_R_usa <- post_usa$hr_SxDetS2S1_R * rSxDetS1_R_usa
post_usa$pSxDetS2_P <- 1 - exp(-rSxDetS2_P_usa)
post_usa$pSxDetS2_D <- 1 - exp(-rSxDetS2_D_usa)
post_usa$pSxDetS2_R <- 1 - exp(-rSxDetS2_R_usa)

rSxDetS3_P_usa <- post_usa$hr_SxDetS3S2_P * rSxDetS2_P_usa
rSxDetS3_D_usa <- post_usa$hr_SxDetS3S2_D * rSxDetS2_D_usa
rSxDetS3_R_usa <- post_usa$hr_SxDetS3S2_R * rSxDetS2_R_usa
post_usa$pSxDetS3_P <- 1 - exp(-rSxDetS3_P_usa)
post_usa$pSxDetS3_D <- 1 - exp(-rSxDetS3_D_usa)
post_usa$pSxDetS3_R <- 1 - exp(-rSxDetS3_R_usa)

rSxDetS4_P_usa <- post_usa$hr_SxDetS4S3_P * rSxDetS3_P_usa
rSxDetS4_D_usa <- post_usa$hr_SxDetS4S3_D * rSxDetS3_D_usa
rSxDetS4_R_usa <- post_usa$hr_SxDetS4S3_R * rSxDetS3_R_usa
post_usa$pSxDetS4_P <- 1 - exp(-rSxDetS4_P_usa)
post_usa$pSxDetS4_D <- 1 - exp(-rSxDetS4_D_usa)
post_usa$pSxDetS4_R <- 1 - exp(-rSxDetS4_R_usa)

# 8b. Build tidy posterior data frame ------------------------------------------

df_post_chile <- post_chile %>%
  pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
  mutate(country = "Chile")

df_post_usa <- post_usa %>%
  pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
  mutate(country = "USA")

df_post <- bind_rows(df_post_chile, df_post_usa)

# Assign same parameter groups as priors (same order as case_when in section 3)
df_post$group <- case_when(
  grepl("IndividuakRiskMult|alpha_lesion|beta_age|Gaussian_Variance", df_post$parameter) ~ "Risk & Onset",
  grepl("AdGrowth",           df_post$parameter) ~ "Adenoma Growth",
  grepl("CancerOnset",        df_post$parameter) ~ "Cancer Onset (Gompertz)",
  grepl("PreclinCancerProg",  df_post$parameter) ~ "Cancer Progression",
  grepl("pSxDetS1",            df_post$parameter) ~ "Sx Detection Prob (Stage 1)",
  grepl("pSxDetS2",            df_post$parameter) ~ "Sx Detection Prob (Stage 2)",
  grepl("pSxDetS3",            df_post$parameter) ~ "Sx Detection Prob (Stage 3)",
  grepl("pSxDetS4",            df_post$parameter) ~ "Sx Detection Prob (Stage 4)",
  grepl("hr_SxDet",            df_post$parameter) ~ "Symptom Detection HR",
  TRUE                                            ~ "Other"
)

df_post$group <- factor(df_post$group, levels = group_order)

# 9. Order posteriors logically (posterior CSV uses different names) ------------

post_param_order <- c(
  # Risk & Onset
  "alpha_lesion_adenoma", "beta_age",
  "IndividuakRiskMultVariance",
  "AdNaturalHistoryPropensity_Gaussian_Variance",
  # Adenoma Growth
  "AdGrowth_Exp_Rate_1_to_6_P", "AdGrowth_Exp_Rate_1_to_6_D", "AdGrowth_Exp_Rate_1_to_6_R",
  "AdGrowth_Exp_Rate_6_to_10_P", "AdGrowth_Exp_Rate_6_to_10_D", "AdGrowth_Exp_Rate_6_to_10_R",
  # Cancer Onset (Gompertz)
  "CancerOnset_Gompertz_Rate_P", "CancerOnset_Gompertz_Rate_D", "CancerOnset_Gompertz_Rate_R",
  "CancerOnset_Gompertz_Shape_P", "CancerOnset_Gompertz_Shape_D", "CancerOnset_Gompertz_Shape_R",
  # Cancer Progression
  "PreclinCancerProg_Exp_Rate_S1S2_P", "PreclinCancerProg_Exp_Rate_S1S2_D", "PreclinCancerProg_Exp_Rate_S1S2_R",
  "PreclinCancerProg_Exp_Rate_S2S3_P", "PreclinCancerProg_Exp_Rate_S2S3_D", "PreclinCancerProg_Exp_Rate_S2S3_R",
  "PreclinCancerProg_Exp_Rate_S3S4_P", "PreclinCancerProg_Exp_Rate_S3S4_D", "PreclinCancerProg_Exp_Rate_S3S4_R",
  # Symptom Detection
  "pSxDetS1_P", "pSxDetS1_D", "pSxDetS1_R",
  # Symptom Detection HR
  "hr_SxDetS2S1_P", "hr_SxDetS2S1_D", "hr_SxDetS2S1_R",
  "hr_SxDetS3S2_P", "hr_SxDetS3S2_D", "hr_SxDetS3S2_R",
  "hr_SxDetS4S3_P", "hr_SxDetS4S3_D", "hr_SxDetS4S3_R",
  # Symptom Detection Probabilities (derived)
  "pSxDetS2_P", "pSxDetS2_D", "pSxDetS2_R",
  "pSxDetS3_P", "pSxDetS3_D", "pSxDetS3_R",
  "pSxDetS4_P", "pSxDetS4_D", "pSxDetS4_R"
)

# Keep only parameters present, add any extras at the end
post_param_order <- post_param_order[post_param_order %in% unique(df_post$parameter)]
post_param_order <- c(post_param_order,
                      setdiff(unique(df_post$parameter), post_param_order))
df_post$parameter <- factor(df_post$parameter, levels = post_param_order)

# Split Stage 4 probabilities by country (magnitude too different to share axis)
s4_params <- c("pSxDetS4_P", "pSxDetS4_D", "pSxDetS4_R")
df_post_s4 <- df_post[df_post$parameter %in% s4_params, ]
df_post    <- df_post[!df_post$parameter %in% s4_params, ]

df_post_s4$parameter <- paste0(df_post_s4$parameter, " (", df_post_s4$country, ")")

# Update factor levels: replace pSxDetS4_X with country-split versions
s4_split_levels <- c(
  "pSxDetS4_P (Chile)", "pSxDetS4_P (USA)",
  "pSxDetS4_D (Chile)", "pSxDetS4_D (USA)",
  "pSxDetS4_R (Chile)", "pSxDetS4_R (USA)"
)
new_levels <- c(setdiff(levels(df_post$parameter), s4_params), s4_split_levels)
df_post <- bind_rows(df_post, df_post_s4)
df_post$parameter <- factor(df_post$parameter, levels = new_levels)

# 9b. Plot posteriors — density curves per parameter ---------------------------

p2 <- ggplot(df_post, aes(x = value, color = country, fill = country)) +
  geom_density(alpha = 0.2, linewidth = 0.8) +
  facet_wrap(~ parameter, scales = "free", ncol = 8) +
  scale_color_manual(values = c("Chile" = "#E41A1C", "USA" = "#377EB8")) +
  scale_fill_manual(values  = c("Chile" = "#E41A1C", "USA" = "#377EB8")) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(size = 16, face = "bold"),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 8),
    legend.position  = "top",
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Posterior Distributions: Chile vs USA (Adenoma F, v0.13.0)",
    x     = "Parameter value",
    y     = "Density",
    color = "Country",
    fill  = "Country"
  )

p2

# 10. Save posteriors ----------------------------------------------------------

ggsave(p2,
       filename = file.path(output_folder, "fig_posteriors_Chile_vs_USA_v0.13.0_Adenoma_F.png"),
       width = 28, height = 18, dpi = 300)

cat("Saved posteriors to:", file.path(output_folder, "fig_posteriors_Chile_vs_USA_v0.13.0_Adenoma_F.png"), "\n")

