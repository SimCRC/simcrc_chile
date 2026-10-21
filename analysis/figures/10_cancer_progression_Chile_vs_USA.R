
# Cancer progression parameters: Chile vs USA  --------------------------------
# Converts calibrated exponential rates to annual transition probabilities
# and summarizes detection probabilities for both countries.

library(dplyr)
library(tidyr)

# 1. Load posteriors -----------------------------------------------------------

output_folder <- "outputs/BayCANN_versions/Chile/Adenoma/F/v0.13.0/v0.13.0.20260406.1214"

post_chile <- read.csv(file.path(output_folder,
                                 "dt_calibrated_posteriors_SimCRC_v0.13.0.20260406.1214_Adenoma_F.csv"))
post_usa   <- read.csv("outputs/BayCANN_versions/USA/Adenoma/F/v0.13.0/v0.13.0.20260313.2334/dt_calibrated_posteriors_SimCRC_v0.13.0.20260313.2334_Adenoma_F.csv")

# 2. Compute annual transition probabilities from rates ------------------------

# Cancer progression: p = 1 - exp(-rate)
for (loc in c("P", "D", "R")) {
  for (trans in c("S1S2", "S2S3", "S3S4")) {
    rate_col <- paste0("PreclinCancerProg_Exp_Rate_", trans, "_", loc)
    prob_col <- paste0("pCancerProg_", trans, "_", loc)
    post_chile[[prob_col]] <- 1 - exp(-post_chile[[rate_col]])
    post_usa[[prob_col]]   <- 1 - exp(-post_usa[[rate_col]])
  }
}

# 3. Compute symptom detection probabilities (S2-S4 from HR × rate) ------------

# Chile
rS1_P_chile <- -log(1 - post_chile$pSxDetS1_P)
rS1_D_chile <- -log(1 - post_chile$pSxDetS1_D)
rS1_R_chile <- -log(1 - post_chile$pSxDetS1_R)

rS2_P_chile <- post_chile$hr_SxDetS2S1_P * rS1_P_chile
rS2_D_chile <- post_chile$hr_SxDetS2S1_D * rS1_D_chile
rS2_R_chile <- post_chile$hr_SxDetS2S1_R * rS1_R_chile

rS3_P_chile <- post_chile$hr_SxDetS3S2_P * rS2_P_chile
rS3_D_chile <- post_chile$hr_SxDetS3S2_D * rS2_D_chile
rS3_R_chile <- post_chile$hr_SxDetS3S2_R * rS2_R_chile

rS4_P_chile <- post_chile$hr_SxDetS4S3_P * rS3_P_chile
rS4_D_chile <- post_chile$hr_SxDetS4S3_D * rS3_D_chile
rS4_R_chile <- post_chile$hr_SxDetS4S3_R * rS3_R_chile

post_chile$pSxDetS2_P <- 1 - exp(-rS2_P_chile)
post_chile$pSxDetS2_D <- 1 - exp(-rS2_D_chile)
post_chile$pSxDetS2_R <- 1 - exp(-rS2_R_chile)
post_chile$pSxDetS3_P <- 1 - exp(-rS3_P_chile)
post_chile$pSxDetS3_D <- 1 - exp(-rS3_D_chile)
post_chile$pSxDetS3_R <- 1 - exp(-rS3_R_chile)
post_chile$pSxDetS4_P <- 1 - exp(-rS4_P_chile)
post_chile$pSxDetS4_D <- 1 - exp(-rS4_D_chile)
post_chile$pSxDetS4_R <- 1 - exp(-rS4_R_chile)

# USA
rS1_P_usa <- -log(1 - post_usa$pSxDetS1_P)
rS1_D_usa <- -log(1 - post_usa$pSxDetS1_D)
rS1_R_usa <- -log(1 - post_usa$pSxDetS1_R)

rS2_P_usa <- post_usa$hr_SxDetS2S1_P * rS1_P_usa
rS2_D_usa <- post_usa$hr_SxDetS2S1_D * rS1_D_usa
rS2_R_usa <- post_usa$hr_SxDetS2S1_R * rS1_R_usa

rS3_P_usa <- post_usa$hr_SxDetS3S2_P * rS2_P_usa
rS3_D_usa <- post_usa$hr_SxDetS3S2_D * rS2_D_usa
rS3_R_usa <- post_usa$hr_SxDetS3S2_R * rS2_R_usa

rS4_P_usa <- post_usa$hr_SxDetS4S3_P * rS3_P_usa
rS4_D_usa <- post_usa$hr_SxDetS4S3_D * rS3_D_usa
rS4_R_usa <- post_usa$hr_SxDetS4S3_R * rS3_R_usa

post_usa$pSxDetS2_P <- 1 - exp(-rS2_P_usa)
post_usa$pSxDetS2_D <- 1 - exp(-rS2_D_usa)
post_usa$pSxDetS2_R <- 1 - exp(-rS2_R_usa)
post_usa$pSxDetS3_P <- 1 - exp(-rS3_P_usa)
post_usa$pSxDetS3_D <- 1 - exp(-rS3_D_usa)
post_usa$pSxDetS3_R <- 1 - exp(-rS3_R_usa)
post_usa$pSxDetS4_P <- 1 - exp(-rS4_P_usa)
post_usa$pSxDetS4_D <- 1 - exp(-rS4_D_usa)
post_usa$pSxDetS4_R <- 1 - exp(-rS4_R_usa)

# 4. Build summary table -------------------------------------------------------

# Parameters to summarize
params_to_summarize <- c(
  # Cancer progression probabilities
  paste0("pCancerProg_", rep(c("S1S2", "S2S3", "S3S4"), each = 3), "_", c("P", "D", "R")),
  # Symptom detection probabilities
  paste0("pSxDetS1_", c("P", "D", "R")),
  paste0("pSxDetS2_", c("P", "D", "R")),
  paste0("pSxDetS3_", c("P", "D", "R")),
  paste0("pSxDetS4_", c("P", "D", "R"))
)

# Load Min Absolute Error parameter sets
load(file.path(output_folder,
               "l_params_calibrated_sets_SimCRC_v0.13.0.20260406.1214_Adenoma_F.RData"))
mae_chile <- l_params_calibrated_sets$Min_AbsolutErr

load("outputs/BayCANN_versions/USA/Adenoma/F/v0.13.0/v0.13.0.20260313.2334/l_params_calibrated_sets_SimCRC_v0.13.0.20260313.2334_Adenoma_F.RData")
mae_usa <- l_params_calibrated_sets$Min_AbsolutErr
rm(l_params_calibrated_sets)

# Compute derived probabilities from MAE parameter sets
compute_mae_probs <- function(mae) {
  probs <- list()
  for (loc in c("P", "D", "R")) {
    for (trans in c("S1S2", "S2S3", "S3S4")) {
      rate_name <- paste0("PreclinCancerProg_Exp_Rate_", trans, "_", loc)
      prob_name <- paste0("pCancerProg_", trans, "_", loc)
      probs[[prob_name]] <- 1 - exp(-mae[[rate_name]])
    }
    # Detection probabilities
    rS1 <- -log(1 - mae[[paste0("pSxDetS1_", loc)]])
    probs[[paste0("pSxDetS1_", loc)]] <- mae[[paste0("pSxDetS1_", loc)]]

    rS2 <- mae[[paste0("hr_SxDetS2S1_", loc)]] * rS1
    probs[[paste0("pSxDetS2_", loc)]] <- 1 - exp(-rS2)

    rS3 <- mae[[paste0("hr_SxDetS3S2_", loc)]] * rS2
    probs[[paste0("pSxDetS3_", loc)]] <- 1 - exp(-rS3)

    rS4 <- mae[[paste0("hr_SxDetS4S3_", loc)]] * rS3
    probs[[paste0("pSxDetS4_", loc)]] <- 1 - exp(-rS4)
  }
  probs
}

mae_probs_chile <- compute_mae_probs(mae_chile)
mae_probs_usa   <- compute_mae_probs(mae_usa)

summarize_params <- function(post_df, country, mae_probs) {
  results <- lapply(params_to_summarize, function(p) {
    vals <- post_df[[p]]
    data.frame(
      parameter  = p,
      country    = country,
      mae        = mae_probs[[p]],
      median     = median(vals),
      mean       = mean(vals),
      lb_95      = quantile(vals, 0.025),
      ub_95      = quantile(vals, 0.975),
      lb_50      = quantile(vals, 0.25),
      ub_50      = quantile(vals, 0.75),
      row.names  = NULL
    )
  })
  do.call(rbind, results)
}

df_summary <- rbind(
  summarize_params(post_chile, "Chile", mae_probs_chile),
  summarize_params(post_usa,   "USA",   mae_probs_usa)
)

# Add readable labels
df_summary$type <- case_when(
  grepl("pCancerProg", df_summary$parameter) ~ "Cancer Progression (annual prob)",
  grepl("pSxDet",      df_summary$parameter) ~ "Symptom Detection (prob)"
)

df_summary$transition <- case_when(
  grepl("S1S2|S1_", df_summary$parameter) ~ "Stage 1 → 2",
  grepl("S2S3|S2_", df_summary$parameter) ~ "Stage 2 → 3",
  grepl("S3S4|S3_", df_summary$parameter) ~ "Stage 3 → 4",
  grepl("S4_",      df_summary$parameter) ~ "Stage 4 (detection)"
)

df_summary$location <- case_when(
  grepl("_P$", df_summary$parameter) ~ "Proximal",
  grepl("_D$", df_summary$parameter) ~ "Distal",
  grepl("_R$", df_summary$parameter) ~ "Rectal"
)

# 5. Print summary table -------------------------------------------------------

cat("\n========================================================================\n")
cat("  Cancer Progression & Detection Parameters: Chile vs USA\n")
cat("  Annual transition probabilities from calibrated posteriors\n")
cat("========================================================================\n\n")

# Format for display
df_display <- df_summary %>%
  mutate(
    value_str = sprintf("MAE=%.4f | Med=%.4f [%.4f, %.4f]", mae, median, lb_95, ub_95)
  ) %>%
  select(type, transition, location, country, value_str) %>%
  pivot_wider(names_from = country, values_from = value_str) %>%
  arrange(type, transition, location)

print(as.data.frame(df_display), row.names = FALSE)

# 6. Save to CSV ---------------------------------------------------------------

write.csv(df_summary,
          file = file.path(output_folder, "cancer_progression_params_Chile_vs_USA.csv"),
          row.names = FALSE)

cat("\nSaved to:", file.path(output_folder, "cancer_progression_params_Chile_vs_USA.csv"), "\n")

# 7. Also export the full posterior draws for both countries --------------------

# Select only the relevant derived columns + original rates
cols_export <- c(
  paste0("PreclinCancerProg_Exp_Rate_", rep(c("S1S2", "S2S3", "S3S4"), each = 3), "_", c("P", "D", "R")),
  params_to_summarize
)

post_chile_export <- post_chile[, intersect(cols_export, colnames(post_chile))]
post_chile_export$country <- "Chile"

post_usa_export <- post_usa[, intersect(cols_export, colnames(post_usa))]
post_usa_export$country <- "USA"

df_draws <- rbind(post_chile_export, post_usa_export)

write.csv(df_draws,
          file = file.path(output_folder, "cancer_progression_draws_Chile_vs_USA.csv"),
          row.names = FALSE)

cat("Saved full draws to:", file.path(output_folder, "cancer_progression_draws_Chile_vs_USA.csv"), "\n")

# 8. Export Min Absolute Error calibrated parameter sets -----------------------

# Already loaded as mae_chile / mae_usa — pivot to long format
df_mae_chile_long <- as.data.frame(mae_chile) %>%
  pivot_longer(everything(), names_to = "parameter", values_to = "Chile")

df_mae_usa_long <- as.data.frame(mae_usa) %>%
  pivot_longer(everything(), names_to = "parameter", values_to = "USA")

df_mae <- inner_join(df_mae_chile_long, df_mae_usa_long, by = "parameter")

write.csv(df_mae,
          file = file.path(output_folder, "calibrated_params_MinAbsErr_Chile_vs_USA.csv"),
          row.names = FALSE)

cat("Saved Min Absolute Error params to:",
    file.path(output_folder, "calibrated_params_MinAbsErr_Chile_vs_USA.csv"), "\n")
