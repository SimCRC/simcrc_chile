save_plot_safe <- function(plot, filename, width = 10, height = 6, ...) {
  
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  
  ext  <- tools::file_ext(filename)
  base <- sub(paste0("\\.", ext, "$"), "", basename(filename))
  dirn <- dirname(filename)
  
  out <- file.path(dirn, paste0(base, ".", ext))
  
  i <- 1L
  while (file.exists(out)) {
    out <- file.path(dirn, paste0(base, "_v", i, ".", ext))
    i <- i + 1L
  }
  
  ggplot2::ggsave(
    filename = out,
    plot     = plot,
    width    = width,
    height   = height,
    ...
  )
  
  message("Saved plot to: ", out)
  invisible(out)
}




# Script to set the calibration parameters for each model
library(simcrc)
library(dplyr)
library(tidyr)
#========================================================================#
# 1. Parameter priors ----
#========================================================================#

# Adenoma-only female
l_params_priors_female_adenoma <- list(
  names = c(
    "IndividuakRiskMultVariance",
    "alpha_lesion_adenoma",
    "beta_age",
    "AdNaturalHistoryPropensity_Gaussian_Variance",
    "AdGrowth_Exp_Rate_1_to_6_P",
    "AdGrowth_Exp_Rate_1_to_6_D",
    "AdGrowth_Exp_Rate_1_to_6_R",
    "AdGrowth_Exp_Rate_6_to_10_P",
    "AdGrowth_Exp_Rate_6_to_10_D",
    "AdGrowth_Exp_Rate_6_to_10_R", 
    "PreclinCancerProg_Exp_Rate_S1S2_P",
    "PreclinCancerProg_Exp_Rate_S1S2_D",
    "PreclinCancerProg_Exp_Rate_S1S2_R",
    "PreclinCancerProg_Exp_Rate_S2S3_P",
    "PreclinCancerProg_Exp_Rate_S2S3_D",
    "PreclinCancerProg_Exp_Rate_S2S3_R",
    "PreclinCancerProg_Exp_Rate_S3S4_P",
    "PreclinCancerProg_Exp_Rate_S3S4_D",
    "PreclinCancerProg_Exp_Rate_S3S4_R",
    "CancerOnset_Gompertz_Shape_P",
    "CancerOnset_Gompertz_Shape_D",
    "CancerOnset_Gompertz_Shape_R",
    "CancerOnset_Gompertz_Rate_P",
    "CancerOnset_Gompertz_Rate_D",
    "CancerOnset_Gompertz_Rate_R",
    "pSxDetS1_P",
    "pSxDetS1_D",
    "pSxDetS1_R",
    "hr_SxDetS2S1_P",
    "hr_SxDetS2S1_D",
    "hr_SxDetS2S1_R",
    "hr_SxDetS3S2_P",
    "hr_SxDetS3S2_D",
    "hr_SxDetS3S2_R",
    "hr_SxDetS4S3_P",
    "hr_SxDetS4S3_D",
    "hr_SxDetS4S3_R"      
  ),
  lb = c(
    IndividuakRiskMultVariance = .01,
    alpha_lesion_adenoma = -8,
    beta_age = 0.02,
    AdNaturalHistoryPropensity_Gaussian_Variance = 0.30,
    AdGrowth_Exp_Rate_1_to_6_P = 0.01, 
    AdGrowth_Exp_Rate_1_to_6_D = 0.01,
    AdGrowth_Exp_Rate_1_to_6_R =  0.01, 
    AdGrowth_Exp_Rate_6_to_10_P = 0.005,
    AdGrowth_Exp_Rate_6_to_10_D = 0.02,
    AdGrowth_Exp_Rate_6_to_10_R = 0.01, 
    PreclinCancerProg_Exp_Rate_S1S2_P = 0.20, 
    PreclinCancerProg_Exp_Rate_S1S2_D = 0.278333*0.85,
    PreclinCancerProg_Exp_Rate_S1S2_R = 0.350170*0.85,
    PreclinCancerProg_Exp_Rate_S2S3_P = 0.246778*0.85, 
    PreclinCancerProg_Exp_Rate_S2S3_D = 0.276626*0.85, 
    PreclinCancerProg_Exp_Rate_S2S3_R = 0.324122*0.85, 
    PreclinCancerProg_Exp_Rate_S3S4_P = 0.394293*0.85, 
    PreclinCancerProg_Exp_Rate_S3S4_D = 0.399649*0.85, 
    PreclinCancerProg_Exp_Rate_S3S4_R = 0.419631*0.85, 
    CancerOnset_Gompertz_Shape_P = 0.00001,
    CancerOnset_Gompertz_Shape_D = 0.001,
    CancerOnset_Gompertz_Shape_R = 0.001,
    CancerOnset_Gompertz_Rate_P = 0.0001,
    CancerOnset_Gompertz_Rate_D = 0.0005,
    CancerOnset_Gompertz_Rate_R = 0.0005,
    pSxDetS1_P = 0.001,
    pSxDetS1_D = 0.001,
    pSxDetS1_R = 0.001,
    hr_SxDetS2S1_P = 1,
    hr_SxDetS2S1_D = 1,
    hr_SxDetS2S1_R = 1,
    hr_SxDetS3S2_P = 1,
    hr_SxDetS3S2_D = 1,
    hr_SxDetS3S2_R = 1,
    hr_SxDetS4S3_P = 1,
    hr_SxDetS4S3_D = 1,
    hr_SxDetS4S3_R = 1
  ),
  ub = c(
    IndividuakRiskMultVariance = .8,
    alpha_lesion_adenoma = -5,
    beta_age = 0.06,
    AdNaturalHistoryPropensity_Gaussian_Variance = 0.5,
    AdGrowth_Exp_Rate_1_to_6_P = 0.23, 
    AdGrowth_Exp_Rate_1_to_6_D = 0.2,
    AdGrowth_Exp_Rate_1_to_6_R = 0.48210154*1.5, 
    AdGrowth_Exp_Rate_6_to_10_P = 0.13,
    AdGrowth_Exp_Rate_6_to_10_D = 0.8,
    AdGrowth_Exp_Rate_6_to_10_R = 0.08189034*5.5,
    PreclinCancerProg_Exp_Rate_S1S2_P = 0.265947*1.15, 
    PreclinCancerProg_Exp_Rate_S1S2_D = 0.278333*1.15, 
    PreclinCancerProg_Exp_Rate_S1S2_R = 0.350170*1.15, 
    PreclinCancerProg_Exp_Rate_S2S3_P = 0.246778*1.15, 
    PreclinCancerProg_Exp_Rate_S2S3_D = 0.4, 
    PreclinCancerProg_Exp_Rate_S2S3_R = 0.5, 
    PreclinCancerProg_Exp_Rate_S3S4_P = 0.52, 
    PreclinCancerProg_Exp_Rate_S3S4_D = 0.7, 
    PreclinCancerProg_Exp_Rate_S3S4_R = 0.9, 
    CancerOnset_Gompertz_Shape_P = 0.0001,   #change for approach with exponential
    CancerOnset_Gompertz_Shape_D = 0.02,  #change for approach with exponential
    CancerOnset_Gompertz_Shape_R = 0.025,  #change for approach with exponential
    CancerOnset_Gompertz_Rate_P = 0.08, #change for approach with exponential
    CancerOnset_Gompertz_Rate_D = 0.04, #0.1 change for approach with exponential
    CancerOnset_Gompertz_Rate_R = 0.18, #change for approach with exponential
    pSxDetS1_P = 0.192,
    pSxDetS1_D = 0.192,
    pSxDetS1_R = 0.192,
    hr_SxDetS2S1_P = 24,
    hr_SxDetS2S1_D = 24,
    hr_SxDetS2S1_R = 24,
    hr_SxDetS3S2_P = 7,
    hr_SxDetS3S2_D = 7,
    hr_SxDetS3S2_R = 7,
    hr_SxDetS4S3_P = 5,
    hr_SxDetS4S3_D = 5,
    hr_SxDetS4S3_R = 5
  )
)

# Adenoma-only male
l_params_priors_male_adenoma <- list(
  names = c(
    "IndividuakRiskMultVariance",
    "alpha_lesion_adenoma",
    "beta_age",
    "AdNaturalHistoryPropensity_Gaussian_Variance",
    "AdGrowth_Exp_Rate_1_to_6_P",
    "AdGrowth_Exp_Rate_1_to_6_D",
    "AdGrowth_Exp_Rate_1_to_6_R",
    "AdGrowth_Exp_Rate_6_to_10_P",
    "AdGrowth_Exp_Rate_6_to_10_D",
    "AdGrowth_Exp_Rate_6_to_10_R", 
    "PreclinCancerProg_Exp_Rate_S1S2_P",
    "PreclinCancerProg_Exp_Rate_S1S2_D",
    "PreclinCancerProg_Exp_Rate_S1S2_R",
    "PreclinCancerProg_Exp_Rate_S2S3_P",
    "PreclinCancerProg_Exp_Rate_S2S3_D",
    "PreclinCancerProg_Exp_Rate_S2S3_R",
    "PreclinCancerProg_Exp_Rate_S3S4_P",
    "PreclinCancerProg_Exp_Rate_S3S4_D",
    "PreclinCancerProg_Exp_Rate_S3S4_R",
    "CancerOnset_Gompertz_Shape_P",
    "CancerOnset_Gompertz_Shape_D",
    "CancerOnset_Gompertz_Shape_R",
    "CancerOnset_Gompertz_Rate_P",
    "CancerOnset_Gompertz_Rate_D",
    "CancerOnset_Gompertz_Rate_R",
    "pSxDetS1_P",
    "pSxDetS1_D",
    "pSxDetS1_R",
    "hr_SxDetS2S1_P",
    "hr_SxDetS2S1_D",
    "hr_SxDetS2S1_R",
    "hr_SxDetS3S2_P",
    "hr_SxDetS3S2_D",
    "hr_SxDetS3S2_R",
    "hr_SxDetS4S3_P",
    "hr_SxDetS4S3_D",
    "hr_SxDetS4S3_R"      
  ),
  lb = c(
    IndividuakRiskMultVariance = .01,
    alpha_lesion_adenoma = -8,
    beta_age = 0.02,
    AdNaturalHistoryPropensity_Gaussian_Variance = 0.30,
    AdGrowth_Exp_Rate_1_to_6_P = 0.01, 
    AdGrowth_Exp_Rate_1_to_6_D = 0.01,
    AdGrowth_Exp_Rate_1_to_6_R =  0.01, 
    AdGrowth_Exp_Rate_6_to_10_P = 0.005,
    AdGrowth_Exp_Rate_6_to_10_D = 0.02,
    AdGrowth_Exp_Rate_6_to_10_R = 0.01, 
    PreclinCancerProg_Exp_Rate_S1S2_P = 0.20, 
    PreclinCancerProg_Exp_Rate_S1S2_D = 0.278333*0.85,
    PreclinCancerProg_Exp_Rate_S1S2_R = 0.350170*0.85,
    PreclinCancerProg_Exp_Rate_S2S3_P = 0.246778*0.85, 
    PreclinCancerProg_Exp_Rate_S2S3_D = 0.276626*0.85, 
    PreclinCancerProg_Exp_Rate_S2S3_R = 0.324122*0.85, 
    PreclinCancerProg_Exp_Rate_S3S4_P = 0.394293*0.85, 
    PreclinCancerProg_Exp_Rate_S3S4_D = 0.399649*0.85, 
    PreclinCancerProg_Exp_Rate_S3S4_R = 0.419631*0.85, 
    CancerOnset_Gompertz_Shape_P = 0.00001,
    CancerOnset_Gompertz_Shape_D = 0.001,
    CancerOnset_Gompertz_Shape_R = 0.001,
    CancerOnset_Gompertz_Rate_P = 0.0001,
    CancerOnset_Gompertz_Rate_D = 0.0005,
    CancerOnset_Gompertz_Rate_R = 0.0005,
    pSxDetS1_P = 0.001,
    pSxDetS1_D = 0.001,
    pSxDetS1_R = 0.001,
    hr_SxDetS2S1_P = 1,
    hr_SxDetS2S1_D = 1,
    hr_SxDetS2S1_R = 1,
    hr_SxDetS3S2_P = 1,
    hr_SxDetS3S2_D = 1,
    hr_SxDetS3S2_R = 1,
    hr_SxDetS4S3_P = 1,
    hr_SxDetS4S3_D = 1,
    hr_SxDetS4S3_R = 1
  ),
  ub = c(
    IndividuakRiskMultVariance = .8,
    alpha_lesion_adenoma = -5,
    beta_age = 0.06,
    AdNaturalHistoryPropensity_Gaussian_Variance = 0.5,
    AdGrowth_Exp_Rate_1_to_6_P = 0.23, 
    AdGrowth_Exp_Rate_1_to_6_D = 0.2,
    AdGrowth_Exp_Rate_1_to_6_R = 0.48210154*1.5, 
    AdGrowth_Exp_Rate_6_to_10_P = 0.13,
    AdGrowth_Exp_Rate_6_to_10_D = 0.8,
    AdGrowth_Exp_Rate_6_to_10_R = 0.08189034*5.5,
    PreclinCancerProg_Exp_Rate_S1S2_P = 0.265947*1.15, 
    PreclinCancerProg_Exp_Rate_S1S2_D = 0.278333*1.15, 
    PreclinCancerProg_Exp_Rate_S1S2_R = 0.350170*1.15, 
    PreclinCancerProg_Exp_Rate_S2S3_P = 0.246778*1.15, 
    PreclinCancerProg_Exp_Rate_S2S3_D = 0.4, 
    PreclinCancerProg_Exp_Rate_S2S3_R = 0.5, 
    PreclinCancerProg_Exp_Rate_S3S4_P = 0.52, 
    PreclinCancerProg_Exp_Rate_S3S4_D = 0.7, 
    PreclinCancerProg_Exp_Rate_S3S4_R = 0.9, 
    CancerOnset_Gompertz_Shape_P = 0.0001,   #change for approach with exponential
    CancerOnset_Gompertz_Shape_D = 0.02,  #change for approach with exponential
    CancerOnset_Gompertz_Shape_R = 0.025,  #change for approach with exponential
    CancerOnset_Gompertz_Rate_P = 0.08, #change for approach with exponential
    CancerOnset_Gompertz_Rate_D = 0.04, #0.1 change for approach with exponential
    CancerOnset_Gompertz_Rate_R = 0.18, #change for approach with exponential
    pSxDetS1_P = 0.192,
    pSxDetS1_D = 0.192,
    pSxDetS1_R = 0.192,
    hr_SxDetS2S1_P = 24,
    hr_SxDetS2S1_D = 24,
    hr_SxDetS2S1_R = 24,
    hr_SxDetS3S2_P = 7,
    hr_SxDetS3S2_D = 7,
    hr_SxDetS3S2_R = 7,
    hr_SxDetS4S3_P = 5,
    hr_SxDetS4S3_D = 5,
    hr_SxDetS4S3_R = 5
  )
)

# female both
l_params_priors_female_both <- list(
  names = c(
    "IndividuakRiskMultVariance", 
    "alpha_lesion_adenoma",
    "beta_age",
    "AdNaturalHistoryPropensity_Gaussian_Variance",
    "AdGrowth_Exp_Rate_1_to_6_P",
    "AdGrowth_Exp_Rate_1_to_6_D",
    "AdGrowth_Exp_Rate_1_to_6_R",
    "AdGrowth_Exp_Rate_6_to_10_P",
    "AdGrowth_Exp_Rate_6_to_10_D",
    "AdGrowth_Exp_Rate_6_to_10_R", 
    "PreclinCancerProg_Exp_Rate_S1S2_P",
    "PreclinCancerProg_Exp_Rate_S1S2_D",
    "PreclinCancerProg_Exp_Rate_S1S2_R",
    "PreclinCancerProg_Exp_Rate_S2S3_P",
    "PreclinCancerProg_Exp_Rate_S2S3_D",
    "PreclinCancerProg_Exp_Rate_S2S3_R",
    "PreclinCancerProg_Exp_Rate_S3S4_P",
    "PreclinCancerProg_Exp_Rate_S3S4_D",
    "PreclinCancerProg_Exp_Rate_S3S4_R",
    "CancerOnset_Gompertz_Shape_P",
    "CancerOnset_Gompertz_Shape_D",
    "CancerOnset_Gompertz_Shape_R",
    "CancerOnset_Gompertz_Rate_P",
    "CancerOnset_Gompertz_Rate_D",
    "CancerOnset_Gompertz_Rate_R",
    "pSxDetS1_P",
    "pSxDetS1_D",
    "pSxDetS1_R",
    "hr_SxDetS2S1_P",
    "hr_SxDetS2S1_D",
    "hr_SxDetS2S1_R",
    "hr_SxDetS3S2_P",
    "hr_SxDetS3S2_D",
    "hr_SxDetS3S2_R",
    "hr_SxDetS4S3_P",
    "hr_SxDetS4S3_D",
    "hr_SxDetS4S3_R",           
    "alpha_lesion_ssp", #SSP 
    "beta_age_ssp",  # SSP
    "hr_SSPGrowth_1_to_6_P",  ##SSP
    "hr_SSPGrowth_6_to_10_P",  ##SSP
    "hr_SSPGrowth_1_to_6_D",  ##SSP
    "hr_SSPGrowth_6_to_10_D",  ##SSP
    "hr_SSPGrowth_1_to_6_R",  ##SSP
    "hr_SSPGrowth_6_to_10_R",  ##SSP
    "hr_SSPCancerOnset_P"  ##SSP
  ),
  lb = c(
    IndividuakRiskMultVariance = 0.01, 
    alpha_lesion_adenoma = -8, 
    beta_age =  0.02, 
    AdNaturalHistoryPropensity_Gaussian_Variance = 0.30,
    AdGrowth_Exp_Rate_1_to_6_P = 0.001, 
    AdGrowth_Exp_Rate_1_to_6_D = 0.001,
    AdGrowth_Exp_Rate_1_to_6_R = 0.001, 
    AdGrowth_Exp_Rate_6_to_10_P = 0.001,
    AdGrowth_Exp_Rate_6_to_10_D = 0.001,
    AdGrowth_Exp_Rate_6_to_10_R = 0.001, 
    PreclinCancerProg_Exp_Rate_S1S2_P = 0.20, #0.228536835,
    PreclinCancerProg_Exp_Rate_S1S2_D = 0.20, #0.270014625,
    PreclinCancerProg_Exp_Rate_S1S2_R = 0.30, #0.260442195,
    PreclinCancerProg_Exp_Rate_S2S3_P = 0.20, #0.0000001,
    PreclinCancerProg_Exp_Rate_S2S3_D = 0.20, #0.0000001,
    PreclinCancerProg_Exp_Rate_S2S3_R = 0.20, #0.0000001,
    PreclinCancerProg_Exp_Rate_S3S4_P = 0.30, #0.0000001,
    PreclinCancerProg_Exp_Rate_S3S4_D = 0.30, #0.0000001,
    PreclinCancerProg_Exp_Rate_S3S4_R = 0.30, #0.0000001,
    CancerOnset_Gompertz_Shape_P = 0.00001,
    CancerOnset_Gompertz_Shape_D = 0.0001,
    CancerOnset_Gompertz_Shape_R = 0.0001,
    CancerOnset_Gompertz_Rate_P = 0.001,
    CancerOnset_Gompertz_Rate_D = 0.001,
    CancerOnset_Gompertz_Rate_R = 0.001,
    pSxDetS1_P = 0.001,
    pSxDetS1_D = 0.001,
    pSxDetS1_R = 0.001,
    hr_SxDetS2S1_P = 1,
    hr_SxDetS2S1_D = 1,
    hr_SxDetS2S1_R = 1,
    hr_SxDetS3S2_P = 1,
    hr_SxDetS3S2_D = 1,
    hr_SxDetS3S2_R = 1,
    hr_SxDetS4S3_P = 1,
    hr_SxDetS4S3_D = 1,
    hr_SxDetS4S3_R = 1,
    alpha_lesion_ssp = -10, #SSP 
    beta_age_ssp  = -0.02,  # SSP 
    hr_SSPGrowth_1_to_6_P = 0.01,  ##SSP
    hr_SSPGrowth_6_to_10_P = 0.01,  ##SSP
    hr_SSPGrowth_1_to_6_D = 0.01,  ##SSP
    hr_SSPGrowth_6_to_10_D = 0.01,  ##SSP
    hr_SSPGrowth_1_to_6_R = 0.01,  ##SSP
    hr_SSPGrowth_6_to_10_R = 0.01,  ##SSP
    hr_SSPCancerOnset_P = 0.8  ##SSP
  ),
  ub = c(
    IndividuakRiskMultVariance = 0.8, 
    alpha_lesion_adenoma = -5, 
    beta_age =  0.06, 
    AdNaturalHistoryPropensity_Gaussian_Variance = 0.5,
    AdGrowth_Exp_Rate_1_to_6_P = 0.1, 
    AdGrowth_Exp_Rate_1_to_6_D = 0.2,
    AdGrowth_Exp_Rate_1_to_6_R = 0.3, 
    AdGrowth_Exp_Rate_6_to_10_P = 0.06,
    AdGrowth_Exp_Rate_6_to_10_D = 0.2,
    AdGrowth_Exp_Rate_6_to_10_R = 0.45,
    PreclinCancerProg_Exp_Rate_S1S2_P = 0.35, 
    PreclinCancerProg_Exp_Rate_S1S2_D = 0.35, 
    PreclinCancerProg_Exp_Rate_S1S2_R = 0.40, 
    PreclinCancerProg_Exp_Rate_S2S3_P = 0.30, 
    PreclinCancerProg_Exp_Rate_S2S3_D = 0.40, 
    PreclinCancerProg_Exp_Rate_S2S3_R = 0.50, 
    PreclinCancerProg_Exp_Rate_S3S4_P = 0.52, 
    PreclinCancerProg_Exp_Rate_S3S4_D = 0.7, 
    PreclinCancerProg_Exp_Rate_S3S4_R = 0.9, 
    CancerOnset_Gompertz_Shape_P = 0.0001,   #change for approach with exponential
    CancerOnset_Gompertz_Shape_D = 0.02,  #change for approach with exponential
    CancerOnset_Gompertz_Shape_R = 0.025,  #change for approach with exponential
    CancerOnset_Gompertz_Rate_P = 0.1, #change for approach with exponential
    CancerOnset_Gompertz_Rate_D = 0.04, #0.1 change for approach with exponential
    CancerOnset_Gompertz_Rate_R = 0.2, #change for approach with exponential
    pSxDetS1_P = 0.192,
    pSxDetS1_D = 0.192,
    pSxDetS1_R = 0.192,
    hr_SxDetS2S1_P = 24,
    hr_SxDetS2S1_D = 24,
    hr_SxDetS2S1_R = 24,
    hr_SxDetS3S2_P = 7,
    hr_SxDetS3S2_D = 7,
    hr_SxDetS3S2_R = 7,
    hr_SxDetS4S3_P = 7,
    hr_SxDetS4S3_D = 7,
    hr_SxDetS4S3_R = 7,
    alpha_lesion_ssp = -5, #SSP 
    beta_age_ssp  = 0.06,  # SSP 
    hr_SSPGrowth_1_to_6_P = 2.5,  ##SSP
    hr_SSPGrowth_6_to_10_P = 2.5,  ##SSP
    hr_SSPGrowth_1_to_6_D = 2.5,  ##SSP
    hr_SSPGrowth_6_to_10_D = 2.5,  ##SSP
    hr_SSPGrowth_1_to_6_R = 2.5,  ##SSP
    hr_SSPGrowth_6_to_10_R = 2.5,  ##SSP
    hr_SSPCancerOnset_P = 2.5  ##SSP
  ) 
)


dt_calibrated_posteriors_SimCRC_v0_12_0_1_Ad_F <- readr::read_csv("data-raw/dt_calibrated_posteriors_SimCRC_v0.12.0.1_Ad_F.csv")

#Get min and max from the calibrated posteriors to use as bounds in the calibration

df_params_US_posteriors <- dt_calibrated_posteriors_SimCRC_v0_12_0_1_Ad_F %>%
  pivot_longer(cols = everything(),
               names_to = "param",
               values_to = "value") %>%
  dplyr::group_by(param) %>%
  dplyr::summarise(
    lb = quantile(value, probs = 0.025, na.rm = TRUE),
    ub = quantile(value, probs = 0.975, na.rm = TRUE)
  ) %>%
  ungroup()

l_params_priors_adenoma_Chile_v2 <- list(
  names = df_params_US_posteriors$param,
  lb = setNames(df_params_US_posteriors$lb, df_params_US_posteriors$param),
  ub = setNames(df_params_US_posteriors$ub, df_params_US_posteriors$param)
)

# Modify a single parameter's bounds
l_params_priors_adenoma_Chile_v2$lb["alpha_lesion_adenoma"] <- -10
l_params_priors_adenoma_Chile_v2$ub["alpha_lesion_adenoma"] <- -5


l_params_priors_adenoma_Chile_v2$lb["hr_SxDetS2S1_P"] <- 1.2
l_params_priors_adenoma_Chile_v2$lb["hr_SxDetS2S1_D"] <- 1.2
l_params_priors_adenoma_Chile_v2$lb["hr_SxDetS2S1_R"] <- 1.2
l_params_priors_adenoma_Chile_v2$lb["hr_SxDetS3S2_P"] <- 4
l_params_priors_adenoma_Chile_v2$lb["hr_SxDetS3S2_D"] <- 4
l_params_priors_adenoma_Chile_v2$lb["hr_SxDetS3S2_R"] <- 4
l_params_priors_adenoma_Chile_v2$lb["hr_SxDetS4S3_P"] <- 8
l_params_priors_adenoma_Chile_v2$lb["hr_SxDetS4S3_D"] <- 8
l_params_priors_adenoma_Chile_v2$lb["hr_SxDetS4S3_R"] <- 8
  
l_params_priors_adenoma_Chile_v2$ub["hr_SxDetS2S1_P"] <- 8
l_params_priors_adenoma_Chile_v2$ub["hr_SxDetS2S1_D"] <- 8
l_params_priors_adenoma_Chile_v2$ub["hr_SxDetS2S1_R"] <- 8
l_params_priors_adenoma_Chile_v2$ub["hr_SxDetS3S2_P"] <- 8.5
l_params_priors_adenoma_Chile_v2$ub["hr_SxDetS3S2_D"] <- 8.5
l_params_priors_adenoma_Chile_v2$ub["hr_SxDetS3S2_R"] <- 8.5
l_params_priors_adenoma_Chile_v2$ub["hr_SxDetS4S3_P"] <- 12
l_params_priors_adenoma_Chile_v2$ub["hr_SxDetS4S3_D"] <- 12
l_params_priors_adenoma_Chile_v2$ub["hr_SxDetS4S3_R"] <- 12


l_params_priors_adenoma_Chile_v2$lb["pSxDetS1_P"]  <-  0.0001
l_params_priors_adenoma_Chile_v2$lb["pSxDetS1_D"]  <-  0.0001
l_params_priors_adenoma_Chile_v2$lb["pSxDetS1_R"]  <-  0.0001


l_params_priors_adenoma_Chile_v2$ub["pSxDetS1_P"]  <-  0.192*.2
l_params_priors_adenoma_Chile_v2$ub["pSxDetS1_D"]  <-  0.192*.2
l_params_priors_adenoma_Chile_v2$ub["pSxDetS1_R"]  <-  0.192*.2


l_params_priors_adenoma_Chile_v2$ub["beta_age"]  <-  0.03474471*2

l_params_priors_adenoma_Chile <- l_params_priors_adenoma_Chile_v2


l_params_priors_adenoma_Chile$names <- l_params_priors_adenoma_Chile$names[l_params_priors_adenoma_Chile$names != "id_draw"]
l_params_priors_adenoma_Chile$lb <- l_params_priors_adenoma_Chile$lb[names(l_params_priors_adenoma_Chile$lb) != "id_draw"]
l_params_priors_adenoma_Chile$ub <- l_params_priors_adenoma_Chile$ub[names(l_params_priors_adenoma_Chile$ub) != "id_draw"]






# Adenoma-only Chile
l_params_priors_adenoma_Chile <- list(
  names = c(
    "IndividuakRiskMultVariance",
    "alpha_lesion_adenoma",
    "beta_age",
    "AdNaturalHistoryPropensity_Gaussian_Variance",
    "AdGrowth_Exp_Rate_1_to_6_P",
    "AdGrowth_Exp_Rate_1_to_6_D",
    "AdGrowth_Exp_Rate_1_to_6_R",
    "AdGrowth_Exp_Rate_6_to_10_P",
    "AdGrowth_Exp_Rate_6_to_10_D",
    "AdGrowth_Exp_Rate_6_to_10_R", 
    "PreclinCancerProg_Exp_Rate_S1S2_P",
    "PreclinCancerProg_Exp_Rate_S1S2_D",
    "PreclinCancerProg_Exp_Rate_S1S2_R",
    "PreclinCancerProg_Exp_Rate_S2S3_P",
    "PreclinCancerProg_Exp_Rate_S2S3_D",
    "PreclinCancerProg_Exp_Rate_S2S3_R",
    "PreclinCancerProg_Exp_Rate_S3S4_P",
    "PreclinCancerProg_Exp_Rate_S3S4_D",
    "PreclinCancerProg_Exp_Rate_S3S4_R",
    "CancerOnset_Gompertz_Shape_P",
    "CancerOnset_Gompertz_Shape_D",
    "CancerOnset_Gompertz_Shape_R",
    "CancerOnset_Gompertz_Rate_P",
    "CancerOnset_Gompertz_Rate_D",
    "CancerOnset_Gompertz_Rate_R",
    "pSxDetS1_P",
    "pSxDetS1_D",
    "pSxDetS1_R",
    "hr_SxDetS2S1_P",
    "hr_SxDetS2S1_D",
    "hr_SxDetS2S1_R",
    "hr_SxDetS3S2_P",
    "hr_SxDetS3S2_D",
    "hr_SxDetS3S2_R",
    "hr_SxDetS4S3_P",
    "hr_SxDetS4S3_D",
    "hr_SxDetS4S3_R"      
  ),
  lb = c(
    IndividuakRiskMultVariance = .01,
    alpha_lesion_adenoma = -10,
    beta_age = 0.01,
    AdNaturalHistoryPropensity_Gaussian_Variance = 0.30,
    AdGrowth_Exp_Rate_1_to_6_P = 0.001, 
    AdGrowth_Exp_Rate_1_to_6_D = 0.001,
    AdGrowth_Exp_Rate_1_to_6_R =  0.001, 
    AdGrowth_Exp_Rate_6_to_10_P = 0.005,
    AdGrowth_Exp_Rate_6_to_10_D = 0.002,
    AdGrowth_Exp_Rate_6_to_10_R = 0.001, 
    PreclinCancerProg_Exp_Rate_S1S2_P = 0.20, 
    PreclinCancerProg_Exp_Rate_S1S2_D = 0.278333*0.1,
    PreclinCancerProg_Exp_Rate_S1S2_R = 0.350170*0.1,
    PreclinCancerProg_Exp_Rate_S2S3_P = 0.246778*0.1, 
    PreclinCancerProg_Exp_Rate_S2S3_D = 0.276626*0.1, 
    PreclinCancerProg_Exp_Rate_S2S3_R = 0.324122*0.1, 
    PreclinCancerProg_Exp_Rate_S3S4_P = 0.394293*0.1, 
    PreclinCancerProg_Exp_Rate_S3S4_D = 0.399649*0.1, 
    PreclinCancerProg_Exp_Rate_S3S4_R = 0.419631*0.1, 
    CancerOnset_Gompertz_Shape_P = 0.00001,
    CancerOnset_Gompertz_Shape_D = 0.00001,
    CancerOnset_Gompertz_Shape_R = 0.00001,
    CancerOnset_Gompertz_Rate_P = 0.00001,
    CancerOnset_Gompertz_Rate_D = 0.00001,
    CancerOnset_Gompertz_Rate_R = 0.00001,
    pSxDetS1_P = 0.0001,
    pSxDetS1_D = 0.0001,
    pSxDetS1_R = 0.0001,
    hr_SxDetS2S1_P = 1.2,
    hr_SxDetS2S1_D = 1.2,
    hr_SxDetS2S1_R = 1.2,
    hr_SxDetS3S2_P = 4.0,
    hr_SxDetS3S2_D = 4.0,
    hr_SxDetS3S2_R = 4.0,
    hr_SxDetS4S3_P = 8.0,
    hr_SxDetS4S3_D = 8.0,
    hr_SxDetS4S3_R = 8.0
  ),
  ub = c(
    IndividuakRiskMultVariance = .8,
    alpha_lesion_adenoma = -5,
    beta_age = 0.06,
    AdNaturalHistoryPropensity_Gaussian_Variance = 0.5,
    AdGrowth_Exp_Rate_1_to_6_P = 0.23*.8, 
    AdGrowth_Exp_Rate_1_to_6_D = 0.2*.8,
    AdGrowth_Exp_Rate_1_to_6_R = 0.72315231*.8, 
    AdGrowth_Exp_Rate_6_to_10_P = 0.13*.8,
    AdGrowth_Exp_Rate_6_to_10_D = 0.8*.8,
    AdGrowth_Exp_Rate_6_to_10_R = 0.45039687*.8,
    PreclinCancerProg_Exp_Rate_S1S2_P = 0.265947*1.15, 
    PreclinCancerProg_Exp_Rate_S1S2_D = 0.278333*1.15, 
    PreclinCancerProg_Exp_Rate_S1S2_R = 0.350170*1.15, 
    PreclinCancerProg_Exp_Rate_S2S3_P = 0.246778*1.15, 
    PreclinCancerProg_Exp_Rate_S2S3_D = 0.4, 
    PreclinCancerProg_Exp_Rate_S2S3_R = 0.5, 
    PreclinCancerProg_Exp_Rate_S3S4_P = 0.52, 
    PreclinCancerProg_Exp_Rate_S3S4_D = 0.7, 
    PreclinCancerProg_Exp_Rate_S3S4_R = 0.9, 
    CancerOnset_Gompertz_Shape_P = 0.0001*.8,   #change for approach with exponential
    CancerOnset_Gompertz_Shape_D = 0.02*.8,  #change for approach with exponential
    CancerOnset_Gompertz_Shape_R = 0.025*.8,  #change for approach with exponential
    CancerOnset_Gompertz_Rate_P = 0.08*.8, #change for approach with exponential
    CancerOnset_Gompertz_Rate_D = 0.04*.8, #0.1 change for approach with exponential
    CancerOnset_Gompertz_Rate_R = 0.18*.8, #change for approach with exponential
    pSxDetS1_P = 0.192*.2,
    pSxDetS1_D = 0.192*.2,
    pSxDetS1_R = 0.192*.2,
    hr_SxDetS2S1_P = 6,
    hr_SxDetS2S1_D = 6,
    hr_SxDetS2S1_R = 6,
    hr_SxDetS3S2_P = 6.5,
    hr_SxDetS3S2_D = 6.5,
    hr_SxDetS3S2_R = 6.5,
    hr_SxDetS4S3_P = 12,
    hr_SxDetS4S3_D = 12,
    hr_SxDetS4S3_R = 12
  )
)






# Adenoma-only female
l_params_priors_adenoma_Chile <- list(
  names = c(
    "IndividuakRiskMultVariance",
    "alpha_lesion_adenoma",
    "beta_age",
    "AdNaturalHistoryPropensity_Gaussian_Variance",
    "AdGrowth_Exp_Rate_1_to_6_P",
    "AdGrowth_Exp_Rate_1_to_6_D",
    "AdGrowth_Exp_Rate_1_to_6_R",
    "AdGrowth_Exp_Rate_6_to_10_P",
    "AdGrowth_Exp_Rate_6_to_10_D",
    "AdGrowth_Exp_Rate_6_to_10_R", 
    "PreclinCancerProg_Exp_Rate_S1S2_P",
    "PreclinCancerProg_Exp_Rate_S1S2_D",
    "PreclinCancerProg_Exp_Rate_S1S2_R",
    "PreclinCancerProg_Exp_Rate_S2S3_P",
    "PreclinCancerProg_Exp_Rate_S2S3_D",
    "PreclinCancerProg_Exp_Rate_S2S3_R",
    "PreclinCancerProg_Exp_Rate_S3S4_P",
    "PreclinCancerProg_Exp_Rate_S3S4_D",
    "PreclinCancerProg_Exp_Rate_S3S4_R",
    "CancerOnset_Gompertz_Shape_P",
    "CancerOnset_Gompertz_Shape_D",
    "CancerOnset_Gompertz_Shape_R",
    "CancerOnset_Gompertz_Rate_P",
    "CancerOnset_Gompertz_Rate_D",
    "CancerOnset_Gompertz_Rate_R",
    "pSxDetS1_P",
    "pSxDetS1_D",
    "pSxDetS1_R",
    "hr_SxDetS2S1_P",
    "hr_SxDetS2S1_D",
    "hr_SxDetS2S1_R",
    "hr_SxDetS3S2_P",
    "hr_SxDetS3S2_D",
    "hr_SxDetS3S2_R",
    "hr_SxDetS4S3_P",
    "hr_SxDetS4S3_D",
    "hr_SxDetS4S3_R"      
  ),
  lb = c(
    IndividuakRiskMultVariance = .01,
    alpha_lesion_adenoma = -8,
    beta_age = 0.02,
    AdNaturalHistoryPropensity_Gaussian_Variance = 0.30,
    AdGrowth_Exp_Rate_1_to_6_P = 0.01, 
    AdGrowth_Exp_Rate_1_to_6_D = 0.01,
    AdGrowth_Exp_Rate_1_to_6_R =  0.01, 
    AdGrowth_Exp_Rate_6_to_10_P = 0.005,
    AdGrowth_Exp_Rate_6_to_10_D = 0.02,
    AdGrowth_Exp_Rate_6_to_10_R = 0.01, 
    PreclinCancerProg_Exp_Rate_S1S2_P = 0.20, 
    PreclinCancerProg_Exp_Rate_S1S2_D = 0.278333*0.85,
    PreclinCancerProg_Exp_Rate_S1S2_R = 0.350170*0.85,
    PreclinCancerProg_Exp_Rate_S2S3_P = 0.246778*0.85, 
    PreclinCancerProg_Exp_Rate_S2S3_D = 0.276626*0.85, 
    PreclinCancerProg_Exp_Rate_S2S3_R = 0.324122*0.85, 
    PreclinCancerProg_Exp_Rate_S3S4_P = 0.394293*0.85, 
    PreclinCancerProg_Exp_Rate_S3S4_D = 0.399649*0.85, 
    PreclinCancerProg_Exp_Rate_S3S4_R = 0.419631*0.85, 
    CancerOnset_Gompertz_Shape_P = 0.00001,
    CancerOnset_Gompertz_Shape_D = 0.001,
    CancerOnset_Gompertz_Shape_R = 0.001,
    CancerOnset_Gompertz_Rate_P = 0.0001,
    CancerOnset_Gompertz_Rate_D = 0.0005,
    CancerOnset_Gompertz_Rate_R = 0.0005,
    pSxDetS1_P = 0.001,
    pSxDetS1_D = 0.001,
    pSxDetS1_R = 0.001,
    hr_SxDetS2S1_P = 1,
    hr_SxDetS2S1_D = 1,
    hr_SxDetS2S1_R = 1,
    hr_SxDetS3S2_P = 1,
    hr_SxDetS3S2_D = 1,
    hr_SxDetS3S2_R = 1,
    hr_SxDetS4S3_P = 1,
    hr_SxDetS4S3_D = 1,
    hr_SxDetS4S3_R = 1
  ),
  ub = c(
    IndividuakRiskMultVariance = .8,
    alpha_lesion_adenoma = -5,
    beta_age = 0.06,
    AdNaturalHistoryPropensity_Gaussian_Variance = 0.5,
    AdGrowth_Exp_Rate_1_to_6_P = 0.23, 
    AdGrowth_Exp_Rate_1_to_6_D = 0.2,
    AdGrowth_Exp_Rate_1_to_6_R = 0.48210154*1.5, 
    AdGrowth_Exp_Rate_6_to_10_P = 0.13,
    AdGrowth_Exp_Rate_6_to_10_D = 0.8,
    AdGrowth_Exp_Rate_6_to_10_R = 0.08189034*5.5,
    PreclinCancerProg_Exp_Rate_S1S2_P = 0.265947*1.15, 
    PreclinCancerProg_Exp_Rate_S1S2_D = 0.278333*1.15, 
    PreclinCancerProg_Exp_Rate_S1S2_R = 0.350170*1.15, 
    PreclinCancerProg_Exp_Rate_S2S3_P = 0.246778*1.15, 
    PreclinCancerProg_Exp_Rate_S2S3_D = 0.4, 
    PreclinCancerProg_Exp_Rate_S2S3_R = 0.5, 
    PreclinCancerProg_Exp_Rate_S3S4_P = 0.52, 
    PreclinCancerProg_Exp_Rate_S3S4_D = 0.7, 
    PreclinCancerProg_Exp_Rate_S3S4_R = 0.9, 
    CancerOnset_Gompertz_Shape_P = 0.0001,   #change for approach with exponential
    CancerOnset_Gompertz_Shape_D = 0.02,  #change for approach with exponential
    CancerOnset_Gompertz_Shape_R = 0.025,  #change for approach with exponential
    CancerOnset_Gompertz_Rate_P = 0.08, #change for approach with exponential
    CancerOnset_Gompertz_Rate_D = 0.04, #0.1 change for approach with exponential
    CancerOnset_Gompertz_Rate_R = 0.18, #change for approach with exponential
    pSxDetS1_P = 0.192,
    pSxDetS1_D = 0.192,
    pSxDetS1_R = 0.192,
    hr_SxDetS2S1_P = 24,
    hr_SxDetS2S1_D = 24,
    hr_SxDetS2S1_R = 24,
    hr_SxDetS3S2_P = 7,
    hr_SxDetS3S2_D = 7,
    hr_SxDetS3S2_R = 7,
    hr_SxDetS4S3_P = 5,
    hr_SxDetS4S3_D = 5,
    hr_SxDetS4S3_R = 5
  )
)


# Modify a single parameter's bounds
l_params_priors_adenoma_Chile$lb["alpha_lesion_adenoma"] <- -10
l_params_priors_adenoma_Chile$ub["alpha_lesion_adenoma"] <- -5


l_params_priors_adenoma_Chile$lb["hr_SxDetS2S1_P"] <- 1.2
l_params_priors_adenoma_Chile$lb["hr_SxDetS2S1_D"] <- 1.2
l_params_priors_adenoma_Chile$lb["hr_SxDetS2S1_R"] <- 1.2
l_params_priors_adenoma_Chile$lb["hr_SxDetS3S2_P"] <- 4
l_params_priors_adenoma_Chile$lb["hr_SxDetS3S2_D"] <- 4
l_params_priors_adenoma_Chile$lb["hr_SxDetS3S2_R"] <- 4
l_params_priors_adenoma_Chile$lb["hr_SxDetS4S3_P"] <- 8
l_params_priors_adenoma_Chile$lb["hr_SxDetS4S3_D"] <- 8
l_params_priors_adenoma_Chile$lb["hr_SxDetS4S3_R"] <- 8

l_params_priors_adenoma_Chile$ub["hr_SxDetS2S1_P"] <- 9
l_params_priors_adenoma_Chile$ub["hr_SxDetS2S1_D"] <- 9
l_params_priors_adenoma_Chile$ub["hr_SxDetS2S1_R"] <- 9
l_params_priors_adenoma_Chile$ub["hr_SxDetS3S2_P"] <- 10.5
l_params_priors_adenoma_Chile$ub["hr_SxDetS3S2_D"] <- 10.5
l_params_priors_adenoma_Chile$ub["hr_SxDetS3S2_R"] <- 10.5
l_params_priors_adenoma_Chile$ub["hr_SxDetS4S3_P"] <- 11
l_params_priors_adenoma_Chile$ub["hr_SxDetS4S3_D"] <- 11
l_params_priors_adenoma_Chile$ub["hr_SxDetS4S3_R"] <- 11


l_params_priors_adenoma_Chile$lb["pSxDetS1_P"]  <-  0.0001
l_params_priors_adenoma_Chile$lb["pSxDetS1_D"]  <-  0.0001
l_params_priors_adenoma_Chile$lb["pSxDetS1_R"]  <-  0.0001


l_params_priors_adenoma_Chile$ub["pSxDetS1_P"]  <-  0.192*.2
l_params_priors_adenoma_Chile$ub["pSxDetS1_D"]  <-  0.192*.2
l_params_priors_adenoma_Chile$ub["pSxDetS1_R"]  <-  0.192*.2


l_params_priors_adenoma_Chile$ub["beta_age"]  <-  0.08



#Claudia's suggestion to increase the upper bound for the variance of the individual risk multiplier to allow for more variability in the calibration process

l_params_priors_adenoma_Chile$ub["IndividuakRiskMultVariance"]  <-  1.2

l_params_priors_adenoma_Chile$lb["AdNaturalHistoryPropensity_Gaussian_Variance"]  <-  .2
l_params_priors_adenoma_Chile$ub["AdNaturalHistoryPropensity_Gaussian_Variance"]  <-  .4


#Posteriors from the US calibration suggest that the alpha_lesion_adenoma parameter could be in a narrower range, so we can adjust the bounds accordingly to focus the calibration on that range.
l_params_priors_adenoma_Chile$ub["alpha_lesion_adenoma"]  <-  -6.5
l_params_priors_adenoma_Chile$lb["alpha_lesion_adenoma"]  <-  -7.6

# Similarly, if the posteriors suggest that the beta_age parameter is likely in a narrower range, we can adjust its bounds as well.
l_params_priors_adenoma_Chile$ub["beta_age"]  <-  0.052
l_params_priors_adenoma_Chile$lb["beta_age"]  <-  0.039



l_params_priors_adenoma_Chile$ub["CancerOnset_Gompertz_Rate_P"]  <-  0.04
l_params_priors_adenoma_Chile$ub["CancerOnset_Gompertz_Rate_R"]  <- 0.08
l_params_priors_adenoma_Chile$ub["CancerOnset_Gompertz_Rate_D"]  <- 0.02




#========================================================================#
# 1. Calibration targets ----
#========================================================================#

l_selected_simcrc_targets_both <- c(
  "Prev0AdOrPreClin_age27",#
  "Prev0AdOrPreClin_age32",
  "Prev0AdOrPreClin_age37",#new
  "Prev0AdOrPreClin_age42",#new
  "Prev0AdOrPreClin_age52",
  "Prev0AdOrPreClin_age62",
  "Prev0AdOrPreClin_age67",
  "Prev0AdOrPreClin_age77",
  "Prev1AdOrPreClin_age27",#new
  "Prev1AdOrPreClin_age32",#new
  "Prev1AdOrPreClin_age37",#new
  "Prev1AdOrPreClin_age42",#new
  "Prev1AdOrPreClin_age52",
  "Prev1AdOrPreClin_age62",
  "Prev1AdOrPreClin_age67",
  "Prev1AdOrPreClin_age77",
  "Prev2AdOrPreClin_age27",#new
  "Prev2AdOrPreClin_age32",#new
  "Prev2AdOrPreClin_age37",#new
  "Prev2AdOrPreClin_age42",
  "Prev2AdOrPreClin_age52",
  "Prev2AdOrPreClin_age62",
  "Prev2AdOrPreClin_age67",
  "Prev2AdOrPreClin_age77",
  "Prev3AdOrPreClin_age27",#new
  "Prev3AdOrPreClin_age32",#new
  "Prev3AdOrPreClin_age37",#new   
  "Prev3AdOrPreClin_age42",#new   
  "Prev3AdOrPreClin_age52",
  "Prev3AdOrPreClin_age62",
  "Prev3AdOrPreClin_age67",
  "Prev3AdOrPreClin_age77",
  "SizeLRGivenAdInP_wtd",
  "SizeMRGivenAdInP_wtd",
  "SizeHRGivenAdInP_wtd",
  "SizeLRGivenAdInD_wtd",
  "SizeMRGivenAdInD_wtd",
  "SizeHRGivenAdInD_wtd",
  "SizeLRGivenAdInR_wtd",
  "SizeMRGivenAdInR_wtd",
  "SizeHRGivenAdInR_wtd",
  "PrevPreclinical_ages50_59",
  "PrevPreclinical_ages60_69",
  "PrevPreclinical_ages70_79",
  "PrevPreclinical_ages80_89",
  "CRCincPer100K_P_ages50_59",
  "CRCincPer100K_P_ages60_69",
  "CRCincPer100K_P_ages70_79",
  "CRCincPer100K_P_ages80_99",
  "CRCincPer100K_D_ages50_59",
  "CRCincPer100K_D_ages60_69",
  "CRCincPer100K_D_ages70_79",
  "CRCincPer100K_D_ages80_99",
  "CRCincPer100K_R_ages50_59",
  "CRCincPer100K_R_ages60_69",
  "CRCincPer100K_R_ages70_79",
  "CRCincPer100K_R_ages80_99",
  "CRCincPer100K_ages50_59",
  "CRCincPer100K_ages60_69",
  "CRCincPer100K_ages70_79",
  "CRCincPer100K_ages80_99",
  "StageDistribution_P_S1",
  "StageDistribution_P_S2",
  "StageDistribution_P_S3",
  "StageDistribution_P_S4",
  "StageDistribution_D_S1",
  "StageDistribution_D_S2",
  "StageDistribution_D_S3",
  "StageDistribution_D_S4",
  "StageDistribution_R_S1",
  "StageDistribution_R_S2",
  "StageDistribution_R_S3",
  "StageDistribution_R_S4",
  "prop_SSP",
  "Prev0SSPOrPreClin_age32",
  "Prev0SSPOrPreClin_age37",
  "Prev0SSPOrPreClin_age42",
  "Prev0SSPOrPreClin_age47",
  "Prev0SSPOrPreClin_age52",
  "Prev0SSPOrPreClin_age57",
  "Prev0SSPOrPreClin_age62",
  "Prev0SSPOrPreClin_age67",
  "Prev0SSPOrPreClin_age72",
  "Prev0SSPOrPreClin_age77",
  "Prev1SSPOrPreClin_age32",
  "Prev1SSPOrPreClin_age37",
  "Prev1SSPOrPreClin_age42",
  "Prev1SSPOrPreClin_age47",
  "Prev1SSPOrPreClin_age52",
  "Prev1SSPOrPreClin_age57",
  "Prev1SSPOrPreClin_age62",
  "Prev1SSPOrPreClin_age67",
  "Prev1SSPOrPreClin_age72",
  "Prev1SSPOrPreClin_age77",
  "Prev2SSPOrPreClin_age32",
  "Prev2SSPOrPreClin_age37",
  "Prev2SSPOrPreClin_age42",
  "Prev2SSPOrPreClin_age47",
  "Prev2SSPOrPreClin_age52",
  "Prev2SSPOrPreClin_age57",
  "Prev2SSPOrPreClin_age62",
  "Prev2SSPOrPreClin_age67",
  "Prev2SSPOrPreClin_age72",
  "Prev2SSPOrPreClin_age77",
  "Prev3SSPOrPreClin_age32",
  "Prev3SSPOrPreClin_age37",
  "Prev3SSPOrPreClin_age42",
  "Prev3SSPOrPreClin_age47",
  "Prev3SSPOrPreClin_age52",
  "Prev3SSPOrPreClin_age57",
  "Prev3SSPOrPreClin_age62",
  "Prev3SSPOrPreClin_age67",
  "Prev3SSPOrPreClin_age72",
  "Prev3SSPOrPreClin_age77",
  "SizeLRGivenSSPInP_wtd",
  "SizeMRGivenSSPInP_wtd",
  "SizeHRGivenSSPInP_wtd",
  "SizeLRGivenSSPInD_wtd",
  "SizeMRGivenSSPInD_wtd",
  "SizeHRGivenSSPInD_wtd",
  "SizeLRGivenSSPInR_wtd",
  "SizeMRGivenSSPInR_wtd",
  "SizeHRGivenSSPInR_wtd" )

# #Selected targets
l_selected_simcrc_targets_adenoma <- c("Prev0AdOrPreClin_age27",
                                     "Prev0AdOrPreClin_age32",
                                     "Prev0AdOrPreClin_age37",
                                     "Prev0AdOrPreClin_age42",
                                     "Prev0AdOrPreClin_age47",
                                     "Prev0AdOrPreClin_age52",
                                     "Prev0AdOrPreClin_age57",
                                     "Prev0AdOrPreClin_age62",
                                     "Prev0AdOrPreClin_age67",
                                     "Prev0AdOrPreClin_age72",
                                     "Prev0AdOrPreClin_age77",
                                     "Prev0AdOrPreClin_age82",
                                     "Prev0AdOrPreClin_age87",
                                     "Prev0AdOrPreClin_age92",
                                     "Prev0AdOrPreClin_age97",
                                     "Prev1AdOrPreClin_age27",   
                                     "Prev1AdOrPreClin_age32",
                                     "Prev1AdOrPreClin_age37",
                                     "Prev1AdOrPreClin_age42",
                                     "Prev1AdOrPreClin_age47",   
                                     "Prev1AdOrPreClin_age52",
                                     "Prev1AdOrPreClin_age57",
                                     "Prev1AdOrPreClin_age62",
                                     "Prev1AdOrPreClin_age67",   
                                     "Prev1AdOrPreClin_age72",    
                                     "Prev1AdOrPreClin_age77",   
                                     "Prev1AdOrPreClin_age82",    
                                     "Prev1AdOrPreClin_age87",   
                                     "Prev1AdOrPreClin_age92",    
                                     "Prev1AdOrPreClin_age97",    
                                     "Prev2AdOrPreClin_age27",    
                                     "Prev2AdOrPreClin_age32",   
                                     "Prev2AdOrPreClin_age37",    
                                     "Prev2AdOrPreClin_age42",    
                                     "Prev2AdOrPreClin_age47",    
                                     "Prev2AdOrPreClin_age52",   
                                     "Prev2AdOrPreClin_age57",    
                                     "Prev2AdOrPreClin_age62",    
                                     "Prev2AdOrPreClin_age67",    
                                     "Prev2AdOrPreClin_age72",   
                                     "Prev2AdOrPreClin_age77",    
                                     "Prev2AdOrPreClin_age82",
                                     "Prev2AdOrPreClin_age87",
                                     "Prev2AdOrPreClin_age92",
                                     "Prev2AdOrPreClin_age97",
                                     "Prev3AdOrPreClin_age27",    
                                     "Prev3AdOrPreClin_age32",    
                                     "Prev3AdOrPreClin_age37",   
                                     "Prev3AdOrPreClin_age42",   
                                     "Prev3AdOrPreClin_age47",   
                                     "Prev3AdOrPreClin_age52",   
                                     "Prev3AdOrPreClin_age57",   
                                     "Prev3AdOrPreClin_age62",    
                                     "Prev3AdOrPreClin_age67",    
                                     "Prev3AdOrPreClin_age72",    
                                     "Prev3AdOrPreClin_age77",   
                                     "Prev3AdOrPreClin_age82",
                                     "Prev3AdOrPreClin_age87",
                                     "Prev3AdOrPreClin_age92",
                                     "Prev3AdOrPreClin_age97",
                                     "SizeLRGivenAdInP_wtd",         
                                     "SizeMRGivenAdInP_wtd",         
                                     "SizeHRGivenAdInP_wtd",         
                                     "SizeLRGivenAdInD_wtd",         
                                     "SizeMRGivenAdInD_wtd",         
                                     "SizeHRGivenAdInD_wtd",         
                                     "SizeLRGivenAdInR_wtd",         
                                     "SizeMRGivenAdInR_wtd",         
                                     "SizeHRGivenAdInR_wtd", 
                                     "PrevPreclinical_ages40_49",
                                     "PrevPreclinical_ages50_59",
                                     "PrevPreclinical_ages60_69",
                                     "PrevPreclinical_ages70_79",
                                     "PrevPreclinical_ages80_89",
                                     "CRCincPer100K_P_ages20_39", 
                                     "CRCincPer100K_P_ages40_49", 
                                     "CRCincPer100K_P_ages50_59", 
                                     "CRCincPer100K_P_ages60_69",
                                     "CRCincPer100K_P_ages70_79", 
                                     "CRCincPer100K_P_ages80_99", 
                                     "CRCincPer100K_D_ages20_39", 
                                     "CRCincPer100K_D_ages40_49",
                                     "CRCincPer100K_D_ages50_59", 
                                     "CRCincPer100K_D_ages60_69", 
                                     "CRCincPer100K_D_ages70_79", 
                                     "CRCincPer100K_D_ages80_99",
                                     "CRCincPer100K_R_ages20_39", 
                                     "CRCincPer100K_R_ages40_49", 
                                     "CRCincPer100K_R_ages50_59", 
                                     "CRCincPer100K_R_ages60_69",
                                     "CRCincPer100K_R_ages70_79", 
                                     "CRCincPer100K_R_ages80_99", 
                                     "CRCincPer100K_ages20_39",   
                                     "CRCincPer100K_ages40_49",  
                                     "CRCincPer100K_ages50_59",   
                                     "CRCincPer100K_ages60_69",   
                                     "CRCincPer100K_ages70_79",   
                                     "CRCincPer100K_ages80_99",  
                                     "StageDistribution_P_S1" ,   
                                     "StageDistribution_P_S2" ,   
                                     "StageDistribution_P_S3" ,   
                                     "StageDistribution_P_S4" ,  
                                     "StageDistribution_D_S1" ,   
                                     "StageDistribution_D_S2" ,   
                                     "StageDistribution_D_S3" ,   
                                     "StageDistribution_D_S4" ,  
                                     "StageDistribution_R_S1" ,   
                                     "StageDistribution_R_S2" ,   
                                     "StageDistribution_R_S3" ,   
                                     "StageDistribution_R_S4" )


#Selected targets
l_selected_simcrc_targets_mexico <- c(
  
  # Mexico Targets
  "PrevAnyAdOrPreClin_age45",
  "PrevAnyAdOrPreClin_age55",
  "PrevAnyAdOrPreClin_age65",
  "PrevAnyAdOrPreClin_age75",
  "PrevAnyAdOrPreClin_age85",
  "CRCincPer100K_ages20_24",
  "CRCincPer100K_ages25_29",
  "CRCincPer100K_ages30_34",
  "CRCincPer100K_ages35_39",
  "CRCincPer100K_ages40_44",
  "CRCincPer100K_ages45_49",
  "CRCincPer100K_ages50_54",
  "CRCincPer100K_ages55_59",
  "CRCincPer100K_ages60_64",
  "CRCincPer100K_ages65_69",
  "CRCincPer100K_ages70_74",
  "CRCincPer100K_ages75_79",
  "CRCincPer100K_ages80_84",
  "CRCincPer100K_ages85_89",
  "CRCincPer100K_ages90_94",
  "StageDistribution_S1",
  "StageDistribution_S2",
  "StageDistribution_S3",
  "StageDistribution_S4")



# #Selected targets
l_selected_simcrc_targets_Chile <- c("Prev0AdOrPreClin_age27",
                                       "Prev0AdOrPreClin_age32",
                                       "Prev0AdOrPreClin_age37",
                                       "Prev0AdOrPreClin_age42",
                                       "Prev0AdOrPreClin_age47",
                                       "Prev0AdOrPreClin_age52",
                                       "Prev0AdOrPreClin_age57",
                                       "Prev0AdOrPreClin_age62",
                                       "Prev0AdOrPreClin_age67",
                                       "Prev0AdOrPreClin_age72",
                                       "Prev0AdOrPreClin_age77",
                                       "Prev0AdOrPreClin_age82",
                                       "Prev0AdOrPreClin_age87",
                                       "Prev0AdOrPreClin_age92",
                                       "Prev0AdOrPreClin_age97",
                                       "Prev1AdOrPreClin_age27",   
                                       "Prev1AdOrPreClin_age32",
                                       "Prev1AdOrPreClin_age37",
                                       "Prev1AdOrPreClin_age42",
                                       "Prev1AdOrPreClin_age47",   
                                       "Prev1AdOrPreClin_age52",
                                       "Prev1AdOrPreClin_age57",
                                       "Prev1AdOrPreClin_age62",
                                       "Prev1AdOrPreClin_age67",   
                                       "Prev1AdOrPreClin_age72",    
                                       "Prev1AdOrPreClin_age77",   
                                       "Prev1AdOrPreClin_age82",    
                                       "Prev1AdOrPreClin_age87",   
                                       "Prev1AdOrPreClin_age92",    
                                       "Prev1AdOrPreClin_age97",    
                                       "Prev2AdOrPreClin_age27",    
                                       "Prev2AdOrPreClin_age32",   
                                       "Prev2AdOrPreClin_age37",    
                                       "Prev2AdOrPreClin_age42",    
                                       "Prev2AdOrPreClin_age47",    
                                       "Prev2AdOrPreClin_age52",   
                                       "Prev2AdOrPreClin_age57",    
                                       "Prev2AdOrPreClin_age62",    
                                       "Prev2AdOrPreClin_age67",    
                                       "Prev2AdOrPreClin_age72",   
                                       "Prev2AdOrPreClin_age77",    
                                       "Prev2AdOrPreClin_age82",
                                       "Prev2AdOrPreClin_age87",
                                       "Prev2AdOrPreClin_age92",
                                       "Prev2AdOrPreClin_age97",
                                       "Prev3AdOrPreClin_age27",    
                                       "Prev3AdOrPreClin_age32",    
                                       "Prev3AdOrPreClin_age37",   
                                       "Prev3AdOrPreClin_age42",   
                                       "Prev3AdOrPreClin_age47",   
                                       "Prev3AdOrPreClin_age52",   
                                       "Prev3AdOrPreClin_age57",   
                                       "Prev3AdOrPreClin_age62",    
                                       "Prev3AdOrPreClin_age67",    
                                       "Prev3AdOrPreClin_age72",    
                                       "Prev3AdOrPreClin_age77",   
                                       "Prev3AdOrPreClin_age82",
                                       "Prev3AdOrPreClin_age87",
                                       "Prev3AdOrPreClin_age92",
                                       "Prev3AdOrPreClin_age97",
                                       "SizeLRGivenAdInP_wtd",         
                                       "SizeMRGivenAdInP_wtd",         
                                       "SizeHRGivenAdInP_wtd",         
                                       "SizeLRGivenAdInD_wtd",         
                                       "SizeMRGivenAdInD_wtd",         
                                       "SizeHRGivenAdInD_wtd",         
                                       "SizeLRGivenAdInR_wtd",         
                                       "SizeMRGivenAdInR_wtd",         
                                       "SizeHRGivenAdInR_wtd", 
                                       "PrevPreclinical_ages40_49",
                                       "PrevPreclinical_ages50_59",
                                       "PrevPreclinical_ages60_69",
                                       "PrevPreclinical_ages70_79",
                                       "PrevPreclinical_ages80_89",
                                       "CRCincPer100K_P_ages20_39", 
                                       "CRCincPer100K_P_ages40_49", 
                                       "CRCincPer100K_P_ages50_59", 
                                       "CRCincPer100K_P_ages60_69",
                                       "CRCincPer100K_P_ages70_79", 
                                       "CRCincPer100K_P_ages80_99", 
                                       "CRCincPer100K_D_ages20_39", 
                                       "CRCincPer100K_D_ages40_49",
                                       "CRCincPer100K_D_ages50_59", 
                                       "CRCincPer100K_D_ages60_69", 
                                       "CRCincPer100K_D_ages70_79", 
                                       "CRCincPer100K_D_ages80_99",
                                       "CRCincPer100K_R_ages20_39", 
                                       "CRCincPer100K_R_ages40_49", 
                                       "CRCincPer100K_R_ages50_59", 
                                       "CRCincPer100K_R_ages60_69",
                                       "CRCincPer100K_R_ages70_79", 
                                       "CRCincPer100K_R_ages80_99", 
                                       "CRCincPer100K_ages20_39",   
                                       "CRCincPer100K_ages40_49",  
                                       "CRCincPer100K_ages50_59",   
                                       "CRCincPer100K_ages60_69",   
                                       "CRCincPer100K_ages70_79",   
                                       "CRCincPer100K_ages80_99",  
                                       "StageDistribution_P_S1" ,   
                                       "StageDistribution_P_S2" ,   
                                       "StageDistribution_P_S3" ,   
                                       "StageDistribution_P_S4" ,  
                                       "StageDistribution_D_S1" ,   
                                       "StageDistribution_D_S2" ,   
                                       "StageDistribution_D_S3" ,   
                                       "StageDistribution_D_S4" ,  
                                       "StageDistribution_R_S1" ,   
                                       "StageDistribution_R_S2" ,   
                                       "StageDistribution_R_S3" ,   
                                       "StageDistribution_R_S4" )





l_calibration_targets_files <- list(
  "Both_female" = "data-raw/20250626_simcrc_targets_ssp_unadjusted.csv",
  "Both_male"  = "data-raw/20250626_simcrc_targets_ssp_unadjusted.csv", 
  "Adenoma_female" = "data-raw/20220909_simcrc_targets_F.csv",
  "Adenoma_male" = "data-raw/20250904_simcrc_targets_M.csv",
  "Mexico" = "data-raw/20251016_simcrc_mx_targets.csv",
  "Chile" = "data-raw/true_target_simcrcRvCH.csv"
)


#========================================================================#
# 2. Calibration parameters ----
#========================================================================#

#ADD LIST WITH LB AND UB FOR EACH MODEL

l_model_female_adenoma <- list(
  lesion_type        = "Adenoma",
  population_type    = "F" ,
  params_priors      = l_params_priors_female_adenoma,
  selected_targets   = l_selected_simcrc_targets_adenoma,
  targest_file       = "data-raw/20220909_simcrc_targets_F.csv",
  SSP_pathway        = FALSE,
  life_table_fem     = get_life_table(sex="female", year = 1977),
  life_table_male    = NULL,
  cohort_year        = 1980,
  birth_year         = 1980,
  p_female           = 1,
  p_white            = 0.8,
  n_pop              = 1e5,
  l_params_init     = simcrc::load_params_init(fromFile = TRUE, filename = simcrc::l_calibrated_params$female$Min_AbsolutErr),
  min_age_lesion_onset = 10
  
)

l_model_male_adenoma <- list(
  lesion_type        = "Adenoma",
  population_type    = "M" ,
  params_priors      = l_params_priors_male_adenoma,
  selected_targets   = l_selected_simcrc_targets_adenoma,
  targest_file       = "data-raw/20250904_simcrc_targets_M.csv",
  SSP_pathway        = FALSE,
  life_table_fem     = NULL,
  life_table_male    = get_life_table(sex="male", year = 1977),
  cohort_year        = 1980,
  birth_year         = 1980,
  p_female           = 0,
  p_white            = 0.8,
  n_pop              = 1e5,
  l_params_init      = simcrc::load_params_init(fromFile = TRUE, filename = simcrc::l_calibrated_params$female$Min_AbsolutErr),
  min_age_lesion_onset = 10
)

l_model_female_both <- list(
  lesion_type        = "Both",
  population_type    = "F" ,
  params_priors      = l_params_priors_female_both,
  selected_targets   = l_selected_simcrc_targets_both,
  targest_file       = "data-raw/20250626_simcrc_targets_ssp_unadjusted.csv",
  SSP_pathway        = TRUE,
  life_table_fem     = get_life_table(sex="female", year = 1977),
  life_table_male    = NULL,
  cohort_year        = 1980,
  birth_year         = 1980,
  p_female           = 1,
  p_white            = 0.8,
  n_pop              = 1e5,
  l_params_init     = simcrc::load_params_init(fromFile = TRUE, filename = simcrc::l_calibrated_params$female$Min_AbsolutErr),
  min_age_lesion_onset = 10
  
)

l_model_male_both <- list(
  lesion_type        = "Both",
  population_type    = "M" ,
  params_priors      = l_params_priors_female_both,
  selected_targets   = l_selected_simcrc_targets_both,
  targest_file       = "data-raw/20250626_simcrc_targets_ssp_unadjusted.csv",
  SSP_pathway        = TRUE,
  life_table_fem     = NULL,
  life_table_male    = get_life_table(sex="male", year = 1977),
  cohort_year        = 1980,
  birth_year         = 1980,
  p_female           = 0,
  p_white            = 0.8,
  n_pop              = 1e5,
  l_params_init     = simcrc::load_params_init(fromFile = TRUE, filename = simcrc::l_calibrated_params$female$Min_AbsolutErr),
  min_age_lesion_onset = 10
)



df_life_table_2017_CH <- readRDS("~/Documents/GitHub/simcrc_chile/data-raw/df_life_table_2017_CH.rds")

#read csv

df_life_table_2017_CH <- read.csv("~/Documents/GitHub/simcrc_chile/data-raw/df_lifetable_2017_CH.csv")


#rename age as Age

colnames(df_life_table_2017_CH)[colnames(df_life_table_2017_CH) == "age"] <- "Age"
colnames(df_life_table_2017_CH)[colnames(df_life_table_2017_CH) == "mortality_rate"] <- "mortality.rates"

l_model_adenoma_Chile <- list(
  lesion_type        = "Adenoma",
  population_type    = "F" ,
  params_priors      = l_params_priors_adenoma_Chile,
  selected_targets   = l_selected_simcrc_targets_Chile,
  targest_file       = "data-raw/true_target_simcrcRvCH.csv",
  SSP_pathway        = FALSE,
  life_table_fem     = df_life_table_2017_CH,
  life_table_male   = NULL,
  cohort_year        = 2017,
  birth_year         = 2017,
  p_female           = 1,
  p_white            = 0.8,
  n_pop              = 1e5,
  l_params_init     = simcrc::load_params_init(fromFile = TRUE, filename = simcrc::l_calibrated_params$female$Min_AbsolutErr),
  min_age_lesion_onset = 10)
  


