# *****************************************************************************
#
# Script: wrangle_fonasa_population.R
#
# Purpose: Summarize the national FONASA-beneficiary export (Nicolas, Tableau
#          public dashboard, Sep 2026) into a national population-by-age-group
#          table for the colonoscopy demand model (D = sum(P_a * k_a)).
#
#          GAP: FONASA only publishes 5-year age bands (EDAD_TRAMO), not
#          single-year-of-age counts. The manuscript's P_a term calls for
#          single year of age, ages 40-100 (Chile CRC Screening CEA - JR
#          tracked v6.docx, demand section). This script's output is an
#          interim 5-year-band proxy; follow up with Nicolas to request a
#          FONASA breakdown by single year of age before finalizing P_a.
#
# Author: Jorge Roa
#
# Email: jorgeroa@stanford.edu
#
# Date Created: 24 Sep 2026
#
# *****************************************************************************
#
# Notes:
#   Source: data-raw/fonasa_beneficiarios_2025.csv (raw export, one row per
#   region x comuna x sexo x edad_tramo x nacionalidad x titular_carga x
#   tramo). TRAMO (A-D) is FONASA's income-tier code, unrelated to age.
#
# *****************************************************************************

library(readr)
library(dplyr)

path_fonasa_raw <- "data-raw/fonasa_beneficiarios_2025.csv"
path_fonasa_out <- "data-raw/fonasa_population_by_age_group_2025.csv"

v_age_group_levels <- c(
  "00 a 02 años", "03 a 04 años", "05 a 09 años", "10 a 14 años",
  "15 a 19 años", "20 a 24 años", "25 a 29 años", "30 a 34 años",
  "35 a 39 años", "40 a 44 años", "45 a 49 años", "50 a 54 años",
  "55 a 59 años", "60 a 64 años", "65 a 69 años", "70 a 74 años",
  "75 a 79 años", "80 a 84 años", "85 a 89 años", "90 a 94 años",
  "95 a 99 años", "Más de 99 años", "S.I."
)

df_fonasa_raw <- read_csv(path_fonasa_raw, show_col_types = FALSE)

df_population_by_age <- df_fonasa_raw |>
  group_by(grupo_edad = factor(EDAD_TRAMO, levels = v_age_group_levels)) |>
  summarise(n_beneficiarios = sum(BENEFICIARIOS), .groups = "drop") |>
  arrange(grupo_edad)

write_csv(df_population_by_age, path_fonasa_out)
