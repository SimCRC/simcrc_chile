# *****************************************************************************
#
# Script: build_final_deck_tables.R
#
# Purpose: Markdown tables for the final Chile deck, all computed from files in
#          this repo: the strategy list, the 68-strategy ICER tables, the model
#          run data, the cost inputs and the FONASA population by age group.
#          Nothing on a slide is typed by hand when it can be read from here.
#
# Author: Jorge Roa
#
# Email: jorgeroa@stanford.edu
#
# ISS: #25   Version: #24
#
# Date Created: 05 October 2026
#
# *****************************************************************************
#
# Notes:
#   Run from the simcrc_chile root:
#     Rscript analysis/build_final_deck_tables.R [dir_tab]
#   Run: calibration v0.15.0.20260921.1018, Min_AbsolutErr, 68 strategies,
#   1,000,000 followed from age 40, willingness to pay 16M CLP per QALY.
#
# *****************************************************************************

v_args  <- commandArgs(trailingOnly = TRUE)
dir_tab <- if (length(v_args) >= 1) v_args[1] else "SimCRC slides/final_deck/tables"
dir.create(dir_tab, recursive = TRUE, showWarnings = FALSE)

n_wtp    <- 16e6
path_run <- file.path("output/2026Chile/Base/ce_results",
                      "2026-09-21 230512.351928_model_run_data_Base.xlsx")

df_all <- read.csv("ce_results/df_icer_all_strategies.csv")
df_fit <- read.csv("ce_results/df_icer_FIT.csv")
df_str <- read.csv("data-raw/df_CH_2026_strategies.csv")
df_run <- readxl::read_excel(path_run, sheet = "model_run_data")

df_all$label <- ifelse(df_all$Strategy == "NoScreening", "No screening",
  ifelse(grepl("_100$", df_all$Strategy),
         sub("^([A-Z]+)_([0-9]+)_.*$", "\\1 \\2 once", df_all$Strategy),
         sub("^([A-Z]+)_([0-9]+)_([0-9]+)_([0-9]+)$", "\\1 \\2-\\3 q\\4",
             df_all$Strategy)))
df_fit$label <- df_all$label[match(df_fit$Strategy, df_all$Strategy)]
df_run$label <- df_all$label[match(df_run$Strategy, df_all$Strategy)]

stopifnot(nrow(df_all) == 68, nrow(df_str) == 78, !anyNA(df_fit$label))

df_str$group <- ifelse(df_str$modality == "NoScreening", "No screening",
  ifelse(df_str$modality == "FIT", "FIT",
         ifelse(df_str$frequency_screening == 100, "Colonoscopy, once",
                "Colonoscopy, repeated")))
df_all$group <- ifelse(df_all$Strategy == "NoScreening", "No screening",
  ifelse(grepl("^FIT", df_all$Strategy), "FIT",
         ifelse(grepl("_100$", df_all$Strategy), "Colonoscopy, once",
                "Colonoscopy, repeated")))

v_grp <- c("No screening", "Colonoscopy, repeated", "Colonoscopy, once", "FIT")
l_row <- vector("list", length(v_grp))
for (i in seq_along(v_grp)) {
  df_g <- df_str[df_str$group == v_grp[i], ]
  v_int <- sort(unique(df_g$frequency_screening))
  l_row[[i]] <- data.frame(
    group = v_grp[i],
    interval = if (v_grp[i] == "No screening") "-"
               else if (v_grp[i] == "Colonoscopy, once") "once"
               else paste(v_int, collapse = ", "),
    begin = if (v_grp[i] == "No screening") "-"
            else paste(sort(unique(df_g$age_to_begin_screening)), collapse = ", "),
    stop = if (v_grp[i] %in% c("No screening", "Colonoscopy, once")) "-"
           else paste(sort(unique(df_g$age_to_end_screening)), collapse = ", "),
    n_raw = nrow(df_g),
    n_unique = sum(df_all$group == v_grp[i])
  )
}
df_tab1 <- do.call(rbind, l_row)
stopifnot(sum(df_tab1$n_raw) == 78, sum(df_tab1$n_unique) == 68)
writeLines(c(
  "| Strategy group | Interval, years | Begin at | Stop at | Strategies | Unique |",
  "| :--- | :--- | :--- | :--- | ---: | ---: |",
  sprintf("| %s | %s | %s | %s | %d | %d |", df_tab1$group, df_tab1$interval,
          df_tab1$begin, df_tab1$stop, df_tab1$n_raw, df_tab1$n_unique),
  sprintf("| **Total** | | | | **%d** | **%d** |", sum(df_tab1$n_raw),
          sum(df_tab1$n_unique))
), file.path(dir_tab, "strategies.md"))

df_cand_all <- df_all[df_all$Status == "ND" & !is.na(df_all$ICER) &
                        df_all$ICER <= n_wtp, ]
df_cand_fit <- df_fit[df_fit$Status == "ND" & !is.na(df_fit$ICER) &
                        df_fit$ICER <= n_wtp, ]
df_best <- rbind(
  cbind(set = "All 68 strategies", df_cand_all[which.max(df_cand_all$Effect),
        c("label", "Cost", "Effect", "ICER")]),
  cbind(set = "FIT strategies only", df_cand_fit[which.max(df_cand_fit$Effect),
        c("label", "Cost", "Effect", "ICER")])
)
writeLines(c(
  "| Strategy set | Choice at 16M per QALY | Cost per 1,000 | QALYs gained per 1,000 | ICER per QALY |",
  "| :--- | :--- | ---: | ---: | ---: |",
  sprintf("| %s | **%s** | %.0fM | %.1f | %.1fM |", df_best$set, df_best$label,
          df_best$Cost / 1e6, df_best$Effect, df_best$ICER / 1e6)
), file.path(dir_tab, "best_choice.md"))

df_nd_all <- df_all[df_all$Status == "ND", ]
df_nd_fit <- df_fit[df_fit$Status == "ND", ]
writeLines(c(
  "| Strategy | Cost per 1,000 | QALYs gained | ICER per QALY |",
  "| :---------------------- | -----: | -----: | --------: |",
  sprintf("| %s%s%s | %.0fM | %.1f | %s |",
          ifelse(df_nd_all$label == df_best$label[1], "**", ""),
          df_nd_all$label,
          ifelse(df_nd_all$label == df_best$label[1], "**", ""),
          df_nd_all$Cost / 1e6, df_nd_all$Effect,
          ifelse(is.na(df_nd_all$ICER), "reference",
                 sprintf("%.1fM", df_nd_all$ICER / 1e6)))
), file.path(dir_tab, "frontier_all.md"))
writeLines(c(
  "| Strategy | Cost per 1,000 | QALYs gained | ICER per QALY |",
  "| :---------------------- | -----: | -----: | --------: |",
  sprintf("| %s%s%s | %.0fM | %.1f | %s |",
          ifelse(df_nd_fit$label == df_best$label[2], "**", ""),
          df_nd_fit$label,
          ifelse(df_nd_fit$label == df_best$label[2], "**", ""),
          df_nd_fit$Cost / 1e6, df_nd_fit$Effect,
          ifelse(is.na(df_nd_fit$ICER), "reference",
                 sprintf("%.1fM", df_nd_fit$ICER / 1e6)))
), file.path(dir_tab, "frontier_fit.md"))

df_bur <- df_run[match(df_nd_all$Strategy, df_run$Strategy), ]
writeLines(c(
  "| Strategy | CRC cases | Cases averted | Deaths averted | Mortality reduction | Screening COL | Follow-up COL | Surveillance COL |",
  "| :------------- | -----: | -----: | -----: | -----: | -----: | -----: | -----: |",
  sprintf("| %s | %.1f | %.1f | %.1f | %.1f%% | %s | %s | %s |", df_bur$label,
          df_bur$CRCCasesper1000, df_bur$CRCCasesAvertedper1000,
          df_bur$CRCDeathsAvertedper1000, 100 * df_bur$MortalityReduction,
          format(round(df_bur$ScreeningColsper1000), big.mark = ",", trim = TRUE),
          format(round(df_bur$DiagnosticFUper1000), big.mark = ",", trim = TRUE),
          format(round(df_bur$SurveillanceColsper1000), big.mark = ",", trim = TRUE))
), file.path(dir_tab, "burden.md"))

writeLines(c(
  "| Strategy | Cost per 1,000 | QALYs gained | ICER per QALY | CRC cases | Cases averted | Deaths averted | Mortality reduction | Screening COL | Follow-up COL | Surveillance COL |",
  "| :------------- | -----: | -----: | -----: | -----: | -----: | -----: | -----: | -----: | -----: | -----: |",
  sprintf("| %s%s%s | %.0fM | %.1f | %s | %.1f | %.1f | %.1f | %.1f%% | %s | %s | %s |",
          ifelse(df_nd_all$label == df_best$label[1], "**", ""),
          df_nd_all$label,
          ifelse(df_nd_all$label == df_best$label[1], "**", ""),
          df_nd_all$Cost / 1e6, df_nd_all$Effect,
          ifelse(is.na(df_nd_all$ICER), "reference",
                 sprintf("%.1fM", df_nd_all$ICER / 1e6)),
          df_bur$CRCCasesper1000, df_bur$CRCCasesAvertedper1000,
          df_bur$CRCDeathsAvertedper1000, 100 * df_bur$MortalityReduction,
          format(round(df_bur$ScreeningColsper1000), big.mark = ",", trim = TRUE),
          format(round(df_bur$DiagnosticFUper1000), big.mark = ",", trim = TRUE),
          format(round(df_bur$SurveillanceColsper1000), big.mark = ",", trim = TRUE))
), file.path(dir_tab, "frontier_outcomes.md"))

df_cost <- read.csv("data-raw/cea_inputs/screen_costs_v3.csv", row.names = 1)
stopifnot(all(df_cost[[1]] == df_cost[[2]]))
v_item <- c(FIT = "FIT",
            `Screening colonoscopy without polypectomy` = "Screening colonoscopy, no polypectomy",
            `Diagnostic colonoscopy without polypectomy` = "Follow-up colonoscopy, no polypectomy",
            `Surveillance colonoscopy without polypectomy` = "Surveillance colonoscopy, no polypectomy",
            `Colonoscopy with polypectomy` = "Colonoscopy with polypectomy",
            `Colonoscopy for diagnosis of a cancer by symptoms` = "Colonoscopy for symptomatic CRC",
            `Fatal perforation` = "Fatal perforation",
            `Serious GI complication` = "Serious GI complication",
            `Other GI complication` = "Other GI complication",
            `Cardiovascular complication` = "Cardiovascular complication")
writeLines(c(
  "| Item | Cost, CLP |",
  "| :--- | ---: |",
  sprintf("| %s | %s |", unname(v_item),
          format(df_cost[names(v_item), 1], big.mark = ",", trim = TRUE))
), file.path(dir_tab, "costs_procedures.md"))

df_care <- read.csv("data-raw/cea_inputs/crc_care_costs.csv", row.names = 1)
stopifnot(all(as.matrix(df_care[1:4, ]) == as.matrix(df_care[5:8, ])))
v_phase <- c("Initial phase", "Continuing phase", "Terminal, died of CRC",
             "Terminal, died of other causes")
writeLines(c(
  "| Phase | Stage I | Stage II | Stage III | Stage IV |",
  "| :--- | ---: | ---: | ---: | ---: |",
  paste0("| ", v_phase, " | ", apply(format(as.matrix(df_care[1:4, ]),
                                           big.mark = ",", trim = TRUE), 1,
                                    paste, collapse = " | "), " |")
), file.path(dir_tab, "costs_cancer_care.md"))

df_pop <- read.csv("data-raw/fonasa_population_by_age_group_2025.csv")
df_old <- df_pop[df_pop$grupo_edad >= "40 a 44 años" &
                   !df_pop$grupo_edad %in% c("S.I.") , ]
df_old <- df_old[order(df_old$grupo_edad), ]
df_old$band <- ifelse(df_old$grupo_edad == "Más de 99 años", "100 or older",
                      sub(" a ", "-", sub(" años", "", df_old$grupo_edad)))
n_old   <- sum(df_old$n_beneficiarios)
n_total <- sum(df_pop$n_beneficiarios)
writeLines(c(
  "| Age band | FONASA beneficiaries, 2025 |",
  "| :--- | ---: |",
  sprintf("| %s | %s |", df_old$band, format(df_old$n_beneficiarios, big.mark = ",", trim = TRUE)),
  sprintf("| **Age 40 and older** | **%s** |", format(n_old, big.mark = ",", trim = TRUE)),
  sprintf("| **All ages** | **%s** |", format(n_total, big.mark = ",", trim = TRUE))
), file.path(dir_tab, "fonasa_bands.md"))

cat("\nbest choices\n"); print(df_best, row.names = FALSE)
cat("\nnon-dominated: all", nrow(df_nd_all), ", FIT incl. no screening",
    nrow(df_nd_fit), "\n")
cat("FONASA 40+ bands:", nrow(df_old), " total 40+:", n_old, " all ages:", n_total, "\n")
cat("Tables written to", dir_tab, "\n")
