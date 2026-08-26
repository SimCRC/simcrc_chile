# *****************************************************************************
#
# Script: benchmark_parallel_cea.R
#
# Purpose: Compare doSNOW/foreach vs mirai for a subset of CEA strategies.
#          Measures wall-clock time and peak system-wide RAM (all R processes)
#          for each approach using identical model inputs and the same n_cores.
#
# Output: ce_results/benchmark_parallel_cea.html
#         ce_results/benchmark_ram_traces.png
#
# First-run setup (uncomment if needed):
#   install.packages(c("mirai", "callr", "ggplot2", "gt", "patchwork"))
#
# Run from the simcrc_chile project root.
#
# *****************************************************************************

remove(list = ls())
gc()


# *****************************************************************************
#### 0. Config ####
# *****************************************************************************

# Number of strategies to benchmark. Set to Inf to run all strategies after
# deduplication (production scale). Use 12 for a quick smoke test.
N_BENCH <- Inf

# Population size for the benchmark run. Use 1e5 for a quick smoke test,
# 1e6 to benchmark at publication scale (slower but realistic).
N_POP_BENCH <- 1e6

# Number of parallel workers (identical for both approaches).
N_CORES <- ceiling(parallel::detectCores() / 2)

# RAM sampling interval in seconds.
RAM_POLL_S <- 0.5

# Output paths (written relative to project root).
OUT_HTML <- "ce_results/benchmark_parallel_cea.html"
OUT_PNG  <- "ce_results/benchmark_ram_traces.png"


# *****************************************************************************
#### 1. Load packages ####
# *****************************************************************************

library(data.table)
library(simcrc)
library(dplyr)
library(readr)
library(CRCmortality)
library(ceacrc)
library(stringr)

# Parallelization
library(parallel)
library(doSNOW)
library(foreach)
library(mirai)

# Benchmark tooling
library(callr)    # background R processes for RAM monitoring

# Reporting
library(ggplot2)
library(gt)
library(patchwork)


# *****************************************************************************
#### 2. Natural history (shared — run once) ####
# *****************************************************************************

simcrc_model_version <- paste0("SimCRC v", as.character(packageVersion("simcrc")))

cat(sprintf("[setup] SimCRC version: %s\n", simcrc_model_version))
cat(sprintf("[setup] Benchmark: %s strategies, n_pop = %s, %d cores\n",
            if (is.infinite(N_BENCH)) "all" else as.character(N_BENCH),
            formatC(N_POP_BENCH, format = "d", big.mark = ","), N_CORES))

# Load calibrated parameters
load("outputs/BayCANN_versions/Chile/Adenoma/F/v0.13.0/v0.13.0.20260406.1214/l_params_calibrated_sets_SimCRC_v0.13.0.20260406.1214_Adenoma_F.RData")
l_params_Min_MSE <- l_params_calibrated_sets$Min_MSE
l_params_all     <- load_params_init(fromFile = TRUE, filename = l_params_Min_MSE)

l_params_all$min_age_lesion_onset <- 10
l_params_all$mort_by_race         <- FALSE

cohort_age <- 40

df_lt_chile <- read.csv("data-raw/df_lifetable_2017_CH.csv")
colnames(df_lt_chile)[colnames(df_lt_chile) == "age"]           <- "Age"
colnames(df_lt_chile)[colnames(df_lt_chile) == "mortality_rate"] <- "mortality.rates"

cat(sprintf("[setup] Building population (n = %s)...\n",
            formatC(N_POP_BENCH, format = "d", big.mark = ",")))

dt_pop <- simcrc::get_dt_population(
  year          = 1980,
  byear         = 1980,
  p_female      = 1,
  p_white       = 0.8,
  n_pop         = N_POP_BENCH,
  dt_life_table_F = df_lt_chile,
  dt_life_table_M = NULL
)

cat("[setup] Running natural history...\n")
l_out_simcrc <- simcr_nathist_ssp_DES(
  l_params_all  = l_params_all,
  dt_pop        = dt_pop,
  SSP_pathway   = FALSE
)
dt_crc_pop <- l_out_simcrc$dt_crc_pop
cat(sprintf("[setup] dt_crc_pop: %.1f MB\n",
            object.size(dt_crc_pop) / 1024^2))


# *****************************************************************************
#### 3. Strategy subset ####
# *****************************************************************************

df_strategies_all <- read_csv("data-raw/df_CH_2026_strategies.csv",
                               show_col_types = FALSE)

l_to_remove <- ceacrc:::identify_strategies_to_remove(
  modality   = "COL",
  start_ages = c(45, 50, 55),
  stop_ages  = c(70, 75, 80, 85),
  intervals  = c(5, 10, 15)
)
df_strategies_all <- df_strategies_all |>
  filter(!strategy %in% l_to_remove)

# Use all strategies when N_BENCH = Inf; otherwise sample a stratified subset.
if (is.infinite(N_BENCH)) {
  set.seed(42)
  df_bench <- df_strategies_all |> slice_sample(prop = 1)
} else {
  n_fit <- round(N_BENCH * 0.6)
  n_col <- N_BENCH - n_fit
  set.seed(42)
  df_bench <- bind_rows(
    df_strategies_all |> filter(modality == "FIT") |> slice_sample(n = n_fit),
    df_strategies_all |> filter(modality == "COL") |> slice_sample(n = n_col)
  ) |> slice_sample(prop = 1)
}

n_fit_bench <- sum(df_bench$modality == "FIT")
n_col_bench <- sum(df_bench$modality == "COL")
cat(sprintf("[setup] Benchmark strategies: %d total (%d FIT, %d COL)\n",
            nrow(df_bench), n_fit_bench, n_col_bench))

# Output folder for benchmark CSVs (temp, cleaned between runs)
bench_out_dir <- normalizePath(
  "output/benchmark_run",
  mustWork = FALSE
)


# *****************************************************************************
#### 4. RAM monitor helpers ####
# *****************************************************************************
# Uses a callr background process to sample total RSS of all R sessions
# (main process + all workers) every RAM_POLL_S seconds.
# macOS: ps reports RSS in KB.

start_ram_monitor <- function(poll_s = 0.5) {
  f <- tempfile(fileext = ".csv")
  bg <- callr::r_bg(
    function(f, poll_s) {
      while (TRUE) {
        rss_kb <- tryCatch(
          as.numeric(system(
            "ps -A -o rss | awk 'NR>1 && $1+0>0 {sum+=$1} END{print sum}'",
            intern = TRUE
          )),
          error = function(e) NA_real_
        )
        cat(
          sprintf("%.6f,%.2f\n", as.numeric(Sys.time()), rss_kb / 1024),
          file = f, append = TRUE
        )
        Sys.sleep(poll_s)
      }
    },
    args      = list(f = f, poll_s = poll_s),
    supervise = TRUE
  )
  list(bg = bg, file = f, start_ts = as.numeric(Sys.time()))
}

stop_ram_monitor <- function(mon) {
  mon$bg$kill()
  Sys.sleep(0.2)   # allow final write
  if (!file.exists(mon$file) || file.size(mon$file) == 0) {
    return(list(peak_mb = NA_real_, baseline_mb = NA_real_, trace = NULL))
  }
  d <- tryCatch(
    read.csv(mon$file, header = FALSE,
             col.names = c("ts", "ram_mb"), stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(d) || nrow(d) == 0) {
    return(list(peak_mb = NA_real_, baseline_mb = NA_real_, trace = NULL))
  }
  d$ram_mb    <- suppressWarnings(as.numeric(d$ram_mb))
  d$elapsed_s <- d$ts - mon$start_ts
  list(
    peak_mb     = max(d$ram_mb,  na.rm = TRUE),
    baseline_mb = min(d$ram_mb,  na.rm = TRUE),
    delta_mb    = max(d$ram_mb,  na.rm = TRUE) - min(d$ram_mb, na.rm = TRUE),
    trace       = d
  )
}


# *****************************************************************************
#### 5. Benchmark: doSNOW / foreach ####
# *****************************************************************************

run_bench_dosnow <- function(df_bench, dt_crc_pop, cohort_age,
                              simcrc_model_version, out_dir, n_cores,
                              ram_poll_s = 0.5) {

  unlink(out_dir, recursive = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cat(sprintf("\n[doSNOW] Starting — %d strategies on %d cores\n",
              nrow(df_bench), n_cores))

  n_ids <- nrow(df_bench)
  mon   <- start_ram_monitor(ram_poll_s)

  cl  <- makeCluster(n_cores)
  registerDoSNOW(cl)
  pb  <- txtProgressBar(min = 0, max = n_ids, style = 3, file = stderr())
  opts <- list(progress = function(n) setTxtProgressBar(pb, n))

  t0 <- proc.time()

  foreach(
    i             = seq_len(n_ids),
    .packages     = c("data.table", "simcrc", "ceacrc", "dplyr"),
    .options.snow = opts
  ) %dopar% {

    run_id             <- df_bench[i, ]
    strategy_name      <- run_id$strategy
    screening_modality <- run_id$modality
    follow_up          <- run_id$follow_up
    screening_years    <- seq(run_id$age_to_begin_screening,
                              run_id$age_to_end_screening,
                              run_id$frequency_screening)

    capture.output({

      res_scr <- screening_detection(
        dt_crc_pop               = dt_crc_pop,
        p_coverage               = run_id$p_coverage,
        p_adherence              = run_id$p_adherence,
        age_to_begin_screening   = run_id$age_to_begin_screening,
        age_to_end_screening     = run_id$age_to_end_screening,
        frequency_screening      = run_id$frequency_screening,
        sens_small_adenoma       = run_id$sens_small_adenoma,
        sens_medium_adenoma      = run_id$sens_medium_adenoma,
        sens_large_adenoma       = run_id$sens_large_adenoma,
        sens_crc                 = run_id$sens_crc,
        sens_by_ad               = run_id$sens_by_ad,
        spec                     = run_id$spec,
        screening_reach          = run_id$screening_reach,
        p_reach_cecum            = run_id$p_reach_cecum,
        p_reach_ascending        = run_id$p_reach_ascending,
        p_reach_transverse       = run_id$p_reach_transverse,
        p_reach_descending       = run_id$p_reach_descending,
        p_reach_sigmoid          = run_id$p_reach_sigmoid,
        p_reach_rectum           = run_id$p_reach_rectum,
        p_death_scr              = run_id$p_death_scr,
        surveillance             = run_id$surveillance,
        confirmation             = run_id$follow_up,
        p_adherence_confirmation = run_id$p_adherence_confirmation,
        sens_small_adenoma_col   = run_id$sens_small_adenoma_col,
        sens_medium_adenoma_col  = run_id$sens_medium_adenoma_col,
        sens_large_adenoma_col   = run_id$sens_large_adenoma_col,
        sens_crc_col             = run_id$sens_crc_col,
        sens_by_ad_col           = run_id$sens_by_ad_col,
        spec_col                 = run_id$spec_col,
        confirmation_reach       = run_id$confirmation_reach,
        p_reach_cecum_conf       = run_id$p_reach_cecum_conf,
        p_reach_ascending_conf   = run_id$p_reach_ascending_conf,
        p_reach_transverse_conf  = run_id$p_reach_transverse_conf,
        p_reach_descending_conf  = run_id$p_reach_descending_conf,
        p_reach_sigmoid_conf     = run_id$p_reach_sigmoid_conf,
        p_reach_rectum_conf      = run_id$p_reach_rectum_conf,
        p_death_conf             = run_id$p_death_conf,
        optimize_memory          = FALSE
      )

      res_surv <- surveillance_detection(
        dt_pop_scr              = res_scr$dt_pop_screening,
        p_adherence             = run_id$p_adherence_surv,
        sens_small_adenoma      = run_id$sens_small_adenoma_surv,
        sens_medium_adenoma     = run_id$sens_medium_adenoma_surv,
        sens_large_adenoma      = run_id$sens_large_adenoma_surv,
        sens_crc                = run_id$sens_crc_surv,
        sens_by_ad              = run_id$sens_by_ad_surv,
        spec                    = run_id$spec_surv,
        surveillance_reach      = run_id$surveillance_reach,
        p_reach_cecum_surv      = run_id$p_reach_cecum_surv,
        p_reach_ascending_surv  = run_id$p_reach_ascending_surv,
        p_reach_transverse_surv = run_id$p_reach_transverse_surv,
        p_reach_descending_surv = run_id$p_reach_descending_surv,
        p_reach_sigmoid_surv    = run_id$p_reach_sigmoid_surv,
        p_reach_rectum_surv     = run_id$p_reach_rectum_surv,
        p_death_surv            = run_id$p_death_surv,
        optimize_memory         = TRUE
      )

      dt_export <- uspstf_summary(
        datatable          = res_surv,
        screening_modality = screening_modality,
        follow_up          = follow_up,
        screening_years    = screening_years,
        min_age            = cohort_age,
        max_age            = 100
      )
    })

    # Write output (same workload as production)
    f <- file.path(out_dir, paste0(strategy_name, ".csv"))
    write.csv(dt_export, f, row.names = FALSE)

    NULL
  }

  elapsed_s <- (proc.time() - t0)[["elapsed"]]
  close(pb)
  cat("\n")
  stopCluster(cl)

  ram <- stop_ram_monitor(mon)
  cat(sprintf("[doSNOW] Done in %.1f s (%.1f min) | peak RAM %.0f MB | delta %.0f MB\n",
              elapsed_s, elapsed_s / 60, ram$peak_mb, ram$delta_mb))

  list(approach = "doSNOW / foreach", elapsed_s = elapsed_s,
       peak_mb = ram$peak_mb, baseline_mb = ram$baseline_mb,
       delta_mb = ram$delta_mb, trace = ram$trace)
}


# *****************************************************************************
#### 6. Benchmark: mirai ####
# *****************************************************************************

run_bench_mirai <- function(df_bench, dt_crc_pop, cohort_age,
                             simcrc_model_version, out_dir, n_cores,
                             ram_poll_s = 0.5) {

  unlink(out_dir, recursive = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cat(sprintf("\n[mirai]  Starting — %d strategies on %d daemons\n",
              nrow(df_bench), n_cores))

  l_strats <- split(df_bench, seq_len(nrow(df_bench)))
  mon      <- start_ram_monitor(ram_poll_s)
  t0       <- proc.time()

  daemons(n_cores, seed = 3)

  everywhere(
    {
      library(data.table)
      library(simcrc)
      library(ceacrc)
      library(dplyr)
    },
    dt_crc_pop           = dt_crc_pop,
    cohort_age           = cohort_age,
    simcrc_model_version = simcrc_model_version,
    out_dir              = out_dir
  )

  mirai_map(
    l_strats,
    function(run_id) {

      strategy_name      <- run_id$strategy
      screening_modality <- run_id$modality
      follow_up          <- run_id$follow_up
      screening_years    <- seq(run_id$age_to_begin_screening,
                                run_id$age_to_end_screening,
                                run_id$frequency_screening)

      capture.output({

        res_scr <- screening_detection(
          dt_crc_pop               = dt_crc_pop,   # from everywhere()
          p_coverage               = run_id$p_coverage,
          p_adherence              = run_id$p_adherence,
          age_to_begin_screening   = run_id$age_to_begin_screening,
          age_to_end_screening     = run_id$age_to_end_screening,
          frequency_screening      = run_id$frequency_screening,
          sens_small_adenoma       = run_id$sens_small_adenoma,
          sens_medium_adenoma      = run_id$sens_medium_adenoma,
          sens_large_adenoma       = run_id$sens_large_adenoma,
          sens_crc                 = run_id$sens_crc,
          sens_by_ad               = run_id$sens_by_ad,
          spec                     = run_id$spec,
          screening_reach          = run_id$screening_reach,
          p_reach_cecum            = run_id$p_reach_cecum,
          p_reach_ascending        = run_id$p_reach_ascending,
          p_reach_transverse       = run_id$p_reach_transverse,
          p_reach_descending       = run_id$p_reach_descending,
          p_reach_sigmoid          = run_id$p_reach_sigmoid,
          p_reach_rectum           = run_id$p_reach_rectum,
          p_death_scr              = run_id$p_death_scr,
          surveillance             = run_id$surveillance,
          confirmation             = run_id$follow_up,
          p_adherence_confirmation = run_id$p_adherence_confirmation,
          sens_small_adenoma_col   = run_id$sens_small_adenoma_col,
          sens_medium_adenoma_col  = run_id$sens_medium_adenoma_col,
          sens_large_adenoma_col   = run_id$sens_large_adenoma_col,
          sens_crc_col             = run_id$sens_crc_col,
          sens_by_ad_col           = run_id$sens_by_ad_col,
          spec_col                 = run_id$spec_col,
          confirmation_reach       = run_id$confirmation_reach,
          p_reach_cecum_conf       = run_id$p_reach_cecum_conf,
          p_reach_ascending_conf   = run_id$p_reach_ascending_conf,
          p_reach_transverse_conf  = run_id$p_reach_transverse_conf,
          p_reach_descending_conf  = run_id$p_reach_descending_conf,
          p_reach_sigmoid_conf     = run_id$p_reach_sigmoid_conf,
          p_reach_rectum_conf      = run_id$p_reach_rectum_conf,
          p_death_conf             = run_id$p_death_conf,
          optimize_memory          = FALSE
        )

        res_surv <- surveillance_detection(
          dt_pop_scr              = res_scr$dt_pop_screening,
          p_adherence             = run_id$p_adherence_surv,
          sens_small_adenoma      = run_id$sens_small_adenoma_surv,
          sens_medium_adenoma     = run_id$sens_medium_adenoma_surv,
          sens_large_adenoma      = run_id$sens_large_adenoma_surv,
          sens_crc                = run_id$sens_crc_surv,
          sens_by_ad              = run_id$sens_by_ad_surv,
          spec                    = run_id$spec_surv,
          surveillance_reach      = run_id$surveillance_reach,
          p_reach_cecum_surv      = run_id$p_reach_cecum_surv,
          p_reach_ascending_surv  = run_id$p_reach_ascending_surv,
          p_reach_transverse_surv = run_id$p_reach_transverse_surv,
          p_reach_descending_surv = run_id$p_reach_descending_surv,
          p_reach_sigmoid_surv    = run_id$p_reach_sigmoid_surv,
          p_reach_rectum_surv     = run_id$p_reach_rectum_surv,
          p_death_surv            = run_id$p_death_surv,
          optimize_memory         = TRUE
        )

        dt_export <- uspstf_summary(
          datatable          = res_surv,
          screening_modality = screening_modality,
          follow_up          = follow_up,
          screening_years    = screening_years,
          min_age            = cohort_age,
          max_age            = 100
        )
      })

      f <- file.path(out_dir, paste0(strategy_name, ".csv"))
      write.csv(dt_export, f, row.names = FALSE)

      NULL
    }
  )[.stop, .progress]

  daemons(0)

  elapsed_s <- (proc.time() - t0)[["elapsed"]]
  ram <- stop_ram_monitor(mon)
  cat(sprintf("\n[mirai]  Done in %.1f s (%.1f min) | peak RAM %.0f MB | delta %.0f MB\n",
              elapsed_s, elapsed_s / 60, ram$peak_mb, ram$delta_mb))

  list(approach = "mirai — mirai_map()", elapsed_s = elapsed_s,
       peak_mb = ram$peak_mb, baseline_mb = ram$baseline_mb,
       delta_mb = ram$delta_mb, trace = ram$trace)
}


# *****************************************************************************
#### 7. Run benchmarks ####
# *****************************************************************************

cat("\n", strrep("=", 60), "\n")
cat("BENCHMARK START\n")
cat(strrep("=", 60), "\n")

r_dosnow <- run_bench_dosnow(
  df_bench             = df_bench,
  dt_crc_pop           = dt_crc_pop,
  cohort_age           = cohort_age,
  simcrc_model_version = simcrc_model_version,
  out_dir              = bench_out_dir,
  n_cores              = N_CORES,
  ram_poll_s           = RAM_POLL_S
)

gc()   # release worker memory before mirai run

r_mirai <- run_bench_mirai(
  df_bench             = df_bench,
  dt_crc_pop           = dt_crc_pop,
  cohort_age           = cohort_age,
  simcrc_model_version = simcrc_model_version,
  out_dir              = bench_out_dir,
  n_cores              = N_CORES,
  ram_poll_s           = RAM_POLL_S
)

# Tidy up temp output folder
unlink(bench_out_dir, recursive = TRUE)


# *****************************************************************************
#### 8. Summary table ####
# *****************************************************************************

df_results <- tibble::tibble(
  Approach    = c(r_dosnow$approach, r_mirai$approach),
  `Time (s)`  = c(r_dosnow$elapsed_s, r_mirai$elapsed_s),
  `Time (min)`= c(r_dosnow$elapsed_s, r_mirai$elapsed_s) / 60,
  `Peak RAM (MB)` = c(r_dosnow$peak_mb, r_mirai$peak_mb),
  `RAM delta (MB)` = c(r_dosnow$delta_mb, r_mirai$delta_mb)
) |>
  mutate(
    speedup    = `Time (s)`[1]    / `Time (s)`,
    ram_change = `Peak RAM (MB)` / `Peak RAM (MB)`[1]
  )

cat("\n", strrep("=", 60), "\n")
cat("RESULTS SUMMARY\n")
cat(strrep("=", 60), "\n")
print(df_results |> select(Approach, `Time (s)`, `Time (min)`, `Peak RAM (MB)`, `RAM delta (MB)`))
cat(sprintf("\nSpeedup (mirai vs doSNOW): %.2fx\n",
            df_results$speedup[df_results$Approach == "mirai — mirai_map()"]))
cat(sprintf("RAM delta change:           %.0f MB → %.0f MB (%.1f%%)\n",
            df_results$`RAM delta (MB)`[1],
            df_results$`RAM delta (MB)`[2],
            (df_results$`RAM delta (MB)`[2] / df_results$`RAM delta (MB)`[1] - 1) * 100))


# *****************************************************************************
#### 9. RAM trace plot ####
# *****************************************************************************

build_trace <- function(res, label) {
  if (is.null(res$trace)) return(NULL)
  res$trace |>
    mutate(approach = label, elapsed_s = elapsed_s)
}

df_traces <- bind_rows(
  build_trace(r_dosnow, "doSNOW / foreach"),
  build_trace(r_mirai,  "mirai — mirai_map()")
) |>
  filter(!is.na(ram_mb), ram_mb > 0)

pal <- c("doSNOW / foreach"      = "#F97316",
         "mirai — mirai_map()"   = "#3B82F6")

p_ram <- ggplot(df_traces, aes(x = elapsed_s, y = ram_mb, colour = approach)) +
  geom_line(linewidth = 0.9, alpha = 0.85) +
  scale_colour_manual(values = pal) +
  scale_x_continuous(labels = scales::label_number(suffix = "s")) +
  scale_y_continuous(labels = scales::label_number(suffix = " MB")) +
  labs(
    title    = sprintf("System-wide R RAM — %s strategies, n_pop = %s, %d cores",
                       nrow(df_bench),
                       formatC(N_POP_BENCH, format = "d", big.mark = ","),
                       N_CORES),
    subtitle = "All R processes combined (main + workers), sampled every 0.5 s",
    x        = "Elapsed time (s)",
    y        = "Total RSS (MB)",
    colour   = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position    = "bottom",
    panel.grid.minor   = element_blank(),
    plot.title         = element_text(face = "bold", size = 11),
    plot.subtitle      = element_text(colour = "#64748b", size = 9)
  )

ggsave(OUT_PNG, p_ram, width = 8, height = 4, units = "in", dpi = 300)
cat(sprintf("Saved: %s\n", OUT_PNG))


# *****************************************************************************
#### 10. HTML report ####
# *****************************************************************************

tbl_results <- df_results |>
  mutate(speedup = sprintf("%.2fx", speedup)) |>
  select(
    Approach,
    `Time (s)`,
    `Time (min)`,
    `Peak RAM (MB)`,
    `RAM delta (MB)`,
    speedup
  ) |>
  gt(rowname_col = "Approach") |>
  tab_header(
    title    = "Parallelization benchmark — CEA Chile",
    subtitle = sprintf("%d strategies · n_pop = %s · %d cores · SimCRC %s",
                       nrow(df_bench),
                       formatC(N_POP_BENCH, format = "d", big.mark = ","),
                       N_CORES,
                       simcrc_model_version)
  ) |>
  fmt_number(columns = `Time (s)`,       decimals = 1) |>
  fmt_number(columns = `Time (min)`,     decimals = 2) |>
  fmt_number(columns = `Peak RAM (MB)`,  decimals = 0) |>
  fmt_number(columns = `RAM delta (MB)`, decimals = 0) |>
  cols_label(
    `Time (s)`       = "Wall time (s)",
    `Time (min)`     = "Wall time (min)",
    `Peak RAM (MB)`  = "Peak RAM (MB)",
    `RAM delta (MB)` = "RAM delta (MB)",
    speedup          = "Speedup vs doSNOW"
  ) |>
  tab_style(
    style     = cell_fill(color = "#DBEAFE"),
    locations = cells_body(rows = Approach == "mirai — mirai_map()")
  ) |>
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(columns = speedup,
                           rows    = Approach == "mirai — mirai_map()")
  ) |>
  tab_source_note(md(paste0(
    "**RAM delta** = peak − baseline (all R processes). ",
    "Baseline includes the main process with dt_crc_pop already in memory. ",
    "For doSNOW, workers serialize dt_crc_pop per task, increasing delta. ",
    "For mirai, `everywhere()` loads it once per daemon — delta reflects only compute overhead."
  ))) |>
  tab_options(
    table.font.size        = px(13),
    heading.align          = "left",
    column_labels.font.weight = "bold"
  )

gtsave(tbl_results, OUT_HTML)
cat(sprintf("Saved: %s\n", OUT_HTML))

cat("\nDone.\n")
