# *****************************************************************************
#
# Script: 14_calibration2_slides.R
#
# Purpose: Figures and tables for the Chilean calibration results deck
#          (slides/calibration2_results_2026-11-11.qmd).
#            - Setup sections: Chilean targets vs the US v0.15.0 model run in
#              the Chilean population, with small prior panels.
#            - Results: dwell and sojourn of both best sets in the Chilean
#              population, and symptom detection rates of both best sets.
#
# Run from the simcrc_chile root:
#   Rscript analysis/figures/14_calibration2_slides.R
#
# Inputs:
#   Chilean calibration: outputs/BayCANN_versions/Chile/Adenoma/F/v0.15.0.901/
#     v0.15.0.901.20261111.2051 (Min_MSE set = Min_AbsolutErr set, draw 693)
#   US v0.15.0.20261026.130537: data-raw/ posteriors, priors and calibrated
#     sets (Min_AbsolutErr; the US run saved no Min_MSE set)
#
# *****************************************************************************

suppressMessages({
  library(tidyverse)
  library(simcrc)
  library(data.table)
})
suppressMessages({
  source("analysis/calibration/01_calibration_setup.R")
  source("R/01_model_input_functions.R")
  source("R/02_calibration_functions.R")
})

dir_out <- "slides/calibration2_2026-11-11"
dir_fig <- file.path(dir_out, "figures")
dir_tab <- file.path(dir_out, "tables")
dir.create(dir_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_tab, recursive = TRUE, showWarnings = FALSE)

cal_run <- "v0.15.0.901.20261111.2051"  # second attempt at calibrating SimCRC v0.15
dir_run <- file.path("outputs/BayCANN_versions/Chile/Adenoma/F/v0.15.0.901", cal_run)
us_tag  <- "v0.15.0.20261026.130537"

col_us <- "#2a78d6"
col_ch <- "#e34948"
col_tg <- "#6b6a66"

ld_sets <- function(f) { e <- new.env(); load(f, e); e$l_params_calibrated_sets }
us_sets <- ld_sets(sprintf("data-raw/l_params_calibrated_sets_SimCRC_%s_Adenoma_F.RData", us_tag))
ch_sets <- ld_sets(file.path(dir_run, sprintf("l_params_calibrated_sets_SimCRC_%s_Adenoma_F.RData", cal_run)))
us_best <- us_sets$Min_AbsolutErr
ch_best <- ch_sets$Min_MSE
stopifnot(identical(ch_sets$Min_MSE$id_draw, ch_sets$Min_AbsolutErr$id_draw))

e_pr <- new.env()
load(file.path(dir_run, sprintf("l_params_priors_SimCRC_SimCRC_%s_Adenoma_F.RData", cal_run)), e_pr)
ch_prior <- get(ls(e_pr)[1], e_pr)
us_post  <- read.csv(sprintf("data-raw/dt_calibrated_posteriors_SimCRC_%s_Adenoma_F.csv", us_tag))
v_par    <- ch_prior$names

short_name <- function(x) x %>%
  str_replace("IndividuakRiskMultVariance", "Individual risk variance") %>%
  str_replace("AdNaturalHistoryPropensity_Gaussian_Variance", "Adenoma propensity variance") %>%
  str_replace("alpha_lesion_adenoma", "Adenoma onset: alpha") %>%
  str_replace("beta_age", "Adenoma onset: beta (age)") %>%
  str_replace("AdGrowth_Exp_Rate_1_to_6_(.)", "Adenoma growth 1-6 mm, \\1") %>%
  str_replace("AdGrowth_Exp_Rate_6_to_10_(.)", "Adenoma growth 6-10 mm, \\1") %>%
  str_replace("PreclinCancerProg_Exp_Rate_S(\\d)S(\\d)_(.)", "Preclinical progression S\\1-S\\2, \\3") %>%
  str_replace("CancerOnset_Gompertz_Shape_(.)", "Cancer onset shape, \\1") %>%
  str_replace("CancerOnset_Gompertz_Rate_(.)", "Cancer onset rate, \\1") %>%
  str_replace("pSxDetS1_(.)", "Symptom detection prob. S1, \\1") %>%
  str_replace("hr_SxDetS(\\d)S(\\d)_(.)", "Symptom detection HR S\\1/S\\2, \\3")

# ---- 2. Best sets in the Chilean population --------------------------------
su <- l_model_adenoma_Chile
make_pop <- function() {
  set.seed(20220906)
  dt <- get_dt_population(year = su$cohort_year, byear = su$birth_year, p_female = su$p_female,
                          p_white = su$p_white, n_pop = su$n_pop,
                          dt_life_table_F = su$life_table_fem, dt_life_table_M = su$life_table_male)
  dt
}
init_for <- function(set) {
  l <- if (!is.null(set$l_params_init)) set$l_params_init else su$l_params_init
  l$min_age_lesion_onset <- 15
  l$mort_by_race <- FALSE
  l
}
run_set <- function(set, dt_pop) {
  set.seed(1)
  calibration_out(v_params_calib = unlist(set$Min[v_par]), l_params_all = init_for(set),
                  dt_pop = dt_pop, id_draw = 1)
}
l_run <- list(US = us_sets, Chile = c(ch_sets, list(l_params_init = NULL)))
l_run$US$Min    <- us_best
l_run$Chile$Min <- ch_best

f_runs <- file.path(dir_tab, "best_set_runs_chile_pop.rds")
if (file.exists(f_runs) && !identical(Sys.getenv("RERUN_SIMS"), "1")) {
  # Reuse saved simulations; set RERUN_SIMS=1 to simulate again.
  out_oc <- readRDS(f_runs)$out_oc
} else {
  dt_pop_oc <- make_pop()
  out_oc    <- lapply(l_run, run_set, dt_pop = dt_pop_oc)
  saveRDS(list(out_oc = out_oc), f_runs)
}

# ---- 2b. Best-set outputs vs Chilean targets -------------------------------
tg <- read.csv(su$targest_file) %>% filter(target_names %in% su$selected_targets)
df_val <- bind_rows(lapply(names(out_oc), function(k) {
  o <- as.data.frame(out_oc[[k]])
  tibble(target_names = tg$target_names, model = as.numeric(o[1, tg$target_names]), set = k)
})) %>% left_join(tg, by = "target_names")
lab_set <- c(Chile = "Chilean best set", US = "US best set (v0.15.0), run in Chile")
df_val$set_lab <- factor(lab_set[df_val$set], levels = lab_set)

grp_title <- c(CRCInc_Overall = "CRC incidence, overall", CRCInc_P = "CRC incidence, proximal",
               CRCInc_D = "CRC incidence, distal", CRCInc_R = "CRC incidence, rectal",
               Prev0Ad = "No adenoma", Prev1Ad = "1 adenoma", Prev2Ad = "2 adenomas",
               Prev3PlusAd = "3+ adenomas", PrevPreclin = "Preclinical cancer")
cat_title <- c(StageDist_P = "Stage at diagnosis, proximal", StageDist_D = "Stage at diagnosis, distal",
               StageDist_R = "Stage at diagnosis, rectal", Size_P = "Adenoma size, proximal",
               Size_D = "Adenoma size, distal", Size_R = "Adenoma size, rectal")
d_cat <- df_val %>% filter(target_groups %in% names(cat_title)) %>%
  mutate(panel = factor(cat_title[target_groups], levels = cat_title),
         cat = case_when(str_detect(target_names, "^StageDistribution") ~ paste("Stage", str_sub(target_names, -1)),
                         str_detect(target_names, "SizeLR") ~ "Low risk",
                         str_detect(target_names, "SizeMR") ~ "Medium risk",
                         str_detect(target_names, "SizeHR") ~ "High risk"),
         cat = factor(cat, levels = c(paste("Stage", 1:4), "Low risk", "Medium risk", "High risk")))
# ---- 3. Dwell and sojourn ----------------------------------------------------
df_dwell <- bind_rows(lapply(c("US", "Chile"), function(k) {
  tibble(set = k, Adenoma_dwell = out_oc[[k]]$Adenoma_dwell, Sojourn_time = out_oc[[k]]$Sojourn_time,
         Total_dwell = out_oc[[k]]$Total_dwell)
}))
write.csv(df_dwell, file.path(dir_tab, "dwell_sojourn_best_sets.csv"), row.names = FALSE)

# Slide table: both best sets in the same Chilean 2017 population (Chilean life
# table, normal other-cause death), so the difference reflects the natural history.
us_d <- filter(df_dwell, set == "US"); ch_d <- filter(df_dwell, set == "Chile")
writeLines(c(
  "| Mean time (years) | US best set | Chile best set | |",
  "| :--- | ---: | ---: | ---: |",
  sprintf("| %s | %.1f | %.1f | %+.1f |", c("Adenoma dwell (adenoma onset to preclinical cancer)",
                                           "Sojourn (preclinical cancer to symptom detection)",
                                           "Total dwell (adenoma onset to symptom detection)"),
          c(us_d$Adenoma_dwell, us_d$Sojourn_time, us_d$Total_dwell),
          c(ch_d$Adenoma_dwell, ch_d$Sojourn_time, ch_d$Total_dwell),
          c(ch_d$Adenoma_dwell - us_d$Adenoma_dwell, ch_d$Sojourn_time - us_d$Sojourn_time,
            ch_d$Total_dwell - us_d$Total_dwell))
), file.path(dir_tab, "dwell_sojourn_table.md"))

# ---- 4. Symptom detection rates ---------------------------------------------
v_loc <- c(P = "Proximal", D = "Distal", R = "Rectal")
df_rate <- bind_rows(lapply(c(US = "US", Chile = "Chile"), function(k) {
  l_par <- if (k == "US") us_best else ch_best
  bind_rows(lapply(names(v_loc), function(loc) {
    v_hr   <- c(l_par[[paste0("hr_SxDetS2S1_", loc)]], l_par[[paste0("hr_SxDetS3S2_", loc)]],
                l_par[[paste0("hr_SxDetS4S3_", loc)]])
    v_prog <- c(l_par[[paste0("PreclinCancerProg_Exp_Rate_S1S2_", loc)]],
                l_par[[paste0("PreclinCancerProg_Exp_Rate_S2S3_", loc)]],
                l_par[[paste0("PreclinCancerProg_Exp_Rate_S3S4_", loc)]], 0)
    v_det  <- cumprod(c(-log(1 - l_par[[paste0("pSxDetS1_", loc)]]), v_hr))
    m_q <- matrix(0, 8, 8)
    m_q[cbind(1:4, 5:8)] <- v_det
    m_q[cbind(1:3, 2:4)] <- v_prog[1:3]
    diag(m_q) <- -rowSums(m_q)
    p_dx <- rowSums(as.matrix(Matrix::expm(m_q))[1:4, 5:8])
    tibble(set = k, loc = loc, stage = 1:4, r_det = v_det, r_prog = v_prog, p_dx = p_dx)
  }))
}))
stopifnot(all(df_rate$p_dx >= 0 & df_rate$p_dx <= 1))
write.csv(df_rate, file.path(dir_tab, "sx_detection_rates_best_sets.csv"), row.names = FALSE)

writeLines(c(paste("Chile:", cal_run, "Min_MSE (= Min_AbsolutErr), draw", ch_best$id_draw),
             paste("US:", us_tag, "Min_AbsolutErr (no Min_MSE set saved)"),
             "Population: Chilean 2017 cohort, 100,000 women, seed 20220906; min lesion onset age 15.",
             paste("Installed simcrc:", packageVersion("simcrc"))),
           file.path(dir_tab, "run_details.txt"))
cat("Done:", dir_out, "\n")

# ---- 5. Setup-section figures: Chilean data vs the US model run in Chile ----
# Mexico City deck style: each section compares the local targets with
# base-case (US) SimCRC before deciding whether to change anything.
# US posterior run in Chile: 200 draws from the US v0.15.0 posterior, each
# simulated in the same Chilean population, to show the spread of the US model.
f_post <- file.path(dir_tab, "us_posterior_runs_chile_pop.rds")
if (file.exists(f_post) && !identical(Sys.getenv("RERUN_SIMS"), "1")) {
  df_us_post_out <- readRDS(f_post)
} else {
  if (!exists("dt_pop_oc")) dt_pop_oc <- make_pop()
  set.seed(20261111)
  v_draw <- sample(nrow(us_post), 200)
  l_init_us <- init_for(us_sets)
  l_out <- parallel::mclapply(v_draw, function(i) {
    set.seed(i)
    o <- tryCatch(calibration_out(v_params_calib = unlist(us_post[i, v_par]), l_params_all = l_init_us,
                                  dt_pop = dt_pop_oc, id_draw = i), error = function(e) NULL)
    if (is.null(o)) return(NULL)
    o <- as.data.frame(o)
    data.frame(draw = i, target_names = tg$target_names, model = as.numeric(o[1, tg$target_names]))
  }, mc.cores = 8)
  df_us_post_out <- bind_rows(l_out)
  saveRDS(df_us_post_out, f_post)
}
cat("US posterior draws simulated:", length(unique(df_us_post_out$draw)), "\n")
# Band = full range of the simulated US posterior draws plus the US best set
# (a posterior draw in the tails for some outputs), so it always contains the bar.
best_us_out <- as.data.frame(out_oc$US)
df_band_src <- bind_rows(df_us_post_out,
                         data.frame(draw = 0, target_names = tg$target_names, model = as.numeric(best_us_out[1, tg$target_names])))
band_us <- df_band_src %>% group_by(target_names) %>%
  summarise(lo95 = min(model, na.rm = TRUE), hi95 = max(model, na.rm = TRUE), .groups = "drop")

lab_us_in_ch <- "US model (v0.15.0 best set) run in Chile"
theme_sec <- theme_bw(base_size = 13) +
  theme(legend.position = "top", legend.justification = "left", panel.grid.minor = element_blank(),
        strip.background = element_blank(), strip.text = element_text(hjust = 0, face = "bold"),
        plot.title = element_text(face = "bold"))
d_us <- df_val %>% filter(set == "US")

sec_num <- function(groups, file, title, ncol, h) {
  d <- d_us %>% filter(target_groups %in% groups) %>% mutate(panel = factor(grp_title[target_groups], levels = grp_title[groups])) %>%
    left_join(band_us, by = "target_names")
  g <- ggplot(d, aes(x = age)) +
    geom_ribbon(aes(ymin = lo95, ymax = hi95, fill = "US posterior range"), alpha = 0.25) +
    geom_errorbar(aes(ymin = stopping_lower_bounds, ymax = stopping_upper_bounds, colour = "Chilean target"), width = 1.8, linewidth = 0.8, alpha = 0.6) +
    geom_point(aes(y = targets, colour = "Chilean target"), size = 2.2) +
    geom_line(aes(y = model, colour = lab_us_in_ch), linewidth = 1.1) +
    facet_wrap(~panel, scales = "free_y", ncol = ncol) +
    scale_colour_manual(NULL, values = setNames(c(col_ch, col_us), c("Chilean target", lab_us_in_ch)),
                        guide = guide_legend(override.aes = list(linetype = c(0, 1), shape = c(16, NA)))) +
    scale_fill_manual(NULL, values = c("US posterior range" = col_us)) +
    labs(x = "age", y = NULL, title = title) + theme_sec
  ggsave(file.path(dir_fig, file), g, width = 13, height = h, dpi = 200)
}
sec_num(c("Prev0Ad", "Prev1Ad", "Prev2Ad", "Prev3PlusAd"), "sec_adenoma_prevalence.png",
        "Adenoma prevalence: Chilean targets vs the US model in Chile", 2, 7.4)
sec_num(c("CRCInc_Overall", "CRCInc_P", "CRCInc_D", "CRCInc_R"), "sec_crc_incidence.png",
        "CRC incidence per 100,000: Chilean targets vs the US model in Chile", 2, 7.4)

d_st <- d_cat %>% filter(set == "US", str_detect(target_groups, "^StageDist")) %>% left_join(band_us, by = "target_names")
# Slide table: stage at diagnosis, Chilean targets vs the US best set run in Chile
st_tab <- d_st %>% transmute(loc = str_match(target_names, "_([PDR])_S")[, 2], stage = str_sub(target_names, -1),
                             ch = targets, us = model)
cell_st <- function(l, k, w) sprintf("%.0f%%", 100 * st_tab[[w]][st_tab$loc == l & st_tab$stage == k])
writeLines(c(
  "| Stage at diagnosis | Proximal: Chile | Proximal: US model | Distal: Chile | Distal: US model | Rectal: Chile | Rectal: US model |",
  "| :--- | ---: | ---: | ---: | ---: | ---: | ---: |",
  sapply(1:4, function(k) paste0("| Stage ", k, " | ", paste(sapply(c("P", "D", "R"), function(l)
    paste(cell_st(l, k, "ch"), cell_st(l, k, "us"), sep = " | ")), collapse = " | "), " |"))
), file.path(dir_tab, "stage_distribution_table.md"))

cat("Section figures done\n")

# ---- 6. Small prior panels for the setup sections ---------------------------
# Chilean priors from 01_calibration_setup.R (second attempt) vs the middle 95%
# of the US posterior.
pr3 <- l_model_adenoma_Chile$params_priors
mini_name <- function(x) short_name(x) %>% str_replace("Symptom detection prob\\. ", "Detection prob. ") %>%
  str_replace("Symptom detection HR ", "Detection HR ") %>% str_replace("Preclinical progression ", "Progression ") %>%
  str_replace("Adenoma onset: ", "Onset ")
mini_ranges <- function(params, file, title) {
  d <- tibble(param = params, lb = pr3$lb[params], ub = pr3$ub[params],
              plo = sapply(params, function(p) quantile(us_post[[p]], 0.025)),
              phi = sapply(params, function(p) quantile(us_post[[p]], 0.975))) %>%
    mutate(same = abs(lb - plo) < 1e-9 * pmax(1, abs(plo)) & abs(ub - phi) < 1e-9 * pmax(1, abs(phi)),
           status = ifelse(same, "Same as US", "Changed"),
           lo = pmin(lb, plo), hi = pmax(ub, phi), w = hi - lo,
           x1 = (lb - lo) / w, x2 = (ub - lo) / w, p1 = (plo - lo) / w, p2 = (phi - lo) / w,
           label = factor(mini_name(param), levels = rev(mini_name(params))))
  g <- ggplot(d, aes(y = label)) +
    geom_segment(aes(x = x1, xend = x2, yend = label, colour = "Chilean prior"), linewidth = 5) +
    geom_segment(aes(x = p1, xend = p2, yend = label, colour = "US posterior 95%"), linewidth = 1.8, lineend = "round") +
    geom_text(aes(x = 1.06, label = status, fontface = ifelse(same, "plain", "bold")), hjust = 0, size = 3.4,
              colour = ifelse(d$same, "#55534f", "#b3312f")) +
    scale_colour_manual(NULL, values = c("Chilean prior" = "#f2b8b7", "US posterior 95%" = col_us)) +
    scale_x_continuous(limits = c(0, 1.5), breaks = NULL, expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "top", legend.justification = "left", panel.grid = element_blank(),
          plot.title = element_text(face = "bold", size = 12), plot.background = element_rect(fill = "white", colour = NA))
  ggsave(file.path(dir_fig, file), g, width = 5.4, height = 1.3 + 0.32 * length(params), dpi = 200)
}
mini_ranges(c("alpha_lesion_adenoma", "beta_age", "IndividuakRiskMultVariance", "AdNaturalHistoryPropensity_Gaussian_Variance",
              paste0("AdGrowth_Exp_Rate_1_to_6_", c("P", "D", "R")), paste0("AdGrowth_Exp_Rate_6_to_10_", c("P", "D", "R"))),
            "mini_priors_adenoma.png", "Adenoma onset and growth priors")
# ---- 7. Adenoma size, pooled over location ----------------------------------
# The Chilean size target is one value for the whole colon, so the model is
# compared pooled over proximal, distal and rectal adenomas, using the same
# size function and age weights as calibration_out().
pooled_size <- function(v_par_i, l_init, dt_pop, seed) {
  set.seed(seed)
  l_all <- update_param_list(l_params_all = l_init, params_updated = v_par_i)
  d <- simcrc::simcr_nathist_ssp_DES(l_params_all = l_all, dt_pop = dt_pop, SSP_pathway = FALSE)$dt_crc_pop
  lp <- get_living_pop(dt_crc_pop = d, min_age = 0, max_age = 100)
  z <- simcrc::get_size_distribution_ssp(dt_crc_pop = d, dt_living_pop = lp, ad_location = c("P", "D", "R"),
                                         v_age = c(50:75), lesion_pathway = "Adenoma", age_groups = FALSE)
  c(low = z[, prop_S], medium = z[, prop_M], high = z[, prop_L])
}
f_pool <- file.path(dir_tab, "us_pooled_size_chile_pop.rds")
if (file.exists(f_pool) && !identical(Sys.getenv("RERUN_SIMS"), "1")) {
  l_pool <- readRDS(f_pool)
} else {
  if (!exists("dt_pop_oc")) dt_pop_oc <- make_pop()
  l_init_us <- init_for(us_sets)
  best_pool <- pooled_size(unlist(us_best[v_par]), l_init_us, dt_pop_oc, 1)
  v_draw_pool <- unique(df_us_post_out$draw)
  post_pool <- do.call(rbind, parallel::mclapply(v_draw_pool, function(i)
    tryCatch(pooled_size(unlist(us_post[i, v_par]), l_init_us, dt_pop_oc, i), error = function(e) rep(NA_real_, 3)),
    mc.cores = 8))
  l_pool <- list(best = best_pool, post = post_pool)
  saveRDS(l_pool, f_pool)
}
tg_size <- tg %>% filter(target_groups == "Size_P") %>%
  mutate(cat = case_when(str_detect(target_names, "SizeLR") ~ "low", str_detect(target_names, "SizeMR") ~ "medium", TRUE ~ "high"))
d_pool <- tibble(cat = c("low", "medium", "high"),
                 model = as.numeric(l_pool$best[c("low", "medium", "high")]),
                 lo95 = apply(rbind(l_pool$post, l_pool$best)[, c("low", "medium", "high")], 2, min, na.rm = TRUE),
                 hi95 = apply(rbind(l_pool$post, l_pool$best)[, c("low", "medium", "high")], 2, max, na.rm = TRUE)) %>%
  left_join(tg_size %>% select(cat, targets, stopping_lower_bounds, stopping_upper_bounds), by = "cat") %>%
  mutate(cat = factor(c(low = "Low risk", medium = "Medium risk", high = "High risk")[cat], levels = c("Low risk", "Medium risk", "High risk")))
write.csv(d_pool, file.path(dir_tab, "adenoma_size_pooled_us_in_chile.csv"), row.names = FALSE)
g_sz <- ggplot(d_pool, aes(x = cat)) +
  geom_col(aes(y = targets, fill = "Chilean target"), width = 0.36, position = position_nudge(x = -0.19)) +
  geom_errorbar(aes(ymin = stopping_lower_bounds, ymax = stopping_upper_bounds), position = position_nudge(x = -0.19),
                width = 0.14, linewidth = 0.6, colour = "#8a2b2a", alpha = 0.7) +
  geom_col(aes(y = model, fill = lab_us_in_ch), width = 0.36, position = position_nudge(x = 0.19)) +
  geom_errorbar(aes(ymin = lo95, ymax = hi95), position = position_nudge(x = 0.19), width = 0.14, linewidth = 0.7, colour = "#17406f") +
  geom_text(aes(y = stopping_upper_bounds, label = sprintf("%.0f%%", 100 * targets)), position = position_nudge(x = -0.19), vjust = -0.5, size = 4.2) +
  geom_text(aes(y = pmax(model, hi95), label = sprintf("%.0f%%", 100 * model)), position = position_nudge(x = 0.19), vjust = -0.5, size = 4.2) +
  scale_fill_manual(NULL, values = setNames(c(col_ch, col_us), c("Chilean target", lab_us_in_ch))) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "proportion of adenomas, all locations", title = "Adenoma size, whole colon: Chilean target vs the US model in Chile",
       subtitle = "Red whisker: Chilean fitting bounds. Blue whisker: range of 200 US posterior draws and the best set.") +
  theme_sec + theme(panel.grid.major.x = element_blank())
ggsave(file.path(dir_fig, "sec_adenoma_size.png"), g_sz, width = 11, height = 6.2, dpi = 200)
mini_ranges(c(paste0("AdGrowth_Exp_Rate_1_to_6_", c("P", "D", "R")), paste0("AdGrowth_Exp_Rate_6_to_10_", c("P", "D", "R"))),
            "mini_priors_growth.png", "Adenoma growth priors")
cat("Size figure done\n")

# ---- 8. Cancer onset hazard by age: US posterior vs Chilean prior -----------
# Yearly hazard of cancer onset for an adenoma of 10 mm or more:
# rate * exp(shape * age) (simcrc CRC_onset_risk, adenoma risk multiplier 1).
# Chilean prior range: hazard at the lowest and highest corners of the box.
v_age_h <- 40:90
d_h <- bind_rows(lapply(c(P = "Proximal", D = "Distal", R = "Rectal") %>% names(), function(l) {
  rn <- paste0("CancerOnset_Gompertz_Rate_", l); sn <- paste0("CancerOnset_Gompertz_Shape_", l)
  hu <- sapply(v_age_h, function(a) quantile(us_post[[rn]] * exp(us_post[[sn]] * a), c(0.025, 0.975)))
  tibble(loc = c(P = "Proximal", D = "Distal", R = "Rectal")[l], age = v_age_h, us_lo = hu[1, ], us_hi = hu[2, ],
         ch_lo = pr3$lb[rn] * exp(pr3$lb[sn] * v_age_h), ch_hi = pr3$ub[rn] * exp(pr3$ub[sn] * v_age_h))
})) %>% mutate(loc = factor(loc, levels = c("Proximal", "Distal", "Rectal")))
write.csv(d_h, file.path(dir_tab, "cancer_onset_hazard_us_vs_chile_prior.csv"), row.names = FALSE)
g_h <- ggplot(d_h, aes(age)) +
  geom_ribbon(aes(ymin = ch_lo, ymax = ch_hi, fill = "Chilean prior range", colour = "Chilean prior range"), alpha = 0.35, linewidth = 0.5) +
  geom_ribbon(aes(ymin = us_lo, ymax = us_hi, fill = "US posterior 95%", colour = "US posterior 95%"), alpha = 0.45, linewidth = 0.5) +
  facet_wrap(~loc, ncol = 1) +
  scale_fill_manual(NULL, values = c("Chilean prior range" = col_ch, "US posterior 95%" = col_us)) +
  scale_colour_manual(NULL, values = c("Chilean prior range" = col_ch, "US posterior 95%" = col_us)) +
  scale_y_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 0.5), labels = c("0.0001", "0.001", "0.01", "0.1", "0.5")) +
  labs(title = "Cancer onset hazard by age", subtitle = "Yearly chance an adenoma of 10 mm+ becomes cancer (log scale)", x = "age", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top", legend.justification = "left", panel.grid.minor = element_blank(),
        strip.text = element_text(hjust = 0, face = "bold"), plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11), plot.background = element_rect(fill = "white", colour = NA))
ggsave(file.path(dir_fig, "mini_onset_hazard.png"), g_h, width = 5.4, height = 7, dpi = 200)
cat("Hazard panel done\n")

# ---- 9. Stage slide: progression (same as US) and symptom detection priors --
mini_grouped <- function(groups, file, title) {
  d <- bind_rows(lapply(names(groups), function(gname) {
    params <- groups[[gname]]
    tibble(group = gname, param = params, lb = pr3$lb[params], ub = pr3$ub[params],
           plo = sapply(params, function(p) quantile(us_post[[p]], 0.025)),
           phi = sapply(params, function(p) quantile(us_post[[p]], 0.975)))
  })) %>%
    mutate(same = abs(lb - plo) < 1e-9 * pmax(1, abs(plo)) & abs(ub - phi) < 1e-9 * pmax(1, abs(phi)),
           status = ifelse(same, "Same as US", "Changed"),
           lo = pmin(lb, plo), hi = pmax(ub, phi), w = hi - lo,
           x1 = (lb - lo) / w, x2 = (ub - lo) / w, p1 = (plo - lo) / w, p2 = (phi - lo) / w,
           label = factor(mini_name(param), levels = rev(mini_name(unlist(groups)))),
           group = factor(group, levels = names(groups)))
  g <- ggplot(d, aes(y = label)) +
    geom_segment(aes(x = x1, xend = x2, yend = label, colour = "Chilean prior"), linewidth = 4) +
    geom_segment(aes(x = p1, xend = p2, yend = label, colour = "US posterior 95%"), linewidth = 1.5, lineend = "round") +
    geom_text(aes(x = 1.06, label = status, fontface = ifelse(same, "plain", "bold")), hjust = 0, size = 3.2,
              colour = ifelse(d$same, "#55534f", "#b3312f")) +
    facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    scale_colour_manual(NULL, values = c("Chilean prior" = "#f2b8b7", "US posterior 95%" = col_us)) +
    scale_x_continuous(limits = c(0, 1.5), breaks = NULL, expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 10.5) +
    theme(legend.position = "top", legend.justification = "left", panel.grid = element_blank(),
          strip.text.y = element_blank(), panel.spacing = unit(0.9, "lines"),
          plot.title = element_text(face = "bold", size = 12), plot.background = element_rect(fill = "white", colour = NA))
  g <- g + geom_text(data = distinct(d, group), aes(x = 0, y = Inf, label = group), inherit.aes = FALSE,
                     hjust = 0, vjust = -0.4, fontface = "bold", size = 3.6)
  ggsave(file.path(dir_fig, file), g, width = 5.4, height = 1.6 + 0.26 * length(unlist(groups)), dpi = 200)
}
mini_grouped(list(
  "Preclinical progression" = c(paste0("PreclinCancerProg_Exp_Rate_S1S2_", c("P", "D", "R")), paste0("PreclinCancerProg_Exp_Rate_S2S3_", c("P", "D", "R")),
                                paste0("PreclinCancerProg_Exp_Rate_S3S4_", c("P", "D", "R"))),
  "Symptom detection" = c(paste0("pSxDetS1_", c("P", "D", "R")), paste0("hr_SxDetS2S1_", c("P", "D", "R")),
                          paste0("hr_SxDetS3S2_", c("P", "D", "R")), paste0("hr_SxDetS4S3_", c("P", "D", "R")))),
  "mini_priors_stage.png", "Progression and symptom detection priors")
cat("Stage prior panel done\n")

# ---- 10. Symptom detection rates as a dot plot ------------------------------
d_dot <- df_rate %>% mutate(loc = factor(c(P = "Proximal", D = "Distal", R = "Rectal")[loc], levels = c("Proximal", "Distal", "Rectal")),
                            who = factor(ifelse(set == "US", "US best set", "Chilean best set"), levels = c("US best set", "Chilean best set")),
                            years = 1 / r_det)
d_seg <- d_dot %>% select(loc, stage, who, r_det) %>% pivot_wider(names_from = who, values_from = r_det)
g_dot <- ggplot(d_dot, aes(x = factor(stage), y = r_det)) +
  geom_segment(data = d_seg, aes(x = factor(stage), xend = factor(stage), y = `US best set`, yend = `Chilean best set`),
               inherit.aes = FALSE, colour = "#b9b7b1", linewidth = 1.2) +
  geom_point(aes(colour = who), size = 4.2) +
  geom_text(aes(label = ifelse(years >= 10, sprintf("%.0f y", years), sprintf("%.1f y", years)), colour = who),
            nudge_x = ifelse(d_dot$who == "US best set", -0.28, 0.28), size = 3.6, show.legend = FALSE) +
  facet_wrap(~loc, nrow = 1) +
  scale_y_log10(breaks = c(0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5), labels = c("0.02", "0.05", "0.1", "0.2", "0.5", "1", "2", "5")) +
  scale_x_discrete(labels = paste("Stage", 1:4)) +
  scale_colour_manual(NULL, values = c("US best set" = col_us, "Chilean best set" = col_ch)) +
  labs(x = "stage at diagnosis", y = "annual symptom detection rate (log scale)",
       title = "Annual symptom detection rate: US vs Chile",
       subtitle = "Labels: average years until detection at that rate (1 / rate), ignoring progression and other-cause death.") +
  theme_bw(base_size = 13) +
  theme(legend.position = "top", legend.justification = "left", panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        strip.background = element_blank(), strip.text = element_text(hjust = 0, face = "bold"), plot.title = element_text(face = "bold"))
ggsave(file.path(dir_fig, "sx_detection_rates_dots.png"), g_dot, width = 13, height = 6, dpi = 200)
cat("Detection dot plot done\n")
