# Fernando's no-recalibration control, decomposed.
#
# Three arms, each 68 or 60 strategies at n_pop = 1e6:
#   A  old code, old params   v0.14.0.2 run, committed at 58389af
#   B  new code, old params   the control: v0.14.0.2 parameters under v0.15.0
#   C  new code, new params   the current v0.15.0 run
#
# B - A isolates the model code change, C - B isolates the prior retune, and
# C - A is the headline difference the deck currently reports as one number.
#
# A and C are read from the deck build cache so only B is post-processed here.
#
# Run from the repo root:  Rscript analysis/compare_code_vs_priors.R

suppressMessages({ library(dplyr); library(dampack) })

CACHE <- paste0("/Users/jorgeroa/Documents/GitHub/presentations/simcrc/chile",
                "/resources/tables/.cea_frames.rds")
CTL_CACHE <- "output/2026Chile/ControlAugParams/.control_frame.rds"
WTP   <- 16e6
TAB   <- paste0("/Users/jorgeroa/Documents/GitHub/presentations/simcrc/chile",
                "/resources/tables")

process_folder <- function(analysis_folder) {
  for (f in c("R/crc_allocation_time.R", "R/uspstf_summary.R",
              "R/process_uspstf_output.R")) source(f)
  out <- ProcessUSPSTFOutput(
    analysis_folder = analysis_folder,
    input_folder = "data-raw/cea_inputs", input_prefix = "2026Chile",
    first_age_of_interest = 40, col_infl_rate = 1.05, discount_rate = 0.03,
    col_spec_adj = 0.86,
    crc_care_costs_file = "crc_care_costs.csv",
    screen_costs_file = "screen_costs_v3.csv",
    crc_care_disutility_file = "crc_care_utility_loss.csv",
    screen_disutility_file = "screen_utility_loss_WithStoolTestValues.csv",
    general_health_utility_weights_file = "GeneralHealthUtilityWeightsByAge.csv",
    selected_outcomes_for_model_data_file =
      "model_data_outcomes_boolean_addDscQALY.csv",
    model_run_data_tag = "_ctl", folder_for_output = "ce_results_control",
    include_SimCRC = FALSE, include_SimCRC_R = TRUE,
    include_MISCAN = FALSE, include_CRCSPIN = FALSE,
    output_template_year = 2028)
  out$Strategy <- gsub("2026Chile_", "", out$Strategy)
  out
}

cf <- readRDS(CACHE)
A  <- cf$old   # old code, old params
C  <- cf$new   # new code, new params

if (file.exists(CTL_CACHE)) {
  message("reusing cached control frame")
  B <- readRDS(CTL_CACHE)
} else {
  message("post-processing the control arm ...")
  B <- process_folder("output/2026Chile/ControlAugParams")
  saveRDS(B, CTL_CACHE)
}
message(sprintf("arms: A=%d  B=%d  C=%d strategies", nrow(A), nrow(B), nrow(C)))

shared <- Reduce(intersect, list(A$Strategy, B$Strategy, C$Strategy))
message("shared strategies: ", length(shared))
pick <- function(d) d[match(shared, d$Strategy), ]
A <- pick(A); B <- pick(B); C <- pick(C)

md_table <- function(df, path, align = NULL) {
  if (is.null(align)) align <- c("l", rep("r", ncol(df) - 1))
  sep <- vapply(align, function(a)
    switch(a, l = ":---", r = "---:", "---"), character(1))
  writeLines(c(paste0("| ", paste(names(df), collapse = " | "), " |"),
               paste0("| ", paste(sep, collapse = " | "), " |"),
               apply(df, 1, function(r)
                 paste0("| ", paste(r, collapse = " | "), " |"))), path)
}

# ---- 1. mean outcome per arm, and the two effects ---------------------------
metrics <- list(
  list(v = "CRCCasesper1000",            lab = "CRC cases per 1,000",
       f = function(x) sprintf("%.1f", x)),
  list(v = "CRCandComplDeathsper1000",   lab = "CRC deaths per 1,000",
       f = function(x) sprintf("%.2f", x)),
  list(v = "DiscountedQALYGainedper1000", lab = "QALYs gained, discounted",
       f = function(x) sprintf("%.1f", x)),
  list(v = "CostsCancerCareper1000",     lab = "cancer-care cost per 1,000",
       f = function(x) sprintf("%.0fM", x / 1e6)),
  list(v = "DiscountedTotalCostsper1000", lab = "discounted cost per 1,000",
       f = function(x) sprintf("%.0fM", x / 1e6)))
pct <- function(x, y) sprintf("%+.1f %%", 100 * (y - x) / x)
rows <- lapply(metrics, function(mm) {
  m <- mm$v
  a <- mean(A[[m]]); b <- mean(B[[m]]); cc <- mean(C[[m]])
  data.frame(outcome = mm$lab,
             `A old / old` = mm$f(a),
             `B new / old` = mm$f(b),
             `C new / new` = mm$f(cc),
             `code effect` = pct(a, b),
             `prior effect` = pct(b, cc),
             check.names = FALSE)
})
tb <- do.call(rbind, rows)
print(tb, row.names = FALSE)
md_table(tb, file.path(TAB, "control_decomposition.md"),
         align = c("l", "r", "r", "r", "r", "r"))

# ---- 2. what each arm would recommend at the threshold ----------------------
decide <- function(d, label) {
  ic <- calculate_icers(cost = d$DiscountedTotalCostsper1000,
                        effect = d$DiscountedQALYGainedper1000,
                        strategies = d$Strategy)
  f <- ic[ic$Status == "ND", ]
  f <- f[order(f$Effect), ]
  v <- f$ICER; v[is.na(v)] <- 0
  i <- max(which(cumsum(v > WTP) == 0))
  nxt <- if (i < nrow(f)) sprintf("%.1fM", v[i + 1] / 1e6) else "-"
  data.frame(arm = label, `frontier size` = nrow(f),
             `chosen at 16M` = {
               q <- strsplit(as.character(f$Strategy[i]), "_")[[1]]
               if (length(q) == 4L) sprintf("%s %s-%s q%s", q[1], q[2], q[3], q[4])
               else as.character(f$Strategy[i])
             },
             `its ICER` = sprintf("%.1fM", v[i] / 1e6),
             `next step` = nxt, check.names = FALSE)
}
dec <- rbind(decide(A, "A  old code, old priors"),
             decide(B, "B  new code, old priors"),
             decide(C, "C  new code, new priors"))
cat("\n")
print(dec, row.names = FALSE)
md_table(dec, file.path(TAB, "control_decision.md"),
         align = c("l", "r", "l", "r", "r"))

# ---- 3. the same two results, said plainly -----------------------------------
# The full tables above are for the appendix. These two are what goes on a
# slide: one says the recommendation only moves when the priors move, the
# other says every outcome agrees.

rec <- dec$`chosen at 16M`
verdict <- data.frame(
  arm = c("A &nbsp; August baseline", "B &nbsp; **the control**",
          "C &nbsp; September"),
  `what changed` = c("&mdash;", "**new model code**",
                     "new model code **+ new priors**"),
  `recommends at 16M` = c(rec[1],
                          paste0(rec[2], " &nbsp; *unchanged*"),
                          paste0("**", rec[3], "** &nbsp; *flips*")),
  check.names = FALSE)
md_table(verdict, file.path(TAB, "control_verdict.md"),
         align = c("l", "l", "l"))

effects <- data.frame(
  outcome = vapply(metrics, function(m) m$lab, ""),
  `changing the code` = vapply(metrics, function(m)
    pct(mean(A[[m$v]]), mean(B[[m$v]])), ""),
  `changing the priors` = vapply(metrics, function(m)
    pct(mean(B[[m$v]]), mean(C[[m$v]])), ""),
  check.names = FALSE)
md_table(effects, file.path(TAB, "control_effects.md"),
         align = c("l", "r", "r"))

cat("\nwrote control_decomposition.md, control_decision.md,",
    "control_verdict.md, control_effects.md\n")
