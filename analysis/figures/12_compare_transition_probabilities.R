# Build the old/new USA comparison used by the focused transition slide deck.
# Run from the repository root:
#   Rscript analysis/figures/12_compare_transition_probabilities.R
# The original 8-state probability calculation stays in plot_sx_competing_risk.R.

suppressPackageStartupMessages(library(ggplot2))
out <- "outputs/comparison/USA_20261026_Chile_20260921/transition_probabilities"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
usa_file <- "data-raw/l_params_calibrated_sets_SimCRC_v0.15.0.20261026.130537_Adenoma_F.RData"
script <- "analysis/figures/plot_sx_competing_risk.R"

for (version in c("previous", "updated")) {
  args <- c(script, file.path(out, version, "figures"), file.path(out, version, "tables"))
  if (version == "updated") args <- c(args, usa_file)
  status <- system2(file.path(R.home("bin"), "Rscript"), vapply(args, shQuote, character(1)))
  if (status != 0) stop("Probability calculation failed for ", version)
}

previous <- read.csv(file.path(out, "previous/tables/diagnosis_probabilities_US_Chile.csv"))
updated <- read.csv(file.path(out, "updated/tables/diagnosis_probabilities_US_Chile.csv"))
stopifnot(identical(previous[c("source", "loc", "stage")], updated[c("source", "loc", "stage")]))
chile <- updated[updated$source == "Chile", ]
stopifnot(identical(previous[previous$source == "Chile", ], chile))
us_old <- previous[previous$source == "United States", ]
us_new <- updated[updated$source == "United States", ]
comparison <- data.frame(
  location = c(P = "Proximal", D = "Distal", R = "Rectum")[us_new$loc],
  stage = us_new$stage, USA_previous = us_old$p_detect,
  USA_updated = us_new$p_detect, Chile = chile$p_detect,
  USA_change_pp = 100 * (us_new$p_detect - us_old$p_detect)
)
write.csv(comparison, file.path(out, "diagnosis_probability_comparison.csv"), row.names = FALSE)

fig_delta <- ggplot(comparison, aes(stage, USA_change_pp)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = .4) +
  geom_col(fill = "#2a78d6", width = .55) +
  geom_text(aes(label = sprintf("%+.2f", USA_change_pp),
                vjust = ifelse(USA_change_pp >= 0, -.5, 1.4)), size = 4) +
  facet_wrap(~factor(location, levels = c("Proximal", "Distal", "Rectum"))) +
  scale_x_continuous(breaks = 1:4) +
  scale_y_continuous(limits = c(-27, 17), breaks = seq(-20, 10, 10)) +
  labs(title = "Change in USA diagnosis probability",
       subtitle = "Updated October 26 set minus the previous USA set. Chile is unchanged.",
       x = "Starting undiagnosed cancer stage", y = "Percentage points") +
  theme_bw(base_size = 13) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
ggsave(file.path(out, "USA_probability_change.png"), fig_delta,
       width = 10, height = 4.6, dpi = 300, bg = "white", device = ragg::agg_png)

md_table <- function(d, path) {
  lines <- c(paste0("| ", paste(names(d), collapse = " | "), " |"),
             paste0("| ", paste(rep("---", ncol(d)), collapse = " | "), " |"),
             apply(d, 1, function(row) paste0("| ", paste(row, collapse = " | "), " |")))
  writeLines(lines, path)
}
display <- data.frame(
  Site = comparison$location, Stage = comparison$stage,
  `USA previous` = sprintf("%.1f%%", 100 * comparison$USA_previous),
  `USA updated` = sprintf("%.1f%%", 100 * comparison$USA_updated),
  Chile = sprintf("%.1f%%", 100 * comparison$Chile),
  `USA change (pp)` = sprintf("%+.2f", comparison$USA_change_pp), check.names = FALSE
)
md_table(display, file.path(out, "diagnosis_probability_comparison.md"))

# A concrete example of the joint effect of changing both exits from stage 3.
rectal_old <- us_old[us_old$loc == "R" & us_old$stage == 3, ]
rectal_new <- us_new[us_new$loc == "R" & us_new$stage == 3, ]
rectal_chile <- chile[chile$loc == "R" & chile$stage == 3, ]
mechanism <- data.frame(
  Quantity = c("Diagnosis rate, per year", "Progression rate, per year", "Diagnosis within one year"),
  `USA previous` = c(sprintf("%.3f", rectal_old$r_det), sprintf("%.3f", rectal_old$r_prog),
                     sprintf("%.1f%%", 100 * rectal_old$p_detect)),
  `USA updated` = c(sprintf("%.3f", rectal_new$r_det), sprintf("%.3f", rectal_new$r_prog),
                    sprintf("%.1f%%", 100 * rectal_new$p_detect)),
  Chile = c(sprintf("%.3f", rectal_chile$r_det), sprintf("%.3f", rectal_chile$r_prog),
            sprintf("%.1f%%", 100 * rectal_chile$p_detect)), check.names = FALSE
)
md_table(mechanism, file.path(out, "rectal_stage_3_example.md"))

writeLines(c(
  "USA previous: packaged SimCRC female v0.13.0.1 Min_AbsolutErr.",
  paste("USA updated:", usa_file),
  "Chile: v0.15.0.20260921.1018 Min_AbsolutErr, identical in both comparisons.",
  "Probability definition: diagnosis at any clinical stage within one year, conditional on the initial undiagnosed stage.",
  "Progression followed by diagnosis within the year is included. Other-cause mortality and screening are excluded.",
  "All probabilities use the same original 8-state rate-matrix calculation.",
  "Percentage-point change = 100 times (updated probability minus previous probability).",
  "Selected parameter sets only; calibration uncertainty is not propagated.",
  "Stage-4 values rounded to 100.0% are near certainty rather than exactly one."
), file.path(out, "comparison_details.txt"))
print(comparison)
message("Saved transition comparison to ", out)
