# Build readable slide assets from the completed, verified comparisons.
# Run from the repo root: Rscript analysis/figures/13_calibration_update_slides.R
# No simulations or calibration are performed here.
suppressPackageStartupMessages({ library(ggplot2); library(data.table) })
root <- "outputs/comparison/USA_20261026_Chile_20260921"
out <- file.path(root, "update_slides")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
colours <- c(Chile = "#F94144", USA = "#F8961E")
stats <- read.csv(file.path(root, "cheat_sheet/natural_history_summary.csv"), check.names = FALSE)
phases <- read.csv(file.path(root, "cheat_sheet/phase_duration_means.csv"))
stopifnot(identical(phases$Phase, c("Small adenoma (<6 mm)", "Medium adenoma (6-10 mm)",
                                  "Large adenoma before cancer", "Cancer before diagnosis (sojourn)")))
hist <- as.data.frame(readRDS(file.path(root, "cheat_sheet/natural_history_plot_data.rds")))
stopifnot(setequal(as.character(unique(hist$type)), c("R v0.15.0-C MAE (Chile)", "R v0.15.0-C MAE (USA)")))
hist$country <- factor(ifelse(grepl("(Chile)", hist$type, fixed = TRUE), "Chile", "USA"), levels = names(colours))
get_stat <- function(stage, country, field = "Mean") {
  stats[[sprintf("%s (R_R v0.15.0-C MAE (%s))", field, country)]][match(stage, stats$stage)]
}
counts <- c(Chile = get_stat("dwell_time", "Chile", "N"), USA = get_stat("dwell_time", "USA", "N"))
for (country in names(colours)) {
  stopifnot(all(vapply(c("dw_t_sm", "dw_t_md", "dw_t_lg", "sojourn_time"),
                      get_stat, numeric(1), country = country, field = "N") == counts[[country]]))
  stopifnot(abs(sum(phases[[country]]) - get_stat("dwell_time", country)) < .025)
}

durations <- hist[hist$stage %in% c("dwell_time", "sojourn_time"), ]
# The saved age field contains one-year duration bins for these two stages.
durations$stage <- factor(durations$stage, levels = c("dwell_time", "sojourn_time"))
labels <- c(
  dwell_time = sprintf("Total dwell time\nChile %.2f years | USA %.2f years", get_stat("dwell_time", "Chile"), get_stat("dwell_time", "USA")),
  sojourn_time = sprintf("Sojourn time\nChile %.2f years | USA %.2f years", get_stat("sojourn_time", "Chile"), get_stat("sojourn_time", "USA")))
p <- ggplot(durations, aes(age, n_ade / 1000, colour = country)) +
  geom_line(linewidth = 1) + facet_wrap(~stage, scales = "free", labeller = labeller(stage = labels)) +
  scale_colour_manual(values = colours, name = NULL) +
  labs(x = "Duration (years)", y = "Diagnosed cases (thousands)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "top", strip.background = element_blank(),
        strip.text = element_text(face = "bold"), panel.grid.minor = element_blank())
ggsave(file.path(out, "dwell_sojourn.png"), p, width = 10, height = 4.6,
       dpi = 300, bg = "white", device = ragg::agg_png)

# Preserve the probability calculation and change only the slide presentation.
for (version in c("previous", "updated")) {
  d <- read.csv(file.path(root, "transition_probabilities", version,
                          "tables/diagnosis_probabilities_US_Chile.csv"))
  us <- d[d$source == "United States", ]; ch <- d[d$source == "Chile", ]
  stopifnot(all(us$loc == ch$loc), all(us$stage == ch$stage))
  d$above <- c(us$p_detect > ch$p_detect, ch$p_detect >= us$p_detect)
  d$country <- factor(ifelse(d$source == "Chile", "Chile", "USA"), levels = names(colours))
  d$loc <- factor(d$loc, levels = c("P", "D", "R"))
  p <- ggplot(d, aes(stage, p_detect, colour = country, group = country)) +
    geom_line(linewidth = 1) + geom_point(size = 3) +
    geom_text(aes(label = sprintf("%.3f", p_detect), vjust = ifelse(above, -.9, 1.9)),
              size = 4, show.legend = FALSE) +
    facet_wrap(~loc, labeller = labeller(loc = c(P = "Proximal colon", D = "Distal colon", R = "Rectum"))) +
    scale_y_continuous(limits = c(-.05, 1.12), breaks = seq(0, 1, .2)) +
    scale_x_continuous(breaks = 1:4, expand = expansion(add = .55)) +
    scale_colour_manual(values = colours, name = NULL) +
    labs(x = "Starting undiagnosed cancer stage", y = "Probability of diagnosis within one year") +
    theme_bw(base_size = 14) +
    theme(legend.position = "top", strip.background = element_blank(),
          strip.text = element_text(face = "bold"), panel.grid.minor = element_blank())
  ggsave(file.path(out, paste0("diagnosis_", version, ".png")), p,
         width = 10, height = 4.6, dpi = 300, bg = "white", device = ragg::agg_png)
}

d <- read.csv(file.path(root, "transition_probabilities/diagnosis_probability_comparison.csv"))
p <- ggplot(d, aes(stage, USA_change_pp)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = .4) +
  geom_col(fill = colours[["USA"]], width = .55) +
  geom_text(aes(label = sprintf("%+.2f", USA_change_pp), vjust = ifelse(USA_change_pp >= 0, -.5, 1.4)), size = 4) +
  facet_wrap(~factor(location, levels = c("Proximal", "Distal", "Rectum"))) +
  scale_x_continuous(breaks = 1:4) +
  scale_y_continuous(limits = c(-27, 17), breaks = seq(-20, 10, 10)) +
  labs(x = "Starting undiagnosed cancer stage", y = "USA change (percentage points)") +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())
ggsave(file.path(out, "diagnosis_change.png"), p, width = 10, height = 4.6,
       dpi = 300, bg = "white", device = ragg::agg_png)

table <- data.frame(
  Phase = c("Small adenoma", "Medium adenoma", "Large adenoma", "Sojourn", "Total dwell"),
  Interval = c("Polyp onset to 6 mm", "6 mm to 10 mm", "10 mm to cancer onset", "Cancer onset to diagnosis", "Polyp onset to diagnosis"),
  `Chile, years` = sprintf("%.2f", c(phases$Chile, get_stat("dwell_time", "Chile"))),
  `USA, years` = sprintf("%.2f", c(phases$USA, get_stat("dwell_time", "USA"))), check.names = FALSE)
writeLines(c(
  paste0("| ", paste(names(table), collapse = " | "), " |"),
  paste0("| ", paste(rep("---", ncol(table)), collapse = " | "), " |"),
  apply(table, 1, function(row) paste0("| ", paste(row, collapse = " | "), " |"))
), file.path(out, "phase_duration_table.md"))
write.csv(table, file.path(out, "phase_duration_table.csv"), row.names = FALSE)
cat("Saved slide assets to", out, "\n")
