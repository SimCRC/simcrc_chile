# *****************************************************************************
#
# Script: plot_sx_competing_risk.R
#
# Purpose: Annual probability of symptomatic detection from each preclinical
#          CRC stage under competing risks, for the final Chile deck. A person
#          in preclinical stage k (Pk) is diagnosed at rate r_det, which moves
#          them to clinical stage k (Ck), or moves on to P(k + 1) at rate r_prog.
#          Method (Fernando): put these rates in one 8 x 8 rate matrix over P1
#          to P4 and C1 to C4, take its matrix exponential for one year, and
#          sum the probabilities of ending in C1 to C4 starting from Pk. This
#          counts a person who moves on and is diagnosed later in the same
#          year. Stage 4 has no progression exit.
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
#     Rscript analysis/plot_sx_competing_risk.R [dir_fig] [dir_tab]
#   Reads the Min_AbsolutErr set of the Chile calibration v0.15.0.20260921.1018
#   and of the US female adenoma model (simcrc::l_calibrated_params). Other
#   cause death is not an exit. Writes one figure (p_sx_us_vs_chile.png), one
#   table (sx_rates_probs_us_vs_chile.md) and the rates under the model's own
#   names (preclinical_rates_US_Chile.csv): PreclinCancerProg_Exp_Rate_SkS(k+1)
#   and rSxDetSk, as in simcrc R/05_simcrc_ssp_des.R.
#
# *****************************************************************************

suppressMessages({
  library(simcrc)
  library(ggplot2)
})

v_args  <- commandArgs(trailingOnly = TRUE)
dir_fig <- if (length(v_args) >= 1) v_args[1] else "SimCRC slides/final_deck/figures"
dir_tab <- if (length(v_args) >= 2) v_args[2] else "SimCRC slides/final_deck/tables"
dir.create(dir_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_tab, recursive = TRUE, showWarnings = FALSE)

cal_run <- "v0.15.0.20260921.1018"
dir_run <- file.path("outputs/BayCANN_versions/Chile/Adenoma/F/v0.15.0", cal_run)

e_cal <- new.env()
load(file.path(dir_run, sprintf(
  "l_params_calibrated_sets_SimCRC_%s_Adenoma_F.RData", cal_run)), envir = e_cal)
l_set <- list(`United States` = simcrc::l_calibrated_params$female$Min_AbsolutErr,
              `Chile` = e_cal$l_params_calibrated_sets$Min_AbsolutErr)
v_loc <- c(P = "Proximal colon", D = "Distal colon", R = "Rectum")

df_rate <- data.frame()
for (src in names(l_set)) {
  for (loc in names(v_loc)) {
    l_par <- l_set[[src]]
    p_det_s1 <- l_par[[paste0("pSxDetS1_", loc)]]
    v_hr_det <- c(l_par[[paste0("hr_SxDetS2S1_", loc)]],
                  l_par[[paste0("hr_SxDetS3S2_", loc)]],
                  l_par[[paste0("hr_SxDetS4S3_", loc)]])
    v_r_prog <- c(l_par[[paste0("PreclinCancerProg_Exp_Rate_S1S2_", loc)]],
                  l_par[[paste0("PreclinCancerProg_Exp_Rate_S2S3_", loc)]],
                  l_par[[paste0("PreclinCancerProg_Exp_Rate_S3S4_", loc)]], 0)
    v_r_det  <- cumprod(c(-log(1 - p_det_s1), v_hr_det))

    m_q <- matrix(0, nrow = 8, ncol = 8)
    m_q[cbind(1:4, 5:8)] <- v_r_det
    m_q[cbind(1:3, 2:4)] <- v_r_prog[1:3]
    diag(m_q) <- -rowSums(m_q)
    p_detect <- rowSums(as.matrix(Matrix::expm(m_q))[1:4, 5:8])

    df_rate <- rbind(df_rate, data.frame(source = src, loc = loc, stage = 1:4,
                                         r_det = v_r_det, r_prog = v_r_prog,
                                         p_detect = p_detect))
  }
}
df_rate$source <- factor(df_rate$source, levels = names(l_set))
df_rate$loc    <- factor(df_rate$loc, levels = names(v_loc))

df_us <- df_rate[df_rate$source == "United States", ]
df_ch <- df_rate[df_rate$source == "Chile", ]
stopifnot(all(df_us$loc == df_ch$loc), all(df_us$stage == df_ch$stage))

df_csv <- rbind(
  data.frame(parameter = paste0("PreclinCancerProg_Exp_Rate_S", df_us$stage, "S",
                                df_us$stage + 1, "_", df_us$loc),
             US = df_us$r_prog, Chile = df_ch$r_prog,
             stage = df_us$stage, loc = df_us$loc)[df_us$stage < 4, ],
  data.frame(parameter = paste0("rSxDetS", df_us$stage, "_", df_us$loc),
             US = df_us$r_det, Chile = df_ch$r_det,
             stage = df_us$stage + 4, loc = df_us$loc)
)
df_csv <- df_csv[order(df_csv$stage, df_csv$loc), c("parameter", "US", "Chile")]
write.csv(df_csv, file.path(dir_tab, "preclinical_rates_US_Chile.csv"),
          row.names = FALSE)

df_rate$above <- c(df_us$p_detect > df_ch$p_detect, df_ch$p_detect >= df_us$p_detect)

fig_us_ch <- ggplot(df_rate, aes(stage, p_detect, colour = source, group = source)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.8) +
  geom_text(aes(label = sprintf("%.3f", p_detect),
                vjust = ifelse(above, -0.9, 1.9)),
            size = 3.5, show.legend = FALSE) +
  facet_wrap(~ loc, labeller = labeller(loc = v_loc)) +
  scale_y_continuous(limits = c(-0.05, 1.12), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(breaks = 1:4, expand = expansion(add = 0.55)) +
  scale_colour_manual(values = c(`United States` = "#2a78d6",
                                 `Chile` = "#eb6834"), name = NULL) +
  labs(title = "Probability of diagnosis within a year: United States and Chile",
       subtitle = "Probability of diagnosis in any clinical stage (C1 to C4) within the year, same calculation for both. Stage 4 has no move-on exit.",
       x = "preclinical CRC stage", y = "annual probability of diagnosis") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top",
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0, face = "bold"),
        plot.title = element_text(face = "bold"))
ggsave(file.path(dir_fig, "p_sx_us_vs_chile.png"), fig_us_ch,
       width = 10, height = 4.6, dpi = 300, bg = "white")

v_loc_short <- c(P = "Proximal", D = "Distal", R = "Rectum")
writeLines(c(
  "| location | stage | diagnosis rate US | diagnosis rate Chile | Chile / US | move-on rate US | move-on rate Chile | P(diagnosis) US | P(diagnosis) Chile |",
  "| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
  sprintf("| %s | %d | %.3f | %.3f | %.2f | %s | %s | %.3f | %.3f |",
          v_loc_short[as.character(df_us$loc)], df_us$stage,
          df_us$r_det, df_ch$r_det, df_ch$r_det / df_us$r_det,
          ifelse(df_us$stage == 4, "-", sprintf("%.3f", df_us$r_prog)),
          ifelse(df_us$stage == 4, "-", sprintf("%.3f", df_ch$r_prog)),
          df_us$p_detect, df_ch$p_detect)
), file.path(dir_tab, "sx_rates_probs_us_vs_chile.md"))

cat("Outputs written to", dir_fig, "and", dir_tab, "\n")
