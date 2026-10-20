# *****************************************************************************
#
# Script: plot_frontier_claudia_style.R
#
# Purpose: Figure 2 (Chile CE frontiers: all strategies, FIT only) and
#          Figure 3 (lifetime colonoscopies vs QALYs against capacity) in the
#          Seguin et al./SDA iconography with ggpop markers. Shape = ages to
#          begin and end screening; colour = modality and interval (Claudia
#          Seguin, 23 Sep 2026). Legend positions are in grid cells (one
#          row = 1) and text sizes in mm.
#
# Author: Jorge Roa
#
# Email: jorgeroa@stanford.edu
#
# ISS: #25   Version: #24
#
# Date Created: 23 September 2026
#
# *****************************************************************************

library(ggplot2)
library(ggpop)
library(patchwork)
library(scales)

n_wtp <- 16e6
# PLACEHOLDER lifetime colonoscopies per 1,000 until the infrastructure study
n_cap <- 1500

# sda2028 SetShapesColorsAndFills(): start age = outline, stop age = fill
v_outline <- c(`45` = "square", `50` = "circle", `55` = "diamond")
v_fill    <- c(`70` = "inset", `75` = "hollow", `80` = "cross", `85` = "solid")
v_once    <- c(`45` = "square-solid", `50` = "circle-solid",
               `55` = "diamond-solid", `60` = "plus-bold",
               `65` = "triangle-down")

# Dataviz validator, all pairs on white: worst dE 18.8 normal, 10.3 colour-blind
v_group <- c(`One-time COL`       = "#08B8A2", `COL every 15 years` = "#1090F3",
             `COL every 10 years` = "#2D3CF0", `COL every 5 years`  = "#574280",
             `FIT every 3 years`  = "#C99D0F", `FIT every 2 years`  = "#DF4601",
             `FIT every year`     = "#980124", `No screening`       = "black")

df_ce <- rbind(
  cbind(analysis = "all", read.csv("ce_results/df_icer_all_strategies.csv")),
  cbind(analysis = "fit", read.csv("ce_results/df_icer_FIT.csv"))
)
# Strategy is MODALITY_start_stop_interval (interval 100 = once) or NoScreening
df_ce <- cbind(df_ce, read.table(text = df_ce$Strategy, sep = "_", fill = TRUE,
                                 colClasses = "character",
                                 col.names = c("modality", "start", "stop",
                                               "interval")))
df_ce <- within(df_ce, {
  none   <- modality == "NoScreening"
  once   <- interval == "100"
  cost_m <- Cost / 1e6
  group  <- ifelse(none, "No screening",
            ifelse(once, "One-time COL",
            ifelse(interval == "1", "FIT every year",
                   paste(modality, "every", interval, "years"))))
  icon   <- ifelse(none, "circle-solid",
            ifelse(once, v_once[start],
                   paste0(v_outline[start], "-", v_fill[stop])))
  label  <- ifelse(none, "No screening",
            ifelse(once, paste0("COL ", start, ", once"),
                   paste0(modality, " ", start, "-", stop, ", ", interval)))
  tag    <- ifelse(is.na(ICER), label, sprintf("%s (%.1fM)", label, ICER / 1e6))
  ages   <- ifelse(once, paste0(start, "y"), paste0(start, "-", stop, "y"))
  below_wtp <- Status == "ND" & ICER <= n_wtp
})

df_all <- df_ce[df_ce$analysis == "all", ]
df_fit <- df_ce[df_ce$analysis == "fit", ]
# Optimal: the most QALYs among efficient strategies with an ICER under WTP
df_all$optimal <- df_all$Effect ==
  max(df_all$Effect[df_all$below_wtp], na.rm = TRUE)
df_fit$optimal <- df_fit$Effect ==
  max(df_fit$Effect[df_fit$below_wtp], na.rm = TRUE)

plot_frontier <- function(df) {
  df_nd <- df[df$Status == "ND", ]
  ggplot(df, aes(Effect, cost_m)) +
    geom_line(data = df_nd, linewidth = 0.6) +
    geom_icon_point(aes(icon = icon, colour = group), size = 1.4, dpi = 400,
                    show.legend = FALSE) +
    geom_text(aes(label = tag, fontface = ifelse(optimal, "bold", "plain")),
              data = df_nd, hjust = 0, vjust = 1, size = 3.1,
              nudge_x = 0.012 * diff(range(df$Effect)),
              nudge_y = -0.018 * diff(range(df$cost_m))) +
    scale_colour_manual(values = v_group) +
    scale_x_continuous(breaks = breaks_pretty(8),
                       expand = expansion(mult = c(0.03, 0.16))) +
    scale_y_continuous(labels = label_number(big.mark = ",", suffix = "M"),
                       breaks = breaks_pretty(8)) +
    labs(x = "Discounted QALYs gained per 1,000",
         y = "Discounted total costs per 1,000 (CLP millions)") +
    theme_classic(base_size = 13) +
    theme(axis.text = element_text(colour = "black"))
}

gg_all <- plot_frontier(df_all)
gg_fit <- plot_frontier(df_fit)

# icon_grid() lays out the plotted shapes; once-only COL goes in row 4
df_leg_grid <- rbind(
  icon_grid(df_fit[df_fit$modality == "FIT", ], icon = "icon", label = "ages",
            row = "start", col = "stop"),
  transform(icon_grid(df_all[df_all$once, ], icon = "icon", label = "ages",
                      row = "modality", col = "start"), row = 4)
)
# The white key makes "Last row: once" a text-only row
df_leg_symbol <- data.frame(section = "symbol", type = c("line", "point"),
                            label = c("Efficient frontier", "Last row: once"),
                            color = c("black", "white"),
                            icon = NA, row = NA, col = NA)
df_key_col <- data.frame(type = "swatch", label = names(v_group)[1:4],
                         color = unname(v_group)[1:4])
df_key_fit <- data.frame(type = "swatch", label = names(v_group)[5:7],
                         color = unname(v_group)[5:7])

# Colour key: COL and FIT columns left of the shape grid, which starts at x = 0
leg_all <- legend_canvas(
  rbind(df_leg_grid, df_leg_symbol),
  grid_title = "Age to begin-age to end screening", symbol_section = "symbol",
  col_spacing = 0.986, row_spacing = 1, label_gap = 0.167,
  symbol_right_gap = 0.875, symbol_key_width = 0.567, symbol_label_gap = 0.06,
  marker_size = 3.717, label_size = 3.634, dpi = 400, x_margin = c(5.106, 2.137)
) +
  annotate("text", x = -2.533, y = 0.85, label = "Modality and interval",
           size = 3.634) +
  key_legend(df_key_col, x = -4.512, y_start = 0, key_width = 0.2375,
             label_gap = 0.1187, label_size = 3.634, swatch_height = 0.56) +
  key_legend(df_key_fit, x = -2.296, y_start = 0, key_width = 0.2375,
             label_gap = 0.1187, label_size = 3.634, swatch_height = 0.56) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))
# legend_box() measures the canvas unexpanded, hence the zero expand above
leg_all <- legend_box(leg_all, width = 10.06, height = 1.52,
                      padding = c(0.018, 0.09), colour = "#231F20")

leg_fit <- legend_canvas(
  rbind(df_leg_grid[df_leg_grid$row < 4, ], df_leg_symbol[1, ]),
  grid_title = "Age to begin-age to end screening", symbol_section = "symbol",
  col_spacing = 0.986, row_spacing = 1, label_gap = 0.167,
  symbol_right_gap = 0.875, symbol_key_width = 0.567, symbol_label_gap = 0.06,
  marker_size = 3.717, label_size = 3.634, dpi = 400, x_margin = c(4.077, 3.167)
) +
  annotate("text", x = -1.583, y = 0.85, label = "Modality and interval",
           size = 3.634) +
  key_legend(df_key_fit, x = -2.454, y_start = 0, key_width = 0.2375,
             label_gap = 0.1187, label_size = 3.634, swatch_height = 0.56) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))
leg_fit <- legend_box(leg_fit, width = 10.06, height = 1.27,
                      padding = c(0.018, 0.09), colour = "#231F20")

fig_all <- gg_all / leg_all + plot_layout(heights = c(1, 0.24))
ggsave("ce_results/frontier_all_claudia.png", fig_all,
       width = 11, height = 8.6, dpi = 300, bg = "white")

fig_fit <- gg_fit / leg_fit + plot_layout(heights = c(1, 0.2))
ggsave("ce_results/frontier_fit_claudia.png", fig_fit,
       width = 11, height = 8.3, dpi = 300, bg = "white")

fig2 <- gg_all / leg_all / gg_fit / leg_fit +
  plot_layout(heights = c(1, 0.24, 1, 0.2)) +
  plot_annotation(tag_levels = list(c("A", "", "B", "")))
ggsave("ce_results/fig2_frontiers_claudia.png", fig2,
       width = 11, height = 16.5, dpi = 300, bg = "white")

df_run <- readxl::read_excel(
  file.path("output/2026Chile/Base/ce_results",
            "2026-09-21 230512.351928_model_run_data_Base.xlsx"),
  sheet = "model_run_data"
)
# Capacity is net of symptomatic diagnosis (Appendix), so demand excludes it
df_run$cols <- df_run$ScreeningColsper1000 + df_run$DiagnosticFUper1000 +
  df_run$SurveillanceColsper1000

df_cap <- merge(df_all[!df_all$none, ], df_run[c("Strategy", "cols")])
df_cap$optimal <- df_cap$optimal |
  df_cap$Strategy %in% df_fit$Strategy[df_fit$optimal]
df_cap <- df_cap[order(df_cap$cols), ]
# On the frontier when no strategy using fewer colonoscopies gains more QALYs
df_cap$on_frontier <- df_cap$Effect == cummax(df_cap$Effect)

gg_cap <- ggplot(df_cap, aes(Effect, cols)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = n_cap,
           fill = "grey93") +
  geom_hline(yintercept = n_cap, linetype = "dashed", colour = "grey30",
             linewidth = 0.6) +
  geom_line(data = df_cap[df_cap$on_frontier, ], linewidth = 0.6) +
  geom_icon_point(aes(icon = icon, colour = group), size = 1.4, dpi = 400,
                  show.legend = FALSE) +
  ggrepel::geom_text_repel(
    aes(label = label, fontface = ifelse(optimal, "bold", "plain")),
    data = df_cap[df_cap$on_frontier & df_cap$modality == "FIT" &
                    df_cap$interval != "1", ],
    size = 3.1, nudge_y = -420, direction = "both", force = 2,
    ylim = c(120, NA), min.segment.length = 0, segment.colour = "grey55",
    box.padding = 0.3, max.overlaps = Inf, seed = 1
  ) +
  ggrepel::geom_text_repel(
    aes(label = label, fontface = ifelse(optimal, "bold", "plain")),
    data = df_cap[df_cap$on_frontier & (df_cap$modality == "COL" |
                                          df_cap$interval == "1"), ],
    size = 3.1, nudge_x = 1.6, direction = "y", hjust = 0,
    min.segment.length = 0, segment.colour = "grey55", box.padding = 0.3,
    max.overlaps = Inf, seed = 1
  ) +
  annotate("text", x = min(df_cap$Effect), y = n_cap - 60, vjust = 1,
           hjust = 0, size = 3.1, colour = "grey30", lineheight = 0.9,
           label = "National colonoscopy\ncapacity (illustrative)") +
  annotate("text", x = min(df_cap$Effect), y = n_cap / 2, hjust = 0,
           size = 3.4, fontface = "italic", colour = "grey30",
           label = "Deliverable") +
  scale_colour_manual(values = v_group) +
  scale_x_continuous(breaks = breaks_pretty(8),
                     expand = expansion(mult = c(0.03, 0.08))) +
  scale_y_continuous(labels = label_number(big.mark = ","),
                     breaks = breaks_pretty(8),
                     expand = expansion(mult = c(0, 0.04))) +
  labs(x = "Discounted QALYs gained per 1,000",
       y = "Lifetime colonoscopies per 1,000") +
  theme_classic(base_size = 13) +
  theme(axis.text = element_text(colour = "black"))

# National capacity takes the third row of the symbol column (x = 4.819)
leg_cap <- leg_all +
  annotate("segment", x = 4.819, xend = 5.386, y = -1, yend = -1,
           linetype = "dashed", colour = "grey30", linewidth = 0.6) +
  annotate("text", x = 5.446, y = -1, label = "National capacity", hjust = 0,
           size = 3.634)

fig3 <- gg_cap / leg_cap + plot_layout(heights = c(1, 0.24))
ggsave("ce_results/fig3_capacity_frontier.png", fig3,
       width = 11, height = 8.6, dpi = 300, bg = "white")
