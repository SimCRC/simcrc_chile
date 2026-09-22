# Cost-effectiveness plane drawn with ggpop, replacing dampack:::plot.icers().
#
# Each strategy is a Font Awesome icon rather than a dot, so modality is carried by
# shape as well as hue and survives greyscale printing and colour-vision deficiency.
# Dominated strategies stay visible but recessive; the frontier is connected and
# labelled; the strategy that is optimal at a given willingness to pay is ringed.
#
# The WTP mark is deliberately the optimal frontier point and NOT a ray from the
# no-screening corner: a ray reads as AVERAGE cost-effectiveness, under which every
# Chile strategy clears 16M per QALY while incremental ICERs run past 400M.

plot_ce_frontier_pop <- function(icers,
                                 wtp = NULL,
                                 title = NULL,
                                 subtitle = NULL,
                                 icons = c(Colonoscopy = "stethoscope",
                                           FIT = "vial",
                                           `No screening` = "ban"),
                                 label_frontier = TRUE) {
  stopifnot(all(c("Strategy", "Cost", "Effect", "Status") %in% names(icers)))

  d <- data.frame(
    Strategy = as.character(icers$Strategy),
    cost     = as.numeric(icers$Cost) / 1e6,
    effect   = as.numeric(icers$Effect),
    status   = as.character(icers$Status),
    stringsAsFactors = FALSE
  )
  d$modality <- factor(
    ifelse(grepl("NoScreening", d$Strategy), "No screening",
           ifelse(grepl("^FIT", d$Strategy), "FIT", "Colonoscopy")),
    levels = names(icons)
  )
  d$efficient <- d$status == "ND"

  # MOD_start_end_interval -> "COL 45-70 q15". Frontier labels are the only text
  # in the panel, so they carry the age range and interval rather than a code the
  # reader has to decode.
  d$label <- vapply(d$Strategy, function(s) {
    if (grepl("NoScreening", s)) return("No screening")
    p <- strsplit(s, "_")[[1]]
    if (length(p) == 4L) sprintf("%s %s-%s q%s", p[1], p[2], p[3], p[4]) else s
  }, character(1), USE.NAMES = FALSE)

  frontier <- d[d$efficient, ]
  frontier <- frontier[order(frontier$effect), ]

  optimal <- NULL
  if (!is.null(wtp) && nrow(frontier) > 0) {
    icer_f <- suppressWarnings(as.numeric(icers$ICER[match(frontier$Strategy,
                                                           icers$Strategy)]))
    icer_f[is.na(icer_f)] <- 0
    ok <- which(cumsum(icer_f > wtp) == 0)
    if (length(ok)) optimal <- frontier[max(ok), ]
  }

  pal <- c(Colonoscopy = "#2a78d6", FIT = "#eb6834", `No screening` = "#1baf7a")
  present <- names(pal)[names(pal) %in% as.character(unique(d$modality))]

  p <- ggplot2::ggplot(d, ggplot2::aes(x = effect, y = cost)) +
    ggplot2::geom_line(data = frontier, colour = "#a01b1b", linewidth = 0.7)

  # one layer per modality: geom_icon_point() takes a single icon per layer, and
  # splitting this way keeps each modality's icon and colour locked together
  for (m in levels(d$modality)) {
    dom <- d[d$modality == m & !d$efficient, ]
    eff <- d[d$modality == m & d$efficient, ]
    if (nrow(dom)) {
      p <- p + ggpop::geom_icon_point(data = dom, icon = unname(icons[[m]]),
                                      colour = pal[[m]], size = 1.3,
                                      alpha = 0.40, show.legend = FALSE)
    }
    if (nrow(eff)) {
      p <- p + ggpop::geom_icon_point(data = eff, icon = unname(icons[[m]]),
                                      colour = pal[[m]], size = 2.3,
                                      show.legend = FALSE)
    }
  }

  # the WTP ring goes on last so it reads above the icons rather than behind them
  if (!is.null(optimal)) {
    p <- p + ggplot2::geom_point(data = optimal, shape = 21, size = 11,
                                 stroke = 1.1, colour = "#26261f", fill = NA)
  }

  # geom_icon_point() emits no legend key, so identity comes from an invisible
  # (alpha 0) dummy layer whose key is forced opaque via override.aes
  legend_df <- data.frame(effect = d$effect[1], cost = d$cost[1],
                          modality = factor(present, levels = levels(d$modality)))
  p <- p +
    ggplot2::geom_point(data = legend_df, ggplot2::aes(colour = modality),
                        alpha = 0, size = 2.6, show.legend = TRUE) +
    ggplot2::scale_colour_manual(values = pal, limits = present, name = NULL) +
    ggplot2::guides(colour = ggplot2::guide_legend(
      override.aes = list(alpha = 1, size = 3))) +
    ggpop::scale_legend_icon(size = 6) +
    ggplot2::scale_x_continuous(labels = scales::label_number(accuracy = 1),
                                breaks = scales::breaks_pretty(n = 10),
                                minor_breaks = NULL) +
    ggplot2::scale_y_continuous(labels = scales::label_number(accuracy = 1,
                                                              big.mark = ",",
                                                              suffix = "M"),
                                breaks = scales::breaks_pretty(n = 10),
                                minor_breaks = NULL) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      title    = title,
      subtitle = subtitle,
      x        = "Discounted QALYs gained per 1,000",
      y        = "Discounted total costs per 1,000 (CLP millions)",
      caption  = paste0(
        "Icons on the frontier are solid; dominated strategies are faded. ",
        paste0(paste(sprintf("%s = %s", unname(icons[present]), tolower(present)),
                     collapse = ", "), "."),
        if (is.null(wtp)) "" else paste0(
          "\nOpen ring = optimal at a willingness to pay of ",
          formatC(wtp / 1e6, format = "f", digits = 1, big.mark = ","),
          "M per QALY",
          if (is.null(optimal)) "." else paste0(" (", optimal$label, ")."),
          " Chosen on the incremental ICER, not the ratio to no screening."))
    )

  # Labels sit on a near-opaque surface plaque: the icon field is dense enough
  # that bare text was unreadable where the frontier passes through it.
  if (label_frontier && nrow(frontier) > 0) {
    p <- p + ggrepel::geom_label_repel(
      data = frontier, ggplot2::aes(label = label),
      size = 2.6, colour = "#26261f",
      fill = scales::alpha("#fcfcfb", 0.92),
      label.size = 0.18, label.r = grid::unit(0.1, "lines"),
      label.padding = grid::unit(0.16, "lines"),
      segment.colour = "#57574f", segment.size = 0.35,
      min.segment.length = 0.2, box.padding = 0.55, point.padding = 0.45,
      max.overlaps = Inf, seed = 1
    )
  }

  # theme_pop() blanks axes, ticks and grid — right for a pictogram, wrong for a
  # cost-effectiveness plane, where both scales carry the result. Keep its
  # typography and legend styling, put the axes back.
  p +
    ggpop::theme_pop(base_size = 10) +
    ggplot2::theme(
      axis.title       = ggplot2::element_text(colour = "#57574f", size = 9.5),
      axis.text        = ggplot2::element_text(colour = "#77776e", size = 8.5),
      axis.ticks       = ggplot2::element_line(colour = "#c9c9c2",
                                               linewidth = 0.3),
      axis.line        = ggplot2::element_line(colour = "#c9c9c2",
                                               linewidth = 0.3),
      panel.grid.major = ggplot2::element_line(colour = "#ececE8",
                                               linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      plot.caption     = ggplot2::element_text(colour = "#77776e", hjust = 0,
                                               size = 7),
      plot.caption.position = "plot"
    )
}
