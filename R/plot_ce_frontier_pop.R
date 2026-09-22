# Cost-effectiveness plane drawn with ggpop, replacing dampack:::plot.icers().
#
# Each strategy is a Font Awesome icon rather than a dot, so modality is carried by
# shape as well as hue and survives greyscale printing and colour-vision deficiency.
# Dominated strategies stay visible but recessive; the frontier is connected and
# labelled; the strategy that is optimal at a given willingness to pay is ringed and
# annotated with its cost, QALYs gained and incremental ICER.
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
  d$icon_name <- unname(icons[as.character(d$modality)])
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
  optimal_icer <- NA_real_
  if (!is.null(wtp) && nrow(frontier) > 0) {
    icer_f <- suppressWarnings(as.numeric(icers$ICER[match(frontier$Strategy,
                                                           icers$Strategy)]))
    icer_f[is.na(icer_f)] <- 0
    ok <- which(cumsum(icer_f > wtp) == 0)
    if (length(ok)) {
      optimal <- frontier[max(ok), ]
      optimal_icer <- icer_f[max(ok)]
    }
  }

  # the selected strategy carries its numbers on the plaque, so the reader does
  # not have to cross-reference the ICER table to see what it costs and buys
  if (!is.null(optimal)) {
    frontier$label[frontier$Strategy == optimal$Strategy] <- sprintf(
      "%s\n%s  -  %.1f QALYs\n%s",
      optimal$label,
      paste0(formatC(optimal$cost, format = "f", digits = 0, big.mark = ","), "M"),
      optimal$effect,
      if (optimal_icer > 0)
        paste0("ICER ", formatC(optimal_icer / 1e6, format = "f", digits = 1),
               "M/QALY") else "reference")
  }

  pal <- c(Colonoscopy = "#2a78d6", FIT = "#eb6834", `No screening` = "#1baf7a")
  present <- names(pal)[names(pal) %in% as.character(unique(d$modality))]

  p <- ggplot2::ggplot(d, ggplot2::aes(x = effect, y = cost)) +
    ggplot2::geom_line(data = frontier, colour = "#a01b1b", linewidth = 0.7)

  # One layer for the dominated set and one for the frontier, each mapping icon in
  # aes() rather than fixing it per layer: ggplot draws EVERY show.legend layer's
  # glyph into EVERY key, so one icon per layer would overlay a stethoscope and a
  # vial in both keys. Only the frontier layer carries the legend.
  dom <- d[!d$efficient, ]
  if (nrow(dom)) {
    p <- p + ggpop::geom_icon_point(
      data = dom, ggplot2::aes(colour = modality, icon = icon_name),
      size = 1.3, alpha = 0.40, show.legend = FALSE)
  }
  p <- p + ggpop::geom_icon_point(
    data = frontier, ggplot2::aes(colour = modality, icon = icon_name),
    size = 2.3, show.legend = TRUE, legend_icons = TRUE)

  # the WTP ring goes on after the icons so it reads above them
  if (!is.null(optimal)) {
    p <- p + ggplot2::geom_point(data = optimal, shape = 21, size = 11,
                                 stroke = 1.1, colour = "#26261f", fill = NA)
  }

  p <- p +
    ggplot2::scale_colour_manual(values = pal, limits = present, name = NULL) +
    ggpop::scale_legend_icon(size = 7) +
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
        "Icons on the frontier are solid; dominated strategies are faded.",
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
    frontier$is_opt <- !is.null(optimal) & frontier$Strategy %in% optimal$Strategy
    lab_args <- list(size = 2.6, colour = "#26261f", lineheight = 1.05,
      fill = scales::alpha("#fcfcfb", 0.92),
      label.size = 0.18, label.r = grid::unit(0.1, "lines"),
      label.padding = grid::unit(0.18, "lines"),
      segment.colour = "#57574f", segment.size = 0.35,
      min.segment.length = 0.2, box.padding = 0.55, max.overlaps = Inf,
      seed = 1)
    # One call, not two: nudge_x/nudge_y are vectorised, so the selected strategy
    # can be pushed clear of its ring while every label still repels every other.
    # Two separate calls cannot see each other and collide.
    rng_x <- diff(range(d$effect)); rng_y <- diff(range(d$cost))
    p <- p + do.call(ggrepel::geom_label_repel, c(list(
      data    = frontier,
      mapping = ggplot2::aes(label = label,
                             fontface = ifelse(is_opt, "bold", "plain")),
      point.padding = 0.55,
      nudge_x = ifelse(frontier$is_opt,  0.085 * rng_x, 0),
      nudge_y = ifelse(frontier$is_opt, -0.16 * rng_y, 0)), lab_args))
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
