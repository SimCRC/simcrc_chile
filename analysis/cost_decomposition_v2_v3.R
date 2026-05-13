# *****************************************************************************
#
# Script: cost_decomposition_v2_v3.R
#
# Purpose: Dissect how total strategy cost is composed and verify it matches
#          the dampack ICER table numbers.
#
#   Cost components (all undiscounted, per 1,000 people):
#     FIT_costsper1000 + COL_costsper1000 + Complication_costsper1000
#     + CostsCancerCareper1000 = TotalCostsper1000          ← sum visible by eye
#
#   Dampack then discounts TotalCostsper1000 → DiscountedTotalCostsper1000
#   which is the "Cost" column in the ICER table.
#
#   Source files matched to ICER CSVs by DiscountedTotalCostsper1000:
#     v2 → 2026-04-08 141620.826308_model_run_data_Base.xlsx
#     v3 → 2026-05-09 142958.018178_model_run_data_Base.xlsx
#
#   Output:
#     table_all_strategies_v2_v3.html  — COL + FIT ND strategies
#     table_fit_only_v2_v3.html        — FIT ND strategies only
#
# Author: Jorge Roa
# Date Created: 12-May-2026
#
# *****************************************************************************

remove(list = ls())
gc()

library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(gt)
library(stringr)
library(ggplot2)


# *****************************************************************************
#### 1.0 Load model_run_data (cost components) ####
# *****************************************************************************

FILE_V2 <- "ce_results/2026-04-08 141620.826308_model_run_data_Base.xlsx"
FILE_V3 <- "ce_results/2026-05-09 142958.018178_model_run_data_Base.xlsx"

v_cost_cols <- c(
  "Strategy",
  "FIT_costsper1000",
  "COL_costsper1000",
  "Complication_costsper1000",
  "CostsCancerCareper1000",
  "TotalCostsper1000",
  "DiscountedTotalCostsper1000"
)

read_mrd <- function(path, version) {
  read_excel(path, sheet = "model_run_data") |>
    select(all_of(v_cost_cols)) |>
    mutate(
      version  = version,
      Strategy = str_remove(Strategy, "^2026Chile_")
    )
}

df_mrd_v2 <- read_mrd(FILE_V2, "v2")
df_mrd_v3 <- read_mrd(FILE_V3, "v3")


# *****************************************************************************
#### 2.0 Load dampack ICER variables ####
# *****************************************************************************

# All-strategies ICER tables
df_icer_all_v2 <- read_csv("ce_results/df_icer_all_strategies_v2.csv",
                            show_col_types = FALSE)
df_icer_all_v3 <- read_csv("ce_results/df_icer_all_strategies.csv",
                            show_col_types = FALSE)

# FIT-only ICER tables
df_icer_fit_v2 <- read_csv("ce_results/df_icer_FIT_v2.csv",
                            show_col_types = FALSE)
df_icer_fit_v3 <- read_csv("ce_results/df_icer_FIT.csv",
                            show_col_types = FALSE)


# *****************************************************************************
#### 3.0 Helper: build wide comparison table for a strategy set ####
# *****************************************************************************

#' Join model_run_data components with dampack ICER variables for two versions
#'
#' @param v_strategies Character vector of strategy names to include
#' @param icer_v2      Dampack ICER data frame for v2
#' @param icer_v3      Dampack ICER data frame for v3
#' @param order_by_col Column name to sort rows (default: DiscountedTotalCostsper1000 v3)
build_wide <- function(v_strategies, icer_v2, icer_v3,
                       order_by_col = "DiscountedTotalCostsper1000_v3") {

  icer_v2_sel <- icer_v2 |>
    filter(Strategy %in% v_strategies) |>
    select(Strategy, Effect,
           Inc_Cost_v2 = Inc_Cost, ICER_v2 = ICER, Status_v2 = Status)

  icer_v3_sel <- icer_v3 |>
    filter(Strategy %in% v_strategies) |>
    select(Strategy,
           Inc_Cost_v3 = Inc_Cost, ICER_v3 = ICER, Status_v3 = Status)

  mrd_v2 <- df_mrd_v2 |>
    filter(Strategy %in% v_strategies) |>
    rename_with(~ paste0(.x, "_v2"), where(is.numeric))

  mrd_v3 <- df_mrd_v3 |>
    filter(Strategy %in% v_strategies) |>
    rename_with(~ paste0(.x, "_v2"), where(is.numeric)) |>  # temp rename
    rename_with(~ str_replace(.x, "_v2$", "_v3"), everything())

  mrd_wide <- full_join(
    mrd_v2 |> select(-version),
    mrd_v3 |> select(-version),
    by = "Strategy"
  )

  df <- mrd_wide |>
    left_join(icer_v2_sel, by = "Strategy") |>
    left_join(icer_v3_sel, by = "Strategy") |>
    mutate(across(c(Inc_Cost_v2, Inc_Cost_v3, ICER_v2, ICER_v3), ~ . / 1e6),
           across(where(is.numeric), ~ . / 1e6)) |>
    arrange(.data[[order_by_col]])

  df
}


# *****************************************************************************
#### 4.0 Helper: render gt table ####
# *****************************************************************************

#' Build a gt table showing components (v2 & v3) + dampack ICER variables
#'
#' Columns per version:
#'   FIT | COL | Complications | Cancer care | Total (undiscounted) ← sum by eye
#'   Discounted total | Inc. cost | ICER | Status                  ← dampack
make_gt <- function(df, title_text, subtitle_text, wtp_m = 16) {

  df |>
    select(
      Strategy, Effect,
      FIT_costsper1000_v2, COL_costsper1000_v2,
      Complication_costsper1000_v2, CostsCancerCareper1000_v2,
      TotalCostsper1000_v2,
      DiscountedTotalCostsper1000_v2, Inc_Cost_v2, ICER_v2, Status_v2,
      FIT_costsper1000_v3, COL_costsper1000_v3,
      Complication_costsper1000_v3, CostsCancerCareper1000_v3,
      TotalCostsper1000_v3,
      DiscountedTotalCostsper1000_v3, Inc_Cost_v3, ICER_v3, Status_v3
    ) |>
    gt(rowname_col = "Strategy") |>
    tab_header(title = title_text, subtitle = subtitle_text) |>
    # ── v2 spanners ──
    tab_spanner(
      label   = "v2 — cost components (undiscounted M CLP per 1,000)",
      columns = c(FIT_costsper1000_v2, COL_costsper1000_v2,
                  Complication_costsper1000_v2, CostsCancerCareper1000_v2,
                  TotalCostsper1000_v2)
    ) |>
    tab_spanner(
      label   = "v2 — dampack (screen_costs_v2)",
      columns = c(DiscountedTotalCostsper1000_v2, Inc_Cost_v2,
                  ICER_v2, Status_v2)
    ) |>
    # ── v3 spanners ──
    tab_spanner(
      label   = "v3 — cost components (undiscounted M CLP per 1,000)",
      columns = c(FIT_costsper1000_v3, COL_costsper1000_v3,
                  Complication_costsper1000_v3, CostsCancerCareper1000_v3,
                  TotalCostsper1000_v3)
    ) |>
    tab_spanner(
      label   = "v3 — dampack (screen_costs_v3)",
      columns = c(DiscountedTotalCostsper1000_v3, Inc_Cost_v3,
                  ICER_v3, Status_v3)
    ) |>
    # ── column labels ──
    cols_label(
      Effect                         = "QALYs gained",
      FIT_costsper1000_v2            = "FIT",
      COL_costsper1000_v2            = "COL",
      Complication_costsper1000_v2   = "Complications",
      CostsCancerCareper1000_v2      = "Cancer care",
      TotalCostsper1000_v2           = "Total",
      DiscountedTotalCostsper1000_v2 = "Disc. total",
      Inc_Cost_v2                    = "Inc. cost",
      ICER_v2                        = "ICER",
      Status_v2                      = "Status",
      FIT_costsper1000_v3            = "FIT",
      COL_costsper1000_v3            = "COL",
      Complication_costsper1000_v3   = "Complications",
      CostsCancerCareper1000_v3      = "Cancer care",
      TotalCostsper1000_v3           = "Total",
      DiscountedTotalCostsper1000_v3 = "Disc. total",
      Inc_Cost_v3                    = "Inc. cost",
      ICER_v3                        = "ICER",
      Status_v3                      = "Status"
    ) |>
    fmt_number(columns = Effect, decimals = 1) |>
    fmt_number(
      columns  = c(FIT_costsper1000_v2, COL_costsper1000_v2,
                   Complication_costsper1000_v2, CostsCancerCareper1000_v2,
                   TotalCostsper1000_v2, DiscountedTotalCostsper1000_v2,
                   Inc_Cost_v2,
                   FIT_costsper1000_v3, COL_costsper1000_v3,
                   Complication_costsper1000_v3, CostsCancerCareper1000_v3,
                   TotalCostsper1000_v3, DiscountedTotalCostsper1000_v3,
                   Inc_Cost_v3),
      decimals = 1
    ) |>
    fmt_number(columns = c(ICER_v2, ICER_v3), decimals = 2) |>
    sub_missing(missing_text = "—") |>
    # Bold Total and Disc. total (the two "sum" columns)
    tab_style(
      style     = cell_text(weight = "bold"),
      locations = cells_body(
        columns = c(TotalCostsper1000_v2, DiscountedTotalCostsper1000_v2,
                    TotalCostsper1000_v3, DiscountedTotalCostsper1000_v3)
      )
    ) |>
    # Cancer care in light blue (dominant, not driven by screening prices)
    tab_style(
      style     = cell_fill(color = "#EFF6FF"),
      locations = cells_body(
        columns = c(CostsCancerCareper1000_v2, CostsCancerCareper1000_v3)
      )
    ) |>
    # Highlight cost-effective rows under WTP
    tab_style(
      style     = cell_fill(color = "#DBEAFE"),
      locations = cells_body(
        rows    = !is.na(ICER_v2) & ICER_v2 <= wtp_m,
        columns = c(DiscountedTotalCostsper1000_v2, ICER_v2, Status_v2)
      )
    ) |>
    tab_style(
      style     = cell_fill(color = "#D1FAE5"),
      locations = cells_body(
        rows    = !is.na(ICER_v3) & ICER_v3 <= wtp_m,
        columns = c(DiscountedTotalCostsper1000_v3, ICER_v3, Status_v3)
      )
    ) |>
    tab_source_note(md(
      "**Components** (undiscounted): FIT + COL + Complications + Cancer care = **Total** — verify by eye.
       Dampack discounts Total → **Disc. total** (`DiscountedTotalCostsper1000`) used for ICER.
       `Inc. cost` = `Inc_Cost` | `Status`: ND non-dominated, ED extended dominance, D dominated.
       Blue = cost-effective v2 | Green = cost-effective v3 | WTP = $16M CLP/QALY."
    )) |>
    tab_options(
      table.font.size           = px(11),
      heading.align             = "left",
      column_labels.font.weight = "bold"
    )
}


# *****************************************************************************
#### 5.0 Table 1: all strategies (COL + FIT) ####
# *****************************************************************************

v_nd_all <- union(
  df_icer_all_v2$Strategy[df_icer_all_v2$Status == "ND"],
  df_icer_all_v3$Strategy[df_icer_all_v3$Status == "ND"]
)

df_all <- build_wide(v_nd_all, df_icer_all_v2, df_icer_all_v3)

tbl_all <- make_gt(
  df_all,
  title_text    = "Cost breakdown — all strategies (COL + FIT), ND frontier",
  subtitle_text = "M CLP per 1,000 people | v2 vs v3 screening costs"
)

gtsave(tbl_all, "ce_results/table_all_strategies_v2_v3.html")
cat("Saved: ce_results/table_all_strategies_v2_v3.html\n")


# *****************************************************************************
#### 6.0 Table 2: FIT strategies only ####
# *****************************************************************************

v_nd_fit <- union(
  df_icer_fit_v2$Strategy[df_icer_fit_v2$Status == "ND"],
  df_icer_fit_v3$Strategy[df_icer_fit_v3$Status == "ND"]
)

df_fit <- build_wide(v_nd_fit, df_icer_fit_v2, df_icer_fit_v3)

tbl_fit <- make_gt(
  df_fit,
  title_text    = "Cost breakdown — FIT strategies only, ND frontier",
  subtitle_text = "M CLP per 1,000 people | v2 vs v3 screening costs"
)

gtsave(tbl_fit, "ce_results/table_fit_only_v2_v3.html")
cat("Saved: ce_results/table_fit_only_v2_v3.html\n")


# *****************************************************************************
#### 7.0 Driver analysis: what explains the v2 → v3 cost increase? ####
#
# For each ND FIT strategy:
#   Δ Total = Δ FIT + Δ COL + Δ Compl + Δ Care
#
# Key finding: COL is the primary driver (57–64% of Δ Total).
# Within COL: ~40% of colonoscopies yield polypectomy (priced +196%),
#   ~60% do not (+49%). Weighted effective COL price: +116%.
# Cancer care contributes minimally (0–4%) despite being the dominant
#   absolute cost, because screening prices do not affect it.
# *****************************************************************************

v_nd_fit <- union(
  df_icer_fit_v2$Strategy[df_icer_fit_v2$Status == "ND"],
  df_icer_fit_v3$Strategy[df_icer_fit_v3$Status == "ND"]
)

v_mrd_driver_cols <- c(
  "Strategy",
  "FITper1000", "Colonoscopiesper1000",
  "FIT_costsper1000", "COL_costsper1000",
  "Complication_costsper1000", "CostsCancerCareper1000",
  "TotalCostsper1000"
)

mrd_v2_d <- read_excel(FILE_V2, sheet = "model_run_data") |>
  mutate(Strategy = str_remove(Strategy, "^2026Chile_")) |>
  filter(Strategy %in% v_nd_fit) |>
  select(all_of(v_mrd_driver_cols))

mrd_v3_d <- read_excel(FILE_V3, sheet = "model_run_data") |>
  mutate(Strategy = str_remove(Strategy, "^2026Chile_")) |>
  filter(Strategy %in% v_nd_fit) |>
  select(all_of(v_mrd_driver_cols))

df_driver <- left_join(mrd_v2_d, mrd_v3_d, by = "Strategy", suffix = c("_v2", "_v3")) |>
  mutate(
    d_FIT   = (FIT_costsper1000_v3            - FIT_costsper1000_v2)            / 1e6,
    d_COL   = (COL_costsper1000_v3            - COL_costsper1000_v2)            / 1e6,
    d_Compl = (Complication_costsper1000_v3   - Complication_costsper1000_v2)   / 1e6,
    d_Care  = (CostsCancerCareper1000_v3      - CostsCancerCareper1000_v2)      / 1e6,
    d_Total = (TotalCostsper1000_v3           - TotalCostsper1000_v2)           / 1e6,
    pct_FIT   = d_FIT   / d_Total * 100,
    pct_COL   = d_COL   / d_Total * 100,
    pct_Compl = d_Compl / d_Total * 100,
    pct_Care  = d_Care  / d_Total * 100,
    # Effective price per event
    eff_FIT_v2  = FIT_costsper1000_v2 / FITper1000_v2,
    eff_FIT_v3  = FIT_costsper1000_v3 / FITper1000_v3,
    eff_COL_v2  = COL_costsper1000_v2 / Colonoscopiesper1000_v2,
    eff_COL_v3  = COL_costsper1000_v3 / Colonoscopiesper1000_v3,
    pct_eff_FIT = (eff_FIT_v3 - eff_FIT_v2) / eff_FIT_v2 * 100,
    pct_eff_COL = (eff_COL_v3 - eff_COL_v2) / eff_COL_v2 * 100,
    Strategy = factor(Strategy, levels = v_nd_fit)
  ) |>
  arrange(d_Total)


# 7a. Table: absolute Δ + % contribution per component -------------------------

tbl_driver <- df_driver |>
  select(Strategy, d_FIT, pct_FIT, d_COL, pct_COL,
         d_Compl, pct_Compl, d_Care, pct_Care, d_Total) |>
  gt(rowname_col = "Strategy") |>
  tab_header(
    title    = "What drives the v2 → v3 cost increase?",
    subtitle = "Δ cost per component (M CLP, undiscounted) and % of total Δ | ND FIT strategies"
  ) |>
  tab_spanner("FIT tests (+51%)",          columns = c(d_FIT,   pct_FIT))   |>
  tab_spanner("Colonoscopies (+116% eff.)", columns = c(d_COL,   pct_COL))   |>
  tab_spanner("Complications (+196%)",      columns = c(d_Compl, pct_Compl)) |>
  tab_spanner("Cancer care (unchanged)",    columns = c(d_Care,  pct_Care))  |>
  cols_label(
    d_FIT = "Δ (M$)", pct_FIT = "%",
    d_COL = "Δ (M$)", pct_COL = "%",
    d_Compl = "Δ (M$)", pct_Compl = "%",
    d_Care  = "Δ (M$)", pct_Care  = "%",
    d_Total = "Δ Total (M$)"
  ) |>
  fmt_number(columns = c(d_FIT, d_COL, d_Compl, d_Care, d_Total), decimals = 1) |>
  fmt_number(columns = c(pct_FIT, pct_COL, pct_Compl, pct_Care),
             decimals = 1, pattern = "{x}%") |>
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_body(columns = c(d_COL, pct_COL, d_Total))) |>
  tab_style(style = cell_fill(color = "#FEF3C7"),
            locations = cells_body(columns = c(d_COL, pct_COL))) |>
  tab_style(style = cell_fill(color = "#EFF6FF"),
            locations = cells_body(columns = c(d_Care, pct_Care))) |>
  tab_source_note(md(
    "**Colonoscopy effective price**: 124K → 279K CLP (+116%) because ~40% of colonoscopies
     yield polypectomy (priced at +196%) and ~60% do not (+49%).
     **Cancer care** is unchanged by screening unit prices — it depends on stage at detection."
  )) |>
  tab_options(table.font.size = px(11), heading.align = "left",
              column_labels.font.weight = "bold")

gtsave(tbl_driver, "ce_results/table_driver_analysis_v2_v3.html")
cat("Saved: ce_results/table_driver_analysis_v2_v3.html\n")


# 7b. Plot: % contribution of each component to Δ Total -----------------------

df_driver_long <- df_driver |>
  select(Strategy, pct_FIT, pct_COL, pct_Compl, pct_Care) |>
  pivot_longer(-Strategy, names_to = "component", values_to = "pct") |>
  mutate(
    component = recode(component,
      pct_FIT   = "FIT tests (+51%)",
      pct_COL   = "Colonoscopies (+116% eff.)",
      pct_Compl = "Complications (+196%)",
      pct_Care  = "Cancer care (unchanged)"
    ),
    component = factor(component, levels = c(
      "Cancer care (unchanged)",
      "Complications (+196%)",
      "FIT tests (+51%)",
      "Colonoscopies (+116% eff.)"
    ))
  )

pal_driver <- c(
  "Colonoscopies (+116% eff.)"  = "#F97316",
  "FIT tests (+51%)"            = "#4E9AF1",
  "Complications (+196%)"       = "#E76F51",
  "Cancer care (unchanged)"     = "#94A3B8"
)

p_driver <- ggplot(df_driver_long,
                   aes(x = Strategy, y = pct, fill = component)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = c(25, 50, 75), linetype = "dashed",
             colour = "white", linewidth = 0.4) +
  scale_fill_manual(values = pal_driver) +
  scale_y_continuous(labels = scales::label_percent(scale = 1),
                     breaks = c(0, 25, 50, 75, 100)) +
  labs(
    title    = "What drives the v2 → v3 cost increase?",
    subtitle = "% of total Δ cost attributable to each component | ND FIT strategies",
    x        = NULL,
    y        = "% of total cost increase",
    fill     = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position    = "bottom",
    axis.text.x        = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank()
  ) +
  guides(fill = guide_legend(nrow = 2))

ggsave("ce_results/plot_driver_analysis_v2_v3.svg",
       p_driver, width = 9, height = 6, units = "in")
cat("Saved: ce_results/plot_driver_analysis_v2_v3.svg\n")
