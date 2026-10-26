# *****************************************************************************
#
# Script: apply_pop_weights.R
#
# Purpose: Reweight the uspstf_summary outputs from the simulated cohort to the
#          actual Chilean (FONASA) population, by single year of age.
#
# Author: Jorge Roa and Claudia Seguin
#
# Email: jorgeroa@stanford.edu/cseguin@stanford.edu
#
# Date Created: 26-Oct-2026
# *****************************************************************************

## Apply population weights ---------------------------------------------------
#'
#' \code{apply_pop_weights} reads every uspstf_summary file in
#' \code{uspstf_folder}, multiplies all outcome columns (everything except age)
#' by a single-age weight, and writes the result, with the same file name and
#' header, to \code{out_folder}.
#'
#'   weight(a) = population(a) / n_alive_NoScreening(a)
#'
#' The denominator is n_alive from the NoScreening file of the same run, so
#' every strategy gets the same weight at a given age and differences between
#' strategies (e.g., CRC deaths averted) are kept. After weighting, NoScreening
#' n_alive equals the target population at each age.
#'
#' Note: with age-specific weights the files no longer describe one closed
#' cohort (deaths by 100 no longer add up to n_alive at 40), so they are meant
#' for population-level totals (cases, deaths, colonoscopies), not for the
#' per-1,000 CEA in ProcessUSPSTFOutput.
#'
#' @param uspstf_folder Folder with the uspstf_summary CSVs
#' (RawModelOutput_SimCRC_R).
#' @param out_folder Folder to write the weighted CSVs to. Must not be inside
#' \code{uspstf_folder}, because ProcessUSPSTFOutput searches it recursively.
#' @param dt_pop Data table with columns \code{age} and
#' \code{n_beneficiarios} (target population by single year of age).
#' @param input_prefix Only files whose name contains this prefix are used.
#' @param population_label Value written to the "population," header line.
#' @return
#' Data table of the weights by age (invisibly).
#' @export
#'
apply_pop_weights <- function(uspstf_folder,
                              out_folder,
                              dt_pop,
                              input_prefix = "",
                              population_label = "FONASA2025") {

  # Compare absolute paths. out_folder may not exist yet, so resolve its
  # nearest existing parent and add back the part that doesn't exist
  uspstf_folder <- normalizePath(uspstf_folder, mustWork = TRUE)
  if (!grepl("^(/|~|[A-Za-z]:)", out_folder)) out_folder <- file.path(getwd(), out_folder)
  out_existing <- out_folder
  while (!dir.exists(out_existing)) out_existing <- dirname(out_existing)
  out_abs <- paste0(normalizePath(out_existing),
                    substring(out_folder, nchar(out_existing) + 1))
  if (startsWith(paste0(out_abs, "/"), paste0(uspstf_folder, "/"))) {
    stop("out_folder must not be inside uspstf_folder: ProcessUSPSTFOutput ",
         "searches that folder recursively and would pick up the weighted files.")
  }
  if (!dir.exists(out_folder)) {
    dir.create(out_folder, recursive = TRUE)
  }

  v_files <- list.files(uspstf_folder,
                        pattern = sprintf(".*%s.*\\.csv$", input_prefix),
                        full.names = TRUE)
  v_noscr <- grep("NoScr", basename(v_files))
  if (length(v_noscr) != 1) {
    stop("Expected exactly one NoScreening file (contains 'NoScr') in ",
         uspstf_folder, ", found ", length(v_noscr), ".")
  }

  # *****************************************************************************
  #### 1.0 Weights from the NoScreening n_alive ####
  # *****************************************************************************

  dt_noscr <- fread(v_files[v_noscr], skip = "age,n_alive")
  dt_weights <- merge(dt_noscr[, .(age, n_alive)],
                      as.data.table(dt_pop)[, .(age, n_beneficiarios)],
                      by = "age", all.x = TRUE)
  if (anyNA(dt_weights$n_beneficiarios)) {
    stop("dt_pop is missing ages: ",
         paste(dt_weights[is.na(n_beneficiarios), age], collapse = ", "))
  }
  dt_weights[, weight := n_beneficiarios / n_alive]

  # *****************************************************************************
  #### 2.0 Weight every file and write it with the same header ####
  # *****************************************************************************

  for (file_in in v_files) {
    v_lines <- readLines(file_in)
    data_start <- grep("^age,n_alive", v_lines)
    stopifnot(length(data_start) == 1)

    # Keep the header lines (including the blank line before the data)
    # verbatim, so ProcessUSPSTFOutput's fixed skip counts still work
    v_header <- v_lines[seq_len(data_start - 1)]
    v_header <- sub("^population,.*$", paste0("population,", population_label),
                    v_header)

    dt_strategy <- fread(file_in, skip = "age,n_alive")
    stopifnot(identical(dt_strategy$age, dt_weights$age))

    v_outcomes <- setdiff(names(dt_strategy), "age")
    dt_strategy[, (v_outcomes) := lapply(.SD, function(x) x * dt_weights$weight),
                .SDcols = v_outcomes]

    file_out <- file.path(out_folder, basename(file_in))
    writeLines(v_header, file_out)
    # Appending column names after the header is intended; silence that warning
    withCallingHandlers(
      write.table(dt_strategy,
                  file = file_out,
                  col.names = TRUE, row.names = FALSE, quote = FALSE, append = TRUE, sep = ","),
      warning = function(w) {
        if (grepl("appending column names", conditionMessage(w))) invokeRestart("muffleWarning")
      }
    )
  }

  invisible(dt_weights)
}
