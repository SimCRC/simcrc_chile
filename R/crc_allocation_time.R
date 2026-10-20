# *****************************************************************************
#
# Script: 12_crc_allocation_time.R
#
# Purpose: Allocate survival time for CRC cases between Initial, Continuing and Terminal phase
#
# Author: Jorge Roa
#
# Email: jorgeroa@stanford.edu
#
# Date Created: 07-Sept-2024
# *****************************************************************************
#
# Notes: 
#   
#
# *****************************************************************************
# *****************************************************************************
# *****************************************************************************

## 1. CRC allocation time ------------------------------------------------------
#'
#'
#' \code{crc_allocation_time} Allocation of survival time for CRC cases between Initial, Continuing and Terminal phase.
#'
#' @param stages Stages of the CRC.
#' @param death_cause CRC or other cause of death.
#' @param datatable Data table at the individual level with age of death, age of diagnosis, stage at diagnosis and cause of death.
#' @param min_age Minimum age for the allocation.
#' @param max_age Maximum age for the allocation.
#' @return 
#' Matrix with the allocation of years by phases \code{m_allocation_time}.
#' @export
crc_allocation_time <- function(stages, death_cause, datatable, min_age, max_age) {
  ## Precondition: refuse impossible survival ------------------------------
  # A row whose diagnosis is dated after its death gives a negative survival_time,
  # and the `survival_time <= 1` branch below adds it straight into the terminal-CRC
  # column with no sign check. By branch exhaustion that is the only path that can
  # emit a negative, which is why every negative cell landed in termcrc. Refuse the
  # input instead of allocating it. Mirrors the loop's own row selection, so it
  # cannot fire on rows the loop would have filtered out anyway.
  if (all(c("age_dx", "age_death", "stage_at_dx", "death_cause")
          %in% names(datatable))) {
    # the parameter and the column are both called death_cause, and inside `[` the
    # column wins; narrowed to the two causes the loop below can actually handle
    v_causes <- intersect(death_cause, c("crc", "oc"))
    dt_seen  <- datatable[!is.na(age_death) & age_death >= min_age &
                            !is.na(stage_at_dx) & stage_at_dx %in% stages &
                            death_cause %in% v_causes, ]

    n_na <- sum(is.na(dt_seen$age_dx))
    if (n_na > 0) {
      stop("crc_allocation_time: ", n_na, " row(s) reach the allocation with ",
           "age_dx missing but stage_at_dx set. Clear stage_at_dx alongside ",
           "age_dx, or these fail later as \"NA/NaN argument\".")
    }

    v_offset <- dt_seen$age_dx - dt_seen$age_death
    n_bad    <- sum(v_offset > 0, na.rm = TRUE)
    if (n_bad > 0) {
      stop("crc_allocation_time: ", n_bad, " row(s) have age_dx > age_death, ",
           "which would be allocated as negative life-years. Worst offset ",
           round(max(v_offset, na.rm = TRUE), 2), " years. ",
           "Reconcile diagnosis and death ages upstream before summarising.")
    }
  }


  
  
  ly_initial_labels <- c('ly_crcI_initial', 'ly_crcII_initial', 'ly_crcIII_initial', 'ly_crcIV_initial')
  ly_contin_labels <-  c('ly_crcI_contin', 'ly_crcII_contin', 'ly_crcIII_contin', 'ly_crcIV_contin')
  ly_termcrc_labels <-  c('ly_crcI_termcrc', 'ly_crcII_termcrc', 'ly_crcIII_termcrc', 'ly_crcIV_termcrc')
  ly_termoc_labels <-  c('ly_crcI_termoc', 'ly_crcII_termoc', 'ly_crcIII_termoc', 'ly_crcIV_termoc')
  
  all_labels <- c(ly_initial_labels, ly_contin_labels, ly_termcrc_labels, ly_termoc_labels)
  
  v_age_full <- min_age:max_age
  
  # Initialize the final matrix with the required number of rows and columns for each stage
  m_allocated_time <- matrix(0, nrow = length(v_age_full), ncol = length(all_labels) + 1)
  colnames(m_allocated_time) <- c("age", all_labels)
  m_allocated_time[, "age"] <- v_age_full
  
  # Process all cases for both survival time conditions
  for (stage in stages) {
    
    for (cause in death_cause) {
      # Subset based on death cause and stage
      if (cause == "crc") {
        subset_stage <- datatable[!is.na(age_death) & age_death >= min_age & !is.na(stage_at_dx) & stage_at_dx == stage & death_cause == "crc", ]
        terminal_col <- which(colnames(m_allocated_time) == paste0('ly_crc', stage, '_termcrc'))
      } else {
        subset_stage <- datatable[!is.na(age_death) & age_death >= min_age & !is.na(stage_at_dx) & stage_at_dx == stage & death_cause == "oc", ]
        terminal_col <- which(colnames(m_allocated_time) == paste0('ly_crc', stage, '_termoc'))
      }
      
      # Check if the subset is empty
      if (nrow(subset_stage) == 0) {
        next  # Skip to the next iteration if no cases are found
      }
      
      initial_col <- which(colnames(m_allocated_time) == paste0('ly_crc', stage, '_initial'))
      continuing_col <- which(colnames(m_allocated_time) == paste0('ly_crc', stage, '_contin'))
      
      # Iterate over all rows of the subset dataset
      for (i in 1:nrow(subset_stage)) {
        # Extract the necessary variables for the current row
        age_dx <- floor(subset_stage$age_dx[i])
        age_death <- floor(subset_stage$age_death[i])
        survival_time <- subset_stage$age_death[i] - subset_stage$age_dx[i]
        
        # Calculate age range for current row
        age_range <- age_dx:age_death
        
        if (survival_time<=1){
          
          initial_fraction <- survival_time - 1
          # Calculate age range for current row
          age_range <- age_dx:age_death
          
          # Iterate over the calculated age range
          for (age in age_range) {
            j <- which(v_age_full == age)
            
            if (age == age_death) { # First year
              m_allocated_time[j, terminal_col] <- m_allocated_time[j, terminal_col] + survival_time
              
            }
          }
        }
        
        # Differentiate based on survival time
        else if (survival_time > 1 && survival_time <= 2) {
          # Handle cases where survival time is between 1 and 2 years
          initial_fraction <- survival_time - 1
          # Calculate age range for current row
          age_range <- age_dx:age_death
          
          # Iterate over the calculated age range
          for (age in age_range) {
            j <- which(v_age_full == age)
            
            if (age == age_death) { # First year
              m_allocated_time[j, terminal_col] <- m_allocated_time[j, terminal_col] + (subset_stage$age_death[i] - age_death)
              
            } else if (age == age_death - 1) { # Second year
              
              # When survival_time is less than 2
              m_allocated_time[j, terminal_col] <- m_allocated_time[j, terminal_col] + (1 - (subset_stage$age_death[i] - age_death))
              
              m_allocated_time[j, initial_col] <- m_allocated_time[j, initial_col] + (initial_fraction)
            }
          }
          
        } else if (survival_time > 2) {
          # Handle cases where survival time is greater than 2 years
          initial_fraction <- subset_stage$age_dx[i] - age_dx
          survival_time <- subset_stage$age_death[i] - subset_stage$age_dx[i]
          
          #Calculate age range for current row
          age_range <- age_dx:age_death
          
          #Iterate over the calculated age range
          for (age in age_range) {
            j <- which(v_age_full == age)
            
            if (age == age_dx) { # First year
              m_allocated_time[j, initial_col] <- m_allocated_time[j, initial_col] + 1 - initial_fraction
              
            } else if (age == age_dx + 1) { # Second year
              
              # Handle when the age_dx is integer.
              if (initial_fraction == 0) {
                
                if (survival_time > 2 && survival_time < 3) {
                  partial_survival <- 3 - survival_time
                  m_allocated_time[j, continuing_col] <- m_allocated_time[j, continuing_col] + (1 - partial_survival)
                  m_allocated_time[j, terminal_col] <- m_allocated_time[j, terminal_col] + partial_survival
                  
                } else if (survival_time >= 3) {
                  m_allocated_time[j, continuing_col] <- m_allocated_time[j, continuing_col] + 1
                  
                } else { 
                  # When survival_time is less than 2
                  m_allocated_time[j, initial_col] <- m_allocated_time[j, initial_col] + initial_fraction
                  m_allocated_time[j, continuing_col] <- m_allocated_time[j, continuing_col] + (1 - initial_fraction)
                }
                
              } else {
                # General case for second year
                m_allocated_time[j, initial_col] <- m_allocated_time[j, initial_col] + initial_fraction
                m_allocated_time[j, continuing_col] <- m_allocated_time[j, continuing_col] + (1 - initial_fraction)
              }
              
            } else if (age < age_death - 1) { # Intermediate years
              m_allocated_time[j, continuing_col] <- m_allocated_time[j, continuing_col] + 1
              
            } else if (age == age_death - 1) { # Second to last year
              continuing_fraction <- subset_stage$age_death[i] - age_death
              m_allocated_time[j, continuing_col] <- m_allocated_time[j, continuing_col] + continuing_fraction
              m_allocated_time[j, terminal_col] <- m_allocated_time[j, terminal_col] + (1 - continuing_fraction)
              
            } else if (age == age_death) { # Last year
              m_allocated_time[j, terminal_col] <- m_allocated_time[j, terminal_col] + (subset_stage$age_death[i] - age_death)
            }
          }
        }
      }
    }
  }
  
  return(m_allocated_time)
}
