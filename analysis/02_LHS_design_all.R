

print(paste0("initializing LHS..."))


# 04 Generate LHS design --------------------------------------------------

SSP_cal <- calibration_setup$SSP_pathway #TRUE if we want to include SSP pathway


l_params_lhs <- sample_uniform_lhs(n_samp = n_sim, 
                                     calibration_setup$params_priors)

m_params_lhs <- l_params_lhs$m_param_samp

set.seed(NULL)

df_params_lhs <- as.data.frame(m_params_lhs) %>% 
  bind_cols(id_draw = 1:n_sim) %>% 
  relocate(id_draw)

char_date   <- format(Sys.time(), "%Y%m%d_%H%M")

m_params_lhs <- as.matrix(df_params_lhs)

#Test one model output

v_params_lhs <- m_params_lhs[2, ]

if (SSP_cal) {
  print("Testing calibration output with SSP pathway")
  l_calib_out <- tryCatch(
    calibration_out_ssp(
      v_params_calib = v_params_lhs,
      l_params_all = l_params_init,
      dt_pop = dt_pop,
      id_draw = 1,
    ),
    error = function(e) rep(NA, 10)  # fallback in case of early error
  )
} else {
  print("Testing calibration output without SSP pathway")
  l_calib_out <- tryCatch(
    calibration_out(
      v_params_calib = v_params_lhs,
      l_params_all = l_params_init,
      dt_pop = dt_pop,
      id_draw = 1,
    ),
    error = function(e) rep(NA, 10)  # fallback in case of early error
  )
}

v_age_names <- colnames(l_calib_out)

# Empty dataframe to be used in case of errors
df_empty <- as.data.frame(matrix(NA, nrow = 1, ncol = length(v_age_names)))
colnames(df_empty) <- v_age_names

# If "dt_empty" has no columns, send stop error
if (ncol(df_empty) == 0) {
  stop("Error: df_empty has no columns. Check calibration output function.")
}

# Matrix for model outputs
df_simcrc_outputs <- data.frame(
  matrix(data = NA,
         nrow = n_sim,
         ncol = length(v_age_names))
)

colnames(df_simcrc_outputs) <- v_age_names


# 06 Initialize calibration rounds ----------------------------------------
## 06.01 Serialized calibration rounds ------------------------------------

if(!parallel)  { 
  
  n_time_init_calib <- Sys.time()
  for (i in 1:n_sim) { # i <- 1
    
    v_params_lhs <- m_params_lhs[i, ]
    
    
    if(SSP_cal) {
      dt_calib_out <- tryCatch(
        {
          dt_calib_out <- calibration_out_ssp(v_params_calib = v_params_lhs,
                                          l_params_all   = l_params_init,
                                          dt_pop         = dt_pop,
                                          id_draw = v_params_lhs[which(names(v_params_lhs)=="id_draw")]  # Passing id_draw for tracking
          )
          
          df_simcrc_outputs[i, ] <- c(i, dt_calib_out)
        },
        error=function(e) {
          message(paste0("error iteration ", i, " error: ", e))
          # Choose a return value in case of error
          df_simcrc_outputs[i, ]  <- c(i, df_empty)
        }
      ) 
    } else {
      dt_calib_out <- tryCatch(
        {
          dt_calib_out <- calibration_out(v_params_calib = v_params_lhs,
                                          l_params_all   = l_params_init,
                                          dt_pop         = dt_pop,
                                          id_draw = v_params_lhs[which(names(v_params_lhs)=="id_draw")]  # Passing id_draw for tracking
          )
          
          df_simcrc_outputs[i, ] <- c(i, dt_calib_out)
        },
        error=function(e) {
          message(paste0("error iteration ", i, " error: ", e))
          # Choose a return value in case of error
          df_simcrc_outputs[i, ]  <- c(i, df_empty)
        }
      ) 
    }
    
    # Display simulation progress every 100 iterations
    if (i %% 100 == 0) {
      print(paste0("Finished iteration ", i, " of ", n_sim))
    }

  }
  n_time_end_calib <- Sys.time()
  n_time_end_calib - n_time_init_calib
  beepr::beep("mario")
}



# 06.02 Parallel calibration rounds --------------------------------------
n_time_init_par <- Sys.time()

if (parallel)  {
  os <- get_os()
  print(paste0("Parallelized PSA on ", os, " using ", no_cores, " cores."))

  # Setup SOCK cluster (no loading full snow package)
  cl_2 <- snow::makeSOCKcluster(no_cores)
  doSNOW::registerDoSNOW(cl_2)

  # Setup progress bar and single-line iteration status
  pb <- txtProgressBar(min = 1, max = n_sim, style = 3)
  progress <- function(n) {
    setTxtProgressBar(pb, n)
    cat(sprintf("[LHS %d/%d]\r", n, n_sim))  # overwrite line in console
    flush.console()
  }
  opts <- list(progress = progress)

  start_time <- Sys.time()
  n_time_init_par <- Sys.time()

  # Initialize cluster object
  #cl <- parallel::makeCluster(no_cores)
  # doParallel::registerDoParallel(cl)
  #* Define combine function
  #* To add elements in an output list see:
  #* https://stackoverflow.com/questions/27279164/output-list-of-two-rbinded-data-frames-with-foreach-in-r
  comb <- function(x, ...) {
    mapply(rbind, x, ..., SIMPLIFY = FALSE)
  }

  df_simcrc_outputs <- foreach::foreach(i = 1:n_sim,
                                        .combine = rbind,
                                        .packages = "simcrc",
                                        .export = NULL,
                                        .options.snow = opts
  ) %dopar% {

    dt_calib_out <- tryCatch(
      {

        set.seed(i)
        v_params_lhs <- m_params_lhs[i, ]
        
        if(SSP_cal) {
          dt_calib_out <- calibration_out_ssp(
            v_params_calib = v_params_lhs,
            l_params_all   = l_params_init,
            dt_pop         = dt_pop,
            id_draw = v_params_lhs[which(names(v_params_lhs)=="id_draw")]  # Passing id_draw for tracking
          )
          
        } else  {
          dt_calib_out <- calibration_out(
            v_params_calib = v_params_lhs,
            l_params_all   = l_params_init,
            dt_pop         = dt_pop,
            id_draw = v_params_lhs[which(names(v_params_lhs)=="id_draw")]  # Passing id_draw for tracking
          )
          
        }

        df_simcrc_outputs <-  dt_calib_out
      },
      error = function(e) {
        message(paste0("error iteration ", i, " error: ", e))
        df_simcrc_outputs <- df_empty
        df_simcrc_outputs$id_draw <- i
      }
    )
    dt_calib_out
  }

  # Transforming the parallel output
  df_simcrc_outputs <- as.data.frame(df_simcrc_outputs)

  for (i in 1:length(df_simcrc_outputs)) {
    df_simcrc_outputs[, i] <- as.numeric(df_simcrc_outputs[, i])
  }
  # Stop clusters
  snow::stopCluster(cl_2)
  close(pb)
  n_time_end_par <- Sys.time()
  n_time_end_par - n_time_init_par

  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "hours"))
  print(paste("Finished:", elapsed, "hours"))

}


beepr::beep("mario")
print(paste0("LHS finished"))


##### NEW VERSION
# 
# if (parallel)  {
# 
#   # Set up a parallel plan (multisession works on Windows/macOS/Linux)
#   # Use most available cores but leave 1 free.
#   library(parallelly)
#   library(future)
#   library(future.apply)
#   library(progressr)
# 
#   handlers("txtprogressbar"); handlers(global = TRUE)
# 
#   n_workers <- max(1, availableCores() - 2)
# 
#   n_workers <- no_cores
#   
#   print(paste("Parallel using:", n_workers, " cores"))
# 
#   plan(multisession, workers = n_workers)
# 
#   B <- n_sim/4  # number of rows (replicates) to generate
# 
#   # Run in parallel.
#   # future.seed = <number> gives reproducible, independent RNG streams per task.
#   # future.scheduling > 1 reduces overhead by chunking tasks per worker.
# 
#   #initialize time
#   n_time_init_par <- Sys.time()
#   n_time_batch_1 <- Sys.time()
# 
#   for (b in 1:4) {
# 
#     message(paste0("Starting batch ", b, " of 4"))
# 
#     with_progress({
#       p <- progressor(along = seq_len(B))
#       res_list <- future_lapply(
#         X = seq_len(B),
#         FUN = function(i)  {
#           on.exit(p(), add = TRUE)
#           tryCatch(
#             {
# 
#               cont <- (b - 1) * B + i
# 
#               set.seed(cont)
# 
#               v_params_lhs <- m_params_lhs[cont, ]
# 
#               dt_calib_out <- calibration_out_ssp(
#                 v_params_calib = v_params_lhs,
#                 l_params_all   = l_params_init,
#                 dt_pop         = dt_pop,
#                 id_draw       = cont
#               )
#             },
#             error = function(e) {
#               message(paste0("error iteration ", cont, " error: ", e))
#               dt_calib_out <- c(i,df_empty)
#             })
#         },
#         future.seed = 123,
#         future.scheduling = 1
#       ) })
# 
#     # Combine the 1-row data.frames into a single data.frame
#     if (b == 1) {
#       df_simcrc_outputs_1 <- do.call(rbind, res_list)
#       total_time_batch_1 <- Sys.time() - n_time_batch_1
#       print(paste0("Batch 1 time: ", round(total_time_batch_1, 2), " ", attr(total_time_batch_1, "units")))
#       n_time_batch_2 <- Sys.time()
#     }
#     if(b ==2) {
#       df_simcrc_outputs_2 <- do.call(rbind, res_list)
#       total_time_batch_2 <- Sys.time() - n_time_batch_2
#       print(paste0("Batch 2 time: ", round(total_time_batch_2, 2), " ", attr(total_time_batch_2, "units")))
#       n_time_batch_3 <- Sys.time()
#     }
#     if(b ==3) {
#       df_simcrc_outputs_3 <- do.call(rbind, res_list)
#       total_time_batch_3 <- Sys.time() - n_time_batch_3
#       print(paste0("Batch 3 time: ", round(total_time_batch_3, 2), " ", attr(total_time_batch_3, "units")))
#       n_time_batch_4 <- Sys.time()
#     }
#     if(b ==4) {
#       df_simcrc_outputs_4 <- do.call(rbind, res_list)
#       total_time_batch_4 <- Sys.time() - n_time_batch_4
#       print(paste0("Batch 4 time: ", round(total_time_batch_4, 2), " ", attr(total_time_batch_4, "units")))
#     }
# 
#   }
#   df_simcrc_outputs <- rbind(df_simcrc_outputs_1,
#                              df_simcrc_outputs_2,
#                              df_simcrc_outputs_3,
#                              df_simcrc_outputs_4)
# 
# 
#   plan(sequential)
# 
# }
# 
# n_time_end_par <- Sys.time()
# 
# # Print time
# lhs_time <- n_time_end_par - n_time_init_par
# print(paste0("LHS time: ", round(lhs_time, 2), " ", attr(lhs_time, "units")))
