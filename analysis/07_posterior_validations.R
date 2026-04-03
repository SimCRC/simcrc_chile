

# 05 Load posterior distributions --------------------------------------------

df_posterior_param <- read.csv(paths_calibration$path_posteriors)

#Remove lp__ and chain columns if present
if("lp__" %in% colnames(df_posterior_param)) df_posterior_param$lp__ <- NULL

m_posterior_param <- as.matrix(df_posterior_param)

n_sim <- dim(m_posterior_param)[1]

set.seed(NULL)

char_date    <- format(Sys.time(), "%Y%m%d_%H%M")

#Test one model output

v_params_lhs <- m_posterior_param[100, ]

   tryCatch(
  {
    if(calibration_setup$SSP_pathway){
      
    l_calib_out <-  calibration_out_ssp (
                    v_params_calib = v_params_lhs, 
                    l_params_all   = l_params_init, 
                    dt_pop         = dt_pop,
                    id_draw        = v_params_lhs[which(names(v_params_lhs)=="id_draw")])  # Passing id_draw for tracking
    }
    if(!calibration_setup$SSP_pathway){
      
     l_calib_out <- calibration_out (
                       v_params_calib = v_params_lhs, 
                       l_params_all   = l_params_init, 
                       dt_pop         = dt_pop,
                       id_draw       = v_params_lhs[which(names(v_params_lhs)=="id_draw")])  # Passing id_draw for tracking
    }
    
  },
  error=function(e) {
    message(paste0("error in iteration error: ", e))
  }
)

output_names <- colnames(l_calib_out) # Get the names of the outputs

#if outputs name are null stop the scrip
if(is.null(output_names)) stop("The outputs of the calibration_out function must have column names")

# Matrix for model outputs
df_simcrc_outputs <- data.frame(
  matrix(data = NA,
         nrow = dim(df_posterior_param)[1],
         ncol = length(output_names)+1 )) # +1 for chain index

colnames(df_simcrc_outputs) <- c(output_names,"chain")

# Empty dataframe to be used in case of errors
df_empty <- as.data.frame(matrix(NA, nrow = 1, ncol = length(output_names)+1))
colnames(df_empty) <- c(output_names,"chain")

### Serialization has not been implemented yet for calibration_out_ssp function
if(!parallel)  { 
  
  n_time_init_calib <- Sys.time()
  for (i in 1:n_sim) { # i <- 1
    
    v_params_lhs <- m_posterior_param[i, ]
    
       tryCatch(
      {

        if(calibration_setup$SSP_pathway){
          
          dt_calib_out <-  calibration_out_ssp (v_params_calib = v_params_lhs, 
                               l_params_all   = l_params_init, 
                               dt_pop         = dt_pop,
                               id_draw       = v_params_lhs[which(names(v_params_lhs)=="id_draw")])  # Passing id_draw for tracking
        }
        if(!calibration_setup$SSP_pathway){
          
          dt_calib_out <- calibration_out (v_params_calib = v_params_lhs, 
                           l_params_all   = l_params_init, 
                           dt_pop         = dt_pop,
                           id_draw       = v_params_lhs[which(names(v_params_lhs)=="id_draw")])  # Passing id_draw for tracking
        }
        
        df_simcrc_outputs[i, ] <- c(i, dt_calib_out)
      },
      error=function(e) {
        message(paste0("error iteration ", i, " error: ", e))
        df_simcrc_outputs <- c(id_draw = i, rep(NA, length(output_names)))
      }
    )    
    
    # Display simulation progress
    if (i/(n_sim/100) == round(i/(n_sim/100), 0)) { # display progress every 2%
      cat('\r', paste(i/n_sim * 100, "% done", sep = " "))
    }
  }
  n_time_end_calib <- Sys.time()
  n_time_end_calib - n_time_init_calib
  beepr::beep("mario")
}

## 06.02 Parallel calibration rounds --------------------------------------

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
  
  #* Define combine function
  #* To add elements in an output list see:
  #* https://stackoverflow.com/questions/27279164/output-list-of-two-rbinded-data-frames-with-foreach-in-r
  comb <- function(x, ...) {  
    mapply(rbind, x, ..., SIMPLIFY = FALSE)
  }
  n_time_init_par <- Sys.time()
  df_simcrc_outputs <- foreach::foreach(i = 1:n_sim,
                                        .combine = rbind,
                                        .packages = "simcrc",
                                        .export = NULL,
                                        .options.snow = opts
  ) %dopar% {
    
    dt_calib_out <- tryCatch(
      {
        v_params_lhs <- m_posterior_param[i, ]
        
        if(calibration_setup$SSP_pathway){
          
          dt_calib_out <-  calibration_out_ssp (v_params_calib = v_params_lhs, 
                                                l_params_all   = l_params_init, 
                                                dt_pop         = dt_pop,
                                                id_draw       = v_params_lhs[which(names(v_params_lhs)=="id_draw")])  # Passing id_draw for tracking
          dt_calib_out$chain <- v_params_lhs[which(names(v_params_lhs)=="chain")]
          }
        if(!calibration_setup$SSP_pathway){
          
          dt_calib_out <- calibration_out (v_params_calib = v_params_lhs, 
                                           l_params_all   = l_params_init, 
                                           dt_pop         = dt_pop,
                                           id_draw       = v_params_lhs[which(names(v_params_lhs)=="id_draw")])  # Passing id_draw for tracking
          dt_calib_out$chain <- v_params_lhs[which(names(v_params_lhs)=="chain")]
          }
        
        df_simcrc_outputs  <- dt_calib_out
      },
      error=function(e) {
        
        message(paste0("error iteration ", i, " error: ", e))
        df_simcrc_outputs <- df_empty
        df_simcrc_outputs$id_draw <- i
      }
    )    
  }
  
  # Transforming the parallel output
  df_simcrc_outputs <- as.data.frame(df_simcrc_outputs)
  # df_simcrc_outputs <- ldply(df_simcrc_outputs, data.frame)
  # df_simcrc_outputs <- t(df_simcrc_outputs)
  # colnames(df_simcrc_outputs)  <- df_simcrc_outputs[1,]
  # df_simcrc_outputs <- data.frame(df_simcrc_outputs[2:dim(df_simcrc_outputs)[1],])
  for (i in 1:length(df_simcrc_outputs)) {
    df_simcrc_outputs[,i] <- as.numeric(df_simcrc_outputs[,i] )
  }
  
  # Stop clusters
  snow::stopCluster(cl_2)
  close(pb)
  n_time_end_par <- Sys.time()
  n_time_end_par - n_time_init_par
  
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "hours"))
  print(paste("Finished:", elapsed, "hours"))

  beepr::beep("mario")
}

BayCANN_version <- param_BayCANN$BayCANN_version

if (save_in_BayCANN ){
  
  path_output_BayCANN_filename <- paths_calibration$path_posterior_outputs
  save(df_simcrc_outputs, file =  path_output_BayCANN_filename ) 
  print(paste0("Saved on: ",path_output_BayCANN_filename))
  
}



