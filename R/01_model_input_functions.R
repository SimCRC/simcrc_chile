



#' Update parameters
#'
#' \code{update_param_list} is used to update list of all parameters with new 
#' values for specific parameters.
#'
#' @param l_params_all List with all parameters of decision model
#' @param params_updated Parameters for which values need to be updated
#' @return 
#' A list with all parameters updated.
#' @export
update_param_list <- function(l_params_all, params_updated){
  
  if (typeof(params_updated) != "list"){
    params_updated <- split(unname(params_updated),names(params_updated)) #convert the named vector to a list
  }
  l_params_all <- modifyList(l_params_all, params_updated) #update the values
  return(l_params_all)
}
