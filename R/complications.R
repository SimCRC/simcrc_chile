#' Cardiovascular complications
#' 
#' \code{cardio} is a function that applies the cardiovascular complications to the input data.
#'
#' @param x A matrix of input data
#' @return A vector of cardiovascular complications
#' 
#' @export
cardio <- function(x){
  sapply(1:nrow(x), function(i) (1/(exp(9.09053 - 0.07056 * x[i,]) + 1) - 1/(exp(9.38297  - 0.07056 * x[i,]) + 1)))
}


#' Serious GI complications
#' 
#' \code{serious_GI} is a function that applies the serious GI complications to the input data.
#' 
#' @param x A matrix of input data
#' @return A vector of serious GI complications
#' 
#' @export
serious_GI <- function(x){
  sapply(1:nrow(x), function(i) (1/(exp(9.27953 - 0.06105 * (x[i,])) + 1) - 1/(exp(10.78719 - 0.06105 * (x[i,])) + 1)))
}


#' Other GI complications
#' 
#' \code{other_GI} is a function that applies the other GI complications to the input data.
#' 
#' @param x A matrix of input data
#' @return A vector of other GI complications
#' 
#' @export
other_GI <- function(x){
  sapply(1:nrow(x), function(i) (1/(exp(8.81404 - 0.05903 * x[i,]) + 1) - 1/(exp(9.61197  - 0.05903 * x[i,]) + 1)))
}

#' Are equal
#' 
#' \code{are.equal} is a function that checks if two numbers are equal.
#' 
#' @param a A number
#' @param b A number
#' @return A boolean
#' 
#' @export
are.equal <- function(a,b){
  if(abs(a-b)<= 0.00001){
    are.equal <- TRUE
  }else{are.equal <- FALSE}
}

