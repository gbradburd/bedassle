#' The 'bedassle' package.
#' 
#' @description A method for modeling pairwise genetic distance as a function of geographic 
#'	and ecological or environmental distance, for estimating the relative contributions 
#'	of these distances to genetic differentiation, and for statistically comparing models 
#'	with different distance predictors (e.g., geographic distance alone vs. geographic AND 
#'	ecological distance). This package contains code for running analyses (which are 
#'	implemented in the modeling language 'rstan') and visualizing and interpreting output. 
#' 
#' @docType package
#' @name bedassle-package
#' @aliases bedassle
#' @useDynLib bedassle, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import rstantools
#' @importFrom rstan sampling
#' 
#' @references 
#' Stan Development Team (2018). RStan: the R interface to Stan. R package version 2.18.2. http://mc-stan.org
#' 
NULL
