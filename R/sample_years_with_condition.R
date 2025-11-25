
#' Sample Years Under Temporal Constraints
#'
#' Randomly selects a subset of years while enforcing constraints on minimum
#' sample size, maximum sample size, and the maximum allowable temporal gap
#' between consecutive selected years. Useful for generating realistic temporal
#' sampling schemes in synthetic or constrained ecological datasets.
#'
#' @title Constrained Year Subsampling
#'
#' @param years A numeric vector of years to sample from.
#' @param min_n Minimum number of years required in the sampled subset.
#' @param max_n Maximum number of years allowed in the sampled subset.
#' @param max_gap Maximum allowable gap (in years) between consecutive
#' selected years.
#' @param max_tries Maximum number of attempts before returning `NULL`
#' if no valid combination is found.
#'
#' @return A sorted numeric vector of sampled years satisfying the constraints,
#' or `NULL` if no valid combination is found after `max_tries` attempts.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Ensures sampled years fall within the specified range of sample sizes
#'   \item Enforces a temporal continuity rule (`diff(years) <= max_gap`)
#'   \item Randomly tries up to `max_tries` combinations before giving up
#' }
#'
#' @author Julie
#'
#' @export
#' 
  sample_years_with_condition <- function(years, min_n = 2, max_n = 15, max_gap = 5, max_tries = 100) {
  years <- sort(unique(years))
  n_years <- length(years)
  
  if (n_years < min_n) return(NULL)
  
  tries <- 0
  while (tries < max_tries) {
    n_pick <- sample(min_n:min(max_n, n_years), 1)
    candidate <- sort(sample(years, n_pick))
    
    if (all(diff(candidate) <= max_gap)) {
      return(candidate)
    }
    
    tries <- tries + 1
  }
  
  return(NULL)  # If no valid combo found after max_tries
}
