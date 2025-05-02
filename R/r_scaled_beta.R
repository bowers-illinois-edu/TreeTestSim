#' Sample from a Scaled Beta Distribution
#'
#' @description
#' Generates random samples from a Beta distribution and then linearly transforms
#' them to lie in the interval \eqn{[a, 1]}.
#'
#' @param a Numeric. The lower bound of the output interval. Must be between 0 and 1 (exclusive).
#' @param shape1 Numeric. The first shape parameter of the Beta distribution.
#' @param shape2 Numeric. The second shape parameter of the Beta distribution.
#' @param n Integer. The number of random samples to generate (default is 1).
#'
#' @return A numeric vector of length n containing the transformed samples.
#' @examples
#' r_scaled_beta(0.2, 2, 5, n = 10)
#' @export
r_scaled_beta <- function(a, shape1, shape2, n = 1) {
  if (!is.numeric(a) || a < 0 || a >= 1) {
    stop("Parameter 'a' must be numeric and in the interval [0, 1).")
  }
  beta_samples <- rbeta(n, shape1, shape2)
  scaled_samples <- a + (1 - a) * beta_samples
  return(scaled_samples)
}
