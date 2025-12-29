#' Legendre-Matern Covariance
#'
#' The truncated Legendre-Matern covariance for the 2-sphere. The truncation is automatically obtained from the polynomial list.
#'
#' @param nu fractal parameter
#' @param alpha length scale parameter
#' @param P_l pre-computed Legendre polynomials using `legendre_Pl_list()`. List of \eqn{N\times N} matrices of length `Lmax`.
#' @returns a covariance matrix of dimension, \eqn{N\times N}. Internal use.
#' @author Aritra Halder <aritra.halder@drexel.edu>
#' @export
legendreM <- function(nu = 1.5, # smoothness
                      alpha = 5,  # scale
                      P_l = NULL){ # precomputed
  Lmax = length(P_l) - 1

  # compute spectral weights
  ell = 0:Lmax
  a_raw = (ell * (ell + 1) + alpha^2)^(-nu)

  a_l = a_raw / sum(a_raw)

  C = 0

  for (l in 0:Lmax) {
    C = C + a_l[l + 1] * P_l[[l + 1]]
  }

  # safety
  C = C + 1e-6 * diag(nrow(C))

  # return correlation matrix
  return(C)
}
