#' Second derivative of Legendre-Matern Covariance
#'
#' Second derivative of Legendre-Matern Covariance used for the variance for posterior inference on gradients.
#'
#' @param Lmax the truncation used for the Legendre-Matern covariance
#' @param nu optional; smoothness or fractal parameter, defaults to 1.5
#' @param alpha length scale parameter
#' @returns a variance matrix of dimension, \eqn{1\times 1}. Internal use.
#' @author Aritra Halder <aritra.halder@drexel.edu>
#' @export
d2legendreM0 <- function(Lmax = 20, # truncation level; increase for more detail (but O(n^2 * L_max))
                         nu = 1.5, # smoothness
                         alpha = 5){ # scale

  # compute spectral weights a_l ~ (l(l+1) + alpha^2)^(-nu)
  ell = 0:Lmax
  a_raw = (ell * (ell + 1) + alpha^2)^(-nu) # actually: (ell * (ell + 1) + alpha^2)^(-nu)

  # normalize so variance at zero is 1:
  # For gamma = 0, cos_gamma = 1, P_l(1) = 1.
  # So C(0) = sum_l a_l, so we want sum_l a_l = 1.
  a_l = a_raw / sum(a_raw)

  C = 0

  for (l in 1:Lmax) {
    C = C + a_l[l] * l * (l + 1)
  }
  return(C/2)
}
