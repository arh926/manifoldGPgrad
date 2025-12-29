#' Gradients of Legendre-Matern Covariance
#'
#' Gradients of Legendre-Matern Covariance used for the cross covariance for posterior inference on gradients.
#'
#' @param nu optional; smoothness or fractal parameter, defaults to 1.5
#' @param alpha length scale parameter
#' @param Pterm the precomputed list of polynomial terms from `dPl()`
#' @returns a cross-covariance matrix of dimension, \eqn{N\times N_G}. Internal use.
#' @author Aritra Halder <aritra.halder@drexel.edu>
#' @export
d1legendreM <- function(nu = 1.5, # smoothness
                        alpha = 5, # scale
                        Pterm = NULL){
  Lmax = length(Pterm)

  ell = 0:Lmax
  a_raw = (ell * (ell + 1) + alpha^2)^(-nu)

  a_l = a_raw / sum(a_raw)
  C = 0

  for (l in 1:Lmax) {
    C = C + a_l[l] * Pterm[[l]] # (l/sint) * (cost * P_l[[l]] - P_lm1[[l]])
  }
  return(C)
}
