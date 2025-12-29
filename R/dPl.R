#' Derivative of Legendre Polynomials
#'
#' Function that pre-computes to derivatives of Legendre-Polynomials
#'
#' @param d a matrix of \eqn{N\times N_G} of distances between observed locations and grid locations
#' @param Lmax optional; truncation used for Legendre-Matern for original process, defaults to 30
#' @importFrom gsl legendre_Pl
#' @returns a list of length `Lmax` containing the polynomial expression term
#' @author Aritra Halder <aritra.halder@drexel.edu>
#' @export
dPl <- function(d = NULL,
                Lmax = 30){

  # d = crossprod(x1, x2)
  # d = ifelse(d>1, 1, ifelse(d < -1, -1, d))

  t = acos(d)
  cost = cos(t)
  sint = sin(t)

  P_term = list()
  for (l in 1:Lmax) {
    # evaluate P_l at all entries of cos_gamma
    P_l = legendre_Pl(l, cost)  # same shape as cos_gamma
    P_lm1 = legendre_Pl(l - 1, cost)  # same shape as cos_gamma
    P_term[[l]] = (l/sint) * (cost * P_l - P_lm1)
  }
  return(P_term)
}
