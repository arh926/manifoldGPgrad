#' Legendre Polynomials
#'
#' The Legendre polynomial list truncated at `Lmax`.
#'
#' @param coords coordinates on the 2-sphere
#' @param Lmax truncation
#' @returns a list of \eqn{N\times N} matrices of length `Lmax`. Internal use.
#' @importFrom gsl legendre_Pl
#' @author Aritra Halder <aritra.halder@drexel.edu>
#' @export
legendre_Pl_list <- function(coords = NULL,
                             Lmax = 20){
  # Dot products between all pairs
  dotprod = tcrossprod(coords, coords)

  # safety clamp
  dotprod[dotprod > 1]  = 1
  dotprod[dotprod < -1] = -1

  gamma = acos(dotprod)   # angular distance in [0, pi]
  cos_gamma = cos(gamma)  # effectively = dotprod

  P_l_list = list()

  for(l in 0:Lmax){
    P_l_list[[(l + 1)]] = legendre_Pl(l, cos_gamma)
  }
  return(P_l_list)
}
