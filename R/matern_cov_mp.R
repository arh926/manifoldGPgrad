#' Matern Covariance for compact Riemanian manifolds
#'
#' Computes the Matern covariance matrix.
#'
#' @param mesh a mesh in polygon file format imported using `vcgPlyRead()` of class `mesh3d`.
#' @param samp list containing barycentric coordinates of \eqn{N} sampled irregular locations, IDs for triangle vertices
#' @param lambda truncated vector of eigen-values
#' @param phi truncated \eqn{N\times T} matrix of eigen functions
#' @param nu fractal parameter controlling process smoothness
#' @param alpha length scale parameter
#' @param d optional; manifold dimension, defaults to 2
#' @param D.out logical, returns the mean of the scaling used, \eqn{C_{\nu,\alpha}(x)}
#' @param cor logical, if TRUE returns a correlation matrix
#' @returns A \eqn{N\times N} covariance or correlation matrix for the realization
#' @author Aritra Halder <aritra.halder@drexel.edu>
#' @export
matern_cov_mp <- function(mesh = NULL,
                          samp = NULL,
                          lambda = NULL,
                          phi = NULL,
                          nu = NULL,
                          alpha = NULL,
                          d = 2,
                          D.out = FALSE,
                          cor = TRUE){
  # triangular faces
  F = t(mesh$it)

  # barycentric coordinates
  u = samp$bary[,1]
  v = samp$bary[,2]
  w = samp$bary[,3]

  # vertices
  v1_idx = F[samp$tri_idx, 1]
  v2_idx = F[samp$tri_idx, 2]
  v3_idx = F[samp$tri_idx, 3]

  # barycentric interpolation
  # Phi_pts: n_points x k_eig
  Phi_pts = u * phi[v1_idx, , drop = FALSE] +
    v * phi[v2_idx, , drop = FALSE] +
    w * phi[v3_idx, , drop = FALSE]

  # spectral weights
  expo   = nu + d/2
  coeffs = (alpha^2 + lambda)^(-expo)

  # K = Phi diag(coeffs) Phi^T
  # avoid forming full diag when possible
  PhiW = t(t(Phi_pts) * sqrt(coeffs))  # N x k, scale each column
  R = PhiW %*% t(PhiW)              # N x N covariance matrix

  if(D.out) cat("Estimated C_nu_alpha:", mean(diag(R)), "\n")

  if(cor){
    # scaling: use when fitting the process
    # otherwise posterior estimates are whack
    diag_R = diag(R)
    D_inv  = 1 / sqrt(diag_R)
    R_star = (D_inv * R) * D_inv
    R_star = (R_star + t(R_star))/2

    return(R_star)
  }else{
    return(R)
  }
}
