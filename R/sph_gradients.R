#' Posterior inference on gradients
#'
#' Posterior Markov Chain Monte Carlo (MCMC) samples for gradients along a pre-supplied vector field 1-1 posterior samples of \eqn{\sigma^2} and \eqn{\alpha}.
#'
#' @param grid_xyz grid locations matrix of dimension \eqn{N_G\times 3} on the 2-sphere
#' @param model a list containing posterior post burn-in samples of \eqn{\sigma^2} and \eqn{\alpha}
#' @param Lmax optional; defaults to 30
#' @returns A list of the posterior samples for gradients at the grid locations.
#' @author Aritra Halder <aritra.halder@drexel.edu>
#' @export
sph_gradients <- function(model = NULL,
                          Lmax = 30,
                          grid_xyz = NULL){
  nmcmc = length(model$sig2)
  coords = model$coords

  ngrid = nrow(grid_xyz)
  nu = model$nu

  d = sapply(1:nrow(grid_xyz), function(y) apply(coords, 1, function(x){
    t = crossprod(x, grid_xyz[y,])
    t = ifelse(t>1, 1, ifelse(t < -1, -1, t))
    t
  }))

  tmp1 = sapply(1:nrow(grid_xyz), function(y) apply(coords, 1, function(x){
    -x[1] * grid_xyz[y, 2] + x[2] * grid_xyz[y, 1]
  }))

  tmp2 = sapply(1:nrow(grid_xyz), function(y) sum(grid_xyz[y, 1:2]^2))

  grad.samples.grid = matrix(NA, nrow = nmcmc, ncol = ngrid)

  # Compute polynomials once
  Pl = legendre_Pl_list(coords = coords, Lmax = Lmax)
  Pterm = dPl(d = d, Lmax = Lmax)

  for(i in 1:nmcmc){
    sig2.x = model$sig2[i]
    alpha.x = model$alpha[i]
    z.x = model$z[i,]

    C = legendreM(nu = nu, alpha = alpha.x, P_l = Pl)
    C.inv = chol2inv(chol(C))

    for(j in 1:ngrid){
      C1 = d1legendreM(nu = nu, alpha = alpha.x,
                       Pterm = lapply(Pterm, function(x) x[, j])) * tmp1[, j]/sqrt(1 - d[, j]^2)

      mu.grad = -crossprod(t(crossprod(C1, C.inv)), z.x)
      Sigma.grad = sig2.x * (d2legendreM0(Lmax = 20, nu = nu, alpha = alpha.x) * tmp2[j] + crossprod(t(crossprod(C1, C.inv)), C1))

      grad.samples.grid[i, j] = rnorm(1, mu.grad, sqrt(Sigma.grad))
    }
    if(i %% 100 == 0){
      if(i == 100) cat("Iteration ", i, "...", "\t")
      else cat(i, "...", "\t")
    }
  }
  return(grad.samples.grid)
}

