#' Posterior samples of sptaial effects and intercept
#'
#' Uses posterior Markov Chain Monte Carlo (MCMC) samples of \eqn{\sigma^2}, \eqn{\tau^2} and \eqn{\alpha} to obtain samples for \eqn{Z(x)} and \eqn{\beta_0}.
#'
#' @param mesh a mesh in polygon file format imported using `vcgPlyRead()` of class `mesh3d`.
#' @param y a \eqn{N\times 1} vector of the observed realizations
#' @param samp list containing barycentric coordinates of \eqn{N} sampled irregular locations, IDs for triangle vertices
#' @param lambda truncated vector of eigen-values
#' @param phi truncated \eqn{N\times T} matrix of eigen functions
#' @param alpha posterior MCMC samples for \eqn{\alpha} from fitting `bayes_sp_manifold()`
#' @param sig2 posterior MCMC samples for \eqn{\sigma^2} from fitting `bayes_sp_manifold()`
#' @param tau2 posterior MCMC samples for \eqn{\tau^2} from fitting `bayes_sp_manifold()`
#' @param nu value of the smoothness/fractal parameter, should be same as that used for `bayes_sp_manifold()`
#' @param silent logical; if TRUE does not print inference
#' @returns A list of posterior samples for \eqn{Z(x)} and \eqn{\beta_0}.
#' @importFrom parallel detectCores
#' @importFrom future plan multisession
#' @importFrom future.apply future_lapply
#' @author Aritra Halder <aritra.halder@drexel.edu>
#' @export
spsample <- function(mesh = NULL,
                     lambda = NULL,
                     phi = NULL,
                     y = NULL,
                     samp = NULL,
                     alpha = NULL,
                     sig2 = NULL,
                     tau2 = NULL,
                     nu = NULL,
                     silent = TRUE){
  N = length(y)

  niter = length(alpha)

  ncores <- detectCores() - 1
  samp.list <- split(1:niter, ceiling(seq_along(1:niter)/(niter/ncores)))
  # parallel.index <- 1:ncores
  # browser()

  plan(multisession, workers = ncores)
  z.beta.list = future_lapply(samp.list, function(x){
    id.x = x # samp.list[[x]]
    n.x = length(id.x)

    z.mcmc = matrix(0, ncol = N, nrow = n.x)
    beta0.mcmc = rep(0, n.x)

    sig2.x = sig2[id.x]
    alpha.x = alpha[id.x]
    tau2.x = tau2[id.x]

    for(i in 1:n.x){
      R.Z = matern_cov_mp(mesh = mesh, samp = samp,
                          lambda = lambda, phi = phi,
                          nu = nu, alpha = alpha.x[i],
                          d = 2, cor = TRUE)
      R.inv.Z = try(chol2inv(chol(R.Z)), silent = TRUE)
      if("try-error" %in% class(R.inv.Z)) R.inv.Z = chol2inv(chol(R.Z + diag(N)))

      ######################
      # Gibbs Update for Z #
      ######################
      Sig.Z.in = R.inv.Z/sig2.x[i] + diag(N)/tau2.x[i]
      chol.Sig.Z.in = try(chol(Sig.Z.in), silent = TRUE)
      if("try-error" %in% class(chol.Sig.Z.in)){
        cat("Bad Iteration::", i, "\n")
        if(i == 1) z.mcmc[i,] = 0
        else if(i == 2) z.mcmc[i,] = z.mcmc[(i - 1),]
        else if(i > 0) z.mcmc[i,] = apply(z.mcmc[1:(i - 1),], 2, median)
      }else{
        Sig.Z = chol2inv(chol.Sig.Z.in)
        mu.Z = crossprod(Sig.Z, y)/tau2.x[i]
        z.mcmc[i,] = z = as.vector(crossprod(chol(Sig.Z), rnorm(N)) + mu.Z)
      }

      #######################################
      # Intercept = mean of spatial effects #
      #######################################
      beta0.mcmc[i] = beta0 = mean(z)
    }
    return(list(z = z.mcmc,
                beta0 = beta0.mcmc))
  }, future.seed = TRUE)

  z.mcmc = do.call(rbind, lapply(z.beta.list, function(x) x$z))
  beta0.mcmc = unlist(lapply(z.beta.list, function(x) x$beta0))

  return(list(z = z.mcmc,
              beta0 = beta0.mcmc))
}
