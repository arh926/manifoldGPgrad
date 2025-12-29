#' Hierarchial spatial model for the 2-sphere
#'
#' Function that fits the hierarchical maodel, \eqn{Y(x)=\beta_0+Z(x)+\epsilon(x)}, \eqn{Z(x)\sim GP(0, K(\cdot;\sigma^2,\alpha))}, and \eqn{\epsilon(x)\sim N(0,\tau^2)}.
#'
#' @param y observed process
#' @param coords observed locations
#' @param niter optional; number of MCMC iterations, defaults to 5000
#' @param nburn optional; number of MCMC burn-in, defaults to 2500
#' @param report optional; MCMC batch length, defaults to 100
#' @param a.v optional; hyperparameter for the shape of the inverse-gamma prior on total variance, defaults to 2
#' @param b.v optional; hyperparameter for the scale of the inverse-gamma prior on total variance, defaults to 2
#' @param lower.alpha optional; hyperparameter for the lower bound for the uniform prior on the length scale parameter, defaults to 0.01
#' @param upper.alpha optional; hyperparameter for the upper bound for the uniform prior on the length scale parameter, defaults to 30
#' @param nu optional; value of the smoothness/fractal parameter, defaults to 1.5
#' @param Lmax optinal; truncation for the Legendre-Matern covariance, defaults to 30
#' @param verbose logical; if TRUE prints inference
#' @param digits optional; precision of print output, defaults to 3 digits
#' @returns A list of the observed process, sampled locations, fractal parameter, posterior post burn-in samples for \eqn{\sigma^2, \tau^2, \alpha}.
#' @importFrom stats median quantile rnorm
#' @author Aritra Halder <aritra.halder@drexel.edu>
#' @export
##########################################################
### Collapsed MH Sampler for Spatial Process on Sphere ###
##########################################################
# Re-parameterized version due to identifiability issues with sig2 and alpha
sphere_bayes_sp <- function(y = NULL,
                            coords = NULL, # (x, y, z) sphere points
                            niter = NULL, nburn = NULL, report = NULL,
                            a.v = NULL, b.v = NULL,
                            lower.alpha = NULL, upper.alpha = NULL,
                            nu = 1.5, # smoothness parameter is fixed for now
                            Lmax = 30,
                            verbose = TRUE,
                            digits = 3){
  # some housekeeping
  N = length(y)
  # MCMC parameters
  if(is.null(niter)){
    niter = 5e3
    nburn = niter/2
    report = 1e2
  }

  if(is.null(nu)){
    nu = 1.5
    print("Fitting with nu = 1.5...")
  }

  # learning rates
  e.v = e.rho = e.alpha = 1e-1

  # Default hyper-parameter settings
  if(is.null(a.v)) a.v = 2
  if(is.null(b.v)) b.v = 2
  if(is.null(lower.alpha)) lower.alpha = 1e-2
  if(is.null(upper.alpha)) upper.alpha = 30
  lower.rho = 0
  upper.rho = 1

  # initial values
  alpha = v = 1
  rho = 0.5

  # batch acceptance probabilities
  accept.p = rep(0, 3)

  # initialize storage
  res_v = res_alpha = res_rho = rep(0, niter)
  accept_m = c()

  # Compute polynomials once
  Pl = legendre_Pl_list(coords = coords, Lmax = Lmax)

  In = diag(N)
  R = legendreM(nu = nu, alpha = alpha, P_l = Pl)

  Sigma = v * (rho * R + (1 - rho) * In)
  chol.Sigma = chol(Sigma)
  inv.Sigma = chol2inv(chol.Sigma)

  # if(trgtFn.compute){
  #   trgtFn.init = -(a.sigma + 1) * log(sig2) - b.sigma/sig2 -
  #     -(a.tau + 1) * log(tau2) - b.tau/tau2 -
  #     sum(log(diag(chol.Sigma)))/4 - t(y - Xb) %*% inv.Sigma %*% (y - Xb)/2
  # }

  for(i in 1:niter){

    ############
    # update v #
    ############
    lv.draw = log(v) + e.v * rnorm(1) # log-normal proposal

    Sigma.draw = exp(lv.draw) * (rho * R + (1 - rho) * In)
    chol.Sigma.draw = chol(Sigma.draw)
    inv.Sigma.draw = chol2inv(chol.Sigma.draw)

    r = -b.v * (exp(-lv.draw) - 1/v) - (a.v + 1) * (lv.draw - log(v)) -
      crossprod(t(crossprod(y, (inv.Sigma.draw - inv.Sigma))), y)/2 -
      sum(log(diag(chol.Sigma.draw)/diag(chol.Sigma)))

    accept.prob = min(r, 0)

    if(log(runif(1)) < accept.prob){
      v = exp(lv.draw)
      res_v[i] = v

      Sigma = Sigma.draw
      chol.Sigma = chol.Sigma.draw
      inv.Sigma = inv.Sigma.draw

      accept.p[1] = accept.p[1] + 1
    }else{
      res_v[i] = v
    }

    ##############
    # update rho #
    ##############
    lrho.draw = log(rho) + e.rho * rnorm(1) # log-normal proposal

    if((exp(lrho.draw) <= lower.rho) | (exp(lrho.draw) >= upper.rho)) accept.prob = -Inf
    else{

      Sigma.draw = v * (exp(lrho.draw) * R + (1 - exp(lrho.draw)) * In)
      chol.Sigma.draw = chol(Sigma.draw)
      inv.Sigma.draw = chol2inv(chol.Sigma.draw)

      r = (lrho.draw - log(rho)) -
        crossprod(t(crossprod(y, (inv.Sigma.draw - inv.Sigma))), y)/2 -
        sum(log(diag(chol.Sigma.draw)/diag(chol.Sigma)))

      accept.prob = min(r, 0)

    }


    if(log(runif(1)) < accept.prob){
      rho = exp(lrho.draw)
      res_rho[i] = rho

      Sigma = Sigma.draw
      chol.Sigma = chol.Sigma.draw
      inv.Sigma = inv.Sigma.draw

      accept.p[2] = accept.p[2] + 1
    }else{
      res_rho[i] = rho
    }


    ################
    # update alpha #
    ################

    lalpha.draw =  log(alpha) + e.alpha * rnorm(1)

    if((exp(lalpha.draw) < lower.alpha) | (exp(lalpha.draw) > upper.alpha)) accept.prob = -Inf
    else{

      R.draw = legendreM(nu = nu, alpha = exp(lalpha.draw), P_l = Pl)

      Sigma.draw = v * (rho * R.draw + (1 - rho) * In)
      chol.Sigma.draw = chol(Sigma.draw)
      inv.Sigma.draw = chol2inv(chol.Sigma.draw)

      r = (lalpha.draw - log(alpha)) - crossprod(t(crossprod(y, (inv.Sigma.draw - inv.Sigma))), y)/2 -
        sum(log(diag(chol.Sigma.draw)/diag(chol.Sigma)))

      accept.prob = min(r, 0)

    }


    if(log(runif(1)) < accept.prob){
      alpha = exp(lalpha.draw)
      res_alpha[i] = alpha

      R = R.draw
      chol.Sigma = chol.Sigma.draw
      inv.Sigma = inv.Sigma.draw

      accept.p[3] = accept.p[3] + 1
    }else{
      res_alpha[i] <- alpha
    }


    # if(trgtFn.compute){
    #   trgtFn[i] = -(a.sigma + 1) * log(sig2) - b.sigma/sig2 -
    #     -(a.tau + 1) * log(tau2) - b.tau/tau2 -
    #     sum(log(diag(chol.Sigma)))/4 - t(y - Xb) %*% inv.Sigma %*% (y - Xb)/2
    # }

    ############@@@@@@@@@@@@@@@##############
    # Adaptive Scaling of Proposal Variance #
    #     MH: optimal scales (33%)          #
    ############@@@@@@@@@@@@@@@##############
    if(i %% report == 0){
      accept.p = accept.p/report
      accept_m = rbind(accept_m, accept.p)
      if(verbose){
        # browser()
        if(i > nburn){
          cat("Iteration::", i, "Acceptance:", accept.p, "Tuning:", round(c(e.v, e.rho, e.alpha), digits), "\n",
              "=+++++++++++++++++++++++++++++++++++++++++++++++++=", "\n",
              "V::", round(median(res_v[nburn:i]), digits = (digits - 1)), "(", round(quantile(res_v[nburn:i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_v[nburn:i], probs = 0.975), digits = (digits - 1)), ")", "\t",
              "Rho::", round(median(res_rho[nburn:i]), digits = (digits - 1)), "(", round(quantile(res_rho[nburn:i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_rho[nburn:i], probs = 0.975), digits = (digits - 1)), ")", "\t",
              "Alpha::", round(median(res_alpha[nburn:i]), digits = (digits - 1)), "(", round(quantile(res_alpha[nburn:i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_alpha[nburn:i], probs = 0.975), digits = (digits - 1)), ")", "\n",
              "=+++++++++++++++++++++++++++++++++++++++++++++++++=", "\n")
        }else{
          cat("Iteration::", i, "Acceptance:", accept.p, "Tuning:",round(c(e.v, e.rho, e.alpha), digits), "\n",
              "=+++++++++++++++++++++++++++++++++++++++++++++++++=", "\n",
              "V::", round(median(res_v[(i - report + 1):i]), digits = (digits - 1)), "(", round(quantile(res_v[(i - report + 1):i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_v[(i - report + 1):i], probs = 0.975), digits = (digits - 1)), ")", "\t",
              "Rho::", round(median(res_rho[(i - report + 1):i]), digits = (digits - 1)), "(", round(quantile(res_rho[(i - report + 1):i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_rho[(i - report + 1):i], probs = 0.975), digits = (digits - 1)), ")", "\t",
              "Alpha::", round(median(res_alpha[(i - report + 1):i]), digits = (digits - 1)), "(", round(quantile(res_alpha[(i - report + 1):i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_alpha[(i - report + 1):i], probs = 0.975), digits = (digits - 1)), ")", "\n",
              "=+++++++++++++++++++++++++++++++++++++++++++++++++=", "\n")
        }
      }
      ####################
      # Proposal Scaling #
      ####################
      step = sapply(accept.p, function(x){
        out = 1
        x = max(0.17, min(x, 0.75))
        if(x > 0.50) out = x/0.50
        else if(x < 0.25) out = x/0.25
        out
      })
      e.v = e.v * step[1]
      e.rho = e.rho * step[2]
      e.alpha = e.alpha * step[3]

      accept.p = rep(0, 3)
    }
  }

  return(list(y = y,
              coords = coords,
              # trgt_fn = c(trgtFn.init, trgtFn),
              nu = nu,
              sig2 = res_v[(nburn + 1):niter] * res_rho[(nburn + 1):niter],
              tau2 =  res_v[(nburn + 1):niter] * (1 - res_rho[(nburn + 1):niter]),
              alpha = res_alpha[(nburn + 1):niter]))
}
