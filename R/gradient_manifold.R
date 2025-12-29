#' Posterior inference on gradients
#'
#' Posterior Markov Chain Monte Carlo (MCMC) samples for gradients along a pre-supplied vector field 1-1 posterior samples of \eqn{\sigma^2} and \eqn{\alpha}. If vector field is not supplied it uses the local \eqn{x}-axis for vertices to generate a smooth vector field.
#'
#' @param mesh a mesh in polygon file format imported using `vcgPlyRead()` of class `mesh3d`.
#' @param coords list containing barycentric coordinates of \eqn{N} sampled irregular locations, IDs for triangle vertices
#' @param grid grid locations generated using `farthest_point_subsample()`
#' @param lambda truncated vector of eigen-values
#' @param phi truncated \eqn{N\times T} matrix of eigen functions
#' @param model a list containing posterior post burn-in samples of \eqn{\sigma^2} and \eqn{\alpha}
#' @param d optional; manifold dimension, defaults to 2
#' @param Vfield optional; vector field on the mesh
#' @returns A list of the posterior samples for the gradients and a list containing: grid points, IDs of triangle vertices and barycentric coordinates and poosterior estimate of gradients.
#' @importFrom Rvcg vcgUpdateNormals
#' @author Aritra Halder <aritra.halder@drexel.edu>
#' @export
gradient_manifold<- function(mesh = NULL,
                             coords = NULL,
                             grid = NULL,
                             model = NULL,
                             lambda = NULL,
                             phi = NULL,
                             d = 2, # manifold dimension
                             Vfield = NULL){ # should be n_triangles x 3
  if(is.null(Vfield)){
    cat("Using x-axis as reference to generate vector field", "\n")
    Vcoords = t(mesh$vb[1:3, ])
    F = t(mesh$it)
    ref = c(1, 0, 0) # x-axis

    # Tangent vector field
    mesh_norm  = vcgUpdateNormals(mesh)
    V_raw = t(mesh_norm$normals)[, 1:3]
    V_vert = matrix(0, nrow(V_raw), 3)

    for(i in seq_len(nrow(Vcoords))){
      n = V_raw[i,]

      v_tan = ref - sum(ref * n) * n
      v_tan = v_tan/sqrt(sum(v_tan^2))

      V_vert[i,] = v_tan
    }
  }else{
    V_vert = Vfield
  }

  # Note:: replace this later
  # barycentric coordinates for the sample
  uX = coords$bary[, 1]
  vX = coords$bary[, 2]
  wX = coords$bary[, 3]

  # vertices for sample
  v1_idx_X = F[coords$tri_idx, 1]
  v2_idx_X = F[coords$tri_idx, 2]
  v3_idx_X = F[coords$tri_idx, 3]

  # eigen-functions at vertices:: f_l(v_k)
  phi_v1_X = phi[v1_idx_X, , drop = FALSE]  # N x k
  phi_v2_X = phi[v2_idx_X, , drop = FALSE]
  phi_v3_X = phi[v3_idx_X, , drop = FALSE]

  # interpolate eigen-functions:: f_l(x_i)
  Phi_X = uX * phi_v1_X + vX * phi_v2_X + wX * phi_v3_X  # N x k

  N = nrow(coords$bary) # sample size
  k = length(lambda) # discrete approx.
  nmcmc = length(model$alpha) # number of samples
  ngrid = nrow(grid$points)

  # posterior estimates of parameters
  alpha.mc = model$alpha
  z.mc = model$z
  sig2.mc = model$sig2
  nu = model$nu

  grad.samples.grid = matrix(NA, nmcmc, ngrid)

  for(i.mcmc in 1:nmcmc){
    R = matern_cov_mp(mesh = mesh, samp = coords,
                      lambda = lambda, phi = phi,
                      nu = nu, alpha = alpha.mc[i.mcmc], d = d,
                      cor = FALSE)

    R = R + diag(nrow(R)) # affected due to scaling
    R.inv = chol2inv(chol(R))

    for (j in 1:ngrid){
      idx_f = F[grid$tri_idx[j], ]
      p1 = Vcoords[idx_f[1], ]
      p2 = Vcoords[idx_f[2], ]
      p3 = Vcoords[idx_f[3], ]

      e1   = p2 - p1
      e2   = p3 - p1
      nvec = cross(e1, e2)
      nvec = nvec / sqrt(sum(nvec^2))  # unit normal

      area2 = sqrt(sum(cross(e1, e2)^2))

      # standard FEM barycentric gradients:
      g1 = cross(nvec, p3 - p2) / area2
      g2 = cross(nvec, p1 - p3) / area2
      g3 = cross(nvec, p2 - p1) / area2

      # barycentric coords at x'
      u_p = grid$bary[j, 1]
      v_p = grid$bary[j, 2]
      w_p = grid$bary[j, 3]

      v1_idx_p = idx_f[1]
      v2_idx_p = idx_f[2]
      v3_idx_p = idx_f[3]

      # vector field at x' (vertex normals interpolated)
      Vxprime = u_p * V_vert[v1_idx_p, ] +
        v_p * V_vert[v2_idx_p, ] +
        w_p * V_vert[v3_idx_p, ]

      # directional derivatives of bary basis along V at x'
      a1 = sum(g1 * Vxprime)   # D_V phi1
      a2 = sum(g2 * Vxprime)   # D_V phi2
      a3 = sum(g3 * Vxprime)   # D_V phi3

      # eigenfunctions at triangle vertices (for all k modes)
      phi_v1_p = phi[v1_idx_p, ]  # length k
      phi_v2_p = phi[v2_idx_p, ]
      phi_v3_p = phi[v3_idx_p, ]

      # interpolated:: D_V f_l(x')
      D_V_f_xprime = a1 * phi_v1_p + a2 * phi_v2_p + a3 * phi_v3_p  # length k

      ## 3) spectral weights
      expo = nu + d/2
      w_spec = (alpha.mc[i.mcmc]^2 + lambda)^(-expo)   # length k

      ## no scaling
      wD  = w_spec * D_V_f_xprime              # length k
      cov_vec_raw = as.numeric(Phi_X %*% wD)    # length N

      KV_xx = sum(w_spec * D_V_f_xprime^2)

      # no scaling anywhere
      mu.grad = as.numeric(crossprod(cov_vec_raw, R.inv %*% z.mc[i.mcmc,]))
      Sigma.grad = sig2.mc[i.mcmc] * (KV_xx - crossprod(cov_vec_raw,
                                                        R.inv %*% cov_vec_raw))

      if (Sigma.grad < 0) {
        cat("Bad Iteration...", "\t")
        grad.samples.grid[i.mcmc, j] = rnorm(1, mu.grad, 1)
      } else {
        grad.samples.grid[i.mcmc, j] = rnorm(1, mu.grad, sqrt(Sigma.grad))
      }
    }
    if(i.mcmc %% 100 == 0){
      if(i.mcmc == 100) cat("Iteration", i.mcmc, "...", "\t")
      else cat(i.mcmc, "...", "\t")
    }
  }
  # Posterior inference
  grad_ci =  t(apply(grad.samples.grid, 2, quantile, probs = c(0.5, 0.025, 0.975)))

  # Gradient Object
  gp_grad_est_obj = list(points = grid$points,
                         values = grad_ci[,1],
                         tri_idx = grid$tri_idx,
                         bary = grid$bary)

  return(list(grad.mc = grad.samples.grid,
              grad_ci = grad_ci,
              gp_grad_est_obj = gp_grad_est_obj))
}
