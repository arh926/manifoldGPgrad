#' Interpolate on compact Riemanian manifolds
#'
#' Function to smoothly interpolate realizations over the mesh using either heat kernel or surface splines.
#' @param mesh a mesh in polygon file format imported using `vcgPlyRead()` of class `mesh3d`.
#' @param gp_sample a list containing funciton values, IDs of triangle vertices for observed locations, barycentric coordinates (obtained from `sample_points_mesh()`)
#' @param tau_noise regularization parameter for heat-kernel (should be small defaults to 0.01)
#' @param lambda_smooth penalty parameter for splinel (should be small defaults to 0.001)
#' @param method choice between: heat kernel and spline
#' @returns a vector of interpolated realizations on the mesh
#' @keywords interpolate_mesh
#' @author Aritra Halder <aritra.halder@drexel.edu>
#' @export
interpolate_mesh <- function(mesh = NULL,
                             gp_sample = NULL,
                             tau_noise = 0.01, # set for heat kernel
                             lambda_smooth = 1e-2, # set for spline
                             method = c("heat kernel", "spline")){
  F = t(mesh$it)

  L_cot = laplacian_cotangent(mesh)
  L = (L_cot + t(L_cot)) / 2

  n = nrow(L)

  v1_idx = F[gp_sample$tri_idx,1]
  v2_idx = F[gp_sample$tri_idx,2]
  v3_idx = F[gp_sample$tri_idx,3]

  # These scatter values to mesh vertices
  f_vert = numeric(n)
  wt_vert = numeric(n)

  for (i in seq_along(gp_sample$values)) {
    f_val = gp_sample$values[i]
    u = gp_sample$bary[i,1]
    v = gp_sample$bary[i,2]
    w = gp_sample$bary[i,3]

    f_vert[v1_idx[i]] = f_vert[v1_idx[i]] + u * f_val
    f_vert[v2_idx[i]] = f_vert[v2_idx[i]] + v * f_val
    f_vert[v3_idx[i]] = f_vert[v3_idx[i]] + w * f_val

    wt_vert[v1_idx[i]] = wt_vert[v1_idx[i]] + u
    wt_vert[v2_idx[i]] = wt_vert[v2_idx[i]] + v
    wt_vert[v3_idx[i]] = wt_vert[v3_idx[i]] + w
  }

  # Avoid division by zero
  obs_mask = wt_vert > 0
  f_vert[obs_mask] = f_vert[obs_mask] / wt_vert[obs_mask]

  if(method == "heat kernel"){
    u_hat = solve(diag(n) + tau_noise * L, f_vert)
  }else if(method == "spline"){
    W = diag(wt_vert)
    A = lambda_smooth * L + W
    rhs = W %*% f_vert

    u_hat = solve(A, rhs)
  }
  return(u_hat = u_hat)
}
