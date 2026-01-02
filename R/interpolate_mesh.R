#' Surface registration on a compact Riemannian manifold from scattered data
#'
#' Function to smoothly interpolate realizations over the mesh using either heat kernel or surface splines. The interpolation is done for scattered or patially observed data on the manifold.
#'
#' @param mesh a mesh in polygon file format imported using `vcgPlyRead()` of class `mesh3d`.
#' @param sample a list containing function values, IDs of triangle vertices for observed locations, barycentric coordinates (obtained from `sample_points_mesh()`)
#' @param L_sym optional; the cotangent Laplacian, if not provided then, is calculated in the function but, takes longer.
#' @param tau_noise optional; regularization parameter for heat-kernel (should be small defaults to 0.01)
#' @param lambda_smooth optional; penalty parameter for spline (should be small defaults to 0.001)
#' @param method choice between: heat kernel and spline
#' @importFrom Matrix Diagonal Cholesky solve
#' @returns a vector of interpolated realizations on the mesh
#' @keywords interpolate_mesh
#' @author Aritra Halder <aritra.halder@drexel.edu>
#' @export
interpolate_mesh <- function(mesh = NULL, L_sym = NULL,
                             sample = NULL,
                             tau_noise = 0.01, # set for heat kernel
                             lambda_smooth = 1e-2, # set for spline
                             method = c("heat kernel", "spline")){
  F = t(mesh$it)

  if(is.null(L_sym)){
    L_cot = laplacian_cotangent(mesh)
    L_sym = (L_cot + t(L_cot)) / 2

    n = nrow(L_sym)
  }

  n = nrow(L_sym)

  v1_idx = F[sample$tri_idx,1]
  v2_idx = F[sample$tri_idx,2]
  v3_idx = F[sample$tri_idx,3]

  # These scatter values to mesh vertices
  f_vert = numeric(n)
  wt_vert = numeric(n)

  for (i in seq_along(sample$values)) {
    f_val = sample$values[i]
    u = sample$bary[i,1]
    v = sample$bary[i,2]
    w = sample$bary[i,3]

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
    A = Diagonal(n) + tau_noise * L_sym
    R = Cholesky(A)
    u_hat = solve(R, f_vert)
  }else if(method == "spline"){
    W = Diagonal(x = wt_vert)
    A = lambda_smooth * L_sym + W
    rhs = W %*% f_vert
    R = Cholesky(A)

    u_hat = solve(R, rhs)
  }
  return(u_hat = u_hat)
}
