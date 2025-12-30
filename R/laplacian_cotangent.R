#' Cotangent Laplacian
#'
#' Function to compute a sparse cotangent Laplacian that approximates the Laplace-Beltrami operator for a compact Riemannian manifold.
#'
#' @param mesh a mesh in polygon file format imported using `vcgPlyRead()` of class `mesh3d`.
#' @returns a sparse cotangent Laplacian matrix.
#' @importFrom Matrix sparseMatrix Diagonal t rowSums
#' @keywords laplacian_cotangent
#' @author Aritra Halder <aritra.halder@drexel.edu>
#' @export
laplacian_cotangent <- function(mesh = NULL){
  V = t(mesh$vb[1:3, ])
  F = t(mesh$it)              # m x 3, 1-based indices
  n = nrow(V)
  m = nrow(F)

  # we'll accumulate off-diagonal weights in sparse form;
  # duplicates (same (i,j)) get summed automatically
  I = integer(0)
  J = integer(0)
  X = numeric(0)

  for (f in seq_len(m)) {
    idx = F[f, ]
    i = idx[1]; j = idx[2]; k = idx[3]

    vi = V[i, ]; vj = V[j, ]; vk = V[k, ]

    # angles at each vertex of the triangle
    alpha = angle_at(vj, vi, vk)  # angle at i
    beta  = angle_at(vi, vj, vk)  # angle at j
    gamma = angle_at(vi, vk, vj)  # angle at k

    # cotangents
    c_alpha = cot(alpha)
    c_beta  = cot(beta)
    c_gamma = cot(gamma)

    # Each triangle contributes cot(angle_opposite_edge)/2 to that edge
    # contributions to edge (j,k) from angle at i:
    w_jk = 0.5 * c_alpha
    # edge (i,k) from angle at j:
    w_ik = 0.5 * c_beta
    # edge (i,j) from angle at k:
    w_ij = 0.5 * c_gamma

    # For symmetry, add both (p,q) and (q,p) entries
    I = c(I, j, k, k, i, i, j)
    J = c(J, k, j, i, k, j, i)
    X = c(X, w_jk, w_jk, w_ik, w_ik, w_ij, w_ij)
  }

  W = sparseMatrix(i = I, j = J, x = X, dims = c(n, n))
  # ensure symmetry
  W = 0.5 * (W + t(W))

  d = rowSums(W)
  D = Diagonal(n, d)

  L = D - W
  L
}
