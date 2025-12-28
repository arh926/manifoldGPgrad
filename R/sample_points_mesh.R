#' Scattered sample locations on a mesh
#'
#' Obtains a sample of scattered locations on a triangulear mesh.
#'
#' @param mesh a mesh in polygon file format imported using `vcgPlyRead()` of class `mesh3d`.
#' @param npoints number of samples points.
#' @returns A list containing a matrix of npoints x 3 points, IDs of vertices forming triangles the barycentric coordinates
#' @importFrom stats runif
#' @author Aritra Halder <aritra.halder@drexel.edu>
#' @export

sample_points_mesh <- function(mesh = NULL,
                               npoints = NULL) {
  V = t(mesh$vb[1:3, ])
  F = t(mesh$it)

  areas = triangle_areas(V, F)
  probs = areas / sum(areas)

  # sample triangles by area
  tri_idx = sample.int(nrow(F), size = npoints, replace = TRUE, prob = probs)

  # barycentric sampling for each chosen triangle
  u = runif(npoints)
  v = runif(npoints)
  # ensure u+v <= 1
  mask = (u + v > 1)
  u[mask] = 1 - u[mask]
  v[mask] = 1 - v[mask]
  w = 1 - u - v

  # coordinates of sampled points
  v1 = V[F[tri_idx, 1], ]
  v2 = V[F[tri_idx, 2], ]
  v3 = V[F[tri_idx, 3], ]

  pts = u * v1 + v * v2 + w * v3  # n_points x 3

  list(points   = pts,
       tri_idx  = tri_idx,
       bary     = cbind(u, v, w))
}
