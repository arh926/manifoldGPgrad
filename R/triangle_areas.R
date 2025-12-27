#' Area of a triangle
#'
#' Function to compute the area of a triangle in the mesh
#'
#' @param V \eqn{K \times 3} matrix containing homologous vertex coordinates
#' @param F \eqn{N_T \times 3} matrix containing vertex indices forming triangular faces
#' @returns a scalar for the area of a triangle. For internal use.
#' @keywords triangle_areas
#' @export
triangle_areas <- function(V, F) {
  # V: n x 3, F: m x 3
  v1 = V[F[,1], ]
  v2 = V[F[,2], ]
  v3 = V[F[,3], ]

  # cross product norm / 2
  a  = v2 - v1
  b  = v3 - v1
  cross = cross(a, b)

  0.5 * sqrt(rowSums(cross^2))
}
