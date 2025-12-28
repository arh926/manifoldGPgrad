#' Angle at Vertex
#'
#' Compute angle at vertex b given triangle (a, b, c) with 3D coords
#'
#' @param a vertex a
#' @param b vertex b
#' @param c vertex c
#' @returns angle at vertex b given triangle (a, b, c) with 3D coords. For internal use only.
#' @keywords angle_at
#' @author Aritra Halder <aritra.halder@drexel.edu>
#' @export
angle_at <- function(a = NULL,
                     b = NULL,
                     c = NULL) {
  v1 = a - b
  v2 = c - b
  cos_theta = sum(v1 * v2) / (sqrt(sum(v1^2)) * sqrt(sum(v2^2)))
  # numerical clamp
  cos_theta = max(min(cos_theta, 1), -1)

  acos(cos_theta)
}
