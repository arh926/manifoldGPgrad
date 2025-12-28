#' Cross Product of Vectors
#'
#' Function for computing cross-product of two vectors.
#'
#' @param a \eqn{3 \times 1} vector.
#' @param b \eqn{3 \times 1} vector.
#' @returns a \eqn{3 \times 1} vector. For internal use only.
#' @keywords cross
#' @author Aritra Halder <aritra.halder@drexel.edu>
#' @export
cross <- function(a = NULL,
                  b = NULL) {
  c(a[2] * b[3] - a[3] * b[2],
    a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1])
}
