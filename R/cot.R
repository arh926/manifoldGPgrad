#' Cotangent Function
#'
#' Function to calculate cotangent of an angle
#'
#' @param theta angle
#' @returns cotangent of the angle. For internal use only.
#' @keywords cot
#' @author Aritra Halder <aritra.halder@drexel.edu>
#' @export
cot <- function(theta = NULL) {
  cos(theta) / sin(theta)
}
