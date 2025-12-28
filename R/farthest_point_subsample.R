#' Farthest Point Sampling for Point Clouds
#'
#' Generates a point cloud resembling a grid for the mesh.
#'
#' @param coords a list containing a dense sample of points, IDs for triangle vertices and barycentric coordinates using `sample_points_mesh()`
#' @param nkeep number of points required in the grid
#' @returns a list of sample of points, IDs for triangle vertices and barycentric coordinates of length nkeep
#' @author Aritra Halder <aritra.halder@drexel.edu>
#' @export
farthest_point_subsample <- function(coords, nkeep) {

  n = nrow(coords$points)
  if (nkeep >= n) return(coords$points)

  # start with a random point
  idx_sel = integer(nkeep)
  idx_sel[1] = sample.int(n, 1)

  # squared distances to first selected point
  diff = sweep(coords$points, 2, coords$points[idx_sel[1],], "-")
  d2   = rowSums(diff^2)

  for (i in 2:nkeep) {
    # pick point farthest from current selected set
    idx_sel[i] = which.max(d2)

    # update distances: min dist to any selected so far
    diff = sweep(coords$points, 2, coords$points[idx_sel[i],], "-")
    d2   = pmin(d2, rowSums(diff^2))
  }

  # points[idx_sel, , drop = FALSE]

  return(list(points   = coords$points[idx_sel, , drop = FALSE],
              tri_idx  = coords$tri_idx[idx_sel, drop = FALSE],
              bary     = coords$bary[idx_sel, , drop = FALSE]))
}
