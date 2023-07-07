#' @title Apollonius diagram and Apollonius graph
#' @description Computation of the Apollonius diagram and the Apollonius graph
#'   of some weighted 2D points. The Apollonius graph is the dual of the
#'   Apollonius diagram. It is also called the additively weighted Voronoï
#'   diagram.
#'
#' @param sites the 2D points, a numeric matrix with two columns (one point
#'   per row)
#' @param radii the weights, a numeric vector of length equal to the number of
#'   points (i.e. the number of rows of \code{sites})
#' @param t0 a positive parameter to increase in case if some infinite edges
#'   of the Apollonius graph are in the wrong direction
#' @param tmax a positive number passed to \code{\link[gyro]{gyroray}},
#'   controlling the length of the infinite edges (i.e. the hyperbolic rays)
#'   of the Apollonius graph
#' @param nsegs a positive integer, the desired number of points of each
#'   finite edge of the Apollonius graph
#' @param nrays a positive integer, the desired number of points of each
#'   infinite edge of the Apollonius graph
#'
#' @return A list with two fields \code{diagram} and \code{graph}. The
#'   \code{diagram} field is a list providing the sites and the faces of the
#'   Apollonius diagram. The \code{graph} field is a list providing the sites
#'   and the edges of the Apollonius graph.
#' @export
#'
#' @importFrom gyro gyrosegment gyroray gyroABt
#' @importFrom abind abind
#'
#' @examples
#' library(Apollonius)
#' sites <- rbind(
#'   c(0, 0),
#'   c(4, 1),
#'   c(2, 4),
#'   c(7, 4),
#'   c(8, 0),
#'   c(5, -2),
#'   c(-4, 4),
#'   c(-2, -1),
#'   c(11, 4),
#'   c(11, 0)
#' )
#' radii <- c(1, 1.5, 1.25, 2, 1.75, 0.5, 0.4, 0.6, 0.7, 0.3)
#' apo <- Apollonius(sites, radii)
#' opar <- par(mar = c(4, 4, 1, 1))
#' plotApolloniusGraph(apo, xlab = "x", ylab = "y")
#' par(opar)
Apollonius <- function(
    sites, radii, t0 = 2, tmax = 10, nsegs = 100L, nrays = 300L
) {
  stopifnot(
    is.numeric(sites), is.matrix(sites), ncol(sites) == 2L, nrow(sites) >= 3L
  )
  storage.mode(sites) <- "double"
  if(anyNA(sites)) {
    stop("Missing values are not allowed.")
  }
  if(anyDuplicated(sites)) {
    stop("Found duplicated sites.")
  }
  stopifnot(is.numeric(radii), length(radii) == nrow(sites), all(radii != 0))
  storage.mode(radii) <- "double"
  if(anyNA(radii)) {
    stop("Found missing value(s) in `radii`.")
  }
  stopifnot(t0 > 1)
  stopifnot(tmax > 1)
  #
  stuff <- ApolloniusCpp(sites, radii)
  neighbors <- stuff[["neighbors"]]
  if(nrow(neighbors) == 0L) {
    stop("The Apollonius diagram is empty.")
  }
  #
  Neighs <- vector("list", nrow(neighbors))
  for(i in 1L:nrow(neighbors)) {
    Neighs_i <- 1L:3L
    neighs_i <- neighbors[i, ]
    remove <- which(is.na(neighs_i))
    for(j in seq_len(i-1L)) {
      if(j %in% neighs_i && i %in% neighbors[j, ]) {
        remove <- c(remove, which(neighs_i == j))
      }
    }
    if(length(remove) > 0L) {
      Neighs_i <- Neighs_i[-remove]
    }
    Neighs[[i]] <- Neighs_i
  }
  #
  commonVertices <-
    abind(stuff[["cvertex1"]], stuff[["cvertex2"]], along = 3L)
  #
  dpoints <- stuff[["dpoints"]]
  vertices <- stuff[["vertices"]]
  #
  hsegments <- vector("list", sum(lengths(Neighs)))
  h <- 1L
  for(i in 1L:nrow(dpoints)) {
    P1 <- dpoints[i, ]
    infinite <- !is.na(P1[3L])
    Neighs_i <- Neighs[[i]]
    vert_i <- vertices[[i]]
    for(k in seq_along(Neighs_i)) {
      vs <- commonVertices[i, Neighs_i[k], ]
      A  <- vert_i[vs[1L], 1L:2L]
      rA <- vert_i[vs[1L], 3L]
      B  <- vert_i[vs[2L], 1L:2L]
      rB <- vert_i[vs[2L], 3L]
      ctr <- (A + B)/2
      P2 <- dpoints[neighbors[i, Neighs_i[k]], c(1L, 2L)]
      if(infinite) {
        AB <- sqrt(c(crossprod(B-A)))
        u <- (B-A) / AB
        P <- A + (rA + (AB - (rA + rB))/2) * u
      } else {
        P <- P1[c(1L, 2L)]
      }
      if(rA != rB) {
        f <- function(log_s) {
          d <- ctr + gyromidpoint(P-ctr, P2-ctr, exp(log_s))
          sqrt(c(crossprod(d-A))) - sqrt(c(crossprod(d-B))) - (rA - rB)
        }
        ur <- uniroot(f, lower = -5, upper = 2, extendInt = "yes")
        s <- exp(ur[["root"]])
      }
      if(infinite) {
        P_is_up <- P1[1L]*P[1L] + P1[2L]*P[2L] > P1[3L]
        if(rA == rB) {
          X <- P2 + t0 * (P - P2)
        } else {
          X <- ctr + gyroABt(P2-ctr, P-ctr, t = t0, s = s)
        }
        X_is_up <- P1[1L]*X[1L] + P1[2L]*X[2L] > P1[3L]
        reverse <- P_is_up != X_is_up
        if(rA == rB) {
          hsegments[[h]] <- ray(P2, P, OtoA = !reverse, tmax = tmax, n = nrays)
        } else {
          hsegments[[h]] <-
            t(ctr + t(gyroray(
              P2-ctr, P-ctr, s = s, OtoA = !reverse, tmax = tmax, n = nrays
            )))
        }
      } else {
        if(rA == rB) {
          hsegments[[h]] <- segment(P, P2, n = nsegs)
        } else {
          hsegments[[h]] <- t(ctr + t(gyrosegment(
            P-ctr, P2-ctr, s = s, n = nsegs
          )))
        }
      }
      h <- h + 1L
    }
  }
  wsites <- cbind(sites, radii)
  colnames(wsites) <- c("x", "y", "weight")
  dsites <- dpoints[is.na(dpoints[, 3L]), ]
  list(
    "diagram" = list("sites" = wsites, "faces" = stuff[["faces"]]),
    "graph"   = list("sites" = dsites, "edges" = hsegments)
  )
}