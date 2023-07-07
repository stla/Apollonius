library(plotrix)
library(gyro)

# gyroray <- function(O, A, s, n = 300, tmax = 20, OtoA = TRUE) {
#   if(OtoA) {
#     t_ <- seq(0, tmax, length.out = n)
#   } else {
#     t_ <- seq(-tmax, 0, length.out = n)
#   }
#   t(vapply(t_, function(t) {
#     gyro:::UgyroABt(O, A, t, s)
#   }, numeric(length(A))))
# }

segment <- function(A, B, n) {
  t_ <- seq(0, 1, length.out = n)
  t(vapply(t_, function(t) {A + t*(B-A)}, numeric(2L)))
}

ray <- function(O, A, n, tmax, OtoA) {
  if(OtoA) {
    t_ <- seq(0, tmax, length.out = n)
  } else {
    t_ <- seq(-tmax, 0, length.out = n)
  }
  t(vapply(t_, function(t) {O + t*(A-O)}, numeric(2L)))
}

p1  <- c(0, 0)
p2  <- c(4, 1)
p3  <- c(2, 4)
p4  <- c(7, 4)
p5  <- c(8, 0)
p6  <- c(5, -2)
p7  <- c(-4, 4)
p8  <- c(-2, -1)
p9  <- c(11, 4)
p10 <- c(11, 0)

sites <- rbind(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10)

sites <- rbind(
  c(0, 0),
  c(4, 1),
  c(2, 4),
  c(7, 4),
  c(8, 0),
  c(5, -2),
  c(-4, 4),
  c(-2, -1),
  c(11, 4),
  c(11, 0)
)
radii <- c(1, 1.5, 1.25, 2, 1.75, 0.5, 0.4, 0.6, 0.7, 0.3)


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
#' opar <- par(mar = c(3, 3, 1, 1))
#' plotApolloniusGraph(apo)
#' par(opar)
Apollonius <- function(
  sites, radii, t0 = 2, tmax = 10, nsegs = 100L, nrays = 300L
) {
  stuff <- Apollonius:::ApolloniusCpp(sites, radii)
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
    abind::abind(stuff[["cvertex1"]], stuff[["cvertex2"]], along = 3L)
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
      f <- function(log_s) {
        d <- ctr + gyromidpoint(P-ctr, P2-ctr, exp(log_s))
        sqrt(c(crossprod(d-A))) - sqrt(c(crossprod(d-B))) - (rA - rB)
      }
      ur <- uniroot(f, lower = -5, upper = 2, extendInt = "yes")
      s <- exp(ur[["root"]])
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
          hsegments[[h]] <- segment(P, P2,n = nsegs)
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


clrs <- randomcoloR::distinctColorPalette(10L)

plotApolloniusGraph <- function(apo, limits = NULL) {
  sites  <- apo[["diagram"]][["sites"]]
  nsites <- nrow(sites)
  radii  <- sites[, "weight"]
  dsites <- apo[["graph"]][["sites"]]
  edges  <- apo[["graph"]][["edges"]]
  #
  #clrs <- randomcoloR::distinctColorPalette(nsites)
  #
  if(is.null(limits)) {
    x <- grDevices::extendrange(sites[, "x"])
    y <- grDevices::extendrange(sites[, "y"])
    limits <- c(min(x[1L], y[1L]), max(x[2L], y[2L]))
  }
  #
  plot(NULL, xlim = limits, ylim = limits, asp = 1, xlab = "x", ylab = "y")
  for(i in 1L:nsites) {
    draw.circle(
      sites[i, "x"], sites[i, "y"], radius = radii[i],
      border = clrs[i], col = clrs[i]
    )
  }
  points(dsites, pch = 19)
  for(i in 1L:length(edges)) {
    lines(edges[[i]], col="black", lwd = 2)
  }
  invisible()
}

# sites <- rbind(
#   c(0, 0),
#   c(10.5, 0),
#   c(20, 0),
#   c(0, 9.5),
#   c(10.5, 9.5),
#   c(20, 9.5),
#   c(0, 16),
#   c(10.5, 16),
#   c(20, 16)
# )
# radii <- seq(5, 1, by = -0.5)

apo <- Apollonius(sites, radii)

# svg("x.svg", width = 8, height = 4)
opar <- par(mar = c(3, 3, 1, 1))
plotApolloniusGraph(apo, limits = c(-5, 25))
par(opar)
# dev.off()
#
# rsvg::rsvg_png("x.svg", "inst/screenshots/agraph02.png", width = 512, height = 256)
# file.remove("x.svg")


for(i in 1L:10L) {
  apo <- Apollonius(sites, radii)
  svg("x.svg", width = 8, height = 4)
  opar <- par(mar = c(3, 3, 1, 1))
  plotApolloniusGraph(apo)#, limits = c(-5, 25))
  par(opar)
  dev.off()
  rsvg::rsvg_png("x.svg", sprintf("zzpic%03d.png", i), width = 512, height = 256)
  radii <- radii * 0.95
}

library(gifski)
pngs <- Sys.glob("zzpic*.png")
gifski(
  png_files = c(pngs, rev(pngs)),
  gif_file = "agraph03.gif",
  width = 512,
  height = 256,
  delay = 1/8
)
file.remove(pngs)
