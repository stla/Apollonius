library(plotrix)
library(gyro)

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
radii <- c(1, 1.5, 1.25, 2, 1.75, 0.5, 0.4, 0.6, 0.7, 0.3)


Apollonius <- function(sites, radii) {
  stuff <- Apollonius:::test(sites, radii)
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
    Neighs_i <- Neighs[[i]]
    vert_i <- vertices[[i]]
    for(k in seq_along(Neighs_i)) {
      vs <- commonVertices[i, Neighs_i[k], ]
      A  <- vert_i[vs[1L], 1L:2L]
      rA <- vert_i[vs[1L], 3L]
      B  <- vert_i[vs[2L], 1L:2L]
      rB <- vert_i[vs[2L], 3L]
      ctr <- (A + B)/2
      P2 <- dpoints[neighbors[i, Neighs_i[k]], ]
      f <- function(log_s) {
        d <- ctr + gyromidpoint(P1-ctr, P2-ctr, exp(log_s))
        sqrt(c(crossprod(d-A))) - sqrt(c(crossprod(d-B))) - (rA - rB)
      }
      ur <- uniroot(f, lower = -5, upper = 2, extendInt = "yes")
      s <- exp(ur$root)
      hsegments[[h]] <- t(ctr + t(gyrosegment(P1-ctr, P2-ctr, s = s)))
      h <- h + 1L
    }
  }
  wsites <- cbind(sites, radii)
  colnames(wsites) <- c("x", "y", "weight")
  list(
    "diagram" = list("sites" = wsites),
    "graph"   = list("sites" = dpoints, "edges" = hsegments)
  )
}




stuff <- Apollonius(sites, radii)
sites <- stuff[["diagram"]][["sites"]]
radii <- sites[, 3L]
dsites <- stuff[["graph"]][["sites"]]
edges <- stuff[["graph"]][["edges"]]

clrs <- randomcoloR::distinctColorPalette(nrow(sites))

# svg("x.svg", width = 8, height = 4)
opar <- par(mar = c(3, 3, 1, 1))
plot(NULL, xlim = c(-1, 10), ylim = c(-3, 6), asp = 1, xlab = "x", ylab = "y")
for(i in 1L:nrow(sites)) {
  draw.circle(
    sites[i, 1L], sites[i, 2L], radius = radii[i],
    border = clrs[i], col = clrs[i], lwd = 4
  )
}
points(dsites, pch = 19)
for(i in 1L:length(edges)) {
  lines(edges[[i]], col="black", lwd = 2)
}
# dev.off()
#
# rsvg::rsvg_png("x.svg", "inst/screenshots/agraph02.png", width = 512, height = 256)
# file.remove("x.svg")
