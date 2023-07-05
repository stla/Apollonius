library(plotrix)
library(gyro)

p1 <- c(0, 0)
p2 <- c(4, 1)
p3 <- c(2, 4)
p4 <- c(7, 4)
p5 <- c(8, 0)
p6 <- c(5, -2)
p7 <- c(-4, 4)
p8 <- c(-2, -1)
p9 <- c(11, 4)
p10 <- c(11, 0)

sites <- rbind(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10)
radii <- c(1, 1.5, 1.25, 2, 1.75, 0.5, 0.4, 0.6, 0.7, 0.3)

clrs <- randomcoloR::distinctColorPalette(nrow(sites))

opar <- par(mar = c(3, 3, 1, 1))
plot(NULL, xlim = c(-1, 10), ylim = c(-3, 6), asp = 1, xlab = "x", ylab = "y")
for(i in 1L:nrow(sites)) {
  draw.circle(
    sites[i, 1L], sites[i, 2L], radius = radii[i],
    border = clrs[i], lwd = 4
  )
}

stuff <- Apollonius:::test(sites, radii)

vertices <- stuff[["vertices"]]

neighbors <- stuff[["neighbors"]]

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
Neighs

commonVertices <-
  abind::abind(stuff[["cvertex1"]], stuff[["cvertex2"]], along = 3L)
commonVertices[1, 3, ] # common vertices face1 with its neighbor 3

dpoints <- stuff[["dpoints"]]

points(dpoints, pch = 19)


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
    f <- function(s) {
      d <- ctr + gyromidpoint(P1-ctr, P2-ctr, s)
      Ad <- d - A
      Bd <- d - B
      sqrt(c(crossprod(Ad))) - sqrt(c(crossprod(Bd))) - (rA - rB)
    }
    ur <- uniroot(f, lower = 0.01, upper = 5)
    s <- ur$root
    hseg <- t(ctr + t(gyrosegment(P1-ctr, P2-ctr, s = s)))
    lines(hseg, col="black", lwd = 2)
  }
}
