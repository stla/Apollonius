library(plotrix)
library(gyro)

p1 <- c(0, 0)
p2 <- c(4, 1)
p3 <- c(2, 4)
p4 <- c(7, 4)
p5 <- c(8, 0)
p6 <- c(5, -2)

sites <- rbind(p1, p2, p3, p4, p5, p6)
radii <- c(1, 1.5, 1, 2, 1, 0.5)

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
Neighs <- c()
x <- na.omit(neighbors[1L, ])[1L]
Neighs[1L] <- which(neighbors[1L, ] == x[1L])
for(i in 2L:nrow(neighbors)) {
  x <-
    setdiff(na.omit(neighbors[i, ]), neighbors[cbind(1L:(i-1L), Neighs[1L:(i-1L)])])
  Neighs[i] <- which(neighbors[i, ] == x[1L])
}

commonVertices <-
  abind::abind(stuff[["cvertex1"]], stuff[["cvertex2"]], along = 3L)
commonVertices[1, 3, ] # common vertices face1 with its neighbor 3

dpoints <- stuff[["dpoints"]]

points(dpoints, pch = 19)


# check equidistances:
P1 <- dpoints[1L, ]
vs <- commonVertices[1L, Neighs[1L], ]
A <- vertices[[1L]][vs[1L], 1L:2L]
rA <- vertices[[1L]][vs[1L], 3L]
B <- vertices[[1L]][vs[2L], 1L:2L]
rB <- vertices[[1L]][vs[2L], 3L]
sqrt(c(crossprod(P1 - A))) - rA
sqrt(c(crossprod(P1 - B))) - rB

ctr <- (A + B)/2
P2 <- dpoints[neighbors[1L, Neighs[1L]], ]
f <- function(s) {
  d <- ctr + gyromidpoint(P1-ctr, P2-ctr, s)
  Ad <- d - A
  Bd <- d - B
  sqrt(c(crossprod(Ad))) - sqrt(c(crossprod(Bd))) - (rA - rB)
}
ur <- uniroot(f, lower = 0.01, upper = 5)
s <- ur$root
#points(rbind(ctr), col = "blue", pch = 19)
hseg <- t(ctr + t(gyrosegment(P1-ctr, P2-ctr, s = s)))
lines(hseg, col="black", lwd = 2)


for(i in 2L:nrow(dpoints)) {
  P1 <- dpoints[i, ]
  vs <- commonVertices[i, Neighs[i], ]
  A <- vertices[[i]][vs[1L], 1L:2L]
  rA <- vertices[[i]][vs[1L], 3L]
  B <- vertices[[i]][vs[2L], 1L:2L]
  rB <- vertices[[i]][vs[2L], 3L]
  ctr <- (A + B)/2
  P2 <- dpoints[neighbors[i, Neighs[i]], ]
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
