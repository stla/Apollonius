library(plotrix)
library(gyro)

gyroray <- function(A, B, s, n=300, tmax = 50){
  t(vapply(seq(0, tmax, length.out = n), function(t){
    gyro:::UgyroABt(A, B, t, s)
  }, numeric(length(A))))
}
gyroray2 <- function(A, B, s, n=300, tmax = -50){
  t(vapply(seq(tmax, 0, length.out = n), function(t){
    gyro:::UgyroABt(A, B, t, s)
  }, numeric(length(A))))
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
radii <- c(1, 1.5, 1.25, 2, 1.75, 0.5, 0.4, 0.6, 0.7, 0.3)

clrs <- randomcoloR::distinctColorPalette(nrow(sites))

opar <- par(mar = c(3, 3, 1, 1))
plot(NULL, xlim = c(-10, 20), ylim = c(-3, 6), asp = 1, xlab = "x", ylab = "y")
for(i in 1L:nrow(sites)) {
  draw.circle(
    sites[i, 1L], sites[i, 2L], radius = radii[i],
    border = clrs[i], lwd = 4
  )
}

stuff <- Apollonius:::test2(sites, radii)
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
vertices <- stuff[["vertices"]]

#points(dpoints, pch = 19)


for(i in 1L:nrow(dpoints)) {
  P1 <- dpoints[i, ]
  if(is.na(P1[3])) {
    P1 <- P1[1:2]
    points(rbind(P1), pch = 19)
    Neighs_i <- Neighs[[i]]
    vert_i <- vertices[[i]]
    for(k in seq_along(Neighs_i)) {
      vs <- commonVertices[i, Neighs_i[k], ]
      A  <- vert_i[vs[1L], 1L:2L]
      rA <- vert_i[vs[1L], 3L]
      B  <- vert_i[vs[2L], 1L:2L]
      rB <- vert_i[vs[2L], 3L]
      ctr <- (A + B)/2
      P2 <- dpoints[neighbors[i, Neighs_i[k]], 1:2]
      f <- function(s) {
        d <- ctr + gyromidpoint(P1-ctr, P2-ctr, s)
        Ad <- d - A
        Bd <- d - B
        sqrt(c(crossprod(Ad))) - sqrt(c(crossprod(Bd))) - (rA - rB)
      }
      tryCatch({
        ur <- uniroot(f, lower = 0.01, upper = 5)
        s <- ur$root
        hseg <- t(ctr + t(gyrosegment(P1-ctr, P2-ctr, s = s)))
        lines(hseg, col="black", lwd = 2)
      }, error = function(e) {
        print("that should not happen")
      })
    }
  } else {
    clr <- randomcoloR::randomColor(hue = "random", luminosity = "random")
    a <- -P1[3]/P1[2]
    b <- -P1[1]/P1[2]
    abline(a, b, col = clr, lwd = 3, lty="dashed")
    Neighs_i <- Neighs[[i]]
    vert_i <- vertices[[i]]
    for(k in seq_along(Neighs_i)) {
      vs <- commonVertices[i, Neighs_i[k], ]
      A  <- vert_i[vs[1L], 1L:2L]
      rA <- vert_i[vs[1L], 3L]
      B  <- vert_i[vs[2L], 1L:2L]
      rB <- vert_i[vs[2L], 3L]
      AB <- sqrt(c(crossprod(B-A)))
      u <- (B-A) / AB
      P <- A + (rA + (AB - (rA + rB))/2) * u
      ctr <- (A + B)/2
      j <- neighbors[i, Neighs_i[k]]
      P2 <- dpoints[j, 1:2]
      f <- function(s) {
        d <- ctr + gyromidpoint(P-ctr, P2-ctr, s)
        Ad <- d - A
        Bd <- d - B
        sqrt(c(crossprod(Ad))) - sqrt(c(crossprod(Bd))) - (rA - rB)
      }
      tryCatch({
        ur <- uniroot(f, lower = 0.001, upper = 15)
        s <- ur$root
        message("i = ", i, "; j = ", j, "; s = ", s)
        face <- stuff[["vertices"]][[i]]
        weights <- face[, 3L]
        v <- face[which(weights != 0)[1L],]
        v_is_up <- P1[1]*v[1] + P1[2]*v[2] > P1[3]
        x <- ctr + gyroABt(P2-ctr, P-ctr, t = 2, s = s)
        print(x)
        P2_is_up <- P1[1]*x[1] + P1[2]*x[2] > P1[3]
        reverse <- v_is_up != P2_is_up
        if(reverse) {
          hseg <- t(ctr + t(gyroray2(P2-ctr, P-ctr, s = s)))
        } else {
          hseg <- t(ctr + t(gyroray(P2-ctr, P-ctr, s = s)))
        }
        lines(hseg, col=clr, lwd = 2)
      }, error = function(e) {
        print("that should not happen")
      })
    }
  }
}
