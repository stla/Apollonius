library(plotrix)
library(gyro)

p1 <- c(0, 0)
p2 <- c(4, 1)
p3 <- c(2, 4)
p4 <- c(7, 4)
p5 <- c(8, 0)

plot(NULL, xlim = c(-3, 6), ylim = c(-3, 6), asp = 1)
draw.circle(p1[1L], p1[2L], radius = 1, border = "red")
draw.circle(p2[1L], p2[2L], radius = 1.5, border = "green")
draw.circle(p3[1L], p3[2L], radius = 1, border = "blue")
draw.circle(p4[1L], p4[2L], radius = 2, border = "gray")
draw.circle(p5[1L], p5[2L], radius = 1, border = "yellow")

Apollonius:::test()
P1 <- c(1.35416, 1.82292)
P2 <- c(3.99138, 3.5431)
P3 <- c(6.4288, 1.12421)
P4 <- c(4, -15.75)
points(rbind(P1), pch = 19)
points(rbind(P2), pch = 19)
points(rbind(P3), pch = 19)
points(rbind(P4), pch = 19)

# check equidistances:
sqrt(c(crossprod(P1 - p3))) - 1
sqrt(c(crossprod(P1 - p2))) - 1.5
ctr <- (p2 + p3)/2
f <- function(s) {
  d <- ctr + gyromidpoint(P1-ctr, P2-ctr, s)
  p3d <- d - p3
  p2d <- d - p2
  sqrt(c(crossprod(p3d))) - sqrt(c(crossprod(p2d))) + 0.5
}
ur <- uniroot(f, lower = 0.1, upper = 5)
s <- ur$root
#points(rbind(ctr), col = "blue", pch = 19)
hseg <- t(ctr + t(gyrosegment(P1-ctr, P2-ctr, s = s)))
lines(hseg, col="navy")

ctr <- (p2 + p4)/2
f <- function(s) {
  d <- ctr + gyromidpoint(P2-ctr, P3-ctr, s)
  p4d <- d - p4
  p2d <- d - p2
  sqrt(c(crossprod(p4d))) - sqrt(c(crossprod(p2d))) - 0.5
}
ur <- uniroot(f, lower = 0.1, upper = 5)
( s <- ur$root )
#points(rbind(ctr), col = "blue", pch = 19)
hseg <- t(ctr + t(gyrosegment(P2-ctr, P3-ctr, s = s)))
lines(hseg, col="navy")

ctr <- (p2 + p1)/2
f <- function(s) {
  d <- ctr + gyromidpoint(P4-ctr, P1-ctr, s)
  p1d <- d - p1
  p2d <- d - p2
  sqrt(c(crossprod(p1d))) - sqrt(c(crossprod(p2d))) + 0.5
}
ur <- uniroot(f, lower = 0.1, upper = 5)
( s <- ur$root )
#points(rbind(ctr), col = "blue", pch = 19)
hseg <- t(ctr + t(gyrosegment(P1-ctr, P4-ctr, s = s)))
lines(hseg, col="navy")

ctr <- (p2 + p5)/2
f <- function(s) {
  d <- ctr + gyromidpoint(P4-ctr, P3-ctr, s)
  p5d <- d - p5
  p2d <- d - p2
  sqrt(c(crossprod(p5d))) - sqrt(c(crossprod(p2d))) + 0.5
}
ur <- uniroot(f, lower = 0.1, upper = 5)
( s <- ur$root )
#points(rbind(ctr), col = "blue", pch = 19)
hseg <- t(ctr + t(gyrosegment(P3-ctr, P4-ctr, s = s)))
lines(hseg, col="navy")






################################################################################
d <- hseg[90, ]
pd <- d - p3
sqrt(c(crossprod(pd))) - 1
pd <- d - p2
sqrt(c(crossprod(pd))) - 1.5


f0 <- ctr
f1 <-

P <- c(1.67381, 1.6631)
a <- (sqrt(c(crossprod(P-p2))) - sqrt(c(crossprod(P-p3)))) / 2
c <- sqrt(c(crossprod(p2 - p3)))/2
b2 <- c^2 - a^2
sqrt(1 + b2/a^2)

vertex <- ctr + a * (p3 - p2)/(2*c) -> d
axis <- p3 - p2
tg <- c(-axis[2], axis[1])

f0 <- ctr
f1 <- vertex - f0
f2 <- tg * 0.4968725

t <- 1
d <- f0 + f1*cosh(t) + f2*sinh(t)

pd <- d - p3
sqrt(c(crossprod(pd))) - 1
pd <- d - p2
sqrt(c(crossprod(pd))) - 1.5

f <- function(h) {
  f2 <- h * tg
  t <- 1
  d <- f0 + f1*cosh(t) + f2*sinh(t)
  p3d <- d - p3
  p2d <- d - p2
  sqrt(c(crossprod(p3d))) - sqrt(c(crossprod(p2d))) + 0.5
}

uniroot(f, lower = 0.1, upper = 5)
h <- 0.4968725

f <- function(s) {
  d <- ctr + gyromidpoint(P1-ctr, P2-ctr, s)
  p3d <- d - p3
  p2d <- d - p2
  sqrt(c(crossprod(p3d))) - sqrt(c(crossprod(p2d))) + 0.5
}

uniroot(f, lower = 0.1, upper = 5)

s <- 2.22207
