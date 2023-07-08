library(Apollonius)

sites <- rbind(
  c(-1, -1),
  c(-1, 1),
  c(1, 1),
  c(1, -1),
  c(0, 0)
)
angle_ <- seq(0, 2*pi, length.out = 13L)[-1L]
circle <- cbind(2 * cos(angle_), 2 * sin(angle_))
sites <- rbind(sites, circle)
radii <- c(rep(2, 5), rep(1, 12))

apo <- Apollonius(sites, radii)

opar <- par(mar = c(2, 2, 1, 1))
plotApolloniusGraph(apo, circles = FALSE, color = "random")
par(opar)
