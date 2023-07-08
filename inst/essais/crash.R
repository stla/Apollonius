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

library(plotrix)
plot(sites, type = "p", pch = 19, asp = 1, ylim = c(-3, 3))
for(i in 1L:nrow(sites)) {
  draw.circle(sites[i,1L], sites[i,2L], radius = radii[i], border = "red", lwd = 2)
}


apo <- Apollonius(sites, radii)

opar <- par(mar = c(2, 2, 1, 1))
plotApolloniusGraph(apo, circles = FALSE, color = "random")
par(opar)
