library(Apollonius)

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
apo <- Apollonius(sites, radii, tmax = 30)
opar <- par(mar = c(4, 4, 1, 1))
plotApolloniusGraph(apo, xlab = "x", ylab = "y")
par(opar)
