library(Apollonius)

phi <- pi*(3 - sqrt(5))
n <- c(0, 50:150)
theta <- n * phi
r <- sqrt(n) * theta
x <- r * cos(theta)
y <- r * sin(theta)
sites <- cbind(x, y) / 1000
radii <- seq(1, 2, length.out = length(n))

apo <- Apollonius(sites, radii, tmax = 80)

opar <- par(mar = c(2, 2, 1, 1))
plotApolloniusGraph(apo, circles = FALSE, color = "random")
par(opar)
