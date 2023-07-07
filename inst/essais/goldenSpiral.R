library(Apollonius)

phi <- (1+sqrt(5)) / 2
theta <- seq(1, 30, length.out = 100)
x <- cos(theta) * phi^(theta/pi)
y <- sin(theta) * phi^(theta/pi)
sites <- cbind(x, y)
radii <- seq(1, 5, length.out = 100L)

apo <- Apollonius(sites, radii)

plotApolloniusGraph(apo, circles = FALSE)
