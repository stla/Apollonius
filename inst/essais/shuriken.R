library(Apollonius)

nsides <- 6L
angles <- seq(0, 2*pi, length.out = nsides+1L)[-1L]
outer_points <- cbind(cos(angles), sin(angles))
inner_points <- outer_points / 4
nsides <- 12L
angles <- seq(0, 2*pi, length.out = nsides+1L)[-1L]
middle_points <- cbind(cos(angles), sin(angles)) / 2
points <- rbind(outer_points, inner_points, middle_points)
angles <- angles + pi/24
middle_points <- cbind(cos(angles), sin(angles)) / 3
points <- rbind(points, middle_points)
middle_points <- cbind(cos(angles), sin(angles)) / 1.5
points <- rbind(points, middle_points)*1

R <- 0.005
r <- 0.001
radii <- c(rep(R, 12), rep(r, 36))
radii <- runif(48, 0.1, 0.3)
apo <- Apollonius(points, radii)

plotApolloniusGraph(apo)
