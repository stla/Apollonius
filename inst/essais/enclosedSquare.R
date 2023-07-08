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

R_ <- c(1, 1.2, 1.4, 1.6)
i <- 4
R <- R_[i]
radii <- c(rep(R, 5), rep(1, 12))
apo <- Apollonius(sites, radii)
png(sprintf("zzpic%03d.png", i), width = 350, height = 350)
opar <- par(mar = c(2, 2, 1, 1))
plotApolloniusGraph(apo, circles = TRUE, color = "red", fill = FALSE,
                    main = paste0("radius = ", R))
par(opar)
dev.off()

