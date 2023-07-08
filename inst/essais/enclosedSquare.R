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


# anim
R_ <- seq(0.5, 1.6, length.out = 40)
for(i in seq_along(R_)) {
  radii <- c(rep(R_[i], 5), rep(1, 12))
  apo <- Apollonius(sites, radii, tmax = 80)
  svg("x.svg", width = 8, height = 8)
  opar <- par(mar = c(2, 2, 1, 1))
  plotApolloniusGraph(
    apo, circles = TRUE, colors = "red", fill = FALSE, limits = c(-3, 3)
  )
  par(opar)
  dev.off()
  rsvg::rsvg_png(
    "x.svg", sprintf("zzpic%03d.png", i), width = 512, height = 512
  )
}

library(gifski)
pngs <- Sys.glob("zzpic*.png")
gifski(
  png_files = c(pngs, rev(pngs)),
  gif_file = "enclosedSquare2.gif",
  width = 512,
  height = 512,
  delay = 1/9
)
file.remove(pngs)
