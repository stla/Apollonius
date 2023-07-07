library(Apollonius)

theta <- seq(0, 25, length.out = 30L)
x <- sqrt(theta) * cos(theta)
y <- sqrt(theta) * sin(theta)
sites <- cbind(x, y)
radii <- seq(1, 4, length.out = 30L)

apo <- Apollonius(sites, radii, t0 = 15)

opar <- par(mar = c(2, 2, 1, 1))
plotApolloniusGraph(apo, circles = FALSE, color = "random")
par(opar)

fplot <- function(rmax) {
  radii <- seq(1, rmax, length.out = 30L)
  apo <- Apollonius(sites, radii, t0 = 20, tmax = 80)
  opar <- par(mar = c(2, 2, 3, 1))
  plotApolloniusGraph(apo, circles = FALSE, color = "red")
  par(opar)
}


rmax_ <- seq(1, 4, length.out = 20)
for(i in seq_along(rmax_)) {
  svg("x.svg", width = 8, height = 8)
  fplot(rmax_[i])
  dev.off()
  rsvg::rsvg_png("x.svg", sprintf("zzpic%03d.png", i), width = 512, height = 512)
}

library(gifski)
pngs <- Sys.glob("zzpic*.png")
gifski(
  png_files = c(pngs, rev(pngs)),
  gif_file = "FermatSpiral.gif",
  width = 512,
  height = 512,
  delay = 1/8
)
file.remove(pngs)
