library(Apollonius)

phi <- (1+sqrt(5)) / 2
theta <- seq(1, 30, length.out = 100)
x <- cos(theta) * phi^(theta/pi)
y <- sin(theta) * phi^(theta/pi)
sites <- cbind(x, y)
radii <- seq(1, 100, length.out = 100L)

apo <- Apollonius(sites, radii)

opar <- par(mar = c(2, 2, 1, 1))
plotApolloniusGraph(apo, circles = FALSE, color = "red")
par(opar)

fplot <- function(rmax) {
  radii <- seq(1, rmax, length.out = 100L)
  apo <- Apollonius(sites, radii, t0 = 3, tmax = 6)
  opar <- par(mar = c(2, 2, 1, 1))
  plotApolloniusGraph(apo, circles = FALSE, color = "red")
  par(opar)
}

fplot(rmax_[12])


rmax_ <- seq(1, 100, length.out = 15)
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
  gif_file = "goldenSpiral.gif",
  width = 512,
  height = 512,
  delay = 1/8
)
file.remove(pngs)
