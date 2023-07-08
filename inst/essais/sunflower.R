library(Apollonius)

phi <- pi*(3 - sqrt(5))
n <- c(0, 50:100)
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


fplot <- function(rmax) {
  radii <- seq(1, rmax, length.out = length(n))
  apo <- Apollonius(sites, radii, tmax = 80)
  opar <- par(mar = c(2, 2, 3, 1))
  plotApolloniusGraph(apo, circles = FALSE, color = "darkred")
  par(opar)
}


rmax_ <- seq(1, 2, length.out = 20)
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
  gif_file = "sunflower.gif",
  width = 512,
  height = 512,
  delay = 1/8
)
file.remove(pngs)
