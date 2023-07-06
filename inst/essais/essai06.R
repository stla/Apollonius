library(plotrix)
library(gyro)

p1  <- c(0, 0)
p2  <- c(4, 1)
p3  <- c(2, 4)
p4  <- c(7, 4)
p5  <- c(8, 0)
p6  <- c(5, -2)
p7  <- c(-4, 4)
p8  <- c(-2, -1)
p9  <- c(11, 4)
p10 <- c(11, 0)

sites <- rbind(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10)
radii <- c(1, 1.5, 1.25, 2, 1.75, 0.5, 0.4, 0.6, 0.7, 0.3)


Apollonius <- function(sites, radii) {
  stuff <- Apollonius:::ApolloniusCpp(sites, radii)
  neighbors <- stuff[["neighbors"]]
  if(nrow(neighbors) == 0L) {
    stop("The Apollonius diagram is empty.")
  }
  #
  Neighs <- vector("list", nrow(neighbors))
  for(i in 1L:nrow(neighbors)) {
    Neighs_i <- 1L:3L
    neighs_i <- neighbors[i, ]
    remove <- which(is.na(neighs_i))
    for(j in seq_len(i-1L)) {
      if(j %in% neighs_i && i %in% neighbors[j, ]) {
        remove <- c(remove, which(neighs_i == j))
      }
    }
    if(length(remove) > 0L) {
      Neighs_i <- Neighs_i[-remove]
    }
    Neighs[[i]] <- Neighs_i
  }
  #
  commonVertices <-
    abind::abind(stuff[["cvertex1"]], stuff[["cvertex2"]], along = 3L)
  #
  dpoints <- stuff[["dpoints"]]
  vertices <- stuff[["vertices"]]
  #
  hsegments <- vector("list", sum(lengths(Neighs)))
  h <- 1L
  for(i in 1L:nrow(dpoints)) {
    P1 <- dpoints[i, ]
    Neighs_i <- Neighs[[i]]
    vert_i <- vertices[[i]]
    for(k in seq_along(Neighs_i)) {
      vs <- commonVertices[i, Neighs_i[k], ]
      A  <- vert_i[vs[1L], 1L:2L]
      rA <- vert_i[vs[1L], 3L]
      B  <- vert_i[vs[2L], 1L:2L]
      rB <- vert_i[vs[2L], 3L]
      ctr <- (A + B)/2
      P2 <- dpoints[neighbors[i, Neighs_i[k]], ]
      f <- function(log_s) {
        d <- ctr + gyromidpoint(P1-ctr, P2-ctr, exp(log_s))
        sqrt(c(crossprod(d-A))) - sqrt(c(crossprod(d-B))) - (rA - rB)
      }
      ur <- uniroot(f, lower = -5, upper = 2, extendInt = "yes")
      s <- exp(ur$root)
      hsegments[[h]] <- t(ctr + t(gyrosegment(P1-ctr, P2-ctr, s = s)))
      h <- h + 1L
    }
  }
  wsites <- cbind(sites, radii)
  colnames(wsites) <- c("x", "y", "weight")
  list(
    "diagram" = list("sites" = wsites, "faces" = stuff[["faces"]]),
    "graph"   = list("sites" = dpoints, "edges" = hsegments)
  )
}


#clrs <- randomcoloR::distinctColorPalette(10L)

plotApolloniusGraph <- function(apo, limits = NULL) {
  sites  <- apo[["diagram"]][["sites"]]
  nsites <- nrow(sites)
  radii  <- sites[, "weight"]
  dsites <- apo[["graph"]][["sites"]]
  edges  <- apo[["graph"]][["edges"]]
  #
  clrs <- randomcoloR::distinctColorPalette(nsites)
  #
  if(is.null(limits)) {
    x <- grDevices::extendrange(sites[, "x"])
    y <- grDevices::extendrange(sites[, "y"])
    limits <- c(min(x[1L], y[1L]), max(x[2L], y[2L]))
  }
  #
  plot(NULL, xlim = limits, ylim = limits, asp = 1, xlab = "x", ylab = "y")
  for(i in 1L:nsites) {
    draw.circle(
      sites[i, "x"], sites[i, "y"], radius = radii[i],
      border = clrs[i], col = clrs[i]
    )
  }
  points(dsites, pch = 19)
  for(i in 1L:length(edges)) {
    lines(edges[[i]], col="black", lwd = 2)
  }
  invisible()
}

sites <- rbind(
  c(0, 0),
  c(10.5, 0),
  c(20, 0),
  c(0, 9.5),
  #c(10.5, 9.5),
  c(20, 9.5),
  c(0, 16),
  c(10.5, 16),
  c(20, 16)
)
radii <- seq(5, 1, by = -0.5)[-5]

apo <- Apollonius(sites, radii/5)

# svg("x.svg", width = 8, height = 4)
opar <- par(mar = c(4, 4, 1, 1))
plotApolloniusGraph(apo)
par(opar)
# dev.off()
#
# rsvg::rsvg_png("x.svg", "inst/screenshots/agraph02.png", width = 512, height = 256)
# file.remove("x.svg")


for(i in 1L:10L) {
  apo <- Apollonius(sites, radii)
  svg("x.svg", width = 8, height = 4)
  opar <- par(mar = c(3, 3, 1, 1))
  plotApolloniusGraph(apo)
  par(opar)
  dev.off()
  rsvg::rsvg_png("x.svg", sprintf("zzpic%03d.png", i), width = 512, height = 256)
  radii <- radii * 0.95
}

library(gifski)
pngs <- Sys.glob("zzpic*.png")
gifski(
  png_files = c(pngs, rev(pngs)),
  gif_file = "agraph03.gif",
  width = 512,
  height = 256,
  delay = 1/8
)
file.remove(pngs)
