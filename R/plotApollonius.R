#' @title Plot Apollonius graph
#' @description Plot an Apollonius graph.
#'
#' @param apo an output of \code{\link{Apollonius}}
#' @param limits either \code{NULL} or a vector of length two passed to the
#'   arguments \code{xlim} and \code{ylim} of \code{\link[graphics]{plot}};
#'   if \code{NULL}, automatic limits are calculated
#' @param circles Boolean, whether to plot the original sites as circles with
#'   the given radii
#' @param colors a character string controlling the colors of the sites;
#'   either \code{"distinct"} for random distinct colors, \code{"random"} for
#'   random colors controlled by the arguments \code{hue} and
#'   \code{luminosity}, or a color name or hex code
#' @param hue,luminosity if \code{colors="random"}, these arguments are passed
#'   to \code{\link[randomcoloR]{randomColor}}
#' @param ... arguments passed to \code{\link[graphics]{plot}}, such as
#'   \code{xlab} and \code{ylab}
#'
#' @return No returned value, called for plotting.
#' @export
#'
#' @importFrom randomcoloR distinctColorPalette randomColor
#' @importFrom grDevices extendrange
#' @importFrom graphics plot points lines
#' @importFrom plotrix draw.circle
#'
#' @examples
#' library(Apollonius)
#' sites <- rbind(
#'   c(0, 0),
#'   c(4, 1),
#'   c(2, 4),
#'   c(7, 4),
#'   c(8, 0),
#'   c(5, -2),
#'   c(-4, 4),
#'   c(-2, -1),
#'   c(11, 4),
#'   c(11, 0)
#' )
#' radii <- c(1, 1.5, 1.25, 2, 1.75, 0.5, 0.4, 0.6, 0.7, 0.3)
#' apo <- Apollonius(sites, radii)
#' opar <- par(mar = c(3, 3, 1, 1))
#' plotApolloniusGraph(apo, colors = "random", xlab = NA, ylab = NA)
#' par(opar)
plotApolloniusGraph <- function(
    apo, limits = NULL, circles = TRUE,
    colors = "distinct", hue = "random", luminosity = "dark", ...
) {
  sites  <- apo[["diagram"]][["sites"]]
  nsites <- nrow(sites)
  radii  <- sites[, "weight"]
  dsites <- apo[["graph"]][["sites"]]
  edges  <- apo[["graph"]][["edges"]]
  #
  if(colors == "distinct") {
    clrs <- distinctColorPalette(nsites)
  } else if(colors == "random") {
    clrs <- randomColor(nsites, hue = hue, luminosity = luminosity)
  } else {
    clrs <- rep(colors, nsites)
  }
  #
  if(is.null(limits)) {
    x <- extendrange(sites[, "x"])
    y <- extendrange(sites[, "y"])
    limits <- c(min(x[1L], y[1L]), max(x[2L], y[2L]))
  }
  #
  plot(NULL, xlim = limits, ylim = limits, asp = 1, ...)
  if(circles) {
    for(i in 1L:nsites) {
      draw.circle(
        sites[i, "x"], sites[i, "y"], radius = radii[i],
        border = clrs[i], col = clrs[i]
      )
    }
  } else {
    for(i in 1L:nsites) {
      points(
        sites[i, "x"], sites[i, "y"], pch = 19, col = clrs[i]
      )
    }
  }
  points(dsites, pch = 19)
  for(i in 1L:length(edges)) {
    lines(edges[[i]], col="black", lwd = 2)
  }
  invisible()
}
