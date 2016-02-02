#' Plotting the empirical Palm intensity
#'
#' Plots the Palm intensity from empirical data.
#'
#' @return A plot showing the empirical Palm intensity.
#' 
#' @inheritParams fit.ns
#' @param breaks The (approximate) number of points plotted.
#' @param xlim The x-axis limits for the plot.
#' @param add Logical, if \code{TRUE} then the line is added to a
#' plot.
#' @param ... Graphical parameters (e.g., to be passed to
#' \link{par}().
#' 
#' @export
empirical.palm <- function(points, lims, breaks = 50, xlim = NULL, add = FALSE, ...){
    error.dims(points, lims)
    n.dims <- ncol(points)
    dists <- pbc_distances(points = points, lims = lims)
    n.points <- nrow(points)
    if (is.null(xlim)){
        xlim <- 0.5*min(apply(lims, 1, diff))
    }
    midpoints <- seq(0, max(xlim), length.out = breaks)
    midpoints <- midpoints[-length(midpoints)]
    h <- diff(midpoints[c(1, 2)])
    midpoints[1] <- midpoints[1] + h/2
    intensities <- numeric(length(midpoints))
    for (i in 1:length(midpoints)){
        halfwidth <- ifelse(i == 1, 0.25*h, 0.5*h)
        n.interval <- sum(dists <= (midpoints[i] + halfwidth)) -
            sum(dists <= (midpoints[i] - halfwidth))
        area <- Vd(midpoints[i] + halfwidth, n.dims) -  Vd(midpoints[i] - halfwidth, n.dims)
        intensities[i] <- n.interval/(n.points*area)
    }
    if (!add){
        par(xaxs = "i")
        par(yaxs = "i")
        plot.new()
        plot.window(xlim = c(0, midpoints[length(midpoints)]), ylim = c(0, max(intensities)*1.04))
        box()
        axis(1)
        axis(2)
    }
    lines(midpoints, intensities, ...)
}

#' Plotting an estimated Palm intensity function.
#'
#' Plots a fitted Palm intensity function from an object returned by
#' \link{fit.ns}().
#'
#' @param x A fitted model from \link{fit.ns}().
#' @param plot.empirical Logical, if \code{TRUE} then the empirical
#' Palm intensity is also plotted.
#' @inheritParams empirical.palm
#' @param ... Graphical parameters (e.g., to be passed to
#' \link{par}().
#'
#' @method plot nspp
#'
#' @export
plot.nspp <- function(x, plot.empirical = FALSE, breaks = NULL, ...){
    Dc <- x$pars["Dc"]
    nu <- x$pars["nu"]
    child.disp <- x$pars["child.disp"]
    R <- x$args$R
    analytic.palm(Dc, nu, child.disp, ncol(x$args$points), c(0, R), ...)
    if (plot.empirical){
        empirical.palm(x$args$points, x$args$lims,
                       breaks = breaks, add = TRUE)
    }
}

## Plots the analytic Palm intensity.
analytic.palm <- function(Dc, nu, child.disp, n.dims, xlim = c(0, 1),dispersion,add = FALSE, ...){
    xx <- seq(xlim[1], xlim[2], length.out = 500)
    yy <- palm.intensity(xx, Dc, nu, child.disp, n.dims,dispersion=dispersion)
    if (!add){
        par(xaxs = "i")
        plot.new()
        plot.window(xlim = xlim, ylim = c(0, max(yy)))
        box()
        abline(h = 0, col = "lightgrey")
        axis(1)
        axis(2)
    }
    lines(xx, yy, ...)
}
