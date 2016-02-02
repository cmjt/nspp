#' Bootstrapping a Neymann-Scott point process model
#'
#' Carries out a bootstrap for Neymann-Scott point process models
#' fitted by \link{fit.ns}().
#'
#' @return The first argument, with added information from the
#' bootstrap procedure.
#'
#' @param fit A fitted object from \link{fit.ns}().
#' @param N The number of bootstrap resamples.
#' @param prog Logical, if \code{TRUE}, a progress bar is printed to
#' the console.
#' @inheritParams sim.ns
#' 
#' @export
boot.ns <- function(fit, rchild, N, prog = TRUE){
    ## Extracting information.
    args <- fit$args    
    pars <- fit$pars
    n.pars <- length(pars)
    lims <- args$lims
    ## Error for fits with known (non-)siblings.
    if (!is.null(args$siblings)){
        stop("Bootstrapping not implemented for models with known (non-)siblings.")
    }
    boots <- matrix(0, nrow = N, ncol = n.pars)
    ## Setting up progress bar.
    if (prog){
        pb <- txtProgressBar(min = 0, max = N, style = 3)
    }
    for (i in 1:N){
        args$points <- tryCatch(sim.ns(pars = pars[c("D", "child.disp", "child.par")], lims = lims,
                              rchild = rchild),error=function(e) e)
        if(inherits(args$points,"error")) next
        args$child.disp.sv <- pars["child.disp"]
        args$child.dist$sv <- pars["child.par"]
        args$trace <- FALSE
        fit.boot <- tryCatch(do.call("fit.ns", args),error=function(e) e)
        if(inherits(fit.boot,"error")) next
        boots[i, ] <- fit.boot$pars
        ## Updating progress bar.
        if (prog){
            setTxtProgressBar(pb, i)
        }
    }
    if (prog){
        close(pb)
    }
    colnames(boots) <- names(pars)
    fit$boots <- boots
    fit$se <- apply(boots, 2, sd)
    class(fit) <- c("boot.nspp", class(fit))
    fit
}
