## Taking points and moving those outside the limits back into the limits.
pbc.fix <- function(points, lims){
    ## Errors for inconsistent dimensions.
    if (!is.matrix(points)){
        points <- matrix(points, ncol = 1)
    }
    if (!is.matrix(lims)){
        lims <- matrix(lims, nrow = 1)
    }
    if (ncol(points) != nrow(lims)){
        stop("The number of columns in 'points' and 'lims' must both equal the number of dimensions.")
    }
    n.points <- nrow(points)
    n.dims <- nrow(lims)
    lim.diffs <- apply(lims, 1, diff)
    for (i in 1:n.points){
        for (j in 1:n.dims){
            ## Checking if ith point's jth dimension is within the limits.
            in.lims <- FALSE
            while (!in.lims){
                ## If too low, increment by the distance between the limits,
                ## If too high, decrement by the distance between the limits.
                if (points[i, j] < lims[j, 1]){
                    points[i, j] <- points[i, j] + lim.diffs[j]
                } else if (points[i, j] > lims[j, 2]){
                    points[i, j] <- points[i, j] - lim.diffs[j]
                }
                in.lims <- points[i, j] >= lims[j, 1] & points[i, j] <= lims[j, 2]
            }
        }
    }
    points
}

## Analytic value for Dc given child.disp and nu.
analytic.Dc <- function(nu, child.disp, n.dists, n.points, R, d,dispersion){
    (n.dists/n.points - nu*Fd(R, child.disp, d,dispersion))/Vd(R, d)
}

## Analytic value for nu given Dc and child.disp.
analytic.nu <- function(Dc, child.disp, n.dists, n.points, R, d,dispersion){
    (n.dists/n.points - Dc*Vd(R, d))/Fd(R, child.disp, d,dispersion)
}

## Surface area of d-dimensional hypersphere with radius r.
Sd <- function(r, d){
    d*pi^(d/2)*r^(d - 1)/gamma(d/2 + 1)
}

## Volume of d-dimensional hypersphere with radius r.
Vd <- function(r, d){
    pi^(d/2)*r^d/gamma(d/2 + 1)
}

## PDF of between-sibling distances.

fd <-function(r, child.disp, d,dispersion=dispersion){
    if (dispersion=="gaussian"){
        2^(1 - d/2)*r^(d - 1)*exp(-r^2/(4*child.disp^2))/((child.disp*sqrt(2))^d*gamma(d/2))
    } else if (dispersion=="uniform"){
        ifelse(r>=2*child.disp,0,
    ((2*d)/(beta((d/2)+(1/2),1/2)))*(child.disp*hyperg_2F1(1/2,(1/2)-(d/2),3/2,1)-((r/2)*hyperg_2F1(1/2,(1/2)-(d/2),3/2,r^2/(4*child.disp^2))))*((r^(d-1))/(child.disp^(d+1))))
    }
}



## CDF of between-sibling distances. ## hypergeometric representation for the matern case

Fd <- function(r, child.disp, d,dispersion=dispersion){
    if (dispersion=="gaussian"){
        pgamma(r^2/(4*child.disp^2), d/2)
    }else if (dispersion=="uniform"){
        (r^d/child.disp^d) - (pbeta(r^2/(4*child.disp^2),1/2,(d/2)+(1/2))*r^d)/(child.disp^d) +( 2^d*pbeta(r^2/(4*child.disp^2),(d/2)+(1/2),(d/2)+(1/2))*beta((d/2)+(1/2),(d/2)+(1/2)))/beta(1/2,(d/2)+(1/2))
    }
}

## Note that Dc + nu/Sd(r, d)*fd(r, child.disp, d) is a correct
## formulation, but the below cancels the r^(d - 1) from both Sd(r, d)
## and fd(r, child.disp, d)
## note Charlotte has chaged this from hardcore version for Thomas process to be general,
## however if dists = 0 is unstable.. Bugger should change to hardcore version--have done
 palm.intensity <- function(r, Dc, nu, child.disp, d,dispersion=dispersion,siblings=NULL){
     if(dispersion=="gaussian"){
         Dc + nu/(pi^(d/2)*d/gamma(d/2 + 1))*
        2^(1 - d/2)*exp(-r^2/(4*child.disp^2))/((child.disp*sqrt(2))^d*gamma(d/2))
     } else if (dispersion=="uniform"){
         ifelse(r<=2*child.disp,
      Dc + ((2*nu)/(child.disp^(d+1)))*((gamma(d/2 +1))/(pi^(d/2 + 0.5)))*(child.disp*hyperg_2F1(1/2,(1/2)-(d/2),3/2,1)-((r/2)*hyperg_2F1(1/2,(1/2)-(d/2),3/2,r^2/(4*child.disp^2)))),Dc)
         
     }
 }

## Separate intensity function for known sibling information to
## optimise performance when there isn't any.
palm.intensity.siblings <- function(r, Dc, nu, child.disp, d, dispersion,siblings){
    ns.intensity <- Dc
    s.intensity <- nu/(pi^(d/2)*d/gamma(d/2 + 1))*
        2^(1 - d/2)*exp(-r^2/(4*child.disp^2))/((child.disp*sqrt(2))^d*gamma(d/2))
    siblings$ns.multipliers*ns.intensity +
        siblings$s.multipliers*s.intensity
}

## Takes matrix component of siblings and turns it into a vector that
## matches up to the distances computed by pbc_distance().
vectorise.siblings <- function(siblings){
    if (nrow(siblings$matrix) != ncol(siblings$matrix)){
        stop("Sibling matrix is not square.")
    }
    n.points <- nrow(siblings$matrix)
    n.comparisons <- n.points^2 - n.points
    vec <- numeric(n.comparisons)
    ns.multipliers <- numeric(n.comparisons)
    s.multipliers <- numeric(n.comparisons)
    k <- 1
    for (i in 1:(n.points - 1)){
        for (j in (i +1):n.points){
            vec[k] <- vec[k + 1] <- siblings$matrix[i, j]
            ## Multipliers for TRUE.
            if (!is.na(siblings$matrix[i, j]) & siblings$matrix[i, j]){
                ns.multipliers[k] <- ns.multipliers[k + 1] <- 0
                s.multipliers[k] <- s.multipliers[k + 1] <- siblings$pT
            }
            ## Multipliers for FALSE.
            if (!is.na(siblings$matrix[i, j]) & !siblings$matrix[i, j]){
                ns.multipliers[k] <- ns.multipliers[k + 1] <- siblings$pF
                s.multipliers[k] <- s.multipliers[k + 1] <- 0
            }
            ## Multipliers for NA.
            if (is.na(siblings$matrix[i, j])){
                ns.multipliers[k] <- ns.multipliers[k + 1] <- 1 - siblings$pF
                s.multipliers[k] <- s.multipliers[k + 1] <- 1 - siblings$pT
            }
            k <- k + 2
        }
    }
    list(vector = vec, ns.multipliers = ns.multipliers, s.multipliers = s.multipliers)
}

## Error function for incompatible dimensions.
error.dims <- function(points, lims){
    if (!is.matrix(points)){
        stop("Argument 'points' must be a matrix.")
    }
    if (!is.matrix(lims)){
        stop("Argument 'lims' must be a matrix.")
    }
    if (ncol(points) != nrow(lims)){
        stop("The number of columns in 'points' and 'lims' must both equal the number of dimensions.")
    }
}

## Functions below calculate partial derivatives for model parameters
## for two-dimensional processes with a Poisson distribution for the
## number of children.
dldD <- function(D, nu, child.disp, n.points, dists, R){
    sum(1/(D + exp(-dists^2/(4*child.disp^2))/(4*pi*child.disp^2))) - n.points*pi*nu*R^2
}

dldnu <- function(D, nu, child.disp, n.points, dists, R){
    length(dists)/nu - n.points*pi*D*R^2 - n.points + n.points*exp(-R^2/(4*child.disp^2))
}

dldchild.disp <- function(D, nu, child.disp, n.points, dists, R){
    sum((-n.points*nu*exp(-dists^2/(4*child.disp^2))/(2*pi*child.disp^3) +
             n.points*nu*dists^2*exp(-dists^2/(4*child.disp^2))/(8*pi*child.disp^5))/
        (n.points*D*nu + n.points*nu*exp(-dists^2/(4*child.disp^2))/(4*pi*child.disp^2))) +
        n.points*nu*(R^2)*exp(-R^2/(4*child.disp^2))/(2*child.disp^3)
}


## Charlotte's additions
##############################

## generatres uniformly distributed childern in a hypersphere

unifsphere <- function(n,d,R){
    ## n = number of points to be generated
    ## d = dimension
    ## R = radius of hypersphere
    # polar coordinates
    x <-  matrix(runif(n*d,-pi,pi),ncol=d)
    x[,d] <- x[,d]/2
    # cartesians
    sin.x <- sin(x)
    cos.x <- cos(x)
    cos.x[,d] <- 1  
    y <- sapply(1:d, function(i){
        if(i==1){
          cos.x[,1]
        } else {
          cos.x[,i]*apply(sin.x[,1:(i-1),drop=F],1,prod)
        }
    })*sqrt(runif(n,0,R^2))
    y <-  as.data.frame(
            t(apply(y,1,'+',rep(0,d)))
          )
    y
}



