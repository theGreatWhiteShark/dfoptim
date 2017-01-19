##
##  h o o k e j e e v e s . R  Hooke-Jeeves Minimization Algorithm
##

## I will not use this one because it makes some very bold steps in the beginning considering the
## shape parameter. If not tuned accordingly it is not usable in the GEV fitting. And I do not see
## the point right now in tweaking yet another fitting procedure.
hjk <- function( par = climex::likelihood.initials( x ), fn = climex::likelihood, x = x,
                control = list(), ... ) {
    if (!is.numeric(par))
        stop("Argument 'par' must be a numeric vector.", call. = FALSE)
    n <- length(par)
    if (n == 1)
        stop("For univariate functions use some different method.", call. = FALSE)
   
    #-- Control list handling ----------
    cntrl <- list(tol      = 1.e-06,
                  maxfeval = Inf,       # set to Inf if no limit wanted
                  maximize = FALSE,     # set to TRUE  for maximization
                  target   = Inf,       # set to Inf for no restriction
                  info     = FALSE)     # for printing interim information
    nmsCo <- match.arg(names(control), choices = names(cntrl), several.ok = TRUE)
    if (!is.null(names(control))) cntrl[nmsCo] <- control

    tol      <- cntrl$tol;      
    maxfeval <- cntrl$maxfeval
    maximize <- cntrl$maximize
    target   <- cntrl$target
    info     <- cntrl$info

	scale <- if (maximize) -1 else 1
    fun <- match.fun(fn)
    ## MOD
    f <- function(x) scale * fun( par, x )
    ## /MOD
    
    #-- Setting steps and stepsize -----
    nsteps <- floor(log2(1/tol))        # number of steps
    steps  <- 2^c(-(0:(nsteps-1)))      # decreasing step size
    dir <- diag(1, n, n)                # orthogonal directions

    x.up <- par                            # start point
    fx <- fbest <- f(x.up)                 # smallest value so far
    fcount <- 1                         # counts number of function calls

    if (info) cat("step\tnofc\tfmin\txpar\n")   # info header

    #-- Start the main loop ------------
    ns <- 0
    while (ns < nsteps && fcount < maxfeval && abs(fx) < target) {
        ns <- ns + 1
        hjs    <- .hjsearch(x.up, f, steps[ns], dir, fcount, maxfeval, target)
        browser()
        x.up      <- hjs$x
        fx     <- hjs$fx
        sf     <- hjs$sf
        fcount <- fcount + hjs$finc

        if (info)
            cat(ns, "\t",  fcount, "\t", fx/scale, "\t", x.up[1], "...\n")
    }

    if (fcount > maxfeval) {
        warning("Function evaluation limit exceeded -- may not converge.")
        conv <- 1
    } else if (abs(fx) > target) {
        warning("Function exceeds min/max value -- may not converge.")
        conv <- 1
    } else {
        conv <- 0
    }

    fx <- fx / scale                    # undo scaling
    return(list(par = x.up, value = fx,
                convergence = conv, feval = fcount, niter = ns))
}

##  Search with a single scale -----------------------------
.hjsearch <- function(xb, f, h, dir, fcount, maxfeval, target) {
    x  <- xb
    xc <- x
    sf <- 0
    finc <- 0
    hje  <- .hjexplore(xb, xc, f, h, dir)
    x    <- hje$x
    fx   <- hje$fx
    sf   <- hje$sf
    finc <- finc + hje$numf

    # Pattern move
    while (sf == 1) {
        d  <- x-xb
        xb <- x
        xc <- x+d
        fb <- fx
        hje  <- .hjexplore(xb, xc, f, h, dir, fb)
        x    <- hje$x
        fx   <- hje$fx
        sf   <- hje$sf
        finc <- finc + hje$numf

        if (sf == 0) {  # pattern move failed
           hje  <- .hjexplore(xb, xb, f, h, dir, fb)
           x    <- hje$x
           fx   <- hje$fx
           sf   <- hje$sf
           finc <- finc + hje$numf
        }
        if (fcount + finc > maxfeval || abs(fx) > target) break
    }

    return(list(x = x, fx = fx, sf = sf, finc = finc))
}

##  Exploratory move ---------------------------------------
.hjexplore <- function(xb, xc, f, h, dir, fbold) {
    n <- length(xb)
    x <- xb

    if (missing(fbold)) {
        fb <- f(x)
        numf <- 1
    } else {
        fb <- fbold
        numf <- 0
    }
    
    ## MOD
    ## this is a quite dirty solution but lets keep it for a moment
    if ( is.nan( fb ) )
        fb <- Inf
    ## /MOD
    fx <- fb
    xt <- xc
    sf <- 0                             # do we find a better point ?
    dirh <- h * dir
    fbold <- fx
    for (k in sample.int(n, n)) {       # resample orthogonal directions
        p <- xt + dirh[, k]
        ## MOD
        ## this is a quite dirty solution but lets keep it for a moment
        ft <- f(p)
        if ( is.nan( ft ) )
            ft <- Inf
        ## /MOD
        numf <- numf + 1
        
        if (ft >= fb) {
            p <- xt - dirh[, k]
            ft <- f(p)
            ## MOD
            ## this is a quite dirty solution but lets keep it for a moment
            if ( is.nan( ft ) )
                ft <- Inf
            ## /MOD
            numf <- numf + 1
        }
        if (ft < fb) {
            sf <- 1
            xt <- p
            fb <- ft
        }
    }
    if (sf == 1) {
        x <- xt
        fx <- fb
    }
    return(list(x = x, fx = fx, sf = sf, numf = numf))
}
