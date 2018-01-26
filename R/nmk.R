##' @title dfoptim::nmk which outputs the parameter vector of every
##'   optimization step
##'
##' @details See \code{\link{dfoptim::nmk}}. This version incorporates
##'   a modified step for avoiding the forbidden region. Its activated
##'   by setting MODIFIED to TRUE. Since its always about the
##'   likelihood function its x argument also will be provided
##'   explicitly (much nicer for debugging; modified regions indicated
##'   by MOD)
##'
##' @param par Initial parameter set. Default = call of the
##'   \code{\link{likelihood.initials}} function with the time.series
##'   as argument
##' @param fn Function which is about to optimize. Default =
##'   \code{\link{likelihood}} 
##' @param x Time series of class 'xts' or 'numeric'.
##' @param model Determines if to use the GEV or GP
##'   distribution. Default = "gev". 
##' @param MODIFIED Flag specifying if the modified or the original
##'   algorithm should be used. Default = TRUE.
##' @param WARNINGS Flag if warnings should be displayed or
##'   not. Default = FALSE. 
##' @param control List of options.
##' @param ... Additional input for the function 'fn'. In case of the
##'   likelihood function the time series 'x' must be provided
##'
##' @return List of elements see original function. In addition two
##'   new elements called \emph{x.updates} and \emph{simplex.final}
##'   are provided. The first one contains the parameters at each
##'   optimization step. The latter one contains the a n x
##'   (n+1)-dimensional matrix containing all the coordinates of the
##'   vertices of the final simplex.
##'
##'   Of class c( "list", "climex.gev.fit" )
##' @author Philipp Mueller
##' @export
##' @examples
##' nmk( x = block( anomalies( temp.potsdam ) ) )
nmk <- function ( par = climex::likelihood.initials( x, model = model ),
                 fn = climex::likelihood, x, model = c( "gev", "gpd" ),
                 MODIFIED = TRUE, WARNINGS = FALSE, control = list(), ...) 
{
  ## Initialization
  if ( missing( model ) ){
    model <- "gev"
  }
  model <- match.arg( model )
  ctrl <- list( tol = 1e-06, maxfeval = min( 5000, max( 1500, 20 * length( par )^2 ) ),
               regsimp = TRUE, maximize = FALSE, restarts.max = 3, trace = FALSE )
  namc <- match.arg( names( control ), choices = names( ctrl ), several.ok = TRUE )
  if ( !all( namc %in% names( ctrl ) ) ) 
    stop( "unknown names in control: ", namc[ !( namc %in% names( ctrl ) ) ] )
  if (!is.null( names( control ) ) ) 
    ctrl[ namc ] <- control
  ftol <- ctrl$tol
  maxfeval <- ctrl$maxfeval
  regsimp <- ctrl$regsimp
  restarts.max <- ctrl$restarts.max
  maximize <- ctrl$maximize
  trace <- ctrl$trace
  ## MOD
  if ( maximize ) {
    fnm <- function( par, x ) -fn( par, x, model = model )
  } else
    fnm <- function( par, x ) fn( par, x, model = model )
  ## /MOD

  ## Starting values
  x0 <- as.numeric( par ) # initial parameters
  n <- length( par )
  if ( n == 1 ) 
    stop( call. = FALSE, "Use `optimize' for univariate optimization" )
  if ( n > 30 && WARNINGS ) 
    warning( "Nelder-Mead should not be used for high-dimensional optimization" )
  V <- cbind( rep( 0, n ), diag( n ) ) # matrix containing the n+1 vertices of the simplex
  f <- rep( 0, n + 1 ) # contains function evaluations
  ## MOD
  f[ 1 ] <- fnm( x0, x )
  if ( is.nan( f[ 1 ] ) && WARNINGS )
    warning( "The supplied initial parameter set to nmk.modified can not be evaluated!" )
  nf <- 0
  restarts <- 0
  ## /MOD
  V[ , 1 ] <- x0
  scale <- max( 1, sqrt( sum( x0^2 ) ) )
  ## to obtain the other three initial points for the simplex the provided initial coordinate is shifted by a constant value alpha[ 2 ] and then moved by alpha[ 1 ] in one of the three dimensions (location, scale and shape ) respectively. But this happens to be completely independent of the forbidden region. So maybe one can think of something better for the GEV case.
  if ( regsimp ) {
    ## So here again. the alphas have to be chosen in such a way that there is no evaluation in the forbidden region
    alpha <- scale/( n * sqrt( 2 ) ) * c( sqrt( n + 1 ) + n - 1, sqrt( n + 1 ) - 1 )
    if ( MODIFIED ){
      if ( is.nan( fnm( x0 + alpha[ 2 ], x ) ) ){
        ## The shifted point is in the forbidden region. Since I see no motivation in the literature to have any specific initialization of the simplex I will just figure out another point.
        alpha.2.vec <- seq( -alpha[ 2 ]* .5, alpha[ 2 ]* .5, , 200 )
        alpha.2.v <- as.numeric( lapply( alpha.2.vec, function( y ) fnm( x0 + y, x ) ) )
        if ( all( is.nan( alpha.2.v ) ) ){
          warning( "In the alpha[2] range in nmk.modified no point can be evaluated. Initialization of the simplex failed." )
          return ( list( par = c( NaN, NaN, NaN ), value = NaN, feval = 0,
                        restarts = 0, convergence = conv, message = message,
                        x.updates = NULL ) )
        }
        alpha.2.v[ is.nan( alpha.2.v ) ] <- Inf
        ## We don't want to be at the very bottom of the parabola in the likelihood with respect to alpha.2.vec but on the lower regions.
        alpha.2.min <- alpha.2.vec[ which.min( abs( alpha.2.v - min( alpha.2.v )* 1.25 ) ) ]
        V[ , - 1 ] <- ( x0 + alpha.2.min )
      } else {
        alpha.2.min <- alpha[ 2 ]
        V[ , -1 ] <- ( x0 + alpha.2.min )
      }
      ## spanning the orthogonal directions.
      diag( V[ , - 1 ] ) <- x0[ 1 : n ] + alpha[ 1 ]
      for ( jj in 2 : ncol( V ) )
        f[ jj ] <- fnm( V[ ,jj ], x )
      if ( any( is.nan( f ) ) ){
        ## At least one point in the initial simplex is in the forbidden region.
        ## In general I don't see any problem in the spanned points having different distances alpha[ 1 ] from the shifted point x0 + alpha[ 2 ]. Sometimes only very close ones can be evaluated. Therefore the distance will start at 0 and the most far away value with an acceptable likelihood will be chosen.
        alpha.1.vec <- seq( 0, alpha[ 1 ]* .5, , 120 )
        for ( cc in 1 : ( ncol( V ) - 1 ) ){
          alpha.1.v <- as.numeric( lapply( alpha.1.vec, function( y ){
            x0[ cc ] <- x0[ cc ] + y
            fnm( x0 + alpha.2.min, x ) } ) )
          alpha.1.v[ is.nan( alpha.1.v ) ] <- Inf
          ## Get the most far away value which is lowest
          alpha.1.min <- length( alpha.1.v ) + 1 -
            which.min( rev( abs( alpha.1.v - min( alpha.1.v )* 1.25 ) ) )
          if ( is.nan( alpha.1.v[ alpha.1.min ] ) )
            alpha.1.min <- which.min( alpha.1.v )
          V[ cc, cc + 1 ] <- x0[ cc ] + alpha.2.min + alpha.1.vec[ alpha.1.min ]
        }
        for ( jj in 2 : ncol( V ) )
          f[ jj ] <- fnm( V[ ,jj ], x )
      }
    } else {
      V[ , -1 ] <- ( x0 + alpha[ 2 ] )
      diag( V[ , -1 ] ) <- x0[ 1 : n ] + alpha[ 1 ]
      ## MOD
      for ( jj in 2 : ncol( V ) )
        f[ jj ] <- fnm( V[ , jj ], x )
      ## /MOD
    }
  } else {
    V[ , -1 ] <- x0 + scale * V[ , -1 ]
    ## MOD
    for ( jj in 2 : ncol( V ) )
      f[ jj ] <- fnm( V[ , jj ], x )
    ## /MOD
  }
  f[ is.nan( f ) ] <- Inf
  nf <- n + 1 # number of function evaluation
  ord <- order( f ) # order the sample to replace always the last entry
  f <- f[ ord ]
  V <- V[ , ord ]

  ## Setting up the Nelder-Mead optimization
  rho <- 1 # weight for the reflection step
  gamma <- 0.5 # rho*gamma is the weight for the outer contraction step
  chi <- 2 # rho*chi is the weight for the expansion step
  sigma <- 0.5
  conv <- 1
  oshrink <- 0
  restarts <- 0
  orth <- 0
  dist <- f[ n + 1 ] - f[ 1 ]
  v <- V[ , -1 ] - V[ , 1 ] # matrix containing the n simplex directions
  delf <- f[ -1 ] - f[ 1 ] 
  diam <- sqrt( colSums( v^2 ) )
  sgrad <- c( crossprod( t( v ), delf ) ) # simplex gradient
  alpha <- 1e-04 * max( diam )/ sqrt( sum( sgrad^2 ) )
  simplex.size <- sum( abs( v ) )/ max( 1, sum( abs( V[ ,1 ] ) ) )
  itc <- 0
  conv <- 0
  ## MOD (but even beforehand)
  if ( model == "gev" ){
    x.new.vector <- data.frame( location = x0[ 1 ], scale = x0[ 2 ],
                               shape = x0[ 3 ], step = 1 )
  } else {
    x.new.vector <- data.frame( scale = x0[ 1 ], shape = x[ 2 ], step = 1 )
  }
  xnew <- NULL # to check later on if it was already modified
  ## /MOD
  message <- "Run into error"
  while ( nf < maxfeval & restarts < restarts.max & dist > ftol & simplex.size > 1e-06 ){
    ## Staring the optimization 
    fbc <- mean( f )
    happy <- 0
    itc <- itc + 1
    xbar <- rowMeans( V[ , 1 : n ] ) # centroid of the convex simplex hull
    
    ## Reflection step
    xr <- ( 1 + rho )* xbar - rho* V[ , n + 1 ]
    ## MOD
    fr <- fnm( xr, x )
    ## There is a problem right here. If this evaluation returns NaN the algorithm doesn't know how to handle it. All the following if statements are false and if fr is NaN for the initial value 'xnew' will not be defined and the algorithm throws an error. Therefore a bunch of different rhos will be evaluated
    if ( is.nan( fr ) && MODIFIED ){
      rho.sequence <- seq( .5, 1.5, .1 )
      x.sequence <- lapply( rho.sequence, function( y ) {
        ( 1 + y )* xbar - y* V[ , n + 1 ] } )
      rho.results <- lapply( x.sequence, function( y ) fnm( y, x ) )
      if ( all( is.nan( as.numeric( rho.results ) ) ) ){
        if ( WARNINGS ){
          browser()
          warning( "no allowed points in the reflection step in nmk.modified" )
        }
        ## Well maybe the other methods are of more luck
        xe <- ( 1 + rho* chi )* xbar - rho* chi* V[ , n + 1 ]
        xco <- ( 1 + rho* gamma )* xbar - rho* gamma* V [ , n + 1 ]
        xci <- ( 1 - gamma )* xbar + gamma* V[ , n + 1 ]
        other.results <- as.numeric( lapply( list( xe, xco, xci ),
                                            function( y ) fnm( y, x ) ) )
        if ( all( is.nan( other.results ) ) && is.null( xnew ) ){
          if ( WARNINGS )
            warning( "All other steps couldn't produce a valid step either. The provided starting point!" )
          ## So especially for the parameter region of the shape where the likelihood isn't even defined it makes no sense to force the algorithm to work. Sometimes one has to let things go.
          fr <- xr <- NaN
        } else {
          other.results[ is.nan( other.results ) ] <- Inf
          fr <- other.results[ which.min( other.results ) ] # This is kinda stupid but shouldn't break anything
          xr <- rbind( xe, xco, xci )[ which.min( other.results ), ]
        }
      } else {
        rho.results <- as.numeric( rho.results )
        rho.results[ is.nan( rho.results ) ] <- Inf
        fr <- rho.results[ which.min( rho.results ) ]
        xr <- x.sequence[[ which.min( rho.results ) ]]
      } }
    ## /MOD
    nf <- nf + 1
    if ( is.nan( fr ) && is.null( xnew ) ){
      ## Is only used when xnew is not defined yet
      return ( list( par = c( NaN, NaN, NaN ), value = NaN, feval = nf,
                    restarts = restarts, convergence = conv, message = message,
                    x.updates = x.new.vector ) )
    } else if ( is.nan( fr ) )
      fr <- Inf
    if ( fr >= f[ 1 ] & fr < f[ n ] ){
      ## Reflection is successful and the update is accepted
      happy <- 1
      xnew <- xr
      fnew <- fr
      ## /Reflection step
    } else if ( fr < f[ 1 ] ){
      
      ## Expansion step
      xe <- ( 1 + rho* chi )* xbar - rho* chi* V[ , n + 1 ]
      ## MOD
      fe <- fnm( xe, x )
      ## /MOD
      if ( is.nan( fe ) ) 
        fe <- Inf
      nf <- nf + 1
      if ( fe < fr ){
        ## Expansion step was successful
        xnew <- xe
        fnew <- fe
        happy <- 1
        ## /Expansion step
      } else {
        xnew <- xr
        fnew <- fr
        happy <- 1
      }
    } else if ( fr >= f[ n ] & fr < f[ n + 1 ] ) {
      
      ## Outer contraction step
      xc <- ( 1 + rho* gamma )* xbar - rho* gamma* V [ , n + 1 ]
      ## MOD
      fc <- fnm( xc, x )
      ## /MOD
      if ( is.nan( fc ) ) 
        fc <- Inf
      nf <- nf + 1
      if ( fc <= fr ) {
        xnew <- xc
        fnew <- fc
        happy <- 1
        ## /Outer contraction step
      }
    } else if ( fr >= f[ n + 1 ] ) {

      ## Inner contraction step
      xc <- ( 1 - gamma )* xbar + gamma* V[ , n + 1 ]
      ## MOD
      fc <- fnm( xc, x )
      ## /MOD
      if ( is.nan( fc ) ) 
        fc <- Inf
      nf <- nf + 1
      if ( fc < f[ n + 1 ] ){
        xnew <- xc
        fnew <- fc
        happy <- 1
        ## /Inner contraction step
      } }
    if ( happy == 1 & oshrink == 1 ) {
      ## Some esoteric stuff. oshrink is set to zero and never touched again
      fbt <- mean(c(f[1:n], fnew))
      delfb <- fbt - fbc
      armtst <- alpha * sum(sgrad^2)
      if (delfb > -armtst/n) {
        if (trace) 
          cat("Trouble - restarting: \n")
        restarts <- restarts + 1
        orth <- 1
        diams <- min(diam)
        sx <- sign(0.5 * sign(sgrad))
        happy <- 0
        V[, -1] <- V[, 1]
        diag(V[, -1]) <- diag(V[, -1]) - diams * sx[1:n]
      }
    }

    ## Re-initializing of the optimization 
    if ( happy == 1 ) {
      ## Overwriting the worst point with the newly found one
      V[ , n + 1 ] <- xnew
      f[ n + 1 ] <- fnew
      ord <- order( f )
      V <- V[ , ord ]
      f <- f[ ord ]
    } else if ( happy == 0 & restarts < restarts.max ) {
      ## Hmm. But since oshrink is strictly zero there is no way in updating 'restart'
      if ( orth == 0 )
        ## This variable has no influence at all. Quite some bugs in here.
        orth <- 1
      V[ , -1 ] <- V[ , 1 ] - sigma * ( V[ , -1 ] - V[ , 1 ] )
      ## MOD
      for ( jj in 2 : ncol( V ) )
        f[ jj ] <- fnm( V[ , jj ], x )
      ## /MOD
      nf <- nf + n
      ord <- order( f )
      V <- V[ , ord ]
      f <- f[ ord ] }
    v <- V[ , -1 ] - V[ , 1 ]
    delf <- f[ - 1] - f[ 1 ]
    diam <- sqrt( colSums( v^2 ) )
    simplex.size <- sum( abs( v ) )/ max( 1, sum( abs( V[ , 1 ] ) ) )
    f[ is.nan( f ) ] <- Inf
    dist <- f[ n + 1 ] - f[ 1 ]
    sgrad <- c( crossprod( t( v ), delf ) )
    ## MOD
    x.new.vector[ itc, ] <- c( xnew, itc )
    ## /MOD
    if ( trace & !( itc%%2 ) ) 
      cat( "iter: ", itc, "\n", "value: ", f[ 1 ], "\n")
  }
  if ( dist <= ftol | simplex.size <= 1e-06 ){
    conv <- 0
    message <- "Successful convergence"
  }
  else if ( nf >= maxfeval ) {
    conv <- 1
    message <- "Maximum number of fevals exceeded"
  }
  else if ( restarts >= restarts.max ) {
    conv <- 2
    message <- "Stagnation in Nelder-Mead"
  }

  res <- list( par = V[ , 1 ], value = f[ 1 ]* ( -1 )^maximize,
              feval = nf, restarts = restarts, convergence = conv,
              message = message, x.updates = x.new.vector,
              simplex.final = V )
  class( res ) <- c( "list", "climex.gev.fit" )
  return( res )        
}
