\name{nmk}
\alias{nmk}
\title{
A Nelder-Mead optimziation algorithm for derivative-free optimization 
}
\description{
An implementation of the Nelder-Mead algorithm for derivative-free optimization.  This implementation is based on the Matlab code of Prof. C.T. Kelley (with his permission).  It generally has better convergence than the Nelder-Mead implemenation in \code{optim}. 
Nelder-Mead is not recommended for high-dimensional (> 20 parameters) optimization problems. 
}
\usage{
nmk(par, fn, control = list(), ...)
}
\arguments{
\item{par}{Starting vector of parameter values}
\item{fn}{
Nonlinear objective function that is to be optimized. 
      A scalar function that takes a real vector as argument and 
      returns a scalar that is the value of the function at that point 
      (see details).}
  \item{control}{A list of control parameters.  See *Details* for more information.
}
  \item{\dots}{Additional arguments passed to \code{fn}
}
}
\details{
Argument \code{control} is a list specifing any changes to default values of algorithm control parameters for the outer loop.  Note that the names of these must be specified completely.  Partial matching will not work.  The list items are as follows:

\code{tol} Convergence tolerance.  Iteration is terminated when the absolute difference in function value between successive iteration is below \code{tol}.  Default is 1.e-06.

\code{maxfeval}: Maximum number of objective function evaluations allowed.  Default is min(5000, max(1500, 20*length(par)^2)).

\code{regsimp} A logical variable indicating whether the starting parameter configuration is a regular simplex.  Default is TRUE.

\code{maximize} A logical variable indicating whether the objective function should be maximized.  Default is FALSE.

\code{restarts.max} Maximum number of times the algorithm should be restarted before declaring failure. Default is 3.

\code{trace} A logical variable indicating whether the starting parameter configuration is a regular simplex.  Default is FALSE.
}
\value{
A list with the following components:
\item{par}{Best estimate of the parameter vector found by the algorithm.}

\item{value}{The value of the objective function at termination.}

\item{feval}{The number of times the objective \code{fn} was evaluated.
}

\item{restarts}{The number of times the algorithm had to be restarted when it stagnated.
}

\item{convergence}{An integer code indicating type of convergence.  \code{0} indicates successful convergence. Positive integer codes indicate failure to converge.    
}

\item{message}{Text message indicating the type of convergence or failure.  
}

}
\references{
C.T. Kelley (1999), Iterative Methods for Optimization, SIAM.
}
\author{
Ravi Varadhan <rvaradhan@jhmi.edu>, Johns Hopkins University
URL:http://www.jhsph.edu/agingandhealth/People/Faculty_personal_pages/Varadhan.html
}
\note{
This algorithm is based on the Matlab code of Prof. C.T. Kelley, given in his book "Iterative methods for optimization".
It is implemented here with the permission of Prof. Kelley and SIAM.  However, there are some non-trivial modifications of the algorithm.
}


\seealso{
\code{\link{optim}}
}
\examples{
 rosbkext <- function(x){
# Extended Rosenbrock function
 n <- length(x)
 sum (100*(x[1:(n-1)]^2 - x[2:n])^2 + (x[1:(n-1)] - 1)^2)
 }

np <- 10
set.seed(123)

p0 <- rnorm(np)
xm1 <- nmk(fn=rosbkext, par=p0) # maximum `fevals' is not sufficient to find correct minimum
xm2 <- nmk(fn=rosbkext, par=p0, control=list(maxfeval=5000)) # finds the correct minimum 
xm3 <- nmk(fn=rosbkext, par=p0, control=list(regsimp=FALSE, maxfeval=5000)) # also finds the correct minimum 

ans.optim <- optim(fn=rosbkext, par=p0, method="Nelder-Mead", control=list(maxit=5000))   # terminates with inferior estimates

### A non-smooth problem
nsf <- function(x) {
	f1 <- x[1]^2 + x[2]^2
	f2 <- x[1]^2 + x[2]^2 + 10 * (-4*x[1] - x[2] + 4)
	f3 <- x[1]^2 + x[2]^2 + 10 * (-x[1] - 2*x[2] + 6)
	max(f1, f2, f3)
}

p0 <- rnorm(2)
xm4 <- nmk(fn=nsf, par=p0)

}
\keyword{optimize}
