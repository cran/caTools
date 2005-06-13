\name{sum.exact}
\alias{sum.exact}
\alias{cumsum.exact}
\alias{runsum.exact}
\title{Basic Sum Operations without Round-off Errors}
\description{
  Functions for performing basic sum operations without round-off errors
}
\usage{
  sum.exact(..., na.rm = FALSE)
  cumsum.exact(x)
  runsum.exact(x,k)
}

\arguments{
  \item{x}{numeric vector}
  \item{...}{numeric vector(s), numbers or other objects to be summed}
  \item{na.rm}{logical. Should missing values be removed?}
  \item{k}{width of moving window; must be an odd integer between one and n }
}

\details{
 All three functions use full precision summation using multiple doubles for 
 intermediate values. The sum of numbers x & y is a=x+y with error term 
 b=error(a+b). That way a+b is equal exactly x+y, so sum of 2 numbers is stored
 as 2 or fewer values, which when added would under-flow. By extension sum of n 
 numbers is calculated with intermediate results stored as array of numbers 
 that can not be added without introducing an error. Only final result is 
 converted to a single number
}

\value{
  Function \code{sum.exact} returns single number. Function \code{cumsum.exact} 
  returns vector of the same length as \code{x}. Function \code{runsum.exact} 
  returns vector of length \code{length(x)-k} and attribute "count" containing
  number of finite (as in \code{\link{is.finite}}) elements in each window. 
}

\references{
  Round-off error correction is based on:
  Shewchuk, Jonathan, \emph{Adaptive Precision Floating-Point Arithmetic and  
    Fast Robust Geometric Predicates},  
   \url{http://www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps}
   and its implementation found at:
   \url{http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/393090}  
} 

\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com} based on
  code by Vadim Ogranovich, which is based on algorithms described in 
  references, pointed out by Gabor Grothendieck. 
} 

\seealso{
 \itemize{       
   \item \code{sum.exact} - is equivalent to \code{\link{sum}}
   \item \code{cumsum.exact} - is equivalent to \code{\link{cumsum}}
   \item \code{runsum.exact} - is similar to \code{\link{runmean}
      (x,k,endrule="trim")}
 }
}

\examples{
  x = c(1, 1e20, 1e40, -1e40, -1e20, -1)
  sum(x)
  sum.exact(x)
  cumsum(x)
  cumsum.exact(x)
}
\keyword{ts}
\keyword{smooth}
\keyword{array}
\keyword{utilities}
