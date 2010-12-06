\name{sumexact, cumsumexact}
\alias{sumexact}
\alias{cumsumexact}
\title{Basic Sum Operations without Round-off Errors}
\description{
  Functions for performing basic sum operations without round-off errors
}
\usage{
  sumexact(..., na.rm = FALSE)
  cumsumexact(x)
}

\arguments{
  \item{x}{numeric vector}
  \item{...}{numeric vector(s), numbers or other objects to be summed}
  \item{na.rm}{logical. Should missing values be removed?}
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
  Function \code{sumexact} returns single number. Function \code{cumsumexact} 
  returns vector of the same length as \code{x}. 
}

\references{
  Round-off error correction is based on:
  Shewchuk, Jonathan, \emph{Adaptive Precision Floating-Point Arithmetic and  
    Fast Robust Geometric Predicates},  
   \url{http://www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps}
   and its implementation found at:
   \url{http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/393090} 
   
   McCullough, D.B., (1998) \emph{Assessing the Reliability of Statistical 
     Software, Part I}, The American Statistician, Vol. 52 No 4, 
     \url{http://www.amstat.org/publications/tas/mccull-1.pdf}
   
   McCullough, D.B., (1999) \emph{Assessing the Reliability of Statistical 
     Software, Part II}, The American Statistician, Vol. 53 No 2
     \url{http://www.amstat.org/publications/tas/mccull.pdf}
     
   NIST Statistical Reference Datasets (StRD) website 
   \url{http://www.nist.gov/itl/div898/strd}

} 

\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com} based on
  code by Vadim Ogranovich, which is based on algorithms described in 
  references, pointed out by Gabor Grothendieck. 
} 

\seealso{
 \itemize{       
   \item \code{\link[base]{sum}} is faster but not error-save version of 
         \code{sumexact}
   \item \code{\link[base]{cumsum}} is equivalent to \code{cumsumexact}
 }
}

\examples{
  x = c(1, 1e20, 1e40, -1e40, -1e20, -1)
  a = sum(x);         print(a)
  b = sumexact(x);    print(b)
  stopifnot(b==0)
  a = cumsum(x);      print(a)
  b = cumsumexact(x); print(b)
  stopifnot(b[6]==0)
}
\keyword{ts}
\keyword{smooth}
\keyword{array}
\keyword{utilities}
\concept{sum with no round-off errors}
\concept{cumulative sum with no round-off errors}



