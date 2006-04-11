\name{trapz}
\alias{trapz}
\title{Trapezoid Rule Numerical Integration}
\description{Computes the integral of Y with respect to X using trapezoid rule 
  integration.}
\usage{trapz(x, y)}

\arguments{
  \item{x}{ Sorted vector of x-axis values. } 
  \item{y}{ Vector of y-axis values. }
}

\details{
  The function has only two lines:
  \preformatted{
    idx = 2:length(x)
    return (as.double( (x[idx] - x[idx-1]) \%*\% (y[idx] + y[idx-1])) / 2)
  }
}

\value{Integral of Y with respect to X or area under the Y curve.} 

\note{
  Trapezoid rule is not the most accurate way of calculating integrals (it is 
  exact for linear functions), for example Simpson's rule (exact for linear and 
  quadratic functions) is more accurate.
} 

\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com}} 
\references{ D. Kincaid & W. Chaney (1991), \emph{Numerical Analysis}, p.445 }

\seealso{
  \itemize{
   \item \code{\link[PROcess]{intg}} from \pkg{PROcess} package
   \item \code{\link[ROC]{trapezint}} from \pkg{ROC} package
   \item \code{\link{integrate}}
   \item Matlab's \code{trapz} function (\url{
     http://www.mathworks.com/access/helpdesk/help/techdoc/ref/trapz.html}) 
 }
}

\examples{
  # integral of sine function in [0, pi] range suppose to be exactly 2.
  # lets calculate it using 10 samples:
  x = (1:10)*pi/10
  trapz(x, sin(x))
  # now lets calculate it using 1000 samples:
  x = (1:1000)*pi/1000
  trapz(x, sin(x))
}

\keyword{math}
\concept{integration}
