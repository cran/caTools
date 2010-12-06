\name{combs}
\alias{combs}
\title{All Combinations of k Elements from Vector v}
\description{Finds all unordered combinations of \code{k} elements from vector 
  \code{v}.}
\usage{combs(v,k)}

\arguments{
  \item{v}{Any numeric vector}
  \item{k}{Number of elements to choose from vector \code{v}. Integer smaller 
    or equal than length of \code{v}.  }
}

\value{
  \code{combs(v,k)} (where \code{v} has length \code{n}) creates a matrix with 
  \eqn{\frac{n!}{(n-k)! k!}}{n!/((n-k)! k!)} (\code{n} choose \code{k}) rows
  and \code{k} columns containing all possible combinations of \code{n} elements 
  taken \code{k} at a time.
} 

\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com}} 

\seealso{
  I discovered recently that R packages already have two functions with 
  similar capabilities:
  \code{\link[gtools]{combinations}} from \pkg{gTools} package and
  Also similar to Matlab's \code{nchoosek} function (\url{
    http://www.mathworks.com/access/helpdesk/help/techdoc/ref/nchoosek.html}) 
}

\examples{
  combs(2:5, 3) # display examples
  combs(c("cats", "dogs", "mice"), 2)
  
  a = combs(1:4, 2)
  b = matrix( c(1,1,1,2,2,3,2,3,4,3,4,4), 6, 2)
  stopifnot(a==b)
}

\keyword{models}
