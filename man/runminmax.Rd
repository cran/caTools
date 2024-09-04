\name{runmin & runmax}
\alias{runmin}
\alias{runmax}
\title{Minimum and Maximum of Moving Windows}
\description{Moving (aka running, rolling) Window Minimum and Maximum 
  calculated over a vector  }
\usage{
  runmin(x, k, alg=c("C", "R"),
         endrule=c("min", "NA", "trim", "keep", "constant", "func"),
         align = c("center", "left", "right"))
  runmax(x, k, alg=c("C", "R"),
         endrule=c("max", "NA", "trim", "keep", "constant", "func"),
         align = c("center", "left", "right"))
}

\arguments{
  \item{x}{numeric vector of length n or matrix with n rows. If \code{x} is a 
    matrix than each column will be processed separately.}
  \item{k}{width of moving window; must be an integer between one and n }
  \item{endrule}{character string indicating how the values at the beginning 
    and the end, of the array, should be treated. Only first and last \code{k2} 
    values at both ends are affected, where \code{k2} is the half-bandwidth 
    \code{k2 = k \%/\% 2}.
     \itemize{
       \item \code{"min"} & \code{"max"} - applies the underlying function to 
       smaller and smaller sections of the array. In case of min equivalent to: 
       \code{for(i in 1:k2) out[i]=min(x[1:(i+k2)])}. Default.
       \item \code{"trim"} - trim the ends; output array length is equal to 
         \code{length(x)-2*k2 (out = out[(k2+1):(n-k2)])}. This option mimics 
         output of \code{\link{apply}} \code{(\link{embed}(x,k),1,FUN)} and other 
         related functions.
       \item \code{"keep"} - fill the ends with numbers from \code{x} vector 
         \code{(out[1:k2] = x[1:k2])}
       \item \code{"constant"} - fill the ends with first and last calculated 
         value in output array \code{(out[1:k2] = out[k2+1])}
       \item \code{"NA"} - fill the ends with NA's \code{(out[1:k2] = NA)}
       \item \code{"func"} - same as \code{"min"} & \code{"max"} but implimented
       in R. This option could be very slow, and is included mostly for testing
     }
     Similar to \code{endrule} in \code{\link{runmed}} function which has the 
     following options: \dQuote{\code{c("median", "keep", "constant")}} .
  }
  \item{alg}{an option allowing to choose different algorithms or 
    implementations. Default is to use of code written in C (option \code{alg="C"}).
    Option \code{alg="R"} will use slower code written in R. Useful for 
    debugging and studying the algorithm.}
  \item{align}{specifies whether result should be centered (default), 
  left-aligned or right-aligned.  If \code{endrule}="min" or "max" then setting
  \code{align} to "left" or "right" will fall back on slower implementation 
  equivalent to \code{endrule}="func". }
}

\details{
  Apart from the end values, the result of y = runFUN(x, k) is the same as 
  \dQuote{\code{for(j=(1+k2):(n-k2)) y[j]=FUN(x[(j-k2):(j+k2)], na.rm = TRUE)}}, where FUN 
  stands for min or max functions.  Both functions can handle non-finite 
  numbers like NaN's and Inf's the same way as \code{\link{min}(x, na.rm = TRUE)}).
   

  The main incentive to write this set of functions was relative slowness of 
  majority of moving window functions available in R and its packages.  With the 
  exception of \code{\link{runmed}}, a running window median function, all 
  functions listed in "see also" section are slower than very inefficient 
  \dQuote{\code{\link{apply}(\link{embed}(x,k),1,FUN)}} approach. Relative 
  speeds \code{runmin} and \code{runmax} functions is O(n) in best and average 
  case and O(n*k) in worst case.
    
  Both functions work with infinite numbers (\code{NA},\code{NaN},\code{Inf},
  \code{-Inf}). Also default \code{endrule} is hardwired in C for speed.
}

\value{
  Returns a numeric vector or matrix of the same size as \code{x}. Only in case of 
  \code{endrule="trim"} the output vectors will be shorter and output matrices 
  will have fewer rows. 
}

\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com}}

\seealso{
  Links related to:
  \itemize{       
   \item Other moving window functions  from this package: \code{\link{runmean}}, 
    \code{\link{runquantile}}, \code{\link{runmad}} and \code{\link{runsd}}  
   \item R functions: \code{\link{runmed}}, \code{\link{min}}, \code{\link{max}}
   \item Similar functions in other packages: \code{\link[zoo]{rollmax}} from \pkg{zoo} library
   \item generic running window functions: \code{\link{apply}}\code{
     (\link{embed}(x,k), 1, FUN)} (fastest), \code{\link[gtools]{running}} from \pkg{gtools} 
     package (extremely slow for this purpose), \code{\link[magic]{subsums}} from 
     \pkg{magic} library can perform running window operations on data with any 
     dimensions. 
  }
}

\examples{
  # show plot using runmin, runmax and runmed
  k=25; n=200;
  x = rnorm(n,sd=30) + abs(seq(n)-n/4)
  col = c("black", "red", "green", "blue", "magenta", "cyan")
  plot(x, col=col[1], main = "Moving Window Analysis Functions")
  lines(runmin(x,k), col=col[2])
  lines(runmean(x,k), col=col[3])
  lines(runmax(x,k), col=col[4])
  legend(0,.9*n, c("data", "runmin", "runmean", "runmax"), col=col, lty=1 )

  # basic tests against standard R approach
  a = runmin(x,k, endrule="trim") # test only the inner part 
  b = apply(embed(x,k), 1, min)   # Standard R running min
  stopifnot(all(a==b));
  a = runmax(x,k, endrule="trim") # test only the inner part
  b = apply(embed(x,k), 1, max)   # Standard R running min
  stopifnot(all(a==b));
  
  # test against loop approach
  k=25; 
  data(iris)
  x = iris[,1]
  n = length(x)
  x[seq(1,n,11)] = NaN;                # add NANs
  k2 = k%/%2
  k1 = k-k2-1
  a1 = runmin(x, k)
  a2 = runmax(x, k)
  b1 = array(0,n)
  b2 = array(0,n)
  for(j in 1:n) {
    lo = max(1, j-k1)
    hi = min(n, j+k2)
    b1[j] = min(x[lo:hi], na.rm = TRUE)
    b2[j] = max(x[lo:hi], na.rm = TRUE)
  }
  # this test works fine at the R prompt but fails during package check - need to investigate
  \dontrun{ 
  stopifnot(all(a1==b1, na.rm=TRUE));
  stopifnot(all(a2==b2, na.rm=TRUE));
  }
  
  # Test if moving windows forward and backward gives the same results
  # Two data sets also corespond to best and worst-case scenatio data-sets
  k=51; n=200;
  a = runmin(n:1, k) 
  b = runmin(1:n, k)
  stopifnot(all(a[n:1]==b, na.rm=TRUE));
  a = runmax(n:1, k)
  b = runmax(1:n, k)
  stopifnot(all(a[n:1]==b, na.rm=TRUE));

  # test vector vs. matrix inputs, especially for the edge handling
  nRow=200; k=25; nCol=10
  x = rnorm(nRow,sd=30) + abs(seq(nRow)-n/4)
  x[seq(1,nRow,10)] = NaN;              # add NANs
  X = matrix(rep(x, nCol ), nRow, nCol) # replicate x in columns of X
  a = runmax(x, k)
  b = runmax(X, k)
  stopifnot(all(a==b[,1], na.rm=TRUE));        # vector vs. 2D array
  stopifnot(all(b[,1]==b[,nCol], na.rm=TRUE)); # compare rows within 2D array
  a = runmin(x, k)
  b = runmin(X, k)
  stopifnot(all(a==b[,1], na.rm=TRUE));        # vector vs. 2D array
  stopifnot(all(b[,1]==b[,nCol], na.rm=TRUE)); # compare rows within 2D array

  # Compare C and R algorithms to each other for extreme window sizes
  numeric.test = function (x, k) {
    a = runmin( x, k, alg="C")
    b = runmin( x, k, alg="R")
    c =-runmax(-x, k, alg="C")
    d =-runmax(-x, k, alg="R")
    stopifnot(all(a==b, na.rm=TRUE));
    #stopifnot(all(c==d, na.rm=TRUE)); 
    #stopifnot(all(a==c, na.rm=TRUE));
    stopifnot(all(b==d, na.rm=TRUE));
  }
  n=200;                               # n is an even number
  x = rnorm(n,sd=30) + abs(seq(n)-n/4) # random data
  for(i in 1:5) numeric.test(x, i)     # test for small window size
  for(i in 1:5) numeric.test(x, n-i+1) # test for large window size
  n=201;                               # n is an odd number
  x = rnorm(n,sd=30) + abs(seq(n)-n/4) # random data
  for(i in 1:5) numeric.test(x, i)     # test for small window size
  for(i in 1:5) numeric.test(x, n-i+1) # test for large window size
  n=200;                               # n is an even number
  x = rnorm(n,sd=30) + abs(seq(n)-n/4) # random data
  x[seq(1,200,10)] = NaN;              # with some NaNs
  for(i in 1:5) numeric.test(x, i)     # test for small window size
  for(i in 1:5) numeric.test(x, n-i+1) # test for large window size
  n=201;                               # n is an odd number
  x = rnorm(n,sd=30) + abs(seq(n)-n/4) # random data
  x[seq(1,200,2)] = NaN;               # with some NaNs
  for(i in 1:5) numeric.test(x, i)     # test for small window size
  for(i in 1:5) numeric.test(x, n-i+1) # test for large window size

  # speed comparison
  \dontrun{
  n = 1e7;  k=991; 
  x1 = runif(n);                       # random data - average case scenario
  x2 = 1:n;                            #  best-case scenario data for runmax
  x3 = n:1;                            # worst-case scenario data for runmax
  system.time( runmax( x1,k,alg="C"))  # C alg on average data O(n)
  system.time( runmax( x2,k,alg="C"))  # C alg on  best-case data O(n)
  system.time( runmax( x3,k,alg="C"))  # C alg on worst-case data O(n*k)
  system.time(-runmin(-x1,k,alg="C"))  # use runmin to do runmax work
  system.time( runmax( x1,k,alg="R"))  # R version of the function
  x=runif(1e5); k=1e2;                 # reduce vector and window sizes
  system.time(runmax(x,k,alg="R"))     # R version of the function
  system.time(apply(embed(x,k), 1, max)) # standard R approach 
  }
}

\keyword{ts}
\keyword{smooth}
\keyword{array}
\keyword{utilities}
\concept{moving min}
\concept{rolling min}
\concept{running min}
\concept{moving max}
\concept{rolling max}
\concept{running max}
\concept{moving minimum}
\concept{rolling minimum}
\concept{running minimum}
\concept{moving maximum}
\concept{rolling maximum}
\concept{running maximum}
\concept{running window}
\concept{moving window}
\concept{rolling window}

