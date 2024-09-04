\name{runquantile}
\alias{runquantile}
\title{Quantile of Moving Window}
\description{Moving (aka running, rolling) Window Quantile calculated over a vector}
\usage{
  runquantile(x, k, probs, type=7, 
         endrule=c("quantile", "NA", "trim", "keep", "constant", "func"),
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
       \item \code{"quantile"} - applies the \code{\link{quantile}} 
         function to smaller and smaller sections of the array. Equivalent to: 
         \code{for(i in 1:k2) out[i]=quantile(x[1:(i+k2)])}. 
       \item \code{"trim"} - trim the ends; output array length is equal to 
         \code{length(x)-2*k2 (out = out[(k2+1):(n-k2)])}. This option mimics 
         output of \code{\link{apply}} \code{(\link{embed}(x,k),1,FUN)} and other 
         related functions.
       \item \code{"keep"} - fill the ends with numbers from \code{x} vector 
         \code{(out[1:k2] = x[1:k2])}
       \item \code{"constant"} - fill the ends with first and last calculated 
         value in output array \code{(out[1:k2] = out[k2+1])}
       \item \code{"NA"} - fill the ends with NA's \code{(out[1:k2] = NA)}
       \item \code{"func"} - same as \code{"quantile"} but implimented
       in R. This option could be very slow, and is included mostly for testing
     }
     Similar to \code{endrule} in \code{\link{runmed}} function which has the 
     following options: \dQuote{\code{c("median", "keep", "constant")}} .
  }
  \item{probs}{numeric vector of probabilities with values in [0,1] range 
    used by \code{runquantile}. For example \code{Probs=c(0,0.5,1)} would be 
    equivalent to running \code{runmin}, \code{\link{runmed}} and \code{runmax}.
    Same as \code{probs} in \code{\link{quantile}} function. }
  \item{type}{an integer between 1 and 9 selecting one of the nine quantile 
    algorithms, same as \code{type} in \code{\link{quantile}} function. 
    Another even more readable description of nine ways to calculate quantiles 
    can be found at \url{http://mathworld.wolfram.com/Quantile.html}. }
  \item{align}{specifies whether result should be centered (default), 
  left-aligned or right-aligned.  If \code{endrule}="quantile" then setting
  \code{align} to "left" or "right" will fall back on slower implementation 
  equivalent to \code{endrule}="func". }
}

\details{
  Apart from the end values, the result of y = runquantile(x, k) is the same as 
  \dQuote{\code{for(j=(1+k2):(n-k2)) y[j]=quintile(x[(j-k2):(j+k2)],na.rm = TRUE)}}. It can handle 
non-finite numbers like NaN's and Inf's (like \code{\link{quantile}(x, na.rm = TRUE)}).

  The main incentive to write this set of functions was relative slowness of 
  majority of moving window functions available in R and its packages.  With the 
  exception of \code{\link{runmed}}, a running window median function, all 
  functions listed in "see also" section are slower than very inefficient 
  \dQuote{\code{\link{apply}(\link{embed}(x,k),1,FUN)}} approach. Relative 
  speeds of \code{runquantile} is O(n*k)

  Functions \code{runquantile} and \code{runmad} are using insertion sort to 
  sort the moving window, but gain speed by remembering results of the previous 
  sort. Since each time the window is moved, only one point changes, all but one 
  points in the window are already sorted. Insertion sort can fix that in O(k) 
  time.
}

\value{
  If \code{x} is a matrix than function \code{runquantile} returns a matrix of 
  size [n \eqn{\times}{x} \code{\link{length}}(probs)]. If \code{x} is vactor 
  a than function \code{runquantile} returns a matrix of size 
  [\code{\link{dim}}(x) \eqn{\times}{x} \code{\link{length}}(probs)]. 
  If \code{endrule="trim"} the output will have fewer rows. 
}

\references{
  \itemize{       
    \item About quantiles: Hyndman, R. J. and Fan, Y. (1996) \emph{Sample 
       quantiles in statistical packages, American Statistician}, 50, 361. 
    \item About quantiles: Eric W. Weisstein. \emph{Quantile}. From MathWorld-- 
     A Wolfram Web Resource. \url{http://mathworld.wolfram.com/Quantile.html} 
    
  \item About insertion sort used in \code{runmad} and \code{runquantile}: 
  R. Sedgewick (1988): \emph{Algorithms}. Addison-Wesley (page 99)
  }
} 

\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com}}

\seealso{
  Links related to:
  \itemize{       
   \item Running Quantile - \code{\link{quantile}}, \code{\link{runmed}}, 
     \code{\link{smooth}}, \code{\link[zoo]{rollmedian}} from 
     \pkg{zoo} library
   \item Other moving window functions  from this package: \code{\link{runmin}}, 
     \code{\link{runmax}}, \code{\link{runmean}}, \code{\link{runmad}} and
     \code{\link{runsd}}
   \item Running Minimum - \code{\link{min}}
   \item Running Maximum - \code{\link{max}}, \code{\link[zoo]{rollmax}} from \pkg{zoo} library
   \item generic running window functions: \code{\link{apply}}\code{
     (\link{embed}(x,k), 1, FUN)} (fastest), \code{\link[gtools]{running}} from \pkg{gtools} 
     package (extremely slow for this purpose), \code{\link[magic]{subsums}} from 
     \pkg{magic} library can perform running window operations on data with any 
     dimensions. 
  }
}

\examples{
  # show plot using runquantile
  k=31; n=200;
  x = rnorm(n,sd=30) + abs(seq(n)-n/4)
  y=runquantile(x, k, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
  col = c("black", "red", "green", "blue", "magenta", "cyan")
  plot(x, col=col[1], main = "Moving Window Quantiles")
  lines(y[,1], col=col[2])
  lines(y[,2], col=col[3])
  lines(y[,3], col=col[4])
  lines(y[,4], col=col[5])
  lines(y[,5], col=col[6])
  lab = c("data", "runquantile(.05)", "runquantile(.25)", "runquantile(0.5)", 
          "runquantile(.75)", "runquantile(.95)")
  legend(0,230, lab, col=col, lty=1 )

  # show plot using runquantile
  k=15; n=200;
  x = rnorm(n,sd=30) + abs(seq(n)-n/4)
  y=runquantile(x, k, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
  col = c("black", "red", "green", "blue", "magenta", "cyan")
  plot(x, col=col[1], main = "Moving Window Quantiles (smoothed)")
  lines(runmean(y[,1],k), col=col[2])
  lines(runmean(y[,2],k), col=col[3])
  lines(runmean(y[,3],k), col=col[4])
  lines(runmean(y[,4],k), col=col[5])
  lines(runmean(y[,5],k), col=col[6])
  lab = c("data", "runquantile(.05)", "runquantile(.25)", "runquantile(0.5)", 
          "runquantile(.75)", "runquantile(.95)")
  legend(0,230, lab, col=col, lty=1 )
  
  # basic tests against runmin & runmax
  y = runquantile(x, k, probs=c(0, 1))
  a = runmin(x,k) # test only the inner part 
  stopifnot(all(a==y[,1], na.rm=TRUE));
  a = runmax(x,k) # test only the inner part
  stopifnot(all(a==y[,2], na.rm=TRUE));
  
  # basic tests against runmed, including testing endrules
  a = runquantile(x, k, probs=0.5, endrule="keep")
  b = runmed(x, k, endrule="keep")
  stopifnot(all(a==b, na.rm=TRUE));
  a = runquantile(x, k, probs=0.5, endrule="constant")
  b = runmed(x, k, endrule="constant")
  stopifnot(all(a==b, na.rm=TRUE));

  # basic tests against apply/embed
  a = runquantile(x,k, c(0.3, 0.7), endrule="trim")
  b = t(apply(embed(x,k), 1, quantile, probs = c(0.3, 0.7)))
  eps = .Machine$double.eps ^ 0.5
  stopifnot(all(abs(a-b)<eps));
  
  # test against loop approach
  # this test works fine at the R prompt but fails during package check - need to investigate
  k=25; n=200;
  x = rnorm(n,sd=30) + abs(seq(n)-n/4) # create random data
  x[seq(1,n,11)] = NaN;                # add NANs
  k2 = k%/%2
  k1 = k-k2-1
  a = runquantile(x, k, probs=c(0.3, 0.8) )
  b = matrix(0,n,2);
  for(j in 1:n) {
    lo = max(1, j-k1)
    hi = min(n, j+k2)
    b[j,] = quantile(x[lo:hi], probs=c(0.3, 0.8), na.rm = TRUE)
  }
  #stopifnot(all(abs(a-b)<eps));
  
  # compare calculation of array ends
  a = runquantile(x, k, probs=0.4, endrule="quantile") # fast C code
  b = runquantile(x, k, probs=0.4, endrule="func")     # slow R code
  stopifnot(all(abs(a-b)<eps));
  
  # test if moving windows forward and backward gives the same results
  k=51;
  a = runquantile(x     , k, probs=0.4)
  b = runquantile(x[n:1], k, probs=0.4)
  stopifnot(all(a[n:1]==b, na.rm=TRUE));

  # test vector vs. matrix inputs, especially for the edge handling
  nRow=200; k=25; nCol=10
  x = rnorm(nRow,sd=30) + abs(seq(nRow)-n/4)
  x[seq(1,nRow,10)] = NaN;              # add NANs
  X = matrix(rep(x, nCol ), nRow, nCol) # replicate x in columns of X
  a = runquantile(x, k, probs=0.6)
  b = runquantile(X, k, probs=0.6)
  stopifnot(all(abs(a-b[,1])<eps));        # vector vs. 2D array
  stopifnot(all(abs(b[,1]-b[,nCol])<eps)); # compare rows within 2D array

  # Exhaustive testing of runquantile to standard R approach
  numeric.test = function (x, k) {
    probs=c(1, 25, 50, 75, 99)/100
    a = runquantile(x,k, c(0.3, 0.7), endrule="trim")
    b = t(apply(embed(x,k), 1, quantile, probs = c(0.3, 0.7), na.rm=TRUE))
    eps = .Machine$double.eps ^ 0.5
    stopifnot(all(abs(a-b)<eps));
  }
  n=50;
  x = rnorm(n,sd=30) + abs(seq(n)-n/4) # nice behaving data
  for(i in 2:5) numeric.test(x, i)     # test small window sizes
  for(i in 1:5) numeric.test(x, n-i+1) # test large window size
  x[seq(1,50,10)] = NaN;               # add NANs and repet the test
  for(i in 2:5) numeric.test(x, i)     # test small window sizes
  for(i in 1:5) numeric.test(x, n-i+1) # test large window size
  
  # Speed comparison
  \dontrun{
  x=runif(1e6); k=1e3+1;
  system.time(runquantile(x,k,0.5))    # Speed O(n*k)
  system.time(runmed(x,k))             # Speed O(n * log(k)) 
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
\concept{moving quantile}
\concept{rolling quantile}
\concept{running quantile}
\concept{moving percentile}
\concept{rolling percentile}
\concept{running percentile}
\concept{moving median}
\concept{rolling median}
\concept{running median}
\concept{moving window}
\concept{rolling window}
\concept{running window}

