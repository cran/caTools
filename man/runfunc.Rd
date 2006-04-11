\name{runmean, runmin, runmax & runquantile}
\alias{runmean}
\alias{runmin}
\alias{runmax}
\alias{runmad}
\alias{runquantile}
\alias{EndRule}
\title{Moving Window Analysis of a Vector}
\description{A collection of functions to perform moving window (running, 
  rolling window) analysis of vectors  }
\usage{
  runmean(x, k, alg=c("C", "R", "exact"),
         endrule=c("NA", "trim", "keep", "constant", "func"))
  runmin (x, k, alg=c("C", "R"),
         endrule=c("NA", "trim", "keep", "constant", "func"))
  runmax (x, k, alg=c("C", "R"),
         endrule=c("NA", "trim", "keep", "constant", "func"))
  runmad (x, k, center=runmed(x,k,endrule="keep"), constant=1.4826,  
          endrule=c("NA", "trim", "keep", "constant", "func"))
  runquantile(x, k, probs, type=7, 
          endrule=c("NA", "trim", "keep", "constant", "func"))
  EndRule(x, y, k, 
          endrule=c("NA", "trim", "keep", "constant", "func"), Func, \dots)
}

\arguments{
  \item{x}{numeric vector of length n}
  \item{k}{width of moving window; must be an odd integer between three and n }
  \item{endrule}{character string indicating how the values at the beginning 
    and the end, of the data, should be treated. Only first and last \code{k2} 
    values at both ends are affected, where \code{k2} is the half-bandwidth 
    \code{k2 = k \%/\% 2}.
     \itemize{
       \item \code{"trim"} - trim the ends output array length is equal to 
         \code{length(x)-2*k2 (out = out[(k2+1):(n-k2)])}. This option mimics 
         output of \code{\link{apply}} \code{(\link{embed}(x,k),1,FUN)} and other 
         related functions.
       \item \code{"keep"} - fill the ends with numbers from \code{x} vector 
         \code{(out[1:k2] = x[1:k2])}
       \item \code{"constant"} - fill the ends with first and last calculated 
         value in output array \code{(out[1:k2] = out[k2+1])}
       \item \code{"NA"} - fill the ends with NA's \code{(out[1:k2] = NA)}
       \item \code{"func"} - applies the underlying function to smaller and 
       smaller sections of the array. For example in case of mean: 
       \code{for(i in 1:k2) out[i]=mean(x[1:i])}. This option is not optimized 
       and could be very slow for large windows.
     }
     Similar to \code{endrule} in \code{\link{runmed}} function which has the 
     following options: \dQuote{\code{c("median", "keep", "constant")}} .
  }
  \item{alg}{an option allowing to choose different algorithms or 
    implementations, if provided. Default is to use of code written in C. 
    Option \code{alg="R"} will use slower code written in R. Useful for 
    debugging and allows extensions in the future.}
  \item{center}{moving window center used by \code{runmad} function defaults 
    to running median (\code{\link{runmed}} function). Similar to \code{center}  
    in \code{\link{mad}} function. }
  \item{constant}{scale factor used by \code{runmad}, such that for Gaussian 
    distribution X, \code{\link{mad}}(X) is the same as \code{\link{sd}}(X). 
    Same as \code{constant} in \code{\link{mad}} function.}
  \item{probs}{numeric vector of probabilities with values in [0,1] range 
    used by \code{runquantile}. For example \code{Probs=c(0,0.5,1)} would be 
    equivalent to running \code{runmin}, \code{\link{runmed}} and \code{runmax}.
    Same as \code{probs} in \code{\link{quantile}} function. }
  \item{type}{an integer between 1 and 9 selecting one of the nine quantile 
    algorithms, same as \code{type} in \code{\link{quantile}} function. 
    Another even more readable description of nine ways to calculate quantiles 
    can be found at \url{http://mathworld.wolfram.com/Quantile.html}. }
  \item{y}{numeric vector of length n, which is partially filled output of 
    one of the \code{run} functions. Function \code{EndRule} will fill the 
    remaining beginning and end sections using method chosen by \code{endrule} 
    argument.}
  \item{Func}{Function name that \code{EndRule} will use in case of 
    \code{endrule="func"}.}
  \item{\dots}{Additional parameters to \code{Func} that \code{EndRule} will 
    use in case of \code{endrule="func"}.}
}

\details{
  Apart from the end values, the result of y = runFUN(x, k) is the same as 
  \dQuote{\code{for(j=(1+k2):(n-k2)) y[j]=FUN(x[(j-k2):(j+k2)])}}, where FUN 
  stands for min, max, mean, mad or quantile functions.

  The main incentive to write this set of functions was relative slowness of 
  majority of moving window functions available in R and its packages.  With 
  exception of \code{\link{runmed}}, a running window median function, all 
  functions listed in "see also" section are slower than very inefficient 
  \dQuote{\code{\link{apply}(\link{embed}(x,k),1,FUN)}} approach. Relative 
  speeds of above functions are as follow:
     \itemize{
       \item \code{runmin}, \code{runmax}, \code{runmean} run at O(n)
       \item \code{runmean(..., alg="exact")} can have worst case speed of 
       O(\eqn{n^2}) for some small data vectors, but average case is still 
         close to O(n).
       \item \code{runquantile} and \code{runmad} run at O(n*k)
       \item \code{\link{runmed}} - related R function runs at O(n*log(k)) 
      }
  
  Functions \code{runquantile} and \code{runmad} are using insertion sort to 
  sort the moving window, but gain speed by remembering results of the previous 
  sort. Since each time the window is moved, only one point changes, all but one 
  points in the window are already sorted. Insertion sort can fix that in O(k) 
  time.
  
  Function \code{runquantile} when run in single probability mode automatically
  recognizes probabilities: 0, 1/2, and 1 as special cases and return output 
  from functions: \code{runmin}, \code{\link{runmed}} and \code{runmax} 
  respectively. 

  All \code{run*} functions are written in C, but \code{runmin}, \code{runmax} 
  and \code{runmean} also have fast R code versions (see argument 
  \code{alg="R"}). Those were included for debugging purposes, and as a fallback 
  in hard-to-port situations. See examples.
  
  Function \code{EndRule} applies one of the five methods (see \code{endrule} 
  argument) to process end-points of the input array \code{x}. 
  
  In case of \code{runmean(..., alg="exact")} function a special algorithm is 
  used (see references section) to ensure that round-off errors do not 
  accumulate. As a result \code{runmean} is more accurate than 
  \code{\link{filter}}(x, rep(1/k,k)) and \code{runmean(..., alg="C")} 
  functions.
  
  All of the functions in this section do not work with infinite numbers
  (\code{NA},\code{NaN},\code{Inf},\code{-Inf}) except for 
  \code{runmean(..., alg="exact")} which omits them.
}

\value{
  Functions \code{runmin}, \code{runmax}, \code{runmean} and \code{runmad} 
  return a numeric vector of the same length as \code{x}. 
  Function \code{runquantile} returns a matrix of size [n \eqn{\times}{x} 
  length(probs)]. In addition \code{x} contains \code{attr}ibute \code{k}
  with (the 'oddified') \code{k}.
} 

\references{
  \itemize{       
    \item About quantiles: Hyndman, R. J. and Fan, Y. (1996) \emph{Sample 
       quantiles in statistical packages, American Statistician}, 50, 361. 
    \item About quantiles: Eric W. Weisstein. \emph{Quantile}. From MathWorld-- 
     A Wolfram Web Resource. \url{http://mathworld.wolfram.com/Quantile.html} 
    
  \item About insertion sort used in \code{runmad} and \code{runquantile}: 
  R. Sedgewick (1988): \emph{Algorithms}. Addison-Wesley (page 99)
   
  \item About round-off error correction used in \code{runmean}:
  Shewchuk, Jonathan \emph{Adaptive Precision Floating-Point Arithmetic and Fast 
    Robust Geometric Predicates},  
   \url{http://www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps}
  
   \item More on round-off error correction can be found at:
   \url{http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/393090 }  
  }
} 

\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com}}
\note{  
  Function \code{runmean(..., alg="exact")} is based by code by Vadim Ogranovich,
  which is based on Python code (see last reference), pointed out by Gabor 
  Grothendieck. 
} 

\seealso{
  Links related to each function:
  \itemize{       
   \item \code{runmean} - \code{\link{mean}}, \code{\link{kernapply}}, 
     \code{\link{filter}}, \code{\link{runsum.exact}}, \code{\link{decompose}},
     \code{\link{stl}},
     \code{\link[fSeries]{rollMean}} from \pkg{fSeries} library, 
     \code{\link[zoo]{rollmean}} from \pkg{zoo} library,
     \code{\link[magic]{subsums}} from \pkg{magic} library,
   \item \code{runmin} - \code{\link{min}}, \code{\link[fSeries]{rollMin}} from 
     \pkg{fSeries} library
   \item \code{runmax} - \code{\link{max}}, \code{\link[fSeries]{rollMax}} from 
     \pkg{fSeries} library, \code{\link[zoo]{rollmax}} from 
     \pkg{zoo} library
   \item \code{runquantile} - \code{\link{quantile}}, \code{\link{runmed}}, 
     \code{\link{smooth}}, \code{\link[zoo]{rollmedian}} from 
     \pkg{zoo} library
   \item \code{runmad} - \code{\link{mad}}, \code{\link[fSeries]{rollVar}} from 
     \pkg{fSeries} library
   \item generic running window functions: \code{\link{apply}}\code{
     (\link{embed}(x,k), 1, FUN)} (fastest), \code{\link[fSeries]{rollFun}} 
     from \pkg{fSeries} (slow), \code{\link[gtools]{running}} from \pkg{gtools} 
     package (extremely slow for this purpose), \code{\link[zoo]{rapply}} from 
     \pkg{zoo} library, \code{\link[magic]{subsums}} from 
     \pkg{magic} library can perform running window operations on data with any 
     dimensions. 
   \item \code{EndRule} - \code{\link{smoothEnds}(y,k)} function is similar to 
     \code{EndRule(x,y,k,endrule="func", median)}
  }
}

\examples{
  # test runmin, runmax and runmed
  k=15; n=200;
  x = rnorm(n,sd=30) + abs(seq(n)-n/4)
  col = c("black", "red", "green", "blue", "magenta", "cyan")
  plot(x, col=col[1], main = "Moving Window Analysis Functions")
  lines(runmin(x,k), col=col[2])
  lines(runmed(x,k), col=col[3])
  lines(runmax(x,k), col=col[4])
  legend(0,.9*n, c("data", "runmin", "runmed", "runmax"), col=col, lty=1 )

  #test runmean and runquantile
  y=runquantile(x, k, probs=c(0, 0.5, 1, 0.25, 0.75), endrule="constant")
  plot(x, col=col[1], main = "Moving Window Quantile")
  lines(runmean(y[,1],k), col=col[2])
  lines(y[,2], col=col[3])
  lines(runmean(y[,3],k), col=col[4])
  lines(y[,4], col=col[5])
  lines(y[,5], col=col[6])
  lab = c("data", "runmean(runquantile(0))", "runquantile(0.5)", 
  "runmean(runquantile(1))", "runquantile(.25)", "runquantile(.75)")
  legend(0,0.9*n, lab, col=col, lty=1 )

  #test runmean and runquantile
  k =25
  m=runmed(x, k)
  y=runmad(x, k, center=m)
  plot(x, col=col[1], main = "Moving Window Analysis Functions")
  lines(m    , col=col[2])
  lines(m-y/2, col=col[3])
  lines(m+y/2, col=col[4])
  lab = c("data", "runmed", "runmed-runmad/2", "runmed+runmad/2")
  legend(0,1.8*n, lab, col=col, lty=1 )

  # numeric comparison between different algorithms
  numeric.test = function (n, k) {
    eps = .Machine$double.eps ^ 0.5
    x = rnorm(n,sd=30) + abs(seq(n)-n/4)
    # numeric comparison : runmean
    a = runmean(x,k)
    b = runmean(x,k, alg="R")
    d = runmean(x,k, alg="exact")
    e = filter(x, rep(1/k,k))
    stopifnot(all(abs(a-b)<eps, na.rm=TRUE));
    stopifnot(all(abs(a-d)<eps, na.rm=TRUE));
    stopifnot(all(abs(a-e)<eps, na.rm=TRUE));
    # numeric comparison : runmin
    a = runmin(x,k, endrule="trim")
    b = runmin(x,k, endrule="trim", alg="R")
    c = apply(embed(x,k), 1, min)
    stopifnot(all(a==b, na.rm=TRUE));
    stopifnot(all(a==c, na.rm=TRUE));
    # numeric comparison : runmax
    a = runmax(x,k, endrule="trim")
    b = runmax(x,k, endrule="trim", alg="R")
    c = apply(embed(x,k), 1, max)
    stopifnot(all(a==b, na.rm=TRUE));
    stopifnot(all(a==c, na.rm=TRUE));
    # numeric comparison : runmad
    a = runmad(x,k, endrule="trim")
    b = apply(embed(x,k), 1, mad)
    stopifnot(all(a==b, na.rm=TRUE));
    # numeric comparison : runquantile
    a = runquantile(x,k, c(0.3, 0.7), endrule="trim")
    b = t(apply(embed(x,k), 1, quantile, probs = c(0.3, 0.7)))
    stopifnot(all(abs(a-b)<eps));
  }
  numeric.test(50, 3) # test different window size vs. vector ...
  numeric.test(50,15) # ... length combinations
  numeric.test(50,49)
  numeric.test(49,49)
  
  # speed comparison
  x=runif(100000); k=991;
  system.time(runmean(x,k))
  system.time(runmean(x,k, alg="R"))
  system.time(runmean(x,k, alg="exact"))
  system.time(filter(x, rep(1/k,k), sides=2)) #the fastest alternative I know
  k=91;
  system.time(runmad(x,k))
  system.time(apply(embed(x,k), 1, mad)) #the fastest alternative I know
  
  # numerical comparison of round-off error handling
  test.runmean = function (x, k) {
    a = k*runmean(x,k, alg="exact")  
    b = k*runmean(x,k, alg="C")  
    d = k*runmean(x,k, alg="R")
    e = k*filter(x, rep(1/k,k))
    f = k* c(NA, NA, apply(embed(x,k), 1, mean), NA, NA)
    x = cbind(x, a, b, d, e, f)
    colnames(x) = c("x","runmean(alg=exact)","runmean(alg=C)",
      "runmean(alg=R)","filter","apply")
    return(x)
  }
  a = rep( c(1, 10,    -10,    -1, 0, 0, 0), 3) # nice-behaving array
  b = rep( c(1, 10^20, -10^20, -1, 0, 0, 0), 3) # round-off error prone array
  d = rep( c(1, 10^20, 10^40, -10^40, -10^20, -1,  0), 3) 
  test.runmean(a, 5) #all runmean algorithms give the same result
  test.runmean(b, 5) #runmean(alg=R) gives wrong result
  test.runmean(d, 5) #only runmean(alg=exact) gives correct result
}

\keyword{ts}
\keyword{smooth}
\keyword{array}
\keyword{utilities}
\concept{moving mean}
\concept{rolling mean}
\concept{running mean}
\concept{moving average}
\concept{rolling average}
\concept{running average}
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
\concept{moving mad}
\concept{rolling mad}
\concept{running mad}
\concept{running window}

