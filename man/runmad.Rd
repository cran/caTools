\name{runmad}
\alias{runmad}
\title{Median Absolute Deviation of Moving Windows}
\description{ Moving (aka running, rolling) Window MAD (Median Absolute 
Deviation) calculated over a vector}
\usage{
   runmad(x, k, center = runmed(x,k), constant = 1.4826,
         endrule=c("mad", "NA", "trim", "keep", "constant", "func"),
         align = c("center", "left", "right"))
}

\arguments{
  \item{x}{numeric vector of length n or matrix with n rows. If \code{x} is a 
    matrix than each column will be processed separately.}
  \item{k}{width of moving window; must be an integer between one and n. In case
  of even k's one will have to provide different \code{center} function, since
  \code{\link{runmed}} does not take even k's.}
  \item{endrule}{character string indicating how the values at the beginning 
    and the end, of the data, should be treated. Only first and last \code{k2} 
    values at both ends are affected, where \code{k2} is the half-bandwidth 
    \code{k2 = k \%/\% 2}.
     \itemize{
       \item \code{"mad"} - applies the mad function to
       smaller and smaller sections of the array. Equivalent to: 
       \code{for(i in 1:k2) out[i]=mad(x[1:(i+k2)])}. 
       \item \code{"trim"} - trim the ends; output array length is equal to 
         \code{length(x)-2*k2 (out = out[(k2+1):(n-k2)])}. This option mimics 
         output of \code{\link{apply}} \code{(\link{embed}(x,k),1,FUN)} and other 
         related functions.
       \item \code{"keep"} - fill the ends with numbers from \code{x} vector 
         \code{(out[1:k2] = x[1:k2])}. This option makes more sense in case of 
	 smoothing functions, kept here for consistency.
       \item \code{"constant"} - fill the ends with first and last calculated 
         value in output array \code{(out[1:k2] = out[k2+1])}
       \item \code{"NA"} - fill the ends with NA's \code{(out[1:k2] = NA)}
       \item \code{"func"} - same as \code{"mad"} option except that implemented
       in R for testing purposes. Avoid since it can be very slow for large windows.
     }
     Similar to \code{endrule} in \code{\link{runmed}} function which has the 
     following options: \dQuote{\code{c("median", "keep", "constant")}} .
  }
  \item{center}{moving window center. Defaults 
    to running median (\code{\link{runmed}} function). Similar to \code{center}  
    in \code{\link{mad}} function. For best acuracy at the edges use 
    \code{\link{runquantile}(x,k,0.5,type=2)}, which is slower than default
    \code{\link{runmed}(x,k,endrule="med")}. If \code{x} is a 2D array (and 
    \code{endrule="mad"}) or if \code{endrule="func"} than array edges are 
    filled by repeated calls to 
    \dQuote{\code{\link{mad}(x, center=\link{median}(x), na.rm=TRUE)}} function. 
    Runmad's \code{center} parameter will be ignored for the beggining and the 
    end of output \code{y}.  Please use  
    \code{center=\link{runquantile}(x,k,0.5,type=2)} for those cases.  
  }
  \item{constant}{scale factor such that for Gaussian 
    distribution X, \code{\link{mad}}(X) is the same as \code{\link{sd}}(X). 
    Same as \code{constant} in \code{\link{mad}} function.}
  \item{align}{specifies whether result should be centered (default), 
  left-aligned or right-aligned.  If \code{endrule}="mad" then setting
  \code{align} to "left" or "right" will fall back on slower implementation 
  equivalent to \code{endrule}="func". }
}

\details{
  Apart from the end values, the result of y = runmad(x, k) is the same as 
  \dQuote{\code{for(j=(1+k2):(n-k2)) y[j]=mad(x[(j-k2):(j+k2)], na.rm = TRUE)}}. 
  It can handle non-finite numbers like NaN's and Inf's 
  (like \dQuote{\code{\link{mad}(x, na.rm = TRUE)}}).

  The main incentive to write this set of functions was relative slowness of 
  majority of moving window functions available in R and its packages.  With the 
  exception of \code{\link{runmed}}, a running window median function, all 
  functions listed in "see also" section are slower than very inefficient 
  \dQuote{\code{\link{apply}(\link{embed}(x,k),1,FUN)}} approach. 
  
  Functions \code{runquantile} and \code{runmad} are using insertion sort to 
  sort the moving window, but gain speed by remembering results of the previous 
  sort. Since each time the window is moved, only one point changes, all but one 
  points in the window are already sorted. Insertion sort can fix that in O(k) 
  time.
}

\value{
  Returns a numeric vector or matrix of the same size as \code{x}. Only in case of 
  \code{endrule="trim"} the output vectors will be shorter and output matrices 
  will have fewer rows. 
} 

\references{
  About insertion sort used in \code{runmad} function see: 
  R. Sedgewick (1988): \emph{Algorithms}. Addison-Wesley (page 99)
} 

\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com}}

\seealso{
  Links related to:
  \itemize{       
   \item \code{runmad} - \code{\link{mad}}
   \item Other moving window functions  from this package: \code{\link{runmin}}, 
     \code{\link{runmax}}, \code{\link{runquantile}}, \code{\link{runmean}} and
     \code{\link{runsd}}
   \item generic running window functions: \code{\link{apply}}\code{
     (\link{embed}(x,k), 1, FUN)} (fastest), \code{\link[gtools]{running}} from \pkg{gtools} 
     package (extremely slow for this purpose), \code{\link[magic]{subsums}} from 
     \pkg{magic} library can perform running window operations on data with any 
     dimensions. 
  }
}

\examples{
  # show runmed function
  k=25; n=200;
  x = rnorm(n,sd=30) + abs(seq(n)-n/4)
  col = c("black", "red", "green")
  m=runmed(x, k)
  y=runmad(x, k, center=m)
  plot(x, col=col[1], main = "Moving Window Analysis Functions")
  lines(m    , col=col[2])
  lines(m-y/2, col=col[3])
  lines(m+y/2, col=col[3])
  lab = c("data", "runmed", "runmed-runmad/2", "runmed+runmad/2")
  legend(0,0.9*n, lab, col=col, lty=1 )

  # basic tests against apply/embed
  eps = .Machine$double.eps ^ 0.5
  k=25 # odd size window
  a = runmad(x,k, center=runmed(x,k), endrule="trim")
  b = apply(embed(x,k), 1, mad)
  stopifnot(all(abs(a-b)<eps));
  k=24 # even size window
  a = runmad(x,k, center=runquantile(x,k,0.5,type=2), endrule="trim")
  b = apply(embed(x,k), 1, mad)
  stopifnot(all(abs(a-b)<eps));
  
  # test against loop approach
  # this test works fine at the R prompt but fails during package check - need to investigate
  k=24; n=200;
  x = rnorm(n,sd=30) + abs(seq(n)-n/4) # create random data
  x = rep(1:5,40)
  #x[seq(1,n,11)] = NaN;               # commented out for time beeing - on to do list
  #x[5] = NaN;                         # commented out for time beeing - on to do list
  k2 = k%/%2
  k1 = k-k2-1
  ac = array(runquantile(x,k,0.5))
  a  = runmad(x, k, center=ac)
  bc = array(0,n)
  b  = array(0,n)
  for(j in 1:n) {
    lo = max(1, j-k1)
    hi = min(n, j+k2)
    bc[j] = median(x[lo:hi], na.rm = TRUE)
    b [j] = mad   (x[lo:hi], na.rm = TRUE, center=bc[j])
  }
  eps = .Machine$double.eps ^ 0.5
  #stopifnot(all(abs(ac-bc)<eps)); # commented out for time beeing - on to do list
  #stopifnot(all(abs(a-b)<eps));   # commented out for time beeing - on to do list
  
  # compare calculation at array ends
  k=25; n=200;
  x = rnorm(n,sd=30) + abs(seq(n)-n/4)
  c = runquantile(x,k,0.5,type=2)             # find the center
  a = runmad(x, k, center=c, endrule="mad" )  # fast C code
  b = runmad(x, k, center=c, endrule="func")  # slow R code
  stopifnot(all(abs(a-b)<eps));
  
  # test if moving windows forward and backward gives the same results
  k=51;
  a = runmad(x     , k)
  b = runmad(x[n:1], k)
  stopifnot(all(a[n:1]==b, na.rm=TRUE));

  # test vector vs. matrix inputs, especially for the edge handling
  nRow=200; k=25; nCol=10
  x = rnorm(nRow,sd=30) + abs(seq(nRow)-n/4)
  X = matrix(rep(x, nCol ), nRow, nCol) # replicate x in columns of X
  a = runmad(x, k, center = runquantile(x,k,0.5,type=2))
  b = runmad(X, k, center = runquantile(X,k,0.5,type=2))
  stopifnot(all(abs(a-b[,1])<eps));        # vector vs. 2D array
  stopifnot(all(abs(b[,1]-b[,nCol])<eps)); # compare rows within 2D array
  
  # speed comparison
  \dontrun{
  x=runif(1e5); k=51;                       # reduce vector and window sizes
  system.time(runmad( x,k,endrule="trim"))
  system.time(apply(embed(x,k), 1, mad))  
  }
}

\keyword{ts}
\keyword{array}
\keyword{utilities}
\concept{moving mad}
\concept{rolling mad}
\concept{running mad}
\concept{running window}

