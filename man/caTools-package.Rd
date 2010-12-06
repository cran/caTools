\name{caTools-package}
\alias{caTools-package}
\alias{caTools}
\docType{package}
\title{
Tools: moving window statistics, GIF, Base64, ROC AUC, etc. 
}
\description{
Contains several basic utility functions including: moving
(rolling, running) window statistic functions, read/write for
GIF and ENVI binary files, fast calculation of AUC, LogitBoost
classifier, base64 encoder/decoder, round-off error free sum
and cumsum, etc.
}
\details{
\tabular{ll}{
Package: \tab caTools\cr
Version: \tab 1.6\cr
Date: \tab Apr 11 2006\cr
Depends: \tab R (>= 2.2.0), bitops\cr
Suggests: \tab MASS, rpart\cr
License: \tab The caMassClass Software License, Version 1.0 (See COPYING
file or \url{http://ncicb.nci.nih.gov/download/camassclasslicense.jsp})\cr
URL: \tab http://ncicb.nci.nih.gov/download/index.jsp\cr
Built: \tab R 2.2.1; i386-pc-mingw32; 2006-04-14 10:45:20; windows\cr
}

Index:
\preformatted{
LogitBoost              LogitBoost Classification Algorithm
predict.LogitBoost      Prediction Based on LogitBoost Algorithm
base64encode            Convert R vectors to/from the Base64 format
colAUC                  Column-wise Area Under ROC Curve (AUC)
combs                   All Combinations of k Elements from Vector v
read.ENVI               Read and Write Binary Data in ENVI Format
read.gif                Read and Write Images in GIF format
runmean                 Mean of a Moving Window
runmin                  Minimum and Maximum of Moving Windows
runquantile             Quantile of Moving Window
runmad                  Median Absolute Deviation of Moving Windows
runsd                   Standard Deviation of Moving Windows
sample.split            Split Data into Test and Train Set
sumexact                Basic Sum Operations without Round-off Errors
trapz                   Trapezoid Rule Numerical Integration
}

}
\author{Jarek Tuszynski <jaroslaw.w.tuszynski@saic.com>}
\keyword{ package }
\examples{
  # GIF image read & write
  write.gif( volcano, "volcano.gif", col=terrain.colors, flip=TRUE, 
           scale="always", comment="Maunga Whau Volcano")
  y = read.gif("volcano.gif", verbose=TRUE, flip=TRUE)
  image(y$image, col=y$col, main=y$comment, asp=1)

  # test runmin, runmax and runmed
  k=25; n=200;
  x = rnorm(n,sd=30) + abs(seq(n)-n/4)
  col = c("black", "red", "green", "brown", "blue", "magenta", "cyan")
  plot(x, col=col[1], main = "Moving Window Analysis Functions (window size=25)")
  lines(runmin (x,k), col=col[2])
  lines(runmed (x,k), col=col[3])
  lines(runmean(x,k), col=col[4])
  lines(runmax (x,k), col=col[5])
  legend(0,.9*n, c("data", "runmin", "runmed", "runmean", "runmax"), col=col, lty=1 )

  # sum vs. sumexact
  x = c(1, 1e20, 1e40, -1e40, -1e20, -1)
  a = sum(x);         print(a)
  b = sumexact(x);    print(b)

}
