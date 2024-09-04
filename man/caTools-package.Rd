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

\author{Jarek Tuszynski <jaroslaw.w.tuszynski@saic.com>}
\keyword{package}
\examples{
  # GIF image read & write
  write.gif( volcano, "volcano.gif", col=terrain.colors, flip=TRUE, 
           scale="always", comment="Maunga Whau Volcano")
  y = read.gif("volcano.gif", verbose=TRUE, flip=TRUE)
  image(y$image, col=y$col, main=y$comment, asp=1)
  file.remove("volcano.gif")

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
