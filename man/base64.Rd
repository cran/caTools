\name{base64}
\alias{base64encode}
\alias{base64decode}
\title{Convert R vectors to/from the Base64 format }
\description{
  Convert R vectors of any type to and from the Base64 format for encrypting
  any binary data as string using alphanumeric subset of ASCII character set.
}
\usage{
  z = base64encode(x, \dots)
  x = base64decode(z, what, \dots)
}

\arguments{
  \item{x}{vector or any structure that can be converted to a vector by  
    \code{\link{as.vector}} function. Strings are also allowed.}
  \item{z}{String with Base64 code, using [A-Z,a-z,0-9,+,/,=] subset of 
    characters}
  \item{what}{Either an object whose mode will give the mode of the vector
          to be created, or a character vector of length one describing
          the mode: one of '"numeric", "double", "integer", "int",
          "logical", "complex", "character", "raw". Same as variable 
          \code{what} in \code{\link{readBin}} and \code{\link{base64decode}}
          functions. }
  \item{\ldots}{parameters to be passed to \code{\link{bin2raw}} and 
    \code{\link{raw2bin}} functions.}
}

\details{
The Base64 encoding is designed to encode arbitrary binary information for 
transmission by electronic mail. It is defined by MIME (Multipurpose Internet 
Mail Extensions) specification RFC 1341, RFC 1421, RFC 2045 and others. 
Triplets of 8-bit octets are encoded as groups of four characters, each 
representing 6 bits of the source 24 bits. Only a 65-character subset 
([A-Z,a-z,0-9,+,/,=]) present in all variants of ASCII and EBCDIC is used, 
enabling 6 bits to be represented per printable character
}

\value{
   Function \code{\link{base64encode}} returns a string with Base64 code.
   Function \code{\link{base64decode}} returns vector of appropriate mode 
    and length (see \code{x} above).
}

\references{
 \itemize{
   \item Base64 description in \emph{Connected: An Internet Encyclopedia}
  \url{http://www.freesoft.org/CIE/RFC/1521/7.htm}
  \item MIME RFC 1341 \url{http://www.faqs.org/rfcs/rfc1341.html}
  \item MIME RFC 1421 \url{http://www.faqs.org/rfcs/rfc1421.html}
  \item MIME RFC 2045 \url{http://www.faqs.org/rfcs/rfc2045.html}
  \item Portions of the code are based on Matlab code by Peter Acklam
   \url{http://home.online.no/~pjacklam/matlab/software/util/datautil/}
 }
}

\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com}} 

\seealso{ 
  \itemize{
  \item \code{\link{bin2raw}} and \code{\link{raw2bin}} are being used to
   convert R vectors to and from the raw binary format.
   \item \code{\link[XML]{xmlValue}} from \pkg{XML} package often reads XML code
   which sometimes is encoded in Base64 format.
  }
}

\examples{
   x = (10*runif(10)>5) # logical
   for (i in c(NA, 1, 2, 4)) {
     y = base64encode(x, size=i)
     z = base64decode(y,typeof(x), size=i)
     stopifnot(x==z)
   }
   print("Checked base64 for encode/decode logical type")

   x = as.integer(1:10) # integer
   for (i in c(NA, 1, 2, 4)) {
     y = base64encode(x, size=i)
     z = base64decode(y,typeof(x), size=i)
     stopifnot(x==z)
   }
   print("Checked base64 encode/decode for integer type")
   
   x = (1:10)*pi        # double
   for (i in c(NA, 4, 8)) {
     y = base64encode(x, size=i)
     z = base64decode(y,typeof(x), size=i)
     stopifnot(mean(abs(x-z))<1e-5)
   }
   print("Checked base64 for encode/decode double type")
   
   x = log(as.complex(-(1:10)*pi))        # complex
   y = base64encode(x)
   z = base64decode(y,typeof(x))
   stopifnot(x==z)
   print("Checked base64 for encode/decode complex type")
  
   x = "Chance favors the prepared mind" # character
   y = base64encode(x)
   z = base64decode(y,typeof(x))
   stopifnot(x==z)
   print("Checked base64 for encode/decode character type")
}
\keyword{file}

