\name{base64encode & base64decode}
\alias{base64encode}
\alias{base64decode}
\title{Convert R vectors to/from the Base64 format }
\description{
  Convert R vectors of any type to and from the Base64 format for encrypting
  any binary data as string using alphanumeric subset of ASCII character set.
}
\usage{
  z = base64encode(x, size=NA, endian=.Platform$endian)
  x = base64decode(z, what, size=NA, signed = TRUE, endian=.Platform$endian)
}

\arguments{
  \item{x}{vector or any structure that can be converted to a vector by  
    \code{\link{as.vector}} function. Strings are also allowed.}
  \item{z}{String with Base64 code, using [A-Z,a-z,0-9,+,/,=] subset of 
    characters}
  \item{what}{Either an object whose mode will give the mode of the vector
          to be created, or a character vector of length one describing
          the mode: one of '"numeric", "double", "integer", "int",
          "logical", "complex", "character", "raw".  
          Same as variable \code{what} in \code{\link{readBin}} functions. }
  \item{size}{ integer.  The number of bytes per element in the byte stream 
    stored in \code{r}. The default, '\code{NA}', uses the natural size. 
    Same as variable \code{size} in \code{\link{readBin}} functions. }
  \item{signed}{logical. Only used for integers of sizes 1 and 2, when it
          determines if the quantity stored as raw should be regarded as a
          signed or unsigned integer.
          Same as variable \code{signed} in \code{\link{readBin}} functions. }
   \item{endian}{If provided, can be used to swap endian-ness. Using '"swap"' 
     will force swapping of byte order. Use '"big"' (big-endian, aka IEEE, 
     aka "network") or '"little"' (little-endian, format used on PC/Intel 
     machines) to indicate type of data encoded in "raw" format.
     Same as variable \code{endian} in \code{\link{readBin}} functions.}
}

\details{
The Base64 encoding is designed to encode arbitrary binary information for 
transmission by electronic mail. It is defined by MIME (Multipurpose Internet 
Mail Extensions) specification RFC 1341, RFC 1421, RFC 2045 and others. 
Triplets of 8-bit octets are encoded as groups of four characters, each 
representing 6 bits of the source 24 bits. Only a 65-character subset 
([A-Z,a-z,0-9,+,/,=]) present in all variants of ASCII and EBCDIC is used, 
enabling 6 bits to be represented per printable character.

Default \code{size}s for different types of \code{what}: \code{logical} - 4, 
 \code{integer} - 4, \code{double} - 8 , \code{complex} - 16, 
 \code{character} - 2, \code{raw} - 1.

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
  \code{\link[XML]{xmlValue}} from \pkg{XML} package reads XML code
   which sometimes is encoded in Base64 format.
   
  \code{\link{readBin}}, \code{\link{writeBin}}
}

\examples{
   x = (10*runif(10)>5) # logical
   for (i in c(NA, 1, 2, 4)) {
     y = base64encode(x, size=i)
     z = base64decode(y, typeof(x), size=i)
     stopifnot(x==z)
   }
   print("Checked base64 for encode/decode logical type")

   x = as.integer(1:10) # integer
   for (i in c(NA, 1, 2, 4)) {
     y = base64encode(x, size=i)
     z = base64decode(y, typeof(x), size=i)
     stopifnot(x==z)
   }
   print("Checked base64 encode/decode for integer type")
   
   x = (1:10)*pi        # double
   for (i in c(NA, 4, 8)) {
     y = base64encode(x, size=i)
     z = base64decode(y, typeof(x), size=i)
     stopifnot(mean(abs(x-z))<1e-5)
   }
   print("Checked base64 for encode/decode double type")
   
   x = log(as.complex(-(1:10)*pi))        # complex
   y = base64encode(x)
   z = base64decode(y, typeof(x))
   stopifnot(x==z)
   print("Checked base64 for encode/decode complex type")
  
   x = "Chance favors the prepared mind" # character
   y = base64encode(x)
   z = base64decode(y, typeof(x))
   stopifnot(x==z)
   print("Checked base64 for encode/decode character type")
}
\keyword{file}
\concept{XML}
