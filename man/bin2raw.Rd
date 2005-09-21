\name{bin2raw & raw2bin}
\alias{bin2raw}
\alias{raw2bin}
\title{Convert R vectors to/from the raw binary format}
\description{
  Convert R vectors of any type to and from the raw binary format, stored as 
  vector of type "raw".
}

\usage{
  r = bin2raw(x, size=NA, endian=.Platform$endian)
  x = raw2bin(r, what, size=NA, signed = TRUE, endian=.Platform$endian)
}

\arguments{
  \item{x}{vector or any structure that can be converted to a vector by  
    \code{\link{as.vector}} function. Strings are also allowed.}
  \item{r}{vector of type "raw"}
  \item{what}{Either an object whose mode will give the mode of the vector
          to be created, or a character vector of length one describing
          the mode: one of '"numeric", "double", "integer", "int",
          "logical", "complex", "character", "raw". Same as variable 
          \code{what} in \code{\link{readBin}} and \code{\link{base64decode}}
          functions. }
  \item{size}{ integer.  The number of bytes per element in the byte stream 
    stored in \code{r}. The default, '\code{NA}', uses the natural size. See 
    details. }
  \item{signed}{logical. Only used for integers of sizes 1 and 2, when it
          determines if the quantity stored as raw should be regarded as a
          signed or unsigned integer.}
   \item{endian}{If provided, can be used to swap endian-ness. Using '"swap"' 
     will force swapping of byte order. Use '"big"' (big-endian, aka IEEE, 
     aka "network") or '"little"' (little-endian, format used on PC/Intel 
     machines) to indicate type of data encoded in "raw" format.}
}

\details{
  Quoting from \code{\link{readBin}} documentation:
  
   "If 'size' is specified and not the natural size of the object,
   each element of the vector is coerced to an appropriate type
   before being written or as it is read.  Possible sizes are 1, 2, 4
   and possibly 8 for integer or logical vectors, and 4, 8 and
   possibly 12/16 for numeric vectors.  (Note that coercion occurs as
   signed types except if 'signed = FALSE' when reading integers of
   sizes 1 and 2.) Changing sizes is unlikely to preserve 'NA's, and
   the extended precision sizes are unlikely to be portable across
   platforms."
   
   Default \code{size}s for different types: \code{logical} - 4, 
   \code{integer} - 4, \code{double} - 8 , \code{complex} - 16, 
   \code{character} - 2, \code{raw} - 1.
}

\value{
  Function \code{bin2raw} returns vector of raw values (see \code{r} 
  above), where each 1-byte raw value correspond to 1-byte of 1-byte of the 
  binary form of other types. Length of the vector is going to be "number of 
  bytes of a  single element in array \code{x}" times \code{length(x)}.
  
  Function \code{raw2bin} returns vector of appropriate mode 
  and length (see \code{x} above), where each 1-byte raw value correspond to 
  1-byte of the binary form of other types. Length of the vector is going to 
  be number of bytes per element in array \code{x} times \code{length(x)}. 
  If parameter \code{what} is equal to "character" than a string (of length 1)
  is returned instead of vector f characters.
}

\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com}} 
\note{ 
  At the moment those two functions use calls to \code{\link{writeBin}} 
  and \code{\link{readBin}} functions to do the job. I hope to change it in 
  the future by writing a C code. 
}

\seealso{  
  \code{\link{readBin}}, \code{\link{writeBin}}
}

\examples{
   print(x <- (1:5)*pi)      
   print(y <- bin2raw(x))
   print(z <- raw2bin(y,"double"))
  
   x = (10*runif(10)>5) # logical
   for (i in c(NA, 1, 2, 4)) {
     y = bin2raw(x, size=i)
     z = raw2bin(y,typeof(x), size=i)
     stopifnot(x==z)
   }
   print("Checked bin2raw and raw2bin conversion for logical type")
   
   x = as.integer(1:10) # integer
   for (i in c(NA, 1, 2, 4)) {
     y = bin2raw(x, size=i)
     z = raw2bin(y,typeof(x), size=i)
     stopifnot(x==z)
   }
   print("Checked bin2raw and raw2bin conversion for integer type")
  
   x = (1:10)*pi        # double
   for (i in c(NA, 4, 8)) {
     y = bin2raw(x, size=i)
     z = raw2bin(y,typeof(x), size=i)
     stopifnot(mean(abs(x-z))<1e-5)
   }
   print("Checked bin2raw and raw2bin conversion for double type")
   
   x = log(as.complex(-(1:10)*pi))        # complex
   y = bin2raw(x)
   z = raw2bin(y,typeof(x))
   stopifnot(x==z)
   print("Checked bin2raw and raw2bin conversion for complex type")
   
   x = "Chance favors the prepared mind" # character
   y = bin2raw(x)
   z = raw2bin(y,typeof(x))
   stopifnot(x==z)
   print("Checked bin2raw and raw2bin conversion for character type")
}
\keyword{file}

