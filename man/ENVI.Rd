\name{read.ENVI & write.ENVI}
\alias{read.ENVI}
\alias{write.ENVI}
\title{Read and Write Binary Data in ENVI Format}
\description{Read and write binary data in ENVI format, which is supported by 
  most GIS software.}
\usage{
  X=read.ENVI(filename, headerfile=paste(filename, ".hdr", sep="")) 
  write.ENVI (X, filename, interleave = c("bsq", "bil", "bip")) 
}

\arguments{
  \item{X}{data to be saved in ENVI file. Can be a matrix or 3D array.}
  \item{filename}{character string with name of the file (connection)}
  \item{headerfile}{optional character string with name of the header file}
  \item{interleave}{optional character string specifying interleave to be used}
}

\details{  
  ENVI binary files use a generalized raster data format that consists of two 
  parts: 
  \itemize{
    \item binary file - flat binary file equivalent to memory dump, as produced by 
    \code{\link{writeBin}} in R or \code{fwrite} in C/C++.
    \item header file - small text (ASCII) file containing the metadata 
     associated with the binary file. This file can contain the following 
     fields, followed by equal sign and a variable: 
    \itemize{
      \item \code{samples} - number of columns \cr
      \item \code{lines} - number of rows \cr
      \item \code{bands} - number of bands (channels, planes) \cr
      \item \code{data type} - following types are supported:
      \itemize{
        \item 1 - 1-byte unsigned integer
        \item 2 - 2-byte signed integer
        \item 3 - 4-byte signed integer
        \item 4 - 4-byte float
        \item 5 - 8-byte double
        \item 9 - 2x8-byte complex number made up from 2 doubles
        \item 12 - 2-byte unsigned integer
       }
      \item \code{header offset} -  number of bytes to skip before 
            raster data starts in binary file. 
      \item \code{interleave} - Permutations of dimensions in binary data:
      \itemize{
        \item \code{BSQ} - Band Sequential (X[col,row,band])
        \item \code{BIL} - Band Interleave by Line (X[col,band,row])
        \item \code{BIP} - Band Interleave by Pixel (X[band,col,row]) 
      }
      \item \code{byte order} - the endian-ness of the saved data: 
      \itemize{
        \item 0 - means little-endian byte order, format used on PC/Intel machines
        \item 1 - means big-endian (aka IEEE, aka "network") byte order, format 
          used on UNIX and Macintosh machines
      }
    }
  }
  Fields \code{samples}, \code{lines}, \code{bands}, \code{data type} are 
  required, while \code{header offset}, \code{interleave}, \code{byte order} are
  optional. All of them are in form of integers except \code{interleave} which
  is a string.
  
  This generic format allows reading of many raw file formats, including those 
  with embedded header information. Also it is a handy binary format to 
  exchange data between PC and UNIX/Mac machines, as well as different 
  languages like: C, Fortran, Matlab, etc. Especially since header files are 
  simple enough to edit by hand.
  
  File type supported by most of GIS (geographic information system) software
  including: ENVI software, Freelook (free file viewer by ENVI), ArcGIS, etc. 
}


\value{ Function \code{read.ENVI} returns either a matrix or 3D array. 
  Function \code{write.ENVI} does not return anything.} 


\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com}} 

\seealso{
  Displaying of images can be done through functions: \code{\link[graphics]{image}},
  \code{\link[fields]{image.plot}} and \code{\link[fields]{add.image}} from 
  \pkg{fields} or \code{\link[spatstat]{plot.im}} from \pkg{spatstat}.
  
  ENVI files are practically C-style memory-dumps as performed by 
  \code{\link{readBin}} and \code{\link{writeBin}} functions plus separate 
  meta-data header file.
  
  GIF file formats can also store 3D data (see \code{\link{read.gif}} and 
  \code{\link{write.gif}} functions).
  
  Packages related to GIS data: \pkg{shapefiles}, \pkg{maptools}, \pkg{sp}, 
  \pkg{spdep}, \pkg{adehabitat}, \pkg{GRASS}, \pkg{PBSmapping}. 
}

\examples{
  X = array(1:60, 3:5)
  write.ENVI(X, "temp.nvi")
  Y = read.ENVI("temp.nvi")
  stopifnot(X == Y)
  readLines("temp.nvi.hdr")
  
  d = c(20,30,40)
  X = array(runif(prod(d)), d)
  write.ENVI(X, "temp.nvi", interleave="bil")
  Y = read.ENVI("temp.nvi")
  stopifnot(X == Y)
  readLines("temp.nvi.hdr")
  
  file.remove("temp.nvi")
  file.remove("temp.nvi.hdr")
}

\keyword{file}
\concept{GIS data I/O}
