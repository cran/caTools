#===========================================================================#
# Written by Jarek Tuszynski. Copyright 2001-2003 SAIC.                     #
# Software developed in conjunction with the National Cancer Institute, and #
# distributed under "caBIO Software License" included in "COPYING" file.    #
#===========================================================================#

bin2raw = function(x, size=NA, endian=.Platform$endian) 
{
  warning("function bin2raw become obsolete due to changes in writeBin: ", 
  "use 'writeBin(x, raw(), ...)' instead.")
  writeBin(x, raw(), size=size, endian=endian)
}

raw2bin = function(r, what, size=NA, signed = TRUE, endian=.Platform$endian)
{
  warning("function raw2bin become obsolete due to changes in readBin: ", 
  "use 'readBin(r, what, n=length(r)%/%size, size=size, ...)' instead.")
  TypeList = c("logical", "integer", "double", "complex", "character", "raw", 
               "numeric", "int")
  if (!is.character(what) || length(what) != 1 || !(what %in% TypeList)) 
    what <- typeof(what)
  if (!is.vector(r) || mode(r) == "list") 
     stop("raw2bin: 'r' has to be vector of type 'raw'")
  if (what=="raw") return(r)
  if (!is.na(size)) nBits=size 
  else nBits = switch(match(what, TypeList), 4, 4, 8, 16, 2, 1, 8, 4) 
  n = length(r)
  if (n%%nBits) stop("raw2bin: number of elements in 'r' is not multiple of 'size'")
  x = readBin(r, what, n = n%/%nBits, size=nBits, signed=signed, endian=endian)
  if (what=="character")  x = paste(x, collapse = "") # convert arrays of characters to strings
  return (x)
}

