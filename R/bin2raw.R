#===========================================================================#
# Written by Jarek Tuszynski. Copyright 2001-2003 SAIC.                     #
# Software developed in conjunction with the National Cancer Institute, and #
# distributed under "caBIO Software License" included in "COPYING" file.    #
#===========================================================================#

bin2raw = function(x, size=NA, endian=.Platform$endian)
{
  if (typeof(x)=="raw") return(x)
  TypeList = c("logical", "integer", "double", "complex", "character", "raw")
  if (is.character(x) && length(x)==1) x = strsplit(x, NULL)[[1]]  # convert strings to arrays of characters
  if (!is.vector(x) || mode(x) == "list") 
     stop("bin2raw: can only write vector objects")
  if (!is.na(size)) nBits=size
  else nBits = switch(match(typeof(x), TypeList), 4, 4, 8, 16, 2, 1)
  if (is.double(x)) { # accelerator for most common case
    if (!is.finite(size)) size = 8;
    if (size!=8 && size!=4) 
      stop("bin2raw: doubles of size ",size," are unknown on this machine")
    swap = as.integer(endian!=.Platform$endian)
    r = .Call("real2raw", x ,as.integer(size), swap, PACKAGE="caTools")    
  } else {            # slower code for exotic cases
    n = length(x)
    fname = tempfile() 
    writeBin(x, fname, size=size, endian=endian)
    r = readBin(fname, "raw", n = n*nBits)
    file.remove(fname)
  }
  return (r)
}

#====================================================================

raw2bin = function(r, what, size=NA, signed = TRUE, endian=.Platform$endian)
{
  TypeList = c("logical", "integer", "double", "complex", "character", "raw", 
               "numeric", "int")
  if (!is.character(what) || length(what) != 1 || !(what %in% TypeList)) 
    what <- typeof(what)
  if (!is.vector(r) || mode(r) == "list") 
     stop("raw2bin: 'r' has to be vector of type 'raw'")
  if (what=="raw") return(r)
  if (!is.na(size)) nBits=size 
  else nBits = switch(match(typeof(x), TypeList), 4, 4, 8, 16, 2, 1, 8, 4) 
  n = length(r)
  if (n%%nBits) stop("raw2bin: number of elements in 'r' is not multiple of 'size'")
  if (what=="double") { # accelerator for most common case
    if (nBits!=8 && nBits!=4) 
      stop("raw2bin: doubles of size ",nBits," are unknown on this machine")
    swap = as.integer(endian!=.Platform$endian)
    x = .Call("raw2real", as.raw(r) ,as.integer(nBits), swap, PACKAGE="caTools") 
  } else {              # slower code for exotic cases
    fname = tempfile()      
    writeBin(as.raw(r), fname)
    x = readBin(fname, what, n = n%/%nBits, size=size, signed=signed, endian=endian)
    file.remove(fname)
  }
  if (what=="character")  x = paste(x, collapse = "") # convert arrays of characters to strings
  return (x)
}

