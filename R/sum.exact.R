#===========================================================================#
# caTools - R library                                                       #
# Copyright (C) 2005 Jarek Tuszynski                                        #
# Distributed under GNU General Public License version 3                    #
#===========================================================================#

sum.exact = function(..., na.rm = FALSE)
{
  x = c(...,  recursive=TRUE)
  if (na.rm) x = x[!is.na(x)]
  else if (any(is.na(x))) return(NA)
  n = length(x)
  .C("sum_exact", as.double(x), y<-as.double(0), as.integer(n),  
      NAOK=TRUE, DUP=FALSE, PACKAGE="caTools") 
  return(y)
}

#==============================================================================

cumsum.exact = function(x)
{
  n = length(x)
  .C("cumsum_exact", as.double(x), y<-double(n), as.integer(n),  
      NAOK=TRUE, DUP=FALSE, PACKAGE="caTools") 
  return(y)
}

#==============================================================================

runsum.exact = function(x, k)
{
  n = length(x)
  k = as.integer(k)
  k2 = k%/%2
  if (k2<1) k2 = 1
  if (k >n) k2 = (n-1)%/%2
  if (k!=1+2*k2)  
    warning("'k' must be odd number bigger than 3 and smaller than 'length(x)'.",
    "Changing 'k' to ", k <- as.integer(1 + 2*k2))
  .C("runsum", as.double(x) ,y<-double(n) ,Size<-integer(n), as.integer(n),  
     as.integer(k), NAOK=TRUE, DUP=FALSE, PACKAGE="caTools") 
  idx = (k2+1):(n-k2)
  y = y[idx]
  Size = Size[idx]
  attr(y, 'count') <- Size
  return(y)
}

