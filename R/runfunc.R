#===========================================================================#
# Written by Jarek Tuszynski. Copyright 2001-2003 SAIC.                     #
# Software developed in conjunction with the National Cancer Institute, and #
# distributed under "caBIO Software License" included in "COPYING" file.    #
#===========================================================================#

runmean = function(x, k, alg=c("C", "R", "exact"), 
                   endrule=c("NA", "trim", "keep", "constant", "func"))
{
  alg = match.arg(alg)
  n = length(x)
  k = as.integer(k)
  k2 = k%/%2
  if (k2<1) k2 = 1
  if (k >n) k2 = (n-1)%/%2
  if (k!=1+2*k2)  
    warning("'k' must be odd number bigger than 3 and smaller than 'length(x)'.",
    "Changing 'k' to ", k <- as.integer(1 + 2*k2))
  y=double(n)
  if (n==k) y[k2+1]=sum.exact(x)/n 
  else if (alg=="C" && is.loaded("runmean")) {
    .C("runmean", as.double(x) ,y , as.integer(n), as.integer(k), 
       NAOK=FALSE, DUP=FALSE, PACKAGE="caTools") 
  } else if (alg=="exact" && is.loaded("runsum")) {
    .C("runsum", as.double(x) ,y ,Size<-integer(n), as.integer(n), 
       as.integer(k), NAOK=TRUE, DUP=FALSE, PACKAGE="caTools") 
    y = y/Size
    if (!all(Size)) y[Size==0] = NA
  } else {     # the similar algorithm implemented in R language
    y = c( sum(x[1:k]), diff(x,k) ); # find the first sum and the differences from it
    y = cumsum(y)/k                  # apply precomputed differences 
    y = c(rep(0,k2), y, rep(0,k2))   # make y the same length as x
  }
  y = EndRule(x, y, k, endrule, mean)
  return(y)
}

#==============================================================================

runmin = function(x, k, alg=c("C", "R"), 
                  endrule=c("NA", "trim", "keep", "constant", "func"))
{
  alg = match.arg(alg)
  n = length(x)
  k = as.integer(k)
  k2 = k%/%2
  if (k2<1) k2 = 1
  if (k >n) k2 = (n-1)%/%2
  if (k!=1+2*k2)  
    warning("'k' must be odd number bigger than 3 and smaller than 'length(x)'.",
    "Changing 'k' to ", k <- as.integer(1 + 2*k2))
  y=double(n)
  if (n==k) y[k2+1]=min(x) else
  if (alg=="C" && is.loaded("runquantile")) {
    .C("runquantile", as.double(x) ,y , as.integer(n), as.integer(k), 
       as.double(1), as.integer(1), NAOK=FALSE, DUP=FALSE, PACKAGE="caTools")
  } else { # the same algorithm implemented in R language
    a <- y[k2+1] <- min(x[1:k])
    for (i in (2+k2):(n-k2)) {
      if (a==y[i-1]) y[i] = min(x[(i-k2):(i+k2)]) # calculate min of the window 
      else           y[i] = min(y[i-1], x[i+k2])  # min of the window is =y[i-1]
      a = x[i-k2] # point that will be removed from the window next
    }
  }
  y = EndRule(x, y, k, endrule, min)
  return(y)
}

#==============================================================================

runmax = function(x, k, alg=c("C", "R"),
                  endrule=c("NA", "trim", "keep", "constant", "func"))
{
  alg = match.arg(alg)
  n = length(x)
  k = as.integer(k)
  k2 = k%/%2
  if (k2<1) k2 = 1
  if (k >n) k2 = (n-1)%/%2
  if (k!=1+2*k2)  
    warning("'k' must be odd number bigger than 3 and smaller than 'length(x)'.",
    "Changing 'k' to ", k <- as.integer(1 + 2*k2))
  y=double(n)
  if (n==k) y[k2+1]=max(x) else
  if (alg=="C" && is.loaded("runquantile")) {
    .C("runquantile", as.double(x) ,y , as.integer(n), as.integer(k), 
       as.double(k), as.integer(1), NAOK=FALSE, DUP=FALSE, PACKAGE="caTools")
  } else { # the same algorithm implemented in R language
    k2 = k%/%2
    a <- y[k2+1] <- max(x[1:k])
    for (i in (2+k2):(n-k2)) {
      if (a==y[i-1]) y[i] = max(x[(i-k2):(i+k2)]) # calculate max of the window 
      else           y[i] = max(y[i-1], x[i+k2])  # max of the window is =y[i-1] 
      a = x[i-k2] # point that will be removed from the window next
    }
  }
  y = EndRule(x, y, k, endrule, max)
  return(y)
}

#==============================================================================

runquantile = function(x, k, probs, type=7, 
                       endrule=c("NA", "trim", "keep", "constant", "func"))
{ ## see http://mathworld.wolfram.com/Quantile.html for very clear definition
  ## of different quatile types
  n    = length(x)
  np   = length(probs) 
  k    = as.integer(k)
  type = as.integer(type)
  k2 = k%/%2
  if (k2<1) k2 = 1
  if (k >n) k2 = (n-1)%/%2
  if (k!=1+2*k2)  
    warning("'k' must be odd number bigger than 3 and smaller than 'length(x)'.",
    "Changing 'k' to ", k <- as.integer(1 + 2*k2))
  if (is.na(type) || (type < 1 | type > 9)) 
    warning("'type' outside allowed range [1,9]; changing 'type' to ", type<-7)
  
  if (np==1 && 2*probs==1) { # special case - return runmed it is faster
    erule = endrule
    if (endrule=="func") erule="median"
    if (endrule=="NA" || endrule=="trim") erule="keep"
    y = runmed(x, k, endrule=erule)
    if (endrule=="NA" || endrule=="trim") 
      y = EndRule(x, y, k, endrule, quantile, probs=probs[i], type=type)
    dim(y) =  c(n,1) 
    return(y)
  }
  # the following code is based on code from quantile.default function
  if (type <= 3) {    ## Types 1, 2 and 3 are discontinuous sample qs. 
    if (type == 3) nppm = k * probs - .5 # n * probs + m; m = -0.5 
    else           nppm = k * probs      # m = 0 
    j = floor(nppm) 
    switch(type, 
      h = ifelse(nppm > j, 1, 0),                   # type 1 
      h = ifelse(nppm > j, 1, 0.5),                 # type 2 
      h = ifelse((nppm==j) && ((j%%2) == 0), 0, 1)) # type 3 
  } else {            ## Types 4 through 9 are continuous sample qs. 
    switch(type - 3, 
     {a<-0; b<-1},  # type 4 
      a<-b<-0.5,    # type 5 
      a<-b<-0,      # type 6 
      a<-b<-1,      # type 7 
      a<-b<-1/3,    # type 8 
      a<-b<-3/8)    # type 9 
    nppm = a + probs * (k + 1 - a - b) # n*probs + m 
    fuzz = 4 * .Machine$double.eps
    j = floor(nppm + fuzz)
    h = nppm - j
    h = ifelse(abs(h) < fuzz, 0, h)
  } 
  nppm = j+h
  nppm = ifelse(nppm<1, 1, nppm)
  nppm = ifelse(nppm>k, k, nppm)
     
  y=double(n*np)
  .C("runquantile", as.double(x) ,y , as.integer(n), as.integer(k), 
     as.double(nppm), as.integer(np), NAOK=FALSE, DUP=FALSE, PACKAGE="caTools")
  dim(y) =  c(n,np) 
  if (endrule=="trim") y = y[(k2+1):(n-k2),]
  else {
    for (i in 1:np) 
      y[,i] = EndRule(x, y[,i], k, endrule, quantile, probs=probs[i], type=type)
  }
  attr(y, "k") = k
  return(y)
}

#==============================================================================

runmad = function(x, k, center = runmed(x,k,endrule="keep"), constant = 1.4826, 
                  endrule=c("NA", "trim", "keep", "constant", "func"))
{
  n = length(x)
  k = as.integer(k)
  constant = as.double(constant)
  k2 = k%/%2
  if (k2<1) k2 = 1
  if (k >n) k2 = (n-1)%/%2
  if (k!=1+2*k2)  
    warning("'k' must be odd number bigger than 3 and smaller than 'length(x)'.",
    "Changing 'k' to ", k <- as.integer(1 + 2*k2))
  y = double(n)
  if (n==k) y[k2+1]=mad(x, constant=1) else
  .C("runmad", as.double(x), as.double(center), y, as.integer(n), 
     as.integer(k), NAOK=FALSE, DUP=FALSE, PACKAGE="caTools")
  y = constant*EndRule(x, y, k, endrule, mad, constant)
  return(y)
}

#==============================================================================

EndRule = function(x, y, k, 
             endrule=c("NA", "trim", "keep", "constant", "func"), Func, ...)
{
  n = length(x)
  if (length(y)!=n) stop("vectors 'x' and 'y' have to have the same length.")
  k = as.integer(k)
  k2 = k%/%2
  if (k2<1) k2 = 1
  if (k >n) k2 = (n-1)%/%2
  if (k!=1+2*k2)  
    warning("'k' must be odd number bigger than 3 and smaller than 'length(x)'.",
    "Changing 'k' to ", k <- as.integer(1 + 2*k2))
  idx1 = 1:k2
  idx2 = (n-k2+1):n
  endrule = match.arg(endrule)
  if (endrule=="NA") {
    y[idx1] = NA
    y[idx2] = NA
  } else if (endrule=="keep") {
    y[idx1] = x[idx1]
    y[idx2] = x[idx2]
  } else if (endrule=="constant") {
    y[idx1] = y[k2+1]
    y[idx2] = y[n-k2]
  } else if (endrule=="trim") {
    y = y[(k2+1):(n-k2)]
  } else if (endrule=="func") {
    for (i in idx1) y[i] = Func(x[1:i], ...)
    for (i in idx2) y[i] = Func(x[i:n], ...)
  }
  attr(y, "k") = k
  return(y)
}

