#===========================================================================#
# caTools - R library                                                       #
# Copyright (C) 2005 Jarek Tuszynski                                        #
# Distributed under GNU General Public License version 3                    #
#===========================================================================#
#source('C:/programs/R/R-2.6.0/src/library/caTools/R/runfunc.R')

runmean = function(x, k, alg=c("C", "R", "fast", "exact"), 
                   endrule=c("mean", "NA", "trim", "keep", "constant", "func"))
{
  alg = match.arg(alg)
  endrule = match.arg(endrule)
  n = length(x)
  if (k<=1) return (x)
  if (k >n) k = n
  k2 = k%/%2
  y=double(n)
   
  if (alg=="exact") {
    .C("runmean_exact", x, y , as.integer(n), as.integer(k), 
       NAOK=TRUE, DUP=FALSE, PACKAGE="caTools") 
  } else if (alg=="C") {
    .C("runmean", as.double(x), y , as.integer(n), as.integer(k), 
       NAOK=TRUE, DUP=FALSE, PACKAGE="caTools") 
  } else if (alg=="fast") {
    .C("runmean_lite", as.double(x), y , as.integer(n), as.integer(k), 
       NAOK=TRUE, DUP=FALSE, PACKAGE="caTools") 
  } else {     # the similar algorithm implemented in R language
    k1 = k-k2-1
    y = c( sum(x[1:k]), diff(x,k) ); # find the first sum and the differences from it
    y = cumsum(y)/k                  # apply precomputed differences 
    y = c(rep(0,k1), y, rep(0,k2))   # make y the same length as x
    if (endrule=="mean") endrule="func"
  }
  if (endrule!="mean") y = EndRule(x, y, k, endrule, mean, na.rm=TRUE)
  return(y)
}

#==============================================================================

runmin = function(x, k, alg=c("C", "R"), 
                  endrule=c("min", "NA", "trim", "keep", "constant", "func"))
{
  alg = match.arg(alg)
  endrule = match.arg(endrule)
  n = length(x)
  if (k<=1) return (x)
  if (k >n) k = n
  y=double(n)
  
  if (alg=="C") {
    .C("runmin", as.double(x) ,y , as.integer(n), as.integer(k), 
       NAOK=TRUE, DUP=FALSE, PACKAGE="caTools")
  } else { # the similar algorithm implemented in R language
    k2 = k%/%2
    k1 = k-k2-1
    a <- y[k1+1] <- min(x[1:k], na.rm=TRUE)
    if (k!=n) for (i in (2+k1):(n-k2)) {
      if (a==y[i-1]) # point leaving the window was the min, so ...
        y[i] = min(x[(i-k1):(i+k2)], na.rm=TRUE) # recalculate min of the window 
      else           # min=y[i-1] is still inside the window
        y[i] = min(y[i-1], x[i+k2 ], na.rm=TRUE) # compare it with the new point 
      a = x[i-k1]    # point that will be removed from the window next
      if (!is.finite(a)) a=y[i-1]+1 # this will force the 'else' option
    }
    if (endrule=="min") endrule="func"
  }
  if (endrule!="min") y = EndRule(x, y, k, endrule, min, na.rm=TRUE)
  return(y)
}

#==============================================================================

runmax = function(x, k, alg=c("C", "R"), 
                  endrule=c("max", "NA", "trim", "keep", "constant", "func"))
{
  alg = match.arg(alg)
  endrule = match.arg(endrule)
  n = length(x)
  k = as.integer(k)
  if (k<=1) return (x)
  if (k >n) k = n
  y=double(n)

  if (alg=="C") {
    .C("runmax", as.double(x) ,y , as.integer(n), as.integer(k), 
       NAOK=TRUE, DUP=FALSE, PACKAGE="caTools")
  } else { # the same algorithm implemented in R language
    k2 = k%/%2
    k1 = k-k2-1
    a <- y[k1+1] <- max(x[1:k], na.rm=TRUE)
    if (k!=n) for (i in (2+k1):(n-k2)) {
      if (a==y[i-1]) # point leaving the window was the max, so ...
        y[i] = max(x[(i-k1):(i+k2)], na.rm=TRUE) # recalculate max of the window 
      else           # max=y[i-1] is still inside the window
        y[i] = max(y[i-1], x[i+k2 ], na.rm=TRUE) # compare it with the new point 
      a = x[i-k1]    # point that will be removed from the window next
      if (!is.finite(a)) a=y[i-1]+1 # this will force the 'else' option
    }
    if (endrule=="max") endrule="func"
  } 
  if (endrule!="max") y = EndRule(x, y, k, endrule, max, na.rm=TRUE)
  return(y)
}

#==============================================================================

runquantile = function(x, k, probs, type=7,
                endrule=c("quantile", "NA", "trim", "keep", "constant", "func"))
{ ## see http://mathworld.wolfram.com/Quantile.html for very clear definition
  ## of different quantile types
  endrule = match.arg(endrule)
  n    = length(x)
  np   = length(probs) 
  k    = as.integer(k)
  type = as.integer(type)
  if (k<=1) return (rep(x,n,np))
  if (k >n) k = n
  if (is.na(type) || (type < 1 | type > 9)) 
    warning("'type' outside allowed range [1,9]; changing 'type' to ", type<-7)
  
  y=double(n*np)
  .C("runquantile", as.double(x) ,y , as.integer(n), as.integer(k), 
       as.double(probs), as.integer(np),as.integer(type), 
       NAOK=TRUE, DUP=FALSE, PACKAGE="caTools")
     
  dim(y) =  c(n,np) 
  if (endrule=="trim") {
    yy = double((n-k+1)*np)
    dim(yy) = c(n-k+1,np) 
    for (i in 1:np) # for each percentile
      yy[,i] = EndRule(x, y[,i], k, endrule, quantile, probs=probs[i], type=type, na.rm=TRUE)
    y=yy
  } else if (endrule!="quantile") {
    for (i in 1:np) # for each percentile
      y[,i] = EndRule(x, y[,i], k, endrule, quantile, probs=probs[i], type=type, na.rm=TRUE)
  }
  attr(y, "k") = k
  return(y)
}

#==============================================================================

runmad = function(x, k, center = runmed(x,k), constant = 1.4826,
                  endrule=c("mad", "NA", "trim", "keep", "constant", "func"))
{
  endrule = match.arg(endrule)
  n = length(x)
  if (k<3) stop("'k' must be larger than 2")
  if (k>n) k = n
  y = double(n)
  .C("runmad", as.double(x), as.double(center), y, as.integer(n), 
       as.integer(k), NAOK=TRUE, DUP=FALSE, PACKAGE="caTools")
  if (endrule!="mad") y = EndRule(x, y, k, endrule, mad, constant=1, na.rm=TRUE)
  return(constant*y)
}

#==============================================================================

runsd = function(x, k, center = runmean(x,k), 
                 endrule=c("sd", "NA", "trim", "keep", "constant", "func"))
{
    endrule = match.arg(endrule)
   n = length(x)
   if (k<3) stop("'k' must be larger than 2")
   if (k>n) k = n
   y = double(n)
   .C("runsd", as.double(x), as.double(center), y, as.integer(n), 
        as.integer(k), NAOK=TRUE, DUP=FALSE, PACKAGE="caTools")
   if (endrule!="sd") y = EndRule(x, y, k, endrule, sd, na.rm=TRUE)
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
  k1 = k-k2-1
  idx1 = 1:k1
  idx2 = (n-k2+1):n
  endrule = match.arg(endrule)
  if (endrule=="NA") {
    y[idx1] = NA
    y[idx2] = NA
  } else if (endrule=="keep") {
    y[idx1] = x[idx1]
    y[idx2] = x[idx2]
  } else if (endrule=="constant") {
    y[idx1] = y[k1+1]
    y[idx2] = y[n-k2]
  } else if (endrule=="trim") {
    y = y[(k1+1):(n-k2)]
  } else if (endrule=="func") {
    for (i in idx1) y[i] = Func(x[1:(i+k2)], ...)
    for (i in idx2) y[i] = Func(x[(i-k1):n], ...)
  }
  return(y)
}

