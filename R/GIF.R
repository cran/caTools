#===========================================================================#
# Written by Jarek Tuszynski. (2005)                                        #
# Software was developed in my private time, but it is distributed under    #
# "caBIO Software License", like the rest of the caTools package.           #
#===========================================================================#

write.gif = function(image, filename, col=NULL, scale=c("smart", "never", "always"), 
    transparent=NULL, comment=NULL, delay=0,  flip=FALSE,  interlace=FALSE)
{
  if (!is.character(filename)) stop("write.gif: 'filename' has to be a string")
  if (length(filename)>1) filename = paste(filename, collapse = "")  # combine characters into a string

  #=================================
  # cast x into a proper dimentions
  #=================================
  dm = dim(image)
  if (is.null(dm)) stop("write.gif: input 'x' has to be an matrix or 3D array")
  if (length(dm)<=2) { # this is a 2D matrix or smaller
    image = as.matrix(image)   # cast to 2D matrix
    if (flip) x = image[,dm[2]:1]
    else x = t(image)
  } else {             # 3D data cube or bigger
    dim(image) = c(dm[1], dm[2], prod(dm)/(dm[1]*dm[2])) # cast to 3D
    if (flip) x = image[,dm[2]:1,]
    else x = aperm(image, c(2,1,3))
  }
  image = 0            # release memory
  dm = dim(x)          # save dimentions and ...
  x = as.vector(x)     # convert to 1D vector
  
  #=================================
  # scale x into a proper range
  #=================================
  scale = match.arg(scale)
  if (!is.null(transparent)) 
   if ((transparent<0) || (transparent>255)) 
    stop("write.gif:'transparent' has to be an integer between 0 and 255")
  mask = !is.finite(x)
  xx = 0
  mColor = 255
  if (any(mask)) {  # some non-finite numbers were found
    if (is.null(transparent)) mColor = 254
    xx = x          # save original x
    x  = x[!mask]   # remove non-finite numbers
  }
  minx = min(x)
  maxx = max(x)
  d = mColor/(maxx-minx)
  if (scale=="never") {
    if ((minx<0) || (maxx>mColor)) 
     warning("write.gif: 'x' is not in proper range and 'scale' is set to 'never',",
     " clipping 'x' to proper range ")
    if (minx<0     ) x[x<0     ] = 0 
    if (maxx>mColor) x[x>mColor] = mColor 
  } else
  if (scale=="always") {
    if ((minx>=0) && (maxx<=1)) 
      x  = mColor*x    # doubles between [0 and 1] -> scale them
    else 
      x = (x-minx)*d   # numbers outside allowed range -> scale them
  } else
  if (scale=="smart") {
    if ((minx<0) || (maxx>mColor)) {
      x = (x-minx)*d   # numbers outside allowed range -> scale them
    } else if ((minx>=0) && (maxx<=1)) {
      if (any(x!=as.integer(x))) x = mColor*x    # doubles between [0 and 1] -> scale them
    }
  }
  maxx = max(x)

  if (length(xx)>1) { # some non-finite numbers were found
    if (is.null(transparent)) transparent = maxx+1
    xx[ mask] = transparent
    xx[!mask] = x
    x = xx
  }
  if (is.null(transparent)) transparent = -1
  x = as.integer(round(x))
  
  #=================================
  # format color palette
  #=================================
  if (is.null(col)) { # no color palette we will use gray scale
    a =  as.integer(seq(0,255,length.out=maxx+1))
    crgb = rbind(a,a,a)
  } else crgb  = col2rgb(col)
  Palette = as.integer(c(256^(2:0) %*% crgb)) # convert to internal int format
  nColor = length(Palette)
  if (nColor<maxx) 
    stop("write.gif: not enough colors in color palette 'col'. Has ",nColor,
         " need at least ", maxx)
  if (nColor<256) Palette = c(Palette, rep(0,256-nColor)) # pad it
  
  # format and cast other input variables into proper format
  param = as.integer(c( dm[2], dm[1], prod(dm)/(dm[1]*dm[2]), nColor, transparent, delay, interlace, 0 ))
  if (is.null(comment)) comment = as.character("")
  else comment = as.character(comment)
  # call C++ function
  .C("imwritegif", filename, x, Palette, param, comment,
     NAOK=FALSE, PACKAGE="caTools") 
  if (param[7]<0) stop("write.gif: cannot open the output file (connection)")
  invisible(NULL)
}

#==============================================================================

read.gif = function(filename, frame=0, flip=FALSE, verbose=FALSE)
{
  if (!is.character(filename)) stop("write.gif: 'filename' has to be a string")
  if (length(filename)>1) filename = paste(filename, collapse = "")  # combine characters into a string
  isURL = length(grep("^http://", filename)) | 
          length(grep("^ftp://",  filename)) | 
          length(grep("^file://", filename))
  if(isURL) {
    tf <- tempfile()
    download.file(filename, tf, mode='wb', quiet=TRUE)
    filename = tf
  }

  x = .Call("imreadgif", filename, as.integer(frame), as.integer(verbose), 
       PACKAGE="caTools") 
  comt = as.character(attr(x, 'comm'))
  if (isURL) file.remove(filename)

  nRow    = x[1]
  nCol    = x[2]
  nBand   = x[3]
  tran    = x[4]
  success = x[5]
  nPixel  = nRow*nCol*nBand
  stats = -success
  if (stats>=6)  {
    warning("write.gif: file '", filename, 
      "' contains multiple color-maps. Use 'frame' > 0.") 
    stats = stats-6
  }
  if (nPixel==0) {
    switch (stats,
    stop("write.gif: cannot open the input file: ", filename),
    stop("write.gif: input file '", filename, "' is not a GIF file"),
    stop("write.gif: unexpected end of file: ", filename),
    stop("write.gif: syntax error in file: ", filename) )
  } else {
    switch (stats, , , 
    warning("write.gif: unexpected end of file: ", filename),
    warning("write.gif: syntax error in file: ", filename),
    warning("write.gif: file '", filename,
      "' contains multiple images (frames) of uneven length. Use 'frame' > 0." ))
  }   
  Palette = x[ 10:265 ]
  x       = x[-(1:265)] # delete non image data
  if (nBand>1) { # 3D data cubes
    dim(x)  = c(nCol, nRow, nBand)
    if (flip) x = x[,ncol(x):1,]
    else x = aperm(x, c(2,1,3))
  } else {       # this is a matrix
    dim(x) = c(nCol, nRow)
    if (flip) x = x[,ncol(x):1]
    else x = t(x)
  }
  Palette = Palette[Palette>=0]
  red     = bitAnd(bitShiftR(Palette,16), 255)
  green   = bitAnd(bitShiftR(Palette, 8), 255)
  blue    = bitAnd(          Palette    , 255)
  Palette = rgb (red, green, blue, maxColorValue = 255)
  if (tran==-1) tran = NULL
  return (list(image=x, col=Palette, transparent=tran, comment=comt))
}

# source("c:/programs/R/rw2011/src/library/caTools/R/GIF.R")
