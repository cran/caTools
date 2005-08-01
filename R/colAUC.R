#===========================================================================#
# Written by Jarek Tuszynski. Copyright 2001-2003 SAIC.                     #
# Software developed in conjunction with the National Cancer Institute, and #
# distributed under "caBIO Software License" included in "COPYING" file.    #
#===========================================================================#

colAUC = function (X, y, p.val=FALSE)
{   
  # given:
  #   X - 2D matrix (of feature columns and samples rows)
  #   y - 1D array identifying which class each sample belongs to (csource(lass numbers start from 1)
  # find AUC of each feature
  nR   = nrow(X) 
  nC   = ncol(X)                        # get dimentions of the data set
  uL   = sort(unique(y))                # find all the classes among the labels and or
  nL   = length(uL)                     # number of unique classes
  if (nL<=1) 
    stop("colAUC: List of labels y have to contain at least 2 class labels.")
  L    = matrix(rep(uL,each=nR),nR,nL)  # store vector L as row vector and copy it into nR rows
  per  = combs(1:nL,2)                  # find all possible pairs of L columns
  nP   = nrow(per)                      # how many possible pairs were found?
  Auc  = matrix(0.5,nP,nC)              # Initialize array to store results
  for (j in 1:nC) {                     # for each column representing a feature
    x = sort(X [, j], index=TRUE)       # sort all columns and store each one in x[[1]]. x[[2]] stores original positions
    nunq = which(diff(x$x)==0)          # find non-unique A's in column j (if vector is [1 1] nunq=1
    n = length(nunq)                    # number of non-unique values
    if (n<nR-1) {                       # make sure all numbers in A column are not the same
      idx = y[x$ix]                     # reorder label vector in the same order as b, or associate label with each number in b
      # assign column for each label (class) and for each point add 1 in the column corresponding to its class
      d = ( matrix(rep(idx,nL),nR,nL) == L ) 
      for (i in 1:nL) d[,i] = cumsum(d[,i])  # cumulative sum of d columns
      if (n) d = d[-nunq, ]             # remove non unique rows if any
      d = rbind( matrix(0,1,nL), d )    # append row of zeros at the beggining
      nD = nrow(d)
      # assume that col#1 ploted on x axis is correct clasification and col#2 (y) is false find AUC
      for (i in 1:nP) {                 # go through all permutations of columns in d
        c1 = per[i,1]                   # and identify 2 classes to be compared
        c2 = per[i,2]
        n  = d[nD,c1]*d[nD,c2]          # normalize area to 1 at the maximum
        if (n>0) Auc[i,j] = trapz(d[,c1], d[,c2])/n  # Trapezoidal numerical integration
       # print(sprintf('%i %f %f ', j, Auc[i,j], n))
      } 
    }
  }
  Auc = 0.5 + abs(0.5-Auc)             # if any auc<0.5 than mirror it to the other side of 0.5 auc is a matrix
  if (p.val) {                         # calculate p-values for the test
    n = colSums(matrix(rep(y,nL),nR,nL) == L ) # cound number of elements in each class
    for (i in 1:nP) {                  # go through all permutations of columns in d
      n1 = n[ per[i,1] ]
      n2 = n[ per[i,2] ]
      Auc[i,] = pnorm( n1*n2*(1-Auc[i,]), mean=n1*n2/2, sd=sqrt(n1*n2*(n1+n2+1)/12) )
    }
  }
  # assign dim names
  rLab = NULL
  for (i in 1:nP)                      # go through all permutations of columns in d
    rLab = c(rLab, paste(uL[per[i,1]]," vs. ",uL[per[i,2]], sep=""))
  rownames(Auc) = as.list(rLab)
  colnames(Auc) = colnames(X)
  return (Auc)
}
