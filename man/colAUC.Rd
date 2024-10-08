\name{colAUC}
\alias{colAUC}
\title{Column-wise Area Under ROC Curve (AUC)}
\description{Calculate Area Under the ROC Curve (AUC) for every column of a 
  matrix. Also, can be used to plot the ROC curves.}
\usage{
  colAUC(X, y, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
}

\arguments{
  \item{X}{A matrix or data frame. Rows contain samples 
    and columns contain features/variables.}
  \item{y}{Class labels for the \code{X} data samples. 
    A response vector with one label for each row/component of \code{X}.
    Can be either a factor, string or a numeric vector.}
  \item{plotROC}{Plot ROC curves. Use only for small number of features. 
    If \code{TRUE}, will set \code{alg} to "ROC".}
  \item{alg}{Algorithm to use: "ROC" integrates ROC curves, while "Wilcoxon"
    uses Wilcoxon Rank Sum Test to get the same results. Default "Wilcoxon" is
    faster. This argument is mostly provided for verification.}
}

\details{
  AUC is a very useful measure of similarity between two classes measuring area
  under "Receiver Operating Characteristic" or ROC curve.
  In case of data with no ties all sections of ROC curve are either horizontal
  or vertical, in case of data with ties diagonal 
  sections can also occur. Area under the ROC curve is calculated using 
  \code{\link{trapz}} function. AUC is always in between 0.5 
  (two classes are statistically identical) and 1.0 (there is a threshold value
   that can achieve a perfect separation between the classes).
  
  Area under ROC Curve (AUC) measure is very similar to Wilcoxon Rank Sum Test 
  (see \code{\link{wilcox.test}}) and Mann-Whitney U Test. 
  
  There are numerous other functions for calculating AUC in other packages. 
  Unfortunately none of them had all the properties that were needed for 
  classification preprocessing, to lower the dimensionality of the data (from 
  tens of thousands to hundreds) before applying standard classification 
  algorithms. 
  
  The main properties of this code are: 
  \itemize{
    \item Ability to work with multi-dimensional data (\code{X} can have many 
          columns).
    \item Ability to work with multi-class datasets (\code{y} can have more 
          than 2 different values).
    \item Speed - this code was written to calculate AUC's of large number of 
          features, fast.
    \item Returned AUC is always bigger than 0.5, which is equivalent of 
    testing for each feature \code{colAUC(x,y)} and \code{colAUC(-x,y)} and
    returning the value of the bigger one.
  }
  If those properties do not fit your problem, see "See Also" and "Examples" 
  sections for AUC 
  functions in other packages that might be a better fit for your needs. 
}

\value{
  An output is a single matrix with the same number of columns as \code{X} and 
  "n choose 2" ( \eqn{\frac{n!}{(n-2)! 2!} = n(n-1)/2}{n!/((n-2)! 2!) = n(n-1)/2} )
  number of rows, 
  where n is number of unique labels in \code{y} list. For example, if \code{y} 
  contains only two unique class labels ( \code{length(unique(lab))==2} ) than
  output 
  matrix will  have a single row containing AUC of each column. If more than 
  two unique labels are present than AUC is calculated for every possible 
  pairing of classes ("n choose 2" of them).
  
  For multi-class AUC "Total AUC" as defined by Hand & Till (2001) can be 
  calculated by \code{\link{colMeans}(auc)}.  
} 

\references{
  \itemize{
    \item Mason, S.J. and Graham, N.E. (1982) Areas beneath the relative 
     operating characteristics (ROC) and relative operating levels (ROL) 
     curves: Statistical significance and interpretation,  Q. J. R. 
     Meteorol. Soc. 30 291-303. 
    \item Fawcett, Tom (2004) ROC Graphs: Notes and Practical 
       Considerations for Researchers
    \item  Hand, David and Till, Robert (2001) A Simple 
          Generalization of the Area Under the ROC Curve for Multiple Class 
          Classification Problems;  Machine Learning 45(2): 171-186
    \item See 
       \url{http://www.medicine.mcgill.ca/epidemiology/hanley/software/} 
       to find articles below: 
     \itemize{
       \item Hanley and McNeil (1982), The Meaning and Use of the Area 
         under a Receiver Operating Characteristic (ROC) Curve, 
         Radiology 143: 29-36.
       \item Hanley and McNeil (1983), A Method of Comparing the Areas  
         under ROC curves derived from same cases, Radiology 148: 839-843.
       \item McNeil and Hanley (1984), Statistical Approaches to the  
         Analysis of ROC curves, Medical Decision Making 4(2): 136-149.
     }
  }
} 

\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com}} 

\seealso{
  \itemize{
  \item \code{\link{wilcox.test}} and \code{\link{pwilcox}}
  \item \code{\link[exactRankTests]{wilcox.exact}} from \pkg{exactRankTests} package
  \item \code{\link[coin]{wilcox_test}} from \pkg{coin} package
  \item \code{\link[ROCR]{performance}} from \pkg{ROCR} package
  \item \code{\link[Epi]{ROC}} from \pkg{Epi} package
  \item \code{\link[verification]{roc.area}} from \pkg{verification} package
  \item \code{\link[Hmisc]{rcorr.cens}} from \pkg{Hmisc} package
  }
}

\examples{
# Load MASS library with "cats" data set that have following columns: sex, body
# weight, hart weight. Calculate how good weights are in predicting sex of cats.
# 2 classes; 2 features; 144 samples
library(MASS); data(cats);
colAUC(cats[,2:3], cats[,1], plotROC=TRUE) 

# Load rpart library with "kyphosis" data set that records if kyphosis
# deformation was present after corrective surgery. Calculate how good age, 
# number and position of vertebrae are in predicting successful operation. 
# 2 classes; 3 features; 81 samples
library(rpart); data(kyphosis);
colAUC(kyphosis[,2:4], kyphosis[,1], plotROC=TRUE)

# Example of 3-class 4-feature 150-sample iris data
data(iris)
colAUC(iris[,-5], iris[,5], plotROC=TRUE)
cat("Total AUC: \n"); 
colMeans(colAUC(iris[,-5], iris[,5]))

# Test plots in case of data without column names
Iris = as.matrix(iris[,-5])
dim(Iris) = c(600,1)
dim(Iris) = c(150,4)
colAUC(Iris, iris[,5], plotROC=TRUE)

# Compare calAUC with other functions designed for similar purpose
auc = matrix(NA,12,3)
rownames(auc) = c("colAUC(alg='ROC')", "colAUC(alg='Wilcox')", "sum(rank)",
    "wilcox.test", "wilcox_test", "wilcox.exact", "roc.area", "AUC", 
    "performance", "ROC", "auROC", "rcorr.cens")
colnames(auc) = c("AUC(x)", "AUC(-x)", "AUC(x+noise)")
X = cbind(cats[,2], -cats[,2], cats[,2]+rnorm(nrow(cats)) )
y = ifelse(cats[,1]=='F',0,1)
for (i in 1:3) {
  x = X[,i]
  x1 = x[y==1]; n1 = length(x1);                 # prepare input data ...
  x2 = x[y==0]; n2 = length(x2);                 # ... into required format
  data = data.frame(x=x,y=factor(y))
  auc[1,i] = colAUC(x, y, alg="ROC") 
  auc[2,i] = colAUC(x, y, alg="Wilcox")
  r = rank(c(x1,x2))
  auc[3,i] = (sum(r[1:n1]) - n1*(n1+1)/2) / (n1*n2)
  auc[4,i] = wilcox.test(x1, x2, exact=0)$statistic / (n1*n2) 
  \dontrun{
  if (require("coin"))
    auc[5,i] = statistic(wilcox_test(x~y, data=data)) / (n1*n2) 
  if (require("exactRankTests"))  
    auc[6,i] = wilcox.exact(x, y, exact=0)$statistic / (n1*n2) 
  if (require("verification"))
    auc[7,i] = roc.area(y, x)$A.tilda 
  if (require("ROC")) 
    auc[8,i] = AUC(rocdemo.sca(y, x, dxrule.sca))    
  if (require("ROCR")) 
    auc[9,i] = performance(prediction( x, y),"auc")@y.values[[1]]
  if (require("Epi"))   auc[10,i] = ROC(x,y,grid=0)$AUC
  if (require("limma")) auc[11,i] = auROC(y, x)
  if (require("Hmisc")) auc[12,i] = rcorr.cens(x, y)[1]
  }
}
print(auc)
stopifnot(auc[1, ]==auc[2, ])   # results of 2 alg's in colAUC must be the same
stopifnot(auc[1,1]==auc[3,1])   # compare with wilcox.test results

# time trials
x = matrix(runif(100*1000),100,1000)
y = (runif(100)>0.5)
system.time(colAUC(x,y,alg="ROC"   ))
system.time(colAUC(x,y,alg="Wilcox"))
}

\keyword{univar}
