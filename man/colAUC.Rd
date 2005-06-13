\name{colAUC}
\alias{colAUC}
\title{Columnwise Area Under ROC Curve (AUC)}
\description{Area Under ROC Curve (AUC) calculated for every column of the 
  matrix.}
\synopsis{colAUC(X, y, p.val=FALSE)}
\usage{
  auc = colAUC(X, y)
  p   = colAUC(X, y, p.val=TRUE)
}

\arguments{
  \item{X}{A matrix or data frame. Rows contain samples 
    and columns contain features/variables.}
  \item{y}{Class labels for the \code{X} data samples. 
    A response vector with one label for each row/component of \code{X}.
    Can be either a factor, string or a numeric vector.}
  \item{p.val}{a boolean flag: if set to TRUE than "Wilcoxon rank sum test" 
    p-values (see \code{\link{wilcox.test}}) will be returned instead of AUC 
    values}
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
  
  This measure is very similar to Wilcoxon rank sum test (see 
  \code{\link{wilcox.test}}), which is also called 
  Mann-Whitney test. Wilcoxon-test's p-value can be calculated by 
  \code{p=pnorm( n1*n2*(1-auc), mean=n1*n2/2, sd=sqrt(n1*n2*(n1+n2+1)/12) )} 
  where \code{n1} and \code{n2} are numbers of elements in two classes being 
  compared.
  
  The main purpose of this function was to calculate AUC's of large number of 
  features, fast. It is being used to help with classification of protein mass 
  spectra data 
  that often have up to 50K features, as a fast and dirty way of lowering 
  dimensionality of the data before applying standard classification algorithms 
  like \code{nnet} or \code{svd}.
}

\value{
  An output is a single matrix with the same number of columns as \code{X} and 
  "n choose 2" ( \eqn{\frac{n!}{(n-2)! 2!}}{n!/((n-2)! 2!)} ) number of rows, 
  where n is number of unique labels in \code{y} list. For example, if \code{y} 
  contains only two unique class labels ( \code{length(unique(lab))==2} ) than
  output 
  matrix will  have a single row containing AUC of each column. If more than 
  two unique labels are present than AUC is calculated for every possible 
  pairing of classes ("n choose 2" of them).
} 

\references{
  \itemize{
    \item Mason, S.J. and Graham, N.E. (1982) \emph{Areas beneath the relative 
     operating characteristics (ROC) and relative operating levels (ROL) 
     curves: Statistical significance and interpretation},  Q. J. R. 
     Meteorol. Soc. textbf{30} 291-303. 
     \item See 
       \url{http://www.medicine.mcgill.ca/epidemiology/hanley/software/} 
       to find articles below: 
     \itemize{
       \item Hanley and McNeil (1982), \emph{The Meaning and Use of the Area 
         under a Receiver Operating Characteristic (ROC) Curve}, 
         Radiology 143: 29-36.
       \item Hanley and McNeil (1983), \emph{A Method of Comparing the Areas  
         under ROC curves derived from same cases}, Radiology 148: 839-843.
       \item McNeil and Hanley (1984), \emph{Statistical Approaches to the  
         Analysis of ROC curves}, Medical Decision Making 4(2): 136-149.
     }
  }
} 

\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com}} 

\seealso{
  \code{\link[ROC]{AUC}} from \pkg{ROC} package, 
  \code{\link[verification]{roc.area}} from \pkg{verification} package, 
  \code{\link{wilcox.test}}
}

\examples{
  # load MASS library with "cats" data set that have following columns: sex, 
  # body weight, hart weight
  library(MASS)
  data(cats)
  colAUC(cats[,2:3], cats[,1]) 
  
  # compare with examples from roc.area function: using Data from Mason and Graham (2002).
  a<- (1981:1995)
  b<- c(0,0,0,1,1,1,0,1,1,0,0,0,0,1,1)
  c<- c(.8, .8, 0, 1,1,.6, .4, .8, 0, 0, .2, 0, 0, 1,1)
  d<- c(.928,.576, .008, .944, .832, .816, .136, .584, .032, .016, .28, .024, 0, .984, .952)
  A<- data.frame(a,b,c,d)
  names(A)<- c("year", "event", "p1", "p2")
  if (library(verification, logical.return=TRUE)) {
    roc.area(A$event, A$p1)           # for model with ties
    roc.area(A$event, A$p2)           # for model without ties
  }
  wilcox.test(p2~event, data=A)
  # colAUC output is the same as roc.area's A.tilda values
  colAUC(A[,3:4], A$event) 
  # colAUC output is the same as roc.area's  and wilcox.test's p values
  colAUC(A[,3:4], A$event, p.val=TRUE) 
  
  # example of 3-class data
  data(iris)
  colAUC(iris[,-5], iris[,5])
}

\keyword{univar}
