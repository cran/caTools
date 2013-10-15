\name{LogitBoost}
\alias{LogitBoost}
\title{LogitBoost Classification Algorithm}
\description{Train logitboost classification algorithm using decision 
  stumps (one node decision trees) as weak learners.  }
\usage{LogitBoost(xlearn, ylearn, nIter=ncol(xlearn))}

\arguments{
  \item{xlearn}{A matrix or data frame with training data. Rows contain samples 
    and columns contain features}
  \item{ylearn}{Class labels for the training data samples. 
    A response vector with one label for each row/component of \code{xlearn}.
    Can be either a factor, string or a numeric vector.}
  \item{nIter}{An integer, describing the number of iterations for
     which boosting should be run, or number of decision stumps that will be 
     used.}
}

\details{
  The function was adapted from logitboost.R function written by Marcel 
  Dettling. See references and "See Also" section. The code was modified in 
  order to make it much faster for very large data sets. The speed-up was 
  achieved by implementing a internal version of decision stump classifier 
  instead of using calls to \code{\link[rpart]{rpart}}. That way, some of the most time 
  consuming operations were precomputed once, instead of performing them at 
  each iteration. Another difference is that training and testing phases of the 
  classification process were split into separate functions.
}

\value{
  An object of class "LogitBoost" including components: 
  \item{Stump}{List of decision stumps (one node decision trees) used:
    \itemize{
      \item column 1: feature numbers or each stump, or which column each stump 
       operates on
      \item column 2: threshold to be used for that column
      \item column 3: bigger/smaller info: 1 means that if values in the column 
       are above threshold than corresponding samples will be labeled as 
       \code{lablist[1]}. Value "-1" means the opposite.
    }
    If there are more than two classes, than several "Stumps" will be
    \code{cbind}'ed
   }
  \item{lablist}{names of each class}
} 

\references{
  Dettling and Buhlmann (2002), \emph{Boosting for Tumor Classification of Gene 
    Expression Data}, available on the web page 
    \url{http://stat.ethz.ch/~dettling/boosting.html}. 
  
  \url{http://www.cs.princeton.edu/~schapire/boost.html}
}

\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com}} 

\seealso{
 \itemize{
   \item \code{\link{predict.LogitBoost}} has prediction half of LogitBoost code
   \item \code{logitboost} function from \pkg{logitboost} library (not in CRAN
   or BioConductor but can be found at 
   \url{http://stat.ethz.ch/~dettling/boosting.html}) is very similar but much
   slower on very large datasets. It also perform optional cross-validation.
 }
}

\examples{
  data(iris)
  Data  = iris[,-5]
  Label = iris[, 5]
  
  # basic interface
  model = LogitBoost(Data, Label, nIter=20)
  Lab   = predict(model, Data)
  Prob  = predict(model, Data, type="raw")
  t     = cbind(Lab, Prob)
  t[1:10, ]

  # two alternative call syntax
  p=predict(model,Data)
  q=predict.LogitBoost(model,Data)
  pp=p[!is.na(p)]; qq=q[!is.na(q)]
  stopifnot(pp == qq)

  # accuracy increases with nIter (at least for train set)
  table(predict(model, Data, nIter= 2), Label)
  table(predict(model, Data, nIter=10), Label)
  table(predict(model, Data),           Label)
  
  # example of spliting the data into train and test set
  mask = sample.split(Label)
  model = LogitBoost(Data[mask,], Label[mask], nIter=10)
  table(predict(model, Data[!mask,], nIter=2), Label[!mask])
  table(predict(model, Data[!mask,]),          Label[!mask])
}

\keyword{classif}
