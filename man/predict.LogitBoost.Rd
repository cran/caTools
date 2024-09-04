\name{predict.LogitBoost}
\alias{predict.LogitBoost}
\title{Prediction Based on LogitBoost Classification Algorithm}
\description{Prediction or Testing using logitboost classification algorithm}
\usage{\method{predict}{LogitBoost}(object, xtest, type = c("class", "raw"), nIter=NA, \dots)}

\arguments{
  \item{object}{An object of class "LogitBoost" see "Value" section of 
    \code{\link{LogitBoost}} for details}
  \item{xtest}{A matrix or data frame with test data. Rows contain samples 
    and columns contain features}
  \item{type}{See "Value" section}
  \item{nIter}{An optional integer, used to lower number of iterations 
    (decision stumps) used in the decision making. If not provided than the 
    number will be the same as the one provided in \code{\link{LogitBoost}}. 
    If provided than the results will be the same as running 
    \code{\link{LogitBoost}} with fewer iterations. }
  \item{\dots}{not used but needed for compatibility with generic predict 
    method}
}

\details{
  Logitboost algorithm relies on a voting scheme to make classifications. Many
  (\code{nIter} of them) week classifiers are applied to each sample and their
  findings are used as votes to make the final classification. The class with 
  the most votes "wins". However, with this scheme it is common for two cases 
  have a tie (the same number of votes), especially if number of iterations is 
  even. In that case NA is returned, instead of a label. 
}

\value{
 If type = "class" (default) label of the class with maximal probability is 
 returned for each sample. If type = "raw", the a-posterior probabilities for 
 each class are returned.
} 

\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com}} 
\seealso{\code{\link{LogitBoost}} has training half of LogitBoost code}
\examples{# See LogitBoost example}
\keyword{classif}
