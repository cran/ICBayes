\name{bcdata}
\alias{bcdata}
\docType{data}
\title{Breast Cosmesis Data}
\description{A general interval-censored data set analyzed in Finkelstein and Wolfe (1985) and can be found in Sun (2006, page 7). Early breast cancer patients 
treated with radiotheraph alone or radiotherapy with adjuvent chemotherapy were examined periodically for breast retraction. Time was measured in months.}
\usage{data(bcdata)}
\format{
  A matrix with 96 rows and 4 columns. Each row (L, R, status, x1) corresponds to a patient in the study.
  \describe{
    \item{\code{L}}{a numeric vector of left-points of observed time intervals}
    \item{\code{R}}{a numeric vector of right-points of observed time intervals}
    \item{\code{status}}{a vector of censorship indicators: 0=left-censored, 1=interval-censored, and 2=right-censored}
    \item{\code{x1}}{a vector of treatment indicators: 0=radiotherapy alone, 1=radiotherapy with adjuvent chemotherapy}
  }
}
\source{
Finkelstein, D. M. and Wolfe, R. A. (1985). A semiparametric model for regression analysis of interval-censored failure time data. \emph{
Biometrics} \bold{41} 933-945.

Sun, J. (2006). \emph{The Statistical Analysis of Interval-censored Failure Time Data}. Springer, New York.
}
\examples{
data(bcdata)
}
\keyword{datasets}
