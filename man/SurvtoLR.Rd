\name{SurvtoLR}
\alias{SurvtoLR}
\title{Transform Surv object to data matrix with L and R columns}
\description{Take a \code{\link[survival]{Surv}} object and transforms it into a data matrix with two columns, L and R, 
representing the left and right points of observed time intervals. For right-censored data, R = NA.}
\usage{SurvtoLR(x)}
\arguments{
  \item{x}{a \code{\link[survival]{Surv}} object}
}
\details{The input Surv object should be in the form of \code{Surv(L,R,type='interval2')}, where R = NA for right-censored data.}
\value{A data matrix with two variables:
  \item{L}{left-points of observed time intervals}
  \item{R}{right-points of observed time intervals}
}
\references{
Michael P. Fay, Pamela A. Shaw (2010). Exact and Asymptotic Weighted Logrank Tests for Interval Censored Data: The interval R Package. 
\emph{Journal of Statistical Software}, \bold{36} 1-34.
}
\examples{
library(survival)
L<-c(45,6,0,46)
R<-c(NA,10,7,NA)
y<-Surv(L,R,type='interval2')
SurvtoLR(y)
}