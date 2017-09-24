\name{summary.ICBayes}
\alias{summary.ICBayes}
\title{summary method for class ICBayes}
\description{present output from function \code{\link{ICBayes}}}
\usage{\method{summary}{ICBayes}(object, ...)}

\arguments{
  \item{object}{class ICBayes object}
  \item{\dots}{other arguments if any}
}

\value{
an object of class \code{ICBayes} containing the following elements
 \item{coef}{a named vector of coefficient estimates}
 \item{coef_ssd}{a named vector of sample standard deviations of coefficient estimates}
 \item{coef_ci}{a named matrix of credible intervals for coefficients}
 \item{LPML}{log pseudo marginal likelihood}
}
