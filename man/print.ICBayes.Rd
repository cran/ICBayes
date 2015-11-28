\name{print.ICBayes}
\alias{print.ICBayes}
\title{print method for ICBayes object}
\description{Print output from function \code{\link{ICBayes}} }
\usage{\method{print}{ICBayes}(x, ...)}
\arguments{
  \item{x}{an object of class \code{"ICBayes"}, i.e., a fitted model}
  \item{...}{other arguments if any}
}
\value{
an object of class \code{ICBayes} containing the following elements
 \item{coef}{a named vector of coefficient estimates}
 \item{coef_ssd}{a named vector of sample standard deviations of coefficient estimates}
 \item{coef_ci}{a named matrix of credible intervals for coefficients}
 \item{LPML}{log pseudo marginal likelihood}
}