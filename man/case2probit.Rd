\name{case2probit}
\alias{case2probit}
\title{Probit model for general interval-censored data}
\description{Fit probit model to general interval-censored data. Use MCMC method to estimate regression coefficients, baseline survival, 
and survival function at user-specified covariate values.}
\usage{
case2probit(L, R, status, xcov, x_user, order,
	    v0, a_eta, b_eta, knots, grids, niter)
}
\arguments{
  \item{L}{a numeric vector of left timepoints of observed time intervals.}
  \item{R}{a numeric vector of right timepoints of observed time intervals.}
  \item{status}{a vector of censoring indicators: 0=left-censored, 1=interval-censored, 2=right-censored.}
  \item{xcov}{a matrix of covariates, each column corresponds to one covariate.}
  \item{x_user}{a vector of user specified covariate values.}
  \item{order}{degree of I-splines (\code{b_l}) (see details). Recommended values are 2-4.} 
  \item{v0}{precision of normal prior for \code{gamma_0}.}
  \item{a_eta}{shape parameter of Gamma prior for \code{gamma_l} (see details).}
  \item{b_eta}{rate parameter of Gamma prior for \code{gamma_l} (see details).}
  \item{knots}{a sequence of points to define I-splines.}
  \item{grids}{a sequence of points where baseline survival function is to be estimated. Default is minimum observed time points.}
  \item{niter}{total number of iterations of MCMC chains.}
}
\details{
The baseline function is modeled by a linear combination of I-splines: 

           \code{gamma_0+sum_{l=1}^{k}(gamma_l*b_l)}.

Regression coefficient vector \code{beta} is sampled from a multivariate normal distribution.
For more information, please see reference.
}
\value{
a list containing the following elements:
 \item{parbeta}{a \code{niter} by \code{p} matrix of MCMC draws of \code{beta_r}, r=1, ..., p. }
 \item{parsurv0}{a \code{niter} by \code{length(grids)} matrix, each row contains the baseline survival at \code{grids} from one iteration .}
 \item{grids}{a sequence of points where baseline survival is estimated.}
}
\references{
Lin, X. and Wang, L. (2009). A semiparametric probit model for case 2 interval-censored failure time data. 
\emph{Statistics in Medicine} \bold{29} 972-981. 
}
\author{Lianming Wang and Xiaoyan Lin. R version by Bo Cai.}
