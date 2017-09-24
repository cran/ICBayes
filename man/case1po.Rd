\name{case1po}
\alias{case1po}
\title{PO model for case 1 interval-censored data}
\description{Fit proportional odds model for case 1 interval-censored data. Use MCMC method to estiamte regression coefficients, baseline survival, 
and survival function at user-specified covariate values.}
\usage{
case1po(L, R, status, xcov, x_user, order, sig0, coef_range, 
	a_eta, b_eta, knots, grids, niter, seed)
}
\arguments{
  \item{L}{a numeric vector of left timepoints of observed time intervals.}
  \item{R}{a numeric vector of right timepoints of observed time intervals.}
  \item{status}{a vector of censoring indicators: 1=left-censored, 0=right-censored.}
  \item{xcov}{a matrix of covariates, each column corresponds to one covariate.}
  \item{x_user}{a vector of user specified covariate values.}
  \item{order}{degree of I-splines (\code{b_l}) (see details). Recommended values are 2-4.}
  \item{sig0}{standard deviation of normal prior for each regression coefficient \code{beta_r}.}
  \item{coef_range}{specify support domain of target density for \code{beta_r} sampled by \code{\link[HI]{arms}} (see details).}
   \item{a_eta}{shape parameter of Gamma prior for \code{gamma_l} (see details).}
  \item{b_eta}{rate parameter of Gamma prior for \code{gamma_l} (see details).}
  \item{knots}{a sequence of points to define I-splines.}
  \item{grids}{a sequence of points where baseline survival function is to be estimated.}
  \item{niter}{total number of iterations of MCMC chains.}
  \item{seed}{a user specified random seed, default is NULL.}
}
\details{
The baseline odds function is approximated by a linear combination of I-splines: 

                    \code{sum_{l=1}^{k}(gamma_l*b_l)}.

Function \code{\link[HI]{arms}} is used to sample each regression coefficient \code{beta_r}, and \code{coef_range} specifies the support of 
the \code{indFunc} 
in \code{\link[HI]{arms}}.
}
\value{
a list containing the following elements:
 \item{parbeta}{a \code{niter} by \code{p} matrix of MCMC draws of \code{beta_r}, r=1, ..., p. }
 \item{parsurv0}{a \code{niter} by \code{length(grids)} matrix, each row contains the baseline survival at \code{grids} from one iteration.}
 \item{parsurv}{a \code{niter} by \code{length(grids)*G} matrix, each row contains the survival at \code{grids} from one iteration. 
	G is the number of sets of user-specified covariate values.}
 \item{parfinv}{a \code{niter} by n matrix, each row contains the inverse PDF of observed interval-censored data from one iteration. 
	This is used for computing LPML later.}
}
\references{
Lin, X. and Wang, L. (2011). Bayesian proportional odds model for analyzing current status data: univariate, clustered, and multivariate. 
\emph{Communication in Statistics-Simulation and Computation}, \bold{40} 1171-1181.
}
\author{Xiaoyan Lin}
