\name{case1ph}
\alias{case1ph}
\title{PH model for case 1 interval-censored data}
\description{Fit proportional hazards model for case 1 interval-censored data. Use MCMC method to estimate regression coefficients, baseline survival, 
and survival at user-specified covariate values.}
\usage{
case1ph(L, R, status, xcov, x_user, order, sig0, coef_range, 
	a_eta, b_eta, knots, grids, niter, seed)
}
\arguments{
  \item{L}{a numeric vector of left timepoints of observed time intervals.}
  \item{R}{a numeric vector of right timepoints of observed time intervals.}
  \item{status}{a vector of censoring indicators: 1=left-censored, 0=right-censored.}
  \item{xcov}{a matrix of covariates, each column corresponds to one covariate.}
  \item{x_user}{a user specified vector of covariate values}
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
The baseline cumulative hazard is approximated by a linear combination of I-splines: 

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
 \item{grids}{a sequence of points where baseline survival is estimated.}
}
\references{
Cai, B., Lin, X., and Wang, L. (2011). Bayesian proportional hazards model for current status data with monotone splines. 
\emph{Computational Statistics and Data Analysis}, \bold{55} 2644-2651.
}
\author{Bo Cai}

