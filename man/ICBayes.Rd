\name{ICBayes}
\alias{ICBayes}
\alias{ICBayes.default}
\alias{ICBayes.formula}
\title{PH, PO, and Probit Models for Interval-Censored Data}
\description{Calls the \code{case1ph}, \code{case2ph}, or \code{case2probit} function to fit the corresponding model.
Give point estimates and credible intervals for regression coefficients and estimation and plot of survival functions.}
\usage{
ICBayes(L, ...)

\method{ICBayes}{default}(L, R, model, status, xcov, x_user=NULL, order=2, 
sig0=10, coef_range=5, v0=0.1, a_eta=1, b_eta=1, 
knots=NULL, grids=NULL, conf.int=0.95,
niter=5000, burnin=1000, thin=1, seed=NULL, ...)
\method{ICBayes}{formula}(formula, data, ...)
}
\arguments{
  \item{L}{a column vector of left-points of observed time intervals.}
  \item{R}{a column vector of right-points of observed time intervals. Use NA to denote infinity.}
  \item{model}{a character string specifying the type of model. Possible values are "case1ph", "case2ph", "case2po", and "case2probit".}
  \item{status}{a vector of censoring indicators. If \code{model}="case1ph", then 1=left-censored, 0=right-censored. If \code{model}="case2ph", "case2po", or "case2probit", then 0=left-censored, 1=interval-censored, 2=right-censored.}
  \item{xcov}{a matrix of covariates, each column corresponds to one covariate.}
  \item{x_user}{a vector of covariate values, default is NULL. Need to specify for survival estimation.}
  \item{order}{degree of I-splines (\code{b_l}) (see details). Recommended values are 2-4. Default is 2.}
  \item{sig0}{standard deviation of normal prior for each regression coefficient \code{beta_r}. 
              Used if \code{model}="case1ph", "case1po", or "case2ph". Default is 10.}
  \item{coef_range}{specify support domain of target density for \code{beta_r} using \code{\link[HI]{arms}} (see details). 
        Used if \code{model}="case1ph", "case1po", or "case2ph". Default is 5.}
  \item{v0}{precision of normal prior for \code{gamma_0}. Used if \code{model}="case2po" or "case2probit". Default is 0.1.}
  \item{a_eta}{shape parameter of Gamma prior for \code{gamma_l} (see details). Default is 1.}
  \item{b_eta}{rate parameter of Gamma prior for \code{gamma_l} (see details). Default is 1.}
  \item{knots}{a sequence of points to define I-splines. Default is a sequence of time points from min to max with length=10.}
  \item{grids}{a sequence of points where survival function is to be estimated. Defalult is a sequence of time points from min to max with length=100.}
  \item{conf.int}{level for a two-sided credible interval on coefficient estimate(s). Default is 0.95.}
  \item{niter}{total number of iterations of MCMC chains. Default is 11000.}
  \item{burnin}{number of iterations to discard at the beginning of an MCMC run. Default is 1000.}
  \item{thin}{specify thinning of MCMC draws. Default is 1.}
  \item{seed}{a use-specified random seed. Default is NULL.}
  \item{formula}{a symbolic description of the model to be fit.}
  \item{data}{a data frame containing the variables in the model.}
  \item{\dots}{values passed to other functions.}
}
\details{
For "case1ph", "case1po", and "case2ph" models, function \code{\link[HI]{arms}} is used to sample regression coefficient \code{beta_r}, and \code{coef_range} specifies the support of the \code{indFunc} 
in \code{\link[HI]{arms}}. The baseline cumulative hazard in "case1ph"and "case2ph" models and the baseline odds function in "case1po" are modeled by a linear combination of I-splines: 
                              
                   \code{sum_{l=1}^{k}(gamma_l*b_l)}.

For "case2probit" model, baseline function is modeled by a linear combination of I-splines: 

           \code{gamma_0+sum_{l=1}^{k}(gamma_l*b_l)}.

For "case2probit" model, regression coefficient vector \code{beta} is sampled from a multivariate normal distribution.

For more information, please see reference.
}
\value{
an object of class \code{ICBayes} containing the following elements:
 \item{coef}{a vector of regression coefficient estimates}
 \item{coef_ssd}{a vector of sample standard deviations of regression coefficient estimates}
 \item{coef_ci}{credible intervals for regression coefficients}
 \item{LPML}{log pseudo marginal likelihood for model selection, the larger the better}
 \item{grids}{the sequance of points where baseline survival functions is estimated}
 \item{S0_m}{estimated baseline survival probabilities at \code{grids}}
 \item{S_m}{a \code{length(grids)*G} by 2 matrix that contains estimated survival probabilities at \code{grids} for \code{x_user}, where G is the number of sets of covariate values}
 \item{S_ci}{credible intervals for survival probablities at \code{grids} for \code{x_user}}
 \item{mcmc_beta}{a \code{niter} by p matrix of mcmc chains for regression coefficients, where niter is the number of iterations and p is the number of covariates}
 \item{mcmc_surv}{a \code{niter} by \code{length{grids}*G} matrix of mcmc chains for survival probabilities at \code{grids}, where niter is the number of iterations and G is the number of sets of covariate values}
}
\author{Chun Pan}
\references{
Cai, B., Lin, X., and Wang, L. (2011). Bayesian proportional hazards model for current status data with monotone splines. 
\emph{Computational Statistics and Data Analysis}, \bold{55} 2644-2651.

Lin, X. and Wang, L. (2009). A semiparametric probit model for case 2 interval-censored failure time data. 
\emph{Statistics in Medicine}, \bold{29} 972-981.

Lin, X. and Wang, L. (2011). Bayesian proportional odds model for analyzing current status data: univariate, clustered, and multivariate. 
\emph{Communication in Statistics-Simulation and Computation}, \bold{40} 1171-1181.

Lin, X., Cai, B., Wang, L., and Zhang, Z. (submitted). Bayesian proportional hazards model for general interval-censored data. 
}
\examples{
# To save time in checking package, niter is set to only 500 iterations.
# formula form
data(bcdata)
bcdata<-data.frame(bcdata) # must be a data frame
try<-ICBayes(Surv(L,R,type='interval2')~x1,data=bcdata,
model='case2ph',status=bcdata[,3],x_user=c(0,1),knots=seq(0.1,60.1,length=10),
grids=seq(0.1,60.1,by=1),niter=500,burnin=100)
# general form
try2<-ICBayes(model='case2ph',L=bcdata[,1],R=bcdata[,2],status=bcdata[,3],
xcov=bcdata[,4],x_user=c(0,1),knots=seq(0.1,60.1,length=10),
grids=seq(0.1,60.1,by=1),niter=500,burnin=100)
}
\seealso{\code{\link{case1ph}}, \code{\link{case1po}}, \code{\link{case2ph}}, \code{\link{case2probit}}}
\keyword{case 1 interval-censored data}
\keyword{case 2 interval-censored data}
