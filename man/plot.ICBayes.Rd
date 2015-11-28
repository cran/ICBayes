\name{plot.ICBayes}
\alias{plot.ICBayes}
\title{Plot baseline survival function}
\description{Plot estimated baseline survival function at \code{grids}, which are stored in the \code{ICBayes} object.}
\usage{\method{plot}{ICBayes}(x, y, ...)}
\arguments{
  \item{x}{a sequence of points where baseline survival probabilities are estimated}
  \item{y}{estiamted baseline survival at grids}
  \item{...}{other arguments}
}
\value{A plot of baseline survival function.}
\examples{
data(bcdata)
try<-ICBayes(Surv(L,R,type='interval2')~x1,data=data.frame(bcdata),
model='case2ph',status=bcdata[,3],p=1,x_user=c(1),
knots=seq(0.1,60.1,length=10),grids=seq(0.1,60.1,by=1),niter=500,burnin=100)
plot.ICBayes(try$grids,try$S0_m)
}
