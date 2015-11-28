\name{lungdata}
\alias{lungdata}
\docType{data}
\title{Lung cancer data}
\description{A case 1 interval-censored data set first presented in Hoel and Walberg (1972) and can be found in Sun (2006, page 6). In the study, 144 male RFM mice 
were raised under two conditions: conventional environment (96 mice) and germfree environment (48 mice). 
Each mouse was "sacrificed" at a random time to see if it had lung tumors. Time was measured in days.  
}
\usage{data(lungdata)}
\format{
  A matrix with 144 rows and 4 columns. Each row (L, R, status, treatment) corresponds to a mouse in the study.
  \describe{
    \item{\code{L}}{left-points of observed intervals}
    \item{\code{R}}{right-points of observed intervals}
    \item{\code{status}}{censorship indicator: 1=left-censor and 0=right-censor}
    \item{\code{treatment}}{treatment indicator: 1=conventional environment, 2=germfree environment}
  }
}
\source{
Hoel, D. G. and Walberg, H. E. (1972). Statistical analysis of survival experiments. 
\emph{Journal of the National Cancer Institute} \bold{49} 361-372.

Sun, J. (2006). \emph{The Statistical Analysis of Interval-censored Failure Time Data}. Springer, New York.
}
\examples{
data(lungdata)
}
\keyword{datasets}
