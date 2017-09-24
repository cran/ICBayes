summary.ICBayes <-
function(object, ...){
   cat("Call:\n")
   print(object$call)
   cat("\ncoef #Estimted regression coefficients\n")
   print(object$coef)
   cat("\ncoef_ssd #Sample standard deviation\n")
   print(object$coef_ssd)
   cat("\ncoef_ci #Credible interval of coefficients\n")
   print(object$coef_ci)
   cat("\nLPML #Log pseudo marginal likelihood\n")
   print(object$LPML)
}
