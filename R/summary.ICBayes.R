summary.ICBayes <-
function(object, ...){

   coef.table<-cbind(object$coef, object$coef_ssd, object$coef_ci[,1], object$coef_ci[,2])
   dimnames(coef.table)<-list(names(object$coef),
			      c('Mean','Std. Dev.',
				paste(object$conf.int*100,'%CI-Low',sep=''),
				paste(object$conf.int*100,'%CI-Upp',sep='')))
   cat("Call:\n")
   print(object$call)

   cat("\nPosterior inference of regression coefficients\n")
   print.default(coef.table)


   cat("\nLog pseudo marginal likelihood: LPML=", object$LPML, sep="")
   cat("\nNegative log-likelihood: NLLK=", -object$LPML, sep="")
   cat("\nNumber of subjects: n=", object$n, "\n", sep="")

}
