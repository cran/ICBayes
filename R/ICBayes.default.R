ICBayes.default <-
function(
L,
R,
model,
status,
xcov,
x_user=NULL,
order=2,
sig0=10,
coef_range=5,
m0=0,
v0=0.1,
a_eta=1,
b_eta=1,
knots=NULL,
grids=NULL,
conf.int=0.95,
niter=5000,
burnin=1000,
thin=1,
seed=NULL,
...){

   if (model=='case2ph'){
est<-case2ph(L,R,status,xcov,x_user,order,sig0,coef_range,
a_eta,b_eta,knots,grids,niter,seed)
   } else
   if (model=='case1ph'){
est<-case1ph(L,R,status,xcov,x_user,order,sig0,coef_range,
a_eta,b_eta,knots,grids,niter,seed)
   } else
   if (model=='case1po'){
est<-case1po(L,R,status,xcov,x_user,order,sig0,coef_range,
a_eta,b_eta,knots,grids,niter,seed)
   } else
   if (model=='case2probit'){
est<-case2probit(L,R,status,xcov,x_user,order,m0,
v0,a_eta,b_eta,knots,grids,niter,seed)
   }

   p=ncol(as.matrix(xcov))
   G=length(x_user)/p
   kgrids=length(grids)

   wbeta=as.matrix(est$parbeta[seq((burnin+thin),niter,by=thin),],ncol=p) # thinned beta samples
   wparsurv0=as.matrix(est$parsurv0[seq((burnin+thin),niter,by=thin),],ncol=kgrids)
   wparsurv=as.matrix(est$parsurv[seq((burnin+thin),niter,by=thin),],ncol=kgrids*G)

   if (model=='case1ph' | model=='case2ph'){
   wparhaz0=as.matrix(est$parhaz0[seq((burnin+thin),niter,by=thin),],ncol=kgrids)
   wparhaz=as.matrix(est$parhaz[seq((burnin+thin),niter,by=thin),],ncol=kgrids*G)}

   coef<-apply(wbeta,2,mean)
   names(coef)<-colnames(xcov)
   coef_ssd<-apply(wbeta,2,sd)
   credinterv<-array(rep(0,p*2),dim=c(p,2))
   credinterv_S<-array(rep(0,kgrids*G*2),dim=c(kgrids*G,2))
   credinterv_h<-array(rep(0,kgrids*G*2),dim=c(kgrids*G,2))

   colnames(credinterv)<-c(paste(100*(1-conf.int)/2,"%CI"),
			   paste(100*(0.5+conf.int/2),"%CI"))
   for (r in 1:p){
	credinterv[r,]<-quantile(wbeta[,r],c((1-conf.int)/2,0.5+conf.int/2))
   }

   colnames(credinterv_S)<-c(paste(100*(1-conf.int)/2,"%CI"),
			     paste(100*(0.5+conf.int/2),"%CI"))
   for (h in 1:(kgrids*G)){
	credinterv_S[h,]<-quantile(wparsurv[,h],c((1-conf.int)/2,0.5+conf.int/2))
   }

   if (model=='case1ph' | model=='case2ph'){
   colnames(credinterv_h)<-c(paste(100*(1-conf.int)/2,"%CI"),
			     paste(100*(0.5+conf.int/2),"%CI"))
   for (h in 1:(kgrids*G)){
	credinterv_h[h,]<-quantile(wparhaz[,h],c((1-conf.int)/2,0.5+conf.int/2))}
   }

   CPO=1/apply(est$parfinv[seq((burnin+thin),niter,by=thin),],2,mean)
   LPML=sum(log(CPO))      # bigger is better   


  if (model=='case1ph' | model=='case2ph'){
  output<-list(
  n = length(L),
  coef = coef,
  coef_ssd = coef_ssd,
  conf.int = conf.int,
  coef_ci = credinterv,
  grids = est$grids,
  S0_m = apply(wparsurv0,2,mean),
  S_m = apply(wparsurv,2,mean),
  S_ci = credinterv_S,
  h0_m = apply(wparhaz0,2,mean),
  h_m = apply(wparhaz,2,mean),
  h_ci = credinterv_h,
  LPML = LPML,
  mcmc_beta = mcmc(est$parbeta),
  mcmc_surv = mcmc(est$parsurv))
  } 


  if (model=='case1po' | model=='case2probit'){
  output<-list(
  n = length(L),
  coef = coef,
  coef_ssd = coef_ssd,
  conf.int = conf.int,
  coef_ci = credinterv,
  grids = est$grids,
  S0_m = apply(wparsurv0,2,mean),
  S_m = apply(wparsurv,2,mean),
  S_ci = credinterv_S,
  LPML = LPML,
  mcmc_beta = mcmc(est$parbeta),
  mcmc_surv = mcmc(est$parsurv))
  } 


  output$call<-match.call()
  class(output)<-"ICBayes"
  output
}
