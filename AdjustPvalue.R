#########################################################################
#               Adaptive version of P-value Adjustment                  #
#########################################################################
# - p           the number of biomarkers                                #
# - tau.min    the parameter for adjusting p-values                     #
#########################################################################

adaptive<-function(p){
  tau.min=tau.min
  if(sum(p<0.05)==0){adj.p=1}
  if(sum(p<0.05)!=0){
    tau=seq(0.05,1,0.1)
    Q=rep(0,length(tau))
    for (tau.i in 1:length(tau)) {
      Q[tau.i]=min(1,quantile(p/tau[tau.i],prob=tau[tau.i],na.rm=T))
    }
    adj.p=min(1,(1-log(tau.min))*min(Q))
  }
  return (adj.p)
}
