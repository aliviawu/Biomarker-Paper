#########################################################################
# Simulation function for survival data with biomarkers                 #
#########################################################################
# - n          the number of patients                                   #
# - p          the number of biomarkers                                 #
# - sim        the number of times of simulation                        #
# - B          the number of iterations                                 #
# - t.prob     the treatment assignment probability                     #
# - rho        the mutual correlation between pairwise biomarkers       #
# - h0         the baseline hazard in years                             #
# - alpha0     the treatment effect                                     #
# - active1    the active prognostic biomarkers                         #
# - active2    the active treatment-biomarker interaction               #
# - alpha      the effect of prognostic biomarkers                      #
# - gamma      the effect of treatment-effect modifiers                 #
# - theta0     the parameter for censoring time generation              #
# - tau.min    the parameter for adjusting p-values                     #
# - m          the allocation rate for multi-splitting                  #
#########################################################################
runsims<- function(n = n, p = p, sim = sim, B = B, t.prob = t.prob, rho = rho, h0 = h0,
                   alpha0 = alpha0, active1 = active1, active2 = active2, alpha = alpha, 
                   gamma = gamma, theta0 = theta0, tau.min = tau.min, k = k){
  alt.methods.sum <- list()
  
  resultlist=foreach(j=1:sim) %dopar% { 
  
  set.seed(j*100)
  pvalue=matrix(rep(1,2*p*B),ncol=2*p)###initial pavalue####
  adj.p=rep(0,p)
  
  ####### Data Generation ##########
  
  # NORMAL Biomarkers
  M_alt=M=rmvnorm(n,mean=rep(0,p),sigma=(diag(1,p,p)+(1-diag(1,p,p))*rho))
  colnames(M_alt)=colnames(M)=paste("X",1:p,sep="")
  
  # Treatment indicators
  H=runif(n)<t.prob
  
  # Interaction Biomarkers 
  MH=M*H
  colnames(MH)=paste(colnames(M),'H',sep='')
  
  M=cbind(H,M,MH)
  index=c(1,seq(2,p+1,1),seq(2,p+1,1))
  M=M[,unlist(lapply(1:(p+1),function(x) which(index==x)))]
  
  active.prog<- paste('X', active1, sep='') 
  active.pred<- paste('X', active2, 'H', sep='')
  colnames(pvalue)=colnames(M)[-1]
  
  # True Hazards
  lambda=h0 * exp(alpha0*H + sum(M[,active.prog] %*% as.matrix(alpha)) + sum( M[,active.pred] %*% as.matrix(gamma)))
  u=runif(n)
  surv.T=-log(u)/lambda
  censor.T=runif(n,min=0,max=theta0)
  y=pmin(surv.T,censor.T)
  fail=(surv.T<=censor.T)
  id= 1:n
  mydata <-as.data.frame(cbind(M_alt, y,fail))
  
  
  ################################## Proposed Three-Stage Method #################################
  
  ####### Stage I: feature screening: we only do screening if p is larger than n/logn##########
  if(p<round(n/log(n))) {
    M <- M
  } else {
    beta_hat <- c()
    for(m in 1:p){
      fit.each<-coxph(Surv(y,fail)~mydata[,m],data=mydata)
      beta_hat <- c(beta_hat, as.numeric(fit.each$coef))
    }
    beta_hat_abs <- abs(beta_hat)
    beta_order <- order(beta_hat_abs, decreasing = TRUE)
    x.index <- beta_order[1:round(n/log(n))]
    #names(mydata)[x.index]
    m.index <- which(colnames(M) %in% names(mydata)[x.index])
    M.select <- M[,c(1,m.index, m.index+1)]

    M <-  M.select
  }

  ####### Stage II: variable selection and p-value adjustment ##########
  
  b=1
  repeat{
    #######split data into training and test######
    chosen<-sample(1:n,size=k*n)  ###index for training###
    nchosen<-which((1:n %in% chosen)==0) ###index for testing###
    train.y<-t(t(y[chosen]))
    train.fail<-t(t(fail[chosen]))
    train.X=cbind(M[chosen,])
    colnames(train.y)="time";colnames(train.fail)="status"
    test.y<-t(t(y[nchosen]))
    test.fail<-t(t(fail[nchosen]))
    colnames(test.y)="time";colnames(test.fail)="status"
    test.X<-M[nchosen,]
    data = list(x = as.matrix(train.X), time = train.y, status=train.fail)


    #######group lasso fit######
    group_temp <- ifelse(p<round(n/log(n)), p, round(n/log(n)))
    cvFit=tryCatch(cv.grpsurv(as.matrix(train.X), Surv(train.y,train.fail), group=c(0,rep(2:(group_temp+1),times=2))),error=function(e){0},warning=function(w){0})
    if(is.list(cvFit)==0) {
      next
    }
    est_pred=coef(cvFit)[-1]#leave out H, do not perform test on H
    param=colnames(train.X[,-1])
    active.predictor=param[est_pred!=0]


    #######ML fit######
    if(length(active.predictor)){
      test.X=data.frame(test.X[,c("H",active.predictor)])
      ml.fit<-tryCatch(coxph(Surv(test.y,test.fail)~., data=test.X),error=function(e){0},warning=function(w){0})
      if(is.list(ml.fit)==0) {
        next
      }
      if(length(which(is.na(summary(ml.fit)[[7]][,5][-1])))) next
      pvalue[b,active.predictor]=summary(ml.fit)[[7]][,5][-1]
    }
    b=b+1
    if(b>B) break
  }
  #######Calculate adjusted p-value######
  tau.min=0.05
  adj.p=apply(pvalue,2,adaptive)

  #######false positives and true positives######
  activeset=names(adj.p)[which(adj.p<=0.05)]
  t1e=ifelse(min(activeset %in% c(active.prog, active.pred))==0,1,0)

  ####### Stage III  #######
  
  if(length(activeset)>0){

    activemh=activeset[grep("H", activeset)] ###find predictive markers###
    if(length(activemh)){
      activem=unlist(strsplit(activemh,"H")) ###insert main effect for the predictive ones##
      activeset=unique(c(activem,activeset))
    }
    alldata=cbind(y,fail)
    colnames(alldata)=c("time","status")
    alldata=data.frame(cbind(alldata,M[,c("H",activeset)]))


    final.fit<-tryCatch(coxph(Surv(time,status)~.,data=alldata),error=function(e){0},warning=function(w){0})
    if(is.list(final.fit)==0) {
      est.parm=NA
      parmname=activeset
      se.parm=NA
      est.tb=cbind(data.frame(cbind(est=est.parm,se=se.parm),row.names=NULL),select=parmname, sim=j, censor=mean(1-fail),t1e=t1e)
    } else{
      est.parm=coef(final.fit)
      parmname=names(est.parm)
      se.parm=sqrt(diag(final.fit$var))
      est.tb=cbind(data.frame(cbind(est=est.parm,se=se.parm),row.names=NULL),select=parmname, sim=j, censor=mean(1-fail),t1e=t1e)
    }
  } else {
    est.tb<-cbind(data.frame(cbind(est=NA,se=NA),row.names=NULL),select=NA, sim=j, censor=mean(1-fail),t1e)
  }
  ################################## Part of Proposed Method Ends #################################
  
  
  ################################## Existing Methods #################################
  
  data<- cbind(id, H, y, fail, M_alt)
  data<- as.data.frame(data)
  colnames(data)[2:4]<- c('treat', 'time', 'status')
  data[,"off.clin"] <- rep(0, nrow(data))
  nfold <- 3
  THRESH <- 1e-06
  set.seed(j*100)
  foldid <- sample(1:nfold, size = nrow(data), replace = TRUE)
  simulation=TRUE
  
  ##Penalized regression: full-lasso;
  fulllasso<- fullLasso(data=data, foldid = foldid, sim = simulation)
  ###check the number of selected biomarkers that are false positive;
  t1e.fulllasso=ifelse(min(names(fulllasso)[-1] %in% c(active.prog, active.pred))==0,1,0)
 

  ###Gradient boosting;
  boost<- gboost(data=data, foldid = foldid, sim=simulation)
  t1e.boost=ifelse(min(names(boost) %in% c(active.prog, active.pred))==0,1,0)
  
  ###Penalized regression: adaptive lasso (grouped weights);
  alassoGw<- gwaLasso(data=data, foldid = foldid, sim=simulation)
  t1e.adapativelasso=ifelse(min(names(alassoGw)[-1] %in% c(active.prog, active.pred))==0,1,0)

  ###Dimension reduction: PCA + lasso;
  pcalasso<- pcaLasso(data=data, foldid = foldid, sim=simulation)
  t1e.pcalasso=ifelse(min(names(pcalasso) %in% c(active.prog, active.pred))==0,1,0)

  ###Group Lasso
  grouplasso<- gLasso(data=data, foldid = foldid)
  t1e.grouplasso<- ifelse(min(names(grouplasso)%in% c(active.prog, active.pred))==0,1,0)

  ###results summary
   alt.methods.sum[[j]] <- list(fulllasso, boost, alassoGw, pcalasso, grouplasso)
   alt.methods<- cbind(t1e.new=t1e, t1e.fulllasso=t1e.fulllasso, t1e.boost=t1e.boost, t1e.adapativelassoo=t1e.adapativelasso, t1e.pcalasso=t1e.pcalasso, t1e.grouplasso=t1e.grouplasso, sim=j)
   
   ################################## Part of Existing Methods End #################################
   
   return(list(est.tb, alt.methods,alt.methods.sum)) 
}

return (resultlist)
} 

