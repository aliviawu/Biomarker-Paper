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
# - k          the allocation rate for multi-splitting                  #
#########################################################################
runsims<- function(n = n, p = p, sim = sim, B = B, t.prob = t.prob, rho = rho, h0 = h0,
                   alpha0 = alpha0, active1 = active1, active2 = active2, alpha = alpha, 
                   gamma = gamma, theta0 = theta0, tau.min = tau.min, k = k){
  alt.methods.sum <- list()
  
  resultlist=foreach(j=1:sim) %dopar% { 
  
    set.seed(j*100)
    pvalue=matrix(rep(1,p*B),ncol=p)###initial pavalue####
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
  colnames(pvalue)=colnames(M_alt)  
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
  
 
  
  ####### Stage II: variable selection and p-value adjustment ##########
  
  b=1
  a=1
  repeat{
    #######split data into training and test######
    a<- a+1
    set.seed(a*100+i+1000)
    
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

    ####### Stage I: feature screening using training dataset ##########
    if(p<round(n/log(n))) {
      train.X <- train.X[,c(1, sapply(1:p, function(x) 2*x), sapply(1:p, function(x) 2*x+1))]
    } else {
      train.X<- as.data.frame(train.X)
      pair_pvalue <- c()
      pair_chiq_stats<- c()
      fit_NULL<- coxph(Surv(train.y, train.fail)~train.X[,1], data = train.X)
      for(m in 1:p){
        fit.each<-coxph(Surv(train.y,train.fail)~train.X[,1]+train.X[,2*m]+train.X[,2*m+1],data=train.X)
        joint_test<- anova(fit_NULL, fit.each)
        pair_pvalue <- c(pair_pvalue, joint_test$`P(>|Chi|)`[2])
        pair_chiq_stats<-c(pair_chiq_stats, joint_test$Chisq[2])
      }
      
      #######here, we only do screening if p is larger than n/logn##########
      
      chiq_stats_order<- order(pair_chiq_stats, decreasing = TRUE)
      x.index<- chiq_stats_order[1:round(n/log(n))]
      M.select<- train.X[,c(1, 2*x.index, 2*x.index+1)]
      
      train.X <-  M.select
    }
    #######group lasso fit######
    group_temp <- ifelse(p<round(n/log(n)), p, length(x.index))
    cv.ridgefullmodel <- cv.glmnet(
      x = as.matrix(train.X),
      y = cbind(time=train.y, status=train.fail),
      family = "cox",
      alpha = 0,
      grouped = TRUE,
      standardize = FALSE,
      penalty.factor = c(0, rep(1, dim(train.X)[2]-1)))
    fit.ridgefullmodel <- glmnet(
      x = as.matrix(train.X),
      y = cbind(time=train.y, status=train.fail),
      family = "cox",
      alpha = 0,
      lambda = cv.ridgefullmodel$lambda.min,
      standardize = FALSE,
      penalty.factor=c(0, rep(1, dim(train.X)[2]-1)))
    coef.biombioi <- as.vector(coef(fit.ridgefullmodel))[-1]
    group=c(rep(1:(group_temp),times=2))
    w<- sapply(1:group_temp, function(x) sqrt(sum(coef.biombioi[which(group==x)]^2))^(-1))
    train.X<- train.X[,c(1,as.vector(sapply(2:(group_temp+1), function(x) c(x, x+group_temp))))]
    wi<- c(w)
    # g<- as.factor(c(0,rep(1:(group_temp+1),each=2)))
    g<- c(0,rep(1:(group_temp),each=2))
    cvFit=tryCatch(cv.grpsurv(as.matrix(train.X), Surv(train.y,train.fail), group=g, group.multiplier=wi),error=function(e){0},warning=function(w){0})
    if(is.list(cvFit)==0) {
      next
    }
    
    est_pred=coef(cvFit)[-1]#leave out H, do not perform test on H
    param=colnames(train.X[,-1])
    active.predictor=param[est_pred!=0]

    num.ap<-length(active.predictor)/2 ####|S|#####
    
    #######ML fit######
    if(length(active.predictor)){
      test.X=data.frame(test.X[,c("H",active.predictor)])
      alt.test<- tryCatch(coxph(Surv(test.y, test.fail)~., data=test.X), error=function(e){0}, warning=function(w){0})
      if(is.list(alt.test)==0){
        next
      }
      for(g in 1:num.ap){
        test.X.pair<- test.X[,-c(2*g, 2*g+1)]
        null.test<- tryCatch(coxph(Surv(test.y, test.fail)~., data=test.X.pair), error=function(e){0}, warning=function(w){0})
        if(is.list(null.test)==0){
          next
        }
        lr.test<- anova(null.test, alt.test)
        pvalue[b, active.predictor[2*k-1]]<- min(lr.test$`P(>|Chi|)`[2]*num.ap,1) ### adjust pvalue
      }
      
    }
    b=b+1
    if(b>B) break
  }
  #######Calculate adjusted p-value######
  tau.min=0.05
  adj.p=apply(pvalue,2,adaptive)

  #######false positives and true positives######
  activeset.pairs=names(adj.p)[which(adj.p<=0.05)]

  ####### Stage III  #######
  
  if(length(activeset.pairs)>0){

    activeset<- c(activeset.pairs, paste(activeset.pairs, 'H', sep=''))
    
    alldata=cbind(y,fail)
    colnames(alldata)=c("time","status")
    alldata=data.frame(cbind(alldata,M[,c("H",activeset)]))


    final.fit<-tryCatch(coxph(Surv(time,status)~.,data=alldata),error=function(e){0},warning=function(w){0})
    if(is.list(final.fit)==0) {
      est.parm=NA
      parmname=activeset
      se.parm=NA
      t1e<- NA
      est.tb=cbind(data.frame(cbind(est=est.parm,se=se.parm),row.names=NULL),select=parmname, sim=j, censor=mean(1-fail), t1e)
    } else{
      pv<- summary(final.fit)[[7]][,5][-1]
      sel.parm<- pv[pv<=0.05]
      est.parm=coef(final.fit)[names(sel.parm)]
      parmname=names(est.parm)
      index<- which(pv<=0.05)
      se.parm=sqrt(diag(final.fit$var))[index+1]
      t1e=ifelse(min(parmname %in% c(active.prog, active.pred))==0,1,0)
      
      if(length(se.parm)!=0){
        est.tb=cbind(data.frame(cbind(est=est.parm,se=se.parm, pvalue=as.numeric(pv[parmname])),row.names=NULL),select=parmname, sim=j, censor=mean(1-fail),t1e=t1e)
      }else{est.tb<-cbind(data.frame(cbind(est=NA,se=NA),row.names=NULL),select=NA, sim=j, censor=mean(1-fail),t1e, fdr=NA)}
      
      
    }
  } else {
    t1e<- NA
    est.tb<-cbind(data.frame(cbind(est=NA,se=NA),row.names=NULL),select=NA, sim=j, censor=mean(1-fail),t1e)
  }
  ################################## End of Proposed Method #################################
  
  
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

