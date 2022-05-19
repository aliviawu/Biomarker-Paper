require(survival)
require(SGL)
require(mvtnorm)
require(foreach)
require(grpreg)
require(glmnet)
require(robustbase)
require(mboost)
require(doMC)
library(corpcor)
library(dplyr)
registerDoMC(cores=6)
demo<- read.table('./FinalCodesResults/real_application/application2/data/GSE6532_LUMINAL_demo.txt',header = TRUE)
#R.utils::gunzip("/Users/aliviawu/Downloads/GSE6532_LUMINAL.RData.gz") 
load("./FinalCodesResults/real_application/application2/data/GSE6532_LUMINAL.RData")
adaptive<-function(p){
  if(sum(p<0.05)==0){adj.p=1}
  if(sum(p<0.05)!=0){
    tau=seq(0.05,1,0.01)
    Q=rep(0,length(tau))
    q<- rep(0, length(tau))
    for (tau.i in 1:length(tau)) {
      Q[tau.i]=min(1,quantile(p/tau[tau.i],prob=tau[tau.i],na.rm=T))
    }
    adj.p=min(1,(1-log(0.05))*min(Q))
  }
  return (adj.p)
}

########## Data Processing #############
data_rfs<- demo[!is.na(demo$t.rfs)&!is.na(demo$e.rfs),]
rfs.tam.idx<- match(data_rfs$samplename, rownames(data.tam))
rfs.tam.sub<- data.tam[rfs.tam.idx[!is.na(rfs.tam.idx)],] 
rfs.tam.sub<- as.data.frame(rfs.tam.sub)
rfs.tam.sub$samplename<- rownames(rfs.tam.sub)

rfs.tam.comb<- merge(rfs.tam.sub, data_rfs[data_rfs$treatment=='tamoxifen',c(1,12,13,16)], by='samplename')

untam.idx<- match(data_rfs$samplename, rownames(data.untreated))
rfs.untam.sub<- data.untreated[untam.idx[!is.na(untam.idx)],]
rfs.untam.sub<- as.data.frame(rfs.untam.sub)
rfs.untam.sub$samplename<- rownames(rfs.untam.sub)
rfs.untam.comb<- merge(rfs.untam.sub, data_rfs[data_rfs$treatment=='none',c(1,12,13,16)], by='samplename')
Breast.rfs<- rbind(rfs.tam.comb, rfs.untam.comb) 
colnames(Breast.rfs)[44932]<- 'treat'
colnames(Breast.rfs)[44930]<- 'time'
colnames(Breast.rfs)[44931]<- 'status'
colnames(Breast.rfs)[1]<- 'id'

check2<- lapply(1:87, function(x) which(is.na(Breast.rfs[x,])))
Breast.rfs.sub<- Breast.rfs[,-c(check2[[1]])]
Breast.rfs.sub$treat<- ifelse(Breast.rfs.sub$treat=='none', 0, 1) 


B=50
tau.min=0.05        ####para q for adaptive adjusting pvalue#####

n<- nrow(Breast.rfs.sub) #380 patients
p<- ncol(Breast.rfs.sub)-4 #44916 biomarkers
pvalue=matrix(rep(1,p*B),ncol=p)###initial pavalue####

adj.p<- rep(0,p)
censor.rate<- mean(1-Breast.rfs.sub$status) #0.6463104
#Breast$treat<- ifelse(Breast$treat==-0.5, 0, 1)
M_alt<- Breast.rfs.sub[,c(2:(p+1))]
H<- Breast.rfs.sub$treat

MH=M_alt*H
colnames(MH)=paste(colnames(M_alt),"H",sep="")
M=cbind(H,M_alt,MH)
index=c(1,seq(2,p+1,1),seq(2,p+1,1))
M=M[,unlist(lapply(1:(p+1),function(x) which(index==x)))]
colnames(pvalue)=colnames(M_alt)
y<- Breast.rfs.sub$time
fail<- Breast.rfs.sub$status
mydata <-as.data.frame(cbind(M_alt, y,fail))

B=50
b=1

j<- 1
active.group<- list()
repeat{
  #######split data into training and test######
  j<- j+1
  set.seed(j*100)
  events<- which(fail==1)
  chosen.events<- sample(events, size = 0.7*length(events))
  nonevents<- which(fail==0)
  chosen.nonevents<- sample(nonevents, size=length(nonevents)*0.7)
  # 
  #chosen<-sample(1:n,size=(1/3)*n)  ###index for training###
  #nchosen<-which((1:n %in% chosen)==0) ###index for testing###
  chosen<- c(chosen.events, chosen.nonevents)
  nchosen<-which((1:n %in% chosen)==0)
  
  #nchosen<- which((1:n %in% chosen)==0)
  train.y<-t(t(y[chosen]))
  train.fail<-t(t(fail[chosen]))
  train.X=cbind(M[chosen,])
  colnames(train.y)="time";colnames(train.fail)="status"
  test.y<-t(t(y[nchosen]))
  test.fail<-t(t(fail[nchosen]))
  colnames(test.y)="time";colnames(test.fail)="status"
  test.X<-M[nchosen,]
  data = list(x = as.matrix(train.X), time = train.y, status=train.fail)
  #######stage I: feature screening##########
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
    
    # x.index<- which(pair_pvalue<=0.05)  #####?????Xue, I think this one is too stringent. I think you can use and order chi-square statistics to get the top n/log(n) as we had before. 
    chiq_stats_order<- order(pair_chiq_stats, decreasing = TRUE)
    x.index<- chiq_stats_order[1:round(n/log(n))]
    M.select<- train.X[,c(1, 2*x.index, 2*x.index+1)]
    #M<- M[,c(1, x.index, )]
    #beta_hat_abs <- abs(beta_hat)
    #beta_order <- order(beta_hat_abs, decreasing = TRUE)
    #x.index <- beta_order[1:round(n/log(n))]
    #names(mydata)[x.index]
    #m.index <- which(colnames(M) %in% names(mydata)[x.index])
    #M.select <- M[,c(1,m.index, m.index+1)]
    
    #######stage II: here, we only do screening if p is larger than n/logn##########
    train.X <-  M.select
  }
  #########group lasso fit##############
  group_temp <- ifelse(p<round(n/log(n)), p, round(n/log(n)))
  cv.ridgefullmodel <- cv.glmnet(
    x = as.matrix(train.X),
    y = cbind(time=train.y, status=train.fail),
    family = "cox",
    alpha = 0,
    grouped = TRUE,
    standardize = FALSE,
    penalty.factor = c(0, rep(1, dim(train.X)[2]-1)))
  ### Fit of the model
  fit.ridgefullmodel <- glmnet(
    x = as.matrix(train.X),
    y = cbind(time=train.y, status=train.fail),
    family = "cox",
    alpha = 0,
    lambda = cv.ridgefullmodel$lambda.min,
    standardize = FALSE,
    penalty.factor=c(0, rep(1, dim(train.X)[2]-1)))
  coef.biombioi <- as.vector(coef(fit.ridgefullmodel))[-1]
  ### Grouped weights (Gw)
  #weights <- 1 / abs(coef.biombioi)
  #wRbm <- mean(abs(coef.biombioi[1:length(pos.biom)]))
  #WRbi <- mean(abs(coef.biombioi[length(pos.biom)+(1:length(pos.biom))]))
  #wRgr <- c(rep(1/wRbm, length(pos.biom)), rep(1 / WRbi, length(pos.biom)))
  #wRgr <- wRgr / mean(wRgr)
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
  active.group[[b]]<- active.predictor 
  print('step1')
  num.ap<-length(active.predictor)/2 ####|S|#####
  # coxph(Surv(time,status)~`205539_at`+`205539_at`*H,data=Breast.sub)))
  ##########ML fit###########
  if(length(active.predictor)){
    test.X=data.frame(test.X[,c("H",active.predictor)])
    alt.test<- tryCatch(coxph(Surv(test.y, test.fail)~., data=test.X), error=function(e){0}, warning=function(w){0})
    if(is.list(alt.test)==0){
      next
    }
    
    for(k in 1:num.ap){
      test.X.pair<- test.X[,-c(2*k, 2*k+1)]
      null.test<- tryCatch(coxph(Surv(test.y, test.fail)~., data=test.X.pair), error=function(e){0}, warning=function(w){0})
      if(is.list(null.test)==0){
        next
      }
      lr.test<- anova(null.test, alt.test)
      pvalue[b, active.predictor[2*k-1]]<- min(lr.test$`P(>|Chi|)`[2]*num.ap,1) ### adjust pvalue
    }
    
    # ml.fit<-tryCatch(coxph(Surv(test.y,test.fail)~., data=test.X),error=function(e){0},warning=function(w){0})
    #if(is.list(ml.fit)==0) {
    #  next
    #}
    #if(length(which(is.na(summary(ml.fit)[[7]][,5][-1])))) next
    #pvalue[b,active.predictor]=summary(ml.fit)[[7]][,5][-1] ###mw: first and third choice
  }
  b=b+1
  print(b)
  if(b>B) break
}
#########Calculate adjusted p-value########

adj.p=apply(pvalue,2,adaptive)
activeset.pairs=names(adj.p)[which(adj.p<=0.05)]
filename2<- paste('./FinalCodesResults/real_application/application2/results','activepairs_breast.RData',sep='') 
save(activeset.pairs, file=filename2)

