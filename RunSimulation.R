#########################################################################
# Generation of data sets and Run for the simulation study               #
#########################################################################
#                                                                        #
# The simulations are for our proposed method and the existing methods.  #
# The simulations should contains 30 studies:                            #
# (1): simulation with n = 150, 300, 500 samples                         #
# (2): For n = 150, simulation with p = 30, 50, 70, 100                  #
# (3): For n = 300, simulation with p = 30, 50, 70, 100, 150             #
# (4): For n = 500, simulation with p = 30, 50, 70, 100, 150, 250        #
# (5): For each pair of n and p, simulation with high censor rate of 60% #
# and low censor rate of 40%                                             #
#                                                                        #
# For each study, four scenarios for alpha0 (effect of treatment         #
# indicator H), alpha4 (effect of biomarker X4), alpha5 (effect of       #
# biomarker X5), and gamma5 (effect of interaction between X5 and        #
# treatment indicator H) have been implemented:                          #
# (S1): alpha0=1, alpha4=1, alpha5=1, gamma5=1                           #
# (S2): alpha0=1, alpha4=0.5, alpha5=1, gamma5=1                         #
# (S3): alpha0=1, alpha4=0.5, alpha5=1, gamma5=1.5                       #
# (S4): alpha0=1, alpha4=0.5, alpha5=1.5, gamma5=2.5                     #
#                                                                        #
#########################################################################
setwd('/Users/aliviawu/Downloads/FinalCodesResults')

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
source('./functions/Simulator.R')
source('./functions/ExistMethods.R')
source('./functions/AdjustPvalue.R')
####Parameters set-up; 

sim=10
B=50####repeat p-value 50 times#####

m.n<-rep(c(150,150,150,150,300,300,300,300,300,500,500,500,500,500,500), 2)
##the dimension of the biomarker;
m.p<-rep(c(30,50,70,100,30,50,70,100,150,30,50,70,100,150,250),2)
##baseline hazard which is set-up as a constant for simplicity;
m.h0<-rep(1/100,length(m.n))
##parameter for censoring time generation (uniform distributed); 
m.theta0<-rep(c(150, 270), each=length(m.n)/2)

###Choose the index of parameters for data simulation; 
for (i in 1 : 2){
  
  n = m.n[i]
  p = m.p[i]
  h0 = m.h0[i]
  theta0 = m.theta0[i]
  alpha0 = 1
  active1 = c(4, 5)
  active2 = c(4)
  rho = 0.15
  t.prob = 0.5
  tau.min = 0.05
  k = 0.5
  
  # Scenario 1: alpha0=1, alpha4=1, alpha5=1, gamma5=1
  alpha = c(1, 1)
  gamma = c(1)
  S1.results<- runsims(n = n, p = p, sim = sim, B = B, t.prob = t.prob, rho = rho, h0 = h0,
                       alpha0 = alpha0, active1 = active1, active2 = active2, alpha = alpha, 
                       gamma = gamma, theta0 = theta0, tau.min = tau.min, k = k)
  result.all=do.call("rbind",S1.results)
  result<- do.call('rbind',result.all[,1]) # results of our method
  result.alt.method<- do.call('rbind', result.all[,2]) # FWER of alternative methods
  result.alt.method.sum<- do.call('rbind', result.all[,3]) # all results of alternative methods
  filname=paste('./results/Sce1/',"results_",paste(c(alpha,gamma,n,p),collapse = "_"),".RData",sep="")
  save(result.all, result, result.alt.method, result.alt.method.sum, file = filname)
  
  # Scenario 2: alpha0=1, alpha4=0.5, alpha5=1, gamma5=1
  alpha = c(0.5, 1)
  gamma = c(1)
  S2.results<- runsims(n = n, p = p, sim = sim, B = B, t.prob = t.prob, rho = rho, h0 = h0,
                       alpha0 = alpha0, active1 = active1, active2 = active2, alpha = alpha, 
                       gamma = gamma, theta0 = theta0, tau.min = tau.min, k = k)
  result.all=do.call("rbind",S2.results)
  result<- do.call('rbind',result.all[,1]) # results of our method
  result.alt.method<- do.call('rbind', result.all[,2]) # FWER of alternative methods
  result.alt.method.sum<- do.call('rbind', result.all[,3]) # all results of alternative methods
  filname=paste('./results/Sce2/',"results_",paste(c(alpha,gamma,n,p),collapse = "_"),".RData",sep="")
  save(result.all, result, result.alt.method, result.alt.method.sum, file = filname)
  
  
  # Scenario 3: alpha0=1, alpha4=0.5, alpha5=1.5, gamma5=1
  alpha = c(0.5, 1)
  gamma = c(1.5)
  S3.results<- runsims(n = n, p = p, sim = sim, B = B, t.prob = t.prob, rho = rho, h0 = h0,
                       alpha0 = alpha0, active1 = active1, active2 = active2, alpha = alpha, 
                       gamma = gamma, theta0 = theta0, tau.min = tau.min, k = k)
  result.all=do.call("rbind",S3.results)
  result<- do.call('rbind',result.all[,1]) # results of our method
  result.alt.method<- do.call('rbind', result.all[,2]) # FWER of alternative methods
  result.alt.method.sum<- do.call('rbind', result.all[,3]) # all results of alternative methods
  filname=paste('./results/Sce3/',"results_",paste(c(alpha,gamma,n,p),collapse = "_"),".RData",sep="")
  save(result.all, result, result.alt.method, result.alt.method.sum, file = filname)

  
  # Scenario 4: alpha0=1, alpha4=0.5, alpha5=2.5, gamma5=1.5
  alpha = c(0.5, 1.5)
  gamma = c(2.5)
  S4.results<- runsims(n = n, p = p, sim = sim, B = B, t.prob = t.prob, rho = rho, h0 = h0,
                       alpha0 = alpha0, active1 = active1, active2 = active2, alpha = alpha, 
                       gamma = gamma, theta0 = theta0, tau.min = tau.min, k = k)
  result.all=do.call("rbind",S4.results)
  result<- do.call('rbind',result.all[,1]) # results of our method
  result.alt.method<- do.call('rbind', result.all[,2]) # FWER of alternative methods
  result.alt.method.sum<- do.call('rbind', result.all[,3]) # all results of alternative methods
  filname=paste('./results/Sce4/',"results_",paste(c(alpha,gamma,n,p),collapse = "_"),".RData",sep="")
  save(result.all, result, result.alt.method, result.alt.method.sum, file = filname)
  
}
  