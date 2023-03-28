############### functions for existing methods##################################
################################################################################

################################################################################
### Penalized regression: full-lasso                                           #
################################################################################



fullLasso <- function(data, foldid, sim) {
  ### Data preparation
  pos.biom <- function(data) {
    which(!colnames(data) %in% c('off.clin', 'treat', 'treat2', 'status', 'time','id'))
  }
  
  lmax <- function(x, y, w, off = 0){
    if(length(off) == 1) off <- rep(0, nrow(x))
    glmnet(x = x,
           y = y,
           penalty.factor = w,
           family = "cox",
           alpha = 1,
           offset = off,
           standardize = FALSE,
           dfmax = sum(w == 0))$lambda[1]/0.7
    
  }
  
  Y <- cbind(time = data$time, 
             status = data$status)
  pos.biom <- pos.biom(data)
  data <- cbind(data, data[, pos.biom])
  pos.bioi <- ncol(data) + 1 - length(pos.biom):1
  data[, pos.bioi] <- data[, pos.biom] * data[, 'treat']
  colnames(data)[pos.bioi] <- paste0(colnames(data)[pos.bioi], 'H')
  pos.trea <- which(colnames(data) == 'treat')  
  
  nfold <- length(unique(foldid))
  THRESH <- 1e-06
  
  ### Lambda.max estimation
  w <- c(0, rep(1, 2 * length(pos.biom)))
  lmax <- lmax(
    x = as.matrix(data[, c(pos.trea, pos.biom, pos.bioi)]),
    y = Y,
    off = data[, "off.clin"],
    w = w)
  epsi <- ifelse(nrow(data) > length(w), 0.001, 0.01)
  lmin <- epsi * lmax
  
  ### Cross-validation to estimate the optimal lambda
  cv.flasso <- cv.glmnet(
    x = as.matrix(data[, c(pos.trea, pos.biom, pos.bioi)]),
    y = Y,
    family = "cox",
    alpha = 1,
    lambda=seq(lmin, lmax, epsi),
    foldid = foldid,
    grouped = TRUE,
    standardize = FALSE,
    thresh = THRESH,
    offset = data[, "off.clin"],
    penalty.factor = w)
  
  ### Fit of the final model
  fit.flasso <- glmnet(
    x = as.matrix(data[, c(pos.trea, pos.biom, pos.bioi)]),
    y = Y,
    family = "cox",
    alpha = 1,
    lambda = cv.flasso$lambda.min * ((nfold - 1) / nfold),
    standardize = FALSE,
    thresh = THRESH,
    offset = data[, "off.clin"],
    penalty.factor = w)
  
  ncoef.flasso <- coef(fit.flasso)@Dimnames[[1]][which(coef(fit.flasso)!= 0)]
  coef.flasso <- coef(fit.flasso)[which(coef(fit.flasso) != 0)]
  names(coef.flasso) <- ncoef.flasso
  
  return(coef.flasso)
}


################################################################################
### Gradient boosting                                                          #
################################################################################

gboost <- function(data, foldid, sim){
  
  # Data preparation
  pos.biom <- function(data) {
    which(!colnames(data) %in% c('off.clin', 'treat', 'treat2', 'status', 'time','id'))
  }
  
  Y <- cbind(time = data$time, 
             status = data$status)
  pos.biom <- pos.biom(data)
  data <- cbind(data, data[, pos.biom])
  pos.bioi <- ncol(data) + 1 - length(pos.biom):1
  data[, pos.bioi] <- data[, pos.biom] * data[, 'treat']
  colnames(data)[pos.bioi] <- paste0(colnames(data)[pos.bioi], 'H')
  pos.trea <- which(colnames(data) == 'treat')
  nfold <- length(unique(foldid))
  THRESH <- 1e-06
  
  # Treatment effect (offset)
  off.trt <- data[, "treat"] * coxph(Surv(time, status) ~ treat, data = data)$coefficients
  
  # Gradient boosting
  boost <- glmboost(
    x = as.matrix(cbind(1, data[, c(pos.biom, pos.bioi)])),
    y = Surv(data$time,
             data$status),
    offset = off.trt + data[, 'off.clin'],
    family = CoxPH(),
    center = T,
    control = boost_control(mstop = 2 * length(pos.biom)))
  
  # Cross-validation to find the optimal number of steps
  cv.boo <- cv(
    weights = model.weights(boost),
    type = "kfold", 
    B = nfold)
  
  cv.boost <- cvrisk(
    object = boost,
    folds = cv.boo,
    papply = lapply)
  
  
  
  # Fit of the final model
  fit <- boost[mstop(cv.boost)]
  coef.boost <- coef(fit)[-1]
  
  return(coef.boost) 
}


################################################################################
### Penalized regression: adaptive lasso (grouped weights)                     #
################################################################################
gwaLasso <- function(data, foldid, sim) {
  ### Data preparation
  pos.biom <- function(data) {
    which(!colnames(data) %in% c('off.clin', 'treat', 'treat2', 'status', 'time','id'))
  }
  
  lmax <- function(x, y, w, off = 0){
    if(length(off) == 1) off <- rep(0, nrow(x))
    glmnet(x = x,
           y = y,
           penalty.factor = w,
           family = "cox",
           alpha = 1,
           offset = off,
           standardize = FALSE,
           dfmax = sum(w == 0))$lambda[1]/0.7
  }
  Y <- cbind(time = data$time, 
             status = data$status)
  pos.biom <- pos.biom(data)
  data <- cbind(data, data[, pos.biom])
  pos.bioi <- ncol(data) + 1 - length(pos.biom):1
  data[, pos.bioi] <- data[, pos.biom] * data[, 'treat']
  colnames(data)[pos.bioi] <- paste0(colnames(data)[pos.bioi], 'H')
  pos.trea <- which(colnames(data) == 'treat')
  nfold <- length(unique(foldid))
  THRESH <- 1e-06
  
  ##### Preliminary step (ridge penalty on the full biomarker-by-treatment interaction model)
  ### Cross-validation to estimate the optimal lambda  
  cv.ridgefullmodel <- cv.glmnet(
    x = as.matrix(data[, c(pos.trea, pos.biom, pos.bioi)]),
    y = Y,
    family = "cox",
    alpha = 0,
    foldid = foldid,
    grouped = TRUE,
    standardize = FALSE,
    thresh = THRESH,
    offset = data[, "off.clin"],
    penalty.factor = c(0, rep(1, 2 * length(pos.biom))))
  ### Fit of the model
  fit.ridgefullmodel <- glmnet(
    x = as.matrix(data[, c(pos.trea, pos.biom, pos.bioi)]),
    y = Y,
    family = "cox",
    alpha = 0,
    lambda = cv.ridgefullmodel$lambda.min * ((nfold-1)/nfold),
    standardize = FALSE,
    thresh = THRESH,
    offset = data[, "off.clin"],
    penalty.factor=c(0, rep(1, 2 * length(pos.biom))))
  coef.biombioi <- as.vector(coef(fit.ridgefullmodel))[-1]
  ### Grouped weights (Gw)
  weights <- 1 / abs(coef.biombioi)
  wRbm <- mean(abs(coef.biombioi[1:length(pos.biom)]))
  WRbi <- mean(abs(coef.biombioi[length(pos.biom)+(1:length(pos.biom))]))
  wRgr <- c(rep(1/wRbm, length(pos.biom)), rep(1 / WRbi, length(pos.biom)))
  wRgr <- wRgr / mean(wRgr)
  #####
  
  ### Resampling for new fold assignments
  foldid <- sample(foldid)
  
  ### Lambda.max estimation
  w <- c(0,wRgr)
  lmax <- lmax(
    x = as.matrix(data[, c(pos.trea, pos.biom, pos.bioi)]),
    y = Y,
    off = data[, "off.clin"],
    w = w)
  epsi <- ifelse(nrow(data) > length(w), 0.0001, 0.01)
  lmin <- epsi * lmax
  
  ### Cross-validation to estimate the optimal lambda
  cv.alassoGw <- cv.glmnet(
    x = as.matrix(data[, c(pos.trea, pos.biom, pos.bioi)]),
    y = Y,
    family = "cox",
    alpha = 1,
    lambda = seq(lmin, lmax, epsi),
    foldid = foldid,
    grouped = TRUE,
    standardize = FALSE,
    thresh = THRESH,
    offset = data[, "off.clin"],
    penalty.factor = w)
  
  
  
  ### Fit of the final model
  fit.alassoGw <- glmnet(
    x = as.matrix(data[, c(pos.trea, pos.biom, pos.bioi)]),
    y = Y,
    family = "cox",
    alpha = 1,
    lambda = cv.alassoGw$lambda.min * ((nfold - 1) / nfold),
    standardize = FALSE,
    thresh = THRESH,
    offset = data[, "off.clin"],
    penalty.factor = w)
  
  ncoef.alassoGw <- coef(fit.alassoGw)@Dimnames[[1]][which(coef(fit.alassoGw) != 0)]
  coef.alassoGw <- coef(fit.alassoGw)[which(coef(fit.alassoGw) != 0)]
  names(coef.alassoGw) <- ncoef.alassoGw
  
  return(coef.alassoGw)
  
}

################################################################################
### Dimension reduction: PCA + lasso                                           #
################################################################################
pcaLasso <- function(data, foldid, sim) {
  ### Data preparation
  pos.biom <- function(data) {
    which(!colnames(data) %in% c('off.clin', 'treat', 'treat2', 'status', 'time','id'))
  }
  lmax <- function(x, y, w, off = 0){
    if(length(off) == 1) off <- rep(0, nrow(x))
    glmnet(x = x,
           y = y,
           penalty.factor = w,
           family = "cox",
           alpha = 1,
           offset = off,
           standardize = FALSE,
           dfmax = sum(w == 0))$lambda[1]/0.7
  }  
  Y <- cbind(time = data$time, 
             status = data$status)
  pos.biom <- pos.biom(data)
  data <- cbind(data, data[, pos.biom])
  pos.bioi <- ncol(data) + 1 - length(pos.biom):1
  data[, pos.bioi] <- data[, pos.biom] * data[, 'treat']
  colnames(data)[pos.bioi] <- paste0(colnames(data)[pos.bioi], 'H')
  pos.trea <- which(colnames(data) == 'treat')
  nfold <- length(unique(foldid))
  THRESH <- 1e-06
  
  ##### Preliminary step (dimension reduction of main effects matrix via PCA)
  ### Estimation of the ncompS first principal components
  ncompS <- round(min(length(pos.biom),
                      table(data[which(data$status == 1), "treat"])))
  pcS <- (fast.svd(t(t(data[, pos.biom]) - colMeans(data[, pos.biom]))))
  PCS <- as.matrix(
    pcS$u[, 1:ncompS, drop=FALSE] %*% diag(pcS$d[1:ncompS], ncompS))
  colnames(PCS) <- gsub(" ", "", format(paste0("PCS", 1:ncompS)))
  ### Cross-validation to estimate the optimal number K
  cv.K <- t(sapply(1:max(foldid),FUN = function(i) {
    PCS.T <- PCS[which(foldid != i),]
    sapply(X = 1:ncompS, FUN = function(X){  
      fit <- coxph(Surv(time, status) ~ PCS.T[,1:X],
                   data = data[which(foldid != i),])
      coef.T <- fit$coefficients
      nlogl.T <- fit$loglik[2]
      off.All <- PCS[, 1:X] %*% as.matrix(coef.T) 
      nlogl.All <- coxph(Surv(time, status) ~ offset(off.All),
                         data = data)$loglik
      nlogl <- nlogl.All - nlogl.T
      return(nlogl)
    })
  }))
  K.PCA <- which.max(colSums(cv.K))
  data.PCA <- cbind(data,PCS[, 1:K.PCA])
  pos.PCS <- grep("PCS", names(data.PCA))
  #####
  
  ### Resampling for new fold assignments
  foldid <- sample(foldid)
  
  ### Lambda.max estimation
  w <- c(rep(0, length(pos.PCS) + 1), rep(1, length(pos.bioi)))
  lmax <- lmax(
    x = as.matrix(data.PCA[, c(pos.trea, pos.PCS, pos.bioi)]),
    y = Y,
    w = w,
    off = data[,"off.clin"])
  epsi <- ifelse(nrow(data) > length(w), 0.0001, 0.01)
  lmin <- epsi * lmax
  
  ### Cross-validation to estimate the optimal lambda
  cv.PCA <- cv.glmnet(
    x = as.matrix(data.PCA[, c(pos.trea, pos.PCS, pos.bioi)]),
    y = Y,
    family = "cox",
    alpha = 1,
    lambda = seq(lmin, lmax, epsi),
    foldid = foldid,
    standardize = FALSE,
    thresh = THRESH,
    offset = data[, "off.clin"],
    penalty.factor = w)
  
  
  ### Fit of the final model
  fit.PCA <- glmnet(
    x = as.matrix(data.PCA[, c(pos.trea, pos.PCS, pos.bioi)]),
    y = Y,
    family = "cox",
    alpha = 1,
    lambda = cv.PCA$lambda.min * ((nfold - 1) / nfold),
    standardize = FALSE,
    thresh = THRESH,
    offset = data[, "off.clin"],
    penalty.factor = w)
  
  coef.PCA <- coef(fit.PCA)[
    which(coef(fit.PCA)!=0)][-c(2:(K.PCA+1))]
  names(coef.PCA) <- coef(fit.PCA)@Dimnames[[1]][
    which(coef(fit.PCA)!=0)][-c(2:(K.PCA+1))]
  #names(coef.PCA) <- gsub("[.]I", "", names(coef.PCA))
  
  if(length(names(coef.PCA))>1){
    biomarker<-  gsub("H", "", names(coef.PCA)[-1])
    data.sub<- data %>% dplyr::select(time, status, treat, biomarker, names(coef.PCA)[-1])
    fit.final<- coxph(Surv(time, status)~., data=data.sub)
    pvalue<- summary(fit.final)[[7]][,5][-1]
    coef.PCA<- summary(fit.final)[[7]][,1][-1][pvalue<=0.05]
    
  }
  return(coef.PCA)
}


################################################################################
### Penalized regression: group-lasso                                          #
################################################################################
gLasso <- function(data, foldid, sim) {
  ### Data preparation
  pos.biom <- function(data) {
    which(!colnames(data) %in% c('off.clin', 'treat', 'treat2', 'status', 'time','id'))
  }
  
  Y <- cbind(time = data$time, 
             status = data$status)
  id <- 1:nrow(data)
  #data <- cbind(id,data)
  pos.biom <- pos.biom(data)
  data <- cbind(data, data[, pos.biom])
  pos.bioi <- ncol(data) + 1 - length(pos.biom):1
  data[, pos.bioi] <- data[, pos.biom] * data[, 'treat']
  colnames(data)[pos.bioi] <- paste0(colnames(data)[pos.bioi], 'H')
  pos.trea <- which(colnames(data) == 'treat')    
  data <- cbind(data,foldid)
  nfold <- length(unique(foldid))
  THRESH <- 1e-06
  
  
  group <- c(0, rep(1:length(pos.biom), 2))
  
  ### Cross-validation to estimate the optimal lambda
  
  cv.gLasso<- cv.grpsurv(X=as.matrix(data[, c(pos.trea, pos.biom, pos.bioi)]),
                         y=Y, group = group, penalty='grLasso', 
                         alpha=1,fold=foldid)
  
  
  
  coef.gLasso<- coef(cv.gLasso)[coef(cv.gLasso)!=0]
  if(length(coef.gLasso)==0){
    coef.grplasso=NA
  }else{
    
    data.sub<- data %>% dplyr::select(time, status, names(coef.gLasso))
    fit.final<- coxph(Surv(time, status)~., data=data.sub)
    pvalue<- summary(fit.final)[[7]][,5][-1]
    coef.grplasso<- summary(fit.final)[[7]][,1][-1][pvalue<=0.05]
  }
  return(coef.grplasso)
}
