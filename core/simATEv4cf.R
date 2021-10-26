library(glmnet)
library(ahaz)
source("core/phi.R")
source("core/psi.R")
source("core/phi-cf.R")
source("core/psi-cf.R")
source("core/pred-cf.R")

sim.ATE.v4cf = function(surv, D, Z, 
                        tau = ifelse(sum(surv[D==1,1]==max(surv[,1]))>2*nfold, 
                                     max(surv[,1]), 
                                     quantile(surv[D==1,1],probs = max(.9,1-2*nfold/sum(D)))),# bn = 20, Kg = 100, #
                        tol=1e-7,maxit=100,
                      true.ate,true.beta,true.gamma,true.lam0, oracle = TRUE, 
                      maxATE = log(2)*(.Machine$double.max.exp-1),
                      nfold = 10)
{
  n = nrow(Z)
  p = ncol(Z)
  
  # Nuisance parameters
  #======================
  # Lasso  
  cv.bhat.nopen = tune.ahazpen(surv,as.matrix(cbind(D,Z)),penalty="lasso",
                         penalty.wgt = c(0,rep(1,p)))
  
  full.bhat = ahazpen(surv,as.matrix(cbind(D,Z)),penalty="lasso")
  full.ghat = glmnet(Z,D,family = "binomial")

  theta.nopen =  as.numeric(
    predict(cv.bhat.nopen,type="coef",lambda = cv.bhat.nopen$lambda.min))[1]

  nlambda.beta = length(full.bhat$lambda)
  nlambda.gamma = length(full.ghat$lambda)
  neworder = c(sample(which(D==1&surv[,1]==tau)),
               sample(which(D==1&surv[,1]<tau)),
               sample(which(D==0)))
  foldid = rep(1:nfold,ceiling(n/nfold))[1:n]
  foldid[neworder] = foldid
  fit.fold <-foreach(k = 1:nfold, .packages=c("ahaz","glmnet","survival")) %dopar%
  {
    fitid = (foldid!=k) 
    Dk = D[fitid]
    Zk = Z[fitid,]
    survk = surv[fitid,]
    # Lasso
    fold.bhat = ahazpen(survk,as.matrix(cbind(Dk,Zk)),penalty="lasso", lambda = full.bhat$lambda)
    fold.ghat = glmnet(Zk,Dk,family = "binomial", lambda = full.ghat$lambda)
    list(fold.bhat = fold.bhat, 
         fold.ghat = fold.ghat)
  }
  
  fold1 = rep(1:nfold,nfold:1 -1)
  fold2 = unlist(lapply(2:10, ':',10))
  fit.fold.fold <-foreach(k = 1:(nfold*(nfold-1)/2), .packages=c("ahaz","glmnet","survival")) %dopar%
  {
    fitid = (!is.element(foldid,c(fold1[k],fold2[k]))) #& (foldid.tune!=k)
    Dk = D[fitid]
    Zk = Z[fitid,]
    survk = surv[fitid,]
    # Lasso
    fold.bhat = ahazpen(survk,as.matrix(cbind(Dk,Zk)),penalty="lasso", lambda = full.bhat$lambda)
    fold.ghat = glmnet(Zk,Dk,family = "binomial", lambda = full.ghat$lambda)
    list(fold.bhat = fold.bhat, 
         fold.ghat = fold.ghat)
  }
  
  
  tune.full <-foreach(k = 1:nfold, .packages=c("ahaz","glmnet","survival")) %dopar%
  {
    source("core/measure_v1.R")
    tuneid = (foldid==k)
    
    # cv for full fit
    #----------------------------------------------------------
    Dk = D[tuneid]
    Zk = Z[tuneid,]
    survk = surv[tuneid,]
    Ytauk = (survk[,1] >= tau)
    
    nz.beta = nonzeroCoef(fit.fold[[k]]$fold.bhat$beta)
    tmp.beta =predict(fit.fold[[k]]$fold.bhat, type = "coef",lambda = full.bhat$lambda)[nz.beta,]
    tmpZ = Z[tuneid,nz.beta-1]
    if(nz.beta[1]==1)
    {
      tmpZ = cbind(Dk,tmpZ)
    }
    tmpfit = ahaz(survk, tmpZ)
    dev.beta = apply(tmp.beta,2,ahaz.mse,x=tmpfit)
    mag.beta = apply(tmpZ %*% tmp.beta,2,mag.ahaz.v1, survk[,1])
    
  
    tmp.p =predict(fit.fold[[k]]$fold.ghat,newx = Zk, type = "response", s = full.ghat$lambda)
    dev.gamma = logistic.dev(Dk,tmp.p)
    mag.gamma = mag.logistic.v1(Dk,survk[,1],Ytauk,tmp.p)

    list(dev.beta = dev.beta, dev.gamma = dev.gamma,
         mag.beta = mag.beta, mag.gamma = mag.gamma)
  }
  
  tune.fold <-foreach(k = 1:nfold, .packages=c("ahaz","glmnet","survival")) %dopar%
  {
    cvlist = which(fold1==k | fold2==k)
    test.fold = setdiff(1:nfold,k)
    
    nlambda.fold.beta = length(fit.fold[[k]]$fold.bhat$lambda)
    dev.fold.beta = mag.fold.beta = matrix(0,nfold-1,nlambda.fold.beta)
    nlambda.fold.gamma = length(fit.fold[[k]]$fold.ghat$lambda)
    dev.fold.gamma = mag.fold.gamma = matrix(0,nfold-1,nlambda.fold.gamma)
    
    for(icv in 1:(nfold-1))
    {
      tuneid = foldid==test.fold[icv]
      Dk = D[tuneid]
      Zk = Z[tuneid,]
      survk = surv[tuneid,]
      Ytauk = (survk[,1] >= tau)
      
      
      nz.beta = nonzeroCoef(fit.fold.fold[[cvlist[icv]]]$fold.bhat$beta)
      tmp.beta =predict(fit.fold.fold[[cvlist[icv]]]$fold.bhat, type = "coef",
                        lambda = fit.fold[[k]]$fold.bhat$lambda)[nz.beta,]
      tmpZ = Z[tuneid,nz.beta-1]
      if(nz.beta[1]==1)
      {
        tmpZ = cbind(Dk,tmpZ)
      }
      tmpfit = ahaz(survk, tmpZ)
      dev.fold.beta[icv,] = apply(tmp.beta,2,ahaz.mse,x=tmpfit)
      mag.fold.beta[icv,] = apply(tmpZ %*% tmp.beta,2,mag.ahaz.v1, survk[,1])
      
      
      tmp.p =predict(fit.fold.fold[[cvlist[icv]]]$fold.ghat,newx = Zk, type = "response", 
                     s = fit.fold[[k]]$fold.ghat$lambda)
      dev.fold.gamma[icv,] = logistic.dev(Dk,tmp.p)
      mag.fold.gamma[icv,] = mag.logistic.v1(Dk,survk[,1],Ytauk,tmp.p)
    }
    
    cvm.dev.beta = apply(dev.fold.beta,2,mean)
    cvm.mag.beta = apply(mag.fold.beta,2,mean)
    cvm.dev.gamma = apply(dev.fold.gamma,2,mean)
    cvm.mag.gamma = apply(mag.fold.gamma,2,mean)
    list(
         hthetabeta = drop(predict(fit.fold[[k]]$fold.bhat,type="coef",
                              lambda = fit.fold[[k]]$fold.bhat$lambda[which.min(cvm.dev.beta)])),
         hgr = drop(predict(fit.fold[[k]]$fold.ghat, type = "coef", 
                            s = fit.fold[[k]]$fold.ghat$lambda[which.min(cvm.dev.gamma)])), 
         cvm.dev.beta = cvm.dev.beta,
         cvm.mag.beta = cvm.mag.beta,
         cvm.dev.gamma = cvm.dev.gamma,
         cvm.mag.gamma = cvm.mag.gamma)
  }

  cvm.dev.beta = apply(sapply(tune.full, function(x) x$dev.beta),1,mean)
  cvm.dev.gamma = apply(sapply(tune.full, function(x) x$dev.gamma),1,mean)
  cvm.mag.beta = apply(sapply(tune.full, function(x) x$mag.beta),1,mean)
  cvm.mag.gamma = apply(sapply(tune.full, function(x) x$mag.gamma),1,mean)
  
  # hat beta star
  hthetabeta =  as.numeric(
    predict(full.bhat,type="coef",lambda = full.bhat$lambda[which.min(cvm.dev.beta)]))
  
  # hat gamma
  hgr = as.numeric(
    predict(full.ghat,type="coef",s = full.ghat$lambda[which.min(cvm.dev.gamma)]))
  
  hps = 1/(1+exp(-hgr[1]-drop(Z%*%hgr[-1])))
  theta.ps = coef(ahaz(surv,cbind(D,hps)))[1]
  
  hthetabeta.cf = sapply(tune.fold, function(x) x$hthetabeta)  
  hgr.cf = sapply(tune.fold, function(x) x$hgr)
  
  fit = list(foldid = foldid,
             hthetabeta = hthetabeta, hgr=hgr,
             hthetabeta.cf = hthetabeta.cf, hgr.cf = hgr.cf, 
             fit.bhat = list(beta = full.bhat$beta, lambda = full.bhat$lambda, 
                             dev.cvm = cvm.dev.beta, mag.cvm = cvm.mag.beta), 
             fit.ghat = list(beta = rbind(full.ghat$a0,full.ghat$beta), 
                              lambda = full.ghat$lambda, 
                              dev.cvm = cvm.dev.gamma, mag.cvm = cvm.mag.gamma),
             fit.cf.bhat = vector("list",nfold), 
             fit.cf.ghat = vector("list",nfold),
             true.ate=true.ate,true.beta=true.beta,true.gamma=true.gamma,true.lam0=true.lam0,
             theta.nopen=theta.nopen, theta.ps = theta.ps)
  for(k in 1:nfold)
  {
    fit$fit.cf.bhat[[k]] = list(beta = fit.fold[[k]]$fold.bhat$beta, 
                                lambda = fit.fold[[k]]$fold.bhat$lambda,
                                dev.cvm = tune.fold[[k]]$cvm.dev.beta,
                                mag.cvm = tune.fold[[k]]$cvm.mag.beta)
    fit$fit.cf.ghat[[k]] = list(beta = rbind(fit.fold[[k]]$fold.ghat$a0,
                                             fit.fold[[k]]$fold.ghat$beta), 
                                lambda = fit.fold[[k]]$fold.ghat$lambda,
                                dev.cvm = tune.fold[[k]]$cvm.dev.gamma,
                                mag.cvm = tune.fold[[k]]$cvm.mag.gamma)
  }
  # ATE 
  #===============
  
  # No cross-fitting --------------------------------
  
  ate = list()
  pD = 1/(1+exp(-drop(Z%*%hgr[-1])-
                  hgr[1]))
  hazbZ = drop(Z%*%hthetabeta[-1])
  # hat : phi with lasso beta and lasso gamma
  tmpfit = phi(hazbZ,hthetabeta[1],pD,
               D,surv,maxit,tol)
  ate$hat = tmpfit$ate
  ate$hat.sd = tmpfit$ate.sd
  fit$hLam = tmpfit$Lam
  
  # checkhat : psi with lasso beta and lasso gamma
  tmpfit = psi(hazbZ,pD,
               D,surv)
  ate$checkhat = tmpfit$ate
  ate$checkhat.sd = tmpfit$ate.sd
  #ate$cLamhat = tmpfit$Lam
  
  # Cross-fitting ----------------------------------
  atecf = list()
  pD = ps.cf(foldid, hgr.cf,Z)
  hazbZ = hazbZ.cf(foldid, 
                   hthetabeta.cf[-1,],Z)
  tmpfit = phi.cf(foldid,
                  hthetabeta.cf[1,],hazbZ, pD,
                  D,surv,maxit,tol)
  atecf$hat = tmpfit$ate
  atecf$hat.sd = tmpfit$ate.sd
  
  # check : psi with regularized lasso beta and regularized lasso gamma
  tmpfit = psi.cf(foldid,hazbZ, pD,
                  D,surv)
  atecf$check = tmpfit$ate
  atecf$check.sd = tmpfit$ate.sd
  
  return(list(fit=fit,ate=ate,atecf=atecf))
}
