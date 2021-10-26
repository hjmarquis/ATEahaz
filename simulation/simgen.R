library(survival)
library(Matrix)

conf.gen.mix = function(n,scont,sbin,pcont,pbin)
{
  as.matrix(cbind(
    matrix(runif(n*scont),n),
    matrix(rbinom(n*sbin,1,0.5),n),
    matrix(runif(n*(pcont-scont)),n),
    matrix(rbinom(n*(pbin-sbin),1,0.5),n)
  ))
}

conf.gen.normal = function(n,p,Z.sd,ar.rho,beta,lam.min,Zlink=identity,buffer = 2)
{
  cormat = triu(ar.rho^outer(-1:-p,-1:-p,'-'),k=1)
  diag(cormat) = sqrt(1-apply(cormat^2,2,sum))
  
  nb = round(n*buffer)
  Z =as.matrix( matrix(rnorm(nb*p),nb,p) %*% cormat) * Z.sd
  accept = which(Zlink(Z)%*%beta > -lam.min)
  n.accept = min(n,length(accept))
  n = n - n.accept
  Z = Z[accept[1:n.accept],]
  
  while(n>0)
  {
    new.Z = as.matrix( matrix(rnorm(nb*p),nb,p) %*% cormat)  * Z.sd
    accept = which(new.Z%*%beta > -lam.min)
    n.accept = min(n,length(accept))
    n = n - n.accept
    Z = rbind(Z,new.Z[accept[1:n.accept],])
  }
  return(Z)
}

inv.Lam.const = function(Dth,bZ,lam)
{
  dhaz = Dth+bZ+lam
  if(any(dhaz <0))
    stop("Negative hazard!")
  -log(runif(length(Dth)))/(dhaz)
}

censor.exp = function(rate)
{
  pmin(rexp(length(rate),rate))
}

logitlink = function(x)
{
  ex = exp(x)
  ex/(1+ex)
}

invT.Cox = function(Dth,bZ,lam)
{
  dhaz = Dth-min(Dth)+exp(bZ)*1
  if(any(dhaz <0))
    stop("Negative hazard!")
  -log(runif(length(Dth)))/(dhaz)
}

data.gen = function(theta,beta,gamma,C.coef,Z,
                    inv.Lam=function(d,z) 
                      inv.Lam.const(d,z,1+max(0,-theta-sum(beta),-sum(beta))),
                    Dlink = logitlink,
                    C.gen = function(x) censor.exp(x),
                    tau = 1,
                    ZDlink = identity,
                    ZTlink = identity)
{
  n=nrow(Z)
  pD = Dlink(drop(gamma[1]+ZDlink(Z)%*%gamma[-1]))
  D = rbinom(n,1,pD)
  dhaz = theta*D + drop(Z%*%beta)
  Tevent = inv.Lam(theta*D,  drop(ZTlink(Z)%*%beta))
  C.par = exp(drop(C.coef[1]+Z%*%C.coef[-1]))
  C = C.gen(C.par)
  X = pmin(Tevent,C)
  X = pmin(X,tau)
  surv = Surv(X,as.numeric(Tevent==X))
  
  return(list(surv=surv,D=D,Z=Z[,-1],U=Z[,1],pD = pD))
}