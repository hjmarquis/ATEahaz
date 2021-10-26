source("core/utility.R")

mag.ahaz.v1 = function(hazbZ, X)
{
  n = length(X)
  time.order = order(-X)
  X = X[time.order]
  hazbZ = hazbZ[time.order]
  Y = rep(1,n)
  
  dX = -diff(c(X,0))
  tie.last = tie.first = 1:n
  for(i in 2:n)
  {
    if(dX[i-1]==0)
    {
      tie.first[i] = tie.first[i-1]
      tie.last[tie.first[i]:i] = i
    }
  }
  
  S1 = cumsum(Y*hazbZ)[tie.last]
  S0 = cumsum(Y)[tie.last]
  
  sqrt(mean(hazbZ^2*X)-
         sum(S1^2/S0/n*dX))
  
  # if(length(beta)>1)
  #   return(sqrt(apply((Z%*%beta)^2*X,2,mean)))
  # return(sqrt(mean((Z*beta)^2*X)))
}

mag.logistic.v1 = function(D,X,Ytau,pD)
{
  pnull = mean(D)
 m0 = mean((1-D)*X)*pnull
 m1 = mean(D*Ytau)*(1-pnull)
  if(is.null(dim(pD)))
  {
    mg0 = mean((1-D)*X*pD)
    mg1 = mean(D*Ytau*(1-pD))
  }else
  {
    mg0 = apply((1-D)*X*pD,2,mean)
    mg1 = apply(D*Ytau*(1-pD),2,mean)
  }
  return(1/mg0+1/mg1)
}

"ahaz.mse"<-function(x, beta)
{
  ## Purpose: MSE from ahaz object and given beta
  ## ----------------------------------------------------------------------
  ## Arguments:
  ##   x     : 'ahaz' object
  ##   beta  : beta estimate
  ## ----------------------------------------------------------------------
  ## Author: Anders Gorst-Rasmussen
  if(length(beta)>1)
    return((t(beta) %*% x$D %*% beta - 2 * t(beta) %*% x$d))
  return(x$D*beta^2-2*beta*x$d)
}


logistic.dev = function(D,p)
{
#  return(-apply(D*log(p)+(1-D)*log(1-p),2,mean))
  -(apply(log(p[D==1,]),2,sum)+apply(log(1-p[D==0,]),2,sum))/length(D)
}
  