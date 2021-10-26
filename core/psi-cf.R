psi.cf = function(foldid,hazbZ,pD,D,surv)
{
  fit = list()
  nfold = max(foldid)
  n = length(D)
  
  ate.num=ate.denom = rep(0,nfold)
  sd.ew = sd.DX = rep(0,n)  
  clamk = Xk = vector("list",nfold)
  for(k in 1:nfold)
  {
    foldk = (foldid==k)
    nk = sum(foldk)
    survk = surv[foldk,]
    Dk = D[foldk]
    time.order = order(-survk[,1],-survk[,2])
    X = survk[time.order,1]
    event = survk[time.order,2]
    Y = rep(1,nk)
    trt = (Dk==1)[time.order]
    DX = Dk[time.order]*X
    
    dX = -diff(c(X,0))
    tie.last = tie.first = 1:nk
    for(i in 2:nk)
    {
      if(dX[i-1]==0)
      {
        tie.first[i] = tie.first[i-1]
        tie.last[tie.first[i]:i] = i
      }
    }
    
    w1 = (Dk*(1-pD[foldk]))[time.order]
    w0 = ((1-Dk)*pD[foldk])[time.order]
    haZ = drop(hazbZ[foldk])[time.order]
    
    S1 = cumsum(w1)[tie.last]
    S0 = cumsum(w0)[tie.last]
    S1Z = cumsum(w1*haZ)[tie.last]
    
    posS1 = S1>0
    
    dclam = (w1*event-S1Z*dX)/S1
    dclam[!posS1]=0
    clamk[[k]] = rev(cumsum(rev(dclam)))[tie.first]
    
    Xk[[k]] = X
    
    ate.num[k] = (-sum(w0*(event-X*haZ))-
                    sum((S0*S1Z/S1*dX)[posS1])+
                    sum((event*w1*S0/S1)[posS1]))
    
    ate.denom[k] = sum(w0*X)
    sd.ew[foldk] = event*(w1-w0)^2
    sd.DX[foldk] = DX
  }
  ATE = sum(ate.num)/sum(ate.denom)
  
  fit$ate = ATE
  fit$ate.sd = sqrt(sum(sd.ew*exp(2*ATE*sd.DX))/sum(ate.denom)^2)

  fit$clam = vector("list",nfold)
  for(k in 1:nfold)
  {
    fit$clam[[k]] = stepfun(sort(Xk[[k]]), c(0,sort(clamk[[k]]-ATE*Xk[[k]])))
  }
  return(fit)
}