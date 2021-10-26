psi = function(hazbZ,pD,D,surv)
{
  fit = list()
  
  n = length(D)
  time.order = order(-surv[,1],-surv[,2])
  X = surv[time.order,1]
  event = surv[time.order,2]
  Y = rep(1,n)
  trt = (D==1)[time.order]
  DX = D[time.order]*X
  
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
  
  w1 = (D*(1-pD))[time.order]
  w0 = ((1-D)*pD)[time.order]
  haZ = drop(hazbZ)[time.order]
  
  S1 = cumsum(w1)[tie.last]
  S0 = cumsum(w0)[tie.last]
  S1Z = cumsum(w1*haZ)[tie.last]
  
  posS1 = S1>0
  ate.num = (-sum(w0*(event-X*haZ))-
               sum((S0*S1Z/S1*dX)[posS1])+
               sum((event*w1*S0/S1)[posS1]))
  
  ate.denom = sum(w0*X)
  ATE = ate.num/ate.denom
  
  fit$ate = ATE
  fit$ate.sd = sqrt(sum(event*(w1-w0)^2*exp(2*ATE*DX))/ate.denom^2)
  
  hazl = drop(hazbZ+ATE*D)[time.order]
  S1Z = cumsum(w1*hazl)[tie.last]
  dLam = (w1*event-S1Z*dX)/S1
  dLam[!posS1]=0
  LamX = rev(cumsum(rev(dLam)))[tie.first]
  fit$Lam = stepfun(sort(X),c(0,sort(LamX)))
  
  return(fit)
}