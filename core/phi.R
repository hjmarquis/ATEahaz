phi = function(hazbZ,init.theta,pD,D,surv, maxit, tol)
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
  
  hazl = drop(hazbZ+init.theta*D)[time.order]
  hazb = drop(hazbZ)[time.order]
  
  SY = cumsum(Y)[tie.last]
  SZ = cumsum(hazl)[tie.last]
  dLam = (event-SZ*dX)/SY
  LamX = rev(cumsum(rev(dLam)))[tie.first]
  fit$Lam = stepfun(sort(X),c(0,sort(LamX)))
  
  ATE = init.theta
  for(step in 1:maxit)
  {
    eATEX = drop(exp(ATE*X))
    if(ATE==0)
    {
      eLamX = LamX
      eATEX.ATE = X
      deATEX.ATE = X^2/2
    }else
    {
      deX = -diff(c(eATEX,1))/ATE
      edLam = (event*eATEX-deX*SZ)/SY
      eLamX = rev(cumsum(rev(edLam)))[tie.first]
      eATEX.ATE = (eATEX-1)/ATE
      deATEX.ATE = (eATEX*X-eATEX.ATE)/ATE
    }
    dX2 = -diff(c(deATEX.ATE,0))
    dteLamX = (event*eATEX*X - dX2*SZ)/SY
    teLamX = rev(cumsum(rev(dteLamX)))[tie.first]
    
    score = sum(-w0*(event-hazb*X-LamX)
                +w1*(event*eATEX-eATEX.ATE*(hazb+ATE) - eLamX))
    dscore = sum(w1*(-eATEX.ATE+event*X*eATEX -
                       deATEX.ATE*(hazb+ATE)-teLamX))
    ATE.new = ATE -score/dscore
    if(abs(ATE-ATE.new)<tol)
      break
    ATE = ATE.new
    if(is.infinite(ATE)|exp(ATE*X[1])==Inf)
    {
      ATE = NaN
      break
    }
  }
  
  var.ATE = sum(event*(w1-w0)^2*exp(2*ATE*DX))/sum(w0*X)^2
  
  fit$ate = ATE
  fit$ate.sd = sqrt(var.ATE)
  
  return(fit)
}