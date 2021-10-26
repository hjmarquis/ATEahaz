phi.cf = function(foldid,init,hazbZ,pD,D,surv, maxit,tol)
{
  fit = list()
  nfold = max(foldid)
  fit$hlam = vector("list",nfold)
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
  
  # initialize
  orD = hazb = hazl = Zb = rep(0,n)
  
  # hat theta: phi with hat beta star and hat gamma
  idin = idout = vector("list",nfold)
  for(k in 1:nfold)
  {
    foldk = foldid==k
    idin[[k]] = is.element(time.order,foldk)
    idout[[k]] = !idin[[k]]
  }
  
  init.theta = mean(init)
  
  w1 = (D*(1-pD))[time.order]
  w0 = ((1-D)*pD)[time.order]
  
  hazl = drop(hazbZ+init[foldid]*D)[time.order]
  hazb = drop(hazbZ)[time.order]

  SY = SZ=vector("list",nfold)
  LamX = rep(0,n)
  for(i in 1:nfold)
  {
    SY[[i]] = cumsum(Y*idout[[i]])[tie.last]
    SZ[[i]] = cumsum(hazl*idout[[i]])[tie.last]
    dLam = (event*idout[[i]]-SZ[[i]]*dX)/SY[[i]]
    dLam[is.nan(dLam)] = 0
    LamXi = rev(cumsum(rev(dLam)))[tie.first]
    LamX[idin[[i]]] = LamXi[idin[[i]]]
    fit$hlam[[i]] = stepfun(sort(X),c(0,sort(LamXi)))
  }
  
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
      eLamX = rep(0,n)
      for(i in 1:nfold)
      {
        edLam = (event*idout[[i]]*eATEX-deX*SZ[[i]])/SY[[i]]
        edLam[is.nan(edLam)] = 0
        eLamX[idin[[i]]] = rev(cumsum(rev(edLam)))[tie.first[idin[[i]]]]
      }
      eATEX.ATE = (eATEX-1)/ATE
      deATEX.ATE = (eATEX*X-eATEX.ATE)/ATE
    }
    dX2 = -diff(c(deATEX.ATE,0))
    teLamX = rep(0,n)
    for(i in 1:nfold)
    {
      dteLamX = (event*eATEX*X*idout[[i]] - dX2*SZ[[i]])/SY[[i]]
      dteLamX[is.nan(dteLamX)] = 0
      teLamX[idin[[i]]] = rev(cumsum(rev(dteLamX)))[tie.first[idin[[i]]]]
    }
    
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