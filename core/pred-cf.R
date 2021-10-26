ps.cf = function(foldid,hgr,Z)
{
  nfold = max(foldid)
  n = nrow(Z)
  pD = rep(0,n)
  for(k in 1:nfold)
  {
    foldk = foldid==k
    pD[foldk] = 1/(1+exp(-drop(Z[foldk,]%*%drop(hgr[-1,k]))-hgr[1,k]))
  }
  return(pD)
}

hazbZ.cf = function(foldid,hbeta,Z)
{
  nfold = max(foldid)
  n = nrow(Z)
  hazbZ = rep(0,n)
  for(k in 1:nfold)
  {
    foldk = foldid==k
    hazbZ[foldk] = Z[foldk,] %*% drop(hbeta[,k])
  }
  return(hazbZ)
}