version = "v4"
library(ahaz)

source("simulation/simgen.R")
source("core/psi-cf.R")
source("core/measure_v1.R")
# ATE
true.ATE = -0.25
lam.min = 0
# Dimension
n = c(300,1500)
N =500

bslist = c(1,1,3,2,4,1,5,1,1)
gslist = c(1,3,1,2,1,4,1,5,8)
slabel = c("v","s","m","d","l","q","u","f")
nfold = 10

nsetup = length(bslist)
Dlinklist = list(logitlink,logitlink,logitlink,logitlink,
                 pnorm, logitlink,logitlink, function(x) as.numeric(x>0))
ZDlinklist = list(identity,identity,identity,identity,
                  identity,function(x) x^2,identity,identity)

hat=checkhat=hat.se=checkhat.se=hat.cf=
  check.cf=hatcf.se=checkcf.se=hbeta.ad=
  hgr.ad =hbeta.mse=hbeta.mag= hgr.mse =hgr.mag= 
  hLam.error =matrix(NA, length(n)*nsetup,N)

bslist.Kb = c(4,1,5,1,1)
gslist.Kb = c(1,4,1,5,8)
nKbsetup = length(bslist.Kb)
nKb = 10
Kb.summary = data.frame(bs = rep(bslist.Kb,each = nKb^2*length(n)),
                              gs = rep(bslist.Kb,each = nKb^2*length(n)),
                              n  = rep(n,each = nKb^2,nKbsetup),
                              Kg = 0, bn = 0,
                              bias.lasso = 0, bias = 0,  MSE.lasso = 0, MSE=0)
Kb.check = matrix(0,nKb^2,N)

thetalasso = thetaps =  matrix(NA, length(n)*nsetup,N)
  for(ni in 1:length(n))  
    for(isetup in 1:nsetup)
    {
      bs = bslist[isetup]
      gs = gslist[isetup]
      
      iKb = which(bs==bslist.Kb & gs == gslist.Kb)
      resultfile = paste("result/",version,"/n",n[ni],"b",slabel[bs],"g",slabel[gs],".rda",sep='')
      shrinkfile = paste("result/",version,"/n",n[ni],"b",slabel[bs],"g",slabel[gs],"shrinkage.rda",sep='')
      fitobj = load(resultfile)
      
      
        
        print(paste("n=", n[ni], ", beta ", slabel[bs], ", gamma ",
                    slabel[gs], sep=''))
        for(k in 1:N)
        {

          dataobj = load(paste("simdata/",version,"/n",n[ni],"b",slabel[bs],"g",slabel[gs],"_",k,".rda",sep=''))
          
          # Treatment Effect
          # =============================================================================
          # Treatment Effect without cross-fitting
          hat[(ni-1)*nsetup + isetup,k]=simATE[[k]]$ate$hat
          checkhat[(ni-1)*nsetup + isetup,k]=simATE[[k]]$ate$checkhat
          hat.se[(ni-1)*nsetup + isetup,k]=simATE[[k]]$ate$hat.sd
          checkhat.se[(ni-1)*nsetup + isetup,k]=simATE[[k]]$ate$checkhat.sd
          
          hat.cf[(ni-1)*nsetup + isetup,k]=simATE[[k]]$atecf$hat
          check.cf[(ni-1)*nsetup + isetup,k]=simATE[[k]]$atecf$check
          hatcf.se[(ni-1)*nsetup + isetup,k]=simATE[[k]]$atecf$hat.sd
          checkcf.se[(ni-1)*nsetup + isetup,k]=simATE[[k]]$atecf$check.sd 
          
          
          # Nuisance parameters
          #==============================================================================
          
          # Propensity model 
          if(gs <= 4)
          {
            hgr.ad[(ni-1)*nsetup + isetup,k] = 
              sum(abs(simATE[[k]]$fit$hgr-simATE[[k]]$fit$true.gamma))
          }else{
            hgr.ad[(ni-1)*nsetup + isetup,k] = 
              sum(abs(simATE[[k]]$fit$hgr-simdata$least.false.gamma))
          }
          
          
          
          if(gs <=7)
          {
            true.pr = Dlinklist[[gs]](simdata$gamma[1]+ simdata$U *simdata$gamma[2]+
                                          drop(ZDlinklist[[gs]](simdata$Z)%*%simdata$gamma[-1:-2]))
          }else if(gs == 8)
          {
            true.pr = simdata$D
          }
          hpr =  rep(0, n[ni])
          for(foldk in 1:nfold)
          {
            foldk.id = simATE[[k]]$fit$foldid==foldk
            
            foldgr = drop(simATE[[k]]$fit$hgr.cf[,foldk])
            hpr[foldk.id] = logitlink(foldgr[1]+ 
                                      drop(simdata$Z[foldk.id,]%*%foldgr[-1]))
          }
          hgr.mse[(ni-1)*nsetup + isetup,k] = 
            sqrt(mean((hpr - true.pr)^2))
          
          # Survival model
          
          if(bs <= 4)
          {
            hbeta.ad[(ni-1)*nsetup + isetup,k] = 
              sum(abs(simATE[[k]]$fit$hthetabeta[-1]-simATE[[k]]$fit$true.beta))
            x = unique(knots(simATE[[k]]$fit$hLam))
            x = x[x<=quantile(x,probs = .75)]
            Lam0x = rep(x*lam.min,2)
            dx.min = min(diff(x))/2
            xeval = c(x-dx.min,x+dx.min)
            hLam.error[(ni-1)*nsetup + isetup,k] = max(abs(simATE[[k]]$fit$hLam(xeval)-Lam0x))
          }else{
            hbeta.ad[(ni-1)*nsetup + isetup,k] = 
              sum(abs(simATE[[k]]$fit$hthetabeta[-1]-simdata$least.false.beta[-1]))
          }
          
          hbz = rep(0, n[ni])
            for(foldk in 1:nfold)
            {
              foldk.id = simATE[[k]]$fit$foldid==foldk
              
              foldbeta = simATE[[k]]$fit$hthetabeta.cf[-1,foldk]
              hbz[foldk.id] = drop(simdata$Z[foldk.id,]%*%foldbeta)
            }
          if(bs != 5)
          {
            true.bz = ZDlinklist[[bs]](simdata$Z)%*%simdata$beta[-1]
           }else{
            true.bz = exp(ZDlinklist[[bs]](simdata$Z)%*%simdata$beta[-1])
          } 
            
            hbeta.mse[(ni-1)*nsetup + isetup,k] = 
              sqrt(mean((hbz-true.bz)^2 * simdata$surv[,1]))

          
#-------prediction errors-------- 
          thetalasso[(ni-1)*nsetup + isetup,k] =
           simATE[[k]]$fit$theta.nopen
          
          hor = exp(drop(simdata$Z%*%simATE[[k]]$fit$hgr[-1])+simATE[[k]]$fit$hgr[1])
          hps = 1/(1+1/hor)
          hps[hor==0] = 0
          thetaps[(ni-1)*nsetup + isetup,k] = coef(ahaz(simdata$surv,cbind(simdata$D,hps)))[1]
          
          hazbZ = simdata$Z%*%simATE[[k]]$fit$hthetabeta.cf[-1,]
          pD = 1/(1+exp(-rep(simATE[[k]]$fit$hgr.cf[1,],each = n[ni])-
            simdata$Z%*%simATE[[k]]$fit$hgr.cf[-1,]))
          Kg = bn = rep(0,nfold)
          tau = max(simdata$surv[,1])
          for(ifold in 1:nfold)
          {
            foldk.id = simATE[[k]]$fit$foldid==foldk
            Kg[ifold] = mag.logistic.v1(simdata$D[foldk.id],
                                        simdata$surv[foldk.id,1],
                                        simdata$surv[foldk.id,1]==tau,
                                        pD[foldk.id,ifold])
            bn[ifold] = mag.ahaz.v1(simdata$surv[foldk.id,1],
                                    hazbZ[foldk.id,ifold])
          }
          hbeta.mag[(ni-1)*nsetup + isetup,k] = max(bn)
          hgr.mag[(ni-1)*nsetup + isetup,k] = max(Kg)
          
        }
        
      }

rm(list = c(dataobj,fitobj))
save(list = objects(), file = paste("result/",version,"/summarycf.RData",sep = ''))
