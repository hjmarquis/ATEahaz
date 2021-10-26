#=======================================#
#                                       #
#           Tune Simulation             #
#                                       #
#=======================================#

# Parameters
#===========

 gintercept = Cmax = chat.eps = an = bn = tau = rep(0,nsetup)
 # least.false = vector("list",nsetup)
 
 # maxiter = 100
 # tol = 1e-6
# Simdata
source("simulation/simgen.R")
set.seed(531)
n=10000
for(j in 1:nsetup)
{
  bs = bslist[j]
  gs = gslist[j]
  {
    print(paste("b=",slabel[bs],", g=",slabel[gs],sep=''))
    beta = bscoef[[bs]]
    Z = conf.gen.normal(n,p+1,Z.sd,rho,beta,lam.min+true.ate,ZDlinklist[[bs]])
    gZ = Z%*%gscoef[[gs]]
    if(gs==5)
    {
      Drate = function(x)
      {
        mean(pnorm(x+gZ))-0.5
      }
      gintercept[j] = uniroot(Drate,c(-100,100))$root
      gamma = c(gintercept[j],gscoef[[gs]])
      p1 = pnorm(drop(gZ + gintercept[j]))
      p0=1-p1
    }else if(gs==6)
    {
      gZ = (Z^2)%*%gscoef[[gs]]
      Drate = function(x)
      {
        or = exp(x+gZ)
        mean(or/(1+or))-0.5
      }
      gintercept[j] = uniroot(Drate,c(-100,100))$root
      gamma = c(gintercept[j],gscoef[[gs]])
      or = exp(drop(gZ + gintercept[j]))
      p0 = 1/(1+or)
      p1 = 1-p0
    }else if(gs==8)
    {
      gintercept[j] = - median(gZ)
      gamma = c(gintercept[j],gscoef[[gs]])
      p0 = drop(gZ)< -gintercept[j]
      p1 = 1-p0
    }else{
      Drate = function(x)
      {
        or = exp(x+gZ)
        mean(or/(1+or))-0.5
      }
      gintercept[j] = uniroot(Drate,c(-100,100))$root
      gamma = c(gintercept[j],gscoef[[gs]])
      or = exp(drop(gZ + gintercept[j]))
      p0 = 1/(1+or)
      p1 = 1-p0
    }
    
    if(bs==5)
    {
      haz0 = exp(drop(Z%*%beta))+max(0,-true.ate)
      haz1 = haz0 + true.ate
    }else{
      haz0 = drop(ZDlinklist[[bs]](Z)%*%beta)+lam.min
      haz1 = haz0 + true.ate
    }
     

     Crate = function(x)
     {
       taurate = function(y)
         mean(p1*exp(-haz1*y)*(x-y)/x)-0.1
       tau = uniroot(taurate, c(0.1,x-0.1))$root 
       mean(p1*((1/haz1+exp(-haz1*tau)*(x-tau-1/haz1))/x)
            + p0*((1/haz0+exp(-haz0*tau)*(x-tau-1/haz0))/x))-0.30
     }
     Cmax[j] = uniroot(Crate, c(.1,100))$root
     taurate = function(y)
       mean(p1*exp(-haz1*y)*(Cmax[j]-y)/Cmax[j])-0.1
     tau[j] = uniroot(taurate, c(0.1,100))$root
     Cunif = function(x){runif(length(x),0,Cmax[j])}
    
    simdata = data.gen(true.ate,beta,
                       gamma,Ccoef,Z,
                       inv.Lam=function(d,z) 
                         invLamlist[[bs]](d,z,lam.min),
                       C.gen =Cunif,
                       tau=tau[j],
                       Dlink = Dlinklist[[gs]],
                       ZDlink = ZDlinklist[[gs]],
                       ZTlink = ZDlinklist[[bs]])
    chat.eps[j] = mean((1-simdata$D)*p1*simdata$surv[,1])
    an[j] = sum(abs(gamma))
    bn[j] = sum(abs(beta))
    print(table(simdata$D)/n)
    print(table(simdata$surv[,2])/n)
    print(table(simdata$surv[,1]==tau[j],simdata$D)/n)
  }
}

save(gintercept,Cmax,chat.eps,an,bn,file=paste("simulation/simconst_main_",version,".rda",sep=''))

