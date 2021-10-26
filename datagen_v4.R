#=======================================#
#                                       #
#           Data Simulation             #
#                                       #
#=======================================#

version = "v4"

# Parameters
#===========

source("simulation/simgen.R")
# ATE
true.ate = -0.25

# Dimension
rho = 0
lam.min = 0
Z.sd = 1
strong.b = 1
weak.b = 0.1
strong.g = 1
weak.g = 0.05

# Setup
bslist = c(1,1,3,2,4,1,5,1,1)
gslist = c(1,3,1,2,1,4,1,5,8)
slabel = c("v","s","m","d","l","q","u","f")
nsetup = length(bslist)
Dlinklist = list(logitlink,logitlink,logitlink,logitlink,
                 pnorm, logitlink,logitlink, function(x) as.numeric(x>0))
ZDlinklist = list(identity,identity,identity,identity,
                  identity,function(x) x^2,identity,identity)
invLamlist = list(inv.Lam.const,inv.Lam.const,inv.Lam.const,inv.Lam.const,
                  invT.Cox,inv.Lam.const,inv.Lam.const)

# Tune
p = 100
slistb.strong = c(1,1,2,4)
slistb.weak   = c(1,5,13,26)
slistg.strong = c(1,1,2,4)
slistg.weak   = c(0,2,8,16)

bscoef =  list(
  very.sparse = c(0,rep(c(strong.b,weak.b),c(slistb.strong[1],slistb.weak[1])),
                  rep(0,p-slistb.strong[1]-slistb.weak[1])),
  sparse = c(0,rep(c(strong.b,weak.b),c(slistb.strong[2],slistb.weak[2])),
             rep(0,p-slistb.strong[2]-slistb.weak[2])),
  moderately.sparse = c(0,rep(c(strong.b,weak.b),c(slistb.strong[3],slistb.weak[3])),
                        rep(0,p-slistb.strong[3]-slistb.weak[3])),
  dense = c(0,rep(c(strong.b,weak.b),c(slistb.strong[4],slistb.weak[4])),
            rep(0,p-slistb.strong[4]-slistb.weak[4])),
  misslink = c(0,rep(c(strong.b,weak.b),c(slistb.strong[1],slistb.weak[1])),
               rep(0,p-slistb.strong[1]-slistb.weak[1])),
  missquad = c(0,rep(c(strong.b,weak.b),c(slistb.strong[1],slistb.weak[1])),
               rep(0,p-slistb.strong[1]-slistb.weak[1])),
  missomit = c(1,rep(c(strong.b,weak.b),c(slistb.strong[1],slistb.weak[1])),
               rep(0,p-slistb.strong[1]-slistb.weak[1]))
)

gscoef = list(
  very.sparse = c(0,rep(c(strong.g,weak.g),c(slistb.strong[1],slistb.weak[1])),
                  rep(0,p-slistb.strong[1]-slistb.weak[1])),
  sparse = c(0,rep(c(strong.g,weak.g),c(slistb.strong[2],slistb.weak[2])),
             rep(0,p-slistb.strong[2]-slistb.weak[2])),
  moderately.sparse = c(0,rep(c(strong.g,weak.g),c(slistb.strong[3],slistb.weak[3])),
                        rep(0,p-slistb.strong[3]-slistb.weak[3])),
  dense = c(0,rep(c(strong.g,weak.g),c(slistb.strong[4],slistb.weak[4])),
            rep(0,p-slistb.strong[4]-slistb.weak[4])),
  misslink = c(0,rep(c(strong.g,weak.g),c(slistb.strong[1],slistb.weak[1])),
               rep(0,p-slistb.strong[1]-slistb.weak[1])),
  missquad = c(0,rep(c(strong.g,weak.g),c(slistb.strong[1],slistb.weak[1])),
               rep(0,p-slistb.strong[1]-slistb.weak[1])),
  missomit = c(1,rep(c(strong.g,weak.g),c(slistb.strong[1],slistb.weak[1])),
               rep(0,p-slistb.strong[1]-slistb.weak[1])),
  fixed    = -c(0,rep(c(strong.g,weak.g),c(slistb.strong[1],slistb.weak[1])),
               rep(0,p-slistb.strong[1]-slistb.weak[1]))
)
Ccoef = rep(0,p+2)
source("simulation/tunesim_main.R")

# Repeats
N = 500
n = c(300,1500)

# Simdata
set.seed(531)

for(j in 1:nsetup)
{
  
  bs = bslist[j]
  gs = gslist[j]
  
  for(ni in 1:length(n))
  {
    p = n[ni]
    Ccoef = rep(0,p+2)
    
    print(paste("n=",n[ni],", b=",slabel[bs],", g=",slabel[gs],sep=''))
    for(k in 1:N)
    {
      filename = paste("simdata/",version,"/n",n[ni],"b",slabel[bs],"g",slabel[gs],"_",k,".rda",sep='')
      beta = c(bscoef[[bs]],rep(0,p+1-length(bscoef[[bs]])))
      gamma = c(gintercept[j],gscoef[[gs]],rep(0,p+1-length(gscoef[[gs]])))
      Cunif = function(x){runif(length(x),0,Cmax[j])}
      Z = conf.gen.normal(n[ni],p+1,Z.sd,rho,beta,lam.min+true.ate, ZDlinklist[[bs]])
      simdata = data.gen(true.ate,beta,
                         gamma,Ccoef,Z,
                         inv.Lam=function(d,z) 
                           invLamlist[[bs]](d,z,lam.min),
                         C.gen =Cunif,
                         tau=tau[j],
                         Dlink = Dlinklist[[gs]],
                         ZDlink = ZDlinklist[[gs]],
                         ZTlink = ZDlinklist[[bs]])
      simdata$true.ate = true.ate
      simdata$beta = beta
      simdata$gamma = gamma
      simdata$lam0 = lam.min
      simdata$eps = chat.eps[j]
      if(bs>=5)
      {
        simdata$least.false.beta = rep(0,n[ni]+1)
        bpos = which(simdata$beta[-1] !=0)
        simdata$least.false.beta[c(1,1+bpos)] = coef(ahaz(simdata$surv,cbind(simdata$D, simdata$Z[,bpos])))
      }
      if(gs>=5)
      {
        simdata$least.false.gamma = rep(0,n[ni]+1)
        bpos = which(simdata$gamma[-1:-2] !=0)
        simdata$least.false.gamma[c(1,1+bpos)] = coef(glm(simdata$D~simdata$Z[,bpos],family = binomial))
      }
      
      save(simdata,file = filename)
    }
    
  }
}
