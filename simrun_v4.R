version = "v4"
library(doParallel)
np = detectCores(logical = F)/2
cl = makeCluster(getOption("cl.cores", np))
registerDoParallel(cl)

n = c(300,1500)
N =500
bslist = c(1,1,3,2,4,1,5,1,1)
gslist = c(1,3,1,2,1,4,1,5,8)
slabel = c("v","s","m","d","l","q","u", "f")
nsetup = length(bslist)

source("core/simATEv4cf.R")
set.seed(531)
for(j in 1:nsetup)
  for(ni in 1:length(n))
  {
    bs = bslist[j]
    gs = gslist[j]
      
      print(paste("n=",n[ni],", b=",slabel[bs],", g=",slabel[gs],sep=''))
      resultfile = paste("result/",version,"/n",n[ni],"b",slabel[bs],"g",slabel[gs],".rda",sep='')

      simATE<-foreach(k = 1:N, .packages=c("ahaz","glmnet"),.errorhandling="stop") %dopar%
      {
        load(paste("simdata/",version,"/n",n[ni],"b",slabel[bs],"g",slabel[gs],"_",k,".rda",sep=''))
        yo = sim.ATE.v4cf(simdata$surv, simdata$D, simdata$Z,
                   true.ate = simdata$true.ate,true.beta = simdata$beta[-1],
                   true.gamma = simdata$gamma[-2], oracle = (bs <4 & gs<4),
                   true.lam0 = simdata$lam0, nfold = 10)
        yo
      }
      
      if(bs>4)
      {
        least.false.beta <- foreach(k = 1:N, .combine = rbind) %dopar%
        {
          load(paste("simdata/",version,"/n",n[ni],"b",slabel[bs],"g",slabel[gs],"_",k,".rda",sep=''))
          simdata$least.false.beta
        }
        simATE$least.false.beta = drop(apply(least.false.beta,2,mean))
      }
      if(gs>4)
      {
        least.false.gamma <- foreach(k = 1:N, .combine = rbind) %dopar%
        {
          load(paste("simdata/",version,"/n",n[ni],"b",slabel[bs],"g",slabel[gs],"_",k,".rda",sep=''))
          simdata$least.false.gamma
        }
        simATE$least.false.gamma = drop(apply(least.false.gamma,2,mean))
      }
      save(simATE,file=resultfile)
  }
stopCluster(cl)

