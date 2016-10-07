#v3 environment model submit
library(rjags)
library(foreach)              
library(doParallel)
library(random)
cl <- makePSOCKcluster(3)              
clusterSetRNGStream(cl)
registerDoParallel(cl)

setwd("~/stat/CMC/v3/");
#".RNG.state" <- c(19900, 14957, 25769)

# constant to be fitted according to the goal
nburn <- 3000;
niter <- 10000;
thin <- 1;
calc <- paste(nburn,niter,thin,sep="|");
version <- "v3";
tomonitor <- c("ec.r","ec.logr","ec.loglambda","mu.ec.logr","sd.ec.logr","mu.lam","sig.lam");

filedat <- paste(version,"data","r",sep=".");
filemod <- paste(version,"env.model","jags",sep=".");

#load data
source(filedat);
setwd("~/stat/CMC/v3/");

init <- function( ) {
  return( list( "ec.data"=ec.init,
                .RNG.name = "lecuyer::RngStream", 
                .RNG.seed = randomNumbers( n = 1, min = 1, max = 1e+06, col = 1 ) ) )
}


### POSTERIOR analysis
print(paste("Start time: ",Sys.time(),sep=""))
mcmc.res <- foreach(i=1:3, .packages=c('rjags','random','coda'), .multicombine=TRUE) %dopar% {
  load.module( "lecuyer" )
  model.jags <- jags.model(file="~/stat/CMC/v3/v3.env.model.jags",data=cntdata,n.chains=1,n.adapt=nburn,
                           inits=init,quiet=FALSE);
  update(model.jags,n.burn=nburn, progress.bar="text");
  mcmc.post <- coda.samples(model.jags,tomonitor,n.iter=niter,thin=thin);
  return(mcmc.post)
}
print(paste("End time: ",Sys.time(),sep=""))
stopCluster(cl)

save(mcmc.res,file="~/stat/CMC/v3/output/v3.mcmc.post.res.rda")

mcmc.post<-as.mcmc.list(c(mcmc.res[[1]],mcmc.res[[2]],mcmc.res[[3]]))

save(mcmc.post,file="~/stat/CMC/v3/output/v3.mcmc.post.rda")

pdf(file="~/stat/CMC/v3/output/mcmc_results.pdf")
plot(mcmc.post)
dev.off()