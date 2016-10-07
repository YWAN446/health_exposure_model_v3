#model diagnostics
load("/Users/yukew/stat/CMC/v3/output/archived/v3.mcmc.post.rda")
load("/Users/yukew/stat/CMC/v3/output/archived/v3.mcmc.post.res.rda")
mcmc.pos<-mcmc.post[,701:741]

#traceplot and density plot
plot(mcmc.pos)
#auto correlation plot
autocorr.plot(mcmc.pos)
#running mean plot
library(mcmcplots)
rmeanplot(mcmc.pos)
#variation with-in and between chain
#This will work if we use different seed for three chains;
gelman.diag(mcmc.pos)
gelman.plot(mcmc.pos)
