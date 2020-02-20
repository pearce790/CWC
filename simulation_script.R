source("All_CWC_Functions_20200213.R")

index <- as.numeric(commandArgs(trailingOnly=T)[1])
file <- paste0("data",index,".RData")
load(file)

#create density
hs <- cwc.resample.kernel.ve.weights(X,resamples = rs,
                                     kmax=0,ngrid = 10,bandwidth = "cv")

#get confidence bounds and number of clusters at different simultaneous coverage probs
probs <- c(seq(0.1,0.9,by=0.1),0.95)
output <- matrix(NA,nrow=length(probs),ncol=5)
for(i in 1:length(probs)){
  target.cov.prob <- probs[i]
  empirical.confidence.bounds <- get.empirical.conf.band(hs,target.cov.prob)
  gaussian.confidence.bounds <-cwc.confidence.bounds(hs,target.cov.prob,
                                                     simultaneous = T,band.type = "gaussian")
  poisson.confidence.bounds <-cwc.confidence.bounds(hs,target.cov.prob,
                                                    simultaneous = T,band.type = "poisson")
  bc.confidence.bounds <- get.bc.conf.band(hs,target.cov.prob)
  e_res <- get_results(hs,X,empirical.confidence.bounds,pruning.criterion = "bootstrap.rise")
  g_res <- get_results(hs,X,gaussian.confidence.bounds,pruning.criterion = "bootstrap.rise")
  p_res <- get_results(hs,X,poisson.confidence.bounds,pruning.criterion = "bootstrap.rise")
  b_res <- get_results(hs,X,bc.confidence.bounds,pruning.criterion = "bootstrap.rise")
  
  output[i,] <- c(target.cov.prob,e_res,g_res,p_res,b_res)
}

#show your output
save(output,file=paste0("results_data",index,".RData"))
