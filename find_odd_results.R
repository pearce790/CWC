#### Create Unified Output File ####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
files <- list.files()
files <- files[!grepl("find",files)]
files <- files[grepl("results",files)]
files <- files[grepl("RData",files)]

oddcases <- c()
for(i in 1:length(files)){
  load(files[i])
  if(1 %in% output){oddcases <- c(oddcases,files[i])}
}
oddcases
rm(files)
rm(i)
rm(output)

setwd("/Users/pearce790/Desktop/Werner_Project/normal_n100_nres200_d1_m2_01_81")


source("All_CWC_Functions_20200128.R")

for(i in sample.int(length(oddcases),5)){
  print(i)
  file <- strsplit(oddcases[i],"_")[[1]][2]
  load(file)
  
  #create density
  hs <- cwc.resample.kernel.ve.weights(X,resamples = rs,kmax=0,ngrid = 10,bandwidth = "cv")
  vw.hs <- hs$resample.vertex.weights
  ew.hs <- hs$resample.edge.weights
  
  target.cov.prob <- 0.95
  
  empirical.confidence.bounds <- get.empirical.conf.band(hs,target.cov.prob)
  gaussian.confidence.bounds <-cwc.confidence.bounds(hs,target.cov.prob,simultaneous = T,band.type = "gaussian")
  poisson.confidence.bounds <-cwc.confidence.bounds(hs,target.cov.prob,simultaneous = T,band.type = "poisson")
  bc.confidence.bounds <- get.bc.conf.band(hs,target.cov.prob)
  
  e_res <- get_results(hs,X,empirical.confidence.bounds,pruning.criterion = "bootstrap.rise")
  g_res <- get_results(hs,X,gaussian.confidence.bounds,pruning.criterion = "bootstrap.rise")
  p_res <- get_results(hs,X,poisson.confidence.bounds,pruning.criterion = "bootstrap.rise")
  b_res <- get_results(hs,X,bc.confidence.bounds,pruning.criterion = "bootstrap.rise")
  
  par(mfrow=c(1,2))
  ord_vertex <- order(hs$tg.vertices)
  midpoint_X_of_edges <- (X[hs$tg.edges[,1]] + X[hs$tg.edges[,2]])/2
  ord_edges <- order(midpoint_X_of_edges)
  
  plot(0, 0, xlim = range(hs$tg.vertices), 
       ylim = c(min(0,vw.hs,gaussian.confidence.bounds$lower.vw),
                max(vw.hs,gaussian.confidence.bounds$upper.ew)),
       type = "n",
       xlab="X",ylab="KDE",main=paste0("Gaussian, ",file,"\nResult=",g_res))
  for (i in 1:200) {points(hs$tg.vertices, vw.hs[,i], pch = ".", col = "red",cex=0.01)}
  lines(hs$tg.vertices[ord_vertex], gaussian.confidence.bounds$lower.vw[ord_vertex],
        col = "blue",pch=1,cex=0.5,type="b")
  lines(midpoint_X_of_edges[ord_edges],(gaussian.confidence.bounds$upper.ew)[ord_edges],
        col = "blue",pch=1,cex=0.5,type="p")
  
  plot(0, 0, xlim = range(hs$tg.vertices), 
       ylim = c(min(0,poisson.confidence.bounds$lower.vw,vw.hs),
                max(vw.hs,poisson.confidence.bounds$upper.ew)),
       type = "n",
       xlab="X",ylab="KDE",main=paste0("Poisson, ",file,"\nResult=",p_res))
  for (i in 1:200) {points(hs$tg.vertices, vw.hs[,i], pch = ".", col = "red",cex=0.01)}
  lines(hs$tg.vertices[ord_vertex], poisson.confidence.bounds$lower.vw[ord_vertex],
        col = "blue",pch=1,cex=0.5,type="b")
  lines(midpoint_X_of_edges[ord_edges],(poisson.confidence.bounds$upper.ew)[ord_edges],
        col = "blue",pch=1,cex=0.5,type="p")
  
  # plot(0, 0, xlim = range(hs$tg.vertices), ylim = c(min(0,min(vw.hs)),max(vw.hs)), type = "n",
  #      xlab="X",ylab="KDE",main=paste0("Empirical, ",file,"\nResult=",e_res))
  # for (i in 1:200) {points(hs$tg.vertices, vw.hs[,i], pch = ".", col = "red",cex=0.01)}
  # lines(hs$tg.vertices[ord_vertex], empirical.confidence.bounds$lower.vw[ord_vertex],
  #       col = "blue",pch=1,cex=0.5,type="b")
  # lines(midpoint_X_of_edges[ord_edges],(empirical.confidence.bounds$upper.ew)[ord_edges],
  #       col = "blue",pch=1,cex=0.5,type="p")
  # 
  # plot(0, 0, xlim = range(hs$tg.vertices), ylim = c(min(0,min(vw.hs)),max(vw.hs)), type = "n",
  #      xlab="X",ylab="KDE",main=paste0("Bias-Corrected, ",file,"\nResult=",b_res))
  # for (i in 1:200) {points(hs$tg.vertices, vw.hs[,i], pch = ".", col = "red",cex=0.01)}
  # lines(hs$tg.vertices[ord_vertex], bc.confidence.bounds$lower.vw[ord_vertex],
  #       col = "blue",pch=1,cex=0.5,type="b")
  # lines(midpoint_X_of_edges[ord_edges],(bc.confidence.bounds$upper.ew)[ord_edges],
  #       col = "blue",pch=1,cex=0.5,type="p")
}
