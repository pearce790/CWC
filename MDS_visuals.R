Nvals <- c(0,1,2,4,8,16)
dvals <- c(4,8)

for(d in dvals){
  plots <- list()
  for(i in 1:length(Nvals)){
    N <- Nvals[i]
    load(paste0("~/Desktop/Werner_Project/n200_M2_N",N,"_d",d,"_normal_1/data2.RData"))
    distMatrix <- as.matrix(dist(X,diag = TRUE, upper = TRUE))
    plots[[i]]<-ggplot(data=as.data.frame(cmdscale(d = distMatrix, k = 2, eig = FALSE)),
                       aes(V1,V2))+geom_point()+
                labs(x=paste0("Average Distance=",round(mean(dist(X)),3)),
                     title=paste0("Two Standard Gaussians,\nd=",d,"; N=",N))
  }
  do.call("grid.arrange", c(plots, nrow=2,ncol=3))
  ggsave(paste0("MDS_d",d,".png"),
         plot = do.call("grid.arrange", c(plots, nrow=2,ncol=3)),
         width =8, height = 6,scale = 1.5,dpi = 300)
}