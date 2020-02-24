#### Create Data ####
setwd("~/Desktop/Werner_Project/Files for New Simulation")
source('All_CWC_Functions_20200213.R', echo=FALSE)
X <- generate_simplex_data(n=5000,M=2,N=100,d=5,distribution="normal",scale=1)

#### Plot Data Together
par(mfrow=c(2,3))
for(i in c(3,7,12,22,52,102)){
  signal_only <- X[,1:2]
  pca_transformed <- X[,1:i] %*% prcomp(X[,1:i])$rotation[,1:2]
  procrusted_data <- procrustes(signal_only,pca_transformed)
  plot(procrusted_data$X,pch=16,cex=.2,xlab="",ylab="",
       main=paste0("Centered Signal Data,\nProcrusted Transform of ",
                   i-2,"\nAdded Noise Dimensions"))
  points(procrusted_data$Yrot,col="blue",pch=16,cex=.2)
  legend("topright",legend=c("Signal","Noisy"),col=c("black","blue"),pch=16)
}

#### Plot Data Side-by-Side
par(mfrow=c(3,2))
plot(procrusted_data$X,pch=16,cex=.3,xlab="",ylab="",
     main=paste0("Centered Signal Data"))
for(i in c(3,7,12,52,102)){
  signal_only <- X[,1:2]
  pca_transformed <- X[,1:i] %*% prcomp(X[,1:i])$rotation[,1:2]
  procrusted_data <- procrustes(signal_only,pca_transformed)
  plot(procrusted_data$Yrot,col="blue",pch=16,cex=.3,xlab="",ylab="",
       main=paste0("Procrusted Transform of ",i-2,"\nAdded Noise Dimensions"))
}
