#### CREATE DATA ####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('All_CWC_Functions_20200213.R', echo=FALSE)
dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
dir <- strsplit(dir,"/")[[1]][length(strsplit(dir,"/")[[1]])]
(dir <- strsplit(dir,"_")[[1]])


(n <- as.numeric(gsub("[^\\d]+", "", dir[1], perl=TRUE)))
(M <- as.numeric(gsub("[^\\d]+", "", dir[2], perl=TRUE)))
(N <- as.numeric(gsub("[^\\d]+", "", dir[3], perl=TRUE)))
(d <- as.numeric(gsub("[^\\d]+", "", dir[4], perl=TRUE)))
(distribution <- dir[5])
(scale <- as.numeric(gsub("[^\\d]+", "", dir[6], perl=TRUE)))

nres <- 1000
num_datasets <- 500


f <- function(x){sample.int(n,n/2,replace = FALSE)}
for(index in 1:num_datasets){
  X <- generate_simplex_data(n,M,N,d,distribution,scale)
  rs <- matrix(unlist(lapply(1:nres,f)),nrow=n/2,ncol=nres)
  save(X,rs,file=paste0("data",index,".RData"))
}
