#### Set Null Value for Number of Clusters ####
null_numclust <- 1 # number of clusters that is the null value

#### Create Unified Output File ####
dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
(setwd(dir))
files <- list.files()
files <- files[grepl(".RData",files)]
files <- files[!grepl("output",files)]
files <- files[grepl("results",files)]

for(i in 1:length(files)){
  if(i == 1){
    load(files[i])
    output_all <- cbind(output[,1],output[,2:5]>null_numclust)
  }else{
    load(files[i])
    output_all[,2:5] <- output_all[,2:5] + as.numeric(output[,2:5]>null_numclust)
  }
  #print(output)
}
output_all[,2:5] <- output_all[,2:5] / length(files)
colnames(output_all) <- c("ConfLevel","Empirical","Gaussian","Poisson","BiasCorr")
output_all #proportion of datasets with greater than null_numclust clusters identified

foldername <- tail(strsplit(dirname(rstudioapi::getActiveDocumentContext()$path),"/")[[1]],1)
save(output_all,file=paste0("output_",foldername,".RData"))

#### Create Output Plot ####
library(tidyverse)
output_all_long <- pivot_longer(as.data.frame(output_all),
                                cols=c(Empirical,Gaussian,Poisson,BiasCorr))


ggplot(data=output_all_long,aes(ConfLevel,value))+ 
  ylim(0,1) + xlim(0,1) + 
  geom_point(aes(color = factor(name))) +
  geom_line(aes(color = factor(name))) +
  labs(x = "Simultaneous Coverage (Confidence)",y="Power",
       color = "Conf. Band Type",
       title=paste0("Power vs. Confidence Level\n(",strsplit(dir,split="/")[[1]][length(strsplit(dir,split="/")[[1]])],")"))+
  theme(legend.position = "bottom")

foldername <- tail(strsplit(dirname(rstudioapi::getActiveDocumentContext()$path),"/")[[1]],1)
filename = paste0("outputplot_",foldername,".png")
ggsave(filename, plot = last_plot(),
       width =8, height = 7,
       scale = 1,dpi = 300)
