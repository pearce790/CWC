#### Set Working Directory ####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


#### Create Single, Unified Dataframe ####
(folders <- list.files()[c(3,5,9,11,13,7)])
unified_outputs = list()
for(i in 1:length(folders)){
  folder <- folders[i]
  load(paste0(folder,"/output_",folder,".RData"))
  noise_levels <- as.numeric(substr(strsplit(folder,"_")[[1]][3],2,2))
  #mode_dist <- as.numeric(substr(strsplit(folder,"_")[[1]][4],2,2))
  #mode_dist <- as.numeric(substr(folder,nchar(folder)-1,nchar(folder)-1))
  print(c(folder,noise_levels))
  unified_outputs[[i]] <- cbind(noise_levels,output_all)
  #print(c(folder,mode_dist))
  #unified_outputs[[i]] <- cbind(mode_dist,output_all)
}
unified_outputs[[6]][,1] <- 16
unified_outputs <- do.call(rbind,unified_outputs)


#### Create Unified Output Plots ####
library(tidyverse)
library(gridExtra)

unified_outputs_long <- pivot_longer(as.data.frame(unified_outputs),
                                cols=c(Empirical,Gaussian,Poisson,BiasCorr))


#mode_dists <- sort(unique(unified_outputs_long$mode_dist))
noise_levels <- sort(unique(unified_outputs_long$noise_levels))
conf_levels <- sort(unique(unified_outputs_long$ConfLevel))
types <- unique(unified_outputs_long$name)

## Plots of Power by Confidence Level at Difference Bands
plots <- list()
for(i in 1:length(conf_levels)){
  cl <- conf_levels[i]
  sub_data <- subset(unified_outputs_long,unified_outputs_long$ConfLevel == cl)
  #plots[[i]] <- ggplot(data=sub_data,aes(mode_dist,value))+ 
  plots[[i]] <- ggplot(data=sub_data,aes(noise_levels,value))+ 
  #  ylim(0,1) + xlim(0,max(mode_dists)) + 
  ylim(0,1) + xlim(0,max(noise_levels)) + 
    geom_point(aes(color = factor(name))) +
    geom_line(aes(color = factor(name))) +
    labs(x = "Number of Noise Dimensions",y="Power",
    #labs(x = "Distance Between Modes",y="Power",
         color = "Conf. Band Type",
         title=paste0("Power vs. Number of Noise Dimensions\nat Confidence Level = ",cl)) +
         #title=paste0("Power vs. Distance Between Modes\nat Confidence Level = ",cl)) +
    theme(legend.position = "bottom")
}
#plots
#do.call("grid.arrange", c(plots, nrow=2,ncol=2))
grid.arrange(plots[[1]],plots[[3]],plots[[5]],plots[[7]],
             plots[[9]],plots[[10]],nrow=2)
ggsave("powervsconflevel_normal.png",
       plot = grid.arrange(plots[[1]],plots[[3]],plots[[5]],plots[[7]],
                           plots[[9]],plots[[10]],nrow=2),
       width =11, height = 6,
       scale = 1.5,dpi = 300)

## Plots of Power by Band Type at Difference Confidence Levels
plots <- list()
for(i in 1:length(types)){
  b <- types[i]
  sub_data <- subset(unified_outputs_long,unified_outputs_long$name == b)
  #plots[[i]] <- ggplot(data=sub_data,aes(mode_dist,value))+ 
  #  ylim(0,1) + xlim(0,max(mode_dists)) + 
  plots[[i]] <- ggplot(data=sub_data,aes(noise_levels,value))+ 
    ylim(0,1) + xlim(0,max(noise_levels)) + 
    geom_point(aes(color = factor(ConfLevel))) +
    geom_line(aes(color = factor(ConfLevel))) +
    #labs(x = "Distance Between Modes",y="Power",
    labs(x = "Number of Noise Dimensions",y="Power",
         color = "Confidence Level",
         #title=paste0("Power vs. Distance Between Modes\nfor ",b," Confidence Bands")) +
         title=paste0("Power vs. Number of Noise Dimensions\nfor ",b," Confidence Bands")) +
    theme(legend.position = "bottom")
}
#plots
do.call("grid.arrange", c(plots, nrow=2,ncol=2))
ggsave("powervsband_normal.png", plot = do.call("grid.arrange", c(plots, nrow=2,ncol=2)),
       width =8, height = 6,
       scale = 1.5,dpi = 300)


