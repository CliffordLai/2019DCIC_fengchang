library(lubridate)
library(magrittr)
library(tidyverse)
library(gridExtra)


#setwd("G:/SCADA_RoundBpre")

for(wtid in 1:33){
  #wtid <- 1 
  tc <- proc.time()
  cat("Processing folder",wtid,"...\n")
  if(!dir.exists("vis")){
    dir.create("vis")
  }
  datapath_vis <- paste0("vis/",dir(path = "dataset/")[[wtid]],"/")
  if(!dir.exists(datapath_vis)){
    dir.create(datapath_vis)
  }
  
  #
  datapath <- paste0("addNA/201807_addNA_",wtid,".rds")
  tr <- readRDS(datapath)
  
  tr$ts <- as_datetime(tr$ts)
  varName <- colnames(tr)[-c(1,2)]
  varClass <- unlist(sapply(tr,class)[-c(1,2)])
  fileName <- paste0(datapath_vis,varName,".png")
  
  cat("Processing variable:")
  for(j in 1:68){
    cat(j,"")
    if(j==30){cat("\n")}
    
    NApos <- which(is.na(tr[,j+2]))
    p0 <- tr %>% ggplot(aes(x=ts,y=eval(parse(text=varName[j])))) +
      ggtitle(paste(varName[j],",",varClass[j],",","NA:",length(NApos),sep=""))+
      geom_point( )+ 
      geom_vline(xintercept = tr$ts[NApos],col="red")
    
    
    # p1 <- autoplot(acf(tr[-NApos,j+2],
    #                    plot=FALSE,lag.max = 100),
    #                main=paste("ACF of",varName[j]))
    # 
    # p2 <- autoplot(acf(diff(as.data.frame(tr[-NApos,j+2])[,1]),
    #                    plot=FALSE,lag.max = 100),
    #                main=paste("ACF of diff",varName[j]))
    
    p1 <- autoplot(acf(tr[,j+2],
                       plot=FALSE,lag.max = 100,na.action=na.pass),
                   main=paste("ACF of",varName[j]))

    p2 <- autoplot(acf(diff(as.data.frame(tr[-NApos,j+2])[,1]),
                       plot=FALSE,lag.max = 100,na.action=na.pass),
                   main=paste("ACF of diff",varName[j]))
    
    #
    p <- grid.arrange(p0,                             # First row with one plot spaning over 2 columns
                      arrangeGrob(p1, p2, ncol = 2), # Second row with 2 plots in 2 different columns
                      nrow = 2)                       # Number of rows
    
    ggsave(fileName[j],plot = p,device="png",width = 10,height = 7)
  }
  tc <- proc.time()-tc
  cat("\n Finish wtid: ",wtid, " in ",tc[3],"s",sep="","\n\n")
}

