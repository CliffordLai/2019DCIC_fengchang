library(lubridate)
library(magrittr)
library(tidyverse)
library(zeallot)
library(Rcpp)
require("myCfun")

#
#setwd("G:/SCADA_RoundBpre")
outputPrefix <- "clusterData/"
varcluster <- readRDS(file = "corMatrix/varcluster2.rds")


#
num_wtid <- 33
start_time <- "2018-07-01 00:00:00"

# This requires large memeory and save time reading data
tr_ls <- vector(mode = "list",num_wtid)
time_ls <- vector(mode = "list",num_wtid)
for(wtid in 1:num_wtid){
  cat(wtid,"")
  datapath <- paste0("addNA/201807_addNA_",wtid,".rds")
  tr_ls[[wtid]] <-  as.data.frame(readRDS(datapath),stringsAsFactors=FALSE)
  time_ls[[wtid]] <- tr_ls[[wtid]]$ts
  tr_ls[[wtid]] <- as.matrix(tr_ls[[wtid]][,-c(1,2)])
  if(wtid%%5==0){invisible(gc())}
}
invisible(gc())


# 
timeIndex_ls <- vector(mode = "list",length = num_wtid)
cat("Processing wtid: ")
for(wtid in 1:num_wtid){
  #wtid <- 1
  cat(wtid," ")
  time1 <- as.numeric(as_datetime(time_ls[[wtid]])-as_datetime(start_time),format="second")
  idx_mat <- matrix(NA,length(time_ls[[wtid]]),num_wtid-1)
  h <- 1
  
  for(wtid_other in setdiff(1:num_wtid,wtid)){
    time2 <- as.numeric(as_datetime(time_ls[[wtid_other]])-as_datetime(start_time),format="second")
    idx <- find1in2(time1,time2)
    idxNA <- abs(time1-time2[idx])
    idx[idxNA>5] <- NA #control how accurate the time is matched
    
    idx_mat[,h] <- idx
    h <- 1+h
    if(h%%5==0){invisible(gc())}
    #cat(wtid_other,"")
  }
  timeIndex_ls[[wtid]] <- idx_mat
}
#saveRDS(timeIndex_ls,paste0(outputPrefix,"/timeIndex_ls.rds")) #2.8g large skip saving
invisible(gc())


#timeIndex_ls <- readRDS(paste0(outputPrefix,"/timeIndex_ls.rds"))
genDatabyCluster <- function(clusterset,timeIndex_ls){
  outputpath <- paste0(outputPrefix,"allwtidOnly",paste(clusterset,collapse = "_"),"/")
  tc_var <- proc.time()
  cat("-----------------------------------------------------------------\n")
  cat("Processing wtid: ")
  for(wtid in 1:num_wtid){
    #wtid <- 1
    cat(wtid," ")
    h <- 1
    nc <- length(clusterset)
    tr1 <- matrix(0,nrow(timeIndex_ls[[wtid]]),(num_wtid-1)*nc)
    colnames(tr1) <- paste0(rep(paste0("w",setdiff(1:num_wtid,wtid),"_"),each=length(clusterset)),colnames(tr_ls[[wtid]])[clusterset] )
    
    for(wtid_other in setdiff(1:num_wtid,wtid)){
      tr1[,1:nc+nc*(h-1)] <- tr_ls[[wtid_other]][timeIndex_ls[[wtid]][,h], clusterset,drop=FALSE]
      h <- 1+h
      if(h%%5==0){invisible(gc())}
      #cat(wtid_other,"")
    }
    
    if(!dir.exists(paste0(outputpath))){
      dir.create(paste0(outputpath))
    }
    
    saveRDS(tr1,paste0(outputpath,wtid,".rds"))
    rm(tr1)
    invisible(gc())
  }
  tc_var <- proc.time()-tc_var
  cat("\n")
  cat("Finish clusterset: ",sep="")
  cat(clusterset)
  cat(" in ",tc_var[3],"s",sep="","\n\n")
  return(invisible())
}

tc_all <- proc.time()
#for(i in 1:1){
for(i in 1:length(varcluster)){
  clusterset <- varcluster[[i]]
  genDatabyCluster(clusterset,timeIndex_ls)
}
tc_all <- proc.time()-tc_all
tc_all #17280.15s


