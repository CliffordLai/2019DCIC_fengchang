library(lubridate)
library(magrittr)
library(tidyverse)

#setwd("G:/SCADA_RoundBpre")

cat("Find the range for each variable: ")
exactvarRange <- array(0,dim=c(33,68,2))
for(wtid in 1:33){
  cat(wtid,"")
  data_wtid <- as.data.frame(readRDS(paste0("addNA/201807_addNA_",wtid,".rds"))[,-c(1,2)])
  exactvarRange[wtid,,] <- t(apply(data_wtid,2,range,na.rm=TRUE))
}
saveRDS(exactvarRange,"exactvarRange.rds")

#
col_types <- rep("d",68)
col_types[c(16,20,47,53,66)] <- "i"
col_types = paste0("c","i",paste0(col_types,collapse = ""))

#----------------Add single series method, used as baseline-------------------------#
sub <- read_csv("template_submit_result.csv",
                col_types=col_types)   #
subts <- sub[,c(1,2)] %>% mutate(obs=1:nrow(sub))
n_replace <- 0
for(wtid in 1:33){
  cat(wtid,"")
  fillpath <- paste0("output_step3/",wtid,"/")
  varName <- dir(fillpath)
  varIdx <-  strsplit(varName,"_|\\.")
  varIdx <- as.numeric(sapply(varIdx,FUN = function(x)x[2]))
  
  for(j in 1:length(varName)){
    trfill <- readRDS(paste0(fillpath,varName[j]))
    trfill <- trfill %>% left_join(subts,c("ts","wtid"))
    sub[trfill$obs,varIdx[j]+2] <- trfill$var
    n_replace <- n_replace+nrow(trfill)
  }
}
n_replace
any(is.na(sub))
sub[is.na(sub)] <- 0

write_csv(sub,path=paste0("submit/","autoPredict.csv"),na="")  #0.67929965000
#-------------------------------------------------------------------------#


#-------------------------------------------------------------------------#
#+ 1,5,38 
sub <- read_csv(paste0("submit/","autoPredict.csv"),
                col_types=col_types)   #
subts <- sub[,c(1,2)] %>% mutate(obs=1:nrow(sub))
n_replace <- 0
for(wtid in 1:33){
  fillpath <- paste0("output_step5/fileVar_roundB1_5_38/tuneSet11/",wtid,"/")
  varName <- dir(fillpath)
  if(length(varName)==0){
    cat("Miss ")
    next
  }
  cat(wtid,"")
  varIdx <-  strsplit(varName,"_|\\.")
  varIdx <- as.numeric(sapply(varIdx,FUN = function(x)x[2]))
  
  for(j in 1:length(varName)){
    trfill <- readRDS(paste0(fillpath,varName[j]))
    trfill <- trfill %>% left_join(subts,c("ts","wtid"))
    sub[trfill$obs,varIdx[j]+2] <- trfill$var
    n_replace <- n_replace+nrow(trfill)
  }
}
n_replace


#-------------------------------------------------------------------------#
#+4_27_34_42_43_46 
n_replace <- 0
for(wtid in 1:33){
  fillpath <- paste0("output_step7/fileVar_roundB4_27_34_42_43_46/tuneSet92/",wtid,"/")
  varName <- dir(fillpath)
  if(length(varName)==0){
    cat("Miss ")
    next
  }
  cat(wtid,"")
  varIdx <-  strsplit(varName,"_|\\.")
  varIdx <- as.numeric(sapply(varIdx,FUN = function(x)x[2]))
  
  for(j in 1:length(varName)){
    trfill <- readRDS(paste0(fillpath,varName[j]))
    trfill <- trfill %>% left_join(subts,c("ts","wtid"))
    sub[trfill$obs,varIdx[j]+2] <- trfill$var
    n_replace <- n_replace+nrow(trfill)
  }
}
n_replace


#-------------------------------------------------------------------------#
#+19_53_59_65_68 
n_replace <- 0
for(wtid in 1:33){
  fillpath <- paste0("output_step7/fileVar_roundB19_53_59_65_68/tuneSet97/",wtid,"/")
  varName <- dir(fillpath)
  if(length(varName)==0){
    cat("Miss ")
    next
  }
  cat(wtid,"")
  varIdx <-  strsplit(varName,"_|\\.")
  varIdx <- as.numeric(sapply(varIdx,FUN = function(x)x[2]))
  
  for(j in 1:length(varName)){
    trfill <- readRDS(paste0(fillpath,varName[j]))
    trfill <- trfill %>% left_join(subts,c("ts","wtid"))
    sub[trfill$obs,varIdx[j]+2] <- trfill$var
    n_replace <- n_replace+nrow(trfill)
  }
}
n_replace
# sub[,c(c(16,20,47,53,66)+2)] <- round(sub[,c(c(16,20,47,53,66)+2)])
unique(sub$var053)
sub$var053 <- ifelse(sub$var053<0.5,0,1)


# Fix range
exactvarRange <- readRDS("exactvarRange.rds")
minRange <- apply(exactvarRange[,,1],2,min)
maxRange <- apply(exactvarRange[,,2],2,max)

for(j in 1:68){
  cat("var",j,":",sep="")
  exceed_min <- 0
  exceed_max <- 0
  
  y <- pull(sub,j+2)
  eminIdx <- (y < minRange[j])
  emaxIdx <- (y > maxRange[j])
  
  if(any(eminIdx,na.rm = TRUE)){
    exceed_min <- exceed_min+sum(eminIdx,na.rm = TRUE)
    sub[which(eminIdx),j+2] <- minRange[j]
  }
  
  if(any(emaxIdx,na.rm = TRUE)){
    exceed_max <- exceed_max+sum(emaxIdx,na.rm = TRUE)
    sub[which(emaxIdx),j+2] <- maxRange[j]
  }

  
  cat("exceed min:",exceed_min,", ")
  cat("exceed max:",exceed_max," ")
  cat("\n")
  if(j %in% c(16,20,47,53,66)){
    sub[,j+2] <- round(sub[,j+2])
  }else{
    sub[,j+2] <- round(sub[,j+2]*100)/100
  }
}


write_csv(sub,path=paste0("submit/","sub_part.csv"),na="")  



#-------------------------------------------------------------------------#

#+ all the other clusters 
sub <- read_csv(paste0("submit/","sub_part.csv"),
                col_types=col_types)   #
subts <- sub[,c(1,2)] %>% mutate(obs=1:nrow(sub))

varcluster <- readRDS(file = "corMatrix/varcluster2.rds")
varcluster <- varcluster[-c(1,3,4,14)]

n_replace <- 0
for(i in 1:length(varcluster)){
  clusterset <- c(6,7,11,14,18,22,29,37,51,55,57,62)
  edition <- paste0("roundB",paste(clusterset,collapse = "_"))
  outputDir <- paste0("output_step7/fileVar_",edition,"/")
  
  if(dir.exists(outputDir)){
    cat(edition,":",sep="")
    for(wtid in 1:33){
      fillpath <- paste0(outputDir,"tuneSet120/",wtid,"/")
      varName <- dir(fillpath)
      if(length(varName)==0){
        next
      }
      cat(wtid,"")
      varIdx <-  strsplit(varName,"_|\\.")
      varIdx <- as.numeric(sapply(varIdx,FUN = function(x)x[2]))
      
      for(j in 1:length(varName)){
        trfill <- readRDS(paste0(fillpath,varName[j]))
        trfill <- trfill %>% left_join(subts,c("ts","wtid"))
        sub[trfill$obs,varIdx[j]+2] <- trfill$var
        n_replace <- n_replace+nrow(trfill)
      }
    }
  }
}
n_replace

# Fix range
exactvarRange <- readRDS("exactvarRange.rds")
minRange <- apply(exactvarRange[,,1],2,min)
maxRange <- apply(exactvarRange[,,2],2,max)

for(j in 1:68){
  cat("var",j,":",sep="")
  exceed_min <- 0
  exceed_max <- 0
  
  y <- pull(sub,j+2)
  eminIdx <- (y < minRange[j])
  emaxIdx <- (y > maxRange[j])
  
  if(any(eminIdx,na.rm = TRUE)){
    exceed_min <- exceed_min+sum(eminIdx,na.rm = TRUE)
    sub[which(eminIdx),j+2] <- minRange[j]
  }
  
  if(any(emaxIdx,na.rm = TRUE)){
    exceed_max <- exceed_max+sum(emaxIdx,na.rm = TRUE)
    sub[which(emaxIdx),j+2] <- maxRange[j]
  }
  
  
  cat("exceed min:",exceed_min,", ")
  cat("exceed max:",exceed_max," ")
  cat("\n")
  if(j %in% c(16,20,47,53,66)){
    sub[,j+2] <- round(sub[,j+2])
  }else{
    sub[,j+2] <- round(sub[,j+2]*100)/100
  }
}


write_csv(sub,path=paste0("submit/","sub_fullcluster.csv"),na="")  #