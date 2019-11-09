library(lubridate)
library(magrittr)
library(tidyverse)
library(imputeTS)
library(pracma)
library(zeallot)
library(Rcpp)
library(mgcv)
require(myCfun)  #my own package

#
#setwd("G:/SCADA_RoundBpre")
outputDir <- "output_step3"
logDir <- "log_step3"
logfile <- paste0(logDir,"/autolog1.txt")

# My own functions
source("FE_utilityMar11.R")

#
if(!dir.exists(outputDir)){
  dir.create(outputDir)
}
if(!dir.exists(logDir)){
  dir.create(logDir)
}

sink(logfile,split=TRUE,append=TRUE)


num_wtid <- 33
nearest_K <- c(2,5,7,10,15,20,25,50,100,200,500)
#number of methods
name_method <- c("mode","nearest","linear",
                 paste0(c("WMA","WMO","LWA","LWO"),rep(nearest_K,each=4)))
length(name_method)
score_wtid_var_method <- array(0,dim=c(num_wtid,68,length(name_method)))


set.seed(111)
for(wtid in 1:num_wtid){
  #wtid <- 1
  tc_wtid <- proc.time()
  cat("-----------------------------------------------------------------\n")
  cat("Processing wtid:",wtid,"\n")
  
  datapath <- paste0("addNA/201807_addNA_",wtid,".rds")
  tr_wtid <- readRDS(datapath)
  
  for(var_j in 1:68){
    tc_var <- proc.time()
    #var_j <- 1
    y0 <- pull(tr_wtid, var_j+2)
    n_y0 <- length(y0)
    isna_y0 <- is.na(y0)
    myNArow <- simNArow(isna_y0,n_y0)
    n_myNArow <- length(myNArow)
    
    # Make my own validation set with known missing values
    y1 <- y0
    y1[myNArow] <- NA
    na_pos_y1 <- which(is.na(y1))
    
    # Save the auto result
    autopred_trainval <- matrix(0,nrow=n_myNArow, ncol=length(name_method))
    
    # Mode
    autopred_trainval[,1] <- rep(Mode(y0[!isna_y0]),n_myNArow)
    
    # Use nearest to compare
    x0 <- 1:n_y0
    autopred_trainval[,2] <- interp1(x0[-na_pos_y1], y1[-na_pos_y1], 
                                     xi = x0[myNArow], method = "nearest")
    
    # linear
    autopred_trainval[,3] <- na.interpolation(y1)[myNArow]
    
    # Prepare for WMA WMO LWA LWO
    tr_WM <- trainWMA(y1,na_pos_y1,kmax=max(nearest_K))
    tmp_idx <- na_pos_y1  %in% myNArow
    for(ki in 1:length(nearest_K)){
      autopred_trainval[,4:7+4*(ki-1)] <- predictWM(tr_WM,k = nearest_K[ki],
                                                    WM = c("WMA","WMO","LWA","LWO"))[tmp_idx,]
    }
    
    if(var_j %in% c(16,20,47,53,66)){
      autopred_trainval <- round(autopred_trainval)
    }
    
    #use my missing set as validation set
    autoScore <- apply(autopred_trainval,2,FUN = function(x){
      mean(eLBscore(y0[myNArow],x,var_j))
    })
    
    score_wtid_var_method[wtid,var_j,] <- autoScore
    
    cat("\nScores:","\n")
    print(data.frame(method=name_method,score=autoScore))
    cat("optimal:",name_method[which.max(autoScore)],
        max(autoScore) )
    tc_var <- proc.time()-tc_var
    cat("\nFinish var ",var_j," of witd ",wtid," in ",tc_var[3],"s",sep="","\n\n")
    
    if(var_j %% 10==0){
      rm(autopred_trainval)
      invisible(gc())
    } 
  }
  
  tc_wtid <- proc.time()-tc_wtid
  cat("\nFinish wtid ",wtid," in ",tc_wtid[3],"s",sep="","\n\n")
  
  # Save the result in case crash
  saveRDS(score_wtid_var_method,paste0(outputDir,"/score_wtid_var_method.rds"))
  rm(autopred_trainval)
  invisible(gc())
}


#Best method for each varaible within each wtid
score_wtid_var <- numeric(68)
for(j in 1:68){
  score_wtid_var[j] <- mean(apply(score_wtid_var_method[,j,],1,max))
  cat("--------------\nvar",j,": ", score_wtid_var[j],sep="")
  print(table(name_method[apply(score_wtid_var_method[,j,],1,which.max)]))
  cat("\n")
}
matrix(score_wtid_var,ncol=1)
mean(score_wtid_var)
saveRDS(score_wtid_var,paste0(outputDir,"/score_wtid_var.rds"))

# Best method for each variable over all wtid 
score_var_method <- apply(score_wtid_var_method,c(2,3),mean)

homogenous_var_score <- apply(score_var_method,1,max)
saveRDS(homogenous_var_score,"output_step3/homogenous_var_score.rds")

data.frame(method=matrix(name_method[apply(score_var_method,1,which.max)],ncol=1),score=homogenous_var_score)
mean(homogenous_var_score) #0.6794363

sink()



# post-processing to fill the test set
for(wtid in 1:num_wtid){
  #wtid <- 1
  tc_wtid <- proc.time()
  cat("-----------------------------------------------------------------\n")
  cat("Processing wtid:",wtid,"in ")
  
  datapath <- paste0("addNA/201807_addNA_",wtid,".rds")
  tr_wtid <- readRDS(datapath)
  
  for(var_j in 1:68){
    tc_var <- proc.time()
    #var_j <- 1
    y0 <- pull(tr_wtid, var_j+2)
    n_y0 <- length(y0)
    isna_y0 <- is.na(y0)
    
    choose_method <- which.max(score_var_method[var_j,])
    
    if(choose_method==1){
      nafill <- rep(Mode(y0[!isna_y0]),sum(isna_y0))
    }else if(choose_method==2){
      x0 <- 1:n_y0
      nafill <- interp1(x0[!isna_y0], y0[!isna_y0], 
                        xi = x0[isna_y0], method = "nearest")
    }else if(choose_method==3){
      nafill <- na.interpolation(y0)[isna_y0]
    }else{
      if((choose_method-3)%%4==1){
        choose_k <- nearest_K[ceiling((choose_method-3)/4)]
        tr_WM <- trainWMA1(y0,which(isna_y0),kmax=choose_k)
        nafill <- predictWM(tr_WM,k = choose_k,
                            WM = c("WMA"))
      }else if((choose_method-3)%%4==2){
        choose_k <- nearest_K[ceiling((choose_method-3)/4)]
        tr_WM <- trainWMA1(y0,which(isna_y0),kmax=choose_k)
        nafill <- predictWM(tr_WM,k = choose_k,
                            WM = c("WMO"))
      }else if((choose_method-3)%%4==3){
        choose_k <- nearest_K[ceiling((choose_method-3)/4)]
        tr_WM <- trainWMA2(y0,which(isna_y0),kmax=choose_k)
        nafill <- predictWM(tr_WM,k = choose_k,
                            WM = c("LWA"))
      }else if((choose_method-3)%%4==0){
        choose_k <- nearest_K[ceiling((choose_method-3)/4)]
        tr_WM <- trainWMA2(y0,which(isna_y0),kmax=choose_k)
        nafill <- predictWM(tr_WM,k = choose_k,
                            WM = c("LWO"))
      }
    }
    
    if(var_j %in% c(16,20,47,53,66)){
      nafill <- round(nafill)
    }else{
      nafill <- round(nafill*100)/100
    }
    
    
    # Save result
    predtestAuto <- data.frame(tr_wtid[isna_y0,1:2],var=nafill)
    if(!dir.exists(paste0(outputDir2,"/",wtid))){
      dir.create(paste0(outputDir2,"/",wtid))
    }
    fillvar_name <- paste0(outputDir2,"/",wtid,"/var_",var_j,".rds")
    saveRDS(predtestAuto,file = fillvar_name)
  }
  
  tc_wtid <- proc.time()-tc_wtid
  cat(tc_wtid[3],"s",sep="","\n")
  invisible(gc())
}









