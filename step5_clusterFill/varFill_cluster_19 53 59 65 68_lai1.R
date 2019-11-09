library(lubridate)
library(magrittr)
library(tidyverse)
library(imputeTS)
library(pracma)
library(zeallot)
library(Rcpp)
library(mgcv)
require("myCfun")

#setwd("G:/SCADA_RoundBpre")
num_wtid <- 33
varcluster <- readRDS(file = "corMatrix/varcluster2.rds")
varcluster

past_score <- readRDS("output_step3/score_wtid_var.rds")
#past_score <- readRDS("output_step3/score_wtid_var_method.rds") #If still running
# past_score <- apply(past_score[1,,],1,max)
ignoreset_history <- which(past_score>0.9)
ignoreset_history


# c(19,53, 59, 65, 68)
rm(list=setdiff(ls(),c("varcluster","num_wtid","ignoreset_history")))
invisible(gc())
source("FE_utilityMar11.R")
source("localApproximateV2fix.R")

clusterset <- c(19,53, 59, 65, 68)
edition <- paste0("roundB",paste(clusterset,collapse = "_"))
pathpre <- paste0("clusterData/","allwtidOnly",paste(clusterset,collapse = "_"),"/")
who <- "lai"
trial <- "1"
wtidSet <- 1:33

outputDir <- paste0("output_step5/fileVar_",edition,"/")
outputDir2 <- paste0("output_step5/fileVarAuto_",edition,"/")
logfile <- paste0("log_step5/record_fileVar_",edition,who,trial,".txt")

ignoreset_history <- setdiff(ignoreset_history,clusterset)

sink(logfile,split=TRUE,append=TRUE)
source("step5_clusterFill/clusterFill.R") #Start filling
sink()

