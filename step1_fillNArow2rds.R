library(lubridate)
library(magrittr)
library(tidyverse)

col_types <- rep("d",68)
col_types[c(16,20,47,53,66)] <- "i"
col_types = paste0("c","i",paste0(col_types,collapse = ""))
sub <- read_csv("template_submit_result.csv",
                col_types=col_types)

for(wtid in 1:33){
  #wtid <- 1 
  tc <- proc.time()
  cat("Processing wtid",wtid,"...\n")
  datapath <- paste0("dataset/",dir(path = "dataset/")[[wtid]],"/201807.csv")
  
  if(dir.exists(paste0("addNA"))){
    dir.create(paste0("addNA"))
  }
  datapath_addNA <- paste0("addNA/201807_addNA_",wtid,".rds")
  
  tr <- read_csv(datapath,col_types=col_types) 
  total_wtid_ts <- sort(unique(c(sub$ts[sub$wtid==wtid],tr$ts)))
  tr <- as.tbl(data.frame(ts=total_wtid_ts,wtid=wtid,
                          stringsAsFactors = FALSE)) %>%
    left_join(tr[,-2],"ts")
  saveRDS(tr,file =datapath_addNA)
  tc <- proc.time()-tc
  cat("\n Finish wtid: ",wtid, " in ",tc[3],"s",sep="","\n\n")
}
