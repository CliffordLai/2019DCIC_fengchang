library(lubridate)
library(magrittr)
library(tidyverse)
library(gridExtra)

#setwd("G:/SCADA_RoundBpre")
num_wtid <- 33


cor_wtid <- vector(mode = "list",length = num_wtid)
for(wtid in 1:33){
  #wtid <- 1 
  tc <- proc.time()
  cat("wtid ",wtid,", ",sep="")
  datapath <- paste0("addNA/201807_addNA_",wtid,".rds")
  tr <- as.matrix(readRDS(datapath)[,-c(1,2)])
  cat("ratio",mean(complete.cases(tr)))
  cor_wtid[[wtid]] <- cor(tr,use = "complete.obs")
  #cor_wtid[[wtid]] <- cor(tr,use = "pairwise.complete.obs")

  tc <- proc.time()-tc
  cat(" in ",tc[3],"s",sep="","\n")
}

saveRDS(cor_wtid,"corMatrix/cor_wtid.rds")


cor_mean <- matrix(0,68,68)
for(wtid in 1:33){
  cor_mean_tmp <- cor_wtid[[wtid]]
  na_pos <- is.na(cor_mean_tmp)
  cor_mean_tmp[na_pos] <- 0
  cor_mean <- cor_mean_tmp+cor_mean
}
cor_mean <- cor_mean/num_wtid

saveRDS(cor_mean,"corMatrix/cor_mean.rds")



cor_mean <- readRDS("corMatrix/cor_mean.rds")

# Visualize clusters
vis_cor_mean <- cor_mean
colnames(vis_cor_mean) <- 1:68
rownames(vis_cor_mean) <- 1:68
hc <- hclust(as.dist((1 - abs(vis_cor_mean))))
plot(hc)
abline(h=0.26,col="red",lty=2)

pdf("ClusterByCor.pdf",width = 14,height=6)
plot(hc)
abline(h=0.26,col="red",lty=2)
dev.off()


# Save clusters
varcluster <- cutree(hc, h=0.26)
sort(varcluster)
max(varcluster)

varcluster2 <- lapply(1:max(varcluster),FUN = function(i){as.numeric(names(varcluster)[varcluster==i])})
saveRDS(varcluster2,file = "corMatrix/varcluster2.rds")

varcluster2 <- readRDS(file = "corMatrix/varcluster2.rds")
sink("varcluster2.txt")
print(varcluster2)
sink()
