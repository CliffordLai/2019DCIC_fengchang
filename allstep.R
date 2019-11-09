
# Please modify the working directory here before running the script!!!!
setwd("G:\SCADA_sanbengzi")


#--------------------Step 1-------------------#
rm(list=ls())
source("step1_fillNArow2rds.R")
invisible(gc())


#--------------Step 2 (ignore)----------------#
#rm(list=ls())
#source("step2_VIS_SCADA.R")
#invisible(gc())


#--------------------Step 3-------------------#
rm(list=ls())
source("step3_singleVarImpute.R")
invisible(gc())


#--------------------Step 4-------------------#
rm(list=ls())
source("step4a_clusterVar.r")
invisible(gc())

rm(list=ls())
source("step4b_genClusterData.R")
invisible(gc())


#--------------------Step 5-------------------#
rm(list=ls())
source("varFill_cluster_1_5_38_lai.R")
invisible(gc())

rm(list=ls())
source("varFill_cluster_3_24_35_40_45_52_lai1.R")
invisible(gc())

rm(list=ls())
source("varFill_cluster_4_27_34_42_43_46_lai1.R")
invisible(gc())

rm(list=ls())
source("varFill_cluster_19 53 59 65 68_lai1.R")
invisible(gc())


# Other clusters that do not provide much improvment
rm(list=ls())
source("varFill_cluster_8_25_26_lai1.R")
invisible(gc())

rm(list=ls())
source("varFill_cluster_28_36_67_lai1.R")
invisible(gc())

rm(list=ls())
source("varFill_cluster_Double1.R")
invisible(gc())

rm(list=ls())
source("varFill_cluster_Single1.R")
invisible(gc())

rm(list=ls())
source("varFill_cluster_Single2.R")
invisible(gc())

rm(list=ls())
source("varFill_cluster_Single3.R")
invisible(gc())

rm(list=ls())
source("varFill_cluster_6_7_11_14_18_22_29_37_51_55_57_62_lai1.R")
invisible(gc())


#--------------------Step 6-------------------#
rm(list=ls())
source("step6_merge.R")
invisible(gc())



