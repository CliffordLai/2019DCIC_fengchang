
#options(warn=2) #investigate warnings
#options(warn=1) #restore

global_impute <- FALSE

eps_samewtid <- sort(c(-0.2,-0.15,-0.1,-0.05,-0.01,-0.005,-0.001, 0 ,0.001,0.005,0.01,0.05,0.1,0.15,0.2,0.25))
eps_diffwtid <- sort(c(-0.2,-0.15,-0.1,-0.05,-0.01,-0.005,-0.001, 0 ,0.001,0.005,0.01,0.05,0.1,0.15,0.2,0.25))

L <- 3
E <- 2

#
eps_tuneset <- expand.grid(same=eps_samewtid,diff=eps_diffwtid)
eps_tuneset <- eps_tuneset[eps_tuneset$same<=eps_tuneset$diff,]
eps_tuneset <- eps_tuneset %>% arrange(same,diff)
eps_tuneset
tuneI <- nrow(eps_tuneset)


########
wtid_var_improvment_ls <- lapply(1:num_wtid,FUN = function(x){lapply(1:68,function(x){})})
wtid_var_aimp_ls <- lapply(1:num_wtid,FUN = function(x){lapply(1:68,function(x){})})
wtid_var_aimpp_ls <- lapply(1:num_wtid,FUN = function(x){lapply(1:68,function(x){})})

ignoreset_history_tmp <- setdiff(ignoreset_history,clusterset)
ignore_clustervar_tmp <- intersect(ignoreset_history,clusterset)
if(length(ignore_clustervar_tmp)!=0){
  newclusterIdx <- clusterset %in% ignoreset_history
  clusterset <- clusterset[newclusterIdx]
}

#################
for(wtid in wtidSet){
  #wtid <- 1
  tc_var <- proc.time()
  cat("-----------------------------------------------------------------\n")
  cat("Processing wtid:",wtid,"\n")
  tr_wtid <-  as.data.frame(readRDS(paste0("addNA/201807_addNA_",wtid,".rds")),stringsAsFactors=FALSE)
  tr_othervar <- readRDS(paste0(pathpre,wtid,".rds"))
  if(length(ignore_clustervar_tmp)!=0){
    tr_othervar <- tr_othervar[,rep(newclusterIdx, 32)]
  }
  invisible(gc())
  
  stat_improveScore_wtid <- 0
  stat_aimpp_wtid <- 0
  num_fill_wtid <- 0
  tune_var_improve <- vector(mode = "list",tuneI) 
  
  for(oset in 1:length(clusterset)){
    na_var <- clusterset[oset]
    cat("------------------------------------\n","var:",na_var,", ",sep="")
    y0 <- tr_wtid[,2+na_var]
    fill_var <- (1:68)[-c(na_var,ignoreset_history_tmp)]
    num_withinwtid <- length(fill_var)
    fill_var <- c( fill_var, clusterset+rep(seq(68,68*(num_wtid-1),by=68),each=length(clusterset)))
    cname <- c(colnames(tr_wtid)[-c(1,2,2+na_var,2+ignoreset_history_tmp)],colnames(tr_othervar))
    na_margin_y0 <- findMargin(y0)
   
    # create NA for each region 
    isna_pos_y0 <- (is.na(y0))
    testNArow <- which(isna_pos_y0)
    region_NA <- simNArow5(isna_pos_y0,n=length(isna_pos_y0),L=L,E=E)

    #
    if(global_impute){
      n_region <- 1
    }else{
      region_NA <- region_NA$region_simNA
      n_region <- length(region_NA)
      cat("number of segments:",n_region,"\n\n")
    }
    
    y0_mode <- Mode(y0[!isna_pos_y0])
    seg_result <- vector(mode = "list",length = n_region)
    for(seg in 1:n_region){
      if(global_impute){
        myNArow1 <- region_NA$myNA_rowIdx1
        myNArow2 <- region_NA$myNA_rowIdx2
        testSegNArow <- region_NA$testNArow
        
        myNArow1 <- unique(setdiff(myNArow1,testSegNArow))
        myNArow2 <- unique(setdiff(myNArow2,testSegNArow))
        myNArow <- unique(setdiff(c(myNArow1,myNArow2),testSegNArow))
      }else{
        myNArow <- region_NA[[seg]]$myNA_rowIdx
        myNArow1 <- region_NA[[seg]]$myNA_rowIdx1
        myNArow2 <- region_NA[[seg]]$myNA_rowIdx2
        testSegNArow <- region_NA[[seg]]$testNArow
      }
    
      #test set
      y0_testseg <- testNArow %in% testSegNArow
      na_margin_y0seg <- na_margin_y0[y0_testseg]
      
      #training set
      y1 <- y0
      y1[myNArow1] <- NA
      na_pos_y1 <- which(is.na(y1))
      y1_myNArow <- na_pos_y1 %in% myNArow1
      na_margin_y1 <- findMargin(y1)[y1_myNArow]
    
      #validation set
      y2 <- y0
      y2[myNArow2] <- NA
      na_pos_y2 <- which(is.na(y2))
      y2_myNArow <- na_pos_y2 %in% myNArow2
      na_margin_y2 <- findMargin(y2)[y2_myNArow]
      

      #---- Start auto method -------------------------------------------#
      # Create local series for training auto method
      n_seg1 <- length(myNArow1)
      n_seg2 <- length(myNArow2)
      aug_segRow <- range(myNArow)
      aug_segRow[1] <- aug_segRow[1]-60
      aug_segRow[2] <- aug_segRow[2]+60
      if(aug_segRow[1]<1){
        aug_segRow[1] <- 1
      }
      while(is.na(y0[aug_segRow[1]])){
        aug_segRow[1] <- aug_segRow[1]-60
        if(aug_segRow[1]<1){
          aug_segRow[1] <- 1
        }
      }
      while(is.na(y0[aug_segRow[2]])){
        aug_segRow[2] <- aug_segRow[2]+60
        if(aug_segRow[2]>length(y0)){
          aug_segRow[2] <- length(y0)
        }
      }
      if(aug_segRow[2]>length(y0)){
        aug_segRow[2] <- length(y0)
      }
      aug_y0 <- y0[aug_segRow[1]:aug_segRow[2]]
      na_pos_aug_y0 <- which(is.na(aug_y0))
      aug_myNArow1 <- myNArow1-aug_segRow[1]+1
      aug_myNArow2 <- myNArow2-aug_segRow[1]+1    
      aug_myNArow <- myNArow-aug_segRow[1]+1
      aug_testSegNArow <- testSegNArow-aug_segRow[1]+1    
      
      aug_y0NA <- aug_y0
      aug_y0NA[c(aug_myNArow1,aug_myNArow2)] <- NA
      na_pos_aug_y0NA <- which(is.na(aug_y0NA))
      
      # Save the auto result
      autopred_trainval <- matrix(0,nrow=n_seg1+n_seg2, ncol=9)

      # Use nearest to compare
      aug_x0 <- 1:length(aug_y0)
      autopred_trainval[,1] <- interp1(aug_x0[-na_pos_aug_y0NA], aug_y0NA[-na_pos_aug_y0NA], 
                                       xi = aug_x0[c(aug_myNArow1,aug_myNArow2)], method = "nearest")
      
      # Mode
      autopred_trainval[,2] <- rep(y0_mode,n_seg1+n_seg2)
      
      # linear
      autopred_trainval[,3] <- na.interpolation(aug_y0NA)[c(aug_myNArow1,aug_myNArow2)]
      
      # Prepare for WMA WMO LWA LWO
      tr_WM <- trainWMA1(aug_y0NA,na_pos_aug_y0NA,kmax=50)
      tmp_idx <- na_pos_aug_y0NA  %in% c(aug_myNArow1,aug_myNArow2)
      autopred_trainval[,4:5] <- predictWM(tr_WM,k = 5,WM = c("WMA","WMO"))[tmp_idx,]
      autopred_trainval[,6:7] <- predictWM(tr_WM,k = 15,WM = c("WMA","WMO"))[tmp_idx,]
      autopred_trainval[,8:9] <- predictWM(tr_WM,k = 50,WM = c("WMA","WMO"))[tmp_idx,]
      
      #use the second set as validation set
      autoScore <- apply(autopred_trainval[1:n_seg2+n_seg1,,drop=FALSE],2,FUN = function(x){
        mean(eLBscore(aug_y0[c(aug_myNArow2)],x,na_var))
      })
      
      autobestIndex <- which.max(autoScore)
      autobestScore <- autoScore[autobestIndex]
      eautobestScore <- eLBscore(aug_y0[c(aug_myNArow2)],autopred_trainval[1:n_seg2+n_seg1,autobestIndex],na_var)
      ypred_auto <- autopred_trainval[1:n_seg2+n_seg1,autobestIndex]
      
      if(autobestIndex==1){
        ypred_auto_test <- interp1(aug_x0[-na_pos_aug_y0NA], aug_y0NA[-na_pos_aug_y0NA], 
                                   xi = aug_x0[c(aug_testSegNArow)], method = "nearest")
        
      }else if(autobestIndex==2){
        ypred_auto_test <- rep(y0_mode,length(testSegNArow))
        
      }else if(autobestIndex==3){
        ypred_auto_test <- na.interpolation(aug_y0NA)[aug_testSegNArow]
        
      }else if(autobestIndex==4){
        tr_WM <- trainWMA1(aug_y0,na_pos_aug_y0,kmax=5)
        tmp_idx <- na_pos_aug_y0  %in% aug_testSegNArow
        ypred_auto_test <- predictWM(tr_WM,k = 5,WM = c("WMA"))[tmp_idx,]
        
      }else if(autobestIndex==5){
        tr_WM <- trainWMA1(aug_y0,na_pos_aug_y0,kmax=5)
        tmp_idx <- na_pos_aug_y0  %in% aug_testSegNArow
        ypred_auto_test <- predictWM(tr_WM,k = 5,WM = c("WMO"))[tmp_idx,]
        
      }else if(autobestIndex==6){
        tr_WM <- trainWMA1(aug_y0,na_pos_aug_y0,kmax=15)
        tmp_idx <- na_pos_aug_y0  %in% aug_testSegNArow
        ypred_auto_test <- predictWM(tr_WM,k = 15,WM = c("WMA"))[tmp_idx,]
        
      }else if(autobestIndex==7){
        tr_WM <- trainWMA1(aug_y0,na_pos_aug_y0,kmax=15)
        tmp_idx <- na_pos_aug_y0  %in% aug_testSegNArow
        ypred_auto_test <- predictWM(tr_WM,k = 15,WM = c("WMO"))[tmp_idx,]
      }else if(autobestIndex==8){
        tr_WM <- trainWMA1(aug_y0,na_pos_aug_y0,kmax=50)
        tmp_idx <- na_pos_aug_y0  %in% aug_testSegNArow
        ypred_auto_test <- predictWM(tr_WM,k = 50,WM = c("WMA"))[tmp_idx,]
        
      }else if(autobestIndex==9){
        tr_WM <- trainWMA1(aug_y0,na_pos_aug_y0,kmax=50)
        tmp_idx <- na_pos_aug_y0  %in% aug_testSegNArow
        ypred_auto_test <- predictWM(tr_WM,k = 50,WM = c("WMO"))[tmp_idx,]
      }
      
      # Probably not useful as these validation sets are not right
      # they are only used for comparing with regression models
      predtestAuto <- data.frame(tr_wtid[testSegNArow,1:2],var=ypred_auto_test)
      rm(aug_y0,aug_y0NA,na_pos_aug_y0,aug_myNArow1,aug_myNArow2,aug_myNArow,aug_testSegNArow,na_pos_aug_y0NA)
      invisible(gc())
      
      #----------Start cross method ------------------------------#
      improveScore_ls <- numeric(length(fill_var))
      aimproveScore_ls <- numeric(length(fill_var))
      valScore_ls <- numeric(length(fill_var))
      caseE_ls <- rep(NA,length(fill_var))
      cut_lag <- rep(NA,length(fill_var))
      caseLS_ls <- rep(NA,length(fill_var))
      Er_ls <- numeric(length(fill_var))
      nafill_test_ls <- matrix(0,nrow = length(testSegNArow),ncol = length(fill_var))
      tune_record_indicator <- matrix(FALSE,tuneI,length(fill_var))
      
      for(i in 1:length(fill_var)){
        #linear model to predict
        if(fill_var[i]>68){
          y_fill_var <- tr_othervar[,i-num_withinwtid]
        }else{
          y_fill_var <- tr_wtid[,2+fill_var[i]]
        }
        
        #If no observation of fill_var in the validation/test set
        if(all(is.na(y_fill_var[myNArow1]))){        
          next
        }
        if(all(is.na(y_fill_var[myNArow2]))){        
          next
        }
        if(all(is.na(y_fill_var[testSegNArow]))){
          next
        }
        
        ignore_valobs <- is.na(y_fill_var[myNArow2])
        
        
        # Single method
        score0 <- numeric(3)
        fit1 <- vector(mode = "list",length = 3)
        mymethod <- c("myls","ols","gam")
        for(k in 1:3){
          nafillval <- fittedVal(mymethod[k])
          score0[k] <- LBscore(y0[myNArow2][!ignore_valobs],nafillval[!ignore_valobs],na_var)
          fit1[[k]] <- optimize(f = function(r_lag){
            nafillval_cut <- fittedVal_cut(r_lag, nafillval)
            LBscore(y0[myNArow2][!ignore_valobs],nafillval_cut[!ignore_valobs],na_var)
            mean(eautobestScore[!ignore_valobs])
          },interval = c(0,1),maximum = TRUE)
        }
        
        caseLS <- which.max(c(score0,sapply(fit1,function(x){x[[2]]})))
        if(length(caseLS)==0){next}
        if(caseLS<=3){
          caseLS <- caseLS
          r_lag <- 0
          cut_lag[i] <- 0
        }else{
          caseLS <- caseLS-3
          r_lag <- fit1[[caseLS]]$maximum
          cut_lag[i] <- round( quantile(na_margin_y1,r_lag,na.rm=TRUE) )-1
        }
        
        # Ensemble method
        nafillval <- fittedVal(mymethod[caseLS])
        nafillval <- fittedVal_cut(r_lag, nafillval)
        efit0bestscore <- eLBscore(y0[myNArow2][!ignore_valobs],nafillval[!ignore_valobs],na_var)
        fit0bestscore <- mean(efit0bestscore)
        
        ensemble <- vector(mode = "list",length = 2)
        er_range <- rbind(c(0.04,0.8),c(0.9,0.95))
        for(k in 1:2){
          ensemble[[k]] <- optimize(f = function(r){
            (C <- quantile(abs(nafillval),r,na.rm=TRUE))
            nafillval2 <- nafillval
            nafillval2[which(abs(nafillval)<C)] <- ypred_auto[which(abs(nafillval)<C)]
            LBscore(y0[myNArow2][!ignore_valobs],nafillval2[!ignore_valobs],na_var)
          },interval = er_range[k,],maximum = TRUE)
        }
        
        if(all(is.na(unlist(sapply(ensemble,function(x){x[2]}))))){
          next
        }
        ensemble <- ensemble[[which.max(sapply(ensemble,function(x){x[2]}))]]
        
        caseE <- FALSE
        Er_i <- 0
        Er_C <- 0
        if(fit0bestscore < ensemble$objective){
          caseE <- TRUE
          Er_i <- ensemble$maximum
          Er_C <- quantile(abs(nafillval),Er_i,na.rm = TRUE)
          fit0bestscore <- ensemble$objective
          
          nafillval_tmp <- nafillval
          nafillval_tmp[which(abs(nafillval)<Er_C)] <- ypred_auto[which(abs(nafillval)<Er_C)]
          efit0bestscore <- eLBscore(y0[myNArow2][!ignore_valobs],nafillval_tmp[!ignore_valobs],na_var)
        }
        
        
        # Test if the improvment is significant or not 
        # use paired t test
        sd_score <- sd(efit0bestscore-eautobestScore[!ignore_valobs])/sqrt(length(eautobestScore[!ignore_valobs]))
        autobestScore_fv <- mean(eautobestScore[!ignore_valobs])
        if(is.na(sd_score)){next}
        
        
        if(fill_var[i]>68){
          eps <- eps_tuneset[,2]
        }else{
          eps <- eps_tuneset[,1]
        }
        if(2*sd_score+eps[1]<0){
          epsUsed <- -2*sd_score
        }else{
          epsUsed <- eps[1]
        }
        
        if(fit0bestscore < (autobestScore_fv+2*sd_score+epsUsed)){
          next
        }
        
        # Update the model and predict the test set
        nafill_test <- fittedTest(r_lag, mymethod[caseLS], caseE, Er_C)
        
        if(mean(is.na(nafill_test))==1){
          next
        }
        
        for(tunei in 1:tuneI){
          if(2*sd_score+eps[tunei]<0){
            epsUsed <- -2*sd_score
          }else{
            epsUsed <- eps[tunei]
          }
          if( fit0bestscore >= (autobestScore_fv+2*sd_score+epsUsed) ){
            tune_record_indicator[tunei,i] <- TRUE
          }else{
            tune_record_indicator[tunei,i] <- FALSE
          }
        }
     
        # Record the score if prediction works
        valScore_i <- (fit0bestscore*sum(!ignore_valobs)+sum(eautobestScore[ignore_valobs]))/length(ignore_valobs)
        
        #Save the imputed values and statistics
        caseE_ls[i] <- caseE
        caseLS_ls[i] <- caseLS
        Er_ls[i] <- Er_i
        valScore_ls[i] <- valScore_i
        nafill_test_ls[,i] <- nafill_test
        improveScore_ls[i] <- fit0bestscore-autobestScore_fv
        aimproveScore_ls[i] <- (fit0bestscore-autobestScore_fv)/pmax(autobestScore_fv,10^-15)
      }
      
      invisible(gc())
      
      if(any(tune_record_indicator[1,])){
        # tuneImax <- 0
        # for(tunei in 1:tuneI){
        #   if(any(tune_record_indicator[tunei,])){
        #     tuneImax <- tuneImax+1
        #   }else{
        #     break
        #   }
        # }
        
        tune_pred <- vector(mode="list",length = tuneI)
        tune_num_fill <- numeric(tuneI)
        tuneImp <- numeric(tuneI)
        tuneAImpp <- numeric(tuneI)
        tune_valScore_ls <- matrix(rep(valScore_ls,tuneI),ncol = tuneI)
        for(tunei in 1:tuneI){
          num_fill <- 0
          num_fill_i <- NULL
          n_na <- length(testSegNArow)
          nafill_test_final <- rep(NA,n_na)
          filled_pos <- rep(FALSE,n_na)
          for(i in order(improveScore_ls,decreasing = TRUE)){
            if(!tune_record_indicator[tunei,i]){
              tune_valScore_ls[i,tunei] <- 0
            }else{
              tofill_pos <- !is.na(nafill_test_ls[,i])
              tofill_pos <- tofill_pos & (!filled_pos)
              
              if(sum(tofill_pos)==0){
                tune_valScore_ls[i,tunei] <- 0
                next
              }
              num_fill_i <- c(num_fill_i,sum(tofill_pos))
              num_fill <- num_fill+sum(tofill_pos)
              nafill_test_final[tofill_pos] <- nafill_test_ls[,i][tofill_pos]
              filled_pos <- filled_pos | tofill_pos
            }
          }
          
          #save the result
          tuneImp[tunei] <- sum(improveScore_ls[tune_valScore_ls[,tunei]!=0]*num_fill_i)
          tuneAImpp[tunei]  <- sum(aimproveScore_ls[tune_valScore_ls[,tunei]!=0]*num_fill_i)/sum(num_fill_i)
          tune_num_fill[tunei] <- num_fill
          tune_pred[[tunei]] <- data.frame(tr_wtid[testSegNArow[filled_pos],1:2],var=nafill_test_final[filled_pos])
        }

        
        seg_result[[seg]] <- list(improved=TRUE,
                                  filled_var=cname[tune_valScore_ls[,1]!=0],
                                  cut_lag=cut_lag[tune_valScore_ls[,1]!=0],
                                  LS=table(factor(caseLS_ls[tune_valScore_ls[,1]!=0],levels = 1:3)),
                                  Er=Er_ls[tune_valScore_ls[,1]!=0],
                                  tune_num_fill=tune_num_fill,
                                  tuneImp=tuneImp,
                                  tuneAImpp=tuneAImpp,
                                  autoScore=autoScore,
                                  autobestIndex=autobestIndex,
                                  filled_Score=valScore_ls[tune_valScore_ls[,1]!=0],
                                  tune_pred=tune_pred,
                                  predtestAuto=predtestAuto,
                                  tuneI=tuneI)
        
      }else{
        seg_result[[seg]] <- list(improved=FALSE,
                                  predtestAuto=predtestAuto)
      }
    }
    
    # Finish one variable
    if(!dir.exists(paste0(outputDir))){
      dir.create(paste0(outputDir))
    }  
    
    # if(!dir.exists(paste0(outputDir,wtid))){
    #   dir.create(paste0(outputDir,wtid))
    # }
    
    if(!dir.exists(paste0(outputDir2))){
      dir.create(paste0(outputDir2))
    }  
    
    if(!dir.exists(paste0(outputDir2,wtid))){
      dir.create(paste0(outputDir2,wtid))
    }
    
    tr_na_auto <- NULL
    tune_tr_na <- vector(mode = "list",tuneI) 
    tune_num_fill_wtid_nvar <- vector(mode = "list",tuneI) 
    tune_stat_improveScore_wtid_nvar <- vector(mode = "list",tuneI) 
    tune_stat_aimpp_wtid_nvar <- vector(mode = "list",tuneI) 

    for(seg in 1:n_region){
      tr_na_auto <- rbind(tr_na_auto,seg_result[[seg]]$predtestAuto)
      if(seg_result[[seg]]$improved){
        cat("Segment: ",na_var,"-",seg,"/",n_region,sep="")
        cat("\nfill_var:",seg_result[[seg]]$filled_var,"\n",sep=" ")
        cat("lagcut:",seg_result[[seg]]$cut_lag,"\n",sep=" ")
        cat("LS:",seg_result[[seg]]$LS,"\n",sep=" ")
        cat("Eratio:",seg_result[[seg]]$Er,"\n",sep=" ")
        cat("filled number:",seg_result[[seg]]$tune_num_fill,"\n",sep=" ")
        cat("overall improved:",signif(seg_result[[seg]]$tuneImp,3),"\n",sep=" ")
        cat("Average improved:",signif(seg_result[[seg]]$tuneImp/seg_result[[seg]]$tune_num_fill,3),"\n",sep=" ")
        cat("Average improved percentage:",signif(seg_result[[seg]]$tuneAImpp,3),"\n",sep=" ")
        cat("all score:",signif(seg_result[[seg]]$autoScore,3),"\n",sep=" ")
        cat("best:",c("nearest","mode","linear","WMA1","WMO1","WMA2","WMO2","WMA3","WMO3")[seg_result[[seg]]$autobestIndex],
            seg_result[[seg]]$autoScore[seg_result[[seg]]$autobestIndex],"\n",sep=" ")
        cat("fillScore:",seg_result[[seg]]$filled_Score,"\n\n",sep=" ")   
        
        for(tunei in 1:seg_result[[seg]]$tuneI){
          tune_var_improve[[tunei]] <- unique(c(tune_var_improve[[tunei]],na_var))
          tune_tr_na[[tunei]] <- rbind(tune_tr_na[[tunei]],seg_result[[seg]]$tune_pred[[tunei]])
          
          tune_stat_improveScore_wtid_nvar[[tunei]] <- c(tune_stat_improveScore_wtid_nvar[[tunei]],
                                                         seg_result[[seg]]$tuneImp[tunei])   
          tune_stat_aimpp_wtid_nvar[[tunei]] <- c( tune_stat_aimpp_wtid_nvar[[tunei]],
                                                   seg_result[[seg]]$tuneAImpp[[tunei]] )
          tune_num_fill_wtid_nvar[[tunei]] <- c(tune_num_fill_wtid_nvar[[tunei]],
                                                seg_result[[seg]]$tune_num_fill[[tunei]])
        }
 
      }else{
        # cat(" pass without improvement, \n")
        # cat("nearest,mode,linear,WMA1,WMA2,WMO1,WMO2","\n",seg_result[[seg]]$autoScore,"\n\n",sep=" ")
      }
    }
    
    for(tunei in 1:tuneI){
      if(!is.null(tune_tr_na[[tunei]])){
        if(!dir.exists(paste0(outputDir,"tuneSet",tunei))){
          dir.create(paste0(outputDir,"tuneSet",tunei))
        }
        if(!dir.exists(paste0(outputDir,"tuneSet",tunei,"/",wtid))){
          dir.create(paste0(outputDir,"tuneSet",tunei,"/",wtid))
        }
        fillvar_name <- paste0(outputDir,"tuneSet",tunei,"/",wtid,"/var_",na_var,".rds")
        saveRDS(tune_tr_na[[tunei]],file = fillvar_name)
      }
    }

    
    if(!is.null(tr_na_auto)){
      fillvar_name <- paste0(outputDir2,wtid,"/var_",na_var,".rds")
      saveRDS(tr_na_auto,file = fillvar_name)
    }
    
    
    if(is.null(tune_stat_improveScore_wtid_nvar[[1]])){
      #stat_improveScore_wtid_nvar <- 0
      cat("pass without improvement\n")
      next
    }
      
    
    wtid_var_improvment_ls[[wtid]][[na_var]] <- sapply(1:tuneI,FUN = function(i){
      sum(tune_stat_improveScore_wtid_nvar[[i]])
    })
    
    wtid_var_aimp_ls[[wtid]][[na_var]] <-  sapply(1:tuneI,FUN = function(i){
      sum(tune_stat_improveScore_wtid_nvar[[i]])/sum(tune_num_fill_wtid_nvar[[i]])
    })
    wtid_var_aimp_ls[[wtid]][[na_var]][is.nan(wtid_var_aimp_ls[[wtid]][[na_var]])] <- 0
    
    wtid_var_aimpp_ls[[wtid]][[na_var]] <- sapply(1:tuneI,FUN = function(i){
      sum(tune_stat_aimpp_wtid_nvar[[i]]*tune_num_fill_wtid_nvar[[i]])/sum(tune_num_fill_wtid_nvar[[i]])
    })
    wtid_var_aimpp_ls[[wtid]][[na_var]][is.nan(wtid_var_aimpp_ls[[wtid]][[na_var]])] <- 0
    
    cat("Overall improvement for var ",na_var," of wtid ",wtid,": ",sep="")
    cat(signif(wtid_var_improvment_ls[[wtid]][[na_var]],4),"\n")
    cat("Average improvement for var ",na_var," of wtid ",wtid,": ",sep="")
    cat(signif(wtid_var_aimp_ls[[wtid]][[na_var]],4),"\n")
    cat("Average improvement percentage for var ",na_var," of wtid ",wtid,": ",sep="")
    cat(signif(wtid_var_aimpp_ls[[wtid]][[na_var]]),"\n\n")
    
    cat("Improvement for var ",na_var," of wtid ",wtid,":\n",sep="")
    print(  cbind(eps_tuneset,
                  imp=signif(wtid_var_improvment_ls[[wtid]][[na_var]],6),
                  aimp=signif(wtid_var_aimp_ls[[wtid]][[na_var]],4),
                  aimpp=signif(wtid_var_aimpp_ls[[wtid]][[na_var]],4)
                  ))
    cat("\n")
    
    
    stat_improveScore_wtid <- stat_improveScore_wtid + wtid_var_improvment_ls[[wtid]][[na_var]]
    stat_aimpp_wtid <- stat_aimpp_wtid + sapply(1:tuneI,FUN = function(i){
      sum(tune_stat_aimpp_wtid_nvar[[i]]*tune_num_fill_wtid_nvar[[i]])
    })
    num_fill_wtid <- num_fill_wtid + sapply(1:tuneI,FUN = function(i){
      sum(tune_num_fill_wtid_nvar[[i]])
    })
    
    rm(tune_tr_na,tr_na_auto)
    invisible(gc())
  }
  
  tc_var <- proc.time()-tc_var
  cat("\n","Finish",wtid," in ",tc_var[3],"s",sep="","\n")
  cat("improvement for wtid ",wtid,":\n",sep="")
  
  if(length(clusterset)>1){
    print(  cbind(eps_tuneset,
                  imp=signif(stat_improveScore_wtid,6),
                  aimp=signif(stat_improveScore_wtid/num_fill_wtid,4),
                  aimpp=signif(stat_aimpp_wtid/num_fill_wtid,4),
                  impvar=sapply(tune_var_improve,FUN = function(x){paste(x,collapse = ",")})))
  }

  cat("\n")

  
  saveRDS(list(imp=wtid_var_improvment_ls,
               aimp=wtid_var_aimp_ls,
               aimpp=wtid_var_aimpp_ls),paste0(outputDir,"wtid_var_imp_ls_",edition,who,trial,".rds") )
  
}




