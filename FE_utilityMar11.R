#library(zeallot)



# Simulation method for auto 
simNArow <- function(isna_pos,
                     n,
                     shape=1.0645002499,
                     rate=0.0025810162,
                     N=23){
  NA_rowIdx <- (1:n)[isna_pos]
  nonNA_rowIdx <- (1:n)[!isna_pos]
  missing_center <- sample(nonNA_rowIdx,size = N,replace = FALSE)
  missing_center <- sort(missing_center)
  missing_halflength <- round(rgamma(N,shape,rate)/2)
  
  missing_location <- cbind(missing_center-missing_halflength,
                            missing_center+missing_halflength)
  myNA_rowIdx <- NULL
  for(i in 1:N){
    myNA_rowIdx <- c(myNA_rowIdx,seq(missing_location[i,1],missing_location[i,2]))
  }
  myNA_rowIdx <- setdiff(unique(myNA_rowIdx),NA_rowIdx)
  myNA_rowIdx <- myNA_rowIdx[myNA_rowIdx>1 & myNA_rowIdx<n]
  return(myNA_rowIdx)
}

# (Local) Simulation method for auto
simNArow2 <- function(isna_pos,n,L=30,E=10){
  NA_rowIdx <- (1:n)[isna_pos]
  nonNA_rowIdx <- (1:n)[!isna_pos]
  
  NA_rowIdx_endpoint <- NA_rowIdx[conEndpoint(NA_rowIdx)]
  dn_region <- length(NA_rowIdx_endpoint)
  
  region_length <- NA_rowIdx_endpoint[seq(2,dn_region,2)]- NA_rowIdx_endpoint[seq(1,dn_region,2)]
  
  
  missing_location1 <- cbind(NA_rowIdx_endpoint[seq(1,dn_region,2)]-region_length-L-E,
                             NA_rowIdx_endpoint[seq(2,dn_region,2)]-region_length-L)
  
  missing_location2 <- cbind(NA_rowIdx_endpoint[seq(1,dn_region,2)]+region_length+L,
                             NA_rowIdx_endpoint[seq(2,dn_region,2)]+region_length+L+E)
  
  region_simNA <- vector(mode = "list",length = dn_region/2)
  myNA_rowIdx <- NULL
  myNA_rowIdx1 <- NULL
  myNA_rowIdx2 <- NULL
  testNArow <- NULL
  for(i in 1:(dn_region/2)){
    new_rowIdx1 <- c(missing_location1[i,1]:missing_location1[i,2])
    new_rowIdx2 <- c(missing_location2[i,1]:missing_location2[i,2])
    new_rowIdx1 <- setdiff(new_rowIdx1,NA_rowIdx)
    new_rowIdx1 <- new_rowIdx1[new_rowIdx1>1 & new_rowIdx1<n]
    new_rowIdx2 <- setdiff(new_rowIdx2,NA_rowIdx)
    new_rowIdx2 <- new_rowIdx2[new_rowIdx2>1 & new_rowIdx2<n]
    
    new_rowIdx <- c(new_rowIdx1,new_rowIdx2)
    testSegNArow <- NA_rowIdx_endpoint[1+2*(i-1)]:NA_rowIdx_endpoint[2+2*(i-1)]
    region_simNA[[i]] <- list(testNArow=testSegNArow,
                              myNA_rowIdx=new_rowIdx,
                              myNA_rowIdx1=new_rowIdx1,
                              myNA_rowIdx2=new_rowIdx2)
    myNA_rowIdx <- c(myNA_rowIdx,new_rowIdx)
    myNA_rowIdx1 <- c(myNA_rowIdx1,new_rowIdx1)
    myNA_rowIdx2 <- c(myNA_rowIdx2,new_rowIdx2)
    testNArow <- c(testNArow,testSegNArow)
  }
  
  myNA_rowIdx <- setdiff(sort(unique(myNA_rowIdx)),NA_rowIdx)
  myNA_rowIdx <- myNA_rowIdx[myNA_rowIdx>1 & myNA_rowIdx<n]
  return(myNA_rowIdx)
}

# (Local) Simulation method for auto Version 2
simNArow3 <- function(isna_pos,n,L=30,E=10){
  NA_rowIdx <- (1:n)[isna_pos]
  nonNA_rowIdx <- (1:n)[!isna_pos]
  
  NA_rowIdx_endpoint <- NA_rowIdx[conEndpoint(NA_rowIdx)]
  dn_region <- length(NA_rowIdx_endpoint)
  
  region_length <- round((NA_rowIdx_endpoint[seq(2,dn_region,2)]- NA_rowIdx_endpoint[seq(1,dn_region,2)])/2) #change this
  
  
  missing_location1 <- cbind(NA_rowIdx_endpoint[seq(1,dn_region,2)]-region_length-L-E,
                             NA_rowIdx_endpoint[seq(2,dn_region,2)]-region_length-L)
  
  missing_location2 <- cbind(NA_rowIdx_endpoint[seq(1,dn_region,2)]+region_length+L,
                             NA_rowIdx_endpoint[seq(2,dn_region,2)]+region_length+L+E)
  
  region_simNA <- vector(mode = "list",length = dn_region/2)
  myNA_rowIdx <- NULL
  myNA_rowIdx1 <- NULL
  myNA_rowIdx2 <- NULL
  testNArow <- NULL
  for(i in 1:(dn_region/2)){
    new_rowIdx1 <- c(missing_location1[i,1]:missing_location1[i,2])
    new_rowIdx2 <- c(missing_location2[i,1]:missing_location2[i,2])
    new_rowIdx1 <- setdiff(new_rowIdx1,NA_rowIdx)
    new_rowIdx1 <- new_rowIdx1[new_rowIdx1>1 & new_rowIdx1<n]
    new_rowIdx2 <- setdiff(new_rowIdx2,NA_rowIdx)
    new_rowIdx2 <- new_rowIdx2[new_rowIdx2>1 & new_rowIdx2<n]
    
    new_rowIdx <- c(new_rowIdx1,new_rowIdx2)
    testSegNArow <- NA_rowIdx_endpoint[1+2*(i-1)]:NA_rowIdx_endpoint[2+2*(i-1)]
    region_simNA[[i]] <- list(testNArow=testSegNArow,
                              myNA_rowIdx=new_rowIdx,
                              myNA_rowIdx1=new_rowIdx1,
                              myNA_rowIdx2=new_rowIdx2)
    myNA_rowIdx <- c(myNA_rowIdx,new_rowIdx)
    myNA_rowIdx1 <- c(myNA_rowIdx1,new_rowIdx1)
    myNA_rowIdx2 <- c(myNA_rowIdx2,new_rowIdx2)
    testNArow <- c(testNArow,testSegNArow)
  }
  
  myNA_rowIdx <- setdiff(sort(unique(myNA_rowIdx)),NA_rowIdx)
  myNA_rowIdx <- myNA_rowIdx[myNA_rowIdx>1 & myNA_rowIdx<n]
  return(myNA_rowIdx)
}

# (Local) Simulation method for cross
simNArow5 <- function(isna_pos,n,L=30,E=10){
  NA_rowIdx <- (1:n)[isna_pos]
  nonNA_rowIdx <- (1:n)[!isna_pos]
  
  NA_rowIdx_endpoint <- NA_rowIdx[conEndpoint(NA_rowIdx)]
  dn_region <- length(NA_rowIdx_endpoint)
  
  region_length <- NA_rowIdx_endpoint[seq(2,dn_region,2)]- NA_rowIdx_endpoint[seq(1,dn_region,2)]
  
  
  missing_location1 <- cbind(NA_rowIdx_endpoint[seq(1,dn_region,2)]-region_length-L-E,
                             NA_rowIdx_endpoint[seq(2,dn_region,2)]-region_length-L)
  
  missing_location2 <- cbind(NA_rowIdx_endpoint[seq(1,dn_region,2)]+region_length+L,
                             NA_rowIdx_endpoint[seq(2,dn_region,2)]+region_length+L+E)
  
  region_simNA <- vector(mode = "list",length = dn_region/2)
  #myNA_rowIdx <- NULL
  myNA_rowIdx1 <- NULL
  myNA_rowIdx2 <- NULL
  testNArow <- NULL
  for(i in 1:(dn_region/2)){
    new_rowIdx1 <- c(missing_location1[i,1]:missing_location1[i,2])
    new_rowIdx2 <- c(missing_location2[i,1]:missing_location2[i,2])
    new_rowIdx1 <- setdiff(new_rowIdx1,NA_rowIdx)
    new_rowIdx1 <- new_rowIdx1[new_rowIdx1>1 & new_rowIdx1<n]
    new_rowIdx2 <- setdiff(new_rowIdx2,NA_rowIdx)
    new_rowIdx2 <- new_rowIdx2[new_rowIdx2>1 & new_rowIdx2<n]
    
    new_rowIdx <- c(new_rowIdx1,new_rowIdx2)
    testSegNArow <- NA_rowIdx_endpoint[1+2*(i-1)]:NA_rowIdx_endpoint[2+2*(i-1)]
    region_simNA[[i]] <- list(testNArow=testSegNArow,
                              myNA_rowIdx=new_rowIdx,
                              myNA_rowIdx1=new_rowIdx1,
                              myNA_rowIdx2=new_rowIdx2)
    #myNA_rowIdx <- c(myNA_rowIdx,new_rowIdx)
    myNA_rowIdx1 <- c(myNA_rowIdx1,new_rowIdx1)
    myNA_rowIdx2 <- c(myNA_rowIdx2,new_rowIdx2)
    testNArow <- c(testNArow,testSegNArow)
  }
  
  #myNA_rowIdx <- sort(unique(myNA_rowIdx))
  
  return(list(region_simNA=region_simNA,
              #myNA_rowIdx=myNA_rowIdx,
              myNA_rowIdx1=myNA_rowIdx1,
              myNA_rowIdx2=myNA_rowIdx2,
              testNArow=testNArow))
}







#---------------------------------------- Evaluation Method---------------------------------------- #
rmse <- function(y,yhat){
  sqrt(mean((y-yhat)^2))
}

LBscore <- function(y,ypred,var_j){
  if(var_j %in% c(53,66)){
    score <- mean(y==ypred)
  }else{
    score <- mean(exp(-100*abs(y-ypred)/pmax(10^-15,abs(y))))
  }
  score
}

eLBscore <- function(y,ypred,var_j){
  if(var_j %in% c(53,66)){
    score <- (y==ypred)
  }else{
    score <- (exp(-100*abs(y-ypred)/pmax(10^-15,abs(y))))
  }
  score
}

Mode <- function(x) {
  ux <- setdiff(unique(x),NA)
  ux[which.max(tabulate(match(x, ux)))]
}

findMargin <- function(x){
  napos <- which(is.na(x))
  diff_napos <- diff(napos)
  left_margin <- locateNNA(diff_napos)
  right_margin <- rev(locateNNA(rev(diff_napos)))
  lr_margin <- pmin(left_margin,right_margin)
  lr_margin
}

nacluster <- function(x){
  napos <- which(is.na(x))
  diffnapos <- diff(napos)
  rle(diffnapos)
}

#R version, depreciated, replaced by Rcpp
# wMode <- function(x,w){
#   weigtedCount <- aggregate(w,by=list(valueGroup=x),FUN=sum)
#   my_wmode <- weigtedCount$valueGroup[which.max(weigtedCount$x)]
#   return(my_wmode)
# }


#------------------- Method of Moving Average/Mode ---------------------------------------#
# y <- c(1,2,3,4,NA,NA,NA,NA,3,2)
# napos <- which(is.na(y))
# nearestFE(y,napos,10,sidediff = "Both")
# tr_WM <- trainWMA(y,napos,10)
# predictWM(tr_WM,5)

nearestFE <- function(y,napos,kmax,sidediff=c("No","Yes","Both")){
  sidediff <- match.arg(sidediff)
  n <- length(y)
  diff_napos <- diff(napos)
  left_margin <- locateNNA(diff_napos)
  right_margin <- rev(locateNNA(rev(diff_napos)))
  n_na <- length(napos)
  
  Margin <- Mat <- SLD <- NULL
  if(sidediff=="No" | sidediff=="Both" ){
    # Find the margin of the nearest 2*kmax neighours within two sides
    nonNA_margin <- locateKNNA(left_margin,right_margin,2*kmax)
    nearestKNNA_Margin <- nonNA_margin[,1:(2*kmax),drop=FALSE]
    nearestKNNA_location <- nearestKNNA_Margin
    nearestKNNA_location <- napos+nearestKNNA_location
    nearestKNNA_location[nearestKNNA_location<=0] <- 1
    nearestKNNA_location[nearestKNNA_location>=n] <- n
    nearestKNNA_Mat <- matrix(y[nearestKNNA_location],ncol=(2*kmax))

    Margin <- abs(nearestKNNA_Margin)
    Mat <- nearestKNNA_Mat
    SLD <- nonNA_margin[,((2*kmax)+1):((2*kmax)+3)]
  }


  yna_leftright_margin_Mat <- yna_leftright_Mat <- NULL
  if(sidediff=="Yes" | sidediff=="Both" ){
    # Find the margin of the nearest kmax neighours for each side
    yna_left_margin_Mat <- matrix(rep(left_margin,kmax),ncol=kmax)+
      matrix(1:kmax-1,nrow=n_na,ncol=kmax,byrow = TRUE)
    yna_right_margin_Mat <- matrix(rep(right_margin,kmax),ncol=kmax)+
      matrix(1:kmax-1,nrow=n_na,ncol=kmax,byrow = TRUE)
    yna_leftright_margin_Mat <- cbind(yna_left_margin_Mat,yna_right_margin_Mat)
    yna_leftright_pos_Mat <- cbind(napos-yna_left_margin_Mat,napos+yna_right_margin_Mat)
    yna_leftright_pos_Mat[yna_leftright_pos_Mat<1] <- 1
    yna_leftright_pos_Mat[yna_leftright_pos_Mat>n] <- n

    yna_leftright_Mat <- matrix(y[yna_leftright_pos_Mat],ncol=2*kmax)
    isna_yna_leftright_Mat <- is.na(yna_leftright_Mat)
    yna_leftright_Mat[isna_yna_leftright_Mat] <- 0
  }

  return(list(Margin=Margin,
              Mat=Mat,
              SLD=SLD,
              LRMargin=yna_leftright_margin_Mat,
              LRMat=yna_leftright_Mat))
}

# train
trainWMA <- function(y,napos,kmax){
  #Currently y_SLD is not used
  c(y_Margin,y_Mat,y_SLD,LRMargin,LRMat) %<-% nearestFE(y,napos,kmax,"Both")

  # Prepare, Clean possible NA
  isna_y_Mat <- is.na(y_Mat)
  y_Mat[isna_y_Mat] <- 0
  UNset <- lapply(1:nrow(y_Mat),function(i){unique(y_Mat[i,])})

  isna_LRMat <- is.na(LRMat)
  LRMat[isna_LRMat] <- 0
  LRUNset <- lapply(1:nrow(LRMat),function(i){unique(LRMat[i,])})

  return(list(y_Margin=y_Margin,
              y_Mat=y_Mat,
              y_SLD=y_SLD,
              isna_y_Mat=isna_y_Mat,
              UNset=UNset,
              LRMargin=LRMargin,
              LRMat=LRMat,
              isna_LRMat=isna_LRMat,
              LRUNset=LRUNset,
              kmax=kmax))
}

# For WMA WMO
trainWMA1 <- function(y,napos,kmax){
  #Currently y_SLD is not used
  c(y_Margin,y_Mat,y_SLD,LRMargin,LRMat) %<-% nearestFE(y,napos,kmax,"No")
  
  # Prepare, Clean possible NA
  isna_y_Mat <- is.na(y_Mat)
  y_Mat[isna_y_Mat] <- 0
  UNset <- lapply(1:nrow(y_Mat),function(i){unique(y_Mat[i,])})
  
  return(list(y_Margin=y_Margin,
              y_Mat=y_Mat,
              y_SLD=y_SLD,
              isna_y_Mat=isna_y_Mat,
              UNset=UNset,
              kmax=kmax))
}


# For LWA LWO
trainWMA2 <- function(y,napos,kmax){
  #Currently y_SLD is not used
  c(y_Margin,y_Mat,y_SLD,LRMargin,LRMat) %<-% nearestFE(y,napos,kmax,"Yes")
  
  # Prepare, Clean possible NA
  isna_LRMat <- is.na(LRMat)
  LRMat[isna_LRMat] <- 0
  LRUNset <- lapply(1:nrow(LRMat),function(i){unique(LRMat[i,])})
  
  return(list(LRMargin=LRMargin,
              LRMat=LRMat,
              isna_LRMat=isna_LRMat,
              LRUNset=LRUNset,
              kmax=kmax))
}

# predict
predictWM <- function(tr_WMA,k,WM=c("WMA","WMO","LWA","LWO")){
  kmax <- tr_WMA$kmax
  if(k>kmax){stop("exceed kmax!")}
  alpha <- 1
  if(!is.null(tr_WMA$y_Mat)){
    nafill_mat <- matrix(0,nrow(tr_WMA$y_Mat),length(WM))
  }else{
    nafill_mat <- matrix(0,nrow(tr_WMA$LRMat),length(WM))
  }
  j  <- 1
  
  if(any(c("WMA","WMO") %in% WM)){
    linear_weight <-   1/(tr_WMA$y_Margin[,1:(2*k),drop=FALSE]+1)^alpha
    linear_weight[tr_WMA$isna_y_Mat[,1:(2*k),drop=FALSE]] <- 0
    linear_weight <- linear_weight/rowSums(linear_weight)

    if("WMA" %in% WM){
      nafill_mat[,j] <- rowSums(tr_WMA$y_Mat[,1:(2*k),drop=FALSE]*linear_weight)
      j <- j+1
    }

    if("WMO" %in% WM){
      nafill_mat[,j] <- wModeFillMat(tr_WMA$UNset,tr_WMA$y_Mat[,1:(2*k),drop=FALSE],linear_weight) #Rcpp version
      j <- j+1
    }
  }

  if(any(c("LWA","LWO") %in% WM) ){
    linear_weight <-   1/(tr_WMA$LRMargin[,c(1:k,kmax+1:k),drop=FALSE]+1)^alpha
    linear_weight[tr_WMA$isna_LRMat[,c(1:k,kmax+1:k),drop=FALSE]] <- 0
    linear_weight <- linear_weight/rowSums(linear_weight)

    if("LWA" %in% WM){
      nafill_mat[,j] <- rowSums(tr_WMA$LRMat[,c(1:k,kmax+1:k),drop=FALSE]*linear_weight)
      j <- j+1
    }

    if("LWO" %in% WM){
      nafill_mat[,j] <- wModeFillMat(tr_WMA$LRUNset,tr_WMA$LRMat[,c(1:k,kmax+1:k),drop=FALSE],linear_weight) #Rcpp version
      j <- j+1
    }
  }

  return(nafill_mat)
}
