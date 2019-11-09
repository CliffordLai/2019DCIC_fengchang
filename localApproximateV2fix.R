# My own LM with MAPE
mylm <- function(y,x){
  nonNaXY <- (!is.na(y)) & (!is.na(x))
  y <- y[nonNaXY]
  x <- x[nonNaXY]
  negscore <- function(b){
    -LBscore(y,b[1]+b[2]*x,na_var)
  }
  ans1 <- optim(par = c(0,1),fn = negscore)
  ans1
}


# val LS
fittedVal <- function(method=c("myls","ols","gam")){
  if(method=="myls"){
    fitval <- mylm(y0[myNArow1],y_fill_var[myNArow1])
    nafillval <- fitval$par[1]+fitval$par[2]*y_fill_var[myNArow2]
    
  }else if(method=="ols"){
    fitval <- lm(y0[myNArow1]~y_fill_var[myNArow1])
    if(!any(is.na(fitval$coefficients))){
      nafillval <- fitval$coefficients[1]+fitval$coefficients[2]*y_fill_var[myNArow2]
    }else{
      nafillval <- rep(0,length(myNArow2))
    }
    
  }else if(method=="gam"){
    fitval <-  try(gam(y1~s(y_fill_var), data=data.frame(y1=y0[myNArow1],y_fill_var=y_fill_var[myNArow1])),silent = TRUE)
    if(!class(fitval)[1]=="try-error"){
      nafillval <- predict(fitval,newdata = data.frame(y_fill_var=y_fill_var[myNArow2]))
    }else{
      nafillval <- rep(0,length(myNArow2))
    }
  }
  
  return(nafillval)
}

fittedVal_cut <- function(r_lag, nafillval){
  fill_lagcut <- round( quantile(na_margin_y1,r_lag,na.rm = TRUE) )
  lagmissVal_val <- na_margin_y2 < fill_lagcut
  nafillval[lagmissVal_val] <- ypred_auto[lagmissVal_val]
  return(nafillval)
}


# test LS
fittedTest <- function(r_lag, method=c("myls","ols","gam"), caseE, Er_C){
  fill_lagcut <- round( quantile(na_margin_y1,r_lag,na.rm = TRUE) )
  lagmissVal_val <- na_margin_y0seg < fill_lagcut
  
  if(method=="myls"){
    fitval <- mylm(y0[myNArow],y_fill_var[myNArow])
    nafilltest <- fitval$par[1]+fitval$par[2]*y_fill_var[testSegNArow]
    
  }else if(method=="ols"){
    fitval <- lm(y0[myNArow]~y_fill_var[myNArow])
    if(!any(is.na(fitval$coefficients))){
      nafilltest <- fitval$coefficients[1]+fitval$coefficients[2]*y_fill_var[testSegNArow]
    }else{
      nafilltest <- rep(0,length(testSegNArow))
    }
    
  }else if(method=="gam"){
    fitval <-  try(gam(y1~s(y_fill_var), data=data.frame(y1=y0[myNArow],y_fill_var=y_fill_var[myNArow])),silent = TRUE)
    if(!class(fitval)[1]=="try-error"){
      nafilltest <- predict(fitval,newdata = data.frame(y_fill_var=y_fill_var[testSegNArow]))
    }else{
      nafilltest <- rep(0,length(testSegNArow))
    }
  }
  
  # cut
  nafilltest[lagmissVal_val] <- NA
  
  # ensemble
  if(caseE){
    nafilltest[which(abs(nafilltest)<Er_C)] <- NA 
  }
 
  return(nafilltest)
}