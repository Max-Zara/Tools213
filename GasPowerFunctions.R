library(lubridate)
library(data.table)

#Load Spot Pamars
NBP.Spot.Params <- c("NBP","UK","NBP","WD","D")
API2.Spot.Params <- c("API2","NL","API2","BOM","D")
EUA.Spot.Params <- c("EUA","EU","EUA","WD","D")
TTF.Spot.Params <- c("TTF","NL","TTF","WD","D")
Spot.Params <- rbind(NBP.Spot.Params,API2.Spot.Params, EUA.Spot.Params,TTF.Spot.Params)
colnames(Spot.Params) <- c("Product","Region","ContractID2","ContractID3","Duration")

get.Storage.Limits <- function(Storage.Noms.Merged, to.Chart = FALSE, add.width = 0){
  Storage.Noms.Merged$DOY <- yday(Storage.Noms.Merged$DDate)
  
  ncol(Storage.Noms.Merged); colnames(Storage.Noms.Merged);
  
  field.ranges <- data.frame(DOY = seq(1,366,1))
  for(icol in 3:25){
    f <- as.formula(paste0(colnames(Storage.Noms.Merged)[icol]," ~ DOY"))
    temp.range <- aggregate(f,Storage.Noms.Merged,range)
    #if(icol == 25){field.ranges <- rbind(field.ranges,temp.range)}else{
      field.ranges <- merge(field.ranges,temp.range,by="DOY", all.x = TRUE, all.y = TRUE)
    #}
    print(paste0(icol,"  Rows: ", nrow(temp.range)))
  }
  head(field.ranges)
  field.ranges <- as.matrix(field.ranges)
  field.ranges[field.ranges==0] <- NA
  ##replace nulls with absolute max
  for(icol in seq(2,46,2)){
    series.range <- range(na.omit(field.ranges[,c(icol,icol+1)]))
    if(!is.finite(series.range[1])){series.range <- c(0,0)}
    field.ranges[is.na(field.ranges[,icol]),icol] <- series.range[1] + add.width*mean(series.range[2]-series.range[1]) 
    field.ranges[is.na(field.ranges[,icol+1]),icol+1] <- series.range[2] - add.width*mean(series.range[2]-series.range[1]) 
    print(paste0("Done:",icol, " Range: ",series.range))
  }
  #field.ranges <- na.locf(field.ranges)
 
  if(to.Chart){
    iprint=0
    par(mfrow=c(3,3))
    for(icol in seq(2,46,2)){
      if(sum(is.na(field.ranges[,c(icol,icol+1)]))<300){
        if(iprint==9){readline("Press Enter for Next Charts....")}
        plot(field.ranges[,icol],main=colnames(field.ranges)[icol],ylim=range(na.omit(field.ranges[,c(icol,icol+1)])),type='l',col='blue')
        lines(field.ranges[,icol+1],type='l',col='red')
        iprint = iprint+1;
        }
    }
  reset_par();
  }
  
  return(field.ranges)
}


#start.date <- as.Date("2001-01-01")
#end.date <- as.Date("2010-12-31")
library(lubridate)
construct.min.max.limits <- function(start.date, end.date, field.ranges){
  temp.ts <- data.frame(DDate = seq(start.date,end.date,1), DOY = yday(seq(start.date,end.date,1)))
  temp.merged <- merge(temp.ts,field.ranges,by="DOY")
  temp.merged <- temp.merged[order(temp.merged$DDate),]
  temp.merged <- ts(temp.merged[,-c(1,2)],start = decimal_date(temp.merged[1,]$DDate),frequency=365)
  #plot(temp.merged[,3],type='l')
  return(temp.merged)
}

library(dlm)
library(zoo)
get.Kalman.Smoothed.Max.Min <- function(temp.merged, get.param=FALSE, param.loc ="Model//Limits//StorageLimitsMLE.csv", i.start.col=1, to.Chart = FALSE){
  
  col.num <- ncol(temp.merged)
  
  buildHarm <- function(u){
    poly.mod <- dlmModPoly(1, dV = u[1], dW=u[2])
    harmon.mod <- dlmModTrig(q=2,tau=365,dV = u[3],dW=c(u[4],u[5],u[6],u[7]))
  }
  
  if(get.param){
    param.MLE <- double()
    print(paste0("Getting MLE Params to ",colnames(temp.merged)[icol]," Writing to:",param.loc))
    for(icol in 1:col.num){
      temp.ts <- temp.merged[,icol]
      
      if(max(temp.ts)-min(temp.ts)!=0){
      
      init <- rep(1,7)
      dlmHarmMLE <- dlmMLE(temp.ts,init,buildHarm)
      param.MLE <- rbind(param.MLE,temp=dlmHarmMLE$par)
      rownames(param.MLE)[nrow(param.MLE)]=colnames(temp.merged)[icol]
      print(paste0("Done: ",colnames(temp.merged)[icol]))
      }else{
        print(paste0("ONLY 0s: ",colnames(temp.merged)[icol]))
      }
    }
    colnames(param.MLE ) <- c("Poly.V","Poly.W","Harmon.V","Harmon.W1","Harmon.W2","Harmon.W3","Harmon.W4")
    write.csv(param.MLE,file=param.loc)
  }
  
  ##Starting Creating the Smoothed Max Mins
  print("Starting Creating the Smoothed Max Mins")
  param.MLE <- read.csv(file=param.loc)
  
  if(i.start.col == 1){
    maxmins <- NULL;maxmins$DOY <- seq(1,366,1)}else{
  maxmins <- read.csv("Model//Limits//StorageSmoothedLimitsMLE.csv")}
  
  for(icol in i.start.col:col.num){
    temp.ts <- temp.merged[,icol]
    
    #get the MLE stats
    temp.param <- unlist(param.MLE[param.MLE$X == colnames(temp.merged)[icol],])[-1]
    if(length(temp.param)>0){
      print(paste0(colnames(temp.merged)[icol],": Found in CSV"))
      dlmHarm <- buildHarm(temp.param)
      temp.smooth.maxmin <- dlmSmooth(temp.ts,dlmHarm)
      #Now Extract the Max Min Level
      #get 2nd year only
      maxmin.year <- window(temp.smooth.maxmin$s,2003,2004)
      #plot(temp.smooth.maxmin$s)
      temp.level <- data.frame(DOY = yday(as.Date(as.numeric((time(maxmin.year)-1990)*365.25),origin="1990-01-01")),Temp= rowSums(maxmin.year))
      colnames(temp.level)[2] <-colnames(temp.merged)[icol]
      #plot(temp.level[1:365,2],type='l',col="red",ylim=range(c(temp.level[1:365,2],temp.ts[1:365]))); lines(temp.ts[1:365])
      maxmins <- merge(maxmins,temp.level,by="DOY")
      #restrict to only 1 value
      maxmins <- maxmins[match(seq(1,365,1),maxmins$DOY),]
      write.csv(maxmins,"Model//Limits//StorageSmoothedLimitsMLE.csv")
    }else{
      print(paste0(colnames(temp.merged)[icol],"Missing in CSV - writing 0s"))
      temp.level <- data.frame(DOY = yday(as.Date(as.numeric((time(maxmin.year)-1990)*365.25),origin="1990-01-01")),Temp= 0)
      colnames(temp.level)[2] <-colnames(temp.merged)[icol]
      maxmins <- merge(maxmins,temp.level,by="DOY")
    }
    
  }
  
  #Make sure maxmins are clean
  maxmins2 <- maxmins; maxmins3 <- maxmins;
  for(icol in seq(2,ncol(maxmins),by=2)){
    maxmins2[,icol] <- unlist(apply(maxmins[,c(icol,icol+1)],1,min))
    maxmins2[,icol+1] <- unlist(apply(maxmins[,c(icol,icol+1)],1,max))
  }
  maxmins <- maxmins2
  
  #number of harmonics q = 3, period is tau
  if(to.Chart){
    reset_par(); par(mfrow=c(3,3))
    iprint=0
    for(i in seq(1,col.num,2)){
      iprint=iprint+1
      plot(maxmins[,i+1],ylim=range(maxmins[,c(i+1,i+2)]),type='l',col="red",main=colnames(maxmins)[i+1]); lines(maxmins[,i+2],type='l',col="darkred")
      print(i)
      if(iprint==9){
        iprint=0
        readline("Press Enter for Next Charts...")
      }
    }
    reset_par()
  }
  
  return(maxmins)
}

get.Seasonal.Cat <- function(x, seasonal.cat="M") {
  if(seasonal.cat =="M"){
    return(month(x))
  }else if(seasonal.cat == "W"){
    return(week(x))
  }else if(seasonal.cat == "D"){
    return(yday(x))
  }else{return(null)}
}

get.Seasonal.Temps <- function(mydb,seasonal.cat = "M", HDD.Temp.Break = 18, to.chart=FALSE){
  Cat1 = "Temperature"; Cat2 = "Seasonal Normal"; Cat3 = "D"
  dSeas <- Get.UK.Gas.Data(mydb,Cat1,Cat2,Cat3)
  
  if(seasonal.cat =="M"){
  dSeas$Seasonal.Cat <- month(dSeas$DDate)
  }else if(seasonal.cat == "W"){
    dSeas$Seasonal.Cat <- week(dSeas$DDate)
  }else if(seasonal.cat == "D"){
    dSeas$Seasonal.Cat <- yday(dSeas$DDate)
  }else{return(null)}
  
  Aggregated.Seasonal <- aggregate(Value ~ Seasonal.Cat,dSeas,mean)
  if(to.chart)
    plot(Aggregated.Seasonal)
  
  Aggregated.Seasonal$HDD <- ifelse(HDD.Temp.Break-Aggregated.Seasonal$Value<0,0,HDD.Temp.Break-Aggregated.Seasonal$Value)
  Aggregated.Seasonal$CDD <- ifelse(Aggregated.Seasonal$Value-HDD.Temp.Break<0,0,Aggregated.Seasonal$Value-HDD.Temp.Break)
  return(Aggregated.Seasonal)
}

get.Actual.Temps <- function(mydb, Cat1, Cat3, HDD.Temp.Break = 18){
  #Get Actual Temperatures
  dTempActual <- Get.UK.Gas.Data(mydb, Cat1,"Forecast" ,Cat3)
  dTempActual$HDD <- ifelse(HDD.Temp.Break-dTempActual$Value<0,0,HDD.Temp.Break-dTempActual$Value)
  dTempActual$CDD <- ifelse(dTempActual$Value-HDD.Temp.Break<0,0,dTempActual$Value-HDD.Temp.Break)
  head(dTempActual)
  return(dTempActual)
}

