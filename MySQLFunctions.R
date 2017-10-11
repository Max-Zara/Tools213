##MYSQL Queries for My DB




##Gets Data from the Natrgid_Data of MYSQL
Get.UK.Gas.Data <- function(mydb, Cat1=0,Cat2=0,Cat3=0){
  
  if(Cat1 == 0 & Cat2 == 0 & Cat3 == 0){
    Sql.Query = paste0("SELECT * FROM natgrid_data ;")
  }  else if(Cat2 == 0 & Cat3 == 0){
    Sql.Query = paste0("SELECT * FROM natgrid_data WHERE `Cat1`='",Cat1,"' ;")    
  }else if(Cat3 == 0){
    Sql.Query = paste0("SELECT * FROM natgrid_data WHERE `Cat1`='",Cat1,"' And `Cat2`='",Cat2,"' ;")
  }else{
    Sql.Query = paste0("SELECT * FROM natgrid_data WHERE `Cat1`='",Cat1,"' And `Cat2`='",Cat2,"' and `Cat3` = '",Cat3,"';")
  }
  
  Sql.result <- dbGetQuery(mydb,Sql.Query)
  Sql.result$DDate <- as.Date(Sql.result$DDate)
  Sql.result$PDate <- as.Date(Sql.result$PDate)
  return(Sql.result)
}

Get.NatGrid.Data <- function(mydb, Cat1 = "0", Cat2 = "0", Cat3 = "0", Cat4 = "0",
                             Cat1.In = FALSE, Cat2.In = FALSE, Cat3.In = FALSE, Cat4.In = FALSE,
                             DDate.Limit = as.Date(c(0,0),origin="1990-01-01")){
  
  date.lim <- ""; if(max(DDate.Limit)!="1990-01-01"){date.lim <- paste0(" And DDate >= ",DDate.Limit[1]," And DDate <= ",DDate.Limit[2]) }
  cat1.lim <- "";  if(Cat1 != "0"){cat1.lim <- paste0(" And Cat1 ",if(Cat1.In){paste0(" IN ",Cat1)}else{paste0("='",Cat1,"'")})};
  cat2.lim <- "";  if(Cat2 != "0"){cat2.lim <- paste0(" And Cat2 ",if(Cat2.In){paste0(" IN ",Cat2)}else{paste0("='",Cat2,"'")})};
  cat3.lim <- ""; if(Cat3 != "0"){cat3.lim <- paste0(" And Cat3 ",if(Cat3.In){paste0(" IN ",Cat3)}else{paste0("='",Cat3,"'")})};
  cat4.lim <- ""; if(Cat4 != "0"){cat4.lim <- paste0(" And Cat4 ",if(Cat3.In){paste0(" IN ",Cat4)}else{paste0("='",Cat4,"'")})};
  
  Sql.Query <- paste0("SELECT * FROM products.natgrid_data WHERE DDATE >= 0",date.lim,cat1.lim,cat2.lim,cat3.lim,cat4.lim)
  print(Sql.Query)
  Sql.result <- dbGetQuery(mydb,Sql.Query)
  Sql.result$DDate <- as.Date(Sql.result$DDate)
  #Sql.result$PDate <- as.Date(Sql.result$PDate)
  return(Sql.result)
}

##Gets Forward Data from Products.Forwards
Get.Fwd.Data <- function(mydb, temp.Region=0,temp.ContractID2=0,temp.ContractID3=0,temp.Duration=0,temp.contractID=0){
  
  Sql.Query = paste0("SELECT * FROM products.forwards WHERE `Region`='",temp.Region
                     ,"' And `ContractID2`='",temp.ContractID2,"'",
                     if(temp.ContractID3 !=0){paste0(" And `ContractID3`='",temp.ContractID3,"'")},
                     if(temp.Duration != 0){paste0(" and `ContractDuration` = '",temp.Duration,"'")},
                      if(temp.contractID != 0){paste0(" and ContractID like ('",temp.contractID,"')")}
                        ,";")
  
  print(Sql.Query)  
  Sql.result <- dbGetQuery(mydb,Sql.Query)
  Sql.result$DDate <- as.Date(Sql.result$DDate)
  #Sql.result$PDate <- as.Date(Sql.result$PDate)
  return(Sql.result)
}