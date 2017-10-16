#Function Tools List

library(tm)
library(SnowballC)
library(data.table)
library(ggplot2)

#replaces punctuation and adds a space " ", this way "hello...Max" becomes "hello Max", and not "helloMax"

myscree <- function(eigs,x=0.8,y=0.1,just=c("right","bottom")){
  vp <- viewport(x=x,y=y,width = 0.2, height=0.2,just=just)
  sp <- qplot(factor(1:length(eigs)),eigs,geom="bar",stat="identity") + labs(x=NULL,y=NULL)
  print(sp,vp=vp) }

movsum <- function(x,n=25){filter(x,rep(1,n), sides=1)}


plot.fit.corr <- function(xvar,yvar,xname,yname,ititle) 
{
  chartdata <- data.frame(xname = xvar, yname = yvar)
  
  if(min(yvar)<0 & max(yvar)>0)
  {
    ggplot(chartdata,aes(x=xname, y=yname)) +
      geom_point(shape=1) +    # Use hollow circles
      geom_smooth(method=lm)+
    geom_hline(aes(yintercept=0)) +labs(x=xname,y=yname) +
    ggtitle(paste(ititle,"| Cor of ",xname," vs ",yname,": ",round(cor(xvar,yvar),2))) 
  }   else   {
    ggplot(chartdata,aes(x=xname, y=yname)) +
      geom_point(shape=1) +    # Use hollow circles
      geom_smooth(method=lm)+ labs(x=xname,y=yname) +
    ggtitle(paste(ititle,"| Cor of ",xname," vs ",yname,": ",round(cor(xvar,yvar),2))) 
  }
}


library(reshape2)

panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}


plot.analyse<- function(x, xname, YLIM=0.48) {
  par(mfrow=c(3,3))
  
  plot(x, main = paste("a.Series(",xname,")"),col="blue",lwd=1, type = 'l')
  grid(lty="dotted",col="gray75")
  
  plot(density(x), main = paste("b.Density(",xname,")"),col="blue",lwd=2)
  
  #renorm the x
  newx <- (x-mean(x))/sd(x)
  
  ran = rnorm(1000000)
  plot(density(ran), ylim=c(0,YLIM), main = paste("c.Density(",xname," vs Normal)"),
       xlim=c(-4,4))
  polygon(density(ran),col="burlywood")
  lines(density(newx),col="blue",lwd=2)
  
  qqnorm(x, main = paste("d. QQPlot vs Normal of ",xname),col="blue",ylab = "tip quantiles")
  qqline(x, col="burlywood", lwd = 2)
  grid(lty="dotted",col="gray75")
  
  logx <- log10(x)
  qqnorm(logx,main=paste("e. QQPlot vs Normal of Log10(",xname,")"),col="blue4")
  qqline(logx,col="burlywood3",lwd=2)
  grid(lty="dotted",col="gray75")
  
  aboveX <- 1/(x)
  qqnorm(aboveX,main=paste("f. QQPlot vs Normal of 1/(",xname,")"),col="blue4")
  qqline(aboveX,col="burlywood3",lwd=2)
  grid(lty="dotted",col="gray75")

  logAbove <- log10(aboveX)
  qqnorm(logAbove,main=paste("g. QQPlot vs Normal of Log10(1/(",xname,"))"),col="blue4")
  qqline(logAbove,col="burlywood3",lwd=2)
  grid(lty="dotted",col="gray75")
  
  sqX <- sqrt(x)
  qqnorm(sqX,main=paste("h. QQPlot vs Normal of Sq.Rt(",xname,")"),col="blue4")
  qqline(sqX,col="burlywood3",lwd=2)
  grid(lty="dotted",col="gray75")
  
  logSq <- log10(sqX)
  qqnorm(logSq,main=paste("i. QQPlot vs Normal of Log10(Sq.Rt(",xname,"))"),col="blue4")
  qqline(logSq,col="burlywood3",lwd=2)
  grid(lty="dotted",col="gray75")
}

plot.dens <- function (x,YLIM=0.48){
  #Plot Density vs Normal Distribution
  
  #renorm the x
  newx <- (x-mean(x))/sd(x)
  
  ran = rnorm(1000000)
  plot(density(ran), ylim=c(0,YLIM), main = "Density of Residuals vs Normal Distribution",
       xlim=c(-4,4))
  polygon(density(ran),col="burlywood")
  lines(density(newx),col="blue",lwd=2)
}

MAE <- function(actual, predicted) {mean(abs(actual - predicted))}

Get.Unique.Country.Combo <- function(all.data){
unique(as.data.frame(all.data[2])[,c("Type","Country")])}

plot.qq <- function(x,Header="QQ plot", logtransform = F) {
 #Plot the Quartiles to see how it compares to the normal distribution 
  if(logtransform)
    x = log10(x)
  
  qqnorm(x, main = Header,col="blue",ylab = "tip quantiles")
  qqline(x,col="burlywood",lwd=2)
  grid(lty="dotted",col="gray75")
}


getname <- function(x) { deparse(substitute(x))}


plot.regress <- function (x, var1, var2, plotresiduals = FALSE){
  #This function plots two variables in a dataset and plots residuals vs the Best Fit line
  attach(x)
  plot(x[,var1], x[,var2], main = paste("Regress",colnames(x)[var1]," vs ",colnames(x)[var2]), pch=20,col="deepskyblue", xlab = colnames(x)[var1], ylab = colnames(x)[var2])
  abline(lm(x[,var2]~x[,var1]),col="dodgerblue4",lty=1,lwd=2) #writes over last plot
  grid(col="gray70")
  
  if(plotresiduals)
  {
    mod1 <- lm(x[,var2]~x[,var1])
    res <- signif(residuals(mod1))
    pre <- predict(mod1) #plot distances between points and the regression line
    segments(x[,var1],x[,var2],x[,var1],pre,col="red")
    library(calibrate)
    textxy(x[,var1],x[,var2],res,cx=0.7)
  }
  detach(x)
}

#Cleans the text by removing Numbers, Punctuation and makes it lower case
CleanStemFunction <- function(x, removeNumbers = FALSE) {
  x <- tm_map(x,content_transformer(tolower))
  
  #Some numbers provide useful information, however most are likely to be unique: UNLESS accompanied by specific ID
  if(removeNumbers)
  x <- tm_map(x,removeNumbers)
  
  #remove words like ours, we, my, me
  x <- tm_map(x,removeWords,stopwords())
  
  x <- tm_map(x,replacePunctuation)
  
  #combine words like learns, learned, learning
  x <- tm_map(x,stemDocument)
  
  x <- tm_map(x,stripWhitespace)
}
#Run this with the Clean Stem Function
replacePunctuation <- function(x) {
  gsub("[[:punct:]]+"," ",x)
}
#stopwords = function(x) {removeWords(x,stopwords())}


MAE <- function(actual, predicted) {mean(abs(actual - predicted))}




normalize <- function(x, optional.rescale = 0) { #optional rescale is how much you want to widen the range on each side
  return((x-min(x)+optional.rescale*(max(x)-min(x)))/((max(x)-min(x))*(1+optional.rescale*2)))
}

#x is the data source, Y is the original data series used for normalizing
renormalize <- function (x, y, optional.rescale = 0){ #optional rescale is how much you want to widen the range on each side
   return (x*(max(y)-min(y))*(1+optional.rescale*2) + min(y)-(optional.rescale*(max(y)-min(y))))
}

#ANALYSIS

panel.cor <- function(x, y, digits = 2, cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  # correlation coefficient
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste("r= ", txt, sep = "")
  text(0.5, 0.6, txt)
  
  # p-value calculation
  p <- cor.test(x, y)$p.value
  txt2 <- format(c(p, 0.123456789), digits = digits)[1]
  txt2 <- paste("p= ", txt2, sep = "")
  if(p<0.01) txt2 <- paste("p= ", "<0.01", sep = "")
  text(0.5, 0.4, txt2)
}



reg <- function(y,x) {
  x <- as.matrix(x)
  x <- cbind(Intercept = 1, x)
  b <- solve(t(x) %*% x) %*% t(x) %*% y
  colnames(b) <- "estimate"
  print(b)
}


plot2dens <- function(x1,x2,heading="Density Comparison",label1="First Series",label2 = "Second Series"){
  max1 <- max(density(x1)$y)
  max2 <- max(density(x2)$y)
  
  plot(density(as.vector(x2)),main=heading)
  polygon(density(x1),col="burlywood")
  lines(density(as.vector(x2)),col="blue",lwd=2)
  legend("topleft",c(label2,label1),col=c("blue","burlywood"),bty="n",pch=c(21,21), pt.bg=c("blue","burlywood"))
}



# Function has.interaction checks whether x is part of a term in terms
# terms is a vector with names of terms from a model

has.interaction <- function(x,terms){
  out <- sapply(terms,function(i){
    sum(1-(strsplit(x,":")[[1]] %in% strsplit(i,":")[[1]]))==0
  })
  return(sum(out)>0)
}


# Function Model.select
# model is the lm object of the full model
# keep is a list of model terms to keep in the model at all times
# sig gives the significance for removal of a variable. Can be 0.1 too (see SPSS)
# verbose=T gives the F-tests, dropped var and resulting model after 
model.select <- function(model,keep,sig=0.05,verbose=F){
  counter=1
  # check input
  if(!is(model,"lm")) stop(paste(deparse(substitute(model)),"is not an lm object\n"))
  # calculate scope for drop1 function
  terms <- attr(model$terms,"term.labels")
  if(missing(keep)){ # set scopevars to all terms
    scopevars <- terms
  } else{            # select the scopevars if keep is used
    index <- match(keep,terms)
    # check if all is specified correctly
    if(sum(is.na(index))>0){
      novar <- keep[is.na(index)]
      warning(paste(
        c(novar,"cannot be found in the model",
          "\nThese terms are ignored in the model selection."),
        collapse=" "))
      index <- as.vector(na.omit(index))
    }
    scopevars <- terms[-index]
  }
  
  # Backward model selection : 
  
  while(T){
    # extract the test statistics from drop.
    test <- drop1(model, scope=scopevars,test="F")
    
    if(verbose){
      cat("-------------STEP ",counter,"-------------\n",
          "The drop statistics : \n")
      print(test)
    }
    
    pval <- test[,dim(test)[2]]
    
    names(pval) <- rownames(test)
    pval <- sort(pval,decreasing=T)
    
    if(sum(is.na(pval))>0) stop(paste("Model",
                                      deparse(substitute(model)),"is invalid. Check if all coefficients are estimated."))
    
    # check if all significant
    if(pval[1]<sig) break # stops the loop if all remaining vars are sign.
    
    # select var to drop
    i=1
    while(T){
      dropvar <- names(pval)[i]
      check.terms <- terms[-match(dropvar,terms)]
      x <- has.interaction(dropvar,check.terms)
      if(x){i=i+1;next} else {break}              
    } # end while(T) drop var
    
    if(pval[i]<sig) break # stops the loop if var to remove is significant
    
    if(verbose){
      cat("\n--------\nTerm dropped in step",counter,":",dropvar,"\n--------\n\n")              
    }
    
    #update terms, scopevars and model
    scopevars <- scopevars[-match(dropvar,scopevars)]
    terms <- terms[-match(dropvar,terms)]
    
    formul <- as.formula(paste(".~.-",dropvar))
    model <- update(model,formul)
    
    if(length(scopevars)==0) {
      warning("All variables are thrown out of the model.\n",
              "No model could be specified.")
      return()
    }
    counter=counter+1
  } # end while(T) main loop
  return(model)
}


data.loading <- function(tickers, start.date, end.date)
{
  # Change the locale
  sl <- Sys.setlocale(locale="US")
  
  # Create the universe of dates
  all.dates <- seq(as.Date(start.date), as.Date(end.date), by="day")
  all.dates <- subset(all.dates,weekdays(all.dates) != "Sunday" & weekdays(all.dates) != "Saturday")
  all.dates.char <- as.matrix(as.character(all.dates))
  
  # Create sparse matrices
  open <- matrix(NA, NROW(all.dates.char), length(tickers))
  hi <- open
  low <- open
  close <- open
  volume <- open
  adj.close <- open
  
  # Name the rows correctly
  rownames(open) <- all.dates.char
  rownames(hi) <- all.dates.char
  rownames(low) <- all.dates.char
  rownames(close) <- all.dates.char
  rownames(volume) <- all.dates.char
  rownames(adj.close) <- all.dates.char
  
  # Split the start and end dates to be used in the ULR later on
  splt <- unlist(strsplit(start.date, "-"))
  a <- as.character(as.numeric(splt[2])-1)
  b <- splt[3]
  c <- splt[1]
  
  splt <- unlist(strsplit(end.date, "-"))
  d <- as.character(as.numeric(splt[2])-1)
  e <- splt[3]
  f <- splt[1]
  
  # Create the two out of the three basic components for the URL loading
  str1 <- "http://ichart.finance.yahoo.com/table.csv?s="
  str3 <- paste("&a=", a, "&b=", b, "&c=", c, "&d=", d, "&e=", e, "&f=", f, "&g=d&ignore=.csv", sep="")
  
  # Main loop for all assets
  for (i in seq(1,length(tickers),1))
  {
    str2 <- tickers[i]
    strx <- paste(str1,str2,str3,sep="")
    x <- read.csv(strx)
    
    datess <- as.matrix(x[1])
    
    replacing <- match(datess, all.dates.char)
    open[replacing,i] <- as.matrix(x[2])
    hi[replacing,i] <- as.matrix(x[3])
    low[replacing,i] <- as.matrix(x[4])
    close[replacing,i] <- as.matrix(x[5])
    volume[replacing,i] <- as.matrix(x[6])
    adj.close[replacing,i] <- as.matrix(x[7])
  }
  
  # Name the cols correctly
  colnames(open) <- tickers
  colnames(hi) <- tickers
  colnames(low) <- tickers
  colnames(close) <- tickers
  colnames(volume) <- tickers
  colnames(adj.close) <- tickers
  
  # Return the ouput
  return(list(open=open, high=hi, low=low, close=close, volume=volume, adj.close=adj.close))
}

GMMGaussian <- function(x){
  
  T <- length(x)
  
  mu_hat <- (1/T)*sum(x)
  sigma_hat <- sqrt((1/T)* sum(x^2)-mu_hat^2)
  
  mom <- c(x - mu_hat, x^2 - mu_hat^2 - sigma_hat^2)
  
  f <- matrix(mom, nrow = 2) 
  d <- matrix(c(-1,  0,  -2*mu_hat,  -2*sigma_hat), nrow =2, ncol =2)
  S <- (1/T)*(f%*%t(f))
  
  V <- solve(t(d) %*% solve(S) %*% d)
  
  result <- list(mu_hat = mu_hat, sigma_hat = sigma_hat,SE_mu = sqrt((1/T)*V[1,1]),SE_sigma = sqrt((1/T)*V[2,2]))
  
  return(result)
}

#Statistic for Bootstrap of Regressors
boot.regress.stat <- function(formula, data, indices){
  d <- data[indices,] #allows boot to select sample
  fit <- lm(formula, data=d)
  return(coef(fit))
}

#In -sample difference in means test by bootstrap
diff.means <- function(d,f){
  n <- nrow(d)
  idx1 <- 1:table(as.numeric(d$fomc))[2]
  idx2 <- seq(length(idx1)+1,n)
  m1 <- sum(d[idx1,1]*f[idx1])/sum(f[idx1])
  m2 <- sum(d[idx2,1]*f[idx2])/sum(f[idx2])
  
  ss1 <- sum(d[idx1,1]^2 * f[idx1])- (m1^2 * sum(f[idx1]))
  ss2 <- sum(d[idx2,1]^2 * f[idx2])- (m2^2 * sum(f[idx2]))
  
  c(m1 - m2, (ss1+ss2)/(sum(f)-2))
  
}

split.data.random <- function(data,p=0.7,s=666) {
  set.seed(s)
  index = sample(1:dim(data)[1])
  train <- data[index[1:floor(dim(data)[1]*p)],]
  test <- data[index[((ceiling(dim(data)[1]*p))+1):dim(data)[1]],]
  return(list(train=train,test=test))
}

split.data.ordered <- function(data,p=0.7) {
  index = (1:dim(data)[1])
  train <- data[index[1:floor(dim(data)[1]*p)],]
  test <- data[index[((ceiling(dim(data)[1]*p))+1):dim(data)[1]],]
  return(list(train=train,test=test))
}

reset_par <- function(){
  op <- structure(list(xlog = FALSE, ylog = FALSE, adj = 0.5, ann = TRUE,
                       ask = FALSE, bg = "transparent", bty = "o", cex = 1, cex.axis = 1,
                       cex.lab = 1, cex.main = 1.2, cex.sub = 1, col = "black",
                       col.axis = "black", col.lab = "black", col.main = "black",
                       col.sub = "black", crt = 0, err = 0L, family = "", fg = "black",
                       fig = c(0, 1, 0, 1), fin = c(6.99999895833333, 6.99999895833333
                       ), font = 1L, font.axis = 1L, font.lab = 1L, font.main = 2L,
                       font.sub = 1L, lab = c(5L, 5L, 7L), las = 0L, lend = "round",
                       lheight = 1, ljoin = "round", lmitre = 10, lty = "solid",
                       lwd = 1, mai = c(1.02, 0.82, 0.82, 0.42), mar = c(5.1, 4.1,
                                                                         4.1, 2.1), mex = 1, mfcol = c(1L, 1L), mfg = c(1L, 1L, 1L,
                                                                                                                        1L), mfrow = c(1L, 1L), mgp = c(3, 1, 0), mkh = 0.001, new = FALSE,
                       oma = c(0, 0, 0, 0), omd = c(0, 1, 0, 1), omi = c(0, 0, 0,
                                                                         0), pch = 1L, pin = c(5.75999895833333, 5.15999895833333),
                       plt = c(0.117142874574832, 0.939999991071427, 0.145714307397962,
                               0.882857125425167), ps = 12L, pty = "m", smo = 1, srt = 0,
                       tck = NA_real_, tcl = -0.5, usr = c(0.568, 1.432, 0.568,
                                                           1.432), xaxp = c(0.6, 1.4, 4), xaxs = "r", xaxt = "s", xpd = FALSE,
                       yaxp = c(0.6, 1.4, 4), yaxs = "r", yaxt = "s", ylbias = 0.2), .Names = c("xlog",
                                                                                                "ylog", "adj", "ann", "ask", "bg", "bty", "cex", "cex.axis",
                                                                                                "cex.lab", "cex.main", "cex.sub", "col", "col.axis", "col.lab",
                                                                                                "col.main", "col.sub", "crt", "err", "family", "fg", "fig", "fin",
                                                                                                "font", "font.axis", "font.lab", "font.main", "font.sub", "lab",
                                                                                                "las", "lend", "lheight", "ljoin", "lmitre", "lty", "lwd", "mai",
                                                                                                "mar", "mex", "mfcol", "mfg", "mfrow", "mgp", "mkh", "new", "oma",
                                                                                                "omd", "omi", "pch", "pin", "plt", "ps", "pty", "smo", "srt",
                                                                                                "tck", "tcl", "usr", "xaxp", "xaxs", "xaxt", "xpd", "yaxp", "yaxs",
                                                                                                "yaxt", "ylbias"))
  par(op)
}