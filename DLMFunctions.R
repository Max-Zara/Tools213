require(dlm)

Get.Gibbs.Params <- function(y,mod.structure,n.sample){
Gibbs.Superiore <- dlmGibbsDIG(y,mod=mod.structure,a.y=1,b.y=1,a.theta=1,b.theta=1,n.sample=n.sample,thin=1,save.states=TRUE)
burn <- 0.2*n.sample;
return(mcmcMeans(cbind(Gibbs.Superiore$dV[-(1:burn)], Gibbs.Superiore$dW[-(1:burn),])))
}

Get.Regression.Filter <- function(y,x , V=0, W = c(0,0)){
  
  if(is.null(ncol(x))){x.ncols <- 1}else{  x.ncols <- ncol(x)
  if(W==c(0,0)){W = rep(0,1+x.ncols)}
  }
  
  ##dlmModReg assumes regression coefficients follow a random walk
  buildReg <- function(u){
    dlmModReg(x, dV = exp(u[1]),dW = exp(u[2:(2+x.ncols)]))
  }
  if(V== 0 & W==c(0,0)){
  outMLE <- dlmMLE(y, parm=rep(0,2+x.ncols), buildReg)
  exp(outMLE$par)
  mod <- buildReg(outMLE$par)}else{
    mod <- buildReg(c(V,W))
  }
  outFilter <- dlmFilter(y,mod)
  outS <- dlmSmooth(y, mod)
  #plot the intercept and Beta and how they change over time
  plot(dropFirst(outS$s ))
  return(dlmFilter(y,mod))
}


chart.Filter <- function(y,dlm.filt){
  fitted.filter <- dlmFilter(y,dlm.filt)
  sdev <- residuals(fitted.filter)$sd
  out <- residuals(fitted.filter)
  beta <- 0 + cumsum(out$res^2)/2
  alpha <- 0 + (1:length(y))/2
  ncol.var <- length(fitted.filter$U.C[[1]])
  tt <- qt(0.95, df=2*alpha) #student t density with degrees of freeedom - as more obs come in - the higher the significance 
  #Ctilde <- unlist(dlmSvd2var(fitted.filter$U.C,fitted.filter$D.C))
  #seq.var <- seq(1,length(Ctilde),by=ncol.var)
  #Var.Val <- Ctilde[seq.var]
  lower <- dropFirst(fitted.filter$f)- tt*sd    #sqrt(dlm.filt$W[1,1]+dlm.filt$W[3,3])  #m is the series of state vectors
  upper <- dropFirst(fitted.filter$f)+ tt*sd    #  sqrt(dlm.filt$W[1,1]+dlm.filt$W[3,3])

  par(mfrow=c(1,1))
  plot(y, ylim=c(min(y), max(y)),
       type='o', col = "darkgrey")
  lines(dropFirst(fitted.filter$f), type='o', lty="431313", pch=3, col="blue")
  lines(dropFirst(fitted.filter$m[,1]+fitted.filter$m[,3]), type = "o", lty="dashed",pch=4, col="lightblue")
  
  lines(lower, lty = 2, lwd = 2, col="red")
  lines(upper, lty = 2, lwd = 2, col="red")
  legend("topleft", bty = "n",
         legend=c("observed", "one-step-ahead forecast - model",
                  "state estimate","bottom 5%","top 5%"),
         lty=c("solid", "431313", "dashed","dashed","dashed"), pch=c(1,3,4,1,1),
         col=c("darkgrey", rep("black", 2),rep("red",2)), inset = 0.05)
  par(mfrow=c(1,1))
}

build.MLE.AR <- function(y, AR.lags, to.chart = FALSE){
  
  level0 <- y[1]
  slope0 <- mean(diff(y))
  
  #Function for MLE to optimise - feed intial values of Level and Slope. Leave DW u[1:2] and AR u[4:..] and Sigma u[3] to solve
  buildGap <- function(u) {
    trend <- dlmModPoly(dV = 1e-7, dW = exp(u[1:2]),
                        m0 = c(level0, slope0), C0 = 2 * diag(2))
    #ARtransPars ensures you have AR < 1
    gap <- dlmModARMA(ar = ARtransPars(u[4:(3+AR.lags)]), sigma2 = exp(u[3]))
    return(trend + gap)}
  
  #Initial Values
  init <- c(-3, -1, -3, rep(.2,AR.lags))
  outMLE <- dlmMLE(y, init, buildGap)
  print(outMLE$message)
  if(outMLE$message == "ERROR: ABNORMAL_TERMINATION_IN_LNSRCH"){
    init <- c(-2, -2, -2, rep(.15,AR.lags))
    outMLE <- dlmMLE(y, init, buildGap)
    print(outMLE$message)
    if(outMLE$message == "ERROR: ABNORMAL_TERMINATION_IN_LNSRCH"){
      init <- c(-1, -3, -1, rep(.1,AR.lags))
      outMLE <- dlmMLE(y, init, buildGap)
      print(outMLE$message)
    }
  }
  dlmGap <- buildGap(outMLE$par)
  stable.check <- GG(dlmGap)[3:(AR.lags+2), 3]
  print(paste("Stability Check:",abs(sum(stable.check))<1, " Coefs: ",c(stable.check)))
  
  if(to.chart){
      chart.Filter(y,dlmGap)
  }
  
  return(dlmGap)
}



Get.Poly.Params.Gibbs <- function(y,poly.n = 2,no.runs = 12000, to.chart= FALSE){
  MCMC <- no.runs
  gibbsOut <- dlmGibbsDIG(y, mod = dlmModPoly(poly.n), a.y = 1, b.y = 1000,
                          a.theta = 10, b.theta = 1000, n.sample = MCMC,
                          thin = 1, save.states = FALSE)
  
  ## output analysis
  burn <- 2000
  
  ## Convergence assessment and autocorrelations
  if(to.chart){
    par(mfrow=c(3,3), mar=c(3.1,2.1,2.1,1.1))
    plot(gibbsOut$dV[-(1:burn)], cex=0.5, xlab="", ylab="", main=expression(V))
    for(i.W in 1:poly.n){
      plot(gibbsOut$dW[-(1:burn),i.W], cex=0.5, xlab="", ylab="", main=paste0("W",i.W))
    }
    #plot(gibbsOut$dW[-(1:burn),2], cex=0.5, xlab="", ylab="", main=expression(W[22]))
    use <- MCMC - burn
    from <- 0.05 * use
    plot(ergMean(gibbsOut$dV[-(1:burn)], from), type="l",
         xaxt="n",xlab="", ylab="")
    at <- pretty(c(0,use),n=3); at <- at[at>=from]
    axis(1, at=at-from, labels=format(at))
    for(i.W in 1:poly.n){
      plot(ergMean(gibbsOut$dW[-(1:burn),i.W], from), type="l",
           xaxt="n",xlab="", ylab="")
      at <- pretty(c(0,use),n=3); at <- at[at>=from]
      axis(1, at=at-from, labels=format(at))
    }
    acf(gibbsOut$dV[-(1:burn)], ci=0, ylim=c(0,1), ylab="", main="")
    for(i.W in 1:poly.n){
      acf(gibbsOut$dW[-(1:burn),i.W], ci=0, ylim=c(0,1), ylab="", main="")
    }
    
    readline(prompt="Press [enter] to continue")
    
    par(mfrow=c(1,3), mar=c(5.1,4.1,1,1))
    for(i.W in 1:poly.n){
      plot(gibbsOut$dV[-(1:burn)], gibbsOut$dW[-(1:burn),i.W],
         xlab=expression(V), ylab=paste0("W",i.W), pch='.')
      if(i.W > 1){
      plot(gibbsOut$dW[-(1:burn),1], gibbsOut$dW[-(1:burn),i.W],
         xlab="W1", ylab=paste0("W",i.W), pch='.')
      }
    }
  }
  
  ## posterior estimates of variances, with estimated MC standard error
  Gibbs.Results <- mcmcMeans(cbind(gibbsOut$dV[-(1:burn)], gibbsOut$dW[-(1:burn),]))
  
  return(Gibbs.Results)
}


Get.Linear.Growth.Params.Gibbs <- function(y,no.runs = 12000, to.chart= FALSE){
  MCMC <- no.runs
  gibbsOut <- dlmGibbsDIG(y, mod = dlmModPoly(2), a.y = 1, b.y = 1000,
                          a.theta = 10, b.theta = 1000, n.sample = MCMC,
                          thin = 1, save.states = FALSE)
  
  ## output analysis
  burn <- 2000
  
  ## Convergence assessment and autocorrelations
  if(to.chart){
    par(mfrow=c(3,3), mar=c(3.1,2.1,2.1,1.1))
    plot(gibbsOut$dV[-(1:burn)], cex=0.5, xlab="", ylab="", main=expression(V))
    plot(gibbsOut$dW[-(1:burn),1], cex=0.5, xlab="", ylab="", main=expression(W[11]))
    plot(gibbsOut$dW[-(1:burn),2], cex=0.5, xlab="", ylab="", main=expression(W[22]))
    use <- MCMC - burn
    from <- 0.05 * use
    plot(ergMean(gibbsOut$dV[-(1:burn)], from), type="l",
         xaxt="n",xlab="", ylab="")
    at <- pretty(c(0,use),n=3); at <- at[at>=from]
    axis(1, at=at-from, labels=format(at))
    plot(ergMean(gibbsOut$dW[-(1:burn),1], from), type="l",
         xaxt="n",xlab="", ylab="")
    at <- pretty(c(0,use),n=3); at <- at[at>=from]
    axis(1, at=at-from, labels=format(at))
    plot(ergMean(gibbsOut$dW[-(1:burn),2], from), type="l",
         xaxt="n",xlab="", ylab="")
    at <- pretty(c(0,use),n=3); at <- at[at>=from]
    axis(1, at=at-from, labels=format(at))
    acf(gibbsOut$dV[-(1:burn)], ci=0, ylim=c(0,1), ylab="", main="")
    acf(gibbsOut$dW[-(1:burn),1], ci=0, ylim=c(0,1), ylab="", main="")
    acf(gibbsOut$dW[-(1:burn),2], ci=0, ylim=c(0,1), ylab="", main="")
   
    readline(prompt="Press [enter] to continue")
     
    par(mfrow=c(1,3), mar=c(5.1,4.1,1,1))
    plot(gibbsOut$dV[-(1:burn)], gibbsOut$dW[-(1:burn),1],
         xlab=expression(V), ylab=expression(W[11]), pch='.')
    plot(gibbsOut$dV[-(1:burn)], gibbsOut$dW[-(1:burn),2],
         xlab=expression(V), ylab=expression(W[22]), pch='.')
    plot(gibbsOut$dW[-(1:burn),1], gibbsOut$dW[-(1:burn),2],
         xlab=expression(W[11]), ylab=expression(W[22]), pch='.')
  }
    
    ## posterior estimates of variances, with estimated MC standard error
  Gibbs.Results <- mcmcMeans(cbind(gibbsOut$dV[-(1:burn)], gibbsOut$dW[-(1:burn),]))
  
  return(Gibbs.Results)
}


Linear.Growth.Model <- function(y,n.poly=2,dV1,dW1,n.ahead=5,to.chart=FALSE){

mod1 <- dlmModPoly(n.poly,dV = dV1, dW = dW1)
mod1Filt <- dlmFilter(y, mod1)
fut1 <- dlmForecast(mod1Filt, n=n.ahead)

print(paste("MAE Model1: ",mean(abs(mod1Filt$f - y)))) 
print(paste("MSE Model1: ",mean((mod1Filt$f - y)^2)))
print(paste("MAE % Model1: ",mean(abs(mod1Filt$f - y) / y)))

sqrt(sum((mod1Filt$f - mod1Filt$y)[-(1:5)]^2) /
       sum(diff(mod1Filt$y[-(1:4)])^2))

  if(to.chart){
  par(mfrow=c(1,1))
  par(mar=c(2,3,1,0) + 0.1, cex=0.7)
  plot(y, ylim=c(min(y), max(y)),
        type='o', col = "darkgrey")
  lines(dropFirst(mod1Filt$f), type='o', lty="431313", pch=3)
  lines(fut1$f, lty="dashed")
  legend("topleft", bty = "n",
         legend=c("observed", "one-step-ahead forecast - model",
                   "5 steps forecast - model"),
         lty=c("solid", "431313", "dashed"), pch=c(1,3,4),
         col=c("darkgrey", rep("black", 2)), inset = 0.05)
  par(mfrow=c(1,1))
  }

return(mod1Filt)
}


Random.Walk.DF <- function(y, disc.factor = 0.7, alpha.prior = 0, beta.prior = -2, to.chart= FALSE, prob  = 0.95){
  mod <- dlmModPoly(1, dV=1)
  modFilt <- dlmFilterDF(y,mod,DF=disc.factor) #discount factor of 0.9
  
  beta0 <- beta.prior; alpha0 <- alpha.prior;
  
  #Filtering Estimates
  out <- residuals(modFilt)
  beta <- beta0 + cumsum(out$res^2)/2
  alpha <- alpha0 + (1:length(y))/2
  Ctilde <- unlist(dlmSvd2var(modFilt$U.C,modFilt$D.C))[-1]
  
  tt <- qt(prob, df=2*alpha) #student t density with degrees of freeedom - as more obs come in - the higher the significance 
  lower <- dropFirst(modFilt$m)- tt*sqrt(Ctilde * beta/alpha) #m is the series of state vectors
  upper <- dropFirst(modFilt$m)+ tt*sqrt(Ctilde * beta/alpha)
  
  if(to.chart){
    par(mfrow=c(1,1))
    plot(y, ylab = "Filtering level estimates", type = "o", 
         ylim = c(min(y), max(y)), col = "darkgray")
    lines(dropFirst(modFilt$m), type = "o")
    lines(lower, lty = 2, lwd = 2)
    lines(upper, lty = 2, lwd = 2)
  }
  return(list(modFilt))
}


#One Step ahead forecasts
Random.Walk.Forecast <- function(modFilt, alpha.prior, beta.prior, to.chart = FALSE, prob = 0.95){
  beta0 <- beta.prior; alpha0 <- alpha.prior;
  alpha <- alpha0 + (1:length(modFilt$y))/2
  out <- residuals(modFilt)
  beta <- beta0 + cumsum(out$res^2)/2
  sigma2 <- c(beta0/(alpha0-1),beta/(alpha-1)) #Original Sigma from pg151 E(sigma^2|y(1:t))=Beta/Alpha-1; Var(Sigma|y(1:t))= Beta^2/(A-1)^2(A-2)
  Qt <- out$sd^2 * sigma2[-length(sigma2)]
  alpha0T <- c(alpha0,alpha)
  tt <- qt(prob, df=2 * alpha0T[-length(alpha0T)])
  parf <- c(beta0/alpha0, beta/alpha)
  parf <- parf[-length(parf)]*out$sd^2
  
  lower <- dropFirst(modFilt$f) - tt * sqrt(parf)
  upper <- dropFirst(modFilt$f) + tt * sqrt(parf)
  
  if(to.chart){
  par(mfrow=c(1,1))
  plot(modFilt$y, ylab = "One-step-ahead forecasts", type = "o",
       ylim = c(min(modFilt$y), max(modFilt$y)), col = "darkgray")
  lines(window(modFilt$f),  type="o") 
  lines(lower, lty = 2, lwd = 2)
  lines(upper, lty = 2, lwd = 2)
  }
  return(list(Filtered.Mean=modFilt$f,Filtered.Lower = lower,Filtered.Upper = upper))
}

Random.Walk.Smoothed <- function(modFilt, alpha.prior, beta.prior, to.chart = FALSE, prob = 0.95){
  beta0 <- beta.prior; alpha0 <- alpha.prior;
  alpha <- alpha0 + (1:length(modFilt$y))/2
  out <- residuals(modFilt)
  beta <- beta0 + cumsum(out$res^2)/2

  modFilt$mod$JW <- matrix(1)
  X <- unlist(dlmSvd2var(modFilt$U.W,modFilt$D.W))[-1]
  modFilt$mod$X <- matrix(X)
  modSmooth <- dlmSmooth(modFilt)
  Stildelist <- dlmSvd2var(modSmooth$U.S,modSmooth$D.S)
  TT <- length(modFilt$y)#last value of Y
  pars <- unlist(Stildelist) * (beta[TT]/alpha[TT])
  tt <- qt(prob, df = 2 * alpha[TT])
  if(to.chart){
    par(mfrow=c(1,1))
    plot(modFilt$y, ylab = "Smoothing level estimates", 
         type = "o", ylim = c(min(modFilt$y), max(modFilt$y)), col = "darkgray")
    lines(dropFirst(modSmooth$s),  type = "o")
    lower <-dropFirst(modSmooth$s - tt * sqrt(pars))
    upper <-dropFirst(modSmooth$s + tt * sqrt(pars)) 
    lines(lower, lty = 3, lwd = 2)
    lines(upper, lty = 3, lwd = 2)
  }
  return(list(Smoothed.Mean=modSmooth$s,Smoothed.Lower = lower, Smoothed.Upper= upper ))
}

Random.Walk.Choose.DF <- function(y,vecDF,to.chart=FALSE){
  mat <- matrix(0, nrow=length(vecDF), ncol=5,dimnames = 
                  list(rep("",length(vecDF)), c("DF", "MAPE", "MAD","MSE","SD")))
  par(mar = c(2, 4, 1, 1) + 0.1, cex = 0.6)
  par(mfrow=c(2,2)) 
  i.printed=0;
  for (i in 1:length(vecDF)){
    mod <- dlmModPoly(1,dV=1)
    modFilt <- dlmFilterDF(y, mod, DF=vecDF[i])
    out <- residuals(modFilt)
    beta0 <- 20
    alpha0 <-2
    beta <- beta0 + 1/2 * cumsum(out$res^2)
    alpha <- alpha0 + (1: length(y))*(1/2)
    sigmaqhat <- beta/(alpha-1) 
    sigmaqhat <- beta/(alpha-1) 
    # Plot the estimates of sigma^2
    # ts.plot(sigmaqhat, main=paste("DF=", as.character(vecDF[i])), ylab="")  
    # Plot the one-step-ahead forecasts   
    if(to.chart){
      plot(window(y), ylab="", ylim = c(min(y,na.rm=T),max(y,na.rm=T)), lty=1, 
           main=paste("DF=", as.character(vecDF[i])), col="darkgray") 
      lines(window(modFilt$f), lty=1)  
      i.printed=i.printed+1;
      if(i.printed==4){
        readline(prompt="Press [enter] to continue")
        i.printed=0;
      }
    }
    mat[i,1] <- vecDF[i]
    mat[i,2] <- mean(abs(modFilt$f - modFilt$y)/modFilt$y,na.rm=T)
    mat[i,3] <- mean(abs(modFilt$f - modFilt$y))
    mat[i,4] <- mean((modFilt$f - modFilt$y)^2)
    mat[i,5] <- sigmaqhat[length(modFilt$y)]
  }
  
  print(round(mat,4))
  return(round(mat,4))
  par(mfrow=c(1,1))
}

dlmFilterDF <- function (y, mod, simplify = FALSE, DF) 
{
  ## storage.mode(y) <- "double"
  mod1 <- mod
  yAttr <- attributes(y)
  ytsp <- tsp(y)
  y <- as.matrix(y)
  timeNames <- dimnames(y)[[1]]
  stateNames <- names(mod$m0)
  m <- rbind(mod$m0, matrix(0, nr = nrow(y), nc = length(mod$m0)))
  a <- matrix(0, nr = nrow(y), nc = length(mod$m0))
  f <- matrix(0, nr = nrow(y), nc = ncol(y))
  U.C <- vector(1 + nrow(y), mode = "list")
  D.C <- matrix(0, 1 + nrow(y), length(mod$m0))
  U.R <- vector(nrow(y), mode = "list")
  D.R <- matrix(0, nrow(y), length(mod$m0))
  U.W <- vector(nrow(y), mode = "list")
  D.W <- matrix(0, nrow(y), length(mod$m0))
  Wliste <- vector(nrow(y), mode = "list")
  P <- vector(nrow(y), mode = "list")
  tmp <- La.svd(mod$V, nu = 0)
  Uv <- t(tmp$vt)
  Dv <- sqrt(tmp$d)
  Dv.inv <- 1/Dv
  Dv.inv[abs(Dv.inv) == Inf] <- 0
  sqrtVinv <- Dv.inv * t(Uv)
  sqrtV <- Dv * Uv
  tmp <- La.svd(mod$C0, nu = 0)
  U.C[[1]] <- t(tmp$vt)
  D.C[1, ] <- sqrt(tmp$d)
  for (i in seq(length = nrow(y))) {
    tF.Vinv <- t(mod$FF) %*% crossprod(sqrtVinv)
    a[i, ] <- mod$GG %*% m[i, ]
    P[[i]] <- mod$GG %*% crossprod(D.C[i,] * t(U.C[[i]])) %*% t(mod$GG)
    Wliste[[i]] <- P[[i]]* ((1-DF)/DF)
    svdW <- La.svd( Wliste[[i]] , nu = 0)
    sqrtW <- sqrt(svdW$d) * svdW$vt
    U.W[[i]] <- t(svdW$vt)
    D.W[i, ] <- sqrt(svdW$d)
    tmp <- La.svd(rbind(D.C[i, ] * t(mod$GG %*% U.C[[i]]), 
                        sqrtW), nu = 0)
    U.R[[i]] <- t(tmp$vt)
    D.R[i, ] <- tmp$d
    f[i, ] <- mod$FF %*% a[i, ]
    D.Rinv <- 1/D.R[i, ]
    D.Rinv[abs(D.Rinv) == Inf] <- 0
    tmp <- La.svd(rbind(sqrtVinv %*% mod$FF %*% U.R[[i]], 
                        diag(x = D.Rinv, nrow = length(D.Rinv))), nu = 0)
    U.C[[i + 1]] <- U.R[[i]] %*% t(tmp$vt)
    foo <- 1/tmp$d
    foo[abs(foo) == Inf] <- 0
    D.C[i + 1, ] <- foo
    m[i + 1, ] <- a[i, ] + crossprod(D.C[i + 1, ] * t(U.C[[i
                                                           + 1]])) %*% tF.Vinv %*% as.matrix(y[i, ] - f[i,])
  }        
  m <- drop(m)
  a <- drop(a)
  f <- drop(f)
  attributes(f) <- yAttr
  ans <- list(m = m, U.C = U.C, D.C = D.C, a = a, U.R = U.R, 
              D.R = D.R, f = f, U.W=U.W, D.W=D.W)
  ans$m <- drop(ans$m)
  ans$a <- drop(ans$a)
  ans$f <- drop(ans$f)
  attributes(ans$f) <- yAttr
  if (!is.null(ytsp)) {
    tsp(ans$a) <- ytsp
    tsp(ans$m) <- c(ytsp[1] - 1/ytsp[3], ytsp[2:3])
    class(ans$a) <- class(ans$m) <- if (length(mod$m0) > 
                                        1) 
      c("mts", "ts")
    else "ts"
  }
  if (!(is.null(timeNames) && is.null(stateNames))) {
    dimnames(ans$a) <- list(timeNames, stateNames)
    dimnames(ans$m) <- list(if (is.null(timeNames)) NULL else c("", 
                                                                timeNames), stateNames)
  }
  if (simplify) 
    return(c(mod = list(mod1), ans))
  else {
    attributes(y) <- yAttr
    ans <- c(y = list(y), mod = list(mod1), ans)
    class(ans) <- "dlmFiltered"
    return(ans)
  }
}

#################################################################
###                                                           ###
###       Unequal variances - t-distributed innovations       ###
###                                                           ###
#################################################################

dlmGibbsDIGt <- function(y, mod, A_y, B_y, A_theta = A_y, B_theta = A_y,
                         nuRange = c(1 : 10, seq(20, 100, by = 10)),
                         alpha_y = 1 / length(nuRange),
                         alpha_theta,
                         n.sample = 1,
                         thin = 0, ind, save.states = FALSE,
                         progressBar = interactive())
#################################################################
#################################################################
###                                                           ###
### y      : data (univariate time series)                    ###
### mod    : model skeleton                                   ###    
### A_y    : upper bound for a_y (prior)                      ###
### B_y    : upper bound for b_y (prior)                      ###
### A_theta: upper bounds for the components of a_theta,      ###
###            recycled if needed (prior)                     ###
### B_theta: upper bounds for the components of b_theta,      ###
###            recycled if needed (prior)                     ###
### nuRange: set of possible values for the degree-of-freedom ###
###            parameter                                      ###
### alpha_y: parameter of the Dirichlet distribution of pi_y  ###
###            (prior)                                        ### 
### alpha_theta: vector parameters of the Dirichlet           ###
###            distributions of the pi_theta's, stored as     ###
###            columns of a matrix; recycled if vector of     ###
###            length 1 or length(nuRange) (prior)            ###
### n.sample: number of MCMC to output                        ###
### thin  : number of MCMC iteration to skip between two      ###        
###           consecutived stored draws                       ###
### ind   : vector of integers specifying the position of the ###
###           unknown, nonconstant variances. If not          ###
###           specified, all are considered unknown           ###
### save.states: logical; should the simulated states be      ###
###           included in the output?                         ###
### progressBar: logical                                      ###
###                                                           ###
#################################################################
#################################################################
{
    m <- NCOL(y)
    nobs <- NROW(y)
    mod$JW <- matrix(0, nrow = ncol(mod$FF), ncol = ncol(mod$FF))
    r <- ncol(mod$FF)
    if ( hasArg(ind) ) {
        ind <- sort(unique(as.integer(ind)))
        s <- 1 : r
        perm <- s[c(ind, s[ !(s %in% ind)])]
        FF(mod) <- mod$FF[, perm, drop = FALSE]
        GG(mod) <- mod$GG[perm, perm, drop = FALSE]
        mod$W <- mod$W[perm, perm, drop = FALSE]
        p <- length(ind)
    }
    else {
        perm <- 1 : r
        p <- r
    }
    diag(mod$JW)[ 1 : p ] <- 1 : p
    mod$JV <- matrix(p + 1)
    mod$X <- matrix(0, nrow = nobs, ncol = p + 1)
    K <- length(nuRange)
    if ( is.numeric(thin) && (thin <- as.integer(thin)) >= 0 ) 
    {
        every <- thin + 1
        mcmc <- n.sample * every
    }
    else
        stop("\"thin\" must be a nonnegative integer")
    if (!all(c(length(A_theta), length(B_theta)) %in% c(1,p)))
        warning("Unexpected length of \"A_theta\" and/or \"B_theta\"")
    A_theta <- rep(A_theta, length.out = p)
    B_theta <- rep(B_theta, length.out = p)
    ablim_theta <- cbind(A_theta, B_theta)
    ablim_y <- c(A_y, B_y)
    if ( !(length(alpha_y) %in% c(1,K)) )
        warning("Unexpected length of \"alpha_y\"")
    alpha_y <- rep(alpha_y, length.out = K)
    if ( hasArg(alpha_theta) ) {
        if ( is.matrix(alpha_theta) ) {
            if ( nrow(alpha_theta) != K || ncol(alpha_theta) != p )
                stop("Wrong dimension of \"alpha_theta\"")
        } else {
            if ( !(length(alpha_theta) %in% c(1,K)) )
                warning("Unexpected length of \"alpha_theta\"")
            alpha_theta <- matrix(rep(alpha_theta, length.out = K*p), nr = K)
        }
    }
    else 
        alpha_theta <- matrix(rep(alpha_y, length.out = K*p), nr = K)
    theta <- matrix(0, nobs + 1, nrow(mod$W))
    ## initialize
    omega_y <- rep(1, nobs)
    nu_y <- rep(100, nobs)
    lambda_y <- 1
    ab_y <- 0.5 * ablim_y
    pi_y <- alpha_y / sum(alpha_y)
    ##
    omega_theta <- matrix(1, nrow = nobs, ncol = p)
    nu_theta <- matrix(100, nrow = nobs, ncol = p)
    lambda_theta <- rep(1, p)
    ab_theta <- 0.5 * ablim_theta
    pi_theta <- alpha_theta / rep(colSums(alpha_theta), each = K) 

    ## memory allocation
    if ( save.states ) 
        gibbsTheta <- array(0, dim = c(dim(theta), n.sample))
    gibbsOmega_y <- matrix(0, nrow = n.sample, ncol = nobs)
    gibbsNu_y <- matrix(0, nrow = n.sample, ncol = nobs)
    gibbsLambda_y <- numeric(n.sample)
    gibbsAB_y <- matrix(0, nrow = n.sample, ncol = 2)
    gibbsPi_y <- matrix(0, nrow = n.sample, ncol = length(nuRange))
    ##
    gibbsOmega_theta <- array(0, dim = c(nobs, p, n.sample))
    gibbsNu_theta <- array(0, dim = c(nobs, p, n.sample))
    gibbsLambda_theta <- matrix(0, nrow = n.sample, ncol = p)
    gibbsAB_theta <- array(0, dim = c(p, 2, n.sample))
    gibbsPi_theta <- array(0, dim = c(K, p, n.sample))

    ## log target and support for use with `arms'
    ldens.ab <- function(x, lambda, ...) {
        rate <- x[1] / x[2]
        dgamma(lambda, rate * x[1], rate, log = TRUE)
    }
    ind.ab <- function(x, ablim, ...) {
        all( x > 1e-6 & x < ablim )
    }

    ## draw from dirichlet distribution, from package gtools
    rdirichlet <- function (n, alpha) 
    {
        l <- length(alpha)
        x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
        sm <- x %*% rep(1, l)
        x/as.vector(sm)
    }

    it.save <- 0
    if ( progressBar )
        pb <- txtProgressBar(0, mcmc, style = 3)
    for (it in 1 : mcmc)
    {
        if ( progressBar )
            setTxtProgressBar(pb, it)
        mod$X[] <- 1 / cbind( omega_theta * rep(lambda_theta, each = nobs),
                             omega_y * lambda_y)

        ## generate states - FFBS
        modFilt <- dlmFilter(y, mod, simplify=TRUE)
        theta[] <- dlmBSample(modFilt)

        ## generate omega_y
        Sq.y.center <- (y - tcrossprod(theta[-1,],mod$FF))^2
        omega_y <- rgamma(nobs, 0.5 * (1 + nu_y),
                                           0.5 * (lambda_y * Sq.y.center + nu_y))
            
        ## generate nu_y
        for (i in 1:nobs) {
            probs <- dgamma(omega_y[i], 0.5 * nuRange, 0.5 * nuRange) * pi_y
            nu_y[i] <- sample(nuRange, size = 1, prob = probs)
        }

        ## generate pi_y
        nuTable <- table(factor(nu_y, levels=nuRange))
        pi_y <- rdirichlet(1, nuTable + alpha_y)
        
        ## generate lambda_y
        SSy <- crossprod(Sq.y.center, omega_y)
        u <- ab_y[1] / ab_y[2]
        lambda_y <- rgamma(1, u * ab_y[1] + 0.5 * nobs,
                           u + 0.5 * SSy)
        
        ## generate a_y and b_y
        ab_y <- arms(ab_y, myldens = ldens.ab, indFunc = ind.ab, n = 1,
                     lambda = lambda_y, ablim = ablim_y)
        
        ## same story for the 'theta' parameters
        ## omega_theta_i
        Sq.theta.center <- (theta[-1,1:p] - tcrossprod(theta[-(nobs + 1), ],
                                                       mod$GG)[,1:p])^2
        omega_theta[] <- rgamma(nobs * p, 0.5 * (1 + nu_theta),
                                0.5 * (rep(lambda_theta, each=nobs) * Sq.theta.center +
                                       nu_theta))
        ## nu_theta_i
        for (j in 1 : nobs)
            for (i in 1 : p)
            {
                probs <- dgamma(omega_theta[j,i], 0.5 * nuRange, 0.5 * nuRange) *
                    pi_theta[, i]
                nu_theta[j, i] <- sample(nuRange, size = 1, prob = probs)
            }

        ## pi_theta_i
        for (i in 1 : p)
        {
            nuTable <- table(factor(nu_theta[, i], levels = nuRange))
            pi_theta[, i] <- rdirichlet(1, nuTable + alpha_theta[, i])
        }

        ## lambda_theta_i
        SStheta <- colSums(Sq.theta.center * omega_theta)
        u <- ab_theta[, 1] / ab_theta[, 2]
        lambda_theta <- rgamma(p, u * ab_theta + 0.5 * nobs,
                               u + 0.5 * SStheta)

        ## a_theta_i & b_theta_i
        for (i in 1 : p)
            ab_theta[i, ] <- arms(ab_theta[i, ], myldens = ldens.ab,
                                 indFunc = ind.ab, n = 1,
                                 lambda = lambda_theta[i], ablim = ablim_theta[i, ])
        
        ## save
        if ( !(it %% every) )
        {
            it.save <- it.save + 1
            if ( save.states )
                gibbsTheta[,,it.save] <- theta
            gibbsOmega_y[it.save,] <- omega_y
            gibbsNu_y[it.save,] <- nu_y
            gibbsPi_y[it.save,] <- pi_y
            gibbsLambda_y[it.save] <- lambda_y
            gibbsAB_y[it.save,] <- ab_y
            ##
            gibbsOmega_theta[,,it.save] <- omega_theta
            gibbsNu_theta[,,it.save] <- nu_theta
            gibbsPi_theta[,,it.save] <- pi_theta
            gibbsLambda_theta[it.save,] <- lambda_theta
            gibbsAB_theta[,,it.save] <- ab_theta
       }
    }

    if ( progressBar )
        close(pb)
    if ( save.states )
        return(list(omega_y = gibbsOmega_y, nu_y = gibbsNu_y,
                    pi_y = gibbsPi_y, lambda_y = gibbsLambda_y, ab_y = gibbsAB_y,
                    omega_theta = gibbsOmega_theta, nu_theta = gibbsNu_theta,
                    pi_theta = gibbsPi_theta, lambda_theta = gibbsLambda_theta,
                    ab_theta = gibbsAB_theta,
                    theta = gibbsTheta[, order(perm), , drop = FALSE]))
    else
        return(list(omega_y = gibbsOmega_y, nu_y = gibbsNu_y,
                    pi_y = gibbsPi_y, lambda_y = gibbsLambda_y, ab_y = gibbsAB_y,
                    omega_theta = gibbsOmega_theta, nu_theta = gibbsNu_theta,
                    pi_theta = gibbsPi_theta, lambda_theta = gibbsLambda_theta,
                    ab_theta = gibbsAB_theta))

#### end dlmGibbsDIGt ####
}
