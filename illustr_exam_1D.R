setwd("~/Documents/Recherche/Article_high-dim_CGP/R codes/1D cases")
source('all_base_functions_1D.R')
source('all_models_1D.R')
library(ggplot2) ## plot zoom inset
library(grid) ## viewport function
library(DiceDesign) ## lhsDesign

## choice of the numerical example
## put 'yes' between quotations
multi_monot_convex_bd_syn = 'yes' # for multiple linear inequalities (monotonicity, convexity and boundedness)
multi_monot_bd_syn = '' # for multiple linear inequalities (monotonicity and constant upper bound)
mult_lin_syn = '' # for multiple linear inequalities: monotonicity and functional lower and upper bounds
mon_syn = '' # for synthetic monotone nondecreasing functions with fixed hyperparameter (LS-ESS and HMC)
bound_syn = '' # for boundedness constraints [a;b] using LS-ESS and HMC
comp_mon_HMC_LS_WC = '' # for synthetic monotone nonincreasing functions with fixed hyperparameter (LS-ESS, LS-WC-FFT and HMC)
runtime_mon_HMC_LS = '' # run-time/iteration
runtime_bound_HMC_LS = '' # run-time/iteration
decreasing_syn = '' # for monotone decreasing toy examples
####################################################


if (multi_monot_convex_bd_syn == 'yes'){
  ###### true monotone and convex function 
  f <- function(x){
    x^2 #  inc and convex fct and bounded between 0 & 1
    # 1/(1+exp(-4*x)) # inc and concave fct
    # exp(-x) # dec & convex fct
    # -x^2 # dec & concave
  }
  ## data
  set.seed(12345)
  n <- 100
  sigN <- 0.1 # noise sd
  xtr <- # lhsDesign(n=n,dimension=1,seed=12345)$design
    runif(n,0,1)
  ytr <- f(xtr) + rnorm(n, 0, sd=sigN)
  N1 <- 9#10# # size of the 1st subdomain
  M <- 5#10# # nb of subdomains
  N <- N1*M # nb of knot points
  eta <- 5000 # smooth approximate parameter
  nsim <- 5000 # nb of mcmc retained samples
  brn <- 1000 # burn in
  thin <- 1
  nu <- 2.5 # smoothness kernel parameter
  l <- round(l_est(nu,c(0,1),0.05),2) # length-scale
  upper <- 1 # upper bound for boundedness constr
  
  ## LS-ESS with only monotonicity
  post_lin_LS <- linCGP.ESS(y=ytr,x=xtr,N1=N1,M=M,nu=nu,l=l,eta=eta,nsim=nsim,burn.in=brn,thin=thin,
                            sig.in=sigN^2,tau.in=1,constrType=c('increasing','boundedness','convex'),
                            prior='Fast.LS', return.plot = T, tol = 0, sseed=12345, lower=0, upper = 1)
  lines(sort(xtr),f(sort(xtr)),type='l',lwd=2)
}




if (multi_monot_bd_syn == 'yes'){
  ###### true function 
  f <- function(x){
    1/(1+exp(-6*x))
  }
  ## data
  set.seed(12345)
  n <- 100
  sigN <- 0.2 # noise sd
  xtr <-  # lhsDesign(n=n,dimension=1,seed=12345)$design
    runif(n,0,1)
  ytr <- f(xtr) + rnorm(n,0,sd=sigN)
  N1 <- 9#10# # size of the 1st subdomain
  M <- 3#10# # nb of subdomains
  N <- N1*M # nb of knot points
  eta <- 200 # smooth approximate parameter
  nsim <- 5000 # nb of mcmc retained samples
  brn <- 1000 # burn in
  thin <- 1
  nu <- 2.5 # smoothness kernel parameter
  l <- round(l_est(nu,c(0,1),0.05),2) # length-scale
  tol <- 1e-9
  upper <- 1 # upper bound for boundedness constr
  
  ## LS-ESS with only monotonicity 
  post_mon_LS <- linCGP.ESS(y=ytr,x=xtr,N1=N1,M=M,nu=nu,l=l,eta=eta,nsim=nsim,burn.in=brn,thin=thin,
                            sig.in=sigN^2,tau.in=1,constrType='increasing',
                            prior='Fast.LS',return.plot=F,tol=tol,sseed=12345)
  ## LS-ESS with monotonicity and boundedness
  post_lin_LS <- linCGP.ESS(y=ytr,x=xtr,N1=N1,M=M,nu=nu,l=l,eta=eta,nsim=nsim,burn.in=brn,thin=thin,
                            sig.in=sigN^2,tau.in=1,constrType=c('increasing','boundedness'),
                            prior='Fast.LS',return.plot=F,tol=tol,sseed=12345,upper=upper)
  
  post_lin_HMC <- linCGP.HMC(y=ytr,x=xtr,N=N,nu=nu,l=l,nsim=nsim,burn.in=brn,thin=thin,
                             sig.in=sigN^2,tau.in=1,constrType=c('increasing','boundedness'),
                             return.plot=F,tol=tol,sseed=12345,upper=upper)
  ## computing some terms for the illustrations 
  f_low_mon <- post_mon_LS$f_low
  f_upp_mon <- post_mon_LS$f_upp
  f_low <- post_lin_LS$f_low
  f_upp <- post_lin_LS$f_upp
  f_low_HMC <- post_lin_HMC$f_low
  f_upp_HMC <- post_lin_HMC$f_upp
  
  fmean <- post_lin_LS$fmean
  MAP <- post_lin_LS$MAP
  fmean_HMC <- post_lin_HMC$fmean
  MAP_HMC <- post_lin_HMC$MAP
  fmean_mon <- post_mon_LS$fmean
  MAP_mon <- post_mon_LS$MAP
  ub <- max(f_upp,fmean,MAP,fmean_mon,MAP_mon,f_upp_mon,
            f_upp_HMC,fmean_HMC,MAP_HMC)
  lb <- min(f_low,fmean,MAP,fmean_mon,MAP_mon,f_low_mon,
            f_low_HMC,fmean_HMC,MAP_HMC)
  ytr <- ytr[order(xtr)]
  xtr <- sort(xtr)
  
  par(mfrow=c(1,1))
  par(mar=c(2.1,2.1,1.7,1.1)) # adapt margins
  ## Illustration only monotonicity
  plot(xtr,ytr,pch='*',lwd=2,lty=1,col='black',
       ylim=range(ub,lb,ytr),xlab='',ylab='')
  polygon(c(xtr,rev(xtr)),y=c(f_low_mon, rev(f_upp_mon)),border=F,col='gray')
  lines(xtr,fmean_mon,type='l',lty=4,lwd=2,col='blue')
  lines(xtr,MAP_mon,type='l',lty=2,lwd=2,col='red')
  points(xtr,ytr,pch='*')
  mtext(text =  'LS-ESS approach', side = 3, line = 0.9, cex = 0.8)
  mtext(text =  'only monotonicity constraints', side = 3, line = 0.1, cex = 0.8)
  x <- seq(0,1,length=100)
  lines(x,f(x),type='l',lwd=2)
  abline(h=1,lty=2,lwd=2)
  legend(0.45,0.8,
         c("true function","MAP","mAP"),
         col = c("black",'red','blue'), 
         text.col = "black",lty=c(1,2,4),
         lwd = c(2,2,2), text.font=1,box.lty=0,cex=0.8,
         bg='transparent')
  ## Illustration monotonicity and boudedness (LS-ESS approach)
  plot(xtr,ytr,pch='*',lwd=2,lty=1,col='black',
       ylim=range(ub,lb,ytr),xlab='',ylab='')
  polygon(c(xtr,rev(xtr)),y=c(f_low, rev(f_upp)),border=F,col='gray')
  lines(xtr,fmean,type='l',lty=4,lwd=2,col='blue')
  lines(xtr,MAP,type='l',lty=2,lwd=2,col='red')
  points(xtr,ytr,pch='*')
  mtext(text =  'LS-ESS approach', side = 3, line = 0.9, cex = 0.8)
  mtext(text =  'monotonicity and boundedness constraints', side = 3, line = 0.1, cex = 0.8)
  x <- seq(0,1,length=100)
  lines(x,f(x),type='l',lwd=2)
  abline(h=1,lty=2,lwd=2)
  legend(0.45,0.8,
         c("true function","MAP","mAP"),
         col = c("black",'red','blue'), 
         text.col = "black",lty=c(1,2,4),
         lwd = c(2,2,2), text.font=1,box.lty=0,cex=0.8,
         bg='transparent')
  
  ## Illustration monotonicity and boudedness (HMC sampler)
  plot(xtr,ytr,pch='*',lwd=2,lty=1,col='black',
       ylim=range(ub,lb,ytr),xlab='',ylab='')
  polygon(c(xtr,rev(xtr)),y=c(f_low_HMC, rev(f_upp_HMC)),border=F,col='gray')
  lines(xtr,fmean_HMC,type='l',lty=4,lwd=2,col='blue')
  lines(xtr,MAP_HMC,type='l',lty=2,lwd=2,col='red')
  points(xtr,ytr,pch='*')
  mtext(text =  'HMC approach', side = 3, line = 0.9, cex = 0.8)
  mtext(text =  'monotonicity and boundedness constraints', side = 3, line = 0.1, cex = 0.8)
  x <- seq(0,1,length=100)
  lines(x,f(x),type='l',lwd=2)
  abline(h=1,lty=2,lwd=2)
  legend(0.45,0.8,
         c("true function","MAP","mAP"),
         col = c("black",'red','blue'), 
         text.col = "black",lty=c(1,2,4),
         lwd = c(2,2,2), text.font=1,box.lty=0,cex=0.8,
         bg='transparent')
}
####################################################






if (mult_lin_syn == 'yes'){
  ###### true function 
  f <- function(x){
    5.6*sqrt(x)+10
  }
  ## lower and upper functions
  f_l <- function(x) 4*x^2+10
  f_u <- function(x) 4*x+12
  ## data
  set.seed(12345)
  n <- 50
  xtr <- runif(n,0,1)
  ytr <- f(xtr) + rnorm(n,0,sd=0.5)
  N1 <- 9#10 # size of the 1st subdomain
  M <- 3#5# # nb of subdomains
  N <- N1*M # nb of knot points
  eta <- 1000 # smooth approximate parameter
  nsim <- 5000 # nb of mcmc retained samples
  brn <- 1000 # burn in
  thin <- 1
  nu <- 2.5 # smoothness kernel parameter
  l <- round(l_est(nu,c(0,1),0.05),2) # length-scale
  sigN <- sqrt(0.5) # noise sd
  delta <- 1/(N-1)
  my_knots <- seq(0,1,by=delta)
  lower <- f_l(my_knots)
  upper <- f_u(my_knots)
  ## LS-ESS with only monotonicity 
  post_mon_LS <- linCGP.ESS(y=ytr,x=xtr,N1=N1,M=M,nu=nu,l=l,eta=eta,nsim=nsim,burn.in=brn,thin=thin,
                            sig.in=sigN^2,tau.in=1,constrType='increasing',
                            prior='Fast.LS',return.plot=F,tol=0,sseed=12345)
  ## LS-ESS with monotonicity and boundedness
  post_lin_LS <- linCGP.ESS(y=ytr,x=xtr,N1=N1,M=M,nu=nu,l=l,eta=eta,nsim=nsim,burn.in=brn,thin=thin,
                            sig.in=sigN^2,tau.in=1,constrType=c('increasing','boundedness'),
                            lower=lower,upper=upper,prior='Fast.LS',return.plot=F,tol=0,sseed=12345)
  
  ## computing some terms for illustration
  f_low_mon <- post_mon_LS$f_low
  f_upp_mon <- post_mon_LS$f_upp
  f_low <- post_lin_LS$f_low
  f_upp <- post_lin_LS$f_upp
  fmean <- post_lin_LS$fmean
  MAP <- post_lin_LS$MAP
  fmean_mon <- post_mon_LS$fmean
  MAP_mon <- post_mon_LS$MAP
  ub <- max(f_upp,fmean,MAP,fmean_mon,MAP_mon,f_upp_mon,upper)
  lb <- min(f_low,fmean,MAP,fmean_mon,MAP_mon,f_low_mon,lower)
  ytr <- ytr[order(xtr)]
  xtr <- sort(xtr)
  
  par(mfrow=c(1,1))
  par(mar=c(2.1,2.1,2.1,1.1)) # adapt margins
  ## Only monotonicity
  plot(xtr,ytr,pch='*',lwd=2,lty=1,col='black',
       ylim=range(ub,lb,ytr),xlab='',ylab='')
  polygon(c(xtr,rev(xtr)),y=c(f_low_mon, rev(f_upp_mon)),border=F,col='gray')
  lines(xtr,fmean_mon,type='l',lty=4,lwd=2,col='blue')
  lines(xtr,MAP_mon,type='l',lty=2,lwd=2,col='red')
  points(xtr,ytr,pch='*')
  mtext(text =  'LS-ESS approach', side = 3, line = 0.9, cex = 0.8)
  mtext(text =  'only monotonicity constraints', side = 3, line = 0.1, cex = 0.8)
  x <- seq(0,1,length=100)
  lines(x,f(x),type='l',lwd=2)
  lines(x,f_l(x),lty=2,lwd=2)
  lines(x,f_u(x),lty=2,lwd=2)
  legend(0,16.5,
         c("true function","MAP","mAP"),
         col = c("black",'red','blue'), 
         text.col = "black",lty=c(1,2,4),
         lwd = c(2,2,2), text.font=1,box.lty=0,cex=0.8,
         bg='transparent')
  ## Monotonicity and boudedness
  plot(xtr,ytr,pch='*',lwd=2,lty=1,col='black',
       ylim=range(ub,lb,ytr),xlab='',ylab='')
  polygon(c(xtr,rev(xtr)),y=c(f_low, rev(f_upp)),border=F,col='gray')
  lines(xtr,fmean,type='l',lty=4,lwd=2,col='blue')
  lines(xtr,MAP,type='l',lty=2,lwd=2,col='red')
  points(xtr,ytr,pch='*')
  mtext(text =  'LS-ESS approach', side = 3, line = 0.9, cex = 0.8)
  mtext(text =  'monotonicity and boundedness constraints', side = 3, line = 0.1, cex = 0.8)
  x <- seq(0,1,length=100)
  lines(x,f(x),type='l',lwd=2)
  lines(x,f_l(x),lty=2,lwd=2)
  lines(x,f_u(x),lty=2,lwd=2)
  legend(0,16.5,
         c("true function","MAP","mAP"),
         col = c("black",'red','blue'), 
         text.col = "black",lty=c(1,2,4),
         lwd = c(2,2,2), text.font=1,box.lty=0,cex=0.8,
         bg='transparent')
  
}
####################################################



if (mon_syn == 'yes'){
  ###### true monotone functions 
  f <- function(x,idx){
    if (idx == 1)
      5*x^2
    else if (idx == 2)
      3/(1+exp(-10*x+2.1))
    else if (idx == 3)
      ifelse(x>=0.6 & x<=1,(5*x-3)^3,0) # flat monotone fct
    else if (idx == 4){
      sum <- 0
      for (l in 1 : 100){
        sum <- sum+(l^(-1.7)*sin(l)*cos(pi*(l-0.5)*(1-x)))
      }
      return(sqrt(2)*sum)
    }
    else if (idx == 5)
      0.32*(10*x + sin(10*x)) # sinusoidal function
    else if (idx == 6)
      -2/(1 + exp(-12*x + 3))
  }
  ## synthetic data
  idx <- 6 # choice of the monotone function
  set.seed(12345)
  n <- 100
  xtr <- #lhsDesign(n=n,dimension=1,seed=12345)$design
    runif(n, 0, 1)
  ytr <- sapply(xtr, function(x) f(x,idx)) + rnorm(n,0,sd=0.5)
  N1 <- 10# # size of the 1st subdomain
  M <- 5# # nb of subdomains
  N <- N1*M # nb of knot points
  eta <- 100 # smooth approximate parameter
  nsim <- 5000#5000 # nb of mcmc retained samples
  brn <- 1000#1000 # burn in
  thin <- 1
  nu <- 2.5 # smoothness kernel parameter
  l <- round(l_est(nu,c(0,1),0.05),2) # length-scale
  sigN <- sqrt(0.5) # noise sd
  tol <- 1e-11
  ## Illustration LS-ESS
  post_mon_LS = linCGP.ESS(y=ytr,x=xtr,N1=N1,M=M,nu=nu,l=l,eta=eta,nsim=nsim,burn.in=brn,thin=thin,
                           sig.in=sigN^2,tau.in=1,constrType='decreasing',
                           prior='Fast.LS',return.plot=T,tol=0,sseed=12345)
  
  mtext(text = 'LS-ESS with monotonicity', side = 3, line = 0.8, cex = 0.8)
  mtext(text =  paste("Running Time (s) = ", round(post_mon_LS$time[3],2), "and dimension N = ",N), side = 3, line = 0.1, cex = 0.8)
  x <- seq(0,1,length=100)
  lines(x,f(x,idx=idx),type='l',lwd=2)
  ## Errors LS-ESS approach
  xtest <- seq(0,1,length=50)
  ytest <- f(x=xtest,idx=idx)
  y_pred_MAP <- fcth(xtest,u=post_mon_LS$knots,N=N)%*%post_mon_LS$z_star
  y_pred_mAP <- fcth(xtest,u=post_mon_LS$knots,N=N)%*%post_mon_LS$z_mean
  errors <- data.frame(MAP=err.meas.report(ytest=ytest,ytr=ytr,mu=y_pred_MAP,out=post_mon_LS,type="all"),
                       mAP=err.meas.report(ytest=ytest,ytr=ytr,mu=y_pred_mAP,out=post_mon_LS,type="all"))
  
  ## Illustration HMC sampler
  post_mon_HMC <- linCGP.HMC(y=ytr,x=xtr,N=N,nu=nu,l=l,nsim=nsim,burn.in=brn,thin=thin,
                             sig.in=sigN^2,tau.in=1,constrType='decreasing',
                             return.plot=T,sseed=12345,tol=tol)
  mtext(text = 'HMC with monotonicity', side = 3, line = 0.8, cex = 0.8)
  mtext(text =  paste("Running Time (s) = ", round(post_mon_HMC$time[3],2), "and dimension N = ",N), side = 3, line = 0.1, cex = 0.8)
  x <- seq(0,1,length=100)
  lines(x,f(x,idx=idx),type='l',lwd=2)
  ## Errors HMC sampler
  y_pred_MAP_HMC <- fcth(xtest,u=post_mon_HMC$knots,N=N)%*%post_mon_HMC$z_star
  y_pred_mAP_HMC <- fcth(xtest,u=post_mon_HMC$knots,N=N)%*%post_mon_HMC$z_mean
  errors_HMC <- data.frame(MAP=err.meas.report(ytest=ytest,ytr=ytr,mu=y_pred_MAP_HMC,out=post_mon_HMC,type="all"),
                           mAP=err.meas.report(ytest=ytest,ytr=ytr,mu=y_pred_mAP_HMC,out=post_mon_HMC,type="all"))
  
}
####################################################








####################################################
############ Boundedness constraints ###############
####################################################


if (bound_syn == 'yes'){
  ###### bounded synthetic data ##########
  f <- function(x){
    # cos(10*x)
    ## fct used for generalization of KW correspondence with noise
    # ifelse(x>=0 & x<= 2/3, cos(2*pi*x+pi/3),0.5) # [-1;0.5]
    ## Andrew Pensoneault (Nonnegativity enforced GPR 2020)
    1/(1+(10*x)^4)+0.5*exp(-100*(x-0.5)^2) # [0;1]
  }
  ## data
  set.seed(12345)
  sigN <- 0.1#0.3# # sd noise
  n <- 100 # nb of training data
  xtr <- lhsDesign(n=n,dimension=1,seed=12345)$design
  # runif(n,0,1)
  ytr <- f(xtr)+rnorm(n,0,sd=sigN)
  lower <- 0#-1#    # lower bound
  upper <- 1#0.5#   # upper bound 
  N1 <- 10#9#   # size of the 1st subdomain
  M <- 10#3#     # nb of subdomains
  N <- N1*M
  eta <- 100 # pdf approximation parameter
  nsim <- 5000 # nb of simulation
  brn <- 1000 # burn in 
  thin <- 1
  nu <- 1.5 # smoothness kernel parameter
  l <- round(l_est(nu,c(0,1),0.05),2) # length-scale
  tol <- 1e-11
  ## LS-ESS approach
  post_bound_LS <- linCGP.ESS(y=ytr,x=xtr,N1=N1,M=M,nu=nu,l=l,eta=eta,nsim=nsim,burn.in=brn,thin=thin,
                              sig.in=sigN^2,tau.in=1,constrType='boundedness',lower=lower,upper=upper,
                              prior='Fast.LS',return.plot=T,tol=tol,sseed=12345)
  abline(h=c(lower,upper),lty=2,lwd=2)
  x <- seq(0,1,length=100)
  lines(x,f(x),type='l',lwd=2)
  mtext(text =  'Large-scale ESS', side = 3, line = 0.8, cex = 0.8)
  mtext(text =  paste("Running Time (s) = ", round(post_bound_LS$time[3],2), "and dimension N = ",N), side = 3, line = 0.1, cex = 0.8)
  legend(0.55,1,
         c("true function","MAP","mAP"),
         col = c('black','red','blue'), 
         text.col = "black",lty=c(1,2,4),
         lwd = c(2,2,2), text.font=1,box.lty=0,cex=0.8,
         bg='transparent')
  ## Illustration FFT-WC
  post_bound_WC <- linCGP.WC.ESS(y=ytr,x=xtr,N=N,nu=nu,l=l,eta=eta,nsim=nsim,burn.in=brn,thin=thin,
                                 sig.in=sigN^2,tau.in=1,constrType='boundedness',
                                 return.plot=T,sseed=12345,lower=lower,upper=upper)
  abline(h=c(upper,lower),lty=2,lwd=2)
  x <- seq(0,1,length=100)
  lines(x,f(x),type='l',lwd=2)
  mtext(text =  'FFT-WC ESS', side = 3, line = 0.8, cex = 0.8)
  mtext(text =  paste("Running Time (s) = ", round(post_bound_WC$time[3],2), "and dimension N = ",N), side = 3, line = 0.1, cex = 0.8)
  legend(0.55,1,
         c("true function","MAP","mAP"),
         col = c('black','red','blue'), 
         text.col = "black",lty=c(1,2,4),
         lwd = c(2,2,2), text.font=1,box.lty=0,cex=0.8,
         bg='transparent')
  
  ## Illustration HMC sampler
  post_bound_HMC <- linCGP.HMC(y=ytr,x=xtr,N=N,nu=nu,l=l,nsim=nsim,burn.in=brn,thin=thin,
                               sig.in=sigN^2,tau.in=1,lower=lower,upper=upper,constrType='boundedness',
                               return.plot=T,sseed=12345,tol=tol)
  x <- seq(0,1,length=100)
  lines(x,f(x),type='l',lwd=2)
  abline(h=c(lower,upper),lty=2,lwd=2)
  mtext(text =  'Hamiltonian Monte Carlo', side = 3, line = 0.8, cex = 0.8)
  mtext(text =  paste("Running Time (s) = ", round(post_bound_HMC$time[3],2), "and dimension N = ",N), side = 3, line = 0.1, cex = 0.8)
  legend(0.55,1,
         c("true function","MAP","mAP"),
         col = c('black','red','blue'), 
         text.col = "black",lty=c(1,2,4),
         lwd = c(2,2,2), text.font=1,box.lty=0,cex=0.8,
         bg='transparent')
}
####################################################


####################################################
#### comparison of behavior of HMC and LS-ESS ######
####################################################
######## comparison HMC vs LS-ESS ##################
if (comp_mon_HMC_LS_WC == 'yes'){
  ### monotone functions
  f <- function(x, idx){
    if (idx == 1)
      5*x^2
    else if (idx == 2)
      # 3/(1+exp(-10*x+2.1))
      - 2/(1 + exp(-12*x + 3))
    else if (idx == 3)
      ifelse(x>=0.6 & x<=1,(5*x-3)^3,0) # flat monotone fct
    else if (idx == 4){
      sum <- 0
      for (l in 1 : 100){
        sum <- sum+(l^(-1.7)*sin(l)*cos(pi*(l-0.5)*(1-x)))
      }
      return(sqrt(2)*sum)
    }
  }
  ## synthetic data
  idx <- 2 # choice of a monotone function
  set.seed(12345)
  n <- 100
  xtr <- lhsDesign(n=n,dimension=1,seed=12345)$design
  # runif(n,0,1)
  y.true <- sapply(xtr,FUN=function(x) f(x,idx=idx)) 
  # sigN <- 0.2 * diff(range(y.true))
  sigN <- 0.5 # noise sd
  ytr <- y.true + rnorm(n,0,sd=sigN)
  
  N1 <- 9#10# # size of the 1st subdomain
  M <- 3#10# # nb of subdomains
  N <- N1*M # nb of knot points
  eta <- 100 # smooth approximate parameter
  nsim <- 5000 # nb of mcmc retained samples
  brn <- 1000 # burn in
  thin <- 1
  nu <- 1.5 # smoothness kernel parameter
  l <- round(l_est(nu,c(0,1),0.05),2) # length-scale
  tol <- 1e-9
  ## Illustration LS-ESS
  post_mon_LS <- linCGP.ESS(y=ytr,x=xtr,N1=N1,M=M,nu=nu,l=l,eta=eta,nsim=nsim,burn.in=brn,thin=thin,
                            sig.in=sigN^2,tau.in=1,constrType='decreasing',
                            prior='Fast.LS',return.plot=T,tol=0,sseed=12345)
  x <- seq(0,1,length=100)
  lines(x,f(x,idx=idx),type='l',lwd=2)
  mtext(text =  'Large-scale ESS', side = 3, line = 0.8, cex = 0.8)
  mtext(text =  paste("Running Time (s) = ", round(post_mon_LS$time[3],2), "and dimension N =",N), side = 3, line = 0.1, cex = 0.8)
  # mtext(text =  paste("RMSE mAP = ", round(rmse_mAP_LS,2), "and RMSE MAP = ", round(rmse_MAP_LS,2)), side = 3, line = 0.1, cex = 0.8)
  legend(0.5,1,
         c("true function","MAP","mAP"),
         col = c('black','red','blue'), 
         text.col = "black",lty=c(1,2,4),
         lwd = c(2,2,2), text.font=1,box.lty=0,cex=0.8,
         bg='transparent')
  ## Illustration WC.ESS
  post_mon_WC <- linCGP.WC.ESS(y=ytr,x=xtr,N=N,nu=nu,l=l,eta=eta,nsim=nsim,burn.in=brn,thin=thin,
                               sig.in=sigN^2,tau.in=1,constrType='decreasing',
                               return.plot=T,sseed=12345,tol=0)
  x <- seq(0,1,length=100)
  lines(x,f(x,idx=idx),type='l',lwd=2)
  mtext(text =  'FFT-WC ESS', side = 3, line = 0.8, cex = 0.8)
  mtext(text =  paste("Running Time (s) = ", round(post_mon_WC$time[3],2), "and dimension N = ",N), side = 3, line = 0.1, cex = 0.8)
  legend(0.5,1,
         c("true function","MAP","mAP"),
         col = c('black','red','blue'), 
         text.col = "black",lty=c(1,2,4),
         lwd = c(2,2,2), text.font=1,box.lty=0,cex=0.8,
         bg='transparent')
  
  ## Illustration HMC
  post_mon_HMC <- linCGP.HMC(y=ytr,x=xtr,N=N,nu=nu,l=l,nsim=nsim,
                             burn.in=brn,thin=thin,sig.in=sigN^2,tau.in=1,
                             return.plot=T,sseed=12345,tol=tol,constrType='decreasing')
  x <- seq(0,1,length=100)
  lines(x,f(x,idx=idx),type='l',lwd=2)
  mtext(text =  'Hamiltonian Monte Carlo', side = 3, line = 0.8, cex = 0.8)
  mtext(text =  paste("Running Time (s) = ", round(post_mon_HMC$time[3],2), "and dimension N = ",N), side = 3, line = 0.1, cex = 0.8)
  legend(0.5,1,
         c("true function","MAP","mAP"),
         col = c('black','red','blue'), 
         text.col = "black",lty=c(1,2,4),
         lwd = c(2,2,2), text.font=1,box.lty=0,cex=0.8,
         bg='transparent')
} 
##################################################





if (runtime_mon_HMC_LS == 'yes'){
  n <- seq(500,5000,by=500) # sample size
  N1 <- 25#9#50
  trial <- 10#
  sigN <- 0.5 # sd noise
  nsim <- 1 # mcmc iterations
  brn <- 0#1000 # burn in
  tol <- 1e-5#1e-11
  nu <- 1.5#2.5
  eta <- 100
  l <- l_est(nu,c(0,1),0.1) # length-scale
  ### monotone functions
  f <- function(x, idx){
    if (idx == 1)
      5*x^2
    else if (idx == 2)
      # 3/(1+exp(-10*x+2.1))
      -2/(1+exp(-12*x+3))
    else if (idx == 3)
      ifelse(x>=0.6 & x<=1,(5*x-3)^3,0) # flat monotone fct
    else if (idx == 4){
      sum <- 0
      for (l in 1 : 100){
        sum <- sum+(l^(-1.7)*sin(l)*cos(pi*(l-0.5)*(1-x)))
      }
      return(sqrt(2)*sum)
    }
  }
  
  timeLS <- matrix(NA,trial,length(n))
  timeHMC <- matrix(NA,trial,length(n))
  for (i in 1 : length(n)){
    print(i)
    M <- ceiling(n[i]/200)
    N <- N1*M
    for (Q in 1 : trial){
      
      ## split samples
      set.seed(123)
      idx <- 2 # choice of a monotone function
      x <- lhsDesign(n=n[i],dimension=1,seed=12345)$design
      # runif(n[i],0,1)
      ytrue <- sapply(x,FUN=function(x) f(x,idx=idx))
      y <- ytrue + rnorm(n[i], 0, sd=sigN)
      
      post_mon_syn_LS <- linCGP.ESS(y,x,N1=N1,M=M,nu=nu,l=l,eta=eta,nsim=nsim,
                                    burn.in=brn,thin=1,sig.in=sigN^2,tau.in=1,constrType='decreasing',
                                    prior='LS.KLE',return.plot=F,sseed=Q,tol=tol)
      post_mon_syn_HMC <- linCGP.HMC(y,x,N=N,nu=nu,l=l,nsim=nsim,constrType='decreasing',
                                     burn.in=brn,thin=1,sig.in=sigN^2,tau.in=1,
                                     sseed=Q,tol=tol,return.plot=F)
      timeLS[Q,i] <- post_mon_syn_LS$time[3]
      timeHMC[Q,i] <- post_mon_syn_HMC$time[3]
    }
  }
  ## Illustration
  par(mfrow=c(1,1))
  par(mar=c(3.1,3.1,2.1,1.1)) # adapt margins
  plot(n,colMeans(timeLS),type='l',lwd=2,
       ylim=range(colMeans(timeLS),colMeans(timeHMC)),
       ylab='',xlab='')
  lines(n,colMeans(timeHMC),type='l',lwd=2,lty=2)
  title(xlab='sample size (n)',ylab='Run-time/iteration (s)',line=2)
  mtext(text='monotonicity constraints',line=0.3)
  legend(1000,0.9*max(colMeans(timeHMC)),c('HMC','LS-ESS'),col='black',
         lwd=2,lty=c(2,1),text.font=1,box.lty=0,cex=1,
         bg='transparent')
  ##################################################
  ## plot zoom inset
  # par(mar=c(0.1,0.1,0.1,0.1)) # adapt margins
  # df = data.frame(n,colMeans(timeLS))
  # G <- ggplot(df,aes(x=n,y=colMeans(timeLS)),color=variable)+
  #   geom_line(aes(x=n, y=colMeans(timeLS),xlab=''), color="black")+
  #   geom_line(aes(x=n,y=colMeans(timeHMC)),col='black',lty=2)
  # 
  # G <- G  + theme_bw()
  # Glabs <- G+labs(title='monotonicity constraints',
  #                 x='sample size (n)',
  #                 y='Run-time/iteration (s)')+
  #   theme(plot.title = element_text(hjust = 0.5, size=10),
  #         plot.subtitle = element_text(hjust = 0.5,size=10))
  # Glabs
  # Gz <- G+coord_cartesian(xlim=c(500,1500),ylim = c(0,0.05),
  #                         expand = T)+
  #   theme(axis.text.x = element_text(size=8),
  #         axis.text.y = element_text(size=8),
  #         axis.title.x = element_text(size=0),
  #         axis.title.y = element_text(size=0))
  # Gzlabs <- Gz+labs(x="", y="")
  # v <- viewport(x = 0.4, y = 0.65,
  #               width = 0.45, height = 0.45)
  # print(Gzlabs, vp=v)
}
####################################################





if (runtime_bound_HMC_LS == 'yes'){
  n <- seq(500,5000,by=500) # sample size
  N1 <- 25#9#50
  trial <- 1#0 # average trial
  sigN <- 0.1#0.3# # sd noise
  nsim <- 1#4000 # mcmc iterations
  brn <- 0#1000 # burn in
  tol <- 1e-11
  nu <- 1.5
  eta <- 100
  l <- l_est(nu,c(0,1),0.05) # length-scale
  f <- function(x){
    # cos(10*x)
    ## fct used for generalization of KW correspondence with noise
    # ifelse(x>=0 & x<= 2/3, cos(2*pi*x+pi/3),0.5) # [-1;0.5]
    ## Andrew Pensoneault (Nonnegativity enforced GPR 2020)
    1/(1+(10*x)^4)+0.5*exp(-100*(x-0.5)^2) # [0;1]
  }
  a <- 0#-1# # lower bound
  b <- 1#0.5# # upper bound
  timeLS <- matrix(NA,trial,length(n))
  timeHMC <- matrix(NA,trial,length(n))
  for (i in 1 : length(n)){
    print(i)
    M <- ceiling(n[i]/200)
    N <- N1*M
    for (Q in 1 : trial){
      
      ## split samples
      set.seed(123)
      x <- runif(n[i],0,1)
      y <- f(x) + rnorm(n[i],0,sd=sigN)
      
      post_bound_syn_LS <- linCGP.ESS(y,x,N1=N1,M=M,nu=nu,l=l,eta=eta,nsim=nsim,
                                      burn.in=brn,thin=1,sig.in=sigN^2,tau.in=1,return.plot=F,
                                      prior='Fast.LS',sseed=Q,tol=0,lower=a,upper=b,
                                      constrType='boundedness')
      post_bound_syn_HMC <- linCGP.HMC(y,x,N=N,nu=nu,l=l,nsim=nsim,return.plot = F,
                                       burn.in=brn,thin=1,sig.in=sigN^2,tau.in=1,
                                       sseed=Q,tol=tol,lower=a,upper=b,constrType='boundedness')
      timeLS[Q,i] <- post_bound_syn_LS$time[3]
      timeHMC[Q,i] <- post_bound_syn_HMC$time[3]
    }
  }
  ## Illustration
  par(mfrow=c(1,1))
  par(mar=c(3.1,3.1,2.1,1.1)) # adapt margins
  plot(n,colMeans(timeLS),type='l',lwd=2,
       ylim=range(colMeans(timeLS),colMeans(timeHMC)),
       ylab='',xlab='')
  # lines(n,colMeans(timeWC),type='l',lwd=2,lty=5)
  lines(n,colMeans(timeHMC),type='l',lwd=2,lty=2)
  title(xlab='sample size (n)',ylab='Run-time/iteration (s)',line=2)
  mtext(text='boundedness constraints',line=0.3)
  legend(1000,0.9*max(colMeans(timeHMC)),c('HMC','LS-ESS'),col='black',
         lwd=2,lty=c(2,1),text.font=1,box.lty=0,cex=1,
         bg='transparent')
}










if (decreasing_syn == 'yes'){
  ###### decreasing example
  f <- function(x,idx){
    if (idx == 1)
      cos(x)
    else if (idx == 2)
      exp(-6*x)
  }
  
  ## synthetic data
  idx <- 2 # choice of the monotone nonincreasing function
  set.seed(12345)
  n <- 100
  xtr <- runif(n,0,1)
  ytr <- sapply(xtr,function(x) f(x,idx))+rnorm(n,0,sd=0.1)
  N1 <- 9# # size of the 1st subdomain
  M <- 3# # nb of subdomains
  N <- N1*M # nb of knot points
  eta <- 100 # smooth approximate parameter
  nsim <- 5000 # nb of mcmc retained samples
  brn <- 1000 # burn in
  thin <- 1
  nu <- 2.5 # smoothness kernel parameter
  l <- round(l_est(nu,c(0,1),0.05),2) # length-scale
  sigN <- sqrt(0.5) # noise sd
  tol <- 1e-11
  ## Illustration LS-ESS
  post_mon_LS = linCGP.ESS(y=ytr,x=xtr,N1=N1,M=M,nu=nu,l=l,eta=eta,nsim=nsim,burn.in=brn,thin=thin,
                           sig.in=sigN^2,tau.in=1,constrType='decreasing',
                           prior='Fast.LS',return.plot=T,tol=tol,sseed=12345)
  mtext(text = 'LS-ESS with monotonicity', side = 3, line = 0.8, cex = 0.8)
  mtext(text =  paste("Running Time (s) = ", round(post_mon_LS$time[3],2), "and dimension N = ",N), side = 3, line = 0.1, cex = 0.8)
  x <- seq(0,1,length=100)
  lines(x,f(x,idx=idx),type='l',lwd=2)
}
############################################


## end