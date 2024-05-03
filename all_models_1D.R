############################################################
#### Function for MCMC samples using ESS and LS algo #######
############################################################

source('all_base_functions_1D.R')
library(Matrix) # crossprod
library(nloptr) # for resolving optim problem
## to install the old version of 'tmg' package
# library(remotes)
# install_version("tmg", "0.3")
library(tmg)




####################################################
########## Linear constraints Axi+B>0 ##############
####################################################
## Function for drawing posterior samples using ESS and LS with fixed hyperparameters \nu and \ell:
## For increasing,decreasing and boundedness functions estimation 
linCGP.ESS <- function(y, x, N1, M, nu, l, est.l = F, eta, nsim, burn.in, thin, tau.in, sig.in, xi.in, lower = -Inf, upper = Inf,
                       constrType, prior, tau.fix, sig.fix, sseed, verbose, return.plot, tol) {
  # y:Response variable; x: vector to form design matrix X (n x N)
  # N1: number of knots first subdomain; M: nb of subdomain
  # nu:smoothness parameter of Matern; l:length-scale parameter of Matern
  # est.l logical if TRUE, the estimated value of length-scale para "l" will be employed 
  # eta: parameter of the approximation function of the indicator functions
  # nsim, burn.in, thin : nsim samples, burning and thinning for MCMC
  # tau.in, sig.in, xi.in : initial values (supplied by user or use the default values)
  # tau.fix,sig.fix : if fixed values of the parameters are to use
  # verbose : logical; if TRUE, prints current status; default is TRUE
  # return.plot : logical; if true a plot of estimate with 95% CI is returned; default is TRUE
  # tol: tolerance for numerical stability
  # prior : method for sampling prior (LS.KLE or Fast.LS or Cholesky decomposition from mvtnorm package); default 'Fast.LS'
  # constrType: type of inequality constraints ('increasing','decreasing','convex','concave','boundedness')
  # lower and upper: options for only boundedness constr
  
  # OUTPUT: Posterior samples on xi,tau,sig and fhat with posterior mean, 95% CI of fhat, mode of posterior distribution z_star
  
  if (length(y) != length(x))
    stop("y and x should be of same length!")
  y <- y[order(x)]
  x <- sort(x)
  n <- length(y)
  if (missing(M))
    M <- 1
  if (missing(N1))
    N1 <- 5
  N <- N1 * M
  if (N <= 1)
    stop('N = N1 x M should be greater than or equal to 2')
  delta <- 1/(N - 1)
  my_knots <- seq(0, 1, by = delta)
  X <- fcth(x, u = my_knots, N)
  
  if (missing(nu))
    stop("nu needs to be supplied")
  if (nu == 0)
    stop("nu cannot be zero")
  if (missing(constrType))
    stop('constrType should be specified and must be one or any logical combination of the following: \'increasing\',
         \'decreasing\', \'convex\',\'concave\',\'boundedness\'')
  if (!missing(l) && l <= 0)
    stop("l should be positive")
  if (missing(return.plot))
    return.plot <- TRUE
  if (!missing(sseed))
    set.seed(sseed)
  if (missing(sseed))
    set.seed(Sys.Date())
  if (missing(verbose))
    verbose <- TRUE
  if (missing(prior))
    prior <- 'Fast.LS'
  if (missing(tol))
    tol <- 1e-6
  if (missing(eta))
    eta <- 200
  if (missing(nsim))
    nsim <- 5000
  if (missing(burn.in))
    burn.in <- 1000
  if (missing(thin))
    thin <- 1
  if (!missing(tau.fix))
    tau.in <- tau.fix
  if (!missing(sig.fix))
    sig.in <- sig.fix
  if (missing(tau.fix) && missing(tau.in))
    tau.in <- 1
  if (missing(sig.fix) && missing(sig.in))
    sig.in <- 0.1 * diff(range(y))
  
  ## Matrix and vector of linear constr
  A_b <- c() # for boundedness constr
  B_b <- c() # for boundedness constr
  A <- c() # for other constr
  B <- c() # for other constr
  if (any(constrType == 'boundedness')) {
    if (missing(lower) & missing(upper))
      stop('Error: lower or upper should be provided')
    if (any(lower >= upper)) 
      stop("The elements from \"upper\" has to be greater than the elements from \"lower\"")
    if (length(lower) == 1)
      lower <- rep(lower, N)
    if (length(upper) == 1)
      upper <- rep(upper, N)
    A_b <- rbind(A,constrSys(N = N, type = 'boundedness', lower = lower, upper = upper)$A)
    B_b <- c(B,constrSys(N = N, type = 'boundedness', lower = lower, upper = upper)$B)
  }
  
  if (any(constrType == 'increasing')) {
    if (any(constrType == 'decreasing'))
      stop('Error: \'increasing\' and \'decreasing\' are not compabitle together')
    A <- rbind(A,constrSys(N = N, type = 'increasing')$A)
    B <- c(B,constrSys(N = N, type = 'increasing')$B)
  }
  if (any(constrType == 'decreasing')) {
    if (any(constrType == 'increasing'))
      stop('Error: \'increasing\' and \'decreasing\' are not compabitle together')
    A <- rbind(A,constrSys(N = N, type = 'decreasing')$A)
    B <- c(B,constrSys(N = N, type = 'decreasing')$B)
  }
  if (any(constrType == 'convex')) {
    if (any(constrType == 'concave'))
      stop('Error: \'convex\' and \'concave\' are not compabitle together')
    A <- rbind(A,constrSys(N = N, type = 'convex')$A)
    B <- c(B,constrSys(N = N, type = 'convex')$B)
  }
  if (any(constrType == 'concave')) {
    if (any(constrType == 'convex'))
      stop('Error: \'convex\' and \'concave\' are not compabitle together')
    A <- rbind(A,constrSys(N = N, type = 'concave')$A)
    B <- c(B,constrSys(N = N, type = 'concave')$B)
  }
  
  ## for the estimation of the length-scale para
  fzetoil <- function(l) {
    K_inv <- tinv(covmat(knot = my_knots, nu = nu, l = l))
    XXK <- crossprod(X) + K_inv
    solve.QP(XXK, dvec = as.vector(t(X)%*%y),
             Amat = t(rbind(A_b, A)), bvec = -c(B_b, B), meq = 0)$solution
  }
  fMAP <- function(l) {
    X %*% fzetoil(l)
  }
  MSPE <- function(l) {
    return(mean((y - fMAP(l))^2))
  }
  if (missing(l)) {
    if (est.l == T) {
      opts <- list("algorithm"="NLOPT_LN_COBYLA",
                   "xtol_rel" = 1e-5)
      l <- nloptr(l_est(nu, range = range(my_knots), 0.05), MSPE, lb = 0.1, ub = 1, opts = opts)$solution
    }
    else if (est.l == F)
      l <- l_est(nu, range = range(my_knots), 0.05)
  }
  if (missing(xi.in)) {
    xi.in <- fzetoil(l)
  }
  
  ## Inverse of prior cov matrix (Toeplitz)
  K_inv <- tinv(covmat(knot = my_knots, nu = nu, l = l))
  
  tau <- tau.in
  sig <- sig.in
  xi_in <- xi.in
  
  em <- nsim + burn.in
  ef <- nsim/thin
  xi_sam <- matrix(NA, nrow = N, ncol = ef)
  tau_sam <- rep(NA, ef)
  sig_sam <- rep(NA, ef)
  fhat_sam <- matrix(NA, nrow = n, ncol = ef)
  
  if (verbose)
    print("MCMC sample draws:")
  
  ptm <- proc.time()
  for (i in 1 : em) {
    # sampling Xi:
    if (prior == 'Fast.LS') {
      nu.ess <- Fast.LS(u = my_knots, M = M, N1 = N1, nu = nu, l = l, tausq = tau, tol = tol, sseedLS = i)
    }
    else if (prior == 'LS.KLE') {
      nu.ess <- LS.KLE(u = my_knots, M = M, N1 = N1, nu = nu, l = l, tausq = tau, tol = tol, sseedLS = i)
    }
    else if (prior == 'chol') {
      nu.ess <- chol_prior(u = my_knots, N = N, nu = nu, l = l, tausq = tau, tol = tol, seedchol = i)
    }
    
    xi_out <- ESS.linear(y = y, X = X, beta = xi_in, nu_ess = nu.ess, sigsq = sig, eta = eta, A = A, B = B, lower = lower, upper = upper, constrType = constrType, seeds = i)
    
    if (length(constrType) == 1) {
      if (constrType == 'increasing') {
        xi_out <- increasing_vector(xi_out)
      }
      else if (constrType == 'decreasing') {
        xi_out <- decreasing_vector(xi_out)
      }
      else if (constrType == 'convex') {
        xi_out <- convex_vector(xi_out)
      }
      else if (constrType == 'concave') {
        xi_out <- concave_vector(xi_out)
      }
      else if (constrType == 'boundedness') {
        xi_out <- pmin(pmax(lower,xi_out),upper)
      }
    }
    else if (length(constrType) == 2) {
      if (any(constrType == 'boundedness') & any(constrType == 'increasing')) {
        xi_out <- increasing_vector(xi_out)
        xi_out <- pmin(pmax(lower, xi_out), upper)
      }
      else if (any(constrType == 'boundedness') & any(constrType == 'decreasing')) {
        xi_out <- decreasing_vector(xi_out)
        xi_out <- pmin(pmax(lower, xi_out), upper)
      }
      else if (any(constrType == 'boundedness') & any(constrType == 'convex')) {
        xi_out <- convex_vector(xi_out)
        xi_out <- pmin(pmax(lower, xi_out), upper)
      }
      else if (any(constrType == 'boundedness') & any(constrType == 'concave')) {
        xi_out <- concave_vector(xi_out)
        xi_out <- pmin(pmax(lower, xi_out), upper)
      }
      else if (any(constrType == 'increasing') & any(constrType == 'convex'))
        xi_out <- inc_conv_vector(xi_out)
      else if (any(constrType == 'increasing') & any(constrType == 'concave'))
        xi_out <- inc_conc_vector(xi_out)
      else if (any(constrType == 'decreasing') & any(constrType == 'convex'))
        xi_out <- dec_conv_vector(xi_out)
      else if (any(constrType == 'decreasing') & any(constrType == 'concave'))
        xi_out <- dec_conc_vector(xi_out)
    }
    else if (length(constrType == 3)) {
      if (any(constrType == 'boundedness') & any(constrType == 'increasing') & any(constrType == 'convex')) {
        xi_out <- inc_conv_vector(xi_out)
        xi_out <- pmin(pmax(lower, xi_out), upper)
      }
      else if (any(constrType == 'boundedness') & any(constrType == 'increasing') & any(constrType == 'concave')) {
        xi_out <- inc_conc_vector(xi_out)
        xi_out <- pmin(pmax(lower, xi_out), upper)
      }
      else if (any(constrType == 'boundedness') & any(constrType == 'decreasing') & any(constrType == 'convex')) {
        xi_out <- dec_conv_vector(xi_out)
        xi_out <- pmin(pmax(lower, xi_out), upper)
      }
      else if (any(constrType == 'boundedness') & any(constrType == 'decreasing') & any(constrType == 'concave')) {
        xi_out <- dec_conc_vector(xi_out)
        xi_out <- pmin(pmax(lower, xi_out), upper)
      }
    }
    set.seed(2 * i)
    # sampling \sigma^2:
    Xxi <- as.vector(X %*% xi_out)
    if (missing(sig.fix)) {
      y_star <- y - Xxi
      sig <- 1/rgamma(1, shape = n/2, rate = sum(y_star^2)/2)
    }
    # sampling \tau^2:
    if (missing(tau.fix))
      tau <- 1/rgamma(1, shape = N/2, rate = (t(xi_out) %*% K_inv %*% xi_out)/2)
    
    # storing MCMC samples:
    if (i > burn.in && i%%thin == 0) {
      xi_sam[,(i-burn.in)/thin] <- xi_out
      sig_sam[(i-burn.in)/thin] <- sig
      tau_sam[(i-burn.in)/thin] <- tau
      fhat_sam[,(i-burn.in)/thin] <- Xxi
    }
    
    if (i %% 1000 == 0 && verbose) {
      print(i)
    }
    
    # renewing the intial value:
    xi_in <- xi_out
  }
  tm <- proc.time()-ptm
  
  ## posterior Mode
  XXK <- crossprod(X)/mean(sig_sam) + K_inv/mean(tau_sam)
  z_star <- solve.QP(XXK, dvec = as.vector(t(X)%*%y)/mean(sig_sam),
                     Amat=t(rbind(A_b, A)), bvec = -c(B_b, B), meq = 0)$solution
  MAP <- X %*% z_star 
  ## mAP estimate
  z_mean <- rowMeans(xi_sam)
  fmean <- rowMeans(fhat_sam) 
  qnt <- apply(fhat_sam, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
  f_low <- qnt[1, ]
  f_upp <- qnt[2, ]
  ub <- max(f_low, f_upp, fmean,MAP)
  lb <- min(f_low, f_upp, fmean,MAP)
  
  if (return.plot) {
    par(mfrow=c(1, 1))
    par(mar=c(2.1, 2.1, 2.1, 1.1)) # adapt margins
    plot(x, y, pch = '*', lwd = 2, lty = 1, col = 'black',
         ylim = range(ub, lb, y), xlab = '', ylab = '')
    polygon(c(x, rev(x)), y = c(f_low, rev(f_upp)), border = F, col = 'gray')
    lines(x, fmean, type = 'l', lty = 4, lwd = 2, col = 'blue')
    lines(x, MAP, type = 'l', lty = 2, lwd = 2, col = 'red')
    points(x, y, pch = '*')
  }
  
  return(list("time" = tm, "xi_sam" = xi_sam, "sig_sam" = sig_sam, "tau_sam" = tau_sam,
              "fhat_sam" = fhat_sam,"fmean" = fmean,
              "f_low" = f_low, "f_upp" = f_upp, "z_star" = z_star,
              "MAP" = MAP, "knots" = my_knots, "z_mean" = z_mean, "length-scale" = l))
}
###################################################






### Function for drawing posterior samples using ESS and FFT-WC with fixed hyperparameters \nu and \ell:
### For increasing,decreasing and boundedness functions estimation 
linCGP.WC.ESS <- function(y,x,N,nu,l,est.l=F,eta,nsim,burn.in,thin,tau.in,sig.in,xi.in,lower=-Inf,upper=Inf,
                          constrType,tau.fix,sig.fix,sseed,verbose,return.plot,tol) {
  # y:Response variable; x: vector to form design matrix X (n x N)
  # N: number of knots 
  # nu:smoothness parameter of Matern; l:length-scale parameter of Matern
  # est.l logical if TRUE, the estimated value of length-scale para "l" will be employed 
  # eta:parameter of the approximation function of the indicator functions
  # nsim, burn.in, thin : nsim samples, burning and thinning for MCMC
  # tau.in, sig.in, xi.in : initial values (supplied by user or use the default values)
  # tau.fix,sig.fix : if fixed values of the parameters are to use
  # verbose : logical; if TRUE, prints current status; default is TRUE
  # return.plot : logical; if true a plot of estimate with 95% CI is returned; default is TRUE
  # tol: tolerance for numerical stability
  # constrType: type of inequality constraints ('increasing','decreasing','boundedness')
  # lower and upper: options for only boundedness constr
  
  # OUTPUT: Posterior samples on xi,tau,sig and fhat with posterior mean, 95% CI of fhat, mode of posterior distribution z_star
  
  if (length(y) != length(x))
    stop("y and x should be of same length!")
  y <- y[order(x)]
  x <- sort(x)
  n <- length(y)
  delta <- 1/(N-1)
  my_knots <- seq(0,1,by=delta)
  X <- fcth(x,u=my_knots,N)
  
  if (missing(nu))
    stop("nu needs to be supplied")
  if (nu == 0)
    stop("nu cannot be zero")
  if (missing(constrType))
    stop('constrType should be specified and must be one of the following: \'increasing\',
         \'decreasing\', \'boundedness\'')
  if (length(constrType)>=4)
    stop('You cannot select more than three shape constraints simultaneously')
  if (!missing(l) && l==0)
    stop("l cannot be zero")
  if (missing(return.plot))
    return.plot <- TRUE
  if (!missing(sseed))
    set.seed(sseed)
  if (missing(sseed))
    set.seed(Sys.Date())
  if (missing(verbose))
    verbose <- TRUE
  # if (missing(l))
  #   l <- l_est(nu,range=range(my_knots),0.05)
  
  # prior covariance K:
  # K <- covmat(my_knots,nu,l)
  # # prior precision:
  # if (min(eigen(K, symmetric = TRUE)$values) <= 0) # numerical stability
  #   K <- K + tol*diag(nrow(K))
  # K_inv <- tinv(K)
  
  if (missing(eta))
    eta <- 200
  if (missing(tol))
    tol <- 1e-6
  if (missing(nsim))
    nsim <- 5000
  if (missing(burn.in))
    burn.in <- 1000
  if (missing(thin))
    thin <- 1
  
  ## Matrix and vector of linear constr
  A_b <- c() # for boundedness constr
  B_b <- c() # for boundedness constr
  A <- c() # for other constr
  B <- c() # for other constr
  if (any(constrType == 'boundedness')) {
    if (missing(lower) & missing(upper))
      stop('Error: lower or upper should be provided')
    if (any(lower >= upper)) 
      stop("The elements from \"upper\" has to be greater than the elements from \"lower\"")
    if (length(lower) == 1)
      lower <- rep(lower,N)
    if (length(upper) == 1)
      upper <- rep(upper,N)
    A_b <- rbind(A,constrSys(N=N,type='boundedness',lower=lower,upper=upper)$A)
    B_b <- c(B,constrSys(N=N,type='boundedness',lower=lower,upper=upper)$B)
  }
  
  if (any(constrType == 'increasing')) {
    A <- rbind(A,constrSys(N=N,type='increasing')$A)
    B <- c(B,constrSys(N=N,type='increasing')$B)
  }
  if (any(constrType == 'decreasing')) {
    A <- rbind(A,constrSys(N=N,type='decreasing')$A)
    B <- c(B,constrSys(N=N,type='decreasing')$B)
  }
  
  if (!missing(tau.fix))
    tau.in <- tau.fix
  if (!missing(sig.fix))
    sig.in <- sig.fix
  if (missing(tau.fix) && missing(tau.in))
    tau.in <- 1
  if (missing(sig.fix) && missing(sig.in))
    sig.in <- 1
  
  ## for the estimation of the length-scale para
  fzetoil <- function(l) {
    K_inv <- tinv(covmat(knot = my_knots, nu = nu, l = l))
    XXK <- crossprod(X) + K_inv
    solve.QP(XXK,dvec=as.vector(t(X)%*%y),
             Amat=t(rbind(A_b, A)),bvec=-c(B_b, B),meq=0)$solution
  }
  fMAP <- function(l) {
    X %*% fzetoil(l)
  }
  MSPE <- function(l) {
    return(mean((y - fMAP(l))^2))
  }
  if (missing(l)) {
    if (est.l == T) {
      opts <- list("algorithm"="NLOPT_LN_COBYLA",
                   "xtol_rel" = 1e-5)
      l <- nloptr(l_est(nu, range = range(my_knots), 0.05), MSPE, lb = 0.1, ub = 1, opts = opts)$solution
    }
    else if (est.l == F)
      l <- l_est(nu,range=range(my_knots),0.05)
  }
  if (missing(xi.in)) {
    xi.in <- fzetoil(l)
  }
  
  K_inv <- tinv(covmat(knot = my_knots, nu = nu, l = l))
  
  
  tau <- tau.in
  sig <- sig.in
  xi_in <- xi.in
  
  em <- nsim + burn.in
  ef <- nsim/thin
  xi_sam <- matrix(NA, nrow = N, ncol = ef)
  tau_sam <- rep(NA, ef)
  sig_sam <- rep(NA, ef)
  fhat_sam <- matrix(NA, nrow = n, ncol = ef)
  
  if (verbose)
    print("MCMC sample draws:")
  
  ptm <- proc.time()
  for (i in 1 : em) {
    # sampling Xi:
    nu.ess <- samp.WC(knot=my_knots, nu=nu, l=l, tausq=tau, sseedWC=i)
    xi_out <- ESS.linear(y=y, X=X, beta=xi_in, nu_ess=nu.ess, sigsq=sig,
                         eta=eta, A=A, B=B, lower=lower, upper=upper, constrType=constrType, seeds=i)
    if (any(constrType == 'increasing')) {
      xi_out <- increasing_vector(xi_out)
    }
    if (any(constrType == 'decreasing')) {
      xi_out <- decreasing_vector(xi_out)
    }
    if (any(constrType == 'boundedness')) {
      xi_out <- pmin(pmax(lower, xi_out), upper)
    }
    set.seed(2*i)
    # sampling \sigma^2:
    Xxi <- as.vector(X %*% xi_out)
    y_star <- y - Xxi
    if (missing(sig.fix))
      sig <- 1/rgamma(1, shape = n/2, rate = sum(y_star^2)/2)
    
    # sampling \tau^2:
    if (missing(tau.fix))
      tau <- 1/rgamma(1, shape = N/2, rate = (t(xi_out) %*% K_inv %*% xi_out)/2)
    
    # storing MCMC samples:
    if (i > burn.in && i%%thin == 0) {
      xi_sam[,(i-burn.in)/thin] <- xi_out
      sig_sam[(i-burn.in)/thin] <- sig
      tau_sam[(i-burn.in)/thin] <- tau
      fhat_sam[,(i-burn.in)/thin] <- Xxi
    }
    
    if (i %% 1000 == 0 && verbose) {
      print(i)
    }
    
    # renewing the intial value:
    xi_in <- xi_out
  }; tm <- proc.time()-ptm
  
  ## posterior Mode
  XXK <- crossprod(X)/mean(sig_sam)+K_inv/mean(tau_sam)
  z_star <- solve.QP(XXK,dvec=as.vector(t(X)%*%y)/mean(sig_sam),
                     Amat=t(rbind(A_b, A)),bvec=-c(B_b, B),meq=0)$solution
  MAP <- X%*%z_star 
  ## mAP estimate
  z_mean <- rowMeans(xi_sam)
  fmean <- rowMeans(fhat_sam) # mAP estimate
  qnt <- apply(fhat_sam,1,function(x) quantile(x,c(0.025,0.975),na.rm=TRUE))
  f_low <- qnt[1,]
  f_upp <- qnt[2,]
  ub <- max(f_low,f_upp,fmean,MAP)
  lb <- min(f_low,f_upp,fmean,MAP)
  
  if (return.plot) {
    par(mfrow=c(1, 1))
    par(mar=c(2.1, 2.1, 2.1, 1.1)) # adapt margins
    plot(x, y, pch = '*', lwd = 2, lty = 1, col = 'black',
         ylim = range(ub, lb, y), xlab = '', ylab = '')
    polygon(c(x, rev(x)), y = c(f_low, rev(f_upp)), border = F, col = 'gray')
    lines(x, fmean, type = 'l', lty = 4, lwd = 2, col = 'blue')
    lines(x, MAP, type = 'l', lty = 2, lwd = 2, col = 'red')
    points(x, y, pch = '*')
  }
  
  return(list("time" = tm, "xi_sam" = xi_sam, "sig_sam" = sig_sam, "tau_sam" = tau_sam,
              "fhat_sam" = fhat_sam, "fmean" = fmean,
              "f_low" = f_low, "f_upp" = f_upp, "z_star" = z_star,
              "MAP" = MAP, "knots" = my_knots, "z_mean" = z_mean, "length-scale" = l))
}
###################################################



### Function for drawing posterior samples using HMC with fixed hyperparameters \nu and \ell:
### For increasing,decreasing and boundedness functions estimation 
linCGP.HMC <- function(y, x, N, nu, l, est.l = F, nsim, burn.in, thin, tau.in, sig.in, xi.in, lower = -Inf, upper = Inf,
                       constrType, tau.fix, sig.fix, sseed, verbose, return.plot, tol) {
  # y:Response variable; x: vector to form design matrix X (n x N)
  # N: number of knots 
  # nu:smoothness parameter of Matern; l:length-scale parameter of Matern
  # est.l logical if TRUE, the estimated value of length-scale para "l" will be employed 
  # nsim, burn.in, thin : nsim samples, burning and thinning for MCMC
  # tau.in, sig.in, xi.in : initial values (supplied by user or use the default values)
  # tau.fix,sig.fix : if fixed values of the parameters are to use
  # verbose : logical; if TRUE, prints current status; default is TRUE
  # return.plot : logical; if true a plot of estimate with 95% CI is returned; default is TRUE
  # tol: tolerance for numerical stability
  # constrType: type of inequality constraints ('increasing','decreasing','boundedness')
  # lower and upper: bounds for only boundedness constr
  
  # OUTPUT: Posterior samples on xi,tau,sig, fhat with posterior mean, 95% CI of fhat, mode of posterior distribution z_star
  
  if (length(y) != length(x))
    stop("y and x should be of same length!")
  y <- y[order(x)]
  x <- sort(x)
  n <- length(y)
  delta <- 1/(N-1)
  my_knots <- seq(0, 1, by = delta)
  X <- fcth(x, u = my_knots, N = N)
  
  if (missing(nu))
    stop("nu needs to be supplied")
  if (nu == 0)
    stop("nu cannot be zero")
  if (missing(constrType))
    stop('constrType should be specified and must be one of the following: \'increasing\',
         \'decreasing\', \'boundedness\'')
  if (!missing(l) && l == 0)
    stop("l cannot be zero")
  if (missing(return.plot))
    return.plot <- TRUE
  if (!missing(sseed))
    set.seed(sseed)
  if (missing(sseed))
    set.seed(Sys.Date())
  if (missing(verbose))
    verbose <- TRUE
  if (missing(tol))
    tol <- 1e-6
  if (missing(nsim))
    nsim <- 5000
  if (missing(burn.in))
    burn.in <- 1000
  if (missing(thin))
    thin <- 1
  
  ## Matrix and vector of linear constr
  A <- c()
  B <- c()
  if (any(constrType == 'boundedness')) {
    if (missing(lower) & missing(upper))
      stop('Error: lower or upper should be provided')
    if (any(lower >= upper)) 
      stop("The elements from \"upper\" has to be greater than the elements from \"lower\"")
    if (length(lower) == 1)
      lower <- rep(lower, N)
    if (length(upper) == 1)
      upper <- rep(upper, N)
    A <- rbind(A, constrSys(N = N, type = 'boundedness', lower = lower, upper = upper)$A)
    B <- c(B, constrSys(N = N, type = 'boundedness', lower = lower, upper = upper)$B)
  }
  
  if (any(constrType == 'increasing')) {
    A <- rbind(A, constrSys(N = N, type = 'increasing')$A)
    B <- c(B, constrSys(N = N, type = 'increasing')$B)
  }
  if (any(constrType == 'decreasing')) {
    A <- rbind(A, constrSys(N = N, type = 'decreasing')$A)
    B <- c(B, constrSys(N = N, type = 'decreasing')$B)
  }
  
  if (!missing(tau.fix))
    tau.in <- tau.fix
  if (!missing(sig.fix))
    sig.in <- sig.fix
  if (missing(tau.fix) && missing(tau.in))
    tau.in <- 1
  if (missing(sig.fix) && missing(sig.in))
    sig.in <- 1
  
  
  ## for the estimation of the length-scale para
  fzetoil <- function(l) {
    K_inv <- tinv(covmat(knot = my_knots, nu = nu, l = l))
    XXK <- crossprod(X) + K_inv
    solve.QP(Dmat = XXK, dvec = as.vector(t(X) %*% y),
             Amat = t(A), bvec = -B, meq = 0)$solution
  }
  fMAP <- function(l) {
    X %*% fzetoil(l)
  }
  MSPE <- function(l) {
    return(mean((y - fMAP(l))^2))
  }
  if (missing(l)) {
    if (est.l == T) {
      opts <- list("algorithm"="NLOPT_LN_COBYLA",
                   "xtol_rel" = 1e-5)
      l <- nloptr(l_est(nu, range = range(my_knots), 0.05), MSPE, lb = 0.1, ub = 1, opts = opts)$solution
    }
    else if (est.l == F)
      l <- l_est(nu, range = range(my_knots), 0.05)
  }
  if (missing(xi.in)) {
    xi.in <- fzetoil(l)
  }
  
  K_inv <- tinv(covmat(knot = my_knots, nu = nu, l = l))
  
  
  tau <- tau.in
  sig <- sig.in
  xi_in <- xi.in
  
  em <- nsim + burn.in
  ef <- nsim/thin
  xi_sam <- matrix(NA, nrow = N, ncol = ef)
  tau_sam <- rep(NA, ef)
  sig_sam <- rep(NA, ef)
  fhat_sam <- matrix(NA, nrow = n, ncol = ef)
  
  if (verbose)
    print("MCMC sample draws:")
  
  ptm <- proc.time()
  for (i in 1 : em) {
    # sampling Xi:
    M <- crossprod(X)/sig + K_inv/tau
    r <- as.vector(t(X) %*% y)/sig
    xi_out <- as.vector(rtmg(n = 1, M = M, r = r, initial = xi_in, f = A, g = B + tol, burn.in = 0))
    set.seed(2 * i)
    # sampling \sigma^2:
    Xxi <- as.vector(X %*% xi_out)
    y_star <- y - Xxi
    if (missing(sig.fix))
      sig <- 1/rgamma(1, shape = n/2, rate = sum(y_star^2)/2)
    
    # sampling \tau^2:
    if (missing(tau.fix))
      tau <- 1/rgamma(1, shape = N/2, rate = (t(xi_out) %*% K_inv %*% xi_out)/2)
    
    # storing MCMC samples:
    if (i > burn.in && i %% thin == 0) {
      xi_sam[, (i-burn.in)/thin] <- xi_out
      sig_sam[(i-burn.in)/thin] <- sig
      tau_sam[(i-burn.in)/thin] <- tau
      fhat_sam[, (i-burn.in)/thin] <- Xxi
    }
    
    if (i %% 1000 == 0 && verbose) {
      print(i)
    }
    
    # renewing the intial value:
    xi_in <- xi_out
  } 
  tm <- proc.time() - ptm
  
  ## posterior Mode
  XXK <- crossprod(X)/mean(sig_sam) + K_inv/mean(tau_sam)
  z_star <- solve.QP(Dmat = XXK, dvec = as.vector(t(X) %*% y)/mean(sig_sam),
                     Amat = t(A), bvec = -B, meq = 0)$solution
  MAP <- X %*% z_star 
  ## mAP estimate
  z_mean <- rowMeans(xi_sam)
  fmean <- rowMeans(fhat_sam) # mAP estimate
  qnt <- apply(fhat_sam, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
  f_low <- qnt[1, ]
  f_upp <- qnt[2, ]
  ub <- max(f_low, f_upp, fmean, MAP)
  lb <- min(f_low, f_upp, fmean, MAP)
  
  if (return.plot) {
    par(mfrow = c(1, 1))
    par(mar = c(2.1, 2.1, 2.1, 1.1)) # adapt margins
    plot(x, y, pch = '*', lwd = 2, lty = 1, col = 'black',
         ylim = range(ub, lb, y), xlab = '', ylab = '')
    polygon(c(x, rev(x)), y = c(f_low, rev(f_upp)), border = F,col = 'gray')
    lines(x, fmean, type = 'l', lty = 4, lwd = 2, col = 'blue')
    lines(x, MAP, type = 'l', lty = 2, lwd = 2, col = 'red')
    points(x, y, pch = '*')
  }
  
  return(list("time" = tm, "xi_sam" = xi_sam, "sig_sam" = sig_sam, "tau_sam" = tau_sam,
              "fhat_sam" = fhat_sam, "fmean" = fmean,
              "f_low" = f_low, "f_upp" = f_upp, "z_star" = z_star,
              "MAP" = MAP, "knots" = my_knots, "z_mean" = z_mean, "length-scale" = l))
}
###################################################



## end               
