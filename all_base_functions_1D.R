### All required functions for using ESS
### Functions related to Wood and Chan algorithm of drawing samples
### MH for sampling from \nu and \ell
### Covariance matrix and design matrix (using basis function) are also defined
### And all related and dependant functions are here

### Required libraries:
library(FastGP) # for tinv and rcpp_rmvnorm function
library(rSPDE) # Matern cov fct
library(mvtnorm);library(MASS)
library(quadprog) # for the MAP
library(Rfast) # matrnorm



## 1D Gaussian covariance function with length-scale l
kGaus <- function(h, l) {
  exp(-0.5 * (h/l))
}

# Covariance matrix for Gaussian kernels
covmat_Gaus <- function(knot, l) {
  kGaus(outer(knot, knot, '-'), l = l) 
}

####################################################
######### Mat\'ern family covariance kernels #######
####################################################
#Given a \nu (smoothness parameter of matern kernel) finding a value of 
# l (length-scale parameter) such that the correlation between the 
# maximum seperation is some small value, say 0.05

# 1D Matern kernel with smoothness nu and length-scale l:
k <- function(h, nu, l) {
  matern.covariance(h = h, kappa = sqrt(2*nu)/l, nu = nu, sigma = 1)
}

## function for uniroot:
fl <- function(l, para) { 
  #para[1]=x, para[2]=y and para[3]=nu of MK : Matern kernel function;
  #para[4]=pre-specified value of the correlation
  a <- k(h = abs(para[1]-para[2]), nu = para[3], l = l)
  return(a - para[4])
}

## function for estimating l:
l_est <- function(nu, range, val) {
  # nu : smoothness; range : c(min, max) of the range of variable
  # val : pre-specified value of the correlation between the maximum seperation
  para <- c(range[1], range[2], nu, val)
  rl <- uniroot(f = fl, interval = c(0.000001, 100000), para)
  return(rl$root)
}

# Covariance matrix for Matern family of kernels
covmat <- function(knot, nu, l) {
  k(outer(knot, knot, '-'), nu = nu, l = l) 
}



### Define design matrix ###
### The basis functions ###

####################################################
######## Maatouk & Bay2017 Basis functions #########
####################################################
## hat basis functions
h <- function(x) {
  ifelse(x >= -1 & x <= 1, 1 - abs(x), 0)
}
hi <- function(x, u, i) {
  delta <- (max(u) - min(u)) / (length(u) - 1)
  h((x - u[i]) / delta)
}


####################################################
########## function of design matrix ###############
####################################################
## fct design matrix (hat basis function)
fcth <- function(x, u, N) {
  n <- length(x)
  h <- matrix(NA, nrow = n, ncol = N)
  for (j in 1 : N) {
    h[, j] <- hi(x, u, i = j)
  }
  return(h)
}
####################################################



## Order of the circulant matrix:
## minimum value of g and m so that G can be embedded into C
min_g <- function(knot) {
  N <- length(knot)
  g <- ceiling(log(2 * N, 2))   #m=2^g and m>=2(n-1) : Wood & Chan notation; 
  #since we are going upto n and not stopping at (n-1), the condition is modified!
  return("g" = g)
}

## forming the circulant matrix:
circulant <- function(x) {
  n <- length(x)
  mat <- matrix(0, nrow = n, ncol = n)
  for (j in 1 : n) {
    mat[j, ] <- c(x[-(1 : (n+1-j))], x[1 : (n+1-j)])
  }
  return(mat)
}

## Function for forming the vector of circulant matrix:
circ_vec <- function(knot, g, nu, l, tausq) {
  delta_N <- 1 / (length(knot) - 1)
  m <- 2**g
  cj <- integer()
  for (j in 1 : m) {
    if (j <= (m/2))
      cj[j] <- (j-1) * delta_N
    else
      cj[j] <- (m-(j-1)) * delta_N
  }
  x <- (tausq * k(cj, nu, l))
  return(x)
}


## Function for finding a g such that C is nnd:
## without forming the circulant matrix and without computing eigen values:
C.eval <- function(knot, g, nu, l, tausq) {
  vec <- circ_vec(knot, g, nu, l, tausq)
  val <- fft(vec) # eigenvalues will be real as the circulant matrix formed by the 
  # vector is by construction is symmetric!
  ev <- min(Re(val))
  return(list("vec" = vec, "min.eig.val" = ev))
}


nnd_C <- function(knot, g, nu, l, tausq) {
  C.vec <- C.eval(knot, g, nu, l, tausq)$vec
  eval <- C.eval(knot, g, nu, l, tausq)$min.eig.val
  if (eval > 0)
    return(list("cj" = C.vec, "g" = g))
  else {
    g <- g + 1
    nnd_C(knot, g, nu, l, tausq)
  }
}

## computing the eigen values of C using FFT:
eigval <- function(knot, nu, l, tausq) {
  g <- min_g(knot)
  c.j <- nnd_C(knot, g, nu, l, tausq)$cj
  lambda <- Re(fft(c.j))
  if (min(lambda) > 0)
    return(lambda)
  else
    stop("nnd condition is NOT satisfied!!!")
}


####################################################
### Samples drawn using Wood and Chan Algorithm ####
####################################################
samp.WC <- function(knot, nu, l, tausq, sseedWC = 1) {
  N <- length(knot)
  lambda <- eigval(knot, nu, l, tausq)
  m <- length(lambda)
  samp.vec <- rep(0, N)
  set.seed(sseedWC)
  a <- rep(0, m)
  a[1] <- sqrt(lambda[1]) * rnorm(1) / sqrt(m)
  a[(m/2)+1] <- sqrt(lambda[(m/2)+1]) * rnorm(1) / sqrt(m)
  i <- sqrt(as.complex(-1))
  for (j in 2 : (m/2)) {
    uj <- rnorm(1) 
    vj <- rnorm(1)
    a[j] <- (sqrt(lambda[j]) * (uj + i * vj)) / (sqrt(2 * m))
    a[m+2-j] <- (sqrt(lambda[j]) * (uj - i * vj)) / (sqrt(2 * m))
  }
  samp <- fft(a)
  samp.vec <- Re(samp[1 : N])
  return(samp.vec)
}
####################################################



####################################################
############## Functions for using ESS #############
####################################################
## ESS for linear constraints AX+B>=0_m
ESS.linear <- function(y, X, beta, nu_ess, eta, sigsq, A, B, seeds = 1, constrType = c('increasing', 'decreasing', 'convex', 'concave', 'boundedness'), lower = -Inf, upper = Inf) {
  
  set.seed(seeds)
  u <- runif(n = 1)
  logy <- loglik_linear(y = y, X = X, eta = eta, beta = beta, sigsq = sigsq, A = A, B = B, lower = lower, upper = upper, constrType = constrType) + log(u)
  
  theta <- runif(n = 1, min = 0, max = 2 * pi) 
  thetamin <- theta - 2 * pi
  thetamax <- theta + 2 * pi
  beta_star <- beta * cos(theta) + nu_ess * sin(theta)
  
  while (loglik_linear(y = y, X = X, eta = eta, beta = beta_star, sigsq = sigsq, A = A, B = B, lower = lower, upper = upper, constrType = constrType) <= logy) {
    if (theta < 0)
      thetamin <- theta
    else
      thetamax <- theta
    theta <- runif(n = 1, min = thetamin, max = thetamax)
    beta_star <- beta * cos(theta) + nu_ess * sin(theta)
  }
  return(beta_star)       
}
####################################################





####################################################
## Defining the loglik function to be used in ESS ## 
### loglik calculates the log of the likelihood: ###
####################################################

## loglik linear constr AX+B>=0_m
loglik_linear <- function(y, X, sigsq, eta, beta, A, B, lower = -Inf, upper = Inf, constrType = c('increasing', 'decreasing', 'convex', 'concave', 'boundedness')) {
  mu <- y - (X %*% beta)
  J_eta <- c()
  if (any(constrType == 'boundedness')) {
    J_eta <- c(J_eta, -log(1+sum(exp(-eta*(upper-beta))+exp(-eta*(beta-lower))+exp(-eta*(upper-lower)))))
  }
  if (any(constrType == 'increasing') || any(constrType == 'decreasing') || any(constrType == 'convex') || any(constrType == 'concave'))
    J_eta <- c(J_eta, -log(1+sum(exp(-eta*(A%*%beta+B)))))
  
  val <- sum(J_eta) - sum(mu^2) / (2 * sigsq)
  return(val)
}
####################################################



####################################################
############## Fast Large-scale ####################
####################################################

### This function generates a multivariate normal (MVN) distribution with a covariance matrix extracted from a Mat\'ern kernel (MK)
### M represents the number of subdomains (M>=1)
### N1 represents the number of equally-spaced points of the 1st subdomain
### u: grid vector, should be of length M times N1
### nu represents the smoothness parameter of the MK (nu>0)
### l represents the length-scale parameter of MK (l>0)
### tol: tolerance (relative to largest variance) for numerical lack of a positive-definiteness problem

library(mvtnorm) # when using `rmvnorm' function
# library(FastGP) # when using `rcppeigen_get_chol' for Cholesky factorization

Fast.LS <- function(u, M, N1, nu, l, tausq, tol, sseedLS = 1) {
  if (missing(tol)) {
    tol <- 1e-5
  }
  if (missing(M)) {
    M <-1
  }
  if (M == 0)
    stop("M cannot be zero")
  if (N1 == 0)
    stop("N1 cannot be zero")
  if (length(u) != M*N1)
    stop("The length of the vector u must be M times N1")
  u1 <- u[1 : N1]
  u2 <- u[(N1+1) : (2*N1)]
  Gamma11 <- tausq * covmat(u1, nu = nu, l = l) #+ tol * diag(N1)
  if (min(eigen(Gamma11, symmetric = TRUE)$values) <= 0) # numerical stability
    Gamma11 <- Gamma11 + tol * diag(N1)
  
  set.seed(sseedLS)
  if (M == 1) {
    # return(as.vector(rcpp_rmvnorm(n = 1, S = Gamma11, mu = rep(0, N1))))
    return(as.vector(Rfast::rmvnorm(n = 1, mu = rep(0, N1), sigma = Gamma11)))
  }
  else 
    Gamma12 <- tausq * k(outer(u2, u1, '-'), nu = nu, l = l)
  Ktilde <- Gamma12 %*% tinv(Gamma11) # coupling matrix
  L <- t(chol(Gamma11 - Gamma12 %*% t(Ktilde) + tol * diag(N1))) %*% solve(t(chol(Gamma11)))
  eta <- Rfast::rmvnorm(n = M, mu = rep(0, N1), sigma = Gamma11)
    # rcpp_rmvnorm(n = M, S = Gamma11, mu = rep(0,N1))
  etaT <- matrix(NA, nrow = M, ncol = N1)
  etaT[1, ] <- eta[1, ]
  for (i in 2 : M) {
    etaT[i, ] <- Ktilde %*% (etaT[i-1, ]) + L %*% (eta[i, ])
  }
  return(as.vector(t(etaT)))
}


## Combination Fast Large-scale and FFT from WC

Fast.LS.WC <- function(u, M, N1, nu, l, tausq, tol, sseedWC = 1) {
  if (missing(tol)) {
    tol <- 1e-6
  }
  if (missing(M)) {
    M <- 1
  }
  if (M == 0)
    stop("M cannot be zero")
  if (N1 == 0)
    stop("N1 cannot be zero")
  if (length(u) != M*N1)
    stop("The length of the vector u must be M times N1")
  u1 <- u[1 : N1]
  u2 <- u[(N1+1) : (2*N1)]
  Gamma12 <- tausq * k(outer(u2, u1, '-'), nu = nu, l = l)
  Gamma11 <- tausq * covmat(u1, nu = nu, l = l) + tol * diag(N1)
  if (M == 1) {
    return(samp.WC(u1, nu = nu, l = l, tausq = tausq, sseedWC))
  }
  else 
    # K <- Gamma12
    Ktilde <- Gamma12 %*% tinv(Gamma11)
  L <- t(chol(Gamma11 - Gamma12 %*% t(Ktilde))) %*% solve(t(chol(Gamma11)))
  eta <- matrix(NA, nrow = M, ncol = N1)
  eta[1, ] <- samp.WC(u1, nu = nu, l = l, tausq = tausq, sseedWC)
  etaT <- matrix(NA, nrow = M, ncol = N1)
  etaT[1, ] <- eta[1, ]
  for (i in 2 : M) {
    eta[i, ] <- samp.WC(u1, nu = nu, l = l, tausq = tausq, sseedWC)
    etaT[i, ] <- Ktilde %*% (etaT[i-1, ]) + L %*% (eta[i, ])
  }
  return(as.vector(t(etaT)))
}
####################################################

### Fast Large-scale for more than one sample
Fast.LS_v <- function(nsim, u, M, N1, nu, l, tausq, tol, sseedLS = 1) {
  if (missing(tol)) {
    tol <- 1e-6
  }
  if (missing(M)) {
    M <- 1
  }
  if (M == 0)
    stop("M cannot be zero")
  if (N1 == 0)
    stop("N1 cannot be zero")
  if (length(u) != M*N1)
    stop("The length of the vector u must be M times N1")
  u1 <- u[1 : N1]
  u2 <- u[(N1+1) : (2*N1)]
  Gamma11 <- tausq * covmat(u1, nu = nu, l = l) + tol * diag(N1)
  set.seed(sseedLS)
  if (M == 1) {
    return(Rfast::rmvnorm(n = nsim, mu = rep(0, N1), sigma = Gamma11))
  }
  else 
    Gamma12 <- tausq * k(outer(u2, u1, '-'), nu = nu, l = l)
  Ktilde <- Gamma12 %*% tinv(Gamma11) # coupling matrix
  L <- t(chol(Gamma11 - Gamma12 %*% t(Ktilde))) %*% solve(t(chol(Gamma11)))
  eta <- t(Rfast::rmvnorm(n = nsim, mu = rep(0, N1), sigma = Gamma11))
  etaT <- list()
  etaT[[1]] <- eta
  for (i in 2 : M) {
    eta <- t(Rfast::rmvnorm(n = nsim, mu = rep(0, N1), sigma = Gamma11))
    etaT[[i]] <- Ktilde %*% (etaT[[i-1]]) + L %*% eta
  }
  return(do.call(rbind,etaT))
}
####################################################




####################################################
################ LS KLE one sample #################
####################################################

### This function generates a multivariate normal (MVN) distribution with a covariance matrix extracted from a Mat\'ern kernel (MK)
### M represents the number of subdomains (M >= 1)
### N1 represents the number of equally-spaced points of the 1st subdomain
### p represents the truncated parameter of the KLE approach
### u: grid vector, should be of length M times N1
### nu represents the smoothness parameter of the MK (nu > 0)
### l represents the length-scale parameter of MK (l > 0)
### tol: tolerance (relative to largest variance) for numerical lack of a positive-definiteness problem

LS.KLE <- function(u, N1, p, M, nu, l, tausq, tol, sseedLS = 1) {
  if (missing(tol)) {
    tol <- 1e-8
  }
  if (missing(p)) {
    p <- N1
  }
  if (missing(M)) {
    M <- 1
  }
  if (M == 0)
    stop("M cannot be zero")
  if (N1 == 0)
    stop("N1 cannot be zero")
  if (length(u) != M*N1)
    stop("The length of the vector u must be M times N1")
  u1 <- u[1 : N1]
  u2 <- u[(N1+1) : (2*N1)]
  Gamma11 <- tausq * covmat(u1, nu, l) + tol * diag(N1)
  set.seed(sseedLS)
  if (M == 1) {
    # return(as.vector(rcpp_rmvnorm(n=1,S=Gamma11,mu=rep(0,N1))))
    return(Rfast::rmvnorm(n = 1, mu = rep(0, N1), sigma = Gamma11))
  }
  else 
    Gamma12 <- tausq * k(outer(u1, u2, '-'), nu, l)
  eig11 <- eigen(Gamma11)
  value11 <- eig11$values[1 : p]
  vector11 <- eig11$vectors[, 1 : p]
  K12 <- (((t(vector11) %*% (Gamma12)) %*% (vector11))/
            sqrt(tcrossprod(value11)))
  L12 <- t(chol(diag(p)-crossprod(K12)))
  eta <- matrnorm(M, p)
  etaT <- matrix(NA, nrow = M, ncol = p)
  f <- matrix(NA, nrow = M, ncol = N1)
  f[1, ] <- vector11 %*% (sqrt(value11) * eta[1, ])  
  etaT[1, ] <- eta[1, ]
  for (i in 2 : M) {
    etaT[i, ] <- t(K12) %*% etaT[(i-1), ] + L12 %*% eta[i, ]
    f[i, ] <- vector11 %*% (sqrt(value11) * etaT[i, ])  
  }
  return(as.vector(t(f)))
}
####################################################



###################################################
########## LS KLE more than one sample ############
###################################################

LS.KLE_v <- function(nsim, u, N1, p, M, nu, l, tausq, tol, sseedLS = 1) {
  if (missing(nsim)) {
    nsim <- 10
  }
  if (missing(tol)) {
    tol <- 1e-8
  }
  if (missing(p)) {
    p <- N1
  }
  if (missing(M)) {
    M <- 1
  }
  if (M == 0)
    stop("M cannot be zero")
  if (N1 == 0)
    stop("N1 cannot be zero")
  if (length(u) != M*N1)
    stop("The length of the vector u must be M times N1")
  u1 <- u[1 : N1]
  u2 <- u[(N1+1) : (2*N1)]
  Gamma11 <- tausq * covmat(u1, nu, l) + tol * diag(N1)
  if (M == 1) {
    return(t(Rfast::rmvnorm(n = 1, mu = rep(0, N1), sigma = Gamma11)))
  }
  else 
    Gamma12 <- tausq * k(outer(u1, u2, '-'), nu, l)
  eig11 <- eigen(Gamma11)
  value11 <- eig11$values[1 : p]
  vector11 <- eig11$vectors[, 1 : p]
  K12 <- (((t(vector11) %*% (Gamma12)) %*% (vector11))/
            sqrt(tcrossprod(value11)))
  L12 <- t(chol(diag(p) - crossprod(K12)))
  eta <- matrnorm(n = p, p = nsim)
  etaT <- list()
  etaT[[1]] <- eta
  f <- list()
  f[[1]] <- vector11 %*% (sqrt(value11) * eta)
  for (i in 2 : M) {
    eta <- matrnorm(n = p, p = nsim)
    etaT[[i]] <- t(K12) %*% etaT[[i-1]] + L12 %*% eta
    f[[i]] <- vector11 %*% (sqrt(value11) * etaT[[i]])
  }
  return(do.call(rbind, f))
}
###################################################



####################################################
############## Cholesky one sample #################
####################################################

### This function generates a multivariate normal (MVN) distribution with a covariance matrix extracted from a Mat\'ern kernel (MK) using Cholesky decomposition
### N represents the size of the MVN
### u: grid vector, should be of length M times N1
### nu represents the smoothness parameter of the MK (nu > 0)
### l represents the length-scale parameter of MK (l > 0)
### tol: tolerance (relative to largest variance) for numerical lack of a positive-definiteness problem


chol_prior <- function(u, N, nu, l, tol = 1e-6, tausq, seedchol = 1) {
  
  Sigma <- tausq * covmat(knot = u, nu = nu, l = l)
  if (min(eigen(Sigma, symmetric = TRUE)$values) <= 0) # numerical stability
    Sigma <- Sigma + tol * diag(N)
  set.seed(seedchol)
  return(as.vector(mvtnorm::rmvnorm(1, mean = rep(0, N), sigma = Sigma, method = "chol")))
}




####################################################
###### maximum between consecutive increasing ######
####################################################
increasing_vector <- function(input_vector) {
  max_values <- pmax(input_vector[-length(input_vector)], input_vector[-1])
  result_vector <- cummax(c(input_vector[1], max_values))
  return(result_vector)
}



####################################################
##### minimum between consecutive decreasing #######
####################################################
decreasing_vector <- function(input_vector) {
  max_values <- pmin(input_vector[-length(input_vector)], input_vector[-1])
  result_vector <- cummin(c(input_vector[1], max_values))
  return(result_vector)
}



####################################################
#### return vector s.t. xi[i+1] >= 2xi[i]-xi[1] ####
####################################################
convex_vector <- function(input) {
  n <- length(input)
  if (n < 3) {
    return(input)
  }
  # output <- input
  for (i in 3 : n) {
    input[i] <- pmax(input[i], 2 * input[i-1] - input[i-2])
  }
  return(input)
}
####################################################


####################################################
#### return vector s.t. xi[i+1] <= 2xi[i]-xi[1] ####
####################################################
concave_vector <- function(input) {
  n <- length(input)
  if (n < 3) {
    return(input)
  }
  for (i in 3 : n) {
    input[i] <- pmin(input[i], 2 * input[i-1] - input[i-2])
  }
  return(input)
}
####################################################




####################################################
## vect s.t. xi[i+1]>= 2xi[i]-xi[1] & xi[i+1]>=xi[i]
####################################################
inc_conv_vector <- function(input) {
  n <- length(input)
  if (n < 2) {
    return(input)
  }
  else if (n == 2) {
    input[2] <- pmax(intput[1], input[2])
    return(input)
  }
  else {
    output <- input
    for (i in 3 : n) {
      output[i] <- pmax(output[i], 2 * output[i-1] - output[i-2], output[i-1])
    }
    return(output)
  }
}
####################################################




####################################################
## vect s.t. xi[i+1]<= 2xi[i]-xi[1] & xi[i+1]>=xi[i]
####################################################
inc_conc_vector <- function(input) {
  n <- length(input)
  if (n < 2) {
    return(input)
  }
  else if (n == 2) {
    input[2] <- pmax(intput[1], input[2])
    return(input)
  }
  else {
    output <- input
    for (i in 3 : n) {
      output[i] <- pmax(pmin(output[i], 2*output[i-1]-output[i-2]), output[i-1]) 
    }
    return(output)
  }
}
####################################################



####################################################
## vect s.t. xi[i+1]>= 2xi[i]-xi[1] & xi[i+1]<=xi[i]
####################################################
dec_conv_vector <- function(input) {
  n <- length(input)
  if (n < 2) {
    return(input)
  }
  else if (n == 2) {
    input[2] <- pmin(intput[1], input[2])
    return(input)
  }
  else {
    output <- input
    for (i in 3 : n) {
      output[i] <- pmin(pmax(output[i], 2 * output[i-1] - output[i-2]), output[i-1])
    }
    return(output)
  }
}
####################################################



####################################################
## vect s.t. xi[i+1]<= 2xi[i]-xi[1] & xi[i+1]<=xi[i]
####################################################
dec_conc_vector <- function(input) {
  n <- length(input)
  if (n < 2) {
    return(input)
  }
  else if (n == 2) {
    input[2] <- pmin(intput[1], input[2])
    return(input)
  }
  else {
    output <- input
    for (i in 3 : n) {
      output[i] <- pmin(output[i], 2 * output[i-1] - output[i-2])
    }
    return(output)
  }
}
####################################################






####################################################
##### matrix for nondecreasing constr in 1D ########
####################################################
Amat1D <- function(N) {
  A <- -diag(N)
  A[cbind(1 : (N-1), 2 : N)] <- 1
  if (N == 2)
    return(t(as.matrix(A[-N, ])))
  else 
    return(A[-N, ])
}
####################################################



####################################################
######## matrix for convex constr in 1D ############
####################################################
AC1D <- function(N) {
  if (N <= 2)
    stop('Error: N should be greater than 3')
  else if (N == 3)
    return(t(as.matrix(c(1, -2, 1))))
  
  else {
    A <- diag(N)
    A[cbind(1 : (N-2), 2 : (N-1))] <- -2
    A[cbind(1 : (N-2), 3 : N)] <- 1
    return(A[-c((N-1), N), ])
  }
}
####################################################






####################################################
############ Error measure calculations ############
####################################################

err.meas.report <- function(ytest, ytr, mu, out, type = "all") {
  ## ytr: training data
  ## ytest: testing data
  ## mu: estimate at xtest
  ## out: output of the model (sig_sam,...)
  
  # precomputing terms
  n <- length(ytr)
  sig_sam <- out$sig_sam
  fhat_sam <- out$fhat_sam
  nsim <- length(sig_sam)
  varsigma <- mean(sig_sam)
  
  # computing the errors
  # WAIC error
  likval <- matrix(0,nrow = n, ncol = nsim)
  for (i in 1 : n) {
    for (j in 1 : nsim) {
      likval[i,j] <- dnorm(ytr[i], mean = fhat_sam[i,j], sd = sqrt(sig_sam[j]))
    }
  }
  lppd <- sum(log(rowMeans(likval)))
  p_waic1 <- 2 * sum(log(rowMeans(likval)) - rowMeans(log(likval)))
  waic1 <- -2 * (lppd - p_waic1)
  #############
  error_abs <- abs(ytest - mu)
  error_sq <- (ytest - mu)^2
  std_error_sq <- error_sq/var(ytest)
  Q2 <- 1 - sum(error_sq)/sum((ytest - mean(ytest))^2)
  pva <- error_sq/varsigma
  # sserror <- ((ytest - mu)^2)/(2*varsigma)
  # logprob <- 0.5*log(2*pi*varsigma) + sserror -
  #   0.5*log(2*pi*var(y)) - ((ytest - mean(y))^2)/(2*var(y))
  # logprob <- logprob[complete.cases(logprob)]
  # cia <- ytest >= (mu-control$nsigma*sqrt(varsigma)) &
  #   ytest <= (mu+control$nsigma*sqrt(varsigma))
  
  error <- c(mae = mean(error_abs), mse = mean(error_sq),
             rmse = sqrt(mean(error_sq)),
             smse = mean(std_error_sq), #msll = mean(logprob),
             Q2 = Q2, pva = abs(log(mean(pva))), WAIC=waic1)#, cia = mean(cia))
  switch(type,
         mae = {return(error["mae"])},
         mse = {return(error["mse"])},
         rmse = {return(error["rmse"])},
         smse = {return(error["smse"])},
         # msll = {return(error["msll"])},
         Q2 = {return(error["Q2"])},
         pva = {return(error["pva"])},
         WAIC = {return(error["WAIC"])},
         # cia = {return(error["cia"])},
         {return(as.matrix(error, nrow = 1))}
  )
}
####################################################



### WAIC criterion function
fun.waic <- function(y, fhat_sam, sig_sam) {
  n <- length(y)
  nsim <- length(sig_sam)
  
  likval <- matrix(0, nrow = n, ncol = nsim)
  for (i in 1 : n) {
    for (j in 1 : nsim) {
      likval[i,j] <- dnorm(ytr[i], mean = fhat_sam[i, j], sd = sqrt(sig_sam[j]))
    }
  }
  lppd <- sum(log(rowMeans(likval)))
  p_waic1 <- 2 * sum(log(rowMeans(likval)) - rowMeans(log(likval)))
  waic1 <- -2 * (lppd - p_waic1)
  return(waic1)
}






constrSys <- function(N, type = c('increasing', 'decreasing', 'convex', 'concave', 'boundedness'), lower = -Inf, upper = Inf) {
  # type <- match.arg(type,type,several.ok=T)
  type <- match.arg(type, several.ok = T) # allow using abbreviation
  switch(type,
         increasing = {
           result <- list('A' = Amat1D(N), 'B' = rep(0, N-1))
         },
         decreasing = {
           result <- list('A' = -Amat1D(N), 'B' = rep(0, N-1))
         },
         convex = {
           result <- list('A' = AC1D(N), 'B' = rep(0, N-2))
         },
         concave = {
           result <- list('A' = -AC1D(N), 'B' = rep(0, N-2))
         },
         boundedness = {
           if (missing(lower) & missing(upper))
             stop('Error: lower or upper should be provided')
           A <- rbind(diag(N), -diag(N))
           B <- c(-lower, upper)
           if (any(B == Inf)) {
             idx <- which(B == Inf)
             B <- B[-idx]
             A <- A[-idx, ]
             if (length(B) == 0) stop("bounds are not defined")
           }
           result <- list('A' = A, 'B' = B)
         }
  )
  return(result)
  
}



## end
