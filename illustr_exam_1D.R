source('all_base_functions_1D.R')
source('all_models_1D.R')
library(ggplot2) ## plot zoom inset
library(grid) ## viewport function
library(DiceDesign) ## lhsDesign

## choice of the numerical example
## put 'yes' between quotations
multi_monot_convex_bd_syn = 'yes' # for multiple shape constraints (monotonicity, convexity and boundedness)
bound_syn = '' # for boundedness constraints [a;b] using LS-ESS and HMC
comp_mon_HMC_LS_WC = '' # for synthetic monotone nonincreasing functions (comparison between LS-ESS, LS-WC-FFT and HMC)
####################################################


if (multi_monot_convex_bd_syn == 'yes') {
  ###### true monotone, bounded and convex function 
  f <- function(x) {
    x^2 #  inc, convex and bounded function between 0 & 1
   }
  ## data
  set.seed(12345)
  n <- 100
  sigN <- 0.1 # noise sd
  xtr <- runif(n, min = 0, max = 1)
  ytr <- f(xtr) + rnorm(n, 0, sd = sigN)
  N1 <- 9 # size of the 1st subdomain
  M <- 5 # nb of subdomains
  N <- N1 * M # nb of knot points
  eta <- 5000 # smooth approximate parameter
  nsim <- 5000 # nb of mcmc retained samples
  brn <- 1000 # burn in
  thin <- 1
  nu <- 2.5 # smoothness kernel parameter
  l <- round(l_est(nu, c(0, 1), 0.05), 2) # length-scale
  upper <- 1 # upper bound for boundedness constr
  
  ## LS-ESS with multiple shape constraints (increasing, boundedness and convexity)
  post_lin_LS <- linCGP.ESS(y = ytr, x = xtr, N1 = N1, M = M, nu = nu, l = l, eta = eta, nsim = nsim, burn.in = brn, thin = thin,
                            sig.in = sigN^2, tau.in = 1, constrType = c('increasing', 'boundedness', 'convex'),
                            prior = 'Fast.LS', return.plot = T, tol = 0, sseed = 12345, lower = 0, upper = 1)
  lines(sort(xtr), f(sort(xtr)), type = 'l', lwd = 2)
}
####################################################









####################################################
############ Boundedness constraints ###############
####################################################


if (bound_syn == 'yes') {
  ###### bounded synthetic data ##########
  f <- function(x) {
    ## Andrew Pensoneault (Nonnegativity enforced GPR 2020)
    1/(1+(10*x)^4)+0.5*exp(-100*(x-0.5)^2) # [0;1]
  }
  ## data
  set.seed(12345)
  sigN <- 0.1 # sd noise
  n <- 100 # nb of training data
  xtr <- lhsDesign(n = n, dimension = 1, seed = 12345)$design
  ytr <- f(xtr) + rnorm(n, 0, sd = sigN)
  lower <- 0    # lower bound
  upper <- 1    # upper bound 
  N1 <- 9   # size of the 1st subdomain
  M <- 3     # nb of subdomains
  N <- N1 * M
  eta <- 100 # pdf approximation parameter
  nsim <- 5000 # nb of simulation
  brn <- 1000 # burn in 
  thin <- 1
  nu <- 1.5 # smoothness kernel parameter
  l <- round(l_est(nu, c(0, 1), 0.05), 2) # length-scale
  tol <- 1e-11
  ## LS-ESS approach
  post_bound_LS <- linCGP.ESS(y = ytr, x = xtr, N1 = N1, M = M, nu = nu, l = l, eta = eta, nsim = nsim, burn.in = brn, thin = thin,
                              sig.in = sigN^2, tau.in = 1, constrType = 'boundedness', lower = lower, upper = upper,
                              prior = 'Fast.LS', return.plot = T, tol = tol, sseed = 12345)
  abline(h = c(lower, upper), lty = 2, lwd = 2)
  x <- seq(from = 0, to = 1, length = 100)
  lines(x, f(x), type = 'l', lwd = 2)
  mtext(text =  'Large-scale ESS', side = 3, line = 0.8, cex = 0.8)
  mtext(text =  paste("Running Time (s) = ", round(post_bound_LS$time[3],2), "and dimension N = ",N), side = 3, line = 0.1, cex = 0.8)
  legend(0.55, 1,
         c("true function", "MAP", "mAP"),
         col = c('black', 'red', 'blue'), 
         text.col = "black", lty = c(1, 2, 4),
         lwd = c(2, 2, 2), text.font = 1, box.lty = 0, cex = 0.8,
         bg = 'transparent')
  ## Illustration FFT-WC
  post_bound_WC <- linCGP.WC.ESS(y = ytr, x = xtr, N = N, nu = nu, l = l, eta = eta, nsim = nsim, burn.in = brn, thin = thin,
                                 sig.in = sigN^2, tau.in = 1, constrType = 'boundedness',
                                 return.plot = T, sseed = 12345, lower = lower, upper = upper)
  abline(h = c(upper, lower), lty = 2, lwd = 2)
  x <- seq(from = 0, to = 1, length = 100)
  lines(x, f(x), type = 'l', lwd = 2)
  mtext(text =  'FFT-WC ESS', side = 3, line = 0.8, cex = 0.8)
  mtext(text =  paste("Running Time (s) = ", round(post_bound_WC$time[3],2), "and dimension N = ",N), side = 3, line = 0.1, cex = 0.8)
  legend(0.55, 1,
         c("true function", "MAP", "mAP"),
         col = c('black', 'red', 'blue'), 
         text.col = "black", lty = c(1, 2, 4),
         lwd = c(2, 2, 2), text.font = 1, box.lty = 0, cex = 0.8,
         bg='transparent')
}
####################################################









####################################################
#### comparison of behavior of HMC and LS-ESS ######
####################################################
######## comparison HMC vs LS-ESS ##################
if (comp_mon_HMC_LS_WC == 'yes') {
  ### monotone functions
  f <- function(x, idx) {
    if (idx == 1)
      5*x^2
    else if (idx == 2)
      # 3/(1+exp(-10*x+2.1))
      - 2/(1 + exp(-12 * x + 3))
    else if (idx == 3)
      ifelse (x >= 0.6 & x <= 1, (5*x-3)^3, 0) # flat monotone fct
    else if (idx == 4) {
      sum <- 0
      for (l in 1 : 100) {
        sum <- sum+(l^(-1.7)*sin(l)*cos(pi*(l-0.5)*(1-x)))
      }
      return(sqrt(2)*sum)
    }
  }
  ## synthetic data
  idx <- 2 # choice of a monotone function
  set.seed(12345)
  n <- 100
  xtr <- lhsDesign(n = n, dimension = 1, seed = 12345)$design
  # runif(n, min = 0, max = 1)
  y.true <- sapply(xtr, FUN = function(x) f(x, idx = idx)) 
  # sigN <- 0.2 * diff(range(y.true))
  sigN <- 0.5 # noise sd
  ytr <- y.true + rnorm(n, 0, sd = sigN)
  
  N1 <- 9#10# # size of the 1st subdomain
  M <- 3#10# # nb of subdomains
  N <- N1*M # nb of knot points
  eta <- 100 # smooth approximate parameter
  nsim <- 5000 # nb of mcmc retained samples
  brn <- 1000 # burn in
  thin <- 1
  nu <- 1.5 # smoothness kernel parameter
  l <- round(l_est(nu, c(0, 1), 0.05), 2) # length-scale
  tol <- 1e-9
  ## Illustration LS-ESS
  post_mon_LS <- linCGP.ESS(y = ytr, x = xtr, N1 = N1, M = M, nu = nu, l = l, eta = eta, nsim = nsim, burn.in = brn, thin = thin,
                            sig.in = sigN^2, tau.in = 1, constrType = 'decreasing',
                            prior = 'Fast.LS', return.plot = T, tol = 0, sseed = 12345)
  x <- seq(from = 0, to = 1, length = 100)
  lines(x, f(x, idx = idx), type = 'l', lwd = 2)
  mtext(text =  'Large-scale ESS', side = 3, line = 0.8, cex = 0.8)
  mtext(text =  paste("Running Time (s) = ", round(post_mon_LS$time[3],2), "and dimension N =",N), side = 3, line = 0.1, cex = 0.8)
  # mtext(text =  paste("RMSE mAP = ", round(rmse_mAP_LS,2), "and RMSE MAP = ", round(rmse_MAP_LS,2)), side = 3, line = 0.1, cex = 0.8)
  legend(0.5, 1,
         c("true function", "MAP", "mAP"),
         col = c('black', 'red', 'blue'), 
         text.col = "black", lty = c(1, 2, 4),
         lwd = c(2, 2, 2), text.font = 1, box.lty = 0, cex = 0.8,
         bg = 'transparent')
  ## Illustration WC.ESS
  post_mon_WC <- linCGP.WC.ESS(y = ytr, x = xtr, N = N, nu = nu, l = l, eta = eta, nsim = nsim, burn.in = brn, thin = thin,
                               sig.in = sigN^2, tau.in = 1, constrType = 'decreasing',
                               return.plot = T, sseed = 12345, tol = 0)
  x <- seq(from = 0, to = 1, length = 100)
  lines(x, f(x, idx = idx), type = 'l', lwd = 2)
  mtext(text =  'FFT-WC ESS', side = 3, line = 0.8, cex = 0.8)
  mtext(text =  paste("Running Time (s) = ", round(post_mon_WC$time[3],2), "and dimension N = ",N), side = 3, line = 0.1, cex = 0.8)
  legend(0.5, 1,
         c("true function", "MAP", "mAP"),
         col = c('black', 'red', 'blue'), 
         text.col = "black", lty = c(1, 2, 4),
         lwd = c(2, 2, 2), text.font = 1, box.lty = 0, cex = 0.8,
         bg = 'transparent')
  
  ## Illustration HMC
  post_mon_HMC <- linCGP.HMC(y = ytr, x = xtr, N = N, nu = nu, l = l, nsim = nsim,
                             burn.in = brn, thin = thin, sig.in = sigN^2, tau.in = 1,
                             return.plot = T, sseed = 12345, tol = tol, constrType = 'decreasing')
  x <- seq(from = 0, to = 1, length = 100)
  lines(x, f(x, idx = idx), type = 'l', lwd = 2)
  mtext(text =  'Hamiltonian Monte Carlo', side = 3, line = 0.8, cex = 0.8)
  mtext(text =  paste("Running Time (s) = ", round(post_mon_HMC$time[3], 2), "and dimension N = ", N), side = 3, line = 0.1, cex = 0.8)
  legend(0.5, 1,
         c("true function", "MAP", "mAP"),
         col = c('black', 'red', 'blue'), 
         text.col = "black", lty=c(1, 2, 4),
         lwd = c(2, 2, 2), text.font = 1, box.lty = 0, cex = 0.8,
         bg = 'transparent')
} 
##################################################





## end
