# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Title: R code accompanying "Extended Joint Models for Longitudinal and       #
#        Time-to-Event Data: A Tutorial"                                       #
# Authors: Pedro Miranda-Afonso and Dimitris Rizopoulos                        #
# Contact: p.mirandaafonso@erasmusmc.nl                                        #
# File: Functions used for data simulation                                     #
# Repository: https://github.com/pedromafonso/tutorial_jmbayes2                #
# Requirements: R and the required R packages, including JMbayes2.             #
# JMbayes2 can be installed from CRAN using install.packages("JMbayes2").      #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Functions ====================================================================
invS <- function(t, i, u_i, b1_i, b2_i, f_i,
                 h0, gammas, alpha1, alpha2, alphaf, betas1, betas2, 
                 long, surv,
                 f_FE1, f_FE2, f_RE1, f_RE2, f_S,
                 tstart = 0, inv_link1 = NULL, inv_link2 = NULL) {
  
  # Baseline covariates
  W <- model.matrix(f_S, data = surv)
  W <- W[, names(gammas), drop = FALSE]
  eta_S <- as.vector(W[i, , drop = FALSE] %*% gammas)
  h <- function(s) {
    time_s <- s + tstart
    data_i <- long[rep(which(long$id == i)[1], length(s)), , drop = FALSE]
    data_i$time <- time_s
    # Longitudinal outcome 1
    X1 <- model.matrix(f_FE1, data = data_i)
    X1 <- X1[, names(betas1), drop = FALSE]
    Z1 <- model.matrix(f_RE1, data = data_i)
    b1 <- matrix(b1_i, nrow = length(s), ncol = length(b1_i), byrow = TRUE)
    eta1 <- as.vector(X1 %*% betas1 + rowSums(Z1 * b1))
    mu1 <- if(is.null(inv_link1)) eta1 else inv_link1(eta1) 
    # Longitudinal outcome 2
    X2 <- model.matrix(f_FE2, data = data_i)
    X2 <- X2[, names(betas2), drop = FALSE]
    Z2 <- model.matrix(f_RE2, data = data_i)
    b2 <- matrix(b2_i, nrow = length(s), ncol = length(b2_i), byrow = TRUE)
    eta2 <- as.vector(X2 %*% betas2 + rowSums(Z2 * b2))
    mu2 <- if(is.null(inv_link2)) eta2 else inv_link2(eta2) 
    # Event hazard
    h0 * exp(eta_S + alpha1 * mu1 + alpha2 * mu2 + alphaf * f_i)
  }
  integrate(h, lower = 0, upper = t)$value + log(u_i)
}

gen_mu <- function(data, betas, b, f_FE, f_RE, inv_link = NULL) {
  X <- model.matrix(f_FE, data = data)
  X <- X[, names(betas), drop = FALSE]
  Z <- model.matrix(f_RE, data = data)
  eta <- as.vector(X %*% betas + rowSums(Z * b[data$id, , drop = FALSE]))
  mu <- if (is.null(inv_link)) eta else inv_link(eta)
  mu
}