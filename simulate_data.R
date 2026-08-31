rm(list = ls())
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Packages =====================================================================

library("survival")
library("MASS")
library("nlme")
library("GLMMadaptive")
# remotes::install_github("drizopoulos/jmbayes2")
library("JMbayes2")

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

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Simulation settings ==========================================================

seed <- 2026
set.seed(seed)
n <- 500L # number of individuals
n_i <- 18L
t_max <- 10

long0 <- expand.grid(time = seq(from = 0, to = t_max, length.out = n_i),
                     id = seq_len(n))
long0 <- round(long0[, 2:1], 2)

surv0 <- data.frame(id = seq_len(n))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Parameters ===================================================================
# Parameter values are chosen for tutorial illustration and are not intended to
# reproduce the natural history of a specific cystic fibrosis cohort.

## Longitudinal submodels ======================================================
# The longitudinal parameters are chosen to generate clinically plausible
# trajectories while retaining sufficient between-subject heterogeneity for
# illustration. They are loosely informed by patterns observed in CF registry
# data but are not calibrated to reproduce a particular cohort.

### lf
# The continuous marker represents lung function on an approximately scaled
# ppFEV1 scale. Lung function is relatively preserved at baseline and declines
# gradually during follow-up. Females are assigned a slightly lower average
# level and a modestly faster decline. Later diagnosis is associated with poorer 
# lung function, reflecting a longer period during which the underlying 
# congenital disease may have progressed without appropriate CF-specific 
# monitoring and treatment. Random intercepts and slopes allow some individual 
# variation around this average trajectory.

f_FE_lf <- ~ time + sex + time:sex + ageD
f_RE_lf <- ~ time

betas_lf <- c("(Intercept)" =  0.850,
              "time"        = -0.030,
              "sex"         = -0.030,
              "time:sex"    = -0.003,
              "ageD"        = -0.015)

sd_b_lf <- c("(Intercept)" = 0.180,
             "time"        = 0.015)

sigma_lf <- 0.060

### pa
# The binary marker represents the presence of Pseudomonas aeruginosa at a
# clinical assessment. Positivity is relatively uncommon at the beginning of
# follow-up but becomes progressively more likely over time. Female sex is
# associated with a small increase in the probability of positivity. Later 
# diagnosis is associated with a higher probability of positivity, reflecting 
# the potential consequences of a longer period without targeted management of 
# CF respiratory disease. Subject-specific random intercepts and slopes allow
# clinically plausible heterogeneity in baseline probability and its evolution.

f_FE_pa <- ~ time + sex + ageD
f_RE_pa <- ~ time

betas_pa <- c("(Intercept)" = -2.150,
              "time"        =  0.170,
              "sex"         =  0.150,
              "ageD"        =  0.060)

sd_b_pa <- c("(Intercept)" = 0.650,
             "time"        = 0.100)

### Random effects
# For simplicity, the random effects used to generate the two longitudinal
# outcomes are assumed independent, both within and between outcomes. Moderate
# random-intercept and random-slope variances are used to generate heterogeneous
# individual trajectories while avoiding excessively extreme simulated values.
# The binary-outcome random-slope variance is kept smaller than its
# random-intercept variance to retain trajectory heterogeneity without
# generating unnecessarily extreme subject-specific trends.

D <- diag(c(sd_b_lf, sd_b_pa)^2)
b <- mvrnorm(n, mu = rep(0, ncol(D)), Sigma = D)
b_lf <- b[, seq_along(sd_b_lf), drop = FALSE]
b_pa <- b[, length(sd_b_lf) + seq_along(sd_b_pa), drop = FALSE]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## Event-time submodels ========================================================

### alv -> dth
# The risk of death increases as lung function deteriorates and as the predicted
# probability of Pseudomonas aeruginosa positivity increases. Female sex is
# assigned a small increase in risk. The prognostic contribution of later
# diagnosis operates both through the longitudinal markers and through a small
# direct hazard effect. A subject-specific frailty contribution represents
# residual mortality heterogeneity not explained by the observed markers.
f_alv_dth <- ~ sex + ageD

h0_alv_dth <- 0.062

gammas_alv_dth <- c("sex"  =  0.100, 
                    "ageD" =  0.020)

alphas_alv_dth <- c("lf" = -0.900, 
                    "pa" =  0.800)

alphaf_alv_dth <- 0.100

### alv -> tx
# Lung transplantation is more likely among individuals with poorer lung function
# and a higher predicted probability of Pseudomonas aeruginosa positivity. These
# associations are slightly stronger than for death, reflecting the importance of
# advanced pulmonary disease in transplant eligibility. The baseline
# hazard is set slightly below that for death; the stronger marker associations
# are intended to result in transplantation occurring slightly more frequently.
# Age at diagnosis has a small direct effect, while no additional frailty
# contribution is included for transplantation.

f_alv_tx <- ~ sex + ageD

h0_alv_tx <- 0.060

gammas_alv_tx <- c("sex"  =  0.100, 
                   "ageD" =  0.020)

alphas_alv_tx <- c("lf" = -1.000, 
                   "pa" =  1.300)

alphaf_alv_tx <- 0.000

### tx -> dth
# Individuals remain at risk of death after transplantation. For illustration,
# poorer lung function and a higher probability of Pseudomonas aeruginosa
# positivity are assumed to remain markers of poorer prognosis after
# transplantation, with attenuated associations. The same underlying
# longitudinal trajectories are allowed to continue after transplantation as a
# simplifying modeling assumption rather than as a literal representation of
# post-transplant physiology. A comparatively higher baseline hazard ensures an
# adequate number of post-transplant deaths. Age at diagnosis has a small direct
# effect; the frailty has no additional direct effect on this transition.

f_tx_dth <- ~ sex + ageD

h0_tx_dth <- 0.150

gammas_tx_dth <- c("sex"  =  0.050, 
                   "ageD" =  0.020)

alphas_tx_dth <- c("lf" = -0.800, 
                   "pa" =  0.600)

alphaf_tx_dth <- 0

# ### pex
f_pex <- ~ sex + ageD

h0_pex <- 0.25

gammas_pex <- c("sex"  = 0.150,
                "ageD" = 0.020)

alphas_pex <- c("lf" = -0.800,
                "pa" =  0.000)

alphaf_pex <- 1

dur_pex <- 0

### frailty
# A normally distributed subject-specific frailty represents unmeasured
# characteristics affecting event risk independently of the longitudinal random
# effects. The frailty enters the recurrent pulmonary-exacerbation intensity
# directly, with coefficient fixed to one, and is weakly associated with the
# terminal death process through alphaf_alv_dth. This induces dependence between
# recurrent exacerbations and mortality beyond that explained by the observed
# longitudinal markers.
sigma_f <- 0.5
f <- rnorm(n, mean = 0, sd = sigma_f)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Generate data ================================================================
## Baseline covariates =========================================================

surv0$sex <- rbinom(n, 1, 0.50) # 0 = male, 1 = female
long0$sex <- surv0$sex[long0$id]
sh_age <- 2
sc_age <- 1.2
#curve(dgamma(x, shape = sh_age, scale = sc_age))
surv0$ageD <- pmin(rgamma(n, shape = sh_age, scale = sc_age), t_max)
long0$ageD <- surv0$ageD[long0$id]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## Longitudinal outcomes =======================================================

### lf
mu_lf <- gen_mu(data = long0,
                betas = betas_lf,
                b = b_lf,
                f_FE = f_FE_lf,
                f_RE = f_RE_lf)

long0$lf <- rnorm(nrow(long0), mean = mu_lf, sd = sigma_lf)

### pa
prob_pa <- gen_mu(data = long0,
                  betas = betas_pa,
                  b = b_pa,
                  f_FE = f_FE_pa,
                  f_RE = f_RE_pa,
                  inv_link = plogis)

long0$pa <- rbinom(nrow(long0),
                   size = 1,
                   prob = prob_pa)

# ### Check simulated longitudinal outcomes
# #### Continuous outcome
# m_lf <- lme(fixed = update(f_FE_lf, lf ~ .),
#             random = as.formula(paste(deparse(f_RE_lf), "| id")),
#             data = long0)
# round(cbind(fixef(m_lf), betas_lf), 2) # FE
# round(cbind(sigma(m_lf), sigma_lf), 2) # Residual sd
# round(as.matrix(getVarCov(m_lf, type = "random.effects")), 2) # D
# ind_lf <- seq_along(sd_b_lf)
# round(D[ind_lf, ind_lf], 2)
# remove(ind_lf)
# #### Binary outcome
# m_pa <- mixed_model(fixed = update(f_FE_pa, pa ~ .),
#                     random = as.formula(paste(deparse(f_RE_pa), "| id")),
#                     family = binomial(link = "logit"),
#                     data = long0)
# round(cbind(fixef(m_pa), betas_pa), 2) # FE
# round(m_pa$D, 2) # D
# ind_pa <- length(sd_b_lf) + seq_along(sd_b_pa)
# round(D[ind_pa, ind_pa], 2)
# remove(ind_pa)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## Event times =================================================================

time_alv_dth <- rep(Inf, n)
time_alv_tx  <- rep(Inf, n)
time_tx_dth  <- rep(Inf, n)
time_pex     <- vector("list", n)

for (i in seq_len(n)) {
  
  max_gap <- t_max
  tstart_i <- 0
  
  # alv -> dth
  u <- runif(1)
  root_tmax <- invS(t = max_gap,
                    i = i,
                    u_i = u,
                    b1_i = b_lf[i, ],
                    b2_i = b_pa[i, ],
                    f_i = f[i],
                    h0 = h0_alv_dth,
                    gammas = gammas_alv_dth,
                    alpha1 = alphas_alv_dth["lf"],
                    alpha2 = alphas_alv_dth["pa"],
                    alphaf = alphaf_alv_dth,
                    betas1 = betas_lf,
                    betas2 = betas_pa,
                    long = long0,
                    surv = surv0,
                    f_FE1 = f_FE_lf,
                    f_FE2 = f_FE_pa,
                    f_RE1 = f_RE_lf,
                    f_RE2 = f_RE_pa,
                    f_S = f_alv_dth,
                    tstart = tstart_i,
                    inv_link1 = NULL,
                    inv_link2 = plogis)
  
  if (root_tmax > 0) {
    
    time_alv_dth[i] <- uniroot(invS,
                               interval = c(0, max_gap),
                               i = i,
                               u_i = u,
                               b1_i = b_lf[i, ],
                               b2_i = b_pa[i, ],
                               f_i = f[i],
                               h0 = h0_alv_dth,
                               gammas = gammas_alv_dth,
                               alpha1 = alphas_alv_dth["lf"],
                               alpha2 = alphas_alv_dth["pa"],
                               alphaf = alphaf_alv_dth,
                               betas1 = betas_lf,
                               betas2 = betas_pa,
                               long = long0,
                               surv = surv0,
                               f_FE1 = f_FE_lf,
                               f_FE2 = f_FE_pa,
                               f_RE1 = f_RE_lf,
                               f_RE2 = f_RE_pa,
                               f_S = f_alv_dth,
                               tstart = tstart_i,
                               inv_link1 = NULL,
                               inv_link2 = plogis)$root
  }
  
  # alv -> tx
  u <- runif(1)
  root_tmax <- invS(t = max_gap,
                    i = i,
                    u_i = u,
                    b1_i = b_lf[i, ],
                    b2_i = b_pa[i, ],
                    f_i = f[i],
                    h0 = h0_alv_tx,
                    gammas = gammas_alv_tx,
                    alpha1 = alphas_alv_tx["lf"],
                    alpha2 = alphas_alv_tx["pa"],
                    alphaf = alphaf_alv_tx,
                    betas1 = betas_lf,
                    betas2 = betas_pa,
                    long = long0,
                    surv = surv0,
                    f_FE1 = f_FE_lf,
                    f_FE2 = f_FE_pa,
                    f_RE1 = f_RE_lf,
                    f_RE2 = f_RE_pa,
                    f_S = f_alv_tx,
                    tstart = tstart_i,
                    inv_link1 = NULL,
                    inv_link2 = plogis)
  
  if (root_tmax > 0) {
    
    time_alv_tx[i] <- uniroot(invS,
                              interval = c(0, max_gap),
                              i = i,
                              u_i = u,
                              b1_i = b_lf[i, ],
                              b2_i = b_pa[i, ],
                              f_i = f[i],
                              h0 = h0_alv_tx,
                              gammas = gammas_alv_tx,
                              alpha1 = alphas_alv_tx["lf"],
                              alpha2 = alphas_alv_tx["pa"],
                              alphaf = alphaf_alv_tx,
                              betas1 = betas_lf,
                              betas2 = betas_pa,
                              long = long0,
                              surv = surv0,
                              f_FE1 = f_FE_lf,
                              f_FE2 = f_FE_pa,
                              f_RE1 = f_RE_lf,
                              f_RE2 = f_RE_pa,
                              f_S = f_alv_tx,
                              tstart = tstart_i,
                              inv_link1 = NULL,
                              inv_link2 = plogis)$root
  }
  
  # tx -> dth
  if (time_alv_tx[i] < time_alv_dth[i]) {
    
    tstart_i <- time_alv_tx[i]
    max_gap <- t_max - tstart_i
    
    u <- runif(1)
    root_tmax <- invS(t = max_gap,
                      i = i,
                      u_i = u,
                      b1_i = b_lf[i, ],
                      b2_i = b_pa[i, ],
                      f_i = f[i],
                      h0 = h0_tx_dth,
                      gammas = gammas_tx_dth,
                      alpha1 = alphas_tx_dth["lf"],
                      alpha2 = alphas_tx_dth["pa"],
                      alphaf = alphaf_tx_dth,
                      betas1 = betas_lf,
                      betas2 = betas_pa,
                      long = long0,
                      surv = surv0,
                      f_FE1 = f_FE_lf,
                      f_FE2 = f_FE_pa,
                      f_RE1 = f_RE_lf,
                      f_RE2 = f_RE_pa,
                      f_S = f_tx_dth,
                      tstart = tstart_i,
                      inv_link1 = NULL,
                      inv_link2 = plogis)
    
    if (root_tmax > 0) {
      
      gap <- uniroot(invS,
                     interval = c(0, max_gap),
                     i = i,
                     u_i = u,
                     b1_i = b_lf[i, ],
                     b2_i = b_pa[i, ],
                     f_i = f[i],
                     h0 = h0_tx_dth,
                     gammas = gammas_tx_dth,
                     alpha1 = alphas_tx_dth["lf"],
                     alpha2 = alphas_tx_dth["pa"],
                     alphaf = alphaf_tx_dth,
                     betas1 = betas_lf,
                     betas2 = betas_pa,
                     long = long0,
                     surv = surv0,
                     f_FE1 = f_FE_lf,
                     f_FE2 = f_FE_pa,
                     f_RE1 = f_RE_lf,
                     f_RE2 = f_RE_pa,
                     f_S = f_tx_dth,
                     tstart = tstart_i,
                     inv_link1 = NULL,
                     inv_link2 = plogis)$root
      
      time_tx_dth[i] <- tstart_i + gap
    }
  }
  
  # pex
  end_i <- min(time_alv_dth[i], time_alv_tx[i], t_max)
  tstart_i <- 0
  
  while (tstart_i < end_i) {
    
    max_gap <- end_i - tstart_i
    u <- runif(1)
    
    root_tmax <- invS(t = max_gap,
                      i = i,
                      u_i = u,
                      b1_i = b_lf[i, ],
                      b2_i = b_pa[i, ],
                      f_i = f[i],
                      h0 = h0_pex,
                      gammas = gammas_pex,
                      alpha1 = alphas_pex["lf"],
                      alpha2 = alphas_pex["pa"],
                      alphaf = alphaf_pex,
                      betas1 = betas_lf,
                      betas2 = betas_pa,
                      long = long0,
                      surv = surv0,
                      f_FE1 = f_FE_lf,
                      f_FE2 = f_FE_pa,
                      f_RE1 = f_RE_lf,
                      f_RE2 = f_RE_pa,
                      f_S = f_pex,
                      tstart = tstart_i,
                      inv_link1 = NULL,
                      inv_link2 = plogis)
    
    if (root_tmax <= 0) break
    
    gap <- uniroot(invS,
                   interval = c(0, max_gap),
                   i = i,
                   u_i = u,
                   b1_i = b_lf[i, ],
                   b2_i = b_pa[i, ],
                   f_i = f[i],
                   h0 = h0_pex,
                   gammas = gammas_pex,
                   alpha1 = alphas_pex["lf"],
                   alpha2 = alphas_pex["pa"],
                   alphaf = alphaf_pex,
                   betas1 = betas_lf,
                   betas2 = betas_pa,
                   long = long0,
                   surv = surv0,
                   f_FE1 = f_FE_lf,
                   f_FE2 = f_FE_pa,
                   f_RE1 = f_RE_lf,
                   f_RE2 = f_RE_pa,
                   f_S = f_pex,
                   tstart = tstart_i,
                   inv_link1 = NULL,
                   inv_link2 = plogis)$root
    
    time_event <- tstart_i + gap
    
    time_pex[[i]] <- c(time_pex[[i]], time_event)
    
    tstart_i <- time_event + dur_pex
    
  }
}

ind_tx <- time_alv_tx < time_alv_dth

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## Datasets ====================================================================
### Basic JM

surv <- data.frame(id = surv0$id,
                   stop = pmin(time_alv_dth, t_max),
                   status = is.finite(time_alv_dth),
                   sex = surv0$sex,
                   ageD = surv0$ageD) 

long <- long0[long0$time <= surv$stop[long0$id], ]

### Competing risks JM
time_comb <- pmin(time_alv_dth, time_alv_tx, t_max)
surv_cr0 <- data.frame(id = surv0$id,
                       stop = time_comb,
                       status = "alv",
                       sex = surv0$sex,
                       ageD = surv0$ageD)
surv_cr0$status[time_alv_tx < time_alv_dth] <- "tx"
surv_cr0$status[time_alv_tx > time_alv_dth] <- "dth"
surv_cr0$status <- factor(surv_cr0$status, levels = c("alv", "tx", "dth"))

long_cr <- long0[long0$time <= surv_cr0$stop[long0$id], ]

### Multistate JM
surv_alv_tx <- data.frame(id = surv0$id,
                          start = 0,
                          stop = time_comb,
                          status = time_alv_tx < time_alv_dth,
                          strat = "alv-tx",
                          sex = surv0$sex,
                          ageD = surv0$ageD)
surv_alv_dth <- data.frame(id = surv0$id,
                           start = 0,
                           stop = time_comb,
                           status = time_alv_dth < time_alv_tx,
                           strat = "alv-dth",
                           sex = surv0$sex,
                           ageD = surv0$ageD)
surv_tx_dth <- data.frame(id = surv0$id[ind_tx],
                          start = time_alv_tx[ind_tx],
                          stop = pmin(time_tx_dth[ind_tx], t_max),
                          status = is.finite(time_tx_dth[ind_tx]),
                          strat = "tx-dth",
                          sex = surv0$sex[ind_tx],
                          ageD = surv0$ageD[ind_tx])

surv_ms <- rbind(surv_alv_tx, surv_alv_dth, surv_tx_dth)
surv_ms$strat <- factor(surv_ms$strat, levels = c("alv-tx", "alv-dth", "tx-dth"))
surv_ms <- surv_ms[order(surv_ms$id, surv_ms$start, surv_ms$strat), ]
rownames(surv_ms) <- NULL

stop_max <- tapply(surv_ms$stop, surv_ms$id, max)
long_ms <- long0[long0$time <= stop_max[long0$id], ]
remove(stop_max)

### Recurrent events JM
n_pex <- lengths(time_pex)

ids_rc <- rep(seq_len(n), n_pex + 1L)

surv_rc <- data.frame(id = ids_rc,
                      sex = surv0$sex[ids_rc],
                      ageD = surv0$ageD[ids_rc])
surv_rc$start <- unlist(Map(function(x) c(0, x), time_pex), use.names = FALSE)
surv_rc$stop  <- unlist(Map(function(x, end) c(x, end), time_pex, surv_cr0$stop), 
                        use.names = FALSE)
surv_rc$status <- as.integer(sequence(n_pex + 1L) <= rep(n_pex, times = n_pex + 1L))
rownames(surv_rc) <- NULL

long_rc <- long0[long0$time <= surv_cr0$stop[long0$id], ]

### Export datasets
saveRDS(long, "long.rds")
saveRDS(surv, "surv.rds")

saveRDS(long_cr,  "long_cr.rds")
saveRDS(surv_cr0, "surv_cr0.rds")

saveRDS(long_ms, "long_ms.rds")
saveRDS(surv_ms, "surv_ms.rds")

saveRDS(long_rc, "long_rc.rds")
saveRDS(surv_rc, "surv_rc.rds")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Data descriptives ============================================================
## Baseline characteristics ====================================================

round(c(female = mean(surv0$sex), 
        mean_ageD = mean(surv0$ageD), 
        sd_ageD = sd(surv0$ageD),
        median_ageD = median(surv0$ageD),
        q25_ageD = quantile(surv0$ageD, 0.25),
        q75_ageD = quantile(surv0$ageD, 0.75)), 2)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## Longitudinal outcomes =======================================================

desc_long <- aggregate(cbind(lf, pa) ~ time, data = long0, FUN = mean)
names(desc_long)[2:3] <- c("mean_lf", "mean_pa")
desc_long$sd_lf <- aggregate(lf ~ time, data = long0, FUN = sd)$lf
desc_long$q025_lf <- aggregate(lf ~ time, data = long0, FUN = quantile, 
                               probs = 0.025)$lf
desc_long$q975_lf <- aggregate(lf ~ time, data = long0, FUN = quantile,
                               probs = 0.975)$lf
desc_long$sd_pa <- aggregate(pa ~ time, data = long0, FUN = sd)$pa
round(desc_long[, c(1, 2, 4:6, 3, 7)], 3)

round(summary(long0$lf), 3)
round(sd(long0$lf), 3)
mean(long0$lf < 0)

round(summary(prob_pa), 3)
round(mean(long0$pa), 3)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## Event-time outcomes =========================================================

### Basic JM
round(c(n = nrow(surv),
        deaths = sum(surv$status),
        death_prop = mean(surv$status),
        median_death_time = median(surv$stop[surv$status]),
        median_followup = median(surv$stop)), 3)

### CR JM: transplantation and death
round(c(finite_alv_tx = mean(is.finite(time_alv_tx)),
        finite_alv_dth = mean(is.finite(time_alv_dth))), 3)

tab_cr <- table(surv_cr0$status)
desc_cr <- data.frame(status = names(tab_cr),
                      n = as.vector(tab_cr),
                      proportion = as.vector(prop.table(tab_cr)))
desc_cr$median_time <- sapply(desc_cr$status, function(x) 
  median(surv_cr0$stop[surv_cr0$status == x]))
desc_cr[, -1] <- round(desc_cr[, -1], 3)
desc_cr

c(tx = sum(surv_cr0$status == "tx"),
  dth = sum(surv_cr0$status == "dth"),
  tx_minus_dth = sum(surv_cr0$status == "tx") - sum(surv_cr0$status == "dth"))

### Mulsitate JM
lev_ms <- levels(surv_ms$strat)
desc_ms <- data.frame(transition = lev_ms)
desc_ms$n_risk <- sapply(lev_ms, function(x) sum(surv_ms$strat == x))
desc_ms$events <- sapply(lev_ms, function(x) sum(surv_ms$status[surv_ms$strat == x]))
desc_ms$person_time <- sapply(lev_ms, function(x) sum(surv_ms$stop[surv_ms$strat == x] - 
                                                        surv_ms$start[surv_ms$strat == x]))
desc_ms$event_prop <- desc_ms$events / desc_ms$n_risk
desc_ms$rate_100py <- 100 * desc_ms$events / desc_ms$person_time
desc_ms$median_event_time <- sapply(lev_ms, function(x) { 
  ind <- surv_ms$strat == x & surv_ms$status == 1
  median(surv_ms$stop[ind])
})
desc_ms$event_prop <- desc_ms$event_prop
desc_ms$person_time <- desc_ms$person_time
desc_ms$rate_100py <- desc_ms$rate_100py
desc_ms$median_event_time <- desc_ms$median_event_time
desc_ms[, -1] <- round(desc_ms[, -1], 3)
desc_ms

#### tx -> dth
ind_tx_dth <- ind_tx & is.finite(time_tx_dth)
round(c(transplanted = sum(ind_tx),
        post_tx_deaths = sum(ind_tx_dth),
        post_tx_death_prop = sum(ind_tx_dth) / sum(ind_tx),
        median_tx_time = median(time_alv_tx[ind_tx]),
        median_post_tx_death_time = median(time_tx_dth[ind_tx_dth]),
        median_tx_to_death = median(time_tx_dth[ind_tx_dth] - time_alv_tx[ind_tx_dth])
), 2)

#### pex

round(c(total_pex = sum(n_pex),
        subjects_with_pex = sum(n_pex > 0),
        prop_with_pex = mean(n_pex > 0),
        mean_pex = mean(n_pex),
        sd_pex = sd(n_pex),
        median_pex = median(n_pex),
        q25_pex = quantile(n_pex, 0.25),
        q75_pex = quantile(n_pex, 0.75),
        q90_pex = quantile(n_pex, 0.90),
        max_pex = max(n_pex)), 2)

table(n_pex)

gap_pex <- with(surv_rc[surv_rc$status == 1, ],
                stop - start)

round(summary(gap_pex), 2)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Listings =====================================================================
## Listing A ===================================================================

fit_long1 <- lme(fixed = lf ~ time * sex + ageD,
                 random =~ time | id,
                 data = long)

fit_surv <- coxph(Surv(stop, status) ~ sex + ageD,
                  data = surv)

fit_jm <- jm(Surv_object = fit_surv,
             Mixed_objects = fit_long1,
             time_var = "time")

summary(fit_jm)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## Listing B ===================================================================

fit_long1_b <- update(fit_long1, fixed = lf ~ time * sex + ns(ageD, 2))

fit_jm_b <- jm(Surv_object = fit_surv,
               Mixed_objects = fit_long1_b,
               time_var = "time")

summary(fit_jm_b)
compare_jm(fit_jm, fit_jm_b, type = "marginal")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## Listing C ===================================================================

fit_long2 <- mixed_model(fixed = pa ~ time + sex + ageD,
                         random =~ time | id,
                         family = binomial(link = "logit"),
                         data = long)

fit_jm2 <- jm(Surv_object = fit_surv,
              Mixed_objects = list(fit_long1, fit_long2),
              time_var = "time")

summary(fit_jm2)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## Listing D ===================================================================

fit_jm3 <- update(fit_jm2,
                  functional_forms =~ value(lf) + slope(lf) + vexpit(value(pa)))

summary(fit_jm3)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## Listing E ===================================================================

fit_jm4 <- update(fit_jm2,
                  functional_forms =~ area(lf, time_window = 0.5) +
                    Delta(lf, time_window = 0.5, standardise = TRUE) +
                    vexpit(value(pa)),
                  n_iter = 7000L, n_burnin = 1000L, n_thin = 2L)
summary(fit_jm4)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## Listing F ===================================================================

surv_rc[1:6, c("id", "start", "stop", "status")]

fit_surv_rc <- coxph(Surv(start, stop, status) ~ sex + ageD,
                     data = surv_rc)

fit_jm_rc <- jm(Surv_object = fit_surv_rc, 
                Mixed_objects = fit_long1, 
                time_var = "time", 
                recurrent = "gap")

summary(fit_jm_rc)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## Listing G ===================================================================

surv_cr0[1:3, c("id", "stop", "status")]

surv_cr <- crisk_setup(surv_cr0,
                       statusVar = "status",
                       censLevel = "alv",
                       nameStrata = "strat")
surv_cr[1:6, c("id", "stop", "status2", "strat")]

fit_long1_cr <- update(fit_long1, data = long_cr)

fit_long2_cr <- update(fit_long2, data = long_cr)

fit_surv_cr <- coxph(Surv(stop, status2) ~ (sex + ageD):strata(strat),
                     data = surv_cr)

fit_jm_cr <- jm(Surv_object = fit_surv_cr,
                Mixed_objects = list(fit_long1_cr, fit_long2_cr),
                time_var = "time",
                functional_forms =~ (value(lf) + vexpit(value(pa))):strat,
                n_iter = 7000L, n_burnin = 1000L, n_thin = 2L)

summary(fit_jm_cr)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## Listing H ===================================================================

surv_ms[1:4, c("id", "start", "stop", "status", "strat")]

fit_long1_ms <- update(fit_long1, data = long_ms)

fit_long2_ms <- update(fit_long2, data = long_ms)

fit_surv_ms <- coxph(Surv(start, stop, status) ~ (sex + ageD):strata(strat),
                     data = surv_ms)

fit_jm_ms <- jm(Surv_object = fit_surv_ms,
                Mixed_objects = list(fit_long1_ms, fit_long2_ms),
                time_var = "time",
                functional_forms =~ (value(lf) + vexpit(value(pa))):strat,
                n_iter = 14000L, n_burnin = 2000L, n_thin = 4L)

summary(fit_jm_ms)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## Listing I ===================================================================

surv_cr0[1:3, c("id", "stop", "status")]

surv_rc[1:6, c("id", "start", "stop", "status")]

surv_comb <- rc_setup(rc_data = surv_rc, trm_data = surv_cr0,
                      idVar = "id", statusVar = "status",
                      startVar = "start", stopVar = "stop",
                      trm_censLevel = "alv",
                      nameStrata = "strat", nameStatus = "status")

surv_comb[1:12, c("id", "start", "stop", "status", "strat")]

fit_surv_comb <- coxph(Surv(start, stop, status) ~ (sex + ageD):strata(strat),
                       data = surv_comb)

fit_jm_comb <- jm(Surv_object = fit_surv_comb, 
                  Mixed_objects = list(fit_long1_cr, fit_long2_cr), 
                  time_var = "time", recurrent = "gap",
                  functional_forms =~ (value(lf) + vexpit(value(pa))):strat,
                  n_iter = 14000L, n_burnin = 2000L, n_thin = 4L)

summary(fit_jm_comb)
