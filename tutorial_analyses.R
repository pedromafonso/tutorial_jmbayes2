# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Title: R code accompanying "Extended Joint Models for Longitudinal and       #
#        Time-to-Event Data: A Tutorial"                                       #
# Authors: Pedro Miranda-Afonso and Dimitris Rizopoulos                        #
# Contact: p.mirandaafonso@erasmusmc.nl                                        #
# File: Analyses and code listings presented in the tutorial                   #
# Repository: https://github.com/pedromafonso/tutorial_jmbayes2                #
# Requirements: R and the required R packages, including JMbayes2.             #
# JMbayes2 can be installed from CRAN using install.packages("JMbayes2").      #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Load data ====================================================================

long <- saveRDS(, "Data/long.rds")
surv <- saveRDS("Data/surv.rds")

long_cr <- saveRDS("Data/long_cr.rds")
surv_cr0 <- saveRDS("Data/surv_cr0.rds")

long_ms <- saveRDS("Data/long_ms.rds")
surv_ms <- saveRDS("Data/surv_ms.rds")

long_rc <- saveRDS("Data/long_rc.rds")
surv_rc <- saveRDS("Data/surv_rc.rds")

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
