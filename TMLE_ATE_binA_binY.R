## ---------------------------
##
## Script name: TMLE_ATE_binA_binY
##
## Purpose of script: Obtain TMLE estimate for ATE of a binary treatment and binary outcome
##
## Author: Kat Hoffman
##
## Date Created: 2020-03-20
##
## ---------------------------

library(tidyverse)
library(SuperLearner)
set.seed(7)

# Superlearner functions from Ivan
source("sl.r")

# using data generating code from Miguel Angel Luque Fernandez' tutorial
generate_data <- function(n){
    w1 <- rbinom(n, size=1, prob=0.5)
    w2 <- rbinom(n, size=1, prob=0.65)
    w3 <- round(runif(n, min=0, max=4), digits=3)
    w4 <- round(runif(n, min=0, max=5), digits=3)
    A  <- rbinom(n, size=1, prob= plogis(-0.4 + 0.2*w2 + 0.15*w3 + 0.2*w4 + 0.15*w2*w4))
    # counterfactual
    Y_1 <- rbinom(n, size=1, prob= plogis(-1 + 1 -0.1*w1 + 0.3*w2 + 0.25*w3 + 0.2*w4 + 0.15*w2*w4))
    Y_0 <- rbinom(n, size=1, prob= plogis(-1 + 0 -0.1*w1 + 0.3*w2 + 0.25*w3 + 0.2*w4 + 0.15*w2*w4))
    # Observed outcome
    Y <- Y_1*A + Y_0*(1 - A)
    # return data.frame
    tibble(w1, w2, w3, w4, A, Y, Y_1, Y_0)
}

# observations N
n <- 10000

# full data set, including Y0 and Y0
dat_full <- generate_data(n)

# calculate the true psi if we saw both outcomes
true_psi <- mean(dat_full$Y_1 - dat_full$Y_0)

# make a data set with observed data only
dat_obs <- dat_full %>%
  select(-Y_1,-Y_0)
 

# set SL libraries
lib <- c('SL.speedglm', # faster glm
         'SL.glmnet', # lasso
         'SL.ranger', # random forest
         'SL.earth')  #


# estimate Q --------------------------------------------------------------

Y <- dat_obs$Y
X_Y <- dat_obs %>% select(-Y)

Q <- SuperLearner(Y = Y, X = X_Y,
                    family=binomial(),
                    SL.library=lib,
                    # metalearner = NNloglik for binary outcomes, NNLS for continuous
                    method = method.NNloglik,
                    # 5-fold cross validation
                    cvControl = list(V = 5))

# estimate g --------------------------------------------------------------

A <- dat_obs$A
X_A <- dat_obs %>% select(-Y, -A)

g <- SuperLearner(Y = A, X = X_A,
                    family=binomial(),
                    SL.library=lib,
                    # metalearner = NNloglik for binary outcomes, NNLS for continuous
                    method = method.NNloglik,
                    # 5-fold cross validation
                    cvControl = list(V = 5))


# predictions for Q -------------------------------------------------------

Q_A <- predict(Q)$pred

X_Y_A1 <- X_Y %>% mutate(A = 1) 
Q_1 <- predict(Q, newdata = as.data.frame(X_Y_A1))$pred

X_Y_A0 <- X_Y %>% mutate(A = 0)
Q_0 <- predict(Q, newdata = X_Y_A0)$pred


# create clever covariate -------------------------------------------------

H_1 <- 1/predict(g)$pred
H_0 <- -1/(1-predict(g)$pred)

# prep data to fit the clever covariate
dat_cc <-
  dat_obs %>%
  mutate(Q_A = as.vector(Q_A),
         H = case_when(A == 1 ~ H_1,
                       A == 0 ~ H_0)) %>%
  select(Y, A, Q_A, H)

# fit parametric working model
glm_fit <- glm(Y ~ -1 + offset(plogis(Q_A)) + H, data=dat_cc, family=binomial)
# get epsilon
eps <- coef(glm_fit)
# also get H as a vector (for Q_A_update)
H <- dat_cc$H

# update expected outcome estimates (Q star) ------------------------------

# update
Q_1_update <- plogis(qlogis(Q_1) + eps*H_1)
Q_0_update <- plogis(qlogis(Q_0) + eps*H_0)
Q_A_update <- plogis(qlogis(Q_A) + eps*H)

# Estimate ATE ------------------------------------------------------------

# get the ATE TMLE
tmle_psi <- mean(Q_1_update - Q_0_update)

# compare to true parameter
tmle_psi
true_psi 


# Calculate SEs -----------------------------------------------------------

ic <- (dat_obs$Y - Q_A_update) * H + Q_1_update + Q_0_update
tmle_se <- sqrt(var(ic)/nrow(dat_obs))

ci_lo <- tmle_psi - 1.96*tmle_se
ci_hi <- tmle_psi + 1.96*tmle_se

pval <- 2 * (1 - pnorm(abs(tmle_psi / tmle_se)))


# Compare with TMLE package -----------------------------------------------------------

tmle_fit <- tmle::tmle(Y, A, X_A,
                       gbound = .000000001, # trying 
                       Q.SL.library = lib, g.SL.library = lib)
tmle_fit$epsilon # mine are different :(
tmle_fit$estimates$ATE
      