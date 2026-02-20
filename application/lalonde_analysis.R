############################################################
# LaLonde Job Training Program Analysis
# RCTE, ATO, ATE, ATT
############################################################

############################
# 1. Packages
############################

packages <- c("dplyr","ranger")

installed <- rownames(installed.packages())
for(p in packages){
  if(!(p %in% installed)) install.packages(p)
}

library(dplyr)
library(ranger)

set.seed(123)

############################
# 2. Load Data
############################

trt <- read.table(
  "https://users.nber.org/~rdehejia/data/nswre74_treated.txt")

colnames(trt) <- c("treat","age","educ","black","hisp",
                   "married","nodegr","re74","re75","re78")

ctrl <- read.table(
  "https://users.nber.org/~rdehejia/data/cps_controls.txt")

colnames(ctrl) <- colnames(trt)

lalonde <- rbind(trt, ctrl)

df <- lalonde %>%
  mutate(T = treat,
         Y = re78)

############################
# 3. Nuisance Models
############################

# Propensity score
ps_model <- glm(
  T ~ age + educ + black + hisp + married +
    nodegr + re74 + re75,
  family = binomial(),
  data = df
)

df$e_hat <- predict(ps_model, type="response")

# Outcome models
form_y <- Y ~ age + educ + black + hisp +
  married + nodegr + re74 + re75

mu0_mod <- ranger(form_y,
                  data = df %>% filter(T==0),
                  num.trees = 200,
                  seed = 123)

mu1_mod <- ranger(form_y,
                  data = df %>% filter(T==1),
                  num.trees = 200,
                  seed = 123)

df$mu0_hat <- predict(mu0_mod, data=df)$predictions
df$mu1_hat <- predict(mu1_mod, data=df)$predictions

############################
# 4. Point Estimators
############################

# ---------- RCTE ----------
df$w_rcte <- 1 - df$mu0_hat

tau_rcte_pi <- mean(df$w_rcte *
                    (df$mu1_hat - df$mu0_hat))

psi_rcte <- df$w_rcte * (
  (df$T/df$e_hat)*(df$Y - df$mu1_hat) -
    ((1-df$T)/(1-df$e_hat))*(df$Y - df$mu0_hat) +
    (df$mu1_hat - df$mu0_hat)
)

tau_rcte_dr <- mean(psi_rcte) /
  mean(df$w_rcte)

# ---------- OWATE ----------
df$w_ow <- df$e_hat*(1-df$e_hat)

tau_owate_pi <- mean(df$w_ow *
                     (df$mu1_hat - df$mu0_hat)) /
  mean(df$w_ow)

psi_ow <- df$w_ow * (
  (df$T/df$e_hat)*(df$Y - df$mu1_hat) -
    ((1-df$T)/(1-df$e_hat))*(df$Y - df$mu0_hat) +
    (df$mu1_hat - df$mu0_hat)
)

tau_owate_dr <- mean(psi_ow) /
  mean(df$w_ow)

# ---------- ATE (AIPW) ----------
aipw <- (df$T/df$e_hat)*(df$Y - df$mu1_hat) -
  ((1-df$T)/(1-df$e_hat))*(df$Y - df$mu0_hat) +
  (df$mu1_hat - df$mu0_hat)

tau_ate_dr <- mean(aipw)

# ---------- ATT ----------
att_i <- df$T*(df$Y - df$mu0_hat)/mean(df$T) +
  ((df$e_hat*(1-df$T))/(1-df$e_hat)) *
  (df$mu1_hat - df$mu0_hat)/mean(df$T)

tau_att <- mean(att_i)

############################
# 5. Print Results
############################

cat("====================================\n")
cat("LaLonde Estimates\n")
cat("====================================\n")

cat("RCTE (Plug-in)  =", tau_rcte_pi, "\n")
cat("RCTE (DR)       =", tau_rcte_dr, "\n")
cat("OWATE (Plug-in) =", tau_owate_pi, "\n")
cat("OWATE (DR)      =", tau_owate_dr, "\n")
cat("ATE (DR)        =", tau_ate_dr, "\n")
cat("ATT (Augmented) =", tau_att, "\n")

############################################################
# 6. Bootstrap (Optional)
############################################################

B <- 1000
set.seed(123)

res_boot <- replicate(B, {

  idx <- sample(nrow(df), replace=TRUE)
  d2  <- df[idx, ]

  ps2 <- glm(T ~ age + educ + black + hisp +
               married + nodegr + re74 + re75,
             family=binomial(),
             data=d2)

  e2 <- predict(ps2, type="response")

  mu0m <- ranger(form_y,
                 data=d2 %>% filter(T==0),
                 num.trees=200)

  mu1m <- ranger(form_y,
                 data=d2 %>% filter(T==1),
                 num.trees=200)

  mu0h <- predict(mu0m, data=d2)$predictions
  mu1h <- predict(mu1m, data=d2)$predictions

  w_r <- 1 - mu0h
  tau_rcte_dr_b <- mean(
    w_r * (
      (d2$T/e2)*(d2$Y - mu1h) -
        ((1-d2$T)/(1-e2))*(d2$Y - mu0h) +
        (mu1h - mu0h)
    )
  ) / mean(w_r)

  tau_rcte_dr_b

}, simplify="matrix")

cat("Bootstrap SE (RCTE DR) =",
    sd(res_boot), "\n")
