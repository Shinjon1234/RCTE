############################################################
# 401(k) Eligibility and Net Financial Assets
# RCTE, ATO, ATE, ATT
############################################################

############################
# 1. Packages
############################

packages <- c("DoubleML","dplyr","ranger")

installed <- rownames(installed.packages())
for(p in packages){
  if(!(p %in% installed)) install.packages(p)
}

library(DoubleML)
library(dplyr)
library(ranger)

set.seed(123)

############################################################
# 2. Load 401(k) Data
############################################################

dml_data <- fetch_401k(
  return_type = "DoubleMLData",
  polynomial_features = FALSE,
  instrument = FALSE
)

df <- dml_data$data

df <- df %>%
  mutate(T = e401,
         Y = net_tfa)

############################################################
# 3. Nuisance Estimation
############################################################

# Propensity score (logit)
ps_model <- glm(
  T ~ age + inc + educ + fsize + marr +
    twoearn + db + pira + hown,
  family = binomial(),
  data = df
)

df$e_hat <- predict(ps_model, type = "response")

# Outcome models (Random Forest)
form_y <- Y ~ age + inc + educ + fsize +
  marr + twoearn + db + pira + hown

mu0_mod <- ranger(
  form_y,
  data = df %>% filter(T == 0),
  num.trees = 500,
  seed = 123
)

mu1_mod <- ranger(
  form_y,
  data = df %>% filter(T == 1),
  num.trees = 500,
  seed = 123
)

df$mu0_hat <- predict(mu0_mod, data = df)$predictions
df$mu1_hat <- predict(mu1_mod, data = df)$predictions

############################################################
# 4. Point Estimators
############################################################

# ---------- RCTE ----------
df$w_rcte <- 1 - df$mu0_hat

tau_rcte_pi <- mean(df$w_rcte *
                    (df$mu1_hat - df$mu0_hat)) /
  mean(df$w_rcte)

psi_rcte <- df$w_rcte * (
  (df$T / df$e_hat) * (df$Y - df$mu1_hat) -
    ((1 - df$T) / (1 - df$e_hat)) * (df$Y - df$mu0_hat) +
    (df$mu1_hat - df$mu0_hat)
)

tau_rcte_dr <- mean(psi_rcte) /
  mean(df$w_rcte)

# ---------- OWATE ----------
df$w_ow <- df$e_hat * (1 - df$e_hat)

tau_owate_pi <- mean(df$w_ow *
                     (df$mu1_hat - df$mu0_hat)) /
  mean(df$w_ow)

psi_ow <- df$w_ow * (
  (df$T / df$e_hat) * (df$Y - df$mu1_hat) -
    ((1 - df$T) / (1 - df$e_hat)) * (df$Y - df$mu0_hat) +
    (df$mu1_hat - df$mu0_hat)
)

tau_owate_dr <- mean(psi_ow) /
  mean(df$w_ow)

# ---------- ATE (AIPW) ----------
aipw_i <- (df$T / df$e_hat) * (df$Y - df$mu1_hat) -
  ((1 - df$T) / (1 - df$e_hat)) * (df$Y - df$mu0_hat) +
  (df$mu1_hat - df$mu0_hat)

tau_ate_dr <- mean(aipw_i)

# ---------- ATT (Augmented) ----------
att_i <- df$T * (df$Y - df$mu0_hat) / mean(df$T) +
  ((df$e_hat * (1 - df$T)) / (1 - df$e_hat)) *
  (df$mu1_hat - df$mu0_hat) / mean(df$T)

tau_att_aug <- mean(att_i)

############################################################
# 5. Print Results
############################################################

cat("============================================\n")
cat("Effect of 401(k) Eligibility on Net Assets\n")
cat("============================================\n")

cat("RCTE (PI)   =", tau_rcte_pi, "\n")
cat("RCTE (DR)   =", tau_rcte_dr, "\n")
cat("OWATE (PI)  =", tau_owate_pi, "\n")
cat("OWATE (DR)  =", tau_owate_dr, "\n")
cat("ATE (DR)    =", tau_ate_dr, "\n")
cat("ATT (Aug)   =", tau_att_aug, "\n")

############################################################
# 6. Bootstrap Inference
############################################################

B <- 500
set.seed(123)

res_boot <- replicate(B, {

  idx <- sample(nrow(df), replace = TRUE)
  d2  <- df[idx, ]

  ps2 <- glm(
    T ~ age + inc + educ + fsize + marr +
      twoearn + db + pira + hown,
    family = binomial(),
    data = d2
  )

  e2 <- predict(ps2, type = "response")

  mu0m <- ranger(form_y,
                 data = d2 %>% filter(T == 0),
                 num.trees = 300)

  mu1m <- ranger(form_y,
                 data = d2 %>% filter(T == 1),
                 num.trees = 300)

  mu0h <- predict(mu0m, data = d2)$predictions
  mu1h <- predict(mu1m, data = d2)$predictions

  # RCTE
  w_rcte <- 1 - mu0h
  tau_rcte_dr_b <- mean(
    w_rcte * (
      (d2$T / e2) * (d2$Y - mu1h) -
        ((1 - d2$T) / (1 - e2)) * (d2$Y - mu0h) +
        (mu1h - mu0h)
    )
  ) / mean(w_rcte)

  tau_rcte_dr_b

}, simplify = "matrix")

cat("Bootstrap SE (RCTE DR) =",
    sd(res_boot), "\n")
