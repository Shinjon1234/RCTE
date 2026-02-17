set.seed(123)

############################
## 1. Reference DGP
############################

generate_data_ref <- function(n,
                              cate_type = c("zero", "binary", "linear", "nonlinear"),
                              ps_type   = c("random", "linear", "interaction", "nonlinear"),
                              p = 10) {
  
  cate_type <- match.arg(cate_type)
  ps_type   <- match.arg(ps_type)
  
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("X", 1:p)
  X <- as.data.frame(X)
  
  mu0 <- with(X, X1 + X5 + X4 * X5)
  
  if (cate_type == "zero") {
    tau <- rep(0, n)
  }
  if (cate_type == "binary") {
    tau <- ifelse(X$X2 >= 0, 2, -1)
  }
  if (cate_type == "linear") {
    tau <- X$X1 + X$X2
  }
  if (cate_type == "nonlinear") {
    tau <- 2 * sin(X$X1 + 0.5 * X$X2 + 0.33 * X$X3) + 
      1.5 * cos(X$X10)
  }
  
  if (ps_type == "random") {
    eX <- rep(0.3, n)
  }
  if (ps_type == "linear") {
    aX <- with(X, X2 + X3 + X5 - X8)
    eX <- pnorm(scale(aX))
  }
  if (ps_type == "interaction") {
    b <- 1 / (1:p)
    aX <- as.vector(as.matrix(X) %*% b + with(X, X2 + X5 + X3 * X8))
    eX <- pnorm(scale(aX))
  }
  if (ps_type == "nonlinear") {
    b <- 1 / (1:p)
    aX <- as.vector(as.matrix(X) %*% b + with(X,
                                              X2 + 1.5 * cos(X4 * X8) + 2 * sin(X5)
    ))
    eX <- pnorm(scale(aX))
  }
  
  A <- rbinom(n, 1, eX)
  U <- rnorm(n)
  
  Y  <- mu0 + tau * A + U
  Y0 <- mu0 + rnorm(n)
  Y1 <- mu0 + tau + rnorm(n)
  
  data.frame(X, A, Y, Y0, Y1, mu0, tau, eX)
}

############################
## 2. Estimators (PI + DR)
############################

estimators_all <- function(dat) {
  X <- dat[, grep("^X", names(dat))]
  A <- dat$A
  Y <- dat$Y
  
  ps_fit <- glm(A ~ ., family = binomial(), data = data.frame(A, X))
  e_hat  <- pmin(pmax(predict(ps_fit, type = "response"), 0.01), 0.99)
  
  mu0_fit <- lm(Y ~ ., data = data.frame(Y, X)[A == 0, ])
  mu1_fit <- lm(Y ~ ., data = data.frame(Y, X)[A == 1, ])
  
  mu0_hat <- predict(mu0_fit, newdata = X)
  mu1_hat <- predict(mu1_fit, newdata = X)
  
  w_rcte_hat <- 1 - mu0_hat
  w_ow_hat   <- pmax(e_hat * (1 - e_hat), 0.01)
  
  ## Plug-in
  ATE_PI   <- mean(mu1_hat - mu0_hat)
  ATT_PI   <- mean(mu1_hat[A == 1] - mu0_hat[A == 1])
  RCTE_PI  <- sum(w_rcte_hat * (mu1_hat - mu0_hat)) / sum(w_rcte_hat)
  OWATE_PI <- sum(w_ow_hat * (mu1_hat - mu0_hat)) / sum(w_ow_hat)
  
  ## DR ATE
  phi_ate <- A * (Y - mu1_hat) / e_hat -
    (1 - A) * (Y - mu0_hat) / (1 - e_hat) +
    (mu1_hat - mu0_hat)
  ATE_DR <- mean(phi_ate)
  
  ## DR ATT (correct EIF)
  pA <- mean(A)
  ATT_DR <- mean(
    A * (Y - mu1_hat) / pA -
      (1 - A) * e_hat * (Y - mu0_hat) / (pA * (1 - e_hat)) +
      A * (mu1_hat - mu0_hat) / pA
  )
  
  ## DR RCTE
  phi_rcte <- w_rcte_hat * (
    A * (Y - mu1_hat) / e_hat -
      (1 - A) * (Y - mu0_hat) / (1 - e_hat) +
      (mu1_hat - mu0_hat)
  )
  RCTE_DR <- sum(phi_rcte) / sum(w_rcte_hat)
  
  ## DR OWATE
  phi_ow <- w_ow_hat * (
    A * (Y - mu1_hat) / e_hat -
      (1 - A) * (Y - mu0_hat) / (1 - e_hat) +
      (mu1_hat - mu0_hat)
  )
  OWATE_DR <- sum(phi_ow) / sum(w_ow_hat)
  
  c(ATE_PI, ATE_DR, ATT_PI, ATT_DR, RCTE_PI, RCTE_DR, OWATE_PI, OWATE_DR)
}

############################
## 3. True Targets
############################

true_targets_all <- function(dat) {
  tau_i <- dat$Y1 - dat$Y0
  eX    <- dat$eX
  mu0   <- dat$mu0
  
  ATE  <- mean(tau_i)
  ATT  <- mean(tau_i[dat$A == 1])
  
  w_rcte <- 1 - mu0
  RCTE   <- mean(w_rcte * tau_i) / mean(w_rcte)
  
  w_ow   <- eX * (1 - eX)
  OWATE  <- mean(w_ow * tau_i) / mean(w_ow)
  
  c(ATE, ATE, ATT, ATT, RCTE, RCTE, OWATE, OWATE)
}

############################
## 4. Monte Carlo
############################

run_mc <- function(R = 300, n = 1000, cate_type, ps_type) {
  est_mat  <- matrix(NA, R, 8)
  true_mat <- matrix(NA, R, 8)
  
  for (r in 1:R) {
    dat <- generate_data_ref(n, cate_type, ps_type)
    est_mat[r, ]  <- estimators_all(dat)
    true_mat[r, ] <- true_targets_all(dat)
  }
  
  colnames(est_mat) <- colnames(true_mat) <-
    c("ATE_PI", "ATE_DR", "ATT_PI", "ATT_DR",
      "RCTE_PI", "RCTE_DR", "OWATE_PI", "OWATE_DR")
  
  list(est = est_mat, true = true_mat)
}

############################
## 5. Summary Table
############################

make_mc_table <- function(est_mat, true_mat) {
  data.frame(
    Estimator = colnames(est_mat),
    True      = round(colMeans(true_mat), 3),
    Mean      = round(colMeans(est_mat), 3),
    Bias      = round(colMeans(est_mat - true_mat), 3),
    Variance  = round(apply(est_mat, 2, var), 3),
    RMSE      = round(sqrt(colMeans((est_mat - true_mat)^2)), 3)
  )
}

############################
## 6. Run All 16 Scenarios
############################

cate_grid <- c("zero", "binary", "linear", "nonlinear")
ps_grid   <- c("random", "linear", "interaction", "nonlinear")

results_all <- list()
k <- 1

for (ct in cate_grid) {
  for (pt in ps_grid) {
    cat("Running:", ct, pt, "\n")
    out <- run_mc(R = 300, n = 1000, cate_type = ct, ps_type = pt)
    tab <- make_mc_table(out$est, out$true)
    tab$CATE_type <- ct
    tab$PS_type   <- pt
    results_all[[k]] <- tab
    k <- k + 1
  }
}

final_results <- do.call(rbind, results_all)
final_results

############################
## 7. Example
############################

cate_grid <- c("zero", "binary", "linear", "nonlinear")
ps_grid   <- c("random", "linear", "interaction", "nonlinear")

scenario_tables <- list()

for (ct in cate_grid) {
  for (pt in ps_grid) {
    key <- paste0("CATE_", ct, "_PS_", pt)
    scenario_tables[[key]] <- subset(final_results,
                                     CATE_type == ct & PS_type == pt)
  }
}

## Example

scenario_tables[["CATE_zero_PS_random"]]
scenario_tables[["CATE_zero_PS_linear"]]
scenario_tables[["CATE_zero_PS_interaction"]]
scenario_tables[["CATE_zero_PS_nonlinear"]]

scenario_tables[["CATE_binary_PS_random"]]
scenario_tables[["CATE_binary_PS_linear"]]
scenario_tables[["CATE_binary_PS_interaction"]]
scenario_tables[["CATE_binary_PS_nonlinear"]]

scenario_tables[["CATE_linear_PS_random"]]
scenario_tables[["CATE_linear_PS_linear"]]
scenario_tables[["CATE_linear_PS_interaction"]]
scenario_tables[["CATE_linear_PS_nonlinear"]]

scenario_tables[["CATE_nonlinear_PS_random"]]
scenario_tables[["CATE_nonlinear_PS_linear"]]
scenario_tables[["CATE_nonlinear_PS_interaction"]]
scenario_tables[["CATE_nonlinear_PS_nonlinear"]]

