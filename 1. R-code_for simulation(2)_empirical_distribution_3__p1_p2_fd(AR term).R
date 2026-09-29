
# =====================================================================
# MONTE CARLO SIMULATION: EMPIRICAL DISTRIBUTIONS FOR KPSS-TYPE TESTS
# Includes: KPSS, FKPSS, KPSS-SB, FKPSS-SB, KPSS-YOFA (fd)
# =====================================================================

# Optional: Set your Working Directory (WD) where the graphs will be saved
setwd("C:/R/arnn") 

# ---------------------------------------------------------------------
# PART 1: Data Generating Process (DGP) & Long-Run Variance
# ---------------------------------------------------------------------
# Generates a trend-stationary or difference-stationary process with AR(1) errors
simulate_kpss_dgp <- function(n, lambda, rho = 0, mu = 0, beta = 0) {
  v <- rnorm(n, mean = 0, sd = 1)
  eps <- numeric(n)
  eps[1] <- v[1]
  if (n > 1) {
    for (t in 2:n) {
      eps[t] <- rho * eps[t-1] + v[t]
    }
  }
  
  var_eps <- 1 / (1 - rho^2) 
  var_u <- lambda * var_eps
  sd_u <- sqrt(var_u)
  
  u <- rnorm(n, mean = 0, sd = sd_u)
  r <- cumsum(u)
  
  t_seq <- 1:n
  y <- mu + beta * t_seq + r + eps
  return(y)
}

# Calculates Long-Run Variance using Bartlett Kernel
calc_lrv <- function(e, l) {
  n <- length(e)
  var_e <- sum(e^2) / n
  if (l > 0) {
    cov_e <- 0
    for (s in 1:l) {
      weight <- 1 - (s / (l + 1))
      cov_s <- sum(e[(s + 1):n] * e[1:(n - s)]) / n
      cov_e <- cov_e + weight * cov_s
    }
    lrv <- var_e + 2 * cov_e
  } else {
    lrv <- var_e
  }
  return(lrv)
}

# ---------------------------------------------------------------------
# PART 2: Test Statistics
# ---------------------------------------------------------------------
# 1. KPSS Test
kpss_test <- function(y, l) {
  t_seq <- 1:length(y)
  mod <- lm(y ~ t_seq)
  e <- residuals(mod)
  eta <- sum(cumsum(e)^2) / (length(e)^2 * calc_lrv(e, l))
  return(eta)
}

# 2. Fourier KPSS (FKPSS) Test
fkpss_test <- function(y, k, l) {
  n <- length(y)
  t_seq <- 1:n
  sin_t <- sin(2 * pi * k * t_seq / n)
  cos_t <- cos(2 * pi * k * t_seq / n)
  
  mod <- lm(y ~ t_seq + sin_t + cos_t)
  e <- residuals(mod)
  eta <- sum(cumsum(e)^2) / (n^2 * calc_lrv(e, l))
  return(eta)
}

# 3. KPSS with Structural Break (KPSS-SB) Test
kpss_sb_test <- function(y, omega, l) {
  n <- length(y)
  t_seq <- 1:n
  n_B <- floor(omega * n)
  
  DU <- ifelse(t_seq > n_B, 1, 0)
  DT_B <- ifelse(t_seq == n_B, 1, 0)
  
  mod <- lm(y ~ t_seq + DU + DT_B)
  e <- residuals(mod)
  eta <- sum(cumsum(e)^2) / (n^2 * calc_lrv(e, l))
  return(eta)
}

# 4. Fourier KPSS with Structural Break (FKPSS-SB) Test
fkpss_sb_test <- function(y, k, omega, l) {
  n <- length(y)
  t_seq <- 1:n
  n_B <- floor(omega * n)
  
  sin_t <- sin(2 * pi * k * t_seq / n)
  cos_t <- cos(2 * pi * k * t_seq / n)
  DU <- ifelse(t_seq > n_B, 1, 0)
  DT_B <- ifelse(t_seq == n_B, 1, 0)
  
  mod <- lm(y ~ t_seq + sin_t + cos_t + DU + DT_B)
  e <- residuals(mod)
  eta <- sum(cumsum(e)^2) / (n^2 * calc_lrv(e, l))
  return(eta)
}

# 5. Original ARNN-KPSS Test (Autoregressive Level Input: p=1 and p=2)
arnn_kpss_test <- function(y, l, p = 1) {
  n <- length(y)
  
  if (p == 1) {
    y_curr <- y[2:n]
    t_trend <- 2:n
    y_lag1 <- y[1:(n - 1)]
    mod <- lm(y_curr ~ t_trend + y_lag1 + I(y_lag1^2) + I(y_lag1^3))
    
  } else if (p == 2) {
    y_curr <- y[3:n]
    t_trend <- 3:n
    y_lag1 <- y[2:(n - 1)]
    y_lag2 <- y[1:(n - 2)]
    mod <- lm(y_curr ~ t_trend + y_lag1 + y_lag2 + 
                I(y_lag1^2) + I(y_lag2^2) + I(y_lag1 * y_lag2) + 
                I(y_lag1^3) + I(y_lag2^3) + I((y_lag1^2) * y_lag2) + I(y_lag1 * (y_lag2^2)))
  }
  
  e <- residuals(mod)
  n_adj <- length(e)
  eta <- sum(cumsum(e)^2) / (n_adj^2 * calc_lrv(e, l))
  return(eta)
}

# 6. ARNN-KPSS-t Test (Trend-Only Input)
arnn_kpss_t_test <- function(y, l) {
  n <- length(y)
  tau <- (1:n) / n # Normalized time trend
  
  mod <- lm(y ~ tau + I(tau^2) + I(tau^3))
  
  e <- residuals(mod)
  eta <- sum(cumsum(e)^2) / (n^2 * calc_lrv(e, l))
  return(eta)
}

# 7. NEW: ARNN-KPSS-fd Test (First-Differenced Autoregressive Input)
arnn_kpss_fd_test <- function(y, l) {
  n <- length(y)
  
  # Need to start at t=3 because we need a lag (t-1) AND a difference (t-1 minus t-2)
  y_curr <- y[3:n]
  t_trend <- 3:n
  
  # dy_{t-1} = y_{t-1} - y_{t-2}
  dy_lag1 <- y[2:(n - 1)] - y[1:(n - 2)]
  
  # Regression on first-differenced lags up to third-order Taylor expansion
  mod <- lm(y_curr ~ t_trend + dy_lag1 + I(dy_lag1^2) + I(dy_lag1^3))
  
  e <- residuals(mod)
  n_adj <- length(e)
  eta <- sum(cumsum(e)^2) / (n_adj^2 * calc_lrv(e, l))
  return(eta)
}

# ---------------------------------------------------------------------
# PART 3: General Empirical Distribution Estimator (All 8 Models)
# ---------------------------------------------------------------------
estimate_empirical_cvs <- function(n_sim = 5000, n_obs = 100, k = 1, omega = 0.125) {
  
  stats <- list(
    KPSS     = numeric(n_sim),
    FKPSS    = numeric(n_sim),
    KPSS_SB  = numeric(n_sim),
    FKPSS_SB = numeric(n_sim),
    ARNN_p1  = numeric(n_sim),
    ARNN_p2  = numeric(n_sim),
    ARNN_t   = numeric(n_sim),
    ARNN_fd  = numeric(n_sim)
  )
  
  l0 <- 0 # L0 bandwidth for pure empirical critical values 
  
  for (i in 1:n_sim) {
    # Generate pure stationary series (H0)
    y <- simulate_kpss_dgp(n_obs, lambda = 0, rho = 0)
    
    # Store calculated statistic for all variants
    stats$KPSS[i]     <- kpss_test(y, l0)
    stats$FKPSS[i]    <- fkpss_test(y, k, l0)
    stats$KPSS_SB[i]  <- kpss_sb_test(y, omega, l0)
    stats$FKPSS_SB[i] <- fkpss_sb_test(y, k, omega, l0)
    stats$ARNN_p1[i]  <- arnn_kpss_test(y, l0, p = 1)
    stats$ARNN_p2[i]  <- arnn_kpss_test(y, l0, p = 2)
    stats$ARNN_t[i]   <- arnn_kpss_t_test(y, l0)
    stats$ARNN_fd[i]  <- arnn_kpss_fd_test(y, l0)
  }
  
  percentiles <- c(0.90, 0.95, 0.99)
  
  cv_matrix <- rbind(
    "KPSS"             = quantile(stats$KPSS, percentiles),
    "FKPSS"            = quantile(stats$FKPSS, percentiles),
    "KPSS_SB"          = quantile(stats$KPSS_SB, percentiles),
    "FKPSS_SB"         = quantile(stats$FKPSS_SB, percentiles),
    "ARNN_KPSS (p=1)"  = quantile(stats$ARNN_p1, percentiles),
    "ARNN_KPSS (p=2)"  = quantile(stats$ARNN_p2, percentiles),
    "ARNN_KPSS-t"      = quantile(stats$ARNN_t, percentiles),
    "ARNN_KPSS-fd"     = quantile(stats$ARNN_fd, percentiles)
  )
  colnames(cv_matrix) <- c("10%", "5%", "1%")
  
  return(cv_matrix)
}

# ---------------------------------------------------------------------
# PART 4: Execution Setup
# ---------------------------------------------------------------------
set.seed(123)

# Define grid of parameters
n_values     <- c(30, 60, 90, 300, 600, 900)
k_values     <- c(1, 2, 3)
omega_values <- c(0.125, 0.375, 0.625, 0.875)
n_simulations <- 5000

cat("========================================================\n")
cat(" ESTIMATING EMPIRICAL CRITICAL VALUES \n")
cat("========================================================\n")

# Collect every (n, k, omega) combination's results here
all_results <- list()
counter <- 1

for (n in n_values) {
  for (k in k_values) {
    for (omega in omega_values) {

      cat(sprintf("\n--- Critical Values for n = %d, k = %d, omega = %.3f ---\n",
                  n, k, omega))

      cvs <- estimate_empirical_cvs(n_sim = n_simulations,
                                     n_obs = n,
                                     k = k,
                                     omega = omega)
      print(round(cvs, 3))

      # Turn the matrix into a tidy data frame and tag it with n, k, omega
      df <- as.data.frame(cvs)
      df$test  <- rownames(cvs)
      df$n     <- n
      df$k     <- k
      df$omega <- omega
      rownames(df) <- NULL

      all_results[[counter]] <- df
      counter <- counter + 1
    }
  }
}

# Combine all combinations into one long-format table
final_results <- do.call(rbind, all_results)

# Reorder columns nicely: n, k, omega, test, then the percentiles
final_results <- final_results[, c("n", "k", "omega", "test", "10%", "5%", "1%")]

# Save the full grid to CSV
write.csv(final_results, "C:/R/arnn/empirical_distribution2.csv", row.names = FALSE)

cat("\Done. Full results (n x k x omega x test) saved to Empirical_distribution2.csv\n")write.csv(final_results, "C:/R/arnn/empirical_distribution2.csv", row.names = FALSE)

cat("\nDone. Full results (n x k x omega x test) saved to Empirical_distribution2.csv\n")