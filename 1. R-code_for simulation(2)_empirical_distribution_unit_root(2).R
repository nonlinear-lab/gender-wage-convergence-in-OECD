# =====================================================================
# MONTE CARLO SIMULATION: EMPIRICAL DISTRIBUTIONS FOR ADF-TYPE TESTS
# Includes: ADF, FADF, ADF-SB, FADF-SB, ARNN-ADF (Variants 1, 2, 3)
# =====================================================================

library(MASS)

# Optional: Set your Working Directory (WD) where the CSV will be saved
setwd("C:/R/arnn") 

# ---------------------------------------------------------------------
# PART 1: Data Generating Process (DGP) for Unit Root Null Hypothesis
# ---------------------------------------------------------------------
# Generates a pure random walk: y_t = y_{t-1} + u_t
simulate_unit_root_dgp <- function(n, rho = 1) {
  u <- rnorm(n, mean = 0, sd = 1)
  y <- numeric(n)
  y[1] <- u[1]
  for (t in 2:n) {
    y[t] <- rho * y[t-1] + u[t]
  }
  return(y)
}

# ---------------------------------------------------------------------
# PART 2: Unit Root Test Statistic Calculators (Fixed k and omega)
# ---------------------------------------------------------------------

get_adf_stats <- function(y, k, omega, LL_fix = 1) {
  n <- length(y)
  dy <- diff(y)
  
  const <- rep(1, n)
  tt <- 1:n
  idx_reg <- (LL_fix + 2):n
  y_dep <- dy[(LL_fix + 1):(n - 1)]
  
  # Calculate exact break index based on omega
  ddd <- floor(omega * n)
  
  # --- 1. Standard ADF ---
  X_adf <- cbind(const[idx_reg], tt[idx_reg], y[(LL_fix + 1):(n - 1)], dy[LL_fix:(n - 2)])
  invx <- tryCatch(solve(t(X_adf) %*% X_adf), error = function(e) MASS::ginv(t(X_adf) %*% X_adf))
  beta <- invx %*% t(X_adf) %*% y_dep
  e <- y_dep - X_adf %*% beta
  ee1 <- as.numeric(t(e) %*% e)
  rss <- ee1 / (length(y_dep) - ncol(X_adf))
  sd_err <- sqrt(diag(rss * invx))
  t_adf <- (beta / sd_err)[3, 1]
  
  # --- 2. Fourier ADF (FADF) for fixed k ---
  sinf <- sin((2 * pi * tt * k) / n)
  cosf <- cos((2 * pi * tt * k) / n)
  X_fadf <- cbind(const[idx_reg], tt[idx_reg], sinf[idx_reg], cosf[idx_reg], 
                  y[(LL_fix + 1):(n - 1)], dy[LL_fix:(n - 2)])
  invx_fadf <- tryCatch(solve(t(X_fadf) %*% X_fadf), error = function(e) MASS::ginv(t(X_fadf) %*% X_fadf))
  beta_fadf <- invx_fadf %*% t(X_fadf) %*% y_dep
  e_fadf <- y_dep - X_fadf %*% beta_fadf
  ee2 <- as.numeric(t(e_fadf) %*% e_fadf)
  rss_fadf <- ee2 / (length(y_dep) - ncol(X_fadf))
  sd_fadf <- sqrt(diag(rss_fadf * invx_fadf))
  t_fadf <- (beta_fadf / sd_fadf)[5, 1]
  
  # --- 3. ADF with Structural Break (ADF-SB) for fixed omega ---
  dbp <- rep(0, n)
  if (ddd >= 1 && ddd <= n) dbp[ddd] <- 1
  
  du <- rep(0, n)
  if (ddd >= 1 && ddd < n) du[(ddd + 1):n] <- 1
  
  X_sb <- cbind(const[idx_reg], tt[idx_reg], dbp[idx_reg], du[idx_reg], 
                y[(LL_fix + 1):(n - 1)], dy[LL_fix:(n - 2)])
  invx_sb <- tryCatch(solve(t(X_sb) %*% X_sb), error = function(e) MASS::ginv(t(X_sb) %*% X_sb))
  beta_sb <- invx_sb %*% t(X_sb) %*% y_dep
  e_sb <- y_dep - X_sb %*% beta_sb
  rss_sb <- as.numeric(t(e_sb) %*% e_sb) / (length(y_dep) - ncol(X_sb))
  sd_sb <- sqrt(diag(rss_sb * invx_sb))
  t_adf_sb <- (beta_sb / sd_sb)[5, 1]
  
  # --- 4. Fourier ADF with Structural Break (FADF-SB) for fixed k and omega ---
  X_fsb <- cbind(const[idx_reg], tt[idx_reg], dbp[idx_reg], du[idx_reg], 
                 sinf[idx_reg], cosf[idx_reg], 
                 y[(LL_fix + 1):(n - 1)], dy[LL_fix:(n - 2)])
  invx_fsb <- tryCatch(solve(t(X_fsb) %*% X_fsb), error = function(e) MASS::ginv(t(X_fsb) %*% X_fsb))
  beta_fsb <- invx_fsb %*% t(X_fsb) %*% y_dep
  e_fsb <- y_dep - X_fsb %*% beta_fsb
  rss_fsb <- as.numeric(t(e_fsb) %*% e_fsb) / (length(y_dep) - ncol(X_fsb))
  sd_fsb <- sqrt(diag(rss_fsb * invx_fsb))
  t_fadf_sb <- (beta_fsb / sd_fsb)[7, 1]
  
  # --- 5. ARNN-ADF Variants (Independent of k and omega) ---
  y_c  <- y[(LL_fix + 2):n]
  y1_c <- y[(LL_fix + 1):(n - 1)]
  y2_c <- y[LL_fix:(n - 2)]
  dyv  <- y[(LL_fix + 1):n] - y[LL_fix:(n - 1)]
  dy1v <- dyv[LL_fix:(n - 2)]
  dep_arnn <- dyv[(LL_fix + 1):(n - 1)]
  
  y1y1 <- y1_c * y1_c
  y1y1y1 <- y1y1 * y1_c
  
  # Variant 1
  yy <- y_c * y_c; yy1 <- y_c * y1_c; yyy <- yy * y_c; yyy1 <- yy * y1_c; yy1y1 <- yy1 * y1_c
  X_a1 <- cbind(const[idx_reg], tt[idx_reg], y1_c, yy, yy1, y1y1, yyy, yyy1, yy1y1, y1y1y1, dy1v)
  invx_a1 <- tryCatch(solve(t(X_a1) %*% X_a1), error = function(e) MASS::ginv(t(X_a1) %*% X_a1))
  beta_a1 <- invx_a1 %*% t(X_a1) %*% dep_arnn
  e_a1 <- dep_arnn - X_a1 %*% beta_a1
  rss_a1 <- as.numeric(t(e_a1) %*% e_a1) / (length(dep_arnn) - ncol(X_a1))
  sd_a1 <- sqrt(diag(rss_a1 * invx_a1))
  t_arnn1 <- (beta_a1 / sd_a1)[3, 1]
  
  # Variant 2
  y1y2 <- y1_c * y2_c; y2y2 <- y2_c * y2_c; y1y1y2 <- y1y1 * y2_c; y1y2y2 <- y1y2 * y2_c; y2y2y2 <- y2y2 * y2_c
  X_a2 <- cbind(const[idx_reg], tt[idx_reg], y1_c, y1y1, y1y2, y2y2, y1y1y1, y1y1y2, y1y2y2, y2y2y2, dy1v)
  invx_a2 <- tryCatch(solve(t(X_a2) %*% X_a2), error = function(e) MASS::ginv(t(X_a2) %*% X_a2))
  beta_a2 <- invx_a2 %*% t(X_a2) %*% dep_arnn
  e_a2 <- dep_arnn - X_a2 %*% beta_a2
  rss_a2 <- as.numeric(t(e_a2) %*% e_a2) / (length(dep_arnn) - ncol(X_a2))
  sd_a2 <- sqrt(diag(rss_a2 * invx_a2))
  t_arnn2 <- (beta_a2 / sd_a2)[3, 1]
  
  # Variant 3
  dyv3 <- y[(LL_fix + 2):n] - y[(LL_fix + 1):(n - 1)]
  dy1v3 <- y[(LL_fix + 1):(n - 1)] - y[LL_fix:(n - 2)]
  X_a3 <- cbind(const[idx_reg], tt[idx_reg], y1_c, y1y1, y1y1y1, dy1v3)
  invx_a3 <- tryCatch(solve(t(X_a3) %*% X_a3), error = function(e) MASS::ginv(t(X_a3) %*% X_a3))
  beta_a3 <- invx_a3 %*% t(X_a3) %*% dyv3
  e_a3 <- dyv3 - X_a3 %*% beta_a3
  rss_a3 <- as.numeric(t(e_a3) %*% e_a3) / (length(dyv3) - ncol(X_a3))
  sd_a3 <- sqrt(diag(rss_a3 * invx_a3))
  t_arnn3 <- (beta_a3 / sd_a3)[3, 1]
  
  return(c(t_adf, t_fadf, t_adf_sb, t_fadf_sb, t_arnn1, t_arnn2, t_arnn3))
}

# ---------------------------------------------------------------------
# PART 3: General Empirical Distribution Estimator
# ---------------------------------------------------------------------
estimate_empirical_cvs <- function(n_sim = 5000, n_obs = 100, k = 1, omega = 0.125) {
  
  # Matrix to store results for 7 tests
  sim_results <- matrix(NA, nrow = n_sim, ncol = 7)
  
  for (i in 1:n_sim) {
    # Generate unit root series (H0)
    y <- simulate_unit_root_dgp(n_obs, rho = 1)
    
    # Calculate all ADF statistics using the specific k and omega
    sim_results[i, ] <- get_adf_stats(y, k = k, omega = omega, LL_fix = 1)
  }
  
  # For ADF tests, we look at the LEFT tail (1%, 5%, 10%)
  percentiles <- c(0.01, 0.05, 0.10)
  
  cv_matrix <- apply(sim_results, 2, quantile, probs = percentiles, na.rm = TRUE)
  rownames(cv_matrix) <- c("1%", "5%", "10%")
  colnames(cv_matrix) <- c("ADF", "FADF", "ADF_SB", "FADF_SB", "ARNN_ADF_1", "ARNN_ADF_2", "ARNN_ADF_3")
  
  return(t(cv_matrix)) # Transpose for nicer formatting
}

# ---------------------------------------------------------------------
# PART 4: Execution Setup
# ---------------------------------------------------------------------
set.seed(123)

# Define grid of parameters
n_values      <- c(50, 100, 500, 1000)
k_values      <- c(1, 2, 3)
omega_values  <- c(0.125, 0.375, 0.625, 0.875)
n_simulations <- 5000

cat("========================================================\n")
cat(" ESTIMATING EMPIRICAL CRITICAL VALUES FOR ADF TESTS \n")
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
      
      # Format to data frame
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
final_results <- final_results[, c("n", "k", "omega", "test", "1%", "5%", "10%")]

# Save the full grid to CSV
write.csv(final_results, "C:/R/arnn/empirical_distribution_ADF.csv", row.names = FALSE)

cat("\nDone. Full results (n x k x omega x test) saved to empirical_distribution_ADF.csv\n")