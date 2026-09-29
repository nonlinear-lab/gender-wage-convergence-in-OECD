#############################################################################
# unit_root_analysis4.R
# Self-contained script to run KPSS, FKPSS, KPSS-SB, FKPSS-SB, and ARNN-KPSS
#############################################################################

library(MASS)

# 1. Set Working Directory
setwd("C:/R/arnn")

# --- Helper Function: Frequency Zero Spectrum (FZS) ---
# Replicates the manual summation lag loops in the GAUSS code for LL <= 6
calc_fzs <- function(e, LL) {
  n <- length(e)
  e3 <- sum(e^2) / n
  
  if (LL == 0) return(e3)
  
  gam <- numeric(LL)
  for (j in 1:LL) {
    sum_lag <- 0
    for (i in (j + 1):n) {
      sum_lag <- sum_lag + e[i] * e[i - j]
    }
    gam[j] <- sum_lag / n
  }
  
  w <- (LL + 1 - 1:LL) / (LL + 1)
  fzs <- e3 + 2 * sum(w * gam)
  
  return(fzs)
}

# 2. Define the core estimation function
run_kpss_family <- function(data_matrix, LL_fix = 1, years = 1997:2023) {
  drr <- nrow(data_matrix)
  dcc <- ncol(data_matrix)
  
  # Output matrix for the 5 tests + optimal breaks/frequencies
  mout <- matrix(0, nrow = dcc, ncol = 10)
  colnames(mout) <- c("KPSS", "FKPSS", "kk_fkpss", "KPSS_SB", "break_year1", 
                      "FKPSS_SB", "kk_fkpss_sb", "break_year2", "ARNN_KPSS", "break_pct")
  
  for (nc in 1:dcc) {
    cat("Processing Series", nc, "...\n")
    
    y <- data_matrix[, nc]
    n0 <- length(y)
    n <- n0
    
    const <- rep(1, n)
    tt <- 1:n
    trend <- tt
    
    # =====================================================================
    # 1. Standard KPSS
    # =====================================================================
    X_kpss <- cbind(const, trend)
    invx <- solve(t(X_kpss) %*% X_kpss)
    beta <- invx %*% t(X_kpss) %*% y
    e0 <- y - X_kpss %*% beta
    
    # Calculate LM statistic
    S <- cumsum(e0)
    e2 <- sum(S^2)
    fzs <- calc_fzs(e0, LL_fix)
    lm_kpss <- e2 / ((n0^2) * fzs)
    
    mout[nc, 1] <- round(lm_kpss, 4)
    
    # =====================================================================
    # 2. Fourier KPSS (FKPSS)
    # =====================================================================
    ssr_f_mins <- numeric(2)
    e_f_list <- list()
    
    for (kk in 1:2) {
      sinf <- sin((2 * pi * tt * kk) / n)
      cosf <- cos((2 * pi * tt * kk) / n)
      X_f <- cbind(const, trend, sinf, cosf)
      invx_f <- solve(t(X_f) %*% X_f)
      beta_f <- invx_f %*% t(X_f) %*% y
      ef <- y - X_f %*% beta_f
      
      ssr_f_mins[kk] <- as.numeric(t(ef) %*% ef)
      e_f_list[[kk]] <- ef
    }
    
    kk_opt <- if (ssr_f_mins[1] < ssr_f_mins[2]) 1 else 2
    ef0 <- e_f_list[[kk_opt]]
    
    S_f <- cumsum(ef0)
    e2_f <- sum(S_f^2)
    fzs_f <- calc_fzs(ef0, LL_fix)
    lm_fkpss <- e2_f / ((n0^2) * fzs_f)
    
    mout[nc, 2] <- round(lm_fkpss, 4)
    mout[nc, 3] <- kk_opt
    
    # =====================================================================
    # 3. KPSS with Structural Break (KPSS-SB)
    # =====================================================================
    ssr1_min <- 1000
    ddd1_min <- 0
    esb0 <- numeric(n)
    
    for (ddd in (LL_fix + 2):(n - 1)) {
      dbp <- rep(0, n); dbp[ddd] <- 1
      du  <- rep(0, n); if (ddd < n) du[(ddd + 1):n] <- 1
      
      X_sb <- cbind(const, trend, dbp, du)
      XtX <- t(X_sb) %*% X_sb
      invx_sb <- tryCatch(solve(XtX), error = function(e) MASS::ginv(XtX))
      
      beta_sb <- invx_sb %*% t(X_sb) %*% y
      e_sb <- y - X_sb %*% beta_sb
      ssr <- as.numeric(t(e_sb) %*% e_sb) / (n - ncol(X_sb))
      
      if (ssr < ssr1_min) {
        ssr1_min <- ssr
        ddd1_min <- ddd
        esb0 <- e_sb
      }
    }
    
    S_sb <- cumsum(esb0)
    e2_sb <- sum(S_sb^2)
    fzs_sb <- calc_fzs(esb0, LL_fix)
    lm_kpss_sb <- e2_sb / ((n0^2) * fzs_sb)
    
    mout[nc, 4] <- round(lm_kpss_sb, 4)
    mout[nc, 5] <- years[ddd1_min]
    
    # =====================================================================
    # 4. Fourier KPSS with Structural Break (FKPSS-SB)
    # =====================================================================
    ssr_bp1_min <- 1000
    ddd2_min <- 0
    fkk_bp <- 0
    ef_bp0 <- numeric(n)
    
    for (kk in 1:2) {
      sinf <- sin((2 * pi * tt * kk) / n)
      cosf <- cos((2 * pi * tt * kk) / n)
      
      for (ddd in (LL_fix + 2):(n - 1)) {
        dbp <- rep(0, n); dbp[ddd] <- 1
        du  <- rep(0, n); if (ddd < n) du[(ddd + 1):n] <- 1
        
        X_fsb <- cbind(const, trend, dbp, du, sinf, cosf)
        XtX <- t(X_fsb) %*% X_fsb
        invx_fsb <- tryCatch(solve(XtX), error = function(e) MASS::ginv(XtX))
        
        beta_fsb <- invx_fsb %*% t(X_fsb) %*% y
        ef_bp <- y - X_fsb %*% beta_fsb
        ssr_bp <- as.numeric(t(ef_bp) %*% ef_bp) / (n - ncol(X_fsb))
        
        if (ssr_bp < ssr_bp1_min) {
          ssr_bp1_min <- ssr_bp
          ddd2_min <- ddd
          fkk_bp <- kk
          ef_bp0 <- ef_bp
        }
      }
    }
    
    S_fsb <- cumsum(ef_bp0)
    e2_fsb <- sum(S_fsb^2)
    fzs_fsb <- calc_fzs(ef_bp0, LL_fix)
    lm_fkpss_sb <- e2_fsb / ((n0^2) * fzs_fsb)
    
    mout[nc, 6] <- round(lm_fkpss_sb, 4)
    mout[nc, 7] <- fkk_bp
    mout[nc, 8] <- years[ddd2_min]
    mout[nc, 10] <- round((ddd2_min / n0) * 100, 2)
    
    # =====================================================================
    # 5. ARNN-KPSS (Non-Linear Autoregressive Neural Network Augmented)
    # =====================================================================
    yv <- y[2:n]
    y1 <- y[1:(n - 1)]
    y1y1 <- y1 * y1
    y1y1y1 <- y1y1 * y1
    
    const1 <- const[1:(n - 1)]
    trend1 <- trend[1:(n - 1)]
    
    X_nn <- cbind(const1, trend1, y1, y1y1, y1y1y1)
    XtX_nn <- t(X_nn) %*% X_nn
    invxx <- tryCatch(solve(XtX_nn), error = function(e) MASS::ginv(XtX_nn))
    
    beta_nn <- invxx %*% t(X_nn) %*% yv
    enn0 <- yv - X_nn %*% beta_nn
    
    S_nn <- cumsum(enn0)
    e2_nn <- sum(S_nn^2)
    fzs_nn <- calc_fzs(enn0, LL_fix)
    # Note: GAUSS ARNN-KPSS divides by (n0^2), the original sample length.
    lm_arnn_kpss <- e2_nn / ((n0^2) * fzs_nn) 
    
    mout[nc, 9] <- round(lm_arnn_kpss, 4)
  }
  
  return(as.data.frame(mout))
}

# 3. Read dataset and execute analysis
data_path <- if (file.exists("data1_2.txt")) "data1_2.txt" else "data1.txt"
data_mat  <- as.matrix(read.table(data_path, header = FALSE))

results <- run_kpss_family(data_mat, LL_fix = 1, years = 1997:2023)

# 4. Format, Display, and Save Results
series_labels <- paste0("Series_", 1:ncol(data_mat))
final_output <- cbind(Series = series_labels, results)

print(final_output)
write.csv(final_output, "kpss_family_results.csv", row.names = FALSE)