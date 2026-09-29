#############################################################################
# unit_root_analysis_3.R
# Self-contained script to run ADF, FADF, ADF-SB, FADF-SB, and ARNN-ADF
#############################################################################

library(MASS)

# 1. Set Working Directory
setwd("C:/R/arnn")

# 2. Define the core estimation function
run_all_tests <- function(data_matrix, LL_fix = 1, years = 1997:2023) {
  drr <- nrow(data_matrix)
  dcc <- ncol(data_matrix)
  
  # Expanded output matrix to include the 3 ARNN-ADF variants
  mout <- matrix(0, nrow = dcc, ncol = 11)
  colnames(mout) <- c("ADF", "FADF", "kk_fadf", "ADF_BP", "break_year1", 
                      "FADF_BP", "kk_fadf_bp", "break_year2",
                      "ARNN_ADF_1", "ARNN_ADF_2", "ARNN_ADF_3")
  
  for (nc in 1:dcc) {
    cat("Processing Series", nc, "...\n")
    
    y <- data_matrix[, nc]
    n <- length(y)
    dy <- diff(y)
    
    const <- rep(1, n)
    tt <- 1:n
    idx_reg <- (LL_fix + 2):n
    y_dep <- dy[(LL_fix + 1):(n - 1)]
    
    # =====================================================================
    # 1. Standard ADF
    # =====================================================================
    X_adf <- cbind(const[idx_reg], tt[idx_reg], y[(LL_fix + 1):(n - 1)], dy[LL_fix:(n - 2)])
    invx <- solve(t(X_adf) %*% X_adf)
    beta <- invx %*% t(X_adf) %*% y_dep
    e <- y_dep - X_adf %*% beta
    ee1 <- as.numeric(t(e) %*% e)
    rss <- ee1 / (length(y_dep) - ncol(X_adf))
    sd_err <- sqrt(diag(rss * invx))
    adf <- (beta / sd_err)[3, 1]
    
    mout[nc, 1] <- round(adf, 4)
    
    # =====================================================================
    # 2. Fourier ADF
    # =====================================================================
    rss_mins <- numeric(2)
    for (kk in 1:2) {
      sinf <- sin((2 * pi * tt * kk) / n)
      cosf <- cos((2 * pi * tt * kk) / n)
      X_f <- cbind(const[idx_reg], tt[idx_reg], sinf[idx_reg], cosf[idx_reg], 
                   y[(LL_fix + 1):(n - 1)], dy[LL_fix:(n - 2)])
      invx_f <- solve(t(X_f) %*% X_f)
      beta_f <- invx_f %*% t(X_f) %*% y_dep
      e_f <- y_dep - X_f %*% beta_f
      rss_mins[kk] <- as.numeric(t(e_f) %*% e_f)
    }
    
    kk_opt <- if (rss_mins[1] < rss_mins[2]) 1 else 2
    
    sinf <- sin((2 * pi * tt * kk_opt) / n)
    cosf <- cos((2 * pi * tt * kk_opt) / n)
    X_fadf <- cbind(const[idx_reg], tt[idx_reg], sinf[idx_reg], cosf[idx_reg], 
                    y[(LL_fix + 1):(n - 1)], dy[LL_fix:(n - 2)])
    invx_fadf <- solve(t(X_fadf) %*% X_fadf)
    beta_fadf <- invx_fadf %*% t(X_fadf) %*% y_dep
    e_fadf <- y_dep - X_fadf %*% beta_fadf
    ee2 <- as.numeric(t(e_fadf) %*% e_fadf)
    rss_fadf <- ee2 / (length(y_dep) - ncol(X_fadf))
    sd_fadf <- sqrt(diag(rss_fadf * invx_fadf))
    fadf <- beta_fadf[5, 1] / sd_fadf[5]
    
    mout[nc, 2] <- round(fadf, 4)
    mout[nc, 3] <- kk_opt
    
    # =====================================================================
    # 3. ADF with Structural Break
    # =====================================================================
    tadf1_min <- 10
    ddd1_min <- 0
    
    for (ddd in (LL_fix + 2):(drr - 1)) {
      dbp <- rep(0, drr); dbp[ddd] <- 1
      du  <- rep(0, drr); if (ddd < drr) du[(ddd + 1):drr] <- 1
      
      X_sb <- cbind(const[idx_reg], tt[idx_reg], dbp[idx_reg], du[idx_reg], 
                    y[(LL_fix + 1):(n - 1)], dy[LL_fix:(n - 2)])
      
      XtX <- t(X_sb) %*% X_sb
      invx_sb <- tryCatch(solve(XtX), error = function(e) MASS::ginv(XtX))
      
      beta_sb <- invx_sb %*% t(X_sb) %*% y_dep
      e_sb <- y_dep - X_sb %*% beta_sb
      ee3 <- as.numeric(t(e_sb) %*% e_sb)
      rss_sb <- ee3 / (length(y_dep) - ncol(X_sb))
      sd_sb <- sqrt(diag(rss_sb * invx_sb))
      adf_bp <- beta_sb[5, 1] / sd_sb[5]
      
      if (!is.na(adf_bp) && adf_bp < tadf1_min) {
        tadf1_min <- adf_bp
        ddd1_min <- ddd
      }
    }
    
    mout[nc, 4] <- round(tadf1_min, 4)
    mout[nc, 5] <- years[ddd1_min]
    
    # =====================================================================
    # 4. Fourier ADF with Structural Break
    # =====================================================================
    tadf1_min <- 10
    ddd1_min <- 0
    kk_min <- 0
    
    for (kk in 1:2) {
      sinf <- sin((2 * pi * tt * kk) / n)
      cosf <- cos((2 * pi * tt * kk) / n)
      
      for (ddd in (LL_fix + 2):(drr - 1)) {
        dbp <- rep(0, drr); dbp[ddd] <- 1
        du  <- rep(0, drr); if (ddd < drr) du[(ddd + 1):drr] <- 1
        
        X_fsb <- cbind(const[idx_reg], tt[idx_reg], dbp[idx_reg], du[idx_reg], 
                       sinf[idx_reg], cosf[idx_reg], 
                       y[(LL_fix + 1):(n - 1)], dy[LL_fix:(n - 2)])
        
        XtX <- t(X_fsb) %*% X_fsb
        invx_fsb <- tryCatch(solve(XtX), error = function(e) MASS::ginv(XtX))
        
        beta_fsb <- invx_fsb %*% t(X_fsb) %*% y_dep
        e_fsb <- y_dep - X_fsb %*% beta_fsb
        ee4 <- as.numeric(t(e_fsb) %*% e_fsb)
        rss_fsb <- ee4 / (length(y_dep) - ncol(X_fsb))
        sd_fsb <- sqrt(diag(rss_fsb * invx_fsb))
        fadf_bp <- beta_fsb[7, 1] / sd_fsb[7]
        
        if (!is.na(fadf_bp) && fadf_bp < tadf1_min) {
          tadf1_min <- fadf_bp
          ddd1_min <- ddd
          kk_min <- kk
        }
      }
    }
    
    mout[nc, 6] <- round(tadf1_min, 4)
    mout[nc, 7] <- kk_min
    mout[nc, 8] <- years[ddd1_min]
    
    # =====================================================================
    # 5. ARNN-ADF Tests (Variants 1, 2, & 3)
    # =====================================================================
    # Core variables
    y_c  <- y[(LL_fix + 2):n]
    y1_c <- y[(LL_fix + 1):(n - 1)]
    y2_c <- y[LL_fix:(n - 2)]
    
    dyv  <- y[(LL_fix + 1):n] - y[LL_fix:(n - 1)]
    dy1v <- dyv[LL_fix:(n - 2)]
    dep_arnn <- dyv[(LL_fix + 1):(n - 1)]
    
    # Polynomial expansions
    y1y1 <- y1_c * y1_c
    y1y1y1 <- y1y1 * y1_c
    
    # --- Variant 1 ---
    yy <- y_c * y_c
    yy1 <- y_c * y1_c
    yyy <- yy * y_c
    yyy1 <- yy * y1_c
    yy1y1 <- yy1 * y1_c
    
    X_a1 <- cbind(const[idx_reg], tt[idx_reg], y1_c, yy, yy1, y1y1, yyy, yyy1, yy1y1, y1y1y1, dy1v)
    invx_a1 <- tryCatch(solve(t(X_a1) %*% X_a1), error = function(e) MASS::ginv(t(X_a1) %*% X_a1))
    beta_a1 <- invx_a1 %*% t(X_a1) %*% dep_arnn
    e_a1 <- dep_arnn - X_a1 %*% beta_a1
    rss_a1 <- as.numeric(t(e_a1) %*% e_a1) / (length(dep_arnn) - ncol(X_a1))
    sd_a1 <- sqrt(diag(rss_a1 * invx_a1))
    mout[nc, 9] <- round((beta_a1 / sd_a1)[3, 1], 4)
    
    # --- Variant 2 ---
    y1y2 <- y1_c * y2_c
    y2y2 <- y2_c * y2_c
    y1y1y2 <- y1y1 * y2_c
    y1y2y2 <- y1y2 * y2_c
    y2y2y2 <- y2y2 * y2_c
    
    X_a2 <- cbind(const[idx_reg], tt[idx_reg], y1_c, y1y1, y1y2, y2y2, y1y1y1, y1y1y2, y1y2y2, y2y2y2, dy1v)
    invx_a2 <- tryCatch(solve(t(X_a2) %*% X_a2), error = function(e) MASS::ginv(t(X_a2) %*% X_a2))
    beta_a2 <- invx_a2 %*% t(X_a2) %*% dep_arnn
    e_a2 <- dep_arnn - X_a2 %*% beta_a2
    rss_a2 <- as.numeric(t(e_a2) %*% e_a2) / (length(dep_arnn) - ncol(X_a2))
    sd_a2 <- sqrt(diag(rss_a2 * invx_a2))
    mout[nc, 10] <- round((beta_a2 / sd_a2)[3, 1], 4)
    
    # --- Variant 3 ---
    dyv3 <- y[(LL_fix + 2):n] - y[(LL_fix + 1):(n - 1)]
    dy1v3 <- y[(LL_fix + 1):(n - 1)] - y[LL_fix:(n - 2)]
    
    X_a3 <- cbind(const[idx_reg], tt[idx_reg], y1_c, y1y1, y1y1y1, dy1v3)
    invx_a3 <- tryCatch(solve(t(X_a3) %*% X_a3), error = function(e) MASS::ginv(t(X_a3) %*% X_a3))
    beta_a3 <- invx_a3 %*% t(X_a3) %*% dyv3
    e_a3 <- dyv3 - X_a3 %*% beta_a3
    rss_a3 <- as.numeric(t(e_a3) %*% e_a3) / (length(dyv3) - ncol(X_a3))
    sd_a3 <- sqrt(diag(rss_a3 * invx_a3))
    mout[nc, 11] <- round((beta_a3 / sd_a3)[3, 1], 4)
  }
  
  return(as.data.frame(mout))
}

# 3. Read dataset and execute analysis
data_path <- if (file.exists("data1_2.txt")) "data1_2.txt" else "data1.txt"
data_mat  <- as.matrix(read.table(data_path, header = FALSE))

results <- run_all_tests(data_mat, LL_fix = 1, years = 1997:2023)

# 4. Format, Display, and Save Results
series_labels <- paste0("Series_", 1:ncol(data_mat))
final_output <- cbind(Series = series_labels, results)

print(final_output)
write.csv(final_output, "all_unit_root_results.csv", row.names = FALSE)