# ============================================================
# Script: SNP-based benchmark models via 5% random subsampling
# Models: LASSO_SNP, GBLUP, GBLUP_AD
# Validation: 10-fold CV, LOMO-CV, and independent test set
# Iterations: 50 random seeds (101-150)
# ============================================================

library(vroom)
library(data.table)
library(glmnet)
library(rrBLUP)
library(caret)
library(magrittr)

run_date <- Sys.Date()
output_dir <- paste0("SNP_Benchmark_Sensitivity_", run_date)
if (!dir.exists(output_dir)) dir.create(output_dir)
cat("Output directory:", output_dir, "\n")

# -- Parameter settings --------------------------------------
SNP_FILE       <- "hybrid_SNP_500_raw.csv"
PHENO_TRAIN    <- "phenotype_training_336.csv"
PHENO_TEST     <- "phenotype_test_92.csv"
SEEDS          <- 101:150
SNP_RATIO      <- 0.05
MIN_REGRESSOR  <- 1
N_FOLDS        <- 10
TRAITS         <- c("FLW_HPH", "FSTW_HPH", "FSHW_HPH")
meta_cols      <- c("hybrid_id", "cms", "restorer")

# -- Read phenotype data -------------------------------------
cat("== Reading phenotype data ==\n")
pheno_train <- fread(PHENO_TRAIN, encoding = "UTF-8")
pheno_test  <- fread(PHENO_TEST, encoding = "UTF-8")

clean_id <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("[^A-Za-z0-9_]", "", x)
  return(x)
}
pheno_train[, hybrid_id := clean_id(hybrid_id)]
pheno_test[, hybrid_id := clean_id(hybrid_id)]
cat("Training samples:", nrow(pheno_train), "\n")
cat("Testing samples:", nrow(pheno_test), "\n\n")

# -- Read CSV column names -----------------------------------
cat("== Reading CSV column names ==\n")
col_names <- names(vroom(SNP_FILE, n_max = 0, show_col_types = FALSE))
snp_cols_all <- setdiff(col_names, meta_cols)
total_snp <- length(snp_cols_all)
cat("Total SNPs:", total_snp, "\n")
cat("Approximate SNPs after 5% sampling:", round(total_snp * SNP_RATIO), "\n\n")

# -- Helper function: Compute dominance matrix ---------------
compute_dominance_matrix <- function(X) {
  D <- 1 - abs(X - 1)
  D[D < 0] <- 0
  return(D)
}

# -- LASSO (Two-step) model function -------------------------
run_lasso_ols <- function(X_train, y_train, X_test, y_test, n_folds = 10, min_regressor = 1, seed = 42) {
  valid_train <- !is.na(y_train)
  X_tr <- X_train[valid_train, , drop = FALSE]
  y_tr <- y_train[valid_train]
  n_tr <- length(y_tr)
  
  valid_test <- !is.na(y_test)
  X_te <- X_test[valid_test, , drop = FALSE]
  y_te <- y_test[valid_test]
  
  set.seed(seed)
  folds <- createFolds(y_tr, k = n_folds, returnTrain = FALSE)
  
  all_pred_cv <- numeric(n_tr)
  all_true_cv <- numeric(n_tr)
  n_features_vec <- numeric(n_folds)
  
  for (i in 1:n_folds) {
    test_idx  <- folds[[i]]
    train_idx <- setdiff(1:n_tr, test_idx)
    
    X_fold_train <- X_tr[train_idx, , drop = FALSE]
    y_fold_train <- y_tr[train_idx]
    X_fold_test  <- X_tr[test_idx, , drop = FALSE]
    y_fold_test  <- y_tr[test_idx]
    
    m <- tryCatch(
      cv.glmnet(x = X_fold_train, y = y_fold_train, nfolds = 10),
      error = function(e) NULL
    )
    if (is.null(m)) {
      all_pred_cv[test_idx] <- NA
      all_true_cv[test_idx] <- y_fold_test
      n_features_vec[i] <- 0
      next
    }
    
    lambda <- m$lambda.min
    j <- sum(coef(m, s = "lambda.min") != 0) - 1
    if (j < min_regressor) {
      idx <- which(m$nzero >= min_regressor)
      if (length(idx) == 0) {
        all_pred_cv[test_idx] <- NA
        all_true_cv[test_idx] <- y_fold_test
        n_features_vec[i] <- 0
        next
      }
      lambda <- m$lambda[max(idx)]
    }
    
    tmp <- coef(m, s = lambda)
    coef_vec <- as.vector(tmp)[-1]
    selected <- which(coef_vec != 0)
    if (length(selected) == 0) {
      all_pred_cv[test_idx] <- NA
      all_true_cv[test_idx] <- y_fold_test
      n_features_vec[i] <- 0
      next
    }
    
    n_features_vec[i] <- length(selected)
    
    X_fold_train_sel <- X_fold_train[, selected, drop = FALSE]
    df_train <- data.frame(outcome = y_fold_train, X_fold_train_sel)
    m_ols <- lm(outcome ~ ., data = df_train)
    
    X_fold_test_sel <- X_fold_test[, selected, drop = FALSE]
    df_test <- data.frame(X_fold_test_sel)
    pred_test <- predict(m_ols, newdata = df_test)
    
    all_pred_cv[test_idx] <- pred_test
    all_true_cv[test_idx] <- y_fold_test
  }
  
  valid_pred_cv <- !is.na(all_pred_cv)
  r_cv <- if (sum(valid_pred_cv) >= 10) {
    cor(all_true_cv[valid_pred_cv], all_pred_cv[valid_pred_cv], use = "complete.obs")
  } else NA_real_
  
  r_test <- NA_real_
  train_r2 <- NA_real_
  test_rmse <- NA_real_
  pred_test_full <- NULL
  m_full <- tryCatch(
    cv.glmnet(x = X_tr, y = y_tr, nfolds = 10),
    error = function(e) NULL
  )
  
  if (!is.null(m_full) && length(y_te) > 5) {
    lambda <- m_full$lambda.min
    j <- sum(coef(m_full, s = "lambda.min") != 0) - 1
    if (j < min_regressor) {
      idx <- which(m_full$nzero >= min_regressor)
      if (length(idx) > 0) lambda <- m_full$lambda[max(idx)]
    }
    tmp <- coef(m_full, s = lambda)
    coef_vec <- as.vector(tmp)[-1]
    selected <- which(coef_vec != 0)
    if (length(selected) > 0) {
      X_tr_sel <- X_tr[, selected, drop = FALSE]
      df_train <- data.frame(outcome = y_tr, X_tr_sel)
      m_ols <- lm(outcome ~ ., data = df_train)
      
      train_pred <- predict(m_ols, newdata = df_train)
      valid_train_full <- !is.na(y_tr) & !is.na(train_pred)
      if (sum(valid_train_full) >= 2) {
        train_r2 <- cor(train_pred[valid_train_full], y_tr[valid_train_full])^2
      }
      
      X_te_sel <- X_te[, selected, drop = FALSE]
      df_test <- data.frame(X_te_sel)
      pred_test_full <- predict(m_ols, newdata = df_test)
      
      if (sum(!is.na(pred_test_full)) >= 5) {
        r_test <- cor(y_te, pred_test_full, use = "complete.obs")
        test_rmse <- sqrt(mean((pred_test_full - y_te)^2, na.rm = TRUE))
      }
    }
  }
  
  return(list(
    r_cv = round(r_cv, 4),
    r_test = round(r_test, 4),
    n_features = round(mean(n_features_vec, na.rm = TRUE), 1),
    train_r2 = round(train_r2, 4),
    test_rmse = round(test_rmse, 4),
    pred_test = pred_test_full,
    y_test = y_te
  ))
}

# -- GBLUP model function ------------------------------------
run_gblup <- function(X_train, y_train, X_test, y_test, n_folds = 10, seed = 42) {
  valid_train <- !is.na(y_train)
  X_tr <- X_train[valid_train, , drop = FALSE]
  y_tr <- y_train[valid_train]
  n_tr <- length(y_tr)
  
  valid_test <- !is.na(y_test)
  X_te <- X_test[valid_test, , drop = FALSE]
  y_te <- y_test[valid_test]
  
  set.seed(seed)
  folds <- createFolds(y_tr, k = n_folds, returnTrain = FALSE)
  
  all_pred_cv <- numeric(n_tr)
  all_true_cv <- numeric(n_tr)
  
  for (i in 1:n_folds) {
    test_idx  <- folds[[i]]
    train_idx <- setdiff(1:n_tr, test_idx)
    
    X_fold_train <- X_tr[train_idx, , drop = FALSE]
    y_fold_train <- y_tr[train_idx]
    X_fold_test  <- X_tr[test_idx, , drop = FALSE]
    y_fold_test  <- y_tr[test_idx]
    
    X_scaled_train <- scale(X_fold_train, center = TRUE, scale = TRUE)
    X_scaled_train[is.nan(X_scaled_train)] <- 0
    train_center <- attr(X_scaled_train, "scaled:center")
    train_scale  <- attr(X_scaled_train, "scaled:scale")
    train_scale[train_scale == 0] <- 1
    
    G_train <- tcrossprod(X_scaled_train) / ncol(X_scaled_train)
    fit <- tryCatch(mixed.solve(y = y_fold_train, K = G_train), error = function(e) NULL)
    if (is.null(fit)) {
      all_pred_cv[test_idx] <- NA
      all_true_cv[test_idx] <- y_fold_test
      next
    }
    
    X_scaled_test <- sweep(X_fold_test, 2, train_center, "-")
    X_scaled_test <- sweep(X_scaled_test, 2, train_scale, "/")
    X_scaled_test[is.nan(X_scaled_test)] <- 0
    G_test_train <- tcrossprod(X_scaled_test, X_scaled_train) / ncol(X_scaled_train)
    G_train_inv <- tryCatch(solve(G_train + diag(1e-6, nrow(G_train))), error = function(e) NULL)
    if (is.null(G_train_inv)) {
      all_pred_cv[test_idx] <- NA
      all_true_cv[test_idx] <- y_fold_test
      next
    }
    
    pred_test <- as.numeric(fit$beta) + G_test_train %*% G_train_inv %*% fit$u
    all_pred_cv[test_idx] <- as.vector(pred_test)
    all_true_cv[test_idx] <- y_fold_test
  }
  
  valid_pred_cv <- !is.na(all_pred_cv)
  r_cv <- if (sum(valid_pred_cv) >= 10) {
    cor(all_true_cv[valid_pred_cv], all_pred_cv[valid_pred_cv], use = "complete.obs")
  } else NA_real_
  
  r_test <- NA_real_
  train_r2 <- NA_real_
  test_rmse <- NA_real_
  pred_test_full <- NULL
  if (length(y_te) > 5) {
    X_scaled_full <- scale(X_tr, center = TRUE, scale = TRUE)
    X_scaled_full[is.nan(X_scaled_full)] <- 0
    full_center <- attr(X_scaled_full, "scaled:center")
    full_scale  <- attr(X_scaled_full, "scaled:scale")
    full_scale[full_scale == 0] <- 1
    
    G_full <- tcrossprod(X_scaled_full) / ncol(X_scaled_full)
    fit_full <- tryCatch(mixed.solve(y = y_tr, K = G_full), error = function(e) NULL)
    if (!is.null(fit_full)) {
      G_full_inv <- solve(G_full + diag(1e-6, nrow(G_full)))
      train_pred <- as.numeric(fit_full$beta) + G_full %*% G_full_inv %*% fit_full$u
      valid_train_full <- !is.na(y_tr) & !is.na(train_pred)
      if (sum(valid_train_full) >= 2) {
        train_r2 <- cor(train_pred[valid_train_full], y_tr[valid_train_full])^2
      }
      
      X_te_scaled <- sweep(X_te, 2, full_center, "-")
      X_te_scaled <- sweep(X_te_scaled, 2, full_scale, "/")
      X_te_scaled[is.nan(X_te_scaled)] <- 0
      G_test_full <- tcrossprod(X_te_scaled, X_scaled_full) / ncol(X_scaled_full)
      pred_test_full <- as.numeric(fit_full$beta) + G_test_full %*% G_full_inv %*% fit_full$u
      valid_test_full <- !is.na(y_te) & !is.na(pred_test_full)
      if (sum(valid_test_full) >= 5) {
        r_test <- cor(pred_test_full[valid_test_full], y_te[valid_test_full], use = "complete.obs")
        test_rmse <- sqrt(mean((pred_test_full[valid_test_full] - y_te[valid_test_full])^2))
      }
    }
  }
  
  return(list(r_cv = round(r_cv, 4), r_test = round(r_test, 4),
              train_r2 = round(train_r2, 4), test_rmse = round(test_rmse, 4),
              pred_test = pred_test_full, y_test = y_te))
}

# -- GBLUP_AD model function ---------------------------------
run_gblup_ad <- function(X_train, y_train, X_test, y_test, n_folds = 10, seed = 42) {
  valid_train <- !is.na(y_train)
  X_tr <- X_train[valid_train, , drop = FALSE]
  y_tr <- y_train[valid_train]
  n_tr <- length(y_tr)
  D_tr <- compute_dominance_matrix(X_tr)
  
  valid_test <- !is.na(y_test)
  X_te <- X_test[valid_test, , drop = FALSE]
  y_te <- y_test[valid_test]
  D_te <- compute_dominance_matrix(X_te)
  
  set.seed(seed)
  folds <- createFolds(y_tr, k = n_folds, returnTrain = FALSE)
  
  all_pred_cv <- numeric(n_tr)
  all_true_cv <- numeric(n_tr)
  
  for (i in 1:n_folds) {
    test_idx  <- folds[[i]]
    train_idx <- setdiff(1:n_tr, test_idx)
    
    X_fold_train <- X_tr[train_idx, , drop = FALSE]
    D_fold_train <- D_tr[train_idx, , drop = FALSE]
    y_fold_train <- y_tr[train_idx]
    X_fold_test  <- X_tr[test_idx, , drop = FALSE]
    D_fold_test  <- D_tr[test_idx, , drop = FALSE]
    y_fold_test  <- y_tr[test_idx]
    
    X_scaled_train <- scale(X_fold_train, center = TRUE, scale = TRUE)
    X_scaled_train[is.nan(X_scaled_train)] <- 0
    center_A <- attr(X_scaled_train, "scaled:center")
    scale_A  <- attr(X_scaled_train, "scaled:scale")
    scale_A[scale_A == 0] <- 1
    
    D_scaled_train <- scale(D_fold_train, center = TRUE, scale = TRUE)
    D_scaled_train[is.nan(D_scaled_train)] <- 0
    center_D <- attr(D_scaled_train, "scaled:center")
    scale_D  <- attr(D_scaled_train, "scaled:scale")
    scale_D[scale_D == 0] <- 1
    
    G_A_train <- tcrossprod(X_scaled_train) / ncol(X_scaled_train)
    G_D_train <- tcrossprod(D_scaled_train) / ncol(D_scaled_train)
    G_train <- G_A_train + G_D_train
    
    fit <- tryCatch(mixed.solve(y = y_fold_train, K = G_train), error = function(e) NULL)
    if (is.null(fit)) {
      all_pred_cv[test_idx] <- NA
      all_true_cv[test_idx] <- y_fold_test
      next
    }
    
    X_scaled_test <- sweep(X_fold_test, 2, center_A, "-")
    X_scaled_test <- sweep(X_scaled_test, 2, scale_A, "/")
    X_scaled_test[is.nan(X_scaled_test)] <- 0
    D_scaled_test <- sweep(D_fold_test, 2, center_D, "-")
    D_scaled_test <- sweep(D_scaled_test, 2, scale_D, "/")
    D_scaled_test[is.nan(D_scaled_test)] <- 0
    
    G_A_test <- tcrossprod(X_scaled_test, X_scaled_train) / ncol(X_scaled_train)
    G_D_test <- tcrossprod(D_scaled_test, D_scaled_train) / ncol(D_scaled_train)
    G_test_train <- G_A_test + G_D_test
    
    G_train_inv <- tryCatch(solve(G_train + diag(1e-6, nrow(G_train))), error = function(e) NULL)
    if (is.null(G_train_inv)) {
      all_pred_cv[test_idx] <- NA
      all_true_cv[test_idx] <- y_fold_test
      next
    }
    
    pred_test <- as.numeric(fit$beta) + G_test_train %*% G_train_inv %*% fit$u
    all_pred_cv[test_idx] <- as.vector(pred_test)
    all_true_cv[test_idx] <- y_fold_test
  }
  
  valid_pred_cv <- !is.na(all_pred_cv)
  r_cv <- if (sum(valid_pred_cv) >= 10) {
    cor(all_true_cv[valid_pred_cv], all_pred_cv[valid_pred_cv], use = "complete.obs")
  } else NA_real_
  
  r_test <- NA_real_
  train_r2 <- NA_real_
  test_rmse <- NA_real_
  pred_test_full <- NULL
  if (length(y_te) > 5) {
    X_scaled_full <- scale(X_tr, center = TRUE, scale = TRUE)
    X_scaled_full[is.nan(X_scaled_full)] <- 0
    full_center_A <- attr(X_scaled_full, "scaled:center")
    full_scale_A  <- attr(X_scaled_full, "scaled:scale")
    full_scale_A[full_scale_A == 0] <- 1
    
    D_scaled_full <- scale(D_tr, center = TRUE, scale = TRUE)
    D_scaled_full[is.nan(D_scaled_full)] <- 0
    full_center_D <- attr(D_scaled_full, "scaled:center")
    full_scale_D  <- attr(D_scaled_full, "scaled:scale")
    full_scale_D[full_scale_D == 0] <- 1
    
    G_A_full <- tcrossprod(X_scaled_full) / ncol(X_scaled_full)
    G_D_full <- tcrossprod(D_scaled_full) / ncol(D_scaled_full)
    G_full <- G_A_full + G_D_full
    
    fit_full <- tryCatch(mixed.solve(y = y_tr, K = G_full), error = function(e) NULL)
    if (!is.null(fit_full)) {
      G_full_inv <- solve(G_full + diag(1e-6, nrow(G_full)))
      train_pred <- as.numeric(fit_full$beta) + G_full %*% G_full_inv %*% fit_full$u
      valid_train_full <- !is.na(y_tr) & !is.na(train_pred)
      if (sum(valid_train_full) >= 2) {
        train_r2 <- cor(train_pred[valid_train_full], y_tr[valid_train_full])^2
      }
      
      X_te_scaled <- sweep(X_te, 2, full_center_A, "-")
      X_te_scaled <- sweep(X_te_scaled, 2, full_scale_A, "/")
      X_te_scaled[is.nan(X_te_scaled)] <- 0
      D_te_scaled <- sweep(D_te, 2, full_center_D, "-")
      D_te_scaled <- sweep(D_te_scaled, 2, full_scale_D, "/")
      D_te_scaled[is.nan(D_te_scaled)] <- 0
      
      G_A_te <- tcrossprod(X_te_scaled, X_scaled_full) / ncol(X_scaled_full)
      G_D_te <- tcrossprod(D_te_scaled, D_scaled_full) / ncol(D_scaled_full)
      G_test_full <- G_A_te + G_D_te
      
      pred_test_full <- as.numeric(fit_full$beta) + G_test_full %*% G_full_inv %*% fit_full$u
      valid_test_full <- !is.na(y_te) & !is.na(pred_test_full)
      if (sum(valid_test_full) >= 5) {
        r_test <- cor(pred_test_full[valid_test_full], y_te[valid_test_full], use = "complete.obs")
        test_rmse <- sqrt(mean((pred_test_full[valid_test_full] - y_te[valid_test_full])^2))
      }
    }
  }
  
  return(list(r_cv = round(r_cv, 4), r_test = round(r_test, 4),
              train_r2 = round(train_r2, 4), test_rmse = round(test_rmse, 4),
              pred_test = pred_test_full, y_test = y_te))
}

# -- Prediction functions (For LOMO-CV) ----------------------
predict_lasso <- function(X_train, y_train, X_test, min_regressor = 1) {
  valid <- !is.na(y_train)
  X_tr <- X_train[valid, , drop = FALSE]
  y_tr <- y_train[valid]
  if (nrow(X_tr) < 2 || ncol(X_tr) == 0) return(rep(NA, nrow(X_test)))
  
  m <- tryCatch(cv.glmnet(x = X_tr, y = y_tr, nfolds = 10), error = function(e) NULL)
  if (is.null(m)) return(rep(NA, nrow(X_test)))
  
  lambda <- m$lambda.min
  j <- sum(coef(m, s = "lambda.min") != 0) - 1
  if (j < min_regressor) {
    idx <- which(m$nzero >= min_regressor)
    if (length(idx) == 0) return(rep(NA, nrow(X_test)))
    lambda <- m$lambda[max(idx)]
  }
  
  tmp <- coef(m, s = lambda)
  coef_vec <- as.vector(tmp)[-1]
  selected <- which(coef_vec != 0)
  if (length(selected) == 0) return(rep(NA, nrow(X_test)))
  
  X_tr_sel <- X_tr[, selected, drop = FALSE]
  df_train <- data.frame(outcome = y_tr, X_tr_sel)
  m_ols <- lm(outcome ~ ., data = df_train)
  
  X_te_sel <- X_test[, selected, drop = FALSE]
  df_test <- data.frame(X_te_sel)
  pred <- predict(m_ols, newdata = df_test)
  return(as.vector(pred))
}

predict_gblup <- function(X_train, y_train, X_test) {
  valid <- !is.na(y_train)
  X_tr <- X_train[valid, , drop = FALSE]
  y_tr <- y_train[valid]
  if (nrow(X_tr) < 2) return(rep(NA, nrow(X_test)))
  
  X_scaled <- scale(X_tr, center = TRUE, scale = TRUE)
  X_scaled[is.nan(X_scaled)] <- 0
  center <- attr(X_scaled, "scaled:center")
  scale_vec <- attr(X_scaled, "scaled:scale")
  scale_vec[scale_vec == 0] <- 1
  
  G <- tcrossprod(X_scaled) / ncol(X_scaled)
  fit <- tryCatch(mixed.solve(y = y_tr, K = G), error = function(e) NULL)
  if (is.null(fit)) return(rep(NA, nrow(X_test)))
  
  X_te_scaled <- sweep(X_test, 2, center, "-")
  X_te_scaled <- sweep(X_te_scaled, 2, scale_vec, "/")
  X_te_scaled[is.nan(X_te_scaled)] <- 0
  G_te <- tcrossprod(X_te_scaled, X_scaled) / ncol(X_scaled)
  G_inv <- tryCatch(solve(G + diag(1e-6, nrow(G))), error = function(e) NULL)
  if (is.null(G_inv)) return(rep(NA, nrow(X_test)))
  
  pred <- as.numeric(fit$beta) + G_te %*% G_inv %*% fit$u
  return(as.vector(pred))
}

predict_gblup_ad <- function(X_train, y_train, X_test) {
  valid <- !is.na(y_train)
  X_tr <- X_train[valid, , drop = FALSE]
  y_tr <- y_train[valid]
  if (nrow(X_tr) < 2) return(rep(NA, nrow(X_test)))
  
  D_tr <- compute_dominance_matrix(X_tr)
  D_te <- compute_dominance_matrix(X_test)
  
  X_scaled <- scale(X_tr, center = TRUE, scale = TRUE)
  X_scaled[is.nan(X_scaled)] <- 0
  center_A <- attr(X_scaled, "scaled:center")
  scale_A <- attr(X_scaled, "scaled:scale")
  scale_A[scale_A == 0] <- 1
  
  D_scaled <- scale(D_tr, center = TRUE, scale = TRUE)
  D_scaled[is.nan(D_scaled)] <- 0
  center_D <- attr(D_scaled, "scaled:center")
  scale_D <- attr(D_scaled, "scaled:scale")
  scale_D[scale_D == 0] <- 1
  
  G_A <- tcrossprod(X_scaled) / ncol(X_scaled)
  G_D <- tcrossprod(D_scaled) / ncol(D_scaled)
  G <- G_A + G_D
  
  fit <- tryCatch(mixed.solve(y = y_tr, K = G), error = function(e) NULL)
  if (is.null(fit)) return(rep(NA, nrow(X_test)))
  
  X_te_scaled <- sweep(X_test, 2, center_A, "-")
  X_te_scaled <- sweep(X_te_scaled, 2, scale_A, "/")
  X_te_scaled[is.nan(X_te_scaled)] <- 0
  D_te_scaled <- sweep(D_te, 2, center_D, "-")
  D_te_scaled <- sweep(D_te_scaled, 2, scale_D, "/")
  D_te_scaled[is.nan(D_te_scaled)] <- 0
  
  G_A_te <- tcrossprod(X_te_scaled, X_scaled) / ncol(X_scaled)
  G_D_te <- tcrossprod(D_te_scaled, D_scaled) / ncol(D_scaled)
  G_te <- G_A_te + G_D_te
  
  G_inv <- tryCatch(solve(G + diag(1e-6, nrow(G))), error = function(e) NULL)
  if (is.null(G_inv)) return(rep(NA, nrow(X_test)))
  
  pred <- as.numeric(fit$beta) + G_te %*% G_inv %*% fit$u
  return(as.vector(pred))
}

# -- Single seed processing function -------------------------
run_single_seed <- function(seed, snp_cols_all, pheno_train, pheno_test,
                            TRAITS, SNP_RATIO, N_FOLDS, MIN_REGRESSOR) {
  
  set.seed(seed)
  n_select <- round(length(snp_cols_all) * SNP_RATIO)
  selected_snps <- sample(snp_cols_all, size = n_select)
  cols_to_select <- c(meta_cols, selected_snps)
  
  snp_subset <- vroom(SNP_FILE, 
                      col_select = all_of(cols_to_select),
                      show_col_types = FALSE,
                      progress = FALSE) %>%
    as.data.table()
  
  snp_subset[, hybrid_id := clean_id(hybrid_id)]
  
  train_ids <- intersect(snp_subset$hybrid_id, pheno_train$hybrid_id)
  snp_train <- snp_subset[hybrid_id %in% train_ids]
  pheno_train_sub <- pheno_train[hybrid_id %in% train_ids]
  setkey(snp_train, hybrid_id)
  setkey(pheno_train_sub, hybrid_id)
  snp_train <- snp_train[pheno_train_sub$hybrid_id]
  X_train <- as.matrix(snp_train[, ..selected_snps])
  X_train[is.na(X_train)] <- 0
  
  cms_train <- snp_train$cms
  restorer_train <- pheno_train_sub$restorer
  
  test_ids <- intersect(snp_subset$hybrid_id, pheno_test$hybrid_id)
  if (length(test_ids) > 0) {
    snp_test <- snp_subset[hybrid_id %in% test_ids]
    pheno_test_sub <- pheno_test[hybrid_id %in% test_ids]
    setkey(snp_test, hybrid_id)
    setkey(pheno_test_sub, hybrid_id)
    snp_test <- snp_test[pheno_test_sub$hybrid_id]
    X_test <- as.matrix(snp_test[, ..selected_snps])
    X_test[is.na(X_test)] <- 0
  } else {
    X_test <- NULL
    pheno_test_sub <- NULL
  }
  
  results <- list()
  
  for (trait in TRAITS) {
    y_train <- pheno_train_sub[[trait]]
    y_test  <- if (!is.null(pheno_test_sub)) pheno_test_sub[[trait]] else NULL
    
    lasso_res <- tryCatch(
      run_lasso_ols(X_train, y_train, X_test, y_test,
                    n_folds = N_FOLDS, min_regressor = MIN_REGRESSOR, seed = seed),
      error = function(e) list(r_cv = NA_real_, r_test = NA_real_, n_features = 0,
                               train_r2 = NA_real_, test_rmse = NA_real_,
                               pred_test = NULL, y_test = NULL)
    )
    
    gblup_res <- tryCatch(
      run_gblup(X_train, y_train, X_test, y_test, n_folds = N_FOLDS, seed = seed),
      error = function(e) list(r_cv = NA_real_, r_test = NA_real_,
                               train_r2 = NA_real_, test_rmse = NA_real_,
                               pred_test = NULL, y_test = NULL)
    )
    
    gblup_ad_res <- tryCatch(
      run_gblup_ad(X_train, y_train, X_test, y_test, n_folds = N_FOLDS, seed = seed),
      error = function(e) list(r_cv = NA_real_, r_test = NA_real_,
                               train_r2 = NA_real_, test_rmse = NA_real_,
                               pred_test = NULL, y_test = NULL)
    )
    
    unique_cms <- unique(cms_train)
    pred_lasso_all <- rep(NA, length(y_train))
    pred_gblup_all <- rep(NA, length(y_train))
    pred_gblup_ad_all <- rep(NA, length(y_train))
    
    for (cc in unique_cms) {
      val_idx <- which(cms_train == cc)
      train_idx <- setdiff(seq_along(y_train), val_idx)
      if (length(train_idx) < 2 || length(val_idx) == 0) next
      
      X_tr <- X_train[train_idx, , drop = FALSE]
      y_tr <- y_train[train_idx]
      X_val <- X_train[val_idx, , drop = FALSE]
      
      pred_lasso <- tryCatch(predict_lasso(X_tr, y_tr, X_val, min_regressor = MIN_REGRESSOR),
                             error = function(e) rep(NA, length(val_idx)))
      pred_gblup <- tryCatch(predict_gblup(X_tr, y_tr, X_val),
                             error = function(e) rep(NA, length(val_idx)))
      pred_gblup_ad <- tryCatch(predict_gblup_ad(X_tr, y_tr, X_val),
                                error = function(e) rep(NA, length(val_idx)))
      
      pred_lasso_all[val_idx] <- pred_lasso
      pred_gblup_all[val_idx] <- pred_gblup
      pred_gblup_ad_all[val_idx] <- pred_gblup_ad
    }
    
    valid_lasso <- !is.na(pred_lasso_all)
    valid_gblup <- !is.na(pred_gblup_all)
    valid_gblup_ad <- !is.na(pred_gblup_ad_all)
    
    r_lomo_lasso <- if (sum(valid_lasso) >= 10) cor(y_train[valid_lasso], pred_lasso_all[valid_lasso], use = "complete.obs") else NA_real_
    r_lomo_gblup <- if (sum(valid_gblup) >= 10) cor(y_train[valid_gblup], pred_gblup_all[valid_gblup], use = "complete.obs") else NA_real_
    r_lomo_gblup_ad <- if (sum(valid_gblup_ad) >= 10) cor(y_train[valid_gblup_ad], pred_gblup_ad_all[valid_gblup_ad], use = "complete.obs") else NA_real_
    
    results[[length(results) + 1]] <- data.frame(
      seed = seed,
      trait = trait,
      LASSO_cv = lasso_res$r_cv,
      LASSO_test = lasso_res$r_test,
      LASSO_lomo = r_lomo_lasso,
      LASSO_trainR2 = lasso_res$train_r2,
      LASSO_testRMSE = lasso_res$test_rmse,
      GBLUP_cv = gblup_res$r_cv,
      GBLUP_test = gblup_res$r_test,
      GBLUP_lomo = r_lomo_gblup,
      GBLUP_trainR2 = gblup_res$train_r2,
      GBLUP_testRMSE = gblup_res$test_rmse,
      GBLUP_AD_cv = gblup_ad_res$r_cv,
      GBLUP_AD_test = gblup_ad_res$r_test,
      GBLUP_AD_lomo = r_lomo_gblup_ad,
      GBLUP_AD_trainR2 = gblup_ad_res$train_r2,
      GBLUP_AD_testRMSE = gblup_ad_res$test_rmse,
      n_features = lasso_res$n_features,
      stringsAsFactors = FALSE
    )
  }
  return(do.call(rbind, results))
}

# -- Main loop: Evaluate 3 models across 50 seeds -------------
cat("== SNP Random Subsampling Sensitivity Analysis (LOMO-CV) ==\n")
cat("Method: Random 5% SNP sampling per seed\n")
cat("Validation: 10-fold CV + LOMO-CV + Independent test set (n=92)\n")
cat("Random seeds: 101-150 (50 iterations)\n")
cat("Date:", run_date, "\n")
cat("Output directory:", output_dir, "\n\n")

SEEDS <- 101:150

start_time <- Sys.time()
results_list <- list()
for (i in seq_along(SEEDS)) {
  seed <- SEEDS[i]
  cat("Processing seed", i, "/", length(SEEDS), " (seed =", seed, ")\n")
  res <- run_single_seed(seed, snp_cols_all, pheno_train, pheno_test,
                         TRAITS, SNP_RATIO, N_FOLDS, MIN_REGRESSOR)
  results_list[[i]] <- as.data.table(res)
  gc()
}
end_time <- Sys.time()
cat("\nComputation completed! Total time:", round(difftime(end_time, start_time, units = "mins"), 2), "mins\n\n")

# -- Summarize results ---------------------------------------
results <- rbindlist(results_list, fill = TRUE)

summary_dt <- results[, .(
  LASSO_cv_mean = round(mean(LASSO_cv, na.rm = TRUE), 4),
  LASSO_cv_sd = round(sd(LASSO_cv, na.rm = TRUE), 4),
  LASSO_test_mean = round(mean(LASSO_test, na.rm = TRUE), 4),
  LASSO_test_sd = round(sd(LASSO_test, na.rm = TRUE), 4),
  LASSO_lomo_mean = round(mean(LASSO_lomo, na.rm = TRUE), 4),
  LASSO_lomo_sd = round(sd(LASSO_lomo, na.rm = TRUE), 4),
  LASSO_trainR2_mean = round(mean(LASSO_trainR2, na.rm = TRUE), 4),
  LASSO_trainR2_sd = round(sd(LASSO_trainR2, na.rm = TRUE), 4),
  LASSO_testRMSE_mean = round(mean(LASSO_testRMSE, na.rm = TRUE), 4),
  LASSO_testRMSE_sd = round(sd(LASSO_testRMSE, na.rm = TRUE), 4),
  GBLUP_cv_mean = round(mean(GBLUP_cv, na.rm = TRUE), 4),
  GBLUP_cv_sd = round(sd(GBLUP_cv, na.rm = TRUE), 4),
  GBLUP_test_mean = round(mean(GBLUP_test, na.rm = TRUE), 4),
  GBLUP_test_sd = round(sd(GBLUP_test, na.rm = TRUE), 4),
  GBLUP_lomo_mean = round(mean(GBLUP_lomo, na.rm = TRUE), 4),
  GBLUP_lomo_sd = round(sd(GBLUP_lomo, na.rm = TRUE), 4),
  GBLUP_trainR2_mean = round(mean(GBLUP_trainR2, na.rm = TRUE), 4),
  GBLUP_trainR2_sd = round(sd(GBLUP_trainR2, na.rm = TRUE), 4),
  GBLUP_testRMSE_mean = round(mean(GBLUP_testRMSE, na.rm = TRUE), 4),
  GBLUP_testRMSE_sd = round(sd(GBLUP_testRMSE, na.rm = TRUE), 4),
  GBLUP_AD_cv_mean = round(mean(GBLUP_AD_cv, na.rm = TRUE), 4),
  GBLUP_AD_cv_sd = round(sd(GBLUP_AD_cv, na.rm = TRUE), 4),
  GBLUP_AD_test_mean = round(mean(GBLUP_AD_test, na.rm = TRUE), 4),
  GBLUP_AD_test_sd = round(sd(GBLUP_AD_test, na.rm = TRUE), 4),
  GBLUP_AD_lomo_mean = round(mean(GBLUP_AD_lomo, na.rm = TRUE), 4),
  GBLUP_AD_lomo_sd = round(sd(GBLUP_AD_lomo, na.rm = TRUE), 4),
  GBLUP_AD_trainR2_mean = round(mean(GBLUP_AD_trainR2, na.rm = TRUE), 4),
  GBLUP_AD_trainR2_sd = round(sd(GBLUP_AD_trainR2, na.rm = TRUE), 4),
  GBLUP_AD_testRMSE_mean = round(mean(GBLUP_AD_testRMSE, na.rm = TRUE), 4),
  GBLUP_AD_testRMSE_sd = round(sd(GBLUP_AD_testRMSE, na.rm = TRUE), 4),
  n_features_mean = round(mean(n_features, na.rm = TRUE), 1),
  n_seeds = .N
), by = trait]

cat(strrep("=", 80), "\n")
cat("  Summary Results (", length(SEEDS), " seeds)\n", sep="")
cat(strrep("=", 80), "\n\n")

cat("[10-fold Cross-Validation Results (r_cv)]\n")
cv_table <- summary_dt[, .(
  trait,
  `LASSO_SNP` = sprintf("%.3f±%.3f", LASSO_cv_mean, LASSO_cv_sd),
  `GBLUP` = sprintf("%.3f±%.3f", GBLUP_cv_mean, GBLUP_cv_sd),
  `GBLUP_AD` = sprintf("%.3f±%.3f", GBLUP_AD_cv_mean, GBLUP_AD_cv_sd)
)]
print(cv_table)

cat("\n[LOMO-CV Results (r_lomo)]\n")
lomo_table <- summary_dt[, .(
  trait,
  `LASSO_SNP` = sprintf("%.3f±%.3f", LASSO_lomo_mean, LASSO_lomo_sd),
  `GBLUP` = sprintf("%.3f±%.3f", GBLUP_lomo_mean, GBLUP_lomo_sd),
  `GBLUP_AD` = sprintf("%.3f±%.3f", GBLUP_AD_lomo_mean, GBLUP_AD_lomo_sd)
)]
print(lomo_table)

cat("\n[Independent Test Set Results (r_test, n=92)]\n")
test_table <- summary_dt[, .(
  trait,
  `LASSO_SNP` = sprintf("%.3f±%.3f", LASSO_test_mean, LASSO_test_sd),
  `GBLUP` = sprintf("%.3f±%.3f", GBLUP_test_mean, GBLUP_test_sd),
  `GBLUP_AD` = sprintf("%.3f±%.3f", GBLUP_AD_test_mean, GBLUP_AD_test_sd)
)]
print(test_table)

summary_wide <- data.table(
  Model = c("LASSO_SNP", "GBLUP", "GBLUP_AD")
)
for (tr in TRAITS) {
  sub_dt <- summary_dt[trait == tr]
  summary_wide[[paste0(tr, "_cv")]] <- c(
    sprintf("%.3f±%.3f", sub_dt$LASSO_cv_mean, sub_dt$LASSO_cv_sd),
    sprintf("%.3f±%.3f", sub_dt$GBLUP_cv_mean, sub_dt$GBLUP_cv_sd),
    sprintf("%.3f±%.3f", sub_dt$GBLUP_AD_cv_mean, sub_dt$GBLUP_AD_cv_sd)
  )
  summary_wide[[paste0(tr, "_lomo")]] <- c(
    sprintf("%.3f±%.3f", sub_dt$LASSO_lomo_mean, sub_dt$LASSO_lomo_sd),
    sprintf("%.3f±%.3f", sub_dt$GBLUP_lomo_mean, sub_dt$GBLUP_lomo_sd),
    sprintf("%.3f±%.3f", sub_dt$GBLUP_AD_lomo_mean, sub_dt$GBLUP_AD_lomo_sd)
  )
  summary_wide[[paste0(tr, "_test")]] <- c(
    sprintf("%.3f±%.3f", sub_dt$LASSO_test_mean, sub_dt$LASSO_test_sd),
    sprintf("%.3f±%.3f", sub_dt$GBLUP_test_mean, sub_dt$GBLUP_test_sd),
    sprintf("%.3f±%.3f", sub_dt$GBLUP_AD_test_mean, sub_dt$GBLUP_AD_test_sd)
  )
}
cat("\n[Publication Format Table (CV + LOMO + External Test)]\n")
print(summary_wide)

fwrite(results, file.path(output_dir, paste0("snp_benchmark_raw_", run_date, ".csv")))
fwrite(summary_dt, file.path(output_dir, paste0("snp_benchmark_summary_", run_date, ".csv")))
fwrite(cv_table, file.path(output_dir, paste0("snp_benchmark_cv_table_", run_date, ".csv")))
fwrite(lomo_table, file.path(output_dir, paste0("snp_benchmark_lomo_table_", run_date, ".csv")))
fwrite(test_table, file.path(output_dir, paste0("snp_benchmark_test_table_", run_date, ".csv")))
fwrite(summary_wide, file.path(output_dir, paste0("snp_benchmark_wide_", run_date, ".csv")))

cat("\nResults saved to:", output_dir, "\n")
cat("\n== Analysis completed ==\n")