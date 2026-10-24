# ============================================================
# Script: PGSI-based LASSO models (Additive, Dominance, Combined)
# Features: 10-fold outer CV, LOMO/LORO CV, external test, and 
#           100-iteration Bootstrap stability analysis.
# ============================================================

library(data.table)
library(glmnet)

run_date <- Sys.Date()

# -- Parameter settings --------------------------------------
USE_FEATURE_SCREENING <- TRUE
USE_LAMBDA_1SE        <- FALSE
MIN_REGRESSOR         <- 1
DO_BOOTSTRAP          <- TRUE
BOOTSTRAP_REPS        <- 100
ROBUST_THRESHOLD      <- 0.7
TOP_PCT               <- 0.1
BOTTOM_PCT            <- 0.1

modes <- c("add_only", "dom_only", "add_dom")

pheno_train <- fread("phenotype_training_336.csv")
pheno_test  <- fread("phenotype_test_92.csv")

required_cols <- c("hybrid_id", "restorer", "FLW_HPH", "FSTW_HPH", "FSHW_HPH")
if (!all(required_cols %in% names(pheno_train))) stop("Missing columns in training phenotype")
if (!all(required_cols %in% names(pheno_test))) stop("Missing columns in test phenotype")

clean_id <- function(x) gsub("[^A-Za-z0-9_]", "", trimws(as.character(x)))
pheno_train[, hybrid_id := clean_id(hybrid_id)]
pheno_test[, hybrid_id := clean_id(hybrid_id)]

# Note: Ensure these names match your phenotype CSV exactly. 
# If you have updated them to Ogu-R3 and Ogu-R8, change them here.
test_restorers <- c("R53", "R58") 
seeds <- 101:150
repeats <- 50
k_folds <- 10
nfolds_inner <- 10
traits <- c("FLW_HPH", "FSTW_HPH", "FSHW_HPH")

# -- Helper functions ----------------------------------------
predict_two_step <- function(X_tr, y_tr, X_te, min_nz = 1, use_lambda_1se = FALSE) {
  keep <- apply(X_tr, 2, var) > 0
  if (sum(keep) == 0) return(rep(mean(y_tr, na.rm = TRUE), nrow(X_te)))
  X_tr <- X_tr[, keep, drop = FALSE]
  X_te <- X_te[, keep, drop = FALSE]
  fit <- tryCatch(cv.glmnet(X_tr, y_tr, alpha = 1, nfolds = nfolds_inner), error = function(e) NULL)
  if (is.null(fit)) return(rep(mean(y_tr, na.rm = TRUE), nrow(X_te)))
  lambda_use <- if (use_lambda_1se) fit$lambda.1se else fit$lambda.min
  nz <- sum(coef(fit, s = lambda_use) != 0)
  if (nz < min_nz) {
    idx <- which(fit$nzero >= min_nz)
    if (length(idx) > 0) lambda_use <- fit$lambda[max(idx)]
  }
  coefs <- coef(fit, s = lambda_use)
  selected <- rownames(coefs)[which(coefs != 0)][-1]
  if (length(selected) == 0) return(rep(mean(y_tr, na.rm = TRUE), nrow(X_te)))
  df_tr <- data.frame(y = y_tr, X_tr[, selected, drop = FALSE])
  lm_fit <- tryCatch(lm(y ~ ., data = df_tr), error = function(e) NULL)
  if (is.null(lm_fit)) return(rep(mean(y_tr, na.rm = TRUE), nrow(X_te)))
  df_te <- as.data.frame(X_te[, selected, drop = FALSE])
  as.vector(predict(lm_fit, newdata = df_te))
}

calc_train_r2 <- function(X, y, min_nz = 1, use_lambda_1se = FALSE) {
  keep <- apply(X, 2, var) > 0
  if (sum(keep) == 0) return(0)
  X <- X[, keep, drop = FALSE]
  fit <- tryCatch(cv.glmnet(X, y, alpha = 1, nfolds = nfolds_inner), error = function(e) NULL)
  if (is.null(fit)) return(0)
  lambda_use <- if (use_lambda_1se) fit$lambda.1se else fit$lambda.min
  nz <- sum(coef(fit, s = lambda_use) != 0)
  if (nz < min_nz) {
    idx <- which(fit$nzero >= min_nz)
    if (length(idx) > 0) lambda_use <- fit$lambda[max(idx)]
  }
  coefs <- coef(fit, s = lambda_use)
  selected <- rownames(coefs)[which(coefs != 0)][-1]
  if (length(selected) == 0) return(0)
  df <- data.frame(y = y, X[, selected, drop = FALSE])
  lm_fit <- tryCatch(lm(y ~ ., data = df), error = function(e) NULL)
  if (is.null(lm_fit)) return(0)
  summary(lm_fit)$r.squared
}

get_lasso_coefs <- function(X, y, min_nz = 1, use_lambda_1se = FALSE, feature_names = NULL) {
  keep <- apply(X, 2, var) > 0
  if (sum(keep) == 0) return(NULL)
  X <- X[, keep, drop = FALSE]
  if (!is.null(feature_names)) feature_names <- feature_names[keep]
  fit <- tryCatch(cv.glmnet(X, y, alpha = 1, nfolds = nfolds_inner), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  lambda_use <- if (use_lambda_1se) fit$lambda.1se else fit$lambda.min
  nz <- sum(coef(fit, s = lambda_use) != 0)
  if (nz < min_nz) {
    idx <- which(fit$nzero >= min_nz)
    if (length(idx) > 0) lambda_use <- fit$lambda[max(idx)]
  }
  coef_vec <- as.vector(coef(fit, s = lambda_use))[-1]
  selected <- which(coef_vec != 0)
  if (length(selected) == 0) return(NULL)
  data.table(feature = if (is.null(feature_names)) colnames(X)[selected] else feature_names[selected],
             coefficient = coef_vec[selected])
}

read_pgsi_mode <- function(file_add, file_dom, mode) {
  if (mode == "add_only") {
    pgsi <- fread(file_add)
  } else if (mode == "dom_only") {
    pgsi <- fread(file_dom)
    setnames(pgsi, gsub("_dom$", "", names(pgsi)))
  } else { # add_dom
    pgsi_add <- fread(file_add)
    pgsi_dom <- fread(file_dom)
    setkey(pgsi_add, hybrid_id); setkey(pgsi_dom, hybrid_id)
    pgsi <- merge(pgsi_add, pgsi_dom, by = c("hybrid_id", "cms", "restorer"), suffixes = c("", "_dom"))
  }
  pgsi[, hybrid_id := clean_id(hybrid_id)]
  pgsi
}

get_important_windows <- function(X_train, y_train, use_lambda_1se, min_regressor) {
  fit <- cv.glmnet(X_train, y_train, alpha = 1, nfolds = nfolds_inner)
  lambda <- if (use_lambda_1se) fit$lambda.1se else fit$lambda.min
  nz <- sum(coef(fit, s = lambda) != 0)
  if (nz < min_regressor) {
    idx <- which(fit$nzero >= min_regressor)
    if (length(idx) > 0) lambda <- fit$lambda[max(idx)]
  }
  which(coef(fit, s = lambda)[-1] != 0)
}

window_configs <- list(
  list(name = "100Kb", add_file = "/biodata/zhy/pgsi_ad/100kb_simi_add.csv",
       dom_file = "/biodata/zhy/pgsi_ad/100kb_simi_dom.csv"),
  list(name = "500Kb", add_file = "/biodata/zhy/pgsi_ad/500kb_simi_add.csv",
       dom_file = "/biodata/zhy/pgsi_ad/500kb_simi_dom.csv"),
  list(name = "1Mb",   add_file = "/biodata/zhy/pgsi_ad/1000kb_simi_add.csv",
       dom_file = "/biodata/zhy/pgsi_ad/1000kb_simi_dom.csv")
)

# -- Main Loop -----------------------------------------------
for (mode in modes) {
  output_dir <- paste0("PGSI_LOMO_LORO_results_", mode, "_", run_date)
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  all_results_total <- list()
  
  for (win in window_configs) {
    window_name <- win$name
    add_file <- win$add_file
    dom_file <- win$dom_file
    
    pgsi <- read_pgsi_mode(add_file, dom_file, mode)
    
    train_restorers <- setdiff(unique(pgsi$restorer), test_restorers)
    pgsi_train <- pgsi[restorer %in% train_restorers]
    pgsi_test  <- pgsi[restorer %in% test_restorers]
    
    matched_train <- intersect(pgsi_train$hybrid_id, pheno_train$hybrid_id)
    pgsi_tr <- pgsi_train[hybrid_id %in% matched_train]
    pheno_tr <- pheno_train[hybrid_id %in% matched_train]
    setkey(pgsi_tr, hybrid_id); setkey(pheno_tr, hybrid_id)
    pgsi_tr <- pgsi_tr[pheno_tr[, .(hybrid_id)]]
    
    matched_test <- intersect(pgsi_test$hybrid_id, pheno_test$hybrid_id)
    if (length(matched_test) > 0) {
      pgsi_te <- pgsi_test[hybrid_id %in% matched_test]
      pheno_te <- pheno_test[hybrid_id %in% matched_test]
      setkey(pgsi_te, hybrid_id); setkey(pheno_te, hybrid_id)
      pgsi_te <- pgsi_te[pheno_te[, .(hybrid_id)]]
    } else {
      pgsi_te <- NULL; pheno_te <- NULL
    }
    
    feat_cols <- setdiff(names(pgsi_tr), c("hybrid_id", "cms", "restorer"))
    X_train_full <- as.matrix(pgsi_tr[, ..feat_cols]); X_train_full[is.na(X_train_full)] <- 0
    cms_vec <- pgsi_tr$cms
    restorer_vec <- pgsi_tr$restorer
    
    if (!is.null(pgsi_te)) {
      X_test_full <- as.matrix(pgsi_te[, ..feat_cols]); X_test_full[is.na(X_test_full)] <- 0
    } else { X_test_full <- NULL }
    
    # Feature screening (only on training set)
    important_windows_list <- list()
    if (USE_FEATURE_SCREENING) {
      for (trait in traits) {
        y_full <- pheno_tr[[trait]]
        valid <- !is.na(y_full)
        X_valid <- X_train_full[valid, , drop = FALSE]
        y_valid <- y_full[valid]
        if (nrow(X_valid) > 10) {
          imp <- get_important_windows(X_valid, y_valid, USE_LAMBDA_1SE, MIN_REGRESSOR)
          important_windows_list[[trait]] <- imp
        } else {
          important_windows_list[[trait]] <- seq_len(ncol(X_train_full))
        }
      }
    }
    
    window_results <- list()
    for (trait in traits) {
      y_train <- pheno_tr[[trait]]
      valid_train <- !is.na(y_train)
      X_train <- X_train_full[valid_train, , drop = FALSE]
      y_train <- y_train[valid_train]
      cms_train <- cms_vec[valid_train]
      restorer_train <- restorer_vec[valid_train]
      
      if (USE_FEATURE_SCREENING && length(important_windows_list[[trait]]) > 0) {
        X_train <- X_train[, important_windows_list[[trait]], drop = FALSE]
      }
      
      if (!is.null(X_test_full)) {
        y_test <- pheno_te[[trait]]
        valid_test <- !is.na(y_test)
        X_test <- X_test_full[valid_test, , drop = FALSE]
        y_test <- y_test[valid_test]
        if (USE_FEATURE_SCREENING && length(important_windows_list[[trait]]) > 0) {
          X_test <- X_test[, important_windows_list[[trait]], drop = FALSE]
        }
      } else {
        X_test <- NULL; y_test <- NULL
      }
      
      # 1. Random 10-fold CV (50 seeds)
      cv_cors <- numeric(repeats)
      train_r2 <- numeric(repeats)
      for (r in 1:repeats) {
        seed <- seeds[r]
        set.seed(seed)
        folds <- sample(rep(1:k_folds, length.out = nrow(X_train)))
        preds <- numeric(nrow(X_train))
        for (fold in 1:k_folds) {
          tr <- folds != fold
          te <- folds == fold
          preds[te] <- predict_two_step(X_train[tr, , drop = FALSE], y_train[tr],
                                        X_train[te, , drop = FALSE],
                                        MIN_REGRESSOR, USE_LAMBDA_1SE)
        }
        valid_pred <- !is.na(y_train) & !is.na(preds)
        cv_cors[r] <- cor(preds[valid_pred], y_train[valid_pred], use = "complete.obs")
        train_r2[r] <- calc_train_r2(X_train, y_train, MIN_REGRESSOR, USE_LAMBDA_1SE)
        
        coefs <- get_lasso_coefs(X_train, y_train, MIN_REGRESSOR, USE_LAMBDA_1SE, colnames(X_train))
        if (!is.null(coefs) && nrow(coefs) > 0) {
          coef_file <- file.path(output_dir, paste0("lasso_coefs_", window_name, "_", trait, "_seed", seed, "_", run_date, ".csv"))
          fwrite(coefs, coef_file)
        }
      }
      mean_cv <- mean(cv_cors, na.rm = TRUE)
      mean_train_r2 <- mean(train_r2, na.rm = TRUE)
      
      # 2. LOMO-CV
      cms_levels <- unique(cms_train)
      lomo_preds <- numeric(nrow(X_train))
      for (cc in cms_levels) {
        te <- cms_train == cc
        tr <- !te
        if (sum(tr) == 0) next
        lomo_preds[te] <- predict_two_step(X_train[tr, , drop = FALSE], y_train[tr],
                                           X_train[te, , drop = FALSE],
                                           MIN_REGRESSOR, USE_LAMBDA_1SE)
      }
      valid_lomo <- !is.na(y_train) & !is.na(lomo_preds)
      r_lomo <- cor(lomo_preds[valid_lomo], y_train[valid_lomo], use = "complete.obs")
      
      # 3. LORO-CV
      restorer_levels <- unique(restorer_train)
      loro_preds <- numeric(nrow(X_train))
      for (rr in restorer_levels) {
        te <- restorer_train == rr
        tr <- !te
        if (sum(tr) == 0) next
        loro_preds[te] <- predict_two_step(X_train[tr, , drop = FALSE], y_train[tr],
                                           X_train[te, , drop = FALSE],
                                           MIN_REGRESSOR, USE_LAMBDA_1SE)
      }
      valid_loro <- !is.na(y_train) & !is.na(loro_preds)
      r_loro <- cor(loro_preds[valid_loro], y_train[valid_loro], use = "complete.obs")
      
      # 4. External test
      if (!is.null(X_test) && length(y_test) >= 3) {
        pred_test <- predict_two_step(X_train, y_train, X_test, MIN_REGRESSOR, USE_LAMBDA_1SE)
        valid_ext <- !is.na(y_test) & !is.na(pred_test)
        r_ext <- cor(pred_test[valid_ext], y_test[valid_ext], use = "complete.obs")
        pred_df <- data.table(hybrid_id = pheno_te$hybrid_id[valid_test],
                              restorer = pheno_te$restorer[valid_test],
                              true_value = y_test, pred_value = pred_test)
        fwrite(pred_df, file.path(output_dir, paste0("test_predictions_", window_name, "_", trait, "_", run_date, ".csv")))
      } else { r_ext <- NA }
      
      # Save results
      res <- data.table(window = window_name, trait = trait,
                        random_cv_r = mean_cv, train_r2 = mean_train_r2,
                        lomo_cv_r = r_lomo, loro_cv_r = r_loro, external_r = r_ext,
                        n_train = nrow(X_train), n_test = ifelse(is.null(X_test), 0, nrow(X_test)))
      window_results[[trait]] <- res
      
      detail <- data.table(window = window_name, trait = trait, seed = seeds[1:length(cv_cors)],
                           cv_r = cv_cors, train_r2 = train_r2)
      fwrite(detail, file.path(output_dir, paste0("lasso_detail_", window_name, "_", trait, "_", run_date, ".csv")))
      fwrite(res, file.path(output_dir, paste0("lasso_summary_", window_name, "_", trait, "_", run_date, ".csv")))
    }
    if (length(window_results) > 0) {
      window_summary <- rbindlist(window_results)
      fwrite(window_summary, file.path(output_dir, paste0("validation_", window_name, "_", run_date, ".csv")))
      all_results_total[[window_name]] <- window_summary
    }
  }
  
  final_all <- rbindlist(all_results_total)
  fwrite(final_all, file.path(output_dir, paste0("validation_all_windows_", run_date, ".csv")))
  
  # ============================================================
  # 5. Bootstrap Stability Analysis (100 iterations)
  # ============================================================
  if (DO_BOOTSTRAP) {
    for (win in window_configs) {
      window_name <- win$name
      add_file <- win$add_file
      dom_file <- win$dom_file
      pgsi <- read_pgsi_mode(add_file, dom_file, mode)
      
      train_restorers <- setdiff(unique(pgsi$restorer), test_restorers)
      pgsi_train <- pgsi[restorer %in% train_restorers]
      matched_train <- intersect(pgsi_train$hybrid_id, pheno_train$hybrid_id)
      pgsi_tr <- pgsi_train[hybrid_id %in% matched_train]
      pheno_tr <- pheno_train[hybrid_id %in% matched_train]
      setkey(pgsi_tr, hybrid_id); setkey(pheno_tr, hybrid_id)
      pgsi_tr <- pgsi_tr[pheno_tr[, .(hybrid_id)]]
      
      feat_cols <- setdiff(names(pgsi_tr), c("hybrid_id", "cms", "restorer"))
      X_full <- as.matrix(pgsi_tr[, ..feat_cols]); X_full[is.na(X_full)] <- 0
      
      important_windows_list <- list()
      if (USE_FEATURE_SCREENING) {
        for (trait in traits) {
          y_full <- pheno_tr[[trait]]
          valid <- !is.na(y_full)
          X_valid <- X_full[valid, , drop = FALSE]
          y_valid <- y_full[valid]
          if (nrow(X_valid) > 10) {
            imp <- get_important_windows(X_valid, y_valid, USE_LAMBDA_1SE, MIN_REGRESSOR)
            important_windows_list[[trait]] <- imp
          } else {
            important_windows_list[[trait]] <- seq_len(ncol(X_full))
          }
        }
      }
      
      for (trait in traits) {
        y_full <- pheno_tr[[trait]]
        valid <- !is.na(y_full)
        X <- X_full[valid, , drop = FALSE]
        y <- y_full[valid]
        
        if (USE_FEATURE_SCREENING && length(important_windows_list[[trait]]) > 0) {
          X <- X[, important_windows_list[[trait]], drop = FALSE]
        }
        
        all_windows <- colnames(X)
        n_windows <- ncol(X)
        n_samples <- nrow(X) # <--- Fixed: Defined 'n' variable here
        
        coef_matrix <- matrix(NA, nrow = BOOTSTRAP_REPS, ncol = n_windows)
        colnames(coef_matrix) <- all_windows
        high_impact_counts <- setNames(integer(n_windows), all_windows)
        
        for (b in 1:BOOTSTRAP_REPS) {
          set.seed(b)
          # Resample hybrids with replacement
          idx <- sample(1:n_samples, size = n_samples, replace = TRUE)
          X_boot <- X[idx, , drop = FALSE]
          y_boot <- y[idx]
          
          keep <- apply(X_boot, 2, var) > 0
          if (sum(keep) == 0) next
          X_boot <- X_boot[, keep, drop = FALSE]
          
          # 5-fold CV for tuning lambda in each bootstrap sample
          fit <- tryCatch(cv.glmnet(X_boot, y_boot, alpha = 1, nfolds = 5), error = function(e) NULL)
          if (is.null(fit)) next
          
          lambda <- if (USE_LAMBDA_1SE) fit$lambda.1se else fit$lambda.min
          nz <- sum(coef(fit, s = lambda) != 0)
          if (nz < MIN_REGRESSOR) {
            idx_lambda <- which(fit$nzero >= MIN_REGRESSOR)
            if (length(idx_lambda) > 0) lambda <- fit$lambda[max(idx_lambda)]
          }
          
          coef_vec <- as.vector(coef(fit, s = lambda))[-1]
          names(coef_vec) <- colnames(X_boot)
          
          full_coef <- rep(0, n_windows)
          names(full_coef) <- all_windows
          common <- intersect(names(coef_vec), all_windows)
          full_coef[common] <- coef_vec[common]
          coef_matrix[b, ] <- full_coef
          
          abs_coef <- abs(full_coef)
          if (length(unique(abs_coef)) == 1) next
          
          # Top/bottom 10% criterion
          q_top <- quantile(abs_coef, probs = 1 - TOP_PCT, na.rm = TRUE)
          q_bottom <- quantile(abs_coef, probs = BOTTOM_PCT, na.rm = TRUE)
          high_impact <- which(abs_coef >= q_top | abs_coef <= q_bottom)
          
          for (idx_win in high_impact) {
            high_impact_counts[all_windows[idx_win]] <- high_impact_counts[all_windows[idx_win]] + 1
          }
        }
        
        freq <- high_impact_counts / BOOTSTRAP_REPS
        freq_dt <- data.table(feature = names(freq), freq = as.numeric(freq))
        setorder(freq_dt, -freq)
        
        coef_mean <- apply(coef_matrix, 2, mean, na.rm = TRUE)
        coef_sd <- apply(coef_matrix, 2, sd, na.rm = TRUE)
        coef_median <- apply(coef_matrix, 2, median, na.rm = TRUE)
        coef_025 <- apply(coef_matrix, 2, quantile, probs = 0.025, na.rm = TRUE)
        coef_975 <- apply(coef_matrix, 2, quantile, probs = 0.975, na.rm = TRUE)
        
        coef_stats <- data.table(feature = all_windows, mean_coef = coef_mean, sd_coef = coef_sd,
                                 median_coef = coef_median, ci_lower = coef_025, ci_upper = coef_975)
        
        final_dt <- merge(freq_dt, coef_stats, by = "feature", all = TRUE)
        setorder(final_dt, -freq)
        fwrite(final_dt, file.path(output_dir, paste0("bootstrap_stats_quantile_", window_name, "_", trait, "_", run_date, ".csv")))
        
        robust <- final_dt[freq >= ROBUST_THRESHOLD]
        if (nrow(robust) > 0) {
          fwrite(robust, file.path(output_dir, paste0("bootstrap_robust_quantile_", window_name, "_", trait, "_", run_date, ".csv")))
        }
      }
    }
  }
}