

library(vroom)
library(data.table)
library(glmnet)
library(lme4)
library(caret)

SNP_FILE <- "hybrid_SNP_500_raw.csv"
PHENO_TRAIN <- "phenotype_training_336.csv"
PHENO_TEST <- "phenotype_test_92.csv"
SEEDS <- 101:150
SNP_RATIO <- 0.05
MIN_REGRESSOR <- 50
ALPHA <- 1
USE_LAMBDA_1SE <- TRUE
N_FOLDS <- 10
TRAITS <- c("FLW_HPH", "FSTW_HPH", "FSHW_HPH")
meta_cols <- c("hybrid_id", "cms", "restorer")

pheno_train <- fread(PHENO_TRAIN)
pheno_test <- fread(PHENO_TEST)
clean_id <- function(x) gsub("[^A-Za-z0-9_]", "", trimws(as.character(x)))
pheno_train[, hybrid_id := clean_id(hybrid_id)]
pheno_test[, hybrid_id := clean_id(hybrid_id)]

col_names <- names(vroom(SNP_FILE, n_max = 0, show_col_types = FALSE))
snp_cols_all <- setdiff(col_names, meta_cols)

# Elastic-net for external test
run_enet <- function(X_tr, y_tr, X_te, y_te, seed) {
  valid_tr <- !is.na(y_tr)
  X_tr <- X_tr[valid_tr, , drop = FALSE]
  y_tr <- y_tr[valid_tr]
  valid_te <- !is.na(y_te)
  X_te <- X_te[valid_te, , drop = FALSE]
  y_te <- y_te[valid_te]
  
  set.seed(seed)
  folds <- createFolds(y_tr, k = N_FOLDS, returnTrain = FALSE)
  pred_cv <- numeric(length(y_tr))
  for (i in seq_along(folds)) {
    test_idx <- folds[[i]]
    train_idx <- setdiff(seq_along(y_tr), test_idx)
    m <- cv.glmnet(X_tr[train_idx, ], y_tr[train_idx], alpha = ALPHA, nfolds = 5)
    lam <- if (USE_LAMBDA_1SE) m$lambda.1se else m$lambda.min
    if (sum(coef(m, s = lam) != 0) - 1 < MIN_REGRESSOR) {
      idx <- which(m$nzero >= MIN_REGRESSOR)
      if (length(idx)) lam <- m$lambda[max(idx)]
    }
    pred_cv[test_idx] <- predict(m, newx = X_tr[test_idx, ], s = lam)
  }
  r_cv <- cor(y_tr, pred_cv, use = "complete.obs")
  
  m_full <- cv.glmnet(X_tr, y_tr, alpha = ALPHA, nfolds = 5)
  lam <- if (USE_LAMBDA_1SE) m_full$lambda.1se else m_full$lambda.min
  if (sum(coef(m_full, s = lam) != 0) - 1 < MIN_REGRESSOR) {
    idx <- which(m_full$nzero >= MIN_REGRESSOR)
    if (length(idx)) lam <- m_full$lambda[max(idx)]
  }
  pred_te <- predict(m_full, newx = X_te, s = lam)
  r_test <- cor(y_te, pred_te, use = "complete.obs")
  list(r_cv = r_cv, r_test = r_test)
}

# LMM + LASSO for a given grouping factor
predict_lasso_lmm <- function(X_tr, y_tr, group_tr, X_val, group_val) {
  if (length(unique(group_tr)) < 2) return(rep(mean(y_tr, na.rm = TRUE), nrow(X_val)))
  lmm <- lmer(y_tr ~ 1 + (1 | group_tr), REML = FALSE)
  resid <- residuals(lmm)
  zero_var <- apply(X_tr, 2, var) == 0
  if (all(zero_var)) return(rep(fixef(lmm)[[1]], nrow(X_val)))
  X_tr <- X_tr[, !zero_var, drop = FALSE]
  X_val <- X_val[, !zero_var, drop = FALSE]
  m <- cv.glmnet(X_tr, resid, alpha = ALPHA, nfolds = 5)
  lam <- if (USE_LAMBDA_1SE) m$lambda.1se else m$lambda.min
  if (sum(coef(m, s = lam) != 0) - 1 < MIN_REGRESSOR) {
    idx <- which(m$nzero >= MIN_REGRESSOR)
    if (length(idx)) lam <- m$lambda[max(idx)]
  }
  pred_resid <- predict(m, newx = X_val, s = lam)
  fixef(lmm)[[1]] + as.vector(pred_resid)
}

# Main loop over seeds
all_res <- list()
for (seed in SEEDS) {
  set.seed(seed)
  selected <- sample(snp_cols_all, size = round(length(snp_cols_all) * SNP_RATIO))
  cols <- c(meta_cols, selected)
  snp <- vroom(SNP_FILE, col_select = all_of(cols), show_col_types = FALSE) %>% as.data.table()
  snp[, hybrid_id := clean_id(hybrid_id)]
  
  # training data
  tr_ids <- intersect(snp$hybrid_id, pheno_train$hybrid_id)
  snp_tr <- snp[hybrid_id %in% tr_ids]
  pheno_tr <- pheno_train[hybrid_id %in% tr_ids]
  setkey(snp_tr, hybrid_id); setkey(pheno_tr, hybrid_id)
  snp_tr <- snp_tr[pheno_tr$hybrid_id]
  X_tr_full <- as.matrix(snp_tr[, ..selected]); X_tr_full[is.na(X_tr_full)] <- 0
  cms_tr <- snp_tr$cms; restorer_tr <- snp_tr$restorer
  
  # test data
  te_ids <- intersect(snp$hybrid_id, pheno_test$hybrid_id)
  if (length(te_ids)) {
    snp_te <- snp[hybrid_id %in% te_ids]
    pheno_te <- pheno_test[hybrid_id %in% te_ids]
    setkey(snp_te, hybrid_id); setkey(pheno_te, hybrid_id)
    snp_te <- snp_te[pheno_te$hybrid_id]
    X_te_full <- as.matrix(snp_te[, ..selected]); X_te_full[is.na(X_te_full)] <- 0
  } else X_te_full <- NULL
  
  for (trait in TRAITS) {
    y_tr <- pheno_tr[[trait]]
    valid <- !is.na(y_tr)
    X_tr <- X_tr_full[valid, , drop = FALSE]
    y_tr <- y_tr[valid]
    cms <- cms_tr[valid]; restorer <- restorer_tr[valid]
    
    # LOMO
    pred_lomo <- rep(NA, length(y_tr))
    for (g in unique(cms)) {
      val <- which(cms == g)
      if (length(val) == 0 || length(setdiff(seq_along(y_tr), val)) < 2) next
      pred_lomo[val] <- predict_lasso_lmm(X_tr[-val, , drop = FALSE], y_tr[-val], cms[-val],
                                          X_tr[val, , drop = FALSE], cms[val])
    }
    r_lomo <- cor(y_tr, pred_lomo, use = "complete.obs")
    
    # LORO
    pred_loro <- rep(NA, length(y_tr))
    for (g in unique(restorer)) {
      val <- which(restorer == g)
      if (length(val) == 0 || length(setdiff(seq_along(y_tr), val)) < 2) next
      pred_loro[val] <- predict_lasso_lmm(X_tr[-val, , drop = FALSE], y_tr[-val], restorer[-val],
                                          X_tr[val, , drop = FALSE], restorer[val])
    }
    r_loro <- cor(y_tr, pred_loro, use = "complete.obs")
    
    # external test
    if (!is.null(X_te_full) && !is.null(pheno_te[[trait]])) {
      y_te <- pheno_te[[trait]]
      valid_te <- !is.na(y_te)
      X_te <- X_te_full[valid_te, , drop = FALSE]
      y_te <- y_te[valid_te]
      enet <- run_enet(X_tr, y_tr, X_te, y_te, seed)
      r_cv <- enet$r_cv
      r_test <- enet$r_test
    } else {
      r_cv <- r_test <- NA
    }
    
    all_res[[length(all_res) + 1]] <- data.table(seed, trait, r_lomo, r_loro, r_cv, r_test)
  }
}

# Aggregate results
res_all <- rbindlist(all_res)
summary <- res_all[, .(
  LOMO_mean = mean(r_lomo, na.rm = TRUE),
  LOMO_sd = sd(r_lomo, na.rm = TRUE),
  LORO_mean = mean(r_loro, na.rm = TRUE),
  LORO_sd = sd(r_loro, na.rm = TRUE),
  CV_mean = mean(r_cv, na.rm = TRUE),
  CV_sd = sd(r_cv, na.rm = TRUE),
  Test_mean = mean(r_test, na.rm = TRUE),
  Test_sd = sd(r_test, na.rm = TRUE)
), by = trait]

print(summary)