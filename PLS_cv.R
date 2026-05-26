
library(vroom)
library(data.table)
library(pls)
library(caret)

SNP_FILE <- "hybrid_SNP_500_raw.csv"
PHENO_TRAIN <- "phenotype_training_336.csv"
PHENO_TEST <- "phenotype_test_92.csv"
SEEDS <- 101:150
SNP_RATIO <- 0.05
N_FOLDS <- 3
TRAITS <- c("FLW_HPH", "FSTW_HPH", "FSHW_HPH")
meta_cols <- c("hybrid_id", "cms", "restorer")
NCOMP <- 3

pheno_train <- fread(PHENO_TRAIN)
pheno_test <- fread(PHENO_TEST)
clean_id <- function(x) gsub("[^A-Za-z0-9_]", "", trimws(as.character(x)))
pheno_train[, hybrid_id := clean_id(hybrid_id)]
pheno_test[, hybrid_id := clean_id(hybrid_id)]

col_names <- names(vroom(SNP_FILE, n_max = 0, show_col_types = FALSE))
snp_cols_all <- setdiff(col_names, meta_cols)

run_pls <- function(X_tr, y_tr, X_te, y_te, seed) {
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
    X_tr_fold <- scale(X_tr[train_idx, ], center = TRUE, scale = TRUE)
    X_tr_fold[is.nan(X_tr_fold)] <- 0
    center <- attr(X_tr_fold, "scaled:center")
    scale_vec <- attr(X_tr_fold, "scaled:scale")
    scale_vec[scale_vec == 0] <- 1
    df <- data.frame(y = y_tr[train_idx], X_tr_fold)
    m <- plsr(y ~ ., data = df, ncomp = NCOMP, validation = "none")
    X_te_fold <- scale(X_tr[test_idx, ], center = center, scale = scale_vec)
    X_te_fold[is.nan(X_te_fold)] <- 0
    pred_cv[test_idx] <- predict(m, newdata = data.frame(X_te_fold), ncomp = NCOMP)
  }
  r_cv <- cor(y_tr, pred_cv, use = "complete.obs")
  
  X_tr_all <- scale(X_tr, center = TRUE, scale = TRUE)
  X_tr_all[is.nan(X_tr_all)] <- 0
  center_all <- attr(X_tr_all, "scaled:center")
  scale_all <- attr(X_tr_all, "scaled:scale")
  scale_all[scale_all == 0] <- 1
  m_full <- plsr(y_tr ~ ., data = data.frame(y = y_tr, X_tr_all), ncomp = NCOMP)
  X_te_all <- scale(X_te, center = center_all, scale = scale_all)
  X_te_all[is.nan(X_te_all)] <- 0
  pred_te <- predict(m_full, newdata = data.frame(X_te_all), ncomp = NCOMP)
  r_test <- cor(y_te, pred_te, use = "complete.obs")
  list(r_cv = r_cv, r_test = r_test)
}

predict_pls <- function(X_tr, y_tr, X_te) {
  valid <- !is.na(y_tr)
  X_tr <- X_tr[valid, , drop = FALSE]
  y_tr <- y_tr[valid]
  if (nrow(X_tr) < 2) return(rep(NA, nrow(X_te)))
  X_tr <- scale(X_tr, center = TRUE, scale = TRUE)
  X_tr[is.nan(X_tr)] <- 0
  center <- attr(X_tr, "scaled:center")
  scale_vec <- attr(X_tr, "scaled:scale")
  scale_vec[scale_vec == 0] <- 1
  m <- plsr(y_tr ~ ., data = data.frame(y = y_tr, X_tr), ncomp = NCOMP)
  X_te <- scale(X_te, center = center, scale = scale_vec)
  X_te[is.nan(X_te)] <- 0
  as.vector(predict(m, newdata = data.frame(X_te), ncomp = NCOMP))
}

all_res <- list()
for (seed in SEEDS) {
  set.seed(seed)
  selected <- sample(snp_cols_all, size = round(length(snp_cols_all) * SNP_RATIO))
  cols <- c(meta_cols, selected)
  snp <- vroom(SNP_FILE, col_select = all_of(cols), show_col_types = FALSE) %>% as.data.table()
  snp[, hybrid_id := clean_id(hybrid_id)]
  
  tr_ids <- intersect(snp$hybrid_id, pheno_train$hybrid_id)
  snp_tr <- snp[hybrid_id %in% tr_ids]
  pheno_tr <- pheno_train[hybrid_id %in% tr_ids]
  setkey(snp_tr, hybrid_id); setkey(pheno_tr, hybrid_id)
  snp_tr <- snp_tr[pheno_tr$hybrid_id]
  X_tr_full <- as.matrix(snp_tr[, ..selected]); X_tr_full[is.na(X_tr_full)] <- 0
  cms_tr <- snp_tr$cms; restorer_tr <- snp_tr$restorer
  
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
    
    pred_lomo <- rep(NA, length(y_tr))
    for (g in unique(cms)) {
      val <- which(cms == g)
      if (length(val) == 0 || length(setdiff(seq_along(y_tr), val)) < 2) next
      pred_lomo[val] <- predict_pls(X_tr[-val, , drop = FALSE], y_tr[-val], X_tr[val, , drop = FALSE])
    }
    r_lomo <- cor(y_tr, pred_lomo, use = "complete.obs")
    
    pred_loro <- rep(NA, length(y_tr))
    for (g in unique(restorer)) {
      val <- which(restorer == g)
      if (length(val) == 0 || length(setdiff(seq_along(y_tr), val)) < 2) next
      pred_loro[val] <- predict_pls(X_tr[-val, , drop = FALSE], y_tr[-val], X_tr[val, , drop = FALSE])
    }
    r_loro <- cor(y_tr, pred_loro, use = "complete.obs")
    
    if (!is.null(X_te_full) && !is.null(pheno_te[[trait]])) {
      y_te <- pheno_te[[trait]]
      valid_te <- !is.na(y_te)
      X_te <- X_te_full[valid_te, , drop = FALSE]
      y_te <- y_te[valid_te]
      pls_res <- run_pls(X_tr, y_tr, X_te, y_te, seed)
      r_cv <- pls_res$r_cv
      r_test <- pls_res$r_test
    } else {
      r_cv <- r_test <- NA
    }
    
    all_res[[length(all_res) + 1]] <- data.table(seed, trait, r_lomo, r_loro, r_cv, r_test)
  }
}

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