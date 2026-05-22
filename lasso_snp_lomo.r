# ============================================================
# SNP 随机抽样敏感性分析 - 混合模型 + LASSO（LOMO-CV）
# 功能：对随机抽取的0.5% SNP子集，先拟合母本随机效应（LMM），
#       再用弹性网预测残差，从而提升留一母本（LOMO）泛化能力
# 验证策略：LOMO-CV（留一母本）+ 外部测试集 + 训练集10折CV
# 循环次数：51次（种子150-200）
# 输出文件夹：SNP_sensitivity_LMM_LASSO_LOMO_YYYY-MM-DD
# 每个种子单独保存结果文件
# ============================================================

# 加载必要的包
library(vroom)
library(data.table)
library(glmnet)
library(lme4)          # 线性混合模型
library(caret)
library(magrittr)

# 获取当前日期
run_date <- Sys.Date()
output_dir <- paste0("SNP_sensitivity_LMM_LASSO_LOMO_", run_date)
if (!dir.exists(output_dir)) dir.create(output_dir)
cat("输出文件夹：", output_dir, "\n")

# ── 参数设置 ────────────────────────────────────────────────
SNP_FILE       <- "hybrid_SNP_500_raw.csv"   # 500行 × 630万列的超大CSV
PHENO_TRAIN    <- "phenotype_training_336.csv"
PHENO_TEST     <- "phenotype_test_92.csv"
SEEDS          <- 150:200                    # 51个随机种子
SNP_RATIO      <- 0.02                       # 随机抽0.5%
MIN_REGRESSOR  <- 50                         # LASSO最小变量数（提高稳定性）
ALPHA          <- 1                          # LASSO（alpha=1）
USE_LAMBDA_1SE <- TRUE                       # 使用更稳定的 lambda.1se
N_FOLDS        <- 10                         # 10折交叉验证（用于外部测试集评估）
TRAITS         <- c("FLW_HPH", "FSTW_HPH", "FSHW_HPH")
meta_cols      <- c("hybrid_id", "cms", "restorer")

# ── 读取表型数据（固定）──────────────────────────────────────
cat("══ 读取表型数据 ══\n")
pheno_train <- fread(PHENO_TRAIN, encoding = "UTF-8")
pheno_test  <- fread(PHENO_TEST, encoding = "UTF-8")

clean_id <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("[^A-Za-z0-9_]", "", x)
  return(x)
}
pheno_train[, hybrid_id := clean_id(hybrid_id)]
pheno_test[, hybrid_id := clean_id(hybrid_id)]
cat("训练集样本数：", nrow(pheno_train), "\n")
cat("测试集样本数：", nrow(pheno_test), "\n\n")

# ── 获取CSV文件的列名（只读表头）──────────────────────────────
cat("══ 读取CSV列名 ══\n")
col_names <- names(vroom(SNP_FILE, n_max = 0, show_col_types = FALSE))
snp_cols_all <- setdiff(col_names, meta_cols)
total_snp <- length(snp_cols_all)
cat("总SNP数：", total_snp, "\n")
cat("0.5%抽样后约：", round(total_snp * SNP_RATIO), "个SNP\n\n")

# ============================================================
# 弹性网直接预测（用于外部测试集，不使用混合模型）
# ============================================================
run_lasso_enet <- function(X_train, y_train, X_test, y_test,
                           n_folds = 10, min_regressor = 20, seed = 42,
                           alpha = 0.5, use_lambda_1se = TRUE) {
  
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
    
    m <- tryCatch(
      cv.glmnet(x = X_fold_train, y = y_fold_train, nfolds = 5, alpha = alpha),
      error = function(e) NULL
    )
    if (is.null(m)) {
      all_pred_cv[test_idx] <- NA
      all_true_cv[test_idx] <- y_tr[test_idx]
      n_features_vec[i] <- 0
      next
    }
    
    lambda <- if (use_lambda_1se) m$lambda.1se else m$lambda.min
    nz <- sum(coef(m, s = lambda) != 0) - 1
    if (nz < min_regressor) {
      idx <- which(m$nzero >= min_regressor)
      if (length(idx) == 0) {
        all_pred_cv[test_idx] <- NA
        all_true_cv[test_idx] <- y_tr[test_idx]
        n_features_vec[i] <- 0
        next
      }
      lambda <- m$lambda[max(idx)]
    }
    
    pred_test <- as.vector(predict(m, newx = X_fold_test, s = lambda))
    all_pred_cv[test_idx] <- pred_test
    all_true_cv[test_idx] <- y_tr[test_idx]
    n_features_vec[i] <- sum(coef(m, s = lambda) != 0) - 1
  }
  
  valid_pred_cv <- !is.na(all_pred_cv)
  r_cv <- if (sum(valid_pred_cv) >= 10) {
    cor(all_true_cv[valid_pred_cv], all_pred_cv[valid_pred_cv], use = "complete.obs")
  } else NA_real_
  
  # 全量训练 + 外部测试
  r_test <- NA_real_
  train_r2 <- NA_real_
  test_rmse <- NA_real_
  pred_test_full <- NULL
  
  if (length(y_te) > 5) {
    m_full <- tryCatch(
      cv.glmnet(x = X_tr, y = y_tr, nfolds = 5, alpha = alpha),
      error = function(e) NULL
    )
    if (!is.null(m_full)) {
      lambda <- if (use_lambda_1se) m_full$lambda.1se else m_full$lambda.min
      nz <- sum(coef(m_full, s = lambda) != 0) - 1
      if (nz < min_regressor) {
        idx <- which(m_full$nzero >= min_regressor)
        if (length(idx) > 0) lambda <- m_full$lambda[max(idx)]
      }
      
      train_pred <- as.vector(predict(m_full, newx = X_tr, s = lambda))
      valid_train_full <- !is.na(y_tr) & !is.na(train_pred)
      if (sum(valid_train_full) >= 2) {
        train_r2 <- cor(train_pred[valid_train_full], y_tr[valid_train_full])^2
      }
      
      pred_test_full <- as.vector(predict(m_full, newx = X_te, s = lambda))
      valid_test_full <- !is.na(y_te) & !is.na(pred_test_full)
      if (sum(valid_test_full) >= 5) {
        r_test <- cor(pred_test_full[valid_test_full], y_te[valid_test_full])
        test_rmse <- sqrt(mean((pred_test_full[valid_test_full] - y_te[valid_test_full])^2))
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

# ============================================================
# 混合模型 + LASSO 预测函数（用于 LOMO-CV）
# 原理：先用 LMM 分离母本随机效应，再用弹性网预测残差
# 对新母本，随机效应为 0，因此预测 = 总体均值 + 残差预测
# ============================================================
predict_lmm_lasso <- function(X_tr, y_tr, cms_tr, X_val, cms_val,
                              min_regressor = 20, alpha = 0.5, use_lambda_1se = TRUE) {
  
  # 1. 检查输入
  if (length(unique(cms_tr)) < 2) {
    # 如果训练集中母本水平少于2，无法拟合随机效应，退化为简单均值预测
    mu <- mean(y_tr, na.rm = TRUE)
    return(rep(mu, nrow(X_val)))
  }
  
  # 2. 拟合线性混合模型（母本随机截距）
  df_tr <- data.frame(y = y_tr, cms = cms_tr)
  # 使用 lmer，固定效应只有截距
  lmm_fit <- tryCatch(
    lmer(y ~ 1 + (1 | cms), data = df_tr, REML = FALSE),
    error = function(e) NULL
  )
  if (is.null(lmm_fit)) {
    # 若拟合失败，退化为均值预测
    mu <- mean(y_tr, na.rm = TRUE)
    return(rep(mu, nrow(X_val)))
  }
  
  # 3. 计算训练集残差（原始 y - 拟合值，拟合值包含总体均值和母本 BLUP）
  y_hat <- fitted(lmm_fit)
  resid_tr <- y_tr - y_hat
  
  # 4. 用弹性网建模残差 ~ SNP
  #    去除方差为零的列
  zero_var <- apply(X_tr, 2, var) == 0
  if (all(zero_var)) {
    mu <- fixef(lmm_fit)[["(Intercept)"]]
    return(rep(mu, nrow(X_val)))
  }
  X_tr_filt <- X_tr[, !zero_var, drop = FALSE]
  
  m_lasso <- tryCatch(
    cv.glmnet(x = X_tr_filt, y = resid_tr, alpha = alpha, nfolds = 5),
    error = function(e) NULL
  )
  if (is.null(m_lasso)) {
    mu <- fixef(lmm_fit)[["(Intercept)"]]
    return(rep(mu, nrow(X_val)))
  }
  
  lambda <- if (use_lambda_1se) m_lasso$lambda.1se else m_lasso$lambda.min
  nz <- sum(coef(m_lasso, s = lambda) != 0)
  if (nz < min_regressor) {
    idx <- which(m_lasso$nzero >= min_regressor)
    if (length(idx) > 0) lambda <- m_lasso$lambda[max(idx)]
  }
  
  # 5. 预测验证集残差
  X_val_filt <- X_val[, !zero_var, drop = FALSE]
  if (ncol(X_val_filt) == 0) {
    pred_resid <- rep(0, nrow(X_val))
  } else {
    pred_resid <- as.vector(predict(m_lasso, newx = X_val_filt, s = lambda))
  }
  
  # 6. 预测验证集表型 = 总体均值 + 残差预测
  mu <- fixef(lmm_fit)[["(Intercept)"]]
  y_pred <- mu + pred_resid
  return(as.vector(y_pred))
}

# ============================================================
# 单次种子处理函数（仅 LASSO，包含 LOMO-CV 和外部测试）
# 每个种子完成后返回结果 data.table
# ============================================================
run_single_seed <- function(seed, snp_cols_all, pheno_train, pheno_test,
                            TRAITS, SNP_RATIO, MIN_REGRESSOR,
                            ALPHA, USE_LAMBDA_1SE) {
  
  set.seed(seed)
  n_select <- round(length(snp_cols_all) * SNP_RATIO)
  selected_snps <- sample(snp_cols_all, size = n_select)
  cols_to_select <- c(meta_cols, selected_snps)
  
  # 读取 SNP 子集
  snp_subset <- vroom(SNP_FILE, 
                      col_select = all_of(cols_to_select),
                      show_col_types = FALSE,
                      progress = FALSE) %>%
    as.data.table()
  snp_subset[, hybrid_id := clean_id(hybrid_id)]
  
  # 对齐训练集
  train_ids <- intersect(snp_subset$hybrid_id, pheno_train$hybrid_id)
  snp_train <- snp_subset[hybrid_id %in% train_ids]
  pheno_train_sub <- pheno_train[hybrid_id %in% train_ids]
  setkey(snp_train, hybrid_id)
  setkey(pheno_train_sub, hybrid_id)
  snp_train <- snp_train[pheno_train_sub$hybrid_id]
  X_train_full <- as.matrix(snp_train[, ..selected_snps])
  X_train_full[is.na(X_train_full)] <- 0
  cms_train <- snp_train$cms          # 用于 LOMO
  
  # 对齐测试集（外部测试）
  test_ids <- intersect(snp_subset$hybrid_id, pheno_test$hybrid_id)
  if (length(test_ids) > 0) {
    snp_test <- snp_subset[hybrid_id %in% test_ids]
    pheno_test_sub <- pheno_test[hybrid_id %in% test_ids]
    setkey(snp_test, hybrid_id)
    setkey(pheno_test_sub, hybrid_id)
    snp_test <- snp_test[pheno_test_sub$hybrid_id]
    X_test_full <- as.matrix(snp_test[, ..selected_snps])
    X_test_full[is.na(X_test_full)] <- 0
  } else {
    X_test_full <- NULL
    pheno_test_sub <- NULL
  }
  
  results <- list()
  
  for (trait in TRAITS) {
    y_train <- pheno_train_sub[[trait]]
    y_test  <- if (!is.null(pheno_test_sub)) pheno_test_sub[[trait]] else NULL
    
    # 剔除训练集中 y 为 NA 的样本
    valid_train <- !is.na(y_train)
    X_train <- X_train_full[valid_train, , drop = FALSE]
    y_train <- y_train[valid_train]
    cms_train_trait <- cms_train[valid_train]
    
    # 外部测试集表型有效样本
    if (!is.null(X_test_full) && !is.null(y_test)) {
      valid_test <- !is.na(y_test)
      X_test <- X_test_full[valid_test, , drop = FALSE]
      y_test <- y_test[valid_test]
    } else {
      X_test <- NULL
      y_test <- NULL
    }
    
    # ------------------- 外部测试集预测（使用常规弹性网，无混合模型） -------------------
    if (!is.null(X_test) && length(y_test) > 0) {
      lasso_ext <- tryCatch(
        run_lasso_enet(X_train, y_train, X_test, y_test,
                       n_folds = N_FOLDS, min_regressor = MIN_REGRESSOR, seed = seed,
                       alpha = ALPHA, use_lambda_1se = USE_LAMBDA_1SE),
        error = function(e) list(r_cv = NA_real_, r_test = NA_real_, n_features = 0,
                                 train_r2 = NA_real_, test_rmse = NA_real_,
                                 pred_test = NULL, y_test = NULL)
      )
      r_cv_ext <- lasso_ext$r_cv
      r_test_ext <- lasso_ext$r_test
      train_r2_ext <- lasso_ext$train_r2
      test_rmse_ext <- lasso_ext$test_rmse
      n_features_ext <- lasso_ext$n_features
    } else {
      r_cv_ext <- NA_real_
      r_test_ext <- NA_real_
      train_r2_ext <- NA_real_
      test_rmse_ext <- NA_real_
      n_features_ext <- 0
    }
    
    # ------------------- LOMO-CV（混合模型 + LASSO，留一母本） -------------------
    unique_cms <- unique(cms_train_trait)
    pred_lasso_all <- rep(NA, length(y_train))
    
    for (cc in unique_cms) {
      val_idx <- which(cms_train_trait == cc)
      train_idx <- setdiff(seq_along(y_train), val_idx)
      if (length(train_idx) < 2 || length(val_idx) == 0) next
      
      X_tr_snp <- X_train[train_idx, , drop = FALSE]
      y_tr <- y_train[train_idx]
      cms_tr <- cms_train_trait[train_idx]
      
      X_val_snp <- X_train[val_idx, , drop = FALSE]
      cms_val <- cms_train_trait[val_idx]   # 全为同一个母本
      
      pred <- tryCatch(
        predict_lmm_lasso(X_tr_snp, y_tr, cms_tr, X_val_snp, cms_val,
                          min_regressor = MIN_REGRESSOR, alpha = ALPHA,
                          use_lambda_1se = USE_LAMBDA_1SE),
        error = function(e) rep(NA, length(val_idx))
      )
      pred_lasso_all[val_idx] <- pred
    }
    
    valid_lomo <- !is.na(pred_lasso_all)
    r_lomo <- if (sum(valid_lomo) >= 10) {
      cor(y_train[valid_lomo], pred_lasso_all[valid_lomo], use = "complete.obs")
    } else NA_real_
    
    # 收集结果
    results[[length(results) + 1]] <- data.table(
      seed = seed,
      trait = trait,
      LASSO_lomo = r_lomo,
      LASSO_cv = r_cv_ext,
      LASSO_test = r_test_ext,
      LASSO_trainR2 = train_r2_ext,
      LASSO_testRMSE = test_rmse_ext,
      n_features = n_features_ext
    )
  }
  return(rbindlist(results))
}

# ============================================================
# 主循环（串行执行，每个种子单独保存结果）
# ============================================================
cat("══ SNP 随机抽样敏感性分析（混合模型 + LASSO，LOMO-CV） ══\n")
cat("方法：每次随机抽取0.5% SNP列，先用 LMM 分离母本随机效应，再用 LASSO 预测残差\n")
cat("验证：LOMO-CV（留一母本）+ 外部测试集(92样本) + 训练集10折CV\n")
cat("参数：LASSO alpha=", ALPHA, ", min_regressor=", MIN_REGRESSOR, ", lambda.1se=", USE_LAMBDA_1SE, "\n")
cat("随机种子：150-200（共51次）\n")
cat("运行日期：", run_date, "\n")
cat("输出文件夹：", output_dir, "\n")
cat("每个种子结果将单独保存为 seed_XXX.csv\n\n")

start_time <- Sys.time()
SEEDS_VEC <- 150:200

all_results_list <- list()

for (i in seq_along(SEEDS_VEC)) {
  seed <- SEEDS_VEC[i]
  cat(sprintf("[%s] 处理种子 %d / %d (seed=%d)\n", 
              format(Sys.time(), "%H:%M:%S"), i, length(SEEDS_VEC), seed))
  
  res <- run_single_seed(seed, snp_cols_all, pheno_train, pheno_test,
                         TRAITS, SNP_RATIO, MIN_REGRESSOR,
                         ALPHA, USE_LAMBDA_1SE)
  
  # 保存该种子的结果到单独文件
  seed_file <- file.path(output_dir, paste0("seed_", seed, ".csv"))
  fwrite(res, seed_file)
  cat("  已保存: ", seed_file, "\n")
  
  all_results_list[[i]] <- res
  
  if (i %% 10 == 0) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    avg_per_seed <- elapsed / i
    remaining <- avg_per_seed * (length(SEEDS_VEC) - i)
    cat(sprintf("  已运行 %.1f 分钟，预计剩余 %.1f 分钟\n", elapsed, remaining))
  }
  
  gc()
}

end_time <- Sys.time()
cat("\n计算完成！总耗时：", round(difftime(end_time, start_time, units = "mins"), 2), "分钟\n\n")

# 合并所有种子结果并保存总汇总
all_results <- rbindlist(all_results_list, fill = TRUE)
summary_dt <- all_results[, .(
  LASSO_lomo_mean = round(mean(LASSO_lomo, na.rm = TRUE), 4),
  LASSO_lomo_sd = round(sd(LASSO_lomo, na.rm = TRUE), 4),
  LASSO_cv_mean = round(mean(LASSO_cv, na.rm = TRUE), 4),
  LASSO_cv_sd = round(sd(LASSO_cv, na.rm = TRUE), 4),
  LASSO_test_mean = round(mean(LASSO_test, na.rm = TRUE), 4),
  LASSO_test_sd = round(sd(LASSO_test, na.rm = TRUE), 4),
  LASSO_trainR2_mean = round(mean(LASSO_trainR2, na.rm = TRUE), 4),
  LASSO_trainR2_sd = round(sd(LASSO_trainR2, na.rm = TRUE), 4),
  LASSO_testRMSE_mean = round(mean(LASSO_testRMSE, na.rm = TRUE), 4),
  LASSO_testRMSE_sd = round(sd(LASSO_testRMSE, na.rm = TRUE), 4),
  n_features_mean = round(mean(n_features, na.rm = TRUE), 1),
  n_seeds = .N
), by = trait]

# 保存总汇总文件
fwrite(all_results, file.path(output_dir, paste0("all_results_", run_date, ".csv")))
fwrite(summary_dt, file.path(output_dir, paste0("summary_", run_date, ".csv")))

# 打印汇总结果
cat(strrep("═", 80), "\n")
cat("  汇总结果（", length(SEEDS_VEC), "个种子）\n", sep="")
cat(strrep("═", 80), "\n\n")

cat("【LOMO-CV结果 (r_lomo)】\n")
print(summary_dt[, .(trait, `LASSO_lomo` = sprintf("%.3f±%.3f", LASSO_lomo_mean, LASSO_lomo_sd))])

cat("\n【训练集10折交叉验证结果 (r_cv)】\n")
print(summary_dt[, .(trait, `LASSO_cv` = sprintf("%.3f±%.3f", LASSO_cv_mean, LASSO_cv_sd))])

cat("\n【外部测试集结果 (r_test, n≈92)】\n")
print(summary_dt[, .(trait, `LASSO_test` = sprintf("%.3f±%.3f", LASSO_test_mean, LASSO_test_sd))])

cat("\n【训练集R2及测试RMSE】\n")
print(summary_dt[, .(trait,
                     `Train R²` = sprintf("%.3f±%.3f", LASSO_trainR2_mean, LASSO_trainR2_sd),
                     `Test RMSE` = sprintf("%.3f±%.3f", LASSO_testRMSE_mean, LASSO_testRMSE_sd),
                     `#Features` = n_features_mean)])

cat("\n✓ 所有结果已保存至：", output_dir, "\n")
cat("  包括每个种子的单独文件（seed_XXX.csv）及总汇总文件。\n")
cat("\n══ 分析完成 ══\n")