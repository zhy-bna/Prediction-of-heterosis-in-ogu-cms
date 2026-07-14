# ============================================================
# Additive model prediction for 49,550 combinations
# Logistic transformation with HPH constrained to (-100%, +100%)
# ============================================================

library(data.table)
library(glmnet)

# ---------- User settings ----------
N_REPEATS <- 100
RANDOM_SEEDS <- 101:200
OUTPUT_TOP_N <- 500

train_add_file <- "/biodata/zhy/pgsi_ad/wang_500kb_simi.csv"
pred_add_file <- "/biodata/zhy/pgsi_ad/991_Restore.50_similarity_index.csv"
pheno_train_file <- "phenotype_training_336.csv"
pheno_test_file <- "phenotype_test_92.csv"

output_dir <- "External_Pred_500Kb_100rep_additive_logit_constrained"
if (!dir.exists(output_dir)) dir.create(output_dir)

# ---------- Logit transformation parameters ----------
L <- -1
U <- 1

logit_transform <- function(y, L, U) {
  p <- (y - L) / (U - L)
  p <- pmax(pmin(p, 1 - 1e-8), 1e-8)
  log(p / (1 - p))
}

inv_logit_transform <- function(eta, L, U) {
  p <- 1 / (1 + exp(-eta))
  L + (U - L) * p
}

# ---------- 1. Load and combine phenotypes ----------
cat("══ Load and combine phenotypes (336 + 92 = 428) ══\n")
pheno_336 <- fread(pheno_train_file)
pheno_336[, hybrid_id := trimws(hybrid_id)]
pheno_92 <- fread(pheno_test_file)
pheno_92[, hybrid_id := trimws(hybrid_id)]
pheno_all <- rbind(pheno_336, pheno_92, fill = TRUE)
setkey(pheno_all, hybrid_id)

# ---------- 2. Load training additive PGSI ----------
cat("Load training additive PGSI...\n")
train_add <- fread(train_add_file)
train_add[, hybrid_id := trimws(hybrid_id)]
setkey(train_add, hybrid_id)

feat_cols <- setdiff(names(train_add), c("hybrid_id", "cms", "restorer"))
cat("Number of features:", length(feat_cols), "\n")

# ---------- 3. Match training phenotypes ----------
common <- intersect(train_add$hybrid_id, pheno_all$hybrid_id)
cat("Matched hybrids:", length(common), "\n")
train_add <- train_add[hybrid_id %in% common]
pheno_train <- pheno_all[hybrid_id %in% common]
setkey(train_add, hybrid_id)
setkey(pheno_train, hybrid_id)

y_train_HPH <- pheno_train$FSHW_HPH
valid <- !is.na(y_train_HPH)
train_add <- train_add[hybrid_id %in% pheno_train$hybrid_id[valid]]
X_train <- as.matrix(train_add[, ..feat_cols])
X_train[is.na(X_train)] <- 0
y_train_HPH <- y_train_HPH[valid]
cat("Training effective samples:", length(y_train_HPH), "\n")
cat("Raw HPH range in training:", range(y_train_HPH), "\n")

y_train_logit <- logit_transform(y_train_HPH, L, U)
cat("Logit range in training:", range(y_train_logit), "\n")

# ---------- 4. Load prediction additive PGSI ----------
cat("Load prediction additive PGSI...\n")
pred_add <- fread(pred_add_file)
pred_add[, hybrid_id := trimws(hybrid_id)]
setkey(pred_add, hybrid_id)

common_feat <- intersect(feat_cols, names(pred_add))
cat("Effective features in prediction set:", length(common_feat), "\n")
if (length(common_feat) == 0) stop("No valid features in prediction set!")

X_pred <- as.matrix(pred_add[, ..common_feat])
X_pred[is.na(X_pred)] <- 0
hybrid_ids_pred <- pred_add$hybrid_id
cms_pred <- pred_add$cms
restorer_pred <- pred_add$restorer
cat("Number of prediction combinations:", nrow(X_pred), "\n")

X_train <- X_train[, match(common_feat, feat_cols), drop = FALSE]
cat("Final training features:", ncol(X_train), "\n")

# ---------- 5. Additive model prediction (100 repeats) ----------
cat("\n══════════════════════════════════════════════════════════════\n")
cat("Start additive model prediction (100 repeats) on logit(HPH)\n")
cat("══════════════════════════════════════════════════════════════\n")

pred_matrix_logit <- matrix(NA, nrow = nrow(X_pred), ncol = N_REPEATS)
rownames(pred_matrix_logit) <- hybrid_ids_pred

for (rep in 1:N_REPEATS) {
  set.seed(RANDOM_SEEDS[rep])
  if (rep %% 10 == 0) cat("Repeat", rep, "/", N_REPEATS, "\n")
  
  fit_cv <- tryCatch(
    cv.glmnet(X_train, y_train_logit, alpha = 1, nfolds = 10),
    error = function(e) NULL
  )
  if (is.null(fit_cv)) next
  
  lambda <- fit_cv$lambda.min
  if (sum(coef(fit_cv, s = lambda) != 0) < 1) {
    idx <- which(fit_cv$nzero >= 1)
    if (length(idx) > 0) lambda <- fit_cv$lambda[max(idx)]
  }
  
  coefs <- coef(fit_cv, s = lambda)
  selected_vars <- rownames(coefs)[which(coefs != 0)][-1]
  
  if (length(selected_vars) == 0) {
    pred_matrix_logit[, rep] <- mean(y_train_logit)
    next
  }
  
  df_train <- data.frame(y = y_train_logit, X_train[, selected_vars, drop = FALSE])
  lm_fit <- tryCatch(lm(y ~ ., data = df_train), error = function(e) NULL)
  if (is.null(lm_fit)) next
  
  df_pred <- as.data.frame(X_pred[, selected_vars, drop = FALSE])
  pred_matrix_logit[, rep] <- predict(lm_fit, newdata = df_pred)
}

pred_mean_logit <- rowMeans(pred_matrix_logit, na.rm = TRUE)
pred_sd_logit <- apply(pred_matrix_logit, 1, sd, na.rm = TRUE)

pred_mean_HPH <- inv_logit_transform(pred_mean_logit, L, U)

results <- data.table(
  hybrid_id = hybrid_ids_pred,
  cms = cms_pred,
  restorer = restorer_pred,
  pred_mean_HPH = pred_mean_HPH,
  pred_mean_logit = pred_mean_logit,
  pred_sd_logit = pred_sd_logit
)

setorder(results, -pred_mean_HPH)
results[, rank := .I]

fwrite(results, file.path(output_dir, "pred_all_49550_100rep_logit_constrained.csv"))
top <- results[1:min(OUTPUT_TOP_N, nrow(results))]
fwrite(top, file.path(output_dir, paste0("top_", OUTPUT_TOP_N, "_logit_constrained.csv")))

cat("Constrained HPH prediction range:", round(range(pred_mean_HPH, na.rm = TRUE), 4), "\n")
cat("Original logit prediction range:", round(range(pred_mean_logit, na.rm = TRUE), 4), "\n")
cat("Done!\n")