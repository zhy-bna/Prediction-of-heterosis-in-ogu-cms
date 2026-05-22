# Script: Calculate additive and dominance PGSI per window

library(data.table)
library(vcfR)
library(parallel)

# defined parameters

vcf_file        <- "/biodata/zhy/pgsi_ad/Brapa_napus.raw.filter.snp.anno.vcf"
ref_genome      <- "/biodata/zhy/pgsi_ad/Brassica_napus.Darmor.v4.1.genome.fa"
samtools_path   <- "/database/biosoft/samtools1.9/samtools"
sample_info_file <- "parent_info_utf8.txt"

window_files <- c(
  "/biodata/zhy/pgsi_ad/100kb_label.csv",
  "/biodata/zhy/pgsi_ad/500kb_label.csv",
  "/biodata/zhy/pgsi_ad/1mb_label.csv"
)
window_sizes <- c(100000, 500000, 1000000)
output_prefix   <- "PGSI_dual_method"
use_parallel    <- TRUE
n_cores         <- 36

keep_chromosomes <- c(
  paste0("chrA", sprintf("%02d", 1:10)),
  paste0("chrC", sprintf("%02d", 1:9)),
  paste0("chrA", sprintf("%02d", 1:10), "_random"),
  paste0("chrC", sprintf("%02d", 1:9), "_random"),
  "chrAnn_random", "chrCnn_random", "chrUnn_random"
)

# read window files

read_window_file <- function(file_path, keep_chr) {
  dt <- fread(file_path, header = FALSE, col.names = "window_id", fill = TRUE)
  dt <- dt[!is.na(window_id) & window_id != ""]
  dt[, order_idx := .I]
  dt[, n_parts := lengths(strsplit(window_id, "_", fixed = TRUE))]
  
  dt[n_parts == 3, c("chr", "start_str", "end_str") := {
    parts <- tstrsplit(window_id, "_", fixed = TRUE)
    list(parts[[1]], parts[[2]], parts[[3]])
  }]
  
  dt[n_parts == 4, c("chr", "start_str", "end_str") := {
    parts <- tstrsplit(window_id, "_", fixed = TRUE)
    list(paste0(parts[[1]], "_", parts[[2]]), parts[[3]], parts[[4]])
  }]
  
  dt[, start := as.numeric(start_str)]
  dt[, end   := as.numeric(end_str)]
  dt <- dt[!is.na(start) & !is.na(end) & !is.na(chr)]
  dt <- dt[chr %in% keep_chr]
  setorder(dt, order_idx)
  windows <- dt[, .(chr, start, end, window_id, order_idx)]
  return(windows)
}


# Assign SNPs to windows

assign_snps_to_windows <- function(snp_pos, windows) {
  snp_pos_copy <- copy(snp_pos)
  snp_pos_copy[, pos_end := pos]
  setkey(snp_pos_copy, chr, pos, pos_end)
  setkey(windows, chr, start, end)
  foverlaps(snp_pos_copy, windows,
            by.x = c("chr", "pos", "pos_end"),
            by.y = c("chr", "start", "end"),
            type = "within")
}


#compute additive and dominance PGSI per window

compute_window_stats <- function(wid, snp_in_window, gt_mat, cms_ids, restorer_ids) {
  snp_idx <- snp_in_window[window_id == wid, snp_id]
  if (length(snp_idx) == 0) return(NULL)
  
  gt_win  <- gt_mat[snp_idx, , drop = FALSE]
  n_sites <- nrow(gt_win)
  
  n_cms  <- length(cms_ids)
  n_rest <- length(restorer_ids)
  
  cms_cols  <- match(cms_ids, colnames(gt_win))
  rest_cols <- match(restorer_ids, colnames(gt_win))
  if (any(is.na(cms_cols)) || any(is.na(rest_cols))) return(NULL)
  
  count_00_mat  <- matrix(0L, nrow = n_cms, ncol = n_rest)
  count_11_mat  <- matrix(0L, nrow = n_cms, ncol = n_rest)
  count_dom_mat <- matrix(0L, nrow = n_cms, ncol = n_rest)
  count_valid_mat <- matrix(0L, nrow = n_cms, ncol = n_rest)
  
  gt_int <- as.matrix(gt_win)
  
  for (i in seq_len(n_sites)) {
    cms_vals <- gt_int[i, cms_cols]
    rest_vals <- gt_int[i, rest_cols]
    
    cms_valid <- !is.na(cms_vals)
    rest_valid <- !is.na(rest_vals)
    valid_mat <- outer(cms_valid, rest_valid, `&`)
    
    cms_is_0 <- cms_vals == 0L
    cms_is_1 <- cms_vals == 1L
    cms_is_2 <- cms_vals == 2L
    rest_is_0 <- rest_vals == 0L
    rest_is_1 <- rest_vals == 1L
    rest_is_2 <- rest_vals == 2L
    
    cms_is_0[is.na(cms_is_0)] <- FALSE
    cms_is_1[is.na(cms_is_1)] <- FALSE
    cms_is_2[is.na(cms_is_2)] <- FALSE
    rest_is_0[is.na(rest_is_0)] <- FALSE
    rest_is_1[is.na(rest_is_1)] <- FALSE
    rest_is_2[is.na(rest_is_2)] <- FALSE
    
    mat_00 <- outer(cms_is_0, rest_is_0, `&`) & valid_mat
    mat_11 <- outer(cms_is_2, rest_is_2, `&`) & valid_mat
    mat_dom <- (
      outer(cms_is_0, rest_is_1, `&`) |
        outer(cms_is_1, rest_is_0, `&`) |
        outer(cms_is_1, rest_is_2, `&`) |
        outer(cms_is_2, rest_is_1, `&`)
    ) & valid_mat
    
    count_00_mat   <- count_00_mat + mat_00
    count_11_mat   <- count_11_mat + mat_11
    count_dom_mat  <- count_dom_mat + mat_dom
    count_valid_mat <- count_valid_mat + valid_mat
  }
  
  res <- data.table(
    window_id = wid,
    cms       = rep(cms_ids, each = n_rest),
    restorer  = rep(restorer_ids, times = n_cms),
    add_PGSI  = as.vector((count_00_mat + count_11_mat) / count_valid_mat),
    dom_PGSI  = as.vector(count_dom_mat / count_valid_mat)
  )
  
  res[is.infinite(add_PGSI) | is.nan(add_PGSI), add_PGSI := NA_real_]
  res[is.infinite(dom_PGSI) | is.nan(dom_PGSI), dom_PGSI := NA_real_]
  return(res)
}


# Create .fai index if needed
fai_file <- paste0(ref_genome, ".fai")
if (!file.exists(fai_file)) {
  system(paste(samtools_path, "faidx", ref_genome))
}

# Read sample info
sample_info <- fread(sample_info_file, encoding = "UTF-8")
cms_ids      <- sample_info[Type == "CMS", SampleID]
restorer_ids <- sample_info[Type == "Restorer", SampleID]

# Read VCF and extract genotypes
vcf    <- read.vcfR(vcf_file, verbose = FALSE)
gt_raw <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
snp_pos <- data.table(
  chr    = vcf@fix[, "CHROM"],
  pos    = as.integer(vcf@fix[, "POS"]),
  snp_id = seq_len(nrow(vcf))
)

# Convert to numeric coding
gt_num <- apply(gt_raw, 2, function(x) {
  x <- gsub("\\|", "/", x)
  code <- rep(NA_integer_, length(x))
  code[x == "0/0"] <- 0L
  code[x == "0/1" | x == "1/0"] <- 1L
  code[x == "1/1"] <- 2L
  code
})
rownames(gt_num) <- snp_pos$snp_id

# Process each window size
for (idx in seq_along(window_files)) {
  win_file <- window_files[idx]
  ws <- window_sizes[idx]
  cat("Processing window size:", ws, "bp\n")
  
  windows <- read_window_file(win_file, keep_chromosomes)
  if (nrow(windows) == 0) next
  
  # Filter windows by chromosomes present in VCF
  windows <- windows[chr %in% unique(snp_pos$chr)]
  if (nrow(windows) == 0) next
  
  snp_in_window <- assign_snps_to_windows(snp_pos, windows)
  windows_with_snps <- unique(snp_in_window$window_id[!is.na(snp_in_window$window_id)])
  if (length(windows_with_snps) == 0) next
  
  # Parallel computation
  if (use_parallel) {
    cl <- makeCluster(n_cores)
    clusterExport(cl, c("snp_in_window", "gt_num", "cms_ids", "restorer_ids", "compute_window_stats"))
    clusterEvalQ(cl, library(data.table))
    result_list <- parLapply(cl, windows_with_snps, function(wid) {
      compute_window_stats(wid, snp_in_window, gt_num, cms_ids, restorer_ids)
    })
    stopCluster(cl)
  } else {
    result_list <- lapply(windows_with_snps, function(wid) {
      compute_window_stats(wid, snp_in_window, gt_num, cms_ids, restorer_ids)
    })
  }
  
  result_list <- result_list[!sapply(result_list, is.null)]
  if (length(result_list) == 0) next
  all_stats <- rbindlist(result_list)
  
  # Convert to wide format
  all_stats[, hybrid_id := paste(cms, restorer, sep = "_")]
  add_wide <- dcast(all_stats, hybrid_id ~ window_id, value.var = "add_PGSI", fill = NA)
  dom_wide <- dcast(all_stats, hybrid_id ~ window_id, value.var = "dom_PGSI", fill = NA)
  
  # Reorder columns to match original window order
  win_order <- windows$window_id[windows$window_id %in% names(add_wide)]
  setcolorder(add_wide, c("hybrid_id", win_order))
  setcolorder(dom_wide, c("hybrid_id", win_order))
  
  # Add cms and restorer columns
  add_wide[, c("cms", "restorer") := tstrsplit(hybrid_id, "_", fixed = TRUE)]
  dom_wide[, c("cms", "restorer") := tstrsplit(hybrid_id, "_", fixed = TRUE)]
  setcolorder(add_wide, c("hybrid_id", "cms", "restorer", win_order))
  setcolorder(dom_wide, c("hybrid_id", "cms", "restorer", win_order))
  
  # Write output
  out_add <- paste0(output_prefix, "_", ws/1000, "Kb_additive.csv")
  out_dom <- paste0(output_prefix, "_", ws/1000, "Kb_dominance.csv")
  fwrite(add_wide, out_add)
  fwrite(dom_wide, out_dom)
  cat("Saved:", out_add, "and", out_dom, "\n")
}

cat("Done.\n")