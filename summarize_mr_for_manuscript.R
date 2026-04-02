options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x) || !nzchar(as.character(x))) y else x

args <- commandArgs(trailingOnly = TRUE)
kv <- list()
if (length(args) > 0) {
  for (a in args) {
    if (!grepl("=", a, fixed = TRUE)) next
    sp <- strsplit(a, "=", fixed = TRUE)[[1]]
    if (length(sp) < 2) next
    kv[[trimws(sp[1])]] <- trimws(paste(sp[-1], collapse = "="))
  }
}

in_dir <- kv$in_dir %||% file.path(getwd(), "Output_Tables", "MR_Bidirectional_Bulk")
res_file <- file.path(in_dir, "MR_Bidirectional_Results.csv")
het_file <- file.path(in_dir, "MR_Heterogeneity.csv")
pleio_file <- file.path(in_dir, "MR_Pleiotropy_EggerIntercept.csv")
status_file <- file.path(in_dir, "MR_Run_Status.csv")

if (!file.exists(res_file)) stop("缺少文件: ", res_file)

read_csv_flex <- function(path) {
  if (!file.exists(path)) return(NULL)
  if (requireNamespace("data.table", quietly = TRUE)) return(data.table::fread(path, data.table = FALSE))
  utils::read.csv(path, check.names = FALSE)
}

res <- read_csv_flex(res_file)
het <- read_csv_flex(het_file)
pleio <- read_csv_flex(pleio_file)
status <- read_csv_flex(status_file)

res <- res[!is.na(res$method), , drop = FALSE]
res$method <- as.character(res$method)

ivw <- res[res$method == "Inverse variance weighted", , drop = FALSE]
egger <- res[res$method == "MR Egger", , drop = FALSE]
wm <- res[res$method == "Weighted median", , drop = FALSE]

pair_key <- function(df) paste(df$direction, df$exposure_trait, df$outcome_trait, sep = "||")
ivw$key <- pair_key(ivw)
egger$key <- pair_key(egger)
wm$key <- pair_key(wm)

sign_label <- function(p, fdr) {
  if (is.na(p)) return("failed")
  if (!is.na(fdr) && fdr < 0.05) return("significant_fdr")
  if (p < 0.05) return("nominal")
  if (p < 0.10) return("borderline")
  "ns"
}

all_keys <- unique(ivw$key)
rows <- lapply(all_keys, function(k) {
  i <- ivw[ivw$key == k, , drop = FALSE]
  e <- egger[egger$key == k, , drop = FALSE]
  w <- wm[wm$key == k, , drop = FALSE]
  bvals <- c(i$b[1], e$b[1], w$b[1])
  bvals <- bvals[is.finite(bvals)]
  same_dir <- if (length(bvals) >= 2) {
    s <- sign(bvals)
    sum(s > 0) >= 2 || sum(s < 0) >= 2
  } else {
    NA
  }
  het_p <- NA_real_
  if (!is.null(het) && nrow(het) > 0) {
    h <- het[het$direction == i$direction[1] & het$exposure_trait == i$exposure_trait[1] & het$outcome_trait == i$outcome_trait[1] & het$method == "Inverse variance weighted", , drop = FALSE]
    if (nrow(h) > 0) het_p <- as.numeric(h$Q_pval[1])
  }
  egger_int_p <- NA_real_
  if (!is.null(pleio) && nrow(pleio) > 0) {
    p <- pleio[pleio$direction == i$direction[1] & pleio$exposure_trait == i$exposure_trait[1] & pleio$outcome_trait == i$outcome_trait[1], , drop = FALSE]
    if (nrow(p) > 0) egger_int_p <- as.numeric(p$pval[1])
  }
  sensitivity <- if (is.na(het_p) || is.na(egger_int_p)) {
    "incomplete"
  } else if (het_p > 0.05 && egger_int_p > 0.05) {
    "pass"
  } else {
    "caution"
  }
  data.frame(
    direction = i$direction[1],
    exposure_trait = i$exposure_trait[1],
    outcome_trait = i$outcome_trait[1],
    nsnp_ivw = as.numeric(i$nsnp[1]),
    ivw_b = as.numeric(i$b[1]),
    ivw_se = as.numeric(i$se[1]),
    ivw_p = as.numeric(i$pval[1]),
    ivw_fdr = as.numeric(i$pval_fdr[1]),
    ivw_significance = sign_label(as.numeric(i$pval[1]), as.numeric(i$pval_fdr[1])),
    egger_b = if (nrow(e) > 0) as.numeric(e$b[1]) else NA_real_,
    egger_p = if (nrow(e) > 0) as.numeric(e$pval[1]) else NA_real_,
    weighted_median_b = if (nrow(w) > 0) as.numeric(w$b[1]) else NA_real_,
    weighted_median_p = if (nrow(w) > 0) as.numeric(w$pval[1]) else NA_real_,
    direction_consistent_2of3 = same_dir,
    heterogeneity_q_p = het_p,
    egger_intercept_p = egger_int_p,
    sensitivity_flag = sensitivity,
    stringsAsFactors = FALSE
  )
})

pair_tbl <- do.call(rbind, rows)
pair_tbl <- pair_tbl[order(pair_tbl$ivw_p), , drop = FALSE]

if (!is.null(status) && nrow(status) > 0) {
  ok <- status[status$status == "ok", , drop = FALSE]
  fail <- status[status$status != "ok", , drop = FALSE]
  status_summary <- data.frame(
    total_pairs = nrow(status),
    ok_pairs = nrow(ok),
    failed_pairs = nrow(fail),
    stringsAsFactors = FALSE
  )
  utils::write.csv(status_summary, file.path(in_dir, "MR_Manuscript_Status_Summary.csv"), row.names = FALSE)
}

count_tbl <- data.frame(
  significant_fdr = sum(pair_tbl$ivw_significance == "significant_fdr", na.rm = TRUE),
  nominal = sum(pair_tbl$ivw_significance == "nominal", na.rm = TRUE),
  borderline = sum(pair_tbl$ivw_significance == "borderline", na.rm = TRUE),
  ns = sum(pair_tbl$ivw_significance == "ns", na.rm = TRUE),
  direction_consistent_2of3 = sum(pair_tbl$direction_consistent_2of3 %in% TRUE, na.rm = TRUE),
  sensitivity_pass = sum(pair_tbl$sensitivity_flag == "pass", na.rm = TRUE),
  stringsAsFactors = FALSE
)

utils::write.csv(pair_tbl, file.path(in_dir, "MR_Manuscript_Main_Table.csv"), row.names = FALSE)
utils::write.csv(count_tbl, file.path(in_dir, "MR_Manuscript_Summary_Counts.csv"), row.names = FALSE)

cat(sprintf("pair_rows=%d\n", nrow(pair_tbl)))
cat(sprintf("significant_fdr=%d\n", count_tbl$significant_fdr))
cat(sprintf("nominal=%d\n", count_tbl$nominal))
cat(sprintf("borderline=%d\n", count_tbl$borderline))
