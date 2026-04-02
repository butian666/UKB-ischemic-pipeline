options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x) || !nzchar(as.character(x))) y else x

log_msg <- function(...) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s\n", ts, paste(..., collapse = " ")))
}

parse_kv_args <- function(args) {
  out <- list()
  if (length(args) == 0) return(out)
  for (a in args) {
    if (!grepl("=", a, fixed = TRUE)) next
    kv <- strsplit(a, "=", fixed = TRUE)[[1]]
    if (length(kv) < 2) next
    key <- trimws(kv[1])
    val <- trimws(paste(kv[-1], collapse = "="))
    out[[key]] <- val
  }
  out
}

to_num <- function(x) suppressWarnings(as.numeric(x))

require_pkgs <- function(pkgs) {
  ok <- vapply(pkgs, function(p) requireNamespace(p, quietly = TRUE), logical(1))
  if (!all(ok)) stop("缺少R包: ", paste(pkgs[!ok], collapse = ", "))
  invisible(TRUE)
}

normalize_source <- function(x) {
  x <- trimws(as.character(x))
  if (!nzchar(x)) return(list(type = "none", value = ""))
  if (grepl("^opengwas:", x, ignore.case = TRUE)) {
    return(list(type = "opengwas", value = sub("^opengwas:", "", x, ignore.case = TRUE)))
  }
  if (grepl("^file:", x, ignore.case = TRUE)) {
    p <- sub("^file:", "", x, ignore.case = TRUE)
    return(list(type = "local", value = normalizePath(p, winslash = "/", mustWork = FALSE)))
  }
  if (file.exists(x)) return(list(type = "local", value = normalizePath(x, winslash = "/", mustWork = FALSE)))
  if (grepl("\\.gz$|\\.txt$|\\.tsv$|\\.csv$", x, ignore.case = TRUE) || grepl("/", x, fixed = TRUE)) {
    return(list(type = "local", value = normalizePath(x, winslash = "/", mustWork = FALSE)))
  }
  list(type = "opengwas", value = x)
}

read_table_flexible <- function(path) {
  if (!file.exists(path)) stop("文件不存在: ", path)
  ext <- tolower(tools::file_ext(path))
  if (requireNamespace("data.table", quietly = TRUE)) {
    return(data.table::fread(path, data.table = FALSE))
  }
  if (ext %in% c("csv")) return(utils::read.csv(path, check.names = FALSE))
  if (ext %in% c("tsv", "txt")) return(utils::read.delim(path, check.names = FALSE))
  utils::read.table(path, header = TRUE, sep = "", check.names = FALSE)
}

find_col <- function(nms, candidates) {
  hit <- candidates[candidates %in% nms]
  if (length(hit) > 0) return(hit[[1]])
  lower <- tolower(nms)
  idx <- match(tolower(candidates), lower)
  idx <- idx[!is.na(idx)]
  if (length(idx) > 0) return(nms[idx[[1]]])
  NA_character_
}

standardize_sumstats <- function(df, trait) {
  nms <- names(df)
  c_snp <- find_col(nms, c("SNP", "snp", "rsid", "rsID", "variant_id", "MarkerName"))
  c_beta <- find_col(nms, c("beta", "b", "BETA", "effect", "Beta"))
  c_se <- find_col(nms, c("se", "SE", "stderr", "StdErr"))
  c_ea <- find_col(nms, c("effect_allele", "EA", "A1", "alt", "allele1"))
  c_oa <- find_col(nms, c("other_allele", "OA", "A2", "ref", "allele2"))
  c_eaf <- find_col(nms, c("eaf", "EAF", "effect_allele_freq", "af", "ALT_FREQS"))
  c_p <- find_col(nms, c("pval", "p_value", "P", "PVAL", "p", "PVALUE"))
  c_n <- find_col(nms, c("n", "N", "samplesize", "sample_size", "Neff"))
  c_chr <- find_col(nms, c("chr", "CHR", "chromosome", "#CHROM"))
  c_pos <- find_col(nms, c("pos", "BP", "position", "base_pair_location", "POS"))
  req <- c(c_snp, c_beta, c_se, c_ea, c_oa, c_p)
  if (any(is.na(req))) stop("本地GWAS文件缺少必要列: SNP/beta/se/effect_allele/other_allele/pval")
  out <- data.frame(
    SNP = as.character(df[[c_snp]]),
    beta = to_num(df[[c_beta]]),
    se = to_num(df[[c_se]]),
    effect_allele = as.character(df[[c_ea]]),
    other_allele = as.character(df[[c_oa]]),
    eaf = if (!is.na(c_eaf)) to_num(df[[c_eaf]]) else NA_real_,
    pval = to_num(df[[c_p]]),
    samplesize = if (!is.na(c_n)) to_num(df[[c_n]]) else NA_real_,
    chr = if (!is.na(c_chr)) as.character(df[[c_chr]]) else NA_character_,
    pos = if (!is.na(c_pos)) to_num(df[[c_pos]]) else NA_real_,
    trait = trait,
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$SNP) & nzchar(out$SNP), , drop = FALSE]
  out <- out[!is.na(out$beta) & !is.na(out$se) & !is.na(out$pval), , drop = FALSE]
  out
}

format_exposure_local <- function(df, trait) {
  TwoSampleMR::format_data(
    df,
    type = "exposure",
    phenotype_col = "trait",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "eaf",
    pval_col = "pval",
    samplesize_col = "samplesize",
    chr_col = "chr",
    pos_col = "pos"
  )
}

format_outcome_local <- function(df, trait) {
  TwoSampleMR::format_data(
    df,
    type = "outcome",
    phenotype_col = "trait",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "eaf",
    pval_col = "pval",
    samplesize_col = "samplesize",
    chr_col = "chr",
    pos_col = "pos"
  )
}

extract_instruments_from_source <- function(source_spec, trait, p_threshold, clump_r2, clump_kb, clump_pop) {
  if (identical(source_spec$type, "opengwas")) {
    dat <- TwoSampleMR::extract_instruments(
      outcomes = source_spec$value,
      p1 = p_threshold,
      clump = TRUE,
      r2 = clump_r2,
      kb = clump_kb
    )
    return(dat)
  }
  if (identical(source_spec$type, "local")) {
    raw <- read_table_flexible(source_spec$value)
    std <- standardize_sumstats(raw, trait)
    std <- std[!is.na(std$pval) & std$pval <= p_threshold, , drop = FALSE]
    if (nrow(std) == 0) stop("本地数据中无达到工具变量阈值的SNP: ", trait)
    exp_dat <- format_exposure_local(std, trait)
    clumped <- TwoSampleMR::clump_data(
      exp_dat,
      clump_kb = clump_kb,
      clump_r2 = clump_r2,
      pop = clump_pop
    )
    return(clumped)
  }
  stop("无效source: ", source_spec$type)
}

extract_outcome_from_source <- function(snps, source_spec, trait) {
  if (identical(source_spec$type, "opengwas")) {
    dat <- TwoSampleMR::extract_outcome_data(snps = snps, outcomes = source_spec$value)
    return(dat)
  }
  if (identical(source_spec$type, "local")) {
    raw <- read_table_flexible(source_spec$value)
    std <- standardize_sumstats(raw, trait)
    std <- std[std$SNP %in% snps, , drop = FALSE]
    if (nrow(std) == 0) return(std)
    out_dat <- format_outcome_local(std, trait)
    return(out_dat)
  }
  stop("无效source: ", source_spec$type)
}

safe_bind_rows <- function(lst) {
  lst <- lst[!vapply(lst, is.null, logical(1))]
  if (length(lst) == 0) return(NULL)
  dplyr::bind_rows(lst)
}

run_one_pair <- function(exp_trait, exp_source, out_trait, out_source, direction, p_threshold, clump_r2, clump_kb, clump_pop, enable_presso = FALSE, enable_detailed = FALSE) {
  meta <- data.frame(
    direction = direction,
    exposure_trait = exp_trait,
    outcome_trait = out_trait,
    exposure_source = paste0(exp_source$type, ":", exp_source$value),
    outcome_source = paste0(out_source$type, ":", out_source$value),
    stringsAsFactors = FALSE
  )
  inst <- tryCatch(extract_instruments_from_source(exp_source, exp_trait, p_threshold, clump_r2, clump_kb, clump_pop), error = function(e) e)
  if (inherits(inst, "error")) {
    return(list(
      status = cbind(meta, status = "failed_exposure_instruments", message = inst$message, stringsAsFactors = FALSE),
      mr = NULL, heterogeneity = NULL, pleiotropy = NULL, harmonised = NULL, leaveoneout = NULL, singlesnp = NULL, presso = NULL
    ))
  }
  if (is.null(inst) || nrow(inst) == 0) {
    return(list(
      status = cbind(meta, status = "no_instruments", message = "工具变量为空", stringsAsFactors = FALSE),
      mr = NULL, heterogeneity = NULL, pleiotropy = NULL, harmonised = NULL, leaveoneout = NULL, singlesnp = NULL, presso = NULL
    ))
  }
  out_dat <- tryCatch(extract_outcome_from_source(unique(inst$SNP), out_source, out_trait), error = function(e) e)
  if (inherits(out_dat, "error")) {
    return(list(
      status = cbind(meta, status = "failed_outcome_extract", message = out_dat$message, stringsAsFactors = FALSE),
      mr = NULL, heterogeneity = NULL, pleiotropy = NULL, harmonised = NULL, leaveoneout = NULL, singlesnp = NULL, presso = NULL
    ))
  }
  if (is.null(out_dat) || nrow(out_dat) == 0) {
    return(list(
      status = cbind(meta, status = "no_outcome_overlap", message = "结果数据与工具变量无重叠SNP", stringsAsFactors = FALSE),
      mr = NULL, heterogeneity = NULL, pleiotropy = NULL, harmonised = NULL, leaveoneout = NULL, singlesnp = NULL, presso = NULL
    ))
  }
  harm <- tryCatch(TwoSampleMR::harmonise_data(inst, out_dat, action = 2), error = function(e) e)
  if (inherits(harm, "error")) {
    return(list(
      status = cbind(meta, status = "failed_harmonise", message = harm$message, stringsAsFactors = FALSE),
      mr = NULL, heterogeneity = NULL, pleiotropy = NULL, harmonised = NULL, leaveoneout = NULL, singlesnp = NULL, presso = NULL
    ))
  }
  harm <- harm[harm$mr_keep, , drop = FALSE]
  if (nrow(harm) == 0) {
    return(list(
      status = cbind(meta, status = "no_harmonised_snps", message = "协调后无可用SNP", stringsAsFactors = FALSE),
      mr = NULL, heterogeneity = NULL, pleiotropy = NULL, harmonised = NULL, leaveoneout = NULL, singlesnp = NULL, presso = NULL
    ))
  }
  method_list <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median")
  mr_res <- tryCatch(TwoSampleMR::mr(harm, method_list = method_list), error = function(e) e)
  if (inherits(mr_res, "error")) {
    return(list(
      status = cbind(meta, status = "failed_mr", message = mr_res$message, stringsAsFactors = FALSE),
      mr = NULL, heterogeneity = NULL, pleiotropy = NULL, harmonised = harm, leaveoneout = NULL, singlesnp = NULL, presso = NULL
    ))
  }
  hetero <- tryCatch(TwoSampleMR::mr_heterogeneity(harm), error = function(e) NULL)
  pleio <- tryCatch(TwoSampleMR::mr_pleiotropy_test(harm), error = function(e) NULL)
  loo <- if (isTRUE(enable_detailed)) tryCatch(TwoSampleMR::mr_leaveoneout(harm), error = function(e) NULL) else NULL
  snp_res <- if (isTRUE(enable_detailed)) tryCatch(TwoSampleMR::mr_singlesnp(harm), error = function(e) NULL) else NULL
  presso_df <- NULL
  if (isTRUE(enable_presso) && requireNamespace("MRPRESSO", quietly = TRUE) && nrow(harm) >= 4) {
    presso_obj <- tryCatch(
      MRPRESSO::mr_presso(
        BetaOutcome = "beta.outcome",
        BetaExposure = "beta.exposure",
        SdOutcome = "se.outcome",
        SdExposure = "se.exposure",
        OUTLIERtest = TRUE,
        DISTORTIONtest = TRUE,
        data = harm,
        NbDistribution = 1000,
        SignifThreshold = 0.05
      ),
      error = function(e) e
    )
    if (!inherits(presso_obj, "error")) {
      presso_df <- data.frame(
        direction = direction,
        exposure_trait = exp_trait,
        outcome_trait = out_trait,
        nsnp = nrow(harm),
        presso_available = TRUE,
        presso_summary = paste(capture.output(str(presso_obj, max.level = 1)), collapse = " | "),
        stringsAsFactors = FALSE
      )
    } else {
      presso_df <- data.frame(
        direction = direction,
        exposure_trait = exp_trait,
        outcome_trait = out_trait,
        nsnp = nrow(harm),
        presso_available = FALSE,
        presso_summary = presso_obj$message,
        stringsAsFactors = FALSE
      )
    }
  }
  status_ok <- cbind(meta, status = "ok", message = "", stringsAsFactors = FALSE)
  add_meta <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    cbind(meta[rep(1, nrow(df)), , drop = FALSE], df, stringsAsFactors = FALSE)
  }
  list(
    status = status_ok,
    mr = add_meta(mr_res),
    heterogeneity = add_meta(hetero),
    pleiotropy = add_meta(pleio),
    harmonised = add_meta(harm),
    leaveoneout = add_meta(loo),
    singlesnp = add_meta(snp_res),
    presso = presso_df
  )
}

build_brain_catalog <- function(catalog_path) {
  if (!file.exists(catalog_path)) stop("脑IDP目录文件不存在: ", catalog_path)
  cat_df <- read_table_flexible(catalog_path)
  req <- c("trait", "source")
  miss <- setdiff(req, names(cat_df))
  if (length(miss) > 0) stop("脑IDP目录缺少字段: ", paste(miss, collapse = ", "))
  cat_df <- cat_df[, intersect(c("trait", "source", "cohort"), names(cat_df)), drop = FALSE]
  cat_df$trait <- as.character(cat_df$trait)
  cat_df$source <- as.character(cat_df$source)
  if (!("cohort" %in% names(cat_df))) cat_df$cohort <- "UKB_BIG40"
  cat_df <- cat_df[!is.na(cat_df$trait) & nzchar(cat_df$trait) & !is.na(cat_df$source) & nzchar(cat_df$source), , drop = FALSE]
  if (nrow(cat_df) == 0) stop("脑IDP目录为空")
  cat_df
}

args <- parse_kv_args(commandArgs(trailingOnly = TRUE))
out_dir <- normalizePath(args$out_dir %||% file.path(getwd(), "Output_Tables", "MR_Bidirectional"), winslash = "/", mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ihd_finngen <- args$ihd_finngen %||% Sys.getenv("IHD_FINNGEN_SOURCE", "opengwas:finn-b-I9_IHD")
ihd_ukb <- args$ihd_ukb %||% Sys.getenv("IHD_UKB_SOURCE", "")
brain_catalog_path <- args$brain_catalog %||% Sys.getenv("BRAIN_CATALOG_PATH", "")
p_threshold <- to_num(args$p_threshold %||% Sys.getenv("MR_P_THRESHOLD", "5e-8"))
if (is.na(p_threshold) || p_threshold <= 0) p_threshold <- 5e-8
clump_r2 <- to_num(args$clump_r2 %||% Sys.getenv("MR_CLUMP_R2", "0.001"))
if (is.na(clump_r2) || clump_r2 <= 0) clump_r2 <- 0.001
clump_kb <- as.integer(to_num(args$clump_kb %||% Sys.getenv("MR_CLUMP_KB", "10000")))
if (is.na(clump_kb) || clump_kb <= 0) clump_kb <- 10000L
clump_pop <- args$clump_pop %||% Sys.getenv("MR_CLUMP_POP", "EUR")
enable_presso <- isTRUE(as.logical(args$enable_presso %||% Sys.getenv("MR_ENABLE_PRESSO", "FALSE")))
enable_detailed <- isTRUE(as.logical(args$enable_detailed %||% Sys.getenv("MR_ENABLE_DETAILED", "FALSE")))

require_pkgs(c("TwoSampleMR", "ieugwasr", "dplyr"))

if (!nzchar(brain_catalog_path)) stop("请提供 brain_catalog 参数（CSV/TSV，至少含 trait 与 source 列）")

brain_catalog <- build_brain_catalog(brain_catalog_path)
ihd_catalog <- data.frame(
  trait = c("IHD_FinnGen", "IHD_UKB"),
  source = c(ihd_finngen, ihd_ukb),
  cohort = c("FinnGen", "UKB"),
  stringsAsFactors = FALSE
)
ihd_catalog <- ihd_catalog[!is.na(ihd_catalog$source) & nzchar(ihd_catalog$source), , drop = FALSE]
if (nrow(ihd_catalog) == 0) stop("IHD来源为空，请至少提供 ihd_finngen 或 ihd_ukb")

log_msg("开始双向MR")
log_msg("IHD来源数:", nrow(ihd_catalog), "脑IDP来源数:", nrow(brain_catalog))

all_status <- list()
all_mr <- list()
all_hetero <- list()
all_pleio <- list()
all_harm <- list()
all_loo <- list()
all_snp <- list()
all_presso <- list()

for (i in seq_len(nrow(ihd_catalog))) {
  for (j in seq_len(nrow(brain_catalog))) {
    exp_trait <- ihd_catalog$trait[i]
    exp_src <- normalize_source(ihd_catalog$source[i])
    out_trait <- brain_catalog$trait[j]
    out_src <- normalize_source(brain_catalog$source[j])
    log_msg("运行:", exp_trait, "->", out_trait)
    res <- run_one_pair(
      exp_trait = exp_trait,
      exp_source = exp_src,
      out_trait = out_trait,
      out_source = out_src,
      direction = "heart_to_brain",
      p_threshold = p_threshold,
      clump_r2 = clump_r2,
      clump_kb = clump_kb,
      clump_pop = clump_pop,
      enable_presso = enable_presso,
      enable_detailed = enable_detailed
    )
    all_status[[length(all_status) + 1]] <- res$status
    all_mr[[length(all_mr) + 1]] <- res$mr
    all_hetero[[length(all_hetero) + 1]] <- res$heterogeneity
    all_pleio[[length(all_pleio) + 1]] <- res$pleiotropy
    all_harm[[length(all_harm) + 1]] <- res$harmonised
    all_loo[[length(all_loo) + 1]] <- res$leaveoneout
    all_snp[[length(all_snp) + 1]] <- res$singlesnp
    all_presso[[length(all_presso) + 1]] <- res$presso
  }
}

for (i in seq_len(nrow(brain_catalog))) {
  for (j in seq_len(nrow(ihd_catalog))) {
    exp_trait <- brain_catalog$trait[i]
    exp_src <- normalize_source(brain_catalog$source[i])
    out_trait <- ihd_catalog$trait[j]
    out_src <- normalize_source(ihd_catalog$source[j])
    log_msg("运行:", exp_trait, "->", out_trait)
    res <- run_one_pair(
      exp_trait = exp_trait,
      exp_source = exp_src,
      out_trait = out_trait,
      out_source = out_src,
      direction = "brain_to_heart",
      p_threshold = p_threshold,
      clump_r2 = clump_r2,
      clump_kb = clump_kb,
      clump_pop = clump_pop,
      enable_presso = enable_presso,
      enable_detailed = enable_detailed
    )
    all_status[[length(all_status) + 1]] <- res$status
    all_mr[[length(all_mr) + 1]] <- res$mr
    all_hetero[[length(all_hetero) + 1]] <- res$heterogeneity
    all_pleio[[length(all_pleio) + 1]] <- res$pleiotropy
    all_harm[[length(all_harm) + 1]] <- res$harmonised
    all_loo[[length(all_loo) + 1]] <- res$leaveoneout
    all_snp[[length(all_snp) + 1]] <- res$singlesnp
    all_presso[[length(all_presso) + 1]] <- res$presso
  }
}

status_df <- safe_bind_rows(all_status)
mr_df <- safe_bind_rows(all_mr)
hetero_df <- safe_bind_rows(all_hetero)
pleio_df <- safe_bind_rows(all_pleio)
harm_df <- safe_bind_rows(all_harm)
loo_df <- safe_bind_rows(all_loo)
snp_df <- safe_bind_rows(all_snp)
presso_df <- safe_bind_rows(all_presso)

if (!is.null(status_df)) utils::write.csv(status_df, file.path(out_dir, "MR_Run_Status.csv"), row.names = FALSE)
if (!is.null(mr_df)) {
  mr_df$pval_fdr <- stats::p.adjust(mr_df$pval, method = "BH")
  utils::write.csv(mr_df, file.path(out_dir, "MR_Bidirectional_Results.csv"), row.names = FALSE)
}
if (!is.null(hetero_df)) utils::write.csv(hetero_df, file.path(out_dir, "MR_Heterogeneity.csv"), row.names = FALSE)
if (!is.null(pleio_df)) utils::write.csv(pleio_df, file.path(out_dir, "MR_Pleiotropy_EggerIntercept.csv"), row.names = FALSE)
if (!is.null(harm_df)) utils::write.csv(harm_df, file.path(out_dir, "MR_Harmonised_SNPs.csv"), row.names = FALSE)
if (!is.null(loo_df)) utils::write.csv(loo_df, file.path(out_dir, "MR_LeaveOneOut.csv"), row.names = FALSE)
if (!is.null(snp_df)) utils::write.csv(snp_df, file.path(out_dir, "MR_SingleSNP.csv"), row.names = FALSE)
if (!is.null(presso_df)) utils::write.csv(presso_df, file.path(out_dir, "MR_PRESSO_Summary.csv"), row.names = FALSE)

summary_df <- data.frame(
  total_pairs = if (!is.null(status_df)) nrow(status_df) else 0L,
  successful_pairs = if (!is.null(status_df)) sum(status_df$status == "ok", na.rm = TRUE) else 0L,
  mr_rows = if (!is.null(mr_df)) nrow(mr_df) else 0L,
  stringsAsFactors = FALSE
)
utils::write.csv(summary_df, file.path(out_dir, "MR_Analysis_Summary.csv"), row.names = FALSE)
log_msg("完成。输出目录:", out_dir)
