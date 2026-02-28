args <- commandArgs(trailingOnly = TRUE)

excel_path <- if (length(args) >= 1) args[[1]] else "IDPs.xlsx"
sheet <- if (length(args) >= 2 && nzchar(args[[2]])) args[[2]] else NA_character_
default_out <- sub("\\.xlsx$", "_filled.xlsx", excel_path, ignore.case = TRUE)
out_xlsx <- if (length(args) >= 3 && nzchar(args[[3]])) args[[3]] else default_out

read_any_table <- function(path, sheet = NA_character_) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("需要 R 包 readxl 才能读取 Excel：请先安装 readxl，或将文件另存为 CSV。")
    }
    if (is.na(sheet)) readxl::read_excel(path, col_names = FALSE) else readxl::read_excel(path, sheet = sheet, col_names = FALSE)
  } else if (ext %in% c("csv", "tsv")) {
    sep <- if (ext == "tsv") "\t" else ","
    utils::read.csv(path, stringsAsFactors = FALSE, sep = sep, check.names = FALSE)
  } else {
    stop("不支持的输入文件格式：", ext)
  }
}

detect_col <- function(nms, candidates) {
  hit <- candidates[candidates %in% nms]
  if (length(hit) > 0) return(hit[[1]])
  lower <- tolower(nms)
  cand_lower <- tolower(candidates)
  idx <- match(cand_lower, lower)
  idx <- idx[!is.na(idx)]
  if (length(idx) > 0) return(nms[idx[[1]]])
  NA_character_
}

rstrip_dots <- function(x) sub("\\.*$", "", as.character(x))

req <- read_any_table(excel_path, sheet = sheet)
req <- as.data.frame(req, stringsAsFactors = FALSE)
if (nrow(req) == 0) stop("Excel/CSV 里没有任何行：", excel_path)
preview_txt <- paste(as.character(unlist(req[seq_len(min(nrow(req), 20)), , drop = FALSE])), collapse = " ")

normalize_header_row <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  scan_n <- min(nrow(df), 20)
  header_tokens <- c("variable name", "variable_name", "idps", "idp", "idp_root", "%", "se", "z", "p_raw", "p_correct", "n")
  header_row <- NA_integer_
  for (i in seq_len(scan_n)) {
    r <- as.character(df[i, , drop = TRUE])
    rc <- tolower(trimws(r))
    hits <- sum(rc %in% header_tokens, na.rm = TRUE)
    if (hits >= 2) {
      header_row <- i
      break
    }
  }
  if (is.na(header_row)) return(df)

  new_names <- as.character(df[header_row, , drop = TRUE])
  bad <- is.na(new_names) | !nzchar(trimws(new_names))
  if (any(bad)) new_names[bad] <- names(df)[bad]

  out <- df[(header_row + 1):nrow(df), , drop = FALSE]
  names(out) <- make.unique(trimws(new_names))

  keep_cols <- vapply(out, function(col) {
    any(!is.na(col) & nzchar(trimws(as.character(col))))
  }, logical(1))
  out <- out[, keep_cols, drop = FALSE]
  out
}

infer_model_from_text <- function(txt) {
  txt <- tolower(paste(txt, collapse = " "))
  if (grepl("\\bis\\s*vs\\s*control\\b", txt)) return("ischemic_vs_control")
  if (grepl("mi\\s*vs\\s*chronic", txt)) return("mi_vs_chronic")
  if (grepl("chronic\\s*vs\\s*control", txt)) return("chronic_vs_control")
  if (grepl("mi\\s*vs\\s*control", txt)) return("mi_vs_control")
  if (grepl("(is|ischemic)\\s*vs\\s*control", txt)) return("ischemic_vs_control")
  NA_character_
}

req <- normalize_header_row(req)

req_idp_col <- detect_col(names(req), c("idp_root", "variable_name", "variable name", "IDPS", "IDP", "idp", "idp_name", "idp name", "name"))
if (is.na(req_idp_col)) stop("未在 Excel/CSV 中找到 IDP 列（期望列名如：idp_root / variable_name / IDP）。")

req_pct_col <- detect_col(names(req), c("%", "percent", "percentage", "delta_percent_change", "delta", "effect"))
if (is.na(req_pct_col)) stop("未在 Excel/CSV 中找到待填充的 % 列（期望列名如：% / percent）。")

req_model_col <- detect_col(names(req), c("model_comparison", "comparison", "group", "model"))
inferred_model <- infer_model_from_text(c(preview_txt, if (is.na(sheet)) "" else sheet))

req$idp_key <- rstrip_dots(req[[req_idp_col]])

default_dir <- file.path(getwd(), "Output_Tables", "Longitudinal", "Median_DeltaRate")
csv_dir <- if (dir.exists(default_dir)) default_dir else getwd()
csv_files <- list.files(csv_dir, pattern = "^Longitudinal_Median_DeltaRate_AllIDPs_Comparison_.*\\.csv$", full.names = FALSE)
if (length(csv_files) == 0 && !identical(csv_dir, getwd())) {
  csv_dir <- getwd()
  csv_files <- list.files(csv_dir, pattern = "^Longitudinal_Median_DeltaRate_AllIDPs_Comparison_.*\\.csv$", full.names = FALSE)
}
if (length(csv_files) == 0) stop("未找到 Longitudinal_Median_DeltaRate_AllIDPs_Comparison_*.csv（默认目录：Output_Tables/Longitudinal/Median_DeltaRate 或当前目录）")

models <- sub("^Longitudinal_Median_DeltaRate_AllIDPs_Comparison_", "", sub("\\.csv$", "", basename(csv_files)))

expand_requests <- function(req_df, models, model_col, inferred_model) {
  if (!is.na(model_col) && model_col %in% names(req_df)) {
    req_df$model_comparison_request <- as.character(req_df[[model_col]])
    req_df
  } else if (!is.na(inferred_model) && nzchar(inferred_model)) {
    req_df$model_comparison_request <- inferred_model
    req_df
  } else {
    expanded <- merge(
      req_df,
      data.frame(model_comparison_request = models, stringsAsFactors = FALSE),
      by = NULL
    )
    expanded
  }
}

req2 <- expand_requests(req, models, req_model_col, inferred_model)

fill_one_model <- function(req_df, model_label) {
  fp <- file.path(csv_dir, paste0("Longitudinal_Median_DeltaRate_AllIDPs_Comparison_", model_label, ".csv"))
  if (!file.exists(fp)) stop("找不到模型对应的对比CSV：", fp)
  dat <- utils::read.csv(fp, stringsAsFactors = FALSE, check.names = FALSE)
  if (!("idp_root" %in% names(dat))) stop("CSV 缺少 idp_root 列：", fp)
  if (!("median_diff_case_minus_control" %in% names(dat))) stop("CSV 缺少 median_diff_case_minus_control 列：", fp)
  dat$idp_key <- rstrip_dots(dat$idp_root)
  map <- dat[, c("idp_key", "median_diff_case_minus_control"), drop = FALSE]
  map <- map[!duplicated(map$idp_key), , drop = FALSE]
  idx <- match(req_df$idp_key, map$idp_key)
  req_df[[req_pct_col]] <- suppressWarnings(as.numeric(map$median_diff_case_minus_control[idx]))
  req_df
}

filled_list <- list()
for (m in unique(as.character(req2$model_comparison_request))) {
  sub_req <- req2[req2$model_comparison_request == m, , drop = FALSE]
  if (nrow(sub_req) == 0) next
  filled_list[[m]] <- fill_one_model(sub_req, m)
}
filled <- do.call(rbind, filled_list[!vapply(filled_list, is.null, logical(1))])
if (is.null(filled) || nrow(filled) == 0) stop("未能填充任何条目：请检查模型标签与 CSV 文件名是否匹配。")

filled <- filled[, setdiff(names(filled), c("idp_key")), drop = FALSE]

if (tolower(tools::file_ext(out_xlsx)) %in% c("xlsx", "xls")) {
  if (!requireNamespace("writexl", quietly = TRUE)) {
    stop("需要 R 包 writexl 才能写出 Excel：请先安装 writexl，或改用输出为 CSV。")
  }
  sheet_name <- if (is.na(sheet)) "Sheet1" else sheet
  writexl::write_xlsx(stats::setNames(list(filled), sheet_name), out_xlsx)
  cat("Saved: ", out_xlsx, "\n", sep = "")
} else {
  utils::write.csv(filled, out_xlsx, row.names = FALSE)
  cat("Saved: ", out_xlsx, "\n", sep = "")
}
