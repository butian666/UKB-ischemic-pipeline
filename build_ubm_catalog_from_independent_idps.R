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

read_csv_flex <- function(path) {
  if (requireNamespace("data.table", quietly = TRUE)) return(data.table::fread(path, data.table = FALSE))
  utils::read.csv(path, check.names = FALSE)
}

normalize_text <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("[^a-z0-9]+", " ", x)
  x <- trimws(gsub("\\s+", " ", x))
  x
}

extract_tokens <- function(x) {
  stop_words <- c(
    "mean","weighted","of","in","on","the","and","for","with","from","to","a","an",
    "volume","area","thickness","left","right","hemisphere","tract","idp","t1","dmri",
    "fa","md","mo","od","l1","l2","l3","icvf","isovf","grey","gray","white","matter",
    "skeleton","normalised","normalized","whole","brain","gyrus","sulcus","cortex"
  )
  x <- normalize_text(x)
  toks <- unlist(strsplit(x, " ", fixed = TRUE))
  toks <- toks[nchar(toks) >= 2]
  toks <- toks[!toks %in% stop_words]
  unique(toks)
}

has_word <- function(x, w) grepl(paste0("\\b", w, "\\b"), x)

score_pair <- function(pheno_text, trait_text) {
  p_norm <- normalize_text(pheno_text)
  t_norm <- normalize_text(trait_text)
  p_toks <- extract_tokens(pheno_text)
  t_toks <- extract_tokens(trait_text)
  inter <- length(intersect(p_toks, t_toks))
  union_n <- length(unique(c(p_toks, t_toks)))
  jacc <- if (union_n > 0) inter / union_n else 0
  hemi_bonus <- 0
  if (has_word(p_norm, "left") && has_word(t_norm, "left")) hemi_bonus <- hemi_bonus + 0.03
  if (has_word(p_norm, "right") && has_word(t_norm, "right")) hemi_bonus <- hemi_bonus + 0.03
  if (has_word(p_norm, "left") && has_word(t_norm, "right")) hemi_bonus <- hemi_bonus - 0.05
  if (has_word(p_norm, "right") && has_word(t_norm, "left")) hemi_bonus <- hemi_bonus - 0.05
  mode_bonus <- 0
  mode_keys <- c("thickness","volume","area","fa","md","od","mo","icvf","isovf","l1","l2","l3","hippocampus","cingulate","amygdala","thalam")
  for (k in mode_keys) {
    if (grepl(k, p_norm, fixed = TRUE) && grepl(k, t_norm, fixed = TRUE)) mode_bonus <- mode_bonus + 0.01
  }
  jacc + hemi_bonus + mode_bonus
}

query_ubm_catalog <- function(max_id = 4500L, batch_size = 200L) {
  if (!requireNamespace("ieugwasr", quietly = TRUE)) stop("缺少R包: ieugwasr")
  ids <- paste0("ubm-b-", seq_len(max_id))
  out <- list()
  idx <- 1L
  for (st in seq(1L, length(ids), by = batch_size)) {
    ed <- min(st + batch_size - 1L, length(ids))
    chunk <- ids[st:ed]
    g <- tryCatch(ieugwasr::gwasinfo(chunk), error = function(e) NULL)
    if (!is.null(g) && nrow(g) > 0) {
      gd <- as.data.frame(g)
      keep <- intersect(c("id","trait","sample_size"), names(gd))
      out[[idx]] <- gd[, keep, drop = FALSE]
      idx <- idx + 1L
    }
  }
  if (length(out) == 0) return(data.frame(id = character(), trait = character(), sample_size = numeric(), stringsAsFactors = FALSE))
  res <- do.call(rbind, out)
  res <- res[!duplicated(res$id), , drop = FALSE]
  res
}

input_csv <- kv$input %||% file.path(getwd(), "Output_Tables", "Four_Model", "All_Independent_Effect_IDPs_for_Visualization_enhanced_parsed.csv")
output_csv <- kv$output %||% file.path(getwd(), "mr_brain_catalog_bulk_from_independent.csv")
mapping_csv <- kv$mapping %||% file.path(getwd(), "mr_brain_catalog_bulk_mapping_detail.csv")
max_id <- suppressWarnings(as.integer(kv$max_id %||% "4500"))
if (is.na(max_id) || max_id < 100) max_id <- 4500L

df <- read_csv_flex(input_csv)
if (!("variable_type" %in% names(df))) stop("输入文件缺少 variable_type 列")
df <- df[df$variable_type == "Brain_Imaging", , drop = FALSE]
if (nrow(df) == 0) stop("输入文件无 Brain_Imaging 行")

pick_col <- function(nms, cand) {
  hit <- cand[cand %in% nms]
  if (length(hit) == 0) NA_character_ else hit[[1]]
}

col_primary <- pick_col(names(df), c("idp_variable","variable_name","display_name"))
col_display <- pick_col(names(df), c("display_name","variable_name","idp_variable"))
if (is.na(col_primary) || is.na(col_display)) stop("无法识别脑表型名称列")

phenos <- unique(df[, c(col_primary, col_display), drop = FALSE])
names(phenos) <- c("phenotype_key", "phenotype_display")
phenos$phenotype_key <- as.character(phenos$phenotype_key)
phenos$phenotype_display <- as.character(phenos$phenotype_display)
phenos <- phenos[!is.na(phenos$phenotype_key) & nzchar(phenos$phenotype_key), , drop = FALSE]

ubm <- query_ubm_catalog(max_id = max_id, batch_size = 200L)
if (nrow(ubm) == 0) stop("未获取到 ubm-b 目录")
ubm <- ubm[!is.na(ubm$trait) & nzchar(ubm$trait), , drop = FALSE]
ubm <- ubm[grepl("idp|aparc|aseg|hipp|amyg|thalam|dmri|fmri|brain|cort", tolower(ubm$trait)), , drop = FALSE]

best_rows <- lapply(seq_len(nrow(phenos)), function(i) {
  pkey <- phenos$phenotype_key[i]
  pdisp <- phenos$phenotype_display[i]
  ptxt <- paste(pkey, pdisp)
  sc <- vapply(ubm$trait, function(t) score_pair(ptxt, t), numeric(1))
  ord <- order(sc, decreasing = TRUE)
  top <- ord[seq_len(min(5L, length(ord)))]
  data.frame(
    phenotype_key = pkey,
    phenotype_display = pdisp,
    matched_id = ubm$id[top],
    matched_trait = ubm$trait[top],
    score = sc[top],
    rank = seq_along(top),
    stringsAsFactors = FALSE
  )
})
map_df <- do.call(rbind, best_rows)
map_df <- map_df[order(map_df$phenotype_key, map_df$rank), , drop = FALSE]
utils::write.csv(map_df, mapping_csv, row.names = FALSE)

best <- map_df[map_df$rank == 1, , drop = FALSE]
best <- best[best$score >= 0.06, , drop = FALSE]
best <- best[!duplicated(best$phenotype_key), , drop = FALSE]

catalog <- data.frame(
  trait = best$phenotype_key,
  source = paste0("opengwas:", best$matched_id),
  cohort = "UKB_BIG40_UBM",
  stringsAsFactors = FALSE
)
utils::write.csv(catalog, output_csv, row.names = FALSE)

cat(sprintf("phenotypes_total=%d\n", nrow(phenos)))
cat(sprintf("catalog_matched=%d\n", nrow(catalog)))
cat(sprintf("catalog_path=%s\n", output_csv))
cat(sprintf("mapping_path=%s\n", mapping_csv))
