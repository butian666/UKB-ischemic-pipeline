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

parse_source <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  if (grepl("^opengwas:", x)) return(sub("^opengwas:", "", x))
  x
}

clean_lab <- function(x, n = 44L) {
  x <- gsub("\\.+", " ", x)
  x <- gsub("\\s+", " ", trimws(x))
  ifelse(nchar(x) > n, paste0(substr(x, 1, n), "..."), x)
}

in_dir <- kv$in_dir %||% file.path(getwd(), "Output_Tables", "MR_Bidirectional_Bulk")
catalog_csv <- kv$catalog_csv %||% file.path(getwd(), "mr_brain_catalog_bulk_from_independent.csv")
out_dir <- kv$out_dir %||% file.path(in_dir, "Manuscript_and_Figures", "Circulation_Enhanced")
ihd_source <- kv$ihd_source %||% "opengwas:finn-b-I9_IHD"
top_n <- suppressWarnings(as.integer(kv$top_n %||% "12"))
if (is.na(top_n) || top_n <= 0) top_n <- 12L
retry_max <- suppressWarnings(as.integer(kv$retry_max %||% "4"))
if (is.na(retry_max) || retry_max < 1) retry_max <- 4L

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

main_csv <- file.path(in_dir, "MR_Manuscript_Main_Table.csv")
if (!file.exists(main_csv)) stop("Missing file: ", main_csv)
tbl <- read_csv_flex(main_csv)
catalog <- read_csv_flex(catalog_csv)
catalog$source_id <- vapply(catalog$source, parse_source, character(1))
trait_to_id <- stats::setNames(catalog$source_id, catalog$trait)
ihd_id <- parse_source(ihd_source)

tbl <- tbl[tbl$ivw_significance %in% c("nominal", "borderline", "significant_fdr"), , drop = FALSE]
tbl <- tbl[order(tbl$ivw_p, decreasing = FALSE), , drop = FALSE]
if (nrow(tbl) == 0) stop("No nominal/borderline/FDR-significant rows found.")
tbl <- tbl[seq_len(min(top_n, nrow(tbl))), , drop = FALSE]
tbl$pair_id <- sprintf("P%02d", seq_len(nrow(tbl)))
tbl$pair_label <- paste0(tbl$pair_id, ": ", clean_lab(tbl$exposure_trait), " -> ", clean_lab(tbl$outcome_trait))

get_ids_for_pair <- function(direction, exposure_trait, outcome_trait) {
  if (direction == "heart_to_brain") {
    exp_id <- ihd_id
    out_id <- trait_to_id[[outcome_trait]]
  } else {
    exp_id <- trait_to_id[[exposure_trait]]
    out_id <- ihd_id
  }
  list(exposure_id = exp_id, outcome_id = out_id)
}

if (!requireNamespace("TwoSampleMR", quietly = TRUE)) stop("TwoSampleMR not available.")
library(TwoSampleMR)
library(ggplot2)

enrich_harm_for_steiger <- function(harm) {
  for (side in c("exposure", "outcome")) {
    s_col <- paste0("samplesize.", side)
    if (!(s_col %in% names(harm))) harm[[s_col]] <- NA_real_
    r_col <- paste0("r.", side)
    if (!(r_col %in% names(harm))) harm[[r_col]] <- NA_real_
  }

  ids_need <- character()
  if ("id.exposure" %in% names(harm)) ids_need <- c(ids_need, unique(as.character(harm$id.exposure)))
  if ("id.outcome" %in% names(harm)) ids_need <- c(ids_need, unique(as.character(harm$id.outcome)))
  ids_need <- unique(ids_need[!is.na(ids_need) & nzchar(ids_need)])

  gi <- NULL
  if (length(ids_need) > 0 && requireNamespace("ieugwasr", quietly = TRUE)) {
    gi <- tryCatch(as.data.frame(ieugwasr::gwasinfo(ids_need)), error = function(e) NULL)
  }

  if (!is.null(gi) && nrow(gi) > 0) {
    sample_n <- gi$sample_size
    if ("ncase" %in% names(gi) && "ncontrol" %in% names(gi)) {
      sample_n <- ifelse(is.na(sample_n), as.numeric(gi$ncase) + as.numeric(gi$ncontrol), sample_n)
    }
    map_sample <- stats::setNames(sample_n, gi$id)
    map_ncase <- if ("ncase" %in% names(gi)) stats::setNames(gi$ncase, gi$id) else NULL
    map_ncontrol <- if ("ncontrol" %in% names(gi)) stats::setNames(gi$ncontrol, gi$id) else NULL
    map_unit <- if ("unit" %in% names(gi)) stats::setNames(gi$unit, gi$id) else NULL
    for (side in c("exposure", "outcome")) {
      id_col <- paste0("id.", side)
      s_col <- paste0("samplesize.", side)
      if (id_col %in% names(harm)) {
        idx <- is.na(harm[[s_col]]) | harm[[s_col]] <= 0
        h_ids <- as.character(harm[[id_col]])
        harm[[s_col]][idx] <- as.numeric(map_sample[h_ids[idx]])
        ncase_col <- paste0("ncase.", side)
        ncontrol_col <- paste0("ncontrol.", side)
        if (!(ncase_col %in% names(harm))) harm[[ncase_col]] <- NA_real_
        if (!(ncontrol_col %in% names(harm))) harm[[ncontrol_col]] <- NA_real_
        if (!is.null(map_ncase)) harm[[ncase_col]][is.na(harm[[ncase_col]])] <- as.numeric(map_ncase[h_ids[is.na(harm[[ncase_col]])]])
        if (!is.null(map_ncontrol)) harm[[ncontrol_col]][is.na(harm[[ncontrol_col]])] <- as.numeric(map_ncontrol[h_ids[is.na(harm[[ncontrol_col]])]])
        unit_col <- paste0("units.", side)
        if (!(unit_col %in% names(harm))) harm[[unit_col]] <- NA_character_
        if (!is.null(map_unit)) harm[[unit_col]][is.na(harm[[unit_col]]) | !nzchar(as.character(harm[[unit_col]]))] <- as.character(map_unit[h_ids[is.na(harm[[unit_col]]) | !nzchar(as.character(harm[[unit_col]]))]])
      }
    }
  }

  for (side in c("exposure", "outcome")) {
    r_col <- paste0("r.", side)
    b_col <- paste0("beta.", side)
    se_col <- paste0("se.", side)
    n_col <- paste0("samplesize.", side)
    unit_col <- paste0("units.", side)
    ncase_col <- paste0("ncase.", side)
    ncontrol_col <- paste0("ncontrol.", side)
    prev_col <- paste0("prevalence.", side)
    af_col <- if (side == "exposure") "eaf.exposure" else "eaf.outcome"

    is_bin <- FALSE
    if (unit_col %in% names(harm)) {
      is_bin <- any(grepl("log odds|binary|bin", tolower(as.character(harm[[unit_col]]))))
    }

    if (is_bin && all(c(b_col, af_col, ncase_col, ncontrol_col) %in% names(harm))) {
      prev <- if (prev_col %in% names(harm) && any(!is.na(harm[[prev_col]]))) harm[[prev_col]] else rep(0.1, nrow(harm))
      harm[[r_col]] <- tryCatch(
        get_r_from_lor(lor = harm[[b_col]], af = harm[[af_col]], ncase = harm[[ncase_col]], ncontrol = harm[[ncontrol_col]], prevalence = prev, model = "logit", correction = FALSE),
        error = function(e) harm[[r_col]]
      )
    }

    need_bsen <- all(c(b_col, se_col, n_col) %in% names(harm))
    if (need_bsen) {
      idx <- is.na(harm[[r_col]]) | !is.finite(harm[[r_col]])
      if (any(idx)) {
        harm[[r_col]][idx] <- tryCatch(
          get_r_from_bsen(b = harm[[b_col]][idx], se = harm[[se_col]][idx], n = harm[[n_col]][idx]),
          error = function(e) harm[[r_col]][idx]
        )
      }
    }
  }
  harm
}

steiger_rows <- list()
raps_rows <- list()
coloc_template <- list()
steiger_diag_rows <- list()
raps_diag_rows <- list()

process_one_pair <- function(r, retry_max) {
  ids <- get_ids_for_pair(r$direction, r$exposure_trait, r$outcome_trait)
  exp_id <- ids$exposure_id
  out_id <- ids$outcome_id
  if (is.null(exp_id) || is.na(exp_id) || !nzchar(exp_id) || is.null(out_id) || is.na(out_id) || !nzchar(out_id)) return(NULL)

  last_pack <- NULL
  for (attempt in seq_len(retry_max)) {
    harm <- NULL
    err <- NULL
    stg_err <- NULL
    raps_err <- NULL
    steiger_row <- NULL
    raps_row <- NULL

    tryCatch({
      exp_dat <- extract_instruments(exp_id, p1 = 5e-8, clump = TRUE, r2 = 0.001, kb = 10000)
      if (nrow(exp_dat) == 0) stop("no_instruments")
      out_dat <- extract_outcome_data(snps = exp_dat$SNP, outcomes = out_id)
      if (nrow(out_dat) == 0) stop("no_outcome_overlap")
      harm <- harmonise_data(exp_dat, out_dat)
      harm <- harm[harm$mr_keep, , drop = FALSE]
      if (nrow(harm) == 0) stop("no_harmonised_snps")
    }, error = function(e) {
      err <<- as.character(e$message)
    })

    coloc_row <- data.frame(
      pair_id = r$pair_id, direction = r$direction, exposure_trait = r$exposure_trait, outcome_trait = r$outcome_trait,
      exposure_gwas_id = exp_id, outcome_gwas_id = out_id, lead_snp = NA_character_, chr = NA_integer_, start_bp = NA_integer_, end_bp = NA_integer_,
      pp_h0 = NA_real_, pp_h1 = NA_real_, pp_h2 = NA_real_, pp_h3 = NA_real_, pp_h4 = NA_real_,
      status = ifelse(is.null(err), "ready_for_coloc", paste0("mr_data_issue: ", err)), stringsAsFactors = FALSE
    )

    steiger_diag <- data.frame(pair_id = r$pair_id, status = "failed_or_empty", message = "not_run", nsnp = ifelse(is.null(harm), NA, nrow(harm)), stringsAsFactors = FALSE)
    raps_diag <- data.frame(pair_id = r$pair_id, status = "failed_or_empty", message = "not_run", nsnp = ifelse(is.null(harm), NA, nrow(harm)), stringsAsFactors = FALSE)

    if (is.null(err) && !is.null(harm) && nrow(harm) >= 2) {
      harm <- enrich_harm_for_steiger(harm)
      stg <- tryCatch(directionality_test(harm), error = function(e) {
        stg_err <<- as.character(e$message)
        NULL
      })
      if (!is.null(stg) && nrow(stg) > 0) {
        v_or_na <- function(df, nms) {
          for (nm in nms) if (nm %in% names(df) && length(df[[nm]]) >= 1) return(as.numeric(df[[nm]][1]))
          NA_real_
        }
        v_logical_or_na <- function(df, nm) {
          if (!(nm %in% names(df)) || length(df[[nm]]) < 1) return(NA)
          as.logical(df[[nm]][1])
        }
        r2_exp_val <- v_or_na(stg, c("r2.exposure", "snp_r2.exposure"))
        r2_out_val <- v_or_na(stg, c("r2.outcome", "snp_r2.outcome"))
        if (is.na(r2_exp_val) && "r.exposure" %in% names(harm)) r2_exp_val <- sum((as.numeric(harm$r.exposure))^2, na.rm = TRUE)
        if (is.na(r2_out_val) && "r.outcome" %in% names(harm)) r2_out_val <- sum((as.numeric(harm$r.outcome))^2, na.rm = TRUE)
        steiger_row <- data.frame(
          pair_id = r$pair_id, direction = r$direction, exposure_trait = r$exposure_trait, outcome_trait = r$outcome_trait,
          nsnp = nrow(harm), r2_exposure = r2_exp_val, r2_outcome = r2_out_val, steiger_p = v_or_na(stg, c("steiger_pval")),
          correct_causal_direction = v_logical_or_na(stg, "correct_causal_direction"), stringsAsFactors = FALSE
        )
        steiger_diag <- data.frame(pair_id = r$pair_id, status = "ok", message = "", nsnp = nrow(harm), stringsAsFactors = FALSE)
      } else {
        steiger_diag <- data.frame(pair_id = r$pair_id, status = "failed_or_empty", message = ifelse(is.null(stg_err), "no_output", stg_err), nsnp = nrow(harm), stringsAsFactors = FALSE)
      }

      mr_methods <- tryCatch(mr(harm, method_list = c("mr_ivw", "mr_raps")), error = function(e) {
        raps_err <<- as.character(e$message)
        NULL
      })
      if (!is.null(mr_methods) && nrow(mr_methods) > 0) {
        ivw <- mr_methods[mr_methods$method == "Inverse variance weighted", , drop = FALSE]
        raps <- mr_methods[grepl("RAPS", mr_methods$method), , drop = FALSE]
        if (nrow(ivw) > 0 || nrow(raps) > 0) {
          raps_row <- data.frame(
            pair_id = r$pair_id, direction = r$direction, exposure_trait = r$exposure_trait, outcome_trait = r$outcome_trait, nsnp = nrow(harm),
            ivw_b = if (nrow(ivw) > 0) as.numeric(ivw$b[1]) else NA_real_, ivw_se = if (nrow(ivw) > 0) as.numeric(ivw$se[1]) else NA_real_, ivw_p = if (nrow(ivw) > 0) as.numeric(ivw$pval[1]) else NA_real_,
            raps_b = if (nrow(raps) > 0) as.numeric(raps$b[1]) else NA_real_, raps_se = if (nrow(raps) > 0) as.numeric(raps$se[1]) else NA_real_, raps_p = if (nrow(raps) > 0) as.numeric(raps$pval[1]) else NA_real_,
            stringsAsFactors = FALSE
          )
          raps_diag <- data.frame(pair_id = r$pair_id, status = "ok", message = "", nsnp = nrow(harm), stringsAsFactors = FALSE)
        } else {
          raps_diag <- data.frame(pair_id = r$pair_id, status = "failed_or_empty", message = "no_output", nsnp = nrow(harm), stringsAsFactors = FALSE)
        }
      } else {
        raps_diag <- data.frame(pair_id = r$pair_id, status = "failed_or_empty", message = ifelse(is.null(raps_err), "no_output", raps_err), nsnp = nrow(harm), stringsAsFactors = FALSE)
      }
    }

    steiger_diag$attempts_used <- attempt
    raps_diag$attempts_used <- attempt
    last_pack <- list(coloc_row = coloc_row, steiger_row = steiger_row, steiger_diag = steiger_diag, raps_row = raps_row, raps_diag = raps_diag)

    done <- identical(as.character(steiger_diag$status), "ok") && identical(as.character(raps_diag$status), "ok")
    msg_all <- paste(err %||% "", stg_err %||% "", raps_err %||% "", collapse = " | ")
    transient <- grepl("502|503|504|timeout|timed out|proxy|temporar|retry|rate", msg_all, ignore.case = TRUE)
    if (done) break
    if (attempt < retry_max && (transient || !done)) Sys.sleep(min(2^attempt, 8))
  }
  last_pack
}

for (i in seq_len(nrow(tbl))) {
  r <- tbl[i, , drop = FALSE]
  res <- process_one_pair(r, retry_max = retry_max)
  if (is.null(res)) next
  coloc_template[[length(coloc_template) + 1]] <- res$coloc_row
  if (!is.null(res$steiger_row)) steiger_rows[[length(steiger_rows) + 1]] <- res$steiger_row
  if (!is.null(res$raps_row)) raps_rows[[length(raps_rows) + 1]] <- res$raps_row
  steiger_diag_rows[[length(steiger_diag_rows) + 1]] <- res$steiger_diag
  raps_diag_rows[[length(raps_diag_rows) + 1]] <- res$raps_diag
}

steiger_df <- if (length(steiger_rows) > 0) do.call(rbind, steiger_rows) else data.frame(
  pair_id = character(), direction = character(), exposure_trait = character(), outcome_trait = character(),
  nsnp = numeric(), r2_exposure = numeric(), r2_outcome = numeric(), steiger_p = numeric(),
  correct_causal_direction = logical(), stringsAsFactors = FALSE
)
raps_df <- if (length(raps_rows) > 0) do.call(rbind, raps_rows) else data.frame(
  pair_id = character(), direction = character(), exposure_trait = character(), outcome_trait = character(),
  nsnp = numeric(), ivw_b = numeric(), ivw_se = numeric(), ivw_p = numeric(),
  raps_b = numeric(), raps_se = numeric(), raps_p = numeric(), stringsAsFactors = FALSE
)
coloc_df <- if (length(coloc_template) > 0) do.call(rbind, coloc_template) else data.frame()
steiger_diag_df <- if (length(steiger_diag_rows) > 0) do.call(rbind, steiger_diag_rows) else data.frame(pair_id = character(), status = character(), message = character(), nsnp = numeric(), stringsAsFactors = FALSE)
raps_diag_df <- if (length(raps_diag_rows) > 0) do.call(rbind, raps_diag_rows) else data.frame(pair_id = character(), status = character(), message = character(), nsnp = numeric(), stringsAsFactors = FALSE)

utils::write.csv(steiger_df, file.path(out_dir, "Steiger_Directionality_Table.csv"), row.names = FALSE)
utils::write.csv(raps_df, file.path(out_dir, "MR_RAPS_Comparison_Table.csv"), row.names = FALSE)
utils::write.csv(steiger_diag_df, file.path(out_dir, "Steiger_Directionality_Diagnostics.csv"), row.names = FALSE)
utils::write.csv(raps_diag_df, file.path(out_dir, "MR_RAPS_Diagnostics.csv"), row.names = FALSE)
utils::write.csv(coloc_df, file.path(out_dir, "Colocalization_Result_Template.csv"), row.names = FALSE)

coloc_input_template <- data.frame(
  pair_id = coloc_df$pair_id,
  exposure_sumstats_path = "",
  outcome_sumstats_path = "",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "ea",
  other_allele_col = "oa",
  pval_col = "pval",
  chr_col = "chr",
  pos_col = "pos",
  sample_size_col = "n",
  stringsAsFactors = FALSE
)
utils::write.csv(coloc_input_template, file.path(out_dir, "Colocalization_Input_Template.csv"), row.names = FALSE)

if (nrow(steiger_df) > 0) {
  steiger_df$direction_label <- ifelse(steiger_df$direction == "heart_to_brain", "Heart to Brain", "Brain to Heart")
  steiger_df$decision <- ifelse(steiger_df$correct_causal_direction %in% TRUE, "Supported", "Not supported")
  p_s1 <- ggplot(steiger_df, aes(x = decision, fill = direction_label)) +
    geom_bar(position = "dodge") +
    labs(x = "Steiger directionality decision", y = "Number of pairs", fill = "MR direction", title = "Steiger directionality summary") +
    theme_bw(base_size = 12)
  ggsave(file.path(out_dir, "Fig6A_Steiger_Direction_Bar.png"), p_s1, width = 8.5, height = 5.5, dpi = 320)

  p_s2 <- ggplot(steiger_df, aes(x = r2_exposure, y = r2_outcome, color = decision, shape = direction_label)) +
    geom_point(size = 2.8, alpha = 0.9) +
    geom_abline(intercept = 0, slope = 1, linetype = 2, color = "#666666") +
    labs(x = "R² explained in exposure", y = "R² explained in outcome", color = "Steiger decision", shape = "MR direction", title = "Steiger R² comparison") +
    theme_bw(base_size = 12)
  ggsave(file.path(out_dir, "Fig6B_Steiger_R2_Scatter.png"), p_s2, width = 7.5, height = 6, dpi = 320)
} else {
  p_s0 <- ggplot(data.frame(x = 1, y = 1, txt = "No valid Steiger output for current pairs.\nPlease provide per-trait sample size/trait metadata\nor run with precomputed r.exposure/r.outcome."), aes(x = x, y = y)) +
    geom_text(aes(label = txt), size = 4.8, lineheight = 1.15) +
    xlim(0, 2) + ylim(0, 2) +
    labs(title = "Steiger directionality: framework placeholder") +
    theme_void(base_size = 12)
  ggsave(file.path(out_dir, "Fig6A_Steiger_Direction_Bar.png"), p_s0, width = 8.5, height = 5.5, dpi = 320)
  ggsave(file.path(out_dir, "Fig6B_Steiger_R2_Scatter.png"), p_s0, width = 7.5, height = 6, dpi = 320)
}

if (nrow(raps_df) > 0) {
  raps_df$direction_label <- ifelse(raps_df$direction == "heart_to_brain", "Heart to Brain", "Brain to Heart")
  p_r1 <- ggplot(raps_df, aes(x = ivw_b, y = raps_b, color = direction_label)) +
    geom_point(size = 2.8, alpha = 0.9) +
    geom_abline(intercept = 0, slope = 1, linetype = 2, color = "#666666") +
    labs(x = "IVW beta", y = "MR-RAPS beta", color = "MR direction", title = "IVW versus MR-RAPS effect estimates") +
    theme_bw(base_size = 12)
  ggsave(file.path(out_dir, "Fig7A_MR_RAPS_vs_IVW_Scatter.png"), p_r1, width = 7.5, height = 6, dpi = 320)

  long_df <- rbind(
    data.frame(pair_id = raps_df$pair_id, method = "IVW", p_value = raps_df$ivw_p, direction_label = raps_df$direction_label, stringsAsFactors = FALSE),
    data.frame(pair_id = raps_df$pair_id, method = "MR-RAPS", p_value = raps_df$raps_p, direction_label = raps_df$direction_label, stringsAsFactors = FALSE)
  )
  long_df$neglog10p <- -log10(pmax(long_df$p_value, 1e-300))
  p_r2 <- ggplot(long_df, aes(x = method, y = neglog10p, fill = method)) +
    geom_boxplot(alpha = 0.75, outlier.shape = NA) +
    geom_jitter(width = 0.12, alpha = 0.7, size = 1.8) +
    facet_wrap(~direction_label) +
    labs(x = NULL, y = "-log10(P)", title = "Significance distribution: IVW vs MR-RAPS") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")
  ggsave(file.path(out_dir, "Fig7B_MR_RAPS_vs_IVW_Pvalue_Distribution.png"), p_r2, width = 10, height = 5.8, dpi = 320)
} else {
  p_r0 <- ggplot(data.frame(x = 1, y = 1, txt = "No MR-RAPS estimates were returned for current pairs.\nInstall/enable mr.raps backend and rerun for full panel."), aes(x = x, y = y)) +
    geom_text(aes(label = txt), size = 4.8, lineheight = 1.15) +
    xlim(0, 2) + ylim(0, 2) +
    labs(title = "MR-RAPS comparison: framework placeholder") +
    theme_void(base_size = 12)
  ggsave(file.path(out_dir, "Fig7A_MR_RAPS_vs_IVW_Scatter.png"), p_r0, width = 7.5, height = 6, dpi = 320)
  ggsave(file.path(out_dir, "Fig7B_MR_RAPS_vs_IVW_Pvalue_Distribution.png"), p_r0, width = 10, height = 5.8, dpi = 320)
}

if (nrow(coloc_df) > 0) {
  pp_cols <- c("pp_h0", "pp_h1", "pp_h2", "pp_h3", "pp_h4")
  for (cl in pp_cols) coloc_df[[cl]][is.na(coloc_df[[cl]])] <- 0
  coloc_top <- coloc_df[seq_len(min(8, nrow(coloc_df))), , drop = FALSE]
  tmp <- do.call(rbind, lapply(seq_len(nrow(coloc_top)), function(i) {
    data.frame(pair_id = coloc_top$pair_id[i], hypothesis = c("H0", "H1", "H2", "H3", "H4"), posterior = as.numeric(coloc_top[i, pp_cols]), stringsAsFactors = FALSE)
  }))
  p_c1 <- ggplot(tmp, aes(x = hypothesis, y = posterior, fill = hypothesis)) +
    geom_col() +
    facet_wrap(~pair_id, ncol = 4) +
    labs(x = "Colocalization hypothesis", y = "Posterior probability", title = "Colocalization posterior template (replace with computed PP values)") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none")
  ggsave(file.path(out_dir, "Fig8A_Colocalization_Posterior_Template.png"), p_c1, width = 12, height = 7.2, dpi = 320)

  placeholder <- data.frame(
    panel = c("Exposure locus", "Outcome locus", "Regional overlap"),
    x = c(1, 1, 1),
    y = c(1, 1, 1),
    text = c("Insert locus association plot\n(-log10 P vs position)", "Insert locus association plot\n(-log10 P vs position)", "Insert overlap/LD pattern panel\nand lead-SNP annotation"),
    stringsAsFactors = FALSE
  )
  p_c2 <- ggplot(placeholder, aes(x = x, y = y)) +
    geom_rect(aes(xmin = 0, xmax = 2, ymin = 0, ymax = 2), fill = "#f5f5f5", color = "#bdbdbd") +
    geom_text(aes(label = text), size = 4.2, lineheight = 1.1) +
    facet_wrap(~panel, ncol = 3) +
    xlim(0, 2) + ylim(0, 2) +
    labs(title = "Colocalization multi-panel layout template") +
    theme_void(base_size = 12)
  ggsave(file.path(out_dir, "Fig8B_Colocalization_MultiPanel_Layout_Template.png"), p_c2, width = 12, height = 4.6, dpi = 320)
}

cat(sprintf("out_dir=%s\n", out_dir))
cat(sprintf("steiger_rows=%d\n", nrow(steiger_df)))
cat(sprintf("raps_rows=%d\n", nrow(raps_df)))
cat(sprintf("coloc_template_rows=%d\n", nrow(coloc_df)))
