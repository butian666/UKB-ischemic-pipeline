args <- commandArgs(trailingOnly = TRUE)

get_script_path <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", ca, value = TRUE)
  if (length(hit) > 0) {
    return(normalizePath(sub("^--file=", "", hit[[1]]), winslash = "/", mustWork = FALSE))
  }
  frames <- sys.frames()
  ofiles <- vapply(frames, function(f) if (!is.null(f$ofile)) as.character(f$ofile) else NA_character_, character(1))
  ofiles <- ofiles[!is.na(ofiles) & nzchar(ofiles)]
  if (length(ofiles) > 0) {
    return(normalizePath(ofiles[[1]], winslash = "/", mustWork = FALSE))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) NA_character_)
    if (!is.na(p) && nzchar(p)) return(normalizePath(p, winslash = "/", mustWork = FALSE))
  }
  NA_character_
}

find_project_root <- function(start_dir) {
  start_dir <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  if (!dir.exists(start_dir)) return(getwd())
  cur <- start_dir
  for (i in 1:12) {
    candidate <- file.path(cur, "Output_Tables", "Cohort", "Matched")
    if (dir.exists(candidate)) return(cur)
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }
  getwd()
}

script_path <- get_script_path()
script_dir <- if (!is.na(script_path)) dirname(script_path) else getwd()
project_root <- find_project_root(script_dir)

default_cohort <- file.path(project_root, "Output_Tables", "Cohort", "Matched", "Final_Matched_Cohort_ischemic_vs_control.csv")
train_path <- if (length(args) >= 1 && nzchar(args[[1]])) args[[1]] else default_cohort
test_path  <- if (length(args) >= 2 && nzchar(args[[2]])) args[[2]] else default_cohort
comparison <- if (length(args) >= 3 && nzchar(args[[3]])) args[[3]] else "ischemic_vs_control"
out_dir    <- if (length(args) >= 4 && nzchar(args[[4]])) args[[4]] else file.path(project_root, "Output_Tables", "Brain_Age_Analysis_Independent")
n_max_arg  <- if (length(args) >= 5 && nzchar(args[[5]])) args[[5]] else NA_character_
ntree_arg  <- if (length(args) >= 6 && nzchar(args[[6]])) args[[6]] else NA_character_

parse_int_arg <- function(x) {
  if (is.na(x) || !nzchar(x)) return(NA_integer_)
  out <- suppressWarnings(as.integer(x))
  if (is.na(out)) NA_integer_ else out
}

n_max <- parse_int_arg(n_max_arg)
ntree <- parse_int_arg(ntree_arg)
if (is.na(ntree)) ntree <- 500L

log_message <- function(...) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s\n", ts, paste(..., collapse = " ")))
}

require_pkgs <- function(pkgs) {
  ok <- vapply(pkgs, function(p) requireNamespace(p, quietly = TRUE), logical(1))
  if (!all(ok)) stop("缺少R包: ", paste(pkgs[!ok], collapse = ", "))
  invisible(TRUE)
}

require_pkg_optional <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    log_message("未安装R包，跳过相关分析:", pkg)
    return(FALSE)
  }
  TRUE
}

read_any <- function(path, n_max = NA_integer_) {
  if (!file.exists(path)) stop("文件不存在: ", path, " (getwd=", getwd(), ")")
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") return(readRDS(path))
  if (ext %in% c("csv", "tsv")) {
    sep <- if (ext == "tsv") "\t" else ","
    nrows <- if (is.na(n_max)) -1L else n_max
    return(utils::read.csv(path, stringsAsFactors = FALSE, sep = sep, check.names = FALSE, nrows = nrows))
  }
  stop("不支持的文件格式: ", ext)
}

coalesce_cols <- function(df, cols) {
  cols <- cols[cols %in% names(df)]
  if (length(cols) == 0) return(rep(NA, nrow(df)))
  out <- df[[cols[[1]]]]
  if (length(cols) >= 2) {
    for (nm in cols[-1]) {
      idx <- is.na(out) | (is.character(out) & !nzchar(out))
      if (any(idx)) out[idx] <- df[[nm]][idx]
    }
  }
  out
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

ensure_baseline_age <- function(df) {
  if ("baseline_age" %in% names(df)) {
    df$baseline_age <- suppressWarnings(as.numeric(df$baseline_age))
    return(df)
  }
  candidates <- c(
    "Age.when.attended.assessment.centre...Instance.2",
    "Age.when.attended.assessment.centre...Instance.0",
    "Age.at.recruitment",
    "Age.at.recruitment...Instance.0",
    "Age.at.recruitment...Instance.2"
  )
  nm <- detect_col(names(df), candidates)
  if (is.na(nm)) stop("无法推导 baseline_age：未找到年龄列")
  df$baseline_age <- suppressWarnings(as.numeric(df[[nm]]))
  df
}

ensure_followup_age <- function(df) {
  if ("followup_age" %in% names(df)) {
    df$followup_age <- suppressWarnings(as.numeric(df$followup_age))
    return(df)
  }
  candidates <- c(
    "Age.when.attended.assessment.centre...Instance.3",
    "Age.at.recruitment...Instance.3",
    "Age.at.recruitment...Instance.2"
  )
  nm <- detect_col(names(df), candidates)
  if (is.na(nm)) {
    df$followup_age <- NA_real_
    return(df)
  }
  df$followup_age <- suppressWarnings(as.numeric(df[[nm]]))
  df
}

derive_group_label <- function(df, comparison) {
  if ("group_label" %in% names(df)) return(as.character(df$group_label))

  if ("Group" %in% names(df)) {
    g <- as.character(df$Group)
    g_lower <- tolower(trimws(g))
    out <- ifelse(grepl("control", g_lower), "Control",
                  ifelse(grepl("ischemic|ischaemic", g_lower), "Ischemic",
                         ifelse(grepl("\\bmi\\b|myocard", g_lower), "MI",
                                ifelse(grepl("chronic", g_lower), "Chronic", tools::toTitleCase(g)))))
    return(out)
  }

  if ("group" %in% names(df)) {
    g <- as.character(df$group)
    if (all(g %in% c("0", "1"), na.rm = TRUE)) {
      return(ifelse(g == "1", "Case", "Control"))
    }
  }

  subtype_col <- detect_col(names(df), c("disease_subtype", "Disease_subtype", "disease.subtype"))
  if (!is.na(subtype_col)) {
    st <- as.character(df[[subtype_col]])
    st_lower <- tolower(st)
    if (comparison == "ischemic_vs_control") {
      return(ifelse(st_lower %in% c("control"), "Control", "Ischemic"))
    }
    if (comparison == "mi_vs_control") {
      return(ifelse(st_lower %in% c("control"), "Control", ifelse(grepl("mi", st_lower), "MI", "Case")))
    }
    if (comparison == "chronic_vs_control") {
      return(ifelse(st_lower %in% c("control"), "Control", ifelse(grepl("chronic", st_lower), "Chronic", "Case")))
    }
    if (comparison == "mi_vs_chronic") {
      return(ifelse(grepl("chronic", st_lower), "Chronic", ifelse(grepl("mi", st_lower), "MI", "Other")))
    }
  }

  dx_col <- detect_col(names(df), c("diagnosis", "cohort", "case_control"))
  if (!is.na(dx_col)) {
    v <- as.character(df[[dx_col]])
    v_lower <- tolower(v)
    if (all(v %in% c("0", "1"), na.rm = TRUE)) return(ifelse(v == "1", "Case", "Control"))
    out <- ifelse(v_lower %in% c("control", "controls", "ctrl"), "Control", tools::toTitleCase(v))
    return(out)
  }

  rep(NA_character_, nrow(df))
}

attach_group_from_age_trend <- function(df, comparison, force = FALSE) {
  if (!("eid" %in% names(df))) return(df)
  if (!force && "group_label" %in% names(df) && any(!is.na(df$group_label))) return(df)
  df$eid <- suppressWarnings(as.numeric(df$eid))
  if (all(is.na(df$eid))) return(df)
  root_dir <- if (exists("project_root", inherits = TRUE) && is.character(project_root) && dir.exists(project_root)) {
    project_root
  } else {
    getwd()
  }
  raw_dir <- file.path(root_dir, "Age_Trend_Plots", "Age_Trend_Raw_Data")
  if (!dir.exists(raw_dir)) return(df)
  files <- list.files(raw_dir, pattern = paste0("_", comparison, "\\.csv$"), full.names = TRUE)
  if (length(files) == 0) return(df)
  best_map <- NULL
  best_cov <- 0
  n_try <- min(length(files), 50L)
  for (i in seq_len(n_try)) {
    map_raw <- tryCatch(utils::read.csv(files[[i]], stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
    if (is.null(map_raw) || nrow(map_raw) == 0) next
    if (!all(c("eid", "Group") %in% names(map_raw))) next
    map <- map_raw[, c("eid", "Group"), drop = FALSE]
    map$eid <- suppressWarnings(as.numeric(map$eid))
    map <- map[!is.na(map$eid), , drop = FALSE]
    map <- map[!duplicated(map$eid), , drop = FALSE]
    cov <- mean(df$eid %in% map$eid)
    if (is.finite(cov) && cov > best_cov) {
      best_cov <- cov
      best_map <- map
    }
    if (best_cov >= 0.95) break
  }
  if (is.null(best_map) || best_cov <= 0) return(df)
  merged <- merge(df, best_map, by = "eid", all.x = TRUE, sort = FALSE)
  merged$group_label <- derive_group_label(merged, comparison)
  merged$Group <- NULL
  merged
}

derive_covariates <- function(df) {
  sex_raw <- coalesce_cols(df, c("Sex", "sex", "sex_factor"))
  sex_factor <- ifelse(tolower(as.character(sex_raw)) %in% c("male", "m", "1"), "Male",
                       ifelse(tolower(as.character(sex_raw)) %in% c("female", "f", "0"), "Female", as.character(sex_raw)))
  sex_factor <- as.factor(sex_factor)

  eth_raw <- coalesce_cols(df, c("Ethnic.background...Instance.0", "Ethnic.background...Instance.1", "Ethnic.background...Instance.2", "Ethnic.background...Instance.3", "ethnicity_factor"))
  ethnicity_factor <- as.factor(as.character(eth_raw))

  center_raw <- coalesce_cols(df, c("UK.Biobank.assessment.centre...Instance.2", "UK.Biobank.assessment.centre...Instance.3", "imaging_center_factor"))
  imaging_center_factor <- as.factor(as.character(center_raw))

  bmi_raw <- coalesce_cols(df, c("Body.mass.index..BMI....Instance.2", "Body.mass.index..BMI....Instance.0", "bmi_baseline"))
  bmi_baseline <- suppressWarnings(as.numeric(bmi_raw))

  sbp_raw <- coalesce_cols(df, c("Systolic.blood.pressure..automated.reading...Instance.2...Array.0", "Systolic.blood.pressure..automated.reading...Instance.2...Array.1", "systolic_bp_baseline"))
  systolic_bp_baseline <- suppressWarnings(as.numeric(sbp_raw))

  dbp_raw <- coalesce_cols(df, c("Diastolic.blood.pressure..automated.reading...Instance.2...Array.0", "Diastolic.blood.pressure..automated.reading...Instance.2...Array.1", "diastolic_bp_baseline"))
  diastolic_bp_baseline <- suppressWarnings(as.numeric(dbp_raw))

  diab_raw <- coalesce_cols(df, c("Diabetes.diagnosed.by.doctor...Instance.0", "Diabetes.diagnosed.by.doctor...Instance.2", "diabetes_factor_baseline"))
  diabetes_factor_baseline <- as.factor(as.character(diab_raw))

  smoke_raw <- coalesce_cols(df, c("Smoking.status...Instance.0", "Smoking.status...Instance.2", "smoking_factor_baseline"))
  smoking_factor_baseline <- as.factor(as.character(smoke_raw))

  alcohol_raw <- coalesce_cols(df, c("Alcohol.intake.frequency....Instance.0", "Alcohol.intake.frequency....Instance.2", "alcohol_factor_baseline"))
  alcohol_factor_baseline <- as.factor(as.character(alcohol_raw))

  out <- data.frame(
    sex_factor = sex_factor,
    ethnicity_factor = ethnicity_factor,
    imaging_center_factor = imaging_center_factor,
    bmi_baseline = bmi_baseline,
    systolic_bp_baseline = systolic_bp_baseline,
    diastolic_bp_baseline = diastolic_bp_baseline,
    diabetes_factor_baseline = diabetes_factor_baseline,
    smoking_factor_baseline = smoking_factor_baseline,
    alcohol_factor_baseline = alcohol_factor_baseline,
    stringsAsFactors = FALSE
  )
  out
}

theme_pub <- function() ggplot2::theme_minimal(base_size = 12)

save_plot <- function(base_path, plot, width = 8, height = 6, dpi = 600) {
  dir.create(dirname(base_path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(paste0(base_path, ".png"), plot = plot, width = width, height = height, dpi = dpi, units = "in")
  ggplot2::ggsave(paste0(base_path, ".pdf"), plot = plot, width = width, height = height, dpi = dpi, units = "in")
}

plot_bag_diagnostics <- function(df, out_dir, tag = "", title_tag = NULL) {
  if (is.null(title_tag)) title_tag <- tag
  tag2 <- ifelse(nzchar(tag), tag, "")
  p_pred <- ggplot2::ggplot(df, ggplot2::aes(x = baseline_age, y = predicted_brain_age)) +
    ggplot2::geom_point(ggplot2::aes(color = brain_age_gap), alpha = 0.6, size = 2) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey40") +
    viridis::scale_color_viridis(option = "C", direction = 1, name = "BAG") +
    ggplot2::labs(title = paste0("Predicted vs Chronological Age", ifelse(nzchar(title_tag), paste0(" (", title_tag, ")"), "")),
                  x = "Chronological Age", y = "Predicted Brain Age") +
    theme_pub()
  save_plot(file.path(out_dir, paste0("BAG_Predicted_vs_Actual", tag2)), p_pred)

  p_bag <- ggplot2::ggplot(df, ggplot2::aes(x = brain_age_gap)) +
    ggplot2::geom_density(fill = "#3182bd", alpha = 0.3, color = "#08519c") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::labs(title = paste0("Brain Age Gap (BAG) Distribution", ifelse(nzchar(title_tag), paste0(" (", title_tag, ")"), "")),
                  x = "BAG", y = "Density") +
    theme_pub()
  save_plot(file.path(out_dir, paste0("BAG_Distribution", tag2)), p_bag)

  if ("group_label" %in% names(df) && any(!is.na(df$group_label))) {
    p_group <- ggplot2::ggplot(df, ggplot2::aes(x = group_label, y = brain_age_gap, fill = group_label)) +
      ggplot2::geom_violin(trim = FALSE, alpha = 0.3) +
      ggplot2::geom_boxplot(width = 0.2, outlier.shape = NA) +
      ggplot2::stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
      ggplot2::labs(title = paste0("BAG by Group", ifelse(nzchar(title_tag), paste0(" (", title_tag, ")"), "")),
                    x = "Group", y = "BAG") +
      theme_pub() +
      ggplot2::theme(legend.position = "none")
    save_plot(file.path(out_dir, paste0("BAG_By_Group", tag2)), p_group)

    p_dist_g <- ggplot2::ggplot(df, ggplot2::aes(x = brain_age_gap, fill = group_label, color = group_label)) +
      ggplot2::geom_density(alpha = 0.25) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
      ggplot2::labs(title = paste0("BAG Distribution by Group", ifelse(nzchar(title_tag), paste0(" (", title_tag, ")"), "")),
                    x = "BAG", y = "Density") +
      theme_pub()
    save_plot(file.path(out_dir, paste0("BAG_Distribution_By_Group", tag2)), p_dist_g)

    p_pred_g <- ggplot2::ggplot(df, ggplot2::aes(x = baseline_age, y = predicted_brain_age, color = group_label)) +
      ggplot2::geom_point(alpha = 0.6, size = 2) +
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey40") +
      ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 1.1) +
      ggplot2::labs(title = paste0("Predicted vs Chronological Age (by Group)", ifelse(nzchar(title_tag), paste0(" (", title_tag, ")"), "")),
                    x = "Chronological Age", y = "Predicted Brain Age") +
      theme_pub()
    save_plot(file.path(out_dir, paste0("BAG_Predicted_vs_Actual_By_Group", tag2)), p_pred_g)
  }

  p_bias <- ggplot2::ggplot(df, ggplot2::aes(x = baseline_age, y = brain_age_gap)) +
    ggplot2::geom_point(alpha = 0.5, size = 1.8) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 1.1, color = "#08519c") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::labs(title = paste0("Age-Bias Diagnostic (BAG vs Age)", ifelse(nzchar(title_tag), paste0(" (", title_tag, ")"), "")),
                  x = "Chronological Age", y = "BAG") +
    theme_pub()
  save_plot(file.path(out_dir, paste0("BAG_Age_Bias", tag2)), p_bias)
}

plot_longitudinal_diagnostics <- function(df_long, out_dir, tag = "", title_tag = NULL) {
  if (is.null(df_long) || nrow(df_long) == 0) return(invisible(NULL))
  if (is.null(title_tag)) title_tag <- tag
  tag2 <- ifelse(nzchar(tag), tag, "")
  if ("group_label" %in% names(df_long) && any(!is.na(df_long$group_label)) && "annual_delta_BAG" %in% names(df_long)) {
    d1 <- df_long[is.finite(df_long$annual_delta_BAG), , drop = FALSE]
    if (nrow(d1) > 0) {
      p1 <- ggplot2::ggplot(d1, ggplot2::aes(x = group_label, y = annual_delta_BAG, fill = group_label)) +
        ggplot2::geom_violin(trim = FALSE, alpha = 0.25) +
        ggplot2::geom_boxplot(width = 0.2, outlier.shape = NA) +
        ggplot2::stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
        ggplot2::labs(title = paste0("Annualized ΔBAG by Group", ifelse(nzchar(title_tag), paste0(" (", title_tag, ")"), "")),
                      x = "Group", y = "Annualized ΔBAG") +
        theme_pub() +
        ggplot2::theme(legend.position = "none")
      save_plot(file.path(out_dir, paste0("BAG_AnnualDelta_By_Group", tag2)), p1)
    }
  }
  if ("group_label" %in% names(df_long) && any(!is.na(df_long$group_label)) && "annual_delta_BAG_cal" %in% names(df_long)) {
    d2 <- df_long[is.finite(df_long$annual_delta_BAG_cal), , drop = FALSE]
    if (nrow(d2) > 0) {
      p2 <- ggplot2::ggplot(d2, ggplot2::aes(x = group_label, y = annual_delta_BAG_cal, fill = group_label)) +
        ggplot2::geom_violin(trim = FALSE, alpha = 0.25) +
        ggplot2::geom_boxplot(width = 0.2, outlier.shape = NA) +
        ggplot2::stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
        ggplot2::labs(title = paste0("Annualized ΔBAG (Calibrated) by Group", ifelse(nzchar(title_tag), paste0(" (", title_tag, ")"), "")),
                      x = "Group", y = "Annualized ΔBAG (Calibrated)") +
        theme_pub() +
        ggplot2::theme(legend.position = "none")
      save_plot(file.path(out_dir, paste0("BAG_AnnualDelta_By_Group_Calibrated", tag2)), p2)

      p3 <- ggplot2::ggplot(d2, ggplot2::aes(x = annual_delta_BAG_cal, fill = group_label, color = group_label)) +
        ggplot2::geom_density(alpha = 0.25) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
        ggplot2::labs(title = paste0("Annualized ΔBAG (Calibrated) Distribution", ifelse(nzchar(title_tag), paste0(" (", title_tag, ")"), "")),
                      x = "Annualized ΔBAG (Calibrated)", y = "Density") +
        theme_pub()
      save_plot(file.path(out_dir, paste0("BAG_AnnualDelta_Distribution_Calibrated", tag2)), p3)
    }
  }
  invisible(TRUE)
}

write_group_stats <- function(df, value_col, out_csv) {
  if (!("group_label" %in% names(df)) || !any(!is.na(df$group_label))) return(invisible(NULL))
  x <- df[[value_col]]
  stats_df <- aggregate(
    x,
    by = list(group_label = df$group_label),
    FUN = function(v) c(mean = mean(v, na.rm = TRUE), sd = stats::sd(v, na.rm = TRUE), n = sum(is.finite(v)))
  )
  stats_df <- data.frame(
    group_label = stats_df$group_label,
    mean = stats_df$x[, "mean"],
    sd = stats_df$x[, "sd"],
    n = stats_df$x[, "n"],
    stringsAsFactors = FALSE
  )
  utils::write.csv(stats_df, out_csv, row.names = FALSE)
  invisible(stats_df)
}

write_adjusted_effect <- function(df, outcome_col, out_csv, term_label) {
  if (!("group_label" %in% names(df)) || !any(df$group_label %in% c("Control", "Ischemic", "Case", "MI", "Chronic"), na.rm = TRUE)) return(invisible(NULL))
  df_reg <- df[df$group_label %in% c("Control", "Ischemic", "Case", "MI", "Chronic"), , drop = FALSE]
  if (!("baseline_age" %in% names(df_reg))) {
    if ("baseline_age.x" %in% names(df_reg)) df_reg$baseline_age <- df_reg[["baseline_age.x"]]
    if (!("baseline_age" %in% names(df_reg)) && "baseline_age.y" %in% names(df_reg)) df_reg$baseline_age <- df_reg[["baseline_age.y"]]
  }
  df_reg$ischemic_binary <- ifelse(df_reg$group_label %in% c("Ischemic", "Case", "MI", "Chronic"), 1L, 0L)
  for (nm in c("sex_factor", "ethnicity_factor", "imaging_center_factor", "diabetes_factor_baseline", "smoking_factor_baseline", "alcohol_factor_baseline")) {
    if (nm %in% names(df_reg)) df_reg[[nm]] <- as.factor(df_reg[[nm]])
  }
  cov_candidates <- c(
    "baseline_age",
    "sex_factor","ethnicity_factor","imaging_center_factor",
    "bmi_baseline","systolic_bp_baseline","diastolic_bp_baseline",
    "diabetes_factor_baseline","smoking_factor_baseline","alcohol_factor_baseline"
  )
  cov_present <- cov_candidates[cov_candidates %in% names(df_reg)]
  rhs <- paste(c("ischemic_binary", cov_present), collapse = " + ")
  fml <- stats::as.formula(paste0(outcome_col, " ~ ", rhs))
  fit <- stats::lm(fml, data = df_reg)
  summ <- summary(fit)
  coefs <- as.data.frame(coef(summ))
  if (!("ischemic_binary" %in% rownames(coefs))) return(invisible(NULL))
  ci <- tryCatch(stats::confint(fit, "ischemic_binary", level = 0.95), error = function(e) matrix(c(NA_real_, NA_real_), nrow = 1))
  reg_rows <- data.frame(
    term = term_label,
    estimate = coefs["ischemic_binary", "Estimate"],
    std_error = coefs["ischemic_binary", "Std. Error"],
    t_value = coefs["ischemic_binary", "t value"],
    p_value = coefs["ischemic_binary", "Pr(>|t|)"],
    conf_low = ci[1, 1],
    conf_high = ci[1, 2],
    n = nrow(df_reg),
    covariates_included = paste(cov_present, collapse = ";"),
    stringsAsFactors = FALSE
  )
  utils::write.csv(reg_rows, out_csv, row.names = FALSE)
  invisible(reg_rows)
}

calibrate_predicted_age_by_controls <- function(df, pred_col = "predicted_brain_age", age_col = "baseline_age", group_col = "group_label", control_label = "Control") {
  if (!(pred_col %in% names(df)) || !(age_col %in% names(df)) || !(group_col %in% names(df))) return(NULL)
  d0 <- df[is.finite(df[[pred_col]]) & is.finite(df[[age_col]]) & df[[group_col]] == control_label, , drop = FALSE]
  if (nrow(d0) < 30) return(NULL)
  fml <- stats::as.formula(paste0(pred_col, " ~ ", age_col))
  fit <- stats::lm(fml, data = d0)
  co <- stats::coef(fit)
  a <- unname(co[[1]])
  b <- unname(co[[2]])
  if (!is.finite(a) || !is.finite(b) || b == 0) return(NULL)
  out <- df
  out$predicted_brain_age <- (out[[pred_col]] - a) / b
  out$brain_age_gap <- out$predicted_brain_age - out[[age_col]]
  list(df = out, intercept = a, slope = b, n_ctrl = nrow(d0))
}

plot_bag_vs_age_by_group <- function(df, out_dir, tag = "", title_tag = NULL) {
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  if (!("baseline_age" %in% names(df)) || !("brain_age_gap" %in% names(df))) return(invisible(NULL))
  if (is.null(title_tag)) title_tag <- tag
  tag2 <- ifelse(nzchar(tag), tag, "")
  p <- ggplot2::ggplot(df, ggplot2::aes(x = baseline_age, y = brain_age_gap, color = group_label)) +
    ggplot2::geom_point(alpha = 0.55, size = 1.8) +
    ggplot2::geom_smooth(method = "loess", se = TRUE, linewidth = 1.1) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::labs(title = paste0("BAG vs Age (by Group)", ifelse(nzchar(title_tag), paste0(" (", title_tag, ")"), "")),
                  x = "Chronological Age", y = "BAG") +
    theme_pub()
  save_plot(file.path(out_dir, paste0("BAG_vs_Age_By_Group", tag2)), p)
  invisible(TRUE)
}

run_age_interaction_suite <- function(df, model_label, calibrated_flag, out_dir, tag = "") {
  if (is.null(df) || nrow(df) == 0) return(list(linear = NULL, spline = NULL, gam = NULL, strata = NULL))
  if (!("baseline_age" %in% names(df)) || !("brain_age_gap" %in% names(df)) || !("group_label" %in% names(df))) {
    return(list(linear = NULL, spline = NULL, gam = NULL, strata = NULL))
  }
  d <- df[df$group_label %in% c("Control", "Ischemic"), , drop = FALSE]
  d <- d[is.finite(d$baseline_age) & is.finite(d$brain_age_gap), , drop = FALSE]
  if (nrow(d) < 30) return(list(linear = NULL, spline = NULL, gam = NULL, strata = NULL))
  d$ischemic_binary <- ifelse(d$group_label == "Ischemic", 1L, 0L)
  d$group_factor <- factor(d$group_label, levels = c("Control", "Ischemic"))

  plot_bag_vs_age_by_group(d, out_dir = out_dir, tag = tag, title_tag = paste0(model_label, ifelse(calibrated_flag == 1, " (Calibrated)", "")))

  fit_lin <- stats::lm(brain_age_gap ~ baseline_age * ischemic_binary, data = d)
  s_lin <- summary(fit_lin)
  co <- as.data.frame(coef(s_lin))
  term_nm <- "baseline_age:ischemic_binary"
  lin_row <- data.frame(
    model = model_label,
    calibrated = calibrated_flag,
    n = nrow(d),
    interaction_term = term_nm,
    interaction_estimate = ifelse(term_nm %in% rownames(co), co[term_nm, "Estimate"], NA_real_),
    interaction_p_value = ifelse(term_nm %in% rownames(co), co[term_nm, "Pr(>|t|)"], NA_real_),
    aic = stats::AIC(fit_lin),
    adj_r2 = unname(s_lin$adj.r.squared),
    stringsAsFactors = FALSE
  )

  fit_ns0 <- stats::lm(brain_age_gap ~ ischemic_binary + splines::ns(baseline_age, df = 4), data = d)
  fit_ns1 <- stats::lm(brain_age_gap ~ ischemic_binary * splines::ns(baseline_age, df = 4), data = d)
  an_ns <- stats::anova(fit_ns0, fit_ns1)
  p_ns <- ifelse(nrow(an_ns) >= 2, an_ns$`Pr(>F)`[2], NA_real_)
  spline_row <- data.frame(
    model = model_label,
    calibrated = calibrated_flag,
    n = nrow(d),
    method = "cubic_spline_ns_df4",
    p_value_interaction = p_ns,
    aic_reduced = stats::AIC(fit_ns0),
    aic_full = stats::AIC(fit_ns1),
    stringsAsFactors = FALSE
  )

  gam_row <- NULL
  if (require_pkg_optional("mgcv")) {
    g0 <- mgcv::gam(brain_age_gap ~ group_factor + s(baseline_age, k = 6), data = d, method = "REML")
    g1 <- mgcv::gam(brain_age_gap ~ group_factor + s(baseline_age, k = 6) + s(baseline_age, by = group_factor, k = 6),
                    data = d, method = "REML")
    an_g <- mgcv::anova.gam(g0, g1, test = "F")
    p_g <- ifelse(nrow(an_g) >= 2, an_g$`Pr(>F)`[2], NA_real_)
    gam_row <- data.frame(
      model = model_label,
      calibrated = calibrated_flag,
      n = nrow(d),
      method = "gam_by_group",
      p_value_interaction = p_g,
      aic_reduced = stats::AIC(g0),
      aic_full = stats::AIC(g1),
      stringsAsFactors = FALSE
    )
  }

  q <- stats::quantile(d$baseline_age, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE, type = 7)
  q <- unique(as.numeric(q))
  if (length(q) < 3) {
    br <- sort(unique(as.numeric(pretty(d$baseline_age, n = 4))))
    if (length(br) >= 3) q <- br
  }
  strata_df <- NULL
  if (length(q) >= 3) {
    d$age_stratum <- cut(d$baseline_age, breaks = q, include.lowest = TRUE, right = TRUE)
    strata_levels <- levels(d$age_stratum)
    out_rows <- lapply(strata_levels, function(st) {
      ds <- d[d$age_stratum == st, , drop = FALSE]
      n_c <- sum(ds$group_label == "Control")
      n_i <- sum(ds$group_label == "Ischemic")
      m_c <- mean(ds$brain_age_gap[ds$group_label == "Control"], na.rm = TRUE)
      m_i <- mean(ds$brain_age_gap[ds$group_label == "Ischemic"], na.rm = TRUE)
      diff <- m_i - m_c
      p <- NA_real_
      if (n_c >= 10 && n_i >= 10) {
        p <- tryCatch(stats::t.test(brain_age_gap ~ group_factor, data = ds)$p.value, error = function(e) NA_real_)
      }
      data.frame(
        model = model_label,
        calibrated = calibrated_flag,
        age_stratum = as.character(st),
        n_control = n_c,
        n_ischemic = n_i,
        mean_bag_control = m_c,
        mean_bag_ischemic = m_i,
        diff_ischemic_minus_control = diff,
        p_value_ttest = p,
        stringsAsFactors = FALSE
      )
    })
    strata_df <- do.call(rbind, out_rows)
  }

  list(linear = lin_row, spline = spline_row, gam = gam_row, strata = strata_df)
}

select_idps <- function(df_train, df_test, max_n = 128) {
  drop_prefix <- c(
    "Age.", "Date.", "UK.", "Participant", "eid", "group", "group_label",
    "Sex", "Ethnic", "Qualifications", "Smoking", "Alcohol", "Body", "Diabetes", "Medication", "Systolic", "Diastolic",
    "weights", "subclass", "distance"
  )
  nms <- intersect(names(df_train), names(df_test))
  nms <- setdiff(nms, c("baseline_age", "predicted_brain_age", "brain_age_gap"))
  keep <- nms[!vapply(nms, function(x) any(startsWith(x, drop_prefix)), logical(1))]
  keep <- keep[vapply(keep, function(nm) is.numeric(df_train[[nm]]) || is.integer(df_train[[nm]]), logical(1))]
  keep <- keep[vapply(keep, function(nm) is.numeric(df_test[[nm]]) || is.integer(df_test[[nm]]), logical(1))]
  keep <- keep[!grepl("followup|incident|baseline|recruit|medication|diabetes|smoking|alcohol|townsend|education|waist|blood|pressure|centre|center|ethnic|sex", tolower(keep))]
  brain_pat <- paste(c(
    "mean\\.thickness", "grey\\.white\\.contrast", "^volume\\.", "\\bvolume\\b",
    "(weighted\\.mean|mean)\\.fa", "(weighted\\.mean|mean)\\.md",
    "(weighted\\.mean|mean)\\.(icvf|od|isovf)", "diffusiv",
    "rfmri", "connectivity", "amplitude",
    "area\\.", "\\barea\\b"
  ), collapse = "|")
  keep_brain <- keep[grepl(brain_pat, keep, ignore.case = TRUE)]
  chosen <- if (length(keep_brain) >= 5) keep_brain else keep
  if (any(grepl("Instance\\.2", chosen))) {
    chosen2 <- chosen[grepl("Instance\\.2", chosen)]
    if (length(chosen2) >= 5) chosen <- chosen2
  }
  if (length(chosen) > max_n) chosen <- chosen[1:max_n]
  chosen
}

median_impute_by_train <- function(train_df, test_df, cols) {
  cols <- cols[cols %in% names(train_df) & cols %in% names(test_df)]
  if (length(cols) == 0) return(list(train = train_df, test = test_df, cols = character(0), medians = numeric(0)))
  med <- vapply(cols, function(nm) {
    x <- suppressWarnings(as.numeric(train_df[[nm]]))
    if (all(is.na(x))) NA_real_ else stats::median(x, na.rm = TRUE)
  }, numeric(1))
  keep <- cols[!is.na(med)]
  med <- med[keep]
  if (length(keep) == 0) return(list(train = train_df, test = test_df, cols = character(0), medians = numeric(0)))
  for (nm in keep) {
    m <- med[[nm]]
    x_tr <- suppressWarnings(as.numeric(train_df[[nm]]))
    x_te <- suppressWarnings(as.numeric(test_df[[nm]]))
    if (anyNA(x_tr)) x_tr[is.na(x_tr)] <- m
    if (anyNA(x_te)) x_te[is.na(x_te)] <- m
    train_df[[nm]] <- x_tr
    test_df[[nm]] <- x_te
  }
  list(train = train_df, test = test_df, cols = keep, medians = med)
}

median_impute_vec <- function(x, med) {
  if (is.null(med) || length(med) == 0) return(x)
  if (anyNA(x)) x[is.na(x)] <- med
  x
}

fit_and_save <- function(train_df, test_df, comparison, out_dir) {
  require_pkgs(c("randomForest", "ggplot2", "viridis", "dplyr"))

  age_lin_rows <- list()
  spline_rows <- list()
  gam_rows <- list()
  strata_rows <- list()

  test_clean_rf_cal <- NULL
  test_clean_en_cal <- NULL

  train_df <- ensure_baseline_age(train_df)
  test_df <- ensure_baseline_age(test_df)
  test_df <- ensure_followup_age(test_df)

  test_df$group_label <- derive_group_label(test_df, comparison)
  if (!("group_label" %in% names(train_df))) train_df$group_label <- derive_group_label(train_df, comparison)

  if (!any(train_df$group_label == "Control", na.rm = TRUE)) train_df <- attach_group_from_age_trend(train_df, comparison, force = TRUE)
  if (!any(test_df$group_label == "Control", na.rm = TRUE)) test_df <- attach_group_from_age_trend(test_df, comparison, force = TRUE)

  if (any(!is.na(train_df$group_label)) && any(train_df$group_label == "Control", na.rm = TRUE)) {
    train_df <- train_df[train_df$group_label == "Control", , drop = FALSE]
    if (exists("n_max", inherits = TRUE) && is.numeric(n_max) && length(n_max) == 1 && !is.na(n_max) && n_max > 0 && nrow(train_df) > n_max) {
      set.seed(1)
      train_df <- train_df[sample.int(nrow(train_df), n_max), , drop = FALSE]
    }
  } else {
    stop("训练数据无法识别 Control 组：需要包含 `group`/`group_label`/`disease_subtype`/`diagnosis` 等列")
  }

  idps <- select_idps(train_df, test_df, max_n = 128)
  if (length(idps) < 5) stop("可用于建模的IDP列太少：", length(idps))
  log_message("用于脑龄预测的IDP列数量:", length(idps))
  imp <- median_impute_by_train(train_df, test_df, idps)
  train_df <- imp$train
  test_df <- imp$test
  idps <- imp$cols
  if (length(idps) < 5) stop("训练集可用中位数填补后可用IDP列太少：", length(idps))
  log_message("中位数填补后用于脑龄预测的IDP列数量:", length(idps))

  train_vars <- c("baseline_age", idps)
  train_clean <- train_df[, train_vars, drop = FALSE]
  train_clean <- train_clean[is.finite(train_clean$baseline_age), , drop = FALSE]
  mtry_val <- max(1L, floor(sqrt(length(idps))))
  log_message("训练Random Forest...")
  rf <- randomForest::randomForest(
    baseline_age ~ .,
    data = train_clean,
    ntree = ntree,
    mtry = mtry_val,
    importance = TRUE
  )

  test_df$.__row_id__ <- seq_len(nrow(test_df))
  test_vars <- intersect(c("baseline_age", idps, "group_label", "eid", ".__row_id__"), names(test_df))
  test_clean <- test_df[, test_vars, drop = FALSE]
  test_clean <- test_clean[is.finite(test_clean$baseline_age), , drop = FALSE]
  test_clean$predicted_brain_age <- stats::predict(rf, newdata = test_clean)
  test_clean$brain_age_gap <- test_clean$predicted_brain_age - test_clean$baseline_age

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(test_clean, file.path(out_dir, "BAG_Predictions.csv"), row.names = FALSE)
  utils::write.csv(
    data.frame(BAG_mean = mean(test_clean$brain_age_gap, na.rm = TRUE),
               BAG_sd = stats::sd(test_clean$brain_age_gap, na.rm = TRUE),
               n = nrow(test_clean)),
    file.path(out_dir, "BAG_Summary.csv"),
    row.names = FALSE
  )

  if (any(!is.na(test_clean$group_label))) {
    stats_df <- aggregate(
      test_clean$brain_age_gap,
      by = list(group_label = test_clean$group_label),
      FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = stats::sd(x, na.rm = TRUE), n = length(x))
    )
    stats_df <- data.frame(
      group_label = stats_df$group_label,
      mean_BAG = stats_df$x[, "mean"],
      sd_BAG = stats_df$x[, "sd"],
      n = stats_df$x[, "n"],
      stringsAsFactors = FALSE
    )
    utils::write.csv(stats_df, file.path(out_dir, "BAG_Group_Stats.csv"), row.names = FALSE)
  }

  plot_bag_diagnostics(test_clean, out_dir = out_dir, tag = "", title_tag = "Random Forest")

  ar0 <- run_age_interaction_suite(test_clean, model_label = "RandomForest", calibrated_flag = 0, out_dir = out_dir, tag = "_RF")
  if (!is.null(ar0$linear)) age_lin_rows[[length(age_lin_rows) + 1]] <- ar0$linear
  if (!is.null(ar0$spline)) spline_rows[[length(spline_rows) + 1]] <- ar0$spline
  if (!is.null(ar0$gam)) gam_rows[[length(gam_rows) + 1]] <- ar0$gam
  if (!is.null(ar0$strata)) strata_rows[[length(strata_rows) + 1]] <- ar0$strata

  cal_rf <- calibrate_predicted_age_by_controls(test_clean, pred_col = "predicted_brain_age", age_col = "baseline_age", group_col = "group_label", control_label = "Control")
  if (!is.null(cal_rf)) {
    test_clean_rf_cal <- cal_rf$df
    utils::write.csv(test_clean_rf_cal, file.path(out_dir, "BAG_Predictions_Calibrated.csv"), row.names = FALSE)
    utils::write.csv(
      data.frame(BAG_mean = mean(test_clean_rf_cal$brain_age_gap, na.rm = TRUE),
                 BAG_sd = stats::sd(test_clean_rf_cal$brain_age_gap, na.rm = TRUE),
                 n = nrow(test_clean_rf_cal),
                 calib_intercept = cal_rf$intercept,
                 calib_slope = cal_rf$slope,
                 calib_n_ctrl = cal_rf$n_ctrl),
      file.path(out_dir, "BAG_Summary_Calibrated.csv"),
      row.names = FALSE
    )
    write_group_stats(test_clean_rf_cal, "brain_age_gap", file.path(out_dir, "BAG_Group_Stats_Calibrated.csv"))
    plot_bag_diagnostics(test_clean_rf_cal, out_dir = out_dir, tag = "_Calibrated", title_tag = "Random Forest (Calibrated)")

    ar1 <- run_age_interaction_suite(test_clean_rf_cal, model_label = "RandomForest", calibrated_flag = 1, out_dir = out_dir, tag = "_RF_Calibrated")
    if (!is.null(ar1$linear)) age_lin_rows[[length(age_lin_rows) + 1]] <- ar1$linear
    if (!is.null(ar1$spline)) spline_rows[[length(spline_rows) + 1]] <- ar1$spline
    if (!is.null(ar1$gam)) gam_rows[[length(gam_rows) + 1]] <- ar1$gam
    if (!is.null(ar1$strata)) strata_rows[[length(strata_rows) + 1]] <- ar1$strata
  }

  cv_en <- NULL
  test_clean_en <- NULL
  if (require_pkg_optional("glmnet")) {
    log_message("训练Elastic Net...")
    set.seed(1)
    x_train <- as.matrix(train_clean[, idps, drop = FALSE])
    y_train <- suppressWarnings(as.numeric(train_clean$baseline_age))
    cv_en <- glmnet::cv.glmnet(
      x = x_train,
      y = y_train,
      alpha = 0.5,
      nfolds = 10,
      standardize = TRUE
    )
    x_test <- as.matrix(test_clean[, idps, drop = FALSE])
    pred_en <- as.numeric(stats::predict(cv_en, newx = x_test, s = "lambda.1se"))
    test_clean_en <- test_clean
    test_clean_en$predicted_brain_age <- pred_en
    test_clean_en$brain_age_gap <- test_clean_en$predicted_brain_age - test_clean_en$baseline_age

    utils::write.csv(test_clean_en, file.path(out_dir, "BAG_Predictions_ElasticNet.csv"), row.names = FALSE)
    utils::write.csv(
      data.frame(BAG_mean = mean(test_clean_en$brain_age_gap, na.rm = TRUE),
                 BAG_sd = stats::sd(test_clean_en$brain_age_gap, na.rm = TRUE),
                 n = nrow(test_clean_en)),
      file.path(out_dir, "BAG_Summary_ElasticNet.csv"),
      row.names = FALSE
    )
    write_group_stats(test_clean_en, "brain_age_gap", file.path(out_dir, "BAG_Group_Stats_ElasticNet.csv"))
    plot_bag_diagnostics(test_clean_en, out_dir = out_dir, tag = "_ElasticNet", title_tag = "Elastic Net")

    ar2 <- run_age_interaction_suite(test_clean_en, model_label = "ElasticNet", calibrated_flag = 0, out_dir = out_dir, tag = "_ElasticNet")
    if (!is.null(ar2$linear)) age_lin_rows[[length(age_lin_rows) + 1]] <- ar2$linear
    if (!is.null(ar2$spline)) spline_rows[[length(spline_rows) + 1]] <- ar2$spline
    if (!is.null(ar2$gam)) gam_rows[[length(gam_rows) + 1]] <- ar2$gam
    if (!is.null(ar2$strata)) strata_rows[[length(strata_rows) + 1]] <- ar2$strata

    cal_en <- calibrate_predicted_age_by_controls(test_clean_en, pred_col = "predicted_brain_age", age_col = "baseline_age", group_col = "group_label", control_label = "Control")
    if (!is.null(cal_en)) {
      test_clean_en_cal <- cal_en$df
      utils::write.csv(test_clean_en_cal, file.path(out_dir, "BAG_Predictions_Calibrated_ElasticNet.csv"), row.names = FALSE)
      utils::write.csv(
        data.frame(BAG_mean = mean(test_clean_en_cal$brain_age_gap, na.rm = TRUE),
                   BAG_sd = stats::sd(test_clean_en_cal$brain_age_gap, na.rm = TRUE),
                   n = nrow(test_clean_en_cal),
                   calib_intercept = cal_en$intercept,
                   calib_slope = cal_en$slope,
                   calib_n_ctrl = cal_en$n_ctrl),
        file.path(out_dir, "BAG_Summary_Calibrated_ElasticNet.csv"),
        row.names = FALSE
      )
      write_group_stats(test_clean_en_cal, "brain_age_gap", file.path(out_dir, "BAG_Group_Stats_Calibrated_ElasticNet.csv"))
      plot_bag_diagnostics(test_clean_en_cal, out_dir = out_dir, tag = "_Calibrated_ElasticNet", title_tag = "Elastic Net (Calibrated)")

      ar3 <- run_age_interaction_suite(test_clean_en_cal, model_label = "ElasticNet", calibrated_flag = 1, out_dir = out_dir, tag = "_ElasticNet_Calibrated")
      if (!is.null(ar3$linear)) age_lin_rows[[length(age_lin_rows) + 1]] <- ar3$linear
      if (!is.null(ar3$spline)) spline_rows[[length(spline_rows) + 1]] <- ar3$spline
      if (!is.null(ar3$gam)) gam_rows[[length(gam_rows) + 1]] <- ar3$gam
      if (!is.null(ar3$strata)) strata_rows[[length(strata_rows) + 1]] <- ar3$strata
    }

    cf <- as.matrix(stats::coef(cv_en, s = "lambda.1se"))
    coef_df <- data.frame(feature = rownames(cf), coef = as.numeric(cf[, 1]), stringsAsFactors = FALSE)
    coef_df <- coef_df[coef_df$feature != "(Intercept)", , drop = FALSE]
    coef_df$abs_coef <- abs(coef_df$coef)
    coef_df <- coef_df[order(coef_df$abs_coef, decreasing = TRUE), , drop = FALSE]
    top_k <- min(30L, nrow(coef_df))
    if (top_k >= 5L) {
      top_df <- coef_df[seq_len(top_k), , drop = FALSE]
      top_df$feature <- factor(top_df$feature, levels = rev(top_df$feature))
      p_coef <- ggplot2::ggplot(top_df, ggplot2::aes(x = feature, y = coef, fill = coef > 0)) +
        ggplot2::geom_col(show.legend = FALSE) +
        ggplot2::coord_flip() +
        ggplot2::labs(title = "Elastic Net Top Coefficients", x = NULL, y = "Coefficient") +
        theme_pub()
      save_plot(file.path(out_dir, "BAG_ElasticNet_Top_Coefficients"), p_coef, width = 10, height = 8)
    }
  }

  if ("followup_age" %in% names(test_df) && any(is.finite(test_df$followup_age))) {
    idps2 <- idps
    idps3 <- sub("Instance\\.2", "Instance.3", idps2)
    ok_pair <- (idps3 %in% names(test_df)) & (idps3 != idps2)
    idps2 <- idps2[ok_pair]
    idps3 <- idps3[ok_pair]
    if (length(idps2) >= 5) {
      dt_t3 <- test_df[, intersect(c("baseline_age", "followup_age", "group_label", "eid", ".__row_id__", idps3, "followup_years"), names(test_df)), drop = FALSE]
      dt_t3 <- dt_t3[is.finite(dt_t3$baseline_age) & is.finite(dt_t3$followup_age), , drop = FALSE]
      newdata3 <- data.frame(matrix(NA_real_, nrow = nrow(dt_t3), ncol = length(idps2)))
      names(newdata3) <- idps2
      for (j in seq_along(idps2)) {
        newdata3[[idps2[[j]]]] <- suppressWarnings(as.numeric(dt_t3[[idps3[[j]]]]))
      }
      med <- imp$medians[idps2]
      for (nm in idps2) {
        if (!is.na(med[[nm]])) newdata3[[nm]] <- median_impute_vec(newdata3[[nm]], med[[nm]])
      }
      pred3 <- stats::predict(rf, newdata = newdata3)
      bag3 <- pred3 - dt_t3$followup_age
      years <- if ("followup_years" %in% names(dt_t3)) suppressWarnings(as.numeric(dt_t3$followup_years)) else (dt_t3$followup_age - dt_t3$baseline_age)
      years[!is.finite(years) | years <= 0] <- NA_real_
      long_df <- data.frame(
        .__row_id__ = dt_t3$.__row_id__,
        eid = if ("eid" %in% names(dt_t3)) dt_t3$eid else NA,
        group_label = dt_t3$group_label,
        baseline_age = dt_t3$baseline_age,
        followup_age = dt_t3$followup_age,
        followup_years = years,
        predicted_brain_age_t3 = pred3,
        brain_age_gap_t3 = bag3,
        stringsAsFactors = FALSE
      )
      long_df <- merge(long_df, test_clean[, c(".__row_id__", "predicted_brain_age", "brain_age_gap"), drop = FALSE], by = ".__row_id__", all.x = TRUE, sort = FALSE)
      names(long_df)[names(long_df) == "predicted_brain_age"] <- "predicted_brain_age_t2"
      names(long_df)[names(long_df) == "brain_age_gap"] <- "brain_age_gap_t2"
      long_df$delta_BAG <- long_df$brain_age_gap_t3 - long_df$brain_age_gap_t2
      long_df$annual_delta_BAG <- long_df$delta_BAG / long_df$followup_years
      long_df <- long_df[is.finite(long_df$annual_delta_BAG), , drop = FALSE]
      utils::write.csv(long_df, file.path(out_dir, "BAG_Longitudinal.csv"), row.names = FALSE)

      if (!is.null(cv_en) && !is.null(test_clean_en)) {
        pred3_en <- as.numeric(stats::predict(cv_en, newx = as.matrix(newdata3), s = "lambda.1se"))
        bag3_en <- pred3_en - dt_t3$followup_age
        long_df_en <- data.frame(
          .__row_id__ = dt_t3$.__row_id__,
          eid = if ("eid" %in% names(dt_t3)) dt_t3$eid else NA,
          group_label = dt_t3$group_label,
          baseline_age = dt_t3$baseline_age,
          followup_age = dt_t3$followup_age,
          followup_years = years,
          predicted_brain_age_t3 = pred3_en,
          brain_age_gap_t3 = bag3_en,
          stringsAsFactors = FALSE
        )
        long_df_en <- merge(long_df_en, test_clean_en[, c(".__row_id__", "predicted_brain_age", "brain_age_gap"), drop = FALSE], by = ".__row_id__", all.x = TRUE, sort = FALSE)
        names(long_df_en)[names(long_df_en) == "predicted_brain_age"] <- "predicted_brain_age_t2"
        names(long_df_en)[names(long_df_en) == "brain_age_gap"] <- "brain_age_gap_t2"
        long_df_en$delta_BAG <- long_df_en$brain_age_gap_t3 - long_df_en$brain_age_gap_t2
        long_df_en$annual_delta_BAG <- long_df_en$delta_BAG / long_df_en$followup_years
        long_df_en <- long_df_en[is.finite(long_df_en$annual_delta_BAG), , drop = FALSE]
        utils::write.csv(long_df_en, file.path(out_dir, "BAG_Longitudinal_ElasticNet.csv"), row.names = FALSE)

        long_ctrl_en <- long_df_en[long_df_en$group_label == "Control" & is.finite(long_df_en$baseline_age) & is.finite(long_df_en$followup_age), , drop = FALSE]
        if (nrow(long_ctrl_en) >= 30) {
          fit_t2_en <- stats::lm(predicted_brain_age_t2 ~ baseline_age, data = long_ctrl_en)
          fit_t3_en <- stats::lm(predicted_brain_age_t3 ~ followup_age, data = long_ctrl_en)
          co2_en <- stats::coef(fit_t2_en)
          co3_en <- stats::coef(fit_t3_en)
          a2_en <- unname(co2_en[1]); b2_en <- unname(co2_en[2])
          a3_en <- unname(co3_en[1]); b3_en <- unname(co3_en[2])
          ok2_en <- is.finite(a2_en) && is.finite(b2_en) && b2_en != 0
          ok3_en <- is.finite(a3_en) && is.finite(b3_en) && b3_en != 0
          if (ok2_en && ok3_en) {
            long_df_en$predicted_brain_age_t2_cal <- (long_df_en$predicted_brain_age_t2 - a2_en) / b2_en
            long_df_en$predicted_brain_age_t3_cal <- (long_df_en$predicted_brain_age_t3 - a3_en) / b3_en
            long_df_en$brain_age_gap_t2_cal <- long_df_en$predicted_brain_age_t2_cal - long_df_en$baseline_age
            long_df_en$brain_age_gap_t3_cal <- long_df_en$predicted_brain_age_t3_cal - long_df_en$followup_age
            long_df_en$delta_BAG_cal <- long_df_en$brain_age_gap_t3_cal - long_df_en$brain_age_gap_t2_cal
            long_df_en$annual_delta_BAG_cal <- long_df_en$delta_BAG_cal / long_df_en$followup_years
            utils::write.csv(long_df_en, file.path(out_dir, "BAG_Longitudinal_Calibrated_ElasticNet.csv"), row.names = FALSE)
            stats_long_cal_en <- aggregate(
              long_df_en$annual_delta_BAG_cal,
              by = list(group_label = long_df_en$group_label),
              FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = stats::sd(x, na.rm = TRUE), n = sum(is.finite(x)))
            )
            stats_long_cal_en <- data.frame(
              group_label = stats_long_cal_en$group_label,
              mean_annual_delta_BAG_cal = stats_long_cal_en$x[, "mean"],
              sd_annual_delta_BAG_cal = stats_long_cal_en$x[, "sd"],
              n = stats_long_cal_en$x[, "n"],
              stringsAsFactors = FALSE
            )
            utils::write.csv(stats_long_cal_en, file.path(out_dir, "BAG_Longitudinal_Group_Stats_Calibrated_ElasticNet.csv"), row.names = FALSE)
          }
        }
        if (nrow(long_df_en) > 0 && any(!is.na(long_df_en$group_label))) {
          stats_long_en <- aggregate(
            long_df_en$annual_delta_BAG,
            by = list(group_label = long_df_en$group_label),
            FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = stats::sd(x, na.rm = TRUE), n = length(x))
          )
          stats_long_en <- data.frame(
            group_label = stats_long_en$group_label,
            mean_annual_delta_BAG = stats_long_en$x[, "mean"],
            sd_annual_delta_BAG = stats_long_en$x[, "sd"],
            n = stats_long_en$x[, "n"],
            stringsAsFactors = FALSE
          )
          utils::write.csv(stats_long_en, file.path(out_dir, "BAG_Longitudinal_Group_Stats_ElasticNet.csv"), row.names = FALSE)
        }
        plot_longitudinal_diagnostics(long_df_en, out_dir = out_dir, tag = "_ElasticNet", title_tag = "Elastic Net")
      }

      long_ctrl <- long_df[long_df$group_label == "Control" & is.finite(long_df$baseline_age) & is.finite(long_df$followup_age), , drop = FALSE]
      if (nrow(long_ctrl) >= 30) {
        fit_t2 <- stats::lm(predicted_brain_age_t2 ~ baseline_age, data = long_ctrl)
        fit_t3 <- stats::lm(predicted_brain_age_t3 ~ followup_age, data = long_ctrl)
        co2 <- stats::coef(fit_t2)
        co3 <- stats::coef(fit_t3)
        a2 <- unname(co2[1]); b2 <- unname(co2[2])
        a3 <- unname(co3[1]); b3 <- unname(co3[2])
        ok2 <- is.finite(a2) && is.finite(b2) && b2 != 0
        ok3 <- is.finite(a3) && is.finite(b3) && b3 != 0
        if (ok2 && ok3) {
          long_df$predicted_brain_age_t2_cal <- (long_df$predicted_brain_age_t2 - a2) / b2
          long_df$predicted_brain_age_t3_cal <- (long_df$predicted_brain_age_t3 - a3) / b3
          long_df$brain_age_gap_t2_cal <- long_df$predicted_brain_age_t2_cal - long_df$baseline_age
          long_df$brain_age_gap_t3_cal <- long_df$predicted_brain_age_t3_cal - long_df$followup_age
          long_df$delta_BAG_cal <- long_df$brain_age_gap_t3_cal - long_df$brain_age_gap_t2_cal
          long_df$annual_delta_BAG_cal <- long_df$delta_BAG_cal / long_df$followup_years
          utils::write.csv(long_df, file.path(out_dir, "BAG_Longitudinal_Calibrated.csv"), row.names = FALSE)
          stats_long_cal <- aggregate(
            long_df$annual_delta_BAG_cal,
            by = list(group_label = long_df$group_label),
            FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = stats::sd(x, na.rm = TRUE), n = sum(is.finite(x)))
          )
          stats_long_cal <- data.frame(
            group_label = stats_long_cal$group_label,
            mean_annual_delta_BAG_cal = stats_long_cal$x[, "mean"],
            sd_annual_delta_BAG_cal = stats_long_cal$x[, "sd"],
            n = stats_long_cal$x[, "n"],
            stringsAsFactors = FALSE
          )
          utils::write.csv(stats_long_cal, file.path(out_dir, "BAG_Longitudinal_Group_Stats_Calibrated.csv"), row.names = FALSE)
        }
      }
      if (nrow(long_df) > 0 && any(!is.na(long_df$group_label))) {
        stats_long <- aggregate(
          long_df$annual_delta_BAG,
          by = list(group_label = long_df$group_label),
          FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = stats::sd(x, na.rm = TRUE), n = length(x))
        )
        stats_long <- data.frame(
          group_label = stats_long$group_label,
          mean_annual_delta_BAG = stats_long$x[, "mean"],
          sd_annual_delta_BAG = stats_long$x[, "sd"],
          n = stats_long$x[, "n"],
          stringsAsFactors = FALSE
        )
        utils::write.csv(stats_long, file.path(out_dir, "BAG_Longitudinal_Group_Stats.csv"), row.names = FALSE)
      }
      plot_longitudinal_diagnostics(long_df, out_dir = out_dir, tag = "", title_tag = "Random Forest")
      cov_df_long <- derive_covariates(test_df)
      cov_df_long$.__row_id__ <- test_df$.__row_id__
      df_long_reg <- merge(long_df, cov_df_long, by = ".__row_id__", all.x = TRUE, sort = FALSE)
      df_long_reg <- df_long_reg[df_long_reg$group_label %in% c("Control", "Ischemic", "Case", "MI", "Chronic"), , drop = FALSE]
      if (nrow(df_long_reg) >= 20) {
        df_long_reg$ischemic_binary <- ifelse(df_long_reg$group_label %in% c("Ischemic", "Case", "MI", "Chronic"), 1L, 0L)
        cov_candidates2 <- c(
          "baseline_age",
          "sex_factor","ethnicity_factor","imaging_center_factor",
          "bmi_baseline","systolic_bp_baseline","diastolic_bp_baseline",
          "diabetes_factor_baseline","smoking_factor_baseline","alcohol_factor_baseline"
        )
        cov_present2 <- cov_candidates2[cov_candidates2 %in% names(df_long_reg)]
        rhs2 <- paste(c("ischemic_binary", cov_present2), collapse = " + ")
        fml2 <- stats::as.formula(paste0("annual_delta_BAG ~ ", rhs2))
        fit2 <- stats::lm(fml2, data = df_long_reg)
        summ2 <- summary(fit2)
        coefs2 <- as.data.frame(coef(summ2))
        if ("ischemic_binary" %in% rownames(coefs2)) {
          ci2 <- tryCatch(stats::confint(fit2, "ischemic_binary", level = 0.95), error = function(e) matrix(c(NA_real_, NA_real_), nrow = 1))
          reg2 <- data.frame(
            term = "ischemic_binary_on_annual_delta_BAG_adjusted",
            estimate = coefs2["ischemic_binary", "Estimate"],
            std_error = coefs2["ischemic_binary", "Std. Error"],
            t_value = coefs2["ischemic_binary", "t value"],
            p_value = coefs2["ischemic_binary", "Pr(>|t|)"],
            conf_low = ci2[1, 1],
            conf_high = ci2[1, 2],
            n = nrow(df_long_reg),
            covariates_included = paste(cov_present2, collapse = ";"),
            stringsAsFactors = FALSE
          )
          utils::write.csv(reg2, file.path(out_dir, "BAG_Longitudinal_Independent_Effect_Adjusted.csv"), row.names = FALSE)
        }
      }

      if ("annual_delta_BAG_cal" %in% names(df_long_reg) && sum(is.finite(df_long_reg$annual_delta_BAG_cal)) >= 20) {
        cov_candidates3 <- c(
          "baseline_age",
          "sex_factor","ethnicity_factor","imaging_center_factor",
          "bmi_baseline","systolic_bp_baseline","diastolic_bp_baseline",
          "diabetes_factor_baseline","smoking_factor_baseline","alcohol_factor_baseline"
        )
        cov_present3 <- cov_candidates3[cov_candidates3 %in% names(df_long_reg)]
        rhs3 <- paste(c("ischemic_binary", cov_present3), collapse = " + ")
        fml3 <- stats::as.formula(paste0("annual_delta_BAG_cal ~ ", rhs3))
        fit3 <- stats::lm(fml3, data = df_long_reg)
        summ3 <- summary(fit3)
        coefs3 <- as.data.frame(coef(summ3))
        if ("ischemic_binary" %in% rownames(coefs3)) {
          ci3 <- tryCatch(stats::confint(fit3, "ischemic_binary", level = 0.95), error = function(e) matrix(c(NA_real_, NA_real_), nrow = 1))
          reg3 <- data.frame(
            term = "ischemic_binary_on_annual_delta_BAG_cal_adjusted",
            estimate = coefs3["ischemic_binary", "Estimate"],
            std_error = coefs3["ischemic_binary", "Std. Error"],
            t_value = coefs3["ischemic_binary", "t value"],
            p_value = coefs3["ischemic_binary", "Pr(>|t|)"],
            conf_low = ci3[1, 1],
            conf_high = ci3[1, 2],
            n = sum(is.finite(df_long_reg$annual_delta_BAG_cal)),
            covariates_included = paste(cov_present3, collapse = ";"),
            stringsAsFactors = FALSE
          )
          utils::write.csv(reg3, file.path(out_dir, "BAG_Longitudinal_Independent_Effect_Adjusted_Calibrated.csv"), row.names = FALSE)
        }
      }

      if (!is.null(cv_en)) {
        long_en_path <- file.path(out_dir, "BAG_Longitudinal_ElasticNet.csv")
        if (file.exists(long_en_path)) {
          long_df_en2 <- utils::read.csv(long_en_path, stringsAsFactors = FALSE, check.names = FALSE)
          long_df_en2$group_label <- as.character(long_df_en2$group_label)
          df_long_reg_en <- merge(long_df_en2, cov_df_long, by = ".__row_id__", all.x = TRUE, sort = FALSE)
          write_adjusted_effect(df_long_reg_en, "annual_delta_BAG", file.path(out_dir, "BAG_Longitudinal_Independent_Effect_Adjusted_ElasticNet.csv"), "ischemic_binary_on_annual_delta_BAG_adjusted_ElasticNet")
        }

        long_en_cal_path <- file.path(out_dir, "BAG_Longitudinal_Calibrated_ElasticNet.csv")
        if (file.exists(long_en_cal_path)) {
          long_df_en3 <- utils::read.csv(long_en_cal_path, stringsAsFactors = FALSE, check.names = FALSE)
          long_df_en3$group_label <- as.character(long_df_en3$group_label)
          df_long_reg_en_cal <- merge(long_df_en3, cov_df_long, by = ".__row_id__", all.x = TRUE, sort = FALSE)
          write_adjusted_effect(df_long_reg_en_cal, "annual_delta_BAG_cal", file.path(out_dir, "BAG_Longitudinal_Independent_Effect_Adjusted_Calibrated_ElasticNet.csv"), "ischemic_binary_on_annual_delta_BAG_cal_adjusted_ElasticNet")
        }
      }
    }
  }

  cov_df <- derive_covariates(test_df)
  cov_df$.__row_id__ <- test_df$.__row_id__
  pred_with_cov <- merge(test_clean, cov_df, by = ".__row_id__", all.x = TRUE, sort = FALSE)
  pred_with_cov <- pred_with_cov[is.finite(pred_with_cov$brain_age_gap) & is.finite(pred_with_cov$baseline_age), , drop = FALSE]
  reg_rows <- write_adjusted_effect(pred_with_cov, "brain_age_gap", file.path(out_dir, "BAG_Independent_Effect_Adjusted.csv"), "ischemic_binary_adjusted")

  if (!is.null(test_clean_rf_cal)) {
    pred_with_cov_cal <- merge(test_clean_rf_cal, cov_df, by = ".__row_id__", all.x = TRUE, sort = FALSE)
    pred_with_cov_cal <- pred_with_cov_cal[is.finite(pred_with_cov_cal$brain_age_gap) & is.finite(pred_with_cov_cal$baseline_age), , drop = FALSE]
    write_adjusted_effect(pred_with_cov_cal, "brain_age_gap", file.path(out_dir, "BAG_Independent_Effect_Adjusted_Calibrated.csv"), "ischemic_binary_adjusted_Calibrated")
  }

  reg_rows_en <- NULL
  if (!is.null(test_clean_en)) {
    pred_with_cov_en <- merge(test_clean_en, cov_df, by = ".__row_id__", all.x = TRUE, sort = FALSE)
    pred_with_cov_en <- pred_with_cov_en[is.finite(pred_with_cov_en$brain_age_gap) & is.finite(pred_with_cov_en$baseline_age), , drop = FALSE]
    reg_rows_en <- write_adjusted_effect(pred_with_cov_en, "brain_age_gap", file.path(out_dir, "BAG_Independent_Effect_Adjusted_ElasticNet.csv"), "ischemic_binary_adjusted_ElasticNet")
  }

  if (!is.null(test_clean_en_cal)) {
    pred_with_cov_en_cal <- merge(test_clean_en_cal, cov_df, by = ".__row_id__", all.x = TRUE, sort = FALSE)
    pred_with_cov_en_cal <- pred_with_cov_en_cal[is.finite(pred_with_cov_en_cal$brain_age_gap) & is.finite(pred_with_cov_en_cal$baseline_age), , drop = FALSE]
    write_adjusted_effect(pred_with_cov_en_cal, "brain_age_gap", file.path(out_dir, "BAG_Independent_Effect_Adjusted_Calibrated_ElasticNet.csv"), "ischemic_binary_adjusted_Calibrated_ElasticNet")
  }

  age_lin_df <- if (length(age_lin_rows) > 0) do.call(rbind, age_lin_rows) else NULL
  if (!is.null(age_lin_df)) utils::write.csv(age_lin_df, file.path(out_dir, "BAG_Age_Interaction_Test.csv"), row.names = FALSE)

  spline_df <- if (length(spline_rows) > 0) do.call(rbind, spline_rows) else NULL
  if (!is.null(spline_df)) utils::write.csv(spline_df, file.path(out_dir, "BAG_Spline_Interaction_Test.csv"), row.names = FALSE)

  gam_df <- if (length(gam_rows) > 0) do.call(rbind, gam_rows) else NULL
  if (!is.null(gam_df)) utils::write.csv(gam_df, file.path(out_dir, "BAG_GAM_Interaction_Test.csv"), row.names = FALSE)

  strata_df_all <- if (length(strata_rows) > 0) do.call(rbind, strata_rows) else NULL
  if (!is.null(strata_df_all)) utils::write.csv(strata_df_all, file.path(out_dir, "BAG_AgeStrata_Estimates.csv"), row.names = FALSE)

  invisible(list(
    model_rf = rf,
    model_elastic_net = cv_en,
    predictions_rf = test_clean,
    predictions_rf_calibrated = test_clean_rf_cal,
    predictions_elastic_net = test_clean_en,
    predictions_elastic_net_calibrated = test_clean_en_cal,
    idps = idps,
    adjusted_effect_rf = reg_rows,
    adjusted_effect_elastic_net = reg_rows_en,
    age_interaction_linear = age_lin_df,
    age_interaction_spline = spline_df,
    age_interaction_gam = gam_df,
    age_strata_estimates = strata_df_all
  ))
}

log_message("Train:", train_path)
log_message("Test:", test_path)
log_message("Comparison:", comparison)
log_message("Out:", out_dir)
log_message("n_max:", ifelse(is.na(n_max), "NA", n_max))
log_message("ntree:", ntree)

train_df <- read_any(train_path, n_max = n_max)
if (!is.na(n_max) && is.numeric(n_max) && length(n_max) == 1 && n_max > 0) {
  tmp <- train_df
  tmp$group_label <- derive_group_label(tmp, comparison)
  if (!any(tmp$group_label == "Control", na.rm = TRUE)) tmp <- attach_group_from_age_trend(tmp, comparison, force = TRUE)
  if (!any(tmp$group_label == "Control", na.rm = TRUE)) {
    n_try <- n_max
    for (k in 1:6) {
      n_try <- n_try * 5L
      tmp2 <- read_any(train_path, n_max = n_try)
      tmp2$group_label <- derive_group_label(tmp2, comparison)
      if (!any(tmp2$group_label == "Control", na.rm = TRUE)) tmp2 <- attach_group_from_age_trend(tmp2, comparison, force = TRUE)
      if (any(tmp2$group_label == "Control", na.rm = TRUE)) {
        tmp <- tmp2
        break
      }
    }
  }
  train_df <- tmp
}
test_df <- read_any(test_path, n_max = n_max)

res <- fit_and_save(train_df, test_df, comparison = comparison, out_dir = out_dir)
log_message("完成")
