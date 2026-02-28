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

derive_group_label <- function(df, comparison) {
  if ("group_label" %in% names(df)) return(as.character(df$group_label))

  if ("Group" %in% names(df)) {
    g <- as.character(df$Group)
    g_lower <- tolower(g)
    out <- ifelse(g_lower %in% c("control", "controls", "ctrl"), "Control",
                  ifelse(grepl("ischemic", g_lower), "Ischemic",
                         ifelse(grepl("mi", g_lower), "MI",
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
  if (length(chosen) > max_n) chosen <- chosen[1:max_n]
  chosen
}

fit_and_save <- function(train_df, test_df, comparison, out_dir) {
  require_pkgs(c("randomForest", "ggplot2", "viridis", "dplyr"))

  train_df <- ensure_baseline_age(train_df)
  test_df <- ensure_baseline_age(test_df)

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

  train_vars <- c("baseline_age", idps)
  train_clean <- stats::na.omit(train_df[, train_vars, drop = FALSE])
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
  test_clean <- stats::na.omit(test_df[, test_vars, drop = FALSE])
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

  cov_df <- derive_covariates(test_df)
  cov_df$.__row_id__ <- test_df$.__row_id__
  pred_with_cov <- merge(test_clean, cov_df, by = ".__row_id__", all.x = TRUE, sort = FALSE)
  pred_with_cov <- pred_with_cov[is.finite(pred_with_cov$brain_age_gap) & is.finite(pred_with_cov$baseline_age), , drop = FALSE]

  reg_rows <- NULL
  if ("group_label" %in% names(pred_with_cov) && any(pred_with_cov$group_label %in% c("Control", "Ischemic", "Case", "MI", "Chronic"), na.rm = TRUE)) {
    df_reg <- pred_with_cov[pred_with_cov$group_label %in% c("Control", "Ischemic", "Case", "MI", "Chronic"), , drop = FALSE]
    df_reg$ischemic_binary <- ifelse(df_reg$group_label %in% c("Ischemic", "Case", "MI", "Chronic"), 1L, 0L)
    df_reg$sex_factor <- as.factor(df_reg$sex_factor)
    df_reg$ethnicity_factor <- as.factor(df_reg$ethnicity_factor)
    df_reg$imaging_center_factor <- as.factor(df_reg$imaging_center_factor)
    df_reg$diabetes_factor_baseline <- as.factor(df_reg$diabetes_factor_baseline)
    df_reg$smoking_factor_baseline <- as.factor(df_reg$smoking_factor_baseline)
    df_reg$alcohol_factor_baseline <- as.factor(df_reg$alcohol_factor_baseline)

    cov_candidates <- c(
      "baseline_age",
      "sex_factor","ethnicity_factor","imaging_center_factor",
      "bmi_baseline","systolic_bp_baseline","diastolic_bp_baseline",
      "diabetes_factor_baseline","smoking_factor_baseline","alcohol_factor_baseline"
    )
    cov_present <- cov_candidates[cov_candidates %in% names(df_reg)]
    rhs <- paste(c("ischemic_binary", cov_present), collapse = " + ")
    fml <- stats::as.formula(paste0("brain_age_gap ~ ", rhs))
    fit <- stats::lm(fml, data = df_reg)
    summ <- summary(fit)
    coefs <- as.data.frame(coef(summ))
    if ("ischemic_binary" %in% rownames(coefs)) {
      ci <- tryCatch(stats::confint(fit, "ischemic_binary", level = 0.95), error = function(e) matrix(c(NA_real_, NA_real_), nrow = 1))
      reg_rows <- data.frame(
        term = "ischemic_binary_adjusted",
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
      utils::write.csv(reg_rows, file.path(out_dir, "BAG_Independent_Effect_Adjusted.csv"), row.names = FALSE)
    }
  }

  p_pred <- ggplot2::ggplot(test_clean, ggplot2::aes(x = baseline_age, y = predicted_brain_age)) +
    ggplot2::geom_point(ggplot2::aes(color = brain_age_gap), alpha = 0.6, size = 2) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey40") +
    viridis::scale_color_viridis(option = "C", direction = 1, name = "BAG") +
    ggplot2::labs(title = "Predicted vs Chronological Age", x = "Chronological Age", y = "Predicted Brain Age") +
    theme_pub()
  save_plot(file.path(out_dir, "BAG_Predicted_vs_Actual"), p_pred)

  p_bag <- ggplot2::ggplot(test_clean, ggplot2::aes(x = brain_age_gap)) +
    ggplot2::geom_density(fill = "#3182bd", alpha = 0.3, color = "#08519c") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::labs(title = "Brain Age Gap (BAG) Distribution", x = "BAG", y = "Density") +
    theme_pub()
  save_plot(file.path(out_dir, "BAG_Distribution"), p_bag)

  if (any(!is.na(test_clean$group_label))) {
    p_group <- ggplot2::ggplot(test_clean, ggplot2::aes(x = group_label, y = brain_age_gap, fill = group_label)) +
      ggplot2::geom_violin(trim = FALSE, alpha = 0.3) +
      ggplot2::geom_boxplot(width = 0.2, outlier.shape = NA) +
      ggplot2::stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
      ggplot2::labs(title = "BAG by Group", x = "Group", y = "BAG") +
      theme_pub() +
      ggplot2::theme(legend.position = "none")
    save_plot(file.path(out_dir, "BAG_By_Group"), p_group)
  }

  invisible(list(model = rf, predictions = test_clean, idps = idps, adjusted_effect = reg_rows))
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
