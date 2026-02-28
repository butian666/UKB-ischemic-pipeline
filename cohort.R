# =================== 统一合并所有数据并建立队列 ===================
library(tidyverse)
library(lubridate)
library(WeightIt)
library(survey)
library(tableone)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(vroom)

parse_ukb_date <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", "NULL", "null")] <- NA
  x <- ifelse(is.na(x), NA, sub("[;, ].*$", "", x))
  y <- suppressWarnings(lubridate::ymd(x))
  i <- is.na(y) & !is.na(x)
  y[i] <- suppressWarnings(lubridate::dmy(x[i]))
  i2 <- is.na(y) & !is.na(x)
  y[i2] <- suppressWarnings(lubridate::mdy(x[i2]))
  as.Date(y)
}

is_nonempty <- function(x) {
  if (is.character(x)) return(!is.na(x) & nzchar(trimws(x)))
  !is.na(x)
}

extract_first_factor <- function(df, pattern, exclude_pattern = "^Date\\.of\\.") {
  cols <- grep(pattern, names(df), ignore.case = TRUE, value = TRUE)
  cols <- cols[!grepl(exclude_pattern, cols)]
  if (length(cols) > 0) factor(df[[cols[1]]]) else factor(NA)
}

extract_baseline_single_or_mean <- function(df, base_pattern, array_pattern, fallback_candidates = character()) {
  cols_base <- grep(base_pattern, names(df), ignore.case = TRUE, value = TRUE)
  cols_base <- cols_base[!grepl("Array\\.[0-9]+$", cols_base)]
  if (length(cols_base) > 0) {
    return(suppressWarnings(as.numeric(df[[cols_base[1]]])))
  }
  cols_array <- grep(array_pattern, names(df), ignore.case = TRUE, value = TRUE)
  if (length(cols_array) > 0) {
    return(suppressWarnings(rowMeans(as.matrix(df[, cols_array, drop = FALSE]), na.rm = TRUE)))
  }
  for (fc in fallback_candidates) {
    if (fc %in% names(df)) return(suppressWarnings(as.numeric(df[[fc]])))
  }
  return(rep(NA_real_, nrow(df)))
}

unify_joined_columns <- function(df) {
  nms <- names(df)
  suf_pat <- "(\\.x|\\.y)+$"
  idx <- grepl(suf_pat, nms)
  if (!any(idx)) return(df)
  base <- sub(suf_pat, "", nms)
  ycnt <- integer(length(nms))
  ycnt[idx] <- stringr::str_count(nms[idx], "\\.y")
  to_remove <- character(0)
  for (bn in unique(base[idx])) {
    cols <- nms[base == bn & idx]
    ord <- order(-ycnt[match(cols, nms)])
    ordered <- cols[ord]
    classes <- vapply(ordered, function(nm) class(df[[nm]])[1], character(1))
    new <- df[[ordered[1]]]
    if (length(unique(classes)) == 1) {
      new <- Reduce(function(a, nm) dplyr::coalesce(a, df[[nm]]), ordered[-1], init = new)
    } else {
      for (nm in ordered[-1]) {
        na_idx <- is.na(new)
        if (any(na_idx)) new[na_idx] <- df[[nm]][na_idx]
      }
    }
    df[[bn]] <- new
    to_remove <- c(to_remove, ordered)
  }
  keep <- setdiff(nms, to_remove)
  df[keep]
}

fast_full_join_list <- function(dfs, by = "eid") {
  if (length(dfs) == 0) return(NULL)
  if (length(dfs) == 1) return(dfs[[1]])
  cur <- dfs
  while (length(cur) > 1) {
    nxt <- vector("list", ceiling(length(cur) / 2))
    k <- 1
    for (i in seq(1, length(cur), by = 2)) {
      if (i == length(cur)) {
        nxt[[k]] <- cur[[i]]
      } else {
        nxt[[k]] <- dplyr::full_join(cur[[i]], cur[[i + 1]], by = by, suffix = c("", ".y"))
      }
      k <- k + 1
    }
    cur <- nxt
  }
  unify_joined_columns(cur[[1]])
}


infer_idp_pairs <- function(df) {
  nms <- names(df)
  ex <- grepl("^(Date\\.of\\.|Age\\.when\\.attended\\.assessment\\.centre|UK\\.Biobank\\.assessment\\.centre)", nms)
  n2 <- nms[grepl("Instance\\.2", nms) & !ex]
  n3 <- nms[grepl("Instance\\.3", nms) & !ex]
  base2 <- sub("\\.\\.\\.Array\\.[0-9]+$", "", sub("\\.\\.\\.Instance\\.2.*$", "", n2))
  base3 <- sub("\\.\\.\\.Array\\.[0-9]+$", "", sub("\\.\\.\\.Instance\\.3.*$", "", n3))
  common <- intersect(base2, base3)
  if (length(common) == 0) return(data.frame(base = character(0), col2 = character(0), col3 = character(0)))
  col2 <- n2[match(common, base2)]
  col3 <- n3[match(common, base3)]
  data.frame(base = common, col2 = col2, col3 = col3, stringsAsFactors = FALSE)
}

idp_missingness_by_group <- function(df, group_var = "Group") {
  pairs <- infer_idp_pairs(df)
  if (nrow(pairs) == 0) return(tibble::tibble())
  groups <- as.character(df[[group_var]])
  ugs <- unique(groups)
  res <- lapply(ugs, function(g) {
    idx <- which(groups == g)
    subdf <- df[idx, , drop = FALSE]
    tibble::tibble(
      idp = pairs$base,
      group = g,
      n = length(idx),
      present2 = vapply(seq_len(nrow(pairs)), function(i) sum(!is.na(subdf[[pairs$col2[i]]])), integer(1)),
      present3 = vapply(seq_len(nrow(pairs)), function(i) sum(!is.na(subdf[[pairs$col3[i]]])), integer(1)),
      both_present = vapply(seq_len(nrow(pairs)), function(i) sum(!is.na(subdf[[pairs$col2[i]]]) & !is.na(subdf[[pairs$col3[i]]])), integer(1))
    ) %>% dplyr::mutate(
      rate_present2 = present2 / n,
      rate_present3 = present3 / n,
      rate_both_present = both_present / n,
      rate_missing3 = 1 - rate_present3
    )
  })
  dplyr::bind_rows(res)
}

cat("=== 重新开始：统一ID后直接合并所有数据 ===\n")
{
  script_path <- NA_character_
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) script_path <- sub("^--file=", "", file_arg[[1]])
  if (is.na(script_path) || !nzchar(script_path)) {
    script_path <- tryCatch(sys.frame(1)$ofile, error = function(e) NA_character_)
  }
  if (!is.na(script_path) && nzchar(script_path)) {
    script_dir <- dirname(normalizePath(script_path))
    if (dir.exists(script_dir)) setwd(script_dir)
  }
}
# 读取数据
cat("重新读取数据...\n")
data1 <- vroom::vroom("UKB.csv", delim = ",", na = "NA", col_select = -starts_with("..."))
data2 <- vroom::vroom("BrainMRI.csv", delim = ",", na = "NA", col_select = -starts_with("..."))

names(data1) <- make.names(names(data1))
names(data2) <- make.names(names(data2))

if ("Participant.ID" %in% names(data1)) data1 <- dplyr::rename(data1, eid = `Participant.ID`)
if ("Participant.ID" %in% names(data2)) data2 <- dplyr::rename(data2, eid = `Participant.ID`)

data1[["eid"]] <- as.character(data1[["eid"]])
data2[["eid"]] <- as.character(data2[["eid"]])

eligible_ukb <- data1 %>%
  dplyr::filter(
    is_nonempty(`Date.of.attending.assessment.centre...Instance.2`),
    is_nonempty(`Date.of.attending.assessment.centre...Instance.3`)
  )

data2_filtered <- data2 %>%
  dplyr::semi_join(eligible_ukb %>% dplyr::select(eid), by = "eid")

cat("\n按ID半连接后进行左连接...\n")
merged_data_complete <- dplyr::left_join(eligible_ukb, data2_filtered, by = "eid", suffix = c("", ".y"))
merged_data_complete <- unify_joined_columns(merged_data_complete)

cat("完整合并后数据维度:", dim(merged_data_complete), "\n")
cat("合并后样本量:", nrow(merged_data_complete), "\n")
cat("合并后变量数:", ncol(merged_data_complete), "\n")



# 兼容不同列名：疾病日期列可能在不同UKB导出中略有差异
pick_first_col <- function(df, patterns) {
  for (pat in patterns) {
    cols <- grep(pat, names(df), ignore.case = TRUE, value = TRUE)
    if (length(cols) > 0) return(cols[1])
  }
  NA_character_
}
col_date_mi <- pick_first_col(merged_data_complete, c("^Date\\.of\\.myocardial\\.infarction$", "myocardial\\.infarction"))
col_date_chronic <- pick_first_col(merged_data_complete, c("^Date\\.I25\\.first\\.reported", "chronic\\.ischaemic\\.heart\\.disease"))
col_date_angina <- pick_first_col(merged_data_complete, c("^Date\\.I20\\.first\\.reported", "angina\\.pectoris"))
col_date_stroke <- pick_first_col(merged_data_complete, c("^Date\\.of\\.stroke$", "\\bstroke\\b"))
col_date_dementia <- pick_first_col(merged_data_complete, c("^Date\\.of\\.all\\.cause\\.dementia\\.report$", "dementia"))
col_diabetes <- pick_first_col(merged_data_complete, c("^Diabetes\\.diagnosed\\.by\\.doctor.*Instance\\.0$", "diabetes.*diagnosed.*Instance\\.0", "diabetes.*Instance\\.0"))
x_date_mi <- if (!is.na(col_date_mi)) merged_data_complete[[col_date_mi]] else rep(NA, nrow(merged_data_complete))
x_date_chronic <- if (!is.na(col_date_chronic)) merged_data_complete[[col_date_chronic]] else rep(NA, nrow(merged_data_complete))
x_date_angina <- if (!is.na(col_date_angina)) merged_data_complete[[col_date_angina]] else rep(NA, nrow(merged_data_complete))
x_date_stroke <- if (!is.na(col_date_stroke)) merged_data_complete[[col_date_stroke]] else rep(NA, nrow(merged_data_complete))
x_date_dementia <- if (!is.na(col_date_dementia)) merged_data_complete[[col_date_dementia]] else rep(NA, nrow(merged_data_complete))
x_diabetes <- if (!is.na(col_diabetes)) merged_data_complete[[col_diabetes]] else rep(NA, nrow(merged_data_complete))

# 在完整的数据上构建队列
cat("\n=== 基于完整合并数据构建缺血性心脏病队列 ===\n")

ischemic_cohort_complete <- merged_data_complete %>%
  filter(
    is_nonempty(`Date.of.attending.assessment.centre...Instance.2`) &
      is_nonempty(`Date.of.attending.assessment.centre...Instance.3`)
  ) %>%
  mutate(
    date_baseline = parse_ukb_date(`Date.of.attending.assessment.centre...Instance.0`), 
    date_img1 = parse_ukb_date(`Date.of.attending.assessment.centre...Instance.2`),
    date_img2 = parse_ukb_date(`Date.of.attending.assessment.centre...Instance.3`),
    date_MI = parse_ukb_date(x_date_mi),
    date_chronic_ischemic = parse_ukb_date(x_date_chronic),
    date_angina = parse_ukb_date(x_date_angina),
    date_stroke = parse_ukb_date(x_date_stroke),
    date_dementia = parse_ukb_date(x_date_dementia),
    
    followup_years = as.numeric(date_img2 - date_img1) / 365.25,
    baseline_age = `Age.when.attended.assessment.centre...Instance.2`,
    followup_age = `Age.when.attended.assessment.centre...Instance.3`,
    age1 = baseline_age,
    age2 = followup_age
  ) %>%
  filter(!is.na(date_img1), !is.na(date_img2), followup_years > 0) %>%
  mutate(
    # 疾病逻辑
    MI_at_baseline = ifelse(!is.na(date_MI) & date_MI <= date_img1, 1, 0),
    MI_at_followup = ifelse(!is.na(date_MI) & date_MI <= date_img2, 1, 0),
    incident_MI = ifelse(!is.na(date_MI) & date_MI > date_img1 & date_MI <= date_img2, 1, 0),

    chronic_at_baseline = ifelse(!is.na(date_chronic_ischemic) & date_chronic_ischemic <= date_img1, 1, 0),
    chronic_at_followup = ifelse(!is.na(date_chronic_ischemic) & date_chronic_ischemic <= date_img2, 1, 0),
    incident_chronic_ischemic = ifelse(!is.na(date_chronic_ischemic) & date_chronic_ischemic > date_img1 & date_chronic_ischemic <= date_img2, 1, 0),
    # === 排除标准逻辑（统一以第一次影像扫描日期为“基线”）===
    angina_at_baseline = ifelse(!is.na(date_angina) & date_angina <= date_img1, 1, 0),
    stroke_at_baseline = ifelse(!is.na(date_stroke) & date_stroke <= date_img1, 1, 0),
    dementia_at_baseline = ifelse(!is.na(date_dementia) & date_dementia <= date_img1, 1, 0),
    
    
    incident_any_ischemic = ifelse(incident_MI == 1 | incident_chronic_ischemic == 1, 1, 0)
  ) %>%
  filter(MI_at_baseline == 0, 
         chronic_at_baseline == 0,
         angina_at_baseline == 0,       
         stroke_at_baseline == 0,       
         dementia_at_baseline == 0) %>%
  
  # 标准化协变量
  mutate(
    age_at_recruitment = {
      ar <- if ("Age.at.recruitment" %in% names(cur_data_all())) {
        suppressWarnings(as.numeric(cur_data_all()[["Age.at.recruitment"]]))
      } else {
        rep(NA_real_, n())
      }
      a0 <- if ("Age.when.attended.assessment.centre...Instance.0" %in% names(cur_data_all())) {
        suppressWarnings(as.numeric(cur_data_all()[["Age.when.attended.assessment.centre...Instance.0"]]))
      } else {
        rep(NA_real_, n())
      }
      dplyr::coalesce(
        ar,
        a0,
        baseline_age - as.numeric(date_img1 - date_baseline) / 365.25
      )
    },
    sex = dplyr::recode(Sex, Male = "Male", Female = "Female", .default = "Male"),
    sex_factor = factor(sex, levels = c("Female", "Male")),
    
    townsend_index = as.numeric(`Townsend.deprivation.index.at.recruitment`),
    
    education_raw = `Qualifications...Instance.0`,
    education = case_when(
      grepl("College|University|degree", education_raw, ignore.case = TRUE) ~ "Higher",
      grepl("A levels|O levels|GCSEs|NVQ|HND|HNC|Professional", education_raw, ignore.case = TRUE) ~ "Medium",
      TRUE ~ "Lower"
    ),
    education_factor = factor(education, levels = c("Lower", "Medium", "Higher")),
    
    ethnicity_raw = `Ethnic.background...Instance.0`,
    ethnicity = ifelse(grepl("White", ethnicity_raw, ignore.case = TRUE), "White", "Non-White"),
    ethnicity_factor = factor(ethnicity, levels = c("White", "Non-White")),
    
    bmi = as.numeric(`Body.mass.index..BMI....Instance.0`),
    systolic_bp = extract_baseline_single_or_mean(
      cur_data_all(),
      base_pattern = "^Systolic\\.blood\\.pressure.*Instance\\.0($|[^A-Za-z0-9])",
      array_pattern = "Systolic\\.blood\\.pressure.*Instance\\.0.*Array\\.[0-9]+$",
      fallback_candidates = c("Systolic.blood.pressure..automated.reading...Instance.0...Array.0")
    ),
    diastolic_bp = extract_baseline_single_or_mean(
      cur_data_all(),
      base_pattern = "^Diastolic\\.blood\\.pressure.*Instance\\.0($|[^A-Za-z0-9])",
      array_pattern = "Diastolic\\.blood\\.pressure.*Instance\\.0.*Array\\.[0-9]+$",
      fallback_candidates = c("Diastolic.blood.pressure..automated.reading...Instance.0...Array.0")
    ),
    
    diabetes = case_when(
      as.character(x_diabetes) %in% c("Yes", "YES", "yes", "1", "True", "TRUE", "true") ~ "Yes",
      TRUE ~ "No"
    ),
    diabetes_factor = factor(diabetes, levels = c("No", "Yes")),
    
    smoking = case_when(
      `Smoking.status...Instance.0` == "Never" ~ "Never",
      `Smoking.status...Instance.0` == "Current" ~ "Current", 
      TRUE ~ "Former"
    ),
    smoking_factor = factor(smoking, levels = c("Never", "Former", "Current")),
    
    alcohol = case_when(
      `Alcohol.intake.frequency....Instance.0` %in% c("Never", "Special occasions only") ~ "Never",
      `Alcohol.intake.frequency....Instance.0` %in% c("Once or twice a week", "Three or four times a week", "Daily or almost daily") ~ "Current",
      TRUE ~ "Former"
    ),
    alcohol_factor = factor(alcohol, levels = c("Never", "Former", "Current")),
    
    Group = factor(ifelse(incident_any_ischemic == 1, "Ischemic Heart Disease", "Control Group"), 
                   levels = c("Control Group", "Ischemic Heart Disease")),
    
    disease_subtype = case_when(
      incident_MI == 1 & incident_chronic_ischemic == 1 ~ "Both MI & Chronic",
      incident_MI == 1 & incident_chronic_ischemic == 0 ~ "MI Only", 
      incident_MI == 0 & incident_chronic_ischemic == 1 ~ "Chronic Only",
      TRUE ~ "Control"
    ),
    group = ifelse(incident_any_ischemic == 1, 1, 0),
    age_at_scan_1 = baseline_age,
    age_at_scan_2 = {
      yob <- if ("Year.of.birth" %in% names(cur_data_all())) {
        suppressWarnings(as.numeric(cur_data_all()[["Year.of.birth"]]))
      } else {
        rep(NA_real_, n())
      }
      dplyr::coalesce(
        followup_age,
        age_at_scan_1 + as.numeric(date_img2 - date_img1) / 365.25,
        age_at_recruitment + as.numeric(date_img2 - date_baseline) / 365.25,
        as.numeric(format(date_img2, "%Y")) - yob
      )
    },
    bmi_baseline = bmi,
    systolic_bp_baseline = systolic_bp,
    diastolic_bp_baseline = diastolic_bp,
    townsend_index_baseline = townsend_index,
    diabetes_factor_baseline = diabetes_factor,
    smoking_factor_baseline = smoking_factor,
    alcohol_factor_baseline = alcohol_factor,
    imaging_center_factor = extract_first_factor(cur_data_all(), "assessment\\.centre.*Instance\\.2$"),
    waist_circumference_baseline = extract_baseline_single_or_mean(
      cur_data_all(),
      base_pattern = "^Waist\\.circumference.*Instance\\.0($|[^A-Za-z0-9])",
      array_pattern = "Waist\\.circumference.*Instance\\.0.*Array\\.[0-9]+$",
      fallback_candidates = c("Waist.circumference...Instance.0")
    ),
    hip_circumference_baseline = extract_baseline_single_or_mean(
      cur_data_all(),
      base_pattern = "^Hip\\.circumference.*Instance\\.0($|[^A-Za-z0-9])",
      array_pattern = "Hip\\.circumference.*Instance\\.0.*Array\\.[0-9]+$",
      fallback_candidates = c("Hip.circumference...Instance.0")
    ),
    waist_hip_ratio_baseline = {
      w <- suppressWarnings(as.numeric(waist_circumference_baseline))
      h <- suppressWarnings(as.numeric(hip_circumference_baseline))
      ifelse(is.finite(w) & is.finite(h) & h != 0, w / h, NA_real_)
    },
    body_fat_percentage_baseline = extract_baseline_single_or_mean(
      cur_data_all(),
      base_pattern = "^Body\\.fat\\.percentage.*Instance\\.0($|[^A-Za-z0-9])",
      array_pattern = "Body\\.fat\\.percentage.*Instance\\.0.*Array\\.[0-9]+$",
      fallback_candidates = c("Body.fat.percentage...Instance.0", "Body.fat.percentage....Instance.0")
    ),
    whole_body_fat_mass_baseline = extract_baseline_single_or_mean(
      cur_data_all(),
      base_pattern = "^Whole\\.body\\.fat\\.mass.*Instance\\.0($|[^A-Za-z0-9])",
      array_pattern = "Whole\\.body\\.fat\\.mass.*Instance\\.0.*Array\\.[0-9]+$",
      fallback_candidates = c("Whole.body.fat.mass...Instance.0", "Whole.body.fat.mass....Instance.0")
    )
  ) 

fast_impute_median <- function(df, imp_vars) {
  dat <- df
  for (nm in imp_vars) {
    if (!nm %in% names(dat)) next
    v <- dat[[nm]]
    if (is.numeric(v) || is.integer(v)) {
      med <- suppressWarnings(stats::median(v, na.rm = TRUE))
      if (is.finite(med)) v[is.na(v)] <- med
      dat[[nm]] <- as.numeric(v)
    } else if (is.factor(v)) {
      tab <- table(v, useNA = "no")
      if (length(tab) > 0) {
        mode_val <- names(tab)[which.max(tab)]
        v[is.na(v)] <- mode_val
      }
      dat[[nm]] <- v
    } else if (is.character(v)) {
      tab <- table(v, useNA = "no")
      if (length(tab) > 0) {
        mode_val <- names(tab)[which.max(tab)]
        v[is.na(v)] <- mode_val
      }
      dat[[nm]] <- v
    }
  }
  dat
}

compute_pscore_glm <- function(df, covariates) {
  fml <- stats::as.formula(paste("group ~", paste(covariates, collapse = " + ")))
  fit <- stats::glm(fml, data = df, family = stats::binomial())
  stats::predict(fit, type = "response")
}

nearest_match <- function(df, pscore, ratio = 2, caliper = 0.2, replace = FALSE) {
  cases <- which(df$group == 1)
  ctrls <- which(df$group == 0)
  ps_c <- pscore[cases]
  ps_t <- pscore[ctrls]
  used <- rep(FALSE, length(ctrls))
  keep_idx <- integer(0)
  for (i in seq_along(cases)) {
    d <- abs(ps_t - ps_c[i])
    ord <- order(d)
    sel <- integer(0)
    for (j in ord) {
      if (!replace && used[j]) next
      if (d[j] <= caliper) {
        sel <- c(sel, ctrls[j])
        used[j] <- TRUE
        if (length(sel) >= ratio) break
      }
    }
    if (length(sel) > 0) keep_idx <- c(keep_idx, cases[i], sel)
  }
  df[sort(unique(keep_idx)), , drop = FALSE]
}
make_flowchart_tj <- function(initial_total, initial_cases, initial_controls, matched_cases, matched_controls, flow_counts = NULL) {
  box_w <- 28
  if (!is.null(flow_counts)) {
    main <- tibble::tibble(
      x = 0,
      y = seq(0, by = -18, length.out = 6),
      rect_h = 6,
      title = c("Initial merged", "With both imaging dates", "Valid follow-up", "Eligible", "PSM Matched", "Final Analysis"),
      info = c(
        paste0("N = ", flow_counts$initial_total),
        paste0("N = ", flow_counts$with_both_imgs),
        paste0("N = ", flow_counts$followup_valid),
        paste0("N = ", flow_counts$final_cohort_n),
        paste0("Cases: ", matched_cases, "   Controls: ", matched_controls),
        paste0("Final N = ", matched_cases + matched_controls)
      )
    )
    side <- tibble::tibble(
      x = 68,
      y = main$y[3],
      rect_h = 18,
      title = "Exclusions",
      info = paste0(
        "• MI at baseline: ", flow_counts$excluded_MI_baseline, "\n",
        "• Chronic ischemic at baseline: ", flow_counts$excluded_chronic_baseline, "\n",
        "• Angina at baseline: ", flow_counts$excluded_angina_baseline, "\n",
        "• Stroke at baseline: ", flow_counts$excluded_stroke_baseline, "\n",
        "• Dementia at baseline: ", flow_counts$excluded_dementia_baseline
      )
    )
    main$xmin <- main$x - box_w
    main$xmax <- main$x + box_w
    main$ymin <- main$y - main$rect_h
    main$ymax <- main$y + main$rect_h
    main$x_info <- main$xmin + 2
    side$xmin <- side$x - box_w
    side$xmax <- side$x + box_w
    side$ymin <- side$y - side$rect_h
    side$ymax <- side$y + side$rect_h
    side$x_info <- side$xmin + 2
    edges_main <- tibble::tibble(
      x = main$x[-nrow(main)],
      xend = main$x[-1],
      y = main$ymin[-nrow(main)] - 1.5,
      yend = main$ymax[-1] + 1.5
    )
    edge_side <- tibble::tibble(
      x = main$xmax[3] + 6,
      y = side$y,
      xend = side$xmin - 6,
      yend = side$y
    )
    ggplot2::ggplot() +
      ggplot2::geom_rect(data = main, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "#FFFFFF", color = "#000000", size = 1.2) +
      ggplot2::geom_rect(data = side, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "#FFFFFF", color = "#000000", size = 1.2) +
      ggplot2::geom_segment(data = edges_main, ggplot2::aes(x = x, y = y, xend = xend, yend = yend), colour = "#000000", size = 0.9, arrow = ggplot2::arrow(type = "closed", length = grid::unit(0.18, "inches"))) +
      ggplot2::geom_segment(data = edge_side, ggplot2::aes(x = x, y = y, xend = xend, yend = yend), colour = "#000000", size = 0.9, arrow = ggplot2::arrow(type = "closed", length = grid::unit(0.18, "inches"))) +
      ggplot2::geom_text(data = main, ggplot2::aes(x = x, y = y + 2.5, label = title), fontface = "bold", size = 5) +
      ggplot2::geom_text(data = main, ggplot2::aes(x = x_info, y = y - 2.2, label = info), hjust = 0, size = 4, lineheight = 1.12) +
      ggplot2::geom_text(data = side, ggplot2::aes(x = x, y = y + 3.0, label = title), fontface = "bold", size = 5) +
      ggplot2::geom_text(data = side, ggplot2::aes(x = x_info, y = y - 2.8, label = info), hjust = 0, size = 4, lineheight = 1.22) +
      ggplot2::theme_void() +
      ggplot2::theme(plot.margin = grid::unit(c(18, 18, 18, 18), "pt"))
  } else {
    nodes <- tibble::tibble(
      x = 0,
      y = c(0, -18, -36),
      rect_h = 6,
      title = c("Eligible", "PSM Matched", "Final Analysis"),
      info = c(
        paste0("Total: ", initial_total, "\nCases: ", initial_cases, "   Controls: ", initial_controls),
        paste0("Cases: ", matched_cases, "   Controls: ", matched_controls),
        paste0("Final N = ", matched_cases + matched_controls)
      )
    )
    nodes$xmin <- nodes$x - box_w
    nodes$xmax <- nodes$x + box_w
    nodes$ymin <- nodes$y - nodes$rect_h
    nodes$ymax <- nodes$y + nodes$rect_h
    nodes$x_info <- nodes$xmin + 2
    edges <- tibble::tibble(
      x = nodes$x[-nrow(nodes)],
      xend = nodes$x[-1],
      y = nodes$ymin[-1] + 2,
      yend = nodes$ymax[-nrow(nodes)] - 2
    )
    ggplot2::ggplot() +
      ggplot2::geom_rect(data = nodes, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "#FFFFFF", color = "#000000", size = 1.2) +
      ggplot2::geom_text(data = nodes, ggplot2::aes(x = x, y = y + 2.3, label = title), fontface = "bold", size = 5) +
      ggplot2::geom_text(data = nodes, ggplot2::aes(x = x_info, y = y - 2.0, label = info), hjust = 0, size = 4, lineheight = 1.12) +
      ggplot2::geom_segment(data = edges, ggplot2::aes(x = x, y = y, xend = xend, yend = yend), colour = "#000000", size = 0.9, arrow = ggplot2::arrow(type = "closed", length = grid::unit(0.18, "inches"))) +
      ggplot2::theme_void() +
      ggplot2::theme(plot.margin = grid::unit(c(18, 18, 18, 18), "pt"))
  }
}

run_psm_pipeline <- function(dat, covariates, ratio = 2, caliper = 0.2) {
  imp_vars <- unique(c("group", covariates))
  dx <- fast_impute_median(dat, imp_vars)
  ps <- compute_pscore_glm(dx, covariates)
  matched <- nearest_match(dx, ps, ratio = ratio, caliper = caliper, replace = FALSE)
  init_tot <- nrow(dx)
  init_cases <- sum(dx$group == 1)
  init_ctrls <- sum(dx$group == 0)
  m_cases <- sum(matched$group == 1)
  m_ctrls <- sum(matched$group == 0)
  psm_dir <- file.path(getwd(), "PSM_Output")
  if (!dir.exists(psm_dir)) dir.create(psm_dir)
  readr::write_csv(matched, file.path(psm_dir, "psm_matched.csv"))
  tbl <- tableone::CreateTableOne(vars = covariates, strata = "group", data = matched, test = FALSE)
  utils::write.csv(as.data.frame(print(tbl, printToggle = FALSE)), file.path(psm_dir, "PSM_Baseline_Table.csv"), row.names = FALSE)
  df_ps <- data.frame(pscore = ps, group = dx$group)
  p <- ggplot2::ggplot(df_ps, ggplot2::aes(x = pscore, fill = factor(group))) + ggplot2::geom_density(alpha = 0.5) + ggplot2::scale_fill_brewer(palette = "Set1") + ggplot2::theme_minimal()
  ggplot2::ggsave(filename = file.path(psm_dir, "PSM_Propensity_Distribution.png"), plot = p, width = 7, height = 4, dpi = 150)
  flow <- make_flowchart_tj(init_tot, init_cases, init_ctrls, m_cases, m_ctrls, NULL)
  ggplot2::ggsave(filename = file.path(psm_dir, "PSM_Flowchart.pdf"), plot = flow, width = 8, height = 3.5)
  list(matched = matched, propensity = ps)
}

 

build_comparison_df <- function(df, comparison) {
  case_levels <- switch(
    comparison,
    "ischemic_vs_control" = c("MI Only", "Chronic Only", "Both MI & Chronic"),
    "mi_vs_control" = c("MI Only", "Both MI & Chronic"),
    "chronic_vs_control" = c("Chronic Only"),
    "mi_vs_chronic" = c("MI Only")
  )
  control_levels <- switch(
    comparison,
    "mi_vs_chronic" = c("Chronic Only"),
    "ischemic_vs_control" = c("Control"),
    "mi_vs_control" = c("Control"),
    "chronic_vs_control" = c("Control")
  )
  df %>% dplyr::filter(disease_subtype %in% c(case_levels, control_levels)) %>% dplyr::mutate(group = ifelse(disease_subtype %in% case_levels, 1L, 0L))
}

compute_smd_df <- function(df, covariates) {
  res <- lapply(covariates, function(v) {
    x <- df[[v]]
    g <- df$group
    if (is.factor(x)) x <- as.numeric(x)
    xt <- x[g == 1]
    xc <- x[g == 0]
    m1 <- suppressWarnings(mean(xt, na.rm = TRUE))
    m0 <- suppressWarnings(mean(xc, na.rm = TRUE))
    s1 <- suppressWarnings(stats::var(xt, na.rm = TRUE))
    s0 <- suppressWarnings(stats::var(xc, na.rm = TRUE))
    sdpooled <- sqrt((s1 + s0) / 2)
    smd <- if (is.finite(sdpooled) && sdpooled > 0) (m1 - m0) / sdpooled else NA_real_
    data.frame(covariate = v, smd = smd)
  })
  do.call(rbind, res)
}

plot_love <- function(before_df, after_df, covariates, out_file) {
  smd_before <- compute_smd_df(before_df, covariates)
  smd_after <- compute_smd_df(after_df, covariates)
  smd_before$stage <- "Before"
  smd_after$stage <- "After"
  dd <- rbind(smd_before, smd_after)
  ord <- dd %>% dplyr::group_by(covariate) %>% dplyr::summarise(mx = max(abs(smd), na.rm = TRUE)) %>% dplyr::arrange(desc(mx))
  dd$covariate <- factor(dd$covariate, levels = ord$covariate)
  p <- ggplot2::ggplot(dd, ggplot2::aes(x = smd, y = covariate, color = stage)) + ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.5)) + ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "#666666") + ggplot2::theme_minimal()
  ggplot2::ggsave(out_file, plot = p, width = 7, height = 6, dpi = 150)
}

generate_baseline_outputs <- function(df, covariates, out_file_csv, out_file_png = NULL) {
  df1 <- df
  if (!is.factor(df1$group)) df1$group <- factor(df1$group, levels = c(0, 1), labels = c("Control", "Ischemic"))
  label_map <- c(
    age_at_scan_1 = "Age at scan 1",
    age_at_scan_2 = "Age at scan 2",
    age_at_recruitment = "Age at recruitment",
    sex_factor = "Sex (male/female)",
    ethnicity_factor = "Ethnicity (white/non-white)",
    education_factor = "Education (lower/medium/higher)",
    imaging_center_factor = "Location of first imaging center (Instance 2)",
    bmi_baseline = "BMI",
    systolic_bp_baseline = "Systolic blood pressure",
    diastolic_bp_baseline = "Diastolic blood pressure",
    townsend_index_baseline = "Townsend deprivation index",
    smoking_factor_baseline = "Smoking status",
    alcohol_factor_baseline = "Alcohol intake",
    diabetes_factor_baseline = "Diabetes status",
    waist_circumference_baseline = "Waist circumference",
    hip_circumference_baseline = "Hip circumference",
    waist_hip_ratio_baseline = "Waist-hip ratio",
    body_fat_percentage_baseline = "Body fat percentage",
    whole_body_fat_mass_baseline = "Whole body fat mass"
  )
  nm_present <- intersect(names(label_map), names(df1))
  for (nm in nm_present) {
    names(df1)[match(nm, names(df1))] <- label_map[[nm]]
  }
  cov_lab <- ifelse(covariates %in% names(label_map), label_map[covariates], covariates)
  facs <- cov_lab[vapply(cov_lab, function(nm) nm %in% names(df1) && is.factor(df1[[nm]]), logical(1))]
  bt <- tableone::CreateTableOne(vars = cov_lab, strata = "group", data = df1, factorVars = facs)
  mat <- print(bt, printToggle = FALSE, test = TRUE, smd = TRUE, showAllLevels = TRUE, noSpaces = TRUE, quote = FALSE)
  df_bt <- as.data.frame(mat)
  df_bt <- tibble::rownames_to_column(df_bt, var = "Variable")
  utils::write.csv(df_bt, out_file_csv, row.names = FALSE)
}

generate_sci_baseline_csv <- function(df, covariates, out_file_csv, group_col = "group", group_labels = c("Control", "Case")) {
  stopifnot(group_col %in% names(df))
  dat <- df
  dat <- dat[!is.na(dat[[group_col]]) & dat[[group_col]] %in% c(0, 1), , drop = FALSE]
  dat[[group_col]] <- as.integer(dat[[group_col]])
  g0 <- dat[dat[[group_col]] == 0, , drop = FALSE]
  g1 <- dat[dat[[group_col]] == 1, , drop = FALSE]
  covariates <- covariates[covariates %in% names(dat)]
  pv_cont <- function(y_name) {
    fml <- stats::as.formula(sprintf("%s ~ %s", y_name, group_col))
    tryCatch(
      stats::t.test(fml, data = dat)$p.value,
      error = function(e) {
        tryCatch(stats::wilcox.test(fml, data = dat)$p.value, error = function(e2) NA_real_)
      }
    )
  }
  pv_cat <- function(x) {
    tt <- table(dat[[group_col]], x, useNA = "no")
    if (nrow(tt) < 2 || ncol(tt) < 2) return(NA_real_)
    tryCatch(
      suppressWarnings(stats::chisq.test(tt)$p.value),
      error = function(e) suppressWarnings(stats::fisher.test(tt)$p.value)
    )
  }
  fmt_mean_sd_range <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- stats::na.omit(x)
    if (length(x) == 0) return("")
    paste0(sprintf("%.1f", mean(x)), "±", sprintf("%.1f", stats::sd(x)), " (", sprintf("%.1f", min(x)), "–", sprintf("%.1f", max(x)), ")")
  }
  fmt_count_pct <- function(x, level) {
    n <- sum(x == level, na.rm = TRUE)
    N <- sum(!is.na(x))
    pct <- ifelse(N > 0, 100 * n / N, NA_real_)
    paste0(n, " (", sprintf("%.1f", pct), "%)")
  }
  fmt_levels <- function(x) {
    xf <- as.factor(x)
    lv <- levels(xf)
    if (length(lv) == 0) return("")
    parts <- vapply(lv, function(l) fmt_count_pct(xf, l), character(1))
    paste(parts, collapse = "/")
  }
  label_map <- c(
    age_at_scan_1 = "Age at scan 1 (mean±s.d. [range])",
    age_at_scan_2 = "Age at scan 2 (mean±s.d. [range])",
    followup_years = "Years between scans 1 and 2 (mean±s.d. [range])",
    sex_factor = "Sex (male/female)",
    ethnicity_factor = "Ethnicity (white/non-white)",
    education_factor = "Education (lower/medium/higher)",
    imaging_center_factor = "Location of first imaging center (Instance 2)",
    bmi_baseline = "BMI (raw; mean±s.d. [range])",
    systolic_bp_baseline = "Systolic blood pressure (raw; mean±s.d. [range])",
    diastolic_bp_baseline = "Diastolic blood pressure (mmHg; mean±s.d. [range])",
    body_fat_percentage_baseline = "Body fat percentage (mean±s.d. [range])",
    whole_body_fat_mass_baseline = "Whole body fat mass (kg; mean±s.d. [range])",
    townsend_index_baseline = "Townsend deprivation index (raw; mean±s.d. [range])",
    diabetes_factor_baseline = "Diagnosed diabetes (factor yes/no)",
    smoking_factor_baseline = "Tobacco smoking",
    alcohol_factor_baseline = "Alcohol-intake frequency",
    waist_circumference_baseline = "Waist circumference (cm; mean±s.d. [range])",
    hip_circumference_baseline = "Hip circumference (cm; mean±s.d. [range])",
    waist_hip_ratio_baseline = "Waist/hip ratio (raw; mean±s.d. [range])"
  )
  rows <- list()
  rows[["Number of participants"]] <- c(nrow(g1), nrow(g0), "-")
  for (v in covariates) {
    lab <- if (!is.null(label_map[[v]])) label_map[[v]] else v
    x <- dat[[v]]
    if (is.numeric(x) || is.integer(x)) {
      p <- pv_cont(v)
      rows[[lab]] <- c(fmt_mean_sd_range(g1[[v]]), fmt_mean_sd_range(g0[[v]]), ifelse(is.na(p), "", sprintf("%.3f", p)))
    } else {
      xf <- as.factor(x)
      p <- pv_cat(xf)
      if (!(v %in% c("sex_factor", "ethnicity_factor", "education_factor", "imaging_center_factor", "diabetes_factor_baseline"))) {
        if (nlevels(xf) > 2) {
          lab <- paste0(lab, " (factor: ", paste(levels(xf), collapse = "/"), ")")
        } else if (nlevels(xf) == 2) {
          lab <- paste0(lab, " (factor ", paste(levels(xf), collapse = "/"), ")")
        }
      }
      rows[[lab]] <- c(fmt_levels(g1[[v]]), fmt_levels(g0[[v]]), ifelse(is.na(p), "", sprintf("%.3f", p)))
    }
  }
  out <- data.frame(
    Variable = names(rows),
    Case = vapply(rows, function(z) as.character(z[[1]]), character(1)),
    Control = vapply(rows, function(z) as.character(z[[2]]), character(1)),
    P_uncorr = vapply(rows, function(z) as.character(z[[3]]), character(1)),
    stringsAsFactors = FALSE
  )
  readr::write_csv(out, out_file_csv)
  invisible(out)
}

generate_topjournal_baseline_csv <- function(df, covariates, out_file_csv) {
  df1 <- df
  if (!is.factor(df1$group)) df1$group <- factor(df1$group, levels = c(0, 1), labels = c("Control", "Ischemic"))
  label_map <- c(
    age_at_scan_1 = "Age at scan 1",
    age_at_scan_2 = "Age at scan 2",
    followup_years = "Years between scans 1 and 2",
    sex_factor = "Sex (male/female)",
    ethnicity_factor = "Ethnicity (white/non-white)",
    education_factor = "Education (lower/medium/higher)",
    bmi_baseline = "BMI",
    systolic_bp_baseline = "Systolic blood pressure",
    diastolic_bp_baseline = "Diastolic blood pressure",
    diabetes_factor_baseline = "Diagnosed diabetes",
    waist_hip_ratio_baseline = "Waist/hip ratio",
    alcohol_factor_baseline = "Alcohol-intake frequency",
    smoking_factor_baseline = "Tobacco smoking",
    townsend_index_baseline = "Townsend deprivation index"
  )
  nm_present <- intersect(names(label_map), names(df1))
  for (nm in nm_present) {
    names(df1)[match(nm, names(df1))] <- label_map[[nm]]
  }
  cov_lab <- ifelse(covariates %in% names(label_map), label_map[covariates], covariates)
  facs <- cov_lab[vapply(cov_lab, function(nm) nm %in% names(df1) && is.factor(df1[[nm]]), logical(1))]
  bt <- tableone::CreateTableOne(vars = cov_lab, strata = "group", data = df1, factorVars = facs)
  mat <- print(bt, printToggle = FALSE, test = TRUE, smd = FALSE, showAllLevels = TRUE, noSpaces = TRUE, quote = FALSE, minMax = TRUE)
  df_bt <- as.data.frame(mat)
  df_bt <- tibble::rownames_to_column(df_bt, var = "Variable")
  if ("Ischemic" %in% names(df_bt) && "Control" %in% names(df_bt)) {
    df_bt <- df_bt[, c("Variable", "Ischemic", "Control", "p")]
  }
  names(df_bt)[names(df_bt) == "p"] <- "P_uncorr"
  top_row <- data.frame(Variable = "Number of participants", Ischemic = sum(df1$group == "Ischemic", na.rm = TRUE), Control = sum(df1$group == "Control", na.rm = TRUE), P_uncorr = "-")
  df_out <- rbind(top_row, df_bt)
  utils::write.csv(df_out, out_file_csv, row.names = FALSE)
}

run_psm_comparison <- function(df, comparison, psm_covariates, baseline_covariates, m = 5, ratio = 2, caliper = 0.2, flow_counts = NULL) {
  dat <- build_comparison_df(df, comparison)
  imp_vars <- unique(c("group", psm_covariates, baseline_covariates, if ("eid" %in% names(dat)) "eid" else NULL))
  baseline_tables <- list()
  matched_k1 <- NULL
  original_k1 <- NULL
  for (k in seq_len(m)) {
    dk <- fast_impute_median(dat, imp_vars)
    ps <- compute_pscore_glm(dk, psm_covariates)
    matched <- nearest_match(dk, ps, ratio = ratio, caliper = caliper, replace = FALSE)
    if (k == 1) {
      matched_k1 <- matched
      original_k1 <- dk
    }
    tbl1 <- tableone::CreateTableOne(vars = baseline_covariates, strata = "group", data = matched, test = FALSE)
    baseline_tables[[k]] <- as.data.frame(print(tbl1, printToggle = FALSE))
  }
  out_dir <- file.path(getwd(), "PSM_Output", comparison)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  baseline_merged <- if (length(baseline_tables) > 0) dplyr::bind_rows(baseline_tables, .id = "imp") else tibble::tibble()
  utils::write.csv(baseline_merged, file.path(out_dir, paste0("Baseline_Table_", comparison, ".csv")), row.names = FALSE)
  if (!is.null(matched_k1)) {
    utils::write.csv(matched_k1, file.path(out_dir, paste0("Matched_Example_", comparison, ".csv")), row.names = FALSE)
    # 基线三线表（全面协变量）
    generate_baseline_outputs(original_k1, baseline_covariates, file.path(out_dir, paste0("Baseline_Table_Before_", comparison, ".csv")), file.path(out_dir, paste0("Baseline_Table_Plot_Before_", comparison, ".png")))
    generate_baseline_outputs(matched_k1, baseline_covariates, file.path(out_dir, paste0("Baseline_Table_After_", comparison, ".csv")), file.path(out_dir, paste0("Baseline_Table_Plot_After_", comparison, ".png")))
    top_vars <- c("age_at_scan_1","age_at_scan_2","sex_factor","ethnicity_factor","followup_years","systolic_bp_baseline","diastolic_bp_baseline","diabetes_factor_baseline","waist_hip_ratio_baseline","bmi_baseline","alcohol_factor_baseline","smoking_factor_baseline","townsend_index_baseline")
    generate_topjournal_baseline_csv(matched_k1, top_vars, file.path(out_dir, "Baseline_Table_ischemic_vs_control.csv"))
    generate_topjournal_baseline_csv(original_k1, top_vars, file.path(out_dir, "Baseline_Table_ischemic_vs_control_Before.csv"))
    # PSM匹配变量三线表（仅五/六项匹配变量）
    generate_sci_baseline_csv(original_k1, psm_covariates, file.path(out_dir, paste0("Baseline_Table_PSM_Before_", comparison, ".csv")))
    generate_sci_baseline_csv(matched_k1, psm_covariates, file.path(out_dir, paste0("Baseline_Table_PSM_After_", comparison, ".csv")))
    ps0 <- compute_pscore_glm(original_k1, psm_covariates)
    df_ps <- data.frame(pscore = ps0, group = factor(original_k1$group, levels = c(0, 1), labels = c("control", "ischemic")))
    p <- ggplot2::ggplot(df_ps, ggplot2::aes(x = pscore, fill = group)) + ggplot2::geom_density(alpha = 0.5) + ggplot2::scale_fill_brewer(palette = "Set1") + ggplot2::labs(fill = NULL) + ggplot2::theme_minimal()
    ggplot2::ggsave(filename = file.path(out_dir, paste0("Propensity_Distribution_", comparison, ".png")), plot = p, width = 7, height = 4, dpi = 150)
  flow <- make_flowchart_tj(nrow(original_k1), sum(original_k1$group == 1), sum(original_k1$group == 0), sum(matched_k1$group == 1), sum(matched_k1$group == 0), flow_counts)
    ggplot2::ggsave(filename = file.path(out_dir, paste0("PSM_Flowchart_", comparison, ".pdf")), plot = flow, width = 7.8, height = 12)
    # Love Plot 仅针对 PSM 匹配变量
    plot_love(original_k1[, unique(c(psm_covariates, "group")), drop = FALSE], matched_k1[, unique(c(psm_covariates, "group")), drop = FALSE], psm_covariates, file.path(out_dir, paste0("Love_Plot_", comparison, ".png")))
    comps <- c("mi_vs_control", "chronic_vs_control", "mi_vs_chronic")
    for (cmp in comps) {
      outc <- file.path(getwd(), "PSM_Output", cmp)
      if (!dir.exists(outc)) dir.create(outc, recursive = TRUE)
      bf <- build_comparison_df(original_k1, cmp)
      af <- build_comparison_df(matched_k1, cmp)
      # 基线三线表（全面协变量）
      generate_baseline_outputs(bf, baseline_covariates, file.path(outc, paste0("Baseline_Table_Before_", cmp, ".csv")))
      generate_baseline_outputs(af, baseline_covariates, file.path(outc, paste0("Baseline_Table_After_", cmp, ".csv")))
      # PSM匹配变量三线表（仅CSV，不绘图）
      generate_sci_baseline_csv(bf, psm_covariates, file.path(outc, paste0("Baseline_Table_PSM_Before_", cmp, ".csv")))
      generate_sci_baseline_csv(af, psm_covariates, file.path(outc, paste0("Baseline_Table_PSM_After_", cmp, ".csv")))
    }
  }
  list(matched_example = matched_k1)
}

cat("最终队列人数:", nrow(ischemic_cohort_complete), "\n")
cat("- 慢性缺血性日期非NA:", sum(!is.na(ischemic_cohort_complete$date_chronic_ischemic)), "\n")
cat("- 慢性缺血性日期 ≤ 扫描1:", sum(!is.na(ischemic_cohort_complete$date_chronic_ischemic) & ischemic_cohort_complete$date_chronic_ischemic <= ischemic_cohort_complete$date_img1), "\n")
cat("- 慢性缺血性日期 在 扫描1-2之间:", sum(!is.na(ischemic_cohort_complete$date_chronic_ischemic) & ischemic_cohort_complete$date_chronic_ischemic > ischemic_cohort_complete$date_img1 & ischemic_cohort_complete$date_chronic_ischemic <= ischemic_cohort_complete$date_img2), "\n")
cat("- 慢性缺血性日期 > 扫描2:", sum(!is.na(ischemic_cohort_complete$date_chronic_ischemic) & ischemic_cohort_complete$date_chronic_ischemic > ischemic_cohort_complete$date_img2), "\n")
cat("- 新发MI:", sum(ischemic_cohort_complete$incident_MI), "\n")
cat("- 新发慢性缺血性:", sum(ischemic_cohort_complete$incident_chronic_ischemic), "\n")
cat("- 任何新发缺血性:", sum(ischemic_cohort_complete$incident_any_ischemic), "\n")
cat("最终数据维度:", dim(ischemic_cohort_complete), "\n")

# 保存完整的分析数据
write.csv(ischemic_cohort_complete, "final_ischemic_cohort.csv", row.names = FALSE)
saveRDS(ischemic_cohort_complete, "final_ischemic_cohort.rds")

# 简洁队列纳入/排除汇报（最小改动）
report_cohort_flow <- function(merged_all, final_cohort) {
  n_initial <- nrow(merged_all)
  n_with_both_imgs <- merged_all %>%
    dplyr::filter(`Date.of.attending.assessment.centre...Instance.2` != "" &
                    `Date.of.attending.assessment.centre...Instance.3` != "") %>%
    nrow()
  
  pre <- merged_all %>%
    dplyr::filter(`Date.of.attending.assessment.centre...Instance.2` != "" &
                    `Date.of.attending.assessment.centre...Instance.3` != "") %>%
    dplyr::mutate(
      date_baseline = parse_ukb_date(`Date.of.attending.assessment.centre...Instance.0`),
      date_img1 = parse_ukb_date(`Date.of.attending.assessment.centre...Instance.2`),
      date_img2 = parse_ukb_date(`Date.of.attending.assessment.centre...Instance.3`),
      date_MI = parse_ukb_date(`Date.of.myocardial.infarction`),
      date_chronic_ischemic = parse_ukb_date(`Date.I25.first.reported..chronic.ischaemic.heart.disease.`),
      date_angina = parse_ukb_date(`Date.I20.first.reported..angina.pectoris.`),
      date_stroke = parse_ukb_date(`Date.of.stroke`),
      date_dementia = parse_ukb_date(`Date.of.all.cause.dementia.report`),
      followup_years = as.numeric(date_img2 - date_img1) / 365.25
    ) %>%
    dplyr::filter(!is.na(date_img1), !is.na(date_img2), followup_years > 0) %>%
    dplyr::mutate(
      MI_at_baseline = ifelse(!is.na(date_MI) & date_MI <= date_img1, 1, 0),
      chronic_at_baseline = ifelse(!is.na(date_chronic_ischemic) & date_chronic_ischemic <= date_img1, 1, 0),
      angina_at_baseline = ifelse(!is.na(date_angina) & date_angina <= date_img1, 1, 0),
      stroke_at_baseline = ifelse(!is.na(date_stroke) & date_stroke <= date_img1, 1, 0),
      dementia_at_baseline = ifelse(!is.na(date_dementia) & date_dementia <= date_img1, 1, 0)
    )
  
  n_followup_valid <- nrow(pre)
  n_excluded_MI_baseline <- sum(pre$MI_at_baseline == 1, na.rm = TRUE)
  n_excluded_chronic_baseline <- sum(pre$chronic_at_baseline == 1, na.rm = TRUE)
  n_excluded_angina_baseline <- sum(pre$angina_at_baseline == 1, na.rm = TRUE)
  n_excluded_stroke_baseline <- sum(pre$stroke_at_baseline == 1, na.rm = TRUE)
  n_excluded_dementia_baseline <- sum(pre$dementia_at_baseline == 1, na.rm = TRUE)
  n_final <- nrow(final_cohort)
  
  flow_df <- tibble::tibble(
    Stage = c(
      "Initial merged",
      "With both imaging dates",
      "Valid follow-up (>0 years)",
      "Excluded: MI at baseline",
      "Excluded: Chronic ischemic at baseline",
      "Excluded: Angina at baseline",
      "Excluded: Stroke at baseline",
      "Excluded: Dementia at baseline",
      "Final analysis cohort"
    ),
    N = c(
      n_initial,
      n_with_both_imgs,
      n_followup_valid,
      n_excluded_MI_baseline,
      n_excluded_chronic_baseline,
      n_excluded_angina_baseline,
      n_excluded_stroke_baseline,
      n_excluded_dementia_baseline,
      n_final
    )
  )
  
  readr::write_csv(flow_df, file.path(getwd(), "Cohort_Flow_Report.csv"))
  cat("\n=== 队列纳入/排除简洁汇报 ===\n")
  print(flow_df)
}

# 立即输出纳入/排除简报
report_cohort_flow(merged_data_complete, ischemic_cohort_complete)

compute_cohort_flow_counts <- function(merged_all, final_cohort) {
  n_initial <- nrow(merged_all)
  n_with_both_imgs <- merged_all %>%
    dplyr::filter(`Date.of.attending.assessment.centre...Instance.2` != "" &
                    `Date.of.attending.assessment.centre...Instance.3` != "") %>%
    nrow()
  pre <- merged_all %>%
    dplyr::filter(`Date.of.attending.assessment.centre...Instance.2` != "" &
                    `Date.of.attending.assessment.centre...Instance.3` != "") %>%
    dplyr::mutate(
      date_baseline = parse_ukb_date(`Date.of.attending.assessment.centre...Instance.0`),
      date_img1 = parse_ukb_date(`Date.of.attending.assessment.centre...Instance.2`),
      date_img2 = parse_ukb_date(`Date.of.attending.assessment.centre...Instance.3`),
      date_MI = parse_ukb_date(`Date.of.myocardial.infarction`),
      date_chronic_ischemic = parse_ukb_date(`Date.I25.first.reported..chronic.ischaemic.heart.disease.`),
      date_angina = parse_ukb_date(`Date.I20.first.reported..angina.pectoris.`),
      date_stroke = parse_ukb_date(`Date.of.stroke`),
      date_dementia = parse_ukb_date(`Date.of.all.cause.dementia.report`),
      followup_years = as.numeric(date_img2 - date_img1) / 365.25
    ) %>%
    dplyr::filter(!is.na(date_img1), !is.na(date_img2), followup_years > 0) %>%
    dplyr::mutate(
      MI_at_baseline = ifelse(!is.na(date_MI) & date_MI <= date_img1, 1, 0),
      chronic_at_baseline = ifelse(!is.na(date_chronic_ischemic) & date_chronic_ischemic <= date_img1, 1, 0),
      angina_at_baseline = ifelse(!is.na(date_angina) & date_angina <= date_img1, 1, 0),
      stroke_at_baseline = ifelse(!is.na(date_stroke) & date_stroke <= date_img1, 1, 0),
      dementia_at_baseline = ifelse(!is.na(date_dementia) & date_dementia <= date_img1, 1, 0)
    )
  list(
    initial_total = n_initial,
    with_both_imgs = n_with_both_imgs,
    followup_valid = nrow(pre),
    excluded_MI_baseline = sum(pre$MI_at_baseline == 1, na.rm = TRUE),
    excluded_chronic_baseline = sum(pre$chronic_at_baseline == 1, na.rm = TRUE),
    excluded_angina_baseline = sum(pre$angina_at_baseline == 1, na.rm = TRUE),
    excluded_stroke_baseline = sum(pre$stroke_at_baseline == 1, na.rm = TRUE),
    excluded_dementia_baseline = sum(pre$dementia_at_baseline == 1, na.rm = TRUE),
    final_cohort_n = nrow(final_cohort)
  )
}

idp_miss_dir <- file.path(getwd(), "IDP_Missingness")
if (!dir.exists(idp_miss_dir)) dir.create(idp_miss_dir)
miss_Group <- idp_missingness_by_group(ischemic_cohort_complete, group_var = "Group")
miss_subtype <- idp_missingness_by_group(ischemic_cohort_complete, group_var = "disease_subtype")
readr::write_csv(miss_Group, file.path(idp_miss_dir, "IDP_missingness_by_Group.csv"))
readr::write_csv(miss_subtype, file.path(idp_miss_dir, "IDP_missingness_by_disease_subtype.csv"))

psm_covariates <- c(
  "age_at_scan_1",
  "sex_factor",
  "ethnicity_factor",
  "education_factor",
  "imaging_center_factor"
)
baseline_covariates <- c(
  "age_at_scan_1", "age_at_scan_2", "sex_factor", "ethnicity_factor", "education_factor",
  "bmi_baseline", "systolic_bp_baseline", "diastolic_bp_baseline",
  "townsend_index_baseline", "smoking_factor_baseline", "alcohol_factor_baseline",
  "diabetes_factor_baseline", "imaging_center_factor",
  "waist_circumference_baseline", "hip_circumference_baseline", "waist_hip_ratio_baseline",
  "body_fat_percentage_baseline", "whole_body_fat_mass_baseline", "followup_years"
)
flow_counts <- compute_cohort_flow_counts(merged_data_complete, ischemic_cohort_complete)
psm_res <- run_psm_comparison(ischemic_cohort_complete, "ischemic_vs_control", psm_covariates, baseline_covariates, m = 5, ratio = 2, caliper = 0.2, flow_counts = flow_counts)
if (!is.null(psm_res$matched_example)) {
  matched_psm <- psm_res$matched_example
  readr::write_csv(matched_psm, file.path(getwd(), "final_psm_cohort_ischemic_vs_control.csv"))
  saveRDS(matched_psm, file.path(getwd(), "final_psm_cohort_ischemic_vs_control.rds"))
}

cat("\n🎉 队列建立流程成功！现在有了包含所有变量的完整分析队列\n")
cat("进一步基线表生成和后续的认知/脑影像分析！\n")
