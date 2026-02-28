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

get_table_dir <- function(...) {
  dir <- file.path(getwd(), "Output_Tables", ...)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  dir
}

table_out_path <- function(subdir, filename) {
  file.path(get_table_dir(subdir), filename)
}

resolve_existing_file <- function(paths) {
  p <- paths[file.exists(paths)]
  if (length(p) > 0) p[[1]] else NA_character_
}

cohort_rds_path <- function() {
  table_out_path("Cohort", "final_ischemic_cohort.rds")
}

matched_cohort_rds_path <- function(model_label) {
  table_out_path(file.path("Cohort", "Matched"), paste0("Final_Matched_Cohort_", model_label, ".rds"))
}

psm_flow_rds_path <- function(model_label) {
  table_out_path(file.path("Cohort", "PSM"), paste0("PSM_Flow_Data_", model_label, ".rds"))
}

combined_four_model_results_csv <- function() {
  table_out_path("Four_Model", "Combined_Four_Model_Z_Analysis_Results.csv")
}

all_independent_effect_idps_csv <- function() {
  table_out_path("Four_Model", "All_Independent_Effect_IDPs_for_Visualization.csv")
}

clean_combined_results_df <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  if (anyDuplicated(names(df)) > 0) names(df) <- make.unique(names(df))
  if (all(c("z_basic_on_extended_sample", "z_statistic_without_non_imaging") %in% names(df))) {
    zb <- suppressWarnings(as.numeric(df$z_basic_on_extended_sample))
    zw <- suppressWarnings(as.numeric(df$z_statistic_without_non_imaging))
    same <- (is.na(zb) & is.na(zw)) | (is.finite(zb) & is.finite(zw) & abs(zb - zw) < 1e-12)
    if (all(same)) df$z_basic_on_extended_sample <- NULL
  }
  df
}

# =================== 模块化运行控制开关 ===================
# 将以下开关设置为 FALSE 可跳过相应步骤（前提是已存在中间数据文件）
env_flag <- function(name, default = FALSE) {
  v <- Sys.getenv(name, unset = "")
  if (!nzchar(v)) return(default)
  tolower(trimws(v)) %in% c("1", "true", "t", "yes", "y")
}
RUN_DATA_FUSION     <- env_flag("RUN_DATA_FUSION", FALSE)   # 1. 数据读取与合并 (生成 final_ischemic_cohort.rds)
RUN_PSM_MATCHING    <- env_flag("RUN_PSM_MATCHING", TRUE)   # 2. PSM倾向性评分匹配 (生成 Final_Matched_Cohort_*.rds)
RUN_IDP_ANALYSIS    <- env_flag("RUN_IDP_ANALYSIS", TRUE)   # 3. IDP统计分析 (生成 Combined_Four_Model_Z_Analysis_Results.csv)
RUN_VISUALIZATION   <- env_flag("RUN_VISUALIZATION", TRUE)  # 4. 可视化 (生成各类图表)
RUN_BRAIN_MAPPING   <- env_flag("RUN_BRAIN_MAPPING", TRUE)  # 5. 脑图映射 (生成脑区映射图)
RUN_REPLOT_SCAN_RESCAN <- env_flag("RUN_REPLOT_SCAN_RESCAN", FALSE)
RUN_REPLOT_IDP_ZDIST   <- env_flag("RUN_REPLOT_IDP_ZDIST", FALSE)

# 设置工作目录
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

if (RUN_DATA_FUSION) {

cat("=== 重新开始：统一ID后直接合并所有数据 ===\n")

# 读取数据
cat("重新读取数据...\n")
data1 <- vroom::vroom("data1.csv", delim = ",", na = "NA")
data2 <- vroom::vroom("data2.csv", delim = ",", na = "NA")
data3 <- vroom::vroom("combined_brain_imaging_data.csv", delim = ",", na = "NA")
names(data1) <- make.names(names(data1))
names(data2) <- make.names(names(data2))
names(data3) <- make.names(names(data3))


# 统一所有数据集的ID列名
cat("统一ID列名...\n")
if("Participant.ID" %in% names(data1)) {
  data1 <- data1 %>% rename(eid = `Participant.ID`)
  cat("✅ data1: Participant.ID → eid\n")
}
if("Participant.ID" %in% names(data2)) {
  data2 <- data2 %>% rename(eid = `Participant.ID`)
  cat("✅ data2: Participant.ID → eid\n")
}

if("participant.id" %in% names(data3)) {
  data3 <- data3 %>% rename(eid = `participant.id`)
  cat("✅ data3: participant.id → eid\n")
} else if("Participant.ID" %in% names(data3)) {
  data3 <- data3 %>% rename(eid = `Participant.ID`)
  cat("✅ data3: Participant.ID → eid\n")
}

# 检查ID列统一情况
cat("检查ID列统一情况:\n")
cat("- data1样本量:", nrow(data1), "，是否有eid:", "eid" %in% names(data1), "\n")
cat("- data2样本量:", nrow(data2), "，是否有eid:", "eid" %in% names(data2), "\n") 
cat("- data3样本量:", nrow(data3), "，是否有eid:", "eid" %in% names(data3), "\n")

# 直接全连接合并所有数据
cat("\n直接合并所有数据集...\n")
merged_data_complete <- data1 %>%
  full_join(data2, by = "eid") %>%
  full_join(data3, by = "eid")

cat("完整合并后数据维度:", dim(merged_data_complete), "\n")
cat("合并后样本量:", nrow(merged_data_complete), "\n")
cat("合并后变量数:", ncol(merged_data_complete), "\n")

# 检查关键变量可用性
cat("\n检查关键变量可用性:\n")
key_vars <- c(
  "Date.of.attending.assessment.centre...Instance.2",
  "Date.of.attending.assessment.centre...Instance.3",
  "Date.of.myocardial.infarction",
  "Date.I25.first.reported..chronic.ischaemic.heart.disease."
)

for (var in key_vars) {
  if (var %in% names(merged_data_complete)) {
    col_vals <- merged_data_complete[[var]]
    if (is.character(col_vals)) {
      non_empty <- sum(!is.na(col_vals) & trimws(col_vals) != "")
    } else {
      non_empty <- sum(!is.na(col_vals))
    }
    cat("✅", var, ":", non_empty, "\n")
  } else {
    cat("❌", var, ": 不存在\n")
  }
}
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
# 在完整的数据上构建队列
cat("\n=== 基于完整合并数据构建缺血性心脏病队列 ===\n")

ischemic_cohort_complete <- unify_joined_columns(merged_data_complete) %>%
  filter(`Date.of.attending.assessment.centre...Instance.2` != "" &
           `Date.of.attending.assessment.centre...Instance.3` != "") %>%
  mutate(
    date_baseline = parse_ukb_date(`Date.of.attending.assessment.centre...Instance.0`), 
    date_img1 = parse_ukb_date(`Date.of.attending.assessment.centre...Instance.2`),
    date_img2 = parse_ukb_date(`Date.of.attending.assessment.centre...Instance.3`),
    date_MI = parse_ukb_date(`Date.of.myocardial.infarction`),
    date_chronic_ischemic = parse_ukb_date(`Date.I25.first.reported..chronic.ischaemic.heart.disease.`),
    date_angina = parse_ukb_date(`Date.I20.first.reported..angina.pectoris.`),
    date_stroke = parse_ukb_date(`Date.of.stroke`),
    date_dementia = parse_ukb_date(`Date.of.all.cause.dementia.report`),
    
    followup_years = as.numeric(date_img2 - date_img1) / 365.25,
    # 供PSM匹配使用：第一次扫描日期（数值）
    date_img1_numeric = as.numeric(date_img1),
    baseline_age = `Age.when.attended.assessment.centre...Instance.2`,
    followup_age = `Age.when.attended.assessment.centre...Instance.3`
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
    age_at_recruitment = baseline_age,
    
    sex = case_when(
      as.character(Sex) == "Male" ~ "Male", 
      as.character(Sex) == "Female" ~ "Female", 
      TRUE ~ "Male"
    ),
    sex_factor = factor(sex, levels = c("Female", "Male")),
    
    townsend_index = as.numeric(`Townsend.deprivation.index.at.recruitment`),
    
    education_raw = as.character(`Qualifications...Instance.0`),
    education = case_when(
      grepl("College|University|degree", education_raw, ignore.case = TRUE) ~ "Higher",
      grepl("A levels|O levels|GCSEs|NVQ|HND|HNC|Professional", education_raw, ignore.case = TRUE) ~ "Medium",
      TRUE ~ "Lower"
    ),
    education_factor = factor(education, levels = c("Lower", "Medium", "Higher")),
    
    ethnicity_raw = as.character(`Ethnic.background...Instance.0`),
    ethnicity = ifelse(grepl("White", ethnicity_raw, ignore.case = TRUE), "White", "Non-White"),
    ethnicity_factor = factor(ethnicity, levels = c("White", "Non-White")),
    
    bmi = as.numeric(`Body.mass.index..BMI....Instance.0`),
    systolic_bp = {
      sys_cols <- grep("Systolic\\.blood\\.pressure.*Instance\\.0.*Array\\.[0-9]+$", names(cur_data_all()), ignore.case = TRUE, value = TRUE)
      if (length(sys_cols) > 0) {
        suppressWarnings(rowMeans(as.matrix(cur_data_all()[, sys_cols, drop = FALSE]), na.rm = TRUE))
      } else {
        suppressWarnings(as.numeric(`Systolic.blood.pressure..automated.reading...Instance.0...Array.0`))
      }
    },
    diastolic_bp = {
      dia_cols <- grep("Diastolic\\.blood\\.pressure.*Instance\\.0.*Array\\.[0-9]+$", names(cur_data_all()), ignore.case = TRUE, value = TRUE)
      if (length(dia_cols) > 0) {
        suppressWarnings(rowMeans(as.matrix(cur_data_all()[, dia_cols, drop = FALSE]), na.rm = TRUE))
      } else {
        suppressWarnings(as.numeric(`Diastolic.blood.pressure..automated.reading...Instance.0...Array.0`))
      }
    },
    
    diabetes = case_when(
      as.character(`Diabetes.diagnosed.by.doctor...Instance.0`) == "Yes" ~ "Yes",
      TRUE ~ "No"
    ),
    diabetes_factor = factor(diabetes, levels = c("No", "Yes")),
    
    smoking = case_when(
      as.character(`Smoking.status...Instance.0`) == "Never" ~ "Never",
      as.character(`Smoking.status...Instance.0`) == "Current" ~ "Current", 
      TRUE ~ "Former"
    ),
    smoking_factor = factor(smoking, levels = c("Never", "Former", "Current")),
    
    alcohol = case_when(
      as.character(`Alcohol.intake.frequency....Instance.0`) %in% c("Never", "Special occasions only") ~ "Never",
      as.character(`Alcohol.intake.frequency....Instance.0`) %in% c("Once or twice a week", "Three or four times a week", "Daily or almost daily") ~ "Current",
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
    )
  ) %>%
  # 处理缺失值
  mutate(
    bmi = ifelse(is.na(bmi) | !is.finite(bmi), median(bmi, na.rm = TRUE), bmi),
    systolic_bp = ifelse(is.na(systolic_bp) | !is.finite(systolic_bp), median(systolic_bp, na.rm = TRUE), systolic_bp),
    diastolic_bp = ifelse(is.na(diastolic_bp) | !is.finite(diastolic_bp), 
                          if ("diastolic_bp_baseline" %in% names(.)) diastolic_bp_baseline else median(diastolic_bp, na.rm = TRUE),
                          diastolic_bp),
    townsend_index = ifelse(is.na(townsend_index) | !is.finite(townsend_index), median(townsend_index, na.rm = TRUE), townsend_index)
  )

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
write.csv(ischemic_cohort_complete, table_out_path("Cohort", "final_ischemic_cohort.csv"), row.names = FALSE)
saveRDS(ischemic_cohort_complete, cohort_rds_path())

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
  
  readr::write_csv(flow_df, table_out_path("Cohort", "Cohort_Flow_Report.csv"))
  cat("\n=== 队列纳入/排除简洁汇报 ===\n")
  print(flow_df)
}

# 立即输出纳入/排除简报
report_cohort_flow(merged_data_complete, ischemic_cohort_complete)

cat("\n🎉 队列建立流程成功！现在有了包含所有变量的完整分析队列\n")
cat("进一步基线表生成和后续的认知/脑影像分析！\n")

} else {
  cat("=== 跳过数据融合步骤，尝试加载已有队列数据 ===\n")
  cohort_rds <- resolve_existing_file(c(cohort_rds_path(), "final_ischemic_cohort.rds"))
  if (!is.na(cohort_rds) && nzchar(cohort_rds)) {
    ischemic_cohort_complete <- readRDS(cohort_rds)
    cat(sprintf("✅ 已加载: %s\n", cohort_rds))
  } else {
    cat("⚠️ 未找到 final_ischemic_cohort.rds，若后续步骤需要此数据可能会报错。\n")
  }
}

# ===================多重插补 + PSM + 纵向ΔIDP分析 + 可视化=========================

suppressPackageStartupMessages({
  library(mice)
})

# ----- 工具函数 -----
# 统一出版主题（Nature风格）与保存工具
theme_pub <- function(base_size = 12, base_family = "Helvetica") {
  # 为保证全局一致的出版风格，将原通用主题代理为 JAMA 风格
  theme_jama_refined(base_size = base_size, base_family = base_family)
}

pub_group_colors <- c(
  "Control" = "#1F78B4",
  "Controls" = "#1F78B4",
  "Ischemic_Heart_Disease" = "#FF7F00",
  "Cases" = "#FF7F00",
  "Myocardial_Infarction" = "#E31A1C",
  "Chronic_Ischemic" = "#33A02C"
)

# 分组标签与颜色统一辅助函数
factorize_group <- function(group, labels = c("Controls","Cases")) {
  gnum <- suppressWarnings(as.numeric(as.character(group)))
  factor(gnum, levels = c(0, 1), labels = labels)
}

get_group_colors <- function() {
  defaults <- c(
    "Control" = "#1F78B4",
    "Ischemic_Heart_Disease" = "#FF7F00",
    "Myocardial_Infarction" = "#E31A1C",
    "Chronic_Ischemic" = "#33A02C",
    "Controls" = "#1F78B4",
    "Cases" = "#FF7F00"
  )
  if (exists("pub_group_colors") && is.character(pub_group_colors)) {
    for (nm in names(defaults)) {
      if (nm %in% names(pub_group_colors)) defaults[nm] <- pub_group_colors[[nm]]
    }
  }
  defaults
}

# SCI publication theme (refined minimal with clear axes, subdued grid, readable legend)
theme_jama_refined <- function(base_size = 12, base_family = "Helvetica") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "#E5E5E5", linewidth = 0.6),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = NA, fill = NA),
      axis.title = ggplot2::element_text(face = "bold", color = "#111111"),
      axis.text = ggplot2::element_text(color = "#222222"),
      axis.line = ggplot2::element_line(color = "#333333", linewidth = 0.6),
      axis.ticks = ggplot2::element_line(color = "#333333"),
      axis.ticks.length = grid::unit(3, "pt"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = base_size + 2, color = "#111111"),
      plot.subtitle = ggplot2::element_text(hjust = 0, size = base_size, color = "#444444"),
      plot.caption = ggplot2::element_text(hjust = 0, size = base_size - 2, color = "#555555"),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold"),
      legend.key.width = grid::unit(10, "pt"),
      legend.key.height = grid::unit(10, "pt")
    )
}

jama_colors <- list(
  control = "#1F78B4",
  case    = "#E31A1C",
  before  = "#34495E",
  after   = "#27AE60"
)



check_dir_writable <- function(dir_path) {
  tryCatch({
    if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    tf <- tempfile(pattern = ".write_test_", tmpdir = dir_path)
    ok <- suppressWarnings(file.create(tf))
    if (isTRUE(ok)) unlink(tf)
    isTRUE(ok)
  }, error = function(e) FALSE)
}

save_pub <- function(filename, plot, width = 10, height = 6, dpi = 300, formats = c("png","pdf")) {
  base <- sub("\n$", "", filename)
  dir_path <- dirname(base)
  if (!check_dir_writable(dir_path)) {
    warning(sprintf("目录不可写或创建失败：%s", dir_path))
  }
  for (fmt in formats) {
    out <- paste0(tools::file_path_sans_ext(base), ".", fmt)
    tryCatch({
      ggplot2::ggsave(out, plot, width = width, height = height, dpi = dpi, units = "in")
      cat(sprintf("✓ 已保存: %s\n", out))
    }, error = function(e) {
      cat(sprintf("❌ 保存失败: %s — %s\n", out, e$message))
    })
  }
}

# JAMA风格CSV导出：统一四舍五入数值列，便于论文表格使用
jama_round_df <- function(df, digits = 3) {
  if (is.null(df)) return(df)
  # 确保为data.frame，并将行名写入第一列，避免CSV丢失行名
  df <- as.data.frame(df)
  rn <- rownames(df)
  if (!is.null(rn) && any(nchar(rn) > 0)) {
    df <- cbind(Variable = rn, df)
    rownames(df) <- NULL
  }
  dplyr::mutate(df, dplyr::across(where(is.numeric), ~ round(., digits)))
}


theme_brain_pub <- function(){
  ggplot2::theme_void(base_family = "Helvetica") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 16, hjust = 0, face = "bold", color = "#111111"),
      plot.subtitle = ggplot2::element_text(size = 12, hjust = 0, color = "#444444"),
      plot.caption = ggplot2::element_text(size = 10, hjust = 0, color = "#555555"),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold"),
      legend.key.size = grid::unit(5, "mm")
    )
}

# 输出目录写入权限校验与汇总
verify_output_dirs <- function(dirs = c(
  file.path(getwd(), "Brain_MRI_Maps"),
  file.path(getwd(), "Age_Trend_Plots"),
  file.path(getwd(), "Individual_IDP_Plots"),
  file.path(getwd(), "Brain_Changes_Overview")
)) {
  res <- lapply(dirs, function(d) list(dir = d, writable = check_dir_writable(d)))
  df <- data.frame(
    Directory = vapply(res, function(x) x$dir, character(1)),
    Writable  = vapply(res, function(x) x$writable, logical(1)),
    stringsAsFactors = FALSE
  )
  out <- file.path(getwd(), "Visualization_Dir_Write_Check.csv")
  utils::write.csv(df, out, row.names = FALSE)
  cat(sprintf("✓ 目录写入权限检查完成：%s\n", out))
  print(df)
  invisible(df)
}

# 简易绘图调用日志
log_plot_call <- function(name) {
  cat(sprintf("→ 调用绘图函数: %s\n", name))
}

## 通用：PSM匹配封装，供全局与MI模块复用
psm_match_standard <- function(data,
                               treatment_col,
                               covariates,
                               method = "nearest",
                               distance = "glm",
                               ratio = 2,
                               caliper = 0.2,
                               replace = FALSE,
                               estimand = "ATT",
                               seed = 2024) {
  stopifnot(is.character(treatment_col), length(treatment_col) == 1)
  stopifnot(is.character(covariates))
  
  missing_vars <- setdiff(c(treatment_col, covariates), names(data))
  if (length(missing_vars) > 0) {
    stop(sprintf("PSM匹配缺少变量: %s", paste(missing_vars, collapse = ", ")))
  }
  
  terms <- paste(covariates, collapse = " + ")
  if (nchar(terms) == 0) stop("PSM匹配需要至少一个协变量")
  ps_formula <- stats::as.formula(paste(treatment_col, "~", terms))
  
  set.seed(seed)
  mobj <- MatchIt::matchit(
    formula = ps_formula,
    data = data,
    method = method,
    distance = distance,
    ratio = ratio,
    caliper = caliper,
    std.caliper = FALSE,
    replace = replace,
    estimand = estimand
  )
  matched <- MatchIt::match.data(mobj)
  list(mobj = mobj, matched_data = matched)
}

## 诊断：计算匹配前后SMD与基线表汇总
compute_psm_balance <- function(original, matched, covariates, group_col = "group") {
  # 安全过滤：仅保留在original与matched同时存在的协变量，避免未定义列错误
  available_covs <- covariates[covariates %in% names(original) & covariates %in% names(matched)]
  if (length(available_covs) == 0) {
    cat("[PSM-诊断] 协变量在数据集中不存在，跳过该诊断\n")
    empty_df <- data.frame()
    return(list(
      smd_plot_data = data.frame(variable = character(), smd = numeric(), phase = character(), stringsAsFactors = FALSE),
      overall_stats = empty_df,
      table_before = empty_df,
      table_after = empty_df,
      p_values = empty_df
    ))
  }
  factor_vars <- available_covs[sapply(original[, available_covs, drop = FALSE], is.factor)]
  # 匹配前
  tbl_before <- tableone::CreateTableOne(vars = available_covs, factorVars = factor_vars,
                                         strata = group_col, data = original, test = TRUE, smd = TRUE)
  df_before <- as.data.frame(print(tbl_before, printToggle = FALSE, showAllLevels = TRUE, test = TRUE, smd = TRUE))
  df_before$phase <- "before"
  # 匹配后
  tbl_after <- tableone::CreateTableOne(vars = available_covs, factorVars = factor_vars,
                                        strata = group_col, data = matched, test = TRUE, smd = TRUE)
  df_after <- as.data.frame(print(tbl_after, printToggle = FALSE, showAllLevels = TRUE, test = TRUE, smd = TRUE))
  df_after$phase <- "after"
  # 提取并聚合SMD（若含水平行，按变量名聚合取最大值）
  extract_smd <- function(df) {
    nm <- rownames(df)
    col_smd <- dplyr::case_when(
      "SMD" %in% names(df) ~ "SMD",
      "StdDif" %in% names(df) ~ "StdDif",
      TRUE ~ NA_character_
    )
    if (is.na(col_smd)) {
      smd <- rep(NA_real_, length(nm))
    } else {
      smd_raw <- as.character(df[[col_smd]])
      smd_raw <- trimws(smd_raw)
      smd <- suppressWarnings(as.numeric(gsub("[^0-9eE.+-]", "", smd_raw)))
    }
    data.frame(variable = nm, smd = smd, phase = df$phase, stringsAsFactors = FALSE)
  }
  smd_before <- extract_smd(df_before)
  smd_after  <- extract_smd(df_after)
  smd_all <- dplyr::bind_rows(smd_before, smd_after) %>%
    # 统一变量名：去掉因子水平（=...）与连续变量的摘要后缀（例如" (mean (sd))"）
    dplyr::mutate(variable = gsub("=.*$|\\s*\\(.*\\)$", "", variable)) %>%
    dplyr::filter(
      !variable %in% c("N", "n", "Missing"),
      !grepl("^X(\\.\\d+)?$", variable)
    ) %>%
    dplyr::group_by(variable, phase) %>%
    dplyr::summarise(
      smd = suppressWarnings({
        vals <- smd
        vals <- vals[is.finite(vals)]
        if (length(vals) == 0) NA_real_ else vals[which.max(abs(vals))]
      }),
      .groups = "drop"
    ) %>%
    tidyr::complete(variable, phase = c("before", "after"))
  overall <- smd_all %>%
    tidyr::pivot_wider(names_from = phase, values_from = smd) %>%
    dplyr::mutate(
      before = as.numeric(before),
      after = as.numeric(after),
      smd_improvement = before - after,
      well_balanced = dplyr::case_when(
        is.na(after) ~ FALSE,
        TRUE ~ abs(after) < 0.1
      )
    ) %>%
    dplyr::summarise(
      n_variables = dplyr::n(),
      n_balanced = sum(well_balanced, na.rm = TRUE),
      balance_rate = mean(well_balanced, na.rm = TRUE) * 100,
      mean_smd_before = mean(abs(before), na.rm = TRUE),
      mean_smd_after = mean(abs(after), na.rm = TRUE),
      mean_improvement = mean(smd_improvement, na.rm = TRUE),
      .groups = "drop"
    )
  # 提取变量层级的p值（若存在）
  extract_pvals <- function(df) {
    pcol <- if ("p" %in% names(df)) {
      "p"
    } else if ("p.value" %in% names(df)) {
      "p.value"
    } else if ("pValue" %in% names(df)) {
      "pValue"
    } else {
      NA_character_
    }
    v <- data.frame(variable_raw = rownames(df), stringsAsFactors = FALSE)
    v$variable <- gsub("=.*$|\\s*\\(.*\\)$", "", v$variable_raw)
    if (!is.na(pcol)) {
      pv_raw <- as.character(df[[pcol]])
      pv_raw <- trimws(pv_raw)
      pv_num <- suppressWarnings(as.numeric(gsub("[^0-9eE.+-]", "", pv_raw)))
      v$p <- pv_num
    } else {
      v$p <- NA_real_
    }
    v <- v %>% dplyr::filter(
      !variable %in% c("N","n","Missing"),
      !grepl("^X(\\.\\d+)?$", variable)
    )
    dplyr::group_by(v, variable) %>%
      dplyr::summarise(p = {
        vv <- v$p[v$variable == unique(variable)[1]]
        vv <- vv[is.finite(vv)]
        if (length(vv) == 0) NA_real_ else vv[1]
      }, .groups = "drop")
  }
  p_before <- extract_pvals(df_before)
  p_after  <- extract_pvals(df_after)
  p_summary <- suppressMessages(dplyr::full_join(p_before, p_after, by = "variable", suffix = c("_before", "_after")))
  # 规范化并仅保留协变量；添加专业标签以便直接发表
  cov_keys <- canonicalize_covariate_id(available_covs)
  label_map <- c(
    age_at_scan_1   = "Age at scan 1",
    age_at_recruitment = "Age at recruitment",
    sex_factor      = "Sex (male/female)",
    ethnicity_factor= "Ethnicity (white/non-white)",
    education_factor = "Education (lower/medium/higher)",
    imaging_center_factor = "Location of first imaging center (Instance 2)"
  )
  p_summary <- p_summary %>%
    dplyr::mutate(variable_key = canonicalize_covariate_id(variable),
                  variable_label = label_map[variable_key]) %>%
    dplyr::filter(variable_key %in% cov_keys) %>%
    dplyr::select(variable_key, variable_label, p_before, p_after) %>%
    dplyr::group_by(variable_key) %>%
    dplyr::summarise(
      variable_label = dplyr::first(na.omit(variable_label)),
      p_before = dplyr::first(na.omit(p_before)),
      p_after  = dplyr::first(na.omit(p_after)),
      .groups = "drop"
    ) %>%
    dplyr::arrange(match(variable_key, cov_keys))
  list(smd_plot_data = smd_all, overall_stats = overall,
       table_before = df_before, table_after = df_after,
       p_values = p_summary)
}

## 辅助：规范化协变量名称以对齐不同命名（用于 Love Plot 过滤与标注）
canonicalize_covariate_id <- function(x) {
  if (is.null(x) || length(x) == 0) return(character(0))
  # 去除因子水平与摘要后缀，统一为小写+下划线
  base <- tolower(gsub("=.*$|\\s*\\(.*\\)$", "", x))
  base <- gsub("[^a-z0-9_]+", "_", base)      # 将非字母数字统一为下划线（如 ..mean..SD.. -> _mean_SD_）
  base <- gsub("^_+|_+$", "", base)
  base <- gsub("_+", "_", base)
  # 移除常见统计摘要后缀，保留变量根名
  # 例如 baseline_age_mean_sd, followup_years_median_iqr, bmi_range 等
  base <- gsub("_(mean|sd|se|median|iqr|range|min|max|p|pvalue|count|counts|percent|percentage)(_[a-z0-9]+)*$", "", base)
  base <- gsub("_(mean_sd|median_iqr)$", "", base)
  map_one <- function(s) {
    # 允许根名后带有其它尾缀（如 _baseline, _male 等），只要匹配到根名即可
    # 年龄：优先匹配扫描2；扫描1仅在变量名严格匹配时归类
    if (grepl("^(age2|age_at_scan_?2|age_at_instance3)(_|$)", s)) return("age_at_scan_2")
    if (grepl("^(baseline_age|age_at_scan_?1|age_at_instance2)(_|$)|^age$", s)) return("age_at_scan_1")
    if (grepl("^(sex|sex_factor|gender|gender_factor)(_|$)", s)) return("sex_factor")
    if (grepl("^(ethnicity|ethnicity_factor|race|race_factor)(_|$)", s)) return("ethnicity_factor")
    if (grepl("^(followup_years|years_between_scans|scan_interval_years|time_between_scans|follow_up_years|followup|interval_years|elapsed_years)(_|$)", s)) return("followup_years")
    s
  }
  vapply(base, map_one, FUN.VALUE = character(1))
}


# ------------------------------------------------------------
# 高标准SCI版 Love Plot（只显示PSM匹配协变量，干净标签）
# ------------------------------------------------------------
plot_love_sci <- function(smd_plot_data, covariates, labels = NULL, title_suffix = NULL) {
  if (missing(covariates) || is.null(covariates) || length(covariates) == 0) {
    stop("plot_love_sci 需要提供用于PSM的协变量向量")
  }
  # 仅保留PSM协变量，且变量名统一去除水平尾缀（并做命名规范化以对齐）
  df <- smd_plot_data %>%
    dplyr::mutate(
      variable = gsub("=.*$|\\s*\\(.*\\)$", "", variable),
      variable_key = canonicalize_covariate_id(variable)
    )
  cov_keys <- unique(canonicalize_covariate_id(covariates))
  df <- df %>% dplyr::filter(variable_key %in% cov_keys)

  # 若筛选后为空，给出更友好的错误提示
  if (nrow(df) == 0) {
    present_keys <- unique(canonicalize_covariate_id(gsub("=.*$|\\s*\\(.*\\)$", "", smd_plot_data$variable)))
    missing_keys <- setdiff(cov_keys, present_keys)
    stop(sprintf("Love Plot 输入为空：检查协变量是否在诊断结果中。缺失：%s", paste(missing_keys, collapse = ", ")))
  }

  # 构建宽表：变量 × {before, after}
  smd_wide <- df %>%
    tidyr::pivot_wider(names_from = phase, values_from = smd, values_fill = NA) %>%
    {
      # 保险：若某一列不存在（极端情况下），补齐空列以避免 rename 报错
      if (!"before" %in% names(.)) .$before <- NA_real_
      if (!"after" %in% names(.))  .$after  <- NA_real_
      .
    } %>%
    dplyr::rename(smd_before = before, smd_after = after) %>%
    dplyr::mutate(variable_key = canonicalize_covariate_id(variable))

  # 默认标签映射（可用户传入覆盖）
  default_labels <- c(
    age_at_scan_1   = "Age at scan 1",
    age_at_recruitment = "Age at recruitment",
    sex_factor      = "Sex (male/female)",
    ethnicity_factor= "Ethnicity (white/non-white)",
    education_factor = "Education (lower/medium/higher)",
    imaging_center_factor = "Location of first imaging center (Instance 2)"
  )
  if (is.null(labels)) labels <- default_labels

  # 保持协变量的显式顺序
  label_vec <- unname(labels[cov_keys])
  # 若存在未覆盖的协变量键，使用键名作为回退标签
  if (any(is.na(label_vec))) {
    unknown_keys <- cov_keys[is.na(label_vec)]
    if (length(unknown_keys) > 0) {
      tmp_labels <- labels
      tmp_labels[unknown_keys] <- unknown_keys
      labels <- tmp_labels
      label_vec <- unname(labels[cov_keys])
    }
  }
  # 去除 NA 与重复，保持原始顺序
  label_vec <- label_vec[!is.na(label_vec)]
  label_vec <- label_vec[!duplicated(label_vec)]
  label_map <- setNames(label_vec, label_vec)  # 使用标签到标签的映射避免重复键
  # 若映射缺失，回退到 variable_key
  smd_wide <- smd_wide %>% dplyr::mutate(variable_label = {
    vl <- labels[variable_key]
    vl[is.na(vl)] <- variable_key[is.na(vl)]
    unname(vl)
  })
  # 确保因子水平唯一且不含NA
  smd_wide$variable_label <- factor(smd_wide$variable_label, levels = unique(c(label_vec, unique(smd_wide$variable_label))))

  # 转换为长表用于散点绘图
  smd_long <- smd_wide %>%
    dplyr::select(variable_label, smd_before, smd_after) %>%
    tidyr::pivot_longer(cols = c(smd_before, smd_after), names_to = "period", values_to = "smd") %>%
    dplyr::mutate(period = factor(period, levels = c("smd_before", "smd_after"), labels = c("Before Matching", "After Matching")))

  # 动态坐标范围
  max_val <- suppressWarnings(max(abs(smd_long$smd), na.rm = TRUE))
  if (!is.finite(max_val)) max_val <- 0.25
  limit_val <- max(0.25, ifelse(max_val > 1, ceiling(max_val*10)/10, ceiling(max_val*20)/20)) * 1.05

  # 颜色映射（JAMA风格）
  color_map <- c("Before Matching" = "#E64B35", "After Matching" = "#4DBBD5")
  if (exists("jama_colors")) {
    if(!is.null(jama_colors$before)) color_map["Before Matching"] <- jama_colors$before
    if(!is.null(jama_colors$after))  color_map["After Matching"]  <- jama_colors$after
  }

  p <- ggplot2::ggplot(smd_wide, ggplot2::aes(y = variable_label)) +
    ggplot2::annotate("rect", xmin = -0.1, xmax = 0.1, ymin = -Inf, ymax = Inf, fill = "gray92", alpha = 0.5) +
    ggplot2::geom_vline(xintercept = 0, linetype = "solid", color = "gray30", linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "gray50", linewidth = 0.3) +
    ggplot2::geom_segment(ggplot2::aes(x = smd_before, xend = smd_after, yend = variable_label),
                          color = "gray70", linewidth = 0.6, arrow = grid::arrow(length = grid::unit(0.15, "cm"), type = "closed")) +
    ggplot2::geom_point(data = dplyr::filter(smd_long, period == "Before Matching"),
                        ggplot2::aes(x = smd),
                        color = color_map["Before Matching"],
                        fill  = color_map["Before Matching"],
                        shape = 21,
                        size = 3.5, stroke = 0.8) +
    ggplot2::geom_point(data = dplyr::filter(smd_long, period == "After Matching"),
                        ggplot2::aes(x = smd),
                        color = color_map["After Matching"],
                        fill  = color_map["After Matching"],
                        shape = 19,
                        size = 3.5, stroke = 0.8) +
    ggplot2::scale_x_continuous(limits = c(-limit_val, limit_val), breaks = scales::pretty_breaks(n = 7), labels = scales::number_format(accuracy = 0.01)) +
    ggplot2::scale_y_discrete(name = NULL) +
    ggplot2::scale_color_manual(name = "", values = color_map) +
    ggplot2::scale_fill_manual(name = "", values = color_map) +
    ggplot2::scale_shape_manual(name = "", values = c("Before Matching" = 21, "After Matching" = 19)) +
    ggplot2::labs(
    title = paste0("Covariate Balance (PSM)", if (!is.null(title_suffix)) paste0(" ", title_suffix) else ""),
    subtitle = "Love plot showing standardized mean differences (SMD)",
    x = "Standardized Mean Difference",
    caption = paste0("PSM covariates: ", paste(unique(unname(labels[cov_keys])), collapse = ", "))
    ) +
    ggplot2::theme_minimal(base_size = 12, base_family = "sans") +
    ggplot2::theme(
      legend.position = "bottom",
      legend.margin = ggplot2::margin(t = -5),
      legend.text = ggplot2::element_text(size = 10, face = "bold"),
      axis.text.y = ggplot2::element_text(color = "black", size = 10),
      axis.text.x = ggplot2::element_text(color = "black", size = 9),
      axis.title.x = ggplot2::element_text(face = "bold", margin = ggplot2::margin(t = 10)),
      panel.grid.major.y = ggplot2::element_line(color = "gray95", linewidth = 0.3),
      panel.grid.major.x = ggplot2::element_line(color = "gray90", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0),
      plot.subtitle = ggplot2::element_text(size = 10, color = "gray40", margin = ggplot2::margin(b = 15)),
      plot.margin = ggplot2::margin(20, 20, 20, 20)
    )

  p
}

# ------------------------------------------------------------
# Comprehensive Brain Changes Overview (bubble chart)
# - Uses Brain_IDPs_Parsed_for_Mapping.csv to summarize mean z-change
# - Point size = number of significant IDPs per region
# - Color = anatomical structure class
# ------------------------------------------------------------
build_brain_changes_overview_data <- function(mapping_csv = file.path(getwd(), "Brain_MRI_Maps", "Brain_IDPs_Parsed_for_Mapping.csv")) {
  if (!file.exists(mapping_csv)) {
    stop(sprintf("Mapping CSV not found: %s", mapping_csv))
  }
  idp_map <- suppressMessages(readr::read_csv(mapping_csv, show_col_types = FALSE))
  # 统一脑区列，尽量回填避免空白聚合行
  coalesce_cols <- function(df, target, candidates) {
    exist <- intersect(candidates, names(df))
    if (length(exist) == 0) {
      df[[target]] <- NA_character_
      return(df)
    }
    acc <- df[[exist[1]]]
    if (length(exist) > 1) {
      for (nm in exist[-1]) acc <- dplyr::coalesce(acc, df[[nm]])
    }
    df[[target]] <- acc
    df
  }
  idp_map <- coalesce_cols(idp_map, "brain_region_final", c("brain_region_final", "ggseg_region", "brain_region", "anatomical_structure", "region"))
  
  # Normalize structure category
  normalize_structure <- function(atlas_type, structure_type, brain_region) {
    atlas_type <- tolower(as.character(atlas_type))
    structure_type <- tolower(as.character(structure_type))
    region <- tolower(as.character(brain_region))
    dplyr::case_when(
      atlas_type %in% c("white_matter", "wm", "white-matter") ~ "white_matter",
      atlas_type %in% c("cortical") ~ "cortical",
      structure_type %in% c("white matter") ~ "white_matter",
      structure_type %in% c("cortical") ~ "cortical",
      structure_type %in% c("subcortical") ~ "subcortical",
      structure_type %in% c("cerebellum") ~ "cerebellum",
      stringr::str_detect(region, "pons|midbrain|brainstem") ~ "brainstem",
      TRUE ~ "unknown"
    )
  }
  
  df <- idp_map %>%
    dplyr::filter(!is.na(brain_region_final), brain_region_final != "", !is.na(z_change)) %>%
    dplyr::mutate(
      structure = normalize_structure(atlas_type, structure_type, brain_region_final)
    ) %>%
    dplyr::group_by(brain_region = brain_region_final, structure) %>%
    dplyr::summarise(
      mean_z_change = mean(z_change, na.rm = TRUE),
      n_idps = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(n_idps))
  
  total_sig <- nrow(idp_map)
  
  list(df = df, total_sig = total_sig)
}

plot_brain_changes_overview <- function(overview_data) {
  df <- overview_data$df
  total_sig <- overview_data$total_sig
  
  if (nrow(df) == 0) {
    stop("No data available to plot brain changes overview.")
  }
  
  # Order y by structure then by count
  df <- df %>%
    dplyr::mutate(structure = factor(structure, levels = c("brainstem","cerebellum","cortical","subcortical","white_matter","unknown"))) %>%
    dplyr::arrange(structure, dplyr::desc(n_idps), dplyr::desc(abs(mean_z_change)))
  df$brain_region <- factor(df$brain_region, levels = unique(df$brain_region))
  
  # Colors per structure
  structure_colors <- c(
    brainstem = "#7f7f7f",
    cerebellum = "#ff9900",
    cortical = "#d62728",
    subcortical = "#1f77b4",
    white_matter = "#2ca02c",
    unknown = "#8a8a8a"
  )
  
  # Dynamic x-axis limits (in percentage)
  x_vals <- df$mean_z_change * 100
  max_abs <- suppressWarnings(max(abs(x_vals), na.rm = TRUE))
  if (!is.finite(max_abs)) max_abs <- 5
  pad <- max(5, min(20, max_abs * 0.1))
  x_lim <- max(25, min(2000, max_abs + pad))
  
  # Label top regions by count
  count_threshold <- max(30, stats::quantile(df$n_idps, 0.9, na.rm = TRUE))
  df_lab <- df %>% dplyr::filter(n_idps >= count_threshold)
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = mean_z_change * 100, y = brain_region)) +
    # shaded common-support interval around 0%
    ggplot2::annotate("rect", xmin = -5, xmax = 5, ymin = -Inf, ymax = Inf, fill = "#F7FAFC", alpha = 0.7) +
    ggplot2::geom_vline(xintercept = c(-5, 5), linetype = "dashed", color = "grey70", linewidth = 0.8) +
    ggplot2::geom_vline(xintercept = 0, color = "black", linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(size = n_idps, color = structure), alpha = 0.85) +
    ggplot2::scale_color_manual(values = structure_colors, name = "Brain Structure") +
    ggplot2::scale_size(range = c(2.5, 12), breaks = c(40, 80, 120, 160), name = "Number\nof IDPs") +
    ggplot2::scale_x_continuous(
      name = "Mean Z-score Change (%)",
      limits = c(-x_lim, x_lim),
      labels = scales::label_number(accuracy = 1)
    ) +
    ggplot2::scale_y_discrete(name = "Brain Region") +
    ggplot2::labs(
      title = "Comprehensive Brain Changes Overview",
      subtitle = "Point size indicates number of significant IDPs per region",
      caption = sprintf("Based on %s significant IDPs", scales::comma(total_sig))
    ) +
    theme_jama_refined() +
    ggplot2::theme(
      legend.key.height = ggplot2::unit(12, "pt"),
      legend.key.width = ggplot2::unit(16, "pt"),
      panel.grid.major.y = ggplot2::element_blank()
    )
  
  if (nrow(df_lab) > 0) {
    p <- p + ggplot2::geom_label(
      data = df_lab,
      ggplot2::aes(label = n_idps),
      size = 3, fill = "#E0E0E0", color = "black", label.size = 0.2, alpha = 0.85,
      nudge_x = 0
    )
  }
  
  p
}

render_brain_changes_overview <- function(output_dir = getwd()) {
  overview <- build_brain_changes_overview_data()
  p <- plot_brain_changes_overview(overview)
  save_pub(
    filename = file.path(output_dir, "Brain_Changes_Overview", "Comprehensive_Brain_Changes_Overview"),
    plot = p,
    width = 12, height = 9
  )
}

## 可视化：倾向评分分布（匹配前后）
plot_ps_distribution <- function(mobj = NULL, original, matched, group_col = "group", distance_col = "distance", layout = c("horizontal","vertical"), fixed_y = FALSE) {
  log_plot_call("plot_ps_distribution")
  get_group <- function(df) factor(ifelse(df[[group_col]] %in% c(1, "1"), "Case", "Control"),
                                   levels = c("Control", "Case"))
  # 优先使用数据列中的距离，如不存在再尝试从mobj获取
  dist_before <- if (distance_col %in% names(original)) original[[distance_col]] else {
    if (!is.null(mobj) && length(mobj$distance) == nrow(original)) mobj$distance else rep(NA_real_, nrow(original))
  }
  dist_after  <- if (distance_col %in% names(matched)) matched[[distance_col]] else {
    if (!is.null(mobj) && length(mobj$distance) == nrow(matched)) mobj$distance else rep(NA_real_, nrow(matched))
  }
  ps_all <- data.frame(propensity_score = dist_before, group = get_group(original), status = "Before Matching")
  ps_mat <- data.frame(propensity_score = dist_after,  group = get_group(matched),  status = "After Matching")
  ps_df <- dplyr::bind_rows(ps_all, ps_mat) %>% dplyr::filter(!is.na(propensity_score)) %>%
    dplyr::mutate(status = factor(status, levels = c("Before Matching", "After Matching")))
  eps <- suppressWarnings(max(1e-6, min(ps_df$propensity_score[ps_df$propensity_score > 0], na.rm = TRUE) / 10))
  low <- eps
  high <- 1 - eps
  ps_df$ps_plot <- pmin(pmax(ps_df$propensity_score, low), high)
  
  # 计算共同支持区间（每个状态各自计算）
  support_df <- ps_df %>%
    dplyr::group_by(status, group) %>%
    dplyr::summarise(min_ps = min(propensity_score, na.rm = TRUE),
                     max_ps = max(propensity_score, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = group, values_from = c(min_ps, max_ps)) %>%
    dplyr::mutate(
      common_min = pmax(min_ps_Control, min_ps_Case, na.rm = TRUE),
      common_max = pmin(max_ps_Control, max_ps_Case, na.rm = TRUE)
    )
  support_df$common_min_plot <- pmin(pmax(support_df$common_min, low), high)
  support_df$common_max_plot <- pmin(pmax(support_df$common_max, low), high)
  
  lay <- match.arg(layout)
  counts <- ps_df %>% dplyr::group_by(status, group) %>% dplyr::summarise(n = dplyr::n(), .groups = "drop") %>% tidyr::pivot_wider(names_from = group, values_from = n) %>% dplyr::mutate(lbl = sprintf("Control N=%s; Case N=%s", scales::comma(Control), scales::comma(Case)))
  cap <- paste0(
    "Before: ", counts$lbl[counts$status == "Before Matching"], "; ",
    "After: ", counts$lbl[counts$status == "After Matching"],
    "; Shaded area indicates common support"
  )
  ggplot2::ggplot(ps_df, ggplot2::aes(x = ps_plot, fill = group)) +
    # 共同支持区间底色
    ggplot2::geom_rect(data = support_df,
                       ggplot2::aes(xmin = common_min_plot, xmax = common_max_plot, ymin = -Inf, ymax = Inf),
                       inherit.aes = FALSE, fill = "#F0F9FF", alpha = 0.25) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins = 40, alpha = 0.65,
                            position = "identity", color = "white", linewidth = 0.3) +
    ggplot2::geom_density(alpha = 0.35, linewidth = 1) +
    ggplot2::facet_wrap(~ status, nrow = if (lay == "horizontal") 1 else 2, ncol = if (lay == "horizontal") 2 else 1, scales = if (fixed_y) "fixed" else "free_y") +
    ggplot2::scale_fill_manual(values = c("Control" = jama_colors$control, "Case" = jama_colors$case), name = "Group") +
    ggplot2::scale_x_continuous(name = "Propensity Score (logit)", trans = scales::logit_trans(),
                                breaks = c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 0.99), limits = c(low, high)) +
    ggplot2::scale_y_continuous(name = "Density") +
    ggplot2::labs(title = "Figure 2 — Propensity Score Distribution",
                  subtitle = NULL,
                  caption = cap) +
    theme_jama_refined()
}

# 简化流程数据汇总（每组）
summarize_psm_flow <- function(original, matched, group_col = "group") {
  safe_sum <- function(x) sum(x, na.rm = TRUE)
  init_total <- nrow(original)
  init_cases <- safe_sum(original[[group_col]] == 1)
  init_controls <- safe_sum(original[[group_col]] == 0)
  match_total <- nrow(matched)
  match_cases <- safe_sum(matched[[group_col]] == 1)
  match_controls <- safe_sum(matched[[group_col]] == 0)
  data.frame(
    initial_total = init_total,
    initial_cases = init_cases,
    initial_controls = init_controls,
    matched_total = match_total,
    matched_cases = match_cases,
    matched_controls = match_controls,
    match_ratio = ifelse(match_cases > 0, round(match_controls / match_cases, 2), NA_real_)
  )
}

# IDP分类与脑图映射：支持外部字典与内置启发式
idp_mapping_cache <- NULL
load_idp_mapping_dictionary <- function(path = "idp_mapping_dictionary.csv", reload = FALSE) {
  if (!reload && !is.null(idp_mapping_cache)) return(idp_mapping_cache)
  if (!file.exists(path)) return(NULL)
  dict <- tryCatch(utils::read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(dict)) {
    # 标准化列名
    names(dict) <- tolower(names(dict))
  }
  idp_mapping_cache <<- dict
  dict
}

match_dict_category <- function(x) {
  dict <- load_idp_mapping_dictionary()
  if (is.null(dict)) return(NULL)
  for (i in seq_len(nrow(dict))) {
    pat <- dict$pattern[i]
    if (!is.na(pat) && nzchar(pat)) {
      if (grepl(pat, x, perl = TRUE, ignore.case = TRUE)) {
        catg <- dict$category[i]
        if (!is.na(catg) && nzchar(catg)) return(catg)
      }
    }
  }
  NULL
}

# English SCI publication category mapping for IDPs（先字典后启发式）
identify_idp_category_en <- function(x) {
  # 先试字典
  catg <- match_dict_category(x)
  if (!is.null(catg)) return(catg)
  # 启发式兜底
  x2 <- tolower(x)
  return(dplyr::case_when(
    stringr::str_detect(x2, "fluid\\.intelligence|reaction\\.time|trail|symbol\\.digit|pairs|match|memory|cognitive|intelligence|digits\\.remembered|errors\\.|time\\.to\\.complete|duration\\.to\\.complete|mean\\.time") ~ "Cognitive",
    stringr::str_detect(x2, "white\\.?\\s*matter\\.?\\s*hyperintensity|wmh") ~ "White Matter Hyperintensity Volume",
    stringr::str_detect(x2, "^volume\\.") ~ "Regional/Tissue Volume",
    stringr::str_detect(x2, "mean\\.thickness") ~ "Cortical Thickness",
    stringr::str_detect(x2, "grey\\.white\\.contrast") ~ "Gray–White Contrast",
    stringr::str_detect(x2, "intensity") ~ "Regional/Tissue Intensity",
    stringr::str_detect(x2, "(weighted\\.mean|mean)\\.fa") ~ "WM Tract FA",
    stringr::str_detect(x2, "(weighted\\.mean|mean)\\.md|diffusiv") ~ "WM Tract Diffusivity",
    stringr::str_detect(x2, "(weighted\\.mean|mean)\\.(l1|l2|l3)\\b") ~ "WM Tract Diffusivity",
    stringr::str_detect(x2, "(weighted\\.mean|mean)\\.mo\\b") ~ "WM Tract MO",
    stringr::str_detect(x2, "(weighted\\.mean|mean)\\.icvf") ~ "WM Tract ICVF",
    stringr::str_detect(x2, "(weighted\\.mean|mean)\\.od") ~ "WM Tract OD",
    stringr::str_detect(x2, "(weighted\\.mean|mean)\\.isovf") ~ "WM Tract ISOVF",
    stringr::str_detect(x2, "^area\\.") ~ "Cortical Surface Area",
    stringr::str_detect(x2, "rfmri.*amplitude|node.*amplitude") ~ "rfMRI Node Amplitude",
    stringr::str_detect(x2, "rfmri.*connectivity|connectivity") ~ "rfMRI Connectivity",
    TRUE ~ "Other/New IDP"
  ))
}

plot_idp_z_distribution_all_models <- function(comparisons = c("ischemic_vs_control", "mi_vs_control", "chronic_vs_control")) {
  src <- NULL
  if (exists("combined_four_model_results") && !is.null(combined_four_model_results) && nrow(combined_four_model_results) > 0) {
    src <- combined_four_model_results
  } else {
    csv_path <- combined_four_model_results_csv()
    if (file.exists(csv_path)) {
      src <- tryCatch(utils::read.csv(csv_path, stringsAsFactors = FALSE), error = function(e) NULL)
    }
  }
  if (is.null(src) || nrow(src) == 0) {
    cat("⚠️ 跳过IDP Z分布图：未找到 Combined_Four_Model_Z_Analysis_Results.csv 或数据为空。\n")
    return(invisible(NULL))
  }
  df <- clean_combined_results_df(src)
  if (!("model_comparison" %in% names(df)) || !("variable_name" %in% names(df))) {
    cat("⚠️ 跳过IDP Z分布图：缺少必要列(model_comparison/variable_name)。\n")
    return(invisible(NULL))
  }
  if (!("variable_type" %in% names(df))) df$variable_type <- NA_character_
  df <- df[df$model_comparison %in% comparisons, , drop = FALSE]
  df <- df[stringr::str_detect(as.character(df$variable_type), stringr::regex("brain\\s*_?\\s*imaging", ignore_case = TRUE)), , drop = FALSE]
  if (nrow(df) == 0) {
    cat("⚠️ 跳过IDP Z分布图：筛选后无脑影像IDP数据。\n")
    return(invisible(NULL))
  }
  z_without  <- if ("z_statistic_without_non_imaging"  %in% names(df)) suppressWarnings(as.numeric(df$z_statistic_without_non_imaging))  else rep(NA_real_, nrow(df))
  z_with     <- if ("z_statistic_with_non_imaging"     %in% names(df)) suppressWarnings(as.numeric(df$z_statistic_with_non_imaging))     else rep(NA_real_, nrow(df))
  cd_with    <- if ("cohens_d_with_non_imaging"        %in% names(df)) suppressWarnings(as.numeric(df$cohens_d_with_non_imaging))        else rep(NA_real_, nrow(df))
  cd_without <- if ("cohens_d_without_non_imaging"     %in% names(df)) suppressWarnings(as.numeric(df$cohens_d_without_non_imaging))     else rep(NA_real_, nrow(df))
  cd_raw     <- if ("cohens_d"                         %in% names(df)) suppressWarnings(as.numeric(df$cohens_d))                         else rep(NA_real_, nrow(df))
  pooled_z <- dplyr::coalesce(z_without, z_with, cd_with, cd_without, cd_raw)
  plot_df <- tibble::tibble(
    idp_category = identify_idp_category_en(as.character(df$variable_name)),
    pooled_z = pooled_z,
    comparison = as.character(df$model_comparison)
  ) %>% tidyr::drop_na(idp_category, pooled_z, comparison)
  if (!exists("idp_levels_master", inherits = TRUE) || is.null(idp_levels_master)) {
    idp_levels_master <- c(
      "Regional/Tissue Volume",
      "Cortical Thickness",
      "Cortical Surface Area",
      "Regional/Tissue Intensity",
      "Gray–White Contrast",
      "White Matter Hyperintensity Volume",
      "WM Tract FA",
      "WM Tract Diffusivity",
      "WM Tract ICVF",
      "WM Tract MO",
      "WM Tract OD",
      "WM Tract ISOVF",
      "rfMRI Node Amplitude",
      "rfMRI Connectivity",
      "Other/New IDP"
    )
  }
  idp_cat_chr <- as.character(plot_df$idp_category)
  extra_levels <- sort(setdiff(unique(idp_cat_chr), idp_levels_master))
  present_levels <- unique(idp_cat_chr)
  ordered_levels <- c(
    idp_levels_master[idp_levels_master %in% present_levels],
    sort(setdiff(present_levels, idp_levels_master))
  )
  plot_df <- plot_df %>%
    dplyr::mutate(
      idp_category = factor(idp_category, levels = ordered_levels),
      comparison = factor(comparison, levels = comparisons)
    )
  idp_x_labels <- levels(plot_df$idp_category)
  n_idp_levels <- length(idp_x_labels)
  jama_cmp_colors <- c(
    "ischemic_vs_control" = "#0072B2",
    "mi_vs_control"       = "#D55E00",
    "chronic_vs_control"  = "#009E73"
  )
  comparison_labels <- c(
    "ischemic_vs_control" = "Ischemic vs Control",
    "mi_vs_control"       = "Myocardial Infarction vs Control",
    "chronic_vs_control"  = "Chronic Ischemia vs Control"
  )
  jama_cmp_colors <- jama_cmp_colors[names(jama_cmp_colors) %in% comparisons]
  comparison_labels <- comparison_labels[names(comparison_labels) %in% comparisons]
  p_z <- ggplot2::ggplot(plot_df, ggplot2::aes(x = idp_category, y = pooled_z, color = comparison)) +
    ggplot2::geom_point(position = ggplot2::position_jitterdodge(jitter.width = 0.20, jitter.height = 0, dodge.width = 0.70), alpha = 0.85, size = 1.8) +
    ggplot2::geom_vline(xintercept = seq(1.5, max(n_idp_levels - 0.5, 1.5), by = 1), color = "#e6e6e6", linetype = "dotted", linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    ggplot2::scale_color_manual(values = jama_cmp_colors, breaks = names(comparison_labels), labels = unname(comparison_labels), name = "Comparison") +
    ggplot2::scale_x_discrete(limits = idp_x_labels, labels = stringr::str_wrap(idp_x_labels, width = 18)) +
    ggplot2::scale_y_continuous(limits = c(-6, 6), breaks = seq(-6, 6, 2)) +
    ggplot2::labs(title = "Z-Statistic Distribution of Brain Imaging IDPs (Three Comparisons)", x = "IDP Category", y = "Model Z-Statistic") +
    theme_jama_refined() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, vjust = 1))
  save_pub("IDP_Z_Distribution_All_Models", p_z, width = 14, height = 6, dpi = 300)
  invisible(p_z)
}

plot_scan_rescan_reproducibility_by_idp_class <- function() {
  in_csv <- table_out_path("Diagnostics", "Diagnostics_ischemic_vs_control_IDP_Reproducibility.csv")
  if (!file.exists(in_csv)) {
    cat("⚠️ 跳过Scan–Rescan可重复性图：未找到 Diagnostics_ischemic_vs_control_IDP_Reproducibility.csv。\n")
    return(invisible(NULL))
  }
  rep_df <- tryCatch(utils::read.csv(in_csv, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(rep_df) || nrow(rep_df) == 0) {
    cat("⚠️ 跳过Scan–Rescan可重复性图：诊断CSV为空。\n")
    return(invisible(NULL))
  }
  need_cols <- c("variable_type", "base_name", "r_cases", "r_controls")
  if (!all(need_cols %in% names(rep_df))) {
    cat("⚠️ 跳过Scan–Rescan可重复性图：诊断CSV缺少必要列。\n")
    return(invisible(NULL))
  }
  rep_df <- rep_df[rep_df$variable_type == "Brain_Imaging", , drop = FALSE]
  rep_plot_df <- dplyr::transmute(
    rep_df,
    base_name = as.character(base_name),
    idp_category = identify_idp_category_en(as.character(base_name)),
    cases = suppressWarnings(as.numeric(r_cases)),
    controls = suppressWarnings(as.numeric(r_controls))
  ) %>% tidyr::drop_na(idp_category)
  if (nrow(rep_plot_df) == 0) {
    cat("⚠️ 跳过Scan–Rescan可重复性图：无可绘制数据。\n")
    return(invisible(NULL))
  }
  long_rep <- data.frame(
    idp_category = rep(rep_plot_df$idp_category, 2),
    group = rep(c("Cases", "Controls"), each = nrow(rep_plot_df)),
    correlation = c(rep_plot_df$cases, rep_plot_df$controls),
    stringsAsFactors = FALSE
  ) %>% tidyr::drop_na(correlation)
  if (!exists("idp_levels_master", inherits = TRUE) || is.null(idp_levels_master)) {
    idp_levels_master <- c(
      "Regional/Tissue Volume",
      "Cortical Thickness",
      "Cortical Surface Area",
      "Regional/Tissue Intensity",
      "Gray–White Contrast",
      "White Matter Hyperintensity Volume",
      "WM Tract FA",
      "WM Tract Diffusivity",
      "WM Tract ICVF",
      "WM Tract OD",
      "WM Tract ISOVF",
      "rfMRI Node Amplitude",
      "rfMRI Connectivity",
      "Other/New IDP"
    )
  }
  long_rep$idp_category <- factor(long_rep$idp_category, levels = {
    uniq <- unique(as.character(long_rep$idp_category))
    lvl <- idp_levels_master[idp_levels_master %in% uniq]
    if (length(lvl) == 0) idp_levels_master else lvl
  })
  long_rep$group <- factor(long_rep$group, levels = c("Cases", "Controls"))
  grp_offsets <- c("Cases" = -0.125, "Controls" = 0.125)
  long_rep$idp_index <- as.numeric(long_rep$idp_category)
  long_rep$offset <- dplyr::coalesce(unname(grp_offsets[as.character(long_rep$group)]), 0)
  long_rep$xpos <- long_rep$idp_index + long_rep$offset
  cmp_colors <- c("Cases" = "#D55E00", "Controls" = "#0072B2")
  p_rep <- ggplot2::ggplot(long_rep, ggplot2::aes(x = xpos, y = correlation, color = group)) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.06, height = 0), alpha = 0.85, size = 1.8) +
    ggplot2::geom_vline(xintercept = seq(1.5, max(length(levels(long_rep$idp_category)) - 0.5, 1.5), by = 1), color = "#e6e6e6", linetype = "dotted", linewidth = 0.5) +
    ggplot2::scale_color_manual(values = cmp_colors, breaks = c("Cases", "Controls"), labels = c("Cases", "Controls"), name = "Group") +
    ggplot2::scale_x_continuous(breaks = seq_along(levels(long_rep$idp_category)), labels = stringr::str_wrap(levels(long_rep$idp_category), width = 18), expand = ggplot2::expansion(mult = c(0.02, 0.02))) +
    ggplot2::scale_y_continuous(limits = c(-0.2, 1.0), breaks = seq(0, 1, 0.2)) +
    ggplot2::labs(title = "Scan–Rescan Reproducibility of Brain Imaging IDPs", x = "IDP Category", y = "Pearson r (Scan 2 vs Scan 1)") +
    theme_jama_refined() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, vjust = 1))
  save_pub("Scan_Rescan_Reproducibility_By_IDP_Class", p_rep, width = 14, height = 6, dpi = 300)
  invisible(p_rep)
}

if (RUN_REPLOT_IDP_ZDIST) {
  plot_idp_z_distribution_all_models()
  cat("✓ 已重跑并保存 IDP_Z_Distribution_All_Models\n")
}
if (RUN_REPLOT_SCAN_RESCAN) {
  plot_scan_rescan_reproducibility_by_idp_class()
  cat("✓ 已重跑并保存 Scan_Rescan_Reproducibility_By_IDP_Class\n")
}

# 补充脑图属性注解（ROI分组/模态/半球/测量类型/显示名）
annotate_idp_attributes <- function(df, var_col = "variable_name") {
  if (is.null(df) || nrow(df) == 0 || !(var_col %in% names(df))) return(df)
  dict <- load_idp_mapping_dictionary()
  # 内置启发式函数
  infer_roi <- function(x) {
    x2 <- tolower(x)
    dplyr::case_when(
      stringr::str_detect(x2, "tract|fasciculus|radiation|peduncle|forceps|corticospinal|cingulum") ~ "WM Tracts",
      stringr::str_detect(x2, "hippocamp") ~ "Hippocampus",
      stringr::str_detect(x2, "amygdal") ~ "Amygdala",
      stringr::str_detect(x2, "thalam") ~ "Thalamus",
      stringr::str_detect(x2, "putamen") ~ "Putamen",
      stringr::str_detect(x2, "caudate") ~ "Caudate",
      stringr::str_detect(x2, "pallidum|globus") ~ "Pallidum",
      stringr::str_detect(x2, "insula") ~ "Insula",
      stringr::str_detect(x2, "cingulate|acc") ~ "Cingulate",
      stringr::str_detect(x2, "frontal") ~ "Frontal",
      stringr::str_detect(x2, "parietal") ~ "Parietal",
      stringr::str_detect(x2, "temporal") ~ "Temporal",
      stringr::str_detect(x2, "occipital") ~ "Occipital",
      stringr::str_detect(x2, "cerebell") ~ "Cerebellum",
      stringr::str_detect(x2, "brainstem") ~ "Brainstem",
      TRUE ~ NA_character_
    )
  }
  infer_modality <- function(x) {
    x2 <- tolower(x)
    dplyr::case_when(
      stringr::str_detect(x2, "(weighted\\.mean|mean)\\.fa") ~ "dMRI-FA",
      stringr::str_detect(x2, "(weighted\\.mean|mean)\\.md|diffusiv") ~ "dMRI-MD",
      stringr::str_detect(x2, "(weighted\\.mean|mean)\\.(l1|l2|l3)\\b") ~ "dMRI-L1/L2/L3",
      stringr::str_detect(x2, "(weighted\\.mean|mean)\\.mo\\b") ~ "dMRI-MO",
      stringr::str_detect(x2, "(weighted\\.mean|mean)\\.icvf") ~ "dMRI-ICVF",
      stringr::str_detect(x2, "(weighted\\.mean|mean)\\.od") ~ "dMRI-OD",
      stringr::str_detect(x2, "(weighted\\.mean|mean)\\.isovf") ~ "dMRI-ISOVF",
      stringr::str_detect(x2, "mean\\.thickness|^area\\.|^volume\\.|intensity|contrast") ~ "T1-Structural",
      stringr::str_detect(x2, "rfmri") ~ "rfMRI",
      TRUE ~ NA_character_
    )
  }
  infer_measure <- function(x) {
    x2 <- tolower(x)
    dplyr::case_when(
      stringr::str_detect(x2, "fa|md|l1|l2|l3|mo|icvf|od|isovf") ~ "WM Tract",
      stringr::str_detect(x2, "thickness|area|contrast") ~ "Cortical",
      stringr::str_detect(x2, "volume|intensity") ~ "Subcortical/Regional",
      stringr::str_detect(x2, "rfmri") ~ "Functional",
      TRUE ~ NA_character_
    )
  }
  infer_hemi <- function(x) {
    x2 <- tolower(x)
    dplyr::case_when(
      stringr::str_detect(x2, "(^|\\W)left(\\W|$)|\\blh\\b|left hemisphere") ~ "Left",
      stringr::str_detect(x2, "(^|\\W)right(\\W|$)|\\brh\\b|right hemisphere") ~ "Right",
      TRUE ~ NA_character_
    )
  }
  pretty_name <- function(x) {
    s <- gsub("\\.", " ", x)
    s <- gsub("_", " ", s)
    s <- trimws(s)
    paste0(toupper(substr(s,1,1)), substr(s,2,nchar(s)))
  }
  # 字典辅助：若有匹配行，覆盖启发式
  apply_dict <- function(x, field) {
    if (is.null(dict)) return(NA_character_)
    x2 <- tolower(x)
    for (i in seq_len(nrow(dict))) {
      pat <- dict$pattern[i]
      if (!is.na(pat) && nzchar(pat) && grepl(pat, x2, perl = TRUE)) {
        val <- dict[[field]][i]
        if (!is.null(val) && !is.na(val) && nzchar(val)) return(val)
      }
    }
    NA_character_
  }
  v <- df[[var_col]]
  df$idp_category <- vapply(v, identify_idp_category_en, FUN.VALUE = character(1))
  df$roi_group    <- vapply(v, function(x) dplyr::coalesce(apply_dict(x, "roi_group"), infer_roi(x)), FUN.VALUE = character(1))
  df$modality     <- vapply(v, function(x) dplyr::coalesce(apply_dict(x, "modality"), infer_modality(x)), FUN.VALUE = character(1))
  df$measure_type <- vapply(v, function(x) dplyr::coalesce(apply_dict(x, "measure_type"), infer_measure(x)), FUN.VALUE = character(1))
  df$hemisphere   <- vapply(v, function(x) dplyr::coalesce(apply_dict(x, "hemisphere"), infer_hemi(x)), FUN.VALUE = character(1))
  df$display_name <- vapply(v, function(x) dplyr::coalesce(apply_dict(x, "display_name"), pretty_name(x)), FUN.VALUE = character(1))
  df
}

# 流程图（每组）：两节点（初始/匹配后）+ 箭头，出版风格
make_psm_flowchart <- function(flow_df, cmp_label, layout = c("horizontal","vertical")) {
  # 计算排除人数（未匹配/超出 caliper）
  excl_total    <- max(flow_df$initial_total    - flow_df$matched_total,    0)
  excl_cases    <- max(flow_df$initial_cases    - flow_df$matched_cases,    0)
  excl_controls <- max(flow_df$initial_controls - flow_df$matched_controls, 0)
  
  # 多节点标签（出版风格）
  lab_init <- sprintf(
    "Eligible cohort — %s\nTotal: %s\nCases: %s\nControls: %s",
    cmp_label,
    scales::comma(flow_df$initial_total),
    scales::comma(flow_df$initial_cases),
    scales::comma(flow_df$initial_controls)
  )
  lab_excl <- sprintf(
    "Excluded (no match / outside caliper)\nTotal: %s\nCases: %s\nControls: %s",
    scales::comma(excl_total),
    scales::comma(excl_cases),
    scales::comma(excl_controls)
  )
  lab_match_cases <- sprintf(
    "Matched Cases\nN: %s",
    scales::comma(flow_df$matched_cases)
  )
  lab_match_controls <- sprintf(
    "Matched Controls\nN: %s",
    scales::comma(flow_df$matched_controls)
  )
  lab_analysis <- sprintf(
    "Analysis & Visualizations\nLove Plot, PS Dist., Flowchart\nSaved: current directory"
  )
  
  lay <- match.arg(layout)
  if (lay == "horizontal") {
    nodes <- tibble::tibble(
      x = c(0, 0, 13, 13, 24),
      y = c(2, -6, 2, -2, 0),
      label = c(lab_init, lab_excl, lab_match_cases, lab_match_controls, lab_analysis)
    )
    node_w <- 6; node_h <- 3
    anchor_bottom <- function(x,y) c(x, y - node_h/2)
    anchor_top    <- function(x,y) c(x, y + node_h/2)
    anchor_left   <- function(x,y) c(x - node_w/2, y)
    anchor_right  <- function(x,y) c(x + node_w/2, y)
    a_init_bot <- anchor_bottom(nodes$x[1], nodes$y[1])
    a_excl_left <- anchor_left(nodes$x[2], nodes$y[2])
    a_case_top  <- anchor_top(nodes$x[3], nodes$y[3])
    a_ctrl_top  <- anchor_top(nodes$x[4], nodes$y[4])
    a_ana_top   <- anchor_top(nodes$x[5], nodes$y[5])
    a_init_right<- anchor_right(nodes$x[1], nodes$y[1])
    a_case_bot  <- anchor_bottom(nodes$x[3], nodes$y[3])
    a_ctrl_bot  <- anchor_bottom(nodes$x[4], nodes$y[4])
    p <- ggplot2::ggplot(nodes) +
      ggplot2::geom_label(ggplot2::aes(x = x, y = y, label = label),
                          size = 4.2, lineheight = 1.05, hjust = 0.5,
                          label.size = 0.6, label.r = grid::unit(6, "pt"), label.padding = grid::unit(6, "pt"),
                          fill = "#F8FAFC", color = "#374151") +
      ggplot2::geom_segment(x = a_init_bot[1], y = a_init_bot[2], xend = a_case_top[1], yend = a_case_top[2],
                            arrow = grid::arrow(length = grid::unit(0.25, "cm")), color = "#374151", linewidth = 0.9) +
      ggplot2::geom_segment(x = a_init_bot[1], y = a_init_bot[2], xend = a_ctrl_top[1], yend = a_ctrl_top[2],
                            arrow = grid::arrow(length = grid::unit(0.25, "cm")), color = "#374151", linewidth = 0.9) +
      ggplot2::geom_segment(x = a_init_right[1], y = a_init_right[2], xend = a_excl_left[1], yend = a_excl_left[2],
                            arrow = grid::arrow(length = grid::unit(0.25, "cm")), color = "#374151", linewidth = 0.9) +
      ggplot2::geom_segment(x = a_case_bot[1], y = a_case_bot[2], xend = a_ana_top[1], yend = a_ana_top[2],
                            arrow = grid::arrow(length = grid::unit(0.25, "cm")), color = "#374151", linewidth = 0.9) +
      ggplot2::geom_segment(x = a_ctrl_bot[1], y = a_ctrl_bot[2], xend = a_ana_top[1], yend = a_ana_top[2],
                            arrow = grid::arrow(length = grid::unit(0.25, "cm")), color = "#374151", linewidth = 0.9) +
      ggplot2::labs(title = paste("PSM Study Flow —", cmp_label)) +
      theme_pub() +
      ggplot2::coord_equal() +
      ggplot2::xlim(-4, 28) + ggplot2::ylim(-10, 8) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0, face = "bold"),
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank()
      )
  } else {
    nodes <- tibble::tibble(
      x = c(0, 14, -6, 6, 0),
      y = c(12, 12, 0, 0, -10),
      label = c(lab_init, lab_excl, lab_match_cases, lab_match_controls, lab_analysis)
    )
    node_w <- 8; node_h <- 4
    anchor_bottom <- function(x,y) c(x, y - node_h/2)
    anchor_top    <- function(x,y) c(x, y + node_h/2)
    anchor_left   <- function(x,y) c(x - node_w/2, y)
    anchor_right  <- function(x,y) c(x + node_w/2, y)
    a_init_bot <- anchor_bottom(nodes$x[1], nodes$y[1])
    a_init_right<- anchor_right(nodes$x[1], nodes$y[1])
    a_excl_left <- anchor_left(nodes$x[2], nodes$y[2])
    a_case_top  <- anchor_top(nodes$x[3], nodes$y[3])
    a_ctrl_top  <- anchor_top(nodes$x[4], nodes$y[4])
    a_case_bot  <- anchor_bottom(nodes$x[3], nodes$y[3])
    a_ctrl_bot  <- anchor_bottom(nodes$x[4], nodes$y[4])
    a_ana_top   <- anchor_top(nodes$x[5], nodes$y[5])
    p <- ggplot2::ggplot(nodes) +
      ggplot2::geom_label(ggplot2::aes(x = x, y = y, label = label),
                          size = 4.5, lineheight = 1.06, hjust = 0.5,
                          label.size = 0.6, label.r = grid::unit(6, "pt"), label.padding = grid::unit(6, "pt"),
                          fill = "#F8FAFC", color = "#374151") +
      ggplot2::geom_segment(x = a_init_bot[1], y = a_init_bot[2], xend = a_case_top[1], yend = a_case_top[2],
                            arrow = grid::arrow(length = grid::unit(0.25, "cm")), color = "#374151", linewidth = 0.9) +
      ggplot2::geom_segment(x = a_init_bot[1], y = a_init_bot[2], xend = a_ctrl_top[1], yend = a_ctrl_top[2],
                            arrow = grid::arrow(length = grid::unit(0.25, "cm")), color = "#374151", linewidth = 0.9) +
      ggplot2::geom_segment(x = a_init_right[1], y = a_init_right[2], xend = a_excl_left[1], yend = a_excl_left[2],
                            arrow = grid::arrow(length = grid::unit(0.25, "cm")), color = "#374151", linewidth = 0.9) +
      ggplot2::geom_segment(x = a_case_bot[1], y = a_case_bot[2], xend = a_ana_top[1], yend = a_ana_top[2],
                            arrow = grid::arrow(length = grid::unit(0.25, "cm")), color = "#374151", linewidth = 0.9) +
      ggplot2::geom_segment(x = a_ctrl_bot[1], y = a_ctrl_bot[2], xend = a_ana_top[1], yend = a_ana_top[2],
                            arrow = grid::arrow(length = grid::unit(0.25, "cm")), color = "#374151", linewidth = 0.9) +
      ggplot2::labs(title = paste("PSM Study Flow —", cmp_label)) +
      theme_pub() +
      ggplot2::coord_equal() +
      ggplot2::xlim(-14, 18) + ggplot2::ylim(-14, 16) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0, face = "bold"),
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank()
      )
  }
  p
}

# ----- 网络/分区解析：从IDP变量名中提取功能网络或分区标签 -----
parse_idp_network <- function(variable_name) {
  v <- tolower(variable_name)
  dplyr::case_when(
    stringr::str_detect(v, "default|dmn|posterior\\.cingulate|precuneus") ~ "Default Mode",
    stringr::str_detect(v, "frontopariet|fpn|executive|dlpfc|middle\\.frontal|inferior\\.parietal") ~ "Frontoparietal",
    stringr::str_detect(v, "salien|sn|insula|anterior\\.cingulate") ~ "Salience",
    stringr::str_detect(v, "dorsal\\.?attention|dan|superior\\.parietal|intraparietal") ~ "Dorsal Attention",
    stringr::str_detect(v, "ventral\\.?attention|van|temporo\\.parietal|supramarginal") ~ "Ventral Attention",
    stringr::str_detect(v, "visual|occipital|calcarine|cuneus|lingual") ~ "Visual",
    stringr::str_detect(v, "somato|sensorimot|precentral|postcentral|rolandic") ~ "Somatomotor",
    stringr::str_detect(v, "limbic|amygdala|hippocamp|parahip|accumbens|olfactory") ~ "Limbic",
    stringr::str_detect(v, "language|superior\\.temporal|middle\\.temporal|parsopercularis|parsorbitalis|parstriangularis") ~ "Language",
    stringr::str_detect(v, "cerebell") ~ "Cerebellar",
    stringr::str_detect(v, "thalam|caudate|putamen|pallid|subthalam|brainstem") ~ "Subcortical",
    TRUE ~ "Unknown"
  )
}

network_palette <- c(
  "Default Mode" = "#1f77b4",
  "Frontoparietal" = "#ff7f0e",
  "Salience" = "#2ca02c",
  "Dorsal Attention" = "#9467bd",
  "Ventral Attention" = "#8c564b",
  "Visual" = "#e377c2",
  "Somatomotor" = "#7f7f7f",
  "Limbic" = "#bcbd22",
  "Language" = "#17becf",
  "Cerebellar" = "#d62728",
  "Subcortical" = "#6a3d9a",
  "Unknown" = "#aaaaaa"
)

# 对齐到 ggsegYeo2011::yeo7 atlas 的标准网络标签
align_to_yeo7_region <- function(label) {
  lbl <- trimws(label)
  dplyr::case_when(
    lbl %in% c("Default Mode", "Default") ~ "Default",
    lbl %in% c("Frontoparietal", "Control", "FPN", "Frontoparietal Control") ~ "Frontoparietal",
    lbl %in% c("Salience", "Ventral Attention", "SN", "Van") ~ "Ventral Attention",
    lbl %in% c("Dorsal Attention", "Dan") ~ "Dorsal Attention",
    lbl %in% c("Visual", "Vis") ~ "Visual",
    lbl %in% c("Somatomotor", "SomMot", "Sensorimotor") ~ "Somatomotor",
    lbl %in% c("Limbic") ~ "Limbic",
    TRUE ~ NA_character_
  )
}

select_idp_cols <- function(df) {
  cols <- names(df)
  idx <- stringr::str_detect(
    cols,
    "(?i)^(Area|Volume|Mean\\.thickness|Grey\\.white\\.contrast|Weighted\\.mean\\.(FA|MD|L1|L2|L3|ISOVF|ICVF|OD)|Mean\\.(FA|MD|L1|L2|L3|ISOVF|ICVF|OD))"
  )
  cols[idx]
}

rubin_pool <- function(estimates, ses) {
  m <- length(estimates)
  Qbar <- mean(estimates, na.rm = TRUE)
  Ubar <- mean(ses^2, na.rm = TRUE)
  B <- stats::var(estimates, na.rm = TRUE)
  Tvar <- Ubar + (1 + 1/m) * B
  z <- Qbar / sqrt(Tvar)
  list(coef = Qbar, se = sqrt(Tvar), z = z)
}


# 快速插补：数值型用中位数、分类型用众数
fast_impute_median <- function(df, imp_vars) {
  dat <- df[, imp_vars, drop = FALSE]
  for (nm in names(dat)) {
    v <- dat[[nm]]
    if (is.numeric(v) || is.integer(v)) {
      med <- suppressWarnings(stats::median(v, na.rm = TRUE))
      if (is.finite(med)) {
        v[is.na(v)] <- med
      }
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
    } else {
      dat[[nm]] <- v
    }
  }
  dat
}

build_comparison_df <- function(df, comparison) {
  if (!(comparison %in% c("ischemic_vs_control", "mi_vs_control", "chronic_vs_control"))) {
    stop(sprintf("Unsupported comparison: %s", comparison))
  }
  case_levels <- switch(
    comparison,
    "ischemic_vs_control" = c("MI Only", "Chronic Only", "Both MI & Chronic"),
    "mi_vs_control" = c("MI Only", "Both MI & Chronic"),
    "chronic_vs_control" = c("Chronic Only")
  )
  control_levels <- c("Control")
  df %>% 
    dplyr::filter(disease_subtype %in% c(case_levels, control_levels)) %>%
    dplyr::mutate(group = if_else(disease_subtype %in% case_levels, 1L, 0L))
}

  run_psm_mi <- function(dat, covariates, idp_cols, comparison, m = 5, seed = 2025,
                       impute_method = getOption("IMPUTE_METHOD", "mice"),
                       mice_maxit = getOption("MI_MAXIT", 5),
                       mice_print = getOption("MI_PRINT", FALSE),
                       psm_ratio = 2,
                       psm_caliper = 0.2) {
  set.seed(seed)
  # 多重插补
  imp_vars <- unique(c("group", covariates, idp_cols, if ("eid" %in% names(dat)) "eid" else NULL))
  if (identical(impute_method, "mice")) {
    mids <- mice::mice(dat[, imp_vars, drop = FALSE], m = m, seed = seed, printFlag = mice_print, maxit = mice_maxit)
  } else {
    mids <- NULL
    m <- 1
  }
  
  pooled_idp <- list()
  baseline_tables <- list()
  matched_k1 <- NULL
  mobj_k1 <- NULL
  original_k1 <- NULL
  
  for (k in 1:m) {
    if (!is.null(mids)) {
      dk <- mice::complete(mids, k)
    } else {
      dk <- fast_impute_median(dat, imp_vars)
    }
    
    # 倾向评分匹配（仅使用用户指定的核心协变量）
    covs_all <- covariates[covariates %in% names(dk)]
    if (length(covs_all) == 0) stop("PSM协变量在数据中不存在")
    psm_fit <- psm_match_standard(
      data = dk,
      treatment_col = "group",
      covariates = covs_all,
      method = "nearest",
      distance = "glm",
      ratio = psm_ratio,
      caliper = psm_caliper,
      replace = FALSE,
      estimand = "ATT",
      seed = seed
    )
    mobj <- psm_fit$mobj
    matched <- psm_fit$matched_data
    if (k == 1) {
      matched_k1 <- matched
      mobj_k1 <- mobj
      original_k1 <- dk
    }
    
    # Baseline表（匹配后）
    tbl1 <- tableone::CreateTableOne(vars = covariates, strata = "group", data = matched, test = FALSE)
    baseline_tables[[k]] <- as.data.frame(print(tbl1, printToggle = FALSE))
    
    # 每个IDP的Z值（匹配后）
    idp_res_k <- purrr::map_df(idp_cols, function(v) {
      vv <- matched[[v]]
      if (!is.numeric(vv)) vv <- suppressWarnings(as.numeric(vv))
      if (all(is.na(vv))) return(NULL)
      # 标准化IDP以便可比
      vv <- scale(vv)
      df_fit <- data.frame(y = as.numeric(vv), group = matched$group)
      fit <- stats::lm(y ~ group, data = df_fit)
      tidy <- broom::tidy(fit)
      eff <- dplyr::filter(tidy, term == "group") %>% dplyr::slice(1)
      data.frame(
        idp = v,
        estimate = eff$estimate,
        std.error = eff$std.error,
        statistic = eff$statistic,
        p.value = eff$p.value
      )
    })
    idp_res_k$imp <- k
    pooled_idp[[k]] <- idp_res_k
  }
  
  idp_all <- if (length(pooled_idp) > 0) dplyr::bind_rows(pooled_idp) else tibble::tibble()
  # Rubin合并
  pooled <- idp_all %>%
    dplyr::group_by(idp) %>%
    dplyr::summarise(
      m = dplyr::n(),
      pooled_coef = rubin_pool(estimate, std.error)$coef,
      pooled_se   = rubin_pool(estimate, std.error)$se,
      pooled_z    = rubin_pool(estimate, std.error)$z,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      comparison = comparison,
      idp_category = identify_idp_category_en(idp)
    )
  
  # 合并baseline表（取第一个插补的展示；若需要也可平均）
  baseline_merged <- if (length(baseline_tables) > 0) dplyr::bind_rows(baseline_tables, .id = "imp") else tibble::tibble()
  
  list(pooled_results = pooled,
       baseline_table = baseline_merged,
       matched_example = matched_k1,
       mobj_example = mobj_k1,
       original_k1 = original_k1)
}

# 工具函数：纵向ΔIDP变化率计算（自动配对 Instance.2 vs Instance.3）
# 说明：支持列名中在 Instance.* 之后仍带有 "...Array.*" 等后缀的情况；
#      根名称以移除 "...Instance.[2|3]" 后的剩余部分为准，从而实现按数组位等精细配对。
compute_delta_rate_long <- function(dat, target_roots = NULL, pairs_df = NULL) {
  to_numeric_safely <- function(v) {
    if (is.factor(v)) v <- as.character(v)
    if (is.logical(v)) return(suppressWarnings(as.numeric(v)))
    if (is.numeric(v)) return(suppressWarnings(as.numeric(v)))
    v <- trimws(as.character(v))
    v[v %in% c("", "NA", "NaN", "NULL", "null")] <- NA
    suppressWarnings(as.numeric(v))
  }
  winsorize <- function(x, frac) {
    if (is.null(frac) || !is.finite(frac) || frac <= 0) return(x)
    xv <- x[is.finite(x)]
    if (length(xv) == 0) return(x)
    qs <- stats::quantile(xv, probs = c(frac, 1 - frac), na.rm = TRUE, type = 7)
    pmin(pmax(x, qs[1]), qs[2])
  }
  rate_method <- getOption("DELTA_RATE_METHOD", "baseline")
  winsor_frac <- getOption("DELTA_RATE_WINSOR_FRAC", 0)
  if (!is.null(pairs_df) && nrow(pairs_df) > 0) {
    req <- c("idp1_variable","idp2_variable")
    if (all(req %in% names(pairs_df))) {
      pairs_df <- pairs_df[!is.na(pairs_df$idp1_variable) & !is.na(pairs_df$idp2_variable), ]
      pairs_df <- pairs_df[pairs_df$idp1_variable %in% names(dat) & pairs_df$idp2_variable %in% names(dat), ]
      if (nrow(pairs_df) == 0) return(tibble::tibble())
      root_from_name <- function(x) {
        gsub("(\\.{1,}|_)Instance\\.?[23]([._].*)?$", "", x)
      }
      out_list <- lapply(seq_len(nrow(pairs_df)), function(i) {
        c1 <- pairs_df$idp1_variable[i]
        c2 <- pairs_df$idp2_variable[i]
        x1 <- to_numeric_safely(dat[[c1]])
        x2 <- to_numeric_safely(dat[[c2]])
        id_vec <- if ("eid" %in% names(dat)) dat$eid else rep(NA, length(x1))
        if (identical(rate_method, "baseline")) {
          diffv <- suppressWarnings(x2 - x1)
          if ("group" %in% names(dat)) {
            g <- dat$group
            ctrl <- g == 0
            case <- g == 1
            idx_ctrl <- ctrl & is.finite(x1) & is.finite(diffv)
            idx_case <- case & is.finite(x1) & is.finite(diffv)
            mean_ctrl <- base::mean(x1[idx_ctrl], na.rm = TRUE)
            mean_case <- base::mean(x1[idx_case], na.rm = TRUE)
            denom <- ifelse(ctrl, abs(mean_ctrl), ifelse(case, abs(mean_case), NA_real_))
            denom[!is.finite(denom) | denom <= .Machine$double.eps] <- NA_real_
            dr <- ifelse(is.finite(diffv), diffv / denom, NA_real_)
          } else {
            denom <- ifelse(is.finite(x1) & abs(x1) > .Machine$double.eps, abs(x1), NA_real_)
            dr <- ifelse(is.finite(diffv), diffv / denom, NA_real_)
          }
        } else if (identical(rate_method, "symmetric")) {
          denom <- ifelse(is.finite(x1) & is.finite(x2), (abs(x1) + abs(x2)) / 2, NA_real_)
          diffv <- suppressWarnings(x2 - x1)
          dr <- ifelse(is.finite(diffv), diffv / denom, NA_real_)
        } else if (identical(rate_method, "log_ratio")) {
          dr <- ifelse(is.finite(x1) & is.finite(x2) & x1 > 0 & x2 > 0, log(x2) - log(x1), NA_real_)
        } else {
          denom <- ifelse(is.finite(x1) & abs(x1) > .Machine$double.eps, abs(x1), NA_real_)
          diffv <- suppressWarnings(x2 - x1)
          dr <- ifelse(is.finite(diffv), diffv / denom, NA_real_)
        }
        dr <- winsorize(dr, winsor_frac)
        tibble::tibble(eid = id_vec, idp_root = root_from_name(c1), idp1 = x1, idp2 = x2, delta_rate = dr)
      })
      non_null <- out_list[!vapply(out_list, is.null, logical(1))]
      return(if (length(non_null) > 0) dplyr::bind_rows(non_null) else tibble::tibble())
    }
  }
  cols <- names(dat)
  # 匹配任何包含 "...Instance.2"/"...Instance.3" 的列（可以在末尾或继续跟随 "...Array.*" 等后缀）
  base_cols   <- cols[stringr::str_detect(cols, "\\.{3,}Instance\\.2(\\.{1,}[A-Za-z0-9_.]+|$)")]
  follow_cols <- cols[stringr::str_detect(cols, "\\.{3,}Instance\\.3(\\.{1,}[A-Za-z0-9_.]+|$)")]

  # 根名称：移除一次 "...Instance.[2|3]"（不要求在字符串末尾）
  # 根名称：从 Instance.2/3 开始剥离到行尾，避免将 .x/.y 或 Array.* 纳入根
  root_fun <- function(x) gsub("\\.{3,}Instance\\.(2|3).*", "", x)
  base_roots   <- root_fun(base_cols)
  follow_roots <- root_fun(follow_cols)
  roots <- intersect(base_roots, follow_roots)
  if (!is.null(target_roots)) {
    target_clean <- unique(sub("\\.{1,}.*$", "", target_roots))
    roots <- intersect(roots, target_clean)
  }
  # 若无可配对的根名称，直接返回空表，避免空列表 rbind 报错
  if (length(roots) == 0) return(tibble::tibble())

  # 若提供目标根集合，则进一步筛选
  if (!is.null(target_roots)) {
    target_clean <- unique(sub("\\.{1,}.*$", "", target_roots))
    roots <- intersect(roots, target_clean)
  }
  if (length(roots) == 0) return(tibble::tibble())

  safe_first <- function(v) if (length(v) > 0) v[1] else NA_character_

  out_list <- purrr::map(roots, function(r) {
    # 根名称可能含有特殊字符，先转义
    esc_r <- stringr::str_replace_all(r, "([\\[\\]{}()*+?.\\\\^$|])", "\\\\$1")
    # 按根名称+Instance位配对（允许后续继续出现 "..." 后缀）
    pattern_base   <- paste0("^", esc_r, "\\.{3,}Instance\\.2(\\.{1,}[A-Za-z0-9_.]+|$)")
    pattern_follow <- paste0("^", esc_r, "\\.{3,}Instance\\.3(\\.{1,}[A-Za-z0-9_.]+|$)")
    c1 <- safe_first(base_cols[stringr::str_detect(base_cols, pattern_base)])
    c2 <- safe_first(follow_cols[stringr::str_detect(follow_cols, pattern_follow)])
    if (is.na(c1) || is.na(c2)) return(NULL)

    x1 <- to_numeric_safely(dat[[c1]])
    x2 <- to_numeric_safely(dat[[c2]])
    id_vec <- if ("eid" %in% names(dat)) dat$eid else rep(NA, length(x1))

    if (identical(rate_method, "baseline")) {
      diffv <- suppressWarnings(x2 - x1)
      if ("group" %in% names(dat)) {
        g <- dat$group
        ctrl <- g == 0
        case <- g == 1
        idx_ctrl <- ctrl & is.finite(x1) & is.finite(diffv)
        idx_case <- case & is.finite(x1) & is.finite(diffv)
        mean_ctrl <- base::mean(x1[idx_ctrl], na.rm = TRUE)
        mean_case <- base::mean(x1[idx_case], na.rm = TRUE)
        denom <- ifelse(ctrl, abs(mean_ctrl), ifelse(case, abs(mean_case), NA_real_))
        denom[!is.finite(denom) | denom <= .Machine$double.eps] <- NA_real_
        dr <- ifelse(is.finite(diffv), diffv / denom, NA_real_)
      } else {
        denom <- ifelse(is.finite(x1) & abs(x1) > .Machine$double.eps, abs(x1), NA_real_)
        dr <- ifelse(is.finite(diffv), diffv / denom, NA_real_)
      }
    } else if (identical(rate_method, "symmetric")) {
      denom <- ifelse(is.finite(x1) & is.finite(x2), (abs(x1) + abs(x2)) / 2, NA_real_)
      diffv <- suppressWarnings(x2 - x1)
      dr <- ifelse(is.finite(diffv), diffv / denom, NA_real_)
    } else if (identical(rate_method, "log_ratio")) {
      dr <- ifelse(is.finite(x1) & is.finite(x2) & x1 > 0 & x2 > 0, log(x2) - log(x1), NA_real_)
    } else {
      denom <- ifelse(is.finite(x1) & abs(x1) > .Machine$double.eps, abs(x1), NA_real_)
      diffv <- suppressWarnings(x2 - x1)
      dr <- ifelse(is.finite(diffv), diffv / denom, NA_real_)
    }
    dr <- winsorize(dr, winsor_frac)

  tibble::tibble(eid = id_vec, idp_root = r, idp1 = x1, idp2 = x2, delta_rate = dr)
  })

  non_null <- out_list[!vapply(out_list, is.null, logical(1))]
  if (length(non_null) > 0) dplyr::bind_rows(non_null) else tibble::tibble()
}

# 并行版本：按根名称（IDP配对）并行计算纵向Δ变化率
# 使用 future + furrr 并行映射；在缺包或不可用时自动回退为顺序。
compute_delta_rate_long_parallel <- function(dat, target_roots = NULL, workers = NULL, pairs_df = NULL) {
  to_numeric_safely <- function(v) {
    if (is.factor(v)) v <- as.character(v)
    if (is.logical(v)) return(suppressWarnings(as.numeric(v)))
    if (is.numeric(v)) return(suppressWarnings(as.numeric(v)))
    v <- trimws(as.character(v))
    v[v %in% c("", "NA", "NaN", "NULL", "null")] <- NA
    suppressWarnings(as.numeric(v))
  }
  winsorize <- function(x, frac) {
    if (is.null(frac) || !is.finite(frac) || frac <= 0) return(x)
    xv <- x[is.finite(x)]
    if (length(xv) == 0) return(x)
    qs <- stats::quantile(xv, probs = c(frac, 1 - frac), na.rm = TRUE, type = 7)
    pmin(pmax(x, qs[1]), qs[2])
  }
  rate_method <- getOption("DELTA_RATE_METHOD", "baseline")
  winsor_frac <- getOption("DELTA_RATE_WINSOR_FRAC", 0)
  if (!is.null(pairs_df) && nrow(pairs_df) > 0) {
    return(compute_delta_rate_long(dat, target_roots = target_roots, pairs_df = pairs_df))
  }
  cols <- names(dat)
  base_cols   <- cols[stringr::str_detect(cols, "\\.{3,}Instance\\.2(\\.{1,}[A-Za-z0-9_.]+|$)")]
  follow_cols <- cols[stringr::str_detect(cols, "\\.{3,}Instance\\.3(\\.{1,}[A-Za-z0-9_.]+|$)")]

  root_fun <- function(x) gsub("\\.{3,}Instance\\.(2|3).*", "", x)
  base_roots   <- root_fun(base_cols)
  follow_roots <- root_fun(follow_cols)
  roots <- intersect(base_roots, follow_roots)
  if (!is.null(target_roots)) {
    target_clean <- unique(sub("\\.{1,}.*$", "", target_roots))
    roots <- intersect(roots, target_clean)
  }
  if (length(roots) == 0) return(tibble::tibble())

  safe_first <- function(v) if (length(v) > 0) v[1] else NA_character_

  has_future <- suppressWarnings(requireNamespace("future", quietly = TRUE))
  has_furrr  <- suppressWarnings(requireNamespace("furrr", quietly = TRUE))
  if (!(has_future && has_furrr)) {
    warning("未发现 future/furrr 包，compute_delta_rate_long_parallel 改为顺序执行。")
    return(compute_delta_rate_long(dat, target_roots = target_roots, pairs_df = NULL))
  }

  old_plan <- future::plan()
on.exit({ if (exists("old_plan") && !is.null(old_plan)) future::plan(old_plan) }, add = TRUE)
  if (is.null(workers)) {
    workers <- max(1L, (if (suppressWarnings(requireNamespace("parallel", quietly = TRUE))) parallel::detectCores() else 2L) - 1L)
  }
  future::plan(future::multisession, workers = workers)

  worker_fun <- function(r) {
    esc_r <- stringr::str_replace_all(r, "([\\[\\]{}()*+?.\\^$|])", "\\\\$1")
    pattern_base   <- paste0("^", esc_r, "\\.{3,}Instance\\.2(\\.{1,}[A-Za-z0-9_.]+|$)")
    pattern_follow <- paste0("^", esc_r, "\\.{3,}Instance\\.3(\\.{1,}[A-Za-z0-9_.]+|$)")
    c1 <- safe_first(base_cols[stringr::str_detect(base_cols, pattern_base)])
    c2 <- safe_first(follow_cols[stringr::str_detect(follow_cols, pattern_follow)])
    if (is.na(c1) || is.na(c2)) return(NULL)

    x1 <- to_numeric_safely(dat[[c1]])
    x2 <- to_numeric_safely(dat[[c2]])
    id_vec <- if ("eid" %in% names(dat)) dat$eid else rep(NA, length(x1))
    if (identical(rate_method, "baseline")) {
      diffv <- suppressWarnings(x2 - x1)
      if ("group" %in% names(dat)) {
        g <- dat$group
        ctrl <- g == 0
        case <- g == 1
        idx_ctrl <- ctrl & is.finite(x1) & is.finite(diffv)
        idx_case <- case & is.finite(x1) & is.finite(diffv)
        mean_ctrl <- base::mean(x1[idx_ctrl], na.rm = TRUE)
        mean_case <- base::mean(x1[idx_case], na.rm = TRUE)
        denom <- ifelse(ctrl, abs(mean_ctrl), ifelse(case, abs(mean_case), NA_real_))
        denom[!is.finite(denom) | denom <= .Machine$double.eps] <- NA_real_
        dr <- ifelse(is.finite(diffv), diffv / denom, NA_real_)
      } else {
        denom <- ifelse(is.finite(x1) & abs(x1) > .Machine$double.eps, abs(x1), NA_real_)
        dr <- ifelse(is.finite(diffv), diffv / denom, NA_real_)
      }
    } else if (identical(rate_method, "symmetric")) {
      denom <- ifelse(is.finite(x1) & is.finite(x2), (abs(x1) + abs(x2)) / 2, NA_real_)
      diffv <- suppressWarnings(x2 - x1)
      dr <- ifelse(is.finite(diffv), diffv / denom, NA_real_)
    } else if (identical(rate_method, "log_ratio")) {
      dr <- ifelse(is.finite(x1) & is.finite(x2) & x1 > 0 & x2 > 0, log(x2) - log(x1), NA_real_)
    } else {
      denom <- ifelse(is.finite(x1) & abs(x1) > .Machine$double.eps, abs(x1), NA_real_)
      diffv <- suppressWarnings(x2 - x1)
      dr <- ifelse(is.finite(diffv), diffv / denom, NA_real_)
    }
    dr <- winsorize(dr, winsor_frac)
    tibble::tibble(eid = id_vec, idp_root = r, idp1 = x1, idp2 = x2, delta_rate = dr)
  }

  out_list <- furrr::future_map(roots, worker_fun, .options = furrr::furrr_options(seed = NULL, packages = c("dplyr","stringr","tibble")))
  non_null <- out_list[!vapply(out_list, is.null, logical(1))]
  if (length(non_null) > 0) dplyr::bind_rows(non_null) else tibble::tibble()
}

compute_delta_rate_long_summary <- function(dat, target_roots = NULL, pairs_df = NULL) {
  res <- compute_delta_rate_long(dat, target_roots = target_roots, pairs_df = pairs_df)
  if (is.null(res) || nrow(res) == 0) return(tibble::tibble())
  res %>%
    dplyr::group_by(idp_root) %>%
    dplyr::summarise(
      n = sum(is.finite(delta_rate), na.rm = TRUE),
      mean = mean(delta_rate, na.rm = TRUE),
      median = stats::median(delta_rate, na.rm = TRUE),
      sd = stats::sd(delta_rate, na.rm = TRUE),
      q25 = stats::quantile(delta_rate, 0.25, na.rm = TRUE, type = 7),
      q75 = stats::quantile(delta_rate, 0.75, na.rm = TRUE, type = 7),
      .groups = "drop"
    )
}

# ---- IDP Descriptive Analysis (Median, Delta, Rate) ----
compute_idp_group_contrast <- function(dat, idp_cols, group_col = "group", comparisons = c("ischemic_vs_control", "mi_vs_control", "chronic_vs_control")) {
  # 确保数据中有必要的列
  if (!all(c("eid", group_col) %in% names(dat))) return(NULL)
  
  res_list <- list()
  
  for (cmp in comparisons) {
    # 构建比较子集
    sub_dat <- tryCatch(build_comparison_df(dat, cmp), error = function(e) NULL)
    if (is.null(sub_dat)) next
    
    # 遍历IDP
    cmp_res <- purrr::map_df(idp_cols, function(v) {
      if (!(v %in% names(sub_dat))) return(NULL)
      val <- sub_dat[[v]]
      grp <- sub_dat[[group_col]]
      
      # 分组数据
      val0 <- val[grp == 0] # Control
      val1 <- val[grp == 1] # Case
      
      # 移除NA
      val0 <- val0[!is.na(val0)]
      val1 <- val1[!is.na(val1)]
      
      if (length(val0) == 0 || length(val1) == 0) return(NULL)
      
      # 计算中位数
      med0 <- stats::median(val0)
      med1 <- stats::median(val1)
      
      # 计算差异
      delta <- med1 - med0
      delta_rate <- delta
      
      data.frame(
        comparison = cmp,
        idp = v,
        median_control = med0,
        median_case = med1,
        delta = delta,
        delta_change_rate = delta_rate,
        n_control = length(val0),
        n_case = length(val1)
      )
    })
    
    res_list[[cmp]] <- cmp_res
  }
  
  dplyr::bind_rows(res_list)
}

# ---- 统一变量规则与基线三线表工具 ----
standardize_variable_rules <- function(dat) {
  nms <- names(dat)
  # 保证常用协变量为因子类型（若存在）
  to_factor <- c("sex_factor","ethnicity_factor",
                 "smoking_factor_baseline","alcohol_factor_baseline","diabetes_factor_baseline","education_factor","imaging_center_factor")
  for (col in to_factor) {
    if (col %in% nms) dat[[col]] <- as.factor(dat[[col]])
  }
  # 将 Instance.0 结尾的列尽量转换为数值型，便于描述性统计
  i0_cols <- nms[stringr::str_detect(nms, "\\.{3,}Instance\\.0$")]
  for (c in i0_cols) {
    dat[[c]] <- suppressWarnings(as.numeric(dat[[c]]))
  }
  # 一次性：从常见 Instance.0 源列创建规范化的 *_baseline 列（不使用 _locf 后缀）
  mk_baseline_from_i0 <- function(pattern, target, as_factor = FALSE) {
    cand <- nms[grepl(pattern, nms, ignore.case = TRUE)]
    if (length(cand) > 0 && !(target %in% names(dat))) {
      dat[[target]] <- dat[[cand[1]]]
      if (as_factor) dat[[target]] <- as.factor(dat[[target]])
    }
  }
  # BMI（基线）
  mk_baseline_from_i0("(\\bbmi\\b|body\\.?mass\\.?index).*\\.{3,}Instance\\.0$", "bmi_baseline")
  # 收缩压（基线）
  mk_baseline_from_i0("(systolic|blood\\.?pressure).*\\.{3,}Instance\\.0$", "systolic_bp_baseline")
  # Townsend 贫困指数（基线）
  mk_baseline_from_i0("townsend.*(deprivation|index).*\\.{3,}Instance\\.0$", "townsend_index_baseline")
  # 糖尿病诊断（基线，因子）
  mk_baseline_from_i0("diabet.*(diagnosis|status).*\\.{3,}Instance\\.0$", "diabetes_factor_baseline", as_factor = TRUE)
  # 吸烟状态（基线，因子）
  mk_baseline_from_i0("smok.*(status|tobacco).*\\.{3,}Instance\\.0$", "smoking_factor_baseline", as_factor = TRUE)
  # 酒精摄入频率（基线，因子）
  mk_baseline_from_i0("alcohol.*(intake|frequency).*\\.{3,}Instance\\.0$", "alcohol_factor_baseline", as_factor = TRUE)
  # 腰/臀围（便于后续派生腰臀比）
  mk_baseline_from_i0("waist.*circumference.*\\.{3,}Instance\\.0$", "waist_circumference_baseline")
  mk_baseline_from_i0("hip.*circumference.*\\.{3,}Instance\\.0$", "hip_circumference_baseline")
  mk_baseline_from_i0("diastolic.*blood\\.?pressure.*\\.{3,}Instance\\.0$", "diastolic_bp_baseline")
  mk_baseline_from_i0("body.*fat.*percentage.*\\.{3,}Instance\\.0$", "body_fat_percentage_baseline")
  mk_baseline_from_i0("whole.*body.*fat.*mass.*\\.{3,}Instance\\.0$", "whole_body_fat_mass_baseline")
  if (("waist_circumference_baseline" %in% names(dat)) && ("hip_circumference_baseline" %in% names(dat)) && !("waist_hip_ratio_baseline" %in% names(dat))) {
    suppressWarnings(dat$waist_hip_ratio_baseline <- as.numeric(dat$waist_circumference_baseline) / as.numeric(dat$hip_circumference_baseline))
  }
  cand_i2 <- nms[grepl("assessment.*centre.*\\.{3,}Instance\\.2$", nms, ignore.case = TRUE) & !grepl("^Date\\.of\\.", nms)]
  if (length(cand_i2) > 0 && !("imaging_center_factor" %in% names(dat))) {
    dat$imaging_center_factor <- as.factor(dat[[cand_i2[1]]])
  }
  dat
}


make_three_line_table <- function(dat, vars, group_col = "group", group_labels = c("Control","Case"), out_base_name = "Descriptive_Table") {
  vars <- vars[vars %in% names(dat)]
  if (length(vars) == 0) return(invisible(NULL))
  factor_vars <- vars[vapply(vars, function(v) is.factor(dat[[v]]) || is.character(dat[[v]]), logical(1))]
  tbl <- tableone::CreateTableOne(vars = vars, strata = group_col, data = dat, factorVars = factor_vars, test = TRUE)
  df <- as.data.frame(print(tbl, printToggle = FALSE, test = TRUE))
  # 友好标签改名
  cn <- colnames(df)
  cn[cn == "level 0"] <- group_labels[1]
  cn[cn == "level 1"] <- group_labels[2]
  colnames(df) <- cn
  readr::write_csv(df, table_out_path("Baseline", paste0(out_base_name, ".csv")))
  invisible(df)
}

# ------------------------------------------------------------
# 高标准SCI版 三线表（用于基线展示）
# ------------------------------------------------------------
make_three_line_table_sci <- function(dat, group_col = "group", group_labels = c("Control","Case"), out_base_name = "Baseline_SCI_Table", exclude_vars = NULL) {
  stopifnot(group_col %in% names(dat))
  # 统一变量名与类型
  if (!("baseline_age" %in% names(dat))) stop("需要 baseline_age 列")
  if (!("followup_years" %in% names(dat))) stop("需要 followup_years 列")
  if (!("sex_factor" %in% names(dat))) stop("需要 sex_factor 列")
  if (!("ethnicity_factor" %in% names(dat))) stop("需要 ethnicity_factor 列")

  dat <- dat %>% dplyr::mutate(
    group_lbl = factor(.data[[group_col]], levels = c(0,1), labels = group_labels),
    age_scan2 = baseline_age + followup_years
  )

  # 连续型变量的P值（自动在t检验与Wilcoxon间回退），支持动态分组列
  pv_cont <- function(df, y_name, grp_name) {
    fml <- stats::as.formula(sprintf("%s ~ %s", y_name, grp_name))
    tryCatch(
      stats::t.test(fml, data = df)$p.value,
      error = function(e) {
        tryCatch(stats::wilcox.test(fml, data = df)$p.value, error = function(e2) NA_real_)
      }
    )
  }

  fmt_mean_sd_range <- function(x) {
    x <- stats::na.omit(as.numeric(x))
    if (length(x) == 0) return("")
    paste0(sprintf("%.1f", mean(x)), "±", sprintf("%.1f", stats::sd(x)), " (", sprintf("%.1f", min(x)), "–", sprintf("%.1f", max(x)), ")")
  }
  fmt_count_pct <- function(x, level) {
    n <- sum(x == level, na.rm = TRUE)
    N <- sum(!is.na(x))
    pct <- ifelse(N > 0, 100*n/N, NA_real_)
    paste0(n, " (", sprintf("%.1f", pct), "%)")
  }

  # 分组向量
  g0 <- dat %>% dplyr::filter(.data[[group_col]] == 0)
  g1 <- dat %>% dplyr::filter(.data[[group_col]] == 1)

  # 统计量
  row_list <- list()
  # 1) 参与者数量
  row_list[["Number of participants"]] <- c(nrow(g1), nrow(g0), "-")

  # 2) 年龄（扫描1）
  p_age1 <- pv_cont(dat, "baseline_age", group_col)
  row_list[["Age at scan 1 (mean±s.d. [range])"]] <- c(
    fmt_mean_sd_range(g1$baseline_age),
    fmt_mean_sd_range(g0$baseline_age),
    sprintf("%.3f", p_age1)
  )

  # 3) 年龄（扫描2，派生）
  p_age2 <- pv_cont(dat, "age_scan2", group_col)
  row_list[["Age at scan 2 (mean±s.d. [range])"]] <- c(
    fmt_mean_sd_range(g1$age_scan2),
    fmt_mean_sd_range(g0$age_scan2),
    sprintf("%.3f", p_age2)
  )

  # 4) 性别（男/女）
  # 使用男性作为展示，括号中同时列出男/女（若需要更改可调整）
  male_level <- ifelse(any(levels(dat$sex_factor) %in% c("Male","male","M")), levels(dat$sex_factor)[grepl("^M|male|Male$", levels(dat$sex_factor))][1], levels(dat$sex_factor)[1])
  female_level <- setdiff(levels(dat$sex_factor), male_level)[1]
  # P值：二分类卡方检验
  sex_tbl <- table(dat[[group_col]], dat$sex_factor)
  p_sex <- tryCatch({
    suppressWarnings(stats::chisq.test(sex_tbl)$p.value)
  }, error = function(e){
    suppressWarnings(stats::fisher.test(sex_tbl)$p.value)
  })
  row_list[["Sex (male/female)"]] <- c(
    paste0(fmt_count_pct(g1$sex_factor, male_level), "/", fmt_count_pct(g1$sex_factor, female_level)),
    paste0(fmt_count_pct(g0$sex_factor, male_level), "/", fmt_count_pct(g0$sex_factor, female_level)),
    sprintf("%.3f", p_sex)
  )

  # 5) 种族（white/non-white）
  white_level <- ifelse(any(levels(dat$ethnicity_factor) %in% c("White","white")), levels(dat$ethnicity_factor)[grepl("White|white", levels(dat$ethnicity_factor))][1], levels(dat$ethnicity_factor)[1])
  nonwhite_levels <- setdiff(levels(dat$ethnicity_factor), white_level)
  g1_nonwhite <- dat %>% dplyr::filter(.data[[group_col]] == 1) %>% dplyr::pull(ethnicity_factor)
  g0_nonwhite <- dat %>% dplyr::filter(.data[[group_col]] == 0) %>% dplyr::pull(ethnicity_factor)
  p_eth <- tryCatch({
    suppressWarnings(stats::chisq.test(table(dat[[group_col]], dat$ethnicity_factor))$p.value)
  }, error = function(e){
    suppressWarnings(stats::fisher.test(table(dat[[group_col]], dat$ethnicity_factor))$p.value)
  })
  fmt_nonwhite <- function(x) {
    n <- sum(x %in% nonwhite_levels, na.rm = TRUE)
    N <- sum(!is.na(x))
    pct <- ifelse(N > 0, 100*n/N, NA_real_)
    paste0(n, " (", sprintf("%.1f", pct), "%)")
  }
  row_list[["Ethnicity (white/non-white)"]] <- c(
    paste0(fmt_count_pct(g1$ethnicity_factor, white_level), "/", fmt_nonwhite(g1_nonwhite)),
    paste0(fmt_count_pct(g0$ethnicity_factor, white_level), "/", fmt_nonwhite(g0_nonwhite)),
    sprintf("%.3f", p_eth)
  )

  # 6) 两次扫描间隔
  p_years <- pv_cont(dat, "followup_years", group_col)
  row_list[["Years between scans 1 and 2 (mean±s.d. [range])"]] <- c(
    fmt_mean_sd_range(g1$followup_years),
    fmt_mean_sd_range(g0$followup_years),
    sprintf("%.3f", p_years)
  )

  # 7) 收缩压（基线）
  if ("systolic_bp_baseline" %in% names(dat)) {
    p_sys <- pv_cont(dat, "systolic_bp_baseline", group_col)
    row_list[["Systolic blood pressure (mmHg; mean±s.d. [range])"]] <- c(
      fmt_mean_sd_range(g1$systolic_bp_baseline),
      fmt_mean_sd_range(g0$systolic_bp_baseline),
      sprintf("%.3f", p_sys)
    )
  }

  if (!("diastolic_bp_baseline" %in% names(dat))) {
    di_cols <- grep("diastolic.*blood.*pressure.*\\.{3,}Instance\\.0(\\.{3}Array\\.[0-9]+)?$", names(dat), ignore.case = TRUE, value = TRUE)
    if (length(di_cols) > 0) {
      tmp <- do.call(cbind, lapply(di_cols, function(c) suppressWarnings(as.numeric(dat[[c]]))))
      dat$diastolic_bp_baseline <- suppressWarnings(rowMeans(tmp, na.rm = TRUE))
      g1$diastolic_bp_baseline <- suppressWarnings(rowMeans(do.call(cbind, lapply(di_cols, function(c) as.numeric(g1[[c]]))), na.rm = TRUE))
      g0$diastolic_bp_baseline <- suppressWarnings(rowMeans(do.call(cbind, lapply(di_cols, function(c) as.numeric(g0[[c]]))), na.rm = TRUE))
    }
  }
  if ("diastolic_bp_baseline" %in% names(dat)) {
    p_dia <- pv_cont(dat, "diastolic_bp_baseline", group_col)
    row_list[["Diastolic blood pressure (mmHg; mean±s.d. [range])"]] <- c(
      fmt_mean_sd_range(g1$diastolic_bp_baseline),
      fmt_mean_sd_range(g0$diastolic_bp_baseline),
      sprintf("%.3f", p_dia)
    )
  }

  # 8) 诊断糖尿病（是/否）
  if ("diabetes_factor_baseline" %in% names(dat)) {
    yes_level <- ifelse(any(levels(dat$diabetes_factor_baseline) %in% c("Yes","Diabetes","Diagnosed")),
                        levels(dat$diabetes_factor_baseline)[grepl("Yes|Diabetes|Diagnosed", levels(dat$diabetes_factor_baseline))][1],
                        levels(dat$diabetes_factor_baseline)[1])
    no_level <- setNames(setdiff(levels(dat$diabetes_factor_baseline), yes_level), NULL)[1]
    diab_tbl <- table(dat[[group_col]], dat$diabetes_factor_baseline)
    p_diab <- tryCatch({
      suppressWarnings(stats::chisq.test(diab_tbl)$p.value)
    }, error = function(e){
      suppressWarnings(stats::fisher.test(diab_tbl)$p.value)
    })
    row_list[["Diagnosed diabetes (yes/no)"]] <- c(
      paste0(fmt_count_pct(g1$diabetes_factor_baseline, yes_level), "/", fmt_count_pct(g1$diabetes_factor_baseline, no_level)),
      paste0(fmt_count_pct(g0$diabetes_factor_baseline, yes_level), "/", fmt_count_pct(g0$diabetes_factor_baseline, no_level)),
      sprintf("%.3f", p_diab)
    )
  }

  # 9) BMI（基线）
  if ("bmi_baseline" %in% names(dat)) {
    p_bmi <- pv_cont(dat, "bmi_baseline", group_col)
    row_list[["BMI (kg/m²; mean±s.d. [range])"]] <- c(
      fmt_mean_sd_range(g1$bmi_baseline),
      fmt_mean_sd_range(g0$bmi_baseline),
      sprintf("%.3f", p_bmi)
    )
  }

  # 10) 腰围（基线，连续）
  # 优先使用 waist_circumference_baseline；若不存在则尝试从 Instance.0 源列派生
  waist_col <- NULL
  if ("waist_circumference_baseline" %in% names(dat)) {
    waist_col <- "waist_circumference_baseline"
  } else {
    # 兜底：直接查找 UKB 源列并作为腰围（单位：cm）
    i0_waist <- grep("waist.*circumference.*\\.{3,}Instance\\.0$", names(dat), ignore.case = TRUE, value = TRUE)
    if (length(i0_waist) > 0) {
      dat$waist_circumference_baseline <- suppressWarnings(as.numeric(dat[[i0_waist[1]]]))
      waist_col <- "waist_circumference_baseline"
      g1$waist_circumference_baseline <- suppressWarnings(as.numeric(g1[[i0_waist[1]]]))
      g0$waist_circumference_baseline <- suppressWarnings(as.numeric(g0[[i0_waist[1]]]))
    }
  }
  if (!is.null(waist_col)) {
    p_waist <- pv_cont(dat, waist_col, group_col)
    row_list[["Waist circumference (cm; mean±s.d. [range])"]] <- c(
      fmt_mean_sd_range(g1[[waist_col]]),
      fmt_mean_sd_range(g0[[waist_col]]),
      sprintf("%.3f", p_waist)
    )
  }

  if ("waist_hip_ratio_baseline" %in% names(dat)) {
    p_whr <- pv_cont(dat, "waist_hip_ratio_baseline", group_col)
    row_list[["Waist-hip ratio (mean±s.d. [range])"]] <- c(
      fmt_mean_sd_range(g1$waist_hip_ratio_baseline),
      fmt_mean_sd_range(g0$waist_hip_ratio_baseline),
      sprintf("%.3f", p_whr)
    )
  } else if (all(c("waist_circumference_baseline","hip_circumference_baseline") %in% names(dat))) {
    dat$waist_hip_ratio_baseline <- suppressWarnings(as.numeric(dat$waist_circumference_baseline) / as.numeric(dat$hip_circumference_baseline))
    g1$waist_hip_ratio_baseline <- suppressWarnings(as.numeric(g1$waist_circumference_baseline) / as.numeric(g1$hip_circumference_baseline))
    g0$waist_hip_ratio_baseline <- suppressWarnings(as.numeric(g0$waist_circumference_baseline) / as.numeric(g0$hip_circumference_baseline))
    p_whr <- pv_cont(dat, "waist_hip_ratio_baseline", group_col)
    row_list[["Waist-hip ratio (mean±s.d. [range])"]] <- c(
      fmt_mean_sd_range(g1$waist_hip_ratio_baseline),
      fmt_mean_sd_range(g0$waist_hip_ratio_baseline),
      sprintf("%.3f", p_whr)
    )
  }

  if (!("hip_circumference_baseline" %in% names(dat))) {
    i0_hip <- grep("hip.*circumference.*\\.{3,}Instance\\.0$", names(dat), ignore.case = TRUE, value = TRUE)
    if (length(i0_hip) > 0) {
      dat$hip_circumference_baseline <- suppressWarnings(as.numeric(dat[[i0_hip[1]]]))
      g1$hip_circumference_baseline <- suppressWarnings(as.numeric(g1[[i0_hip[1]]]))
      g0$hip_circumference_baseline <- suppressWarnings(as.numeric(g0[[i0_hip[1]]]))
    }
  }

  if (!("body_fat_percentage_baseline" %in% names(dat))) {
    i0_bfp <- grep("body.*fat.*percentage.*\\.{3,}Instance\\.0$", names(dat), ignore.case = TRUE, value = TRUE)
    if (length(i0_bfp) > 0) {
      dat$body_fat_percentage_baseline <- suppressWarnings(as.numeric(dat[[i0_bfp[1]]]))
      g1$body_fat_percentage_baseline <- suppressWarnings(as.numeric(g1[[i0_bfp[1]]]))
      g0$body_fat_percentage_baseline <- suppressWarnings(as.numeric(g0[[i0_bfp[1]]]))
    }
  }
  if ("body_fat_percentage_baseline" %in% names(dat)) {
    p_bfp <- pv_cont(dat, "body_fat_percentage_baseline", group_col)
    row_list[["Body fat percentage (mean±s.d. [range])"]] <- c(
      fmt_mean_sd_range(g1$body_fat_percentage_baseline),
      fmt_mean_sd_range(g0$body_fat_percentage_baseline),
      sprintf("%.3f", p_bfp)
    )
  }

  if (!("whole_body_fat_mass_baseline" %in% names(dat))) {
    i0_wbfm <- grep("whole.*body.*fat.*mass.*\\.{3,}Instance\\.0$", names(dat), ignore.case = TRUE, value = TRUE)
    if (length(i0_wbfm) > 0) {
      dat$whole_body_fat_mass_baseline <- suppressWarnings(as.numeric(dat[[i0_wbfm[1]]]))
      g1$whole_body_fat_mass_baseline <- suppressWarnings(as.numeric(g1[[i0_wbfm[1]]]))
      g0$whole_body_fat_mass_baseline <- suppressWarnings(as.numeric(g0[[i0_wbfm[1]]]))
    }
  }
  if ("whole_body_fat_mass_baseline" %in% names(dat)) {
    p_wbfm <- pv_cont(dat, "whole_body_fat_mass_baseline", group_col)
    row_list[["Whole body fat mass (kg; mean±s.d. [range])"]] <- c(
      fmt_mean_sd_range(g1$whole_body_fat_mass_baseline),
      fmt_mean_sd_range(g0$whole_body_fat_mass_baseline),
      sprintf("%.3f", p_wbfm)
    )
  }

  # 11) 酒精摄入频率（有序编码后以均值±SD展示）
  if ("alcohol_factor_baseline" %in% names(dat)) {
    ord_score <- function(f) {
      lv <- levels(f)
      score_map <- setNames(seq_along(lv), lv)
      # 语义化顺序优化
      score_map[grepl("Never|No", lv, ignore.case = TRUE)] <- 0
      score_map[grepl("Special|rare|occasion", lv, ignore.case = TRUE)] <- 1
      score_map[grepl("month", lv, ignore.case = TRUE)] <- 2
      score_map[grepl("1-2|once|twice.*week|weekly", lv, ignore.case = TRUE)] <- 3
      score_map[grepl("3-4|several.*week", lv, ignore.case = TRUE)] <- 4
      score_map[grepl("daily|almost", lv, ignore.case = TRUE)] <- 5
      as.numeric(score_map[as.character(f)])
    }
    dat$alcohol_score <- ord_score(dat$alcohol_factor_baseline)
    g1$alcohol_score <- ord_score(g1$alcohol_factor_baseline)
    g0$alcohol_score <- ord_score(g0$alcohol_factor_baseline)
    p_alc <- pv_cont(dat, "alcohol_score", group_col)
    row_list[["Alcohol-intake frequency (a.u.; mean±s.d. [range])"]] <- c(
      fmt_mean_sd_range(g1$alcohol_score),
      fmt_mean_sd_range(g0$alcohol_score),
      sprintf("%.3f", p_alc)
    )
  }

  # 12) 吸烟（有序编码：从不/既往/当前）
  if ("smoking_factor_baseline" %in% names(dat)) {
    ord_smoke <- function(f) {
      lv <- levels(f)
      score_map <- setNames(seq_along(lv), lv)
      score_map[grepl("Never|Non|No", lv, ignore.case = TRUE)] <- 0
      score_map[grepl("Former|Previous|Ex", lv, ignore.case = TRUE)] <- 1
      score_map[grepl("Current|Yes", lv, ignore.case = TRUE)] <- 2
      as.numeric(score_map[as.character(f)])
    }
    dat$smoking_score <- ord_smoke(dat$smoking_factor_baseline)
    g1$smoking_score <- ord_smoke(g1$smoking_factor_baseline)
    g0$smoking_score <- ord_smoke(g0$smoking_factor_baseline)
    p_smoke <- pv_cont(dat, "smoking_score", group_col)
    row_list[["Tobacco smoking (a.u.; mean±s.d. [range])"]] <- c(
      fmt_mean_sd_range(g1$smoking_score),
      fmt_mean_sd_range(g0$smoking_score),
      sprintf("%.3f", p_smoke)
    )
  }

  # 13) Townsend 贫困指数（基线，连续）
  if ("townsend_index_baseline" %in% names(dat)) {
    p_town <- pv_cont(dat, "townsend_index_baseline", group_col)
    row_list[["Townsend deprivation index (mean±s.d. [range])"]] <- c(
      fmt_mean_sd_range(g1$townsend_index_baseline),
      fmt_mean_sd_range(g0$townsend_index_baseline),
      sprintf("%.3f", p_town)
    )
  }

  # ============ 补充：Townsend（原始连续）、教育（raw/派生/因子）、生活方式 & 原始连续指标 ============
  # 14) Townsend（原始，连续）
  if ("townsend_index" %in% names(dat)) {
    p_town_raw <- pv_cont(dat, "townsend_index", group_col)
    row_list[["Townsend deprivation index (raw; mean±s.d. [range])"]] <- c(
      fmt_mean_sd_range(g1$townsend_index),
      fmt_mean_sd_range(g0$townsend_index),
      sprintf("%.3f", p_town_raw)
    )
  }

  # 小工具：按给定水平顺序输出“计数(%)”并以“/”拼接
  fmt_multi_levels <- function(x, levels_order) {
    f <- factor(x, levels = levels_order)
    paste(vapply(levels_order, function(lv) fmt_count_pct(f, lv), character(1)), collapse = "/")
  }

  # 15) 教育（因子：Lower/Medium/Higher）
  if ("education_factor" %in% names(dat)) {
    lv_ed <- levels(dat$education_factor)
    if (length(lv_ed) == 3 && all(c("Lower","Medium","Higher") %in% lv_ed)) {
      ed_levels <- c("Lower","Medium","Higher")
    } else {
      ed_levels <- lv_ed
    }
    ed_tbl <- table(dat[[group_col]], dat$education_factor)
    p_ed <- tryCatch({
      suppressWarnings(stats::chisq.test(ed_tbl)$p.value)
    }, error = function(e){
      suppressWarnings(stats::fisher.test(ed_tbl)$p.value)
    })
    row_list[["Education (lower/medium/higher)"]] <- c(
      fmt_multi_levels(g1$education_factor, ed_levels),
      fmt_multi_levels(g0$education_factor, ed_levels),
      sprintf("%.3f", p_ed)
    )
  }

  # 16) 教育（派生：与education_factor一致的三分类）
  if ("education" %in% names(dat)) {
    ed_levels <- c("Lower","Medium","Higher")
    ed_tbl2 <- table(dat[[group_col]], factor(dat$education, levels = ed_levels))
    p_ed2 <- tryCatch({
      suppressWarnings(stats::chisq.test(ed_tbl2)$p.value)
    }, error = function(e){
      suppressWarnings(stats::fisher.test(ed_tbl2)$p.value)
    })
    row_list[["Education (derived; lower/medium/higher)"]] <- c(
      fmt_multi_levels(g1$education, ed_levels),
      fmt_multi_levels(g0$education, ed_levels),
      sprintf("%.3f", p_ed2)
    )
  }

  # 17) 教育（原始 raw；显示Top-3类别）
  if ("education_raw" %in% names(dat)) {
    vec_all <- dat$education_raw
    vec_all <- vec_all[!is.na(vec_all) & vec_all != ""]
    top3 <- names(sort(table(vec_all), decreasing = TRUE))[seq_len(min(3, length(unique(vec_all))))]
    # 若top3为空则跳过
    if (length(top3) > 0) {
      ed_raw_tbl <- table(dat[[group_col]], factor(dat$education_raw))
      p_edraw <- tryCatch({
        suppressWarnings(stats::chisq.test(ed_raw_tbl)$p.value)
      }, error = function(e){
        suppressWarnings(stats::fisher.test(ed_raw_tbl)$p.value)
      })
      fmt_topk <- function(df_group) {
        paste(vapply(top3, function(cat){
          n <- sum(df_group == cat, na.rm = TRUE)
          N <- sum(!is.na(df_group))
          pct <- ifelse(N > 0, 100*n/N, NA_real_)
          paste0(cat, "=", n, " (", sprintf("%.1f", pct), "%)")
        }, character(1)), collapse = "; ")
      }
      row_list[["Education (raw; top-3 categories)"]] <- c(
        fmt_topk(g1$education_raw),
        fmt_topk(g0$education_raw),
        sprintf("%.3f", p_edraw)
      )
    }
  }

  # 18) BMI（原始，连续）
  if ("bmi" %in% names(dat)) {
    p_bmi_raw <- pv_cont(dat, "bmi", group_col)
    row_list[["BMI (raw; mean±s.d. [range])"]] <- c(
      fmt_mean_sd_range(g1$bmi),
      fmt_mean_sd_range(g0$bmi),
      sprintf("%.3f", p_bmi_raw)
    )
  }

  # 19) 收缩压（原始，连续）
  if ("systolic_bp" %in% names(dat)) {
    p_sys_raw <- pv_cont(dat, "systolic_bp", group_col)
    row_list[["Systolic blood pressure (raw; mean±s.d. [range])"]] <- c(
      fmt_mean_sd_range(g1$systolic_bp),
      fmt_mean_sd_range(g0$systolic_bp),
      sprintf("%.3f", p_sys_raw)
    )
  }

  # 20) 糖尿病（原始：Yes/No）
  if ("diabetes" %in% names(dat)) {
    diab_levels <- c("No","Yes")
    diab_tbl_raw <- table(dat[[group_col]], factor(dat$diabetes, levels = diab_levels))
    p_diab_raw <- tryCatch({
      suppressWarnings(stats::chisq.test(diab_tbl_raw)$p.value)
    }, error = function(e){
      suppressWarnings(stats::fisher.test(diab_tbl_raw)$p.value)
    })
    row_list[["Diagnosed diabetes (raw yes/no)"]] <- c(
      fmt_multi_levels(g1$diabetes, diab_levels),
      fmt_multi_levels(g0$diabetes, diab_levels),
      sprintf("%.3f", p_diab_raw)
    )
  }

  # 21) 糖尿病（因子：Yes/No）
  if ("diabetes_factor" %in% names(dat)) {
    lv <- levels(dat$diabetes_factor)
    diab_levels <- if (all(c("No","Yes") %in% lv)) c("No","Yes") else lv
    diab_tbl_fac <- table(dat[[group_col]], dat$diabetes_factor)
    p_diab_fac <- tryCatch({
      suppressWarnings(stats::chisq.test(diab_tbl_fac)$p.value)
    }, error = function(e){
      suppressWarnings(stats::fisher.test(diab_tbl_fac)$p.value)
    })
    row_list[["Diagnosed diabetes (factor yes/no)"]] <- c(
      fmt_multi_levels(g1$diabetes_factor, diab_levels),
      fmt_multi_levels(g0$diabetes_factor, diab_levels),
      sprintf("%.3f", p_diab_fac)
    )
  }

  # 22) 吸烟（原始：Never/Former/Current）
  if ("smoking" %in% names(dat)) {
    smk_levels <- c("Never","Former","Current")
    smk_tbl_raw <- table(dat[[group_col]], factor(dat$smoking, levels = smk_levels))
    p_smk_raw <- tryCatch({
      suppressWarnings(stats::chisq.test(smk_tbl_raw)$p.value)
    }, error = function(e){
      suppressWarnings(stats::fisher.test(smk_tbl_raw)$p.value)
    })
    row_list[["Tobacco smoking (raw: never/former/current)"]] <- c(
      fmt_multi_levels(g1$smoking, smk_levels),
      fmt_multi_levels(g0$smoking, smk_levels),
      sprintf("%.3f", p_smk_raw)
    )
  }

  # 23) 吸烟（因子：Never/Former/Current）
  if ("smoking_factor" %in% names(dat)) {
    lv <- levels(dat$smoking_factor)
    smk_levels <- if (all(c("Never","Former","Current") %in% lv)) c("Never","Former","Current") else lv
    smk_tbl_fac <- table(dat[[group_col]], dat$smoking_factor)
    p_smk_fac <- tryCatch({
      suppressWarnings(stats::chisq.test(smk_tbl_fac)$p.value)
    }, error = function(e){
      suppressWarnings(stats::fisher.test(smk_tbl_fac)$p.value)
    })
    row_list[["Tobacco smoking (factor: never/former/current)"]] <- c(
      fmt_multi_levels(g1$smoking_factor, smk_levels),
      fmt_multi_levels(g0$smoking_factor, smk_levels),
      sprintf("%.3f", p_smk_fac)
    )
  }

  # 24) 酒精摄入（原始：Never/Former/Current）
  if ("alcohol" %in% names(dat)) {
    alc_levels <- c("Never","Former","Current")
    alc_tbl_raw <- table(dat[[group_col]], factor(dat$alcohol, levels = alc_levels))
    p_alc_raw <- tryCatch({
      suppressWarnings(stats::chisq.test(alc_tbl_raw)$p.value)
    }, error = function(e){
      suppressWarnings(stats::fisher.test(alc_tbl_raw)$p.value)
    })
    row_list[["Alcohol-intake frequency (raw: never/former/current)"]] <- c(
      fmt_multi_levels(g1$alcohol, alc_levels),
      fmt_multi_levels(g0$alcohol, alc_levels),
      sprintf("%.3f", p_alc_raw)
    )
  }

  # 25) 酒精摄入（因子：Never/Former/Current）
  if ("alcohol_factor" %in% names(dat)) {
    lv <- levels(dat$alcohol_factor)
    alc_levels <- if (all(c("Never","Former","Current") %in% lv)) c("Never","Former","Current") else lv
    alc_tbl_fac <- table(dat[[group_col]], dat$alcohol_factor)
    p_alc_fac <- tryCatch({
      suppressWarnings(stats::chisq.test(alc_tbl_fac)$p.value)
    }, error = function(e){
      suppressWarnings(stats::fisher.test(alc_tbl_fac)$p.value)
    })
    row_list[["Alcohol-intake frequency (factor: never/former/current)"]] <- c(
      fmt_multi_levels(g1$alcohol_factor, alc_levels),
      fmt_multi_levels(g0$alcohol_factor, alc_levels),
      sprintf("%.3f", p_alc_fac)
    )
  }

  # 组装三线表
  df <- do.call(rbind, lapply(names(row_list), function(nm){
    v <- row_list[[nm]]
    data.frame(
      Variable = nm,
      `Case` = v[1],
      `Control` = v[2],
      P_uncorr = v[3],
      stringsAsFactors = FALSE
    )
  }))

  if (!is.null(exclude_vars)) {
    df <- df[!(df$Variable %in% exclude_vars), , drop = FALSE]
  }
  out_file <- table_out_path("Baseline", paste0(out_base_name, ".csv"))
  readr::write_csv(df, out_file)
  invisible(df)
}

make_topjournal_baseline_table <- function(dat, group_col = "group", out_file, case_label = "Ischemic", control_label = "Control") {
  stopifnot(group_col %in% names(dat))
  dat <- dat[!is.na(dat[[group_col]]) & dat[[group_col]] %in% c(0, 1), , drop = FALSE]
  dat[[group_col]] <- as.integer(dat[[group_col]])
  g0 <- dat[dat[[group_col]] == 0, , drop = FALSE]
  g1 <- dat[dat[[group_col]] == 1, , drop = FALSE]
  pv_cont <- function(df, y_name, grp_name) {
    fml <- stats::as.formula(sprintf("%s ~ %s", y_name, grp_name))
    tryCatch(
      stats::t.test(fml, data = df)$p.value,
      error = function(e) {
        tryCatch(stats::wilcox.test(fml, data = df)$p.value, error = function(e2) NA_real_)
      }
    )
  }
  pv_cat <- function(df, x, grp_name) {
    tt <- table(df[[grp_name]], x, useNA = "no")
    if (nrow(tt) < 2 || ncol(tt) < 2) return(NA_real_)
    tryCatch(
      suppressWarnings(stats::chisq.test(tt)$p.value),
      error = function(e) suppressWarnings(stats::fisher.test(tt)$p.value)
    )
  }
  fmt_mean_sd_range <- function(x) {
    x <- stats::na.omit(suppressWarnings(as.numeric(x)))
    if (length(x) == 0) return("")
    paste0(
      sprintf("%.1f", mean(x)), "±", sprintf("%.1f", stats::sd(x)),
      " (", sprintf("%.1f", min(x)), "–", sprintf("%.1f", max(x)), ")"
    )
  }
  fmt_count_pct <- function(x, level) {
    n <- sum(x == level, na.rm = TRUE)
    N <- sum(!is.na(x))
    pct <- ifelse(N > 0, 100 * n / N, NA_real_)
    paste0(n, " (", sprintf("%.1f", pct), "%)")
  }
  fmt_yes_only <- function(x, yes_level) {
    n <- sum(x == yes_level, na.rm = TRUE)
    N <- sum(!is.na(x))
    pct <- ifelse(N > 0, 100 * n / N, NA_real_)
    paste0(n, " (", sprintf("%.1f", pct), "%)")
  }
  num_col <- function(df, candidates) {
    nm <- candidates[candidates %in% names(df)]
    if (length(nm) == 0) return(rep(NA_real_, nrow(df)))
    suppressWarnings(as.numeric(df[[nm[[1]]]]))
  }
  age1 <- num_col(dat, c("baseline_age", "age_at_instance2", "age_at_scan_1", "age_at_recruitment"))
  age2 <- rep(NA_real_, nrow(dat))
  if ("followup_age" %in% names(dat)) {
    age2 <- suppressWarnings(as.numeric(dat$followup_age))
  } else if ("age_at_instance3" %in% names(dat)) {
    age2 <- suppressWarnings(as.numeric(dat$age_at_instance3))
  } else if (all(c("baseline_age", "followup_years") %in% names(dat))) {
    age2 <- suppressWarnings(as.numeric(dat$baseline_age) + as.numeric(dat$followup_years))
  } else if (all(c("age_at_recruitment", "followup_years") %in% names(dat))) {
    age2 <- suppressWarnings(as.numeric(dat$age_at_recruitment) + as.numeric(dat$followup_years))
  }
  dat$.__age1 <- age1
  dat$.__age2 <- age2
  g0$.__age1 <- age1[dat[[group_col]] == 0]
  g1$.__age1 <- age1[dat[[group_col]] == 1]
  g0$.__age2 <- age2[dat[[group_col]] == 0]
  g1$.__age2 <- age2[dat[[group_col]] == 1]
  sex_col <- if ("sex_factor" %in% names(dat)) "sex_factor" else NULL
  eth_col <- if ("ethnicity_factor" %in% names(dat)) "ethnicity_factor" else NULL
  diab_col <- if ("diabetes_factor_baseline" %in% names(dat)) "diabetes_factor_baseline" else if ("diabetes_factor" %in% names(dat)) "diabetes_factor" else NULL
  sex <- if (!is.null(sex_col)) as.factor(dat[[sex_col]]) else factor(NA)
  eth <- if (!is.null(eth_col)) as.factor(dat[[eth_col]]) else factor(NA)
  diab <- if (!is.null(diab_col)) as.factor(dat[[diab_col]]) else factor(NA)
  sbp <- num_col(dat, c("systolic_bp_baseline", "systolic_bp", "systolic_blood_pressure_baseline"))
  dbp <- num_col(dat, c("diastolic_bp_baseline", "diastolic_bp", "diastolic_blood_pressure_baseline"))
  years <- num_col(dat, c("followup_years"))
  whr <- num_col(dat, c("waist_hip_ratio_baseline", "waist_hip_ratio"))
  if (all(is.na(whr))) {
    waist <- num_col(dat, c(
      "waist_circumference_baseline",
      "waist_circumference",
      "Waist.circumference...Instance.0"
    ))
    hip <- num_col(dat, c(
      "hip_circumference_baseline",
      "hip_circumference",
      "Hip.circumference...Instance.0"
    ))
    whr <- suppressWarnings(waist / hip)
  }
  bmi <- num_col(dat, c("bmi_baseline", "bmi"))
  town <- num_col(dat, c("townsend_index_baseline", "townsend_index"))
  dat$.__sbp <- sbp
  dat$.__dbp <- dbp
  dat$.__years <- years
  dat$.__whr <- whr
  dat$.__bmi <- bmi
  dat$.__town <- town
  g0$.__sbp <- sbp[dat[[group_col]] == 0]
  g1$.__sbp <- sbp[dat[[group_col]] == 1]
  g0$.__dbp <- dbp[dat[[group_col]] == 0]
  g1$.__dbp <- dbp[dat[[group_col]] == 1]
  g0$.__years <- years[dat[[group_col]] == 0]
  g1$.__years <- years[dat[[group_col]] == 1]
  g0$.__whr <- whr[dat[[group_col]] == 0]
  g1$.__whr <- whr[dat[[group_col]] == 1]
  g0$.__bmi <- bmi[dat[[group_col]] == 0]
  g1$.__bmi <- bmi[dat[[group_col]] == 1]
  g0$.__town <- town[dat[[group_col]] == 0]
  g1$.__town <- town[dat[[group_col]] == 1]
  lv_sex <- levels(sex)
  male_idx <- which(tolower(lv_sex) %in% c("male", "m"))
  male_level <- if (length(male_idx) > 0) lv_sex[male_idx[[1]]] else lv_sex[[1]]
  female_level <- setdiff(levels(sex), male_level)[1]
  white_level <- ifelse(any(levels(eth) %in% c("White", "white")), levels(eth)[grepl("White|white", levels(eth))][1], levels(eth)[1])
  nonwhite_levels <- setdiff(levels(eth), white_level)
  yes_level <- if (length(levels(diab)) == 0) NA_character_ else ifelse(any(levels(diab) %in% c("Yes", "YES", "Diabetes", "Diagnosed")), levels(diab)[grepl("Yes|Diabetes|Diagnosed", levels(diab))][1], levels(diab)[length(levels(diab))])
  row_list <- list()
  row_list[["Number of participants"]] <- list(val = c(as.character(nrow(g1)), as.character(nrow(g0))), p = NA_real_)
  p_age1 <- pv_cont(dat, ".__age1", group_col)
  row_list[["Age at scan 1 (mean±s.d. (range))"]] <- list(val = c(fmt_mean_sd_range(g1$.__age1), fmt_mean_sd_range(g0$.__age1)), p = p_age1)
  p_age2 <- pv_cont(dat, ".__age2", group_col)
  row_list[["Age at scan 2 (mean±s.d. (range))"]] <- list(val = c(fmt_mean_sd_range(g1$.__age2), fmt_mean_sd_range(g0$.__age2)), p = p_age2)
  p_sex <- if (!is.null(sex_col)) pv_cat(dat, sex, group_col) else NA_real_
  row_list[["Sex (male/female)"]] <- list(
    val = c(
      if (!is.null(sex_col)) paste0(fmt_count_pct(as.factor(g1[[sex_col]]), male_level), "/", fmt_count_pct(as.factor(g1[[sex_col]]), female_level)) else "",
      if (!is.null(sex_col)) paste0(fmt_count_pct(as.factor(g0[[sex_col]]), male_level), "/", fmt_count_pct(as.factor(g0[[sex_col]]), female_level)) else ""
    ),
    p = p_sex
  )
  p_eth <- if (!is.null(eth_col)) pv_cat(dat, eth, group_col) else NA_real_
  fmt_nonwhite <- function(x) {
    n <- sum(x %in% nonwhite_levels, na.rm = TRUE)
    N <- sum(!is.na(x))
    pct <- ifelse(N > 0, 100 * n / N, NA_real_)
    paste0(n, " (", sprintf("%.1f", pct), "%)")
  }
  row_list[["Ethnicity (white/non-white)"]] <- list(
    val = c(
      if (!is.null(eth_col)) paste0(fmt_count_pct(as.factor(g1[[eth_col]]), white_level), "/", fmt_nonwhite(as.factor(g1[[eth_col]]))) else "",
      if (!is.null(eth_col)) paste0(fmt_count_pct(as.factor(g0[[eth_col]]), white_level), "/", fmt_nonwhite(as.factor(g0[[eth_col]]))) else ""
    ),
    p = p_eth
  )
  p_years <- pv_cont(dat, ".__years", group_col)
  row_list[["Years between scans 1 and 2 (mean±s.d. (range))"]] <- list(val = c(fmt_mean_sd_range(g1$.__years), fmt_mean_sd_range(g0$.__years)), p = p_years)
  p_sbp <- pv_cont(dat, ".__sbp", group_col)
  row_list[["Systolic blood pressure (mmHg)"]] <- list(val = c(fmt_mean_sd_range(g1$.__sbp), fmt_mean_sd_range(g0$.__sbp)), p = p_sbp)
  p_dbp <- pv_cont(dat, ".__dbp", group_col)
  row_list[["Diastolic blood pressure (mmHg)"]] <- list(val = c(fmt_mean_sd_range(g1$.__dbp), fmt_mean_sd_range(g0$.__dbp)), p = p_dbp)
  p_diab <- if (!is.null(diab_col)) pv_cat(dat, diab, group_col) else NA_real_
  row_list[["Diagnosed diabetes"]] <- list(
    val = c(
      if (!is.null(diab_col) && !is.na(yes_level)) fmt_yes_only(as.factor(g1[[diab_col]]), yes_level) else "",
      if (!is.null(diab_col) && !is.na(yes_level)) fmt_yes_only(as.factor(g0[[diab_col]]), yes_level) else ""
    ),
    p = p_diab
  )
  p_whr <- pv_cont(dat, ".__whr", group_col)
  row_list[["Waist/hip ratio"]] <- list(val = c(fmt_mean_sd_range(g1$.__whr), fmt_mean_sd_range(g0$.__whr)), p = p_whr)
  p_bmi <- pv_cont(dat, ".__bmi", group_col)
  row_list[["BMI (kg m−2)"]] <- list(val = c(fmt_mean_sd_range(g1$.__bmi), fmt_mean_sd_range(g0$.__bmi)), p = p_bmi)
  p_town <- pv_cont(dat, ".__town", group_col)
  row_list[["Townsend deprivation index"]] <- list(val = c(fmt_mean_sd_range(g1$.__town), fmt_mean_sd_range(g0$.__town)), p = p_town)
  df_out <- do.call(rbind, lapply(names(row_list), function(nm) {
    vv <- row_list[[nm]]
    data.frame(
      Variable = nm,
      Case = vv$val[[1]],
      Control = vv$val[[2]],
      P_uncorr = vv$p,
      stringsAsFactors = FALSE
    )
  }))
  names(df_out) <- c("Variable", case_label, control_label, "P_uncorr")
  utils::write.csv(df_out, out_file, row.names = FALSE, quote = TRUE, na = "NA")
  invisible(df_out)
}

summarize_age_distributions <- function(dat, comparison, group_col = "group", group_labels = c("Control","Case")) {
  grp <- factor(dat[[group_col]], levels = c(0,1), labels = group_labels)
  age1 <- if ("age_at_instance2" %in% names(dat)) suppressWarnings(as.numeric(dat$age_at_instance2)) else suppressWarnings(as.numeric(dat$baseline_age))
  age2 <- if ("age_at_instance3" %in% names(dat)) suppressWarnings(as.numeric(dat$age_at_instance3)) else if ("followup_age" %in% names(dat)) suppressWarnings(as.numeric(dat$followup_age)) else if (all(c("baseline_age","followup_years") %in% names(dat))) suppressWarnings(as.numeric(dat$baseline_age) + as.numeric(dat$followup_years)) else rep(NA_real_, length(age1))
  df <- data.frame(age1 = age1, age2 = age2, group = grp)
  df <- df[is.finite(df$age1) & is.finite(df$age2) & !is.na(df$group), , drop = FALSE]
  if (nrow(df) == 0) return(invisible(NULL))
  lo <- floor(min(c(df$age1, df$age2), na.rm = TRUE)/5)*5
  hi <- ceiling(max(c(df$age1, df$age2), na.rm = TRUE)/5)*5
  breaks <- seq(lo, hi, by = 5)
  df$band1 <- cut(df$age1, breaks = breaks, include.lowest = TRUE, right = FALSE)
  df$band2 <- cut(df$age2, breaks = breaks, include.lowest = TRUE, right = FALSE)
  tab1 <- df %>% dplyr::group_by(band1, group) %>% dplyr::tally(name = "n") %>% dplyr::mutate(age_band = as.character(band1)) %>% dplyr::select(age_band, group, n) %>% dplyr::arrange(age_band, group)
  tab2 <- df %>% dplyr::group_by(band2, group) %>% dplyr::tally(name = "n") %>% dplyr::mutate(age_band = as.character(band2)) %>% dplyr::select(age_band, group, n) %>% dplyr::arrange(age_band, group)
  readr::write_csv(tab1, table_out_path(file.path("Baseline","Age"), paste0("Age_Band_Distribution_", comparison, "_Scan1.csv")))
  readr::write_csv(tab2, table_out_path(file.path("Baseline","Age"), paste0("Age_Band_Distribution_", comparison, "_Scan2.csv")))
  x1_case <- df$age1[df$group == group_labels[2]]
  x1_ctrl <- df$age1[df$group == group_labels[1]]
  x2_case <- df$age2[df$group == group_labels[2]]
  x2_ctrl <- df$age2[df$group == group_labels[1]]
  ks1_p <- suppressWarnings(tryCatch(stats::ks.test(x1_case, x1_ctrl)$p.value, error = function(e) NA_real_))
  ks2_p <- suppressWarnings(tryCatch(stats::ks.test(x2_case, x2_ctrl)$p.value, error = function(e) NA_real_))
  has_nortest <- requireNamespace("nortest", quietly = TRUE)
  lillie_p <- function(x) {
    if (has_nortest) {
      suppressWarnings(nortest::lillie.test(x)$p.value)
    } else {
      suppressWarnings(stats::shapiro.test(x)$p.value)
    }
  }
  p_l_case1 <- lillie_p(x1_case)
  p_l_ctrl1 <- lillie_p(x1_ctrl)
  p_l_case2 <- lillie_p(x2_case)
  p_l_ctrl2 <- lillie_p(x2_ctrl)
  signif1 <- if (!is.na(ks1_p) && ks1_p < 0.05) "differ significantly" else "do not differ significantly"
  signif2 <- if (!is.na(ks2_p) && ks2_p < 0.05) "differ significantly" else "do not differ significantly"
  method_lab <- if (has_nortest) "Lilliefors" else "Shapiro-Wilk"
  txt <- paste0(
    "Age distributions for ischemic heart desease participants and controls at each time point ",
    if (!is.na(ks1_p) && !is.na(ks2_p) && ks1_p >= 0.05 && ks2_p >= 0.05) "do not differ significantly." else "show differences (see P values below).",
    "\nTwo-sample Kolmogorov-Smirnov was used to compute the P values for age comparisons, since age for each group was not normally distributed (", method_lab, " P values).\n",
    "Scan 1 KS P = ", format(ks1_p, digits = 3), "; Scan 2 KS P = ", format(ks2_p, digits = 3), ".\n",
    "Normality test P values (", method_lab, "): Scan 1 Case=", format(p_l_case1, digits = 3), ", Control=", format(p_l_ctrl1, digits = 3), "; Scan 2 Case=", format(p_l_case2, digits = 3), ", Control=", format(p_l_ctrl2, digits = 3), ".\n"
  )
  con <- file(table_out_path(file.path("Baseline","Age"), paste0("Age_Distribution_Description_", comparison, ".txt")), open = "w", encoding = "UTF-8")
  writeLines(txt, con = con, sep = "\n")
  close(con)
  fill_cols <- c("Control" = "#1F78B4", "Case" = "#FF7F00")
  p1 <- ggplot2::ggplot(df, ggplot2::aes(x = age1, fill = group)) +
    ggplot2::geom_histogram(binwidth = 5, position = "identity", alpha = 0.5, color = NA, boundary = 0) +
    ggplot2::scale_fill_manual(values = fill_cols, name = "Group") +
    ggplot2::labs(title = "Age at Scan 1", x = "Age (years)", y = "Count") +
    theme_jama_refined()
  p2 <- ggplot2::ggplot(df, ggplot2::aes(x = age2, fill = group)) +
    ggplot2::geom_histogram(binwidth = 5, position = "identity", alpha = 0.5, color = NA, boundary = 0) +
    ggplot2::scale_fill_manual(values = fill_cols, name = "Group") +
    ggplot2::labs(title = "Age at Scan 2", x = "Age (years)", y = "Count") +
    theme_jama_refined()
  g <- gridExtra::arrangeGrob(p1, p2, ncol = 1)
  save_pub(paste0("Figure_Age_Distribution_", comparison), g, width = 10, height = 8, dpi = 300)
  invisible(list(ks_scan1 = ks1_p, ks_scan2 = ks2_p))
}

if (RUN_PSM_MATCHING) {
# ---- 主执行：四个比较的匹配、汇总与可视化 ----
## 移除try包装，保证流程通畅，改为直接执行
  # 识别IDP列与PSM协变量
  # 应用统一变量规则，保证类型一致与基线列可用于描述性统计
  ischemic_cohort_complete <- standardize_variable_rules(ischemic_cohort_complete)
  idp_cols <- select_idp_cols(ischemic_cohort_complete)
  # PSM匹配协变量（按最新要求）：Age at scan 1 与 Age at scan 2 独立协变量
  # 同时保留性别、种族与两次扫描间隔
  psm_covariates <- c(
    "age_at_instance2",
    "age_at_recruitment",
    "sex_factor",
    "ethnicity_factor",
    "education_factor",
    "imaging_center_factor"
  )
  # 先不与总表交集，等待派生年龄列后再与主比较数据交集
  covariates <- psm_covariates
  
  # 仅在总队列（ischemic vs control）执行一次PSM
  main_cmp <- "ischemic_vs_control"
  dat_main <- build_comparison_df(ischemic_cohort_complete, main_cmp)
  # 在PSM前派生年龄列，确保 age_at_instance2/3 可用（与 compute_age_instances 逻辑一致）
  {
    n <- nrow(dat_main)
    if (!("baseline_age" %in% names(dat_main)) && ("age_at_recruitment" %in% names(dat_main))) {
      dat_main$baseline_age <- dat_main$age_at_recruitment
    }
    if (!("age_at_recruitment" %in% names(dat_main)) && ("baseline_age" %in% names(dat_main))) {
      dat_main$age_at_recruitment <- dat_main$baseline_age
    }
    if (!("followup_age" %in% names(dat_main)) || length(dat_main$followup_age) != n) {
      dat_main$followup_age <- rep(NA_real_, n)
    }
    if (!("followup_years" %in% names(dat_main)) || length(dat_main$followup_years) != n) {
      dat_main$followup_years <- rep(NA_real_, n)
    }
    # Age at scan 1
    dat_main$age_at_instance2 <- if ("baseline_age" %in% names(dat_main)) dat_main$baseline_age else dat_main$age_at_recruitment
    # Age at scan 2（优先 followup_age，其次 recruitment+followup_years，否则 +2 年）
    dat_main$age_at_instance3 <- dat_main$age_at_recruitment + 2
    has_fup_years <- !is.na(dat_main$followup_years)
    dat_main$age_at_instance3[has_fup_years] <- dat_main$age_at_recruitment[has_fup_years] + dat_main$followup_years[has_fup_years]
    has_fup_age <- !is.na(dat_main$followup_age)
    dat_main$age_at_instance3[has_fup_age] <- dat_main$followup_age[has_fup_age]
    # 补充：在派生完年龄后，重新与主比较数据交集，确保年龄被纳入PSM协变量
    covariates <- psm_covariates %>% intersect(names(dat_main))
  }
  # 非影像协变量（用于四模型与平衡报告，采用LOCF统一命名）
  # 统一使用 *_baseline 别名，避免 _locf 后缀
  base_non_img <- c("townsend_index_baseline","bmi_baseline","systolic_bp_baseline","diastolic_bp_baseline",
                    "diabetes_factor_baseline","smoking_factor_baseline","alcohol_factor_baseline",
                    "waist_circumference_baseline")
  # 根据最新要求：移除与运动/体力活动相关的协变量，不纳入敏感性分析
  extra_exercise <- character(0)
  # 基线（Instance.0）自然状况与血检：按关键词自动抓取
  i0_pattern <- "\\.{3,}Instance\\.0$"
  i0_conditions <- grep(paste0("(",
                               paste(c("blood.pressure","bp","hypertension","diabetes","smoking","alcohol","bmi","townsend"), collapse = ")|("),
                               ")"),
                        names(dat_main), ignore.case = TRUE, value = TRUE)
  i0_conditions <- i0_conditions[grepl(i0_pattern, i0_conditions)]
  i0_labs <- grep(paste0("(",
                         paste(c("calcium","cholesterol","hdl","ldl","triglyceride","glucose","hb?a1c","creatinine","urea","sodium","potassium","alt","ast","alp","bilirubin","crp","vitamin","hematocrit","hemoglobin","platelet","white.blood","red.blood","neutrophil","lymphocyte","monocyte","eosinophil","basophil"), collapse = ")|("),
                         ")"),
                  names(dat_main), ignore.case = TRUE, value = TRUE)
  i0_labs <- i0_labs[grepl(i0_pattern, i0_labs)]
  non_imaging_covariates <- unique(c(base_non_img, i0_conditions, i0_labs)) %>% intersect(names(dat_main))
  res_main <- run_psm_mi(
    dat_main, covariates, idp_cols,
    comparison = main_cmp,
    m = 1,
    impute_method = "median",
    mice_maxit = 0,
    mice_print = FALSE,
    psm_ratio = 2,
    psm_caliper = 0.2
  )
  
  # 诊断与可视化（主比较）
  diag_main <- compute_psm_balance(res_main$original_k1, res_main$matched_example, covariates, group_col = "group")
  fig_love_main <- plot_love_sci(
    diag_main$smd_plot_data,
    covariates = covariates,
    title_suffix = "(Latest PSM covariates)"
  )
  if (exists("fig_love_main") && !is.null(fig_love_main)) {
    save_pub(paste0("Figure_LovePlot_", main_cmp), fig_love_main, width = 10, height = 8, dpi = 300)
  } else {
    warning(sprintf("Love Plot 未生成：%s。请检查协变量与诊断数据。", main_cmp))
  }
  fig_ps_main <- plot_ps_distribution(res_main$mobj_example, res_main$original_k1, res_main$matched_example, group_col = "group", layout = "vertical", fixed_y = TRUE)
  if (exists("fig_ps_main") && !is.null(fig_ps_main)) {
    save_pub(paste0("Figure_PS_Distribution_", main_cmp), fig_ps_main, width = 7, height = 10, dpi = 300)
  }
  rm(fig_love_main, fig_ps_main)
  gc()
  if (!is.null(diag_main$p_values) && nrow(diag_main$p_values) > 0) {
    readr::write_csv(jama_round_df(as.data.frame(diag_main$p_values), digits = 3), table_out_path(file.path("Baseline","PSM"), paste0("PSM_PValues_", main_cmp, ".csv")))
  }
  readr::write_csv(jama_round_df(as.data.frame(diag_main$overall_stats), digits = 3), table_out_path(file.path("Baseline","PSM"), paste0("PSM_Balance_Overall_", main_cmp, ".csv")))
  cat(sprintf("✓ 保存匹配诊断与表格: %s\n", main_cmp))
  
  # 保存全局匹配队列及流程图（主比较）
  out_match_csv_main <- table_out_path(file.path("Cohort","Matched"), paste0("Final_Matched_Cohort_", main_cmp, ".csv"))
  out_match_rds_main <- matched_cohort_rds_path(main_cmp)
  join_cols <- intersect(c("eid","weights","subclass","distance"), names(res_main$matched_example))
  matched_join <- res_main$matched_example[, join_cols, drop = FALSE]
  # 关键修复：主队列的最终保存应为“匹配到的全量列”，
  # 而不是仅保存PSM子集（group+协变量+IDP）。否则会丢失认知/自然状况/血检等列。
  matched_full_main <- dplyr::inner_join(dat_main, matched_join, by = "eid")
  readr::write_csv(matched_full_main, out_match_csv_main)
  saveRDS(matched_full_main, out_match_rds_main)
  make_topjournal_baseline_table(
    standardize_variable_rules(dat_main),
    group_col = "group",
    out_file = table_out_path(file.path("Baseline","PSM"), paste0("Baseline_Table_Before_", main_cmp, ".csv")),
    case_label = "Ischemic",
    control_label = "Control"
  )
  make_topjournal_baseline_table(
    standardize_variable_rules(matched_full_main),
    group_col = "group",
    out_file = table_out_path(file.path("Baseline","PSM"), paste0("Baseline_Table_After_", main_cmp, ".csv")),
    case_label = "Ischemic",
    control_label = "Control"
  )
  # 基线三线表（高标准SCI版，仅显示PSM协变量，匹配后队列）
  make_three_line_table_sci(
    matched_full_main,
    group_col = "group",
    group_labels = c("Control","Case"),
    out_base_name = paste0("Baseline_SCI_Table_", main_cmp),
    exclude_vars = c(
      "Education (lower/medium/higher)",
      "Education (derived; lower/medium/higher)",
      "Education (raw; top-3 categories)",
      "Waist circumference (cm; mean±s.d. [range])",
      "Tobacco smoking (raw: never/former/current)",
      "Alcohol-intake frequency (raw: never/former/current)",
      "Diagnosed diabetes (raw yes/no)"
    )
  )
  summarize_age_distributions(matched_full_main, comparison = main_cmp, group_col = "group", group_labels = c("Control","Case"))
  
  flow_df_main <- summarize_psm_flow(res_main$original_k1, res_main$matched_example, group_col = "group")
  flow_df_main$comparison <- main_cmp
  saveRDS(flow_df_main, psm_flow_rds_path(main_cmp))
  fig_flow_main <- make_psm_flowchart(flow_df_main, cmp_label = main_cmp, layout = "vertical")
  save_pub(paste0("Figure1_PSM_Flowchart_", main_cmp), fig_flow_main, width = 7, height = 10, dpi = 300)
  rm(fig_flow_main)
  gc()
  cat(sprintf("✓ 保存匹配队列与流程数据/流程图: %s\n", main_cmp))
  
  # 汇总主比较的IDP结果
  pooled_all <- list()
  pooled_all[[main_cmp]] <- res_main$pooled_results
  
  # 子队列无需再次PSM：在全局匹配eid上复用并输出诊断和图表
  sub_comparisons <- c("mi_vs_control", "chronic_vs_control")
  pb_sub <- txtProgressBar(min = 1, max = length(sub_comparisons), style = 3)
  idx_sub <- 0
  for (cmp in sub_comparisons) {
    idx_sub <- idx_sub + 1
    setTxtProgressBar(pb_sub, idx_sub)
    
    dat_cmp <- build_comparison_df(ischemic_cohort_complete, cmp)
    # 原始数据添加距离（若存在），用于PS分布展示
    original_with_dist <- dplyr::left_join(dat_cmp, matched_join, by = "eid")
    # 子队列匹配结果：限制到全局匹配的eid
    matched_sub <- dplyr::inner_join(dat_cmp, matched_join, by = "eid")

    case_label <- dplyr::case_when(
      cmp == "mi_vs_control" ~ "MI",
      cmp == "chronic_vs_control" ~ "Chronic",
      TRUE ~ "Case"
    )
    control_label <- "Control"
    diag_sub <- compute_psm_balance(dat_cmp, matched_sub, covariates, group_col = "group")
    if (!is.null(diag_sub$p_values) && nrow(diag_sub$p_values) > 0) {
      readr::write_csv(jama_round_df(as.data.frame(diag_sub$p_values), digits = 3), table_out_path(file.path("Baseline","PSM"), paste0("PSM_PValues_", cmp, ".csv")))
    }
    readr::write_csv(jama_round_df(as.data.frame(diag_sub$overall_stats), digits = 3), table_out_path(file.path("Baseline","PSM"), paste0("PSM_Balance_Overall_", cmp, ".csv")))
    make_topjournal_baseline_table(
      standardize_variable_rules(dat_cmp),
      group_col = "group",
      out_file = table_out_path(file.path("Baseline","PSM"), paste0("Baseline_Table_Before_", cmp, ".csv")),
      case_label = case_label,
      control_label = control_label
    )
    make_topjournal_baseline_table(
      standardize_variable_rules(matched_sub),
      group_col = "group",
      out_file = table_out_path(file.path("Baseline","PSM"), paste0("Baseline_Table_After_", cmp, ".csv")),
      case_label = case_label,
      control_label = control_label
    )
    
  make_three_line_table_sci(
    matched_sub,
    group_col = "group",
    group_labels = c("Control", "Case"),
    out_base_name = paste0("Baseline_SCI_Table_", cmp),
    exclude_vars = c(
      "Education (lower/medium/higher)",
      "Education (derived; lower/medium/higher)",
      "Education (raw; top-3 categories)",
      "Waist circumference (cm; mean±s.d. [range])",
      "Tobacco smoking (raw: never/former/current)",
      "Alcohol-intake frequency (raw: never/former/current)",
      "Diagnosed diabetes (raw yes/no)"
    )
  )
    summarize_age_distributions(matched_sub, comparison = cmp, group_col = "group", group_labels = c("Control","Case"))
    cat(sprintf("✓ 子队列三线表已生成: %s\n", cmp))
    
    # 保存子队列匹配数据与流程图
    out_match_csv <- table_out_path(file.path("Cohort","Matched"), paste0("Final_Matched_Cohort_", cmp, ".csv"))
    out_match_rds <- matched_cohort_rds_path(cmp)
    readr::write_csv(matched_sub, out_match_csv)
    saveRDS(matched_sub, out_match_rds)
    
    # 子比较不生成流程图，也不保存PSM流程数据
    
    # 子队列IDP效应（不重新PSM，仅在matched_sub上估计）
    idp_res <- purrr::map_df(idp_cols, function(v) {
      vv <- matched_sub[[v]]
      if (!is.numeric(vv)) vv <- suppressWarnings(as.numeric(vv))
      if (all(is.na(vv))) return(NULL)
      vv <- scale(vv)
      df_fit <- data.frame(y = as.numeric(vv), group = matched_sub$group)
      fit <- stats::lm(y ~ group, data = df_fit)
      tidy <- broom::tidy(fit)
      eff <- dplyr::filter(tidy, term == "group") %>% dplyr::slice(1)
      data.frame(
        idp = v,
        estimate = eff$estimate,
        std.error = eff$std.error,
        statistic = eff$statistic,
        p.value = eff$p.value
      )
    })
    idp_res$imp <- 1
    pooled <- idp_res %>%
      dplyr::group_by(idp) %>%
      dplyr::summarise(
        m = dplyr::n(),
        pooled_coef = rubin_pool(estimate, std.error)$coef,
        pooled_se   = rubin_pool(estimate, std.error)$se,
        pooled_z    = rubin_pool(estimate, std.error)$z,
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        comparison = cmp,
        idp_category = identify_idp_category_en(idp)
      )
    pooled_all[[cmp]] <- pooled
  }
  close(pb_sub)
  
  # 直接使用四模型综合CSV作为Z分布的数据源，更稳健一致
  csv_path <- combined_four_model_results_csv()
  if (file.exists(csv_path)) {
    fallback <- readr::read_csv(csv_path, show_col_types = FALSE)
    # 稳健派生 pooled_z：优先使用无/有非影像Z；缺失时回退到基础Z与效应量
    z_with     <- if ("z_statistic_with_non_imaging"     %in% names(fallback)) suppressWarnings(as.numeric(fallback$z_statistic_with_non_imaging))     else rep(NA_real_, nrow(fallback))
    z_without  <- if ("z_statistic_without_non_imaging"  %in% names(fallback)) suppressWarnings(as.numeric(fallback$z_statistic_without_non_imaging))  else rep(NA_real_, nrow(fallback))
    cd_with    <- if ("cohens_d_with_non_imaging"        %in% names(fallback)) suppressWarnings(as.numeric(fallback$cohens_d_with_non_imaging))        else rep(NA_real_, nrow(fallback))
    cd_without <- if ("cohens_d_without_non_imaging"     %in% names(fallback)) suppressWarnings(as.numeric(fallback$cohens_d_without_non_imaging))     else rep(NA_real_, nrow(fallback))
    cd_raw     <- if ("cohens_d"                          %in% names(fallback)) suppressWarnings(as.numeric(fallback$cohens_d))                          else rep(NA_real_, nrow(fallback))
    fallback   <- fallback %>% dplyr::mutate(pooled_z = dplyr::coalesce(z_without, z_with, cd_with, cd_without, cd_raw))

    # 统一判断“脑影像”类型；若类型列不匹配，再用变量名模式兜底识别脑IDP
    name_brain_pattern <- stringr::regex(
      paste(c(
        "mean\\.thickness", "^volume\\.", "grey\\.white\\.contrast", "intensity",
        "(weighted\\.mean|mean)\\.fa", "(weighted\\.mean|mean)\\.md|diffusiv",
        "(weighted\\.mean|mean)\\.icvf", "(weighted\\.mean|mean)\\.od",
        "(weighted\\.mean|mean)\\.isovf", "rfmri.*amplitude", "connectivity"
      ), collapse = "|"), ignore_case = TRUE)
    fallback <- fallback %>% dplyr::mutate(
      is_brain_imaging = dplyr::coalesce(
        stringr::str_detect(as.character(variable_type), stringr::regex("brain\\s*_?\\s*imaging", ignore_case = TRUE)),
        FALSE
      )
    )
    fallback <- fallback %>% dplyr::mutate(
      is_brain_imaging = ifelse(is_brain_imaging, TRUE, stringr::str_detect(as.character(variable_name), name_brain_pattern))
    )

    pooled_results_all <- fallback %>%
      dplyr::filter(is_brain_imaging) %>%
      dplyr::transmute(
        idp_category = identify_idp_category_en(variable_name),
        pooled_z = pooled_z,
        comparison = model_comparison
      ) %>%
      tidyr::drop_na(idp_category, comparison, pooled_z)

    cat(sprintf("Z分布绘图数据行数（CSV）：%d\n", nrow(pooled_results_all)))
    rm(fallback); gc()
  } else {
    non_null <- pooled_all[!vapply(pooled_all, is.null, logical(1))]
    pooled_results_all <- if (length(non_null) > 0) dplyr::bind_rows(non_null) else tibble::tibble()
    if (nrow(pooled_results_all) == 0) {
      warning("未找到CSV且计算结果为空：IDP Z分布图可能为空。")
    }
  }
  
  # Z-statistic distribution (four comparisons), SCI/JAMA-style
  # Factorize IDP category with fixed English levels for consistent ordering
  idp_levels_master <- c(
    "Regional/Tissue Volume",
    "Cortical Thickness",
    "Cortical Surface Area",
    "Regional/Tissue Intensity",
    "Gray–White Contrast",
    "White Matter Hyperintensity Volume",
    "WM Tract FA",
    "WM Tract Diffusivity",
    "WM Tract ICVF",
    "WM Tract MO",
    "WM Tract OD",
    "WM Tract ISOVF",
    "rfMRI Node Amplitude",
    "rfMRI Connectivity",
    "Other/New IDP"
  )
  # 仅保留四个比较模型，移除其他来源导致的灰色点；若数据为空则跳过绘图
  if (!exists("pooled_results_all") || is.null(pooled_results_all)) pooled_results_all <- tibble::tibble()
  if (nrow(pooled_results_all) > 0) {
    pooled_results_all <- pooled_results_all %>%
      dplyr::mutate(
        comparison = dplyr::case_when(
          comparison %in% c("ischemic_vs_control","Ischemic_vs_Control","Ischemic vs Control") ~ "ischemic_vs_control",
          comparison %in% c("mi_vs_control","MI_vs_Control","Myocardial Infarction vs Control","Myocardial_Infarction_vs_Control") ~ "mi_vs_control",
          comparison %in% c("chronic_vs_control","Chronic_Ischemia_vs_Control","Chronic Ischemia vs Control") ~ "chronic_vs_control",
          TRUE ~ NA_character_
        )
      ) %>%
      dplyr::filter(!is.na(comparison)) %>%
      dplyr::mutate(
        idp_category = factor(idp_category, levels = {
          present <- unique(as.character(idp_category))
          c(
            idp_levels_master[idp_levels_master %in% present],
            sort(setdiff(present, idp_levels_master))
          )
        }),
        comparison   = factor(comparison, levels = c(
          "ischemic_vs_control",
          "mi_vs_control",
          "chronic_vs_control"
        ))
      )
    idp_x_labels <- levels(pooled_results_all$idp_category)
    n_idp_levels <- length(idp_x_labels)
  } else {
    warning("IDP Z分布绘图数据为空：跳过分布图绘制与保存。")
  }
  
  jama_cmp_colors <- c(
    "ischemic_vs_control" = "#0072B2",
    "mi_vs_control"       = "#D55E00",
    "chronic_vs_control"  = "#009E73"
  )
  comparison_labels <- c(
    "ischemic_vs_control" = "Ischemic vs Control",
    "mi_vs_control"       = "Myocardial Infarction vs Control",
    "chronic_vs_control"  = "Chronic Ischemia vs Control"
  )
  if (nrow(pooled_results_all) > 0) {
  p_z <- ggplot(pooled_results_all,
               aes(x = idp_category, y = pooled_z, color = comparison)) +
    ggplot2::geom_point(
      position = ggplot2::position_jitterdodge(jitter.width = 0.20, jitter.height = 0, dodge.width = 0.70),
      alpha = 0.85, size = 1.8
    ) +
    # 类目间分隔线（更易感知各IDP类目下的四模型并列）
    ggplot2::geom_vline(xintercept = seq(1.5, max(n_idp_levels - 0.5, 1.5), by = 1),
                        color = "#e6e6e6", linetype = "dotted", linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    ggplot2::scale_color_manual(values = jama_cmp_colors,
                                breaks = names(comparison_labels),
                                labels = unname(comparison_labels),
                                name = "Comparison") +
    ggplot2::scale_x_discrete(limits = idp_x_labels, labels = stringr::str_wrap(idp_x_labels, width = 18), drop = FALSE) +
    ggplot2::scale_y_continuous(limits = c(-6, 6), breaks = seq(-6, 6, 2)) +
    ggplot2::labs(
      title = "Z-Statistic Distribution of Brain Imaging IDPs (Three Comparisons)",
      x = "IDP Category",
      y = "Model Z-Statistic"
    ) +
    theme_jama_refined() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, vjust = 1)
    )
  save_pub("IDP_Z_Distribution_All_Models", p_z, width = 14, height = 6, dpi = 300)
  cat("✓ Saved IDP Z-Distribution plot: PNG/PDF\n")
  rm(p_z); gc()
  } else {
    cat("⚠️ 跳过IDP Z分布图：无可用数据。\n")
  }
  
  cat("✓ 总队列一次PSM + 子队列复用全局匹配完成。\n")
}

# =================== 检查并加载分析数据列表（如果跳过了 PSM 步骤） ===================
if (!exists("analysis_data_list") || length(analysis_data_list) == 0) {
  cat("检测到内存中缺少 analysis_data_list，尝试从磁盘加载匹配队列...\n")
  analysis_data_list <- list()
  
  models <- c("ischemic_vs_control", "mi_vs_control", "chronic_vs_control")
  for (model in models) {
    rds_file <- resolve_existing_file(c(
      matched_cohort_rds_path(model),
      paste0("Final_Matched_Cohort_", model, ".rds")
    ))
    if (!is.na(rds_file) && nzchar(rds_file)) {
      cat(sprintf("  加载 %s: %s\n", model, rds_file))
      analysis_data_list[[model]] <- readRDS(rds_file)
    }
  }
  
  if (length(analysis_data_list) == 0) {
    warning("⚠️ 未能加载任何分析数据列表 (analysis_data_list)，后续步骤可能会失败。请检查是否已运行 PSM 匹配步骤或文件是否存在。")
  } else {
    cat(sprintf("✓ 已恢复 %d 个模型的分析数据。\n", length(analysis_data_list)))
  }
}

if (RUN_IDP_ANALYSIS) {
# =================== 按照四模型设计的缺血性心肌病IDP纵向变化Z值分析 ===================
cat("=== 按照四模型设计的缺血性心肌病IDP纵向变化Z值分析 ===\n")

# =================== 第一步：读取数据并识别变量 ===================
cat("第一步：读取数据并识别变量...\n")

if(!exists("matched_cohort")) {
  # 优先使用内存中的主匹配队列（在前面的 PSM 步骤生成）
  if (exists("matched_full_main", inherits = TRUE)) {
    matched_cohort <- matched_full_main
    message("已使用内存中的主匹配队列作为四模型分析基础。")
  } else {
    # 回退：在工作目录下寻找任意已保存的主匹配队列 RDS
    candidates <- c(
      matched_cohort_rds_path("ischemic_vs_control"),
      matched_cohort_rds_path("mi_vs_control"),
      matched_cohort_rds_path("chronic_vs_control"),
      "Final_Matched_Cohort_ischemic_vs_control.rds",
      "Final_Matched_Cohort_mi_vs_control.rds",
      "Final_Matched_Cohort_chronic_vs_control.rds"
    )
    cand_path <- candidates[file.exists(candidates)]
    if (length(cand_path) > 0) {
      matched_cohort <- readRDS(cand_path[[1]])
      message(sprintf("已加载匹配队列：%s", cand_path[[1]]))
    } else if (exists("res_main") && !is.null(res_main$matched_example)) {
      matched_cohort <- res_main$matched_example
      message("已使用 res_main$matched_example 作为临时匹配队列。")
    } else {
      stop("匹配队列不可用，请先运行主比较的PSM以生成结果或确认工作目录。")
    }
  }
}

if (!exists("matched_cohort") || is.null(matched_cohort)) {
  stop("匹配队列不可用，无法继续。")
}
cat("匹配队列样本数:", nrow(matched_cohort), "\n")

# 认知和脑影像变量识别模式
cognitive_patterns <- c(
  "Duration.to.complete", "Total.errors", "symbol.digit", "Mean.time.to.correctly.identify",
  "Fluid.intelligence", "Maximum.digits", "Number.of.correct.matches", 
  "Number.of.incorrect.matches", "Time.to.complete.round", "Pairs.matching",
  "Trail.making", "Tower.rearranging", "Symbol.digit.substitution",
  "Reaction.time", "reasoning", "cognitive", "memory", "attention", "executive",
  "numeric.memory", "prospective.memory", "intelligence.score"
)

brain_patterns <- c(
  "Volume.of","Mean.thickness","Area.of","Mean.FA","Mean.MD","Mean.MO","Mean.OD","Mean.ICVF",
  "Mean.ISOVF","Mean.CBF","Mean.ATT","Mean.L1","Mean.L2","Mean.L3","Mean.absolute","Mean.relative",
  "Median.BOLD","Median.z.statistic","Total.volume","T2.FLAIR","Arterial.spin","Volumetric.scaling",
  "Inverted.contrast","Inverted.signal","Amount.of.warping","T1.surface","T1.structural","Functional.brain",
  "Discrepancy.between","Weighted.mean","X90th.","Grey.white.contrast","percentile.of.BOLD","Median.BOLD",
  "Mean.tfMRI","volume.of.white.matter","volume.of.grey.matter","volume.of.brain"
)

# 模式去重，减少冗余
cognitive_patterns <- unique(cognitive_patterns)
brain_patterns <- unique(brain_patterns)

# 通用识别函数：按包含/排除模式一次性识别并去重
identify_by_patterns <- function(all_var_names, include_patterns, exclude_patterns = character(0), var_type = NULL) {
  if (length(include_patterns) == 0) return(character(0))
  inc_regex <- paste0("(", paste(include_patterns, collapse = ")|("), ")")
  matches <- all_var_names[grepl(inc_regex, all_var_names, ignore.case = TRUE)]
  if (length(exclude_patterns) > 0) {
    exc_regex <- paste0("(", paste(exclude_patterns, collapse = ")|("), ")")
    matches <- matches[!grepl(exc_regex, matches, ignore.case = TRUE)]
  }
  matches <- unique(matches)
  if (!is.null(var_type)) {
    cat(sprintf("%s变量: %d 个\n", var_type, length(matches)))
  }
  matches
}

# 识别变量：统一整合四个比较的最终匹配队列的变量名，减少遗漏
get_all_variable_names <- function(matched_cohort) {
  vec <- names(matched_cohort)
  cmps <- c("ischemic_vs_control", "mi_vs_control", "chronic_vs_control")
  for (cmp in cmps) {
    rds_path <- resolve_existing_file(c(
      matched_cohort_rds_path(cmp),
      paste0("Final_Matched_Cohort_", cmp, ".rds")
    ))
    if (!is.na(rds_path) && nzchar(rds_path)) {
      dat <- tryCatch(readRDS(rds_path), error = function(e) NULL)
      if (!is.null(dat)) {
        vec <- c(vec, names(dat))
      }
    }
  }
  unique(vec)
}
all_variable_names <- get_all_variable_names(matched_cohort)

# 来自“纵向分析2.R”的认知识别增强：纯认知变量（排除脑影像关键词）
identify_cognitive_variables_pure <- function(all_names) {
  cognitive_test_keywords <- c(
    "fluid.intelligence", "intelligence.score",
    "mean.time.to.correctly.identify", "reaction.time",
    "time.to.complete", "pairs.matching",
    "trail.making", "digit.symbol", "verbal", "numeric.memory",
    "processing.speed", "attention", "working.memory",
    "executive.function", "reasoning"
  )
  brain_imaging_keywords_exclude <- c(
    "volume.of", "cortical", "subcortical", "white.matter", "grey.matter",
    "thickness", "area", "FA", "MD", "BOLD", "connectivity", "tfMRI"
  )
  identify_by_patterns(all_names, cognitive_test_keywords, brain_imaging_keywords_exclude)
}

# 动态补充：基于最终队列变量名出现的候选关键词，自动完善模式
original_cognitive_patterns <- cognitive_patterns
original_brain_patterns <- brain_patterns

augment_patterns <- function(patterns, candidates, names_vec) {
  hits <- candidates[sapply(candidates, function(k) any(grepl(k, names_vec, ignore.case = TRUE)))]
  unique(c(patterns, hits))
}

cognitive_candidates <- c(
  "fluid.intelligence","intelligence.score","mean.time.to.correctly.identify","reaction.time",
  "time.to.complete","pairs.matching","trail.making","digit.symbol","verbal","numeric.memory",
  "prospective.memory","maximum.digits","symbol.digit","symbol.digit.substitution",
  "total.errors","duration.to.complete","number.of.correct.matches","number.of.incorrect.matches"
)

brain_candidates <- c(
  "Mean.RD","Mean.AD","Mean.R2s","Mean.FC","R2s","SWI","Total.surface.area",
  "Cortical.thickness","Volume.of.grey.matter","Volume.of.white.matter","Volume.of.brain",
  "Median.z.statistic","percentile.of.BOLD","Mean.tfMRI"
)

cognitive_patterns <- augment_patterns(cognitive_patterns, cognitive_candidates, all_variable_names)
brain_patterns <- augment_patterns(brain_patterns, brain_candidates, all_variable_names)

added_cog <- setdiff(cognitive_patterns, original_cognitive_patterns)
added_brain <- setdiff(brain_patterns, original_brain_patterns)
cat(sprintf("动态补充认知关键词: %d 个\n", length(added_cog)))
if(length(added_cog) > 0) cat(paste0("  + ", paste(added_cog, collapse = ", "), "\n"))
cat(sprintf("动态补充脑影像关键词: %d 个\n", length(added_brain)))
if(length(added_brain) > 0) cat(paste0("  + ", paste(added_brain, collapse = ", "), "\n"))

# 先基于模式识别两类变量
cognitive_vars_pattern <- identify_by_patterns(all_variable_names, cognitive_patterns, var_type = "认知")
brain_imaging_vars <- identify_by_patterns(all_variable_names, brain_patterns, var_type = "脑影像")

# 使用纵向分析2.R的增强策略识别纯认知变量并与模式识别结果合并
cognitive_vars_pure <- identify_cognitive_variables_pure(all_variable_names)
cat(sprintf("纯认知关键词识别: %d 个\n", length(cognitive_vars_pure)))

cognitive_vars <- unique(c(cognitive_vars_pattern, cognitive_vars_pure))
cat(sprintf("合并后认知变量: %d 个\n", length(cognitive_vars)))

# 处理重叠变量
overlapping_vars <- intersect(cognitive_vars, brain_imaging_vars)
if(length(overlapping_vars) > 0) {
  brain_imaging_vars <- setdiff(brain_imaging_vars, overlapping_vars)
}

# =================== 第二步：创建Instance 2-3变量配对 ===================
cat("\n第二步：创建Instance 2-3变量配对...\n")

create_variable_pairs <- function(cognitive_vars, brain_vars) {
  pairs <- list()
  
  # 认知变量配对
  if(length(cognitive_vars) > 0) {
    cog_inst2 <- grep("Instance(\\.|_|\\s)?2", cognitive_vars, value = TRUE)
    for(var2 in cog_inst2) {
      var2_root <- sub("(\\.x|\\.y)$", "", var2)
      var3_raw <- gsub("Instance(\\.|_|\\s)?2", "Instance.3", var2_root)
      base_name <- gsub("\\.\\.\\.Instance(\\.|_|\\s)?2.*$|Instance(\\.|_|\\s)?2.*$", "", var2_root)
      base_name <- gsub("\\.$", "", base_name)
      pairs[[length(pairs) + 1]] <- list(
        idp1_var = var2_root, idp2_var = var3_raw, base_name = base_name, variable_type = "Cognitive"
      )
    }
  }
  
  # 脑影像变量配对
  if(length(brain_vars) > 0) {
    brain_inst2 <- grep("Instance(\\.|_|\\s)?2", brain_vars, value = TRUE)
    for(var2 in brain_inst2) {
      var2_root <- sub("(\\.x|\\.y)$", "", var2)
      var3_raw <- gsub("Instance(\\.|_|\\s)?2", "Instance.3", var2_root)
      base_name <- gsub("\\.\\.\\.Instance(\\.|_|\\s)?2.*$|Instance(\\.|_|\\s)?2.*$", "", var2_root)
      base_name <- gsub("\\.$", "", base_name)
      pairs[[length(pairs) + 1]] <- list(
        idp1_var = var2_root, idp2_var = var3_raw, base_name = base_name, variable_type = "Brain_Imaging"
      )
    }
  }
  
  return(pairs)
}

variable_pairs <- create_variable_pairs(cognitive_vars, brain_imaging_vars)
cognitive_pairs_count <- sum(sapply(variable_pairs, function(x) x$variable_type == "Cognitive"))
brain_pairs_count <- sum(sapply(variable_pairs, function(x) x$variable_type == "Brain_Imaging"))

cat("成功配对变量:", length(variable_pairs), "个\n")
cat("- 认知变量配对:", cognitive_pairs_count, "个\n")
cat("- 脑影像变量配对:", brain_pairs_count, "个\n")

# =================== 第三步：准备四种模型的分析数据 ===================
cat("\n第三步：准备四种模型的分析数据...\n")

# 计算年龄变量（健壮处理，避免0长度替换错误）
compute_age_instances <- function(df) {
  n <- nrow(df)
  # 保证关键列存在且长度匹配
  if (!("baseline_age" %in% names(df)) && ("age_at_recruitment" %in% names(df))) {
    df$baseline_age <- df$age_at_recruitment
  }
  if (!("age_at_recruitment" %in% names(df)) && ("baseline_age" %in% names(df))) {
    df$age_at_recruitment <- df$baseline_age
  }
  if (!("followup_age" %in% names(df)) || length(df$followup_age) != n) {
    df$followup_age <- rep(NA_real_, n)
  }
  if (!("followup_years" %in% names(df)) || length(df$followup_years) != n) {
    df$followup_years <- rep(NA_real_, n)
  }
  
  # Age at instance 2：优先baseline_age
  if ("baseline_age" %in% names(df)) {
    df$age_at_instance2 <- df$baseline_age
  } else {
    df$age_at_instance2 <- df$age_at_recruitment
  }
  
  # Age at instance 3：优先followup_age，其次recruitment+followup_years，否则+2年
  df$age_at_instance3 <- df$age_at_recruitment + 2
  has_fup_years <- !is.na(df$followup_years)
  df$age_at_instance3[has_fup_years] <- df$age_at_recruitment[has_fup_years] + df$followup_years[has_fup_years]
  has_fup_age <- !is.na(df$followup_age)
  df$age_at_instance3[has_fup_age] <- df$followup_age[has_fup_age]
  
  df
}

matched_cohort <- compute_age_instances(matched_cohort)

canonicalize_instance_columns <- function(df) {
  n <- nrow(df)
  nm <- names(df)
  norm <- function(s) {
    s <- sub("(\\.x|\\.y)$", "", s)
    s <- gsub("Instance[_ ]([23])", "Instance.\\1", s, perl = TRUE)
    s
  }
  keys <- sapply(nm, norm)
  groups <- split(nm, keys)
  for (b in names(groups)) {
    cols <- groups[[b]]
    v <- rep(NA_real_, n)
    for (c in cols) {
      v <- dplyr::coalesce(v, suppressWarnings(as.numeric(df[[c]])))
    }
    df[[b]] <- v
  }
  df
}

matched_cohort <- canonicalize_instance_columns(matched_cohort)

 dedupe_instance_columns <- function(df) {
   nm <- names(df)
   norm <- function(s) {
     s <- sub("(\\.x|\\.y)$", "", s)
     s <- gsub("Instance[_ ]([23])", "Instance.\\1", s, perl = TRUE)
     s
   }
   keys <- sapply(nm, norm)
   groups <- split(nm, keys)
   drop <- unlist(Map(function(b, cols) setdiff(cols, b), names(groups), groups))
   keep <- setdiff(nm, drop)
   df[, keep, drop = FALSE]
 }

 enforce_two_timepoint_presence <- function(df) {
   nm <- names(df)
   inst2 <- nm[grepl("Instance\\.2$", nm)]
   inst3 <- nm[grepl("Instance\\.3$", nm)]
   if (length(inst2) == 0 || length(inst3) == 0) return(df)
   to_num <- function(x) suppressWarnings(as.numeric(x))
   m2 <- as.data.frame(lapply(df[, inst2, drop = FALSE], to_num))
   m3 <- as.data.frame(lapply(df[, inst3, drop = FALSE], to_num))
   has2 <- rowSums(!is.na(m2)) > 0
   has3 <- rowSums(!is.na(m3)) > 0
   df[has2 & has3, , drop = FALSE]
 }

 matched_cohort <- dedupe_instance_columns(matched_cohort)
 matched_cohort <- enforce_two_timepoint_presence(matched_cohort)

cat("提示：四模型分析数据将于第五步按每个比较的匹配队列重建。\n")

# =================== 第四步：Z统计量分析函数 ===================
cat("\n第四步：定义Z统计量分析函数...\n")

z_statistics_four_model_analysis <- function(data, pairs, model_name,
                                             drop_demographics = FALSE,
                                             non_imaging_covariates = NULL,
                                             parallel = TRUE,
                                             n_cores = NULL) {
  
  results <- list()
  cat("分析模型:", model_name, "- 变量对数:", length(pairs), "\n")
  
  if(length(pairs) == 0) {
    cat("无可分析的变量对\n")
    return(NULL)
  }
  
  # 单个pair分析函数：供串行/并行重用
  analyze_one_pair <- function(pair, data, drop_demographics, non_imaging_covariates, model_name) {
    if(!all(c(pair$idp1_var, pair$idp2_var) %in% names(data))) return(NULL)
    idp1_data <- data[[pair$idp1_var]]
    idp2_data <- data[[pair$idp2_var]]
    if(!is.numeric(idp1_data) || !is.numeric(idp2_data)) return(NULL)
    complete_cases <- !is.na(idp1_data) & !is.na(idp2_data)
    if(sum(complete_cases) < 3) return(NULL)
    analysis_data <- data[complete_cases, ]
    analysis_data$Age2 <- dplyr::coalesce(analysis_data$Age2,
                                          analysis_data$age_at_instance3,
                                          analysis_data$age_at_instance2 + analysis_data$followup_years,
                                          analysis_data$baseline_age + analysis_data$followup_years,
                                          analysis_data$age_at_instance2,
                                          analysis_data$baseline_age,
                                          analysis_data$age_at_recruitment)
    analysis_data$Age2 <- ifelse(is.na(analysis_data$Age2) | !is.finite(analysis_data$Age2), stats::median(analysis_data$Age2, na.rm = TRUE), analysis_data$Age2)
    if (!("case_control_binary" %in% names(analysis_data))) analysis_data$case_control_binary <- NA_real_
    if ("comparison_group" %in% names(analysis_data)) {
      if (exists("cmp_labels") && model_name %in% names(cmp_labels)) {
        pos_lab <- cmp_labels[[model_name]][1]
        analysis_data$case_control_binary <- dplyr::coalesce(analysis_data$case_control_binary, as.numeric(analysis_data$comparison_group == pos_lab))
      } else {
        labs <- unique(analysis_data$comparison_group)
        analysis_data$case_control_binary <- dplyr::coalesce(analysis_data$case_control_binary, as.numeric(analysis_data$comparison_group == labs[1]))
      }
    } else if ("group" %in% names(analysis_data)) {
      analysis_data$case_control_binary <- dplyr::coalesce(analysis_data$case_control_binary, as.numeric(as.character(analysis_data$group)) == 1)
    }
    analysis_data <- analysis_data[!is.na(analysis_data$case_control_binary) & !is.na(analysis_data$Age2), ]
    analysis_data$IDP1 <- analysis_data[[pair$idp1_var]]
    analysis_data$IDP2 <- analysis_data[[pair$idp2_var]]
    case_control_demeaned <- analysis_data$case_control_binary - mean(analysis_data$case_control_binary)
    aging_multiplier <- 10^(analysis_data$Age2 * 0.0524 - 3.27)
    analysis_data$Case_vs_Control <- case_control_demeaned * aging_multiplier
    confounds_basic <- c("age_diff_linear", "age_sq_diff_linear")
    if (!isTRUE(drop_demographics)) {
      if("sex_factor" %in% names(analysis_data)) confounds_basic <- c(confounds_basic, "sex_factor")
      if("ethnicity_factor" %in% names(analysis_data)) confounds_basic <- c(confounds_basic, "ethnicity_factor")
    }
    out_df <- NULL
    tryCatch({
      formula_basic <- "IDP2 ~ Case_vs_Control + IDP1"
      available_confounds <- confounds_basic[confounds_basic %in% names(analysis_data)]
      if(length(available_confounds) > 0) {
        formula_basic <- paste(formula_basic, "+", paste(available_confounds, collapse = " + "))
      }
      model_basic <- lm(as.formula(formula_basic), data = analysis_data)
      coef_table_basic <- coef(summary(model_basic))
      case_control_row <- which(rownames(coef_table_basic) == "Case_vs_Control")
      if(length(case_control_row) > 0) {
        z_without_non_imaging <- coef_table_basic[case_control_row, "t value"]
        p_without_non_imaging <- coef_table_basic[case_control_row, "Pr(>|t|)"]
        beta_without_non_imaging <- coef_table_basic[case_control_row, "Estimate"]
        se_without_non_imaging <- coef_table_basic[case_control_row, "Std. Error"]
        pooled_sd <- sqrt(var(analysis_data$IDP2, na.rm = TRUE))
        cohens_d <- beta_without_non_imaging / pooled_sd
        group_stats <- analysis_data %>%
          dplyr::group_by(case_control_binary) %>%
          dplyr::summarise(
            n = dplyr::n(),
            idp1_mean = mean(IDP1, na.rm = TRUE),
            idp1_sd = sd(IDP1, na.rm = TRUE),
            idp2_mean = mean(IDP2, na.rm = TRUE),
            idp2_sd = sd(IDP2, na.rm = TRUE),
            age_mean = mean(Age2, na.rm = TRUE),
            age_sd = sd(Age2, na.rm = TRUE),
            .groups = "drop"
          )
        significance_level <- dplyr::case_when(
          abs(z_without_non_imaging) >= 3.29 ~ "p<0.001",
          abs(z_without_non_imaging) >= 2.58 ~ "p<0.01",
          abs(z_without_non_imaging) >= 1.96 ~ "p<0.05",
          abs(z_without_non_imaging) >= 1.65 ~ "p<0.10",
          TRUE ~ "p>=0.10"
        )
        independent_effect <- (p_without_non_imaging < 0.05)
        z_with_non_imaging <- NA_real_
        p_with_non_imaging <- NA_real_
        beta_with_non_imaging <- NA_real_
        se_with_non_imaging <- NA_real_
        cohens_d_with_non_imaging <- NA_real_
        available_non_img <- NULL
        if (!is.null(non_imaging_covariates) && length(non_imaging_covariates) > 0) {
          available_non_img <- non_imaging_covariates[non_imaging_covariates %in% names(analysis_data)]
        }
        if (!is.null(available_non_img) && length(available_non_img) > 0) {
          formula_with_nonimg <- paste(formula_basic, "+", paste(available_non_img, collapse = " + "))
          model_nonimg <- tryCatch(lm(as.formula(formula_with_nonimg), data = analysis_data), error = function(e) NULL)
          if (!is.null(model_nonimg)) {
            ct_nonimg <- coef(summary(model_nonimg))
            rr_nonimg <- which(rownames(ct_nonimg) == "Case_vs_Control")
            if (length(rr_nonimg) > 0) {
              z_with_non_imaging <- ct_nonimg[rr_nonimg, "t value"]
              p_with_non_imaging <- ct_nonimg[rr_nonimg, "Pr(>|t|)"]
              beta_with_non_imaging <- ct_nonimg[rr_nonimg, "Estimate"]
              se_with_non_imaging <- ct_nonimg[rr_nonimg, "Std. Error"]
              pooled_sd2 <- sqrt(var(analysis_data$IDP2, na.rm = TRUE))
              cohens_d_with_non_imaging <- beta_with_non_imaging / pooled_sd2
            }
          }
        }
        out_df <- data.frame(
          model_comparison = model_name,
          variable_name = substr(pair$base_name, 1, 100),
          variable_type = pair$variable_type,
          idp1_variable = pair$idp1_var,
          idp2_variable = pair$idp2_var,
          n_total = nrow(analysis_data),
          n_cases = sum(analysis_data$case_control_binary == 1),
          n_controls = sum(analysis_data$case_control_binary == 0),
          case_prevalence = mean(analysis_data$case_control_binary),
          z_statistic_without_non_imaging = round(z_without_non_imaging, 4),
          p_value_without_non_imaging = p_without_non_imaging,
          p_formatted = ifelse(p_without_non_imaging < 0.001, "<0.001", sprintf("%.4f", p_without_non_imaging)),
          significance_level = significance_level,
          z_statistic_with_non_imaging = round(z_with_non_imaging, 4),
          p_value_with_non_imaging = p_with_non_imaging,
          beta_estimate_with_non_imaging = round(beta_with_non_imaging, 6),
          standard_error_with_non_imaging = round(se_with_non_imaging, 6),
          beta_estimate = round(beta_without_non_imaging, 6),
          standard_error = round(se_without_non_imaging, 6),
          cohens_d = round(cohens_d, 4),
          cohens_d_with_non_imaging = round(cohens_d_with_non_imaging, 4),
          effect_size_category = dplyr::case_when(
            abs(cohens_d) < 0.2 ~ "Small",
            abs(cohens_d) < 0.5 ~ "Moderate",
            abs(cohens_d) >= 0.5 ~ "Large"
          ),
          independent_effect = independent_effect,
          independent_of_confounds = independent_effect,
          cases_idp1_mean = ifelse(nrow(group_stats) >= 2, round(group_stats$idp1_mean[group_stats$case_control_binary == 1][1], 4), NA),
          controls_idp1_mean = ifelse(nrow(group_stats) >= 2, round(group_stats$idp1_mean[group_stats$case_control_binary == 0][1], 4), NA),
          cases_idp2_mean = ifelse(nrow(group_stats) >= 2, round(group_stats$idp2_mean[group_stats$case_control_binary == 1][1], 4), NA),
          controls_idp2_mean = ifelse(nrow(group_stats) >= 2, round(group_stats$idp2_mean[group_stats$case_control_binary == 0][1], 4), NA),
          cases_age_mean = ifelse(nrow(group_stats) >= 2, round(group_stats$age_mean[group_stats$case_control_binary == 1][1], 2), NA),
          controls_age_mean = ifelse(nrow(group_stats) >= 2, round(group_stats$age_mean[group_stats$case_control_binary == 0][1], 2), NA),
          model_formula = formula_basic,
          model_formula_with_non_imaging = if (exists("formula_with_nonimg")) formula_with_nonimg else NA_character_,
          stringsAsFactors = FALSE
        )
      }
      out_df
    }, error = function(e) NULL)
  }

  # 决定是否并行：若缺少 foreach/doParallel 则回退为串行
  use_par <- isTRUE(parallel) && requireNamespace("foreach", quietly = TRUE) && requireNamespace("doParallel", quietly = TRUE)
  if (use_par) {
    # 附加 foreach 以启用 %dopar% 操作符
    suppressPackageStartupMessages(library(foreach))
    suppressPackageStartupMessages(library(doParallel))
    cores <- if (is.null(n_cores) || !is.numeric(n_cores) || n_cores < 1) max(1L, parallel::detectCores() - 1L) else as.integer(n_cores)
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    on.exit({ try(parallel::stopCluster(cl), silent = TRUE) }, add = TRUE)
    results <- foreach::foreach(i = seq_along(pairs), .combine = 'c', .multicombine = TRUE,
                                 .packages = c('dplyr','tidyr','stringr'),
                                 .export = c('analyze_one_pair')) %dopar% {
      list(analyze_one_pair(pairs[[i]], data, drop_demographics, non_imaging_covariates, model_name))
    }
  } else {
    pb <- txtProgressBar(min = 0, max = length(pairs), style = 3)
    for (i in seq_along(pairs)) {
      setTxtProgressBar(pb, i)
      results[[i]] <- analyze_one_pair(pairs[[i]], data, drop_demographics, non_imaging_covariates, model_name)
    }
    close(pb)
  }
  
  if(length(results) > 0) {
    valid_results <- results[!sapply(results, is.null)]
    if (length(valid_results) == 0) {
      return(NULL)
    }
    results_df <- do.call(rbind, valid_results)
    if(nrow(results_df) > 0) {
      # 按Z统计量绝对值排序
      results_df <- results_df[order(abs(results_df$z_statistic_without_non_imaging), decreasing = TRUE), ]
      rownames(results_df) <- NULL
    }
    return(results_df)
  }
  return(NULL)
}

# 置换-FWE计算：对同一模型内的所有检验进行 max-|t| 置换，返回每个结果的 FWE p 值
compute_fwe_pvalues_for_results <- function(data_model, results_df, drop_demographics = FALSE, n_perm = 50, seed = 202311) {
  if (is.null(results_df) || nrow(results_df) == 0) return(rep(NA_real_, 0))
  set.seed(seed)
  
  # 全局老化乘子（按行复用，避免反复计算）
  aging_mult_global <- 10^(data_model$Age2 * 0.0524 - 3.27)
  
  # 为每个结果构建最小元数据，避免存储整数据帧
  build_item <- function(row) {
    idp1 <- results_df$idp1_variable[row]
    idp2 <- results_df$idp2_variable[row]
    if (!(idp1 %in% names(data_model)) || !(idp2 %in% names(data_model))) return(NULL)
    idp1_data <- data_model[[idp1]]
    idp2_data <- data_model[[idp2]]
    cc <- !is.na(idp1_data) & !is.na(idp2_data) & !is.na(data_model$case_control_binary) & !is.na(data_model$Age2)
    if (sum(cc) < 3) return(NULL)
    row_idx <- which(cc)
    
    # 可用混杂项
    conf_basic <- c("age_diff_linear", "age_sq_diff_linear")
    if (!isTRUE(drop_demographics)) {
      if ("sex_factor" %in% names(data_model)) conf_basic <- c(conf_basic, "sex_factor")
      if ("ethnicity_factor" %in% names(data_model)) conf_basic <- c(conf_basic, "ethnicity_factor")
    }
    avail_conf <- conf_basic[conf_basic %in% names(data_model)]
    
    # 构建一次性最小数据子集并计算观测 t
    df_min <- data.frame(
      IDP2 = data_model[row_idx, idp2],
      IDP1 = data_model[row_idx, idp1],
      Case_vs_Control = (data_model$case_control_binary[row_idx] - mean(data_model$case_control_binary[row_idx])) * aging_mult_global[row_idx],
      stringsAsFactors = FALSE
    )
    if (length(avail_conf) > 0) {
      df_min[avail_conf] <- data_model[row_idx, avail_conf, drop = FALSE]
    }
    form <- "IDP2 ~ Case_vs_Control + IDP1"
    if (length(avail_conf) > 0) form <- paste0(form, " + ", paste(avail_conf, collapse = " + "))
    
    fit0 <- tryCatch(lm(as.formula(form), data = df_min), error = function(e) NULL)
    if (is.null(fit0)) return(NULL)
    ct <- coef(summary(fit0))
    rr <- which(rownames(ct) == "Case_vs_Control")
    if (length(rr) == 0) return(NULL)
    t_obs <- as.numeric(ct[rr, "t value"])
    
    list(row_idx = row_idx, idp1 = idp1, idp2 = idp2, avail_conf = avail_conf, formula = form, t_obs = t_obs)
  }
  
  items <- lapply(seq_len(nrow(results_df)), build_item)
  keep <- !sapply(items, is.null)
  if (!any(keep)) return(rep(NA_real_, nrow(results_df)))
  items_idx <- which(keep)
  
  # 观测 |t|
  obs_abs_t <- sapply(items[items_idx], function(it) abs(it$t_obs))
  exceed_count <- numeric(length(items_idx))
  
  # 置换进度条
  cat(sprintf("[FWE] compute_fwe_pvalues_for_results: n_perm=%d, items=%d\n", n_perm, length(items_idx)))
  pb_perm <- txtProgressBar(min = 1, max = n_perm, style = 3)

  # 置换循环：按次迭代更新最大 |t|，避免大向量分配
  for (b in seq_len(n_perm)) {
    max_abs_t <- 0
    for (j in seq_along(items_idx)) {
      it <- items[[items_idx[j]]]
      ridx <- it$row_idx
      # 仅构建必要列的最小数据帧
      perm_cc <- sample(data_model$case_control_binary[ridx])
      df_min <- data.frame(
        IDP2 = data_model[ridx, it$idp2],
        IDP1 = data_model[ridx, it$idp1],
        Case_vs_Control = (perm_cc - mean(perm_cc)) * aging_mult_global[ridx],
        stringsAsFactors = FALSE
      )
      if (length(it$avail_conf) > 0) {
        df_min[it$avail_conf] <- data_model[ridx, it$avail_conf, drop = FALSE]
      }
      fitp <- tryCatch(lm(as.formula(it$formula), data = df_min), error = function(e) NULL)
      if (is.null(fitp)) next
      ctp <- coef(summary(fitp))
      rrp <- which(rownames(ctp) == "Case_vs_Control")
      if (length(rrp) > 0) {
        tv <- abs(as.numeric(ctp[rrp, "t value"]))
        if (is.finite(tv) && tv > max_abs_t) max_abs_t <- tv
      }
    }
      if (is.finite(max_abs_t)) {
        exceed_count <- exceed_count + as.numeric(max_abs_t >= obs_abs_t)
      }
      # 可选：释放临时对象
      rm(max_abs_t)
      setTxtProgressBar(pb_perm, b)
    }
  close(pb_perm)
  
  # FWE p 值（加一修正）
  p_fwe <- rep(NA_real_, nrow(results_df))
  for (j in seq_along(items_idx)) {
    idx <- items_idx[j]
    p_fwe[idx] <- (1 + exceed_count[j]) / (n_perm + 1)
  }
  p_fwe
}

# 基于候选非成像变量的Z值衰减分析（逐模型、逐IDP对）
z_attenuation_by_non_imaging <- function(data, pairs, model_name, drop_demographics = FALSE, non_imaging_covariates = NULL) {
  if(length(pairs) == 0) return(NULL)
  out_rows <- list()
  if (!("diastolic_bp" %in% names(data)) && ("diastolic_bp_baseline" %in% names(data))) {
    data$diastolic_bp <- data$diastolic_bp_baseline
  }
  if (!("systolic_bp" %in% names(data)) && ("systolic_bp_baseline" %in% names(data))) {
    data$systolic_bp <- data$systolic_bp_baseline
  }
  # 候选非成像变量：吸烟、饮酒、Townsend、BMI、糖尿病、收缩压，若存在则扩展饮食/运动相关变量
  if (!is.null(non_imaging_covariates) && length(non_imaging_covariates) > 0) {
    non_img_candidates <- non_imaging_covariates[non_imaging_covariates %in% names(data)]
  } else {
    non_img_candidates <- c("smoking_factor", "alcohol_factor", "townsend_index",
                            "bmi", "diabetes_factor", "systolic_bp", "diastolic_bp")
    extra_exercise <- grep("exercise|physical|activity", names(data), ignore.case = TRUE, value = TRUE)
    non_img_candidates <- unique(c(non_img_candidates, extra_exercise))
  }
  
  # 先对所有IDP对计算基础模型未校正p值，并进行FDR校正；后续仅对FDR显著者做衰减分析
  base_pvals <- rep(NA_real_, length(pairs))
  for(i in seq_along(pairs)) {
    pair <- pairs[[i]]
    if(!all(c(pair$idp1_var, pair$idp2_var) %in% names(data))) next
    idp1_data <- data[[pair$idp1_var]]
    idp2_data <- data[[pair$idp2_var]]
    if(!is.numeric(idp1_data) || !is.numeric(idp2_data)) next
    complete_cases <- !is.na(idp1_data) & !is.na(idp2_data)
    if(sum(complete_cases) < 3) next
    analysis_data <- data[complete_cases, ]
    analysis_data$IDP1 <- analysis_data[[pair$idp1_var]]
    analysis_data$IDP2 <- analysis_data[[pair$idp2_var]]
    analysis_data$Age2 <- dplyr::coalesce(analysis_data$Age2,
                                          analysis_data$age_at_instance3,
                                          analysis_data$age_at_instance2 + analysis_data$followup_years,
                                          analysis_data$baseline_age + analysis_data$followup_years,
                                          analysis_data$age_at_instance2,
                                          analysis_data$baseline_age,
                                          analysis_data$age_at_recruitment)
    analysis_data$Age2 <- ifelse(is.na(analysis_data$Age2) | !is.finite(analysis_data$Age2), stats::median(analysis_data$Age2, na.rm = TRUE), analysis_data$Age2)
    if (!("case_control_binary" %in% names(analysis_data))) analysis_data$case_control_binary <- NA_real_
    if ("comparison_group" %in% names(analysis_data)) {
      if (exists("cmp_labels") && model_name %in% names(cmp_labels)) {
        pos_lab <- cmp_labels[[model_name]][1]
        analysis_data$case_control_binary <- dplyr::coalesce(analysis_data$case_control_binary, as.numeric(analysis_data$comparison_group == pos_lab))
      } else {
        labs <- unique(analysis_data$comparison_group)
        analysis_data$case_control_binary <- dplyr::coalesce(analysis_data$case_control_binary, as.numeric(analysis_data$comparison_group == labs[1]))
      }
    } else if ("group" %in% names(analysis_data)) {
      analysis_data$case_control_binary <- dplyr::coalesce(analysis_data$case_control_binary, as.numeric(as.character(analysis_data$group)) == 1)
    }
    analysis_data <- analysis_data[!is.na(analysis_data$case_control_binary) & !is.na(analysis_data$Age2), ]
    case_control_demeaned <- analysis_data$case_control_binary - mean(analysis_data$case_control_binary)
    aging_multiplier <- 10^(analysis_data$Age2 * 0.0524 - 3.27)
    analysis_data$Case_vs_Control <- case_control_demeaned * aging_multiplier
    confounds_basic <- c("age_diff_linear", "age_sq_diff_linear")
    if (!isTRUE(drop_demographics)) {
      if("sex_factor" %in% names(analysis_data)) confounds_basic <- c(confounds_basic, "sex_factor")
      if("ethnicity_factor" %in% names(analysis_data)) confounds_basic <- c(confounds_basic, "ethnicity_factor")
    }
    formula_basic <- "IDP2 ~ Case_vs_Control + IDP1"
    available_confounds <- confounds_basic[confounds_basic %in% names(analysis_data)]
    if(length(available_confounds) > 0) {
      formula_basic <- paste(formula_basic, "+", paste(available_confounds, collapse = " + "))
    }
    base_fit <- tryCatch(lm(as.formula(formula_basic), data = analysis_data), error = function(e) NULL)
    if(is.null(base_fit)) next
    base_coef <- coef(summary(base_fit))
    cc_row <- which(rownames(base_coef) == "Case_vs_Control")
    if(length(cc_row) == 0) next
    base_pvals[i] <- base_coef[cc_row, "Pr(>|t|)"]
  }
  # FDR筛选显著IDP
  p_bh <- stats::p.adjust(base_pvals, method = "BH")
  significant_idx <- which(!is.na(p_bh) & p_bh < 0.05)
  if(length(significant_idx) == 0) return(NULL)
  
  for(i in significant_idx) {
    pair <- pairs[[i]]
    idp1_data <- data[[pair$idp1_var]]
    idp2_data <- data[[pair$idp2_var]]
    complete_cases <- !is.na(idp1_data) & !is.na(idp2_data)
    analysis_data <- data[complete_cases, ]
    analysis_data$IDP1 <- analysis_data[[pair$idp1_var]]
    analysis_data$IDP2 <- analysis_data[[pair$idp2_var]]
    analysis_data$Age2 <- dplyr::coalesce(analysis_data$Age2,
                                          analysis_data$age_at_instance3,
                                          analysis_data$age_at_instance2 + analysis_data$followup_years,
                                          analysis_data$baseline_age + analysis_data$followup_years,
                                          analysis_data$age_at_instance2,
                                          analysis_data$baseline_age,
                                          analysis_data$age_at_recruitment)
    analysis_data$Age2 <- ifelse(is.na(analysis_data$Age2) | !is.finite(analysis_data$Age2), stats::median(analysis_data$Age2, na.rm = TRUE), analysis_data$Age2)
    if (!("case_control_binary" %in% names(analysis_data))) analysis_data$case_control_binary <- NA_real_
    if ("comparison_group" %in% names(analysis_data)) {
      if (exists("cmp_labels") && model_name %in% names(cmp_labels)) {
        pos_lab <- cmp_labels[[model_name]][1]
        analysis_data$case_control_binary <- dplyr::coalesce(analysis_data$case_control_binary, as.numeric(analysis_data$comparison_group == pos_lab))
      } else {
        labs <- unique(analysis_data$comparison_group)
        analysis_data$case_control_binary <- dplyr::coalesce(analysis_data$case_control_binary, as.numeric(analysis_data$comparison_group == labs[1]))
      }
    } else if ("group" %in% names(analysis_data)) {
      analysis_data$case_control_binary <- dplyr::coalesce(analysis_data$case_control_binary, as.numeric(as.character(analysis_data$group)) == 1)
    }
    analysis_data <- analysis_data[!is.na(analysis_data$case_control_binary) & !is.na(analysis_data$Age2), ]
    case_control_demeaned <- analysis_data$case_control_binary - mean(analysis_data$case_control_binary)
    aging_multiplier <- 10^(analysis_data$Age2 * 0.0524 - 3.27)
    analysis_data$Case_vs_Control <- case_control_demeaned * aging_multiplier
    confounds_basic <- c("age_diff_linear", "age_sq_diff_linear")
    if (!isTRUE(drop_demographics)) {
      if("sex_factor" %in% names(analysis_data)) confounds_basic <- c(confounds_basic, "sex_factor")
      if("ethnicity_factor" %in% names(analysis_data)) confounds_basic <- c(confounds_basic, "ethnicity_factor")
    }
    formula_basic <- "IDP2 ~ Case_vs_Control + IDP1"
    available_confounds <- confounds_basic[confounds_basic %in% names(analysis_data)]
    if(length(available_confounds) > 0) {
      formula_basic <- paste(formula_basic, "+", paste(available_confounds, collapse = " + "))
    }
    base_fit <- tryCatch(lm(as.formula(formula_basic), data = analysis_data), error = function(e) NULL)
    if(is.null(base_fit)) next
    base_coef <- coef(summary(base_fit))
    cc_row <- which(rownames(base_coef) == "Case_vs_Control")
    if(length(cc_row) == 0) next
    z0 <- base_coef[cc_row, "t value"]
    beta0 <- base_coef[cc_row, "Estimate"]
    se0 <- base_coef[cc_row, "Std. Error"]
    
    for(var in non_img_candidates) {
      if(!(var %in% names(analysis_data))) next
      f_att <- paste(formula_basic, "+", var)
      fit_att <- tryCatch(lm(as.formula(f_att), data = analysis_data), error = function(e) NULL)
      if(is.null(fit_att)) next
      ct <- coef(summary(fit_att))
      rr <- which(rownames(ct) == "Case_vs_Control")
      if(length(rr) == 0) next
      z1 <- ct[rr, "t value"]
      p1 <- ct[rr, "Pr(>|t|)"]
      beta1 <- ct[rr, "Estimate"]
      se1 <- ct[rr, "Std. Error"]
      z_change_pct <- ((z1 - z0) / z0) * 100
      attenuated_25pct <- (abs(z1) < abs(z0)) && (abs(z_change_pct) >= 25)
      
      out_rows[[length(out_rows) + 1]] <- data.frame(
        model_comparison = model_name,
        variable_name = substr(pair$base_name, 1, 100),
        variable_type = pair$variable_type,
        idp1_variable = pair$idp1_var,
        idp2_variable = pair$idp2_var,
        non_imaging_variable = var,
        z_without_non_imaging = round(z0, 4),
        beta_without_non_imaging = round(beta0, 6),
        se_without_non_imaging = round(se0, 6),
        z_with_non_imaging = round(z1, 4),
        p_value_with_non_imaging = p1,
        beta_with_non_imaging_single = round(beta1, 6),
        se_with_non_imaging_single = round(se1, 6),
        z_change_percent = round(z_change_pct, 2),
        attenuated_25pct = attenuated_25pct,
        model_formula_with_non_imaging = f_att,
        stringsAsFactors = FALSE
      )
    }
  }
  if(length(out_rows) == 0) return(NULL)
  do.call(rbind, out_rows)
}

# 置换-FWE计算（扩展模型）：对同一模型内的所有“加入单个非成像变量”的检验进行 max-|t| 置换，返回每行的 FWE p 值
compute_fwe_pvalues_for_attenuation <- function(data_model, attenuation_df, drop_demographics = FALSE, n_perm = 50, seed = 202311) {
  if (is.null(attenuation_df) || nrow(attenuation_df) == 0) return(rep(NA_real_, 0))
  set.seed(seed)
  
  aging_mult_global <- 10^(data_model$Age2 * 0.0524 - 3.27)
  
  build_item <- function(row) {
    idp1 <- attenuation_df$idp1_variable[row]
    idp2 <- attenuation_df$idp2_variable[row]
    form_att <- attenuation_df$model_formula_with_non_imaging[row]
    if (is.na(form_att) || !nzchar(form_att)) return(NULL)
    if (!(idp1 %in% names(data_model)) || !(idp2 %in% names(data_model))) return(NULL)
    idp1_data <- data_model[[idp1]]
    idp2_data <- data_model[[idp2]]
    cc <- !is.na(idp1_data) & !is.na(idp2_data) & !is.na(data_model$case_control_binary) & !is.na(data_model$Age2)
    if (sum(cc) < 3) return(NULL)
    row_idx <- which(cc)
    
    # 构建一次性最小数据子集并计算观测 t
    df_min <- data.frame(
      IDP2 = data_model[row_idx, idp2],
      IDP1 = data_model[row_idx, idp1],
      Case_vs_Control = (data_model$case_control_binary[row_idx] - mean(data_model$case_control_binary[row_idx])) * aging_mult_global[row_idx],
      stringsAsFactors = FALSE
    )
    # 提取公式涉及的协变量列
    vars_in_form <- setdiff(all.vars(as.formula(form_att)), c("IDP2", "Case_vs_Control", "IDP1"))
    if (length(vars_in_form) > 0) {
      avail <- vars_in_form[vars_in_form %in% names(data_model)]
      if (length(avail) > 0) df_min[avail] <- data_model[row_idx, avail, drop = FALSE]
    }
    
    fit0 <- tryCatch(lm(as.formula(form_att), data = df_min), error = function(e) NULL)
    if (is.null(fit0)) return(NULL)
    ct <- coef(summary(fit0))
    rr <- which(rownames(ct) == "Case_vs_Control")
    if (length(rr) == 0) return(NULL)
    t_obs <- as.numeric(ct[rr, "t value"])
    
    list(row_idx = row_idx, idp1 = idp1, idp2 = idp2, formula = form_att, t_obs = t_obs, extra_vars = vars_in_form)
  }
  
  items <- lapply(seq_len(nrow(attenuation_df)), build_item)
  keep <- !sapply(items, is.null)
  if (!any(keep)) return(rep(NA_real_, nrow(attenuation_df)))
  items_idx <- which(keep)
  
  # 逐次置换，记录每次的最大 |t|（内存友好）
  cat(sprintf("[FWE] compute_fwe_pvalues_for_attenuation: n_perm=%d, items=%d\n", n_perm, length(items_idx)))
  pb_perm <- txtProgressBar(min = 1, max = n_perm, style = 3)
  max_t_perm <- numeric(n_perm)
  for (b in seq_len(n_perm)) {
    max_abs_t <- 0
    for (j in seq_along(items_idx)) {
      it <- items[[items_idx[j]]]
      ridx <- it$row_idx
      perm_cc <- sample(data_model$case_control_binary[ridx])
      df_min <- data.frame(
        IDP2 = data_model[ridx, it$idp2],
        IDP1 = data_model[ridx, it$idp1],
        Case_vs_Control = (perm_cc - mean(perm_cc)) * aging_mult_global[ridx],
        stringsAsFactors = FALSE
      )
      if (length(it$extra_vars) > 0) {
        avail <- it$extra_vars[it$extra_vars %in% names(data_model)]
        if (length(avail) > 0) df_min[avail] <- data_model[ridx, avail, drop = FALSE]
      }
      fitp <- tryCatch(lm(as.formula(it$formula), data = df_min), error = function(e) NULL)
      if (is.null(fitp)) next
      ctp <- coef(summary(fitp))
      rrp <- which(rownames(ctp) == "Case_vs_Control")
      if (length(rrp) > 0) {
        tv <- abs(as.numeric(ctp[rrp, "t value"]))
        if (is.finite(tv) && tv > max_abs_t) max_abs_t <- tv
      }
    }
    max_t_perm[b] <- max_abs_t
    rm(max_abs_t)
    setTxtProgressBar(pb_perm, b)
  }
  close(pb_perm)
  
  p_fwe <- rep(NA_real_, nrow(attenuation_df))
  for (j in seq_along(items_idx)) {
    idx <- items_idx[j]
    t_obs <- abs(items[[idx]]$t_obs)
    p_fwe[idx] <- (1 + sum(max_t_perm >= t_obs)) / (n_perm + 1)
  }
  p_fwe
}

# =================== 第五步：执行四种模型的Z统计量分析 ===================
cat("\n第五步：执行四种模型的Z统计量分析...\n")
# 计数口径开关：限定“稳定/独立”和“中介/受影响”的计数范围
# 可选取值："none"（不限）、"raw"（原始p<0.05）、"fdr"（FDR<0.05）
# 默认改为 "raw"，按用户要求以显著差异为前置门槛
if (!exists("attenuation_count_gate_mode")) {
  attenuation_count_gate_mode <- "raw"
}

# 为确保每个模型严格使用其对应的匹配队列，这里从每个比较组的RDS文件重建分析数据集
## 移除try包装，改为直接执行
  if (!exists("comparisons")) {
    comparisons <- c("ischemic_vs_control", "mi_vs_control", "chronic_vs_control")
  }
  
  # 定义各比较的组标签（用于可视化一致性）
  cmp_labels <- list(
    ischemic_vs_control = c("Ischemic_Heart_Disease", "Control"),
    mi_vs_control       = c("Myocardial_Infarction", "Control"),
    chronic_vs_control  = c("Chronic_Ischemic", "Control")
  )
  
  # ========= 修改：仅读取主队列RDS，并从主队列派生三个子比较 =========
  analysis_data_list <- list()
  main_rds <- resolve_existing_file(c(
    matched_cohort_rds_path("ischemic_vs_control"),
    "Final_Matched_Cohort_ischemic_vs_control.rds"
  ))
  if (is.na(main_rds) || !nzchar(main_rds)) {
    stop("未找到主匹配队列文件。请先生成 ischemic_vs_control 的匹配结果。")
  }
  dat_main <- readRDS(main_rds)
  # 计算年龄变量（健壮处理，避免0长度替换错误）
  dat_main <- compute_age_instances(dat_main)
  dat_main <- canonicalize_instance_columns(dat_main)
  dat_main <- dedupe_instance_columns(dat_main)
  dat_main <- enforce_two_timepoint_presence(dat_main)
  # 主队列的二分类与标签（case: Ischemic_Heart_Disease, control: Control）
  grp_num_main <- suppressWarnings(as.numeric(as.character(dat_main$group)))
  dat_main$case_control_binary <- as.numeric(grp_num_main == 1)
  dat_main$comparison_group <- ifelse(grp_num_main == 1, cmp_labels[["ischemic_vs_control"]][1], cmp_labels[["ischemic_vs_control"]][2])
  # 派生分析所需变量
  dat_main$Age2 <- dplyr::coalesce(dat_main$age_at_instance3,
                                   dat_main$age_at_instance2 + dat_main$followup_years,
                                   dat_main$baseline_age + dat_main$followup_years,
                                   dat_main$age_at_instance2,
                                   dat_main$baseline_age,
                                   dat_main$age_at_recruitment)
  dat_main$Age2 <- ifelse(is.na(dat_main$Age2) | !is.finite(dat_main$Age2), stats::median(dat_main$Age2, na.rm = TRUE), dat_main$Age2)
  dat_main$age_diff <- dat_main$age_at_instance3 - dat_main$age_at_instance2
  dat_main$age_diff_sq <- (dat_main$age_at_instance3 - dat_main$age_at_instance2)^2
  dat_main$age_diff_linear <- dat_main$age_at_instance3 - dat_main$age_at_instance2
  dat_main$age_sq_diff_linear <- dat_main$age_at_instance3^2 - dat_main$age_at_instance2^2
  # 放入主队列
  analysis_data_list[["ischemic_vs_control"]] <- dat_main

  write_instance_name_diagnostics <- function(df) {
    nm <- names(df)
    inst3 <- nm[grepl("Instance(\\.|_|\\s)?3", nm)]
    if (length(inst3) == 0) return(invisible(NULL))
    n <- nrow(df)
    metrics <- lapply(inst3, function(v) {
      x <- suppressWarnings(as.numeric(df[[v]]))
      miss <- sum(is.na(x))
      data.frame(
        timepoint = "Age2",
        variable_name = v,
        n_total = n,
        n_missing = miss,
        missing_rate = round(miss / n, 6),
        stringsAsFactors = FALSE
      )
    })
    out <- do.call(rbind, metrics)
    utils::write.csv(out, table_out_path("Diagnostics", "Instance3_Name_Diagnostics.csv"), row.names = FALSE)
    cat(sprintf("[诊断] 已写出 Instance3_Name_Diagnostics.csv：%d 列\n", nrow(out)))
  }

  write_instance_name_diagnostics(dat_main)

  write_canonical_mapping <- function(df) {
    nm <- names(df)
    if (length(nm) == 0) return(invisible(NULL))
    norm <- function(s) {
      s <- sub("(\\.x|\\.y)$", "", s)
      s <- gsub("Instance[_ ]([23])", "Instance.\\1", s, perl = TRUE)
      s
    }
    base <- sapply(nm, norm)
    keep <- grepl("Instance\\.[23]", base)
    if (!any(keep)) return(invisible(NULL))
    tp <- ifelse(grepl("Instance\\.2", base), "Age1", ifelse(grepl("Instance\\.3", base), "Age2", NA_character_))
    out <- data.frame(
      base_name = base[keep],
      variant_name = nm[keep],
      timepoint = tp[keep],
      stringsAsFactors = FALSE
    )
    out <- out[order(out$base_name, out$variant_name), ]
    utils::write.csv(out, table_out_path("Diagnostics", "Instance_Canonical_Column_Mapping.csv"), row.names = FALSE)
    cat(sprintf("[诊断] 已写出 Instance_Canonical_Column_Mapping.csv：%d 行\n", nrow(out)))
  }

  write_canonical_mapping(dat_main)

  write_instance2_name_diagnostics <- function(df) {
    nm <- names(df)
    inst2 <- nm[grepl("Instance(\\.|_|\\s)?2", nm)]
    if (length(inst2) == 0) return(invisible(NULL))
    n <- nrow(df)
    metrics <- lapply(inst2, function(v) {
      x <- suppressWarnings(as.numeric(df[[v]]))
      miss <- sum(is.na(x))
      data.frame(
        timepoint = "Age1",
        variable_name = v,
        n_total = n,
        n_missing = miss,
        missing_rate = round(miss / n, 6),
        stringsAsFactors = FALSE
      )
    })
    out <- do.call(rbind, metrics)
    utils::write.csv(out, table_out_path("Diagnostics", "Instance2_Name_Diagnostics.csv"), row.names = FALSE)
    cat(sprintf("[诊断] 已写出 Instance2_Name_Diagnostics.csv：%d 列\n", nrow(out)))
  }

  write_instance2_name_diagnostics(dat_main)

  write_pair_root_presence <- function(df) {
    nm <- names(df)
    norm <- function(s) {
      s <- sub("(\\.x|\\.y)$", "", s)
      s <- gsub("Instance[_ ]([23])", "Instance.\\1", s, perl = TRUE)
      s
    }
    bases <- unique(sapply(nm, norm))
    inst2_bases <- bases[grepl("Instance\\.2$", bases)]
    inst3_of <- function(b) sub("Instance\\.2$", "Instance.3", b)
    rows <- lapply(inst2_bases, function(b2) {
      b3 <- inst3_of(b2)
      has2 <- b2 %in% nm
      has3 <- b3 %in% nm
      n <- nrow(df)
      x2 <- if (has2) suppressWarnings(as.numeric(df[[b2]])) else rep(NA_real_, n)
      x3 <- if (has3) suppressWarnings(as.numeric(df[[b3]])) else rep(NA_real_, n)
      n2 <- sum(!is.na(x2))
      n3 <- sum(!is.na(x3))
      n_pair <- sum(!is.na(x2) & !is.na(x3))
      data.frame(
        base_name_age1 = b2,
        base_name_age2 = b3,
        has_age1_col = has2,
        has_age2_col = has3,
        n_total = n,
        n_age1_nonmissing = n2,
        n_age2_nonmissing = n3,
        n_pair_complete = n_pair,
        coverage_pair = round(n_pair / n, 6),
        stringsAsFactors = FALSE
      )
    })
    out <- do.call(rbind, rows)
    utils::write.csv(out, table_out_path("Diagnostics", "Instance_Pair_Root_Presence.csv"), row.names = FALSE)
    cat(sprintf("[诊断] 已写出 Instance_Pair_Root_Presence.csv：%d 行\n", nrow(out)))
  }

  write_pair_root_presence(dat_main)

  write_idp_missingness_pair_report <- function(df) {
    nm <- names(df)
    n <- nrow(df)
    case_ok <- "case_control_binary" %in% names(df)
    base2 <- nm[grepl("Instance(\\.|_|\\s)?2$", nm) & !grepl("(\\.x|\\.y)$", nm)]
    rows <- lapply(base2, function(b2) {
      b3 <- sub("Instance(\\.|_|\\s)?2$", "Instance.3", b2)
      has2 <- b2 %in% nm
      has3 <- b3 %in% nm
      x2 <- if (has2) suppressWarnings(as.numeric(df[[b2]])) else rep(NA_real_, n)
      x3 <- if (has3) suppressWarnings(as.numeric(df[[b3]])) else rep(NA_real_, n)
      miss1 <- sum(is.na(x2))
      miss2 <- sum(is.na(x3))
      pair_complete <- sum(!is.na(x2) & !is.na(x3))
      rate1 <- if (n > 0) round(miss1 / n, 6) else NA_real_
      rate2 <- if (n > 0) round(miss2 / n, 6) else NA_real_
      pair_rate <- if (n > 0) round(pair_complete / n, 6) else NA_real_
      n_pair_cases <- if (case_ok) sum(!is.na(x2) & !is.na(x3) & df$case_control_binary == 1, na.rm = TRUE) else NA_integer_
      n_pair_controls <- if (case_ok) sum(!is.na(x2) & !is.na(x3) & df$case_control_binary == 0, na.rm = TRUE) else NA_integer_
      flag_pair_gap <- (!is.na(rate1) && !is.na(rate2) && !is.na(pair_rate) && rate1 <= 0.2 && rate2 <= 0.2 && pair_rate <= 0.3)
      data.frame(
        base_name_age1 = b2,
        base_name_age2 = b3,
        n_total = n,
        n_age1_missing = miss1,
        age1_missing_rate = rate1,
        n_age2_missing = miss2,
        age2_missing_rate = rate2,
        n_pair_complete = pair_complete,
        pair_complete_rate = pair_rate,
        n_pair_complete_cases = n_pair_cases,
        n_pair_complete_controls = n_pair_controls,
        flag_pair_gap = flag_pair_gap,
        stringsAsFactors = FALSE
      )
    })
    out <- do.call(rbind, rows)
    utils::write.csv(out, table_out_path("Diagnostics", "IDP_Missingness_By_Timepoint_and_Pair.csv"), row.names = FALSE)
    cat(sprintf("[诊断] 已写出 IDP_Missingness_By_Timepoint_and_Pair.csv：%d 行\n", nrow(out)))
  }

  write_idp_missingness_pair_report(dat_main)

  write_descriptive_pair_counts_ischemic <- function(dat, variable_pairs) {
    base_df <- dat
    main_rds <- resolve_existing_file(c(
      matched_cohort_rds_path("ischemic_vs_control"),
      "Final_Matched_Cohort_ischemic_vs_control.rds"
    ))
    if (!is.na(main_rds) && nzchar(main_rds)) {
      base_df <- readRDS(main_rds)
      base_df <- compute_age_instances(base_df)
      base_df <- canonicalize_instance_columns(base_df)
      base_df <- dedupe_instance_columns(base_df)
      if ("group" %in% names(base_df)) {
        grp_num_main <- suppressWarnings(as.numeric(as.character(base_df$group)))
        base_df$case_control_binary <- as.numeric(grp_num_main == 1)
        base_df$comparison_group <- ifelse(grp_num_main == 1, cmp_labels[["ischemic_vs_control"]][1], cmp_labels[["ischemic_vs_control"]][2])
      } else if ("Group" %in% names(base_df)) {
        base_df$case_control_binary <- as.numeric(as.character(base_df$Group) %in% c("Ischemic Heart Disease"))
        base_df$comparison_group <- ifelse(base_df$case_control_binary == 1, cmp_labels[["ischemic_vs_control"]][1], cmp_labels[["ischemic_vs_control"]][2])
      }
    }
    pairs_eval <- Filter(function(p) p$variable_type == "Brain_Imaging", variable_pairs)
    pairs_eval <- Filter(function(p) all(c(p$idp1_var, p$idp2_var) %in% names(base_df)), pairs_eval)
    if (length(pairs_eval) == 0) return(invisible(NULL))
    n <- nrow(base_df)
    case_ok <- "case_control_binary" %in% names(base_df)
    rows <- lapply(pairs_eval, function(p) {
      x1 <- suppressWarnings(as.numeric(base_df[[p$idp1_var]]))
      x2 <- suppressWarnings(as.numeric(base_df[[p$idp2_var]]))
      both <- !is.na(x1) & !is.na(x2)
      inst2_only <- !is.na(x1) & is.na(x2)
      inst3_only <- is.na(x1) & !is.na(x2)
      n_both <- sum(both)
      rate_both <- if (n > 0) round(n_both / n, 6) else NA_real_
      n_both_cases <- if (case_ok) sum(both & base_df$case_control_binary == 1, na.rm = TRUE) else NA_integer_
      n_both_controls <- if (case_ok) sum(both & base_df$case_control_binary == 0, na.rm = TRUE) else NA_integer_
      n_inst2_only_cases <- if (case_ok) sum(inst2_only & base_df$case_control_binary == 1, na.rm = TRUE) else NA_integer_
      n_inst2_only_controls <- if (case_ok) sum(inst2_only & base_df$case_control_binary == 0, na.rm = TRUE) else NA_integer_
      n_inst3_only_cases <- if (case_ok) sum(inst3_only & base_df$case_control_binary == 1, na.rm = TRUE) else NA_integer_
      n_inst3_only_controls <- if (case_ok) sum(inst3_only & base_df$case_control_binary == 0, na.rm = TRUE) else NA_integer_
      data.frame(
        base_name = p$base_name,
        idp1_variable = p$idp1_var,
        idp2_variable = p$idp2_var,
        n_total = n,
        n_both = n_both,
        rate_both = rate_both,
        n_both_cases = n_both_cases,
        n_both_controls = n_both_controls,
        n_inst2_only_cases = n_inst2_only_cases,
        n_inst2_only_controls = n_inst2_only_controls,
        n_inst3_only_cases = n_inst3_only_cases,
        n_inst3_only_controls = n_inst3_only_controls,
        stringsAsFactors = FALSE
      )
    })
    out <- do.call(rbind, rows)
    utils::write.csv(out, table_out_path("Diagnostics", "Descriptive_Pair_Counts_ischemic_vs_control.csv"), row.names = FALSE)
    cat(sprintf("[描述] 已写出 Descriptive_Pair_Counts_ischemic_vs_control.csv：%d 行\n", nrow(out)))
    out
  }

  write_descriptive_pair_counts_ischemic(dat_main, variable_pairs)

  # 从完整合并数据的 disease_subtype 中取 eid，派生子队列（确保变量一致）
  if (!exists("ischemic_cohort_complete")) {
    warning("未找到 ischemic_cohort_complete，无法根据 disease_subtype 派生子队列。仅保留主队列。")
  } else if (!("eid" %in% names(dat_main))) {
    warning("主匹配队列缺少 eid 列，无法派生子队列。仅保留主队列。")
  } else {
    mi_eids_all <- tryCatch({
      subset(ischemic_cohort_complete, disease_subtype %in% c("MI Only", "Both MI & Chronic"))$eid
    }, error = function(e) integer(0))
    mi_only_eids <- tryCatch({
      subset(ischemic_cohort_complete, disease_subtype %in% c("MI Only"))$eid
    }, error = function(e) integer(0))
    chronic_eids <- tryCatch({
      subset(ischemic_cohort_complete, disease_subtype %in% c("Chronic Only"))$eid
    }, error = function(e) integer(0))
    control_eids <- tryCatch({
      subset(ischemic_cohort_complete, disease_subtype %in% c("Control"))$eid
    }, error = function(e) integer(0))

    # mi_vs_control：从主队列取 MI + Control
    mi_mask <- dat_main$eid %in% mi_eids_all
    ctrl_mask <- dat_main$eid %in% control_eids
    mi_ctrl_mask <- mi_mask | ctrl_mask
    if (any(mi_ctrl_mask)) {
      dat_mi_ctrl <- dat_main[mi_ctrl_mask, , drop = FALSE]
      dat_mi_ctrl$comparison_group <- ifelse(dat_mi_ctrl$eid %in% mi_eids_all, cmp_labels[["mi_vs_control"]][1], cmp_labels[["mi_vs_control"]][2])
      dat_mi_ctrl$case_control_binary <- as.numeric(dat_mi_ctrl$eid %in% mi_eids_all)
      analysis_data_list[["mi_vs_control"]] <- dat_mi_ctrl
    } else {
      warning("根据主队列未能派生 mi_vs_control（可能是 eid 集合为空）。")
    }

    # chronic_vs_control：从主队列取 Chronic + Control
    ch_mask <- dat_main$eid %in% chronic_eids
    ch_ctrl_mask <- ch_mask | ctrl_mask
    if (any(ch_ctrl_mask)) {
      dat_ch_ctrl <- dat_main[ch_ctrl_mask, , drop = FALSE]
      dat_ch_ctrl$comparison_group <- ifelse(dat_ch_ctrl$eid %in% chronic_eids, cmp_labels[["chronic_vs_control"]][1], cmp_labels[["chronic_vs_control"]][2])
      dat_ch_ctrl$case_control_binary <- as.numeric(dat_ch_ctrl$eid %in% chronic_eids)
      analysis_data_list[["chronic_vs_control"]] <- dat_ch_ctrl
    } else {
      warning("根据主队列未能派生 chronic_vs_control（可能是 eid 集合为空）。")
    }

  }

  # 样本量统计（从主队列派生）
  cat("三种模型样本量统计（从主队列派生）：\n")
  for (model_name in names(analysis_data_list)) {
    data <- analysis_data_list[[model_name]]
    group_counts <- table(data$comparison_group)
    cat(sprintf("%-20s - 总样本: %5d - %s\n",
                model_name, nrow(data),
                paste(names(group_counts), "=", group_counts, collapse = ", ")))
  }
  diagnose_idp_pair_completeness <- function(data, variable_pairs, model_label) {
    pairs_eval <- Filter(function(p) p$variable_type == "Brain_Imaging", variable_pairs)
    if (length(pairs_eval) == 0) return(NULL)
    n_total <- nrow(data)
    out_rows <- lapply(pairs_eval, function(p) {
      if (!all(c(p$idp1_var, p$idp2_var) %in% names(data))) return(NULL)
      x1 <- suppressWarnings(as.numeric(data[[p$idp1_var]]))
      x2 <- suppressWarnings(as.numeric(data[[p$idp2_var]]))
      pair_ok <- !is.na(x1) & !is.na(x2)
      group_ok <- if ("case_control_binary" %in% names(data)) !is.na(data$case_control_binary) else TRUE
      age_ok <- if ("Age2" %in% names(data)) !is.na(data$Age2) else TRUE
      model_ok <- pair_ok & group_ok & age_ok
      n_cases_pair <- if ("case_control_binary" %in% names(data)) sum(pair_ok & data$case_control_binary == 1, na.rm = TRUE) else NA_integer_
      n_ctrls_pair <- if ("case_control_binary" %in% names(data)) sum(pair_ok & data$case_control_binary == 0, na.rm = TRUE) else NA_integer_
      n_cases_model <- if ("case_control_binary" %in% names(data)) sum(model_ok & data$case_control_binary == 1, na.rm = TRUE) else NA_integer_
      n_ctrls_model <- if ("case_control_binary" %in% names(data)) sum(model_ok & data$case_control_binary == 0, na.rm = TRUE) else NA_integer_
      data.frame(
        model = model_label,
        variable_name = p$base_name,
        idp1_variable = p$idp1_var,
        idp2_variable = p$idp2_var,
        n_total = n_total,
        n_pair_complete = sum(pair_ok),
        n_model_complete = sum(model_ok),
        n_pair_missing_age2 = sum(pair_ok & !age_ok),
        n_pair_missing_group = sum(pair_ok & !group_ok),
        n_cases_pair_complete = n_cases_pair,
        n_controls_pair_complete = n_ctrls_pair,
        n_cases_model_complete = n_cases_model,
        n_controls_model_complete = n_ctrls_model,
        stringsAsFactors = FALSE
      )
    })
    out_df <- do.call(rbind, out_rows[!sapply(out_rows, is.null)])
    if (!is.null(out_df) && nrow(out_df) > 0) {
      utils::write.csv(out_df, table_out_path("Diagnostics", paste0("Diagnostics_IDP_Completeness_", model_label, ".csv")), row.names = FALSE)
    }
    out_df
  }
  for (ml in names(analysis_data_list)) {
    try(diagnose_idp_pair_completeness(analysis_data_list[[ml]], variable_pairs, ml), silent = TRUE)
  }
  # 此处不再进行并集变量补齐：四个模型均源自同一主队列的列集合
  # ============ 新增：完备性与样本量诊断输出（横断-认知） ============
  write_cross_sectional_diagnostics <- function(analysis_data_list, variable_pairs,
                                                basic_threshold = 3, extended_threshold = 25) {
    if (length(variable_pairs) == 0) return(invisible(NULL))
    idp_all <- unique(unlist(lapply(variable_pairs, function(p) c(p$idp1_var, p$idp2_var))))
    for (model_name in names(analysis_data_list)) {
      dat <- analysis_data_list[[model_name]]
      idp_vars <- idp_all[idp_all %in% names(dat)]
      non_missing_rates <- vapply(idp_vars, function(v) mean(!is.na(dat[[v]])), numeric(1))
      idp_rank_df <- data.frame(IDP = idp_vars, NonMissingRate = non_missing_rates, stringsAsFactors = FALSE)
      idp_rank_df <- idp_rank_df[order(idp_rank_df$NonMissingRate, decreasing = TRUE), ]
      base_non_img <- c("townsend_index_baseline","bmi_baseline","systolic_bp_baseline","diastolic_bp_baseline",
                        "diabetes_factor_baseline","smoking_factor_baseline","alcohol_factor_baseline",
                        "waist_circumference_baseline")
      non_img_covars <- base_non_img
      non_img_covars <- intersect(non_img_covars, names(dat))
      dat$Age2 <- dplyr::coalesce(dat$Age2,
                                  dat$age_at_instance3,
                                  dat$age_at_instance2 + dat$followup_years,
                                  dat$baseline_age + dat$followup_years,
                                  dat$age_at_instance2,
                                  dat$baseline_age,
                                  dat$age_at_recruitment)
      dat$Age2 <- ifelse(is.na(dat$Age2) | !is.finite(dat$Age2), stats::median(dat$Age2, na.rm = TRUE), dat$Age2)
      if (!("case_control_binary" %in% names(dat))) dat$case_control_binary <- NA_real_
      if ("comparison_group" %in% names(dat)) {
        if (exists("cmp_labels") && model_name %in% names(cmp_labels)) {
          pos_lab <- cmp_labels[[model_name]][1]
          dat$case_control_binary <- dplyr::coalesce(dat$case_control_binary, as.numeric(dat$comparison_group == pos_lab))
        } else {
          labs <- unique(dat$comparison_group)
          dat$case_control_binary <- dplyr::coalesce(dat$case_control_binary, as.numeric(dat$comparison_group == labs[1]))
        }
      } else if ("group" %in% names(dat)) {
        dat$case_control_binary <- dplyr::coalesce(dat$case_control_binary, as.numeric(as.character(dat$group)) == 1)
      }
      pair_stats <- lapply(variable_pairs, function(p) {
        idp1 <- p$idp1_var; idp2 <- p$idp2_var
        if (!all(c(idp1, idp2) %in% names(dat))) return(NULL)
        basic_cc <- !is.na(dat[[idp1]]) & !is.na(dat[[idp2]]) & !is.na(dat$case_control_binary) & !is.na(dat$Age2)
        n_basic <- sum(basic_cc)
        extended_cc <- basic_cc
        if (length(non_img_covars) > 0) {
          extended_cc <- extended_cc & apply(!is.na(dat[non_img_covars]), 1, all)
        }
        n_extended <- sum(extended_cc)
        data.frame(idp1_var = idp1, idp2_var = idp2, n_basic_complete = n_basic, n_extended_complete = n_extended,
                   pass_basic = n_basic >= basic_threshold, pass_extended = n_extended >= extended_threshold,
                   stringsAsFactors = FALSE)
      })
      pair_stats_valid <- pair_stats[!sapply(pair_stats, is.null)]
      if (length(pair_stats_valid) == 0) next
      pair_stats_df <- do.call(rbind, pair_stats_valid)
      summary_df <- data.frame(Model = model_name, SampleSize = nrow(dat), PairsCount = nrow(pair_stats_df),
                               Proportion_Passing_Basic = mean(pair_stats_df$pass_basic),
                               Proportion_Passing_Extended = mean(pair_stats_df$pass_extended),
                               stringsAsFactors = FALSE)
      excluded_basic <- subset(pair_stats_df, !pass_basic)
      if (nrow(excluded_basic) > 0) excluded_basic$Reason <- sprintf("基础样本不足(<%d)", basic_threshold)
      excluded_extended <- subset(pair_stats_df, pass_basic & !pass_extended)
      if (nrow(excluded_extended) > 0) excluded_extended$Reason <- sprintf("扩展样本不足(<%d)", extended_threshold)
      utils::write.csv(idp_rank_df, table_out_path("Diagnostics", paste0("Diagnostics_", model_name, "_IDP_NonMissing_Ranking.csv")), row.names = FALSE)
      utils::write.csv(pair_stats_df, table_out_path("Diagnostics", paste0("Diagnostics_", model_name, "_Pair_Completeness_Counts.csv")), row.names = FALSE)
      utils::write.csv(summary_df, table_out_path("Diagnostics", paste0("Diagnostics_", model_name, "_Summary.csv")), row.names = FALSE)
      if (nrow(excluded_basic) > 0) utils::write.csv(excluded_basic, table_out_path("Diagnostics", paste0("Diagnostics_", model_name, "_Excluded_Basic.csv")), row.names = FALSE)
      if (nrow(excluded_extended) > 0) utils::write.csv(excluded_extended, table_out_path("Diagnostics", paste0("Diagnostics_", model_name, "_Excluded_Extended.csv")), row.names = FALSE)
      cat(sprintf("[诊断] 已输出横断诊断：%s\n", model_name))
    }
    invisible(TRUE)
  }
  write_cross_sectional_diagnostics(analysis_data_list, variable_pairs)

  # ============ 新增：IDP扫描-重扫描可重复性评估与筛选 ============
  compute_idp_reproducibility <- function(data, variable_pairs, group_var = "comparison_group",
                                          r_threshold = 0.5, min_complete = 20,
                                          restrict_type = c("Brain_Imaging","Cognitive")) {
    if (length(variable_pairs) == 0) return(NULL)
    # 限定评估的IDP类型（默认优先脑影像）
    restrict_type <- match.arg(restrict_type, several.ok = TRUE)
    pairs_eval <- Filter(function(p) p$variable_type %in% restrict_type, variable_pairs)
    if (length(pairs_eval) == 0) return(NULL)
    # 组顺序：优先使用约定标签，保证 r_cases/r_controls 语义一致
    expected_groups <- c("Ischemic_Heart_Disease", "Control")
    present_groups <- unique(data[[group_var]])
    groups <- intersect(expected_groups, present_groups)
    if (length(groups) < length(present_groups)) {
      # 有意外组名时仍按出现顺序追加
      groups <- c(groups, setdiff(present_groups, groups))
    }
    to_numeric_safe <- function(x) {
      if (is.numeric(x)) return(x)
      if (is.factor(x)) return(suppressWarnings(as.numeric(as.character(x))))
      if (is.logical(x) || is.integer(x)) return(as.numeric(x))
      suppressWarnings(as.numeric(x))
    }
    out_rows <- lapply(pairs_eval, function(p) {
      idp1 <- p$idp1_var; idp2 <- p$idp2_var
      if (!all(c(idp1, idp2) %in% names(data))) return(NULL)
      x1 <- to_numeric_safe(data[[idp1]])
      x2 <- to_numeric_safe(data[[idp2]])
      # 按组计算IDP1~IDP2的Pearson相关
      r_values <- vapply(groups, function(g) {
        idx <- !is.na(x1) & !is.na(x2) & data[[group_var]] == g
        if (sum(idx) >= 3) {
          suppressWarnings(stats::cor(x1[idx], x2[idx], use = "complete.obs"))
        } else {
          NA_real_
        }
      }, numeric(1))
      r_mean <- suppressWarnings(mean(r_values, na.rm = TRUE))
      # 完整样本量（两次都不缺失）
      n_complete <- sum(!is.na(x1) & !is.na(x2))
      data.frame(
        variable_type = p$variable_type,
        base_name = p$base_name,
        idp1_var = idp1,
        idp2_var = idp2,
        r_cases = ifelse(length(r_values) >= 1, r_values[1], NA_real_),
        r_controls = ifelse(length(r_values) >= 2, r_values[2], NA_real_),
        r_mean = r_mean,
        n_complete = n_complete,
        pass_r = is.finite(r_mean) && r_mean >= r_threshold,
        pass_n = n_complete >= min_complete,
        pass_both = isTRUE(is.finite(r_mean) && r_mean >= r_threshold) && n_complete >= min_complete,
        stringsAsFactors = FALSE
      )
    })
    out_df <- do.call(rbind, out_rows[!sapply(out_rows, is.null)])
    out_df
  }

  # 在主比较（ischemic_vs_control）上执行可重复性评估并输出
  reproducibility_results <- NULL
  if ("ischemic_vs_control" %in% names(analysis_data_list)) {
    dat_rep <- analysis_data_list[["ischemic_vs_control"]]
    use_full <- isTRUE(as.logical(getOption("viz.use.full.cohort", Sys.getenv("VIZ_USE_FULL_COHORT", "FALSE"))))
    if (use_full && exists("ischemic_cohort_complete")) {
      dat_full <- ischemic_cohort_complete
      dat_full <- compute_age_instances(dat_full)
      dat_full <- canonicalize_instance_columns(dat_full)
      dat_full <- dedupe_instance_columns(dat_full)
      dat_full <- enforce_two_timepoint_presence(dat_full)
      if ("group" %in% names(dat_full)) {
        grp_num_full <- suppressWarnings(as.numeric(as.character(dat_full$group)))
        dat_full$case_control_binary <- as.numeric(grp_num_full == 1)
        dat_full$comparison_group <- ifelse(grp_num_full == 1, cmp_labels[["ischemic_vs_control"]][1], cmp_labels[["ischemic_vs_control"]][2])
      }
      dat_full$Age2 <- dplyr::coalesce(dat_full$age_at_instance3,
                                       dat_full$age_at_instance2 + dat_full$followup_years,
                                       dat_full$baseline_age + dat_full$followup_years,
                                       dat_full$age_at_instance2,
                                       dat_full$baseline_age,
                                       dat_full$age_at_recruitment)
      dat_full$Age2 <- ifelse(is.na(dat_full$Age2) | !is.finite(dat_full$Age2), stats::median(dat_full$Age2, na.rm = TRUE), dat_full$Age2)
      dat_rep <- dat_full
    }
    reproducibility_results <- compute_idp_reproducibility(
      data = dat_rep,
      variable_pairs = variable_pairs,
      group_var = "comparison_group",
      r_threshold = 0.5,
      min_complete = 10,
      restrict_type = c("Brain_Imaging")
    )
    if (!is.null(reproducibility_results) && nrow(reproducibility_results) > 0) {
      utils::write.csv(reproducibility_results, table_out_path("Diagnostics", "Diagnostics_ischemic_vs_control_IDP_Reproducibility.csv"), row.names = FALSE)
      # 输出摘要：初始、r筛选后、样本量筛选后数量
      n_initial <- sum(reproducibility_results$variable_type == "Brain_Imaging")
      n_after_r <- sum(reproducibility_results$variable_type == "Brain_Imaging" & reproducibility_results$pass_r)
      n_after_both <- sum(reproducibility_results$variable_type == "Brain_Imaging" & reproducibility_results$pass_both)
      corr_groups <- suppressWarnings(stats::cor(reproducibility_results$r_cases, reproducibility_results$r_controls, use = "complete.obs"))
      summary_df <- data.frame(
        Model = "ischemic_vs_control",
        Initial_IDPs = n_initial,
        After_r_ge_0_5 = n_after_r,
        After_r_and_n_ge_10 = n_after_both,
        r_vector_correlation_cases_vs_controls = corr_groups,
        stringsAsFactors = FALSE
      )
      utils::write.csv(summary_df, table_out_path("Diagnostics", "Diagnostics_ischemic_vs_control_IDP_Reproducibility_Summary.csv"), row.names = FALSE)
      cat(sprintf("[可重复性] 初始IDP数: %d; r>=0.5后: %d; 再筛n>=10后: %d; r向量相关(cases vs controls): %.3f\n",
                  n_initial, n_after_r, n_after_both, corr_groups))

      plot_scan_rescan_reproducibility_by_idp_class()

      # ============ 新增：依据 ischemic_vs_control 的可重复性结果对脑影像IDP进行全局筛选 ============
      keep_bases <- reproducibility_results$base_name[
        reproducibility_results$variable_type == "Brain_Imaging" & reproducibility_results$pass_both
      ]
      n_brain_before <- sum(vapply(variable_pairs, function(p) p$variable_type == "Brain_Imaging", logical(1)))
      variable_pairs <- c(
        Filter(function(p) p$variable_type == "Cognitive", variable_pairs),
        Filter(function(p) p$variable_type == "Brain_Imaging" && p$base_name %in% keep_bases, variable_pairs)
      )
      n_brain_after <- sum(vapply(variable_pairs, function(p) p$variable_type == "Brain_Imaging", logical(1)))
      utils::write.csv(
        reproducibility_results[
          reproducibility_results$variable_type == "Brain_Imaging" & reproducibility_results$pass_both,
        ],
        table_out_path("Diagnostics", "Filtered_Brain_IDPs_from_Ischemic_vs_Control.csv"), row.names = FALSE
      )
      cat(sprintf("[可重复性筛选] 脑影像IDP配对从 %d → %d（依据 ischemic_vs_control）\n", n_brain_before, n_brain_after))
      # 依据筛选后的集合，重跑横断诊断以反映最终分析集
      write_cross_sectional_diagnostics(analysis_data_list, variable_pairs)

      TRUE
    } else {
      cat("[可重复性] 无可评估的脑影像IDP。\n")
    }
  }

  write_onepage_events_and_longitudinal <- function(analysis_data_list, variable_pairs) {
    new_mi <- if (exists("ischemic_cohort_complete") && ("incident_MI" %in% names(ischemic_cohort_complete))) sum(ischemic_cohort_complete$incident_MI, na.rm = TRUE) else NA_integer_
    new_chronic <- if (exists("ischemic_cohort_complete") && ("incident_chronic_ischemic" %in% names(ischemic_cohort_complete))) sum(ischemic_cohort_complete$incident_chronic_ischemic, na.rm = TRUE) else NA_integer_
    new_any <- if (exists("ischemic_cohort_complete") && ("incident_any_ischemic" %in% names(ischemic_cohort_complete))) sum(ischemic_cohort_complete$incident_any_ischemic, na.rm = TRUE) else NA_integer_
    out_rows <- lapply(names(analysis_data_list), function(model_label) {
      dat <- analysis_data_list[[model_label]]
      pairs_eval <- Filter(function(p) p$variable_type == "Brain_Imaging", variable_pairs)
      pairs_eval <- Filter(function(p) all(c(p$idp1_var, p$idp2_var) %in% names(dat)), pairs_eval)
      n_total <- nrow(dat)
      if (length(pairs_eval) == 0) {
        data.frame(
          Model = model_label,
          n_total = n_total,
          Pairs_Brain = 0,
          Mean_coverage_pair = NA_real_,
          Median_coverage_pair = NA_real_,
          Mean_coverage_model = NA_real_,
          Median_coverage_model = NA_real_,
          New_MI = new_mi,
          New_Chronic_Ischemic = new_chronic,
          New_Any_Ischemic = new_any,
          stringsAsFactors = FALSE
        )
      } else {
        age_ok <- if ("Age2" %in% names(dat)) !is.na(dat$Age2) else rep(TRUE, nrow(dat))
        group_ok <- if ("case_control_binary" %in% names(dat)) !is.na(dat$case_control_binary) else rep(TRUE, nrow(dat))
        n_pair_vec <- vapply(pairs_eval, function(p) {
          x1 <- suppressWarnings(as.numeric(dat[[p$idp1_var]]))
          x2 <- suppressWarnings(as.numeric(dat[[p$idp2_var]]))
          sum(!is.na(x1) & !is.na(x2))
        }, integer(1))
        n_model_vec <- vapply(pairs_eval, function(p) {
          x1 <- suppressWarnings(as.numeric(dat[[p$idp1_var]]))
          x2 <- suppressWarnings(as.numeric(dat[[p$idp2_var]]))
          pair_ok <- !is.na(x1) & !is.na(x2)
          sum(pair_ok & age_ok & group_ok)
        }, integer(1))
        cov_pair <- n_pair_vec / n_total
        cov_model <- n_model_vec / n_total
        data.frame(
          Model = model_label,
          n_total = n_total,
          Pairs_Brain = length(pairs_eval),
          Mean_coverage_pair = round(mean(cov_pair, na.rm = TRUE), 4),
          Median_coverage_pair = round(stats::median(cov_pair, na.rm = TRUE), 4),
          Mean_coverage_model = round(mean(cov_model, na.rm = TRUE), 4),
          Median_coverage_model = round(stats::median(cov_model, na.rm = TRUE), 4),
          New_MI = new_mi,
          New_Chronic_Ischemic = new_chronic,
          New_Any_Ischemic = new_any,
          stringsAsFactors = FALSE
        )
      }
    })
    df <- do.call(rbind, out_rows)
    utils::write.csv(df, table_out_path("Diagnostics", "Diagnostics_OnePage_Events_and_Longitudinal_Sample.csv"), row.names = FALSE)
    df
  }

  build_viz_idps_for_all_models <- function(analysis_data_list, variable_pairs) {
    out_rows <- list()
    for (model_label in names(analysis_data_list)) {
      dat <- analysis_data_list[[model_label]]
      pairs_eval <- Filter(function(p) p$variable_type == "Brain_Imaging", variable_pairs)
      if (length(pairs_eval) == 0) next
      rep_df <- compute_idp_reproducibility(dat, variable_pairs, group_var = "comparison_group", r_threshold = 0.5, min_complete = 10, restrict_type = c("Brain_Imaging"))
      if (is.null(rep_df) || nrow(rep_df) == 0) next
      rep_df <- rep_df[rep_df$variable_type == "Brain_Imaging" & is.finite(rep_df$r_mean), ]
      rep_df <- rep_df[!is.na(rep_df$idp1_var) & !is.na(rep_df$idp2_var), ]
      n_samp <- nrow(dat)
      ratio <- suppressWarnings(as.numeric(getOption("viz.coverage.min", Sys.getenv("VIZ_COVERAGE_MIN", "0.85"))))
      if (!is.finite(ratio) || ratio <= 0 || ratio > 1) ratio <- 0.85
      coverage_min <- floor(ratio * n_samp)
      target_min <- suppressWarnings(as.integer(getOption("viz.target.min", Sys.getenv("VIZ_TARGET_MIN", "24"))))
      if (!is.finite(target_min) || target_min < 0) target_min <- 24
      target_max <- suppressWarnings(as.integer(getOption("viz.target.max", Sys.getenv("VIZ_TARGET_MAX", as.character(nrow(rep_df))))))
      if (!is.finite(target_max) || target_max < 1) target_max <- max(nrow(rep_df), target_min)
      # 仅要求可重复性（pass_both），不再使用硬性覆盖阈值；排序优先高覆盖
      primary <- rep_df[rep_df$pass_both, ]
      primary <- primary[order(-primary$n_complete, -primary$r_mean), ]
      primary <- primary[!duplicated(primary$base_name), ]
      if (nrow(primary) < target_min) {
        need_extra <- target_min - nrow(primary)
        if (need_extra > 0) {
        pool <- rep_df[rep_df$pass_both & !(rep_df$base_name %in% unique(primary$base_name)), ]
        pool <- pool[order(-pool$n_complete, -pool$r_mean), ]
          primary <- rbind(primary, head(pool, need_extra))
          primary <- primary[!duplicated(primary$base_name), ]
        }
        primary <- head(primary, target_max)
      }
      one <- if (nrow(primary) > 0) {
        data.frame(
          model_comparison = rep(model_label, nrow(primary)),
          idp1_variable = primary$idp1_var,
          idp2_variable = primary$idp2_var,
          variable_name = primary$base_name,
          variable_type = "Brain_Imaging",
          n_complete = primary$n_complete,
          n_total = n_samp,
          coverage = round(primary$n_complete / n_samp, 4),
          stringsAsFactors = FALSE
        )
      } else {
        data.frame(
          model_comparison = character(0),
          idp1_variable = character(0),
          idp2_variable = character(0),
          variable_name = character(0),
          variable_type = character(0),
          n_complete = integer(0),
          n_total = integer(0),
          coverage = numeric(0),
          stringsAsFactors = FALSE
        )
      }
      out_rows[[model_label]] <- one
    }
    TRUE
    invisible(TRUE)
  }

  build_viz_idps_for_all_models(analysis_data_list, variable_pairs)
  write_onepage_events_and_longitudinal(analysis_data_list, variable_pairs)

  all_four_model_results <- list()

## 进度条：逐模型Z统计量分析（健壮化：处理模型数量为0或1的情况）
n_models <- length(analysis_data_list)
if (n_models == 0) {
  cat("[Z-分析] 无可分析模型（analysis_data_list 为空），跳过该步骤。\n")
} else {
  pb_models <- txtProgressBar(min = 0, max = n_models, style = 3)
  idx_models <- 0
  for(model_name in names(analysis_data_list)) {
    idx_models <- idx_models + 1
    setTxtProgressBar(pb_models, idx_models)
    cat("\n", rep("=", 50), "\n")
    cat("开始分析模型:", model_name, "\n")
    cat(rep("=", 50), "\n")
  
    # 构建该模型的非影像协变量列表（仅使用可用列；采用LOCF命名）
    data_model <- analysis_data_list[[model_name]]
  base_non_img <- c("townsend_index_baseline","bmi_baseline","systolic_bp_baseline","diastolic_bp_baseline",
                    "diabetes_factor_baseline","smoking_factor_baseline","alcohol_factor_baseline",
                    "waist_circumference_baseline")
  non_img_for_model <- base_non_img %>% intersect(names(data_model))
  
  results_full <- z_statistics_four_model_analysis(
    data = data_model,
    pairs = variable_pairs,
    model_name = model_name,
    drop_demographics = FALSE,
    non_imaging_covariates = non_img_for_model
  )
  
  results_simplified <- z_statistics_four_model_analysis(
    data = data_model,
    pairs = variable_pairs,
    model_name = paste0(model_name, "_no_sex_ethnicity"),
    drop_demographics = TRUE,
    non_imaging_covariates = non_img_for_model
  )
  
  # 基于口径开关构建“显著对”作为衰减分析的前置门槛
  # 确保具备原始p与FDR列
  if (!"p_value_unadjusted" %in% names(results_full)) {
    results_full$p_value_unadjusted <- results_full$p_value_without_non_imaging
  }
  if (!"p_value_fdr" %in% names(results_full)) {
    results_full$p_value_fdr <- stats::p.adjust(results_full$p_value_unadjusted, method = "BH")
  }
  if (!"significant_fdr" %in% names(results_full)) {
    results_full$significant_fdr <- results_full$p_value_fdr < 0.05
  }
  
  if (!"p_value_unadjusted" %in% names(results_simplified)) {
    results_simplified$p_value_unadjusted <- results_simplified$p_value_without_non_imaging
  }
  if (!"p_value_fdr" %in% names(results_simplified)) {
    results_simplified$p_value_fdr <- stats::p.adjust(results_simplified$p_value_unadjusted, method = "BH")
  }
  if (!"significant_fdr" %in% names(results_simplified)) {
    results_simplified$significant_fdr <- results_simplified$p_value_fdr < 0.05
  }
  
  gate_mask_full <- rep(TRUE, nrow(results_full))
  gate_mask_simpl <- rep(TRUE, nrow(results_simplified))
  if (identical(attenuation_count_gate_mode, "raw")) {
    gate_mask_full <- results_full$p_value_unadjusted < 0.05
    gate_mask_simpl <- results_simplified$p_value_unadjusted < 0.05
  } else if (identical(attenuation_count_gate_mode, "fdr")) {
    gate_mask_full <- results_full$significant_fdr
    gate_mask_simpl <- results_simplified$significant_fdr
  }
  signif_pairs_full <- variable_pairs
  signif_pairs_simpl <- variable_pairs
  # 仅保留符合门槛的IDP对
  if (length(variable_pairs) > 0) {
    allowed_keys_full <- paste(results_full$idp1_variable[gate_mask_full], results_full$idp2_variable[gate_mask_full], sep = "||")
    allowed_keys_simpl <- paste(results_simplified$idp1_variable[gate_mask_simpl], results_simplified$idp2_variable[gate_mask_simpl], sep = "||")
    signif_pairs_full <- Filter(function(p) paste(p$idp1_var, p$idp2_var, sep = "||") %in% allowed_keys_full, variable_pairs)
    signif_pairs_simpl <- Filter(function(p) paste(p$idp1_var, p$idp2_var, sep = "||") %in% allowed_keys_simpl, variable_pairs)
  }
  
  # 计算并保存Z值衰减结果
  attenuation_full <- z_attenuation_by_non_imaging(
    data = data_model,
    pairs = signif_pairs_full,
    model_name = model_name,
    drop_demographics = FALSE,
    non_imaging_covariates = non_img_for_model
  )
  if(!is.null(attenuation_full) && nrow(attenuation_full) > 0) {
    # 扩展模型的FDR/FWE与稳定/独立/部分中介标记
    attenuation_full$p_value_with_non_imaging_fdr <- stats::p.adjust(attenuation_full$p_value_with_non_imaging, method = "BH")
    attenuation_full$p_value_with_non_imaging_fwe <- compute_fwe_pvalues_for_attenuation(data_model, attenuation_full, drop_demographics = FALSE, n_perm = 50)
    # 稳定/独立仅依据 |Δz%| < 25% 判定
    attenuation_full$stable_independent <- (abs(attenuation_full$z_change_percent) < 25)
    attenuation_full$partially_mediated <- !attenuation_full$stable_independent
    # 聚合到IDP对层面
    full_pair_flags <- attenuation_full %>%
      dplyr::group_by(model_comparison, idp1_variable, idp2_variable) %>%
      dplyr::summarise(
        independent_of_nonimaging = all(stable_independent),
        partially_mediated_by = paste(non_imaging_variable[!stable_independent], collapse = "; "),
        .groups = "drop"
      )
    # 保存并连接到完整模型结果
    write.csv(attenuation_full, table_out_path(file.path("Four_Model","PerModel"), paste0("Four_Model_Z_Attenuation_", model_name, ".csv")), row.names = FALSE)
    if(!is.null(results_full) && nrow(results_full) > 0) {
      results_full <- results_full %>% dplyr::left_join(full_pair_flags, by = c("model_comparison","idp1_variable","idp2_variable"))
      results_full$partially_mediated_any <- !dplyr::coalesce(results_full$independent_of_nonimaging, FALSE)
    }
  }
  
  attenuation_simplified <- z_attenuation_by_non_imaging(
    data = data_model,
    pairs = signif_pairs_simpl,
    model_name = paste0(model_name, "_no_sex_ethnicity"),
    drop_demographics = TRUE,
    non_imaging_covariates = non_img_for_model
  )
  if(!is.null(attenuation_simplified) && nrow(attenuation_simplified) > 0) {
    # 扩展模型的FDR/FWE与稳定/独立/部分中介标记（精简模型）
    attenuation_simplified$p_value_with_non_imaging_fdr <- stats::p.adjust(attenuation_simplified$p_value_with_non_imaging, method = "BH")
    attenuation_simplified$p_value_with_non_imaging_fwe <- compute_fwe_pvalues_for_attenuation(data_model, attenuation_simplified, drop_demographics = TRUE, n_perm = 50)
    # 稳定/独立仅依据 |Δz%| < 25% 判定
    attenuation_simplified$stable_independent <- (abs(attenuation_simplified$z_change_percent) < 25)
    attenuation_simplified$partially_mediated <- !attenuation_simplified$stable_independent
    # 聚合到IDP对层面
    simple_pair_flags <- attenuation_simplified %>%
      dplyr::group_by(model_comparison, idp1_variable, idp2_variable) %>%
      dplyr::summarise(
        independent_of_nonimaging = all(stable_independent),
        partially_mediated_by = paste(non_imaging_variable[!stable_independent], collapse = "; "),
        .groups = "drop"
      )
    # 保存并连接到精简模型结果
    write.csv(attenuation_simplified, table_out_path(file.path("Four_Model","PerModel"), paste0("Four_Model_Z_Attenuation_", model_name, "_no_sex_ethnicity.csv")), row.names = FALSE)
    if(!is.null(results_simplified) && nrow(results_simplified) > 0) {
      results_simplified <- results_simplified %>% dplyr::left_join(simple_pair_flags, by = c("model_comparison","idp1_variable","idp2_variable"))
      results_simplified$partially_mediated_any <- !dplyr::coalesce(results_simplified$independent_of_nonimaging, FALSE)
    }
  }
  
  # 保存与汇总：完整模型（加入FDR与置换FWE校正）
  if(!is.null(results_full) && nrow(results_full) > 0) {
    # FDR校正在模型内的未校正p值上
    results_full$p_value_unadjusted <- results_full$p_value_without_non_imaging
    results_full$p_value_fdr <- stats::p.adjust(results_full$p_value_unadjusted, method = "BH")
    results_full$significant_fdr <- results_full$p_value_fdr < 0.05
    # 置换-FWE（max-|t|）
    results_full$p_value_fwe <- compute_fwe_pvalues_for_results(data_model, results_full, drop_demographics = FALSE, n_perm = 50)
    results_full$significant_fwe <- results_full$p_value_fwe < 0.05
    # 独立效应：FDR<0.05 或 原始p<0.05
    results_full$independent_effect <- (results_full$significant_fdr | (results_full$p_value_unadjusted < 0.05))
    results_full$independent_of_confounds <- results_full$independent_effect
    # 统一列集：确保所有模型结果拥有相同的衰减衍生列（即便为NA），避免rbind列不一致
    required_cols <- c("independent_of_nonimaging","partially_mediated_by","partially_mediated_any")
    for (cc in required_cols) {
      if (!cc %in% names(results_full)) {
        results_full[[cc]] <- NA
      }
    }
    all_four_model_results[[model_name]] <- results_full
    write.csv(results_full, table_out_path(file.path("Four_Model","PerModel"), paste0("Four_Model_Z_Analysis_", model_name, ".csv")), row.names = FALSE)
    cat("\n模型", model_name, "（含性别/种族）分析完成:\n")
    cat("- 总分析变量:", nrow(results_full), "个\n")
    cat("- 认知变量:", sum(results_full$variable_type == "Cognitive"), "个\n")
    cat("- 脑影像变量:", sum(results_full$variable_type == "Brain_Imaging"), "个\n")
    cat("- 原始显著(raw p<0.05):", sum(results_full$p_value_unadjusted < 0.05, na.rm = TRUE), "个\n")
    cat("- 显著(FDR<0.05):", sum(results_full$significant_fdr, na.rm = TRUE), "个\n")
    cat("- 显著(FWE<0.05):", sum(results_full$significant_fwe, na.rm = TRUE), "个\n")
    cat("- 独立效应IDP:", sum(results_full$independent_effect, na.rm = TRUE), "个\n")
    if("independent_of_nonimaging" %in% names(results_full)) {
      gate_mask_full <- rep(TRUE, nrow(results_full))
      if (identical(attenuation_count_gate_mode, "raw")) {
        gate_mask_full <- results_full$p_value_unadjusted < 0.05
      } else if (identical(attenuation_count_gate_mode, "fdr")) {
        gate_mask_full <- results_full$significant_fdr
      }
      stable_count_full <- sum(dplyr::coalesce(results_full$independent_of_nonimaging, FALSE) & dplyr::coalesce(gate_mask_full, FALSE), na.rm = TRUE)
      mediated_count_full <- sum(dplyr::coalesce(results_full$partially_mediated_any, FALSE) & dplyr::coalesce(gate_mask_full, FALSE), na.rm = TRUE)
      cat("- 稳定/独立(对非成像变量):", stable_count_full, "个\n")
      cat("- 部分中介/受影响:", mediated_count_full, "个\n")
    }
    cat("- 大效应(|d|>0.5):", sum(abs(results_full$cohens_d) > 0.5, na.rm = TRUE), "个\n")
    ranked_full <- results_full %>%
      dplyr::filter(independent_effect == TRUE) %>%
      dplyr::mutate(sort_key = dplyr::coalesce(abs(z_statistic_with_non_imaging), abs(z_statistic_without_non_imaging))) %>%
      dplyr::arrange(dplyr::desc(sort_key)) %>%
      dplyr::mutate(rank_index = dplyr::row_number()) %>%
      dplyr::select(-sort_key)
    write.csv(ranked_full, table_out_path(file.path("Four_Model","PerModel"), paste0("Ranked_Meaningful_IDPs_", model_name, ".csv")), row.names = FALSE)
  } else {
    cat("模型", model_name, "（含性别/种族）无有效结果\n")
  }
  
  # 保存与汇总：精简模型（去除性别/种族，加入FDR与置换FWE校正）
  if(!is.null(results_simplified) && nrow(results_simplified) > 0) {
    results_simplified$p_value_unadjusted <- results_simplified$p_value_without_non_imaging
    results_simplified$p_value_fdr <- stats::p.adjust(results_simplified$p_value_unadjusted, method = "BH")
    results_simplified$significant_fdr <- results_simplified$p_value_fdr < 0.05
    results_simplified$p_value_fwe <- compute_fwe_pvalues_for_results(data_model, results_simplified, drop_demographics = TRUE, n_perm = 50)
    results_simplified$significant_fwe <- results_simplified$p_value_fwe < 0.05
    results_simplified$independent_effect <- (results_simplified$significant_fdr | (results_simplified$p_value_unadjusted < 0.05))
    results_simplified$independent_of_confounds <- results_simplified$independent_effect
    # 统一列集：确保精简模型结果拥有相同的衰减衍生列（即便为NA），避免rbind列不一致
    required_cols <- c("independent_of_nonimaging","partially_mediated_by","partially_mediated_any")
    for (cc in required_cols) {
      if (!cc %in% names(results_simplified)) {
        results_simplified[[cc]] <- NA
      }
    }
    all_four_model_results[[paste0(model_name, "_no_sex_ethnicity")]] <- results_simplified
    write.csv(results_simplified, table_out_path(file.path("Four_Model","PerModel"), paste0("Four_Model_Z_Analysis_", model_name, "_no_sex_ethnicity.csv")), row.names = FALSE)
    cat("\n模型", model_name, "_no_sex_ethnicity（精简）分析完成:\n")
    cat("- 总分析变量:", nrow(results_simplified), "个\n")
    cat("- 认知变量:", sum(results_simplified$variable_type == "Cognitive"), "个\n")
    cat("- 脑影像变量:", sum(results_simplified$variable_type == "Brain_Imaging"), "个\n")
    cat("- 原始显著(raw p<0.05):", sum(results_simplified$p_value_unadjusted < 0.05, na.rm = TRUE), "个\n")
    cat("- 显著(FDR<0.05):", sum(results_simplified$significant_fdr, na.rm = TRUE), "个\n")
    cat("- 显著(FWE<0.05):", sum(results_simplified$significant_fwe, na.rm = TRUE), "个\n")
    cat("- 独立效应IDP:", sum(results_simplified$independent_effect, na.rm = TRUE), "个\n")
    if("independent_of_nonimaging" %in% names(results_simplified)) {
      gate_mask_simpl <- rep(TRUE, nrow(results_simplified))
      if (identical(attenuation_count_gate_mode, "raw")) {
        gate_mask_simpl <- results_simplified$p_value_unadjusted < 0.05
      } else if (identical(attenuation_count_gate_mode, "fdr")) {
        gate_mask_simpl <- results_simplified$significant_fdr
      }
      stable_count_simpl <- sum(dplyr::coalesce(results_simplified$independent_of_nonimaging, FALSE) & dplyr::coalesce(gate_mask_simpl, FALSE), na.rm = TRUE)
      mediated_count_simpl <- sum(dplyr::coalesce(results_simplified$partially_mediated_any, FALSE) & dplyr::coalesce(gate_mask_simpl, FALSE), na.rm = TRUE)
      cat("- 稳定/独立(对非成像变量):", stable_count_simpl, "个\n")
      cat("- 部分中介/受影响:", mediated_count_simpl, "个\n")
    }
    cat("- 大效应(|d|>0.5):", sum(abs(results_simplified$cohens_d) > 0.5, na.rm = TRUE), "个\n")
    ranked_simplified <- results_simplified %>%
      dplyr::filter(independent_effect == TRUE) %>%
      dplyr::mutate(sort_key = dplyr::coalesce(abs(z_statistic_with_non_imaging), abs(z_statistic_without_non_imaging))) %>%
      dplyr::arrange(dplyr::desc(sort_key)) %>%
      dplyr::mutate(rank_index = dplyr::row_number()) %>%
      dplyr::select(-sort_key)
    write.csv(ranked_simplified, table_out_path(file.path("Four_Model","PerModel"), paste0("Ranked_Meaningful_IDPs_", model_name, "_no_sex_ethnicity.csv")), row.names = FALSE)
  } else {
    cat("模型", model_name, "_no_sex_ethnicity（精简）无有效结果\n")
  }
  }
  close(pb_models)
}
## 已在模型循环结束后关闭进度条，避免未定义对象错误

# =================== 第六步：合并结果并筛选独立效应IDP ===================
cat("\n第六步：合并结果并筛选独立效应IDP...\n")

if(length(all_four_model_results) > 0) {
  # 使用按列名对齐的拼接，避免不同模型结果列集不完全一致导致报错
  combined_four_model_results <- dplyr::bind_rows(all_four_model_results)
  combined_four_model_results <- clean_combined_results_df(combined_four_model_results)
  # 注入IDP分类与脑图映射属性
  var_col <- if ("variable_name" %in% names(combined_four_model_results)) "variable_name" else if ("idp2_variable" %in% names(combined_four_model_results)) "idp2_variable" else NA_character_
  if (!is.na(var_col)) {
    combined_four_model_results <- annotate_idp_attributes(combined_four_model_results, var_col = var_col)
  }
  combined_four_model_results <- clean_combined_results_df(combined_four_model_results)
  # 应用计数口径开关，用于综合与分模型统计
  gate_mask_global <- rep(TRUE, nrow(combined_four_model_results))
  if (identical(attenuation_count_gate_mode, "raw")) {
    gate_mask_global <- combined_four_model_results$p_value_without_non_imaging < 0.05
  } else if (identical(attenuation_count_gate_mode, "fdr")) {
    gate_mask_global <- combined_four_model_results$significant_fdr
  }
  write.csv(combined_four_model_results, combined_four_model_results_csv(), row.names = FALSE)
  
  cat("=== 四模型综合结果摘要 ===\n")
  cat("总分析变量数:", nrow(combined_four_model_results), "\n")
  cat("原始p<0.05总数:", sum(combined_four_model_results$p_value_without_non_imaging < 0.05, na.rm = TRUE), "\n")
  cat("FDR<0.05总数:", sum(combined_four_model_results$significant_fdr, na.rm = TRUE), "\n")
  cat("FWE<0.05总数:", sum(combined_four_model_results$significant_fwe, na.rm = TRUE), "\n")
  cat("独立效应IDP总数:", sum(combined_four_model_results$independent_effect, na.rm = TRUE), "\n")
  if("independent_of_nonimaging" %in% names(combined_four_model_results)) {
    stable_total <- sum(dplyr::coalesce(combined_four_model_results$independent_of_nonimaging, FALSE) & dplyr::coalesce(gate_mask_global, FALSE), na.rm = TRUE)
    mediated_total <- sum(dplyr::coalesce(combined_four_model_results$partially_mediated_any, FALSE) & dplyr::coalesce(gate_mask_global, FALSE), na.rm = TRUE)
    cat("稳定/独立(对非成像变量)总数:", stable_total, "\n")
    cat("部分中介/受影响总数:", mediated_total, "\n")
  }
  
  # 按模型统计独立效应IDP
  model_independent_summary <- combined_four_model_results %>%
    dplyr::mutate(
      gate_mask = dplyr::case_when(
        identical(attenuation_count_gate_mode, "raw") ~ (p_value_without_non_imaging < 0.05),
        identical(attenuation_count_gate_mode, "fdr") ~ significant_fdr,
        TRUE ~ TRUE
      )
    ) %>%
    dplyr::group_by(model_comparison, variable_type) %>%
    dplyr::summarise(
      total_vars = dplyr::n(),
      raw_p_lt_05 = sum(p_value_without_non_imaging < 0.05, na.rm = TRUE),
      independent_effects_bh = sum(independent_effect, na.rm = TRUE),
      significant_fdr = sum(significant_fdr, na.rm = TRUE),
      significant_fwe = sum(significant_fwe, na.rm = TRUE),
      independent_of_nonimaging_count = sum(dplyr::coalesce(independent_of_nonimaging, FALSE) & dplyr::coalesce(gate_mask, FALSE), na.rm = TRUE),
      partially_mediated_any_count = sum(dplyr::coalesce(partially_mediated_any, FALSE) & dplyr::coalesce(gate_mask, FALSE), na.rm = TRUE),
      large_effects = sum(abs(cohens_d) > 0.5, na.rm = TRUE),
      large_effects_with_non_imaging = sum(abs(dplyr::coalesce(cohens_d_with_non_imaging, 0)) > 0.5, na.rm = TRUE),
      mean_abs_z = round(mean(abs(z_statistic_without_non_imaging), na.rm = TRUE), 3),
      mean_abs_z_with_non_imaging = round(mean(abs(dplyr::coalesce(z_statistic_with_non_imaging, NA_real_)), na.rm = TRUE), 3),
      .groups = "drop"
    )
  
  print(model_independent_summary)
  write.csv(model_independent_summary, table_out_path("Four_Model", "Four_Model_Independent_Effects_Summary.csv"), row.names = FALSE)
  
  # 筛选独立效应IDP用于可视化（综合集合）
  independent_effect_idps <- combined_four_model_results %>%
    dplyr::filter(independent_effect == TRUE) %>%
    dplyr::arrange(desc(abs(dplyr::coalesce(z_statistic_with_non_imaging, z_statistic_without_non_imaging))))
  # 注解独立效应集合
  if (!is.na(var_col)) {
    independent_effect_idps <- annotate_idp_attributes(independent_effect_idps, var_col = var_col)
  }
  
  if(nrow(independent_effect_idps) > 0) {
    cat("\n发现", nrow(independent_effect_idps), "个独立效应IDP，用于后续可视化分析\n")
    
    # 清理冗余：不再输出认知/脑影像子集CSV，统一用综合CSV
    cat("- 认知独立效应IDP:", sum(independent_effect_idps$variable_type == "Cognitive", na.rm = TRUE), "个\n")
    cat("- 脑影像独立效应IDP:", sum(independent_effect_idps$variable_type == "Brain_Imaging", na.rm = TRUE), "个\n")
    
    write.csv(independent_effect_idps, all_independent_effect_idps_csv(), row.names = FALSE)
    
    # 额外输出严格集合（仅 FDR<0.05），用于对比分析
    strict_fdr_idps <- combined_four_model_results %>%
      dplyr::filter(significant_fdr == TRUE) %>%
      dplyr::arrange(desc(abs(dplyr::coalesce(z_statistic_with_non_imaging, z_statistic_without_non_imaging))))
    if (!is.na(var_col)) {
      strict_fdr_idps <- annotate_idp_attributes(strict_fdr_idps, var_col = var_col)
    }
    if(nrow(strict_fdr_idps) > 0) {
      write.csv(strict_fdr_idps, table_out_path("Four_Model", "Strict_FDR_Independent_IDPs.csv"), row.names = FALSE)
      cat("- 严格集合(FDR<0.05)独立效应IDP:", nrow(strict_fdr_idps), "个\n")
    } else {
      cat("- 严格集合(FDR<0.05)暂无独立效应IDP\n")
    }
    
    
    
    # 显示综合集合的顶级独立效应IDP（前20个）
    cat("\n顶级独立效应IDP（综合，前20个）:\n")
    top_independent <- utils::head(independent_effect_idps, 20)
    for(k in 1:nrow(top_independent)) {
      cat(sprintf("  %d. %s (%s, %s) - |Z|=%.3f, p=%s, d=%.3f\n",
                  k, substr(top_independent$variable_name[k], 1, 35),
                  top_independent$model_comparison[k],
                  top_independent$variable_type[k],
                  abs(dplyr::coalesce(top_independent$z_statistic_with_non_imaging[k], top_independent$z_statistic_without_non_imaging[k])),
                  top_independent$p_formatted[k],
                  dplyr::coalesce(top_independent$cohens_d_with_non_imaging[k], top_independent$cohens_d[k])))
    }

    # 按模型分别选择Top20独立效应IDP并输出与展示
    cat("\n按模型分别选择Top20独立效应IDP并输出...\n")
    base_models <- unique(gsub("_no_sex_ethnicity$", "", combined_four_model_results$model_comparison))
    top20_dir <- get_table_dir("Four_Model", "Top20_Independent_IDPs")

    for (mb in base_models) {
      # Full 模型 Top20（独立效应，按|Z|排序）
      sub_full <- combined_four_model_results %>%
        dplyr::filter(model_comparison == mb, independent_effect == TRUE) %>%
        dplyr::mutate(sort_key = dplyr::coalesce(abs(z_statistic_with_non_imaging), abs(z_statistic_without_non_imaging))) %>%
        dplyr::arrange(dplyr::desc(sort_key)) %>%
        dplyr::mutate(rank_index = dplyr::row_number()) %>%
        dplyr::select(-sort_key)
      if (!is.na(var_col) && nrow(sub_full) > 0) {
        sub_full <- annotate_idp_attributes(sub_full, var_col = var_col)
      }
      top20_full <- utils::head(sub_full, 20)
      full_path <- file.path(top20_dir, paste0("Top20_Independent_IDPs_", mb, ".csv"))
      utils::write.csv(top20_full, full_path, row.names = FALSE)

      # No Sex/Ethnicity 精简模型 Top20
      mb_simpl <- paste0(mb, "_no_sex_ethnicity")
      sub_simpl <- combined_four_model_results %>%
        dplyr::filter(model_comparison == mb_simpl, independent_effect == TRUE) %>%
        dplyr::mutate(sort_key = dplyr::coalesce(abs(z_statistic_with_non_imaging), abs(z_statistic_without_non_imaging))) %>%
        dplyr::arrange(dplyr::desc(sort_key)) %>%
        dplyr::mutate(rank_index = dplyr::row_number()) %>%
        dplyr::select(-sort_key)
      if (!is.na(var_col) && nrow(sub_simpl) > 0) {
        sub_simpl <- annotate_idp_attributes(sub_simpl, var_col = var_col)
      }
      top20_simpl <- utils::head(sub_simpl, 20)
      simpl_path <- file.path(top20_dir, paste0("Top20_Independent_IDPs_", mb, "_no_sex_ethnicity.csv"))
      utils::write.csv(top20_simpl, simpl_path, row.names = FALSE)

      # 控制台展示：每个模型Top20（Full与Simplified）
      cat(sprintf("\n模型 %s - Top20 独立效应IDP（Full）: %d 项\n", mb, nrow(top20_full)))
      if (nrow(top20_full) > 0) {
        for (k in seq_len(nrow(top20_full))) {
          cat(sprintf("  %2d. %s - |Z|=%.3f, p=%s, d=%.3f\n",
                      k,
                      substr(top20_full$variable_name[k], 1, 35),
                      abs(dplyr::coalesce(top20_full$z_statistic_with_non_imaging[k], top20_full$z_statistic_without_non_imaging[k])),
                      dplyr::coalesce(top20_full$p_formatted[k], sprintf("%.3e", top20_full$p_value_without_non_imaging[k])),
                      dplyr::coalesce(top20_full$cohens_d_with_non_imaging[k], top20_full$cohens_d[k])))
        }
      }
      cat(sprintf("模型 %s - Top20 独立效应IDP（No Sex/Ethnicity）: %d 项\n", mb, nrow(top20_simpl)))
      if (nrow(top20_simpl) > 0) {
        for (k in seq_len(nrow(top20_simpl))) {
          cat(sprintf("  %2d. %s - |Z|=%.3f, p=%s, d=%.3f\n",
                      k,
                      substr(top20_simpl$variable_name[k], 1, 35),
                      abs(dplyr::coalesce(top20_simpl$z_statistic_with_non_imaging[k], top20_simpl$z_statistic_without_non_imaging[k])),
                      dplyr::coalesce(top20_simpl$p_formatted[k], sprintf("%.3e", top20_simpl$p_value_without_non_imaging[k])),
                      dplyr::coalesce(top20_simpl$cohens_d_with_non_imaging[k], top20_simpl$cohens_d[k])))
        }
      }
    }

    top20_img_dir <- get_table_dir("Four_Model", "Top20_Brain_Imaging_IDPs")
    safe_abs_col <- function(df, nm) {
      if (is.null(df) || nrow(df) == 0) return(numeric(0))
      if (!(nm %in% names(df))) return(rep(NA_real_, nrow(df)))
      suppressWarnings(abs(as.numeric(df[[nm]])))
    }
    img_all <- combined_four_model_results %>%
      dplyr::filter(variable_type == "Brain_Imaging", !grepl("_no_sex_ethnicity$", model_comparison))
    if (nrow(img_all) > 0) {
      img_all$sort_key <- dplyr::coalesce(
        safe_abs_col(img_all, "z_statistic_with_non_imaging"),
        safe_abs_col(img_all, "z_statistic_without_non_imaging")
      )
      img_all <- img_all[order(-img_all$sort_key), , drop = FALSE]
      img_all$rank_index <- seq_len(nrow(img_all))
      img_all$sort_key <- NULL
    }
    top20_img_all <- utils::head(img_all, 20)
    utils::write.csv(top20_img_all, file.path(top20_img_dir, "Top20_Brain_Imaging_IDPs_AllModels.csv"), row.names = FALSE)

    for (mb in base_models) {
      sub_img_full <- combined_four_model_results %>%
        dplyr::filter(model_comparison == mb, variable_type == "Brain_Imaging")
      if (nrow(sub_img_full) > 0) {
        sub_img_full$sort_key <- dplyr::coalesce(
          safe_abs_col(sub_img_full, "z_statistic_with_non_imaging"),
          safe_abs_col(sub_img_full, "z_statistic_without_non_imaging")
        )
        sub_img_full <- sub_img_full[order(-sub_img_full$sort_key), , drop = FALSE]
        sub_img_full$rank_index <- seq_len(nrow(sub_img_full))
        sub_img_full$sort_key <- NULL
      }
      if (!is.na(var_col) && nrow(sub_img_full) > 0) {
        sub_img_full <- annotate_idp_attributes(sub_img_full, var_col = var_col)
      }
      top20_img_full <- utils::head(sub_img_full, 20)
      utils::write.csv(top20_img_full, file.path(top20_img_dir, paste0("Top20_Brain_Imaging_IDPs_", mb, ".csv")), row.names = FALSE)

      mb_simpl <- paste0(mb, "_no_sex_ethnicity")
      sub_img_simpl <- combined_four_model_results %>%
        dplyr::filter(model_comparison == mb_simpl, variable_type == "Brain_Imaging")
      if (nrow(sub_img_simpl) > 0) {
        sub_img_simpl$sort_key <- dplyr::coalesce(
          safe_abs_col(sub_img_simpl, "z_statistic_with_non_imaging"),
          safe_abs_col(sub_img_simpl, "z_statistic_without_non_imaging")
        )
        sub_img_simpl <- sub_img_simpl[order(-sub_img_simpl$sort_key), , drop = FALSE]
        sub_img_simpl$rank_index <- seq_len(nrow(sub_img_simpl))
        sub_img_simpl$sort_key <- NULL
      }
      if (!is.na(var_col) && nrow(sub_img_simpl) > 0) {
        sub_img_simpl <- annotate_idp_attributes(sub_img_simpl, var_col = var_col)
      }
      top20_img_simpl <- utils::head(sub_img_simpl, 20)
      utils::write.csv(top20_img_simpl, file.path(top20_img_dir, paste0("Top20_Brain_Imaging_IDPs_", mb, "_no_sex_ethnicity.csv")), row.names = FALSE)
    }
    
    
    # =================== 第六.5b步：Age1/Age2横截面GLM（二类模型+置换校正） ===================
    cat("\n第六.5b步：Age1/Age2横截面GLM（Binary与Age-modulated；含FDR与置换FWE）...\n")
    
    # 统一“完整样本”口径：若存在配对时间点，必须两次扫描均非缺失（与纵向一致）
    complete_case_mask_cross <- function(df, v, age_col, require_sex = FALSE) {
      if (!(v %in% names(df))) return(rep(FALSE, nrow(df)))
      base_mask <- !is.na(df[[v]]) & !is.na(df$case_control_binary) & !is.na(df[[age_col]])
      pair_candidates <- character(0)
      if (grepl("Instance(\\.|_|\\s)?2$", v)) {
        vr <- sub("Instance(\\.|_|\\s)?2$", "Instance.3", v)
        pair_candidates <- c(vr, paste0(vr, ".x"), paste0(vr, ".y"))
      } else if (grepl("Instance(\\.|_|\\s)?3(\\.[xy])?$", v)) {
        vr <- sub("Instance(\\.|_|\\s)?3(\\.[xy])?$", "Instance.2", v)
        pair_candidates <- c(vr, paste0(vr, ".x"), paste0(vr, ".y"))
      }
      pair_candidates <- pair_candidates[pair_candidates %in% names(df)]
      if (length(pair_candidates) > 0) {
        ok <- rowSums(!is.na(df[, pair_candidates, drop = FALSE])) > 0
        base_mask <- base_mask & ok
      }
      if (isTRUE(require_sex) && ("sex_factor" %in% names(df))) {
        base_mask <- base_mask & !is.na(df$sex_factor)
      }
      base_mask
    }
    
    # 置换-FWE：横截面GLM（通用，按变量集max-|t|）
compute_fwe_pvalues_for_cross_sectional_glm <- function(df, var_names, instance = 2L, term = c("group", "modulated"), include_sex = TRUE, n_perm = 50, seed = 202401, n_cores = NULL) {
      term <- match.arg(term)
      if (length(var_names) == 0) return(numeric(0))
      set.seed(seed)
      age_col <- if (instance == 2L) "age_at_instance2" else "age_at_instance3"
      has_sex <- include_sex && ("sex_factor" %in% names(df))
      # 观测t值收集
      build_item <- function(v) {
        if (!(v %in% names(df))) return(NULL)
        y <- df[[v]]
        if (!is.numeric(y)) return(NULL)
        cc <- complete_case_mask_cross(df, v, age_col, require_sex = has_sex)
        if (sum(cc) < 3) return(NULL)
        d <- df[cc, c("case_control_binary", age_col, if (has_sex) "sex_factor" else NULL), drop = FALSE]
        d$y <- y[cc]
        d$age <- d[[age_col]]; d$age_sq <- d$age^2
        if (term == "group") {
          if (has_sex) {
            f <- y ~ case_control_binary + age + age_sq + sex_factor
          } else {
            f <- y ~ case_control_binary + age + age_sq
          }
          fit <- tryCatch(stats::lm(f, data = d), error = function(e) NULL)
          if (is.null(fit)) return(NULL)
          ct <- coef(summary(fit)); ridx <- which(rownames(ct) == "case_control_binary")
          if (length(ridx) == 0) return(NULL)
          list(df = d, time = instance, t_obs = as.numeric(ct[ridx, "t value"]))
        } else {
          cc_dm <- d$case_control_binary - mean(d$case_control_binary)
          d$Case_vs_Control_modulated <- cc_dm * (10^(d$age * 0.0524 - 3.27))
          if (has_sex) {
            f <- y ~ Case_vs_Control_modulated + age + age_sq + sex_factor
          } else {
            f <- y ~ Case_vs_Control_modulated + age + age_sq
          }
          fit <- tryCatch(stats::lm(f, data = d), error = function(e) NULL)
          if (is.null(fit)) return(NULL)
          ct <- coef(summary(fit)); ridx <- which(rownames(ct) == "Case_vs_Control_modulated")
          if (length(ridx) == 0) return(NULL)
          list(df = d, time = instance, t_obs = as.numeric(ct[ridx, "t value"]))
        }
      }
      items <- lapply(var_names, build_item)
      keep <- !sapply(items, is.null)
      if (!any(keep)) return(rep(NA_real_, length(var_names)))
      items_idx <- which(keep)
      # 置换：打乱分组，重建回归量，按次计算max |t|并累计超过次数（支持并行）
      exceed_count <- numeric(length(items_idx))
      obs_abs_t <- sapply(items[items_idx], function(it) abs(it$t_obs))
      use_par <- requireNamespace("foreach", quietly = TRUE) && requireNamespace("doParallel", quietly = TRUE)
      cores <- if (is.null(n_cores) || !is.numeric(n_cores) || n_cores < 1) max(1L, parallel::detectCores(logical = TRUE) - 1L) else as.integer(n_cores)
      if (use_par && cores > 1) {
        suppressPackageStartupMessages(library(foreach))
        suppressPackageStartupMessages(library(doParallel))
        cl <- parallel::makeCluster(cores)
        doParallel::registerDoParallel(cl)
        on.exit({ if (exists("cl") && !is.null(cl)) parallel::stopCluster(cl) }, add = TRUE)
        perm_max <- foreach::foreach(b = seq_len(n_perm), .combine = 'c', .packages = c('stats')) %dopar% {
          max_abs_t <- 0
          for (j in seq_along(items_idx)) {
            it <- items[[items_idx[j]]]; d <- it$df
            perm_cc <- sample(d$case_control_binary)
            d$perm_cc <- perm_cc
            d$age <- d$age; d$age_sq <- d$age_sq
            if (term == "group") {
              if (has_sex) {
                f <- y ~ perm_cc + age + age_sq + sex_factor
              } else {
                f <- y ~ perm_cc + age + age_sq
              }
              fitp <- tryCatch(stats::lm(f, data = d), error = function(e) NULL)
              if (!is.null(fitp)) {
                ctp <- coef(summary(fitp)); ridx <- which(rownames(ctp) == "perm_cc")
                if (length(ridx) > 0) {
                  tv <- abs(as.numeric(ctp[ridx, "t value"]))
                  if (is.finite(tv) && tv > max_abs_t) max_abs_t <- tv
                }
              }
            } else {
              cc_dm <- perm_cc - mean(perm_cc)
              d$Case_vs_Control_modulated <- cc_dm * (10^(d$age * 0.0524 - 3.27))
              if (has_sex) {
                f <- y ~ Case_vs_Control_modulated + age + age_sq + sex_factor
              } else {
                f <- y ~ Case_vs_Control_modulated + age + age_sq
              }
              fitp <- tryCatch(stats::lm(f, data = d), error = function(e) NULL)
              if (!is.null(fitp)) {
                ctp <- coef(summary(fitp)); ridx <- which(rownames(ctp) == "Case_vs_Control_modulated")
                if (length(ridx) > 0) {
                  tv <- abs(as.numeric(ctp[ridx, "t value"]))
                  if (is.finite(tv) && tv > max_abs_t) max_abs_t <- tv
                }
              }
            }
          }
          max_abs_t
        }
        exceed_count <- vapply(obs_abs_t, function(obs) sum(perm_max >= obs), integer(1))
      } else {
        for (b in seq_len(n_perm)) {
          max_abs_t <- 0
          for (j in seq_along(items_idx)) {
            it <- items[[items_idx[j]]]; d <- it$df
            # 置换分组标签
            perm_cc <- sample(d$case_control_binary)
            d$perm_cc <- perm_cc
            d$age <- d$age; d$age_sq <- d$age_sq
            if (term == "group") {
              if (has_sex) {
                f <- y ~ perm_cc + age + age_sq + sex_factor
              } else {
                f <- y ~ perm_cc + age + age_sq
              }
              fitp <- tryCatch(stats::lm(f, data = d), error = function(e) NULL)
              if (is.null(fitp)) { next }
              ctp <- coef(summary(fitp)); ridx <- which(rownames(ctp) == "perm_cc")
              if (length(ridx) > 0) {
                tv <- abs(as.numeric(ctp[ridx, "t value"]))
                if (is.finite(tv) && tv > max_abs_t) max_abs_t <- tv
              }
            } else {
              cc_dm <- perm_cc - mean(perm_cc)
              d$Case_vs_Control_modulated <- cc_dm * (10^(d$age * 0.0524 - 3.27))
              if (has_sex) {
                f <- y ~ Case_vs_Control_modulated + age + age_sq + sex_factor
              } else {
                f <- y ~ Case_vs_Control_modulated + age + age_sq
              }
              fitp <- tryCatch(stats::lm(f, data = d), error = function(e) NULL)
              if (is.null(fitp)) { next }
              ctp <- coef(summary(fitp)); ridx <- which(rownames(ctp) == "Case_vs_Control_modulated")
              if (length(ridx) > 0) {
                tv <- abs(as.numeric(ctp[ridx, "t value"]))
                if (is.finite(tv) && tv > max_abs_t) max_abs_t <- tv
              }
            }
          }
          if (is.finite(max_abs_t)) {
            exceed_count <- exceed_count + as.numeric(max_abs_t >= obs_abs_t)
          }
        }
      }
      # FWE p值
      p_fwe <- rep(NA_real_, length(var_names))
      for (j in seq_along(items_idx)) {
        idx <- items_idx[j]
        p_fwe[idx] <- (1 + exceed_count[j]) / (n_perm + 1)
      }
      p_fwe
    }

    # 工具函数：统一扁平化并筛选数值型变量（兼容list/字符向量），并按Instance过滤
    flatten_numeric_vars <- function(df, vars, instance) {
      v_all <- unique(as.character(unlist(vars)))
      inst_regex <- paste0("Instance\\.", instance)
      v <- v_all[grepl(inst_regex, v_all) & v_all %in% names(df)]
      v[vapply(v, function(x) is.numeric(df[[x]]), logical(1))]
    }
    
    # 模型1：Binary组间差异（Age+Age^2+Sex）
cross_sectional_glm_binary <- function(df, vars, instance, type_label, model_label, output_prefix = "CrossSection_GLM_Binary", n_perm = 50, n_cores = NULL) {
      time_tag <- if (instance == 2L) "Age1" else "Age2"
      use_vars <- flatten_numeric_vars(df, vars, instance)
      if (length(use_vars) == 0) { cat(sprintf("[GLM-Binary] %s-%s 无可用变量\n", type_label, time_tag)); return(invisible(NULL)) }
      age_col <- if (instance == 2L) "age_at_instance2" else "age_at_instance3"
      has_sex <- "sex_factor" %in% names(df)
      has_eth <- "ethnicity_factor" %in% names(df)
      compute_row <- function(v) {
        y <- df[[v]]
        cc <- complete_case_mask_cross(df, v, age_col, require_sex = has_sex)
        if (sum(cc) < 3) return(NULL)
        d <- df[cc, c("case_control_binary", age_col, if (has_sex) "sex_factor" else NULL, if (has_eth) "ethnicity_factor" else NULL), drop = FALSE]
        d$y <- y[cc]; d$age <- d[[age_col]]; d$age_sq <- d$age^2
        conf_terms <- paste(c(if (has_sex) "sex_factor" else NULL, if (has_eth) "ethnicity_factor" else NULL), collapse = " + ")
        f <- stats::as.formula(paste("y ~ case_control_binary + age + age_sq", if (nchar(conf_terms) > 0) paste("+", conf_terms) else ""))
        fit <- tryCatch(stats::lm(f, data = d), error = function(e) NULL)
        if (is.null(fit)) return(NULL)
        ct <- coef(summary(fit)); ridx <- which(rownames(ct) == "case_control_binary")
        if (length(ridx) == 0) return(NULL)
        beta <- ct[ridx, "Estimate"]; se <- ct[ridx, "Std. Error"]; tval <- ct[ridx, "t value"]; pval <- ct[ridx, "Pr(>|t|)"]
        sd_y <- sqrt(stats::var(d$y, na.rm = TRUE)); d_es <- if (is.finite(sd_y) && sd_y > 0) beta / sd_y else NA_real_
        data.frame(model_comparison = model_label, variable_type = type_label, timepoint = time_tag,
                   variable_name = v, beta_group = beta, se_group = se, z_statistic = tval, p_value = pval,
                   cohens_d = d_es, n = nrow(d), stringsAsFactors = FALSE)
      }
      use_par <- requireNamespace("foreach", quietly = TRUE) && requireNamespace("doParallel", quietly = TRUE)
      cores <- if (is.null(n_cores) || !is.numeric(n_cores) || n_cores < 1) max(1L, parallel::detectCores(logical = TRUE) - 1L) else as.integer(n_cores)
      if (use_par && cores > 1) {
        suppressPackageStartupMessages(library(foreach))
        suppressPackageStartupMessages(library(doParallel))
        cl <- parallel::makeCluster(cores)
        doParallel::registerDoParallel(cl)
        on.exit({ if (exists("cl") && !is.null(cl)) parallel::stopCluster(cl) }, add = TRUE)
        rows <- foreach::foreach(v = use_vars, .combine = 'c', .multicombine = TRUE, .packages = c('stats'), .export = c('complete_case_mask_cross')) %dopar% {
          list(compute_row(v))
        }
        rows_valid <- rows[!sapply(rows, is.null)]
        if (length(rows_valid) == 0) { cat(sprintf("[GLM-Binary] %s-%s 无有效结果\n", type_label, time_tag)); return(invisible(NULL)) }
        res <- do.call(rbind, rows_valid)
      } else {
        rows <- lapply(use_vars, compute_row)
        rows_valid <- rows[!sapply(rows, is.null)]
        if (length(rows_valid) == 0) { cat(sprintf("[GLM-Binary] %s-%s 无有效结果\n", type_label, time_tag)); return(invisible(NULL)) }
        res <- do.call(rbind, rows_valid)
      }
      if (is.null(res) || nrow(res) == 0) { cat(sprintf("[GLM-Binary] %s-%s 无有效结果\n", type_label, time_tag)); return(invisible(NULL)) }
      res$p_fdr <- stats::p.adjust(res$p_value, method = "BH")
      res$p_fwe <- compute_fwe_pvalues_for_cross_sectional_glm(df, res$variable_name, instance = instance, term = "group", include_sex = has_sex, n_perm = n_perm, n_cores = n_cores)
      res$significant_fdr <- res$p_fdr < 0.05
      res$significant_fwe <- res$p_fwe < 0.05
      out_path <- paste0(output_prefix, "_", type_label, "_", time_tag, "_", model_label, ".csv")
      utils::write.csv(res, out_path, row.names = FALSE)
      cat(sprintf("[GLM-Binary] 已输出 %s\n", out_path))
      invisible(res)
    }
    
    # 模型2：Age-modulated组间差异（demeaned(Group)×10^(Age*0.0524−3.27) + Age + Age^2 + Sex）
cross_sectional_glm_modulated <- function(df, vars, instance, type_label, model_label, output_prefix = "CrossSection_GLM_Modulated", n_perm = 50, n_cores = NULL) {
      time_tag <- if (instance == 2L) "Age1" else "Age2"
      use_vars <- flatten_numeric_vars(df, vars, instance)
      if (length(use_vars) == 0) { cat(sprintf("[GLM-Modulated] %s-%s 无可用变量\n", type_label, time_tag)); return(invisible(NULL)) }
      age_col <- if (instance == 2L) "age_at_instance2" else "age_at_instance3"
      has_sex <- "sex_factor" %in% names(df)
      has_eth <- "ethnicity_factor" %in% names(df)
      compute_row <- function(v) {
        y <- df[[v]]
        cc <- complete_case_mask_cross(df, v, age_col, require_sex = has_sex)
        if (sum(cc) < 3) return(NULL)
        d <- df[cc, c("case_control_binary", age_col, if (has_sex) "sex_factor" else NULL, if (has_eth) "ethnicity_factor" else NULL), drop = FALSE]
        d$y <- y[cc]; d$age <- d[[age_col]]; d$age_sq <- d$age^2
        cc_dm <- d$case_control_binary - mean(d$case_control_binary)
        d$Case_vs_Control_modulated <- cc_dm * (10^(d$age * 0.0524 - 3.27))
        conf_terms <- paste(c(if (has_sex) "sex_factor" else NULL, if (has_eth) "ethnicity_factor" else NULL), collapse = " + ")
        f <- stats::as.formula(paste("y ~ Case_vs_Control_modulated + age + age_sq", if (nchar(conf_terms) > 0) paste("+", conf_terms) else ""))
        fit <- tryCatch(stats::lm(f, data = d), error = function(e) NULL)
        if (is.null(fit)) return(NULL)
        ct <- coef(summary(fit)); ridx <- which(rownames(ct) == "Case_vs_Control_modulated")
        if (length(ridx) == 0) return(NULL)
        beta <- ct[ridx, "Estimate"]; se <- ct[ridx, "Std. Error"]; tval <- ct[ridx, "t value"]; pval <- ct[ridx, "Pr(>|t|)"]
        sd_y <- sqrt(stats::var(d$y, na.rm = TRUE)); d_es <- if (is.finite(sd_y) && sd_y > 0) beta / sd_y else NA_real_
        data.frame(model_comparison = model_label, variable_type = type_label, timepoint = time_tag,
                   variable_name = v, beta_modulated = beta, se_modulated = se, z_statistic = tval, p_value = pval,
                   cohens_d = d_es, n = nrow(d), stringsAsFactors = FALSE)
      }
      use_par <- requireNamespace("foreach", quietly = TRUE) && requireNamespace("doParallel", quietly = TRUE)
      cores <- if (is.null(n_cores) || !is.numeric(n_cores) || n_cores < 1) max(1L, parallel::detectCores(logical = TRUE) - 1L) else as.integer(n_cores)
      if (use_par && cores > 1) {
        suppressPackageStartupMessages(library(foreach))
        suppressPackageStartupMessages(library(doParallel))
        cl <- parallel::makeCluster(cores)
        doParallel::registerDoParallel(cl)
        on.exit({ if (exists("cl") && !is.null(cl)) parallel::stopCluster(cl) }, add = TRUE)
        rows <- foreach::foreach(v = use_vars, .combine = 'c', .multicombine = TRUE, .packages = c('stats'), .export = c('complete_case_mask_cross')) %dopar% {
          list(compute_row(v))
        }
        rows_valid <- rows[!sapply(rows, is.null)]
        if (length(rows_valid) == 0) { cat(sprintf("[GLM-Modulated] %s-%s 无有效结果\n", type_label, time_tag)); return(invisible(NULL)) }
        res <- do.call(rbind, rows_valid)
      } else {
        rows <- lapply(use_vars, compute_row)
        rows_valid <- rows[!sapply(rows, is.null)]
        if (length(rows_valid) == 0) { cat(sprintf("[GLM-Modulated] %s-%s 无有效结果\n", type_label, time_tag)); return(invisible(NULL)) }
        res <- do.call(rbind, rows_valid)
      }
      if (is.null(res) || nrow(res) == 0) { cat(sprintf("[GLM-Modulated] %s-%s 无有效结果\n", type_label, time_tag)); return(invisible(NULL)) }
      res$p_fdr <- stats::p.adjust(res$p_value, method = "BH")
      res$p_fwe <- compute_fwe_pvalues_for_cross_sectional_glm(df, res$variable_name, instance = instance, term = "modulated", include_sex = has_sex, n_perm = n_perm, n_cores = n_cores)
      res$significant_fdr <- res$p_fdr < 0.05
      res$significant_fwe <- res$p_fwe < 0.05
      out_path <- paste0(output_prefix, "_", type_label, "_", time_tag, "_", model_label, ".csv")
      utils::write.csv(res, out_path, row.names = FALSE)
      cat(sprintf("[GLM-Modulated] 已输出 %s\n", out_path))
      invisible(res)
    }
    
run_cross_sectional_glm_for_model <- function(model_label = "ischemic_vs_control", n_perm = 50, n_cores = NULL) {
      if (!exists("analysis_data_list") || is.null(analysis_data_list[[model_label]])) {
        cat(sprintf("[GLM] 模型 %s 数据不存在，跳过\n", model_label)); return(invisible(NULL))
      }
      df <- analysis_data_list[[model_label]]
      if (!exists("cognitive_vars") || !exists("brain_imaging_vars")) {
        cat("[GLM] 变量列表不存在，跳过\n"); return(invisible(NULL))
      }
      # 统一入口：给定变量集与类别标签，运行Age1/Age2两类模型（Binary/Modulated）
      run_cross_section_set <- function(df, vars, type_label, model_label, n_perm = 50, n_cores = NULL) {
        cross_sectional_glm_binary(df, vars, 2L, type_label, model_label, n_perm = n_perm, n_cores = n_cores)
        cross_sectional_glm_binary(df, vars, 3L, type_label, model_label, n_perm = n_perm, n_cores = n_cores)
        cross_sectional_glm_modulated(df, vars, 2L, type_label, model_label, n_perm = n_perm, n_cores = n_cores)
        cross_sectional_glm_modulated(df, vars, 3L, type_label, model_label, n_perm = n_perm, n_cores = n_cores)
      }
      run_cross_section_set(df, cognitive_vars, "Cognitive", model_label, n_perm = n_perm, n_cores = n_cores)
      run_cross_section_set(df, brain_imaging_vars, "Brain_Imaging", model_label, n_perm = n_perm, n_cores = n_cores)
    }
    
    # 运行：缺血性心肌病 vs 对照 的横截面GLM
run_cross_sectional_glm_for_model("ischemic_vs_control", n_perm = 50, n_cores = max(1L, parallel::detectCores(logical = TRUE) - 1L))
    # 运行：其余比较模型的横截面GLM（Binary 与 Age-modulated；含FDR与置换FWE）
run_cross_sectional_glm_for_model("mi_vs_control", n_perm = 50, n_cores = max(1L, parallel::detectCores(logical = TRUE) - 1L))
run_cross_sectional_glm_for_model("chronic_vs_control", n_perm = 50, n_cores = max(1L, parallel::detectCores(logical = TRUE) - 1L))
    
    # =================== 第六.6步：摘要导出（缺失值、显著IDP变化率） ===================
    cat("\n第六.6步：导出摘要（缺失值、显著IDP变化率）...\n")
    
    summarize_missing_values_for_model <- function(model_label = "ischemic_vs_control", force = FALSE) {
      if (!exists("analysis_data_list") || is.null(analysis_data_list[[model_label]])) {
        cat(sprintf("[摘要-缺失值] 模型 %s 数据不存在，跳过\n", model_label)); return(invisible(NULL))
      }
      df <- analysis_data_list[[model_label]]
      if (!exists("cognitive_vars") || !exists("brain_imaging_vars")) {
        cat("[摘要-缺失值] 变量列表不存在，跳过\n"); return(invisible(NULL))
      }
      get_vars <- function(vars, instance) {
        # 兼容字符向量或列表输入，统一扁平化为字符向量
        v_all <- unique(as.character(unlist(vars)))
        rgx <- paste0("Instance\\.", instance)
        v <- v_all[grepl(rgx, v_all) & v_all %in% names(df)]
        v[sapply(v, function(x) is.numeric(df[[x]]))]
      }
      sets <- list(
        list(type = "Cognitive", time = "Age1", vars = get_vars(cognitive_vars, 2)),
        list(type = "Cognitive", time = "Age2", vars = get_vars(cognitive_vars, 3)),
        list(type = "Brain_Imaging", time = "Age1", vars = get_vars(brain_imaging_vars, 2)),
        list(type = "Brain_Imaging", time = "Age2", vars = get_vars(brain_imaging_vars, 3))
      )
      overview_rows <- list()
      for (s in sets) {
        vars <- s$vars
        if (length(vars) == 0) { cat(sprintf("[摘要-缺失值] %s-%s 无可用变量，跳过\n", s$type, s$time)); next }
        # 若已存在且不强制，复用已导出的CSV以避免重复计算
        per_path <- paste0("Missingness_PerVariable_", s$type, "_", s$time, "_", model_label, ".csv")
        if (!force && file.exists(per_path)) {
          cat(sprintf("[摘要-缺失值] 复用已有 %s\n", per_path))
          per_var_cached <- tryCatch(utils::read.csv(per_path, stringsAsFactors = FALSE), error = function(e) NULL)
          if (!is.null(per_var_cached) && nrow(per_var_cached)) {
            overview_rows[[length(overview_rows) + 1]] <- data.frame(
              model_comparison = model_label,
              variable_type = s$type,
              timepoint = s$time,
              var_count = nrow(per_var_cached),
              mean_missing_rate = mean(per_var_cached$missing_rate, na.rm = TRUE),
              median_missing_rate = stats::median(per_var_cached$missing_rate, na.rm = TRUE),
              max_missing_rate = max(per_var_cached$missing_rate, na.rm = TRUE),
              top_missing_var = per_var_cached$variable_name[which.max(per_var_cached$missing_rate)],
              stringsAsFactors = FALSE
            )
            next
          }
        }
        has_group <- "case_control_binary" %in% names(df)
        n_cases <- if (has_group) sum(df$case_control_binary == 1, na.rm = TRUE) else NA_integer_
        n_ctrls <- if (has_group) sum(df$case_control_binary == 0, na.rm = TRUE) else NA_integer_
        total <- nrow(df)
        # 向量化一次性计算缺失率（显著提速）
        M <- is.na(df[, vars, drop = FALSE])
        miss_counts <- colSums(M)
        miss_rate <- if (total > 0) miss_counts / total else rep(NA_real_, length(vars))
        if (has_group) {
          mask_cases <- which(df$case_control_binary == 1)
          mask_ctrls <- which(df$case_control_binary == 0)
          miss_cases <- if (length(mask_cases) > 0) colSums(M[mask_cases, , drop = FALSE]) else rep(NA_integer_, length(vars))
          miss_ctrls <- if (length(mask_ctrls) > 0) colSums(M[mask_ctrls, , drop = FALSE]) else rep(NA_integer_, length(vars))
          rate_cases <- if (!is.na(n_cases) && n_cases > 0) miss_cases / n_cases else rep(NA_real_, length(vars))
          rate_ctrls <- if (!is.na(n_ctrls) && n_ctrls > 0) miss_ctrls / n_ctrls else rep(NA_real_, length(vars))
        } else {
          rate_cases <- rep(NA_real_, length(vars)); rate_ctrls <- rep(NA_real_, length(vars))
        }
        per_var <- data.frame(
          model_comparison = model_label,
          variable_type = s$type,
          timepoint = s$time,
          variable_name = vars,
          n_total = total,
          n_missing = as.integer(miss_counts),
          missing_rate = miss_rate,
          n_cases = n_cases,
          n_controls = n_ctrls,
          missing_rate_cases = rate_cases,
          missing_rate_controls = rate_ctrls,
          stringsAsFactors = FALSE
        )
        per_path <- paste0("Missingness_PerVariable_", s$type, "_", s$time, "_", model_label, ".csv")
        utils::write.csv(per_var, per_path, row.names = FALSE)
        cat(sprintf("[摘要-缺失值] 已输出 %s (%d 变量)\n", per_path, nrow(per_var)))
        overview_rows[[length(overview_rows) + 1]] <- data.frame(
          model_comparison = model_label,
          variable_type = s$type,
          timepoint = s$time,
          var_count = length(vars),
          mean_missing_rate = mean(per_var$missing_rate, na.rm = TRUE),
          median_missing_rate = stats::median(per_var$missing_rate, na.rm = TRUE),
          max_missing_rate = max(per_var$missing_rate, na.rm = TRUE),
          top_missing_var = per_var$variable_name[which.max(per_var$missing_rate)],
          stringsAsFactors = FALSE
        )
      }
      if (length(overview_rows) > 0) {
        overview <- do.call(rbind, overview_rows)
        out_path <- paste0("Missingness_Overview_", model_label, ".csv")
        utils::write.csv(overview, out_path, row.names = FALSE)
        cat(sprintf("[摘要-缺失值] 已输出 %s\n", out_path))
      }
    }
    
    # =================== 第六.5c步：横纵一致性（空间重叠，按IDP根名称） ===================
    
    # 读取指定模型的四模型纵向结果（独立效应）
    read_four_model_independent <- function(model_label) {
      fpath <- file.path(getwd(), paste0("Four_Model_Z_Analysis_", model_label, ".csv"))
      if (!file.exists(fpath)) { cat(sprintf("[重叠] 未找到四模型结果：%s\n", fpath)); return(NULL) }
      df <- tryCatch(utils::read.csv(fpath, stringsAsFactors = FALSE), error = function(e) NULL)
      if (is.null(df) || !nrow(df)) return(NULL)
      # 仅脑影像，且独立效应（并对协变量独立）
      if (!("variable_type" %in% names(df))) df$variable_type <- NA_character_
      if (!("independent_effect" %in% names(df))) df$independent_effect <- FALSE
      if (!("independent_of_confounds" %in% names(df))) df$independent_of_confounds <- FALSE
      df <- df[df$variable_type == "Brain_Imaging" & df$independent_effect == TRUE, , drop = FALSE]
      # 若存在 independent_of_confounds，则进一步过滤
      df <- df[df$independent_of_confounds == TRUE | !("independent_of_confounds" %in% names(df)), , drop = FALSE]
      unique(df$variable_name)
    }
    
    # 计算横断Age-modulated显著与纵向独立效应的空间重叠
    compute_cross_long_overlap_for_model <- function(model_label) {
      # 读取横断 Age-modulated 的脑影像两时间点结果
      files <- c(
        paste0("CrossSection_GLM_Modulated_Brain_Imaging_Age1_", model_label, ".csv"),
        paste0("CrossSection_GLM_Modulated_Brain_Imaging_Age2_", model_label, ".csv")
      )
      cs_list <- list()
      for (fp in files) {
        if (!file.exists(fp)) { cat(sprintf("[重叠] 横断文件缺失：%s\n", fp)); next }
        df <- tryCatch(utils::read.csv(fp, stringsAsFactors = FALSE), error = function(e) NULL)
        if (is.null(df) || !nrow(df)) next
        time_tag <- if (grepl("Age1", fp)) "Age1" else "Age2"
        # 识别显著（FDR或FWE）
        sig <- df[(df$significant_fdr == TRUE) | (df$significant_fwe == TRUE), , drop = FALSE]
        if (!nrow(sig)) next
        sig$variable_root <- gsub("\\.{3,}Instance\\.[0-9]+.*$", "", sig$variable_name)
        sig$timepoint <- time_tag
        cs_list[[length(cs_list) + 1]] <- sig
      }
      if (length(cs_list) == 0) { cat(sprintf("[重叠] 模型 %s 横断Age-modulated无显著结果\n", model_label)); return(invisible(NULL)) }
      cs_all <- do.call(rbind, cs_list)
      # 聚合横断根名称：统计出现的时间点与最小p值
      agg <- lapply(split(cs_all, cs_all$variable_root), function(df) {
        tp <- paste(sort(unique(df$timepoint)), collapse = "+")
        # 使用p_value的最小值（或p_fdr最小值）
        min_p <- suppressWarnings(min(df$p_value, na.rm = TRUE))
        min_p_fdr <- suppressWarnings(min(df$p_fdr, na.rm = TRUE))
        data.frame(variable_root = df$variable_root[1], cross_timepoints = tp, cross_min_p = min_p, cross_min_p_fdr = min_p_fdr, stringsAsFactors = FALSE)
      })
      cross_df <- do.call(rbind, agg)
      
      # 读取纵向独立效应集合
      long_roots <- read_four_model_independent(model_label)
      long_set <- unique(long_roots)
      
      # 生成重叠结果
      cross_df$longitudinal_support <- cross_df$variable_root %in% long_set
      out_path <- paste0("Overlap_CrossSection_Modulated_vs_FourModel_", model_label, ".csv")
      utils::write.csv(cross_df, out_path, row.names = FALSE)
      cat(sprintf("[重叠] 已输出 %s（%d 条目，重叠 %d）\n", out_path, nrow(cross_df), sum(cross_df$longitudinal_support)))
      
      # 摘要
      summary_df <- data.frame(
        model_comparison = model_label,
        cross_significant_roots = nrow(cross_df),
        overlap_with_longitudinal = sum(cross_df$longitudinal_support),
        proportion_overlap = ifelse(nrow(cross_df) > 0, sum(cross_df$longitudinal_support) / nrow(cross_df), NA_real_),
        stringsAsFactors = FALSE
      )
      sum_path <- paste0("Overlap_Summary_", model_label, ".csv")
      utils::write.csv(summary_df, sum_path, row.names = FALSE)
      cat(sprintf("[重叠] 已输出 %s\n", sum_path))
      invisible(list(detail = cross_df, summary = summary_df))
    }
    
    # 对四个比较模型执行空间重叠
    compute_cross_long_overlap_for_model("ischemic_vs_control")
    compute_cross_long_overlap_for_model("mi_vs_control")
    compute_cross_long_overlap_for_model("chronic_vs_control")
    
    
    
    
    
compute_idp_change_rate_summary_A <- function(model_label = "ischemic_vs_control", independent_only = FALSE) {
      comb_csv <- combined_four_model_results_csv()
      if (!file.exists(comb_csv)) { cat("[IDP变化率-A] 未找到 Combined_Four_Model_Z_Analysis_Results.csv，跳过\n"); return(invisible(NULL)) }
      comb_df <- tryCatch(utils::read.csv(comb_csv, stringsAsFactors = FALSE), error = function(e) NULL)
      if (is.null(comb_df) || nrow(comb_df) == 0) { cat("[IDP变化率-A] Combined 结果为空，跳过\n"); return(invisible(NULL)) }
      idps <- comb_df[comb_df$model_comparison == model_label & comb_df$variable_type %in% c("Brain_Imaging", "Cognitive"), , drop = FALSE]
      if (independent_only && ("independent_effect" %in% names(idps))) {
        idps <- idps[isTRUE(idps$independent_effect) | idps$independent_effect == TRUE, , drop = FALSE]
      }
      if (nrow(idps) == 0) { cat(sprintf("[IDP变化率-A] 模型 %s 无可用条目\n", model_label)); return(invisible(NULL)) }

      rds_path <- resolve_existing_file(c(
        matched_cohort_rds_path(model_label),
        paste0("Final_Matched_Cohort_", model_label, ".rds")
      ))
      if (is.na(rds_path) || !nzchar(rds_path)) { cat(sprintf("[IDP变化率-A] 未找到匹配队列: %s\n", model_label)); return(invisible(NULL)) }
      dat <- readRDS(rds_path)
      if (!("group" %in% names(dat))) { cat(sprintf("[IDP变化率-A] 匹配队列缺少 'group': %s\n", rds_path)); return(invisible(NULL)) }

      get_num <- function(df, nm) {
        if (nm %in% names(df)) suppressWarnings(as.numeric(df[[nm]])) else rep(NA_real_, nrow(df))
      }
      get_chr <- function(df, nm) {
        if (nm %in% names(df)) as.character(df[[nm]]) else rep(NA_character_, nrow(df))
      }

      if (independent_only && ("variable_name" %in% names(idps))) {
        z_key <- abs(dplyr::coalesce(
          get_num(idps, "z_statistic_without_non_imaging"),
          get_num(idps, "z_statistic_with_non_imaging")
        ))
        p_key <- dplyr::coalesce(
          get_num(idps, "p_value_without_non_imaging"),
          get_num(idps, "p_value")
        )
        idps <- idps %>%
          dplyr::mutate(.z_key = z_key, .p_key = p_key) %>%
          dplyr::arrange(dplyr::desc(.z_key), .p_key) %>%
          dplyr::group_by(.data$variable_name) %>%
          dplyr::slice(1) %>%
          dplyr::ungroup() %>%
          dplyr::select(-.z_key, -.p_key)
      }

      z_stat <- dplyr::coalesce(
        get_num(idps, "z_statistic_without_non_imaging"),
        get_num(idps, "z_statistic_with_non_imaging")
      )
      p_val <- dplyr::coalesce(
        get_num(idps, "p_value_without_non_imaging"),
        get_num(idps, "p_value")
      )
      cd <- dplyr::coalesce(
        get_num(idps, "cohens_d"),
        get_num(idps, "cohens_d_without_non_imaging"),
        get_num(idps, "cohens_d_with_non_imaging")
      )

      clean_key <- function(x) sub("\\.*$", "", as.character(x))
      idps$idp_key <- clean_key(idps$variable_name)

      pairs_df <- NULL
      if (all(c("idp1_variable","idp2_variable") %in% names(idps))) {
        pairs_df <- unique(idps[, c("idp1_variable","idp2_variable"), drop = FALSE])
      }

      dl <- NULL
      if (!is.null(pairs_df) && nrow(pairs_df) > 0) {
        dl <- compute_delta_rate_long_parallel(dat, pairs_df = pairs_df)
      } else {
        dl <- compute_delta_rate_long_parallel(dat, target_roots = unique(idps$variable_name))
      }
      if (is.null(dl) || nrow(dl) == 0) { cat(sprintf("[IDP变化率-A] 未能计算Δ率: %s\n", model_label)); return(invisible(NULL)) }

      if ("eid" %in% names(dat) && ("eid" %in% names(dl))) {
        grp_df <- dat[, c("eid","group"), drop = FALSE]
        dl <- dplyr::left_join(dl, grp_df, by = "eid")
      } else {
        cat(sprintf("[IDP变化率-A] 无法按 eid 合并 group，跳过: %s\n", model_label))
        return(invisible(NULL))
      }

      dl$idp_key <- clean_key(dl$idp_root)
      dl$raw_change <- suppressWarnings(dl$idp2 - dl$idp1)
      dl$Group <- factorize_group(dl$group, labels = c("Control","Case"))

      sum_df <- dl %>%
        dplyr::filter(!is.na(idp_key), !is.na(Group)) %>%
        dplyr::group_by(idp_key, Group) %>%
        dplyr::summarise(
          n = {
            idx <- is.finite(raw_change) & is.finite(idp1)
            sum(idx, na.rm = TRUE)
          },
          idp1_mean = {
            idx <- is.finite(raw_change) & is.finite(idp1)
            if (sum(idx, na.rm = TRUE) > 0) base::mean(idp1[idx], na.rm = TRUE) else NA_real_
          },
          change_mean = {
            idx <- is.finite(raw_change) & is.finite(idp1)
            if (sum(idx, na.rm = TRUE) > 0) base::mean(raw_change[idx], na.rm = TRUE) else NA_real_
          },
          percent_change_mean = dplyr::if_else(
            is.finite(idp1_mean) & abs(idp1_mean) > .Machine$double.eps,
            (change_mean / abs(idp1_mean)) * 100,
            as.numeric(NA)
          ),
          .groups = "drop"
        )

      ctrl <- sum_df %>%
        dplyr::filter(Group == "Control") %>%
        dplyr::select(idp_key,
                      controls_n = n,
                      controls_change_mean = change_mean,
                      controls_percent_change_mean = percent_change_mean)

      case <- sum_df %>%
        dplyr::filter(Group == "Case") %>%
        dplyr::select(idp_key,
                      cases_n = n,
                      cases_change_mean = change_mean,
                      cases_percent_change_mean = percent_change_mean)

      merged_rates <- dplyr::full_join(ctrl, case, by = "idp_key") %>%
        dplyr::mutate(delta_percent_change = cases_percent_change_mean - controls_percent_change_mean)

      out <- data.frame(
        model_comparison = get_chr(idps, "model_comparison"),
        variable_name = get_chr(idps, "variable_name"),
        variable_type = get_chr(idps, "variable_type"),
        z_statistic = z_stat,
        p_value = p_val,
        cohens_d = cd,
        stringsAsFactors = FALSE
      )
      out$idp_key <- idps$idp_key

      out <- dplyr::left_join(out, merged_rates, by = "idp_key")
      out$idp_key <- NULL
      out
    }

    summarize_significant_idp_change_rates <- function(model_label = "ischemic_vs_control") {
      out <- compute_idp_change_rate_summary_A(model_label = model_label, independent_only = TRUE)
      if (is.null(out) || nrow(out) == 0) return(invisible(NULL))
      out_path <- table_out_path(file.path("Longitudinal", "Summary"), paste0("Significant_IDP_Change_Rate_Summary_", model_label, ".csv"))
      utils::write.csv(out, out_path, row.names = FALSE)
      cat(sprintf("[摘要-IDP变化率] 已输出 %s (%d 个IDP)\n", out_path, nrow(out)))
      invisible(out)
    }

    summarize_all_idp_change_rates_mean <- function(model_label = "ischemic_vs_control") {
      out <- compute_idp_change_rate_summary_A(model_label = model_label, independent_only = FALSE)
      if (is.null(out) || nrow(out) == 0) return(invisible(NULL))
      out_path <- table_out_path(file.path("Longitudinal", "Summary"), paste0("All_IDP_Change_Rate_Summary_", model_label, ".csv"))
      utils::write.csv(out, out_path, row.names = FALSE)
      cat(sprintf("[全量-IDP均值变化率] 已输出 %s (%d 个IDP)\n", out_path, nrow(out)))
      invisible(out)
    }
    
    summarize_delta_rate_median_vs_control_all_idps <- function(model_label = "ischemic_vs_control") {
      rds_path <- resolve_existing_file(c(
        matched_cohort_rds_path(model_label),
        paste0("Final_Matched_Cohort_", model_label, ".rds")
      ))
      if (is.na(rds_path) || !nzchar(rds_path)) { cat(sprintf("[中位数-Δ率-全量] 未找到匹配队列: %s\n", model_label)); return(invisible(NULL)) }
      dat <- readRDS(rds_path)
      if (!("group" %in% names(dat))) {
        warning(sprintf("匹配队列缺少 'group' 列，跳过ΔIDP分析: %s", rds_path))
        return(invisible(NULL))
      }

      pairs_df <- NULL
      meta_df <- NULL
      comb_csv <- combined_four_model_results_csv()
      if (file.exists(comb_csv)) {
        comb_df <- tryCatch(utils::read.csv(comb_csv, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
        if (!is.null(comb_df) && nrow(comb_df)) {
          sub_df <- comb_df[comb_df$model_comparison == model_label & comb_df$variable_type %in% c("Brain_Imaging", "Cognitive"), , drop = FALSE]
          if (nrow(sub_df) && all(c("idp1_variable","idp2_variable") %in% names(sub_df))) {
            pairs_df <- unique(sub_df[, c("idp1_variable","idp2_variable"), drop = FALSE])
            meta_cols <- intersect(c("idp1_variable","variable_name","variable_type","idp_category"), names(sub_df))
            meta_df <- unique(sub_df[, meta_cols, drop = FALSE])
          }
        }
      }

      dl <- compute_delta_rate_long_parallel(dat, pairs_df = pairs_df)
      if (is.null(dl) || nrow(dl) == 0) { cat("[中位数-Δ率-全量] 未能计算Δ率，跳过\n"); return(invisible(NULL)) }

      if ("eid" %in% names(dat) && ("eid" %in% names(dl))) {
        grp_df <- dat[, c("eid","group"), drop = FALSE]
        dl <- dplyr::left_join(dl, grp_df, by = "eid")
      } else {
        dl$group <- rep(dat$group, times = length(unique(dl$idp_root)))
      }

      dl$Group <- factorize_group(dl$group, labels = c("Control","Case"))

      median_df <- dl %>%
        dplyr::filter(is.finite(delta_rate), !is.na(idp_root), !is.na(Group)) %>%
        dplyr::mutate(delta_percent = delta_rate * 100) %>%
        dplyr::group_by(idp_root, Group) %>%
        dplyr::summarise(
          n = dplyr::n(),
          median_delta_percent = stats::median(delta_percent, na.rm = TRUE),
          mean_delta_percent = base::mean(delta_percent, na.rm = TRUE),
          q25 = stats::quantile(delta_percent, 0.25, na.rm = TRUE, type = 7),
          q75 = stats::quantile(delta_percent, 0.75, na.rm = TRUE, type = 7),
          .groups = "drop"
        )

      ctrl <- median_df %>%
        dplyr::filter(Group == "Control") %>%
        dplyr::select(idp_root,
                      control_n = n,
                      control_median = median_delta_percent,
                      control_mean = mean_delta_percent,
                      control_q25 = q25,
                      control_q75 = q75)
      case <- median_df %>%
        dplyr::filter(Group == "Case") %>%
        dplyr::select(idp_root,
                      case_n = n,
                      case_median = median_delta_percent,
                      case_mean = mean_delta_percent,
                      case_q25 = q25,
                      case_q75 = q75)

      merged <- dplyr::full_join(ctrl, case, by = "idp_root") %>%
        dplyr::mutate(
          median_diff_case_minus_control = case_median - control_median,
          mean_diff_case_minus_control = case_mean - control_mean,
          model_comparison = model_label
        )

      if (!is.null(meta_df) && nrow(meta_df)) {
        if ("idp1_variable" %in% names(meta_df)) {
          root_from_name <- function(x) gsub("(\\.{1,}|_)Instance\\.?[23]([._].*)?$", "", x)
          meta_df$idp_root <- root_from_name(meta_df$idp1_variable)
        } else if ("variable_name" %in% names(meta_df)) {
          meta_df$idp_root <- as.character(meta_df$variable_name)
        }
        if ("idp_root" %in% names(meta_df)) {
          meta_df <- meta_df %>%
            dplyr::filter(!is.na(.data$idp_root) & nzchar(as.character(.data$idp_root))) %>%
            dplyr::group_by(.data$idp_root) %>%
            dplyr::slice(1) %>%
            dplyr::ungroup()
          meta_keep <- intersect(c("idp_root","variable_type","idp_category"), names(meta_df))
          merged <- dplyr::left_join(merged, meta_df[, meta_keep, drop = FALSE], by = "idp_root")
        }
      }
      if (!("idp_category" %in% names(merged))) merged$idp_category <- NA_character_
      merged$idp_category <- ifelse(is.na(merged$idp_category), identify_idp_category_en(merged$idp_root), merged$idp_category)
      if (!("variable_type" %in% names(merged))) merged$variable_type <- NA_character_

      out_csv <- table_out_path(file.path("Longitudinal", "Median_DeltaRate"), paste0("Longitudinal_Median_DeltaRate_AllIDPs_Comparison_", model_label, ".csv"))
      utils::write.csv(merged, out_csv, row.names = FALSE)
      cat(sprintf("[中位数-Δ率-全量] 已输出 %s (IDP=%d)\n", out_csv, nrow(merged)))
      invisible(merged)
    }
    
    # 统一执行纵向摘要导出（缺失值、显著IDP变化率、Δ率中位数对比）
    safe_write_csv <- function(df, path) {
      try(utils::write.csv(df, path, row.names = FALSE), silent = TRUE)
    }

    combine_change_rate_summaries <- function(out_path = table_out_path(file.path("Longitudinal", "Summary"), "Longitudinal_Change_Rate_Summary_All_Models.csv")) {
      files <- list.files(get_table_dir(file.path("Longitudinal", "Summary")), pattern = "^Significant_IDP_Change_Rate_Summary_.*\\.csv$", full.names = TRUE)
      if (length(files) == 0) return(invisible(NULL))
      dfs <- lapply(files, function(fp) {
        tryCatch(utils::read.csv(fp, stringsAsFactors = FALSE), error = function(e) NULL)
      })
      dfs <- lapply(dfs, function(x) {
        if (is.null(x) || nrow(x) == 0) return(x)
        if ("model_comparison" %in% names(x)) {
          x <- x[!is.na(x$model_comparison) & nzchar(as.character(x$model_comparison)), , drop = FALSE]
        }
        x
      })
      dfs <- Filter(function(x) !is.null(x) && nrow(x) > 0, dfs)
      if (length(dfs) == 0) return(invisible(NULL))
      all_df <- dplyr::bind_rows(dfs)
      safe_write_csv(all_df, out_path)
      cat(sprintf("[综合-变化率] 已合并输出: %s (n=%d)\n", out_path, nrow(all_df)))
      invisible(all_df)
    }

    combine_all_idp_change_rate_summaries <- function(out_path = table_out_path(file.path("Longitudinal", "Summary"), "Longitudinal_Change_Rate_Summary_All_Models_AllIDPs.csv")) {
      files <- list.files(get_table_dir(file.path("Longitudinal", "Summary")), pattern = "^All_IDP_Change_Rate_Summary_.*\\.csv$", full.names = TRUE)
      if (length(files) == 0) return(invisible(NULL))
      dfs <- lapply(files, function(fp) {
        tryCatch(utils::read.csv(fp, stringsAsFactors = FALSE), error = function(e) NULL)
      })
      dfs <- lapply(dfs, function(x) {
        if (is.null(x) || nrow(x) == 0) return(x)
        if ("model_comparison" %in% names(x)) {
          x <- x[!is.na(x$model_comparison) & nzchar(as.character(x$model_comparison)), , drop = FALSE]
        }
        x
      })
      dfs <- Filter(function(x) !is.null(x) && nrow(x) > 0, dfs)
      if (length(dfs) == 0) return(invisible(NULL))
      all_df <- dplyr::bind_rows(dfs)
      safe_write_csv(all_df, out_path)
      cat(sprintf("[综合-全量均值变化率] 已合并输出: %s (n=%d)\n", out_path, nrow(all_df)))
      invisible(all_df)
    }

    run_longitudinal_summary_all_models <- function(models = c("ischemic_vs_control","mi_vs_control","chronic_vs_control")) {
      # 缺失值摘要仅在主比较生成一次（按当前逻辑）
      if ("ischemic_vs_control" %in% models) {
        summarize_missing_values_for_model("ischemic_vs_control")
      }
      # Δ率中位数对比按模型批量生成
      for (m in models) {
        summarize_delta_rate_median_vs_control_all_idps(m)
        # 同步生成显著IDP变化率摘要（若有独立效应IDP）
        summarize_significant_idp_change_rates(m)
        summarize_all_idp_change_rates_mean(m)
      }
      # 合并所有模型的显著IDP变化率摘要
      combine_change_rate_summaries()
      combine_all_idp_change_rate_summaries()
    }

    run_longitudinal_summary_all_models(c("ischemic_vs_control","mi_vs_control","chronic_vs_control"))

    # =================== 敏感性分析森林图（按模型） ===================
    
    # 按模型生成敏感性森林图（每个模型单独选择Top IDP并绘图）
    create_forest_plots_for_model <- function(model_base_label, combined_results, independent_idps, top_n = 20, overlay_n = 10, include_nonimaging = FALSE, plot_set = NA_character_) {
      if (is.null(combined_results) || nrow(combined_results) == 0) {
        return(NULL)
      }
      model_full <- model_base_label
      model_simpl <- paste0(model_base_label, "_no_sex_ethnicity")
      sub_df <- combined_results %>% dplyr::filter(model_comparison %in% c(model_full, model_simpl))
      if (nrow(sub_df) == 0) return(NULL)
      
      out_dir <- file.path("Forest_Plots", model_base_label)
      if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
      
      # 选择该模型的Top IDP名称：按模型内的独立效应与|Z|降序选TopN
      sub_full <- sub_df %>% dplyr::filter(model_comparison == model_full)
      if ("variable_type" %in% names(sub_full)) {
        sub_full <- sub_full %>% dplyr::filter(variable_type == "Brain_Imaging")
      }
      rank_df <- sub_full %>%
        dplyr::filter(independent_effect == TRUE) %>%
        dplyr::mutate(sort_key = dplyr::coalesce(abs(z_statistic_with_non_imaging), abs(z_statistic_without_non_imaging))) %>%
        dplyr::arrange(dplyr::desc(sort_key))
      if (nrow(rank_df) == 0) {
        # 若该模型无独立效应，则回退使用|beta|进行排序
        rank_df <- sub_full %>% dplyr::mutate(sort_key = abs(beta_estimate)) %>% dplyr::arrange(dplyr::desc(sort_key))
      }
      top_names <- head(unique(rank_df$variable_name), min(top_n, nrow(rank_df)))
      
      # 全模型森林图（按该模型分面）
      full_df <- sub_df %>%
        dplyr::filter(!grepl("_no_sex_ethnicity$", model_comparison), variable_name %in% top_names) %>%
        dplyr::select(model_comparison, variable_name, variable_type, beta_estimate, standard_error, p_formatted) %>%
        dplyr::mutate(lower = beta_estimate - 1.96 * standard_error, upper = beta_estimate + 1.96 * standard_error) %>%
        dplyr::filter(is.finite(beta_estimate), is.finite(standard_error))
      
      if (nrow(full_df) > 0) {
        p_full <- ggplot(full_df, aes(x = beta_estimate, y = reorder(variable_name, abs(beta_estimate)), color = variable_type)) +
          geom_point(size = 2) +
          geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
          facet_wrap(~ model_comparison, scales = "free_x") +
          geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
          scale_color_manual(values = c("Cognitive" = "#1F78B4", "Brain_Imaging" = "#E31A1C")) +
          labs(title = paste0("Forest Plot (Full) - ", model_base_label), x = "Beta Estimate (±95% CI)", y = "IDP (Top)") +
          theme_minimal() +
          theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), legend.position = "bottom", strip.text = element_text(size = 10, face = "bold"))
        ggsave(file.path(out_dir, paste0("Forest_Full_Top", top_n, "_", model_base_label, ".png")), p_full, width = 12, height = 9, dpi = 300)
        ggsave(file.path(out_dir, paste0("Forest_Full_Top", top_n, "_", model_base_label, ".pdf")), p_full, width = 12, height = 9)
      }
      
      # 性别/种族敏感性（全模型 vs 精简）
      sens_names <- head(top_names, min(overlay_n, length(top_names)))
      sens_df <- sub_df %>%
        dplyr::filter(variable_name %in% sens_names) %>%
        dplyr::mutate(analysis = ifelse(grepl("_no_sex_ethnicity$", model_comparison), "No Sex/Ethnicity", "Full"), model_base = sub("_no_sex_ethnicity$", "", model_comparison)) %>%
        dplyr::select(model_base, analysis, variable_name, variable_type, beta_estimate, standard_error) %>%
        dplyr::mutate(lower = beta_estimate - 1.96 * standard_error, upper = beta_estimate + 1.96 * standard_error) %>%
        dplyr::filter(is.finite(beta_estimate), is.finite(standard_error))
      
      if (nrow(sens_df) > 0) {
        dodge <- position_dodge(width = 0.5)
        p_sens <- ggplot(sens_df, aes(x = beta_estimate, y = reorder(variable_name, abs(beta_estimate)), color = analysis)) +
          geom_point(size = 2, position = dodge) +
          geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, position = dodge) +
          facet_wrap(~ model_base, scales = "free_x") +
          geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
          scale_color_manual(values = c("Full" = "#4DAF4A", "No Sex/Ethnicity" = "#984EA3")) +
          labs(title = paste0("Forest Plot (Sensitivity: Sex/Ethnicity) - ", model_base_label), x = "Beta Estimate (±95% CI)", y = "IDP (Top)") +
          theme_minimal() +
          theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), legend.position = "bottom", strip.text = element_text(size = 10, face = "bold"))
        ggsave(file.path(out_dir, paste0("Forest_Sensitivity_Top", min(overlay_n, length(top_names)), "_", model_base_label, ".png")), p_sens, width = 12, height = 9, dpi = 300)
        ggsave(file.path(out_dir, paste0("Forest_Sensitivity_Top", min(overlay_n, length(top_names)), "_", model_base_label, ".pdf")), p_sens, width = 12, height = 9)

        # 额外输出：固定 Top20 的敏感性森林图
        sens_names20 <- head(top_names, min(20, length(top_names)))
        sens_df20 <- sub_df %>%
          dplyr::filter(variable_name %in% sens_names20) %>%
          dplyr::mutate(analysis = ifelse(grepl("_no_sex_ethnicity$", model_comparison), "No Sex/Ethnicity", "Full"), model_base = sub("_no_sex_ethnicity$", "", model_comparison)) %>%
          dplyr::select(model_base, analysis, variable_name, variable_type, beta_estimate, standard_error) %>%
          dplyr::mutate(lower = beta_estimate - 1.96 * standard_error, upper = beta_estimate + 1.96 * standard_error) %>%
          dplyr::filter(is.finite(beta_estimate), is.finite(standard_error))
        if (nrow(sens_df20) > 0) {
          dodge20 <- position_dodge(width = 0.5)
          p_sens20 <- ggplot(sens_df20, aes(x = beta_estimate, y = reorder(variable_name, abs(beta_estimate)), color = analysis)) +
            geom_point(size = 2, position = dodge20) +
            geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, position = dodge20) +
            facet_wrap(~ model_base, scales = "free_x") +
            geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
            scale_color_manual(values = c("Full" = "#4DAF4A", "No Sex/Ethnicity" = "#984EA3")) +
            labs(title = paste0("Forest Plot (Sensitivity: Sex/Ethnicity, Top20) - ", model_base_label), x = "Beta Estimate (±95% CI)", y = "IDP (Top)") +
            theme_minimal() +
            theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), legend.position = "bottom", strip.text = element_text(size = 10, face = "bold"))
          ggsave(file.path(out_dir, paste0("Forest_Sensitivity_Top20_", model_base_label, ".png")), p_sens20, width = 12, height = 9, dpi = 300)
          ggsave(file.path(out_dir, paste0("Forest_Sensitivity_Top20_", model_base_label, ".pdf")), p_sens20, width = 12, height = 9)
        }
        
        # 稳健性数值汇总（Full vs No Sex/Ethnicity）
        sens_full <- sens_df %>% dplyr::filter(analysis == "Full") %>% dplyr::select(variable_name, beta_full = beta_estimate, se_full = standard_error)
        sens_no   <- sens_df %>% dplyr::filter(analysis == "No Sex/Ethnicity") %>% dplyr::select(variable_name, beta_no = beta_estimate, se_no = standard_error)
        sens_comp <- dplyr::inner_join(sens_full, sens_no, by = "variable_name")
        
        # p值一致性：优先使用未校正p值，如不可用则尝试FDR或通用p_value
        full_p <- sub_df %>% dplyr::filter(model_comparison == model_full, variable_name %in% sens_names) %>% dplyr::select(variable_name, p_unadj = p_value_unadjusted, p_fdr = p_value_fdr)
        no_p   <- sub_df %>% dplyr::filter(model_comparison == model_simpl, variable_name %in% sens_names) %>% dplyr::select(variable_name, p_unadj_no = p_value_unadjusted, p_fdr_no = p_value_fdr)
        p_comp <- dplyr::inner_join(full_p, no_p, by = "variable_name")
        prop_sig_unadj <- if (all(c("p_unadj","p_unadj_no") %in% names(p_comp))) mean((p_comp$p_unadj < 0.05) == (p_comp$p_unadj_no < 0.05), na.rm = TRUE) else NA
        prop_sig_fdr   <- if (all(c("p_fdr","p_fdr_no") %in% names(p_comp))) mean((p_comp$p_fdr < 0.05) == (p_comp$p_fdr_no < 0.05), na.rm = TRUE) else NA
        p_any          <- dplyr::coalesce(p_comp$p_unadj, p_comp$p_fdr)
        p_any_no       <- dplyr::coalesce(p_comp$p_unadj_no, p_comp$p_fdr_no)
        prop_sig_any   <- mean((p_any < 0.05) == (p_any_no < 0.05), na.rm = TRUE)
        
        robustness_sexethnicity <- data.frame(
          model_label = model_base_label,
          n_idps_compared = nrow(sens_comp),
          corr_beta_full_vs_no = {
            idx <- is.finite(sens_comp$beta_full) & is.finite(sens_comp$beta_no)
            if (sum(idx, na.rm = TRUE) >= 2) suppressWarnings(stats::cor(sens_comp$beta_full[idx], sens_comp$beta_no[idx], use = "complete.obs")) else NA_real_
          },
          prop_same_sign = mean(sign(sens_comp$beta_full) == sign(sens_comp$beta_no), na.rm = TRUE),
          mean_abs_beta_diff = mean(abs(sens_comp$beta_full - sens_comp$beta_no), na.rm = TRUE),
          median_abs_beta_diff = stats::median(abs(sens_comp$beta_full - sens_comp$beta_no), na.rm = TRUE),
          mean_abs_se_diff = mean(abs(sens_comp$se_full - sens_comp$se_no), na.rm = TRUE),
          median_abs_se_diff = stats::median(abs(sens_comp$se_full - sens_comp$se_no), na.rm = TRUE),
          prop_significance_consistent_unadj = prop_sig_unadj,
          prop_significance_consistent_fdr = prop_sig_fdr,
          prop_significance_consistent_any = prop_sig_any,
          stringsAsFactors = FALSE
        )
        utils::write.csv(robustness_sexethnicity, file.path(out_dir, "SexEthnicity_Sensitivity_Comparison_Summary.csv"), row.names = FALSE)
      }
      
      # 非影像协变量敏感性（With vs Without Non-Imaging；仅对全模型）
      has_nonimg_cols <- all(c("beta_estimate_with_non_imaging", "standard_error_with_non_imaging", "p_value_with_non_imaging", "p_value_without_non_imaging") %in% names(sub_df))
      if (include_nonimaging && has_nonimg_cols) {
        base_df <- sub_df %>%
          dplyr::filter(!grepl("_no_sex_ethnicity$", model_comparison), variable_name %in% sens_names) %>%
          dplyr::select(model_comparison, variable_name, variable_type, beta_estimate, standard_error, beta_estimate_with_non_imaging, standard_error_with_non_imaging, p_value_without_non_imaging, p_value_with_non_imaging) %>%
          dplyr::filter(is.finite(beta_estimate) | is.finite(beta_estimate_with_non_imaging))
        
        without_df <- base_df %>% dplyr::transmute(model_comparison, variable_name, variable_type, analysis = "Without Non-Imaging", beta = beta_estimate, se = standard_error, lower = beta_estimate - 1.96 * standard_error, upper = beta_estimate + 1.96 * standard_error)
        with_df    <- base_df %>% dplyr::transmute(model_comparison, variable_name, variable_type, analysis = "With Non-Imaging",    beta = beta_estimate_with_non_imaging, se = standard_error_with_non_imaging, lower = beta_estimate_with_non_imaging - 1.96 * standard_error_with_non_imaging, upper = beta_estimate_with_non_imaging + 1.96 * standard_error_with_non_imaging)
        
        nonimg_df <- dplyr::bind_rows(without_df, with_df) %>% dplyr::filter(is.finite(beta), is.finite(se))
        if (nrow(nonimg_df) > 0) {
          dodge2 <- position_dodge(width = 0.55)
          p_nonimg <- ggplot(nonimg_df, aes(x = beta, y = reorder(variable_name, abs(beta)), color = analysis)) +
            geom_point(size = 2, position = dodge2) +
            geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, position = dodge2) +
            facet_wrap(~ model_comparison, scales = "free_x") +
            geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
            scale_color_manual(values = c("Without Non-Imaging" = "#FF7F00", "With Non-Imaging" = "#4DAF4A")) +
            labs(title = paste0("Forest Plot (Non-Imaging Sensitivity) - ", model_base_label), x = "Beta Estimate (±95% CI)", y = "IDP (Top)") +
            theme_minimal() +
            theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), legend.position = "bottom", strip.text = element_text(size = 10, face = "bold"))
          ggsave(file.path(out_dir, paste0("Forest_NonImaging_Top", min(overlay_n, length(top_names)), "_", model_base_label, ".png")), p_nonimg, width = 12, height = 9, dpi = 300)
          ggsave(file.path(out_dir, paste0("Forest_NonImaging_Top", min(overlay_n, length(top_names)), "_", model_base_label, ".pdf")), p_nonimg, width = 12, height = 9)
          
          # 稳健性数值汇总（相关、同号率、平均差、显著性一致性）
          robustness_nonimg <- base_df %>%
            dplyr::group_by(model_comparison) %>%
            dplyr::summarise(
              n_idps_compared = dplyr::n(),
              corr_beta_with_vs_without = {
                idx <- is.finite(beta_estimate) & is.finite(beta_estimate_with_non_imaging)
                if (sum(idx, na.rm = TRUE) >= 2) {
                  suppressWarnings(stats::cor(beta_estimate[idx], beta_estimate_with_non_imaging[idx]))
                } else {
                  NA_real_
                }
              },
              prop_same_sign = mean(sign(beta_estimate) == sign(beta_estimate_with_non_imaging), na.rm = TRUE),
              mean_abs_beta_diff = mean(abs(beta_estimate - beta_estimate_with_non_imaging), na.rm = TRUE),
              median_abs_beta_diff = stats::median(abs(beta_estimate - beta_estimate_with_non_imaging), na.rm = TRUE),
              mean_abs_se_diff = mean(abs(standard_error - standard_error_with_non_imaging), na.rm = TRUE),
              median_abs_se_diff = stats::median(abs(standard_error - standard_error_with_non_imaging), na.rm = TRUE),
              prop_significance_consistent = mean((p_value_without_non_imaging < 0.05) == (p_value_with_non_imaging < 0.05), na.rm = TRUE),
              .groups = "drop"
            )
          utils::write.csv(robustness_nonimg, file.path(out_dir, "NonImaging_Sensitivity_Comparison_Summary.csv"), row.names = FALSE)
        }
      }
      
      # 逐变量简化模型：从全非成像协变量中仅移除指定组，输出专属森林图与稳健性CSV
      per_variable_dirs <- c()
      if (include_nonimaging && exists("analysis_data_list") && model_full %in% names(analysis_data_list)) {
        data_model <- analysis_data_list[[model_full]]
        # 确定候选非成像协变量组（存在于当前模型数据中）
        base_non_img <- c("smoking_factor_baseline", "alcohol_factor_baseline",
                           "townsend_index_baseline", "bmi_baseline",
                           "waist_circumference_baseline",
                           "diabetes_factor_baseline", "systolic_bp_baseline", "diastolic_bp_baseline")
        # 按研究要求，移除饮食/运动等扩展变量，仅保留核心非影像协变量
        non_img_groups <- intersect(base_non_img, names(data_model))
        # 仅对Top IDP对应的配对进行评估，提升效率
        if (exists("variable_pairs")) {
          pairs_top <- Filter(function(p) p$base_name %in% sens_names, variable_pairs)
        } else {
          pairs_top <- list()
        }
        if (length(non_img_groups) > 0 && length(pairs_top) > 0) {
          for (g in non_img_groups) {
            dir_g <- file.path(out_dir, "NonImaging_PerVariable", g)
            if (!dir.exists(dir_g)) dir.create(dir_g, recursive = TRUE)
            per_variable_dirs <- unique(c(per_variable_dirs, file.path("NonImaging_PerVariable", g)))
            
            # 使用“全非成像协变量移除g”作为扩展模型；与“全非成像协变量”进行对比
            non_img_minus <- setdiff(non_img_groups, g)
            results_minus <- tryCatch(
              z_statistics_four_model_analysis(
                data_model,
                pairs_top,
                model_name = paste0(model_full, "_minus_", g),
                drop_demographics = FALSE,
                non_imaging_covariates = non_img_minus
              ), error = function(e) NULL
            )
            if (is.null(results_minus) || nrow(results_minus) == 0) next
            
            join_full <- sub_full %>%
              dplyr::filter(variable_name %in% sens_names) %>%
              dplyr::select(variable_name, variable_type, beta_full_ext = beta_estimate_with_non_imaging, se_full_ext = standard_error_with_non_imaging, p_full_ext = p_value_with_non_imaging)
            join_minus <- results_minus %>%
              dplyr::filter(variable_name %in% sens_names) %>%
              dplyr::select(variable_name, variable_type, beta_minus_ext = beta_estimate_with_non_imaging, se_minus_ext = standard_error_with_non_imaging, p_minus_ext = p_value_with_non_imaging)
            join_df <- dplyr::inner_join(join_full, join_minus, by = c("variable_name", "variable_type"))
            if (nrow(join_df) == 0) next
            
            # p值一致性（未校正与BH-FDR）
            prop_sig_unadj_ext <- mean((join_df$p_full_ext < 0.05) == (join_df$p_minus_ext < 0.05), na.rm = TRUE)
            fdr_full <- stats::p.adjust(join_df$p_full_ext, method = "BH")
            fdr_minus <- stats::p.adjust(join_df$p_minus_ext, method = "BH")
            prop_sig_fdr_ext <- mean((fdr_full < 0.05) == (fdr_minus < 0.05), na.rm = TRUE)
            
            # 稳健性汇总
            robustness_group <- data.frame(
              model_label = model_base_label,
              group_removed = g,
              n_idps_compared = nrow(join_df),
              corr_beta_full_vs_minus = {
                idx2 <- is.finite(join_df$beta_full_ext) & is.finite(join_df$beta_minus_ext)
                if (sum(idx2, na.rm = TRUE) >= 2) suppressWarnings(stats::cor(join_df$beta_full_ext[idx2], join_df$beta_minus_ext[idx2], use = "complete.obs")) else NA_real_
              },
              prop_same_sign = mean(sign(join_df$beta_full_ext) == sign(join_df$beta_minus_ext), na.rm = TRUE),
              mean_abs_beta_diff = mean(abs(join_df$beta_full_ext - join_df$beta_minus_ext), na.rm = TRUE),
              median_abs_beta_diff = stats::median(abs(join_df$beta_full_ext - join_df$beta_minus_ext), na.rm = TRUE),
              mean_abs_se_diff = mean(abs(join_df$se_full_ext - join_df$se_minus_ext), na.rm = TRUE),
              median_abs_se_diff = stats::median(abs(join_df$se_full_ext - join_df$se_minus_ext), na.rm = TRUE),
              prop_significance_consistent_unadj = prop_sig_unadj_ext,
              prop_significance_consistent_fdr = prop_sig_fdr_ext,
              stringsAsFactors = FALSE
            )
            utils::write.csv(robustness_group, file.path(dir_g, "NonImaging_PerVariable_Sensitivity_Summary.csv"), row.names = FALSE)
            
            # 对比明细表
            comp_df <- join_df %>% dplyr::arrange(dplyr::desc(abs(beta_full_ext)))
            utils::write.csv(comp_df, file.path(dir_g, "NonImaging_PerVariable_Comparison.csv"), row.names = FALSE)
            
            # 专属森林图（全扩展 vs 去除当前组）
            plot_full <- comp_df %>% dplyr::transmute(variable_name, variable_type, analysis = "Full Extended", beta = beta_full_ext, se = se_full_ext, lower = beta_full_ext - 1.96 * se_full_ext, upper = beta_full_ext + 1.96 * se_full_ext)
            plot_minus <- comp_df %>% dplyr::transmute(variable_name, variable_type, analysis = paste0("Minus ", g), beta = beta_minus_ext, se = se_minus_ext, lower = beta_minus_ext - 1.96 * se_minus_ext, upper = beta_minus_ext + 1.96 * se_minus_ext)
            plot_df <- dplyr::bind_rows(plot_full, plot_minus)
            dodge3 <- position_dodge(width = 0.55)
            p_group <- ggplot(plot_df, aes(x = beta, y = reorder(variable_name, abs(beta)), color = analysis)) +
              geom_point(size = 2, position = dodge3) +
              geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, position = dodge3) +
              facet_wrap(~ variable_type, scales = "free_x") +
              geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
              scale_color_manual(values = setNames(c("#4DAF4A", "#FF7F00"), c("Full Extended", paste0("Minus ", g)))) +
              labs(title = paste0("Forest Plot (Per-Variable Sensitivity: ", g, ") - ", model_base_label), x = "Beta Estimate (±95% CI)", y = "IDP (Top)") +
              theme_minimal() +
              theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), legend.position = "bottom", strip.text = element_text(size = 10, face = "bold"))
          ggsave(file.path(dir_g, paste0("Forest_NonImaging_PerVariable_", g, "_Top", min(overlay_n, length(top_names)), "_", model_base_label, ".png")), p_group, width = 12, height = 9, dpi = 300)
          ggsave(file.path(dir_g, paste0("Forest_NonImaging_PerVariable_", g, "_Top", min(overlay_n, length(top_names)), "_", model_base_label, ".pdf")), p_group, width = 12, height = 9)
          }
        }
      }

      # 新增：按IDP输出协变量影响森林图（Full Extended vs 去除各核心协变量）
      per_idp_dirs <- c()
      if (include_nonimaging && exists("analysis_data_list") && model_full %in% names(analysis_data_list)) {
        data_model <- analysis_data_list[[model_full]]
        base_non_img <- c("smoking_factor_baseline", "alcohol_factor_baseline",
                           "townsend_index_baseline", "bmi_baseline",
                           "waist_circumference_baseline",
                           "diabetes_factor_baseline", "systolic_bp_baseline", "diastolic_bp_baseline")
        non_img_groups <- intersect(base_non_img, names(data_model))

        # 针对每个Top IDP，比较Full Extended与逐一去除协变量后的估计
        if (length(non_img_groups) > 0) {
          per_idp_root_dir <- file.path(out_dir, "NonImaging_PerIDP")
          if (!dir.exists(per_idp_root_dir)) dir.create(per_idp_root_dir, recursive = TRUE)
          per_idp_dirs <- c(per_idp_dirs, "NonImaging_PerIDP")

          # 为快速定位，限定到Top IDP配对集合
          pairs_top_all <- if (exists("variable_pairs")) Filter(function(p) p$base_name %in% sens_names, variable_pairs) else list()
          # 预取Full Extended的IDP效果
          full_ext_df <- sub_full %>%
            dplyr::filter(variable_name %in% sens_names) %>%
            dplyr::select(variable_name, variable_type, beta_full_ext = beta_estimate_with_non_imaging, se_full_ext = standard_error_with_non_imaging, p_full_ext = p_value_with_non_imaging)

          for (idp_nm in sens_names) {
            # 收集该IDP在“去除单一协变量”下的估计
            minus_list <- list()
            for (g in non_img_groups) {
              non_img_minus <- setdiff(non_img_groups, g)
              res_minus <- tryCatch(
                z_statistics_four_model_analysis(
                  data_model,
                  Filter(function(p) p$base_name == idp_nm, pairs_top_all),
                  model_name = paste0(model_full, "_minus_", g),
                  drop_demographics = FALSE,
                  non_imaging_covariates = non_img_minus
                ), error = function(e) NULL
              )
              if (!is.null(res_minus) && nrow(res_minus) > 0) {
                row_i <- res_minus %>% dplyr::filter(variable_name == idp_nm) %>%
                  dplyr::select(variable_name, variable_type, beta_minus_ext = beta_estimate_with_non_imaging,
                                se_minus_ext = standard_error_with_non_imaging, p_minus_ext = p_value_with_non_imaging)
                if (nrow(row_i) > 0) {
                  row_i$analysis <- paste0("Minus ", g)
                  minus_list[[length(minus_list) + 1]] <- row_i
                }
              }
            }
            idp_minus_df <- if (length(minus_list) > 0) dplyr::bind_rows(minus_list) else NULL
            idp_full_row <- full_ext_df %>% dplyr::filter(variable_name == idp_nm)
            if (!is.null(idp_minus_df) && nrow(idp_minus_df) > 0 && nrow(idp_full_row) > 0) {
              # 组合绘图数据
              plot_full <- idp_full_row %>% dplyr::transmute(variable_name, variable_type, analysis = "Full Extended",
                                                             beta = beta_full_ext, se = se_full_ext,
                                                             lower = beta_full_ext - 1.96 * se_full_ext,
                                                             upper = beta_full_ext + 1.96 * se_full_ext)
              plot_minus <- idp_minus_df %>% dplyr::transmute(variable_name, variable_type, analysis,
                                                               beta = beta_minus_ext, se = se_minus_ext,
                                                               lower = beta_minus_ext - 1.96 * se_minus_ext,
                                                               upper = beta_minus_ext + 1.96 * se_minus_ext)
              plot_df <- dplyr::bind_rows(plot_full, plot_minus) %>% dplyr::filter(is.finite(beta), is.finite(se))
              if (nrow(plot_df) > 1) {
                dodge4 <- position_dodge(width = 0.55)
                p_idp <- ggplot(plot_df, aes(x = beta, y = reorder(analysis, abs(beta)), color = analysis)) +
                  geom_point(size = 2, position = dodge4) +
                  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, position = dodge4) +
                  facet_wrap(~ variable_type, scales = "free_x") +
                  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
                  scale_color_manual(values = c("Full Extended" = "#4DAF4A")) +
                  labs(title = paste0("Forest (Per-IDP Non-Imaging): ", idp_nm, " - ", model_base_label),
                       x = "Beta Estimate (±95% CI)", y = "Analysis") +
                  theme_minimal() +
                  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
                        legend.position = "none", strip.text = element_text(size = 10, face = "bold"))
                # 安全的文件名
                safe_idp <- gsub("[^A-Za-z0-9_]+", "_", substr(idp_nm, 1, 60))
                ggsave(file.path(per_idp_root_dir, paste0("Forest_NonImaging_PerIDP_", safe_idp, "_", model_base_label, ".png")), p_idp, width = 10, height = 7, dpi = 300)
                ggsave(file.path(per_idp_root_dir, paste0("Forest_NonImaging_PerIDP_", safe_idp, "_", model_base_label, ".pdf")), p_idp, width = 10, height = 7)
              }
            }
          }
        }
      }
      
      data.frame(
        plot_set = plot_set,
        model_label = model_base_label,
        files_full = ifelse(nrow(full_df) > 0, paste0("Forest_Full_Top", top_n, "_", model_base_label, ".png + .pdf"), "None"),
        files_sensitivity = ifelse(nrow(sens_df) > 0, paste0("Forest_Sensitivity_Top", min(overlay_n, length(top_names)), "_", model_base_label, ".png + .pdf"), "None"),
        files_sensitivity_top20 = ifelse(exists("sens_df20") && !is.null(sens_df20) && nrow(sens_df20) > 0, paste0("Forest_Sensitivity_Top20_", model_base_label, ".png + .pdf"), "None"),
        files_nonimaging = ifelse(include_nonimaging && exists("nonimg_df") && !is.null(nonimg_df) && nrow(nonimg_df) > 0, paste0("Forest_NonImaging_Top", min(overlay_n, length(top_names)), "_", model_base_label, ".png + .pdf"), "None"),
        sexethnicity_summary = ifelse(nrow(sens_df) > 0, "SexEthnicity_Sensitivity_Comparison_Summary.csv", "None"),
        nonimaging_summary = ifelse(include_nonimaging && has_nonimg_cols && exists("nonimg_df") && !is.null(nonimg_df) && nrow(nonimg_df) > 0, "NonImaging_Sensitivity_Comparison_Summary.csv", "None"),
        per_variable_sensitivity_dir = ifelse(length(per_variable_dirs) > 0, paste(unique(per_variable_dirs), collapse = "; "), "None"),
        per_idp_sensitivity_dir = ifelse(length(per_idp_dirs) > 0, paste(unique(per_idp_dirs), collapse = "; "), "None"),
        output_dir = out_dir,
        stringsAsFactors = FALSE
      )
    }
    
    cat("\n=== 分模型敏感性森林图 ===\n")
    # 分模型森林图：逐模型生成独立的敏感性分析和图件
    base_models <- unique(gsub("_no_sex_ethnicity$", "", combined_four_model_results$model_comparison))
    per_model_reports <- lapply(base_models, function(mb) {
      n_indep <- sum(combined_four_model_results$model_comparison == mb & combined_four_model_results$independent_effect == TRUE, na.rm = TRUE)
      # 若该模型下没有“独立效应IDP”，回退到20以保证至少有输出
      n_indep <- ifelse(is.na(n_indep) || n_indep <= 0, 20, n_indep)
      main_rep <- create_forest_plots_for_model(
        mb,
        combined_four_model_results,
        independent_effect_idps,
        top_n = n_indep,
        overlay_n = n_indep,
        include_nonimaging = FALSE,
        plot_set = "TopN_Independent"
      )
      top10_rep <- create_forest_plots_for_model(
        mb,
        combined_four_model_results,
        independent_effect_idps,
        top_n = min(10, n_indep),
        overlay_n = min(10, n_indep),
        include_nonimaging = FALSE,
        plot_set = "Top10"
      )
      dplyr::bind_rows(main_rep, top10_rep)
    })
    non_null_reports <- per_model_reports[!sapply(per_model_reports, is.null)]
    per_model_reports <- if (length(non_null_reports) > 0) dplyr::bind_rows(non_null_reports) else NULL
    if (!is.null(per_model_reports) && nrow(per_model_reports) > 0) {
      utils::write.csv(per_model_reports, table_out_path("Forest_Plots", "Forest_PerModel_Summary.csv"), row.names = FALSE)
      cat("已输出分模型森林图汇总: ", table_out_path("Forest_Plots", "Forest_PerModel_Summary.csv"), "\n")
    }
  }
}

# =================== 第七步：心梗与慢性缺血性Z统计量相关性分析 ===================
cat("\n第七步：心梗与慢性缺血性Z统计量相关性分析...\n")

if("mi_vs_control" %in% names(all_four_model_results) && 
   "chronic_vs_control" %in% names(all_four_model_results)) {
  
  mi_results <- all_four_model_results[["mi_vs_control"]]
  chronic_results <- all_four_model_results[["chronic_vs_control"]]
  
  # 匹配相同的变量
  common_variables <- intersect(mi_results$variable_name, chronic_results$variable_name)
  
  if(length(common_variables) > 0) {
    mi_matched <- mi_results[mi_results$variable_name %in% common_variables, ]
    chronic_matched <- chronic_results[chronic_results$variable_name %in% common_variables, ]
    
    # 确保变量顺序一致
    mi_matched <- mi_matched[order(mi_matched$variable_name), ]
    chronic_matched <- chronic_matched[order(chronic_matched$variable_name), ]
    
    # 计算相关性
    z_correlation <- cor(mi_matched$z_statistic_without_non_imaging, 
                         chronic_matched$z_statistic_without_non_imaging,
                         use = "complete.obs")
    
    cor_test <- cor.test(mi_matched$z_statistic_without_non_imaging,
                         chronic_matched$z_statistic_without_non_imaging)
    
    cat("心梗组与慢性缺血性组Z统计量相关性:\n")
    cat("- 共同变量数:", length(common_variables), "个\n")
    cat("- 相关系数 r =", round(z_correlation, 4), "\n")
    cat("- 95%置信区间: [", round(cor_test$conf.int[1], 4), ", ", round(cor_test$conf.int[2], 4), "]\n")
    cat("- p值:", ifelse(cor_test$p.value < 0.001, "<0.001", sprintf("%.4f", cor_test$p.value)), "\n")
    
    # 保存相关性数据
    correlation_data <- data.frame(
      variable_name = mi_matched$variable_name,
      variable_type = mi_matched$variable_type,
      mi_z_statistic = mi_matched$z_statistic_without_non_imaging,
      chronic_z_statistic = chronic_matched$z_statistic_without_non_imaging,
      mi_independent = mi_matched$independent_effect,
      chronic_independent = chronic_matched$independent_effect,
      both_independent = mi_matched$independent_effect & chronic_matched$independent_effect,
      z_difference = abs(mi_matched$z_statistic_without_non_imaging - chronic_matched$z_statistic_without_non_imaging),
      stringsAsFactors = FALSE
    )
    
    write.csv(correlation_data, table_out_path(file.path("Four_Model", "Diagnostics"), "MI_vs_Chronic_Z_Correlation_Analysis.csv"), row.names = FALSE)
  }
}

# （新增）纵向ΔIDP变化率分析（按匹配队列）
cat("\n（新增）纵向ΔIDP变化率分析（按匹配队列）...\n")

if (!exists("comparisons")) {
  comparisons <- c("ischemic_vs_control", "mi_vs_control", "chronic_vs_control")
}
# 进度条：纵向ΔIDP变化率（各比较）
pb_delta <- txtProgressBar(min = 0, max = length(comparisons), style = 3)
idx_delta <- 0
for (cmp in comparisons) {
  rds_path <- resolve_existing_file(c(
    matched_cohort_rds_path(cmp),
    paste0("Final_Matched_Cohort_", cmp, ".rds")
  ))
  if (is.na(rds_path) || !nzchar(rds_path)) {
    warning(sprintf("未找到匹配队列文件，跳过ΔIDP分析: %s", cmp))
    idx_delta <- idx_delta + 1; setTxtProgressBar(pb_delta, idx_delta)
    next
  }
  dat <- tryCatch(readRDS(rds_path), error = function(e) NULL)
  if (is.null(dat)) {
    warning(sprintf("读取匹配队列失败，跳过ΔIDP分析: %s", rds_path))
    idx_delta <- idx_delta + 1; setTxtProgressBar(pb_delta, idx_delta)
    next
  }
  # 安全检查：必须存在分组列
  if (!("group" %in% names(dat))) {
    warning(sprintf("匹配队列缺少 'group' 列，跳过ΔIDP分析: %s", rds_path))
    idx_delta <- idx_delta + 1; setTxtProgressBar(pb_delta, idx_delta)
    next
  }
  # 基线年龄：优先 baseline_age，其次 age_at_instance2，再次 age_at_recruitment
  if (!("baseline_age" %in% names(dat))) {
    dat$baseline_age <- if ("age_at_instance2" %in% names(dat)) dat$age_at_instance2 else dat$age_at_recruitment
  }
  # 检查是否存在任一配对列（Instance.2 与 Instance.3），不再局限于脑影像IDP，允许认知变量参与
  has_inst2 <- any(grepl("\\.{3,}Instance\\.2(\\.|$)", names(dat)))
  has_inst3 <- any(grepl("\\.{3,}Instance\\.3(\\.|$)", names(dat)))
  if (!(has_inst2 && has_inst3)) {
    warning(sprintf("未发现 Instance.2/Instance.3 配对列，跳过ΔIDP分析: %s", rds_path))
    idx_delta <- idx_delta + 1; setTxtProgressBar(pb_delta, idx_delta)
    next
  }
  # 仅对独立效应的IDP根进行Δ率计算以加速（不做MICE插补，避免耗时与卡顿）
  target_roots <- NULL
  csv_path <- all_independent_effect_idps_csv()
  idps <- tryCatch(utils::read.csv(csv_path, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(idps) && nrow(idps)) {
      idps <- idps[idps$model_comparison == cmp, , drop = FALSE]
      if (nrow(idps)) {
        # 优先使用 idp1_variable（更稳定且包含 Instance.2），否则回退 variable_name
        root_src <- if ("idp1_variable" %in% names(idps)) idps$idp1_variable else idps$variable_name
        target_roots <- unique(gsub("\\.\\.\\.Instance\\.[0-9]+", "", root_src))
        target_roots <- target_roots[!is.na(target_roots) & nzchar(target_roots)]
      }
    }
  # 若该模型无“独立效应IDP”，与箱线图一致地回退到四模型综合结果的Top集合（包含认知与脑影像）
  if (is.null(target_roots) || length(target_roots) == 0) {
    comb_csv <- combined_four_model_results_csv()
    comb_df <- tryCatch(utils::read.csv(comb_csv, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(comb_df) && nrow(comb_df)) {
      sub_df <- comb_df[comb_df$model_comparison == cmp, , drop = FALSE]
      if (nrow(sub_df)) {
        sub_df$sort_key <- dplyr::coalesce(abs(sub_df$z_statistic_with_non_imaging), abs(sub_df$z_statistic_without_non_imaging))
        sub_df <- sub_df[order(-sub_df$sort_key), , drop = FALSE]
        fallback_n <- 24L
        top_df <- head(sub_df, min(fallback_n, nrow(sub_df)))
        root_src <- if ("idp1_variable" %in% names(top_df)) top_df$idp1_variable else top_df$variable_name
        target_roots <- unique(gsub("\\.\\.\\.Instance\\.[0-9]+", "", root_src))
        target_roots <- target_roots[!is.na(target_roots) & nzchar(target_roots)]
        if (length(target_roots) > 0) {
          message(sprintf("[Δ率] %s 使用回退Top%d集合（与箱线图一致）", cmp, length(target_roots)))
        }
      }
    }
  }
  dl <- compute_delta_rate_long_parallel(dat, target_roots = target_roots)
  if (is.null(dl) || nrow(dl) == 0) {
    idx_delta <- idx_delta + 1; setTxtProgressBar(pb_delta, idx_delta)
    next
  }
  cat(sprintf("[Δ率] %s 计算完成：IDP根=%d，记录数=%d\n", cmp, length(unique(dl$idp_root)), nrow(dl)))
  # 将协变量按 Participant.ID 左连接，确保行对齐
  if (("eid" %in% names(dat)) && ("eid" %in% names(dl))) {
    cov_df <- dat[, c("eid","baseline_age","group"), drop = FALSE]
    dl <- dplyr::left_join(dl, cov_df, by = "eid")
  } else {
    # 无 Participant.ID 时：按 idp_root 重复向量长度以匹配行数
    kroots <- length(unique(dl$idp_root))
    dl$baseline_age <- rep(dat$baseline_age, times = kroots)
    dl$group <- rep(dat$group, times = kroots)
  }
  dl$comparison <- cmp
  dl$idp_category <- identify_idp_category_en(dl$idp_root)
  delta_all <- dl
  if (is.null(delta_all) || nrow(delta_all) == 0) {
    idx_delta <- idx_delta + 1; setTxtProgressBar(pb_delta, idx_delta)
    next
  }

  # 年龄分布图（双面板）
  delta_all$delta_percent <- delta_all$delta_rate * 100
  # 统一分组变量命名（数值优先，字符串回退）
  Group_tmp <- factorize_group(delta_all$group, labels = c("Controls","Cases"))
  if (all(is.na(Group_tmp)) && is.character(delta_all$group)) {
    Group_tmp <- ifelse(grepl("control", delta_all$group, ignore.case = TRUE), "Controls",
                        ifelse(grepl("case|ischemic|myocardial|chronic", delta_all$group, ignore.case = TRUE), "Cases", NA))
    Group_tmp <- factor(Group_tmp, levels = c("Controls","Cases"))
  }
  delta_all$Group <- Group_tmp
  idx_delta <- idx_delta + 1; setTxtProgressBar(pb_delta, idx_delta)
}
close(pb_delta)
cat("✓ 纵向ΔIDP变化率分析完成（按匹配队列）。\n")

} else {
  cat("⚠️ 跳过 IDP 统计分析步骤。\n")
}

# =================== 第八步：完善的可视化分析代码 ===================
cat("=== 完善的缺血性心肌病IDP可视化分析 ===\n")

# =================== 公共：DK脑区名称规范化（已在文件顶部定义） ===================
load_brain_mapping_packages_once <- function() {
  pkgs <- c("ggseg","ggsegExtra","ggsegYeo2011","ggsegYeo","ggsegHO","ggsegJHU","ggseg3d","plotly")
  ap <- character(0)
  for (pkg in pkgs) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      suppressMessages(suppressWarnings(library(pkg, character.only = TRUE)))
      ap <- c(ap, pkg)
      cat(sprintf("✓ 已加载脑图包: %s\n", pkg))
    } else {
      cat(sprintf("✗ 未安装脑图包: %s\n", pkg))
    }
  }
  assign("available_packages", ap, envir = .GlobalEnv)
}

# 公共：DK脑区名称规范化（全局可用）
  if (!exists("normalize_dk_region", inherits = TRUE)) {
    normalize_dk_region <- function(s) {
      v <- tolower(s)
      v <- gsub("[_.]+", " ", v)
      v <- gsub("\\n", " ", v)
      v <- gsub("[[:space:]]+", " ", v)
      dplyr::case_when(
        stringr::str_detect(v, "superior[_\\s-]?temporal|sup\\.?\\s*temp") ~ "superiortemporal",
        stringr::str_detect(v, "middle[_\\s-]?temporal|mid\\.?\\s*temp") ~ "middletemporal",
        stringr::str_detect(v, "inferior[_\\s-]?temporal|inf\\.?\\s*temp") ~ "inferiortemporal",
        stringr::str_detect(v, "supramarginal") ~ "supramarginal",
        stringr::str_detect(v, "pars[_\\s-]?opercularis|opercularis") ~ "parsopercularis",
        stringr::str_detect(v, "pars[_\\s-]?triangularis|triangularis") ~ "parstriangularis",
        stringr::str_detect(v, "pars[_\\s-]?orbitalis|orbitalis") ~ "parsorbitalis",
        stringr::str_detect(v, "postcentral") ~ "postcentral",
        stringr::str_detect(v, "precentral") ~ "precentral",
        stringr::str_detect(v, "superior[_\\s-]?frontal|sup\\.?\\s*frontal") ~ "superiorfrontal",
        stringr::str_detect(v, "rostral[_\\s-]?middle[_\\s-]?frontal|rostral\\s*mid\\s*frontal") ~ "rostralmiddlefrontal",
        stringr::str_detect(v, "caudal[_\\s-]?middle[_\\s-]?frontal|caudal\\s*mid\\s*frontal") ~ "caudalmiddlefrontal",
        stringr::str_detect(v, "medial[_\\s-]?orbit[_\\s-]?ofrontal|medial[_\\s-]?orbito[_\\s-]?frontal|medial\\s*orbitofrontal") ~ "medialorbitofrontal",
        stringr::str_detect(v, "lateral[_\\s-]?orbit[_\\s-]?ofrontal|lateral[_\\s-]?orbito[_\\s-]?frontal|lat\\.?\\s*orbitofrontal") ~ "lateralorbitofrontal",
        stringr::str_detect(v, "rostral[_\\s-]?anterior[_\\s-]?cingulate|rostral\\s*ant\\s*cingulate") ~ "rostralanteriorcingulate",
        stringr::str_detect(v, "caudal[_\\s-]?anterior[_\\s-]?cingulate|caudal\\s*ant\\s*cingulate") ~ "caudalanteriorcingulate",
        stringr::str_detect(v, "posterior[_\\s-]?cingulate") ~ "posteriorcingulate",
        stringr::str_detect(v, "lateral[_\\s-]?occipital|lat\\.?\\s*occipital") ~ "lateraloccipital",
        stringr::str_detect(v, "\\bcuneus\\b") ~ "cuneus",
        stringr::str_detect(v, "lingual") ~ "lingual",
        stringr::str_detect(v, "peri[_\\s-]?calcarine|calcarine") ~ "pericalcarine",
        stringr::str_detect(v, "fusiform") ~ "fusiform",
        stringr::str_detect(v, "parahipp") ~ "parahippocampal",
        stringr::str_detect(v, "entorhinal") ~ "entorhinal",
        stringr::str_detect(v, "hippocamp") ~ "hippocampus",
        stringr::str_detect(v, "precuneus") ~ "precuneus",
        stringr::str_detect(v, "superior[_\\s-]?parietal|sup\\.?\\s*parietal") ~ "superiorparietal",
        stringr::str_detect(v, "inferior[_\\s-]?parietal|inf\\.?\\s*parietal") ~ "inferiorparietal",
        stringr::str_detect(v, "banks[_\\s-]?sts|superior[_\\s-]?temporal[_\\s-]?sulcus|sts\\b") ~ "bankssts",
        stringr::str_detect(v, "paracentral") ~ "paracentral",
        TRUE ~ NA_character_
      )
    }
  }
  if (!exists("normalize_wm_tract", inherits = TRUE)) {
    normalize_wm_tract <- function(s) {
      if (is.null(s)) return(NA_character_)
      v <- tolower(s)
      v <- gsub("[_.]+", " ", v)
      v <- gsub("\\n", " ", v)
      v <- gsub("[[:space:]]+", " ", v)
      dplyr::case_when(
        stringr::str_detect(v, "\\bcst\\b|corticospinal") ~ "corticospinal tract",
        stringr::str_detect(v, "\\bslf\\b|superior\\s*longitudinal") ~ "superior longitudinal fasciculus",
        stringr::str_detect(v, "\\bilf\\b|inferior\\s*longitudinal") ~ "inferior longitudinal fasciculus",
        stringr::str_detect(v, "\\baf\\b|arcuate") ~ "arcuate fasciculus",
        stringr::str_detect(v, "\\bifof\\b|inferior\\s*fronto\\s*occipital") ~ "inferior fronto-occipital fasciculus",
        stringr::str_detect(v, "\\bunc\\b|uncinate") ~ "uncinate fasciculus",
        stringr::str_detect(v, "cingulum|\\bcg\\b") ~ "cingulum",
        stringr::str_detect(v, "corpus\\s*callosum|\\bcc\\b") ~ "corpus callosum",
        stringr::str_detect(v, "optic\\s*radiation|\\bor\\b") ~ "optic radiation",
        stringr::str_detect(v, "thalamic\\s*radiation|thalamocortical|\\batr\\b") ~ "thalamic radiation",
        stringr::str_detect(v, "middle\\s*cerebellar\\s*peduncle|\\bmcp\\b") ~ "middle cerebellar peduncle",
        stringr::str_detect(v, "superior\\s*cerebellar\\s*peduncle|\\bscp\\b") ~ "superior cerebellar peduncle",
        stringr::str_detect(v, "inferior\\s*cerebellar\\s*peduncle|\\bicp\\b") ~ "inferior cerebellar peduncle",
        stringr::str_detect(v, "pontine\\s*crossing\\s*tract|\\bpct\\b") ~ "pontine crossing tract",
        TRUE ~ NA_character_
      )
    }
  }
# =================== 1. 增强版脑区映射字典 ===================
create_comprehensive_brain_dictionary <- function() {
  brain_mapping_dict <- tribble(
    ~pattern, ~brain_region, ~structure_type, ~atlas_type, ~ggseg_region, ~priority,
    
    # === 额叶区域 ===
    "superiorfrontal|superior\\.frontal|front\\.sup", "Superior Frontal", "Cortical", "cortical", "superiorfrontal", 1,
    "rostralmiddlefrontal|rostral\\.middle\\.frontal", "Rostral Middle Frontal", "Cortical", "cortical", "rostralmiddlefrontal", 1,
    "caudalmiddlefrontal|caudal\\.middle\\.frontal", "Caudal Middle Frontal", "Cortical", "cortical", "caudalmiddlefrontal", 1,
    "middlefrontal|middle\\.frontal|front\\.middle", "Middle Frontal", "Cortical", "cortical", "rostralmiddlefrontal", 2,
    "lateralorbitofrontal|lateral\\.orbitofrontal", "Lateral Orbitofrontal", "Cortical", "cortical", "lateralorbitofrontal", 1,
    "medialorbitofrontal|medial\\.orbitofrontal", "Medial Orbitofrontal", "Cortical", "cortical", "medialorbitofrontal", 1,
    "parstriangularis|pars\\.triangularis", "Pars Triangularis", "Cortical", "cortical", "parstriangularis", 1,
    "parsopercularis|pars\\.opercularis", "Pars Opercularis", "Cortical", "cortical", "parsopercularis", 1,
    "parsorbitalis|pars\\.orbitalis", "Pars Orbitalis", "Cortical", "cortical", "parsorbitalis", 1,
    "precentral|g\\.precentral", "Precentral", "Cortical", "cortical", "precentral", 1,
    "paracentral|g\\.s\\.paracentral", "Paracentral", "Cortical", "cortical", "paracentral", 1,
    "frontalpole|frontal\\.pole", "Frontal Pole", "Cortical", "cortical", "frontalpole", 1,
    
    # === 顶叶区域 ===
    "superiorparietal|superior\\.parietal", "Superior Parietal", "Cortical", "cortical", "superiorparietal", 1,
    "inferiorparietal|inferior\\.parietal", "Inferior Parietal", "Cortical", "cortical", "inferiorparietal", 1,
    "postcentral|g\\.postcentral", "Postcentral", "Cortical", "cortical", "postcentral", 1,
    "precuneus|g\\.precuneus", "Precuneus", "Cortical", "cortical", "precuneus", 1,
    "supramarginal", "Supramarginal", "Cortical", "cortical", "supramarginal", 1,
    "angular", "Angular", "Cortical", "cortical", "inferiorparietal", 1,
    
    # === 颞叶区域 ===
    "superiortemporal|superior\\.temporal", "Superior Temporal", "Cortical", "cortical", "superiortemporal", 1,
    "middletemporal|middle\\.temporal", "Middle Temporal", "Cortical", "cortical", "middletemporal", 1,
    "inferiortemporal|inferior\\.temporal", "Inferior Temporal", "Cortical", "cortical", "inferiortemporal", 1,
    "transversetemporal|transverse\\.temporal", "Transverse Temporal", "Cortical", "cortical", "transversetemporal", 1,
    "temporalpole|temporal\\.pole|pole\\.temporal", "Temporal Pole", "Cortical", "cortical", "temporalpole", 1,
    "parahippocampal", "Parahippocampal", "Cortical", "cortical", "parahippocampal", 1,
    "entorhinal", "Entorhinal", "Cortical", "cortical", "entorhinal", 1,
    "fusiform", "Fusiform", "Cortical", "cortical", "fusiform", 1,
    "bankssts|banks\\.sts", "Banks STS", "Cortical", "cortical", "bankssts", 1,
    
    # === 枕叶区域 ===
    "lateraloccipital|lateral\\.occipital", "Lateral Occipital", "Cortical", "cortical", "lateraloccipital", 1,
    "pericalcarine", "Pericalcarine", "Cortical", "cortical", "pericalcarine", 1,
    "cuneus|g\\.cuneus", "Cuneus", "Cortical", "cortical", "cuneus", 1,
    "lingual", "Lingual", "Cortical", "cortical", "lingual", 1,
    
    # === 扣带回 ===
    "rostralanteriorcingulate", "Rostral Anterior Cingulate", "Cortical", "cortical", "rostralanteriorcingulate", 1,
    "caudalanteriorcingulate", "Caudal Anterior Cingulate", "Cortical", "cortical", "caudalanteriorcingulate", 1,
    "posteriorcingulate", "Posterior Cingulate", "Cortical", "cortical", "posteriorcingulate", 1,
    "isthmuscingulate", "Isthmus Cingulate", "Cortical", "cortical", "isthmuscingulate", 1,
    
    # === 脑岛 ===
    "insula", "Insula", "Cortical", "cortical", "insula", 1,
    
    # === 皮下核团 ===
    "caudate", "Caudate", "Subcortical", "subcortical", "caudate", 1,
    "putamen", "Putamen", "Subcortical", "subcortical", "putamen", 1,
    "pallidum", "Pallidum", "Subcortical", "subcortical", "pallidum", 1,
    "accumbens", "Nucleus Accumbens", "Subcortical", "subcortical", "accumbens", 1,
    
    # === 海马体和杏仁核 ===
    "whole\\.amygdala|amygdala$", "Amygdala", "Subcortical", "subcortical", "amygdala", 1,
    "whole\\.hippocampus|hippocampus$", "Hippocampus", "Subcortical", "subcortical", "hippocampus", 1,
    "HATA", "Hippocampal Amygdaloid Transition Area", "Subcortical", "subcortical", "amygdala", 1,
    "Corticoamygdaloid.*transition|Cortico-amygdaloid.*transition|Corticoamygdaloid\\.transitio", "Corticoamygdaloid Transition", "Subcortical", "subcortical", "amygdala", 1,
    "GC[-_]?ML[-_]?DG", "Dentate Gyrus (GC-ML-DG)", "Subcortical", "subcortical", "hippocampus", 1,
    "CA1|CA2|CA3|CA4|subiculum|presubiculum|hippocampal\\.head|hippocampal\\.body", "Hippocampal Subfields", "Subcortical", "subcortical", "hippocampus", 2,
    
    # === 丘脑 ===
    "whole\\.thalamus|thalamus$", "Thalamus", "Subcortical", "subcortical", "thalamus", 1,
    "LGN|lateral\\.geniculate", "Lateral Geniculate Nucleus", "Thalamus", "subcortical", "thalamus", 0.5,
    "MGN|medial\\.geniculate", "Medial Geniculate Nucleus", "Thalamus", "subcortical", "thalamus", 0.5,
    "VPL|ventral\\.posterior\\.lateral", "Ventral Posterior Lateral", "Thalamus", "subcortical", "thalamus", 0.5,
    "VPM|ventral\\.posterior\\.medial", "Ventral Posterior Medial", "Thalamus", "subcortical", "thalamus", 0.5,
    "VA\\b|ventral\\.anterior", "Ventral Anterior", "Thalamus", "subcortical", "thalamus", 0.5,
    "VL\\b|ventral\\.lateral", "Ventral Lateral", "Thalamus", "subcortical", "thalamus", 0.5,
    "VLa|VL\\.a|ventral\\.lateral\\.anterior", "Ventral Lateral Anterior", "Thalamus", "subcortical", "thalamus", 0.5,
    "VLp|VL\\.p|ventral\\.lateral\\.posterior", "Ventral Lateral Posterior", "Thalamus", "subcortical", "thalamus", 0.5,
    "MDl|MD\\.l|mediodorsal.*lateral", "Mediodorsal Lateral", "Thalamus", "subcortical", "thalamus", 0.5,
    "VA\\.mc|VAmc|ventral\\.anterior.*magnocellular", "Ventral Anterior Magnocellular", "Thalamus", "subcortical", "thalamus", 0.5,
    "CL\\b|central\\s+lateral", "Central Lateral", "Thalamus", "subcortical", "thalamus", 0.5,
    "CM\\b|centromedian", "Centromedian", "Thalamus", "subcortical", "thalamus", 0.5,
    
    # === 小脑 ===
    "cerebellum|cerebell", "Cerebellum", "Cerebellum", "cerebellum", "", 2,
    
    # === 脑干 ===
    "whole\\.brainstem|brain\\.stem|brainstem", "Brainstem", "Brainstem", "brainstem", "", 1,
    
    # === 白质纤维束 ===
    "corpus\\.callosum", "Corpus Callosum", "White Matter", "white_matter", "", 2,
    "superior\\.longitudinal\\.fasciculus", "Superior Longitudinal Fasciculus", "White Matter", "white_matter", "", 1,
    "inferior\\.longitudinal\\.fasciculus", "Inferior Longitudinal Fasciculus", "White Matter", "white_matter", "", 1,
    "uncinate\\.fasciculus", "Uncinate Fasciculus", "White Matter", "white_matter", "", 1,
    "cingulum", "Cingulum", "White Matter", "white_matter", "", 2,
    "internal\\.capsule", "Internal Capsule", "White Matter", "white_matter", "", 2,
    "corona\\.radiata", "Corona Radiata", "White Matter", "white_matter", "", 2,
    "tract\\.", "White Matter Tract", "White Matter", "white_matter", "", 5,
    
    # === 通用匹配模式 ===
    "frontal|front", "Frontal", "Cortical", "cortical", "superiorfrontal", 3,
    "parietal|pariet", "Parietal", "Cortical", "cortical", "superiorparietal", 3,
    "temporal|temp", "Temporal", "Cortical", "cortical", "middletemporal", 3,
    "occipital|occip", "Occipital", "Cortical", "cortical", "lateraloccipital", 3,
    "cingulat|cingul", "Cingulate", "Cortical", "cortical", "rostralanteriorcingulate", 3,
    "^volume\\.", "Volume Measure", "Unknown", "unknown", "", 6,
    "^area\\.", "Surface Area Measure", "Unknown", "unknown", "", 6,
    "mean\\.thickness", "Cortical Thickness", "Unknown", "unknown", "", 6
  )
  
  measure_type_dict <- tribble(
    ~pattern, ~measure_type,
    "^area\\.", "Surface Area",
    "^volume\\.", "Volume", 
    "mean\\.thickness", "Cortical Thickness",
    "grey\\.white\\.contrast", "Grey-White Contrast",
    "weighted\\.mean\\.(fa|md|l1|l2|l3|isovf|icvf|od)", "DTI Tract Weighted",
    "mean\\.(fa|md|l1|l2|l3|isovf|icvf|od).*on\\.fa\\.skeleton", "DTI Skeleton",
    "mean\\.(fa|md|l1|l2|l3|isovf|icvf|od)", "DTI Measure",
    "median\\.bold|90th\\.percentile.*z\\.statistic|median\\.z\\.statistic", "fMRI Activation",
    "head\\.motion|tfmri\\.head\\.motion", "Head Motion",
    "percentile|bold|motion", "fMRI/Motion"
  )
  
  return(list(
    brain_mapping = brain_mapping_dict,
    measure_types = measure_type_dict
  ))
}

# =================== 2. 字典驱动解析函数 ===================
parse_ukb_idp_comprehensive <- function(idp_names, z_changes = NULL) {
  mapping_dicts <- create_comprehensive_brain_dictionary()
  brain_dict <- mapping_dicts$brain_mapping
  measure_dict <- mapping_dicts$measure_types
  # 读取外部字典（如存在），优先使用以提升覆盖率
  external_dict_path <- file.path("ischemic", "idp_mapping_dictionary.csv")
  external_dict <- NULL
  if (file.exists(external_dict_path)) {
    suppressWarnings({
      external_dict <- utils::read.csv(external_dict_path, stringsAsFactors = FALSE)
    })
  }
  
  if (is.null(z_changes) || length(z_changes) != length(idp_names)) {
    z_changes <- rep(NA_real_, length(idp_names))
  }
  idp_data <- data.frame(
    idp_variable = idp_names,
    z_change = z_changes,
    brain_region = "Unidentified",
    measure_type = "Other", 
    hemisphere = "bilateral",
    anatomical_structure = NA,
    ggseg_region = NA,
    atlas_type = "unknown",
    structure_type = "Other",
    mapping_confidence = 0,
    stringsAsFactors = FALSE
  )
  
  cat("Starting comprehensive IDP parsing...\n")
  cat("Processing", length(idp_names), "IDPs\n")
  
  for(i in 1:length(idp_names)) {
    idp_original <- idp_names[i]
    idp_clean <- tolower(idp_original)
    
    # 检测半球（更健壮：支持 left/right 及 lh/rh，含词边界）
    if (grepl("(^|\\W)(left|lh)(\\W|$)", idp_clean, perl = TRUE)) {
      idp_data$hemisphere[i] <- "left"
    } else if (grepl("(^|\\W)(right|rh)(\\W|$)", idp_clean, perl = TRUE)) {
      idp_data$hemisphere[i] <- "right"
    } else if (grepl("vermis", idp_clean)) {
      idp_data$hemisphere[i] <- "midline"
    } else {
      idp_data$hemisphere[i] <- "bilateral"
    }
    
    # 检测测量类型
    for(j in 1:nrow(measure_dict)) {
      if(grepl(measure_dict$pattern[j], idp_clean)) {
        idp_data$measure_type[i] <- measure_dict$measure_type[j]
        break
      }
    }
    
    # 先尝试外部字典映射（若存在），再回退到内部字典
    external_matched <- FALSE
    if (!is.null(external_dict)) {
      for (j2 in seq_len(nrow(external_dict))) {
        patt <- tolower(external_dict$pattern[j2])
        if (!is.na(patt) && patt != "" && grepl(patt, idp_clean, perl = TRUE)) {
          # 使用外部字典的字段进行填充
          br <- external_dict$display_name[j2]
          rg <- external_dict$roi_group[j2]
          mt <- external_dict$measure_type[j2]
          st <- NA_character_
          # 结构类型：根据 roi_group 或 measure_type 推断
          st <- dplyr::case_when(
            grepl("WM Tracts", rg, ignore.case = TRUE) ~ "White Matter",
            grepl("Cortical", mt, ignore.case = TRUE) ~ "Cortical",
            grepl("Subcortical", mt, ignore.case = TRUE) ~ "Subcortical",
            TRUE ~ idp_data$structure_type[i]
          )
          idp_data$brain_region[i]    <- dplyr::coalesce(br, rg, idp_data$brain_region[i])
          idp_data$measure_type[i]    <- dplyr::coalesce(mt, idp_data$measure_type[i])
          idp_data$structure_type[i]  <- dplyr::coalesce(st, idp_data$structure_type[i])
          idp_data$atlas_type[i]      <- ifelse(idp_data$structure_type[i] == "Cortical", "cortical",
                                                ifelse(idp_data$structure_type[i] == "Subcortical", "subcortical",
                                                       ifelse(idp_data$structure_type[i] == "White Matter", "white_matter", idp_data$atlas_type[i])))
          idp_data$mapping_confidence[i] <- max(idp_data$mapping_confidence[i], 7)
          external_matched <- TRUE
          # 不break，允许后续更具体模式继续覆盖
        }
      }
    }
    
    # 脑区映射 - 按优先级匹配（内部字典）
    best_match <- NULL
    highest_priority <- 999
    for(k in 1:nrow(brain_dict)) {
      pattern <- brain_dict$pattern[k]
      priority <- brain_dict$priority[k]
      
      if(grepl(pattern, idp_clean, perl = TRUE) && priority < highest_priority) {
        best_match <- k
        highest_priority <- priority
      }
    }
    
    # 应用最佳匹配
    if(!is.null(best_match)) {
      idp_data$brain_region[i] <- brain_dict$brain_region[best_match]
      idp_data$structure_type[i] <- brain_dict$structure_type[best_match]
      idp_data$atlas_type[i] <- brain_dict$atlas_type[best_match]
      idp_data$ggseg_region[i] <- brain_dict$ggseg_region[best_match]
      idp_data$mapping_confidence[i] <- max(idp_data$mapping_confidence[i], 10 - highest_priority)
    }
    
    # 设置解剖结构描述
    idp_data$anatomical_structure[i] <- gsub("\\.", " ", idp_original)
  }
  
  cat("IDP parsing completed.\n")
  return(idp_data)
}

# 统一箱线图函数：参数化支持 individual/grid，含标准化选项
create_idp_boxplots_unified <- function(independent_idps, analysis_data_list,
                                        output_dir = "Individual_IDP_Plots",
                                        mode = c("individual", "grid"),
                                        top_n = 12,
                                        delta_normalize_by = c("baseline_sd", "baseline_mean", "none"),
                                        percent_baseline = c("baseline_mean", "per_participant"),
                                        parallel = TRUE,
                                        workers = NULL,
                                        point_downsample_frac = 1,
                                        dpi = 600) {
  mode <- match.arg(mode)
  delta_normalize_by <- match.arg(delta_normalize_by)
  percent_baseline <- match.arg(percent_baseline)

  if (nrow(independent_idps) == 0) {
    cat("无独立效应IDP可供可视化\n")
    return(NULL)
  }
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # 颜色（与全局一致）
  group_colors <- get_group_colors()

  # grid模式：前top_n个IDP组合展示绝对（可标准化）
  if (mode == "grid") {
    top_idps <- utils::head(independent_idps, top_n)
    plot_list <- list()
    for (i in seq_len(nrow(top_idps))) {
      idp_info <- top_idps[i, ]
      model_name <- idp_info$model_comparison
      idp1_var <- idp_info$idp1_variable
      idp2_var <- idp_info$idp2_variable
      if (!(model_name %in% names(analysis_data_list))) next
      data <- analysis_data_list[[model_name]]
      if (!all(c(idp1_var, idp2_var) %in% names(data))) next
      pd <- data %>%
        dplyr::filter(!is.na(.data[[idp1_var]]), !is.na(.data[[idp2_var]])) %>%
        dplyr::mutate(IDP1 = .data[[idp1_var]], IDP2 = .data[[idp2_var]], Raw_Change = IDP2 - IDP1, Group = comparison_group)
      if (nrow(pd) == 0) next
      bmean <- mean(pd$IDP1, na.rm = TRUE)
      bsd <- stats::sd(pd$IDP1, na.rm = TRUE)
      sd_safe <- ifelse(is.na(bsd) || bsd == 0, 1e-8, bsd)
      mean_safe <- ifelse(is.na(bmean) || bmean == 0, 1e-8, abs(bmean))
      pd <- pd %>% dplyr::mutate(IDP_Change_Std = dplyr::case_when(
        delta_normalize_by == "baseline_sd" ~ Raw_Change / sd_safe,
        delta_normalize_by == "baseline_mean" ~ Raw_Change / mean_safe,
        TRUE ~ Raw_Change
      ))
      ylab_abs <- if (delta_normalize_by == "baseline_sd") "IDP Change (SD units)" else if (delta_normalize_by == "baseline_mean") "IDP Change (per baseline mean)" else "IDP Change (Instance3 - Instance2)"
      p <- ggplot2::ggplot(pd, ggplot2::aes(x = Group, y = IDP_Change_Std, fill = Group)) +
        ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA) +
        ggplot2::geom_jitter(width = 0.2, alpha = 0.5, size = 0.8) +
        ggplot2::scale_fill_manual(values = group_colors) +
        ggplot2::labs(title = paste("IDP Change:", substr(idp_info$variable_name, 1, 40)), subtitle = paste("Model:", model_name), x = NULL, y = ylab_abs) +
        theme_pub() + ggplot2::theme(legend.position = "none")
      plot_list[[length(plot_list) + 1]] <- p
    }
    if (length(plot_list) == 0) return(NULL)
    combined_plot <- do.call(gridExtra::grid.arrange, c(plot_list, ncol = 4))
    save_pub(file.path(output_dir, "Independent_IDPs_Boxplots_Grid"), combined_plot, width = 16, height = 12, dpi = 600)
    cat("已保存分组箱线图网格: ", file.path(output_dir, "Independent_IDPs_Boxplots_Grid"), ".png/pdf\n", sep = "")
    return(invisible(combined_plot))
  }

  # individual模式：每个IDP两面板（绝对+百分比），支持标准化
  plot_idps <- if (!is.null(top_n) && is.finite(top_n) && top_n > 0) {
    utils::head(independent_idps, top_n)
  } else {
    independent_idps
  }
  point_downsample_frac <- suppressWarnings(as.numeric(point_downsample_frac))
  if (!is.finite(point_downsample_frac) || point_downsample_frac <= 0) point_downsample_frac <- 1
  if (point_downsample_frac > 1) point_downsample_frac <- 1

  get_sample <- function(df, frac) {
    if (!is.finite(frac) || frac >= 1) return(df)
    n_target <- max(1L, as.integer(ceiling(nrow(df) * frac)))
    suppressWarnings(dplyr::slice_sample(df, n = min(n_target, nrow(df))))
  }

  render_one <- function(i) {
    idp_info <- plot_idps[i, ]
    model_name <- idp_info$model_comparison
    idp1_var <- idp_info$idp1_variable
    idp2_var <- idp_info$idp2_variable
    variable_name <- idp_info$variable_name
    if (!(model_name %in% names(analysis_data_list))) return(NULL)
    data <- analysis_data_list[[model_name]]
    if (!all(c(idp1_var, idp2_var) %in% names(data))) return(NULL)
    pd <- data %>%
      dplyr::filter(!is.na(.data[[idp1_var]]), !is.na(.data[[idp2_var]])) %>%
      dplyr::mutate(IDP1 = .data[[idp1_var]], IDP2 = .data[[idp2_var]], Raw_Change = IDP2 - IDP1, Group = comparison_group)
    if (nrow(pd) == 0) return(NULL)
    bmean <- mean(pd$IDP1, na.rm = TRUE)
    bsd <- stats::sd(pd$IDP1, na.rm = TRUE)
    sd_safe <- ifelse(is.na(bsd) || bsd == 0, 1e-8, bsd)
    mean_safe <- ifelse(is.na(bmean) || bmean == 0, 1e-8, abs(bmean))
    pd <- pd %>% dplyr::mutate(
      IDP_Change = Raw_Change,
      IDP_Change_Std = dplyr::case_when(
        delta_normalize_by == "baseline_sd" ~ Raw_Change / sd_safe,
        delta_normalize_by == "baseline_mean" ~ Raw_Change / mean_safe,
        TRUE ~ Raw_Change
      ),
      IDP_Change_Percent = dplyr::case_when(
        percent_baseline == "per_participant" & !is.na(IDP1) & IDP1 != 0 ~ (Raw_Change) / abs(IDP1) * 100,
        percent_baseline == "baseline_mean" ~ (Raw_Change) / mean_safe * 100,
        TRUE ~ NA_real_
      )
    )
    pd_pts <- get_sample(pd, point_downsample_frac)
    ylab_abs <- if (delta_normalize_by == "baseline_sd") "IDP Change (SD units)" else if (delta_normalize_by == "baseline_mean") "IDP Change (per baseline mean)" else "IDP Change (Instance3 - Instance2)"
    ylab_pct <- if (percent_baseline == "baseline_mean") "IDP Change (% of baseline mean)" else "IDP Change (%)"
    p1 <- ggplot2::ggplot(pd, ggplot2::aes(x = Group, y = IDP_Change_Std, fill = Group)) +
      ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6, color = "black", linewidth = 0.5) +
      ggplot2::geom_jitter(data = pd_pts, width = 0.2, alpha = 0.4, size = 0.8, color = "black") +
      ggplot2::scale_fill_manual(values = group_colors) +
      ggplot2::labs(title = paste("IDP Absolute Change"), subtitle = paste(substr(variable_name, 1, 50)), x = "", y = ylab_abs) +
      theme_pub() + ggplot2::theme(legend.position = "none")
    p2 <- ggplot2::ggplot(pd, ggplot2::aes(x = Group, y = IDP_Change_Percent, fill = Group)) +
      ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6, color = "black", linewidth = 0.5) +
      ggplot2::geom_jitter(data = pd_pts, width = 0.2, alpha = 0.4, size = 0.8, color = "black") +
      ggplot2::scale_fill_manual(values = group_colors) +
      ggplot2::labs(title = paste("IDP Percentage Change"), subtitle = paste("Model:", gsub("_", " ", model_name)), x = "", y = ylab_pct) +
      theme_pub() + ggplot2::theme(legend.position = "none")
    combined_plot <- cowplot::plot_grid(p1, p2, ncol = 2, rel_widths = c(1, 1))
    filename_base <- paste0("IDP_", sprintf("%03d", i), "_", gsub("[^A-Za-z0-9]", "_", substr(variable_name, 1, 30)))
    save_pub(file.path(output_dir, filename_base), combined_plot, width = 12, height = 6, dpi = dpi)
    list(plot = combined_plot, data = pd, info = idp_info)
  }

  results <- NULL
  use_par <- isTRUE(parallel) && requireNamespace("future", quietly = TRUE) && requireNamespace("furrr", quietly = TRUE)
  if (use_par) {
    old_plan <- future::plan()
    on.exit({ if (exists("old_plan") && !is.null(old_plan)) future::plan(old_plan) }, add = TRUE)
    if (is.null(workers)) {
      workers <- max(1L, (if (suppressWarnings(requireNamespace("parallel", quietly = TRUE))) parallel::detectCores() else 2L) - 1L)
    }
    future::plan(future::multisession, workers = workers)
    map_opts <- furrr::furrr_options(seed = NULL, packages = c("ggplot2","dplyr","stringr","cowplot"))
    results <- furrr::future_map(seq_len(nrow(plot_idps)), render_one, .options = map_opts)
  } else {
    results <- lapply(seq_len(nrow(plot_idps)), render_one)
  }

  valid_results <- results[!sapply(results, is.null)]
  cat("已保存", length(valid_results), "个统一箱线图至目录:", output_dir, "\n")
  invisible(valid_results)
}

# =================== 3. 可视化IDP选择（含回退机制） ===================
get_visualization_idps <- function(independent_idps, combined_results = NULL, fallback_n = NULL) {
  combined_src <- NULL
  if (!is.null(combined_results) && nrow(combined_results) > 0) {
    combined_src <- tibble::as_tibble(combined_results)
  } else {
    csv_path <- combined_four_model_results_csv()
    if (file.exists(csv_path)) {
      combined_src <- tryCatch({
        tmp <- utils::read.csv(csv_path, stringsAsFactors = FALSE)
        tmp <- clean_combined_results_df(tmp)
        tibble::as_tibble(tmp)
      }, error = function(e) NULL)
    }
  }

  if (!is.null(combined_src) && nrow(combined_src) > 0) {
    df <- combined_src
    if ("variable_type" %in% names(df)) {
      df <- df[df$variable_type %in% c("Brain_Imaging", "Cognitive"), , drop = FALSE]
    }
    if (all(c("idp1_variable", "idp2_variable") %in% names(df))) {
      df <- df[!is.na(df$idp1_variable) & nzchar(as.character(df$idp1_variable)) &
                 !is.na(df$idp2_variable) & nzchar(as.character(df$idp2_variable)), , drop = FALSE]
    }
    if (all(c("model_comparison", "variable_name") %in% names(df))) {
      df <- dplyr::distinct(df, model_comparison, variable_name, .keep_all = TRUE)
    }
    if (nrow(df) > 0) {
      z_num <- if ("z_statistic_without_non_imaging" %in% names(df)) {
        suppressWarnings(as.numeric(df$z_statistic_without_non_imaging))
      } else if ("z_statistic_with_non_imaging" %in% names(df)) {
        suppressWarnings(as.numeric(df$z_statistic_with_non_imaging))
      } else {
        rep(NA_real_, nrow(df))
      }
      df$sort_key <- abs(z_num)
      df <- df[order(-df$sort_key), , drop = FALSE]
      df$sort_key <- NULL
      if (!is.null(fallback_n) && is.finite(fallback_n) && fallback_n > 0 && nrow(df) > fallback_n) {
        df <- df[seq_len(fallback_n), , drop = FALSE]
      }
      return(tibble::as_tibble(df))
    }
  }

  if (!is.null(independent_idps) && nrow(independent_idps) > 0) {
    return(tibble::as_tibble(independent_idps))
  }

  return(data.frame())
}

# 简单统计：哪些IDP在对应模型数据中具备完整的两个时间点变量
summarize_plot_eligibility <- function(idps, analysis_data_list) {
  if (is.null(idps) || nrow(idps) == 0) return(list(total = 0, eligible = 0))
  ok <- 0
  for (i in seq_len(nrow(idps))) {
    model_name <- idps$model_comparison[i]
    idp1 <- idps$idp1_variable[i]
    idp2 <- idps$idp2_variable[i]
    if (model_name %in% names(analysis_data_list)) {
      dat <- analysis_data_list[[model_name]]
      if (all(c(idp1, idp2) %in% names(dat))) ok <- ok + 1
    }
  }
  list(total = nrow(idps), eligible = ok)
}


# =================== 4. 年龄相关变化和变化率趋势图 ===================

# 统一年龄趋势绘图函数：参数化支持 percentage/absolute/both，默认percentage单面版带散点
create_age_trend_plots <- function(independent_idps,
                                   analysis_data_list,
                                   output_dir = "Age_Trend_Plots",
                                   top_n = 20,
                                   min_group_n = 5,
                                   panel = c("percentage", "absolute", "both"),
                                   include_points = TRUE,
                                   smooth_span = 0.8,
                                   delta_normalize_by = c("baseline_sd", "baseline_mean", "none"),
                                   percent_baseline = c("baseline_mean", "per_participant"),
                                   dpi = 600) {
  panel <- match.arg(panel)
  delta_normalize_by <- match.arg(delta_normalize_by)
  percent_baseline <- match.arg(percent_baseline)
  if (nrow(independent_idps) == 0) return(NULL)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # 选择最显著的IDP
  top_idps <- utils::head(independent_idps, top_n)
  completeness_records <- list()

  # 颜色映射（统一来源）
  group_colors <- get_group_colors()

  for (i in seq_len(nrow(top_idps))) {
    idp_info <- top_idps[i, ]
    model_name <- idp_info$model_comparison
    idp1_var <- idp_info$idp1_variable
    idp2_var <- idp_info$idp2_variable
    variable_name <- idp_info$variable_name

    if (model_name %in% names(analysis_data_list)) {
      data <- analysis_data_list[[model_name]]
      if (all(c(idp1_var, idp2_var) %in% names(data))) {
        plot_data <- data %>%
          dplyr::filter(!is.na(.data[[idp1_var]]), !is.na(.data[[idp2_var]])) %>%
          dplyr::mutate(
            IDP1 = .data[[idp1_var]],
            IDP2 = .data[[idp2_var]],
            Raw_Change = IDP2 - IDP1,
            Group = comparison_group
          ) %>%
          dplyr::filter(!is.na(Group), !is.na(Age2))

        # 基线统计（参与者的首次扫描整体分布）
        baseline_mean <- mean(plot_data$IDP1, na.rm = TRUE)
        baseline_sd <- stats::sd(plot_data$IDP1, na.rm = TRUE)
        sd_safe <- ifelse(is.na(baseline_sd) || baseline_sd == 0, 1e-8, baseline_sd)
        mean_safe <- ifelse(is.na(baseline_mean) || baseline_mean == 0, 1e-8, abs(baseline_mean))

        # 标准化绝对变化 + 百分比变化（基线可选）
        plot_data <- plot_data %>%
          dplyr::mutate(
            IDP_Change = Raw_Change,
            IDP_Change_Std = dplyr::case_when(
              delta_normalize_by == "baseline_sd" ~ Raw_Change / sd_safe,
              delta_normalize_by == "baseline_mean" ~ Raw_Change / mean_safe,
              TRUE ~ Raw_Change
            ),
            IDP_Change_Percent = dplyr::case_when(
              percent_baseline == "per_participant" & !is.na(IDP1) & IDP1 != 0 ~ (Raw_Change) / abs(IDP1) * 100,
              percent_baseline == "baseline_mean" ~ (Raw_Change) / mean_safe * 100,
              TRUE ~ NA_real_
            )
          )

        # 计算两时间点完整样本量（按组），并应用覆盖度门槛
        valid_groups <- if (exists("cmp_labels") && model_name %in% names(cmp_labels)) cmp_labels[[model_name]] else sort(unique(plot_data$Group))
        counts_df <- plot_data %>% dplyr::filter(Group %in% valid_groups) %>% dplyr::group_by(Group) %>% dplyr::summarise(n_complete = dplyr::n(), .groups = "drop")
        for (g in valid_groups) {
          if (!g %in% counts_df$Group) counts_df <- dplyr::bind_rows(counts_df, dplyr::tibble(Group = g, n_complete = 0L))
        }
        counts_df <- counts_df %>% dplyr::arrange(factor(Group, levels = valid_groups))
        eligible_groups <- counts_df %>% dplyr::filter(n_complete >= min_group_n) %>% dplyr::pull(Group)
        status <- if (length(eligible_groups) >= 2) "two_group" else if (length(eligible_groups) == 1) "single_group" else "skipped"

        # 记录完备性
        completeness_records[[length(completeness_records) + 1]] <- dplyr::tibble(
          panel_type = panel,
          model_comparison = model_name,
          variable_name = variable_name,
          group_a = ifelse(length(valid_groups) >= 1, valid_groups[1], NA_character_),
          group_b = ifelse(length(valid_groups) >= 2, valid_groups[2], NA_character_),
          n_a_complete = ifelse(length(valid_groups) >= 1, counts_df$n_complete[counts_df$Group == valid_groups[1]][1], NA_integer_),
          n_b_complete = ifelse(length(valid_groups) >= 2, counts_df$n_complete[counts_df$Group == valid_groups[2]][1], NA_integer_),
          status = status,
          threshold = min_group_n
        )

        # 根据覆盖度门槛筛选可绘制的组
        if (status == "skipped") {
          next
        }
        plot_data_eff <- if (status == "two_group") {
          plot_data %>% dplyr::filter(Group %in% valid_groups)
        } else {
          plot_data %>% dplyr::filter(Group %in% eligible_groups)
        }
        cap_text <- if (status == "single_group") paste0("另一组覆盖度不足(n<", min_group_n, ")，仅显示: ", eligible_groups[1]) else NULL

        # 绘图：根据panel参数
        if (panel == "percentage") {
          p <- ggplot2::ggplot(plot_data_eff, ggplot2::aes(x = Age2, y = IDP_Change_Percent, color = Group)) +
            { if (isTRUE(include_points)) ggplot2::geom_point(alpha = 0.5, size = 0.8) else ggplot2::geom_blank() } +
            ggplot2::geom_smooth(method = "loess", se = TRUE, span = smooth_span, alpha = 0.2, linewidth = 1.2) +
            ggplot2::scale_color_manual(values = group_colors, drop = FALSE) +
            ggplot2::labs(
              title = paste0(substr(variable_name, 1, 40)),
              subtitle = "Percentage Change vs Age",
              x = "Age (years)", y = if (percent_baseline == "baseline_mean") "ΔIDP (% of baseline mean)" else "ΔIDP (%)", color = "Group", caption = cap_text
            ) +
            theme_jama_refined() +
            ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1, size = 2)))
          filename_base <- paste0("Age_Trend_Percentage_", sprintf("%03d", i), "_", gsub("[^A-Za-z0-9]", "_", substr(variable_name, 1, 30)))
          save_pub(file.path(output_dir, filename_base), p, width = 6, height = 6, dpi = dpi)
        } else if (panel == "absolute") {
          p <- ggplot2::ggplot(plot_data_eff, ggplot2::aes(x = Age2, y = IDP_Change_Std, color = Group)) +
            { if (isTRUE(include_points)) ggplot2::geom_point(alpha = 0.5, size = 0.8) else ggplot2::geom_blank() } +
            ggplot2::geom_smooth(method = "loess", se = TRUE, span = smooth_span, alpha = 0.2, linewidth = 1.2) +
            ggplot2::scale_color_manual(values = group_colors, drop = FALSE) +
            ggplot2::labs(
              title = paste0(substr(variable_name, 1, 40)),
              subtitle = paste(
                "Absolute Change vs Age | Z =",
                round(suppressWarnings(as.numeric(if ("z_statistic_without_non_imaging" %in% names(idp_info)) idp_info$z_statistic_without_non_imaging else NA_real_)), 3)
              ),
              x = "Age (years)", y = if (delta_normalize_by == "baseline_sd") "ΔIDP (SD units)" else if (delta_normalize_by == "baseline_mean") "ΔIDP (per baseline mean)" else "Absolute Change", color = "Group", caption = cap_text
            ) +
            theme_jama_refined() +
            ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1, size = 2)))
          filename_base <- paste0("Age_Trend_Absolute_", sprintf("%03d", i), "_", gsub("[^A-Za-z0-9]", "_", substr(variable_name, 1, 30)))
          save_pub(file.path(output_dir, filename_base), p, width = 6, height = 6, dpi = dpi)
        } else {
          p1 <- ggplot2::ggplot(plot_data_eff, ggplot2::aes(x = Age2, y = IDP_Change_Std, color = Group)) +
            { if (isTRUE(include_points)) ggplot2::geom_point(alpha = 0.5, size = 0.8) else ggplot2::geom_blank() } +
            ggplot2::geom_smooth(method = "loess", se = TRUE, span = smooth_span, alpha = 0.2, linewidth = 1.2) +
            ggplot2::scale_color_manual(values = group_colors, drop = FALSE) +
            ggplot2::labs(
              title = paste(substr(variable_name, 1, 40)),
              subtitle = paste(
                "Absolute Change vs Age | Z =",
                round(suppressWarnings(as.numeric(if ("z_statistic_without_non_imaging" %in% names(idp_info)) idp_info$z_statistic_without_non_imaging else NA_real_)), 3)
              ),
              x = "Age (years)", y = if (delta_normalize_by == "baseline_sd") "ΔIDP (SD units)" else if (delta_normalize_by == "baseline_mean") "ΔIDP (per baseline mean)" else "Absolute Change", color = "Group", caption = cap_text
            ) + theme_jama_refined() + ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1, size = 2)))

          p2 <- ggplot2::ggplot(plot_data_eff, ggplot2::aes(x = Age2, y = IDP_Change_Percent, color = Group)) +
            { if (isTRUE(include_points)) ggplot2::geom_point(alpha = 0.5, size = 0.8) else ggplot2::geom_blank() } +
            ggplot2::geom_smooth(method = "loess", se = TRUE, span = smooth_span, alpha = 0.2, linewidth = 1.2) +
            ggplot2::scale_color_manual(values = group_colors, drop = FALSE) +
            ggplot2::labs(
              title = paste(substr(variable_name, 1, 40)),
              subtitle = "Percentage Change vs Age",
              x = "Age (years)", y = if (percent_baseline == "baseline_mean") "ΔIDP (% of baseline mean)" else "ΔIDP (%)", color = "Group", caption = cap_text
            ) + theme_jama_refined() + ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1, size = 2)))

          combined_plot <- cowplot::plot_grid(p1, p2, ncol = 1, rel_heights = c(1, 1))
          filename_base <- paste0("Age_Trend_Both_", sprintf("%03d", i), "_", gsub("[^A-Za-z0-9]", "_", substr(variable_name, 1, 30)))
          save_pub(file.path(output_dir, filename_base), combined_plot, width = 10, height = 12, dpi = dpi)
        }
      }
    }
  }

  completeness_df <- if (length(completeness_records) > 0) dplyr::bind_rows(completeness_records) else NULL
  if (!is.null(completeness_df)) {
    utils::write.csv(completeness_df, file.path(output_dir, "Age_Trend_Completeness.csv"), row.names = FALSE)
    cat("已写入统一趋势完备性报表: ", file.path(output_dir, "Age_Trend_Completeness.csv"), "\n")
  }
  return(completeness_df)
}

# 基于个体级数据的年龄×组别斜率差异检验：验证Case组随年龄变化率高于Control
compute_age_slope_difference_for_idps <- function(independent_idps,
                                                  analysis_data_list,
                                                  used_metric = c("percent_per_participant", "sd_units", "baseline_mean_percent"),
                                                  min_group_n = 5,
                                                  output_path = file.path("Age_Trend_Plots", "Age_Trend_Slope_Difference_AllModels.csv")) {
  used_metric <- match.arg(used_metric)
  if (nrow(independent_idps) == 0) return(invisible(NULL))
  out_rows <- list()

  for (i in seq_len(nrow(independent_idps))) {
    idp_info <- independent_idps[i, ]
    model_name <- idp_info$model_comparison
    idp1_var <- idp_info$idp1_variable
    idp2_var <- idp_info$idp2_variable
    variable_name <- idp_info$variable_name
    if (!(model_name %in% names(analysis_data_list))) next
    data <- analysis_data_list[[model_name]]
    if (!all(c(idp1_var, idp2_var, "Age2", "comparison_group") %in% names(data))) next

    df <- data %>%
      dplyr::filter(!is.na(.data[[idp1_var]]), !is.na(.data[[idp2_var]]), !is.na(Age2), !is.na(comparison_group)) %>%
      dplyr::mutate(IDP1 = .data[[idp1_var]], IDP2 = .data[[idp2_var]], Raw_Change = IDP2 - IDP1)
    if (nrow(df) == 0) next

    # 归一化选择
    bmean <- mean(df$IDP1, na.rm = TRUE)
    bsd   <- stats::sd(df$IDP1, na.rm = TRUE)
    sd_safe <- ifelse(is.na(bsd) || bsd == 0, 1e-8, bsd)
    mean_safe <- ifelse(is.na(bmean) || bmean == 0, 1e-8, abs(bmean))
    df <- df %>% dplyr::mutate(
      y_metric = dplyr::case_when(
        used_metric == "percent_per_participant" & !is.na(IDP1) & IDP1 != 0 ~ (Raw_Change) / abs(IDP1) * 100,
        used_metric == "baseline_mean_percent" ~ (Raw_Change) / mean_safe * 100,
        used_metric == "sd_units" ~ Raw_Change / sd_safe,
        TRUE ~ NA_real_
      )
    )
    if (exists("cmp_labels") && (model_name %in% names(cmp_labels))) {
      pos_lab <- cmp_labels[[model_name]][1]
      neg_lab <- cmp_labels[[model_name]][2]
      df$Group <- dplyr::case_when(
        df$comparison_group == pos_lab ~ "Case",
        df$comparison_group == neg_lab ~ "Control",
        TRUE ~ NA_character_
      )
    } else {
      df$Group <- dplyr::case_when(
        tolower(df$comparison_group) %in% c("control", "healthy") ~ "Control",
        TRUE ~ "Case"
      )
    }
    df$Group <- factor(df$Group, levels = c("Control", "Case"))
    df <- df %>% dplyr::filter(!is.na(Group))
    counts_df <- df %>% dplyr::group_by(Group) %>% dplyr::summarise(n = dplyr::n(), .groups = "drop")
    if (any(is.na(df$y_metric)) || any(!is.finite(df$y_metric))) df <- df %>% dplyr::filter(is.finite(y_metric))

    # 组覆盖度阈值（安全提取，防止NA索引导致长度>1）
    n_case <- counts_df %>% dplyr::filter(Group == "Case") %>% dplyr::pull(n)
    n_ctrl <- counts_df %>% dplyr::filter(Group == "Control") %>% dplyr::pull(n)
    n_case <- if (length(n_case) == 0) 0L else n_case[1]
    n_ctrl <- if (length(n_ctrl) == 0) 0L else n_ctrl[1]
    if (n_case < min_group_n || n_ctrl < min_group_n) {
      out_rows[[length(out_rows) + 1]] <- dplyr::tibble(
        model_comparison = model_name,
        variable_name = variable_name,
        slope_control = NA_real_,
        slope_case = NA_real_,
        slope_diff = NA_real_,
        p_interaction = NA_real_,
        n_control = n_ctrl,
        n_case = n_case,
        used_metric = used_metric,
        status = "insufficient_group"
      )
      next
    }

    # 线性模型：y ~ Age2 * Group (+sex_factor若存在)
    has_sex <- "sex_factor" %in% names(df)
    f_str <- if (has_sex) "y_metric ~ Age2 * Group + sex_factor" else "y_metric ~ Age2 * Group"
    fit <- tryCatch(stats::lm(stats::as.formula(f_str), data = df), error = function(e) NULL)
    if (is.null(fit)) next
    ct <- tryCatch(coef(summary(fit)), error = function(e) NULL)
    if (is.null(ct)) next

    # 提取斜率系数
    b_age <- if ("Age2" %in% rownames(ct)) as.numeric(ct["Age2", "Estimate"]) else NA_real_
    b_int <- if ("Age2:GroupCase" %in% rownames(ct)) as.numeric(ct["Age2:GroupCase", "Estimate"]) else NA_real_
    p_int <- if ("Age2:GroupCase" %in% rownames(ct)) as.numeric(ct["Age2:GroupCase", "Pr(>|t|)"]) else NA_real_

    out_rows[[length(out_rows) + 1]] <- dplyr::tibble(
      model_comparison = model_name,
      variable_name = variable_name,
      slope_control = b_age,
      slope_case = b_age + b_int,
      slope_diff = b_int,
      p_interaction = p_int,
      n_control = n_ctrl,
      n_case = n_case,
      used_metric = used_metric,
      status = "ok"
    )
  }

  out_df <- if (length(out_rows) > 0) dplyr::bind_rows(out_rows) else NULL
  if (!is.null(out_df) && nrow(out_df) > 0) {
    utils::write.csv(out_df, output_path, row.names = FALSE)
    cat("已输出年龄×组别斜率差异检验: ", output_path, "\n")
  }
  invisible(out_df)
}
# 并行+降采样版本：统一年龄趋势绘图，提升绘制效率
create_age_trend_plots_fast <- function(independent_idps,
                                        analysis_data_list,
                                        output_dir = "Age_Trend_Plots",
                                        top_n = NULL,
                                        min_group_n = 5,
                                        panel = c("percentage", "absolute", "both"),
                                        include_points = TRUE,
                                        smooth_span = 0.8,
                                        delta_normalize_by = c("baseline_sd", "baseline_mean", "none"),
                                        percent_baseline = c("baseline_mean", "per_participant"),
                                        parallel = TRUE,
                                        workers = NULL,
                                        point_downsample_frac = 0.5,
                                        smooth_sample_frac = 0.5,
                                        smooth_sample_max = 5000,
                                        se_confint = TRUE,
                                        dpi = 600,
                                        export_raw_data = TRUE) {
  panel <- match.arg(panel)
  delta_normalize_by <- match.arg(delta_normalize_by)
  percent_baseline <- match.arg(percent_baseline)
  if (nrow(independent_idps) == 0) return(NULL)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # 若未指定top_n或无效，则绘制所有独立效应IDP
  if (is.null(top_n) || !is.finite(top_n) || top_n <= 0) {
    top_idps <- independent_idps
  } else {
    top_idps <- utils::head(independent_idps, top_n)
  }
  completeness_records <- list()

  # 原始数据导出目录（可选）
  raw_dir <- file.path(output_dir, "Age_Trend_Raw_Data")
  if (isTRUE(export_raw_data) && !dir.exists(raw_dir)) dir.create(raw_dir, recursive = TRUE)

  group_colors <- get_group_colors()

  # 并行配置
  if (isTRUE(parallel)) {
    has_future <- suppressWarnings(requireNamespace("future", quietly = TRUE))
    has_furrr  <- suppressWarnings(requireNamespace("furrr", quietly = TRUE))
    if (has_future && has_furrr) {
      old_plan <- future::plan()
on.exit({ if (exists("old_plan") && !is.null(old_plan)) future::plan(old_plan) }, add = TRUE)
      if (is.null(workers)) {
        workers <- max(1L, (if (suppressWarnings(requireNamespace("parallel", quietly = TRUE))) parallel::detectCores() else 2L) - 1L)
      }
      future::plan(future::multisession, workers = workers)
      mapper <- furrr::future_map
      map_opts <- furrr::furrr_options(seed = NULL, packages = c("ggplot2","dplyr","stringr","cowplot"))
    } else {
      warning("并行被请求，但未发现 future/furrr 包，改为顺序执行。")
      mapper <- purrr::map
      map_opts <- NULL
    }
  } else {
    mapper <- purrr::map
    map_opts <- NULL
  }

  get_group_sample <- function(df, frac, max_n) {
    if (!is.finite(frac) || frac >= 1) return(df)
    if (!is.finite(max_n)) max_n <- nrow(df)
    n_target <- max(1L, min(nrow(df), as.integer(ceiling(nrow(df) * frac)), max_n))
    suppressWarnings(dplyr::slice_sample(df, n = n_target))
  }

  render_one <- function(i) {
    idp_info <- top_idps[i, ]
    model_name <- idp_info$model_comparison
    idp1_var <- idp_info$idp1_variable
    idp2_var <- idp_info$idp2_variable
    variable_name <- idp_info$variable_name

    out_rec <- NULL
    if (model_name %in% names(analysis_data_list)) {
      data <- analysis_data_list[[model_name]]
      if (all(c(idp1_var, idp2_var) %in% names(data))) {
        plot_data <- data %>%
          dplyr::filter(!is.na(.data[[idp1_var]]), !is.na(.data[[idp2_var]])) %>%
          dplyr::mutate(
            IDP1 = .data[[idp1_var]],
            IDP2 = .data[[idp2_var]],
            Raw_Change = IDP2 - IDP1,
            Group = comparison_group
          ) %>%
          dplyr::filter(!is.na(Group), !is.na(Age2))

        baseline_mean <- mean(plot_data$IDP1, na.rm = TRUE)
        baseline_sd <- stats::sd(plot_data$IDP1, na.rm = TRUE)
        sd_safe <- ifelse(is.na(baseline_sd) || baseline_sd == 0, 1e-8, baseline_sd)
        mean_safe <- ifelse(is.na(baseline_mean) || baseline_mean == 0, 1e-8, abs(baseline_mean))

        plot_data <- plot_data %>%
          dplyr::mutate(
            IDP_Change = Raw_Change,
            IDP_Change_Std = dplyr::case_when(
              delta_normalize_by == "baseline_sd" ~ Raw_Change / sd_safe,
              delta_normalize_by == "baseline_mean" ~ Raw_Change / mean_safe,
              TRUE ~ Raw_Change
            ),
            IDP_Change_Percent = dplyr::case_when(
              percent_baseline == "per_participant" & !is.na(IDP1) & IDP1 != 0 ~ (Raw_Change) / abs(IDP1) * 100,
              percent_baseline == "baseline_mean" ~ (Raw_Change) / mean_safe * 100,
              TRUE ~ NA_real_
            )
          )

        # 导出每个IDP的个体级原始纵向数据（含归一化）
        if (isTRUE(export_raw_data)) {
          safe_name <- gsub("[^A-Za-z0-9]", "_", substr(variable_name, 1, 60))
          # 修复长度不匹配：使用plot_data中的Participant.ID以保证长度与当前数据一致
          id_col <- if ("eid" %in% names(plot_data)) plot_data$eid else rep(NA, nrow(plot_data))
          raw_out <- plot_data %>%
            dplyr::mutate(
              eid = id_col,
              delta_rate = dplyr::case_when(
                !is.na(IDP1) & IDP1 != 0 ~ (Raw_Change) / abs(IDP1),
                TRUE ~ NA_real_
              )
            ) %>%
            dplyr::transmute(
              model_comparison = model_name,
              variable_name = variable_name,
              idp1_variable = idp1_var,
              idp2_variable = idp2_var,
              eid = eid,
              Group = Group,
              Age2 = Age2,
              IDP1 = IDP1,
              IDP2 = IDP2,
              Raw_Change = Raw_Change,
              IDP_Change_Std = IDP_Change_Std,
              IDP_Change_Percent = IDP_Change_Percent,
              delta_rate = delta_rate
            )
          utils::write.csv(raw_out, file.path(raw_dir, paste0("Age_Trend_RawData_", sprintf("%03d", i), "_", safe_name, "_", model_name, ".csv")), row.names = FALSE)
        }

        valid_groups <- if (exists("cmp_labels") && model_name %in% names(cmp_labels)) cmp_labels[[model_name]] else sort(unique(plot_data$Group))
        counts_df <- plot_data %>% dplyr::filter(Group %in% valid_groups) %>% dplyr::group_by(Group) %>% dplyr::summarise(n_complete = dplyr::n(), .groups = "drop")
        for (g in valid_groups) {
          if (!g %in% counts_df$Group) counts_df <- dplyr::bind_rows(counts_df, dplyr::tibble(Group = g, n_complete = 0L))
        }
        counts_df <- counts_df %>% dplyr::arrange(factor(Group, levels = valid_groups))
        eligible_groups <- counts_df %>% dplyr::filter(n_complete >= min_group_n) %>% dplyr::pull(Group)
        status <- if (length(eligible_groups) >= 2) "two_group" else if (length(eligible_groups) == 1) "single_group" else "skipped"

        out_rec <- dplyr::tibble(
          panel_type = panel,
          model_comparison = model_name,
          variable_name = variable_name,
          group_a = ifelse(length(valid_groups) >= 1, valid_groups[1], NA_character_),
          group_b = ifelse(length(valid_groups) >= 2, valid_groups[2], NA_character_),
          n_a_complete = ifelse(length(valid_groups) >= 1, counts_df$n_complete[counts_df$Group == valid_groups[1]][1], NA_integer_),
          n_b_complete = ifelse(length(valid_groups) >= 2, counts_df$n_complete[counts_df$Group == valid_groups[2]][1], NA_integer_),
          status = status,
          threshold = min_group_n
        )

        if (status == "skipped") return(out_rec)
        plot_data_eff <- if (status == "two_group") plot_data %>% dplyr::filter(Group %in% valid_groups) else plot_data %>% dplyr::filter(Group %in% eligible_groups)
        cap_text <- if (status == "single_group") paste0("另一组覆盖度不足(n<", min_group_n, ")，仅显示: ", eligible_groups[1]) else NULL

        point_df <- if (isTRUE(include_points) && point_downsample_frac < 1.0) {
          dplyr::bind_rows(lapply(split(plot_data_eff, plot_data_eff$Group), get_group_sample, frac = point_downsample_frac, max_n = Inf))
        } else plot_data_eff
        smooth_df <- if ((smooth_sample_frac < 1.0) || is.finite(smooth_sample_max)) {
          dplyr::bind_rows(lapply(split(plot_data_eff, plot_data_eff$Group), get_group_sample, frac = smooth_sample_frac, max_n = smooth_sample_max))
        } else plot_data_eff

        if (panel == "percentage") {
          p <- ggplot2::ggplot(plot_data_eff, ggplot2::aes(x = Age2, y = IDP_Change_Percent, color = Group)) +
            { if (isTRUE(include_points)) ggplot2::geom_point(data = point_df, alpha = 0.5, size = 0.8) else ggplot2::geom_blank() } +
            ggplot2::geom_smooth(data = smooth_df, method = "loess", se = se_confint, span = smooth_span, alpha = 0.2, linewidth = 1.2) +
            ggplot2::scale_color_manual(values = group_colors, drop = FALSE) +
            ggplot2::labs(
              title = paste0(substr(variable_name, 1, 40)),
              subtitle = "Percentage Change vs Age",
              x = "Age (years)", y = if (percent_baseline == "baseline_mean") "ΔIDP (% of baseline mean)" else "ΔIDP (%)", color = "Group", caption = cap_text
            ) + theme_jama_refined() + ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1, size = 2)))
          filename_base <- paste0("Age_Trend_Percentage_", sprintf("%03d", i), "_", gsub("[^A-Za-z0-9]", "_", substr(variable_name, 1, 30)))
          save_pub(file.path(output_dir, filename_base), p, width = 6, height = 6, dpi = dpi)
        } else if (panel == "absolute") {
          p <- ggplot2::ggplot(plot_data_eff, ggplot2::aes(x = Age2, y = IDP_Change_Std, color = Group)) +
            { if (isTRUE(include_points)) ggplot2::geom_point(data = point_df, alpha = 0.5, size = 0.8) else ggplot2::geom_blank() } +
            ggplot2::geom_smooth(data = smooth_df, method = "loess", se = se_confint, span = smooth_span, alpha = 0.2, linewidth = 1.2) +
            ggplot2::scale_color_manual(values = group_colors, drop = FALSE) +
            ggplot2::labs(
              title = paste0(substr(variable_name, 1, 40)),
              subtitle = paste(
                "Absolute Change vs Age | Z =",
                round(suppressWarnings(as.numeric(if ("z_statistic_without_non_imaging" %in% names(idp_info)) idp_info$z_statistic_without_non_imaging else NA_real_)), 3)
              ),
              x = "Age (years)", y = if (delta_normalize_by == "baseline_sd") "ΔIDP (SD units)" else if (delta_normalize_by == "baseline_mean") "ΔIDP (per baseline mean)" else "Absolute Change", color = "Group", caption = cap_text
            ) + theme_jama_refined() + ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1, size = 2)))
          filename_base <- paste0("Age_Trend_Absolute_", sprintf("%03d", i), "_", gsub("[^A-Za-z0-9]", "_", substr(variable_name, 1, 30)))
          save_pub(file.path(output_dir, filename_base), p, width = 6, height = 6, dpi = dpi)
        } else {
          p1 <- ggplot2::ggplot(plot_data_eff, ggplot2::aes(x = Age2, y = IDP_Change_Std, color = Group)) +
            { if (isTRUE(include_points)) ggplot2::geom_point(data = point_df, alpha = 0.5, size = 0.8) else ggplot2::geom_blank() } +
            ggplot2::geom_smooth(data = smooth_df, method = "loess", se = se_confint, span = smooth_span, alpha = 0.2, linewidth = 1.2) +
            ggplot2::scale_color_manual(values = group_colors, drop = FALSE) +
            ggplot2::labs(
              title = paste(substr(variable_name, 1, 40)),
              subtitle = paste("Absolute Change vs Age | Z =", round(idp_info$z_statistic_without_non_imaging, 3)),
              x = "Age (years)", y = if (delta_normalize_by == "baseline_sd") "ΔIDP (SD units)" else if (delta_normalize_by == "baseline_mean") "ΔIDP (per baseline mean)" else "Absolute Change", color = "Group", caption = cap_text
            ) + theme_jama_refined() + ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1, size = 2)))

          p2 <- ggplot2::ggplot(plot_data_eff, ggplot2::aes(x = Age2, y = IDP_Change_Percent, color = Group)) +
            { if (isTRUE(include_points)) ggplot2::geom_point(data = point_df, alpha = 0.5, size = 0.8) else ggplot2::geom_blank() } +
            ggplot2::geom_smooth(data = smooth_df, method = "loess", se = se_confint, span = smooth_span, alpha = 0.2, linewidth = 1.2) +
            ggplot2::scale_color_manual(values = group_colors, drop = FALSE) +
            ggplot2::labs(
              title = paste(substr(variable_name, 1, 40)),
              subtitle = "Percentage Change vs Age",
              x = "Age (years)", y = if (percent_baseline == "baseline_mean") "ΔIDP (% of baseline mean)" else "ΔIDP (%)", color = "Group", caption = cap_text
            ) + theme_jama_refined() + ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1, size = 2)))

          combined_plot <- cowplot::plot_grid(p1, p2, ncol = 1, rel_heights = c(1, 1))
          filename_base <- paste0("Age_Trend_Both_", sprintf("%03d", i), "_", gsub("[^A-Za-z0-9]", "_", substr(variable_name, 1, 30)))
          save_pub(file.path(output_dir, filename_base), combined_plot, width = 10, height = 12, dpi = dpi)
        }
      }
    }
    out_rec
  }

  idx <- seq_len(nrow(top_idps))
  recs <- if (!is.null(map_opts)) mapper(idx, render_one, .options = map_opts) else mapper(idx, render_one)
  for (r in recs) { if (!is.null(r)) completeness_records[[length(completeness_records) + 1]] <- r }

  completeness_df <- if (length(completeness_records) > 0) dplyr::bind_rows(completeness_records) else NULL
  if (!is.null(completeness_df)) {
    utils::write.csv(completeness_df, file.path(output_dir, "Age_Trend_Completeness_Fast.csv"), row.names = FALSE)
    cat("已写入统一趋势完备性报表(FAST): ", file.path(output_dir, "Age_Trend_Completeness_Fast.csv"), "\n")
  }
  # 汇总原始个体级数据（可选）：合并本次生成的所有RawData文件
  if (isTRUE(export_raw_data)) {
    raw_files <- list.files(raw_dir, pattern = "^Age_Trend_RawData_.*\\.csv$", full.names = TRUE)
    if (length(raw_files) > 0) {
      raw_dfs <- lapply(raw_files, function(fp) tryCatch(utils::read.csv(fp, stringsAsFactors = FALSE), error = function(e) NULL))
      raw_dfs <- Filter(function(x) !is.null(x) && nrow(x) > 0, raw_dfs)
      if (length(raw_dfs) > 0) {
        raw_all <- dplyr::bind_rows(raw_dfs)
        utils::write.csv(raw_all, file.path(output_dir, "Age_Trend_RawData_AllModels.csv"), row.names = FALSE)
        cat("已合并导出个体级原始纵向数据: ", file.path(output_dir, "Age_Trend_RawData_AllModels.csv"), "\n")
      }
    }
  }
  return(completeness_df)
}

# =================== 5. MRI脑图映射函数 ===================
create_mri_brain_mapping <- function(brain_independent_idps, output_dir = "Brain_MRI_Maps") {
  load_brain_mapping_packages_once()
  
  cat("创建MRI脑图映射...\n")
  
  if(missing(brain_independent_idps) || nrow(brain_independent_idps) == 0) {
    cat("无脑影像独立效应IDP可供映射\n")
    return(NULL)
  }
  
  # 创建输出目录
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 解析IDP脑区信息
  cat("解析脑影像IDP信息...\n")
  z_vec <- if ("z_statistic_without_non_imaging" %in% names(brain_independent_idps)) {
    brain_independent_idps$z_statistic_without_non_imaging
  } else if ("z_change" %in% names(brain_independent_idps)) {
    brain_independent_idps$z_change
  } else if ("z" %in% names(brain_independent_idps)) {
    brain_independent_idps$z
  } else {
    rep(NA_real_, nrow(brain_independent_idps))
  }
  brain_parsed <- parse_ukb_idp_comprehensive(brain_independent_idps$variable_name, z_vec)
  
  # 合并数据
  # 使用列绑定并自动修复重复列名，避免后续filter/mutate因重复列名报错
  if (anyDuplicated(c(names(brain_independent_idps), names(brain_parsed))) > 0) {
    cat("检测到重复列名，正在修复以避免filter报错...\n")
  }
  brain_mapping_data <- dplyr::bind_cols(
      tibble::as_tibble(brain_independent_idps),
      tibble::as_tibble(brain_parsed),
      .name_repair = "unique"
    ) %>%
    dplyr::filter(structure_type %in% c("Cortical", "Subcortical", "White Matter")) %>%
    dplyr::mutate(network = parse_idp_network(variable_name))
  z_src <- if ("z_statistic_without_non_imaging" %in% names(brain_mapping_data)) {
    brain_mapping_data$z_statistic_without_non_imaging
  } else if ("z_change" %in% names(brain_mapping_data)) {
    brain_mapping_data$z_change
  } else if ("z" %in% names(brain_mapping_data)) {
    brain_mapping_data$z
  } else {
    rep(NA_real_, nrow(brain_mapping_data))
  }
  brain_mapping_data$z_change <- suppressWarnings(as.numeric(z_src))
  brain_mapping_data <- brain_mapping_data %>% dplyr::arrange(dplyr::desc(abs(z_change)))

  # 参数化映射逻辑：通过选项或环境变量控制阈值与取值方式
  # BRAIN_MAP_MODE: 'abs_z' 或 'signed_z'（默认 'abs_z'）
  # BRAIN_MAP_Z_THRESHOLD: 数值阈值（默认 0；<=0 表示不过滤）
  # BRAIN_MAP_MIN_SELECTED: 阈值筛选后若选中太少，触发回退（默认 25）
  # BRAIN_MAP_TOP_N: 回退时按 |Z| 选择 Top N（默认 60）
  map_mode <- tolower(trimws(as.character(getOption("brain.map.mode", Sys.getenv("BRAIN_MAP_MODE", "abs_z")))))
  if (!map_mode %in% c("abs_z","signed_z")) map_mode <- "abs_z"
  z_thresh <- suppressWarnings(as.numeric(getOption("brain.map.z_threshold", Sys.getenv("BRAIN_MAP_Z_THRESHOLD", "1.96"))))
  if (!is.finite(z_thresh)) z_thresh <- 0
  min_selected <- suppressWarnings(as.integer(getOption("brain.map.min_selected", Sys.getenv("BRAIN_MAP_MIN_SELECTED", "25"))))
  if (!is.finite(min_selected) || min_selected < 1) min_selected <- 25L
  top_n <- suppressWarnings(as.integer(getOption("brain.map.top_n", Sys.getenv("BRAIN_MAP_TOP_N", "60"))))
  if (!is.finite(top_n) || top_n < 1) top_n <- 60L
  # 计算效应列
  brain_mapping_data$effect_magnitude <- if (map_mode == "abs_z") abs(brain_mapping_data$z_change) else brain_mapping_data$z_change
  # 标记选中行
  sel <- if (z_thresh > 0) {
    abs(brain_mapping_data$z_change) > z_thresh
  } else {
    rep(TRUE, nrow(brain_mapping_data))
  }
  sel[is.na(sel)] <- FALSE
  if (sum(sel) < min_selected) {
    ord <- order(abs(brain_mapping_data$z_change), decreasing = TRUE, na.last = NA)
    k <- min(top_n, length(ord))
    sel <- rep(FALSE, nrow(brain_mapping_data))
    if (k > 0) sel[ord[seq_len(k)]] <- TRUE
  }
  brain_mapping_data$selected_for_mapping <- sel

  # 从根源统一字段：稳健合并 hemisphere 与 measure_type，避免 name_repair 后的 ...XX 后缀导致找不到分组列
  unify_first_nonempty <- function(df, base) {
    cols <- grep(paste0("^", base), names(df), value = TRUE)
    if (length(cols) == 0) {
      df[[base]] <- NA
      return(df)
    }
    out <- df[[cols[1]]]
    if (length(cols) > 1) {
      for (col in cols[-1]) {
        out <- dplyr::coalesce(out, df[[col]])
      }
    }
    df[[base]] <- out
    df
  }
  brain_mapping_data <- unify_first_nonempty(brain_mapping_data, "hemisphere")
  brain_mapping_data <- unify_first_nonempty(brain_mapping_data, "measure_type")
  # 预防类型警告：确保关键文本列为字符类型
  cols_to_char <- intersect(c("ggseg_region","brain_region","anatomical_structure","measure_type","hemisphere"), names(brain_mapping_data))
  if (length(cols_to_char) > 0) {
    brain_mapping_data[cols_to_char] <- lapply(brain_mapping_data[cols_to_char], function(x) as.character(x))
  }
  # 统一并标准化 measure_type 与 ROI 分组
  normalize_measure_type <- function(s, structure_type) {
    s_low <- tolower(trimws(as.character(s)))
    st <- tolower(trimws(as.character(structure_type)))
    dplyr::case_when(
      grepl("surface\\.?\u00A0?area|cortical\\s*surface", s_low) ~ "Cortical Surface Area",
      grepl("thick|cortical\\s*thick", s_low) ~ "Cortical Thickness",
      grepl("volume", s_low) & st == "cortical" ~ "Cortical Volume",
      grepl("volume", s_low) & st == "subcortical" ~ "Subcortical Volume",
      grepl("\\bfa\\b|fractional\\s*anisotropy", s_low) ~ "WM Tract FA",
      grepl("\\bmd\\b|mean\\s*diffusiv", s_low) ~ "WM Tract MD",
      grepl("cbf|perfusion|asl", s_low) ~ "Perfusion (CBF)",
      TRUE ~ ifelse(st == "white matter", "White Matter Metric", "Other")
    )
  }
  classify_roi_group <- function(structure_type, network) {
    st <- tolower(trimws(as.character(structure_type)))
    net <- tolower(trimws(as.character(network)))
    dplyr::case_when(
      st == "cortical" ~ "Cortical",
      st == "subcortical" ~ "Subcortical",
      st == "white matter" ~ "White Matter",
      !is.na(net) & net != "" & net != "unknown" ~ "Network",
      TRUE ~ "Unknown"
    )
  }
  # 规范化半球取值到 left/right/bilateral，并派生稳定绘图列 hemi
  brain_mapping_data <- brain_mapping_data %>%
    dplyr::mutate(hemisphere = tolower(trimws(as.character(hemisphere)))) %>%
    dplyr::mutate(hemisphere = dplyr::case_when(
      hemisphere %in% c("l", "left", "lh") ~ "left",
      hemisphere %in% c("r", "right", "rh") ~ "right",
      hemisphere %in% c("b", "bilateral", "both") ~ "bilateral",
      grepl("left", hemisphere) ~ "left",
      grepl("right", hemisphere) ~ "right",
      TRUE ~ "bilateral"
    )) %>%
    dplyr::mutate(hemi = dplyr::case_when(
      hemisphere %in% c("left", "bilateral") ~ "left",
      hemisphere %in% c("right") ~ "right",
      TRUE ~ "left"
    )) %>%
    dplyr::mutate(
      measure_type_std = normalize_measure_type(measure_type, structure_type),
      roi_group = classify_roi_group(structure_type, network)
    )

  # 结构类型纠偏：当 measure_type 为 Subcortical/Regional 且区域名命中皮层下关键词时，将 structure_type 归一为 Subcortical；
  # 当 measure_type 指示白质度量时归一为 White Matter；其余保持原值。
  subcortical_kw <- "(hippoc|amygdal|thalam|caudate|putamen|pallid|globus|accumbens)"
  cortical_hint_kw <- "(insula|cingulate|frontal|parietal|occipital|temporal|precuneus|fusiform|entorhinal|parahipp|lingual|cuneus|supramarginal|angular)"
  # 映射类别归类：用于稳健驱动各模块的输入筛选
  classify_map_category <- function(structure_type, measure_type_std, region_text, network) {
    st <- tolower(trimws(as.character(structure_type)))
    mt <- tolower(trimws(as.character(measure_type_std)))
    r  <- tolower(trimws(as.character(region_text)))
    net <- tolower(trimws(as.character(network)))
    is_vascular <- grepl("(anterior\\.?cingulate|rostral\\.?anterior\\.?cingulate|caudal\\.?anterior\\.?cingulate|paracentral|medial\\.?frontal|superior\\.?frontal|frontal\\.?medial|rostral\\.?middle\\.?frontal|insula|precentral|postcentral|middle\\.?temporal|superior\\.?temporal|inferior\\.?temporal|supramarginal|angular|inferior\\.?parietal|parsopercularis|parstriangularis|pars\\.?orbitalis|middle\\.?frontal|bank\\.?sts|occipital|calcarine|pericalcarine|cuneus|lingual|lateral\\.?occipital|fusiform|hippocamp|parahip|entorhinal|precuneus|posterior\\.?cingulate|cerebell|brainstem|pons|medulla|vermis)", r)
    dplyr::case_when(
      grepl("(hippocamp|ca1|ca2|ca3|ca4|dentate|dg|subicul|presubicul)", r) ~ "Hippocampus",
      st == "subcortical" & !grepl("hippocamp", r) ~ "Subcortical",
      st == "white matter" | grepl("\\bwm\\b|tract|fa|md", mt) ~ "White Matter",
      !is.na(net) & net != "" & net != "unknown" ~ "Network",
      is_vascular ~ "Vascular",
      st == "cortical" ~ "Cortical",
      TRUE ~ "Unknown"
    )
  }
  brain_mapping_data <- brain_mapping_data %>%
    dplyr::mutate(region_src_low = tolower(trimws(as.character(dplyr::coalesce(ggseg_region, brain_region, anatomical_structure))))) %>%
    dplyr::mutate(structure_type = dplyr::case_when(
      # 白质度量统一到 White Matter
      tolower(measure_type) %in% c("wm tract", "dti skeleton") ~ "White Matter",
      grepl("white\\s*matter|\\bwm\\b|tract|\\bfa\\b|\\bmd\\b", tolower(measure_type_std)) ~ "White Matter",
      # 皮层下区域关键词 + Subcortical/Regional 则归 Subcortical
      grepl(subcortical_kw, region_src_low) & tolower(measure_type) %in% c("subcortical", "subcortical/regional") ~ "Subcortical",
      # 其余 Subcortical/Regional 若明显皮层提示则归 Cortical
      tolower(measure_type) %in% c("subcortical", "subcortical/regional") & !grepl(subcortical_kw, region_src_low) & grepl(cortical_hint_kw, region_src_low) ~ "Cortical",
      TRUE ~ structure_type
    )) %>%
    # Hippocampus 独立分类（包括亚区；覆盖任何误归入 Subcortical 的情况）
    dplyr::mutate(structure_type = dplyr::if_else(
      grepl("(hippocamp|ca1|ca2|ca3|ca4|dentate|dg|subicul|presubicul)", region_src_low),
      "Hippocampus", structure_type
    )) %>%
    dplyr::mutate(roi_group = classify_roi_group(structure_type, network)) %>%
    dplyr::mutate(region_for_cat = tolower(trimws(as.character(dplyr::coalesce(ggseg_region, brain_region, anatomical_structure))))) %>%
    dplyr::mutate(map_category = classify_map_category(structure_type, measure_type_std, region_for_cat, network)) %>%
    dplyr::select(-region_src_low, -region_for_cat)
  
  # Cortical ggseg_region 回填：若缺失则根据通用脑区名映射到DK atlas可识别的region
  fallback_region_map <- c(
    "insula" = "insula",
    "cingulate" = "rostral anterior cingulate",
    "frontal lobe" = "superior frontal",
    "temporal lobe" = "middle temporal",
    "parietal lobe" = "superior parietal",
    "occipital lobe" = "lateral occipital",
    "precuneus" = "precuneus",
    "lingual" = "lingual",
    "cuneus" = "cuneus",
    "fusiform" = "fusiform",
    "parahippocampal" = "parahippocampal",
    "entorhinal" = "entorhinal",
    "banks sts" = "bankssts"
  )
  brain_mapping_data <- brain_mapping_data %>%
    dplyr::mutate(src_text = dplyr::coalesce(ggseg_region, brain_region, anatomical_structure)) %>%
    dplyr::mutate(brain_key = {
      k <- tolower(trimws(as.character(src_text))); k <- gsub("[ ]+", " ", k); k
    }) %>%
    dplyr::mutate(alt_dk = normalize_dk_region(brain_key)) %>%
    dplyr::mutate(
      ggseg_region = dplyr::if_else(
        structure_type == "Cortical" & (is.na(ggseg_region) | ggseg_region == ""),
        dplyr::coalesce(alt_dk, unname(fallback_region_map[brain_key]), ggseg_region),
        ggseg_region
      )
    ) %>%
    dplyr::select(-brain_key, -src_text, -alt_dk)

  if (requireNamespace("ggseg", quietly = TRUE)) {
    if (!exists("dk")) suppressWarnings(data(dk, package = "ggseg"))
    if (exists("dk")) {
      dk_obj <- get("dk", inherits = TRUE)
      dk_df <- if (is.list(dk_obj) && "data" %in% names(dk_obj)) dk_obj$data else dk_obj
      valid_dk_regions <- if (is.data.frame(dk_df) && "region" %in% names(dk_df)) unique(as.character(dk_df$region)) else character(0)
      valid_dk_regions <- valid_dk_regions[!is.na(valid_dk_regions) & nzchar(valid_dk_regions)]
      norm_key <- function(x) {
        x <- tolower(trimws(as.character(x)))
        x <- gsub("[^a-z0-9]+", "", x)
        x
      }
      dk_key <- norm_key(valid_dk_regions)
      dk_map <- stats::setNames(valid_dk_regions, dk_key)
      canon_dk_region <- function(x) {
        if (is.null(x) || is.na(x) || !nzchar(x)) return(NA_character_)
        if (x %in% valid_dk_regions) return(x)
        k <- norm_key(x)
        if (k %in% names(dk_map)) return(unname(dk_map[[k]]))
        NA_character_
      }
      infer_dk_region <- function(txt) {
        if (is.null(txt) || is.na(txt) || !nzchar(txt)) return(NA_character_)
        n1 <- normalize_dk_region(txt)
        c1 <- canon_dk_region(n1)
        if (!is.na(c1)) return(c1)
        c0 <- canon_dk_region(txt)
        if (!is.na(c0)) return(c0)
        t <- tolower(trimws(as.character(txt)))
        if (grepl("\\bcollat\\b|collateral", t)) return("parahippocampal")
        if (grepl("suborbital|orbital", t)) return("medial orbitofrontal")
        if (grepl("lunatus|occip", t) && !grepl("calcarine", t)) return("lateral occipital")
        if (grepl("calcarine", t)) return("pericalcarine")
        if (grepl("lat\\s*fis|lateral\\s*fiss|sylv", t)) return("superior temporal")
        if (grepl("\\bfront\\b|frontal", t) && grepl("\\binf\\b|inferior", t)) return("pars triangularis")
        if (grepl("\\btemp\\b|temporal", t) && grepl("\\bsup\\b|superior", t)) return("superior temporal")
        if (grepl("\\btemp\\b|temporal", t) && grepl("\\bmid\\b|middle", t)) return("middle temporal")
        if (grepl("\\btemp\\b|temporal", t) && grepl("\\binf\\b|inferior", t)) return("inferior temporal")
        if (grepl("\\bpariet\\b|parietal", t) && grepl("\\bsup\\b|superior", t)) return("superior parietal")
        if (grepl("\\bpariet\\b|parietal", t) && grepl("\\binf\\b|inferior", t)) return("inferior parietal")
        tk <- norm_key(t)
        if (!nzchar(tk)) return(NA_character_)
        hit <- which(nzchar(dk_key) & nchar(dk_key) >= 5 & vapply(dk_key, function(p) grepl(p, tk, fixed = TRUE), logical(1)))
        if (length(hit) > 0) {
          best <- hit[which.max(nchar(dk_key[hit]))]
          return(valid_dk_regions[best])
        }
        dmat <- utils::adist(tk, dk_key, partial = TRUE, ignore.case = TRUE)
        d <- if (is.matrix(dmat)) as.numeric(dmat[1, ]) else as.numeric(dmat)
        if (length(d) == 0) return(NA_character_)
        best <- which.min(d)
        if (is.finite(d[best]) && d[best] <= 3) return(valid_dk_regions[best])
        NA_character_
      }
      cortical_idx <- which(brain_mapping_data$structure_type == "Cortical")
      if (length(cortical_idx) > 0) {
        cur <- as.character(brain_mapping_data$ggseg_region[cortical_idx])
        cur2 <- vapply(cur, canon_dk_region, character(1))
        bad <- is.na(cur2) | !nzchar(cur2) | !(cur2 %in% valid_dk_regions)
        if (any(bad)) {
          txt <- dplyr::coalesce(
            brain_mapping_data$display_name[cortical_idx],
            brain_mapping_data$anatomical_structure[cortical_idx],
            brain_mapping_data$brain_region[cortical_idx],
            brain_mapping_data$variable_name[cortical_idx]
          )
          newv <- vapply(txt, infer_dk_region, character(1))
          cur2[bad] <- newv[bad]
          brain_mapping_data$ggseg_region[cortical_idx] <- cur2
        }
      }
    }
  }

  # 统一最终脑区列供下游汇总使用
  brain_mapping_data <- brain_mapping_data %>%
    dplyr::mutate(yeo_region = align_to_yeo7_region(network)) %>%
    dplyr::mutate(wm_tract = dplyr::if_else(
      structure_type == "White Matter",
      normalize_wm_tract(dplyr::coalesce(anatomical_structure, brain_region, ggseg_region)),
      NA_character_
    )) %>%
    dplyr::mutate(brain_region_final = dplyr::coalesce(ggseg_region, brain_region, wm_tract, anatomical_structure, yeo_region))

  # 基本一致性检查与日志输出
  missing_region_rate <- mean(is.na(brain_mapping_data$brain_region_final) | brain_mapping_data$brain_region_final == "")
  cat(sprintf("脑区统一列缺失比例: %.2f%%\n", 100 * missing_region_rate))

  # 精简重复后缀列：若存在同名基础列则移除带 ...NN 后缀列
  drop_duplicate_suffix_columns <- function(df) {
    nm <- names(df)
    suffix_idx <- grep("\\.\\.\\.[0-9]+$", nm)
    if (length(suffix_idx) == 0) return(df)
    keep <- rep(TRUE, length(nm))
    for (i in suffix_idx) {
      base <- sub("\\.\\.\\.[0-9]+$", "", nm[i])
      if (base %in% nm) {
        keep[i] <- FALSE
      }
    }
    df[, keep, drop = FALSE]
  }
  brain_mapping_data <- drop_duplicate_suffix_columns(brain_mapping_data)
  
  # 保存解析结果
  write.csv(brain_mapping_data, file.path(output_dir, "Brain_IDPs_Parsed_for_Mapping.csv"), row.names = FALSE)
  # 兼容旧版Python脚本：在项目根目录生成 All_Independent_Effect_IDPs_for_Visualization.csv
  # 为避免重复存储，皮质脑图的保存移至专用子目录函数（create_cortical_mapping）
  # 此处不再在根目录直接输出皮质图像
  
  # 皮层下（Subcortical）根目录保存已移除，改用海马相关映射专用目录（create_hippocampal_mapping）
  
  # 返回结果对象
  return(list(brain_mapping_data = brain_mapping_data))
}

effect_dir_from_z <- function(z) ifelse(mean(z, na.rm = TRUE) > 0, "Increase", "Decrease")

# =================== 5.1 皮质与皮质下分别映射（独立函数） ===================
create_cortical_mapping <- function(brain_mapping_data, output_dir = file.path("Brain_MRI_Maps", "Cortical")) {
  load_brain_mapping_packages_once()
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if(!requireNamespace("ggseg", quietly = TRUE)) return(NULL)
  if(!exists("dk")) {
    if (requireNamespace("ggseg", quietly = TRUE)) {
      suppressWarnings(data(dk, package = "ggseg"))
    }
  }
  dk_obj <- if (exists("dk", inherits = TRUE)) get("dk", inherits = TRUE) else NULL
  dk_df <- if (!is.null(dk_obj) && is.list(dk_obj) && "data" %in% names(dk_obj)) dk_obj$data else dk_obj
  valid_regions <- if (is.data.frame(dk_df) && "region" %in% names(dk_df)) unique(as.character(dk_df$region)) else character(0)
  valid_regions <- valid_regions[!is.na(valid_regions) & nzchar(valid_regions)]
  norm_key <- function(x) {
    x <- tolower(trimws(as.character(x)))
    x <- gsub("[^a-z0-9]+", "", x)
    x
  }
  dk_key <- norm_key(valid_regions)
  dk_map <- stats::setNames(valid_regions, dk_key)
  canon_region <- function(x) {
    if (is.null(x) || is.na(x) || !nzchar(x)) return(NA_character_)
    if (x %in% valid_regions) return(x)
    k <- norm_key(x)
    if (k %in% names(dk_map)) return(unname(dk_map[[k]]))
    NA_character_
  }
  
  cortical_data <- brain_mapping_data %>%
    dplyr::filter(.data[["selected_for_mapping"]]) %>%
    dplyr::filter(structure_type == "Cortical", ggseg_region != "", !is.na(ggseg_region)) %>%
    dplyr::mutate(ggseg_region = vapply(as.character(ggseg_region), canon_region, character(1))) %>%
    dplyr::filter(!is.na(ggseg_region) & ggseg_region != "") %>%
    dplyr::group_by(ggseg_region, hemi) %>%
    dplyr::summarise(
      avg_z_change = mean(abs(z_change), na.rm = TRUE),
      max_z_change = max(abs(z_change), na.rm = TRUE),
      count = dplyr::n(),
      effect_direction = effect_dir_from_z(z_change),
      .groups = "drop"
    )
  
  if(nrow(cortical_data) == 0) return(NULL)
  
  left_data <- cortical_data %>% 
    dplyr::filter(hemi == "left") %>%
    dplyr::select(region = ggseg_region, z_value = max_z_change, direction = effect_direction) %>%
    dplyr::mutate(hemi = "left")
  
  right_data <- cortical_data %>% 
    dplyr::filter(hemi == "right") %>%
    dplyr::select(region = ggseg_region, z_value = max_z_change, direction = effect_direction) %>%
    dplyr::mutate(hemi = "right")
  
  combined_ggseg_data <- dplyr::bind_rows(left_data, right_data)
  
  p_brain <- combined_ggseg_data %>%
    ggplot2::ggplot(ggplot2::aes(fill = z_value)) +
    ggseg::geom_brain(atlas = dk, position = ggseg::position_brain(hemi ~ side)) +
    { safe_z_max <- suppressWarnings(max(combined_ggseg_data$z_value, na.rm = TRUE)); if (!is.finite(safe_z_max) || safe_z_max <= 0) safe_z_max <- 1e-6; viridis::scale_fill_viridis(option = "C", direction = 1, name = "|Z|", limits = c(0, safe_z_max)) } +
    ggplot2::labs(
      title = "Cortical IDPs Mapping",
      subtitle = paste("Max |Z| per region | n:", nrow(cortical_data)),
      caption = "DK atlas (ggseg)"
    ) +
    theme_brain_pub()
  
  save_pub(file.path(output_dir, "Cortical_IDPs_Mapping"), p_brain, width = 12, height = 8, dpi = 300)
  invisible(p_brain)
}

# Subcortical and WhiteMatterVoxel mapping removed per user request (duplicate with Network fiber bundles)


# 海马相关区域映射（强调认知/记忆脑区）
create_hippocampal_mapping <- function(brain_mapping_data, output_dir = file.path("Brain_MRI_Maps", "Hippocampus")) {
  load_brain_mapping_packages_once()
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  hm_regions <- brain_mapping_data %>%
    dplyr::filter(.data[["selected_for_mapping"]]) %>%
    dplyr::mutate(region_src = dplyr::coalesce(ggseg_region, brain_region),
                  region_src_lower = tolower(region_src)) %>%
    dplyr::filter(!is.na(region_src_lower) & region_src_lower != "") %>%
    dplyr::filter(
      grepl("hippoc", region_src_lower) |
      grepl("parahipp", region_src_lower) |
      grepl("entorh", region_src_lower) |
      grepl("subicul", region_src_lower) |
      grepl("dentate", region_src_lower)
    ) %>%
    dplyr::group_by(region_src, hemi) %>%
    dplyr::summarise(
      z_value = max(abs(z_change), na.rm = TRUE),
      effect_direction = effect_dir_from_z(z_change),
      .groups = "drop"
    )
  if(nrow(hm_regions) == 0) return(NULL)

  # 优先尝试：ggseg + aseg 映射海马整体
  p_map <- NULL
  if (requireNamespace("ggseg", quietly = TRUE)) {
    if(!exists("aseg")) {
      if (requireNamespace("ggseg", quietly = TRUE)) {
        suppressWarnings(data(aseg, package = "ggseg"))
      }
    }
    hipp <- hm_regions %>% dplyr::filter(grepl("hippoc", tolower(region_src)))
    if(nrow(hipp) > 0 && exists("aseg")) {
      hipp <- hipp %>% dplyr::mutate(region = "hippocampus")
      safe_z_max_h <- suppressWarnings(max(hipp$z_value, na.rm = TRUE))
      if (!is.finite(safe_z_max_h) || safe_z_max_h <= 0) safe_z_max_h <- 1e-6
      # 部分 atlas（如 aseg）不包含 side 列，使用默认位置更稳健
      p_map <- ggplot2::ggplot(hipp, ggplot2::aes(fill = z_value)) +
        ggseg::geom_brain(atlas = aseg) +
        viridis::scale_fill_viridis(option = "C", direction = 1, name = "|Z|",
                                    limits = c(0, safe_z_max_h)) +
        ggplot2::labs(title = "Hippocampus Mapping (aseg)", subtitle = paste("Max |Z| per hemisphere | n:", nrow(hipp))) +
        theme_brain_pub()
      save_pub(file.path(output_dir, "Hippocampal_Atlas_Map"), p_map, width = 10, height = 6, dpi = 600)
    }
  }

  # 其次：Harvard–Oxford 半球聚合展示
  p_ho <- NULL
  if (requireNamespace("ggsegHO", quietly = TRUE) && requireNamespace("ggseg", quietly = TRUE)) {
    suppressWarnings(data(ho, package = "ggsegHO"))
    hemi_avg <- hm_regions %>% dplyr::group_by(hemisphere = hemi) %>% dplyr::summarise(stat_avg = mean(z_value, na.rm = TRUE), .groups = "drop")
    ho_df <- data.frame(region = "hippocampus", hemi = hemi_avg$hemisphere, z_value = hemi_avg$stat_avg)
    if(nrow(ho_df) > 0 && exists("ho")) {
      p_ho <- ggplot2::ggplot(ho_df, ggplot2::aes(fill = z_value)) +
        ggseg::geom_brain(atlas = ho, position = ggseg::position_brain(hemi ~ side)) +
        ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, name = "Effect Size") +
        ggplot2::labs(title = "Hippocampus (Harvard-Oxford)", subtitle = "Hemisphere-averaged across hippocampal-related regions") +
        theme_brain_pub()
      save_pub(file.path(output_dir, "Hippocampus_HO_Aggregate_Map"), p_ho, width = 10, height = 6, dpi = 600)
    }
  }

  # 回退：若上述脑图均未生成，则绘制条形图
  p_bar <- NULL
  if (is.null(p_map) && is.null(p_ho)) {
    p_bar <- hm_regions %>%
      dplyr::mutate(region_hemi = paste0(region_src, " (", hemi, ")")) %>%
      ggplot2::ggplot(ggplot2::aes(x = reorder(region_hemi, z_value), y = z_value, fill = z_value)) +
      ggplot2::geom_col() + ggplot2::coord_flip() +
      viridis::scale_fill_viridis(option = "C", direction = 1, name = "|Z|") +
      ggplot2::labs(title = "Hippocampal & MTL IDPs", subtitle = "Max |Z| per region and hemisphere") +
      theme_brain_pub()
    save_pub(file.path(output_dir, "Hippocampal_Memory_Regions_Barplot"), p_bar, width = 10, height = 6, dpi = 600)
  }

  utils::write.csv(hm_regions, file.path(output_dir, "Hippocampal_Regions_Summary.csv"), row.names = FALSE)
  return(list(plot_brain_aseg = p_map, plot_ho = p_ho, plot_bar = p_bar, data = hm_regions))
}

# =================== 5.3 功能网络映射（Yeo 7-network 优先，回退条形图） ===================
create_network_mapping <- function(brain_mapping_data, output_dir = file.path("Brain_MRI_Maps", "Network")) {
  load_brain_mapping_packages_once()
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  net_data <- brain_mapping_data %>%
    dplyr::filter(.data[["selected_for_mapping"]]) %>%
    dplyr::mutate(network_label = parse_idp_network(variable_name)) %>%
    dplyr::filter(!is.na(network_label) & network_label != "Unknown") %>%
    dplyr::mutate(yeo_region = align_to_yeo7_region(network_label)) %>%
    dplyr::filter(!is.na(yeo_region)) %>%
    dplyr::group_by(network_label, hemi) %>%
    dplyr::summarise(
      z_value = max(abs(z_change), na.rm = TRUE),
      count = dplyr::n(),
      effect_direction = effect_dir_from_z(z_change),
      .groups = "drop"
    )

  if(nrow(net_data) == 0) return(NULL)

  p_net_bar <- net_data %>%
    dplyr::mutate(network_hemi = paste0(network_label, " (", hemi, ")")) %>%
    ggplot2::ggplot(ggplot2::aes(x = reorder(network_hemi, z_value), y = z_value, fill = network_label)) +
    ggplot2::geom_col() + ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = network_palette) +
    ggplot2::labs(title = "Yeo 7-Network Mapping (summary)", subtitle = "Max |Z| per network and hemisphere") +
    theme_brain_pub()
  save_pub(file.path(output_dir, "Brain_IDPs_MRI_Mapping_By_Network"), p_net_bar, width = 10, height = 6, dpi = 600)

  if (requireNamespace("ggsegYeo2011", quietly = TRUE) && requireNamespace("ggseg", quietly = TRUE)) {
    try({
      suppressWarnings(data(yeo7, package = "ggsegYeo2011"))
      gg_net <- net_data %>%
        dplyr::mutate(region = align_to_yeo7_region(network_label)) %>%
        dplyr::filter(!is.na(region))
      p_net_map <- ggplot2::ggplot(ggplot2::aes(fill = z_value), data = gg_net) +
        ggseg::geom_brain(atlas = yeo7, position = ggseg::position_brain(hemi ~ side)) +
        viridis::scale_fill_viridis(option = "C", direction = 1, name = "|Z|") +
        ggplot2::labs(title = "Yeo 7-Network Brain Map", subtitle = paste("n:", nrow(gg_net))) +
        theme_brain_pub()
      save_pub(file.path(output_dir, "Brain_IDPs_MRI_Mapping_By_Network_Brain"), p_net_map, width = 12, height = 8, dpi = 300)
    }, silent = TRUE)
  } else if (requireNamespace("ggsegYeo", quietly = TRUE)) {
    try({
      suppressWarnings(data(yeo7, package = "ggsegYeo"))
      gg_net <- net_data %>% dplyr::rename(network = network_label)
      p_net_map <- ggplot2::ggplot(ggplot2::aes(fill = z_value), data = gg_net) +
        ggseg::geom_brain(atlas = yeo7, position = ggseg::position_brain(hemi ~ side)) +
        viridis::scale_fill_viridis(option = "C", direction = 1, name = "|Z|") +
        ggplot2::labs(title = "Yeo 7-Network Brain Map", subtitle = paste("n:", nrow(gg_net))) +
        theme_brain_pub()
      save_pub(file.path(output_dir, "Brain_IDPs_MRI_Mapping_By_Network_Brain"), p_net_map, width = 12, height = 8, dpi = 300)
    }, silent = TRUE)
  }

  invisible(net_data)
}

# 直接从结果表生成 Yeo7 网络级条形图与脑图

# =================== 5.2.1 海马亚区（CA1-4/DG/Subiculum）精细映射 ===================
visualize_hippocampal_subfields <- function(results, stat_column = "cohens_d", output_dir = file.path("Brain_MRI_Maps", "Hippocampus")) {
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  cat("生成海马亚区映射...\n")

  # 选择变量名列以解析（支持 idp_pair / variable_name / IDP）
  var_col <- intersect(c("idp_pair", "variable_name", "IDP"), names(results))
  if(length(var_col) == 0) {
    warning("结果表缺少可解析的变量名列（idp_pair/variable_name/IDP）。")
    return(NULL)
  }
  var_col <- var_col[1]

  # 过滤海马相关条目
  hippo_results <- results %>%
    dplyr::filter(stringr::str_detect(.data[[var_col]], "(?i)hippocamp|ca[1-4]|dentate|subiculum|hippocampal\\.tail|hippocampal\\.fissure"))

  if(nrow(hippo_results) == 0) {
    cat("未找到海马亚区数据\n")
    return(NULL)
  }

  # 选择统计列；若不存在则尝试备用列
  stat_col <- stat_column
  if(!(stat_col %in% names(hippo_results))) {
    fallbacks <- c("z_statistic_without_non_imaging", "z_statistic_with_non_imaging", "effect_size")
    fb <- intersect(fallbacks, names(hippo_results))
    if(length(fb) == 0) stop("未找到效应量/统计列，请指定 stat_column。") else stat_col <- fb[1]
  }

  # 解析半球与亚区标签
  parse_hemi <- function(s) {
    v <- tolower(s)
    dplyr::case_when(
      stringr::str_detect(v, "\\b(left|lh|l\\b|\\.L)\\b") ~ "left",
      stringr::str_detect(v, "\\b(right|rh|r\\b|\\.R)\\b") ~ "right",
      TRUE ~ "left" # 默认留左，便于绘图
    )
  }

  parse_subfield <- function(s) {
    v <- tolower(s)
    dplyr::case_when(
      stringr::str_detect(v, "ca1") ~ "CA1",
      stringr::str_detect(v, "ca2|ca3") ~ "CA2-3",
      stringr::str_detect(v, "ca4|dg|dentate") ~ "CA4/Dentate Gyrus",
      stringr::str_detect(v, "subiculum") ~ "Subiculum",
      stringr::str_detect(v, "hippocampal\\.tail|hippo.*tail") ~ "Hippocampal Tail",
      stringr::str_detect(v, "fissure") ~ "Hippocampal Fissure",
      TRUE ~ NA_character_
    )
  }

  hippo_data <- hippo_results %>%
    dplyr::mutate(raw = .data[[var_col]]) %>%
    dplyr::mutate(hemisphere = parse_hemi(raw), subfield = parse_subfield(raw)) %>%
    dplyr::filter(!is.na(subfield)) %>%
    dplyr::mutate(stat_value = suppressWarnings(as.numeric(.data[[stat_col]])))

  if(nrow(hippo_data) == 0) {
    warning("未能解析到具体海马亚区标签。")
    return(NULL)
  }

  # 亚区角色说明
  roles <- data.frame(
    subfield = c("CA1", "CA2-3", "CA4/Dentate Gyrus", "Subiculum", "Hippocampal Tail", "Hippocampal Fissure"),
    functional_role = c(
      "Pattern separation",
      "Memory consolidation & associative learning",
      "Neurogenesis & pattern completion",
      "Memory output & cortical relay",
      "Spatial memory (posterior extension)",
      "Structural landmark"
    ),
    stringsAsFactors = FALSE
  )

  hippo_summary <- hippo_data %>%
    dplyr::group_by(subfield, hemisphere) %>%
    dplyr::summarise(
      mean_effect = mean(stat_value, na.rm = TRUE),
      se = stats::sd(stat_value, na.rm = TRUE) / sqrt(dplyr::n()),
      n = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::left_join(roles, by = c("subfield" = "subfield"))

  # 条形图
  p_hippo <- ggplot2::ggplot(hippo_summary, ggplot2::aes(x = reorder(subfield, mean_effect), y = mean_effect)) +
    ggplot2::geom_col(ggplot2::aes(fill = hemisphere), position = ggplot2::position_dodge(width = 0.8), color = "black", width = 0.7) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = mean_effect - se, ymax = mean_effect + se, group = hemisphere), position = ggplot2::position_dodge(width = 0.8), width = 0.3, linewidth = 0.8) +
    ggplot2::scale_fill_manual(values = c("left" = "#1f77b4", "right" = "#ff7f0e"), name = "Hemisphere") +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "Hippocampal Subfield", y = sprintf("Effect Size (%s)", stat_col), title = "Hippocampal Subfield Changes", subtitle = "FreeSurfer subfields (aggregated)") +
    theme_brain_pub()
  save_pub(file.path(output_dir, "Hippocampal_Subfields_Effects"), p_hippo, width = 10, height = 6, dpi = 600)

  # 若可用，使用 Harvard-Oxford 映射海马整体（按半球汇总效果作为对照图）
  p_ho <- NULL
  if (requireNamespace("ggsegHO", quietly = TRUE) && requireNamespace("ggseg", quietly = TRUE)) {
    suppressWarnings(data(ho, package = "ggsegHO"))
    hemi_avg <- hippo_summary %>% dplyr::group_by(hemisphere) %>% dplyr::summarise(stat_avg = mean(mean_effect, na.rm = TRUE), .groups = "drop")
    ho_df <- data.frame(region = "hippocampus", hemi = hemi_avg$hemisphere, z_value = hemi_avg$stat_avg)
    p_ho <- ggplot2::ggplot(ho_df, ggplot2::aes(fill = z_value)) +
      ggseg::geom_brain(atlas = ho, position = ggseg::position_brain(hemi ~ side)) +
      ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, name = "Effect Size") +
      ggplot2::labs(title = "Hippocampus (Harvard-Oxford)", subtitle = "Hemisphere-averaged across subfields") +
      theme_brain_pub()
    save_pub(file.path(output_dir, "Hippocampus_HO_Aggregate_Map"), p_ho, width = 10, height = 6, dpi = 600)
  }

  utils::write.csv(hippo_summary, file.path(output_dir, "Hippocampal_Subfields_Summary.csv"), row.names = FALSE)
  cat("✓ 海马亚区映射完成\n")
  return(list(plot_bar = p_hippo, plot_ho = p_ho, data = hippo_summary))
}


# =================== 5.5 血管分布区域分析 ===================
create_vascular_territory_analysis <- function(brain_mapping_data, output_dir = file.path("Brain_MRI_Maps", "Vascular")) {
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  assign_vascular_territory <- function(region_name) {
    v <- tolower(region_name)
    dplyr::case_when(
      # ACA 领地
      stringr::str_detect(v, "(anterior\\.?cingulate|rostral\\.?anterior\\.?cingulate|caudal\\.?anterior\\.?cingulate|paracentral|medial\\.?frontal|superior\\.?frontal|frontal\\.?medial|rostral\\.?middle\\.?frontal)") ~ "ACA",
      # MCA 领地
      stringr::str_detect(v, "(insula|precentral|postcentral|middle\\.?temporal|superior\\.?temporal|inferior\\.?temporal|supramarginal|angular|inferior\\.?parietal|parsopercularis|parstriangularis|pars\\.?orbitalis|middle\\.?frontal|bank\\.?sts)") ~ "MCA",
      # PCA 领地
      stringr::str_detect(v, "(occipital|calcarine|pericalcarine|cuneus|lingual|lateral\\.?occipital|fusiform|inferior\\.?temporal|hippocamp|parahip|entorhinal|precuneus|posterior\\.?cingulate)") ~ "PCA",
      # 椎基底系统
      stringr::str_detect(v, "(cerebell|brainstem|pons|medulla|vermis)") ~ "Vertebrobasilar",
      TRUE ~ "Uncertain"
    )
  }

  vas_df <- brain_mapping_data %>%
    dplyr::filter(.data[["selected_for_mapping"]]) %>%
    dplyr::mutate(region_src = dplyr::coalesce(ggseg_region, brain_region, anatomical_structure)) %>%
    dplyr::filter(!is.na(region_src) & region_src != "") %>%
    dplyr::mutate(territory = assign_vascular_territory(region_src)) %>%
    dplyr::group_by(territory) %>%
    dplyr::summarise(
      max_z = max(abs(z_change), na.rm = TRUE),
      mean_z = mean(abs(z_change), na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(max_z))

  if(nrow(vas_df) == 0) return(NULL)

  p_vas <- ggplot2::ggplot(vas_df, ggplot2::aes(x = reorder(territory, max_z), y = max_z, fill = mean_z)) +
    ggplot2::geom_col() + ggplot2::coord_flip() +
    viridis::scale_fill_viridis(option = "C", direction = 1, name = "Mean |Z|") +
    ggplot2::labs(title = "Vascular Territories Analysis", subtitle = "Max and mean |Z| per territory") +
    theme_brain_pub()
  save_pub(file.path(output_dir, "Vascular_Territories_Analysis"), p_vas, width = 9, height = 6, dpi = 600)

  utils::write.csv(vas_df, file.path(output_dir, "Vascular_Territories_Summary.csv"), row.names = FALSE)
  invisible(vas_df)
}

# =================== 5.5.1 脑血管高风险区映射（ACA/MCA/PCA） ===================

# 脑血管领地的脑图可视化（DK atlas）：将 ACA/MCA/PCA 对应皮层区域在脑图上高亮
visualize_vascular_territories_brain_map <- function(brain_mapping_data, output_dir = file.path("Brain_MRI_Maps", "Vascular")) {
  load_brain_mapping_packages_once()
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (!requireNamespace("ggseg", quietly = TRUE)) return(NULL)
  if(!exists("dk")) {
    suppressWarnings(data(dk, package = "ggseg"))
  }
  dk_obj <- if (exists("dk", inherits = TRUE)) get("dk", inherits = TRUE) else NULL
  dk_df <- if (!is.null(dk_obj) && is.list(dk_obj) && "data" %in% names(dk_obj)) dk_obj$data else dk_obj
  valid_regions <- if (is.data.frame(dk_df) && "region" %in% names(dk_df)) unique(as.character(dk_df$region)) else character(0)
  valid_regions <- valid_regions[!is.na(valid_regions) & nzchar(valid_regions)]
  norm_key <- function(x) {
    x <- tolower(trimws(as.character(x)))
    x <- gsub("[^a-z0-9]+", "", x)
    x
  }
  dk_key <- norm_key(valid_regions)
  dk_map <- stats::setNames(valid_regions, dk_key)
  canon_region <- function(x) {
    if (is.null(x) || is.na(x) || !nzchar(x)) return(NA_character_)
    if (x %in% valid_regions) return(x)
    k <- norm_key(x)
    if (k %in% names(dk_map)) return(unname(dk_map[[k]]))
    NA_character_
  }

  # 只映射皮层区域，使用统一的规范化函数归类到 ACA/MCA/PCA

  vascular_territories <- data.frame(
    region = c(
      # MCA
      "superiortemporal","middletemporal","inferiortemporal","supramarginal","parsopercularis","parstriangularis","postcentral","precentral",
      # ACA
      "superiorfrontal","rostralmiddlefrontal","medialorbitofrontal","rostralanteriorcingulate","caudalanteriorcingulate",
      # PCA
      "lateraloccipital","lingual","cuneus","pericalcarine","fusiform","parahippocampal","entorhinal"
    ),
    artery = c(rep("MCA", 8), rep("ACA", 5), rep("PCA", 7)),
    stringsAsFactors = FALSE
  )
  vascular_territories$region <- vapply(as.character(vascular_territories$region), canon_region, character(1))
  vascular_territories <- vascular_territories[!is.na(vascular_territories$region) & vascular_territories$region != "", , drop = FALSE]

  # 规范化函数就绪检查：优先全局 normalize_dk_region，缺失时使用局部备用实现
  local_normalizer <- if (exists("normalize_dk_region", inherits = TRUE)) {
    normalize_dk_region
  } else {
    function(s) {
      v <- tolower(s)
      v <- gsub("[_.]+", " ", v)
      v <- gsub("\\n", " ", v)
      v <- gsub("[[:space:]]+", " ", v)
      dplyr::case_when(
        stringr::str_detect(v, "superior[_\\s-]?temporal") ~ "superiortemporal",
        stringr::str_detect(v, "middle[_\\s-]?temporal") ~ "middletemporal",
        stringr::str_detect(v, "inferior[_\\s-]?temporal") ~ "inferiortemporal",
        stringr::str_detect(v, "supramarginal") ~ "supramarginal",
        stringr::str_detect(v, "pars[_\\s-]?opercularis|opercularis") ~ "parsopercularis",
        stringr::str_detect(v, "pars[_\\s-]?triangularis|triangularis") ~ "parstriangularis",
        stringr::str_detect(v, "postcentral") ~ "postcentral",
        stringr::str_detect(v, "precentral") ~ "precentral",
        stringr::str_detect(v, "superior[_\\s-]?frontal") ~ "superiorfrontal",
        stringr::str_detect(v, "rostral[_\\s-]?middle[_\\s-]?frontal") ~ "rostralmiddlefrontal",
        stringr::str_detect(v, "medial[_\\s-]?orbit[_\\s-]?ofrontal|medial[_\\s-]?orbito[_\\s-]?frontal") ~ "medialorbitofrontal",
        stringr::str_detect(v, "rostral[_\\s-]?anterior[_\\s-]?cingulate") ~ "rostralanteriorcingulate",
        stringr::str_detect(v, "caudal[_\\s-]?anterior[_\\s-]?cingulate") ~ "caudalanteriorcingulate",
        stringr::str_detect(v, "lateral[_\\s-]?occipital") ~ "lateraloccipital",
        stringr::str_detect(v, "\\bcuneus\\b") ~ "cuneus",
        stringr::str_detect(v, "lingual") ~ "lingual",
        stringr::str_detect(v, "peri[_\\s-]?calcarine|calcarine") ~ "pericalcarine",
        stringr::str_detect(v, "fusiform") ~ "fusiform",
        stringr::str_detect(v, "parahipp") ~ "parahippocampal",
        stringr::str_detect(v, "entorhinal") ~ "entorhinal",
        stringr::str_detect(v, "hippocamp") ~ "hippocampus",
        TRUE ~ NA_character_
      )
    }
  }

  df <- brain_mapping_data %>%
    dplyr::filter(structure_type == "Cortical") %>%
    dplyr::mutate(src = dplyr::coalesce(ggseg_region, brain_region, anatomical_structure)) %>%
    dplyr::mutate(region = local_normalizer(src)) %>%
    dplyr::mutate(region = vapply(as.character(region), canon_region, character(1))) %>%
    dplyr::filter(!is.na(region)) %>%
    dplyr::group_by(region, hemi) %>%
    dplyr::summarise(z_value = max(abs(z_change), na.rm = TRUE), .groups = "drop") %>%
    dplyr::left_join(vascular_territories, by = "region") %>%
    dplyr::filter(!is.na(artery))

  if(nrow(df) == 0) return(NULL)

  # 将不同动脉领地着色（红=高风险MCA，橙=ACA，蓝=低风险PCA）
  artery_colors <- c("MCA" = "#d62728", "ACA" = "#ff7f0e", "PCA" = "#1f77b4")
  p_brain <- ggplot2::ggplot(df, ggplot2::aes(fill = artery)) +
    ggseg::geom_brain(atlas = dk) +
    ggplot2::scale_fill_manual(values = artery_colors, name = "Territory") +
    ggplot2::labs(title = "Vascular Territories (DK atlas)", subtitle = paste("n:", nrow(df))) +
    theme_brain_pub()
  save_pub(file.path(output_dir, "Vascular_Territories_Brain_Map"), p_brain, width = 12, height = 8, dpi = 600)

  invisible(p_brain)
}
# Initialize visualization-related variables to avoid missing object errors during wrapup
vis_idps <- NULL
has_vis <- FALSE

has_idps <- exists("independent_effect_idps") && !is.null(independent_effect_idps) && nrow(independent_effect_idps) > 0
# 通过回退机制获取可视化IDP
vis_idps <- get_visualization_idps(
  if (has_idps) independent_effect_idps else NULL,
  if (exists("combined_four_model_results")) combined_four_model_results else NULL,
  fallback_n = NULL
)
has_vis <- !is.null(vis_idps) && nrow(vis_idps) > 0

if ((RUN_VISUALIZATION || RUN_BRAIN_MAPPING) && has_vis) {
  cat("开始创建完善的可视化图表...\n")
  
  # 统一校验/创建输出目录
  verify_output_dirs()
  
  idps_for_plots <- vis_idps
  csv_indep <- table_out_path("Visualization", "Independent_Brain_IDPs_for_Visualization.csv")
  if (!is.null(idps_for_plots) && nrow(idps_for_plots) > 0) {
    utils::write.csv(as.data.frame(idps_for_plots), csv_indep, row.names = FALSE)
    cat(sprintf("已基于 Combined_Four_Model_Z_Analysis_Results.csv 生成可视化IDP列表：%s\n", csv_indep))
  }
  # 修复重复列名以避免dplyr在filter/mutate中报错
  if (!is.null(idps_for_plots) && nrow(idps_for_plots) > 0) {
    old_names <- names(idps_for_plots)
    if (anyDuplicated(old_names) > 0) {
      dup_count <- sum(duplicated(old_names))
      names(idps_for_plots) <- make.unique(old_names)
      cat(sprintf("检测到可视化IDP数据中有 %d 个重复列名，已自动修复。\n", dup_count))
    }
    idps_for_plots <- tibble::as_tibble(idps_for_plots)
  }
  elig <- summarize_plot_eligibility(idps_for_plots, analysis_data_list)
  cat(sprintf("候选IDP总数: %d; 具备完整变量的可绘制IDP: %d\n", elig$total, elig$eligible))

  completeness_trends <- NULL
  
  if (RUN_VISUALIZATION) {
    # 1. 创建每个独立效应IDP的单独箱线图
    cat("\n=== 创建单独箱线图 ===\n")
    individual_plots <- create_idp_boxplots_unified(
      independent_idps = idps_for_plots,
      analysis_data_list = analysis_data_list,
      output_dir = "Individual_IDP_Plots",
      mode = "individual",
      top_n = min(30, nrow(idps_for_plots)),
      delta_normalize_by = "baseline_sd",
      percent_baseline = "baseline_mean"
    )
    
    # 2. 统一年龄趋势图（默认：percentage单面版带散点）
    cat("\n=== 统一年龄趋势图（默认：percentage单面版带散点） ===\n")
    completeness_trends <- create_age_trend_plots_fast(
      independent_idps = idps_for_plots,
      analysis_data_list = analysis_data_list,
      top_n = NULL,
      min_group_n = 5,
      panel = "percentage",
      include_points = TRUE,
      delta_normalize_by = "baseline_sd",
      percent_baseline = "baseline_mean",
      parallel = TRUE,
      point_downsample_frac = 0.5,
      smooth_sample_frac = 0.5,
      smooth_sample_max = 5000,
      se_confint = TRUE,
      dpi = 600,
      export_raw_data = TRUE
    )
    # 年龄×组别斜率差异检验（与上面的归一化选择保持一致）
    slope_df <- compute_age_slope_difference_for_idps(
      independent_idps = idps_for_plots,
      analysis_data_list = analysis_data_list,
      used_metric = "baseline_mean_percent",
      min_group_n = 5,
      output_path = file.path("Age_Trend_Plots", "Age_Trend_Slope_Difference_AllModels.csv")
    )
  }

  brain_independent <- data.frame()
  
  if (RUN_BRAIN_MAPPING) {
    # 3. 创建MRI脑图映射（仅脑影像IDP）
    cat("\n=== 创建MRI脑图映射 ===\n")
    combined_src_map <- NULL
    if (exists("combined_four_model_results") && !is.null(combined_four_model_results) && nrow(combined_four_model_results) > 0) {
      combined_src_map <- tibble::as_tibble(combined_four_model_results)
    } else {
      csv_path_z <- combined_four_model_results_csv()
      if (file.exists(csv_path_z)) {
        combined_src_map <- tryCatch({
          tmp <- utils::read.csv(csv_path_z, stringsAsFactors = FALSE)
          tmp <- clean_combined_results_df(tmp)
          tibble::as_tibble(tmp)
        }, error = function(e) NULL)
      }
    }
    brain_independent <- if (!is.null(combined_src_map) && nrow(combined_src_map) > 0) {
      df <- combined_src_map
      if ("variable_type" %in% names(df)) df <- df[df$variable_type == "Brain_Imaging", , drop = FALSE]
      if ("independent_effect" %in% names(df)) df <- df[df$independent_effect == TRUE, , drop = FALSE]
      df
    } else {
      dplyr::filter(idps_for_plots, variable_type == "Brain_Imaging")
    }
    
    if (nrow(brain_independent) > 0) {
      needs_z <- (!("z_statistic_without_non_imaging" %in% names(brain_independent)) ||
                    all(is.na(brain_independent$z_statistic_without_non_imaging)))
      if (needs_z) {
        combined_src <- NULL
        if (exists("combined_four_model_results") && !is.null(combined_four_model_results) && nrow(combined_four_model_results) > 0) {
          combined_src <- tibble::as_tibble(combined_four_model_results)
        } else {
          csv_path_z <- combined_four_model_results_csv()
          if (file.exists(csv_path_z)) {
            combined_src <- tryCatch({
              tmp <- utils::read.csv(csv_path_z, stringsAsFactors = FALSE)
              tmp <- clean_combined_results_df(tmp)
              tibble::as_tibble(tmp)
            }, error = function(e) NULL)
          }
        }
        if (!is.null(combined_src) && nrow(combined_src) > 0) {
          z_map <- combined_src[, intersect(c("model_comparison","variable_name","variable_type","z_statistic_without_non_imaging","z_statistic_with_non_imaging","independent_effect"), names(combined_src)), drop = FALSE]
          z_map <- dplyr::distinct(z_map, model_comparison, variable_name, variable_type, .keep_all = TRUE)
          brain_independent <- dplyr::left_join(brain_independent, z_map, by = c("model_comparison","variable_name","variable_type"))
        }
      }
  
      mapping_all <- create_mri_brain_mapping(brain_independent, output_dir = "Brain_MRI_Maps")
      if (!is.null(mapping_all) && !is.null(mapping_all$brain_mapping_data) && nrow(mapping_all$brain_mapping_data) > 0) {
        all_indep <- NULL
        if (!is.null(combined_src_map) && nrow(combined_src_map) > 0) {
          df <- combined_src_map
          if ("independent_effect" %in% names(df)) df <- df[df$independent_effect == TRUE, , drop = FALSE]
          all_indep <- df
        } else if (exists("independent_effect_idps") && !is.null(independent_effect_idps) && nrow(independent_effect_idps) > 0) {
          all_indep <- tibble::as_tibble(independent_effect_idps)
        } else {
          all_indep <- tibble::as_tibble(brain_independent)
        }
        map_df <- tibble::as_tibble(mapping_all$brain_mapping_data)
        key_cols <- intersect(
          c("model_comparison", "variable_name", "variable_type", "idp1_variable", "idp2_variable"),
          intersect(names(all_indep), names(map_df))
        )
        if (length(key_cols) == 0) {
          key_cols <- intersect(c("model_comparison", "variable_name"), intersect(names(all_indep), names(map_df)))
        }
        add_cols <- setdiff(names(map_df), names(all_indep))
        map_join <- if (length(add_cols) > 0) {
          map_df %>%
            dplyr::select(dplyr::all_of(c(key_cols, add_cols))) %>%
            dplyr::distinct(dplyr::across(dplyr::all_of(key_cols)), .keep_all = TRUE)
        } else {
          map_df %>%
            dplyr::select(dplyr::all_of(key_cols)) %>%
            dplyr::distinct(dplyr::across(dplyr::all_of(key_cols)), .keep_all = TRUE)
        }
        unified_ie <- tibble::as_tibble(all_indep) %>% dplyr::left_join(map_join, by = key_cols)
        if (!("selected_for_mapping" %in% names(unified_ie))) unified_ie$selected_for_mapping <- FALSE
        unified_ie$selected_for_mapping <- suppressWarnings(as.logical(unified_ie$selected_for_mapping))
        unified_ie$selected_for_mapping[is.na(unified_ie$selected_for_mapping)] <- FALSE
        if ("variable_type" %in% names(unified_ie)) {
          unified_ie$selected_for_mapping[unified_ie$variable_type != "Brain_Imaging"] <- FALSE
        }
        utils::write.csv(as.data.frame(unified_ie), all_independent_effect_idps_csv(), row.names = FALSE)
      }
  
      # 按比较组分别映射
      if ("model_comparison" %in% names(brain_independent)) {
        comparisons <- unique(brain_independent$model_comparison)
        cat(sprintf("检测到 %d 个比较组，将分别生成脑图...\n", length(comparisons)))
        
        for (comp in comparisons) {
          cat(sprintf("\n--- 处理比较组: %s ---\n", comp))
          sub_data <- brain_independent %>% dplyr::filter(model_comparison == comp)
          comp_dir <- file.path("Brain_MRI_Maps", comp)
          
          mapping_result <- create_mri_brain_mapping(sub_data, output_dir = comp_dir)
          
          if (!is.null(mapping_result) && !is.null(mapping_result$brain_mapping_data)) {
            suppressWarnings(create_cortical_mapping(mapping_result$brain_mapping_data, output_dir = file.path(comp_dir, "Cortical")))
            # create_subcortical_mapping removed per user request
            suppressWarnings(create_hippocampal_mapping(mapping_result$brain_mapping_data, output_dir = file.path(comp_dir, "Hippocampus")))
            suppressWarnings(create_network_mapping(mapping_result$brain_mapping_data, output_dir = file.path(comp_dir, "Network")))
            suppressWarnings(create_vascular_territory_analysis(mapping_result$brain_mapping_data, output_dir = file.path(comp_dir, "Vascular")))
            suppressWarnings(visualize_vascular_territories_brain_map(mapping_result$brain_mapping_data, output_dir = file.path(comp_dir, "Vascular")))
          }
        }
      } else {
        # 兼容旧版无 model_comparison 列的情况
        cat("未检测到 model_comparison 列，进行整体映射...\n")
        mapping_result <- create_mri_brain_mapping(brain_independent)
        if (!is.null(mapping_result) && !is.null(mapping_result$brain_mapping_data)) {
          suppressWarnings(create_cortical_mapping(mapping_result$brain_mapping_data))
          # create_subcortical_mapping removed per user request
          suppressWarnings(create_hippocampal_mapping(mapping_result$brain_mapping_data))
          suppressWarnings(create_network_mapping(mapping_result$brain_mapping_data))
          suppressWarnings(create_vascular_territory_analysis(mapping_result$brain_mapping_data))
          suppressWarnings(visualize_vascular_territories_brain_map(mapping_result$brain_mapping_data))
        }
      }
      
      cat("脑影像IDP映射完成:\n")
      cat("- 解析脑影像IDP:", nrow(brain_independent), "个\n")
      cat("- 映射数据已保存至 Brain_MRI_Maps/<Comparison>/ 目录\n")
      cat("- 已生成皮层、海马亚区、功能网络(Yeo7)、血管分布区域分析\n")
    } else {
      cat("⚠️  无脑影像独立效应IDP可供映射\n")
    }
    
    # 3.1 综合脑变化总览气泡图（使用解析后的映射CSV）
    cat("\n=== 创建脑变化总览气泡图 ===\n")
    try({
      render_brain_changes_overview(output_dir = getwd())
      cat("已生成: Brain_Changes_Overview/Comprehensive_Brain_Changes_Overview.{pdf,png}\n")
    }, silent = TRUE)
  }
  
  # 按用户要求：移除严格集合相关脑图映射与气泡总览
  
  # 4. 创建综合摘要报告
  cat("\n=== 创建综合摘要报告 ===\n")
  # 兜底：若 independent_effect_idps 不在内存，尝试从CSV读取；若仍不可用则置为空数据框
  if (!exists("independent_effect_idps") || is.null(independent_effect_idps)) {
    csv_path_ie <- all_independent_effect_idps_csv()
    if (file.exists(csv_path_ie)) {
      independent_effect_idps <- tryCatch(utils::read.csv(csv_path_ie, stringsAsFactors = FALSE), error = function(e) NULL)
    }
    if (is.null(independent_effect_idps)) independent_effect_idps <- data.frame()
  }
  safe_nrow <- function(x) { if (is.null(x)) 0L else nrow(x) }
  # 聚合完备性计数
  agg_abs_two <- 0
  agg_abs_single <- 0
  agg_abs_skipped <- 0
  agg_pct_two <- if (!is.null(completeness_trends)) sum(completeness_trends$status == "two_group", na.rm = TRUE) else 0
  agg_pct_single <- if (!is.null(completeness_trends)) sum(completeness_trends$status == "single_group", na.rm = TRUE) else 0
  agg_pct_skipped <- if (!is.null(completeness_trends)) sum(completeness_trends$status == "skipped", na.rm = TRUE) else 0
  
  visualization_summary <- data.frame(
    Visualization_Type = c(
      "Individual IDP Boxplots", 
      "Age-Related Trends", 
      "MRI Brain Mapping", 
      "Brain Changes Overview"
    ),
    Files_Generated = c(
      sprintf("%d PNG + PDF files", safe_nrow(independent_effect_idps)),
      sprintf("%d age trend files", min(30, safe_nrow(independent_effect_idps))),
      sprintf("%d brain mapping files", ifelse(nrow(brain_independent) > 0, 6, 0)),
      "Overview bubble chart"
    ),
    Output_Directory = c(
      "Individual_IDP_Plots/",
      "Age_Trend_Plots/", 
      "Brain_MRI_Maps/",
      "Brain_Changes_Overview/"
    ),
    Description = c(
      "Box plots with jittered points for each significant IDP",
      "年龄相关百分比变化趋势（单面板，含散点）", 
      "MRI brain mapping with region analysis",
      "Comprehensive bubble chart summarizing regional brain changes"
    ),
    Min_Group_N_Threshold = c(NA, 5, NA, NA),
    Age_Abs_TwoGroup = c(NA, agg_abs_two, NA, NA),
    Age_Abs_SingleGroup = c(NA, agg_abs_single, NA, NA),
    Age_Abs_Skipped = c(NA, agg_abs_skipped, NA, NA),
    Age_Pct_TwoGroup = c(NA, agg_pct_two, NA, NA),
    Age_Pct_SingleGroup = c(NA, agg_pct_single, NA, NA),
    Age_Pct_Skipped = c(NA, agg_pct_skipped, NA, NA),
    stringsAsFactors = FALSE
  )
  # 详细的完备性报表（按IDP逐项记录）
  try({
    comp_unified_path <- table_out_path(file.path("Visualization", "Age_Trend"), "Age_Trend_Completeness.csv")
    if (!is.null(completeness_trends)) utils::write.csv(completeness_trends, comp_unified_path, row.names = FALSE)
  }, silent = TRUE)
  write.csv(visualization_summary, table_out_path("Visualization", "Visualization_Summary_Report.csv"), row.names = FALSE)
  print(visualization_summary)
}

if (!has_vis) {
  cat("⚠️  无可用于回退的IDP，未生成新的图像\n")
  cat("请检查四模型综合结果或重新运行分析以产生候选IDP\n")
}

 

# =================== 绘图总结报告 ===================
cat("\n", rep("=", 60), "\n")
cat("🎉 完善的缺血性心肌病IDP可视化分析完成！\n")
cat(rep("=", 60), "\n")

cat("生成的可视化内容:\n")
cat("📊 Individual_IDP_Plots/ - 每个独立效应IDP的单独箱线图\n")
cat("📈 Age_Trend_Plots/ - IDP变化随年龄变化的趋势图\n") 
cat("🧠 Brain_MRI_Maps/ - MRI脑图映射和区域分析\n")

cat("📈 Forest_Plots/ - 森林图（四模型与敏感性对比）\n")
cat("📋 Visualization_Summary_Report.csv - 可视化摘要报告\n")

cat("\n可视化特点:\n")
cat("✓ 参考您提供的样式设计箱线图和年龄趋势图\n")
cat("✓ 每个独立效应IDP都有单独的详细图表\n")
cat("✓ 包含绝对变化和百分比变化两种视角\n")
cat("✓ MRI脑图映射基于真实脑区解剖结构\n")
cat("✓ 使用ggseg包创建标准脑图（如可用）\n")
cat("✓ Python脚本提供高级3D可视化选项\n")

cat("\n后续使用建议:\n")
cat("1. 查看 Individual_IDP_Plots/ 中每个IDP的详细图表\n")
cat("2. 分析 Age_Trend_Plots/ 中的年龄相关变化模式\n") 
cat("3. 运行 Python 脚本获取高级脑图映射\n")
cat("4. 基于可视化结果进一步解释生物学意义\n")

cat(rep("=", 60), "\n")
# =================== 总体总结报告 ===================
cat("\n", rep("=", 60), "\n")
cat("🎉 四模型缺血性心肌病IDP纵向变化Z值分析完成！\n")
cat(rep("=", 60), "\n")

if(exists("all_four_model_results") && length(all_four_model_results) > 0) {
  
  cat("分析总览:\n")
  for(model_name in names(all_four_model_results)) {
    results <- all_four_model_results[[model_name]]
    independent_count <- sum(results$independent_effect)
    cat(sprintf("%-25s: %4d变量 (%3d认知, %3d脑影像), %3d独立效应IDP\n",
                model_name,
                nrow(results),
                sum(results$variable_type == "Cognitive"),
                sum(results$variable_type == "Brain_Imaging"),
                independent_count))
  }
  
  if(exists("independent_effect_idps")) {
    cat(sprintf("\n总独立效应IDP: %d 个\n", nrow(independent_effect_idps)))
    cat("- 这些IDP不受混淆变量影响，表明缺血性心肌病在年龄之外对认知和脑影像的独立影响\n")
  }
}

cat("\n四个分析模型:\n")
cat("✓ 缺血性心肌病 vs 对照组\n")
cat("✓ 心梗 vs 对照组\n")
cat("✓ 慢性缺血性 vs 对照组\n")
cat("✓ 心梗 vs 慢性缺血性\n")

cat("\nZ值分析特点:\n")
cat("✓ Case-versus-Control = demeaned(binary) × 10^(Age2×0.0524−3.27)\n")
cat("✓ 混淆变量: 年龄差异、年龄平方差、种族、性别\n")
cat("✓ 关注独立效应IDP（不受混淆变量影响）\n")
cat("✓ 为后续非成像变量分析提供基础\n")

cat("\n生成文件:\n")
cat("📊 Four_Model_Z_Analysis_[model_name].csv (各模型Z统计量结果)\n")
cat("📊 Combined_Four_Model_Z_Analysis_Results.csv (综合Z统计量结果)\n")
cat("📊 Four_Model_Independent_Effects_Summary.csv (独立效应摘要)\n")
cat("📊 All_Independent_Effect_IDPs_for_Visualization.csv (所有独立效应IDP)\n")
cat("📊 Strict_FDR_Independent_IDPs.csv (严格集合，仅FDR<0.05)\n")

cat("📊 MI_vs_Chronic_Z_Correlation_Analysis.csv (心梗vs慢性缺血性相关性)\n")

cat("\n可视化文件:\n")
cat("📈 Independent_IDPs_Boxplots.png/pdf (分组箱线图)\n")
cat("📈 Brain_MRI_Maps/Brain_IDPs_Parsed_for_Mapping.csv (脑影像映射数据)\n")

cat("📈 Forest_Plots/Forest_Full_Top20.png/pdf (四模型Full森林图)\n")
cat("📈 Forest_Plots/Forest_Sensitivity_Top10.png/pdf (敏感性对比森林图)\n")

# 结果文件简洁展示与汇总
report_results_overview <- function() {
  files <- c(
    list.files(get_table_dir("Cohort", "Matched"), pattern = "Final_Matched_Cohort_.*\\.(csv|rds)$", full.names = FALSE),
    list.files(getwd(), pattern = "Final_Matched_Cohort_.*\\.(csv|rds)$", full.names = FALSE),
    list.files(get_table_dir("Four_Model"), pattern = "\\.csv$", full.names = TRUE, recursive = TRUE),
    list.files(file.path(getwd(), "Brain_MRI_Maps"), pattern = "\\.(pdf|png|tiff)$", full.names = FALSE),
    list.files(file.path(getwd(), "Forest_Plots"), pattern = "\\.(pdf|png)$", full.names = FALSE),
    list.files(getwd(), pattern = "Figure.*\\.(pdf|png|tiff)$", full.names = FALSE)
  )
  overview_df <- tibble::tibble(File = files)
  readr::write_csv(overview_df, table_out_path("Overview", "Results_Overview.csv"))
  cat("\n=== 结果文件清单（简洁展示） ===\n")
  cat(paste0("- ", files), sep = "\n")
}

# 执行结果文件汇总展示
report_results_overview()

cat(rep("=", 60), "\n")
cat("分析完成！独立效应IDP已识别，可进行后续研究。\n")

# ============================================================================
# 脑龄预测与脑龄差（BAG）统计与可视化
# ============================================================================

# 轻量日志输出（带时间戳）
log_message <- function(...) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s\n", ts, paste(..., collapse = " ")))
}

# 一次性尝试加载脑龄建模所需包（caret / randomForest / data.table / ggplot2 / viridis）
load_brain_age_packages_once <- local({
  loaded <- FALSE
  function(verbose = TRUE) {
    if (loaded) return(invisible(TRUE))
    pkgs <- c("caret", "randomForest", "data.table", "ggplot2", "viridis")
    ok <- logical(length(pkgs))
    for (i in seq_along(pkgs)) {
      p <- pkgs[i]
      ok[i] <- requireNamespace(p, quietly = TRUE)
      if (verbose) cat(sprintf("- %s: %s\n", p, if (ok[i]) "OK" else "Missing"))
    }
    if (!all(ok)) {
      miss <- paste(pkgs[!ok], collapse = ", ")
      warning(sprintf("脑龄分析所需包缺失: %s。请先安装，例如: install.packages(c(%s))",
                      miss, paste(sprintf("\"%s\"", pkgs[!ok]), collapse = ", ")))
    }
    loaded <- TRUE
    invisible(all(ok))
  }
})

# 自动选择可用于脑龄预测的IDP特征列
# 规则：
# - 优先使用 Independent_Brain_IDPs_for_Visualization.csv 中的 Brain_Imaging 变量
# - 若不存在，则回退到 All_Independent_Effect_IDPs_for_Visualization.csv
# - 再回退到 Brain_IDPs_Parsed_for_Mapping.csv（从解析结果中抽取 topN 的度量列）
get_brain_age_idps <- function(top_n = 64, prefer_csv = TRUE) {
  # 候选CSV路径（相对于工作目录）
  csv_indep <- table_out_path("Visualization", "Independent_Brain_IDPs_for_Visualization.csv")
  csv_all   <- all_independent_effect_idps_csv()
  csv_map   <- file.path(getwd(), "Brain_MRI_Maps", "Brain_IDPs_Parsed_for_Mapping.csv")

  idp_vars <- character(0)
  src <- NULL

  if (prefer_csv && file.exists(csv_indep)) {
    dt <- tryCatch(utils::read.csv(csv_indep, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(dt)) {
      src <- "Independent_Brain_IDPs_for_Visualization.csv"
      cols <- intersect(c("idp1_variable","idp2_variable","variable_name","IDP"), names(dt))
      if (length(cols) > 0) {
        raw <- unique(unlist(dt[cols]))
        raw <- raw[!is.na(raw) & nzchar(raw)]
        idp_vars <- raw
      }
    }
  }

  if (length(idp_vars) == 0 && file.exists(csv_all)) {
    dt <- tryCatch(utils::read.csv(csv_all, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(dt)) {
      src <- "All_Independent_Effect_IDPs_for_Visualization.csv"
      cols <- intersect(c("idp1_variable","idp2_variable","variable_name","IDP"), names(dt))
      if (length(cols) > 0) {
        raw <- unique(unlist(dt[cols]))
        raw <- raw[!is.na(raw) & nzchar(raw)]
        idp_vars <- raw
      }
    }
  }

  if (length(idp_vars) == 0 && file.exists(csv_map)) {
    dt <- tryCatch(utils::read.csv(csv_map, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(dt)) {
      src <- "Brain_IDPs_Parsed_for_Mapping.csv"
      candidate_cols <- c("variable_name","IDP","anatomical_structure","brain_region_final")
      cols <- intersect(candidate_cols, names(dt))
      if (length(cols) > 0) {
        raw <- unique(unlist(dt[cols]))
        raw <- raw[!is.na(raw) & nzchar(raw)]
        idp_vars <- raw
      }
    }
  }

  # 清洗变量名
  idp_vars <- unique(idp_vars)
  idp_vars <- idp_vars[grepl("[A-Za-z]", idp_vars)]
  if (length(idp_vars) > top_n) idp_vars <- idp_vars[1:top_n]
  log_message(sprintf("脑龄预测将使用特征列（源: %s）: %d 个", ifelse(is.null(src), "unknown", src), length(idp_vars)))
  idp_vars
}

# 训练脑龄预测模型并生成脑龄差（BAG）统计与可视化
# @param data       健康对照组数据（用于训练）
# @param test_data  病例/其他组数据（用于评估与可视化）
# @param idp_list   要用于建模的IDP列名向量
# @param output_dir 输出目录
create_brain_age_gap_analysis <- function(data, test_data, idp_list, output_dir = "Brain_Age_Analysis") {
  load_brain_age_packages_once()
  log_message("开始脑龄预测分析...")

  # 转换为 data.table（若可用）
  dt_train <- if (requireNamespace("data.table", quietly = TRUE)) data.table::as.data.table(data) else data
  dt_test  <- if (requireNamespace("data.table", quietly = TRUE)) data.table::as.data.table(test_data) else test_data

  # 输出目录
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # 检查关键列
  if (!("baseline_age" %in% names(dt_train))) stop("训练数据缺少 baseline_age 列")
  if (!("baseline_age" %in% names(dt_test))) stop("测试数据缺少 baseline_age 列")
  
  idp_list <- unique(as.character(idp_list))
  idp_list <- idp_list[!is.na(idp_list) & nzchar(idp_list)]
  if (length(idp_list) == 0) stop("未提供可用于脑龄预测的IDP列")
  
  idp_in_train <- intersect(idp_list, names(dt_train))
  miss_train <- setdiff(idp_list, idp_in_train)
  if (length(miss_train) > 0) {
    log_message(sprintf("训练数据缺少 %d 个IDP列，将跳过: %s", length(miss_train), paste(miss_train, collapse = ", ")))
  }
  idp_in_both <- intersect(idp_in_train, names(dt_test))
  miss_test <- setdiff(idp_in_train, idp_in_both)
  if (length(miss_test) > 0) {
    log_message(sprintf("测试数据缺少 %d 个IDP列，将跳过: %s", length(miss_test), paste(miss_test, collapse = ", ")))
  }
  idp_list <- idp_in_both
  if (length(idp_list) < 5) {
    stop(sprintf("可用于建模的IDP列太少（%d个），无法训练脑龄模型", length(idp_list)))
  }

  # 准备训练数据（仅健康对照）
  train_vars <- c("baseline_age", idp_list)
  # 统一使用 .SDcols 选择列，兼容不同版本 data.table
  dt_train_clean <- if (inherits(dt_train, "data.table")) {
    stats::na.omit(dt_train[, .SD, .SDcols = train_vars])
  } else {
    stats::na.omit(dt_train[, train_vars, drop = FALSE])
  }

  # 训练随机森林模型
  log_message("训练脑龄预测模型（Random Forest）...")
  mtry_val <- max(1L, floor(sqrt(length(idp_list))))
  rf_model <- randomForest::randomForest(
    baseline_age ~ .,
    data = dt_train_clean,
    ntree = 500,
    mtry = mtry_val,
    importance = TRUE
  )

  # 近后期树的 R² 估计（稳健切片）
  if (!is.null(rf_model$rsq) && length(rf_model$rsq) >= 50) {
    tail_idx <- seq(max(1, length(rf_model$rsq) - 99), length(rf_model$rsq))
    r2_mean <- mean(rf_model$rsq[tail_idx], na.rm = TRUE)
    log_message(sprintf("模型 R²（后期树平均）= %.3f", r2_mean))
  }

  # 预测病例组的脑龄
  # 若存在分组列，将其与训练列一起保留
  group_col <- intersect(c("diagnosis","cohort","group","case_control"), names(dt_test))
  # 若不存在分组列，避免将 NA 传入列选择
  extras_keep <- intersect(c("eid"), names(dt_test))
  test_vars <- if (length(group_col) > 0) unique(c(train_vars, group_col[1], extras_keep)) else unique(c(train_vars, extras_keep))
  dt_test_clean <- if (inherits(dt_test, "data.table")) {
    stats::na.omit(dt_test[, .SD, .SDcols = test_vars])
  } else {
    stats::na.omit(dt_test[, test_vars, drop = FALSE])
  }
  dt_test_clean$predicted_brain_age <- stats::predict(rf_model, newdata = dt_test_clean)

  # 计算脑龄差（Brain Age Gap, BAG）
  dt_test_clean$brain_age_gap <- dt_test_clean$predicted_brain_age - dt_test_clean$baseline_age

  # 汇总统计
  bag_mean <- mean(dt_test_clean$brain_age_gap, na.rm = TRUE)
  bag_sd   <- stats::sd(dt_test_clean$brain_age_gap, na.rm = TRUE)
  log_message(sprintf("BAG均值 = %.2f, 标准差 = %.2f", bag_mean, bag_sd))

  # 可视化：预测脑龄 vs 实际年龄（散点 + 对角线）
  p_bag <- ggplot2::ggplot(dt_test_clean, ggplot2::aes(x = baseline_age, y = predicted_brain_age)) +
    ggplot2::geom_point(ggplot2::aes(color = brain_age_gap), alpha = 0.6, size = 2) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey40") +
    viridis::scale_color_viridis(option = "C", direction = 1, name = "BAG") +
    ggplot2::labs(title = "Predicted vs Chronological Age", x = "Chronological Age", y = "Predicted Brain Age") +
    theme_pub()

  # 可视化：BAG 分布（密度图）
  p_bag_density <- ggplot2::ggplot(dt_test_clean, ggplot2::aes(x = brain_age_gap)) +
    ggplot2::geom_density(fill = "#3182bd", alpha = 0.3, color = "#08519c") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::labs(title = "Brain Age Gap (BAG) Distribution", x = "BAG", y = "Density") +
    theme_pub()

  # 可视化：若存在分组列（如 diagnosis/cohort），绘制分组小提琴+箱线图
  p_bag_group <- NULL
  p_bag_density_group <- NULL
  p_pred_group <- NULL
  bag_group_stats <- NULL
  slope_diff_df <- NULL
  if (length(group_col) > 0) {
    gcol <- group_col[1]
    normalize_group <- function(v) {
      vv_raw <- v
      vv_chr <- as.character(vv_raw)
      vv_lower <- tolower(vv_chr)
      if (all(vv_chr %in% c("0","1"), na.rm = TRUE)) {
        return(ifelse(vv_chr == "1", "Ischemic", "Control"))
      }
      dplyr::case_when(
        vv_lower %in% c("control","controls","ctrl") ~ "Control",
        vv_lower %in% c("mi","ischemic","ischemic heart disease","ihd","ischemic_vs_control") ~ "Ischemic",
        vv_lower %in% c("chronic","chronic disease") ~ "Chronic",
        TRUE ~ tools::toTitleCase(vv_chr)
      )
    }
    dt_test_clean$group_label <- normalize_group(dt_test_clean[[gcol]])

    # 分组 BAG 小提琴+箱线图，标注组均值
    p_bag_group <- ggplot2::ggplot(dt_test_clean, ggplot2::aes(x = group_label, y = brain_age_gap, fill = group_label)) +
      ggplot2::geom_violin(trim = FALSE, alpha = 0.3) +
      ggplot2::geom_boxplot(width = 0.2, outlier.shape = NA) +
      ggplot2::stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
      ggplot2::labs(title = "BAG by Group", x = "Group", y = "BAG (Predicted - Chronological)") +
      theme_pub() +
      ggplot2::theme(legend.position = "none")

    # 分组密度叠加，突出分布整体位移
    p_bag_density_group <- ggplot2::ggplot(dt_test_clean, ggplot2::aes(x = brain_age_gap, fill = group_label, color = group_label)) +
      ggplot2::geom_density(alpha = 0.25) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
      ggplot2::labs(title = "BAG Distribution by Group", x = "BAG", y = "Density") +
      theme_pub()

    # 预测脑龄 vs 实际年龄（分组配色 + 组别线性趋势），用于体现加速衰老的斜率差
    p_pred_group <- ggplot2::ggplot(dt_test_clean, ggplot2::aes(x = baseline_age, y = predicted_brain_age, color = group_label)) +
      ggplot2::geom_point(alpha = 0.6, size = 2) +
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey40") +
      ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
      ggplot2::labs(title = "Predicted vs Chronological Age (by Group)", x = "Chronological Age", y = "Predicted Brain Age") +
      theme_pub()

    # 组别 BAG 统计与效应量（相对 Control）
    bag_group_stats <- dt_test_clean %>%
      dplyr::group_by(group_label) %>%
      dplyr::summarise(
        mean_BAG = mean(brain_age_gap, na.rm = TRUE),
        sd_BAG = stats::sd(brain_age_gap, na.rm = TRUE),
        n = dplyr::n(),
        .groups = "drop"
      )
    if (any(bag_group_stats$group_label == "Control")) {
      control_bag <- dt_test_clean$brain_age_gap[dt_test_clean$group_label == "Control"]
      comp_groups <- setdiff(bag_group_stats$group_label, "Control")
      comp_rows <- lapply(comp_groups, function(g) {
        g_bag <- dt_test_clean$brain_age_gap[dt_test_clean$group_label == g]
        if (length(g_bag) >= 3 && length(control_bag) >= 3) {
          tt <- tryCatch(stats::t.test(g_bag, control_bag), error = function(e) NULL)
          d <- tryCatch({
            m1 <- mean(g_bag, na.rm = TRUE); m0 <- mean(control_bag, na.rm = TRUE)
            s1 <- stats::sd(g_bag, na.rm = TRUE); s0 <- stats::sd(control_bag, na.rm = TRUE)
            n1 <- sum(!is.na(g_bag)); n0 <- sum(!is.na(control_bag))
            sp <- sqrt(((n1 - 1)*s1^2 + (n0 - 1)*s0^2)/(n1 + n0 - 2))
            (m1 - m0)/sp
          }, error = function(e) NA_real_)
          data.frame(group_label = g,
                     diff_vs_control = mean(g_bag, na.rm = TRUE) - mean(control_bag, na.rm = TRUE),
                     p_value = if (!is.null(tt)) tt$p.value else NA_real_,
                     cohen_d = d,
                     stringsAsFactors = FALSE)
        } else {
          data.frame(group_label = g, diff_vs_control = NA_real_, p_value = NA_real_, cohen_d = NA_real_)
        }
      })
      bag_group_stats <- bag_group_stats %>% dplyr::left_join(do.call(rbind, comp_rows), by = "group_label")
    }

    # 年龄×组别斜率差（Predicted ~ Age）用于刻画加速衰老
    slope_diff_df <- dt_test_clean %>%
      dplyr::group_by(group_label) %>%
      dplyr::summarise(
        slope = {
          fit <- tryCatch(stats::lm(predicted_brain_age ~ baseline_age, data = cur_data_all()), error = function(e) NULL)
          if (is.null(fit)) NA_real_ else stats::coef(fit)["baseline_age"]
        },
        .groups = "drop"
      )
    if (any(slope_diff_df$group_label == "Control")) {
      s0 <- slope_diff_df$slope[slope_diff_df$group_label == "Control"][1]
      slope_diff_df$diff_vs_control <- slope_diff_df$slope - s0
    }
    bag_group_reg <- NULL
    if (any(dt_test_clean$group_label == "Control")) {
      df_lm <- dt_test_clean[dt_test_clean$group_label %in% c("Control","Ischemic"), , drop = FALSE]
      df_lm <- df_lm[is.finite(df_lm$brain_age_gap) & is.finite(df_lm$baseline_age), , drop = FALSE]
      if (nrow(df_lm) >= 10) {
        df_lm$group_binary <- ifelse(df_lm$group_label == "Ischemic", 1L, 0L)
        fit_lm <- tryCatch(stats::lm(brain_age_gap ~ group_binary + baseline_age, data = df_lm), error = function(e) NULL)
        if (!is.null(fit_lm)) {
          summ <- summary(fit_lm)
          coefs <- as.data.frame(coef(summ))
          if ("group_binary" %in% rownames(coefs)) {
            est <- coefs["group_binary","Estimate"]
            se <- coefs["group_binary","Std. Error"]
            tval <- coefs["group_binary","t value"]
            pval <- coefs["group_binary","Pr(>|t|)"]
            ci <- tryCatch(stats::confint(fit_lm, "group_binary", level = 0.95), error = function(e) matrix(c(NA_real_, NA_real_), nrow = 1))
            bag_group_reg <- data.frame(
              term = "Ischemic_vs_Control",
              estimate = est,
              std_error = se,
              t_value = tval,
              p_value = pval,
              conf_low = ci[1,1],
              conf_high = ci[1,2],
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
    bag_group_reg_multi <- NULL
    if (any(dt_test_clean$group_label == "Control")) {
      df_pred <- dt_test_clean
      df_cov <- if (inherits(dt_test, "data.table")) as.data.frame(dt_test) else dt_test
      key_cols <- intersect(c("eid"), names(df_pred))
      if (length(key_cols) > 0) {
        df_lm_m <- merge(df_pred, df_cov, by = key_cols, all.x = TRUE)
      } else {
        df_pred$.__row_id__ <- seq_len(nrow(df_pred))
        df_cov$.__row_id__ <- seq_len(nrow(df_cov))
        df_lm_m <- merge(df_pred, df_cov, by = ".__row_id__", all.x = TRUE)
      }
      df_lm_m <- df_lm_m[df_lm_m$group_label %in% c("Control","Ischemic"), , drop = FALSE]
      df_lm_m$group_binary <- ifelse(df_lm_m$group_label == "Ischemic", 1L, 0L)
      cov_candidates <- c(
        "baseline_age",
        "sex_factor","ethnicity_factor","imaging_center_factor",
        "bmi_baseline","systolic_bp_baseline","diastolic_bp_baseline",
        "townsend_index_baseline","diabetes_factor_baseline",
        "smoking_factor_baseline","alcohol_factor_baseline"
      )
      cov_present <- intersect(cov_candidates, names(df_lm_m))
      cov_present <- setdiff(cov_present, c("group_binary"))
      if (length(cov_present) > 0 && is.finite(sum(df_lm_m$brain_age_gap, na.rm = TRUE))) {
        for (nm in cov_present) {
          if (is.character(df_lm_m[[nm]]) || is.factor(df_lm_m[[nm]])) {
            df_lm_m[[nm]] <- as.factor(df_lm_m[[nm]])
          } else {
            df_lm_m[[nm]] <- suppressWarnings(as.numeric(df_lm_m[[nm]]))
          }
        }
        rhs <- paste(c("group_binary", cov_present), collapse = " + ")
        fml <- stats::as.formula(paste0("brain_age_gap ~ ", rhs))
        fit_lm2 <- tryCatch(stats::lm(fml, data = df_lm_m), error = function(e) NULL)
        if (!is.null(fit_lm2)) {
          summ2 <- summary(fit_lm2)
          coefs2 <- as.data.frame(coef(summ2))
          if ("group_binary" %in% rownames(coefs2)) {
            est <- coefs2["group_binary","Estimate"]
            se <- coefs2["group_binary","Std. Error"]
            tval <- coefs2["group_binary","t value"]
            pval <- coefs2["group_binary","Pr(>|t|)"]
            ci2 <- tryCatch(stats::confint(fit_lm2, "group_binary", level = 0.95), error = function(e) matrix(c(NA_real_, NA_real_), nrow = 1))
            bag_group_reg_multi <- data.frame(
              term = "Ischemic_vs_Control_adjusted",
              estimate = est,
              std_error = se,
              t_value = tval,
              p_value = pval,
              conf_low = ci2[1,1],
              conf_high = ci2[1,2],
              n = nrow(df_lm_m),
              covariates_included = paste(cov_present, collapse = ";"),
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }

  # 特征重要性（Top 20）
  imp <- tryCatch(randomForest::importance(rf_model), error = function(e) NULL)
  p_feat_imp <- NULL
  if (!is.null(imp)) {
    imp_df <- data.frame(feature = rownames(imp), importance = imp[, 1], stringsAsFactors = FALSE)
    imp_df <- imp_df[order(imp_df$importance, decreasing = TRUE), ]
    top_k <- min(20, nrow(imp_df))
    imp_top <- imp_df[1:top_k, ]
    p_feat_imp <- ggplot2::ggplot(imp_top, ggplot2::aes(x = reorder(feature, importance), y = importance, fill = importance)) +
      ggplot2::geom_col() + ggplot2::coord_flip() +
      viridis::scale_fill_viridis(option = "C", direction = 1, name = "Importance") +
      ggplot2::labs(title = "Feature Importance (Random Forest)", x = "Feature", y = "Importance") +
      theme_pub()
  }

  # 保存图与数据
  if (!is.null(p_bag))       save_pub(file.path(output_dir, "BAG_Predicted_vs_Actual"), p_bag, width = 8, height = 6, dpi = 600)
  if (!is.null(p_bag_density)) save_pub(file.path(output_dir, "BAG_Distribution"), p_bag_density, width = 8, height = 6, dpi = 600)
  if (!is.null(p_bag_group))  save_pub(file.path(output_dir, "BAG_By_Group"), p_bag_group, width = 8, height = 6, dpi = 600)
  if (!is.null(p_bag_density_group)) save_pub(file.path(output_dir, "BAG_Distribution_By_Group"), p_bag_density_group, width = 8, height = 6, dpi = 600)
  if (!is.null(p_pred_group)) save_pub(file.path(output_dir, "BAG_Predicted_vs_Actual_By_Group"), p_pred_group, width = 8, height = 6, dpi = 600)
  if (!is.null(p_feat_imp))   save_pub(file.path(output_dir, "BAG_Feature_Importance"), p_feat_imp, width = 8, height = 6, dpi = 600)

  utils::write.csv(dt_test_clean, file.path(output_dir, "BAG_Predictions.csv"), row.names = FALSE)
  bag_summary <- data.frame(BAG_mean = bag_mean, BAG_sd = bag_sd, n = nrow(dt_test_clean))
  utils::write.csv(bag_summary, file.path(output_dir, "BAG_Summary.csv"), row.names = FALSE)
  if (!is.null(bag_group_stats)) utils::write.csv(bag_group_stats, file.path(output_dir, "BAG_Group_Stats.csv"), row.names = FALSE)
  if (!is.null(slope_diff_df)) utils::write.csv(slope_diff_df, file.path(output_dir, "BAG_Age_Slope_Difference.csv"), row.names = FALSE)
  if (exists("bag_group_reg") && !is.null(bag_group_reg)) utils::write.csv(bag_group_reg, file.path(output_dir, "BAG_Group_Regression.csv"), row.names = FALSE)
  if (exists("bag_group_reg_multi") && !is.null(bag_group_reg_multi)) utils::write.csv(bag_group_reg_multi, file.path(output_dir, "BAG_Group_Regression_Multivariable.csv"), row.names = FALSE)

  log_message("脑龄预测分析完成，结果已保存至:", output_dir)
  return(list(
    model = rf_model,
    predictions = dt_test_clean,
    bag_summary = bag_summary,
    plots = list(pred_vs_actual = p_bag, bag_density = p_bag_density, bag_group = p_bag_group, feat_importance = p_feat_imp)
  ))
}

# 便捷入口：从匹配队列文件加载数据并运行BAG分析
run_brain_age_gap_for <- function(train_cohort_csv, test_cohort_csv, idp_list = NULL, output_dir = "Brain_Age_Analysis") {
  if (!file.exists(train_cohort_csv) || !file.exists(test_cohort_csv)) {
    stop("指定的训练/测试队列CSV不存在")
  }
  train_df <- tryCatch(utils::read.csv(train_cohort_csv, stringsAsFactors = FALSE), error = function(e) NULL)
  test_df  <- tryCatch(utils::read.csv(test_cohort_csv, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(train_df) || is.null(test_df)) stop("无法读取队列CSV")

  # 自动选择IDPs（若未提供）
  if (is.null(idp_list) || length(idp_list) == 0) idp_list <- get_brain_age_idps(top_n = 64)
  if (length(idp_list) == 0) stop("未能解析到可用于脑龄预测的IDP变量名")
  
  idp_list <- unique(as.character(idp_list))
  idp_list <- idp_list[!is.na(idp_list) & nzchar(idp_list)]
  idp_list <- intersect(idp_list, intersect(names(train_df), names(test_df)))
  if (length(idp_list) < 5) stop("训练/测试队列可用的脑龄特征列不足")

  # 仅保留建模所需列（减少内存与NA）
  keep_train <- intersect(c("baseline_age", idp_list), names(train_df))
  cov_candidates <- c(
    "sex_factor","ethnicity_factor","imaging_center_factor",
    "bmi_baseline","systolic_bp_baseline","diastolic_bp_baseline",
    "townsend_index_baseline","diabetes_factor_baseline",
    "smoking_factor_baseline","alcohol_factor_baseline"
  )
  keep_test  <- intersect(
    c("baseline_age", idp_list, "diagnosis","cohort","group","case_control","eid", cov_candidates),
    names(test_df)
  )
  train_df <- train_df[, keep_train, drop = FALSE]
  test_df  <- test_df[, keep_test, drop = FALSE]

  create_brain_age_gap_analysis(train_df, test_df, idp_list, output_dir = output_dir)
}

# 运行脑龄分析
run_brain_age_gap_for(
   train_cohort_csv = table_out_path(file.path("Cohort", "Matched"), "Final_Matched_Cohort_ischemic_vs_control.csv"),
   test_cohort_csv  = table_out_path(file.path("Cohort", "Matched"), "Final_Matched_Cohort_mi_vs_control.csv"),
   idp_list = NULL,
   output_dir = get_table_dir("Brain_Age_Analysis")
 )
