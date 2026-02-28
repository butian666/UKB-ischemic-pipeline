options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
})

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
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = base_size + 2, color = "#111111"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = base_size, color = "#444444"),
      plot.caption = ggplot2::element_text(hjust = 0, size = base_size - 2, color = "#555555"),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold"),
      legend.key.width = grid::unit(10, "pt"),
      legend.key.height = grid::unit(10, "pt")
    )
}

theme_longitudinal_pub <- function(base_size = 12, base_family = "Helvetica") {
  theme_jama_refined(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      legend.position = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.background = ggplot2::element_rect(fill = scales::alpha("white", 0.75), color = "grey80"),
      legend.key = ggplot2::element_rect(fill = NA, color = NA),
      strip.text = ggplot2::element_text(size = base_size - 2, face = "bold"),
      plot.margin = ggplot2::margin(6, 10, 6, 6)
    )
}

check_dir_writable <- function(dir_path) {
  tryCatch({
    if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    tf <- tempfile(pattern = ".write_test_", tmpdir = dir_path)
    ok <- suppressWarnings(file.create(tf))
    if (isTRUE(ok)) unlink(tf)
    isTRUE(ok)
  }, error = function(e) FALSE)
}

save_pub <- function(filename, plot, width = 10, height = 6, dpi = 600, formats = c("png", "pdf")) {
  base <- sub("\n$", "", filename)
  dir_path <- dirname(base)
  if (!check_dir_writable(dir_path)) stop(paste0("output dir not writable: ", dir_path))
  for (fmt in formats) {
    out <- paste0(tools::file_path_sans_ext(base), ".", fmt)
    ggplot2::ggsave(out, plot, width = width, height = height, dpi = dpi, units = "in")
  }
}

get_group_colors <- function() {
  c(
    "Controls" = "#1F78B4",
    "Cases" = "#FF7F00",
    "Control" = "#1F78B4",
    "Case" = "#FF7F00",
    "Ischemic_Heart_Disease" = "#FF7F00",
    "Myocardial_Infarction" = "#E31A1C",
    "Chronic_Ischemic" = "#33A02C"
  )
}

safe_filename <- function(x, max_len = 80) {
  x <- as.character(x)
  x <- ifelse(is.na(x) | !nzchar(x), "NA", x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  ifelse(nchar(x) > max_len, substr(x, 1, max_len), x)
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  if (is.logical(x)) return(isTRUE(x))
  x <- tolower(trimws(as.character(x)))
  if (x %in% c("1", "true", "t", "yes", "y")) return(TRUE)
  if (x %in% c("0", "false", "f", "no", "n")) return(FALSE)
  default
}

parse_args <- function(args) {
  out <- list()
  if (length(args) == 0) return(out)
  for (a in args) {
    if (!startsWith(a, "--")) next
    kv <- sub("^--", "", a)
    key <- sub("=.*$", "", kv)
    val <- sub("^[^=]*=?", "", kv)
    if (identical(key, val)) val <- TRUE
    out[[key]] <- val
  }
  out
}

to_numeric_safely <- function(v) {
  if (is.factor(v)) v <- as.character(v)
  if (is.logical(v)) return(suppressWarnings(as.numeric(v)))
  if (is.numeric(v)) return(suppressWarnings(as.numeric(v)))
  v <- trimws(as.character(v))
  v[v %in% c("", "NA", "NaN", "NULL", "null")] <- NA
  suppressWarnings(as.numeric(v))
}

root_from_name <- function(x) gsub("(\\.{1,}|_)Instance\\.?[23]([._].*)?$", "", x)

find_age_col <- function(dat, preferred = NULL) {
  if (!is.null(preferred) && preferred %in% names(dat)) return(preferred)
  for (cand in c("baseline_age", "Age2", "age2", "age", "Age.when.attended.assessment.centre...Instance.2")) {
    if (cand %in% names(dat)) return(cand)
  }
  NA_character_
}

summarize_percent <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(list(n = NA_integer_, median = NA_real_, mean = NA_real_, q25 = NA_real_, q75 = NA_real_))
  }
  list(
    n = length(x),
    median = stats::median(x, na.rm = TRUE),
    mean = base::mean(x, na.rm = TRUE),
    q25 = as.numeric(stats::quantile(x, 0.25, na.rm = TRUE, type = 7)),
    q75 = as.numeric(stats::quantile(x, 0.75, na.rm = TRUE, type = 7))
  )
}

age_slope <- function(age, pct) {
  idx <- is.finite(age) & is.finite(pct)
  if (sum(idx) < 2) return(NA_real_)
  fit <- stats::lm(pct[idx] ~ age[idx])
  as.numeric(stats::coef(fit)[2])
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

model_label <- if (!is.null(args$model)) as.character(args$model) else "mi_vs_chronic"
combined_csv <- if (!is.null(args$combined)) as.character(args$combined) else file.path("Output_Tables", "Four_Model", "Combined_Four_Model_Z_Analysis_Results.csv")
out_dir <- if (!is.null(args$out)) as.character(args$out) else file.path("Output_Tables", "Longitudinal", "Median_DeltaRate")
plot_dir <- if (!is.null(args$plot_out)) as.character(args$plot_out) else file.path("Output_Tables", "Longitudinal", "Plots_Count_Cognitive", model_label)
plot_dpi <- if (!is.null(args$plot_dpi)) suppressWarnings(as.numeric(args$plot_dpi)) else 600
if (!is.finite(plot_dpi) || plot_dpi < 300) plot_dpi <- 600
make_plots <- as_bool(args$plots, default = FALSE)
idp_filter <- if (!is.null(args$idp)) as.character(args$idp) else NULL
smooth_k <- if (!is.null(args$smooth_k)) suppressWarnings(as.numeric(args$smooth_k)) else 1
if (!is.finite(smooth_k) || smooth_k < 0) smooth_k <- 1

cohort_rds <- if (!is.null(args$cohort)) {
  as.character(args$cohort)
} else {
  file.path("Output_Tables", "Cohort", "Matched", paste0("Final_Matched_Cohort_", model_label, ".rds"))
}

if (!file.exists(combined_csv)) stop(paste0("missing combined csv: ", combined_csv))
if (!file.exists(cohort_rds)) stop(paste0("missing cohort rds: ", cohort_rds))
if (!dir.exists(out_dir)) stop(paste0("missing output dir: ", out_dir))

comb <- utils::read.csv(combined_csv, check.names = FALSE)
sub_df <- comb[comb$model_comparison == model_label & comb$variable_type %in% c("Brain_Imaging", "Cognitive"), , drop = FALSE]
if (nrow(sub_df) == 0) stop(paste0("no rows for model: ", model_label))
if (!all(c("idp1_variable", "idp2_variable") %in% names(sub_df))) stop("combined results missing idp1_variable/idp2_variable")

sub_df$idp_root <- root_from_name(sub_df$idp1_variable)

get_num <- function(df, nm) if (nm %in% names(df)) suppressWarnings(as.numeric(df[[nm]])) else rep(NA_real_, nrow(df))
z_key <- abs(coalesce(get_num(sub_df, "z_statistic_without_non_imaging"), get_num(sub_df, "z_statistic_with_non_imaging")))
p_key <- coalesce(get_num(sub_df, "p_value_without_non_imaging"), get_num(sub_df, "p_value_with_non_imaging"))

meta <- sub_df %>%
  mutate(.z_key = z_key, .p_key = p_key) %>%
  arrange(desc(.z_key), .p_key) %>%
  group_by(.data$idp_root) %>%
  slice(1) %>%
  ungroup() %>%
  select(
    idp_root,
    model_comparison,
    variable_name,
    variable_type,
    idp_category,
    idp1_variable,
    idp2_variable,
    z_statistic_without_non_imaging,
    z_statistic_with_non_imaging,
    p_value_without_non_imaging,
    p_value_with_non_imaging
  )

if (!is.null(idp_filter)) {
  meta <- meta[meta$idp_root %in% idp_filter | meta$variable_name %in% idp_filter, , drop = FALSE]
}

dat <- readRDS(cohort_rds)
if (!("group" %in% names(dat))) stop("cohort missing group")
age_col <- find_age_col(dat, preferred = if (!is.null(args$age_col)) as.character(args$age_col) else NULL)
age_vec <- if (!is.na(age_col)) to_numeric_safely(dat[[age_col]]) else rep(NA_real_, nrow(dat))

is_count_like <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(FALSE)
  idx <- is.finite(x)
  x <- x[idx]
  if (length(x) == 0) return(FALSE)
  is_int_nonneg <- abs(x - round(x)) < 1e-8 & x >= 0
  mean(is_int_nonneg, na.rm = TRUE) >= 0.95
}

res_list <- lapply(seq_len(nrow(meta)), function(i) {
  c1 <- meta$idp1_variable[i]
  c2 <- meta$idp2_variable[i]
  if (!(c1 %in% names(dat)) || !(c2 %in% names(dat))) return(NULL)
  x1 <- to_numeric_safely(dat[[c1]])
  x2 <- to_numeric_safely(dat[[c2]])
  diffv <- suppressWarnings(x2 - x1)
  abs_change <- ifelse(is.finite(diffv), diffv, NA_real_)
  denom_ind <- ifelse(is.finite(x1) & abs(x1) > .Machine$double.eps, x1, NA_real_)
  pct_change <- ifelse(is.finite(diffv) & is.finite(denom_ind), (diffv / denom_ind) * 100, NA_real_)
  count_like <- is_count_like(x1) && is_count_like(x2)
  smooth_denom <- ifelse(count_like & is.finite(x1), x1 + smooth_k, NA_real_)
  smooth_pct_change <- ifelse(count_like & is.finite(diffv) & is.finite(smooth_denom), (diffv / smooth_denom) * 100, NA_real_)

  g <- dat$group
  ctrl_idx <- g == 0
  case_idx <- g == 1

  s_ctrl_pct <- summarize_percent(pct_change[ctrl_idx])
  s_case_pct <- summarize_percent(pct_change[case_idx])
  s_ctrl_abs <- summarize_percent(abs_change[ctrl_idx])
  s_case_abs <- summarize_percent(abs_change[case_idx])
  s_ctrl_smooth <- summarize_percent(smooth_pct_change[ctrl_idx])
  s_case_smooth <- summarize_percent(smooth_pct_change[case_idx])

  slope_ctrl_pct <- age_slope(age_vec[ctrl_idx], pct_change[ctrl_idx])
  slope_case_pct <- age_slope(age_vec[case_idx], pct_change[case_idx])
  slope_ctrl_abs <- age_slope(age_vec[ctrl_idx], abs_change[ctrl_idx])
  slope_case_abs <- age_slope(age_vec[case_idx], abs_change[case_idx])
  slope_ctrl_smooth <- age_slope(age_vec[ctrl_idx], smooth_pct_change[ctrl_idx])
  slope_case_smooth <- age_slope(age_vec[case_idx], smooth_pct_change[case_idx])

  preferred_metric <- if (count_like) "abs_change" else "percent_change"
  if (count_like) {
    pref_ctrl <- s_ctrl_abs
    pref_case <- s_case_abs
    pref_slope_ctrl <- slope_ctrl_abs
    pref_slope_case <- slope_case_abs
  } else {
    pref_ctrl <- s_ctrl_pct
    pref_case <- s_case_pct
    pref_slope_ctrl <- slope_ctrl_pct
    pref_slope_case <- slope_case_pct
  }

  data.frame(
    idp_root = meta$idp_root[i],
    model_comparison = meta$model_comparison[i],
    variable_name = meta$variable_name[i],
    variable_type = meta$variable_type[i],
    idp_category = meta$idp_category[i],
    z_statistic = coalesce(meta$z_statistic_without_non_imaging[i], meta$z_statistic_with_non_imaging[i]),
    p_value = coalesce(meta$p_value_without_non_imaging[i], meta$p_value_with_non_imaging[i]),
    is_count_like = count_like,
    preferred_metric = preferred_metric,
    control_n_percent = s_ctrl_pct$n,
    control_median_percent = s_ctrl_pct$median,
    control_mean_percent = s_ctrl_pct$mean,
    control_q25_percent = s_ctrl_pct$q25,
    control_q75_percent = s_ctrl_pct$q75,
    case_n_percent = s_case_pct$n,
    case_median_percent = s_case_pct$median,
    case_mean_percent = s_case_pct$mean,
    case_q25_percent = s_case_pct$q25,
    case_q75_percent = s_case_pct$q75,
    median_diff_percent_case_minus_control = s_case_pct$median - s_ctrl_pct$median,
    mean_diff_percent_case_minus_control = s_case_pct$mean - s_ctrl_pct$mean,
    control_n_abs = s_ctrl_abs$n,
    control_median_abs = s_ctrl_abs$median,
    control_mean_abs = s_ctrl_abs$mean,
    control_q25_abs = s_ctrl_abs$q25,
    control_q75_abs = s_ctrl_abs$q75,
    case_n_abs = s_case_abs$n,
    case_median_abs = s_case_abs$median,
    case_mean_abs = s_case_abs$mean,
    case_q25_abs = s_case_abs$q25,
    case_q75_abs = s_case_abs$q75,
    median_diff_abs_case_minus_control = s_case_abs$median - s_ctrl_abs$median,
    mean_diff_abs_case_minus_control = s_case_abs$mean - s_ctrl_abs$mean,
    control_n_smooth_percent = s_ctrl_smooth$n,
    control_median_smooth_percent = s_ctrl_smooth$median,
    control_mean_smooth_percent = s_ctrl_smooth$mean,
    control_q25_smooth_percent = s_ctrl_smooth$q25,
    control_q75_smooth_percent = s_ctrl_smooth$q75,
    case_n_smooth_percent = s_case_smooth$n,
    case_median_smooth_percent = s_case_smooth$median,
    case_mean_smooth_percent = s_case_smooth$mean,
    case_q25_smooth_percent = s_case_smooth$q25,
    case_q75_smooth_percent = s_case_smooth$q75,
    median_diff_smooth_percent_case_minus_control = s_case_smooth$median - s_ctrl_smooth$median,
    mean_diff_smooth_percent_case_minus_control = s_case_smooth$mean - s_ctrl_smooth$mean,
    control_age_slope_percent = slope_ctrl_pct,
    case_age_slope_percent = slope_case_pct,
    age_slope_diff_percent_case_minus_control = slope_case_pct - slope_ctrl_pct,
    control_age_slope_abs = slope_ctrl_abs,
    case_age_slope_abs = slope_case_abs,
    age_slope_diff_abs_case_minus_control = slope_case_abs - slope_ctrl_abs,
    control_age_slope_smooth_percent = slope_ctrl_smooth,
    case_age_slope_smooth_percent = slope_case_smooth,
    age_slope_diff_smooth_percent_case_minus_control = slope_case_smooth - slope_ctrl_smooth,
    preferred_control_n = pref_ctrl$n,
    preferred_control_median = pref_ctrl$median,
    preferred_control_mean = pref_ctrl$mean,
    preferred_case_n = pref_case$n,
    preferred_case_median = pref_case$median,
    preferred_case_mean = pref_case$mean,
    preferred_median_diff_case_minus_control = pref_case$median - pref_ctrl$median,
    preferred_mean_diff_case_minus_control = pref_case$mean - pref_ctrl$mean,
    preferred_control_age_slope = pref_slope_ctrl,
    preferred_case_age_slope = pref_slope_case,
    preferred_age_slope_diff_case_minus_control = pref_slope_case - pref_slope_ctrl,
    stringsAsFactors = FALSE
  )
})

res <- bind_rows(res_list)

out_csv <- file.path(out_dir, paste0("Longitudinal_IndivPercentChange_Summary_", model_label, ".csv"))
utils::write.csv(res, out_csv, row.names = FALSE)

if (nrow(res) > 0) {
  top <- res %>% arrange(desc(preferred_age_slope_diff_case_minus_control)) %>% head(10)
  print(top[, c("idp_root", "preferred_metric", "preferred_age_slope_diff_case_minus_control", "preferred_case_age_slope", "preferred_control_age_slope")])
}

if (isTRUE(make_plots)) {
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  meta_cog <- meta %>% filter(.data$variable_type == "Cognitive")
  if (nrow(meta_cog) == 0) {
    message("no cognitive rows for plotting")
  } else {
    group_colors <- get_group_colors()
    per_idp_dir <- file.path(plot_dir, "Per_IDP")
    if (!dir.exists(per_idp_dir)) dir.create(per_idp_dir, recursive = TRUE, showWarnings = FALSE)

    build_plot_df_one <- function(i) {
      c1 <- meta_cog$idp1_variable[i]
      c2 <- meta_cog$idp2_variable[i]
      if (!(c1 %in% names(dat)) || !(c2 %in% names(dat))) return(NULL)
      x1 <- to_numeric_safely(dat[[c1]])
      x2 <- to_numeric_safely(dat[[c2]])
      diffv <- suppressWarnings(x2 - x1)
      count_like <- is_count_like(x1) && is_count_like(x2)
      if (!isTRUE(count_like)) return(NULL)
      abs_change <- ifelse(is.finite(diffv), diffv, NA_real_)
      smooth_denom <- ifelse(is.finite(x1), x1 + smooth_k, NA_real_)
      smooth_pct_change <- ifelse(is.finite(diffv) & is.finite(smooth_denom), (diffv / smooth_denom) * 100, NA_real_)
      g <- dat$group
      grp <- dplyr::case_when(
        g == 0 ~ "Controls",
        g == 1 ~ "Cases",
        TRUE ~ NA_character_
      )
      data.frame(
        model_comparison = model_label,
        idp_root = meta_cog$idp_root[i],
        variable_name = meta_cog$variable_name[i],
        idp_category = meta_cog$idp_category[i],
        Group = factor(grp, levels = c("Controls", "Cases")),
        Age = age_vec,
        abs_change = abs_change,
        smooth_pct_change = smooth_pct_change,
        stringsAsFactors = FALSE
      ) %>%
        filter(!is.na(.data$Group), is.finite(.data$Age))
    }

    plot_df <- dplyr::bind_rows(lapply(seq_len(nrow(meta_cog)), build_plot_df_one))
    plot_df <- plot_df %>%
      mutate(
        variable_name = as.character(.data$variable_name),
        variable_name = ifelse(is.na(variable_name) | !nzchar(variable_name), .data$idp_root, variable_name)
      ) %>%
      filter(is.finite(.data$abs_change) | is.finite(.data$smooth_pct_change))

    if (nrow(plot_df) == 0) {
      message("no count-type cognitive rows with usable values for plotting")
    } else {
      age_min <- suppressWarnings(min(plot_df$Age, na.rm = TRUE))
      age_max <- suppressWarnings(max(plot_df$Age, na.rm = TRUE))
      age_breaks <- if (is.finite(age_min) && is.finite(age_max)) {
        seq(floor(age_min / 5) * 5, ceiling(age_max / 5) * 5, by = 5)
      } else {
        scales::pretty_breaks(n = 6)(plot_df$Age)
      }
      key_ages <- c(60, 70)
      key_ages <- key_ages[key_ages >= min(age_breaks, na.rm = TRUE) & key_ages <= max(age_breaks, na.rm = TRUE)]
      key_df <- if (length(key_ages) > 0) data.frame(Age = key_ages, label = paste0(key_ages, "y")) else data.frame(Age = numeric(0), label = character(0))

      idps_to_plot <- unique(plot_df$idp_root)
      for (idp in idps_to_plot) {
        df_one <- plot_df %>% filter(.data$idp_root == idp)
        if (nrow(df_one) == 0) next
        vtitle <- df_one$variable_name[!is.na(df_one$variable_name) & nzchar(df_one$variable_name)][1]
        if (is.na(vtitle) || !nzchar(vtitle)) vtitle <- idp
        safe_idp <- safe_filename(idp, max_len = 80)
        base_name <- paste0(model_label, "__", safe_idp)

        p_box_abs <- ggplot2::ggplot(df_one %>% filter(is.finite(.data$abs_change)), ggplot2::aes(x = Group, y = abs_change, fill = Group)) +
          ggplot2::geom_boxplot(alpha = 0.75, width = 0.65, color = "black", linewidth = 0.4, outlier.alpha = 0.8, outlier.size = 1.2) +
          ggplot2::geom_jitter(width = 0.18, alpha = 0.25, size = 1.0, color = "black") +
          ggplot2::scale_fill_manual(values = group_colors, name = "") +
          ggplot2::labs(
            title = vtitle,
            subtitle = paste0(model_label, " | Absolute Change (IDP2 - IDP1)"),
            x = NULL,
            y = "Absolute Change"
          ) +
          theme_longitudinal_pub() +
          ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(alpha = 1)))

        save_pub(file.path(per_idp_dir, paste0("Boxplot_Absolute__", base_name)), p_box_abs, width = 6, height = 6, dpi = plot_dpi)

        p_box_smooth <- ggplot2::ggplot(df_one %>% filter(is.finite(.data$smooth_pct_change)), ggplot2::aes(x = Group, y = smooth_pct_change, fill = Group)) +
          ggplot2::geom_boxplot(alpha = 0.75, width = 0.65, color = "black", linewidth = 0.4, outlier.alpha = 0.8, outlier.size = 1.2) +
          ggplot2::geom_jitter(width = 0.18, alpha = 0.25, size = 1.0, color = "black") +
          ggplot2::scale_fill_manual(values = group_colors, name = "") +
          ggplot2::labs(
            title = vtitle,
            subtitle = paste0(model_label, " | Smoothed Relative Change (%) = (IDP2-IDP1)/(IDP1+", smooth_k, ") x 100"),
            x = NULL,
            y = "Smoothed Delta (%)"
          ) +
          theme_longitudinal_pub() +
          ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(alpha = 1)))

        save_pub(file.path(per_idp_dir, paste0("Boxplot_SmoothedPct__", base_name)), p_box_smooth, width = 6, height = 6, dpi = plot_dpi)

        p_trend_abs <- ggplot2::ggplot(df_one %>% filter(is.finite(.data$abs_change)), ggplot2::aes(x = Age, y = abs_change, color = Group)) +
          ggplot2::geom_point(alpha = 0.25, size = 1.0) +
          ggplot2::geom_smooth(method = "loess", se = TRUE, level = 0.95, span = 0.8, alpha = 0.18, linewidth = 1.0) +
          { if (nrow(key_df) > 0) ggplot2::geom_vline(data = key_df, ggplot2::aes(xintercept = Age), inherit.aes = FALSE, linetype = "dashed", color = "grey50", linewidth = 0.35) else ggplot2::geom_blank() } +
          { if (nrow(key_df) > 0) ggplot2::geom_text(data = key_df, ggplot2::aes(x = Age, y = Inf, label = label), inherit.aes = FALSE, vjust = 1.4, size = 3, color = "grey30") else ggplot2::geom_blank() } +
          ggplot2::scale_color_manual(values = group_colors, name = "") +
          ggplot2::scale_x_continuous(breaks = age_breaks) +
          ggplot2::labs(
            title = vtitle,
            subtitle = paste0(model_label, " | Absolute Change vs Age"),
            x = "Age (years)",
            y = "Absolute Change"
          ) +
          theme_longitudinal_pub() +
          ggplot2::coord_cartesian(clip = "off") +
          ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1, linewidth = 1.2)))

        save_pub(file.path(per_idp_dir, paste0("AgeTrend_Absolute__", base_name)), p_trend_abs, width = 6, height = 6, dpi = plot_dpi)

        p_trend_smooth <- ggplot2::ggplot(df_one %>% filter(is.finite(.data$smooth_pct_change)), ggplot2::aes(x = Age, y = smooth_pct_change, color = Group)) +
          ggplot2::geom_point(alpha = 0.25, size = 1.0) +
          ggplot2::geom_smooth(method = "loess", se = TRUE, level = 0.95, span = 0.8, alpha = 0.18, linewidth = 1.0) +
          { if (nrow(key_df) > 0) ggplot2::geom_vline(data = key_df, ggplot2::aes(xintercept = Age), inherit.aes = FALSE, linetype = "dashed", color = "grey50", linewidth = 0.35) else ggplot2::geom_blank() } +
          { if (nrow(key_df) > 0) ggplot2::geom_text(data = key_df, ggplot2::aes(x = Age, y = Inf, label = label), inherit.aes = FALSE, vjust = 1.4, size = 3, color = "grey30") else ggplot2::geom_blank() } +
          ggplot2::scale_color_manual(values = group_colors, name = "") +
          ggplot2::scale_x_continuous(breaks = age_breaks) +
          ggplot2::labs(
            title = vtitle,
            subtitle = paste0(model_label, " | Smoothed Relative Change vs Age (IDP1+", smooth_k, ")"),
            x = "Age (years)",
            y = "Smoothed Delta (%)"
          ) +
          theme_longitudinal_pub() +
          ggplot2::coord_cartesian(clip = "off") +
          ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1, linewidth = 1.2)))

        save_pub(file.path(per_idp_dir, paste0("AgeTrend_SmoothedPct__", base_name)), p_trend_smooth, width = 6, height = 6, dpi = plot_dpi)
      }
    }
  }
}
