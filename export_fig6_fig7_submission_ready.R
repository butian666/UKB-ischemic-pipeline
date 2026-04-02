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

clean_lab <- function(x, n = 52L) {
  x <- gsub("\\.+", " ", x)
  x <- gsub("\\s+", " ", trimws(x))
  ifelse(nchar(x) > n, paste0(substr(x, 1, n), "..."), x)
}

in_dir <- kv$in_dir %||% file.path(getwd(), "Output_Tables", "MR_Bidirectional_Bulk", "Manuscript_and_Figures", "Circulation_Enhanced")
out_dir <- kv$out_dir %||% file.path(in_dir, "Submission_Ready")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

steiger_csv <- file.path(in_dir, "Steiger_Directionality_Table.csv")
raps_csv <- file.path(in_dir, "MR_RAPS_Comparison_Table.csv")
if (!file.exists(steiger_csv) || !file.exists(raps_csv)) stop("缺少Steiger或RAPS结果表")

steiger <- read_csv_flex(steiger_csv)
raps <- read_csv_flex(raps_csv)

library(ggplot2)

theme_pub <- function() {
  theme_bw(base_size = 9, base_family = "Arial") +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(linewidth = 0.5, colour = "black"),
      axis.ticks = element_line(linewidth = 0.5, colour = "black"),
      legend.title = element_text(size = 8, face = "plain"),
      legend.text = element_text(size = 8),
      plot.title = element_text(size = 10, face = "bold", hjust = 0),
      strip.background = element_blank(),
      strip.text = element_text(size = 8, face = "bold")
    )
}

save_tiff600 <- function(plot_obj, filename, width, height) {
  ggsave(filename = file.path(out_dir, filename), plot = plot_obj, width = width, height = height, dpi = 600, units = "in", compression = "lzw")
}

steiger$direction_label <- ifelse(steiger$direction == "heart_to_brain", "Heart to Brain", "Brain to Heart")
steiger$decision <- ifelse(steiger$correct_causal_direction %in% TRUE, "Supported", "Not supported")
pal_dir <- c("Heart to Brain" = "#1F4E79", "Brain to Heart" = "#8C2D04")
pal_dec <- c("Supported" = 16, "Not supported" = 17)

p6a <- ggplot(steiger, aes(x = direction_label, fill = direction_label)) +
  geom_bar(linewidth = 0.3, colour = "black") +
  scale_fill_manual(values = pal_dir) +
  labs(
    title = "Figure 6A. Steiger directionality support across prioritized pairs",
    x = "MR direction",
    y = "Number of prioritized pairs",
    fill = "MR direction"
  ) +
  theme_pub()
save_tiff600(p6a, "Fig6A_Steiger_Direction_Bar_600dpi.tiff", width = 6.6, height = 4.2)

p6b <- ggplot(steiger, aes(x = r2_exposure, y = r2_outcome, color = direction_label, shape = decision)) +
  geom_point(size = 2.6, stroke = 0.45) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, linewidth = 0.4, colour = "#636363") +
  scale_color_manual(values = pal_dir) +
  scale_shape_manual(values = pal_dec) +
  labs(
    title = "Figure 6B. Exposure versus outcome variance explained in Steiger testing",
    x = expression(sum(r[exposure]^2)),
    y = expression(sum(r[outcome]^2)),
    color = "MR direction",
    shape = "Steiger decision"
  ) +
  theme_pub()
save_tiff600(p6b, "Fig6B_Steiger_R2_Scatter_600dpi.tiff", width = 6.6, height = 5.0)

raps$direction_label <- ifelse(raps$direction == "heart_to_brain", "Heart to Brain", "Brain to Heart")
raps$pair_label <- paste0(raps$pair_id, " | ", clean_lab(raps$outcome_trait, 36))
raps$ivw_l <- raps$ivw_b - 1.96 * raps$ivw_se
raps$ivw_u <- raps$ivw_b + 1.96 * raps$ivw_se
raps$raps_l <- raps$raps_b - 1.96 * raps$raps_se
raps$raps_u <- raps$raps_b + 1.96 * raps$raps_se

long_beta <- rbind(
  data.frame(pair_label = raps$pair_label, direction_label = raps$direction_label, method = "IVW", beta = raps$ivw_b, lo = raps$ivw_l, hi = raps$ivw_u, stringsAsFactors = FALSE),
  data.frame(pair_label = raps$pair_label, direction_label = raps$direction_label, method = "MR-RAPS", beta = raps$raps_b, lo = raps$raps_l, hi = raps$raps_u, stringsAsFactors = FALSE)
)
long_beta$pair_label <- factor(long_beta$pair_label, levels = rev(unique(raps$pair_label)))

p7a <- ggplot(long_beta, aes(x = beta, y = pair_label, color = method)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.4, colour = "#636363") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.15, linewidth = 0.45, position = position_dodge(width = 0.55), orientation = "y") +
  geom_point(size = 1.9, position = position_dodge(width = 0.55)) +
  facet_wrap(~direction_label, scales = "free_y", ncol = 1) +
  scale_color_manual(values = c("IVW" = "#1B9E77", "MR-RAPS" = "#D95F02")) +
  labs(
    title = "Figure 7A. Effect-size comparison with 95% CI (IVW vs MR-RAPS)",
    x = "Causal effect estimate (beta, 95% CI)",
    y = "Prioritized pair",
    color = "Estimator"
  ) +
  theme_pub()
save_tiff600(p7a, "Fig7A_Effect_Comparison_with95CI_600dpi.tiff", width = 7.2, height = 7.8)

p7b <- ggplot(raps, aes(x = ivw_b, y = raps_b, color = direction_label)) +
  geom_errorbar(aes(ymin = raps_l, ymax = raps_u), linewidth = 0.4, alpha = 0.75) +
  geom_errorbar(aes(xmin = ivw_l, xmax = ivw_u), linewidth = 0.4, alpha = 0.75, orientation = "y") +
  geom_point(size = 2.3, stroke = 0.35) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, linewidth = 0.45, colour = "#636363") +
  scale_color_manual(values = c("Heart to Brain" = "#4D4D4D", "Brain to Heart" = "#1F78B4")) +
  labs(
    title = "Figure 7B. Agreement between IVW and MR-RAPS estimates",
    x = "IVW beta (95% CI)",
    y = "MR-RAPS beta (95% CI)",
    color = "MR direction"
  ) +
  theme_pub()
save_tiff600(p7b, "Fig7B_IVW_vs_RAPS_with95CI_600dpi.tiff", width = 6.6, height = 5.3)

legend_lines <- c(
  "Figure 6. Steiger directionality analysis in prioritized heart-brain MR pairs.",
  "Figure 6A displays the number of prioritized pairs by MR direction (Heart to Brain vs Brain to Heart).",
  "Figure 6B plots summed SNP-level variance explained in exposure (x-axis) against outcome (y-axis), with color indicating MR direction and point shape indicating Steiger decision. The dashed diagonal indicates equality.",
  "",
  "Figure 7. Robustness of causal effect estimates using IVW and MR-RAPS.",
  "Figure 7A compares pair-specific causal effect estimates (beta) and 95% confidence intervals from IVW and MR-RAPS for prioritized pairs, stratified by MR direction.",
  "Figure 7B presents pair-level agreement between IVW and MR-RAPS estimates with 95% confidence intervals on both axes; the dashed diagonal indicates perfect agreement."
)
writeLines(legend_lines, con = file.path(out_dir, "Figure_Legends_Fig6_Fig7_EN.txt"))

cat(sprintf("out_dir=%s\n", out_dir))
cat(sprintf("steiger_n=%d\n", nrow(steiger)))
cat(sprintf("raps_n=%d\n", nrow(raps)))
