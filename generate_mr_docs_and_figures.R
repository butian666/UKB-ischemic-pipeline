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

in_dir <- kv$in_dir %||% file.path(getwd(), "Output_Tables", "MR_Bidirectional_Bulk")
out_dir <- kv$out_dir %||% file.path(in_dir, "Manuscript_and_Figures")
catalog_csv <- kv$catalog_csv %||% file.path(getwd(), "mr_brain_catalog_bulk_from_independent.csv")
ihd_source <- kv$ihd_source %||% "opengwas:finn-b-I9_IHD"
top_n <- suppressWarnings(as.integer(kv$top_n %||% "30"))
if (is.na(top_n) || top_n <= 0) top_n <- 30L
detail_top_n <- suppressWarnings(as.integer(kv$detail_top_n %||% "6"))
if (is.na(detail_top_n) || detail_top_n <= 0) detail_top_n <- 6L

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "Figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "Figures", "MR_Detail_Plots"), recursive = TRUE, showWarnings = FALSE)

read_csv_flex <- function(path) {
  if (requireNamespace("data.table", quietly = TRUE)) return(data.table::fread(path, data.table = FALSE))
  utils::read.csv(path, check.names = FALSE)
}

sanitize_filename <- function(x) {
  y <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
  y <- gsub("_+", "_", y)
  y <- gsub("^_|_$", "", y)
  if (!nzchar(y)) y <- "NA"
  y
}

parse_source <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  if (grepl("^opengwas:", x)) return(sub("^opengwas:", "", x))
  x
}

main_csv <- file.path(in_dir, "MR_Manuscript_Main_Table.csv")
count_csv <- file.path(in_dir, "MR_Manuscript_Summary_Counts.csv")
status_csv <- file.path(in_dir, "MR_Analysis_Summary.csv")
if (!file.exists(main_csv)) stop("缺少文件: ", main_csv)
if (!file.exists(count_csv)) stop("缺少文件: ", count_csv)
if (!file.exists(status_csv)) stop("缺少文件: ", status_csv)

tbl <- read_csv_flex(main_csv)
cnt <- read_csv_flex(count_csv)
st <- read_csv_flex(status_csv)

tbl$direction <- as.character(tbl$direction)
tbl$ivw_significance <- as.character(tbl$ivw_significance)
tbl$evidence_tier <- ifelse(tbl$ivw_significance == "significant_fdr", "Tier1_FDR",
                     ifelse(tbl$ivw_significance == "nominal" & tbl$sensitivity_flag == "pass", "Tier2_Nominal_SensitivityPass",
                     ifelse(tbl$ivw_significance == "nominal", "Tier3_Nominal",
                     ifelse(tbl$ivw_significance == "borderline", "Tier4_Borderline", "Tier5_Null"))))
tbl$neglog10p <- -log10(pmax(tbl$ivw_p, 1e-300))
tbl$direction_plot <- ifelse(tbl$direction == "heart_to_brain", "Heart to Brain", "Brain to Heart")
tbl$pair_label <- paste(tbl$exposure_trait, "=>", tbl$outcome_trait)

library(ggplot2)

p1 <- ggplot(tbl, aes(x = ivw_b, y = neglog10p, color = ivw_significance, shape = direction_plot)) +
  geom_point(alpha = 0.8, size = 2.2) +
  scale_color_manual(values = c(significant_fdr = "#b2182b", nominal = "#ef8a62", borderline = "#fddbc7", ns = "#4d4d4d", failed = "#969696")) +
  labs(x = "IVW beta", y = "-log10(P)", color = "Significance", shape = "Direction", title = "Bidirectional MR volcano plot") +
  theme_bw(base_size = 12)
ggsave(file.path(out_dir, "Figures", "Fig1_MR_Volcano_AllPairs.png"), p1, width = 9, height = 6, dpi = 320)

ord <- order(tbl$ivw_p, decreasing = FALSE)
top_tbl <- tbl[ord[seq_len(min(top_n, nrow(tbl)))], , drop = FALSE]
top_tbl$pair_short <- paste0(seq_len(nrow(top_tbl)), ". ", substr(top_tbl$outcome_trait, 1, 40))
p2 <- ggplot(top_tbl, aes(x = reorder(pair_short, ivw_p), y = neglog10p, fill = direction_plot)) +
  geom_col() +
  coord_flip() +
  labs(x = NULL, y = "-log10(P)", fill = "Direction", title = paste0("Top ", nrow(top_tbl), " associations by IVW P")) +
  theme_bw(base_size = 11)
ggsave(file.path(out_dir, "Figures", "Fig2_MR_TopAssociations_Bar.png"), p2, width = 10, height = 8, dpi = 320)

heat_tbl <- top_tbl
heat_tbl$exposure_short <- ifelse(nchar(heat_tbl$exposure_trait) > 32, paste0(substr(heat_tbl$exposure_trait, 1, 32), "..."), heat_tbl$exposure_trait)
heat_tbl$outcome_short <- ifelse(nchar(heat_tbl$outcome_trait) > 32, paste0(substr(heat_tbl$outcome_trait, 1, 32), "..."), heat_tbl$outcome_trait)
p3 <- ggplot(heat_tbl, aes(x = outcome_short, y = exposure_short, fill = neglog10p)) +
  geom_tile() +
  scale_fill_gradient(low = "#f7fbff", high = "#08306b") +
  labs(x = "Outcome", y = "Exposure", fill = "-log10(P)", title = "MR signal heatmap (top associations)") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave(file.path(out_dir, "Figures", "Fig3_MR_Heatmap_TopPairs.png"), p3, width = 12, height = 8, dpi = 320)

cons_tbl <- data.frame(
  metric = c("2-of-3 direction consistent", "Sensitivity passed", "Nominal", "Borderline", "FDR significant"),
  value = c(sum(tbl$direction_consistent_2of3 %in% TRUE, na.rm = TRUE), sum(tbl$sensitivity_flag == "pass", na.rm = TRUE), sum(tbl$ivw_significance == "nominal", na.rm = TRUE), sum(tbl$ivw_significance == "borderline", na.rm = TRUE), sum(tbl$ivw_significance == "significant_fdr", na.rm = TRUE))
)
p4 <- ggplot(cons_tbl, aes(x = metric, y = value, fill = metric)) +
  geom_col(show.legend = FALSE) +
  labs(x = NULL, y = "Count", title = "Consistency and sensitivity summary") +
  theme_bw(base_size = 12)
ggsave(file.path(out_dir, "Figures", "Fig4_MR_Consistency_Sensitivity.png"), p4, width = 7.5, height = 5, dpi = 320)

tier_tbl <- aggregate(pair_label ~ evidence_tier, data = tbl, FUN = length)
names(tier_tbl)[2] <- "count"
p5 <- ggplot(tier_tbl, aes(x = evidence_tier, y = count, fill = evidence_tier)) +
  geom_col(show.legend = FALSE) +
  labs(x = "Evidence tier", y = "Count", title = "Evidence tiers for MR results") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(out_dir, "Figures", "Fig5_MR_Evidence_Tier.png"), p5, width = 8.5, height = 5.5, dpi = 320)

catalog <- read_csv_flex(catalog_csv)
catalog$source_id <- vapply(catalog$source, parse_source, character(1))
trait_to_id <- stats::setNames(catalog$source_id, catalog$trait)
ihd_id <- parse_source(ihd_source)

if (requireNamespace("TwoSampleMR", quietly = TRUE)) {
  library(TwoSampleMR)
  pos_tbl <- tbl[tbl$ivw_significance %in% c("significant_fdr", "nominal", "borderline"), , drop = FALSE]
  pos_tbl <- pos_tbl[order(pos_tbl$ivw_p, decreasing = FALSE), , drop = FALSE]
  if (nrow(pos_tbl) > 0) {
    pos_tbl <- pos_tbl[seq_len(min(detail_top_n, nrow(pos_tbl))), , drop = FALSE]
    detail_status <- list()
    for (i in seq_len(nrow(pos_tbl))) {
      r <- pos_tbl[i, , drop = FALSE]
      pair_id <- paste0(sprintf("%02d", i), "_", sanitize_filename(r$direction), "_", sanitize_filename(substr(r$exposure_trait, 1, 40)), "_to_", sanitize_filename(substr(r$outcome_trait, 1, 40)))
      exp_id <- NA_character_
      out_id <- NA_character_
      if (r$direction == "heart_to_brain") {
        exp_id <- ihd_id
        out_id <- trait_to_id[[r$outcome_trait]]
      } else {
        exp_id <- trait_to_id[[r$exposure_trait]]
        out_id <- ihd_id
      }
      if (is.null(exp_id) || is.na(exp_id) || !nzchar(exp_id) || is.null(out_id) || is.na(out_id) || !nzchar(out_id)) {
        detail_status[[length(detail_status) + 1]] <- data.frame(pair = pair_id, status = "skip_missing_id", stringsAsFactors = FALSE)
        next
      }
      ok <- TRUE
      msg <- "ok"
      tryCatch({
        exp_dat <- extract_instruments(exp_id, p1 = 5e-8, clump = TRUE, r2 = 0.001, kb = 10000)
        if (is.null(exp_dat) || nrow(exp_dat) == 0) stop("no_instruments")
        out_dat <- extract_outcome_data(snps = exp_dat$SNP, outcomes = out_id)
        if (is.null(out_dat) || nrow(out_dat) == 0) stop("no_outcome_overlap")
        harm <- harmonise_data(exp_dat, out_dat)
        harm <- harm[harm$mr_keep, , drop = FALSE]
        if (nrow(harm) < 1) stop("no_harmonised_snps")
        mr_res <- mr(harm)
        if (nrow(mr_res) < 1) stop("no_mr_result")

        scat <- mr_scatter_plot(mr_res, harm)
        if (length(scat) > 0) ggsave(file.path(out_dir, "Figures", "MR_Detail_Plots", paste0(pair_id, "_scatter.png")), scat[[1]], width = 6.5, height = 5.2, dpi = 320)

        sing <- tryCatch(mr_singlesnp(harm), error = function(e) NULL)
        if (!is.null(sing) && nrow(sing) > 0) {
          forest <- tryCatch(mr_forest_plot(sing), error = function(e) NULL)
          if (!is.null(forest) && length(forest) > 0) ggsave(file.path(out_dir, "Figures", "MR_Detail_Plots", paste0(pair_id, "_forest.png")), forest[[1]], width = 7.2, height = 8.2, dpi = 320)
          funnel <- tryCatch(mr_funnel_plot(sing), error = function(e) NULL)
          if (!is.null(funnel) && length(funnel) > 0) ggsave(file.path(out_dir, "Figures", "MR_Detail_Plots", paste0(pair_id, "_funnel.png")), funnel[[1]], width = 6.2, height = 5.2, dpi = 320)
        }

        loo <- tryCatch(mr_leaveoneout(harm), error = function(e) NULL)
        if (!is.null(loo) && nrow(loo) > 0) {
          loo_p <- tryCatch(mr_leaveoneout_plot(loo), error = function(e) NULL)
          if (!is.null(loo_p) && length(loo_p) > 0) ggsave(file.path(out_dir, "Figures", "MR_Detail_Plots", paste0(pair_id, "_leaveoneout.png")), loo_p[[1]], width = 7.2, height = 7.8, dpi = 320)
        }
      }, error = function(e) {
        ok <<- FALSE
        msg <<- as.character(e$message)
      })
      detail_status[[length(detail_status) + 1]] <- data.frame(pair = pair_id, direction = r$direction, exposure = r$exposure_trait, outcome = r$outcome_trait, status = ifelse(ok, "ok", "failed"), message = msg, stringsAsFactors = FALSE)
    }
    detail_status_df <- do.call(rbind, detail_status)
    utils::write.csv(detail_status_df, file.path(out_dir, "Figures", "MR_Detail_Plots", "detail_plot_status.csv"), row.names = FALSE)
  }
}

top_hits <- tbl[order(tbl$ivw_p), c("direction","exposure_trait","outcome_trait","nsnp_ivw","ivw_b","ivw_se","ivw_p","ivw_fdr","ivw_significance","direction_consistent_2of3","sensitivity_flag")]
top_hits <- top_hits[seq_len(min(20, nrow(top_hits))), , drop = FALSE]
utils::write.csv(top_hits, file.path(out_dir, "Top20_MR_Associations.csv"), row.names = FALSE)

n_total <- st$total_pairs[1]
n_ok <- st$successful_pairs[1]
n_rows <- st$mr_rows[1]
n_sig <- cnt$significant_fdr[1]
n_nom <- cnt$nominal[1]
n_bor <- cnt$borderline[1]
n_cons <- cnt$direction_consistent_2of3[1]
n_sens <- cnt$sensitivity_pass[1]

cn_lines <- c(
  "# 方法学写法（中文）",
  "",
  "## 研究设计",
  paste0("本研究采用双向两样本孟德尔随机化（Two-sample MR）评估缺血性心肌病（IHD）与脑影像表型（IDPs）的潜在因果关联。IHD暴露来源为 FinnGen（finn-b-I9_IHD），脑结构候选表型来自 OpenGWAS 中可用的 ubm-* 脑影像GWAS条目，并由自动映射流程从独立效应IDP清单中构建。"),
  "",
  "## 工具变量构建与数据处理",
  "暴露端SNP筛选阈值设为 P<5e-8，LD去相关参数设置为 r2<0.001、窗口10,000 kb，EUR参考群体；随后进行结局提取、等位基因协调、回文位点与不一致等位基因处理。",
  "",
  "## 因果估计与稳健性分析",
  "主分析方法为IVW；同时报告MR-Egger与Weighted median。对于具备条件的结果对，进一步评估异质性（Cochran’s Q）与水平多效性（Egger截距）。",
  "",
  "## 多重检验与证据分层",
  "对IVW P值进行FDR校正；结果按FDR显著、名义显著、边缘显著与不显著分层，并结合方向一致性（IVW/Egger/WM至少2/3方向一致）与敏感性通过情况进行证据分级。"
)

en_lines <- c(
  "# Methods (English)",
  "",
  "## Study design",
  "We performed bidirectional two-sample Mendelian randomization (MR) to evaluate potential causal links between ischemic heart disease (IHD) and brain imaging-derived phenotypes (IDPs). IHD exposure was obtained from FinnGen (finn-b-I9_IHD). Brain phenotypes were drawn from OpenGWAS ubm-* neuroimaging GWAS resources and mapped from the independent-effect IDP list via an automated matching workflow.",
  "",
  "## Instrument selection and harmonization",
  "Instrumental SNPs were selected at P<5e-8 with LD clumping (r2<0.001, 10,000 kb, EUR reference). Outcome data were extracted and harmonized with allele alignment, and ambiguous palindromic/incompatible variants were filtered according to standard MR procedures.",
  "",
  "## Causal estimation and robustness",
  "IVW was used as the primary estimator, with MR-Egger and weighted median as complementary methods. Heterogeneity (Cochran's Q) and directional pleiotropy (Egger intercept) were evaluated when data conditions allowed.",
  "",
  "## Multiple testing and evidence grading",
  "IVW P values were FDR-adjusted. Evidence was graded by FDR significance, nominal significance, borderline significance, and null findings, further combined with direction consistency (>=2/3 consistent signs across IVW/Egger/WM) and sensitivity-pass flags."
)

res_cn <- c(
  "# 结果段落 + 补充材料说明（中文）",
  "",
  "## 主结果",
  paste0("本次批量双向MR共测试 ", n_total, " 个暴露-结局组合，其中 ", n_ok, " 对成功完成分析，形成 ", n_rows, " 条方法学结果记录。"),
  paste0("FDR显著关联数为 ", n_sig, "；名义显著（P<0.05）为 ", n_nom, "；边缘显著（0.05<=P<0.10）为 ", n_bor, "。"),
  paste0("在证据稳健性方面，", n_cons, " 对结果满足至少2/3方法方向一致，", n_sens, " 对满足敏感性通过（异质性与Egger截距均未提示明显偏倚）。"),
  "",
  "## 图表与补充材料",
  "主文图建议纳入：全量火山图（Fig1）、Top关联条形图（Fig2）、信号热图（Fig3）、方向一致性与敏感性汇总（Fig4）、证据分层图（Fig5）。",
  "补充材料建议纳入：Top20关联清单、每个阳性/边缘阳性结果的scatter/forest/funnel/leave-one-out图组、详细运行状态表与映射细节表。",
  "",
  "## 顶刊风格扩展图建议",
  "建议补充共定位（colocalization）证据图、Steiger方向检验汇总图、MR-RAPS/CAUSE敏感性一致性图、跨队列重复验证森林图与负对照结果图，以增强因果链可信度。"
)

res_en <- c(
  "# Results + Supplement Description (English)",
  "",
  "## Main findings",
  paste0("A total of ", n_total, " exposure-outcome pairs were tested in bidirectional MR; ", n_ok, " pairs were successfully analyzed, yielding ", n_rows, " method-level result rows."),
  paste0("The number of FDR-significant associations was ", n_sig, "; nominal associations (P<0.05) were ", n_nom, "; borderline associations (0.05<=P<0.10) were ", n_bor, "."),
  paste0("Regarding robustness, ", n_cons, " pairs showed at least 2-of-3 directional consistency across IVW/Egger/weighted median, and ", n_sens, " pairs passed the sensitivity flag."),
  "",
  "## Figures and supplementary outputs",
  "Recommended main figures include: global volcano plot (Fig1), top-association bar chart (Fig2), signal heatmap (Fig3), consistency/sensitivity summary (Fig4), and evidence-tier distribution (Fig5).",
  "Recommended supplements include: Top20 association table, per-pair scatter/forest/funnel/leave-one-out panels for positive/borderline pairs, full run-status tables, and mapping detail tables.",
  "",
  "## Additional top-journal style figures",
  "To further strengthen causal interpretation, add colocalization plots, Steiger direction summary, MR-RAPS/CAUSE sensitivity concordance panels, cross-cohort replication forest plots, and negative-control analyses."
)

writeLines(c(cn_lines, "", "---", "", en_lines), con = file.path(out_dir, "Manuscript_Methods_CN_EN.md"))
writeLines(c(res_cn, "", "---", "", res_en), con = file.path(out_dir, "Manuscript_Results_and_Supp_CN_EN.md"))

cat(sprintf("out_dir=%s\n", out_dir))
cat(sprintf("pairs_total=%d\n", n_total))
cat(sprintf("pairs_success=%d\n", n_ok))
cat(sprintf("mr_rows=%d\n", n_rows))
cat(sprintf("nominal=%d\n", n_nom))
cat(sprintf("borderline=%d\n", n_bor))
