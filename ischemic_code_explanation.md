# 缺血性心脏病-认知/脑影像：技术说明（核心 R 脚本解读 + 结果解读）

本文聚焦两个脚本：
- 主流水线：[ischemic_heart_decease-cognition.R](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R)
- 独立纵向汇总脚本：[longitudinal_individual_delta_rate_summary.R](file:///Users/zbt/Documents/MI-COG2/ischemic/longitudinal_individual_delta_rate_summary.R)

两者的关系可以理解为：
- `ischemic_heart_decease-cognition.R` 负责“从队列 → 匹配 → 四模型纵向统计 → 可视化/脑图映射 → 纵向变化率汇总（多种口径）”的一条龙产物；
- `longitudinal_individual_delta_rate_summary.R` 是“拿四模型合并结果 + 某个比较的匹配队列”，额外计算“个体变化（绝对/百分比/平滑百分比）+ 年龄斜率”的独立脚本，适合单独复跑、做敏感性或补图。

## 0. SCI Methods（可直接用于论文的方法段落 + 公式口径）

### 0.1 研究设计与总体流程

本项目为 UK Biobank（UKB）两时点影像-认知纵向分析。以影像随访的两次扫描为配对时间点：Instance.2（基线，记为 $t_1$）与 Instance.3（随访，记为 $t_2$）。总体流程如下：
1) 三表合并构建分析队列；2) 在主比较（ischemic_vs_control）上执行倾向评分匹配（PSM），并复用匹配到的 `eid` 生成子比较队列；3) 对每个 IDP（认知或脑影像）执行配对的纵向回归，输出 Z 统计与效应量；4) 进行多重比较校正（FDR 与置换 FWE）；5) 对显著 IDP 生成变化率摘要与脑区映射可视化。

本段方法严格对应主脚本 [ischemic_heart_decease-cognition.R](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R) 与两个辅助脚本 [longitudinal_individual_delta_rate_summary.R](file:///Users/zbt/Documents/MI-COG2/ischemic/longitudinal_individual_delta_rate_summary.R)、[extract_idps_from_longitudinal_median_delta.R](file:///Users/zbt/Documents/MI-COG2/ischemic/extract_idps_from_longitudinal_median_delta.R) 的当前实现。

### 0.2 数据来源、时间定义与随访长度

三份输入数据表（`data1.csv`、`data2.csv`、`combined_brain_imaging_data.csv`）以 `eid` 统一并做全连接合并。两次影像日期分别记为 $D_{i,1}$（Instance.2）与 $D_{i,2}$（Instance.3），随访长度（年）定义为：

$$
FU_i = \frac{D_{i,2} - D_{i,1}}{365.25}
$$

纳入分析需满足：两次影像日期均存在且 $FU_i>0$。

### 0.3 暴露定义：纵向“新发”缺血性心脏病及其亚型

脚本以 $t_1$ 为基线，排除基线前既往发生的事件（例如 MI、慢性缺血、心绞痛、卒中、痴呆；判定条件为事件日期 $\le D_{i,1}$）。在此基础上，将 $t_1$ 与 $t_2$ 之间发生的事件定义为“新发（incident）”。以任一事件类型为例，其新发判定为：

$$
\mathrm{incident\_event}_i = \mathbb{1}(D^{event}_i > D_{i,1} \ \land\ D^{event}_i \le D_{i,2})
$$

主比较为 ischemic_vs_control（任何缺血性心脏病新发 vs 对照），并衍生子比较 mi_vs_control、chronic_vs_control、mi_vs_chronic（部分版本已移除该分组，但不影响主比较口径）。

### 0.4 倾向评分匹配（PSM）

在主比较中，使用逻辑回归估计倾向评分（propensity score）：

$$
e_i = P(T_i=1\mid \mathbf{X}_i) = \mathrm{logit}^{-1}\left(\alpha + \sum_{k}\beta_k X_{ik}\right)
$$

其中 $T_i$ 为病例指示变量（1=case，0=control），$\mathbf{X}_i$ 为匹配协变量（当前版本包含扫描年龄、性别、种族、教育、影像中心等；见 [ischemic_heart_decease-cognition.R:L2638-L2645](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L2638-L2645)）。匹配采用最近邻 2:1（control:case）匹配，caliper=0.2，估计目标为 ATT（Average Treatment effect on the Treated），实现封装于 [psm_match_standard](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L585-L622)。

匹配平衡通过标准化均差（SMD）评估，并输出匹配前后基线表与 Love plot（见 [compute_psm_balance](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L625-L737)）。对连续变量的 SMD 形式为：

$$
\mathrm{SMD}=\frac{\bar{x}_1-\bar{x}_0}{s_p},\quad
s_p=\sqrt{\frac{s_1^2+s_0^2}{2}}
$$

脚本在汇总层面以 $|\mathrm{SMD}|<0.1$ 作为“平衡良好”的经验阈值（见 [ischemic_heart_decease-cognition.R:L692-L695](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L692-L695)）。

### 0.5 IDP 纵向配对与主回归模型（Z 统计）

对每个 IDP，脚本将 Instance.2 作为基线测量 $Y_{i1}$、Instance.3 作为随访测量 $Y_{i2}$，并通过变量名规则自动配对（见 [ischemic_heart_decease-cognition.R:L3211-L3243](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L3211-L3243)）。主分析采用 ANCOVA 型纵向回归（控制基线值），核心形式为：

$$
Y_{i2}=\beta_0+\beta_1\cdot G_i+\beta_2\cdot Y_{i1}+\sum_{k}\gamma_k Z_{ik}+\varepsilon_i
$$

其中 $Z_{ik}$ 为混杂项（最小集合为随访年龄差与其二次项差，外加可选人口学变量）：

$$
\Delta A_i=A_{i2}-A_{i1},\quad \Delta A_i^{(2)}=A_{i2}^2-A_{i1}^2
$$

病例效应采用“年龄调制的病例指示变量”（age-modulated case effect），在每个 IDP 的可用样本内先去均值，再按随访年龄 $A_{i2}$ 进行缩放（见 [z_statistics_four_model_analysis](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L3731-L3929)）：

$$
C_i\in\{0,1\},\quad \tilde{C}_i = C_i-\bar{C}
$$

$$
w_i=10^{(0.0524\cdot A_{i2}-3.27)},\quad G_i=\tilde{C}_i\cdot w_i
$$

Z 统计取回归中 $G_i$ 的 $t$ 值（脚本直接使用 `t value` 作为 Z；见 [ischemic_heart_decease-cognition.R:L3793-L3799](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L3793-L3799)）：

$$
z_j = t(\hat{\beta}_{1j})=\frac{\hat{\beta}_{1j}}{\mathrm{SE}(\hat{\beta}_{1j})}
$$

效应量以 Cohen’s $d$ 近似表征（脚本口径为 $\hat{\beta}_1$ 除以 $Y_{2}$ 的总体标准差；见 [ischemic_heart_decease-cognition.R:L3800-L3802](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L3800-L3802)）：

$$
d_j=\frac{\hat{\beta}_{1j}}{\mathrm{SD}(Y_{\cdot 2})}
$$

### 0.6 非成像混杂的敏感性与“独立效应”判定

当提供候选非成像协变量集合（吸烟、饮酒、Townsend、BMI、糖尿病、血压及可选运动变量等）时，脚本拟合扩展模型并记录 $z$ 的变化，用于评估病例效应对非成像混杂的敏感性（见 [z_attenuation_by_non_imaging](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L4038-L4158)）。核心变化百分比定义为：

$$
\Delta z\%=\frac{z_{\mathrm{adj}}-z_{\mathrm{base}}}{|z_{\mathrm{base}}|}\times 100
$$

脚本在汇总层面以 $|\Delta z\%|<25\%$ 作为“相对稳健/独立”的经验阈值（用于 downstream 的独立效应筛选与可视化口径）。

### 0.7 多重比较：FDR（BH）与置换 FWE（Max-|t|）

FDR 采用 Benjamini–Hochberg（BH）方法对同一比较内的多重检验进行校正（见 [z_attenuation_by_non_imaging](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L4112-L4115) 的 `p.adjust(method="BH")` 用法）。

FWE 采用置换检验的 max-|t| 方法（见 [compute_fwe_pvalues_for_results](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L3932-L4036)）：对每次置换 $b=1,\dots,B$，在每个 IDP 的有效样本内对病例指示变量 $C_i$ 随机重排，重新拟合与主模型同形的回归，得到每个 IDP 的 $t_{b,j}$，并记录

$$
M_b=\max_j |t_{b,j}|
$$

对第 $j$ 个 IDP 的 FWE p 值使用加一修正：

$$
p^{\mathrm{FWE}}_j=\frac{1+\sum_{b=1}^{B}\mathbb{1}(M_b\ge |t_{\mathrm{obs},j}|)}{B+1}
$$

### 0.8 变化率与中位数差异（Median Delta Rate）

除回归 Z 统计外，脚本提供 IDP 的变化率摘要，核心来自 [compute_delta_rate_long](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L1836-L1999) 与 [summarize_delta_rate_median_vs_control_all_idps](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L6160-L6260)。

记绝对变化为：

$$
\Delta Y_i = Y_{i2}-Y_{i1}
$$

默认变化率方法为 “baseline-normalized”（可通过 `options(DELTA_RATE_METHOD=...)` 切换），其分母在存在 `group` 列时使用组内基线均值的绝对值，从而避免个体 $Y_{i1}\approx 0$ 导致的比例爆炸：

$$
d_g = \left|\mathrm{mean}(Y_{i1}\mid g)\right|,\quad
r_i = \frac{\Delta Y_i}{d_{g(i)}}
$$

脚本还实现了对称分母（symmetric）与对数比（log_ratio）两种备选方法：

$$
r_i^{\mathrm{symmetric}}=\frac{Y_{i2}-Y_{i1}}{\frac{|Y_{i1}|+|Y_{i2}|}{2}},\quad
r_i^{\mathrm{log\_ratio}}=\log(Y_{i2})-\log(Y_{i1})
$$

在 “Median Delta Rate” 汇总中，先将 $r_i$ 转为百分数 $100r_i$，再分别取病例组与对照组的中位数并做差：

$$
\Delta_{\mathrm{median}} = \mathrm{median}(100r_i\mid \mathrm{case})-\mathrm{median}(100r_i\mid \mathrm{control})
$$

该结果写入 `Output_Tables/Longitudinal/Median_DeltaRate/Longitudinal_Median_DeltaRate_AllIDPs_Comparison_<model>.csv`，并可用脚本 [extract_idps_from_longitudinal_median_delta.R](file:///Users/zbt/Documents/MI-COG2/ischemic/extract_idps_from_longitudinal_median_delta.R) 将表格中的 `median_diff_case_minus_control` 自动回填到用户指定的 Excel/CSV 模板中（以 `idp_root/variable_name` 作为键）。

### 0.9 个体层面变化与年龄斜率（独立纵向脚本）

独立脚本 [longitudinal_individual_delta_rate_summary.R](file:///Users/zbt/Documents/MI-COG2/ischemic/longitudinal_individual_delta_rate_summary.R) 在匹配队列内逐个体计算：

$$
\mathrm{abs\_change}_i = Y_{i2}-Y_{i1},\quad
\mathrm{pct\_change}_i = \frac{Y_{i2}-Y_{i1}}{\mathrm{denom}_i}\times 100
$$

并对“计数型变量”（近似非负整数）优先使用绝对变化，非计数变量使用百分比变化作为 `preferred_metric`（见 [longitudinal_individual_delta_rate_summary.R:L199-L209](file:///Users/zbt/Documents/MI-COG2/ischemic/longitudinal_individual_delta_rate_summary.R#L199-L209)）。年龄斜率通过在病例组与对照组内分别拟合线性模型 $change\sim age$ 获取，并计算斜率差（case−control）。

### 0.10 脑图映射与网络/血管分区汇总

对显著且（按阈值）独立的脑影像 IDP，脚本解析 `variable_name` 以提取结构类型、半球、脑区与测量类型，输出 [Brain_IDPs_Parsed_for_Mapping.csv](file:///Users/zbt/Documents/MI-COG2/ischemic/Brain_MRI_Maps/ischemic_vs_control/Brain_IDPs_Parsed_for_Mapping.csv) 作为脑图绘制的“单一事实来源”（见 [create_mri_brain_mapping](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L7722-L8062)）。皮层脑图基于 DK atlas，网络脑图对齐到 Yeo7（见 [align_to_yeo7_region](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L1646-L1657)），血管分区脑图将 DK 皮层区按字符串规则映射到 ACA/MCA/PCA 并在“区域×半球”层面汇总（见 [visualize_vascular_territories_brain_map](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L8739-L8841)）。

## 1. 主流水线 ischemic_heart_decease-cognition.R：从输入到输出

### 1.1 模块化开关与目录约定

脚本顶部有五个开关（可按需跳过步骤，但要求相应中间文件已存在）：[ischemic_heart_decease-cognition.R:L60-L66](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L60-L66)
- `RUN_DATA_FUSION`：读取/合并原始 CSV，构建最终队列（输出 Cohort RDS/CSV）。
- `RUN_PSM_MATCHING`：在总比较执行一次 PSM，并复用匹配 `eid` 生成子比较队列（输出 Matched RDS/CSV + 诊断图与表）。
- `RUN_IDP_ANALYSIS`：四模型纵向 Z 统计分析（输出 Four_Model 主表/分模型表）。
- `RUN_VISUALIZATION`：统计图、诊断图、箱线图、Z 分布等。
- `RUN_BRAIN_MAPPING`：基于显著/独立效应 IDP 输出脑区映射数据与脑图。

脚本统一把结果放入 `Output_Tables/`，通过 `get_table_dir()/table_out_path()` 保证目录存在：[ischemic_heart_decease-cognition.R:L12-L21](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L12-L21)。

### 1.2 Step 1：数据融合与队列构建（RUN_DATA_FUSION）

入口块： [ischemic_heart_decease-cognition.R:L83-L422](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L83-L422)

**输入**
- `data1.csv`、`data2.csv`、`combined_brain_imaging_data.csv`（用 vroom 读入）。

**核心处理逻辑**
- ID 统一：把不同表的 `Participant.ID/participant.id` 统一为 `eid`，[L97-L114](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L97-L114)。
- 三表 full join：先合并 `data1+data2`，再与 `data3` 合并，[L122-L126](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L122-L126)。
- 日期解析：`parse_ukb_date()` 支持 ymd/dmy/mdy 并清理“分号/逗号后续信息”，用于 UKB 混合格式日期列，[L154-L164](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L154-L164)。
- 合并列去重：`unify_joined_columns()` 将 `.x/.y` 等后缀列合并回基础列（优先 `.y`，再 coalesce），避免后续协变量/IDP 引用冲突，[L166-L194](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L166-L194)。

**纳入/排除标准（以第一次影像扫描 Instance.2 作为“基线”）**
- 需要两次影像日期均非空：`Instance.2` 与 `Instance.3`，[L198-L201](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L198-L201)。
- 随访时间 `followup_years > 0`：[L211-L217](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L211-L217)。
- 排除基线前已发生的：MI、慢性缺血、心绞痛、卒中、痴呆（全部以 `date_* <= date_img1` 判定），[L219-L240](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L219-L240)。

**分组定义（纵向“新发”）**
- `incident_any_ischemic = incident_MI | incident_chronic_ischemic`，[L232-L233](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L232-L233)。
- `Group` 用于人类可读输出：`Ischemic Heart Disease` vs `Control Group`，[L304-L306](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L304-L306)。
- `disease_subtype`：MI only / Chronic only / Both / Control，[L307-L312](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L307-L312)。

**协变量标准化（用于 PSM、基线表与后续敏感性分析）**
- `sex_factor / ethnicity_factor / education_factor` 等按字符串派生并 factor 化，[L241-L265](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L241-L265)。
- BMI、血压、Townsend 缺失用中位数填补（这里属于“就地填补”，不是 mice 的多重插补），[L314-L322](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L314-L322)。

**输出**
- `Output_Tables/Cohort/final_ischemic_cohort.csv`、`Output_Tables/Cohort/final_ischemic_cohort.rds`，[L334-L336](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L334-L336)。
- `Output_Tables/Cohort/Cohort_Flow_Report.csv`（最小纳入/排除简报），[L339-L405](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L339-L405)。

### 1.3 Step 2：PSM 匹配（RUN_PSM_MATCHING）

入口块： [ischemic_heart_decease-cognition.R:L2629-L3015](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L2629-L3015)

**关键设计：只做一次 PSM，然后复用匹配 eid**
- 仅在 `ischemic_vs_control` 上跑 PSM（2:1 最近邻、caliper=0.2），然后把匹配出来的 `eid` 作为“全局匹配集”，再在该集合里切出 `mi_vs_control`、`chronic_vs_control`、`mi_vs_chronic`，[L2649-L2779](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L2649-L2779)。
- 这样可以让子比较在“相同匹配策略与匹配池”下更可比，避免每个子比较单独匹配导致样本结构不一致。

**PSM 协变量（最新版本）**
- `age_at_instance2`（扫描1年龄）、`age_at_recruitment`（基线年龄别名）、`sex_factor`、`ethnicity_factor`、`education_factor`、`imaging_center_factor`，[L2638-L2645](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L2638-L2645)。
- PSM 前显式派生 `age_at_instance2/3`，保证协变量可用，[L2652-L2677](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L2652-L2677)。

**非成像协变量集合（用于后续四模型扩展/衰减分析）**
- 以 `*_baseline` 别名为主（避免 `_locf` 后缀差异），并额外按关键词抓取 Instance.0 的自然状况与血检列，[L2678-L2697](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L2678-L2697)。

**匹配诊断与输出**
- Love plot：`Figure_LovePlot_ischemic_vs_control.*`，[L2710-L2724](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L2710-L2724)。
- PS 分布：`Figure_PS_Distribution_ischemic_vs_control.*`，[L2721-L2724](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L2721-L2724)。
- 基线表（匹配前/后、P 值、整体平衡指标）输出到 `Output_Tables/Baseline/PSM/`，[L2727-L2732](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L2727-L2732)。
- 最重要的修复点：保存“匹配后的全量列”，而不是只保存 PSM 子集列，避免丢失认知/血检等后续分析所需变量，[L2740-L2744](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L2740-L2744)。

**匹配队列输出（后续所有分析的基础数据）**
- `Output_Tables/Cohort/Matched/Final_Matched_Cohort_<model>.csv/.rds`（4 个比较都有），主比较在 [L2736-L2745](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L2736-L2745)，子比较在后续循环中生成。

### 1.4 Step 3：四模型纵向 IDP 分析（RUN_IDP_ANALYSIS）

入口块： [ischemic_heart_decease-cognition.R:L3041-L5705](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L3041-L5705)

该模块的目标是：对每个 IDP（含认知与脑影像）在两个时间点（Instance.2 → Instance.3）的变化，进行“病例 vs 对照”的纵向比较，并用一套四模型设计（含敏感性变体）输出 Z 统计、效应量、多重校正与独立性判定。

#### 1.4.1 变量识别与两时间点配对

**识别规则**
- 认知变量/脑影像变量用关键词模式识别，并做动态补充；同时引入“纯认知”识别策略以防止漏检，[L3083-L3201](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L3083-L3201)。

**配对规则**
- 核心思想：用变量名中的 `Instance.2` 推出对应 `Instance.3`，并生成 `base_name` 作为 IDP 根名，[L3211-L3243](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L3211-L3243)。
- 产物：`variable_pairs` 列表（每项包含 `idp1_var/idp2_var/base_name/variable_type`），并输出配对计数，[L3245-L3251](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L3245-L3251)。

#### 1.4.2 年龄列与 Instance 列的“稳健化”

由于数据来源多表合并，常出现：
- `baseline_age/age_at_recruitment/followup_age` 命名不一致；
- 同一变量存在 `.x/.y` 或 `Instance_2/Instance 2` 等多种写法；
脚本在四模型分析前进行了集中规整：
- `compute_age_instances()` 派生 `age_at_instance2/age_at_instance3`，[L3257-L3288](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L3257-L3288)。
- `canonicalize_instance_columns()` 用 `coalesce` 把同一“规范 key”的多个列合成一个列，[L3292-L3311](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L3292-L3311)。
- `dedupe_instance_columns()` 删除冗余列，[L3315-L3327](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L3315-L3327)。
- `enforce_two_timepoint_presence()` 仅保留“至少一个 Instance.2 与一个 Instance.3 数值存在”的受试者，[L3329-L3340](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L3329-L3340)。

#### 1.4.3 四模型 Z 统计：核心回归与年龄加权

核心函数：`z_statistics_four_model_analysis()`，[ischemic_heart_decease-cognition.R:L3350-L3870](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L3350-L3870)（该范围包含主要建模、提取 t 值为 Z、以及一部分结果拼装）。

**基础回归式（每个 IDP 配对一条回归）**
- 目标：检验 “病例 vs 对照” 是否影响 `IDP2`（随访值），同时控制 `IDP1`（基线值）：
  - `IDP2 ~ Case_vs_Control + IDP1 + confounds`
- confounds 至少包含 `age_diff_linear/age_sq_diff_linear`，并可选加入 `sex_factor/ethnicity_factor`，[L3399-L3403](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L3399-L3403)。

**Case_vs_Control 的定义（年龄加权）**
- `case_control_binary`：从 `group`（0/1）或 `comparison_group` 反推病例标记，[L3381-L3392](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L3381-L3392)。
- 去均值：`case_control_binary - mean(case_control_binary)`，[L3396](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L3396)。
- 年龄权重：`10^(Age2*0.0524 - 3.27)`，[L3397-L3398](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L3397-L3398)。
- 最终：`Case_vs_Control = demeaned * aging_multiplier`，[L3398](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L3398)。

这意味着模型并不是单纯的 `case_control_binary`，而是把病例效应按随访年龄（`Age2`）做了幅度缩放。解释结果时应把它理解为“病例效应随年龄调制（age-modulated case effect）”。

**Z 统计的来源**
- 从 `lm()` 的系数表中取 `Case_vs_Control` 的 `t value` 作为 `z_statistic_without_non_imaging`，[L3411-L3417](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L3411-L3417)。

**效应量（Cohen’s d）**
- `cohens_d = beta / sd(IDP2)`（pooled_sd 用 `IDP2` 的方差开根），[L3418-L3421](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L3418-L3421)。

#### 1.4.4 多重比较：FDR 与 FWE

脚本同时输出两条口径：
- FDR：BH 校正（通常是 `p_value` → `p_fdr`）。
- FWE：置换 max-|t|（Max-t）得到 family-wise 校正 p。

FWE 的核心实现位于 `compute_fwe_pvalues_for_results()`（在脚本中定义并在各比较调用）。由于该函数较长，建议通过 IDE 直接跳转到其定义处并阅读其“置换流程、随机种子、n_perm”的配置。

#### 1.4.5 非成像混杂与衰减（attenuation）分析：独立性判定

当 `non_imaging_covariates` 给定时，脚本会对每个 IDP 做“逐一加入单个非成像变量”的敏感性：
- 观察 `Case_vs_Control` 的 Z 变化百分比：`z_change_percent`
- 判定“稳定独立”：`abs(z_change_percent) < 25`（即 |ΔZ%|<25%）

这套逻辑的用途是：把“显著的纵向关联”进一步分成
- 可能主要由非成像因素解释（加入某个混杂后 Z 明显衰减）；
- 相对稳健（加入候选混杂后 Z 仍相近）。

#### 1.4.6 纵向变化率摘要（在主脚本内部直接生成）

主脚本除了四模型 Z 统计，还会额外生成“变化率”类的结果（偏结果汇总/解释层面）：

1) **显著 IDP 的均值变化率摘要（按 Control/Case 分组）**
- `summarize_significant_idp_change_rates()` 会调用 `compute_idp_change_rate_summary_A(independent_only=TRUE)` 输出 `Output_Tables/Longitudinal/Summary/Significant_IDP_Change_Rate_Summary_<model>.csv`，[L5853-L5860](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L5853-L5860)。
- 该摘要取 `raw_change = IDP2 - IDP1`，并计算 `percent_change_mean = change_mean / |idp1_mean| * 100`，[L5792-L5816](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L5792-L5816)。

2) **全量 IDP 的均值变化率摘要**
- `summarize_all_idp_change_rates_mean()` 输出 `All_IDP_Change_Rate_Summary_<model>.csv`，[L5862-L5869](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L5862-L5869)。

3) **全量 IDP 的“中位数 Δ率对比”（推荐用于稳健描述）**
- `summarize_delta_rate_median_vs_control_all_idps()` 输出：
  - `Output_Tables/Longitudinal/Median_DeltaRate/Longitudinal_Median_DeltaRate_AllIDPs_Comparison_<model>.csv`，[L5968-L5971](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L5968-L5971)。
- 该函数会优先从四模型合并结果里抽取 `idp1_variable/idp2_variable` 的配对集合：
  - `pairs_df <- unique(sub_df[, c("idp1_variable","idp2_variable"), drop = FALSE])`，[L5889-L5892](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L5889-L5892)。
- 这段逻辑保证“纵向 Δ 率的计算口径”与四模型中实际分析的 IDP 配对一致，避免因配对集合不一致导致解释困难。

### 1.5 Step 5：脑图映射（RUN_BRAIN_MAPPING）

核心函数：`create_mri_brain_mapping()`，[ischemic_heart_decease-cognition.R:L7722-L8062](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L7722-L8062)

**输入**
- `brain_independent_idps`：通常来自四模型结果中“脑影像 + 独立效应”的子集（脚本在更靠后的主流程中构建）。

**输出**
- `Brain_MRI_Maps/Brain_IDPs_Parsed_for_Mapping.csv`：把 `variable_name` 解析成结构类型、半球、measure_type、network 等字段并统一成 `brain_region_final`，[L8053](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L8053)。
- 皮层脑图：`Brain_MRI_Maps/Cortical/Cortical_IDPs_Mapping.*`（ggseg + DK atlas），[L8067-L8136](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L8067-L8136)。
- 海马相关脑图/回退条形图：`Brain_MRI_Maps/Hippocampus/*`，[L8142-L8220](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L8142-L8220)。
- 网络映射（Yeo7 优先）：`Brain_MRI_Maps/Network/*`，[L8223-L8279](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L8223-L8279)。

## 2. 独立纵向脚本 longitudinal_individual_delta_rate_summary.R：个体变化与年龄斜率

脚本入口与参数解析： [longitudinal_individual_delta_rate_summary.R:L145-L167](file:///Users/zbt/Documents/MI-COG2/ischemic/longitudinal_individual_delta_rate_summary.R#L145-L167)

### 2.1 输入与输出

**必需输入**
- 四模型合并结果：默认 `Output_Tables/Four_Model/Combined_Four_Model_Z_Analysis_Results.csv`，[L147-L149](file:///Users/zbt/Documents/MI-COG2/ischemic/longitudinal_individual_delta_rate_summary.R#L147-L149)。
- 指定比较的匹配队列 RDS：默认 `Output_Tables/Cohort/Matched/Final_Matched_Cohort_<model>.rds`，[L158-L162](file:///Users/zbt/Documents/MI-COG2/ischemic/longitudinal_individual_delta_rate_summary.R#L158-L162)。

**主要输出**
- `Longitudinal_IndivPercentChange_Summary_<model>.csv` 写入 `out_dir`（默认是 `Output_Tables/Longitudinal/Median_DeltaRate`），在脚本后半段生成。
- 可选图：当 `--plots=true` 时，会对“计数型认知变量”生成按年龄分布的纵向变化图（输出到 `plot_dir`），默认 `Output_Tables/Longitudinal/Plots_Count_Cognitive/<model>/`，[L150-L153](file:///Users/zbt/Documents/MI-COG2/ischemic/longitudinal_individual_delta_rate_summary.R#L150-L153)。

### 2.2 与四模型结果的对齐方式（为什么要先做 meta）

脚本会先从 `Combined_Four_Model_Z_Analysis_Results.csv` 过滤出某个 `model_label` 的 `Brain_Imaging/Cognitive` 行，并要求存在 `idp1_variable/idp2_variable` 两列，[L168-L172](file:///Users/zbt/Documents/MI-COG2/ischemic/longitudinal_individual_delta_rate_summary.R#L168-L172)。

然后按 `idp_root`（从 `idp1_variable` 去掉 `Instance.2/3` 后缀得到）做去重：对同一根名，只保留 |Z| 最大（再按 p 最小）的那条记录作为“代表行”，[L173-L197](file:///Users/zbt/Documents/MI-COG2/ischemic/longitudinal_individual_delta_rate_summary.R#L173-L197)。

这样做的效果是：当同一 IDP 根名因为不同模型变体/字段重名出现多行时，后续的“变化率汇总”不会重复计算同一根名。

### 2.3 个体层面的变化指标与 preferred_metric

对每个 `idp1_variable/idp2_variable`，脚本读取匹配队列中两个时间点的值，构造三种变化：
- `abs_change = IDP2 - IDP1`（绝对变化）
- `pct_change = (IDP2 - IDP1) / denom * 100`（百分比变化；分母用基线值或其它分母策略，代码在后续段落完成）
- `smooth_pct_change`：当变量呈“计数型”时，用 `IDP1 + smooth_k` 作为平滑分母，避免基线=0 导致爆炸，[L155-L156](file:///Users/zbt/Documents/MI-COG2/ischemic/longitudinal_individual_delta_rate_summary.R#L155-L156)。

脚本通过 `is_count_like()` 判断变量是否“近似非负整数计数”（≥95% 值满足“非负整数”），[L208-L216](file:///Users/zbt/Documents/MI-COG2/ischemic/longitudinal_individual_delta_rate_summary.R#L208-L216)。
- 若 `is_count_like=TRUE`：`preferred_metric` 设为 `abs_change`（因为百分比在计数型变量上更易受 0/小分母影响）。
- 否则：`preferred_metric` 设为 `percent_change`（更利于跨 IDP 尺度对比）。

### 2.4 年龄斜率（case/control 分别拟合，再做差）

脚本会为 Control 与 Case 分别拟合 `change ~ age` 的线性回归，提取斜率：
- `control_age_slope_*` 与 `case_age_slope_*`
- 并输出 `age_slope_diff_*_case_minus_control` 用于对比病例组与对照组“变化随年龄的趋势差异”

年龄列用 `find_age_col()` 自动寻找（优先 `baseline_age`），[L116-L122](file:///Users/zbt/Documents/MI-COG2/ischemic/longitudinal_individual_delta_rate_summary.R#L116-L122)。

## 3. 结合结果文件：如何解读关键输出

### 3.1 四模型合并主结果

主结果表：
- `Output_Tables/Four_Model/Combined_Four_Model_Z_Analysis_Results.csv`（由主脚本写出，路径函数见 [ischemic_heart_decease-cognition.R:L39-L41](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L39-L41)）

解读建议优先关注这些字段（字段名以实际 CSV 为准）：
- `model_comparison`：四个比较之一（以及可能的 `_no_sex_ethnicity` 变体）。
- `variable_type`：`Brain_Imaging` / `Cognitive`。
- `variable_name`：IDP 根名（用于跨模块汇总与映射）。
- `idp1_variable / idp2_variable`：Instance.2 与 Instance.3 的具体列名（用于回溯原始数据）。
- `z_statistic_without_non_imaging`、`p_value_without_non_imaging`：基础模型（不加非成像混杂）的统计强度。
- `z_statistic_with_non_imaging`、`p_value_with_non_imaging`：扩展模型（加非成像混杂）的统计强度。
- `p_value_fdr`、`p_value_fwe`：多重比较校正后的显著性（FDR/FWE 两套口径）。
- `cohens_d`：效应量（基于 `IDP2` 的尺度标准化）。
- `independent_effect` / `independent_of_confounds` / `stable_independent`（如存在）：对“显著性 + 稳健性”的二次判定标签。

### 3.2 中位数 Δ率对比（全量 IDP）

文件示例：
- [Longitudinal_Median_DeltaRate_AllIDPs_Comparison_ischemic_vs_control.csv](file:///Users/zbt/Documents/MI-COG2/ischemic/Output_Tables/Longitudinal/Median_DeltaRate/Longitudinal_Median_DeltaRate_AllIDPs_Comparison_ischemic_vs_control.csv)

关键字段解释（见该 CSV 表头）：
- `control_median / case_median`：两组的 `delta_percent = delta_rate * 100` 的中位数。
- `median_diff_case_minus_control`：病例组中位数 − 对照组中位数，适合用作“组间差异”的稳健描述。
- `variable_type / idp_category`：用于分层汇总（如把脑影像按 Surface Area / Thickness / WM Tract FA 等分类）。

### 3.3 个体变化率 + 年龄斜率（独立脚本输出）

当你需要回答“变化随年龄是否更陡”时，`longitudinal_individual_delta_rate_summary.R` 输出的这些字段更直接：
- `preferred_age_slope_diff_case_minus_control`
- `preferred_case_age_slope`、`preferred_control_age_slope`
- `preferred_median_diff_case_minus_control`、`preferred_mean_diff_case_minus_control`

## 4. 推荐复现顺序（最小闭环）

- Step 1：运行主脚本到 PSM（确保生成匹配队列 RDS/CSV）
  - 产物：`Output_Tables/Cohort/Matched/Final_Matched_Cohort_<model>.rds`
- Step 2：运行主脚本到 Four_Model（生成合并结果）
  - 产物：`Output_Tables/Four_Model/Combined_Four_Model_Z_Analysis_Results.csv`
- Step 3：需要更细的“个体变化 + 年龄斜率”时，单跑独立纵向脚本
  - 使用 `--model=<model>` 指向某个比较，复用 Step 1/2 的产物

## 5. 常见问题定位（按现象）

- 现象：某些比较结果条目少/显著数明显偏低
  - 先看该比较匹配队列的有效样本量（`Final_Matched_Cohort_<model>`）
  - 再看该比较可配对的 `idp1_variable/idp2_variable` 是否在队列中真实存在（变量缺失会被跳过）
  - 最后看多重校正阈值（FWE 在置换数不大时更保守）

- 现象：纵向 Δ 率对不上四模型显著结果
  - 确认 Δ 率计算是否使用了与四模型一致的配对集合（主脚本在 [L5889-L5892](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R#L5889-L5892) 已做对齐；若另写脚本需保持同样逻辑）

## 6. 其他脚本：用途、输入输出与与主流水线的关系

本仓库里还有一些“单用途脚本”，它们不一定被主流水线直接调用，但经常用于补充分析、二次加工或可视化。下面按功能分组说明。

### 6.1 cohort.R：早期/独立的队列构建与 PSM 工具脚本

脚本：[cohort.R](file:///Users/zbt/Documents/MI-COG2/ischemic/cohort.R)

**定位**
- `cohort.R` 是一个“从原始表合并 → 构建 incident 队列 → PSM → 输出匹配诊断与基线表”的自包含脚本。
- 当前主流水线 [ischemic_heart_decease-cognition.R](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic_heart_decease-cognition.R) 已把同类逻辑内置并统一输出到 `Output_Tables/`；`cohort.R` 更像“可单独跑的最小队列/匹配模块”。

**关键函数**
- `compute_pscore_glm()`：二项逻辑回归计算 propensity score，[cohort.R:L392](file:///Users/zbt/Documents/MI-COG2/ischemic/cohort.R#L392)
- `nearest_match()`：最近邻匹配（默认 2:1、caliper=0.2、不可重复使用对照），[cohort.R:L398](file:///Users/zbt/Documents/MI-COG2/ischemic/cohort.R#L398)
- `run_psm_comparison()`：对指定比较执行 m 次“中位数插补 + 匹配”，并输出匹配示例与基线表/诊断图，[cohort.R:L665](file:///Users/zbt/Documents/MI-COG2/ischemic/cohort.R#L665)
- `report_cohort_flow()`：生成最小纳入/排除简报 CSV，[cohort.R:L737](file:///Users/zbt/Documents/MI-COG2/ischemic/cohort.R#L737)

**主要输入（脚本内部读取）**
- `data1.csv`、`data2.csv`、`combined_brain_imaging_data.csv`（与主脚本同源；在脚本前半段合并）。

**主要输出（写到当前工作目录 getwd() 下）**
- `final_ischemic_cohort.csv/.rds`：最终队列
- `Cohort_Flow_Report.csv`：纳入/排除简报
- `IDP_Missingness/IDP_missingness_by_*.csv`：按组/亚型的 IDP 缺失率
- `PSM_Output/<comparison>/`：基线表、PS 分布图、Love plot、PSM flowchart、匹配示例等（`comparison` 包含 `ischemic_vs_control` 及其子比较）

### 6.2 brain_age_independent_analysis.R：脑龄（BAG）预测与“独立效应”回归

脚本：[brain_age_independent_analysis.R](file:///Users/zbt/Documents/MI-COG2/ischemic/brain_age_independent_analysis.R)

**定位**
- 从匹配队列里抽取“可用的脑影像 IDP 列”，在 Control 组上训练 Random Forest 预测年龄；
- 在目标数据（通常仍是匹配队列）上输出 `predicted_brain_age` 和 `brain_age_gap (BAG)`；
- 再用线性回归检验 `BAG ~ ischemic_binary + 协变量`，输出缺血相关的调整后效应（作为“独立性”参考）。

**入口参数（命令行 args）**
- `train_path`、`test_path`：训练/测试数据（csv/tsv/rds 均可）；默认指向 `Output_Tables/Cohort/Matched/Final_Matched_Cohort_ischemic_vs_control.csv`，[brain_age_independent_analysis.R:L40-L45](file:///Users/zbt/Documents/MI-COG2/ischemic/brain_age_independent_analysis.R#L40-L45)
- `comparison`：用于派生 `group_label`（Control/MI/Chronic/Ischemic），默认 `ischemic_vs_control`
- `out_dir`：输出目录，默认 `Output_Tables/Brain_Age_Analysis_Independent`
- `n_max`：控制训练样本下采样（只对 Control 训练集生效）
- `ntree`：Random Forest 树数（默认 500）

**关键函数**
- `ensure_baseline_age()`：自动从多候选列推导 `baseline_age`，[brain_age_independent_analysis.R:L105](file:///Users/zbt/Documents/MI-COG2/ischemic/brain_age_independent_analysis.R#L105)
- `select_idps()`：基于列名规则过滤出“像脑影像 IDP 的数值列”，并限制最大列数（默认 128），[brain_age_independent_analysis.R:L265](file:///Users/zbt/Documents/MI-COG2/ischemic/brain_age_independent_analysis.R#L265)
- `fit_and_save()`：训练 RF、生成 BAG、按组汇总与回归检验并输出图表，[brain_age_independent_analysis.R:L290](file:///Users/zbt/Documents/MI-COG2/ischemic/brain_age_independent_analysis.R#L290)

**主要输出（在 out_dir）**
- `BAG_Predictions.csv`：包含 `baseline_age/predicted_brain_age/brain_age_gap`（以及 `group_label/eid` 若存在）
- `BAG_Summary.csv`、`BAG_Group_Stats.csv`
- `BAG_Independent_Effect_Adjusted.csv`：`ischemic_binary` 的调整后回归系数（若组信息可解析）
- `BAG_Predicted_vs_Actual.png/pdf`、`BAG_Distribution.png/pdf` 等图

### 6.3 extract_idps_from_longitudinal_median_delta.R：把 Δ率结果回填到 Excel/模板

脚本：[extract_idps_from_longitudinal_median_delta.R](file:///Users/zbt/Documents/MI-COG2/ischemic/extract_idps_from_longitudinal_median_delta.R)

**定位**
- 给一个 Excel/CSV 模板（包含 IDP 名称列 + 一个待填充的 `%` 列），
- 从 `Longitudinal_Median_DeltaRate_AllIDPs_Comparison_<model>.csv` 里匹配 `idp_root`，
- 把 `median_diff_case_minus_control` 回填到模板的 `%` 列，输出一个“已填充”的新文件。

**输入**
- 模板文件：默认 `IDPs.xlsx`，也可传 CSV/TSV
- Δ率结果目录：优先 `Output_Tables/Longitudinal/Median_DeltaRate/`，找不到则回退到当前目录，[extract_idps_from_longitudinal_median_delta.R:L94-L101](file:///Users/zbt/Documents/MI-COG2/ischemic/extract_idps_from_longitudinal_median_delta.R#L94-L101)

**关键逻辑**
- 自动识别模板里“IDP 列”和“% 列”（基于候选列名集合匹配），[extract_idps_from_longitudinal_median_delta.R:L83-L88](file:///Users/zbt/Documents/MI-COG2/ischemic/extract_idps_from_longitudinal_median_delta.R#L83-L88)
- 自动推断/展开 `model_comparison`：
  - 模板若带 model 列则按行指定；
  - 否则尝试从 sheet 名/前 20 行文本推断；
  - 再否则对所有模型文件展开并逐个回填，[extract_idps_from_longitudinal_median_delta.R:L105-L120](file:///Users/zbt/Documents/MI-COG2/ischemic/extract_idps_from_longitudinal_median_delta.R#L105-L120)
- 每个模型的回填核心：`fill_one_model()`，[extract_idps_from_longitudinal_median_delta.R:L124-L136](file:///Users/zbt/Documents/MI-COG2/ischemic/extract_idps_from_longitudinal_median_delta.R#L124-L136)

**输出**
- 默认输出：把输入文件名改为 `_filled.xlsx`（或输出 CSV）
- 依赖：读 Excel 需要 `readxl`，写 Excel 需要 `writexl`（脚本会在缺包时报错提示安装）

### 6.4 T1WI_mapping.py：变量名增强解析（脑区/半球/测量类型）

脚本：[T1WI_mapping.py](file:///Users/zbt/Documents/MI-COG2/ischemic/T1WI_mapping.py)

**定位**
- 读取 `All_Independent_Effect_IDPs_for_Visualization.csv`（通常来自四模型独立效应汇总），
- 用规则+词典把 `variable_name` 解析成结构化字段（`brain_region/lobe/hemisphere/measure_type/...`），
- 输出“增强解析版 CSV”与“未识别清单/统计摘要”，便于后续脑图映射或人工校对。

**主要入口**
- `process_brain_data_with_enhanced_parser(file_path, save_outputs=True)`，[T1WI_mapping.py:L1267](file:///Users/zbt/Documents/MI-COG2/ischemic/T1WI_mapping.py#L1267)
- 脚本直接运行时默认寻找：
  - `Output_Tables/Four_Model/All_Independent_Effect_IDPs_for_Visualization.csv`
  - 或当前目录下同名文件，[T1WI_mapping.py:L1399-L1412](file:///Users/zbt/Documents/MI-COG2/ischemic/T1WI_mapping.py#L1399-L1412)

**输出（默认同目录/同前缀）**
- `<input>_enhanced_parsed.csv`
- `<input>_unknown_regions.csv`
- `<input>_summary_stats.txt`

### 6.5 effect_brain_mapping.py：基于 nilearn 的 NIfTI/平面图/交互 HTML 脑图

脚本：[effect_brain_mapping.py](file:///Users/zbt/Documents/MI-COG2/ischemic/effect_brain_mapping.py)

**定位**
- 建立一个“IDP → 脑区（atlas label / MNI 坐标）”的解析与映射流水线；
- 将显著/独立 IDP 生成体素空间的统计图（NIfTI），并输出静态平面图、玻璃脑、表面投影以及交互 HTML。

**运行方式**
- 脚本运行时默认读取：
  - `Output_Tables/Four_Model/All_Independent_Effect_IDPs_for_Visualization.csv`
  - 或当前目录下同名文件，[effect_brain_mapping.py:L1211-L1217](file:///Users/zbt/Documents/MI-COG2/ischemic/effect_brain_mapping.py#L1211-L1217)
- 主入口 `main()` 会创建 `IDPBrainMappingPipeline()` 并执行 `process_micog_data()`，[effect_brain_mapping.py:L1208-L1221](file:///Users/zbt/Documents/MI-COG2/ischemic/effect_brain_mapping.py#L1208-L1221)

**典型输出（脚本打印的清单）**
- `parsed_idp_results.csv`
- `idp_brain_map.nii.gz`
- `anatomical_planes.png`
- `glass_brain_3d.png`
- `surface_projections.png`
- `statistical_report.png`
- `interactive_brain_map.html`，[effect_brain_mapping.py:L1195-L1202](file:///Users/zbt/Documents/MI-COG2/ischemic/effect_brain_mapping.py#L1195-L1202)

### 6.6 wm_jhu_overlay.py：白质 JHU 叠加图 + UofM 概率图渲染

脚本：[wm_jhu_overlay.py](file:///Users/zbt/Documents/MI-COG2/ischemic/wm_jhu_overlay.py)

**定位**
- 读取 `Brain_IDPs_Parsed_for_Mapping.csv`（主脚本 R 侧脑图模块会生成），
- 将白质 tract 名称规范化到 JHU label 体系，
- 生成每个 tract 的 voxel overlay（强度通常编码为 max |z|），并输出 montage/coverage 统计；
- 可选渲染 Harvard-Oxford overlays 与 UofM 白质概率图（含 glass brain）。

**默认输入/输出位置**
- 默认输入：`Brain_MRI_Maps/Brain_IDPs_Parsed_for_Mapping.csv`，[wm_jhu_overlay.py:L41-L43](file:///Users/zbt/Documents/MI-COG2/ischemic/wm_jhu_overlay.py#L41-L43)
- 默认输出目录：
  - `Brain_MRI_Maps/WhiteMatterVoxel/`（JHU voxel overlays）
  - 以及 `Brain_MRI_Maps/Cortical/`、`Brain_MRI_Maps/Subcortical/`、`Brain_MRI_Maps/Network/` 的一些补充输出，[wm_jhu_overlay.py:L43-L46](file:///Users/zbt/Documents/MI-COG2/ischemic/wm_jhu_overlay.py#L43-L46)

**命令行入口**
- `python wm_jhu_overlay.py <Brain_IDPs_Parsed_for_Mapping.csv> [--atlas-dir ... --uofm-dir ... --prob-threshold ...]`，[wm_jhu_overlay.py:L1602-L1647](file:///Users/zbt/Documents/MI-COG2/ischemic/wm_jhu_overlay.py#L1602-L1647)

### 6.7 UofM atlas 附带的 3D Slicer 工具脚本（非主流水线）

脚本：[CaptureRotationVideo.py](file:///Users/zbt/Documents/MI-COG2/ischemic/ischemic/atlases/uofm_jhu_atlas_v1/Video_Capture_Scripts/CaptureRotationVideo.py)

**定位**
- 该脚本是给 3D Slicer 的 Python Console 使用的“截图序列采集”工具，用于旋转 3D 视图并逐帧导出 PNG。
- 它包含硬编码的输出目录与依赖（`slicer/vtk` 环境），不属于本仓库主分析流水线的一部分。
