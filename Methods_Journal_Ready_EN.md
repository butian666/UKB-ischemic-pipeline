# Methods (Journal-Ready; aligned with current analysis code)

## Study Design and Data Source

This study uses longitudinal data from the UK Biobank (UKB) with paired brain imaging assessments at two time points. Imaging time points are defined as Instance.2 (baseline; \(t_1\)) and Instance.3 (follow-up; \(t_2\)). All analyses are performed at the participant level and link phenotypes across source tables via the unified participant identifier (\(\mathrm{eid}\)).

## Imaging Time Points and Follow-up Time

For participant \(i\), imaging dates are denoted \(D_{i,1}\) (Instance.2) and \(D_{i,2}\) (Instance.3). Follow-up time in years is computed as:

\[
FU_i=\frac{D_{i,2}-D_{i,1}}{365.25}.
\]

Participants are required to have non-missing \(D_{i,1}\) and \(D_{i,2}\) with \(FU_i>0\).

## Exposure Definition: Incident Ischemic Heart Disease and Subtypes

All exposure definitions are anchored to the baseline imaging date \(D_{i,1}\). For each relevant cardiovascular or neurological event type, prevalent cases are excluded if the corresponding event date satisfies \(D^{event}_i \le D_{i,1}\). Incident events are defined as events occurring after baseline and on or before follow-up:

\[
\mathrm{incident\_event}_i=\mathbb{1}(D^{event}_i>D_{i,1}\ \land\ D^{event}_i\le D_{i,2}).
\]

The primary comparison is incident ischemic heart disease versus controls (ischemic\_vs\_control). Secondary comparisons include incident myocardial infarction versus controls (mi\_vs\_control) and incident chronic ischemic heart disease versus controls (chronic\_vs\_control). Control participants are defined as those without incident events in the corresponding window and meeting the baseline exclusion criteria.

## Propensity Score Matching (PSM)

To mitigate confounding, propensity score matching is performed in the primary comparison. A propensity score is estimated using logistic regression:

\[
e_i=P(T_i=1\mid\mathbf{X}_i)=\mathrm{logit}^{-1}\left(\alpha+\sum_k\beta_k X_{ik}\right),
\]

where \(T_i\) indicates case status (1=case, 0=control) and \(\mathbf{X}_i\) denotes the prespecified matching covariates (e.g., age at scan, sex, ethnicity, education, imaging center; covariate availability is checked within the analysis dataset prior to modeling).

Matching is implemented via nearest-neighbor matching with a 2:1 control-to-case ratio, caliper 0.2, without replacement, targeting the average treatment effect on the treated (ATT).

Covariate balance is assessed using standardized mean differences (SMD). For a continuous covariate \(x\), SMD is computed as:

\[
\mathrm{SMD}=\frac{\bar{x}_1-\bar{x}_0}{s_p},\quad s_p=\sqrt{\frac{s_1^2+s_0^2}{2}},
\]

where subscripts 1 and 0 refer to cases and controls, respectively. Balance is summarized using \(|\mathrm{SMD}|<0.1\) as a conventional threshold indicating acceptable balance.

## Longitudinal Imaging-Derived Phenotypes (IDPs) and Pairing Across Time

For each imaging-derived phenotype (IDP), the baseline measurement is \(Y_{i1}\) (Instance.2) and the follow-up measurement is \(Y_{i2}\) (Instance.3). IDPs are paired across instances using standardized naming rules to ensure that the same phenotype (including array position, where applicable) is compared across \(t_1\) and \(t_2\).

## Primary Longitudinal Model and Z Statistics

The primary analysis uses an ANCOVA-type longitudinal regression that models the follow-up value while adjusting for the baseline value:

\[
Y_{i2}=\beta_0+\beta_1\cdot G_i+\beta_2\cdot Y_{i1}+\sum_k\gamma_k Z_{ik}+\varepsilon_i,
\]

where \(Z_{ik}\) includes follow-up time–related nuisance covariates derived from the two imaging ages:

\[
\Delta A_i=A_{i2}-A_{i1},\quad \Delta A_i^{(2)}=A_{i2}^2-A_{i1}^2.
\]

The case effect is represented by an age-modulated case indicator. Let \(C_i\in\{0,1\}\) denote case status and \(\tilde{C}_i=C_i-\bar{C}\) be the within-sample mean-centered case indicator. A participant-specific modulation factor is defined based on follow-up age \(A_{i2}\):

\[
w_i=10^{(0.0524\cdot A_{i2}-3.27)},\quad G_i=\tilde{C}_i\cdot w_i.
\]

The primary test statistic for each IDP \(j\) is the \(t\)-statistic associated with \(\hat{\beta}_{1j}\), which is used as the Z statistic:

\[
z_j=t(\hat{\beta}_{1j})=\frac{\hat{\beta}_{1j}}{\mathrm{SE}(\hat{\beta}_{1j})}.
\]

Effect size is summarized using a standardized coefficient:

\[
d_j=\frac{\hat{\beta}_{1j}}{\mathrm{SD}(Y_{\cdot 2})},
\]

where \(\mathrm{SD}(Y_{\cdot 2})\) is the standard deviation of the follow-up IDP measurement within the analyzed sample for that IDP.

## Sensitivity to Non-imaging Covariates and “Independent Effect” Criterion

To evaluate robustness to additional non-imaging covariates (e.g., smoking, alcohol, deprivation index, BMI, diabetes, blood pressure, and optional activity-related variables), an expanded model is fitted when these covariates are available. Sensitivity is quantified by the relative change in Z statistic:

\[
\Delta z\%=\frac{z_{\mathrm{adj}}-z_{\mathrm{base}}}{|z_{\mathrm{base}}|}\times 100.
\]

For downstream reporting and visualization, an IDP may be designated as having a relatively robust/independent effect when \(|\Delta z\%|<25\%\) under the expanded adjustment set.

## Multiple Testing Correction

False discovery rate (FDR) is controlled within each comparison using the Benjamini–Hochberg procedure applied to the IDP-wise p-values.

Family-wise error (FWE) control is additionally implemented using a permutation-based max-\(|t|\) approach. For each permutation \(b=1,\dots,B\), the case indicator \(C_i\) is permuted within the valid sample for each IDP, the same longitudinal model is refitted, and the maximum absolute \(t\)-statistic across IDPs is recorded:

\[
M_b=\max_j |t_{b,j}|.
\]

The FWE-adjusted p-value for IDP \(j\) is computed using a plus-one correction:

\[
p^{\mathrm{FWE}}_j=\frac{1+\sum_{b=1}^{B}\mathbb{1}(M_b\ge |t_{\mathrm{obs},j}|)}{B+1}.
\]

## Derived Change Metrics and Brain Mapping Outputs

In addition to regression-based statistics, the analysis produces summary measures of longitudinal change and, for imaging phenotypes passing prespecified criteria, structured outputs for brain mapping. Brain mapping is performed by parsing each selected IDP into atlas-aligned region labels and hemisphere assignment, and then visualizing region-level summaries for cortical, network (Yeo7-aligned), and vascular-territory groupings.

