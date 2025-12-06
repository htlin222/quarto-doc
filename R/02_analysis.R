# ==============================================================================
# 02_analysis.R
# Breast Cancer Exercise Intervention Trial - Main Analysis
# ==============================================================================

# Load required packages
library(tidyverse)
library(survival)
library(survminer)
library(ggsurvfit)
library(gtsummary)
library(broom)
library(nlme)
library(lme4)
library(gt)

# Set working directory context
set.seed(20241206)

# ==============================================================================
# Load Data
# ==============================================================================

breast_data <- readRDS("data/breast_cancer_exercise_trial.rds")
longitudinal_data <- readRDS("data/breast_cancer_longitudinal.rds")

message(sprintf("Loaded %d patients", nrow(breast_data)))
message(sprintf("Loaded %d longitudinal observations", nrow(longitudinal_data)))

# ==============================================================================
# 1. BASELINE CHARACTERISTICS TABLE (Table 1)
# ==============================================================================

message("\n=== Creating Table 1: Baseline Characteristics ===")

# Prepare data for Table 1
table1_data <- breast_data |>
  dplyr::select(
    treatment_group, age, sex, bmi, menopause_status,
    tumor_stage, tumor_grade, hormone_receptor_status, her2_status,
    molecular_subtype, lymph_nodes_positive, surgery_type,
    chemo_regimen, radiation_therapy, endocrine_therapy,
    targeted_therapy, ecog_ps, comorbidities
  ) |>
  dplyr::mutate(
    treatment_group = factor(
      treatment_group,
      levels = c("control", "exercise"),
      labels = c("Control", "Exercise")
    ),
    tumor_stage = factor(tumor_stage, levels = c("I", "II", "III")),
    tumor_grade = factor(tumor_grade),
    ecog_ps = factor(ecog_ps)
  )

# Create Table 1 with gtsummary
table1 <- table1_data |>
  tbl_summary(
    by = treatment_group,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    label = list(
      age ~ "Age (years)",
      sex ~ "Sex",
      bmi ~ "BMI (kg/m²)",
      menopause_status ~ "Menopausal status",
      tumor_stage ~ "Tumor stage",
      tumor_grade ~ "Tumor grade",
      hormone_receptor_status ~ "Hormone receptor status",
      her2_status ~ "HER2 status",
      molecular_subtype ~ "Molecular subtype",
      lymph_nodes_positive ~ "Positive lymph nodes",
      surgery_type ~ "Surgery type",
      chemo_regimen ~ "Chemotherapy regimen",
      radiation_therapy ~ "Radiation therapy",
      endocrine_therapy ~ "Endocrine therapy",
      targeted_therapy ~ "Targeted therapy",
      ecog_ps ~ "ECOG performance status",
      comorbidities ~ "Comorbidities"
    ),
    missing = "no"
  ) |>
  add_p() |>
  add_overall() |>
  modify_header(label ~ "**Characteristic**") |>
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Treatment Group**") |>
  bold_labels()

# Save Table 1
table1 |>
  as_gt() |>
  gt::gtsave("output/tables/table1_baseline_characteristics.html")

# Export to text format
table1_text <- table1 |> as_tibble()
write.csv(table1_text, "output/tables/table1_baseline_characteristics.csv",
          row.names = FALSE)

message("Table 1 saved to output/tables/")

# ==============================================================================
# 2. SURVIVAL ANALYSIS - Disease-Free Survival
# ==============================================================================

message("\n=== Performing Survival Analysis ===")

# Prepare survival data
surv_data <- breast_data |>
  dplyr::mutate(
    treatment_group = factor(
      treatment_group,
      levels = c("control", "exercise"),
      labels = c("Control", "Exercise")
    )
  )

# --- 2a. Kaplan-Meier Analysis for DFS ---

# Fit KM model
km_dfs <- survfit(
  Surv(dfs_time, dfs_event) ~ treatment_group,
  data = surv_data
)

# Create DFS survival plot with ggsurvfit
p_dfs <- km_dfs |>
  ggsurvfit(linewidth = 1) +
  add_confidence_interval() +
  add_risktable(
    risktable_stats = c("n.risk", "cum.event"),
    risktable_group = "risktable_stats",
    stats_label = list(
      n.risk = "At Risk",
      cum.event = "Events"
    )
  ) +
  add_pvalue(location = "annotation", caption = "Log-rank {p.value}") +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8")) +
  scale_x_continuous(breaks = seq(0, 48, 12), limits = c(0, 48)) +
  scale_y_continuous(
    limits = c(0.5, 1),
    labels = scales::percent
  ) +
  labs(
    title = "Disease-Free Survival by Treatment Group",
    x = "Time (months)",
    y = "Disease-Free Survival Probability"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )

ggsave(
  "output/figures/fig1_dfs_kaplan_meier.png",
  p_dfs, width = 10, height = 8, dpi = 300
)
ggsave(
  "output/figures/fig1_dfs_kaplan_meier.pdf",
  p_dfs, width = 10, height = 8
)

message("DFS Kaplan-Meier plot saved")

# --- 2b. Cox Proportional Hazards Model for DFS ---

# Unadjusted Cox model
cox_dfs_unadj <- coxph(
  Surv(dfs_time, dfs_event) ~ treatment_group,
  data = surv_data
)

# Adjusted Cox model (stratified by randomization factors)
cox_dfs_adj <- coxph(
  Surv(dfs_time, dfs_event) ~ treatment_group +
    strata(stratification_stage, stratification_bmi, stratification_menopause),
  data = surv_data
)

# Fully adjusted Cox model
cox_dfs_full <- coxph(
  Surv(dfs_time, dfs_event) ~ treatment_group + age + bmi +
    tumor_stage + tumor_grade + hormone_receptor_status + her2_status +
    lymph_nodes_positive + surgery_type,
  data = surv_data
)

# Create Cox model results table
cox_table_dfs <- tbl_regression(
  cox_dfs_full,
  exponentiate = TRUE,
  label = list(
    treatment_group ~ "Treatment group",
    age ~ "Age (per year)",
    bmi ~ "BMI (per kg/m²)",
    tumor_stage ~ "Tumor stage",
    tumor_grade ~ "Tumor grade",
    hormone_receptor_status ~ "Hormone receptor status",
    her2_status ~ "HER2 status",
    lymph_nodes_positive ~ "Positive lymph nodes",
    surgery_type ~ "Surgery type"
  )
) |>
  bold_p() |>
  bold_labels()

cox_table_dfs |>
  as_gt() |>
  gt::gtsave("output/tables/table2_cox_dfs.html")

# Extract Cox model results with broom
cox_dfs_tidy <- broom::tidy(cox_dfs_full, exponentiate = TRUE, conf.int = TRUE)
write.csv(cox_dfs_tidy, "output/reports/cox_dfs_result.csv", row.names = FALSE)

# Generate text report for DFS analysis
sink("output/reports/survival_analysis_dfs_result.txt")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("DISEASE-FREE SURVIVAL ANALYSIS RESULTS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

cat("1. KAPLAN-MEIER ANALYSIS\n")
cat("-" |> rep(40) |> paste(collapse = ""), "\n")
print(summary(km_dfs, times = c(12, 24, 36, 48)))

cat("\n2. LOG-RANK TEST\n")
cat("-" |> rep(40) |> paste(collapse = ""), "\n")
logrank_test <- survdiff(Surv(dfs_time, dfs_event) ~ treatment_group, data = surv_data)
print(logrank_test)

cat("\n3. COX PROPORTIONAL HAZARDS MODEL (Unadjusted)\n")
cat("-" |> rep(40) |> paste(collapse = ""), "\n")
print(summary(cox_dfs_unadj))

cat("\n4. COX PROPORTIONAL HAZARDS MODEL (Stratified)\n")
cat("-" |> rep(40) |> paste(collapse = ""), "\n")
print(summary(cox_dfs_adj))

cat("\n5. COX PROPORTIONAL HAZARDS MODEL (Fully Adjusted)\n")
cat("-" |> rep(40) |> paste(collapse = ""), "\n")
print(summary(cox_dfs_full))

cat("\n6. HAZARD RATIO SUMMARY\n")
cat("-" |> rep(40) |> paste(collapse = ""), "\n")
hr_unadj <- exp(coef(cox_dfs_unadj))["treatment_groupExercise"]
hr_adj <- exp(coef(cox_dfs_adj))["treatment_groupExercise"]
ci_unadj <- exp(confint(cox_dfs_unadj))["treatment_groupExercise", ]
ci_adj <- exp(confint(cox_dfs_adj))["treatment_groupExercise", ]

cat(sprintf("Unadjusted HR: %.3f (95%% CI: %.3f - %.3f)\n",
            hr_unadj, ci_unadj[1], ci_unadj[2]))
cat(sprintf("Stratified HR: %.3f (95%% CI: %.3f - %.3f)\n",
            hr_adj, ci_adj[1], ci_adj[2]))

sink()

message("DFS analysis report saved to output/reports/survival_analysis_dfs_result.txt")

# ==============================================================================
# 3. SURVIVAL ANALYSIS - Overall Survival
# ==============================================================================

message("\n=== Performing Overall Survival Analysis ===")

# Fit KM model for OS
km_os <- survfit(
  Surv(os_time, os_event) ~ treatment_group,
  data = surv_data
)

# Create OS survival plot
p_os <- km_os |>
  ggsurvfit(linewidth = 1) +
  add_confidence_interval() +
  add_risktable(
    risktable_stats = c("n.risk", "cum.event"),
    risktable_group = "risktable_stats"
  ) +
  add_pvalue(location = "annotation", caption = "Log-rank {p.value}") +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8")) +
  scale_x_continuous(breaks = seq(0, 48, 12), limits = c(0, 48)) +
  scale_y_continuous(
    limits = c(0.7, 1),
    labels = scales::percent
  ) +
  labs(
    title = "Overall Survival by Treatment Group",
    x = "Time (months)",
    y = "Overall Survival Probability"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )

ggsave(
  "output/figures/fig2_os_kaplan_meier.png",
  p_os, width = 10, height = 8, dpi = 300
)
ggsave(
  "output/figures/fig2_os_kaplan_meier.pdf",
  p_os, width = 10, height = 8
)

# Cox model for OS
cox_os <- coxph(
  Surv(os_time, os_event) ~ treatment_group +
    strata(stratification_stage, stratification_bmi, stratification_menopause),
  data = surv_data
)

cox_os_tidy <- broom::tidy(cox_os, exponentiate = TRUE, conf.int = TRUE)
write.csv(cox_os_tidy, "output/reports/cox_os_result.csv", row.names = FALSE)

# Generate OS report
sink("output/reports/survival_analysis_os_result.txt")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("OVERALL SURVIVAL ANALYSIS RESULTS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

cat("1. KAPLAN-MEIER ANALYSIS\n")
cat("-" |> rep(40) |> paste(collapse = ""), "\n")
print(summary(km_os, times = c(12, 24, 36, 48)))

cat("\n2. LOG-RANK TEST\n")
cat("-" |> rep(40) |> paste(collapse = ""), "\n")
logrank_os <- survdiff(Surv(os_time, os_event) ~ treatment_group, data = surv_data)
print(logrank_os)

cat("\n3. COX PROPORTIONAL HAZARDS MODEL (Stratified)\n")
cat("-" |> rep(40) |> paste(collapse = ""), "\n")
print(summary(cox_os))
sink()

message("OS analysis report saved")

# ==============================================================================
# 4. ADHERENCE ANALYSIS
# ==============================================================================

message("\n=== Performing Adherence Analysis ===")

# Filter exercise group
exercise_adherence <- breast_data |>
  dplyr::filter(treatment_group == "exercise") |>
  dplyr::select(
    patient_id, phase1_behavior_sessions, phase1_mandatory_exercise,
    phase1_optional_exercise, phase2_behavior_sessions,
    phase3_behavior_sessions, overall_adherence
  )

# Create adherence summary table
adherence_summary <- exercise_adherence |>
  dplyr::summarise(
    `Phase 1 Behavior Sessions (mean)` = mean(phase1_behavior_sessions, na.rm = TRUE),
    `Phase 1 Behavior Sessions (SD)` = sd(phase1_behavior_sessions, na.rm = TRUE),
    `Phase 1 Mandatory Exercise (mean)` = mean(phase1_mandatory_exercise, na.rm = TRUE),
    `Phase 1 Mandatory Exercise (SD)` = sd(phase1_mandatory_exercise, na.rm = TRUE),
    `Phase 2 Behavior Sessions (mean)` = mean(phase2_behavior_sessions, na.rm = TRUE),
    `Phase 2 Behavior Sessions (SD)` = sd(phase2_behavior_sessions, na.rm = TRUE),
    `Phase 3 Behavior Sessions (mean)` = mean(phase3_behavior_sessions, na.rm = TRUE),
    `Phase 3 Behavior Sessions (SD)` = sd(phase3_behavior_sessions, na.rm = TRUE),
    `Overall Adherence (mean)` = mean(overall_adherence, na.rm = TRUE),
    `Overall Adherence (SD)` = sd(overall_adherence, na.rm = TRUE)
  )

# Save adherence report
sink("output/reports/adherence_analysis_result.txt")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("ADHERENCE ANALYSIS RESULTS (Exercise Group Only)\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

cat("ADHERENCE SUMMARY\n")
cat("-" |> rep(40) |> paste(collapse = ""), "\n")
cat(sprintf("N = %d patients in exercise group\n\n", nrow(exercise_adherence)))

cat("Phase 1 Behavior Sessions (out of 12):\n")
cat(sprintf("  Mean: %.1f (SD: %.1f)\n",
            adherence_summary$`Phase 1 Behavior Sessions (mean)`,
            adherence_summary$`Phase 1 Behavior Sessions (SD)`))
cat(sprintf("  Adherence: %.1f%%\n\n",
            adherence_summary$`Phase 1 Behavior Sessions (mean)` / 12 * 100))

cat("Phase 1 Mandatory Exercise (out of 12):\n")
cat(sprintf("  Mean: %.1f (SD: %.1f)\n",
            adherence_summary$`Phase 1 Mandatory Exercise (mean)`,
            adherence_summary$`Phase 1 Mandatory Exercise (SD)`))
cat(sprintf("  Adherence: %.1f%%\n\n",
            adherence_summary$`Phase 1 Mandatory Exercise (mean)` / 12 * 100))

cat("Phase 2 Behavior Sessions (out of 12):\n")
cat(sprintf("  Mean: %.1f (SD: %.1f)\n",
            adherence_summary$`Phase 2 Behavior Sessions (mean)`,
            adherence_summary$`Phase 2 Behavior Sessions (SD)`))
cat(sprintf("  Adherence: %.1f%%\n\n",
            adherence_summary$`Phase 2 Behavior Sessions (mean)` / 12 * 100))

cat("Phase 3 Behavior Sessions (out of 24):\n")
cat(sprintf("  Mean: %.1f (SD: %.1f)\n",
            adherence_summary$`Phase 3 Behavior Sessions (mean)`,
            adherence_summary$`Phase 3 Behavior Sessions (SD)`))
cat(sprintf("  Adherence: %.1f%%\n\n",
            adherence_summary$`Phase 3 Behavior Sessions (mean)` / 24 * 100))

cat("Overall Adherence:\n")
cat(sprintf("  Mean: %.1f%% (SD: %.1f%%)\n",
            adherence_summary$`Overall Adherence (mean)`,
            adherence_summary$`Overall Adherence (SD)`))
sink()

message("Adherence analysis report saved")

# ==============================================================================
# 5. ADVERSE EVENTS ANALYSIS
# ==============================================================================

message("\n=== Performing Adverse Events Analysis ===")

# Prepare AE data
ae_data <- breast_data |>
  dplyr::mutate(
    treatment_group = factor(
      treatment_group,
      levels = c("control", "exercise"),
      labels = c("Control", "Exercise")
    ),
    any_ae = ae_any > 0,
    musculoskeletal_ae = ae_musculoskeletal > 0,
    cardiovascular_ae = ae_cardiovascular > 0,
    lymphedema_ae = ae_lymphedema > 0,
    grade3plus_ae = ae_any >= 3
  )

# Create AE table
ae_table <- ae_data |>
  dplyr::select(
    treatment_group, any_ae, musculoskeletal_ae, cardiovascular_ae,
    lymphedema_ae, grade3plus_ae
  ) |>
  tbl_summary(
    by = treatment_group,
    statistic = all_categorical() ~ "{n} ({p}%)",
    label = list(
      any_ae ~ "Any adverse event",
      musculoskeletal_ae ~ "Musculoskeletal adverse event",
      cardiovascular_ae ~ "Cardiovascular adverse event",
      lymphedema_ae ~ "Lymphedema",
      grade3plus_ae ~ "Grade 3+ adverse event"
    )
  ) |>
  add_p() |>
  bold_p() |>
  bold_labels()

ae_table |>
  as_gt() |>
  gt::gtsave("output/tables/table3_adverse_events.html")

# Save AE report
sink("output/reports/adverse_events_analysis_result.txt")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("ADVERSE EVENTS ANALYSIS RESULTS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

ae_summary <- ae_data |>
  dplyr::group_by(treatment_group) |>
  dplyr::summarise(
    n = dplyr::n(),
    any_ae_n = sum(any_ae),
    any_ae_pct = mean(any_ae) * 100,
    msk_ae_n = sum(musculoskeletal_ae),
    msk_ae_pct = mean(musculoskeletal_ae) * 100,
    cv_ae_n = sum(cardiovascular_ae),
    cv_ae_pct = mean(cardiovascular_ae) * 100,
    lymph_ae_n = sum(lymphedema_ae),
    lymph_ae_pct = mean(lymphedema_ae) * 100,
    grade3_n = sum(grade3plus_ae),
    grade3_pct = mean(grade3plus_ae) * 100,
    .groups = "drop"
  )

cat("ADVERSE EVENT SUMMARY BY TREATMENT GROUP\n")
cat("-" |> rep(50) |> paste(collapse = ""), "\n\n")

for (grp in unique(ae_summary$treatment_group)) {
  row <- ae_summary |> dplyr::filter(treatment_group == grp)
  cat(sprintf("%s Group (N = %d)\n", grp, row$n))
  cat(sprintf("  Any AE: %d (%.1f%%)\n", row$any_ae_n, row$any_ae_pct))
  cat(sprintf("  Musculoskeletal: %d (%.1f%%)\n", row$msk_ae_n, row$msk_ae_pct))
  cat(sprintf("  Cardiovascular: %d (%.1f%%)\n", row$cv_ae_n, row$cv_ae_pct))
  cat(sprintf("  Lymphedema: %d (%.1f%%)\n", row$lymph_ae_n, row$lymph_ae_pct))
  cat(sprintf("  Grade 3+: %d (%.1f%%)\n\n", row$grade3_n, row$grade3_pct))
}

# Chi-square tests
cat("STATISTICAL COMPARISONS\n")
cat("-" |> rep(50) |> paste(collapse = ""), "\n\n")

test_any <- chisq.test(table(ae_data$treatment_group, ae_data$any_ae))
test_msk <- chisq.test(table(ae_data$treatment_group, ae_data$musculoskeletal_ae))
test_lymph <- chisq.test(table(ae_data$treatment_group, ae_data$lymphedema_ae))

cat(sprintf("Any AE: Chi-square = %.2f, p = %.4f\n",
            test_any$statistic, test_any$p.value))
cat(sprintf("Musculoskeletal AE: Chi-square = %.2f, p = %.4f\n",
            test_msk$statistic, test_msk$p.value))
cat(sprintf("Lymphedema: Chi-square = %.2f, p = %.4f\n",
            test_lymph$statistic, test_lymph$p.value))
sink()

message("Adverse events analysis report saved")

message("\n=== Analysis Complete ===")
message("All tables saved to output/tables/")
message("All figures saved to output/figures/")
message("All text reports saved to output/reports/")
