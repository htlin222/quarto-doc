# ==============================================================================
# 04_subgroup_analysis.R
# Subgroup Analysis with Forest Plots for Breast Cancer Exercise Trial
# ==============================================================================

# Load required packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(survival)
library(broom)
library(purrr)
library(forcats)

# ------------------------------------------------------------------------------
# Load Data
# ------------------------------------------------------------------------------

analysis_data <- readRDS("data/breast_cancer_exercise_trial.rds")
longitudinal_data <- readRDS("data/breast_cancer_longitudinal.rds")

# Get baseline MVPA from longitudinal data
baseline_mvpa_df <- longitudinal_data |>
  dplyr::filter(time_months == 0) |>
  dplyr::select(patient_id, baseline_mvpa)

analysis_data <- analysis_data |>
  dplyr::left_join(baseline_mvpa_df, by = "patient_id")

# Create subgroup variables
analysis_data <- analysis_data |>
  dplyr::mutate(
    # Age groups
    age_group = dplyr::case_when(
      age < 50 ~ "<50 years",
      age >= 50 & age < 65 ~ "50-64 years",
      age >= 65 ~ ">=65 years"
    ),
    age_group = factor(age_group, levels = c("<50 years", "50-64 years", ">=65 years")),
    # BMI groups
    bmi_group = dplyr::if_else(bmi <= 25, "BMI <=25", "BMI >25"),
    bmi_group = factor(bmi_group, levels = c("BMI <=25", "BMI >25")),
    # Baseline MVPA groups (median split)
    mvpa_group = dplyr::if_else(
      baseline_mvpa <= stats::median(baseline_mvpa, na.rm = TRUE),
      "Low baseline MVPA",
      "High baseline MVPA"
    ),
    mvpa_group = factor(mvpa_group, levels = c("Low baseline MVPA", "High baseline MVPA")),
    # Chemotherapy regimen groups
    chemo_group = dplyr::case_when(
      chemo_regimen %in% c("AC-T", "TAC") ~ "Anthracycline-based",
      chemo_regimen == "TC" ~ "TC regimen",
      TRUE ~ "Other"
    ),
    chemo_group = factor(chemo_group, levels = c("Anthracycline-based", "TC regimen", "Other"))
  )

# ------------------------------------------------------------------------------
# Function to Calculate Subgroup HRs
# ------------------------------------------------------------------------------

calculate_subgroup_hr <- function(data, subgroup_var, outcome = "dfs") {
  if (outcome == "dfs") {
    time_var <- "dfs_time"
    event_var <- "dfs_event"
  } else {
    time_var <- "os_time"
    event_var <- "os_event"
  }

  subgroup_levels <- unique(data[[subgroup_var]])
  subgroup_levels <- subgroup_levels[!is.na(subgroup_levels)]

  results <- purrr::map_dfr(subgroup_levels, function(level) {
    subset_data <- data[data[[subgroup_var]] == level, ]

    n_events_exercise <- sum(
      subset_data[[event_var]][subset_data$treatment_group == "exercise"], na.rm = TRUE
    )
    n_events_control <- sum(
      subset_data[[event_var]][subset_data$treatment_group == "control"], na.rm = TRUE
    )
    n_exercise <- sum(subset_data$treatment_group == "exercise", na.rm = TRUE)
    n_control <- sum(subset_data$treatment_group == "control", na.rm = TRUE)

    if (n_events_exercise < 3 || n_events_control < 3) {
      return(data.frame(
        subgroup = subgroup_var, level = as.character(level),
        n_total = nrow(subset_data), n_exercise = n_exercise, n_control = n_control,
        n_events = n_events_exercise + n_events_control,
        hr = NA_real_, hr_lower = NA_real_, hr_upper = NA_real_, p_value = NA_real_
      ))
    }

    formula_str <- paste0("Surv(", time_var, ", ", event_var, ") ~ treatment_group")
    cox_fit <- tryCatch(
      survival::coxph(stats::as.formula(formula_str), data = subset_data),
      error = function(e) NULL
    )

    if (is.null(cox_fit)) {
      return(data.frame(
        subgroup = subgroup_var, level = as.character(level),
        n_total = nrow(subset_data), n_exercise = n_exercise, n_control = n_control,
        n_events = n_events_exercise + n_events_control,
        hr = NA_real_, hr_lower = NA_real_, hr_upper = NA_real_, p_value = NA_real_
      ))
    }

    cox_summary <- summary(cox_fit)
    data.frame(
      subgroup = subgroup_var, level = as.character(level),
      n_total = nrow(subset_data), n_exercise = n_exercise, n_control = n_control,
      n_events = n_events_exercise + n_events_control,
      hr = cox_summary$conf.int[1, 1],
      hr_lower = cox_summary$conf.int[1, 3],
      hr_upper = cox_summary$conf.int[1, 4],
      p_value = cox_summary$coefficients[1, 5]
    )
  })

  results
}

# ------------------------------------------------------------------------------
# Function to Test for Interaction
# ------------------------------------------------------------------------------

test_interaction <- function(data, subgroup_var, outcome = "dfs") {
  if (outcome == "dfs") {
    time_var <- "dfs_time"
    event_var <- "dfs_event"
  } else {
    time_var <- "os_time"
    event_var <- "os_event"
  }

  formula_int <- paste0("Surv(", time_var, ", ", event_var, ") ~ treatment_group * ", subgroup_var)
  formula_main <- paste0("Surv(", time_var, ", ", event_var, ") ~ treatment_group + ", subgroup_var)

  cox_int <- tryCatch(
    survival::coxph(stats::as.formula(formula_int), data = data),
    error = function(e) NULL
  )
  cox_main <- tryCatch(
    survival::coxph(stats::as.formula(formula_main), data = data),
    error = function(e) NULL
  )

  if (is.null(cox_int) || is.null(cox_main)) return(NA_real_)

  lrt <- stats::anova(cox_main, cox_int)
  lrt$`Pr(>|Chi|)`[2]
}

# ------------------------------------------------------------------------------
# Calculate Subgroup Results for DFS
# ------------------------------------------------------------------------------

subgroups <- c(
  "age_group", "menopause_status", "tumor_stage", "hormone_receptor_status",
  "her2_status", "bmi_group", "mvpa_group", "chemo_group"
)

subgroup_labels <- c(
  "age_group" = "Age",
  "menopause_status" = "Menopause Status",
  "tumor_stage" = "Tumor Stage",
  "hormone_receptor_status" = "Hormone Receptor Status",
  "her2_status" = "HER2 Status",
  "bmi_group" = "BMI",
  "mvpa_group" = "Baseline Physical Activity",
  "chemo_group" = "Chemotherapy Regimen"
)

# Calculate HRs for each subgroup
dfs_subgroup_results <- purrr::map_dfr(subgroups, function(sg) {
  calculate_subgroup_hr(analysis_data, sg, outcome = "dfs")
})

# Calculate interaction p-values
interaction_pvalues <- purrr::map_dbl(subgroups, function(sg) {
  test_interaction(analysis_data, sg, outcome = "dfs")
})
names(interaction_pvalues) <- subgroups

# Add interaction p-values and labels
dfs_subgroup_results <- dfs_subgroup_results |>
  dplyr::mutate(
    interaction_p = interaction_pvalues[subgroup],
    subgroup_label = subgroup_labels[subgroup]
  )

# ------------------------------------------------------------------------------
# Overall HR for Reference
# ------------------------------------------------------------------------------

overall_cox <- survival::coxph(Surv(dfs_time, dfs_event) ~ treatment_group, data = analysis_data)
overall_summary <- summary(overall_cox)

overall_result <- data.frame(
  subgroup = "Overall", level = "All Patients",
  n_total = nrow(analysis_data),
  n_exercise = sum(analysis_data$treatment_group == "exercise"),
  n_control = sum(analysis_data$treatment_group == "control"),
  n_events = sum(analysis_data$dfs_event),
  hr = overall_summary$conf.int[1, 1],
  hr_lower = overall_summary$conf.int[1, 3],
  hr_upper = overall_summary$conf.int[1, 4],
  p_value = overall_summary$coefficients[1, 5],
  interaction_p = NA_real_,
  subgroup_label = "Overall"
)

all_results <- dplyr::bind_rows(overall_result, dfs_subgroup_results)

# ------------------------------------------------------------------------------
# Create Forest Plot
# ------------------------------------------------------------------------------

# Prepare data for forest plot
forest_data <- all_results |>
  dplyr::filter(!is.na(hr)) |>
  dplyr::mutate(
    display_label = dplyr::if_else(subgroup == "Overall", "Overall", paste0("  ", level)),
    hr_text = sprintf("%.2f (%.2f-%.2f)", hr, hr_lower, hr_upper)
  )

# Create plotting order
plot_levels <- c()
plot_levels <- c(plot_levels, "Overall")
for (sg in subgroups) {
  plot_levels <- c(plot_levels, subgroup_labels[sg])
  sg_levels <- forest_data |>
    dplyr::filter(subgroup == sg) |>
    dplyr::pull(display_label)
  plot_levels <- c(plot_levels, sg_levels)
}

# Add header rows
header_rows <- data.frame(
  display_label = subgroup_labels,
  hr = NA_real_, hr_lower = NA_real_, hr_upper = NA_real_,
  is_header = TRUE
)

forest_plot_data <- forest_data |>
  dplyr::mutate(is_header = FALSE) |>
  dplyr::bind_rows(header_rows) |>
  dplyr::mutate(display_label = factor(display_label, levels = rev(plot_levels)))

# Create forest plot
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

forest_plot <- ggplot(
  forest_plot_data |> dplyr::filter(!is.na(hr)),
  aes(x = hr, y = display_label)
) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = hr_lower, xmax = hr_upper), height = 0.25, linewidth = 0.7) +
  geom_point(size = 3.5, shape = 18, color = "steelblue") +
  scale_x_log10(
    breaks = c(0.25, 0.5, 0.75, 1, 1.5, 2),
    limits = c(0.2, 2.5),
    labels = c("0.25", "0.5", "0.75", "1", "1.5", "2")
  ) +
  labs(
    x = "\nHazard Ratio (95% CI)",
    y = "",
    title = "Forest Plot: Disease-Free Survival by Subgroup",
    subtitle = "Exercise vs Control"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(hjust = 0),
    plot.title = element_text(face = "bold")
  ) +
  annotate("text", x = 0.3, y = 0.5, label = "Favors Exercise", hjust = 0, size = 3) +
  annotate("text", x = 1.8, y = 0.5, label = "Favors Control", hjust = 1, size = 3)

ggsave("output/figures/forest_plot_dfs.png", forest_plot, width = 10, height = 12, dpi = 300)
ggsave("output/figures/forest_plot_dfs.pdf", forest_plot, width = 10, height = 12)

# ------------------------------------------------------------------------------
# OS Subgroup Analysis
# ------------------------------------------------------------------------------

os_subgroup_results <- purrr::map_dfr(subgroups, function(sg) {
  calculate_subgroup_hr(analysis_data, sg, outcome = "os")
})

os_interaction_pvalues <- purrr::map_dbl(subgroups, function(sg) {
  test_interaction(analysis_data, sg, outcome = "os")
})

os_subgroup_results <- os_subgroup_results |>
  dplyr::mutate(
    interaction_p = os_interaction_pvalues[subgroup],
    subgroup_label = subgroup_labels[subgroup]
  )

# Overall OS
overall_cox_os <- survival::coxph(Surv(os_time, os_event) ~ treatment_group, data = analysis_data)
overall_summary_os <- summary(overall_cox_os)

overall_result_os <- data.frame(
  subgroup = "Overall", level = "All Patients",
  n_total = nrow(analysis_data),
  n_exercise = sum(analysis_data$treatment_group == "exercise"),
  n_control = sum(analysis_data$treatment_group == "control"),
  n_events = sum(analysis_data$os_event),
  hr = overall_summary_os$conf.int[1, 1],
  hr_lower = overall_summary_os$conf.int[1, 3],
  hr_upper = overall_summary_os$conf.int[1, 4],
  p_value = overall_summary_os$coefficients[1, 5],
  interaction_p = NA_real_,
  subgroup_label = "Overall"
)

all_results_os <- dplyr::bind_rows(overall_result_os, os_subgroup_results)

# ------------------------------------------------------------------------------
# Create Subgroup Summary Table
# ------------------------------------------------------------------------------

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

subgroup_table <- all_results |>
  dplyr::filter(!is.na(hr)) |>
  dplyr::mutate(
    hr_ci = sprintf("%.2f (%.2f-%.2f)", hr, hr_lower, hr_upper),
    p_fmt = sprintf("%.3f", p_value),
    int_p_fmt = dplyr::if_else(is.na(interaction_p), "-", sprintf("%.3f", interaction_p))
  ) |>
  dplyr::select(
    Subgroup = subgroup_label, Level = level,
    N = n_total, Events = n_events,
    `HR (95% CI)` = hr_ci, `P-value` = p_fmt, `Interaction P` = int_p_fmt
  )

write.csv(subgroup_table, "output/tables/subgroup_analysis_dfs.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# Generate Text Report
# ------------------------------------------------------------------------------

dir.create("output/reports", recursive = TRUE, showWarnings = FALSE)

report_lines <- c(
  paste(rep("=", 70), collapse = ""),
  "SUBGROUP ANALYSIS RESULTS",
  "Disease-Free Survival by Pre-specified Subgroups",
  paste("Generated:", Sys.time()),
  paste(rep("=", 70), collapse = ""),
  "",
  "METHODS",
  "-------",
  "Cox proportional hazards regression estimated hazard ratios (HR)",
  "and 95% confidence intervals for exercise vs control within subgroups.",
  "Interaction p-values from likelihood ratio tests.",
  "",
  paste(rep("=", 70), collapse = ""),
  "OVERALL RESULT",
  paste(rep("=", 70), collapse = ""),
  sprintf("N = %d patients", overall_result$n_total),
  sprintf("Events = %d", overall_result$n_events),
  sprintf("HR = %.2f (95%% CI: %.2f - %.2f)", overall_result$hr,
          overall_result$hr_lower, overall_result$hr_upper),
  sprintf("P-value = %.4f", overall_result$p_value),
  "",
  paste(rep("=", 70), collapse = ""),
  "SUBGROUP RESULTS",
  paste(rep("=", 70), collapse = "")
)

for (sg in subgroups) {
  sg_results <- all_results |> dplyr::filter(subgroup == sg)
  report_lines <- c(report_lines, "", paste0(subgroup_labels[sg], ":"),
                    paste(rep("-", nchar(subgroup_labels[sg]) + 1), collapse = ""))

  for (i in seq_len(nrow(sg_results))) {
    row <- sg_results[i, ]
    if (!is.na(row$hr)) {
      report_lines <- c(report_lines, sprintf(
        "  %s: HR = %.2f (%.2f - %.2f), N = %d, Events = %d",
        row$level, row$hr, row$hr_lower, row$hr_upper, row$n_total, row$n_events
      ))
    }
  }

  int_p <- interaction_pvalues[sg]
  if (!is.na(int_p)) {
    report_lines <- c(report_lines, sprintf("  Interaction p-value: %.4f", int_p))
  }
}

report_lines <- c(report_lines, "",
  paste(rep("=", 70), collapse = ""),
  "OS SUBGROUP SUMMARY",
  paste(rep("=", 70), collapse = ""),
  sprintf("Overall OS: HR = %.2f (95%% CI: %.2f - %.2f)",
          overall_result_os$hr, overall_result_os$hr_lower, overall_result_os$hr_upper),
  "",
  paste(rep("=", 70), collapse = ""),
  "FILES GENERATED",
  paste(rep("=", 70), collapse = ""),
  "- output/figures/forest_plot_dfs.png",
  "- output/figures/forest_plot_dfs.pdf",
  "- output/tables/subgroup_analysis_dfs.csv",
  paste(rep("=", 70), collapse = "")
)

writeLines(report_lines, "output/reports/subgroup_analysis_result.txt")

# Save results
saveRDS(all_results, "output/subgroup_results_dfs.rds")
saveRDS(all_results_os, "output/subgroup_results_os.rds")

cat("\n=== Subgroup Analysis Complete ===\n")
cat("Generated files:\n")
cat("  - output/figures/forest_plot_dfs.png\n")
cat("  - output/figures/forest_plot_dfs.pdf\n")
cat("  - output/tables/subgroup_analysis_dfs.csv\n")
cat("  - output/reports/subgroup_analysis_result.txt\n")
