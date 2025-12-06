# ==============================================================================
# 03_longitudinal_analysis.R
# Breast Cancer Exercise Trial - Longitudinal & Repeated Measures Analysis
# ==============================================================================

library(tidyverse)
library(nlme)
library(lme4)
library(broom)
library(broom.mixed)
library(gtsummary)
library(ggplot2)
library(patchwork)

set.seed(20241206)

# ==============================================================================
# Load Data
# ==============================================================================

breast_data <- readRDS("data/breast_cancer_exercise_trial.rds")
longitudinal_data <- readRDS("data/breast_cancer_longitudinal.rds")

# Prepare longitudinal data
long_data <- longitudinal_data |>
  dplyr::mutate(
    treatment_group = factor(
      treatment_group,
      levels = c("control", "exercise"),
      labels = c("Control", "Exercise")
    ),
    time_factor = factor(time_months),
    time_years = time_months / 12,
    tumor_stage = factor(tumor_stage)
  )

message(sprintf("Loaded %d longitudinal observations", nrow(long_data)))

# ==============================================================================
# 1. MVPA (Physical Activity) Analysis
# ==============================================================================

message("\n=== Analyzing MVPA (Physical Activity) ===")

# Linear mixed effects model for MVPA
mvpa_model <- lme(
  mvpa_met_hours ~ treatment_group * time_months + baseline_mvpa + age + bmi,
  random = ~ time_months | patient_id,
  data = long_data,
  method = "REML",
  control = lmeControl(opt = "optim")
)

# Extract model summary
mvpa_summary <- summary(mvpa_model)
mvpa_tidy <- broom.mixed::tidy(mvpa_model, effects = "fixed", conf.int = TRUE)

# Save MVPA results
write.csv(mvpa_tidy, "output/reports/mvpa_mixed_model_result.csv", row.names = FALSE)

# Create MVPA trajectory plot
mvpa_means <- long_data |>
  dplyr::group_by(treatment_group, time_months) |>
  dplyr::summarise(
    mean_mvpa = mean(mvpa_met_hours, na.rm = TRUE),
    se_mvpa = sd(mvpa_met_hours, na.rm = TRUE) / sqrt(dplyr::n()),
    n = dplyr::n(),
    .groups = "drop"
  )

p_mvpa <- ggplot(
  mvpa_means,
  aes(x = time_months, y = mean_mvpa, color = treatment_group, fill = treatment_group)
) +
  geom_ribbon(
    aes(ymin = mean_mvpa - 1.96 * se_mvpa, ymax = mean_mvpa + 1.96 * se_mvpa),
    alpha = 0.2, color = NA
  ) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8")) +
  scale_x_continuous(breaks = c(0, 6, 12, 18, 24, 36)) +
  labs(
    title = "MVPA Over Time by Treatment Group",
    x = "Time (months)",
    y = "MVPA (MET-hours/week)",
    color = "Treatment",
    fill = "Treatment"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("output/figures/fig3_mvpa_trajectory.png", p_mvpa,
       width = 8, height = 6, dpi = 300)

# ==============================================================================
# 2. VO2max Analysis
# ==============================================================================

message("\n=== Analyzing VO2max ===")

vo2max_model <- lme(
  vo2max_predicted ~ treatment_group * time_months + baseline_vo2max + age,
  random = ~ time_months | patient_id,
  data = long_data,
  method = "REML",
  control = lmeControl(opt = "optim")
)

vo2max_tidy <- broom.mixed::tidy(vo2max_model, effects = "fixed", conf.int = TRUE)
write.csv(vo2max_tidy, "output/reports/vo2max_mixed_model_result.csv", row.names = FALSE)

# VO2max trajectory plot
vo2max_means <- long_data |>
  dplyr::group_by(treatment_group, time_months) |>
  dplyr::summarise(
    mean_vo2max = mean(vo2max_predicted, na.rm = TRUE),
    se_vo2max = sd(vo2max_predicted, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop"
  )

p_vo2max <- ggplot(
  vo2max_means,
  aes(x = time_months, y = mean_vo2max, color = treatment_group, fill = treatment_group)
) +
  geom_ribbon(
    aes(ymin = mean_vo2max - 1.96 * se_vo2max, ymax = mean_vo2max + 1.96 * se_vo2max),
    alpha = 0.2, color = NA
  ) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8")) +
  scale_x_continuous(breaks = c(0, 6, 12, 18, 24, 36)) +
  labs(
    title = "Predicted VO2max Over Time by Treatment Group",
    x = "Time (months)",
    y = "VO2max (ml/kg/min)",
    color = "Treatment",
    fill = "Treatment"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("output/figures/fig4_vo2max_trajectory.png", p_vo2max,
       width = 8, height = 6, dpi = 300)

# ==============================================================================
# 3. SF-36 Physical Functioning Analysis
# ==============================================================================

message("\n=== Analyzing SF-36 Physical Functioning ===")

sf36_model <- lme(
  sf36_physical_functioning ~ treatment_group * time_months + age + tumor_stage,
  random = ~ time_months | patient_id,
  data = long_data,
  method = "REML",
  control = lmeControl(opt = "optim"),
  na.action = na.omit
)

sf36_tidy <- broom.mixed::tidy(sf36_model, effects = "fixed", conf.int = TRUE)
write.csv(sf36_tidy, "output/reports/sf36_pf_mixed_model_result.csv", row.names = FALSE)

# SF-36 trajectory plot
sf36_means <- long_data |>
  dplyr::filter(!is.na(sf36_physical_functioning)) |>
  dplyr::group_by(treatment_group, time_months) |>
  dplyr::summarise(
    mean_sf36 = mean(sf36_physical_functioning, na.rm = TRUE),
    se_sf36 = sd(sf36_physical_functioning, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop"
  )

p_sf36 <- ggplot(
  sf36_means,
  aes(x = time_months, y = mean_sf36, color = treatment_group, fill = treatment_group)
) +
  geom_ribbon(
    aes(ymin = mean_sf36 - 1.96 * se_sf36, ymax = mean_sf36 + 1.96 * se_sf36),
    alpha = 0.2, color = NA
  ) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8")) +
  scale_x_continuous(breaks = c(0, 6, 12, 18, 24, 36)) +
  scale_y_continuous(limits = c(50, 80)) +
  labs(
    title = "SF-36 Physical Functioning Over Time",
    x = "Time (months)",
    y = "SF-36 Physical Functioning Score",
    color = "Treatment",
    fill = "Treatment"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("output/figures/fig5_sf36_pf_trajectory.png", p_sf36,
       width = 8, height = 6, dpi = 300)

# ==============================================================================
# 4. Fatigue Score Analysis
# ==============================================================================

message("\n=== Analyzing Fatigue Score ===")

fatigue_model <- lme(
  fatigue_score ~ treatment_group * time_months + age + tumor_stage,
  random = ~ time_months | patient_id,
  data = long_data,
  method = "REML",
  control = lmeControl(opt = "optim"),
  na.action = na.omit
)

fatigue_tidy <- broom.mixed::tidy(fatigue_model, effects = "fixed", conf.int = TRUE)
write.csv(fatigue_tidy, "output/reports/fatigue_mixed_model_result.csv", row.names = FALSE)

# Fatigue trajectory plot
fatigue_means <- long_data |>
  dplyr::filter(!is.na(fatigue_score)) |>
  dplyr::group_by(treatment_group, time_months) |>
  dplyr::summarise(
    mean_fatigue = mean(fatigue_score, na.rm = TRUE),
    se_fatigue = sd(fatigue_score, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop"
  )

p_fatigue <- ggplot(
  fatigue_means,
  aes(x = time_months, y = mean_fatigue, color = treatment_group, fill = treatment_group)
) +
  geom_ribbon(
    aes(ymin = mean_fatigue - 1.96 * se_fatigue, ymax = mean_fatigue + 1.96 * se_fatigue),
    alpha = 0.2, color = NA
  ) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8")) +
  scale_x_continuous(breaks = c(0, 6, 12, 18, 24, 36)) +
  labs(
    title = "Fatigue Score Over Time (Higher = Less Fatigue)",
    x = "Time (months)",
    y = "FACIT-Fatigue Score",
    color = "Treatment",
    fill = "Treatment"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("output/figures/fig6_fatigue_trajectory.png", p_fatigue,
       width = 8, height = 6, dpi = 300)

# ==============================================================================
# 5. Combined PRO Panel Figure
# ==============================================================================

message("\n=== Creating Combined PRO Panel ===")

# Quality of Life plot
qol_means <- long_data |>
  dplyr::filter(!is.na(quality_of_life_global)) |>
  dplyr::group_by(treatment_group, time_months) |>
  dplyr::summarise(
    mean_qol = mean(quality_of_life_global, na.rm = TRUE),
    se_qol = sd(quality_of_life_global, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop"
  )

p_qol <- ggplot(
  qol_means,
  aes(x = time_months, y = mean_qol, color = treatment_group, fill = treatment_group)
) +
  geom_ribbon(
    aes(ymin = mean_qol - 1.96 * se_qol, ymax = mean_qol + 1.96 * se_qol),
    alpha = 0.2, color = NA
  ) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8")) +
  scale_x_continuous(breaks = c(0, 12, 24, 36)) +
  labs(title = "Quality of Life", x = "Time (months)", y = "Score") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none")

# Anxiety/Depression plot
hads_means <- long_data |>
  dplyr::filter(!is.na(anxiety_depression)) |>
  dplyr::group_by(treatment_group, time_months) |>
  dplyr::summarise(
    mean_hads = mean(anxiety_depression, na.rm = TRUE),
    se_hads = sd(anxiety_depression, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop"
  )

p_hads <- ggplot(
  hads_means,
  aes(x = time_months, y = mean_hads, color = treatment_group, fill = treatment_group)
) +
  geom_ribbon(
    aes(ymin = mean_hads - 1.96 * se_hads, ymax = mean_hads + 1.96 * se_hads),
    alpha = 0.2, color = NA
  ) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8")) +
  scale_x_continuous(breaks = c(0, 12, 24, 36)) +
  labs(title = "Anxiety/Depression (Lower = Better)", x = "Time (months)", y = "HADS Score") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none")

# Body weight plot
weight_means <- long_data |>
  dplyr::group_by(treatment_group, time_months) |>
  dplyr::summarise(
    mean_weight = mean(body_weight, na.rm = TRUE),
    se_weight = sd(body_weight, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop"
  )

p_weight <- ggplot(
  weight_means,
  aes(x = time_months, y = mean_weight, color = treatment_group, fill = treatment_group)
) +
  geom_ribbon(
    aes(ymin = mean_weight - 1.96 * se_weight, ymax = mean_weight + 1.96 * se_weight),
    alpha = 0.2, color = NA
  ) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8")) +
  scale_x_continuous(breaks = c(0, 12, 24, 36)) +
  labs(title = "Body Weight", x = "Time (months)", y = "Weight (kg)") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none")

# 6-min walk distance plot
walk_means <- long_data |>
  dplyr::group_by(treatment_group, time_months) |>
  dplyr::summarise(
    mean_walk = mean(six_min_walk_distance, na.rm = TRUE),
    se_walk = sd(six_min_walk_distance, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop"
  )

p_walk <- ggplot(
  walk_means,
  aes(x = time_months, y = mean_walk, color = treatment_group, fill = treatment_group)
) +
  geom_ribbon(
    aes(ymin = mean_walk - 1.96 * se_walk, ymax = mean_walk + 1.96 * se_walk),
    alpha = 0.2, color = NA
  ) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8")) +
  scale_x_continuous(breaks = c(0, 12, 24, 36)) +
  labs(title = "6-Minute Walk Distance", x = "Time (months)", y = "Distance (m)") +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom")

# Combine plots
combined_plot <- (p_sf36 + p_fatigue) / (p_qol + p_hads) / (p_weight + p_walk) +
  plot_annotation(
    title = "Patient-Reported Outcomes and Physical Measures Over Time",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

ggsave("output/figures/fig7_pro_panel.png", combined_plot,
       width = 12, height = 14, dpi = 300)
ggsave("output/figures/fig7_pro_panel.pdf", combined_plot,
       width = 12, height = 14)

# ==============================================================================
# 6. Generate Longitudinal Analysis Report
# ==============================================================================

message("\n=== Generating Longitudinal Analysis Report ===")

sink("output/reports/longitudinal_analysis_result.txt")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("LONGITUDINAL / REPEATED MEASURES ANALYSIS RESULTS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

cat("1. MVPA (Moderate-to-Vigorous Physical Activity)\n")
cat("-" |> rep(50) |> paste(collapse = ""), "\n")
cat("Model: mvpa_met_hours ~ treatment_group * time_months + baseline + age + bmi\n")
cat("Random effects: ~ time_months | patient_id\n\n")
print(summary(mvpa_model))

cat("\n\nFixed Effects Coefficients:\n")
print(mvpa_tidy)

cat("\n\n2. VO2max (Cardiorespiratory Fitness)\n")
cat("-" |> rep(50) |> paste(collapse = ""), "\n")
cat("Model: vo2max_predicted ~ treatment_group * time_months + baseline + age\n")
cat("Random effects: ~ time_months | patient_id\n\n")
print(summary(vo2max_model))

cat("\n\nFixed Effects Coefficients:\n")
print(vo2max_tidy)

cat("\n\n3. SF-36 Physical Functioning\n")
cat("-" |> rep(50) |> paste(collapse = ""), "\n")
cat("Model: sf36_physical_functioning ~ treatment_group * time_months + age + stage\n")
cat("Random effects: ~ time_months | patient_id\n\n")
print(summary(sf36_model))

cat("\n\nFixed Effects Coefficients:\n")
print(sf36_tidy)

cat("\n\n4. FACIT-Fatigue Score\n")
cat("-" |> rep(50) |> paste(collapse = ""), "\n")
cat("Model: fatigue_score ~ treatment_group * time_months + age + stage\n")
cat("Random effects: ~ time_months | patient_id\n\n")
print(summary(fatigue_model))

cat("\n\nFixed Effects Coefficients:\n")
print(fatigue_tidy)

cat("\n\n5. INTERPRETATION SUMMARY\n")
cat("-" |> rep(50) |> paste(collapse = ""), "\n")

# Extract interaction effects
mvpa_interaction <- mvpa_tidy |>
  dplyr::filter(term == "treatment_groupExercise:time_months")
vo2max_interaction <- vo2max_tidy |>
  dplyr::filter(term == "treatment_groupExercise:time_months")
sf36_interaction <- sf36_tidy |>
  dplyr::filter(term == "treatment_groupExercise:time_months")
fatigue_interaction <- fatigue_tidy |>
  dplyr::filter(term == "treatment_groupExercise:time_months")

cat("\nTreatment x Time Interaction Effects:\n\n")

cat(sprintf("MVPA: %.3f MET-hours/week per month (95%% CI: %.3f to %.3f, p = %.4f)\n",
            mvpa_interaction$estimate,
            mvpa_interaction$conf.low,
            mvpa_interaction$conf.high,
            mvpa_interaction$p.value))

cat(sprintf("VO2max: %.3f ml/kg/min per month (95%% CI: %.3f to %.3f, p = %.4f)\n",
            vo2max_interaction$estimate,
            vo2max_interaction$conf.low,
            vo2max_interaction$conf.high,
            vo2max_interaction$p.value))

cat(sprintf("SF-36 PF: %.3f points per month (95%% CI: %.3f to %.3f, p = %.4f)\n",
            sf36_interaction$estimate,
            sf36_interaction$conf.low,
            sf36_interaction$conf.high,
            sf36_interaction$p.value))

cat(sprintf("Fatigue: %.3f points per month (95%% CI: %.3f to %.3f, p = %.4f)\n",
            fatigue_interaction$estimate,
            fatigue_interaction$conf.low,
            fatigue_interaction$conf.high,
            fatigue_interaction$p.value))

sink()

message("Longitudinal analysis report saved to output/reports/longitudinal_analysis_result.txt")
message("\n=== Longitudinal Analysis Complete ===")
