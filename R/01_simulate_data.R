# ==============================================================================
# 01_simulate_data.R
# Breast Cancer Exercise Intervention Trial - Data Simulation
# ==============================================================================

# Load required packages
library(tidyverse)
library(lubridate)

# Import .data pronoun for tidy evaluation
utils::globalVariables(".data")

set.seed(20241206)

# ==============================================================================
# Parameters
# ==============================================================================

n_patients <- 1000
n_centers <- 8

# Baseline event rates (3-year DFS)
p_control_3y <- 0.82
p_exercise_3y <- 0.87
hr_target <- 0.75

# Follow-up duration (months)
max_followup <- 48

# ==============================================================================
# Generate Baseline Demographics and Clinical Characteristics
# ==============================================================================

generate_baseline_data <- function(n) {
  patient_id <- sprintf("PT%04d", 1:n)

  # Demographics
  age <- round(rnorm(n, mean = 52, sd = 10))
  age <- pmax(35, pmin(75, age))

  sex <- sample(c("Female", "Male"), n, replace = TRUE, prob = c(0.99, 0.01))

  bmi <- round(rnorm(n, mean = 26, sd = 5), 1)
  bmi <- pmax(18, pmin(45, bmi))

  # Menopause status (correlated with age)
  menopause_prob <- pnorm(age, mean = 50, sd = 3)
  menopause_status <- ifelse(
    runif(n) < menopause_prob,
    "postmenopausal",
    "premenopausal"
  )
  menopause_status[sex == "Male"] <- NA

  # Tumor characteristics
  tumor_stage <- sample(
    c("I", "II", "III"), n,
    replace = TRUE, prob = c(0.40, 0.45, 0.15)
  )
  tumor_grade <- sample(1:3, n, replace = TRUE, prob = c(0.25, 0.50, 0.25))

  # Hormone receptor status
  er_status <- sample(
    c("positive", "negative"), n,
    replace = TRUE, prob = c(0.75, 0.25)
  )
  pr_status <- ifelse(
    er_status == "positive",
    sample(c("positive", "negative"), n, replace = TRUE, prob = c(0.85, 0.15)),
    sample(c("positive", "negative"), n, replace = TRUE, prob = c(0.20, 0.80))
  )

  hormone_receptor_status <- dplyr::case_when(
    er_status == "positive" | pr_status == "positive" ~ "HR+",
    TRUE ~ "HR-"
  )

  her2_status <- sample(
    c("positive", "negative"), n,
    replace = TRUE, prob = c(0.20, 0.80)
  )

  # Molecular subtype
  molecular_subtype <- dplyr::case_when(
    hormone_receptor_status == "HR+" &
      her2_status == "negative" &
      tumor_grade <= 2 ~ "Luminal A",
    hormone_receptor_status == "HR+" &
      her2_status == "negative" &
      tumor_grade == 3 ~ "Luminal B (HER2-)",
    hormone_receptor_status == "HR+" &
      her2_status == "positive" ~ "Luminal B (HER2+)",
    hormone_receptor_status == "HR-" &
      her2_status == "positive" ~ "HER2-enriched",
    TRUE ~ "Triple-negative"
  )

  # Lymph nodes
  lymph_nodes_positive <- rpois(
    n,
    lambda = ifelse(
      tumor_stage == "I", 0.5,
      ifelse(tumor_stage == "II", 2, 5)
    )
  )

  # Surgery type
  surgery_type <- ifelse(
    tumor_stage == "I",
    sample(
      c("lumpectomy", "mastectomy"), n,
      replace = TRUE, prob = c(0.70, 0.30)
    ),
    sample(
      c("lumpectomy", "mastectomy"), n,
      replace = TRUE, prob = c(0.40, 0.60)
    )
  )

  # Treatment
  chemo_regimen <- sample(
    c("AC-T", "TAC", "TC", "Other"), n,
    replace = TRUE, prob = c(0.40, 0.25, 0.25, 0.10)
  )

  radiation_therapy <- ifelse(
    surgery_type == "lumpectomy",
    sample(c("yes", "no"), n, replace = TRUE, prob = c(0.95, 0.05)),
    sample(c("yes", "no"), n, replace = TRUE, prob = c(0.50, 0.50))
  )

  endocrine_therapy <- dplyr::case_when(
    hormone_receptor_status == "HR-" ~ "none",
    menopause_status == "premenopausal" ~ sample(
      c("tamoxifen", "AI", "none"), n,
      replace = TRUE, prob = c(0.80, 0.10, 0.10)
    ),
    TRUE ~ sample(
      c("tamoxifen", "AI", "none"), n,
      replace = TRUE, prob = c(0.30, 0.60, 0.10)
    )
  )

  targeted_therapy <- ifelse(
    her2_status == "positive",
    sample(c("yes", "no"), n, replace = TRUE, prob = c(0.90, 0.10)),
    "no"
  )

  # Performance status and comorbidities
  ecog_ps <- sample(0:1, n, replace = TRUE, prob = c(0.70, 0.30))

  diabetes_flag <- sample(
    c("yes", "no"), n,
    replace = TRUE,
    prob = c(0.12, 0.88)
  )
  hypertension_flag <- sample(
    c("yes", "no"), n,
    replace = TRUE,
    prob = c(0.17, 0.83)
  )

  comorbidities <- dplyr::case_when(
    diabetes_flag == "yes" & hypertension_flag == "yes" ~
      "diabetes, hypertension",
    diabetes_flag == "yes" ~ "diabetes",
    hypertension_flag == "yes" ~ "hypertension",
    TRUE ~ "none"
  )

  # Trial center
  stratification_center <- sample(
    paste0("Center_", LETTERS[1:n_centers]), n,
    replace = TRUE
  )

  data.frame(
    patient_id = patient_id,
    age = age,
    sex = sex,
    bmi = bmi,
    menopause_status = menopause_status,
    tumor_stage = tumor_stage,
    tumor_grade = tumor_grade,
    er_status = er_status,
    pr_status = pr_status,
    hormone_receptor_status = hormone_receptor_status,
    her2_status = her2_status,
    molecular_subtype = molecular_subtype,
    lymph_nodes_positive = lymph_nodes_positive,
    surgery_type = surgery_type,
    chemo_regimen = chemo_regimen,
    radiation_therapy = radiation_therapy,
    endocrine_therapy = endocrine_therapy,
    targeted_therapy = targeted_therapy,
    ecog_ps = ecog_ps,
    comorbidities = comorbidities,
    stratification_center = stratification_center,
    stringsAsFactors = FALSE
  )
}

# ==============================================================================
# Randomization
# ==============================================================================

randomize_patients <- function(df) {
  n <- nrow(df)

  # Stratification variables
  df$stratification_stage <- df$tumor_stage
  df$stratification_bmi <- ifelse(df$bmi <= 25, "BMI<=25", "BMI>25")
  df$stratification_menopause <- ifelse(
    is.na(df$menopause_status), "male",
    df$menopause_status
  )

  # Block randomization within strata
  df$strata <- paste(
    df$stratification_stage, df$stratification_bmi,
    df$stratification_menopause,
    sep = "_"
  )

  df <- df |>
    dplyr::group_by(strata) |>
    dplyr::mutate(
      rand_order = sample(dplyr::n()),
      treatment_group = ifelse(rand_order %% 2 == 0, "exercise", "control")
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-rand_order)

  # Randomization date (spread over 24 months)
  df$randomization_date <- as.Date("2020-01-01") +
    sample(0:730, n, replace = TRUE)

  df
}

# ==============================================================================
# Generate Exercise and Physical Measurements (Longitudinal)
# ==============================================================================
generate_longitudinal_data <- function(df) {
  timepoints <- c(0, 6, 12, 18, 24, 36)
  n_pts <- nrow(df)

  long_data <- df |>
    tidyr::crossing(time_months = timepoints) |>
    dplyr::arrange(patient_id, time_months)

  # Baseline MVPA (MET-hours/week)
  baseline_mvpa_vec <- rnorm(n_pts, mean = 8, sd = 5)
  baseline_mvpa_vec <- pmax(0, baseline_mvpa_vec)

  long_data <- long_data |>
    dplyr::group_by(patient_id) |>
    dplyr::mutate(
      baseline_mvpa = baseline_mvpa_vec[match(patient_id[1], df$patient_id)],
      mvpa_met_hours = dplyr::case_when(
        time_months == 0 ~ baseline_mvpa,
        treatment_group == "exercise" ~ baseline_mvpa +
          rnorm(dplyr::n(), mean = 8 * (1 - exp(-time_months / 12)), sd = 3),
        TRUE ~ baseline_mvpa + rnorm(dplyr::n(), mean = 1, sd = 2)
      ),
      mvpa_met_hours = pmax(0, mvpa_met_hours)
    ) |>
    dplyr::ungroup()

  # Set baseline MVPA properly
  long_data <- long_data |>
    dplyr::group_by(patient_id) |>
    dplyr::mutate(baseline_mvpa = dplyr::first(mvpa_met_hours)) |>
    dplyr::ungroup()

  # VO2max predicted (ml/kg/min)
  long_data <- long_data |>
    dplyr::group_by(patient_id) |>
    dplyr::mutate(
      baseline_vo2max = 25 - 0.2 * (dplyr::first(age) - 50) + rnorm(1, 0, 3),
      vo2max_predicted = dplyr::case_when(
        time_months == 0 ~ baseline_vo2max,
        treatment_group == "exercise" ~ baseline_vo2max +
          rnorm(dplyr::n(), 2.5 * (1 - exp(-time_months / 12)), 1),
        TRUE ~ baseline_vo2max + rnorm(dplyr::n(), mean = 0.3, sd = 0.8)
      )
    ) |>
    dplyr::ungroup()

  # 6-minute walk distance (meters)
  long_data <- long_data |>
    dplyr::group_by(patient_id) |>
    dplyr::mutate(
      baseline_6mwd = 450 - 2 * (dplyr::first(age) - 50) + rnorm(1, 0, 50),
      six_min_walk_distance = dplyr::case_when(
        time_months == 0 ~ baseline_6mwd,
        treatment_group == "exercise" ~ baseline_6mwd +
          rnorm(dplyr::n(), 40 * (1 - exp(-time_months / 12)), 20),
        TRUE ~ baseline_6mwd + rnorm(dplyr::n(), mean = 5, sd = 15)
      )
    ) |>
    dplyr::ungroup()

  # Body weight
  long_data <- long_data |>
    dplyr::group_by(patient_id) |>
    dplyr::mutate(
      baseline_weight = dplyr::first(bmi) * (1.65)^2,
      body_weight = dplyr::case_when(
        time_months == 0 ~ baseline_weight,
        treatment_group == "exercise" ~ baseline_weight -
          rnorm(dplyr::n(), 2 * (1 - exp(-time_months / 18)), 1.5),
        TRUE ~ baseline_weight + rnorm(dplyr::n(), mean = 0.5, sd = 1)
      )
    ) |>
    dplyr::ungroup()

  # Waist and hip circumference
  long_data <- long_data |>
    dplyr::group_by(patient_id) |>
    dplyr::mutate(
      baseline_waist = 70 + dplyr::first(bmi) * 1.2 + rnorm(1, 0, 5),
      waist_circumference = dplyr::case_when(
        time_months == 0 ~ baseline_waist,
        treatment_group == "exercise" ~ baseline_waist -
          rnorm(dplyr::n(), 3 * (1 - exp(-time_months / 12)), 2),
        TRUE ~ baseline_waist + rnorm(dplyr::n(), mean = 0.5, sd = 1.5)
      ),
      baseline_hip = baseline_waist + rnorm(1, 10, 3),
      hip_circumference = dplyr::case_when(
        time_months == 0 ~ baseline_hip,
        treatment_group == "exercise" ~ baseline_hip -
          rnorm(dplyr::n(), 2 * (1 - exp(-time_months / 12)), 1.5),
        TRUE ~ baseline_hip + rnorm(dplyr::n(), mean = 0.3, sd = 1)
      )
    ) |>
    dplyr::ungroup()

  # Grip strength (kg)
  long_data <- long_data |>
    dplyr::group_by(patient_id) |>
    dplyr::mutate(
      baseline_grip = ifelse(dplyr::first(sex) == "Female", 25, 40) -
        0.2 * (dplyr::first(age) - 50) + rnorm(1, 0, 5),
      grip_strength = dplyr::case_when(
        time_months == 0 ~ baseline_grip,
        treatment_group == "exercise" ~ baseline_grip +
          rnorm(dplyr::n(), 3 * (1 - exp(-time_months / 12)), 1.5),
        TRUE ~ baseline_grip + rnorm(dplyr::n(), mean = 0.2, sd = 1)
      )
    ) |>
    dplyr::ungroup()

  # Sit-to-stand test
  long_data <- long_data |>
    dplyr::group_by(patient_id) |>
    dplyr::mutate(
      baseline_sts = 12 - 0.1 * (dplyr::first(age) - 50) + rnorm(1, 0, 2),
      sit_to_stand_test = dplyr::case_when(
        time_months == 0 ~ baseline_sts,
        treatment_group == "exercise" ~ baseline_sts +
          rnorm(dplyr::n(), 3 * (1 - exp(-time_months / 12)), 1),
        TRUE ~ baseline_sts + rnorm(dplyr::n(), mean = 0.3, sd = 0.8)
      ),
      sit_to_stand_test = round(pmax(0, sit_to_stand_test))
    ) |>
    dplyr::ungroup()

  long_data
}

# ==============================================================================
# Generate Adherence Data
# ==============================================================================

generate_adherence_data <- function(df) {
  exercise_pts <- df |>
    dplyr::filter(treatment_group == "exercise") |>
    dplyr::pull(patient_id) |>
    unique()
  n_exercise <- length(exercise_pts)

  # Base adherence tendency per patient
  adherence_tendency <- rbeta(n_exercise, 3, 1.5)

  adherence_df <- data.frame(
    patient_id = exercise_pts,
    adherence_tendency = adherence_tendency
  )

  adherence_df <- adherence_df |>
    dplyr::mutate(
      phase1_behavior_sessions = round(12 * rbeta(
        n_exercise, adherence_tendency * 5 + 1, (1 - adherence_tendency) * 3 + 1
      )),
      phase1_mandatory_exercise = round(12 * rbeta(
        n_exercise, adherence_tendency * 4 + 1, (1 - adherence_tendency) * 3 + 1
      )),
      phase1_optional_exercise = round(12 * rbeta(
        n_exercise, adherence_tendency * 3 + 1, (1 - adherence_tendency) * 4 + 1
      )),
      phase2_behavior_sessions = round(12 * rbeta(
        n_exercise, adherence_tendency * 4 + 1, (1 - adherence_tendency) * 4 + 1
      )),
      phase3_behavior_sessions = round(24 * rbeta(
        n_exercise, adherence_tendency * 3 + 1, (1 - adherence_tendency) * 5 + 1
      )),
      overall_adherence = (
        phase1_behavior_sessions + phase1_mandatory_exercise +
          phase2_behavior_sessions + phase3_behavior_sessions
      ) / 60 * 100
    ) |>
    dplyr::select(-adherence_tendency)

  # Control group gets NA for adherence
  control_pts <- df |>
    dplyr::filter(treatment_group == "control") |>
    dplyr::pull(patient_id) |>
    unique()

  control_adherence <- data.frame(
    patient_id = control_pts,
    phase1_behavior_sessions = NA_integer_,
    phase1_mandatory_exercise = NA_integer_,
    phase1_optional_exercise = NA_integer_,
    phase2_behavior_sessions = NA_integer_,
    phase3_behavior_sessions = NA_integer_,
    overall_adherence = NA_real_
  )

  dplyr::bind_rows(adherence_df, control_adherence)
}

# ==============================================================================
# Generate Patient-Reported Outcomes
# ==============================================================================

generate_pro_data <- function(df) {
  long_data <- df |>
    dplyr::select(patient_id, treatment_group, age, tumor_stage) |>
    dplyr::distinct() |>
    tidyr::crossing(time_months = c(0, 6, 12, 18, 24, 36))

  long_data <- long_data |>
    dplyr::group_by(patient_id) |>
    dplyr::mutate(
      baseline_sf36_pf = 65 - 0.3 * (dplyr::first(age) - 50) + rnorm(1, 0, 10),
      sf36_physical_functioning = dplyr::case_when(
        time_months == 0 ~ baseline_sf36_pf,
        treatment_group == "exercise" ~ pmin(
          100, baseline_sf36_pf +
            rnorm(dplyr::n(), 8 * (1 - exp(-time_months / 12)), 5)
        ),
        TRUE ~ pmin(100, baseline_sf36_pf + rnorm(dplyr::n(), 2, 4))
      ),
      baseline_sf36_pr = 55 + rnorm(1, 0, 15),
      sf36_physical_role = dplyr::case_when(
        time_months == 0 ~ baseline_sf36_pr,
        treatment_group == "exercise" ~ pmin(
          100, baseline_sf36_pr +
            rnorm(dplyr::n(), 10 * (1 - exp(-time_months / 12)), 8)
        ),
        TRUE ~ pmin(100, baseline_sf36_pr + rnorm(dplyr::n(), 3, 6))
      ),
      baseline_sf36_bp = 60 + rnorm(1, 0, 12),
      sf36_bodily_pain = dplyr::case_when(
        time_months == 0 ~ baseline_sf36_bp,
        treatment_group == "exercise" ~ pmin(
          100, baseline_sf36_bp +
            rnorm(dplyr::n(), 6 * (1 - exp(-time_months / 12)), 5)
        ),
        TRUE ~ pmin(100, baseline_sf36_bp + rnorm(dplyr::n(), 1, 4))
      ),
      baseline_fatigue = 35 + rnorm(1, 0, 8),
      fatigue_score = dplyr::case_when(
        time_months == 0 ~ baseline_fatigue,
        treatment_group == "exercise" ~ pmin(
          52, baseline_fatigue +
            rnorm(dplyr::n(), 5 * (1 - exp(-time_months / 12)), 3)
        ),
        TRUE ~ pmin(52, baseline_fatigue + rnorm(dplyr::n(), 1, 2.5))
      ),
      baseline_qol = 60 + rnorm(1, 0, 12),
      quality_of_life_global = dplyr::case_when(
        time_months == 0 ~ baseline_qol,
        treatment_group == "exercise" ~ pmin(
          100, baseline_qol +
            rnorm(dplyr::n(), 7 * (1 - exp(-time_months / 12)), 5)
        ),
        TRUE ~ pmin(100, baseline_qol + rnorm(dplyr::n(), 2, 4))
      ),
      baseline_hads = 12 + rnorm(1, 0, 5),
      anxiety_depression = dplyr::case_when(
        time_months == 0 ~ baseline_hads,
        treatment_group == "exercise" ~ pmax(
          0, baseline_hads -
            rnorm(dplyr::n(), 3 * (1 - exp(-time_months / 12)), 2)
        ),
        TRUE ~ pmax(0, baseline_hads - rnorm(dplyr::n(), 0.5, 1.5))
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::select(
      patient_id, time_months, sf36_physical_functioning,
      sf36_physical_role, sf36_bodily_pain, fatigue_score,
      quality_of_life_global, anxiety_depression
    )

  long_data
}

# ==============================================================================
# Generate Survival Data
# ==============================================================================

generate_survival_data <- function(df) {
  n <- nrow(df)

  # Calculate baseline hazard based on prognostic factors
  df <- df |>
    dplyr::mutate(
      risk_score =
        ifelse(tumor_stage == "II", 0.3, 0) +
        ifelse(tumor_stage == "III", 0.8, 0) +
        ifelse(tumor_grade == 2, 0.2, 0) +
        ifelse(tumor_grade == 3, 0.5, 0) +
        ifelse(hormone_receptor_status == "HR-", 0.4, 0) +
        ifelse(her2_status == "positive", 0.2, 0) +
        0.02 * lymph_nodes_positive +
        ifelse(molecular_subtype == "Triple-negative", 0.5, 0) +
        0.01 * (age - 50),
      treatment_effect = ifelse(treatment_group == "exercise", log(hr_target), 0),
      log_hazard = -4.5 + risk_score + treatment_effect
    )

  # Generate event times using exponential distribution
  df <- df |>
    dplyr::mutate(
      true_event_time = rexp(n, rate = exp(log_hazard)),
      censoring_time = runif(n, min = max_followup - 12, max = max_followup),
      dropout_time = ifelse(
        runif(n) < 0.05,
        runif(n, min = 6, max = max_followup - 6),
        Inf
      ),
      dfs_time = pmin(true_event_time, censoring_time, dropout_time),
      dfs_event = as.integer(
        true_event_time <= pmin(censoring_time, dropout_time)
      ),
      follow_up_time = dfs_time
    )

  # Event types for those with events
  df <- df |>
    dplyr::mutate(
      event_type = dplyr::case_when(
        dfs_event == 0 ~ "censored",
        runif(n) < 0.30 ~ "local_recurrence",
        runif(n) < 0.50 ~ "distant_recurrence",
        runif(n) < 0.10 ~ "second_primary_breast",
        runif(n) < 0.05 ~ "second_primary_other",
        TRUE ~ "death_without_recurrence"
      )
    )

  # Overall survival (death occurs later than recurrence for most)
  df <- df |>
    dplyr::mutate(
      death_time = dplyr::case_when(
        event_type %in% c("distant_recurrence") ~
          dfs_time + rexp(n, rate = 0.08),
        event_type %in% c("local_recurrence", "second_primary_breast") ~
          dfs_time + rexp(n, rate = 0.02),
        event_type == "death_without_recurrence" ~ dfs_time,
        TRUE ~ Inf
      ),
      os_time = pmin(death_time, censoring_time, dropout_time),
      os_event = as.integer(
        death_time <= pmin(censoring_time, dropout_time)
      ),
      death_cause = dplyr::case_when(
        os_event == 0 ~ NA_character_,
        event_type %in% c(
          "distant_recurrence", "local_recurrence",
          "death_without_recurrence"
        ) ~ "breast_cancer",
        event_type == "second_primary_other" ~ "other_cancer",
        runif(n) < 0.15 ~ "cardiovascular",
        TRUE ~ "other"
      )
    )

  df |>
    dplyr::select(
      patient_id, follow_up_time, dfs_event, dfs_time,
      event_type, os_event, os_time, death_cause
    )
}

# ==============================================================================
# Generate Adverse Events
# ==============================================================================

generate_ae_data <- function(df) {
  n <- nrow(df)

  df |>
    dplyr::mutate(
      ae_musculoskeletal = dplyr::case_when(
        treatment_group == "exercise" ~ sample(
          0:2, n, replace = TRUE, prob = c(0.70, 0.25, 0.05)
        ),
        TRUE ~ sample(0:2, n, replace = TRUE, prob = c(0.90, 0.08, 0.02))
      ),
      ae_cardiovascular = sample(
        0:3, n, replace = TRUE, prob = c(0.96, 0.03, 0.008, 0.002)
      ),
      ae_lymphedema = dplyr::case_when(
        surgery_type == "mastectomy" & treatment_group == "exercise" ~ sample(
          0:2, n, replace = TRUE, prob = c(0.80, 0.15, 0.05)
        ),
        surgery_type == "mastectomy" ~ sample(
          0:2, n, replace = TRUE, prob = c(0.82, 0.14, 0.04)
        ),
        treatment_group == "exercise" ~ sample(
          0:2, n, replace = TRUE, prob = c(0.92, 0.06, 0.02)
        ),
        TRUE ~ sample(0:2, n, replace = TRUE, prob = c(0.94, 0.05, 0.01))
      ),
      ae_other = sample(
        0:3, n, replace = TRUE, prob = c(0.85, 0.10, 0.04, 0.01)
      ),
      ae_any = pmax(ae_musculoskeletal, ae_cardiovascular, ae_lymphedema,
                    ae_other),
      ae_relationship = dplyr::case_when(
        ae_any == 0 ~ NA_character_,
        treatment_group == "control" ~ sample(
          c("unrelated", "unlikely"), n, replace = TRUE, prob = c(0.7, 0.3)
        ),
        ae_musculoskeletal > 0 ~ sample(
          c("possibly_related", "probably_related", "definitely_related"),
          n, replace = TRUE, prob = c(0.4, 0.4, 0.2)
        ),
        TRUE ~ sample(
          c("unrelated", "unlikely", "possibly_related"),
          n, replace = TRUE, prob = c(0.5, 0.3, 0.2)
        )
      )
    ) |>
    dplyr::select(
      patient_id, ae_musculoskeletal, ae_cardiovascular,
      ae_lymphedema, ae_other, ae_any, ae_relationship
    )
}

# ==============================================================================
# Main Data Generation Pipeline
# ==============================================================================

message("Generating baseline data...")
baseline_df <- generate_baseline_data(n_patients)

message("Randomizing patients...")
baseline_df <- randomize_patients(baseline_df)

message("Generating longitudinal exercise/fitness data...")
longitudinal_df <- generate_longitudinal_data(baseline_df)

message("Generating adherence data...")
adherence_df <- generate_adherence_data(baseline_df)

message("Generating patient-reported outcomes...")
pro_df <- generate_pro_data(baseline_df)

message("Generating survival data...")
survival_df <- generate_survival_data(baseline_df)

message("Generating adverse events data...")
ae_df <- generate_ae_data(baseline_df)

# ==============================================================================
# Merge All Data
# ==============================================================================

# Baseline dataset (one row per patient)
breast_data <- baseline_df |>
  dplyr::left_join(adherence_df, by = "patient_id") |>
  dplyr::left_join(survival_df, by = "patient_id") |>
  dplyr::left_join(ae_df, by = "patient_id")

# Longitudinal dataset (multiple rows per patient)
longitudinal_full <- longitudinal_df |>
  dplyr::left_join(pro_df, by = c("patient_id", "time_months"))

# ==============================================================================
# Save Data
# ==============================================================================

message("Saving datasets...")

# Main analysis dataset
readr::write_csv(breast_data, "data/breast_cancer_exercise_trial.csv")

# Longitudinal dataset
readr::write_csv(longitudinal_full, "data/breast_cancer_longitudinal.csv")

# Also save as RDS for R analysis
saveRDS(breast_data, "data/breast_cancer_exercise_trial.rds")
saveRDS(longitudinal_full, "data/breast_cancer_longitudinal.rds")

message("Data generation complete!")
message(sprintf("  - %d patients randomized", nrow(breast_data)))
message(sprintf(
  "  - %d in exercise group, %d in control group",
  sum(breast_data$treatment_group == "exercise"),
  sum(breast_data$treatment_group == "control")
))
message(sprintf(
  "  - %d DFS events (%.1f%%)",
  sum(breast_data$dfs_event),
  100 * mean(breast_data$dfs_event)
))
message(sprintf(
  "  - %d deaths (%.1f%%)",
  sum(breast_data$os_event),
  100 * mean(breast_data$os_event)
))
