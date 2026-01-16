# Fatigue Model Verification & Validation (V&V)

This repository contains the baseline fatigue simulation model and a complete set of **Verification** and **Validation** scripts used to evaluate the model’s correctness (verification) and realism/credibility (validation).

---

## Repository Contents

### Core Model
- **`baseline_model.m`**  
  The original (baseline) fatigue simulation model.

---

## Verification Methods (Model Correctness)

Verification checks whether the model logic behaves as intended under controlled test conditions.

### 1) Value Substitution
**Goal:** Substitute known/controlled inputs to confirm the model returns reasonable and expected outputs.

- **Main file:** `verification_value_substitution.m`
- **Typical output:** console results and/or plots (depends on your script)

### 2) Equilibrium State
**Goal:** Run the model long enough to see whether key states stabilise at a sensible steady condition (no exploding/unstable trends).

- **Main file:** `verification_equilibrium_state.m`
- **Typical output:** time-series plots showing steady-state behaviour

### 3) Changes in State
**Goal:** Change inputs over time (e.g., increase stress/training pressure) and confirm the model state changes in the correct direction.

- **Main file:** `verification_changes_in_state.m`
- **Typical output:** scenario-based plots comparing state changes

---

## Validation Methods (Model Credibility)

Validation checks whether the model is realistic and aligns with expected behaviour or real data patterns.

### 1) Input Sensitivity Analysis
**Goal:** Test how sensitive the model output is to **input variables** (e.g., sleep, stress, training intensity) and identify which inputs most influence fatigue.

- **Main file:** `validation_input_sensitivity_analysis.m`
- **Typical output:** sensitivity plots / ranking of influential inputs

### 2) Parameter Sensitivity Analysis
**Goal:** Test how sensitive the model output is to **model parameters** (internal tuning values) to identify which parameters require careful calibration.

- **Main file:** `validation_parameter_sensitivity_analysis.m`
- **Typical output:** parameter influence plots / key parameter list

### 3) Statistical Testing (Questionnaire Validation)
**Goal:** Compare simulated fatigue vs. questionnaire-observed fatigue using a **paired t-test**.

This validation is a 3-step pipeline:

#### Step A — Preprocess questionnaire data (Python)
- **Input:** `Questionnaire_Data.csv`
- **Script:** `process_questionnaire_data.py`
- **Output:** `Questionnaire_Data_Processed.csv`

#### Step B — Generate validation dataset (MATLAB)
- **Script:** `validation_dataset_from_questionnaire.m`
- **Input:** `Questionnaire_Data_Processed.csv`
- **Output (typical):** `validation_data.mat`  
  (contains simulated fatigue, observed fatigue, differences, sample size, etc.)

#### Step C — Run statistical paired t-test (MATLAB)
- **Script:** `validation_statistical_paired_ttest.m`
- **Input:** `validation_data.mat`
- **Output (typical):** statistical results (p-value, confidence interval) + plots and/or a text summary

---

## Method → File Map (Quick Reference)

| Category | Method | Main File(s)  | Data Needed |
|---------|--------|--------------|------------|
| Verification | Value Substitution | `verification_value_substitution.m` | None |
| Verification | Equilibrium State | `verification_equilibrium_state.m`  | None |
| Verification | Changes in State | `verification_changes_in_state.m`  | None |
| Validation | Input Sensitivity | `validation_input_sensitivity_analysis.m`  | None |
| Validation | Parameter Sensitivity | `validation_parameter_sensitivity_analysis.m`  | None |
| Validation | Statistical Testing | `process_questionnaire_data.py`, `validation_dataset_from_questionnaire.m`, `validation_statistical_paired_ttest.m` | `Questionnaire_Data.csv` |

---

## How to Run

### Requirements
- MATLAB (recommended) or GNU Octave (if compatible with your scripts)
- Python 3.x (for questionnaire processing)

### 1) Run Verification
Open MATLAB and run each script (order is flexible):
```matlab
run('verification_value_substitution.m')
run('verification_equilibrium_state.m')
run('verification_changes_in_state.m')
```

### 2) Run Validation (Sensitivity Analyses)
```
run('validation_input_sensitivity_analysis.m')
run('validation_parameter_sensitivity_analysis.m')
```

### 3) Run Validation (Statistical Testing)

A. Python preprocessing
```python process_questionnaire_data.py```

B. MATLAB: generate simulation-vs-observed dataset
```run('validation_dataset_from_questionnaire.m')```

C. MATLAB: run paired t-test
```run('validation_statistical_paired_ttest.m')```