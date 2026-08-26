# Changelog

## 2026-08-25 — v1.1.0: SimCRC v0.13.0 calibration checkpoint

### Added

**SimCRC v0.13.0 calibration with US posteriors as Chile priors** (ISS: #21)

Rewired the Chile prior derivation to start from the US v0.13.0 calibrated posteriors rather than hand-specified bounds, then modified only the Chile-specific parameters (`alpha_lesion_adenoma`, hazard rates, and preclinical-to-clinical detection probabilities) per Nicolás's April review. Produced BayCANN run `v0.13.0.20260630.1105` (Chile, Adenoma, Female): 37 parameters, 110 outputs, 4x360-node tanh emulator, 4 Stan chains at 100k iterations thinned by 100, emulator loss 0.00205 and MAE 0.029. Twelve rounds of coverage checks preceded the accepted design; six calibrated parameter sets were selected (`Max_lp`, `Max_log_likelihood`, `Min_AbsolutErr`, `Min_MSE`, `Post_mean`, `Post_median`).

Files: `analysis/01_calibration_setup.R`, `analysis/12_best_param_set.R`, `analysis/20_SimCRC_calibration_all.R`, `data-raw/inputs_priors_SimCRC_v0.13.0.20260617.182238_Adenoma_F.csv`, `data-raw/dt_calibrated_posteriors_SimCRC_v0.13.0.20260617.182238_Adenoma_F.csv`, `outputs/BayCANN_versions/Chile/Adenoma/F/v0.13.0/v0.13.0.20260630.1105/`

---

### Fixed

**Metadata-column leak inflated the calibration parameter set** (ISS: #21)

`pivot_longer(everything())` in the prior setup pulled non-parameter columns (`id_draw`, `lp__`, `chain`) into the prior set, so `x_names` no longer matched the calibrated `Xq` matrix and all three BayCANN calibration scripts failed at `colnames(df_post_ann) <- x_names`. Generalized the metadata guard across the three scripts so the parameter-count mismatch cannot recur.

Files: `analysis/06_BayCANN_calibration_all.R`, `analysis/06_BayCANN_calibration_all_k3.R`, `analysis/06_BayCANN_calibration_k3.R`

---

### Added

**CEA parallelization benchmark** (ISS: #21)

Benchmarked doSNOW/foreach against mirai at 1M population and 68 strategies on 5 cores: 252 min vs 233 min (1.08x), with mirai using 4.8 GB more peak RAM. Both ship the population table once per worker, so the speedup comes from dynamic scheduling alone — mirai's value is reproducible RNG and error surfacing, not throughput. Documents the three levers for larger runs (`optimize_memory`, Linux fork, Slurm).

Files: `analysis/benchmark_parallel_cea.R`, `analysis/cea_analysis_chile_parallel.R`, `ce_results/benchmark_parallel_cea.html`, `ce_results/benchmark_ram_traces.png`, `cea_parallelization_architecture.html`

---

**2026 Chilean screening cost workbook**

Tracked `screen_costs_CH2026.xlsx`, the source workbook behind the v3 screening cost inputs, so the derivation of the CSV inputs is reproducible from the repo.

Files: `screen_costs_CH2026.xlsx`

---

### Changed

**Raw model output directory renamed**

Renamed `output/2026Chile/Base/RawModelOutput_SimCRC_R` to `RawModelOutput_SimCRC` and regenerated the strategy sweep (68 COL/FIT strategy files) under the new path.

Files: `output/2026Chile/Base/RawModelOutput_SimCRC/`

---

## 2026-05-09 — v1.0.0: Initial release

### Added

**Calibration infrastructure** (ISS: #5)

Built the full BayCANN calibration pipeline for the Chilean SimCRC model (Adenoma, Female). Includes LHS design, coverage analysis, ANN emulator training via keras3/TensorFlow, Bayesian posterior sampling with Stan, posterior validation, and best parameter set selection. Calibrated to Chilean epidemiological targets using priors derived from the US model with Chile-specific adjustments.

Files: `R/01_model_input_functions.R`, `R/02_calibration_functions.R`, `R/02_calibration_functions_ssp.R`, `R/04_BayCANN_setup.R`, `R/05_calibration_targets.R`, `R/06_validation_functions.R`, `analysis/01_calibration_setup.R`, `analysis/01_wrangle_targets.R`, `analysis/02_LHS_design_all.R`, `analysis/03_Coverage_analysis.R`, `analysis/06_BayCANN_calibration_all.R`, `analysis/07_posterior_validations.R`, `analysis/12_best_param_set.R`, `analysis/20_SimCRC_calibration_all.R`

---

**Target analysis — USA vs Chile** (ISS: #1)

Added Chilean calibration target data and comparison analysis against US targets for both numerical (CRC incidence, adenoma prevalence) and categorical (stage/size distribution) targets.

Files: `data-raw/true_target_simcrcRvCH.csv`, `analysis/01_wrangle_targets.R`

---

**Cost-effectiveness analysis** (ISS: #11)

Implemented CEA for CRC screening strategies under Chilean cost and life-table assumptions. Includes single-thread and parallel versions, strategy definition, ICER computation, and CE frontier plots using `ceacrc` and `dampack`.

Files: `analysis/cea_analysis_chile.R`, `analysis/cea_analysis_chile_parallel.R`, `data-raw/cea_inputs/screen_costs.csv`, `data-raw/cea_inputs/screen_costs_v2.csv`, `data-raw/cea_inputs/crc_care_costs.csv`

---

**Chile vs USA prior and posterior comparison** (ISS: #18)

Added scripts to compare calibrated prior ranges and cancer progression parameters between Chile and USA models, validating that Chile-specific adjustments are sensible relative to the established US calibration.

Files: `analysis/09_compare_priors_Chile_vs_USA.R`, `analysis/10_cancer_progression_Chile_vs_USA.R`

---

**Updated Chilean treatment costs (v3)** (ISS: #18)

Implemented 2026 Chilean screening costs (+49–196% vs v2). Higher costs shift the efficient frontier: NoScreening and FIT re-enter the non-dominated set, showing cheaper screening modalities become more competitive under realistic Chilean pricing.

Files: `data-raw/cea_inputs/screen_costs_v3.csv`, `analysis/cea_analysis_chile.R`, `analysis/cea_analysis_chile_parallel.R`

---

**Calibration results slides**

Built Quarto reveal.js presentation with calibration results (priors, posteriors, validation, coverage, ANN performance) and v2-vs-v3 CEA cost comparison. Includes PDF export tooling and SimCRC template library.

Files: `SimCRC slides/calibration_results_Adenoma_F_v0.13.0.qmd`, `SimCRC slides/simcrc-template_library.qmd`, `SimCRC slides/export_pdf.sh`

---

**BayCANN calibration outputs (v0.12.x and v0.13.0)**

Versioned calibration outputs for Chile and USA (Adenoma, Female) across multiple calibration runs, including priors, posteriors, validation plots, and ANN prediction performance.

Files: `outputs/BayCANN_versions/Chile/`, `outputs/BayCANN_versions/USA/`
