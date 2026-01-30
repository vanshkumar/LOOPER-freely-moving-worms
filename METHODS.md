# METHODS.md

Shared methods for the experiments in `EXPERIMENTS.md`.

This is intentionally **research/prototype** oriented:
- prefer “works and informative” over “perfect and general,”
- keep the first implementation minimal,
- save intermediate artifacts so you can iterate.

**Source of truth:** the MATLAB scripts in this repo are authoritative. This document is kept in sync with the current code; if anything diverges, trust the scripts.

---

## 1) How this maps to EXPERIMENTS

`EXPERIMENTS.md` lays out the progression:
1. **Positive control (Kato, immobilized):** confirm LOOPER recovers loop‑like dynamics.
2. **Target test (Atanas, freely moving):** run the same protocol on freely moving worms.
3. **Decision point:** proceed to perturbation only if the loop scaffold is stable.

This document is the **implementation blueprint** for those steps: shared pipeline, metrics, and dataset‑specific settings.

---

## 2) Shared pipeline (both datasets)

**Step A — Load data**
- Each worm provides `RawData` (neurons × time), `uid`, `dt_sec`, `num_neurons`, `T`.

**Step B — Preprocess (minimal)**
- Avoid double‑smoothing or double‑z‑scoring.
- Optional detrending is applied **only on the training half** when enabled.

**Step C — Fit LOOPER**
- **Fidelity:** train on the full trace.
- **Stationarity:** train on the **first half** (pre), then project the **full trace**.

**Step D — Evaluate**
- Use LOOPER’s scaffold (`BestStateMap`, `BestLoopAssignments`) to derive:
  - loop ID α and phase bin θ,
  - distance to scaffold `d(t)`,
  - split‑based drift metrics and phase continuity metrics.

This protocol is identical across Kato and Atanas; only the **dataset‑specific preprocessing/params** differ (see Section 5).

---

## 3) Evaluation logic and metrics

We run two tests for **both datasets**:

### 3.1 Fidelity test (full‑trace fit)
- Train LOOPER on the full trace.
- Use diagnostics + summary metrics to judge whether loops are recovered.

**Primary indicators:**
- `recon_corr_full` (correlation between reconstructed and embedded stream; **not** true R²)
- Stable loop assignment (few switches relative to trace length)
- Smooth phase progression (`phase_frac_small` high, `phase_var` low)

### 3.2 Stationarity test (split‑half)
- Train on the first half, project the full trace.
- Evaluate whether the learned scaffold remains valid in the second half.

**Primary indicators:**
- Small post‑split drift (`mean_post` near 0, low `post_slope`, low `d_peak`)
- Phase continuity remains stable in the post half
- `recon_corr_post` is not drastically worse than pre‑half reconstruction

These are **qualitative criteria** (no fixed thresholds). The stationarity test is intentionally strict; failure here does not necessarily rule out loop‑like structure under trial‑based or windowed validation.

### 3.3 Summary schema (identical across datasets)

Each eval writes a `summary.csv` row with:
- split stability: `mean_pre`, `mean_post`, `post_slope`, `d_peak`, `split_idx`
- loop stability: `unique_loops`, `loop_switches`, `segments`
- phase continuity: `phase_frac_small`, `phase_var`, `phase_speed_bins_per_min`, `median_d`, `mean_segment_len`
- reconstruction fidelity: `recon_corr_full`
- stationarity hold‑out fidelity: `recon_corr_post` (stationarity only)
- diagnostics: `validation_score_mean`, `validation_score_std` (LOOPER internal)

Non‑applicable fields are recorded as `NaN` to keep the schema identical.

---

## 4) Projection + distance‑to‑scaffold details

We project raw data onto the learned scaffold using **empirical state means** from `FinalStream` (OSF‑style):
- For each scaffold bin `(α, θ)`, compute the mean in `FinalStream` for the training points assigned to that bin.
- Use per‑state STD (fallback to global) to compute a **diagonal z‑score distance**.
- Use the same metric for both assignment and `d(t)`.

For split alignment, `compute_delta_d`:
- trims indices by `DelayCount*DelayTime`,
- clamps pre/post indices to `[1, procLen]`,
- enforces a strictly pre‑split baseline window.

---

## 5) Dataset‑specific settings (intentional divergences)

**Atanas (freely moving, baseline):**
- `ZScore = false` (already z‑scored in JSONs)
- Delay embedding scaled to dt (~0.6 s): `DelayTime = 4`, `DelayCount = 5`  
  (chosen via `atanas_autocorr_delay.m`, using the paper’s autocorrelation rule)
- `MinReturnTime = 4`, `MaxCheckTime = 4` (scaled to dt)
- **Detrending disabled** by default

**Kato (immobilized, positive control):**
- `ZScore = true`
- Paper defaults (dt ~0.34 s): `DelayTime = 10`, `DelayCount = 10`
- `MinReturnTime = 10`, `MaxCheckTime = 10`
- **Detrending enabled** (per‑neuron linear detrend)

Both datasets use the same fidelity + stationarity protocol and identical evaluation metrics.

---

## 6) Kato shared‑neuron run (paper‑style fidelity)

The shared run concatenates all 5 Kato worms on the **intersection of identified neuron labels**:
- per‑worm detrend + per‑worm z‑score **before** concatenation
- LOOPER `ZScore` is disabled to avoid double normalization
- `TrialData` encodes worm ID
- Uses `MinimumReturnTime = 50` (matches OSF worm checkpoint) and `MaxCheckTime = 1` to prevent memory blow‑ups on long traces

This is treated as the closest approximation to the paper’s Fig‑5B reconstruction metric.

---

## 7) Outputs & reproducibility

Always save:
- per‑run `*.mat` with `worm` + `saveData`
- per‑eval `summary.csv`
- diagnostic figures in `results/<run>/diagnostics/`

Scripts use fixed filenames (no timestamps) and overwrite old results.
Keep `results/` in sync with the scripts you intend to cite.
