# METHODS.md

Shared methods for the experiments in `EXPERIMENTS.md`.

This is intentionally **research/prototype** oriented:
- prefer “works and informative” over “perfect and general,”
- keep the first implementation minimal,
- save intermediate artifacts so you can iterate.

**Source of truth:** the MATLAB scripts in this repo are authoritative. This document is kept in sync with the current code; if anything diverges, trust the scripts.

---

## 1) LLM quick path (what matters in practice)

1. Load per‑worm data from the **raw dataset** and metadata:
   `X` (time x neurons), `neurons` (optional IDs), `dt`, and `ranges_raw` (pre/post split).
2. Apply minimal preprocessing (do **not** double‑smooth or double‑normalize).
3. Fit LOOPER per worm on the **pre** split and save core outputs:
   `BestModel`, `BestEmission`, `BestStateMap`, `FinalStream`.
4. Use `BestStateMap(:,1)` as **α** (trajectory ID), `BestStateMap(:,2)` as **θ bin**.
5. For distance‑to‑scaffold, use the emission‑space proxy (z‑scored, per‑state STD).
6. Log all params and artifacts so experiments are reproducible.

Current status: **freely moving Atanas baseline does not yield stable LOOPER‑style loops** under this split‑based evaluation; heat‑pulse path is deprecated. See `RESULTS.md`.

Experiment 1 uses the paper windows (`ranges_raw`) for pre/post segments; scripts warn if non‑contiguous, and either skip or error if missing (depending on the run).
For Atanas (dt ~0.6 s), delay embedding uses the **paper’s autocorrelation rule** (minimum across neurons of the autocorrelation time constant), which yields **DelayTime = 4 frames** (~2.4 s) and **DelayCount = 5** (~12 s total). To keep time scales aligned, **MinReturnTime** and **MaxCheckTime** are set to **4 frames** as well.
For Kato 2015 runs, we currently **detrend per neuron across time** (standard) rather than the OSF script’s column‑wise detrend (across neurons per timepoint). The OSF behavior appears to remove a global population mode, but it reduced reconstruction quality in our runs (R² dropped from ~0.64 to ~0.49), so we reverted while awaiting clarification from the authors.

Repo-specific notes:
- Atanas `trace_array` is already z‑scored; do **not** z‑score again.
- Baseline split uses `ranges_raw`; scripts warn if non‑contiguous.
- Both `experiment1_all_baseline_run` and `experiment1_single_baseline_run` keep detrending **off** by default (toggle in script if needed).

Additional takeaways relevant to Atanas runs:
- LOOPER diagnostics (`plotReconstruction`) operate on the **delay‑embedded stream** (`FinalStream`) and report **correlation**, not true R²; treat these as **sanity checks** rather than paper‑style validation metrics.
- If you concatenate segments or use pre/post windows, **mark trial boundaries** (`TrialData`) to avoid artificial transitions across discontinuities.
- `MinimumReturnTime` is a **time‑separation constraint** for neighbor selection; scale it to the dataset’s dt and the timescale of recurring dynamics.

---

## 2) Data assumptions (raw datasets)

Each recording should provide:
- `X`: neural activity matrix of shape `[T, n]` where `n` is the number of recorded neurons.
- `neurons`: neuron ID strings aligned to columns of `X` (when available).
- `is_labeled`: logical vector aligned to columns of `X` (identity known vs unknown).
- `dt`: sampling interval.
- required metadata for this repo’s experiments:
  - **pre/post split** (`ranges_raw`) for Experiment 1.

Notes:
- “Recorded” vs “labeled” are different: `is_labeled` indicates neuron identity, not whether
  the neuron was recorded (all columns in `X` are recorded).

---

## 3) Mask handling (critical)

Whenever you compute distances, embeddings, or projections:
- use the recorded neurons for that worm (all columns in `X`).

If you compare across worms/datasets:
- use the **intersection of neuron IDs** (by name) when available:
  - `idx_common = intersect(neurons_wormA, neurons_wormB, ...)`
- if neuron identity alignment is weak or missing:
  - analyze **within‑worm** first, then compare scaffold‑level statistics across worms.

---

## 4) Preprocessing (minimal)

Default to **no additional normalization** unless you have a specific reason.
Many raw datasets already include resampling/smoothing/normalization steps.
Avoid double‑smoothing or double‑z‑scoring.

If needed:
- per‑neuron z‑score within each worm:
  - `Xz[:,i] = (X[:,i] - mean)/std` over time
- light temporal smoothing only if traces are very noisy (avoid double‑smoothing).

**Detrending policy (Experiment 1):**
- Default is **off**.
- If used, fit a **linear trend on the pre‑split only** and apply it to the full trace.
  This avoids leakage from the evaluation half while keeping train/eval consistent.

---

## 5) LOOPER pipeline (conceptual)

You want a mapping:
- raw activity snapshots \(x(t)\) → latent state \(y(t)\) (e.g., diffusion map)
- latent dynamics → scaffold state \((\alpha(t), \theta(t))\)

Typical steps:
1. (optional but often used) delay embedding of \(Xz\)
2. build a neighborhood graph / local transition model
3. compute diffusion map coordinates \(y(t)\)
4. run LOOPER to extract stable 1D trajectories / strands:
   - assign each time point a trajectory ID \(\alpha\)
   - assign a phase-like coordinate \(\theta\) along that trajectory

### 5.1) LOOPER outputs in this repo (what you actually need)

The MATLAB pipeline writes these fields (names used in the OSF scripts):
- `BestModel`: transition matrix between scaffold states (loop/phase bins)
- `BestEmission`: emission means per scaffold state
- `BestStateMap`: per‑timepoint mapping to (loop ID, phase bin)
- `BestLoopAssignments`: loop cluster membership mapping

Use these for scaffold alignment, phase velocity estimation, reconstruction, and distance proxy.

### 5.2) LOOPER parameters you will likely tune

Important knobs (see `LOOPER.md` / Table 2 of the paper for reference defaults):
- `NearestNeighbors`: local neighborhood size in diffusion map
- `RepopulateDensity`: matrix powering target for connectivity
- `MinReturnTime`: separation for repeated peaks in local neighbor search
- `DistanceMeasure` (for `reduceMatrix`): default `correlation`
- `MaxCheckTime`: time horizon used in KL evaluation
- `TotalStates`: pre‑loop cluster count
- `UseTerminalState`: whether to add terminal state for trial stitching

If you tune these, always log them alongside results.
We also set a wider `PutativeLoopCounts` range (`[8 7 6 5 4 3 2]`) in current
Atanas scripts to reduce the chance that the MDL minimum sits on a boundary
(per the paper’s note on loop/cluster count ranges).

### 5.3) Trial handling

If your data is trialized, LOOPER expects correct trial boundaries:
- `preprocessData.m` builds `processedTrialSwitches` and uses terminal state logic
  to avoid invalid cross‑trial transitions.
- If you don’t have trials, treat the full recording as one continuous stream.

---

## 6) Distance-to-scaffold \(d(t)\) (transverse deviation proxy)

Goal: a scalar measuring how far the current latent point is from the assigned strand.

Minimal implementation:
1. represent each strand as a set of prototype points in latent space:
   - e.g., bin by θ and store mean latent coordinate per bin
2. for time point t:
   - find the prototype point for the assigned \((\alpha(t), \theta(t))\)
   - define \(d(t) = \| y(t) - y_{\alpha,\theta} \|\)

Alternative (simpler):
- \(d(t)\) = distance to nearest point on the strand’s point cloud.

This is the key input to “split stability” analyses. If latent coordinates are not saved,
use a **data‑space proxy**:
- \(d(t) = \| X(t) - \mu_{\alpha,\theta} \|\) where \(\mu_{\alpha,\theta}\) is the
  emission mean for the assigned scaffold bin (`BestEmission`).
- In practice we use the **z‑scored version** of this distance (normalize by
  per‑state STD from `FinalStream`), matching OSF scripts and making \(\Delta d\)
  comparable across states.

### 6.1) Eval projection details (Experiment 1)

The eval scripts project **raw data** onto the learned scaffold with:
- **empirical state means** from `FinalStream` (OSF‑style), computed per
  `(loop ID, phase bin)` using `BestStateMap`/`BestLoopAssignments`.
- per‑state STD from `FinalStream` using the same matching; fallback to a
  global mean/STD if a state has no points.
- distance metric for both assignment and \(d(t)\): diagonal z‑score distance
  (sum of squared z‑scores, then \(d=\sqrt{\cdot}\)).
- a sanity warning is emitted if the processed length differs from
  `T_raw - DelayCount*DelayTime` by more than 1.

For split alignment, `compute_delta_d`:
- trims indices by `DelayCount*DelayTime`,
- clamps pre/post indices to `[1, procLen]`,
- and ensures the baseline window is strictly pre‑split.

### 6.2) Phase continuity diagnostics (loop‑ishness)

To check whether recovered trajectories are genuinely cyclic, diagnostics
summarize **θ continuity within constant‑α segments**:
- fraction of steps with \(|\Delta\theta| \le 1\) bin,
- variance of wrapped \(\Delta\theta\),
- cycles per minute (unwrapped θ change / bins),
- median \(d(t)\) within segments,
- mean segment length (frames).

These are reported in `run_looper_diagnostics` and summarized for
Atanas baseline runs in `summary.csv`.

### 6.3) Summary metrics (Experiment 1)

Each eval writes a `summary.csv` row with:
- split stability: `mean_pre`, `mean_post`, `post_slope`, `d_peak`
- loop stability: `unique_loops`, `loop_switches`, `segments`
- phase continuity: `phase_frac_small`, `phase_var`, `cycles_per_min`,
  `median_d`, `mean_segment_len`

---

## 7) Kato 2015 (positive control) specifics

Kato scripts in `kato_looper/` are used as a positive control:
- **Minimal single‑worm run** uses paper‑style params (Table 2) and
  per‑neuron detrend (across time), then LOOPER z‑scores.
- **Shared‑neuron run** concatenates five worms on the identified‑label intersection:
  per‑worm detrend + per‑worm z‑score **before** concatenation, then LOOPER
  z‑score is disabled to avoid double normalization. `TrialData` encodes worm ID.
- Shared run uses `MinimumReturnTime = 50` and `MaxCheckTime = 1` to limit memory
  use on long concatenated traces.

---

## 8) Outputs & reproducibility (current repo)

Always save:
- per‑run `*.mat` with `worm` + `saveData` (LOOPER outputs)
- per‑eval `summary.csv` (split stability + phase continuity metrics)
- diagnostic figures in `results/<run>/diagnostics/`
- optional: `final_stream_pca.fig/.png` from LOOPER

Notes:
- Scripts use fixed filenames (no timestamps) and overwrite old results.
- Keep `results/` in sync with the scripts you intend to cite in write‑ups.
