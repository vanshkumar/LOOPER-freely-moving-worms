# Methods

This document describes the processing pipeline, evaluation protocol, metric definitions, and dataset-specific parameter choices used in this repo.

For the experimental design and rationale, see [EXPERIMENTS.md](https://github.com/vanshkumar/LOOPER-freely-moving-worms/blob/main/EXPERIMENTS.md). For how to run the scripts, see [RUNBOOK.md](https://github.com/vanshkumar/LOOPER-freely-moving-worms/blob/main/RUNBOOK.md). For results, see [RESULTS.md](https://github.com/vanshkumar/LOOPER-freely-moving-worms/blob/main/RESULTS.md).

**Source of truth:** the MATLAB scripts are authoritative. If anything here diverges from the code, trust the scripts.

---

## Pipeline overview

Every experiment in this repo follows the same sequence:

```
Raw data (neurons × time)
  → Preprocessing (detrend, z-score)
  → LOOPER fitting (delay embed → diffusion map → state reduction → loop discovery → phase binning)
  → Evaluation (reconstruction correlation, phase continuity, stationarity drift)
  → Outputs (summary.csv, diagnostic figures, saved model)
```

The pipeline is implemented in two core scripts:
- **`looper_run_core.m`** handles data loading through LOOPER fitting and saves the model.
- **`looper_eval_core.m`** loads the saved model, runs diagnostics, computes metrics, and writes `summary.csv`.

Each entry-point script (e.g., `atanas_all_fidelity.m`) is a thin wrapper that calls these two functions with the appropriate dataset and mode.

---

## 1. Data loading

### Kato 2015

Loaded by `kato_2015/load_kato_data.m` from `KATO_WT_NoStim.mat`.

Each of the 5 worms provides:
- `RawData`: neurons × time matrix (bleach-corrected deltaF/F)
- `dt_sec`: ~0.344 s (derived from time vector, ~2.9 Hz)
- `neuron_ids`: string labels where identified
- T: 3021-3311 timepoints per worm (~17-19 min)

### Atanas 2023

Loaded by `atanas-data/load_atanas_data.m` from per-worm JSON files in `atanas-data/baseline/`.

Each of the 21 worms provides:
- `RawData`: neurons × time matrix (from `trace_array`, already z-scored)
- `dt_sec`: ~0.60 s (derived from `timestamp_confocal`)
- Behavioral covariates: velocity, angular_velocity, head/body curvature, pumping, reversal events (not used in current analysis)
- T: 1600 timepoints per worm (~16 min)

---

## 2. Preprocessing

### Detrending

**Kato:** Per-neuron linear detrend is applied. For each neuron, a linear fit (slope, intercept) is computed on the training portion of the data and subtracted from the full trace. This removes slow drift (e.g., photobleaching residuals). Implemented in `looper_helpers.detrend_fit` / `looper_helpers.detrend_apply`.

**Atanas:** Detrending is disabled. The Atanas traces are already z-scored in the source JSONs, and applying an additional linear detrend was not found to be necessary.

### Z-scoring

**Kato:** LOOPER's internal z-scoring is enabled (`PreprocessData.ZScore = true`). Each neuron is centered and scaled during LOOPER's preprocessing step.

**Atanas:** LOOPER's internal z-scoring is disabled (`PreprocessData.ZScore = false`) because the traces are already z-scored. Enabling it would double-normalize.

### Delay embedding

LOOPER delay-embeds the neuron × time matrix before computing the diffusion map. This converts a snapshot at time t into a short trajectory segment: `[x(t), x(t-Δ), x(t-2Δ), ..., x(t-kΔ)]` plus the corresponding time derivatives. This captures temporal structure that a single time-slice cannot.

Two parameters control the embedding:
- **DelayTime** (Δ): spacing between delays, in frames.
- **DelayCount** (k): number of delays.

The effective embedding window is `DelayTime × DelayCount` frames. The first `DelayTime × DelayCount` timepoints are trimmed (they have no complete delay history), so the embedded stream is shorter than the original.

**Kato parameters:** `DelayTime = 10`, `DelayCount = 10`. These are the paper's default values for this dataset (dt ~0.34 s, so the embedding window is ~34 s).

**Atanas parameters:** `DelayTime = 4`, `DelayCount = 5`. These were chosen by running `atanas_autocorr_delay.m`, which estimates the autocorrelation decay time following the rule from Brennan & Proekt (eLife 2019), then scales to the slower sampling rate (~0.60 s, so the embedding window is ~12 s).

---

## 3. LOOPER pipeline

LOOPER is run via `LOOPER.m`, a headless wrapper around the original GUI application's pipeline. The steps are:

### 3a. Diffusion map construction

`buildDiffusionMap.m` computes a data-driven transition operator on the embedded neural trajectory:

1. For each time point, find the `NearestNeighbors` (= 8) nearest neighbors in the embedded space.
2. Compute local distances that combine static similarity (activity vectors) and velocity similarity (smoothed derivative cosine distance), via `findBestAxes.m`.
3. Build an asymmetric Markov transition matrix: entry (i, j) is the probability of transitioning from state i to state j. The asymmetry is key — it captures the directionality of neural dynamics.
4. Raise the matrix to successive powers until `RepopulateDensity` (= 0.95) connectivity is achieved. This smooths the transition structure.
5. Exclude self-transitions within `MinimumReturnTime` frames (prevents trivially short loops).

The diffusion map is the most expensive step and is cached to `results/<run>/cache/` for re-use. If you change preprocessing or parameters, clear the cache to avoid stale results (see [RUNBOOK.md](https://github.com/vanshkumar/LOOPER-freely-moving-worms/blob/main/RUNBOOK.md)).

### 3b. State reduction

`reduceMatrix.m` coarse-grains the transition matrix:

1. Hierarchical clustering on the rows of the transition matrix (using correlation distance).
2. For each candidate number of clusters, build a reduced transition matrix and measure how well it approximates the full dynamics (KL divergence over 1 to `MaxCheckTime` time steps).
3. Select the cluster count that minimizes divergence.

The result is a set of ~10-25 coarse states, each representing a cluster of similar time points.

### 3c. Loop discovery

`buildMinimalModelFromMatrix.m` identifies loop-like trajectories:

1. Build transition-from and transition-to matrices from the coarse state time series.
2. For each state, find the shortest loop (closed path) through the transition graph.
3. For each candidate loop count (from `PutativeLoopCounts = [8,7,6,5,4,3,2]`), assign states to loops, discretize phase, build a Markov model on (loop, phase) states, and compute log-likelihood.
4. Select the loop count that maximizes likelihood.

### 3d. Phase binning → scaffold

`buildMinimalModelFromLoops.m` aligns loops to a common phase reference, bins each loop into `TotalStates` (= 25) phase bins, and assigns each time point to a (loop ID α, phase bin θ) state.

The scaffold is LOOPER's core output: a map from each time point to a discrete (α, θ) state, plus an emission model that gives the expected neural activity at each state.

---

## 4. Evaluation protocol

### 4a. Fidelity test

Train LOOPER on the full trace. The scaffold and emission model come directly from the training fit. Evaluate:

- **Reconstruction quality:** How well does the scaffold's emission model reconstruct the original embedded neural stream?
- **Phase continuity:** Is phase progression smooth within each loop?
- **Loop structure:** How many loops, how stable are the assignments?

### 4b. Stationarity test (split-half)

Train LOOPER on the **first half** of the trace (`splitIdx = floor(T/2)`). The scaffold is learned only from the pre-split data. Then:

1. **Project the full trace** onto the learned scaffold. For each time point in the full (embedded) trace, find the nearest scaffold state by z-scored distance to the emission means. This assigns every time point — including the held-out second half — to a (α, θ) state with a distance d.
2. **Compute distance-to-scaffold over time.** The distance d(t) measures how far each time point is from the nearest scaffold state. During the training half, d(t) should be low (the scaffold was fit to this data). During the held-out half, d(t) may increase if the scaffold doesn't generalize.
3. **Compute drift metrics.** Baseline-subtract d(t) using the pre-split median, then measure mean, slope, and peak of the post-split deviation.
4. **Compute held-out reconstruction.** Correlate the scaffold reconstruction with the actual neural stream in the post-split half only.

The projection logic is in `looper_eval_core.m` and uses empirical state means and standard deviations from the training-half `FinalStream`.

---

## 5. Metric definitions

Each evaluation writes a row to `summary.csv` with the following columns. All metrics are computed identically for both datasets.

### Reconstruction quality

| Metric | Definition | Interpretation |
|--------|-----------|----------------|
| `recon_corr_full` | Pearson correlation between the original embedded stream and the scaffold's reconstruction (full trace). Computed by `plotReconstruction.m`. | Higher = better in-sample fit. Note: this is Pearson r, not R-squared. Not directly comparable to the paper's R² values. |
| `recon_corr_post` | Same correlation, computed only on the held-out post-split half. Stationarity runs only. | Positive = scaffold still fits. Near zero or negative = scaffold does not generalize. |

### Loop structure

| Metric | Definition | Interpretation |
|--------|-----------|----------------|
| `unique_loops` | Number of distinct loop IDs (α) assigned by LOOPER. | Typically 2-6. More loops = more complex dynamics. |
| `loop_switches` | Number of times the loop assignment changes between consecutive time points. | High relative to trace length = unstable assignment. |
| `segments` | Number of contiguous runs of the same loop ID. | = loop_switches + 1. |

### Phase continuity

| Metric | Definition | Interpretation |
|--------|-----------|----------------|
| `phase_frac_small` | Fraction of within-loop phase steps where \|Δθ\| ≤ 1 bin (after wrapping). | High (>0.95) = smooth phase progression. Low = jumpy. |
| `phase_var` | Variance of the wrapped per-step Δθ values. | Lower = smoother. Kato median ~0.05, Atanas median ~0.29. |
| `phase_speed_bins_per_min` | Mean rate of phase progression through the scaffold. | Higher = faster cycling. |
| `median_d` | Median z-scored distance from each time point to its assigned scaffold state. | Lower = tighter fit to the scaffold. |
| `mean_segment_len` | Average length (in time points) of contiguous same-loop runs. | Longer = more stable loop assignment. |

### Stationarity drift (stationarity runs only)

| Metric | Definition | Interpretation |
|--------|-----------|----------------|
| `mean_pre` | Mean baseline-subtracted distance-to-scaffold in the pre-split half. | Should be near 0 (scaffold was fit here). |
| `mean_post` | Mean baseline-subtracted distance-to-scaffold in the post-split half. | Large positive = scaffold doesn't fit held-out data. |
| `post_slope` | Linear slope of distance-to-scaffold over time in the post-split half. | Positive slope = progressive drift away from scaffold. |
| `d_peak` | Peak distance-to-scaffold in the post-split half. | Maximum deviation from scaffold. |

### LOOPER internal diagnostics

| Metric | Definition | Interpretation |
|--------|-----------|----------------|
| `validation_score_mean` | Mean of LOOPER's internal segment-wise validation score (from `validateModel.m`). | Higher = better internal consistency. |
| `validation_score_std` | Standard deviation of the validation score. | Lower = more consistent across segments. |

Non-applicable fields (e.g., stationarity metrics in fidelity runs) are recorded as `NaN`.

---

## 6. Dataset-specific parameter summary

| Parameter | Kato | Atanas | Reason for difference |
|-----------|------|--------|----------------------|
| `ZScore` | true | false | Atanas already z-scored in source JSONs |
| `DelayTime` | 10 | 4 | Scaled to sampling rate (~0.34 s vs ~0.60 s) |
| `DelayCount` | 10 | 5 | Autocorrelation-based rule from paper |
| `MinimumReturnTime` | 10 | 4 | Scaled to sampling rate |
| `MaxCheckTime` | 10 | 4 | Scaled to sampling rate |
| `Detrending` | per-neuron linear | disabled | Atanas pre-normalized |
| `NearestNeighbors` | 8 | 8 | Same |
| `RepopulateDensity` | 0.95 | 0.95 | Same |
| `TotalStates` | 25 | 25 | Same |
| `PutativeLoopCounts` | [8 7 6 5 4 3 2] | [8 7 6 5 4 3 2] | Same |

Both datasets use identical evaluation protocols and metrics. The only differences are in preprocessing and LOOPER parameters, driven by the different sampling rates and normalization conventions.

---

## 7. Kato shared-neuron concatenated run

The shared-neuron run is a special case that concatenates all 5 Kato worms using only the neurons identified in all animals. This is the closest approximation to the paper's Fig-5B analysis.

Procedure (`kato_shared_run.m`):
1. Find the intersection of neuron labels across all 5 worms (letter-based matching). Result: 22 shared neurons.
2. Per-worm: detrend each neuron, then z-score.
3. Concatenate the 5 worms into one long trace (15,662 timepoints, ~90 min).
4. Set `TrialData` to encode worm boundaries (so LOOPER knows where one worm ends and the next begins).
5. Disable LOOPER's internal z-scoring (to avoid double-normalizing).
6. Set `MinimumReturnTime = 50` (matching the paper's OSF checkpoint) and `MaxCheckTime = 1` (to prevent memory issues on the long concatenated trace).

This run is fidelity-only (no split-half) because the concatenation creates a multi-worm trace where a temporal split would not be meaningful.

---

## 8. Diagnostic figures

Each experiment run produces diagnostic plots. Here is what they show and where to find representative examples.

### Per-worm diagnostics (in `diagnostics/` subdirectories)

| Plot | Filename pattern | What it shows |
|------|-----------------|---------------|
| Loop structure | `*_plotLoops.png` | 3D rendering of LOOPER's loops as tubes in model space. Color distinguishes loops. |
| Reconstruction | `*_reconstruction_full.png` | Top: original neural activity (neurons × time). Bottom: scaffold reconstruction. Title includes recon_corr value. |
| Phase over time | `*_loop_phase_time.png` | Phase (θ) plotted against time. Smooth ramps = clean cycling. Jumps = phase discontinuities. |
| Validation | `*_validateModel.png` | Trajectory in first 3 PCs, color-coded by state. Shows dynamics, emissions, and transitions. |
| Transition matrix | `*_transition_matrix.png` | State-to-state transition probabilities. Band structure along the diagonal indicates sequential cycling. |

### Stationarity-specific plots (in `plots/` subdirectories)

| Plot | Filename pattern | What it shows |
|------|-----------------|---------------|
| Stationarity recovery | `*_stationarity_recovery.png` | Distance-to-scaffold (Δd) over time, centered on the split point. Flat pre-split + rising post-split = stationarity failure. |
| Loop assignment | `*_loop_assignment.png` | Loop ID (α) over time for the full projected trace. Changes in post-split pattern indicate scaffold mismatch. |

### PCA stream plots (in run root directories)

| Plot | Filename pattern | What it shows |
|------|-----------------|---------------|
| Final dynamics stream | `*_final_stream_pca.png` | The delay-embedded neural trajectory projected onto PC1 vs PC2. Dense looping = cyclic dynamics. Sparse/diffuse = less structured. |

### Key comparisons

- **Kato shared-neuron loops** (`results/kato_shared/diagnostics/kato_shared_plotLoops.png`): Two well-separated loops (blue and orange tubes), consistent with the paper's two-loop structure.
- **Kato stationarity drift** (`results/kato_all/stationarity/plots/worm_1_stationarity_recovery.png`): Flat pre-split, roughly linear rise to Δd ~ 190 post-split.
- **Atanas stationarity drift** (`results/atanas_single/stationarity/worm_2022-06-14-01_stationarity_recovery.png`): Same pattern but more variable, Δd rises gradually to ~50.
- **OSF reference** (`results/osf_wormlooper2_diagnostics/wormlooper2_plotLoops.png`): Paper's own output on the Kato data, showing similar two-loop structure.

---

## 9. Known issues and caveats

- **Reconstruction correlation vs. paper's R².** Our `recon_corr_full` is Pearson r; the LOOPER paper reports R². These are not directly comparable. Our shared-neuron Kato run gives r = 0.661; the paper reports R² ≈ 0.79 on the same dataset. The gap likely reflects differences in neuron set (22 vs. 15 curated), detrending method, and `MinimumReturnTime` parameter.

- **No null model.** We have not run LOOPER on time-shuffled data, and neither does the original paper. A null model would calibrate what `recon_corr` LOOPER produces by chance on data with the same spectral properties but no temporal structure. This is a limitation shared with the published C. elegans results.

- **Behavioral covariates unused.** The Atanas dataset includes velocity, angular velocity, curvature, and other behavioral time series. These have not been correlated with LOOPER phase or used to condition the analysis. Behavior-conditioned LOOPER is a natural next step.

- **Detrending method differs from paper.** This repo uses per-neuron linear detrend across time. The paper's OSF script (`buildWormData.m`) uses column-wise detrending. The effect on results has not been systematically compared.

- **Embedding space incomparability.** Distance-to-scaffold values between Kato and Atanas are not directly comparable because the embedding spaces differ (different neuron counts, delay parameters, preprocessing). Reconstruction correlations are more comparable but still computed in different spaces.
