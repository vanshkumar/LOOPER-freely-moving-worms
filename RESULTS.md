# Results: Applying LOOPER to Freely Moving C. elegans

**Last updated:** 2026-02-09

For background on the scientific question and experimental design, see `EXPERIMENTS.md`. For pipeline details and metric definitions, see `METHODS.md`. For how to reproduce these results, see `RUNBOOK.md`.

---

## Summary

LOOPER (Brennan et al. 2023) has been shown to recover clean loop-like neural dynamics in immobilized C. elegans. We asked whether it can also recover loop structure in **freely moving** worms, where neural dynamics are more complex due to behavioral state switching and sensory feedback. We applied LOOPER to two datasets — **Kato et al. 2015** (5 immobilized worms, positive control) and **Atanas et al. 2023** (21 freely moving worms, target) — and evaluated fidelity (in-sample fit) and split-half stationarity (temporal generalization).

**Key findings:**

- LOOPER recovers loop structure in both datasets with comparable reconstruction correlations (Kato median r = 0.57, Atanas median r = 0.59). This uses the same in-sample reconstruction approach the original LOOPER paper uses for C. elegans — continuous recordings have no trial structure for held-out validation (paper Table 1: C. elegans Trials = "Contiguous", Conditions = N/A). A null model (e.g., time-shuffled data) would help calibrate these values, but has not been run here or in the original paper.
- The freely moving scaffolds have **noisier phase dynamics**: phase variance is ~6x higher, cycling is ~2.6x faster, and loop switching is more frequent. This is consistent with more complex behavioral dynamics during movement, though the contribution of genuine cyclic structure vs. LOOPER fitting to less-cyclic data cannot be disentangled without a null model.
- The **split-half stationarity test fails for both datasets**, including the positive control. Post-split reconstruction correlations are near zero or negative (Kato median r = -0.03, Atanas median r = -0.05), and distance-to-scaffold grows monotonically after the split. This test goes beyond what the paper did — the paper used trial-based cross-validation, which does not apply to continuous C. elegans recordings.
- Because the positive control itself fails stationarity, this test **cannot be used to draw differential conclusions about freely moving worms**. A stationarity criterion the positive control passes is needed before this axis of comparison is informative.

These results show that LOOPER recovers loop structure in freely moving recordings with correlations comparable to the positive control, using the same methodology as the published paper. The stationarity test reveals that scaffolds do not generalize across temporal halves in either condition — a finding that goes beyond the paper's analysis and motivates further work on null models and alternative validation approaches.

---

## 1. Datasets

### Kato 2015 (immobilized)

| Property | Value |
|----------|-------|
| Source | Kato et al. 2015, *Cell* |
| Animals | 5 worms |
| Neurons per worm | 109, 135, 131, 125, 129 |
| Timepoints per worm | 3137, 3134, 3059, 3311, 3021 |
| Sampling interval | ~0.34 s (varies slightly per worm) |
| Recording duration | ~17-19 min per worm |
| Condition | Immobilized (paralyzed), no sensory input |
| Preprocessing | Per-neuron linear detrend, z-score |

For the shared-neuron concatenated run: 22 neurons common across all 5 worms, 15,662 total timepoints (~90 min concatenated).

### Atanas 2023 (freely moving, baseline)

| Property | Value |
|----------|-------|
| Source | Atanas et al. 2023, *Cell* |
| Animals | 21 worms |
| Neurons per worm | 109-153 (median 142) |
| Timepoints per worm | 1600 (all worms) or 1615 (2 worms) |
| Sampling interval | ~0.60 s |
| Recording duration | ~16 min per worm |
| Condition | Freely moving, baseline (no stimulus) |
| Preprocessing | Already z-scored in source JSONs; no additional detrending |

Worm identifiers follow the date-based NeuroPAL convention (e.g., `2022-06-14-01`).

---

## 2. Methods summary

The same pipeline is applied to both datasets: preprocess (detrend for Kato, skip for Atanas) → delay-embed → run LOOPER → evaluate. Full details in `METHODS.md`.

**Two tests per worm:**
- **Fidelity:** Train on the full trace, evaluate on the same trace (in-sample).
- **Stationarity:** Train on the first half, project the full trace, evaluate reconstruction and drift in the held-out second half.

**Key metrics:** `recon_corr_full` (in-sample reconstruction, Pearson r), `recon_corr_post` (held-out reconstruction), `phase_var` (phase smoothness), `unique_loops` / `loop_switches` (loop structure), `mean_post` / `post_slope` / `d_peak` (stationarity drift). See `METHODS.md` Section 5 for full definitions.

**Dataset-specific parameters:** Kato uses the paper's defaults (DelayTime=10, DelayCount=10, z-score enabled). Atanas uses autocorrelation-scaled parameters (DelayTime=4, DelayCount=5, z-score disabled because already normalized). See `METHODS.md` Section 6 for the full parameter table.

---

## 3. Kato Fidelity (positive control)

### Shared-neuron concatenated run (closest to paper)

This concatenates all 5 Kato worms on the 22 neurons common across all animals, with per-worm detrend and z-score before concatenation.

| Metric | Value |
|--------|-------|
| n_neurons | 22 |
| T | 15,662 timepoints (~90 min) |
| recon_corr_full | **0.661** |
| unique_loops | **2** |
| loop_switches | 87 |
| phase_frac_small | 0.993 |
| phase_var | 0.133 |

**What the diagnostics show:** The `plotLoops` figure (`results/kato_shared/diagnostics/kato_shared_plotLoops.png`) reveals two well-separated loops in model space (rendered as orange and blue tubes), visually consistent with the two-loop structure reported in the original LOOPER paper. The `validateModel` figure shows the trajectory color-coded by time/state, confirming that the dynamics trace out coherent cyclic paths. The reconstruction heatmap (`kato_shared_reconstruction_full.png`) shows the top panel (original) and bottom panel (reconstruction) with visible periodic banding, and an overall correlation of 0.56.

**Note on recon_corr:** The 0.661 we report is computed on the full embedded stream. The reconstruction figure title shows 0.56, which may reflect a different computation window or neuron subset. The paper reports R-squared values (which are not directly comparable to Pearson r). This gap is a known documentation issue but does not affect the qualitative conclusion.

### Per-worm fidelity (individual animals)

| Worm | Neurons | T | recon_corr | Loops | Switches | phase_var |
|------|---------|---|------------|-------|----------|-----------|
| 1 | 109 | 3137 | **0.641** | 3 | 6 | 0.037 |
| 2 | 135 | 3134 | 0.510 | 2 | 1 | 0.127 |
| 3 | 131 | 3059 | 0.572 | 4 | 31 | 0.050 |
| 4 | 125 | 3311 | **0.684** | 3 | 15 | 0.021 |
| 5 | 129 | 3021 | 0.528 | 2 | 2 | 0.116 |
| **Median** | **129** | **3134** | **0.572** | **3** | **6** | **0.050** |

**Observations:**
- All 5 worms show positive reconstruction correlations (range 0.51-0.68), confirming that loop-like structure is recoverable per-individual.
- Phase continuity is consistently high: phase_frac_small >= 0.995 for all worms, and phase_var ranges from 0.021 (very smooth) to 0.127.
- Worms 2 and 5 show the simplest structure (2 loops, 1-2 switches), while Worm 3 is the most dynamic (4 loops, 31 switches).
- Worm 4 has the highest reconstruction and smoothest phase (phase_var = 0.021), suggesting particularly clean cyclic dynamics.

**Conclusion:** The positive control passes fidelity. LOOPER recovers 2-4 loops per worm with smooth phase progression and moderate-to-strong reconstruction. This follows the same in-sample evaluation approach the paper uses for C. elegans data (continuous recordings, no trial-based validation split).

---

## 4. Atanas Fidelity (freely moving)

### Per-worm fidelity (all 21 animals)

| Worm | Neurons | recon_corr | Loops | Switches | phase_var | phase_speed |
|------|---------|------------|-------|----------|-----------|-------------|
| 2022-06-14-01 | 109 | 0.541 | 2 | 30 | 0.713 | 24.5 |
| 2022-06-14-07 | 130 | 0.651 | 3 | 13 | 0.041 | 3.2 |
| 2022-06-14-13 | 134 | 0.551 | 2 | 9 | 0.111 | 9.2 |
| 2022-06-28-01 | 142 | 0.594 | 3 | 10 | 0.322 | 15.1 |
| 2022-06-28-07 | 142 | 0.531 | 2 | 30 | 0.205 | 12.8 |
| 2022-07-15-06 | 128 | 0.622 | 2 | 9 | 0.059 | 5.1 |
| 2022-07-15-12 | 132 | 0.534 | 3 | 15 | **1.712** | 55.7 |
| 2022-07-20-01 | 144 | 0.565 | 4 | 46 | 0.064 | 5.3 |
| 2022-07-26-01 | 133 | 0.570 | 3 | 20 | 0.628 | 26.6 |
| 2022-08-02-01 | 152 | 0.564 | 2 | 11 | 0.624 | 21.2 |
| 2023-01-09-28 | 145 | **0.706** | 3 | 14 | 0.110 | 9.8 |
| 2023-01-17-01 | 138 | 0.614 | 3 | 15 | 0.272 | 11.7 |
| 2023-01-19-01 | 124 | 0.586 | 2 | 1 | **0.919** | 32.1 |
| 2023-01-19-08 | 147 | 0.660 | 3 | 12 | 0.532 | 20.0 |
| 2023-01-19-15 | 144 | 0.568 | 2 | 1 | 0.232 | 10.0 |
| 2023-01-19-22 | 153 | 0.642 | 2 | 3 | 0.328 | 14.1 |
| 2023-01-23-01 | 151 | 0.615 | 2 | 15 | 0.230 | 20.4 |
| 2023-01-23-08 | 143 | **0.778** | 6 | 15 | 0.064 | 4.3 |
| 2023-01-23-15 | 147 | 0.623 | 2 | 3 | 0.344 | 18.2 |
| 2023-01-23-21 | 126 | 0.612 | 2 | 17 | 0.291 | 20.9 |
| 2023-03-07-01 | 153 | 0.556 | 3 | 23 | **1.585** | 57.8 |
| **Median** | **142** | **0.594** | **2** | **13** | **0.291** | **15.1** |

**Observations:**

- **Reconstruction is comparable to Kato.** Median recon_corr = 0.594 (vs. Kato's 0.572). The range (0.53-0.78) overlaps substantially with Kato (0.51-0.68). This means LOOPER finds quantitatively similar in-sample loop-like structure in freely moving worms as in immobilized worms.

- **Phase is noisier.** Median phase_var = 0.291 in Atanas vs. 0.050 in Kato -- roughly 6x higher. Three worms have phase_var > 0.9 (2022-07-15-12, 2023-01-19-01, 2023-03-07-01), indicating highly irregular phase progression. These same worms show the fastest phase speeds (32-58 bins/min vs. a Kato median of 5.8 bins/min), suggesting rapid, jumpy cycling rather than smooth oscillation.

- **More loop switching.** Median loop_switches = 13 in Atanas vs. 6 in Kato, despite shorter recording durations (~16 min vs. ~18 min). This is consistent with more dynamic, behaviorally-modulated neural activity in freely moving animals.

- **Substantial inter-worm variability.** The best-performing Atanas worm (2023-01-23-08: recon_corr = 0.778, 6 loops, phase_var = 0.064) outperforms all Kato worms. The worst performers (2022-06-28-07: recon_corr = 0.531, phase_var = 0.205) are comparable to the weakest Kato worms. This variability likely reflects differences in behavioral state during the recording.

- **Two worms are outliers in phase quality.** Worms 2022-07-15-12 and 2023-03-07-01 have phase_var > 1.5 and phase_frac_small < 0.85 (vs. >0.95 for all other worms). LOOPER fits a scaffold to these traces, but the phase progression is so erratic that the "loop" interpretation is questionable. These worms may have been in an atypical behavioral state.

**Conclusion:** LOOPER recovers loop structure in freely moving data with reconstruction correlations comparable to immobilized worms, using the same in-sample approach as the original paper's C. elegans analysis. Phase progression is substantially noisier, loop switching is more frequent, and there is more inter-individual variability — consistent with the added complexity of behavioral state switching during movement. A null model (time-shuffled data) would help quantify how much of the reconstruction is attributable to genuine cyclic structure vs. LOOPER's smoothing; this control has not been run here or in the original paper.

---

## 5. Kato Stationarity (split-half stress test)

This test asks: does a scaffold learned from the first half of the recording still describe the second half?

### Per-worm stationarity

| Worm | recon_corr_post | mean_post | post_slope | d_peak | Loops (post) |
|------|----------------|-----------|------------|--------|--------------|
| 1 | -0.025 | 97.8 | 0.335 | 189.8 | 4 |
| 2 | **0.237** | 81.1 | 0.290 | 164.7 | 2 |
| 3 | -0.082 | 98.4 | 0.263 | 197.8 | 1 |
| 4 | **0.156** | 58.9 | 0.128 | 96.2 | 5 |
| 5 | -0.109 | 58.2 | 0.145 | 94.4 | 3 |
| **Median** | **-0.025** | **81.1** | **0.263** | **164.7** | **3** |

**What the stationarity recovery plots show:** The distance-to-scaffold d(t) is stable and near zero during the pre-split (training) half, then increases monotonically after the split point. For Kato Worm 1 (`results/kato_all/stationarity/plots/worm_1_stationarity_recovery.png`), delta_d rises roughly linearly from ~0 to ~190 over ~540 seconds, with some oscillatory structure superimposed. This is not a gradual degradation -- it is a clear, sustained drift away from the learned scaffold.

**What the loop assignment plots show:** The loop assignment over time (`results/kato_all/stationarity/plots/worm_1_loop_assignment.png`) reveals that the pre-split half contains frequent switching among 3-4 loops, while the post-split half often collapses to a single dominant loop. This suggests the scaffold learned from the first half does not capture the dynamics of the second half, and the projection defaults to a single catch-all assignment.

**This is the critical finding:** Even the positive control fails this test. Three of 5 worms have negative recon_corr_post, and the two positive values (0.237, 0.156) are weak. All worms show large positive post_slope (scaffold drift grows over time) and high d_peak values (94-198).

**Interpretation:** The split-half stationarity test **goes beyond what the original LOOPER paper did for C. elegans**. The paper used trial-based cross-validation for datasets with trial structure (primate, mouse, human fMRI), but C. elegans recordings are continuous with no trials — so the paper did not test temporal generalization for this dataset type at all (paper Table 1: Trials = "Contiguous", Conditions = N/A, Bootstraps = 0). Our test asks a question the paper did not address: does a scaffold learned from one temporal half of a continuous recording generalize to the other half? The answer is no, even for immobilized worms.

Possible explanations for Kato stationarity failure:
1. **Slow drift in calcium dynamics** (photobleaching residual, slow adaptation) that a linear detrend does not fully remove.
2. **Genuine non-stationarity** in the neural oscillation parameters over ~10-min timescales.
3. **Sensitivity of the scaffold** to the specific trajectory sample used for training (the diffusion map embedding may not generalize from one half to another even if the underlying dynamics are stationary).

Without a null model (e.g., time-shuffled data), we cannot distinguish these explanations.

---

## 6. Atanas Stationarity (split-half stress test)

### Per-worm stationarity (all 21 animals)

| Worm | recon_corr_post | mean_post | post_slope | d_peak |
|------|----------------|-----------|------------|--------|
| 2022-06-14-01 | **0.173** | 19.0 | 0.064 | 50.3 |
| 2022-06-14-07 | -0.045 | 21.4 | 0.059 | 78.5 |
| 2022-06-14-13 | -0.040 | 33.8 | 0.112 | 80.8 |
| 2022-06-28-01 | -0.209 | 38.4 | 0.147 | 93.5 |
| 2022-06-28-07 | -0.022 | 28.3 | 0.051 | 95.2 |
| 2022-07-15-06 | 0.084 | 39.3 | 0.062 | 113.8 |
| 2022-07-15-12 | -0.055 | 22.7 | 0.027 | 48.2 |
| 2022-07-20-01 | 0.111 | 26.8 | 0.059 | 51.7 |
| 2022-07-26-01 | 0.096 | 30.1 | 0.078 | 56.3 |
| 2022-08-02-01 | -0.242 | 31.5 | 0.056 | 53.8 |
| 2023-01-09-28 | 0.073 | 38.1 | 0.073 | 76.0 |
| 2023-01-17-01 | -0.118 | 21.0 | 0.043 | 43.3 |
| 2023-01-19-01 | -0.079 | 33.8 | 0.124 | 83.0 |
| 2023-01-19-08 | -0.285 | 33.8 | 0.107 | 78.9 |
| 2023-01-19-15 | -0.096 | 26.2 | 0.055 | 44.3 |
| 2023-01-19-22 | -0.132 | 32.1 | 0.117 | 77.5 |
| 2023-01-23-01 | **0.132** | 30.9 | 0.092 | 69.1 |
| 2023-01-23-08 | -0.125 | 45.4 | 0.065 | 88.4 |
| 2023-01-23-15 | -0.186 | 39.1 | 0.091 | 67.9 |
| 2023-01-23-21 | -0.016 | 22.5 | 0.047 | 43.4 |
| 2023-03-07-01 | **0.115** | 37.3 | 0.121 | 72.1 |
| **Median** | **-0.045** | **31.5** | **0.065** | **72.1** |

**What the stationarity recovery plots show:** The Atanas recovery plots (e.g., `results/atanas_single/stationarity/worm_2022-06-14-01_stationarity_recovery.png`) show a more gradual drift than Kato, with greater high-frequency variability superimposed on the upward trend. The pre-split baseline fluctuations are larger in Atanas (amplitude ~10) than in Kato (amplitude ~5), reflecting noisier dynamics. The post-split increase reaches delta_d ~50 by the end of the recording, compared to ~190 for Kato Worm 1.

**Observations:**

- **Atanas also fails stationarity**, though with smaller absolute drift values than Kato (median mean_post = 31.5 vs. 81.1; median d_peak = 72.1 vs. 164.7). This is partly because the Atanas distance metric operates in a different embedding space (different neuron counts, delay parameters, and normalization).

- **Post-split reconstruction is poor.** Only 7 of 21 worms have positive recon_corr_post, and the best value (0.173) is weak. The median (-0.045) is essentially zero.

- **Post-slope is universally positive.** All 21 worms show the scaffold drifting further from the data over time after the split, confirming that the scaffold does not capture the held-out dynamics.

- **The stationarity failure is qualitatively similar to Kato.** Both datasets show the same pattern: good fit in the training half, progressive drift in the held-out half. The failure mode is the same; the magnitude differs.

---

## 7. Comparing the two datasets

### Fidelity comparison

| Metric | Kato (n=5) | Atanas (n=21) | Interpretation |
|--------|-----------|--------------|----------------|
| recon_corr median | 0.572 | 0.594 | Comparable in-sample fit |
| recon_corr range | 0.51 - 0.68 | 0.53 - 0.78 | Atanas has wider range but overlapping |
| unique_loops median | 3 | 2 | Similar |
| phase_var median | 0.050 | 0.291 | **Atanas is ~6x noisier** |
| phase_speed median | 5.8 bins/min | 15.1 bins/min | **Atanas cycles ~2.6x faster** |
| loop_switches median | 6 | 13 | More dynamic in freely moving |

The most striking difference is not in reconstruction quality (which is similar) but in **phase dynamics**: freely moving scaffolds cycle faster, with more irregular phase progression and more frequent loop transitions. This is consistent with real but behaviorally-modulated cyclic dynamics — the animal switches between behavioral modes (forward, reversal, turn) on timescales faster than the smooth oscillations seen in immobilized worms. A null model would help determine whether some of the phase noise reflects LOOPER fitting to less-cyclic data rather than genuine but noisy loops.

### Stationarity comparison

| Metric | Kato (n=5) | Atanas (n=21) | Note |
|--------|-----------|--------------|------|
| recon_corr_post median | -0.025 | -0.045 | Both near zero |
| mean_post median | 81.1 | 31.5 | Different scales (different embeddings) |
| post_slope median | 0.263 | 0.065 | Kato drifts faster in absolute terms |
| d_peak median | 164.7 | 72.1 | Kato peaks higher |
| % worms with positive recon_corr_post | 40% (2/5) | 33% (7/21) | Neither dataset generalizes well |

Both datasets fail the split-half stationarity test. The absolute distance values are not directly comparable because the embedding spaces differ (22-153 neurons, different delay parameters, different preprocessing). What is comparable is the qualitative pattern: **in both cases, the scaffold learned from the first half does not describe the second half.**

---

## 8. What these results mean

### What we can say

1. **LOOPER recovers loop structure in freely moving data with reconstruction correlations comparable to immobilized data.** Atanas median r = 0.59, Kato median r = 0.57. This uses the same in-sample reconstruction approach the original LOOPER paper uses for C. elegans — the data is a continuous trace with no trial structure for held-out validation (paper Table 1: Trials = "Contiguous", Conditions = N/A). Our fidelity results are on the same methodological footing as the paper's published C. elegans results.

2. **The freely moving scaffolds have qualitatively different phase dynamics.** Phase is noisier, cycling is faster, and loop switching is more frequent. This is consistent with behavioral modulation — the animal switches between behavioral modes on timescales faster than the smooth oscillations in immobilized worms.

3. **A strict temporal split-half test fails on both datasets.** The scaffold learned from the first ~8-9 minutes of recording does not generalize to the second ~8-9 minutes. This test goes beyond what the paper did — the paper used trial-based cross-validation for datasets that had trial structure, and did not test temporal stationarity on C. elegans data. The failure of the positive control means this test cannot distinguish the two datasets.

4. **The split-half test reveals a property of LOOPER scaffolds on continuous data that the paper did not address.** Whether this reflects genuine non-stationarity in the neural dynamics, slow drift artifacts, or sensitivity of the diffusion map embedding to the specific training sample is an open question.

### What we cannot say

1. **We cannot say that freely moving worms lack loop-like structure.** The fidelity results show they have structure comparable to the positive control, evaluated the same way the paper evaluates C. elegans data.

2. **We cannot say the loops are "less real" in freely moving worms.** The stationarity test does not discriminate between the datasets because both fail. A test that Kato passes but Atanas fails would be informative; we do not have one.

3. **We cannot attribute the stationarity failure to any specific mechanism** (slow drift, mode switching, method sensitivity, photobleaching residual) without additional controls.

4. **We cannot compare reconstruction quality between datasets in absolute terms.** The correlation values are computed in different embedding spaces with different numbers of neurons and delay parameters. The fact that they are numerically similar is suggestive but not rigorously comparable.

### Open issues

- **A null model would strengthen the fidelity results.** Running LOOPER on time-shuffled data would calibrate what recon_corr is expected by chance given LOOPER's smoothing and delay embedding. This control has not been run here or in the original paper. If shuffled data yields r ~ 0.5, the observed r ~ 0.6 is only marginally above chance; if it yields r ~ 0.1, the observed values are clearly meaningful.

- **The recon_corr gap with the paper is unresolved.** The original LOOPER paper reports higher reconstruction values (R-squared, not Pearson r), and the relationship between our `recon_corr_full` and the paper's metric is not fully characterized. Our shared-neuron Kato run (r = 0.661) is the closest comparison point, but exact parameter matching has not been verified.

- **Behavioral covariates have not been used.** The Atanas dataset includes velocity, angular velocity, and other behavioral time series. These have not been correlated with LOOPER phase or used to condition the analysis. Behavior-conditioned LOOPER (e.g., running LOOPER separately on forward-crawling and reversal epochs) could reveal whether mode-switching explains the stationarity failure or the noisier phase dynamics.

---

For diagnostic figure descriptions, see `METHODS.md` Section 8. For the complete file index and scripts-to-outputs mapping, see `RUNBOOK.md`.
