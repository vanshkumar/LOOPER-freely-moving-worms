# RESULTS.md

Negative and positive results from the LOOPER experiments in this repo.  
Updated: 2026-01-27

---

## Scope

We evaluated whether LOOPER recovers stable, loop-like dynamics in:

- **Freely moving NeuroPAL worms** (Atanas baseline recordings; train on first half, project full trace).
- **Immobilized Kato 2015 worms** (positive control; expected to match LOOPER paper behavior).

Heat-pulse experiments were deprecated after baseline failed to recover stable loops.

---

## Atanas baseline (freely moving) — **negative result**

**Setup**
- Train LOOPER on **first half** of each baseline recording (paper windows).
- Project full trace onto learned scaffold.
- Evaluate Δd, loop assignment stability, and phase-continuity metrics.

**Aggregate results** (from `results/experiment1_all_baseline/summary.csv`, n=21 worms):

- median `mean_post` ≈ **31.46**
- median `post_slope` ≈ **0.06498**
- median `d_peak` ≈ **72.09**
- median `cycles_per_min` ≈ **0.49** (max ≈ **0.95**)
- median `unique_loops` = **2**
- median `loop_switches` = **22**
- median `median_d` ≈ **36.49**

**Interpretation**
- Post-split distances to the scaffold grow large and continue drifting upward.
- Loop ID (`alpha`) is unstable over time (frequent switches).
- Phase continuity metrics are not strong enough to support stable cyclic dynamics.

**Conclusion**
Under this evaluation regime, **LOOPER does not recover stable loop-like dynamics in freely moving Atanas baseline worms**.

---

## Detrending check (single worm, Atanas baseline)

We tested a simple linear detrend (fit on pre-split, applied to full trace) for one worm (2022-06-14-01).

**No detrend**
- `mean_post` = **19.05**
- `post_slope` = **0.064**
- `d_peak` = **50.30**
- `loop_switches` = **6**

**With detrend**
- `mean_post` = **26.69**
- `post_slope` = **0.102**
- `d_peak` = **64.06**
- `loop_switches` = **59**

**Conclusion**
Detrending **did not improve** generalization or loop stability; it worsened drift and alpha switching.

---

## Kato 2015 (immobilized) — **positive control**

LOOPER successfully recovers stable loops in Kato data, consistent with the LOOPER paper.

**Key numbers (shared‑neuron concatenated run):**
- `summary.csv`: `results/kato_shared/summary.csv`
- `n_neurons = 22`, `T = 15662`, `dt_sec = 0.344`
- `unique_loops = 3`, `loop_switches = 96` (segments = 97, mean_segment_len = 194.0)
- `phase_frac_small = 0.995`, `phase_var = 0.111`
- `median_d = 14.97`
- `validation_score_mean = 10.46` (std = 0.815)
- `cycles_per_min = 0.201` (interpret with caution; see note below)

**Comparability caveat (important):**
The paper’s Fig 5B reports **R² reconstruction** (and equivalent #PCs) and validation
correlations on held‑out trials. Our `summary.csv` is based on LOOPER diagnostics
correlations in embedded space and phase‑continuity metrics, so it is **not directly
comparable** to Fig 5B. The cycles/min estimate is especially sensitive to loop switching
and concatenation across worms, so treat it as **qualitative only**.

---

## Artifacts / where to look

- Atanas baseline summary: `results/experiment1_all_baseline/summary.csv`
- Atanas single baseline: `results/experiment1_single_baseline/summary.csv`
- Kato summaries (to generate):
  - `kato_looper_single_eval` → `results/kato_single/summary.csv`
  - `kato_looper_shared_eval` → `results/kato_shared/summary.csv`
