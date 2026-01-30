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

**Aggregate results** (from `results/atanas_all/stationarity/summary.csv`, n=21 worms):

- median `mean_post` ≈ **31.46**
- median `post_slope` ≈ **0.06498**
- median `d_peak` ≈ **72.09**
- median `phase_speed_bins_per_min` = **(recompute after eval rerun)**
- median `unique_loops` = **2**
- median `loop_switches` = **22**
- median `median_d` ≈ **36.49**

**Interpretation**
- Post-split distances to the scaffold grow large and continue drifting upward.
- Loop ID (`alpha`) is unstable over time (frequent switches).
- Phase continuity metrics are not strong enough to support stable cyclic dynamics.

**Conclusion (nuanced)**
Under this **strict split‑half stationarity test**, we do **not** see evidence that a single
stable scaffold generalizes across the full recording in freely moving Atanas baseline worms.
This does **not** rule out loop‑like structure on shorter windows or under trial‑style /
pseudo‑trial validation schemes.

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

LOOPER recovers **loop‑like structure** in Kato data (high phase continuity, low switching),
consistent with the paper’s qualitative behavior, even though strict split‑half generalization
does not hold (see note below).

**Key numbers (shared‑neuron concatenated run):**
- `summary.csv`: `results/kato_shared/summary.csv`
- `n_neurons = 22`, `T = 15662`, `dt_sec = 0.344`
- `unique_loops = 3`, `loop_switches = 96` (segments = 97, mean_segment_len = 194.0)
- `phase_frac_small = 0.995`, `phase_var = 0.111`
- `median_d = 14.97`
- `validation_score_mean = 10.46` (std = 0.815)
- `phase_speed_bins_per_min = (recompute after eval rerun)`

**Comparability caveat (important):**
The paper’s Fig 5B reports **R² reconstruction** (and equivalent #PCs) and validation
correlations on held‑out trials. Our `summary.csv` is based on LOOPER diagnostics
correlations in embedded space and phase‑continuity metrics, so it is **not directly
comparable** to Fig 5B. The phase‑speed estimate is especially sensitive to loop switching
and concatenation across worms, so treat it as **qualitative only**.

**Split‑half nuance (Atanas vs Kato):**
The split‑half test is a **strong stationarity test**. Kato still shows loop‑like phase behavior
under this test (high phase continuity), but the learned scaffold does **not** generalize across
the full trace. Atanas shows **both** poor generalization **and** weaker loop‑like phase
structure. So Kato remains a stronger positive control even though split‑half generalization
fails there too.

---

## Artifacts / where to look

- Atanas stationarity summary: `results/atanas_all/stationarity/summary.csv`
- Atanas single stationarity: `results/atanas_single/stationarity/summary.csv`
- Kato summaries (to generate):
## Kato single-worm checks
- **Fidelity (full trace):**
  - `kato_looper/kato_single_fidelity.m` → `results/kato_single/fidelity/summary.csv`
- **Stationarity (split‑half):**
  - `kato_looper/kato_single_stationarity.m` → `results/kato_single/stationarity/summary.csv`
- `kato_looper/kato_shared_eval.m` → `results/kato_shared/summary.csv`
