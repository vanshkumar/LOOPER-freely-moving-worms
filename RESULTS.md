# RESULTS.md

Summary of results from the LOOPER experiments in this repo.  
Updated: 2026-01-30

---

## TL;DR (plain language)

LOOPER can reconstruct loop‑like dynamics **within the same trace** (fidelity) for both
Kato (immobilized) and Atanas (freely moving). But the **strict split‑half test (stationarity)**
(train on first half, check second half) fails **even for Kato**, so this is a
stress test that is stricter than the paper’s trial‑style validation.

## Scope

We ran the same **fidelity** and **stationarity** evaluations on two datasets:

- **Kato 2015 (immobilized) – positive control**
- **Atanas NeuroPAL baseline (freely moving) – target dataset**

**Fidelity** = train on full trace, evaluate reconstruction / loop metrics on the same trace.  
**Stationarity** = train on first half, project full trace, evaluate post‑half drift away from the scaffold.

Heat‑pulse experiments are **on hold** after the baseline stationarity stress‑test results (see below).

---

## Kato shared‑neuron concatenated run (paper‑style repro)
`results/kato_shared/summary.csv`

- `n_neurons` = **22**, `T` = **15662**, `dt_sec` = **0.344**
- `recon_corr_full` = **0.661**
- `unique_loops` = **2**, `loop_switches` = **87**
- `phase_frac_small` = **0.993**, `phase_var` = **0.133**

**Setup note:** This concatenates all 5 worms on the intersection of identified neurons (22) with per‑worm detrend+z‑score before concatenation; this is the closest run to the paper’s shared‑neuron setup (though not an exact reproduction).

**Interpretation:** Full‑trace reconstruction is strong and phase continuity is high. Loop switching is higher because this run concatenates multiple worms into a single trace.

---

## Kato 2015 (immobilized) — positive control

### Fidelity (all worms, n=5)
`results/kato_all/fidelity/summary.csv`

- median `recon_corr_full` = **0.572**
- median `unique_loops` = **3**, `loop_switches` = **6**
- median `phase_frac_small` = **0.997**, `phase_var` = **0.050**
- median `phase_speed_bins_per_min` = **5.81**

**Interpretation:** LOOPER recovers stable loop‑like structure in immobilized worms under the full‑trace fidelity test (positive control passes).

### Stationarity (split‑half, n=5) — **stress/generalization test**
`results/kato_all/stationarity/summary.csv`

- median `mean_post` = **81.10**
- median `post_slope` = **0.263**
- median `d_peak` = **164.69**
- median `recon_corr_post` = **−0.025**
- median `phase_frac_small` = **0.994**, `phase_var` = **0.048**

**Interpretation:** The **strict split‑half stationarity test fails** (large drift, poor post‑half reconstruction), even though the phase metrics still look loop‑like. This indicates the split‑half test is a **stress/generalization criterion** that is **stricter** than the paper’s trial‑style validation, and may not align with how persistence was measured in the paper.

---

## Atanas baseline (freely moving) — fails the split‑half stress test

### Fidelity (full trace, n=21)
`results/atanas_all/fidelity/summary.csv`

- median `recon_corr_full` = **0.594**
- median `unique_loops` = **2**, `loop_switches` = **14**
- median `phase_frac_small` = **0.986**, `phase_var` = **0.291**
- median `phase_speed_bins_per_min` = **15.11**

**Interpretation:** Fidelity metrics are **comparable in magnitude to the Kato positive control** (same‑trace reconstruction and high phase continuity), suggesting loop‑like structure is present in‑sample, though phase smoothness is weaker than Kato. Fidelity alone does not establish **temporal stability** or generalization.

### Stationarity (split‑half, n=21) — **stress test**
`results/atanas_all/stationarity/summary.csv`

- median `mean_post` = **31.46**
- median `post_slope` = **0.06498**
- median `d_peak` = **72.09**
- median `recon_corr_post` = **−0.045**
- median `unique_loops` = **2**, `loop_switches` = **22**
- median `phase_frac_small` = **0.973**, `phase_var` = **0.309**

**Interpretation:** The learned scaffold **does not generalize** to the second half of the trace under this strict criterion. Distance‑to‑scaffold grows and loop assignments are unstable. This supports the statement: **freely moving baseline worms fail the split‑half stress test**, but it does **not** by itself rule out loop‑like structure on shorter windows or under trial‑style validation.

Negative `recon_corr_post` means the held‑out reconstruction is **poor** (it does not match the post‑half data).

Note: split‑half drift can reflect either **true scaffold drift** or **mode mixing** (switching between multiple loop‑like regimes). A failure here does not distinguish those cases.

---

## Overall conclusion (relative to EXPERIMENTS)

1. **Positive control succeeds for fidelity** on Kato immobilized worms, consistent with the paper’s reported behavior.
2. **Strict split‑half stationarity fails** even for Kato worms, so this test is a **stress criterion** that is stronger than the paper’s trial‑based validation.
3. **Atanas fidelity shows loop‑like structure** with metrics comparable in magnitude to Kato (in‑sample fit only).
4. **Atanas stationarity fails the split‑half stress test** and shows weaker phase‑continuity metrics than Kato. This is evidence **against** long‑range stationarity under this criterion, not definitive evidence that loops are absent, and may reflect a metric stricter than the paper’s validation regime.
5. Heat‑pulse analyses are **on hold** pending a better‑matched stationarity test or behavior‑segmented analyses.

---

## Next step (stationarity metric)

The current split‑half stationarity test is a **stress test** that even the Kato positive control fails.  
Before drawing stronger conclusions about Atanas stationarity, we should define a **stationarity criterion that Kato reliably passes**, then apply the same criterion to Atanas. Two candidate directions:

- **Paper‑style trial validation** (closer to LOOPER’s reported regime).
- **Windowed stationarity** (shorter contiguous windows rather than full‑trace split‑half).
- **Behavior‑conditioned stationarity** (segment by locomotor mode / behavior state using Atanas behavior time series; fit/validate within‑mode, or allow multiple scaffolds).

This keeps the positive control meaningful and avoids comparing Atanas to an overly strict bar.

---

## Artifacts / where to look

- Atanas fidelity summaries: `results/atanas_all/fidelity/summary.csv`
- Atanas stationarity summaries: `results/atanas_all/stationarity/summary.csv`
- Kato fidelity summaries: `results/kato_all/fidelity/summary.csv`
- Kato stationarity summaries: `results/kato_all/stationarity/summary.csv`
- Kato shared run: `results/kato_shared/summary.csv`
