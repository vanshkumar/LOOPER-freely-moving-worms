# EXPERIMENTS.md

These are **prototype research experiments** designed to test two parts of the hypothesis:

> Neuron-level (gene-shaped) nonlinear feedback controllers yield **conserved macro dynamics** (e.g., locomotion loops / LOOPER scaffolds) with **animal-specific neural coordinates**, and perturbations reveal **state-dependent correction** (“controller fingerprints”).

You’ll run these inside a repo that contains:
- the **LOOPER code** (or a wrapper script that runs the LOOPER pipeline), and
- the **raw dataset files + metadata** (event times, mode labels, etc.).

This file specifies:
1) the two experiments,
2) what to compute,
3) what plots / sanity checks to generate,
4) what constitutes “support” vs “falsification” of the targeted subclaims.

If anything here is unclear for implementation, treat it as a prompt for your coding LLM (“implement X exactly as described here”).

**Source of truth:** the MATLAB scripts in this repo are authoritative. This document is kept in sync with the current code; if anything diverges, trust the scripts.

---

## Quick start

Recommended order:

1. Implement the shared pipeline in `METHODS.md` (load → mask handling → embedding → LOOPER fit/projection → scaffold distance/phase).
2. Run **Experiment 1** on the **Atanas 2023 baseline NeuroPAL dataset** (this repo’s `atanas-data/baseline/`).
3. **Experiment 2 is on hold** until we have a dataset where LOOPER recovers stable loops.

Important data‑availability notes:
- For this repo, E1 uses **Atanas baseline** worms; we use `ranges_raw` to split the trace in half.
- Some recommended datasets are **immobilized**; LOOPER scaffolds there reflect internal state
  dynamics, not necessarily locomotion loops.

Data availability checklist (before you start):
- Do you have `ranges_raw` (pre/post split) for E1?
- (E2 on hold) Do you have behavior/mode labels (roaming/dwelling, fed/starved, wake/sleep)?
- Do you have neuron ID mappings for cross‑worm comparisons (otherwise keep analyses within‑worm)?

LLM‑ready minimum viable run:
1. Use **Atanas baseline** worms in `atanas-data/baseline/`.
2. Confirm `ranges_raw` exists for each worm.
3. Load worms, fit LOOPER per worm, save `BestModel`, `BestEmission`, `BestStateMap`.
4. Compute α/θ from `BestStateMap` and run the metrics described below.

Current scripts in this repo:
- `experiment1_single_baseline_run.m` / `experiment1_single_baseline_eval.m`
- `experiment1_all_baseline_run.m` / `experiment1_all_baseline_eval.m`
- `DEPRECATED_experiment1_single_run.m` / `DEPRECATED_experiment1_single_eval.m`

Current outputs (as implemented):
- `results/experiment1_all_baseline/summary.csv`
- `results/experiment1_all_baseline/plots/*.png`
- `results/experiment1_all_baseline/diagnostics/*.png`
- `results/experiment1_single_baseline/summary.csv`

---

# Experiment 1 — Baseline split stability test (freely moving)

**Status:** Baseline freely moving Atanas worms did **not** yield stable LOOPER‑style loops.
Heat‑pulse plan is **deprecated** in this repo.

### Targeted subclaim
If a stable locomotion loop exists in freely moving worms, LOOPER should:
- recover a small number of stable loops (α),
- yield smooth phase progression (θ) within loops,
- generalize across time (pre → post split) without large drift from the scaffold.

### Dataset used in this repo (fixed for now)
We are **only using the Atanas 2023 baseline NeuroPAL dataset** stored at:
- `atanas-data/baseline/` (per‑worm JSONs)

---

## E1.1 Data selection & split

For each baseline worm recording:
1. Use the **paper windows** stored in `ranges_raw` to split the trace in half.
2. Fit LOOPER on the **pre** half only.
3. Project the **full trace** onto the learned scaffold.

---

## E1.2 Core computations (split‑based)

You need these objects (see `METHODS.md` for definitions):
- LOOPER scaffold giving \((\alpha(t), \theta(t))\) for every time point.
- Distance to scaffold \(d(t)\) (transverse deviation proxy).

Compute the following (per worm):
1. **Baseline‑normalized deviation**:
   - \( \Delta d(t) = d(t) - median(d(t)\in pre) \)
   - Evaluate drift after split (`mean_post`, `post_slope`, `d_peak`).
2. **Loop stability**:
   - `unique_loops`, `loop_switches`, `segments`.
3. **Phase continuity (loop‑ishness)**:
   - fraction of small phase steps, phase variance, cycles/min,
   - median \(d(t)\) and mean segment length.

---

## E1.3 Aggregation & summary

Aggregate across worms using `results/experiment1_all_baseline/summary.csv`.
This repo treats large post‑split drift and unstable α as evidence against stable loops.

---

## E1.4 Controls (future)

If we revisit E1 with a dataset that shows stable loops, add:
- event‑time shuffles,
- neuron‑wise circular shifts,
- (if applicable) stimulus‑label swaps.

---

## E1.5 Outcome (current)

Weakened: LOOPER scaffold is unstable across the split and post‑split drift dominates.
See `RESULTS.md` for the negative result summary.

---

# Experiment 2 — Top-down constraint persists across modes (gain scheduling)

**Status:** **On hold.** We are not attempting E2 until we have a dataset where
LOOPER reliably recovers stable loops. The negative result on freely moving
Atanas baseline makes E2 questionable without immobilized or otherwise stable data.

### High-level idea (brief)
Across distinct internal modes (neuromodulatory or metabolic regimes), the **macro scaffold topology**
(loops / merge–split graph) is conserved, while “controller parameters” change:
- phase velocity profile \(v(\theta)\),
- diffusion/noise \(D(\theta)\),
- occupancy/dwell times,
- branch probabilities.

This is **not** “return to baseline.” It is “same macro object, different operating point.”

---

## Minimal outputs to generate (current)

Experiment 1:
- `results/experiment1_all_baseline/summary.csv`
- `results/experiment1_all_baseline/plots/*.png`

Experiment 2:
- On hold.
