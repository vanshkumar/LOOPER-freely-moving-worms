# EXPERIMENTS.md

These are **prototype research experiments** testing a single core hypothesis:

> Neuron‑level (gene‑shaped) nonlinear feedback controllers yield **conserved macro dynamics** (e.g., locomotion loops / LOOPER scaffolds) with **animal‑specific neural coordinates**, and perturbations reveal **state‑dependent correction** (“controller fingerprints”).

This file is intentionally **high‑level**. For implementation details and exact metrics, see:
- `METHODS.md` (definitions + metrics)
- `RUNBOOK.md` (how to run the scripts)

**Source of truth:** the MATLAB scripts in this repo are authoritative.

---

## Progression (what we are trying to do)

1) **Positive control (Kato, immobilized).**
   Verify that LOOPER + our evaluation metrics recover loop‑like dynamics on a standard dataset.

2) **Target test (Atanas, freely moving baseline).**
   Run the **same fidelity + stationarity protocol** to test whether loops generalize to freely moving worms.

3) **Decision point.**
   If loops are recovered in Atanas under this protocol, proceed to perturbation/heat‑pulse analyses.
   If not, treat as a negative result and pause E2.

---

# Experiment 1 — Baseline split stability test (freely moving)

**Status:** Primary experiment definition for this repo.
The heat‑pulse path is **out of scope** (on hold).

### Targeted subclaim
If a stable locomotion loop exists in freely moving worms, LOOPER should:
- recover a small number of stable loops (α),
- yield smooth phase progression (θ) within loops,
- generalize across time (pre → post split) without large drift from the scaffold.

### Method (high level)
- Split each trace in half (pre/post).
- Fit LOOPER on the pre half.
- Project the full trace and evaluate drift + loop‑ishness metrics.

This split‑half test is intentionally **strict** and is not the same as the paper’s trial‑based validation.

---

# Experiment 2 — Top‑down constraint persists across modes (gain scheduling)

**Status:** **On hold.** E2 is deferred until we have a dataset where LOOPER reliably recovers stable loops.

### High‑level idea (brief)
Across distinct internal modes, the **macro scaffold topology** is conserved while “controller parameters”
change (phase velocity profile, diffusion/noise, branch probabilities). This is **not** “return to baseline.”

---

For details, definitions, and outputs, see `METHODS.md` and `RUNBOOK.md`.
