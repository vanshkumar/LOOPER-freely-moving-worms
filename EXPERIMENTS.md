# Experiments

This document describes the experimental design used in this repo: what we're testing, why, and how the experiments relate to each other.

For implementation details and metric definitions, see `METHODS.md`. For how to run the scripts, see `RUNBOOK.md`. For results, see `RESULTS.md`.

**Source of truth:** the MATLAB scripts are authoritative. If this document diverges from the code, trust the scripts.

---

## The question

LOOPER (Brennan et al. 2023) has been shown to recover clean loop-like trajectories in neural dynamics from **immobilized** C. elegans (Kato et al. 2015). In immobilized worms, neural activity produces relatively stereotyped oscillations — the dynamics are low-dimensional and cyclic, and LOOPER captures them well.

**Freely moving worms are harder.** During locomotion, the animal cycles between behavioral modes (forward crawling, reversals, turns, pauses), proprioceptive feedback modulates neural activity, and the dynamics are expected to be more variable. Can LOOPER still recover meaningful loop-like structure in this setting?

This repo applies the same LOOPER pipeline to both datasets — Kato 2015 (immobilized, positive control) and Atanas et al. 2023 (freely moving, target) — and evaluates two properties:

1. **Fidelity:** Can LOOPER fit loop structure to the data at all? (In-sample test.)
2. **Stationarity:** Does a scaffold learned from one part of the recording generalize to the rest? (Held-out test.)

Together, these tell us whether LOOPER finds structure in freely moving data and whether that structure is temporally stable.

---

## Why these two tests

**Fidelity** is the basic check: does LOOPER find anything? If you train LOOPER on a neural time series and ask it to reconstruct the same time series from its scaffold, how well does it do? High reconstruction correlation means the scaffold captures real variance in the data. Low correlation means LOOPER didn't find meaningful cyclic structure.

This is the same evaluation approach the original LOOPER paper uses for C. elegans data. Because the recordings are continuous (no trial structure), the paper does not use held-out validation for this dataset — it reports in-sample reconstruction on the full concatenated trace (paper Table 1: Trials = "Contiguous", Conditions = N/A). Our fidelity test follows the same methodology.

**Stationarity** asks whether the scaffold generalizes. We train LOOPER on the first half of a recording, then project the full trace onto the learned scaffold. If the scaffold captures genuine, temporally stable dynamics, it should fit the second half reasonably well. If it doesn't, either the dynamics changed (non-stationarity) or the scaffold overfit to the training half.

The split-half stationarity test used here **goes beyond what the paper did for C. elegans**. The paper used trial-based cross-validation for datasets with trial structure (primate, mouse, human fMRI), but C. elegans recordings are continuous — so the paper did not test temporal generalization on worm data at all. Our half-split asks a new question: does a scaffold learned from one temporal half generalize to the other? This is a harder test than anything in the paper because it requires generalization across a long temporal gap (~8-10 minutes), during which slow drift, mode switching, or other non-stationarities can accumulate. Failure on this test does **not** rule out loop-like structure — it means the structure is not stationary on this timescale under this test.

---

## Experimental progression

The experiments follow a controlled progression:

### Step 1: Positive control (Kato, immobilized)

Run LOOPER on the 5 Kato worms. This dataset is well-characterized — the LOOPER paper demonstrated loop recovery on it. If our pipeline doesn't recover loops here, something is wrong with our implementation.

We run both fidelity (full-trace) and stationarity (split-half) on each worm, plus a shared-neuron concatenated run across all 5 worms (the closest approximation to the paper's analysis).

### Step 2: Target test (Atanas, freely moving)

Run the **same pipeline and metrics** on the 21 Atanas baseline worms. The only differences are dataset-specific preprocessing and delay-embedding parameters (scaled to the different sampling rate). The evaluation protocol is identical.

### Step 3: Decision point

Compare the two datasets:

- If LOOPER recovers loops in both datasets (fidelity passes) but the positive control passes stationarity while the target fails, we learn something specific about freely moving dynamics — the loop structure is present but not temporally stable, possibly due to behavioral state switching.
- If both datasets fail stationarity (which is what happened), the test itself is too strict to discriminate. We need a different stationarity criterion — one the positive control reliably passes — before drawing conclusions about the freely moving data.

---

## Experiment 1: Baseline split stability

**Status:** Complete. This is the primary experiment in this repo.

### What we're testing

If LOOPER can recover stable loop structure in freely moving worms, the scaffold should:
- Reconstruct neural activity reasonably well in-sample (fidelity).
- Show smooth, progressive phase cycling within loops.
- Generalize from the first half to the second half of the recording (stationarity).

### Design

For each worm in each dataset:

1. **Fidelity run:** Train LOOPER on the full trace. Evaluate reconstruction correlation and phase continuity on the same trace.
2. **Stationarity run:** Train LOOPER on the first half (pre-split). Project the full trace onto the learned scaffold. Measure reconstruction quality and distance-to-scaffold drift in the second half (post-split).

### What counts as success or failure

There are no fixed thresholds. Instead, we compare the two datasets on the same metrics:

- **Fidelity success:** Positive reconstruction correlations, small number of stable loops, smooth phase progression. Both datasets are evaluated on the same criteria.
- **Stationarity success:** Post-split reconstruction correlation that is positive and not drastically lower than pre-split. Low distance-to-scaffold drift. The scaffold should still "fit" the held-out data.
- **Stationarity failure:** Near-zero or negative post-split reconstruction, monotonically increasing distance-to-scaffold after the split. The scaffold learned from the first half does not describe the second half.

### Why the split-half test is strict

The LOOPER paper validated using trial-based cross-validation (e.g., hold out one trial at a time). This works when there are natural trial boundaries. The Atanas baseline data is a continuous ~16-minute recording with no trial structure, so we use a temporal split instead.

A temporal split-half is a harder test because:
- The held-out data is always temporally distant from the training data (the entire second half of the recording).
- Any slow drift, photobleaching residual, or behavioral state changes that accumulate over minutes will degrade performance.
- The LOOPER scaffold (specifically the diffusion map embedding) is fit to one specific stretch of data and may not extrapolate well even if the underlying dynamics are similar.

Failure on this test is informative (there is temporal non-stationarity at this timescale) but not definitive (it doesn't mean loops are absent).

### Current status and conclusions

**Fidelity:** Both datasets pass. LOOPER finds comparable in-sample loop structure in freely moving and immobilized worms. See `RESULTS.md` Sections 3-4.

**Stationarity:** Both datasets fail, including the positive control. The scaffold does not generalize from the first half to the second half in either case. Because the positive control also fails, we cannot use this test to draw conclusions specific to freely moving worms. See `RESULTS.md` Sections 5-7.

**Implication:** The next step is to find a stationarity criterion that the Kato positive control reliably passes (e.g., trial-style validation, windowed stationarity, or behavior-conditioned scaffolds), then apply that same criterion to the Atanas data.

---

## Experiment 2: Heat-pulse perturbation recovery

**Status: On hold.** This experiment was originally planned to test whether the LOOPER scaffold recovers after a thermal perturbation (using the Atanas heat-pulse recordings). It was deferred for two reasons: (1) the heat pulse in the Atanas data is a nociceptive stimulus that elicits a specific behavioral response, not a clean perturbation to ongoing locomotion dynamics; and (2) the positive control's stationarity failure means we don't yet have a validated baseline to measure "recovery" against.

If a suitable stationarity criterion is established (see above), this experiment could be revisited using either the Atanas heat data or behavior-conditioned analysis of natural perturbation events.

---

## Longer-term direction

The deeper question behind this work is whether the loop-like structure in C. elegans neural dynamics reflects intrinsic attractor dynamics (present even without sensory feedback, as in immobilized worms) or whether it is shaped and stabilized by closed-loop sensory feedback during movement. Addressing this would require perturbation experiments or careful comparison of scaffold structure between conditions — but it first requires establishing that LOOPER recovers temporally stable scaffolds in freely moving data, which is what this repo works toward.
