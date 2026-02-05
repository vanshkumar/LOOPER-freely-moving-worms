# CODEX_SUGGESTIONS

## Repository Understanding (Current State)

The repository is doing a careful LOOPER replication and stress test in two regimes:

- Positive control: Kato 2015 immobilized worms.
- Target: Atanas 2023 freely moving worms (baseline first, heat on hold).

The same two evaluations are applied to both:

- Fidelity: fit and evaluate on the same trace.
- Stationarity stress test: fit first half, evaluate second half.

Current result pattern (from docs and summaries):

- Fidelity is decent in both datasets.
- The strict split-half stationarity test fails in Atanas and also fails in Kato.
- This means the current stationarity metric is likely stricter than the paper-style validation regime.
- The key scientific bottleneck is now metric validity and mode/state mixing, not just "does LOOPER fit in-sample."

This is a strong and honest framing.

## 1) Process Improvements for Co-Working with LLMs Across the Science Cycle

### A. Put a Hard "Claim Contract" in Front of Every LLM Run

Before any coding prompt, force a short structured preamble:

- `Claim`: what statement you want to support or refute.
- `Observable`: exact metric(s) and file(s) that will move.
- `Comparator`: what baseline/control must be passed.
- `Fail condition`: explicit reason to reject the run.

Why this helps here:

- Your current ambiguity is not coding quality. It is claim quality (stationarity definition vs paper claim).
- A contract prevents LLMs from optimizing the wrong target.

### B. Use a Two-LLM Workflow by Default

For non-trivial analyses, split LLM roles:

- `Builder`: writes code/docs.
- `Skeptic`: only attacks assumptions, leakage, and over-interpretation.

Apply skeptic review to:

- split definitions,
- leakage in preprocessing,
- proxy metrics (`recon_corr_full` usage),
- conclusions in `RESULTS.md`.

This directly matches your current failure mode (metric mismatch, not implementation bugs).

### C. Add a "Metric Provenance" Block to Every Result CSV

Write a sidecar metadata file (or extra columns) on each run with:

- commit hash,
- script entrypoint,
- parameter struct hash,
- dataset split policy,
- metric implementation version.

Why:

- LLM-assisted repos drift quickly in definitions even when code still runs.
- You need to know if two CSV rows are scientifically comparable.

### D. Adopt a Small Test Ladder for Metrics

For every new metric, require 3 tests before interpretation:

1. Synthetic sanity test (known ground truth).
2. Positive control pass/fail expectation (Kato).
3. Negative control behavior (time-shuffled surrogate).

This should be mandatory before LLM-generated narrative text is allowed to mention the metric.

### E. Separate "Exploration" and "Decision" Artifacts

Keep two tracks:

- `explore/*`: fast, LLM-heavy iterations.
- `decision/*`: minimal, reviewed, reproducible artifacts only.

Then only `decision/*` can feed `RESULTS.md` conclusions.

This prevents accidental promotion of promising but unstable analyses.

### F. Add a Structured "Evidence Ledger"

Create one table (markdown or CSV) with columns:

- `claim_id`
- `status` (`supported`, `weak`, `refuted`, `open`)
- `evidence_file`
- `metric`
- `control_passed`
- `notes`

This is the single best protection against LLM overconfident synthesis in science projects.

### G. Prompt Templates for Scientific Rigor

Use fixed prompt headers for LLM tasks:

- "List assumptions first."
- "List 3 ways this could be wrong."
- "State what control would falsify your proposal."
- "Do not interpret before reporting raw effect directions."

This improves non-coding LLM outputs (methods, interpretation, writing), not just code generation.

### H. Repo-Level Process Changes Worth Adding

High value files to add:

- `docs/CLAIMS_LEDGER.md`
- `docs/METRIC_REGISTRY.md`
- `docs/ANALYSIS_DECISION_LOG.md`
- `templates/LLM_TASK_CONTRACT.md`

These are lightweight and directly aligned with how this repo is already operating.

## 2) Science Suggestions to Move the Project Forward

### Priority 1: Re-anchor Stationarity to a Passing Positive Control

Do not treat the current split-half criterion as the main gate yet.

Recommendation:

- Define at least one stationarity criterion that Kato passes reliably.
- Then apply the same criterion to Atanas.

Candidate criteria:

- Windowed stationarity (short rolling windows instead of full half split).
- Behavior-matched split (same behavior composition pre and post).
- Trial-like pseudo-epochs in continuous data (event-aligned segments).

This is the most important next scientific step.

### Priority 2: Explicitly Model Behavioral Mode Mixture in Atanas

Given Atanas includes behavior time series and known state-dependent remapping in the paper:

- Segment by behavior state first (forward, reversal, turning, pumping regime, etc.).
- Fit per-mode scaffold(s).
- Evaluate within-mode and cross-mode transfer separately.

Interpretation target:

- Distinguish true scaffold drift from mode-mixing failure.

Right now those are confounded.

### Priority 3: Use Heat Data as a Controlled Perturbation of State

Heat is currently on hold, but it is likely the strongest testbed for your core question.

Run pre-registered tests:

- Geometry stability: does scaffold topology stay similar pre/post heat?
- Occupancy dynamics: do alpha/theta occupancies and transition probabilities shift?
- Recovery dynamics: does model return to baseline scaffold or move to alternate scaffold?

This connects directly to intrinsic vs feedback-stabilized dynamics framing.

### Priority 4: Compare Multiple Generative Families, Not Only LOOPER

For this question, include at least one comparator that supports regime switching:

- switching LDS / HMM with explicit state transitions,
- mode-conditioned LOOPER,
- simple behavior-conditioned baseline.

Decision criterion should include:

- predictive power on held-out windows,
- interpretability of branch/merge structure,
- robustness to split policy.

### Priority 5: Move From Aggregate Medians to Per-Worm Effect Taxonomy

Your summary medians are useful, but scientifically you likely have subtypes.

Classify worms by:

- stable scaffold,
- multi-scaffold switching,
- strong drift,
- low-signal/noisy.

Then ask what predicts class membership:

- neuron count,
- behavioral composition,
- encoding-change indices from Atanas metadata,
- recording quality proxies.

### Priority 6: Test the "Computational Scaffold" Claim Directly

For each fitted scaffold, test whether branch/merge structure predicts behavior decisions:

- Can scaffold state predict upcoming reversal/turn better than matched low-dim PCA baselines?
- Are branch points aligned with behavior transitions above chance under shuffles?
- Are errors structured (specific confusions), not just noisy?

This moves from geometry description to computation claim.

### Priority 7: Cross-Worm Neuron-Identity Anchoring in Atanas

You have NeuroPAL labels in many recordings. Use them.

- Build shared-neuron analyses analogous to Kato shared run.
- Quantify whether scaffold events map to consistent neuron-class motifs.
- Track which neuron classes are most involved in remapping periods.

This will greatly strengthen biological interpretation.

## Suggested 4-Week Execution Plan

Week 1:

- Finalize stationarity metric registry and positive-control pass criteria.
- Implement 2 alternate stationarity metrics.

Week 2:

- Run behavior-conditioned baseline Atanas analyses.
- Produce within-mode vs cross-mode generalization tables.

Week 3:

- Run pre/post heat scaffold shift analysis with pre-registered metrics.
- Add surrogate controls for each metric.

Week 4:

- Consolidate evidence ledger.
- Update `RESULTS.md` with claim-status labels (`supported`, `weak`, `open`).

## Final Note

The repo is already in a good scientific posture: you did not overclaim from in-sample fits, and you used a strong control. The next leap is to tighten metric governance and explicitly model state/mode structure, which is strongly supported by both the Atanas and Kato papers you are building on.
