# CODEX_SUGGESTIONS (Revised After Reading `December 2025 Aliveline`)

This revision is based on the full daily log timeline, not just repo code/docs.

## What The Logs Show (First-Principles Diagnosis)

Your month produced real progress, but the logs show four recurring bottlenecks:

1. Question quality lagged behind execution speed.
2. Dataset assumptions were often discovered late.
3. Evaluation metrics were selected before metric validity was fully established.
4. LLM output quality was high for implementation, mixed for framing/prioritization.

The most important insight from the logs is not "LLMs made mistakes." It is:

`Execution bandwidth > epistemic bandwidth` most of the time.

That is why code advanced quickly while question quality, control design, and interpretation lagged.

## 1) Process Improvements For LLM Co-Working Across The Whole Science Loop

### A. Add A Mandatory Question Gate Before Any Coding

For every experiment idea, require this one-page spec before touching code:

- `Question`: one sentence, decision-relevant.
- `Testable claim`: falsifiable with current data.
- `Observable`: exactly what quantity will move.
- `Comparator`: what baseline/control must be beaten.
- `Non-identifiability`: what this experiment cannot distinguish.
- `Kill condition`: explicit condition to stop.

This directly addresses the Day 31 realization that the heat-recovery question was under-specified.

### B. Add A Dataset Suitability Gate

Before implementing analysis, fill a short table per dataset:

- recording regime (`immobilized`, `freely moving`, `stimulated`)
- trial structure vs continuous stream
- known nonstationarities and inputs
- required metadata available or missing
- whether the target claim is identifiable in this dataset

This would have flagged the Kato `WT_Stim` vs `NoStim` mismatch earlier and saved iteration loops.

### C. Separate Three LLM Roles

Use three explicit roles in sequence:

- `Architect`: refine question and controls only.
- `Builder`: implement code only.
- `Auditor`: attack leakage/confounds and interpretation.

Do not let `Builder` choose the scientific question.

### D. Enforce A Metric Validation Ladder

A metric can enter [RESULTS.md](https://github.com/vanshkumar/LOOPER-freely-moving-worms/blob/main/RESULTS.md) only after passing:

1. synthetic sanity check,
2. positive-control expectation check,
3. surrogate/null behavior check.

The logs show this was learned empirically after split-half stress tests failed even for Kato.

### E. Track Evidence, Not Narrative

Create `docs/CLAIMS_LEDGER.md` with one row per claim:

- `claim_id`
- `status` (`supported`, `open`, `refuted`)
- `best_evidence_file`
- `controls_passed`
- `known confounds`

This prevents accidental over-interpretation during fast LLM-assisted writing.

### F. Introduce A Weekly Human Checkpoint By Design

The Connor email was a high-leverage intervention and came late.

Make this default:

- one external expert check every 7 days,
- with a fixed template: question, current evidence, what could falsify, next decision.

### G. Add A "Thinking Buffer" To The Workflow

Your reflections explicitly note that 20 minutes of unstructured thinking corrected major framing errors.

Operationalize this:

- After any major run, no coding for 20-30 minutes.
- Write: what the run says, what it does not say, next cheapest disambiguating test.

### H. Add Provenance To Every Summary Output

Append metadata sidecars to each result:

- commit hash,
- script entrypoint,
- parameter hash,
- metric version tag,
- split policy,
- dataset id.

This helps when multiple LLM/model/tool versions are involved.

### I. Build A Lightweight "LLM Arbitration" Rule

When two models disagree:

- require concrete reproduction steps from both,
- run smallest discriminating check,
- log winner and reason in `.codex/LEARNINGS.md`.

This addresses the observed "web model found bug, codex missed it" dynamic.

### J. Process Changes Worth Adding To Repo

- `docs/QUESTION_GATE.md`
- `docs/DATASET_GATE.md`
- `docs/METRIC_VALIDATION.md`
- `docs/CLAIMS_LEDGER.md`
- `templates/WEEKLY_EXPERT_CHECK.md`

## 2) Science Suggestions To Move This Project Forward

### A. Reframe The Core Scientific Program

Current high-level hypothesis is ambitious and multi-layered (`genes -> local rules -> scaffold -> behavior`).

With current data, prioritize claims that are identifiable now:

1. Is scaffold geometry stable within behavior modes?
2. Are failures mostly mode-mixing or true drift?
3. Does scaffold state add predictive value for behavior transitions?

Only then climb toward stronger mechanistic claims.

### B. Replace Global Split-Half As Primary Stationarity Test

Global split-half is a stress test, not the main scientific criterion.

Use three stationarity views:

- `windowed contiguous stationarity`
- `behavior-matched stationarity`
- `event-aligned pseudo-trial stationarity`

Require Kato pass on at least one before applying identical criterion to Atanas.

### C. Decompose Geometry Versus Occupancy

For Atanas and heat:

- test geometry invariance (shape/topology),
- separately test occupancy and phase-velocity changes.

Many apparent scaffold changes are likely occupancy/mixture shifts, not geometry collapse.

### D. Use Atanas Metadata More Aggressively

You have behavior series and encoding-change information.

Use them to:

- stratify runs by behavioral regime,
- relate scaffold instability to behavioral composition and remapping indices,
- test whether instability clusters in specific state-change periods.

### E. Make Heat-Pulse Analysis A State-Perturbation Study

Do not ask "does it return to one global scaffold?" first.

Ask:

- does heat induce temporary regime switch?
- does pre-heat model predict post-heat within matched behavior mode?
- does system return to prior mode occupancy profile over time?

This is better aligned with the Atanas paper’s state-dependent remapping results.

### F. Add One Comparator Family Beyond LOOPER

Use one explicit switching model baseline:

- switching LDS or HMM.

Compare against LOOPER on:

- held-out predictive performance,
- robustness to split policy,
- interpretability of regime switches.

This tests whether LOOPER-specific assumptions are driving conclusions.

### G. Upgrade From Median Summaries To Worm-Level Taxonomy

Classify each worm into:

- stable-within-mode,
- mode-switch dominated,
- drift dominated,
- low-signal/ambiguous.

Then ask what predicts class membership.

This is more scientifically informative than aggregate medians alone.

### H. Define Explicit Falsifiers For Near-Term Claims

Example:

- If behavior-conditioned stationarity still fails on Kato and Atanas, then "mode-mixing explains instability" is weakened.
- If scaffold state does not improve transition prediction over simple baselines, "computational scaffold usefulness" is weakened.

## 8-Week Practical Plan

Week 1:

- Implement question gate, dataset gate, metric validation docs.
- Register 2-3 near-term claims in claims ledger.

Week 2:

- Implement behavior-conditioned and windowed stationarity.
- Re-run Kato and Atanas baseline with provenance tags.

Week 3:

- Add geometry-vs-occupancy decomposition analyses.
- Produce worm-level taxonomy table.

Week 4:

- Add switching-model comparator baseline.
- Run matched evaluations.

Week 5:

- Heat-pulse state-perturbation analyses using mode-matched evaluation.

Week 6:

- Cross-worm neuron-identity anchored analyses where labels permit.

Week 7:

- Consolidate claim status and falsifier outcomes.

Week 8:

- Write updated [RESULTS.md](https://github.com/vanshkumar/LOOPER-freely-moving-worms/blob/main/RESULTS.md) around claim statuses, not around one global narrative.

## Closing

The logs show strong execution capability and strong self-correction. The next leap is to harden the epistemic control system so that LLM speed compounds scientific clarity instead of outrunning it.
