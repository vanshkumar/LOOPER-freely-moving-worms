# Suggestions for LOOPER-freely-moving-worms

Prepared after reading the full repository (code, docs, results, LEARNINGS.md),
the LOOPER paper (Brennan et al. 2023, PLoS Comp Bio), and its supplement.

---

## Part 1: Improving the LLM-assisted science process

These observations are about the **workflow** of using LLMs as co-workers for
computational neuroscience, based on what went well and what could be improved
in this repo.

### 1.1 What this repo does well

**Separation of truth from explanation.** The convention that "scripts are the
source of truth" and docs are kept in sync is the right call. Many LLM-assisted
projects accumulate stale prose that actively misleads the next session. This
repo's discipline around that is strong.

**LEARNINGS.md as institutional memory.** The signpost-based log in
`.codex/LEARNINGS.md` is an effective pattern. Each entry records a specific
failure mode and its resolution, which prevents the LLM (or a new human
collaborator) from repeating the same mistake. The entries about
`pcaBasis` caching, `Rcorr` vs `Rsquared`, and `saveData` field reattachment
are particularly high-value because they capture non-obvious MATLAB gotchas
that would otherwise cost hours each time.

**Incremental complexity.** Starting with a single-worm fidelity run, then
all-worms, then shared-neuron concatenation, then stationarity is a sound
progression. It let you validate the pipeline at each step before adding
complexity.

### 1.2 Where the LLM-assisted process can improve

#### A. Pre-registration of analysis decisions before coding

The biggest risk with LLM-assisted analysis is that decisions get baked into
code without a deliberate record of *why*. In this repo, several important
choices were made implicitly:

- **Half-split as the stationarity test** (rather than k-fold, rolling-window,
  or trial-style validation). The half-split was chosen, discovered to be too
  strict, and then documented as a "stress test" after the fact. A better
  workflow: before writing any stationarity code, write a short decision
  document ("which stationarity tests should we run and what would pass/fail
  mean?"), get the LLM to critique it, then implement.

- **Delay embedding parameters** (DelayTime=4, DelayCount=5 for Atanas).
  These were chosen via `atanas_autocorr_delay.m` using a 1/e crossing rule.
  The choice is reasonable, but there is no sensitivity analysis showing
  that results are robust to, say, DelayTime=3 or 6. The LLM could have
  been asked to generate a parameter sweep script before committing to a
  single set.

- **Detrending on vs. off for Atanas.** Detrending is disabled for Atanas but
  enabled for Kato. The rationale (Atanas traces are already z-scored) is
  stated, but it is worth checking whether Atanas traces still have slow
  drift that z-scoring per-neuron does not remove (it does not if the drift
  is shared across neurons).

**Recommendation:** Maintain a `DECISIONS.md` (or a decisions section in
EXPERIMENTS.md) where every non-trivial analysis choice is recorded *before*
implementation, with a one-sentence rationale and the alternatives considered.
Ask the LLM to play devil's advocate on each entry before proceeding.

#### B. Explicit positive and negative controls at every stage

The Kato positive control is a good idea, but the repo currently has no
**negative control** (e.g., shuffled data, phase-randomized surrogates, or
isotropic noise) to establish what LOOPER metrics look like when there is
genuinely no loop structure. Without this, it is hard to interpret numbers
like `recon_corr_full = 0.594` for Atanas: is that high relative to chance,
or does LOOPER always produce moderate correlations due to smoothing and delay
embedding?

**Recommendation:** Add a shuffled-data control (temporally shuffle each neuron
independently, destroying dynamics but preserving marginal statistics) and run
the same pipeline. This gives a floor for every metric and makes the fidelity
results interpretable. This is a one-session task for an LLM.

#### C. Use the LLM for literature search and conceptual critique earlier

The LEARNINGS.md log suggests that most LLM interactions focused on debugging
MATLAB code (path issues, caching, struct fields). That is valuable but
underuses the LLM. Before writing the stationarity test, it would have been
productive to ask the LLM:

- "What are the standard ways to test temporal stationarity of a nonlinear
  dynamical model? Which ones are used in the LOOPER paper and in comparable
  methods (SLDs, HMMs, rSLDS)?"
- "Is a half-split the right granularity, given that Kato recordings are ~90
  minutes and C. elegans locomotor bouts last seconds?"

This kind of conceptual pre-flight check is where LLMs add the most value per
token, and it is currently underrepresented in this repo's workflow.

#### D. Structured parameter sensitivity analysis

The repo runs each analysis with a single parameter set. LOOPER's paper
supplement (Fig. SH) shows that only `NearestNeighbors` has a large effect,
but this was tested on primate data, not worm data. Freely moving worm data
is much noisier and shorter per recording, so the sensitivity profile may
differ.

**Recommendation:** Write a generic parameter sweep wrapper (the LLM can
scaffold this in one pass) that varies NearestNeighbors, TotalStates, and
DelayCount, runs the fidelity pipeline, and saves summary.csv for each
combination. Even a coarse 3x3x3 grid on a single representative worm
would be informative. Document which parameters the results are robust to.

#### E. Better figure archiving and narrative linking

The diagnostic figures are saved as PNGs per worm, but there is no
composite summary figure or notebook that stitches them together with
narrative. For LLM-assisted science, a good pattern is:

1. Run the analysis.
2. Ask the LLM to generate a summary figure script that loads all summary.csv
   files and produces a multi-panel overview (distributions of recon_corr,
   phase_frac_small, etc., comparing Kato vs. Atanas).
3. Save that figure alongside the data.

This makes it much easier for a future LLM session (or a human reader) to
quickly assess the state of the project without re-running diagnostics.

#### F. Version-pinning the LOOPER library

The repo includes two copies of the LOOPER code (`LOOPER_github_2020` and
`LOOPER_OSF_hosted_2022`) but does not record which git commit or OSF version
they correspond to. If the upstream code is updated, there is no way to know
whether this repo's results are reproducible.

**Recommendation:** Add a `VERSIONS.md` or note in README with the exact
commit hash / download date for each external dependency.

---

## Part 2: Suggestions on the science

### 2.1 The core finding is important but under-interpreted

The main result is: **freely moving worms show in-sample loop-like structure
comparable to immobilized worms, but the scaffold does not generalize across
a half-split.** This is a genuinely interesting finding. It means one of:

1. The loop structure is **non-stationary** on the timescale of the full
   recording (drift, slow modulation, or mode-switching).
2. The loop structure is stationary, but the **scaffold learned from half the
   data is too noisy** to generalize (underfitting / insufficient data).
3. The loop structure is stationary, but the **projection/distance metric**
   used for evaluation is too sensitive to small drifts in the embedding space.

The repo currently frames this as "stress test too strict" (possibility 3),
but possibilities 1 and 2 are equally plausible and more scientifically
interesting. The next steps should aim to distinguish these.

### 2.2 Specific scientific recommendations

#### A. Behavior-conditioned LOOPER (highest priority)

The Atanas dataset includes rich behavioral covariates (velocity,
angular_velocity, head_curvature, body_curvature, reversal_events,
forwardness, etc.). The most likely explanation for stationarity failure
in freely moving worms is **behavioral mode-switching**: the neural dynamics
look like different loops during forward crawling vs. reversals vs. turns.

**Concrete next step:** Segment each Atanas trace by behavioral state (e.g.,
forward vs. not-forward using the `forwardness` time series), run LOOPER
separately on each behavioral epoch, and test whether within-mode scaffolds
are more stationary. This directly tests the "mode-switching confound"
hypothesis mentioned in RESULTS.md.

If within-mode scaffolds are stable, the finding becomes: "freely moving
worms have mode-specific loops that are individually stationary," which is
a publishable and mechanistically interpretable result.

If within-mode scaffolds are still non-stationary, that points toward genuine
drift or data-length issues, which is also informative.

#### B. Surrogate / null model controls

As mentioned above, run LOOPER on:

1. **Time-shuffled data** (shuffle each neuron independently): destroys
   temporal dynamics, preserves marginal distributions. LOOPER should find
   no loops.
2. **Phase-randomized surrogates** (Fourier surrogate): preserves power
   spectrum and cross-correlations but destroys phase relationships. Tests
   whether loop structure requires temporal coherence beyond what linear
   correlations provide.
3. **AR(1) surrogates** (fit an autoregressive model per neuron, generate
   synthetic data): preserves first-order temporal correlations. Tests
   whether loop structure requires nonlinear dynamics.

If LOOPER finds comparable `recon_corr_full` on surrogates, the fidelity
results are not meaningful. If surrogates have much lower scores, the
fidelity results are validated as reflecting genuine nonlinear structure.

#### C. Windowed stationarity instead of half-split

The half-split is the coarsest possible stationarity test (2 windows). A more
informative approach:

1. **Sliding-window LOOPER:** Fit LOOPER on a window of length W (e.g., 1/4
   of the trace), slide by W/2, and compute reconstruction quality for each
   window. This gives a time course of scaffold quality that can be correlated
   with behavioral state, recording time, etc.

2. **Cross-validation with multiple splits:** Use 5-fold temporal
   cross-validation (train on 4/5, test on 1/5) and report mean held-out
   reconstruction. This is less sensitive to a single bad split point and
   gives error bars.

3. **Leave-one-bout-out:** If behavioral epochs are segmented, train on all
   bouts except one and test on the held-out bout. This is the most natural
   analog of the paper's trial-based validation for continuous data.

#### D. Direct comparison with Kato paper's validation regime

The LOOPER paper validates on C. elegans by training on a subset of trials
and testing on held-out trials (Fig. 5B, supplement Fig. SD). The Kato data
in this repo is continuous (no trials), but the shared-neuron run uses
`TrialData` to encode worm ID, effectively treating each worm as a "trial."

The paper reports R^2 = 0.79 for C. elegans reconstruction (Fig. 5B).
This repo reports `recon_corr_full = 0.661` for the shared run, which is a
correlation (not R^2). Squaring it gives ~0.44, which is substantially lower
than the paper's number. This gap could be due to:

- Different neuron sets (22 vs. paper's 15 or 18).
- Different detrending (per-neuron vs. column-wise).
- Different MinimumReturnTime (50 in shared run vs. 10 in paper).

**Recommendation:** Systematically vary these to find which parameter explains
the gap. If the gap cannot be closed, document it as a reproduction note.
This is important: if the positive control does not fully reproduce, the
negative result on Atanas is harder to interpret.

#### E. Quantify the timescale of scaffold breakdown

The current `post_slope` metric captures whether d(t) grows linearly after
the split, but it does not tell you *when* the scaffold breaks down. A more
informative analysis:

- Compute the held-out reconstruction correlation as a function of time since
  the split (e.g., in 1-minute bins): `recon_corr(t)`.
- Plot this for all worms and both datasets.
- If Kato scaffolds remain valid for 5 minutes but break at 10, and Atanas
  scaffolds break at 2 minutes, that quantifies the stationarity horizon and
  directly tests whether the freely moving context shortens it.

#### F. Neuron-level analysis of loop participation

The Atanas dataset has neuron identity labels. An underexploited analysis:

- Which neurons contribute most to the loop structure? (Compute
  reconstruction quality when each neuron is dropped.)
- Do the "loop neurons" overlap with known motor command interneurons
  (AVA, AVB, RIB, etc.)?
- Are there neurons that contribute to loops in immobilized worms (Kato)
  but not freely moving worms (Atanas), or vice versa?

This leverages the NeuroPAL identity information that makes the Atanas
dataset unique, and connects the abstract loop structure to circuit-level
neuroscience.

#### G. Phase velocity as a behavioral readout

The LOOPER scaffold assigns each timepoint a phase (theta). The instantaneous
phase velocity (d_theta/dt) should correlate with behavioral variables if
the loop reflects locomotion. Specifically:

- Phase velocity should correlate with crawling speed.
- Phase should reverse during reversal events.
- Phase velocity should vary across behavioral modes.

This is a strong test of whether the recovered loop is locomotion-related
(as expected from Kato et al. 2015's work on the motor command sequence)
or an artifact. If phase velocity does not correlate with any behavioral
variable, the loop structure may be capturing something other than locomotion.

#### H. Consider the data quality difference explicitly

Kato recordings are calcium imaging in immobilized worms in a microfluidic
device with ~130 neurons, ~15 min, dt=0.34s. Atanas recordings are in freely
moving worms with NeuroPAL, variable neuron counts (the JSON files), dt~0.6s.
The freely moving data likely has:

- More motion artifacts.
- Lower and more variable SNR.
- Fewer timepoints per behavioral cycle (if locomotion cycles are ~5s,
  that is ~8 samples/cycle at 0.6s dt).

Before attributing stationarity failure to biology, rule out that it is
a **data quality** artifact. One way: subsample the Kato data to match the
Atanas sampling rate and neuron count, add noise to match estimated Atanas
SNR, and re-run the stationarity test. If Kato also fails under these
degraded conditions, the "freely moving worms lack stationary loops"
conclusion is weakened.

### 2.3 Framing for a potential paper

If the above analyses are done, the story could be:

> "LOOPER recovers loop-like structure in freely moving C. elegans
> neural dynamics, but the scaffold is non-stationary on multi-minute
> timescales. Behavioral segmentation reveals that within-mode
> scaffolds are [more stable / still unstable], suggesting that
> [mode-switching / genuine drift / proprioceptive feedback] modulates
> the locomotion loop. Phase velocity from the LOOPER scaffold
> correlates with [crawling speed / behavioral state], confirming
> that the recovered structure reflects the motor command sequence
> rather than a preprocessing artifact."

This is stronger than "loops fail the stress test" because it:
1. Includes a positive statement (loops exist in-sample).
2. Quantifies the stationarity horizon.
3. Connects to behavior via NeuroPAL identity and behavioral covariates.
4. Is grounded in proper null models.

---

## Summary: priority-ordered action items

| Priority | Action | Why |
|----------|--------|-----|
| 1 | Null model controls (shuffled, surrogates) | Without these, fidelity results are uninterpretable |
| 2 | Behavior-conditioned LOOPER (segment by forwardness) | Most likely explanation for stationarity failure; uses existing data |
| 3 | Windowed stationarity (sliding window or k-fold) | Replaces the too-coarse half-split with a graded measure |
| 4 | Phase velocity vs. behavior correlation | Validates that loops reflect locomotion, not artifacts |
| 5 | Parameter sensitivity sweep | Ensures results are robust; one-session LLM task |
| 6 | Kato reproduction gap analysis (recon_corr vs. paper R^2) | Validates the positive control quantitatively |
| 7 | Timescale of scaffold breakdown | Quantifies the stationarity horizon |
| 8 | SNR-matched Kato degradation test | Rules out data-quality confound |
| 9 | Neuron-level dropout analysis | Connects loops to circuit identity |

---

## Process recommendations summary

| Recommendation | One-liner |
|----------------|-----------|
| Pre-register analysis decisions | Write DECISIONS.md before coding |
| Negative controls at every stage | Shuffled data is mandatory |
| Use LLM for conceptual critique early | Ask "what could go wrong" before implementing |
| Parameter sensitivity analysis | Sweep, don't pick a single point |
| Composite summary figures | Auto-generate overview figures from summary.csv |
| Version-pin external code | Record commit hashes for LOOPER library |
