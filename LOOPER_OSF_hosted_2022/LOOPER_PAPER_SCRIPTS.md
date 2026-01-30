# LOOPER Paper (OSF Scripts) — C. elegans Notes

This file summarizes **C. elegans–specific** scripts in `LOOPER_OSF_hosted_2022/` and what they do. The OSF folder contains many scripts for other datasets (primate, mouse LFP, fMRI, RNN, etc.); this document focuses only on the worm‑specific pieces.

---

## C. elegans–specific scripts found

### 1) `buildWormData.m`

**Purpose**: Builds a shared‑neuron, z‑scored dataset across multiple worm recordings (Kato‑style whole‑brain imaging).

**Input sources**:
- `F:/Dropbox/FullBrain/wbdata/` (hard‑coded path).
- Loads `neuronList18` from the same directory.
- Loads all `TS*.mat` files (each contains a `wbData` struct).

**Key processing steps**:
1. **Load raw traces**:  
   `rawData = wbData.deltaFOverF_bc'`  
   (bleach‑corrected ΔF/F, transposed to neurons × time).
2. **Detrend**:  
   `rawData = detrend(rawData);`  
   Note: `rawData` is neurons × time, so MATLAB detrends column‑wise (across neurons
   at each timepoint). This is closer to removing a global population mode than
   per‑neuron detrending across time. In this repo we **do not** apply this OSF
   detrend for Kato runs because it reduced reconstruction quality (R² ~0.64 → ~0.49);
   we currently use per‑neuron detrending pending author clarification.
3. **Identify worm ID** based on `wbData.FlNm` (the file name):
   - suffix `b/c/d/e/f` maps to worm IDs `5/3/4/1/2`.
4. **Shared neuron alignment**:
   - Uses `neuronList(wormID,:)` to select the same set of neurons across worms.
5. **Per‑neuron z‑score**:
   - `sharedData = zscore(sharedData, 0, 2);` (z‑score across time for each neuron).
6. **Concatenate** all worms into a single long matrix.
7. **PCA basis**:
   - `pcaBasis = pca(concatData);`

**Outputs**:
- The script computes PCA basis and explained variance but does **not** save to disk in this version. It appears intended as a preprocessing step for later LOOPER runs or visualization.

**Implications**:
- This script shows that in the original OSF workflow, **worm data were detrended (across neurons per timepoint), z‑scored per neuron**, and **restricted to a shared neuron set** before LOOPER. We currently diverge on detrending for Kato runs due to the R² drop noted above.
- This is consistent with the paper’s Table 1 note that C. elegans data are z‑scored.

---

### 2) `validateData.m` (generic, used across datasets)

While not worm‑specific, this script is the standard LOOPER diagnostic for loop/phase‑over‑time. It assumes:
- `saveData.BestStateMap(:,1)` is loop ID (α),
- `saveData.BestStateMap(:,2)` is phase bin (θ),
- trial boundaries are encoded via `saveData.RawTrialSwitches`.

For continuous worm data, it effectively plots **loop assignment and phase progression over time**. In this repo, diagnostics may skip `validateData` when trial lengths vary, because the OSF script assumes uniform trial lengths.

---

## What is *not* present in OSF for worms

There is no dedicated “run LOOPER on worm data” script in the OSF folder beyond `buildWormData.m`. That suggests the original LOOPER run for C. elegans may have been:
- executed via the GUI or separate scripts not included here, or
- stored in a `.mat` result file that isn’t part of this repo.

---

## Fig 5 (C. elegans) reconstruction note (inferred from generic LOOPER scripts)

The paper reports **R² = 0.79 (3 PCs)** for the C. elegans panel in Fig 5. In LOOPER’s generic reconstruction code (`plotReconstruction.m`), the reported “R²” is actually a **Pearson correlation** between the observed stream and LOOPER reconstruction (not squared). The “(3 PCs)” is typically computed by asking how many PCA components explain the same fraction of variance as the LOOPER correlation; this logic appears in the `ComapreToMarkov*` scripts.

There is **no worm‑specific Fig 5 plotting script** in the OSF folder. However, `buildWormData.m` concatenates **all worm traces** (after detrend + shared‑neuron selection + z‑score) into a single long matrix before computing PCA. This strongly suggests that the LOOPER model used for Fig 5B likely operated on a **concatenated multi‑worm continuous trace** rather than a single worm.

**Inference (based on OSF + paper):** no worm‑specific validation split is described in the paper or OSF scripts. The Fig‑5 caption references trial boundaries (vertical lines), but the C. elegans panel has none, so the most plausible interpretation is that **Fig 5B reports full‑trace reconstruction correlation on the concatenated shared‑neuron worm dataset**, not a held‑out validation split.

---

## Practical takeaways for this repo

- **Preprocessing** from OSF for worms:
  - detrend neurons × time (column‑wise; removes global population mode)
  - z‑score each neuron across time
  - (this repo currently uses per‑neuron detrending instead; see note above)
  - restrict to a shared neuron list (if cross‑worm comparison)
- **No explicit LOOPER parameter block** for worms appears in OSF scripts; the paper’s Table 1/2 remain the authoritative defaults.

If we later need to match OSF more closely, the main actionable step would be to replicate the **detrend + shared‑neuron selection** before running LOOPER on Kato‑style data.

---

## Supplement (journal.pcbi.1010784.s001.pdf): C. elegans details

The supplement clarifies the **exact neuron subset** used for the Kato data:
- They restrict to **15 neurons** that are “unambiguously identified” in every worm:
  **AIBL, AIBR, ALA, AVAL, AVAR, AVBL, AVER, RID, RIML, RIMR, RMED, RMEL, RMER, VB01, VB02**.

This is smaller than:
- the **raw intersection** across 5 worms in `KATO_WT_NoStim.mat` (46 labels),
- the **lettered intersection** across worms (22 labels),
- and the **OSF checkpoint** (`wormlooper2.mat`), which uses **18 neurons**.

So the supplement implies a **curated 15‑neuron list** for the paper figure, which
does not exactly match the OSF checkpoint’s 18‑neuron run.

### Diagnostics script (in this repo)
You can run LOOPER’s standard diagnostics on the OSF checkpoint via:
- `LOOPER_OSF_hosted_2022/run_wormlooper2_diagnostics.m`
  - Outputs to: `results/osf_wormlooper2_diagnostics/`
## `wormlooper2.mat` (OSF output checkpoint)

The file `LOOPER_OSF_hosted_2022/wormlooper2.mat` appears to be a saved LOOPER
output (`saveData`) for the C. elegans analysis. Inspecting it reveals:

- **RawData shape**: `18 × 15662` (neurons × time), consistent with a **shared
  neuron set of size 18** and **concatenated time** across all worms.
- **TrialData length**: `1 × 15662`, with **TrialSwitches length 5**, indicating
  the concatenated stream includes **per‑worm boundaries** (5 worms).
- **FinalStream shape**: `15152 × 198` (after delay embedding and filtering).
- **BestModel shape**: `23 × 23`, **BestEmission**: `23 × 198`.

### Parameters stored in `saveData` (from the MAT file)
Top‑level LOOPER params:
- `NearestNeighbors = 8`
- `UseLocalDimensions = 1`
- `RepopulateDensity = 0.95`
- `MinimumReturnTime = 50`
- `MaxCheckTime = 10`
- `TotalStates = 25`
- `UseTerminalState = 0`

PreprocessData params:
- `ZScore = 1`
- `Smoothing = 1`
- `DelayTime = 10`
- `DelayCount = 10`
- `InputLambda = 1`, `OutputLambda = 1`, `TrialLambda = 1`
- `MergeStarts = 0`, `MergeEnds = 0`

**Implication:** the OSF output confirms a **shared‑neuron concatenation** across
all 5 worms (with trial boundaries) and uses the paper‑style preprocessing
settings; the minimum return time used here (`50`) is **larger than the default
10** often cited for C. elegans in the paper. In this repo we keep **10** for
single‑worm runs and **50** for the shared‑worm concatenation to match OSF.

---
