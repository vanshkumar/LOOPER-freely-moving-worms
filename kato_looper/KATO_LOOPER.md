# KATO_LOOPER.md

This folder contains **run/eval scripts for LOOPER on the Kato 2015 dataset** and
documents how they map to the original LOOPER paper’s C. elegans analyses.

This is intentionally **standalone** for the `kato_looper/` subfolder.

---

## What the scripts do

### 1) Single‑worm run + eval (implemented)
- **Run:** `kato_looper_single_run.m`
  - Loads the first Kato worm via `load_kato_data`.
  - Applies per‑neuron detrend across time (pending OSF clarification).
  - Runs LOOPER with paper‑style C. elegans parameters.
  - Saves to `results/kato_single/kato_single.mat` plus the final PCA figure.

- **Eval:** `kato_looper_single_eval.m`
  - Loads the saved `.mat` and runs `run_looper_diagnostics`.
  - Saves diagnostics plots under `results/kato_single/diagnostics/`.
  - Writes `results/kato_single/summary.csv`.

### 2) All‑worms run + eval
Goal: run LOOPER independently on all 5 worms and evaluate each result.

- `kato_looper_all_run.m`
  - Loop over all worms from `load_kato_data`.
  - Save one result per worm:
    - `results/kato_all/worm_1.mat` … `worm_5.mat`
  - Save a manifest with params + worm IDs.

- `kato_looper_all_eval.m`
  - Iterate over `results/kato_all/*.mat`.
  - Run diagnostics and save plots to per‑worm subfolders.
  - Writes `results/kato_all/summary.csv`.

### 3) Shared‑neuron, concatenated run + eval
Goal: reproduce the **paper‑style** worm analysis that appears to concatenate
traces across all 5 worms **using a shared neuron set** (per OSF script).

- `kato_looper_shared_run.m`
  - Compute shared neuron IDs across worms.
  - Subset each worm to the shared set:
    - default: intersection of identified labels (IDs containing letters)
    - note: paper curated 15‑neuron list is **not used** here because VB01 is
      missing in worm 5 of `KATO_WT_NoStim.mat`.
  - Detrend + z‑score per neuron.
  - Concatenate across worms and run LOOPER once.
  - Save mapping: worm offsets/lengths + shared IDs.

- `kato_looper_shared_eval.m`
  - Run diagnostics on the concatenated model.
  - Writes `results/kato_shared/summary.csv`.
  - (Optional) split reconstructions back by worm using saved offsets.

---

## How this relates to the LOOPER paper

The LOOPER paper’s C. elegans results (Fig. 5) report high reconstruction
accuracy. The OSF script `buildWormData.m` suggests that their workflow:
- uses **all 5 Kato worms**,
- restricts to a **shared neuron set**,
- detrends (column‑wise on neurons × time),
- z‑scores each neuron across time,
- and concatenates data before downstream analyses.

This folder mirrors that structure:
- **single‑worm** runs for sanity checks,
- **all‑worms** runs for per‑worm comparison,
- **shared‑neuron concatenation** to align most closely with the paper.

From the OSF checkpoint `wormlooper2.mat` (shared‑worm run):
- shared neuron set size appears to be **18** (RawData = 18 × 15662),
- 5 worm boundaries are encoded in `TrialSwitches`,
- `MinimumReturnTime` is **50** (larger than the paper’s default 10).

Current Kato scripts:
- **Single‑worm** runs (`kato_looper_single_run.m`, `kato_looper_all_run.m`) use
  `MinimumReturnTime = 10` (paper/default).
- **Shared‑worm** run (`kato_looper_shared_run.m`) uses `MinimumReturnTime = 50`
  to match the OSF checkpoint.
- **Shared‑worm** run sets `MaxCheckTime = 1` to avoid `N x N x MaxCheckTime`
  allocations in `reduceMatrix` for long concatenated traces.
- All Kato scripts set `PutativeLoopCounts = [8 7 6 5 4 3 2]` to match
  `wormlooper2.mat`.

Intersection note (from this repo’s `KATO_WT_NoStim.mat`):
- Direct intersection across all 5 worms yields **46** shared labels.
- Restricting to labels that contain letters (likely identified neurons) yields **22**:
  AIBL, AIBR, ALA, ASKR, AVAL, AVAR, AVBL, AVBR, AVEL, AVER, RIBL, RID,
  RIML, RIMR, RIVL, RIVR, RMED, RMEL, RMER, SMDDR, SMDVR, VB02.
- The OSF shared set size (18) likely reflects a **curated list** (neuronList18)
  rather than raw intersection.
 - The **supplement** (journal.pcbi.1010784.s001.pdf) states a **15‑neuron** curated
   set for the Kato analysis: AIBL, AIBR, ALA, AVAL, AVAR, AVBL, AVER, RID, RIML,
   RIMR, RMED, RMEL, RMER, VB01, VB02.

---

## Detrend note (current deviation)

The OSF script detrends a neurons × time matrix **column‑wise**, which removes a
global population mode rather than per‑neuron trends. In this repo we currently
use **per‑neuron detrending across time** because the OSF‑style detrend reduced
reconstruction quality (R² ~0.64 → ~0.49). We will revisit once the authors
clarify the intended preprocessing.

(See `METHODS.md` and `LOOPER_OSF_hosted_2022/LOOPER_PAPER_SCRIPTS.md`.)
