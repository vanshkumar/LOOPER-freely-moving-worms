# KATO_LOOPER.md

This folder contains **run/eval scripts for LOOPER on the Kato 2015 dataset** and
documents how they map to the original LOOPER paper’s C. elegans analyses.

This is intentionally **standalone** for the `kato_looper/` subfolder.

---

## What the scripts do

### 1) Single‑worm runs
#### Fidelity (full‑trace reconstruction)
- **Run + eval:** `kato_single_fidelity.m`
  - Loads the first Kato worm via `load_kato_data`.
  - Applies per‑neuron detrend across time (pending OSF clarification).
  - Trains on the **full trace** (paper‑style fidelity).
  - **Goal:** check whether LOOPER can *faithfully reconstruct* an immobilized worm’s
    continuous dynamics (closest to Fig‑5B style fidelity).
  - Saves to `results/kato_single/fidelity/kato_single_fidelity.mat` and writes
    `results/kato_single/fidelity/summary.csv` + diagnostics.
  - Implementation: calls `looper_run_core.m` + `looper_eval_core.m`.

#### Stationarity (split‑half)
- **Run + eval:** `kato_single_stationarity.m`
  - Loads the first Kato worm via `load_kato_data`.
  - Applies per‑neuron detrend across time (pending OSF clarification).
  - Trains on the **first half** and evaluates on the full trace (strict split‑half **stress test**).
  - **Goal:** test whether a scaffold learned early remains valid later in the same recording.
  - Saves to `results/kato_single/stationarity/kato_single_stationarity.mat` and writes
    `results/kato_single/stationarity/summary.csv` + diagnostics.
  - Computes a **held‑out post‑half reconstruction correlation** (`recon_corr_post`) as
    a strict generalization check (not the same as the paper’s trial‑style validation).
  - Implementation: calls `looper_run_core.m` + `looper_eval_core.m`.

### 2) All‑worms runs
Goal: run LOOPER independently on all 5 worms and evaluate each result.

#### Fidelity (full‑trace reconstruction)
- **Run + eval:** `kato_all_fidelity.m`
  - Loop over all worms from `load_kato_data`.
  - Train on the **full trace** for each worm (paper‑style fidelity).
  - **Goal:** check reproducibility of reconstruction fidelity across worms.
  - Save one result per worm (`results/kato_all/fidelity/worm_*.mat`) and a manifest.
  - Writes `results/kato_all/fidelity/summary.csv` + diagnostics.
  - Implementation: calls `looper_run_core.m` + `looper_eval_core.m`.

#### Stationarity (split‑half)
- **Run + eval:** `kato_all_stationarity.m`
  - Loop over all worms from `load_kato_data`.
  - Train on the **first half** for each worm (strict split‑half **stress test**).
  - **Goal:** check scaffold stability across time for each worm.
  - Save one result per worm (`results/kato_all/stationarity/worm_*.mat`) and a manifest.
  - Writes `results/kato_all/stationarity/summary.csv` + diagnostics.
  - Implementation: calls `looper_run_core.m` + `looper_eval_core.m`.

### 3) Shared‑neuron, concatenated run + eval
Goal: reproduce the **paper‑style** worm analysis that appears to concatenate
traces across all 5 worms **using a shared neuron set** (per OSF script).
This is a **fidelity‑only** run; stationarity is not well‑defined once worms
are concatenated.

- `kato_shared_run.m`
  - Compute shared neuron IDs across worms.
  - Subset each worm to the shared set:
    - intersection of identified labels (IDs containing letters)
    - note: the paper’s curated 15‑neuron list is **not used** because worm 5
      lacks **VB01** in `NeuronNames` (confirmed in `KATO_WT_NoStim.mat`).
  - Detrend + z‑score per neuron.
  - Concatenate across worms and run LOOPER once.
  - Save mapping: worm offsets/lengths + shared IDs.
  - **Goal:** approximate the OSF/paper analysis that pools worms into a single
    long continuous trace.

- `kato_shared_eval.m`
  - Run diagnostics on the concatenated model.
  - Writes `results/kato_shared/summary.csv`.
  - (Optional) split reconstructions back by worm using saved offsets.
  - **Goal:** assess the pooled‑worm scaffold and compare to `wormlooper2.mat`.

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

Inference note: the Fig‑5 caption references trial boundaries (vertical lines), but the
C. elegans panel shows none; combined with the lack of worm‑specific validation code
in OSF, we treat **shared‑neuron concatenation + full‑trace reconstruction** as the
closest match to Fig‑5B’s reconstruction metric.

From the OSF checkpoint `wormlooper2.mat` (shared‑worm run):
- shared neuron set size appears to be **18** (RawData = 18 × 15662),
- 5 worm boundaries are encoded in `TrialSwitches`,
- `MinimumReturnTime` is **50** (larger than the paper’s default 10).

Current Kato scripts:
- **Single‑worm** runs (`kato_single_fidelity.m`, `kato_single_stationarity.m`) use
  `MinimumReturnTime = 10` (paper/default).
- **All‑worm** runs (`kato_all_fidelity.m`, `kato_all_stationarity.m`) use
  `MinimumReturnTime = 10` (paper/default).
- **Shared‑worm** run (`kato_shared_run.m`) overrides `MinimumReturnTime = 50`
  (vs paper default 10) to match the OSF checkpoint behavior.
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
