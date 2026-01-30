# LOOPER on freely moving worms (Atanas + Kato)

This repo tests whether **LOOPER** (a method for finding loop‑like neural dynamics)
can recover 1‑D trajectories in **freely moving C. elegans** (Atanas NeuroPAL).
It uses **immobilized Kato 2015** worms as a positive control.

High‑level status:
- **Atanas baseline (freely moving):** LOOPER **fails the strict split‑half stationarity test**
  used here; this is a stress‑test result, not a definitive “no loops” claim.
- **Kato 2015 (immobilized):** LOOPER recovers loop‑like structure in fidelity runs,
  broadly consistent with the paper’s qualitative behavior.

If anything in docs diverges, **scripts are the source of truth**.

---

## Where to start

If you’re new, read in this order:
1. **RESULTS.md** — current findings + key numbers (split‑half stress test vs fidelity).
2. **EXPERIMENTS.md** — the experiment design used here (baseline split stability).
3. **METHODS.md** — processing + evaluation details.
4. **RUNBOOK.md** — how to run the scripts.

---

## Data

- **Atanas NeuroPAL**: `atanas-data/` (per‑worm JSONs; baseline + heat).
  - Stationarity split uses a **half‑split** (same as Kato); `ranges_raw` is only relevant for heat‑pulse runs.
- **Kato 2015**: `kato_2015/` (MAT file; immobilized worms).

Dataset docs:
- `atanas-data/DATASET.md`
- `kato_2015/DATASET.md`

---

## Main scripts

**Terminology** (used below):
- **Fidelity** = train on the full trace, then evaluate on the same trace.
- **Stationarity** = train on the first half, then project the full trace and test whether the second half still matches.

Baseline (freely moving):
- Fidelity: `atanas_single_fidelity.m`, `atanas_all_fidelity.m`
- Stationarity: `atanas_single_stationarity.m`, `atanas_all_stationarity.m`

Kato (positive control):
- Fidelity: `kato_looper/kato_single_fidelity.m`, `kato_looper/kato_all_fidelity.m`
- Stationarity: `kato_looper/kato_single_stationarity.m`, `kato_looper/kato_all_stationarity.m`
- Shared‑neuron concatenation: `kato_looper/kato_shared_run.m` /
  `kato_looper/kato_shared_eval.m`

---

## Key outputs

- Atanas stationarity summary: `results/atanas_all/stationarity/summary.csv`
- Atanas single‑worm stationarity summary: `results/atanas_single/stationarity/summary.csv`
- Kato shared summary: `results/kato_shared/summary.csv`

Diagnostics figures are under each run’s `diagnostics/` folder.

**Note:** large binary artifacts (e.g., `*.mat`) are **not** checked in.
Re‑run the scripts above to regenerate results locally.
