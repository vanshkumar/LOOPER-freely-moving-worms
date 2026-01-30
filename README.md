# LOOPER on freely moving worms (Atanas + Kato)

This repo investigates whether **LOOPER** can recover 1‑D loop dynamics from
**freely moving C. elegans** (Atanas NeuroPAL) and uses **immobilized Kato 2015**
worms as a positive control.

High‑level status:
- **Atanas baseline (freely moving):** LOOPER does **not** recover stable loops
  under the split‑half evaluation used here (negative result).
- **Heat‑pulse path:** deprecated (kept only as reference).
- **Kato 2015 (immobilized):** LOOPER recovers loop‑like structure in shared‑neuron runs,
  broadly consistent with the paper’s qualitative behavior.

If anything in docs diverges, **scripts are the source of truth**.

---

## Where to start

- **RESULTS.md** — current findings + key numbers (Atanas negative result, Kato positive control).
- **EXPERIMENTS.md** — the actual experiment design used here (baseline split stability).
- **METHODS.md** — processing + evaluation details.
- **RUNBOOK.md** — how to run the scripts.

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
