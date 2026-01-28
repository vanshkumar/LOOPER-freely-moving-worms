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
  - Baseline split uses `ranges_raw` from the paper.
- **Kato 2015**: `kato_2015/` (MAT file; immobilized worms).

Dataset docs:
- `atanas-data/DATASET.md`
- `kato_2015/DATASET.md`

---

## Main scripts

Baseline (freely moving):
- Single‑worm: `experiment1_single_baseline_run.m` / `experiment1_single_baseline_eval.m`
- All‑worms: `experiment1_all_baseline_run.m` / `experiment1_all_baseline_eval.m`

Kato (positive control):
- Single‑worm: `kato_looper/kato_looper_single_run.m` / `kato_looper/kato_looper_single_eval.m`
- All‑worms: `kato_looper/kato_looper_all_run.m` / `kato_looper/kato_looper_all_eval.m`
- Shared‑neuron concatenation: `kato_looper/kato_looper_shared_run.m` /
  `kato_looper/kato_looper_shared_eval.m`

Deprecated heat‑pulse path (kept for reference only):
- `DEPRECATED_experiment1_single_run.m` / `DEPRECATED_experiment1_single_eval.m`

---

## Key outputs

- Atanas baseline summary: `results/experiment1_all_baseline/summary.csv`
- Atanas single‑worm summary: `results/experiment1_single_baseline/summary.csv`
- Kato shared summary: `results/kato_shared/summary.csv`

Diagnostics figures are under each run’s `diagnostics/` folder.

**Note:** large binary artifacts (e.g., `*.mat`) are **not** checked in.
Re‑run the scripts above to regenerate results locally.
