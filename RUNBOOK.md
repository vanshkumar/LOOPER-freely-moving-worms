# RUNBOOK.md

Minimal, LLM-friendly runbooks to execute the LOOPER experiments with raw datasets + metadata.

---

## Shared setup

Inputs (required):
- **Raw dataset** files (per dataset; format varies by source)
- **Metadata** with event times (E1) and/or mode labels (E2)
- LOOPER code: `LOOPER_github_2020/Functions/`

Optional inputs:
- Additional logs (stimulus types/strengths, behavior traces, speed, etc.)

Outputs (current scripts):
- per-run `*.mat` (worm + saveData)
- per-eval `summary.csv`
- diagnostic figures in `results/<run>/diagnostics/`

Core LOOPER outputs to use:
- `BestStateMap` → α/θ (`(:,1)` = loop ID, `(:,2)` = phase bin)
- `BestModel` → transition graph
- `BestEmission` → emission means (for distance proxy)

Minimal LOOPER call (MATLAB):
```
saveData = struct;
saveData.RawData = X';          % channels x time
saveData.TrialData = ones(1, size(X,1));  % single continuous trial
params = struct;                % optional overrides
saveData = LOOPER(saveData, true, [], [], [], params);
```

---

## Current repo scripts (run vs eval)

Run scripts (fit LOOPER + save `*.mat`):
- `kato_looper/kato_looper_single_run.m` → `results/kato_single/kato_single.mat`
- `kato_looper/kato_looper_all_run.m` → `results/kato_all/worm_*.mat`
- `kato_looper/kato_looper_shared_run.m` → `results/kato_shared/kato_shared.mat`
- `experiment1_single_baseline_run.m` → `results/experiment1_single_baseline/experiment1_single_baseline.mat`
- `experiment1_all_baseline_run.m` → `results/experiment1_all_baseline/worm_*.mat`
- `DEPRECATED_experiment1_single_run.m` → `results/experiment1_single/experiment1_single.mat`

Eval scripts (diagnostics + plots):
- `kato_looper/kato_looper_single_eval.m`
- `kato_looper/kato_looper_all_eval.m`
- `kato_looper/kato_looper_shared_eval.m`
- `experiment1_single_baseline_eval.m`
- `experiment1_all_baseline_eval.m`
- `DEPRECATED_experiment1_single_eval.m`

These eval scripts **always** load the fixed output filenames above (no per‑worm or dated filenames).

---

## Experiment 1 runbook: Baseline split stability (current)

Prereqs:
- You must have **ranges_raw** (pre/post split) aligned to the raw dataset timebase.

Steps:
1. Load worms from the **Atanas baseline** dataset:
   - Use `load_atanas_data('baseline')`.
2. For each worm:
   - Fit LOOPER on the **pre** split.
   - Save `saveData` + `worm`.
3. Evaluate:
   - Project the full trace onto the scaffold.
   - Compute split‑based summary metrics and phase continuity metrics.
4. Aggregate:
   - Use `results/experiment1_all_baseline/summary.csv`.

Minimal outputs:
- `results/experiment1_all_baseline/summary.csv`
- `results/experiment1_all_baseline/plots/*.png`

---

## Experiment 2 runbook: On hold

E2 is on hold until we have a dataset where LOOPER reliably recovers stable loops.

---

## Notes / gotchas

- If you lack neuron identity alignment, compare **within‑worm** and aggregate scaffold‑level stats.
- If `ranges_raw` is unavailable, you cannot run Experiment 1 as written.
