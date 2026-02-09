# Runbook

How to run the experiments in this repo and what to expect.

For what these experiments test and why, see `EXPERIMENTS.md`. For pipeline and metric details, see `METHODS.md`.

---

## Prerequisites

- **MATLAB R2023b+** (older versions may work but are untested)
- **Statistics and Machine Learning Toolbox** (used by LOOPER for clustering)
- **Data files in place:**
  - `kato_2015/KATO_WT_NoStim.mat` (56 MB) — the Kato dataset
  - `atanas-data/baseline/*.json` (21 files) — the Atanas dataset
- **LOOPER source code** is included in `LOOPER_github_2020/` — no separate download needed

The entry-point scripts add the necessary paths (`LOOPER_github_2020/Functions/`, etc.) automatically. If you run into "function not found" errors, make sure you're running from the repo root.

---

## Running experiments

Each experiment is run by calling a single MATLAB script from the repo root. Each script runs LOOPER and then evaluates the results, writing a `summary.csv` and diagnostic figures.

### Kato 2015 (positive control)

| Script | What it does | Approx. runtime |
|--------|-------------|-----------------|
| `kato_looper/kato_single_fidelity.m` | One worm, fidelity (full-trace fit) | 2-5 min |
| `kato_looper/kato_single_stationarity.m` | One worm, stationarity (split-half) | 2-5 min |
| `kato_looper/kato_all_fidelity.m` | All 5 worms, fidelity | 15-30 min |
| `kato_looper/kato_all_stationarity.m` | All 5 worms, stationarity | 15-30 min |
| `kato_looper/kato_shared_run.m` | Shared-neuron concatenation (fidelity only) — fit step | 10-20 min |
| `kato_looper/kato_shared_eval.m` | Shared-neuron concatenation — evaluation step | 1-2 min |

The shared-neuron run is two scripts because the fit (`kato_shared_run.m`) and evaluation (`kato_shared_eval.m`) are separated. Run them in order.

### Atanas 2023 (freely moving, target)

| Script | What it does | Approx. runtime |
|--------|-------------|-----------------|
| `atanas_single_fidelity.m` | One worm, fidelity | 2-5 min |
| `atanas_single_stationarity.m` | One worm, stationarity | 2-5 min |
| `atanas_all_fidelity.m` | All 21 worms, fidelity | 30-60 min |
| `atanas_all_stationarity.m` | All 21 worms, stationarity | 30-60 min |

### Recommended order

If you're running everything from scratch:

```matlab
% 1. Quick sanity check — does LOOPER run at all?
kato_single_fidelity

% 2. Positive control — does our pipeline reproduce known results?
kato_all_fidelity
kato_shared_run
kato_shared_eval

% 3. Target — does LOOPER find loops in freely moving data?
atanas_all_fidelity

% 4. Stationarity — do scaffolds generalize across time?
kato_all_stationarity
atanas_all_stationarity
```

---

## What each script outputs

Every run produces three types of output:

### 1. Saved model (`.mat` files)

Per-worm LOOPER output containing the `saveData` struct (scaffold, emission model, state assignments, etc.) plus metadata (`worm` struct with `uid`, `dt_sec`, `TrainSplit` info).

- Single-worm: `results/<dataset>_single/<mode>/<dataset>_single_<mode>.mat`
- All-worm: `results/<dataset>_all/<mode>/worm_<uid>.mat` (one per worm)
- Shared-neuron: `results/kato_shared/kato_shared.mat`

### 2. Summary CSV

One row per worm with all evaluation metrics. See `METHODS.md` Section 5 for column definitions.

- `results/<dataset>_<scope>/<mode>/summary.csv`

### 3. Diagnostic figures

PNG files showing loop structure, reconstruction quality, phase continuity, and (for stationarity) distance-to-scaffold drift.

- Per-worm diagnostics: `results/<dataset>_<scope>/<mode>/diagnostics/<worm>/*.png`
- Stationarity recovery plots: `results/<dataset>_<scope>/<mode>/plots/*.png`
- PCA stream plots: `results/<dataset>_<scope>/<mode>/*_final_stream_pca.png`

See `METHODS.md` Section 8 for what each plot type shows.

---

## Output directory structure

```
results/
├── kato_single/
│   ├── fidelity/
│   │   ├── kato_single_fidelity.mat
│   │   ├── summary.csv
│   │   ├── diagnostics/
│   │   └── final_stream_pca.png
│   └── stationarity/
│       ├── kato_single_stationarity.mat
│       ├── summary.csv
│       ├── diagnostics/
│       ├── plots/          (recovery curves, loop assignment over time)
│       └── projections/    (full projection .mat files)
│
├── kato_all/
│   ├── fidelity/
│   │   ├── worm_1.mat ... worm_5.mat
│   │   ├── manifest.mat
│   │   ├── summary.csv
│   │   ├── cache/          (diffusion map caches)
│   │   └── diagnostics/    (per-worm subdirectories)
│   └── stationarity/
│       ├── (same structure as fidelity, plus plots/ and projections/)
│
├── kato_shared/
│   ├── kato_shared.mat
│   ├── diffusion_cache.mat
│   ├── summary.csv
│   └── diagnostics/
│
├── atanas_single/
│   ├── fidelity/ and stationarity/   (same structure as kato_single)
│
├── atanas_all/
│   ├── fidelity/ and stationarity/   (same structure as kato_all, but 21 worms)
│
└── osf_wormlooper2_diagnostics/
    ├── wormlooper2_diag.mat
    └── *.png                          (reference output from paper's analysis)
```

---

## Complete file index

### Summary CSVs

| File | Dataset | Test | Rows |
|------|---------|------|------|
| `results/kato_shared/summary.csv` | Kato shared-neuron | Fidelity | 1 |
| `results/kato_all/fidelity/summary.csv` | Kato per-worm | Fidelity | 5 |
| `results/kato_all/stationarity/summary.csv` | Kato per-worm | Stationarity | 5 |
| `results/kato_single/fidelity/summary.csv` | Kato worm 1 | Fidelity | 1 |
| `results/kato_single/stationarity/summary.csv` | Kato worm 1 | Stationarity | 1 |
| `results/atanas_all/fidelity/summary.csv` | Atanas per-worm | Fidelity | 21 |
| `results/atanas_all/stationarity/summary.csv` | Atanas per-worm | Stationarity | 21 |
| `results/atanas_single/fidelity/summary.csv` | Atanas worm 1 | Fidelity | 1 |
| `results/atanas_single/stationarity/summary.csv` | Atanas worm 1 | Stationarity | 1 |

### Diagnostic image counts

| Directory | PNGs | Content |
|-----------|------|---------|
| `results/kato_shared/diagnostics/` | 4 | Shared-neuron loop structure |
| `results/kato_all/fidelity/diagnostics/` | ~25 | 5 worms × 5 plots |
| `results/kato_all/stationarity/diagnostics/` | ~25 | 5 worms × 5 plots |
| `results/kato_all/stationarity/plots/` | 10 | 5 worms × (recovery + assignment) |
| `results/atanas_all/fidelity/diagnostics/` | ~105 | 21 worms × 5 plots |
| `results/atanas_all/stationarity/diagnostics/` | ~105 | 21 worms × 5 plots |
| `results/atanas_all/stationarity/` | ~42 | PCA stream plots + MATLAB figs |
| `results/osf_wormlooper2_diagnostics/` | 4 | Reference output from paper |

### Scripts → outputs mapping

| Script | What it generates |
|--------|------------------|
| `kato_looper/kato_shared_run.m` + `kato_shared_eval.m` | `results/kato_shared/` |
| `kato_looper/kato_all_fidelity.m` | `results/kato_all/fidelity/` |
| `kato_looper/kato_all_stationarity.m` | `results/kato_all/stationarity/` |
| `kato_looper/kato_single_fidelity.m` | `results/kato_single/fidelity/` |
| `kato_looper/kato_single_stationarity.m` | `results/kato_single/stationarity/` |
| `atanas_all_fidelity.m` | `results/atanas_all/fidelity/` |
| `atanas_all_stationarity.m` | `results/atanas_all/stationarity/` |
| `atanas_single_fidelity.m` | `results/atanas_single/fidelity/` |
| `atanas_single_stationarity.m` | `results/atanas_single/stationarity/` |

---

## Caching

The diffusion map computation (the most expensive step) is cached per-worm to avoid redundant recomputation.

- Cache location: `results/<dataset>_<scope>/<mode>/cache/worm_<uid>_diffusion_cache.mat`
- The shared-neuron run caches to `results/kato_shared/diffusion_cache.mat`

**When to clear the cache:** If you change any of the following, delete the corresponding cache files before re-running:
- Preprocessing (detrending, z-scoring)
- Delay embedding parameters (`DelayTime`, `DelayCount`)
- `NearestNeighbors` or `MinimumReturnTime`
- The input data itself

To clear all caches:
```matlab
% From repo root
delete('results/*/cache/*.mat')
delete('results/*/*/cache/*.mat')
delete('results/kato_shared/diffusion_cache.mat')
```

---

## Troubleshooting

**"Function not found" errors:** Make sure you're running from the repo root. The entry-point scripts call `addpath` for `LOOPER_github_2020/Functions/` and its subdirectories. If running a helper function directly, you may need to add paths manually.

**Out-of-memory on long traces:** The diffusion map step builds an N×N matrix where N is the number of (embedded) time points. For the shared-neuron concatenated run (~15,000 points), this requires ~1.8 GB. If MATLAB runs out of memory, try closing other applications or running single-worm analyses first.

**Results look different after re-running:** Scripts overwrite old results with no timestamps. If you changed parameters between runs, the old results are gone. Consider copying the `results/` directory before re-running with different settings.

**Stale diffusion cache:** If results look wrong after changing parameters, check whether the cache is still present. The cache stores the diffusion map from a previous run and will be re-used even if parameters changed. Delete the cache and re-run.

**LOOPER produces only 1 loop:** This occasionally happens when the dynamics are not strongly cyclic or when `MinimumReturnTime` is too large relative to the natural cycle period. It is not necessarily an error — some worms may genuinely lack multi-loop structure.
