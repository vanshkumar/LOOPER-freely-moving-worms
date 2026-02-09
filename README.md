# LOOPER on Freely Moving C. elegans

LOOPER (Brennan et al. 2023, *PNAS*) extracts loop-like trajectories from high-dimensional neural time series. It has been shown to recover clean cyclic structure in **immobilized** C. elegans (Kato et al. 2015). This repo asks: **can LOOPER also recover loop structure in freely moving worms**, where neural dynamics are more complex due to behavioral state switching, sensory feedback, and proprioception?

We use the Kato 2015 immobilized data as a positive control and apply the same pipeline to the Atanas et al. 2023 NeuroPAL dataset (21 freely moving worms, baseline behavior).

---

## What we found

1. **LOOPER recovers loop structure in freely moving data with reconstruction correlations comparable to immobilized data** (median r = 0.59 vs. 0.57). This uses the same in-sample reconstruction approach the original LOOPER paper uses for C. elegans (continuous traces have no trial structure for held-out validation). A null model (e.g., time-shuffled data) would help calibrate these values but has not been run — this is also the case for the paper's published C. elegans results.

2. **The freely moving scaffolds have noisier phase dynamics.** Phase variance is ~6x higher, cycling is ~2.6x faster, and loop switching is more frequent. This is consistent with more complex behavioral dynamics during movement — but without a null model, we cannot rule out that some of the noisier phase reflects LOOPER fitting to less-cyclic data.

3. **A strict split-half stationarity test fails on both datasets, including the positive control.** A scaffold learned from the first half of a recording does not generalize to the second half. This test goes beyond what the paper did (the paper used trial-based cross-validation, which is not applicable to continuous recordings). The failure of the positive control means this test cannot be used to draw differential conclusions about freely moving worms.

4. **Next steps:** A null model baseline would calibrate the fidelity results, and a stationarity criterion that the positive control passes would enable meaningful comparison of temporal stability between datasets.

See [RESULTS.md](https://github.com/vanshkumar/LOOPER-freely-moving-worms/blob/main/RESULTS.md) for the full analysis with per-worm data tables and diagnostic figure descriptions.

---

## Documentation

Read in this order:

| Document | What it covers |
|----------|---------------|
| **[RESULTS.md](https://github.com/vanshkumar/LOOPER-freely-moving-worms/blob/main/RESULTS.md)** | What we found — per-worm tables, cross-dataset comparisons, interpretation |
| **[EXPERIMENTS.md](https://github.com/vanshkumar/LOOPER-freely-moving-worms/blob/main/EXPERIMENTS.md)** | Why we ran these tests — scientific question, experimental design, decision logic |
| **[METHODS.md](https://github.com/vanshkumar/LOOPER-freely-moving-worms/blob/main/METHODS.md)** | How we computed everything — pipeline, parameters, metrics, diagnostic plots |
| **[RUNBOOK.md](https://github.com/vanshkumar/LOOPER-freely-moving-worms/blob/main/RUNBOOK.md)** | How to reproduce — prerequisites, running scripts, output structure, troubleshooting |

---

## Repo structure

```
LOOPER-freely-moving-worms/
│
├── RESULTS.md / EXPERIMENTS.md / METHODS.md / RUNBOOK.md
│
├── atanas_all_fidelity.m          # Run LOOPER on all 21 Atanas worms (full-trace)
├── atanas_all_stationarity.m      # Run LOOPER on all 21 Atanas worms (split-half)
├── atanas_single_fidelity.m       # Run LOOPER on one Atanas worm (full-trace)
├── atanas_single_stationarity.m   # Run LOOPER on one Atanas worm (split-half)
├── atanas_helpers.m               # Atanas-specific parameters and preprocessing
├── atanas_autocorr_delay.m        # Delay-embedding parameter estimation
│
├── kato_looper/
│   ├── kato_all_fidelity.m        # Run LOOPER on all 5 Kato worms (full-trace)
│   ├── kato_all_stationarity.m    # Run LOOPER on all 5 Kato worms (split-half)
│   ├── kato_single_fidelity.m     # Run LOOPER on one Kato worm (full-trace)
│   ├── kato_single_stationarity.m # Run LOOPER on one Kato worm (split-half)
│   ├── kato_shared_run.m          # Concatenated shared-neuron run (paper-style)
│   ├── kato_shared_eval.m         # Evaluate shared-neuron run
│   └── kato_helpers.m             # Kato-specific parameters
│
├── looper_run_core.m              # Shared runner: data loading → LOOPER fitting → saving
├── looper_eval_core.m             # Shared evaluator: diagnostics → metrics → summary.csv
├── looper_helpers.m               # Shared utilities (split, detrend, projection)
├── run_looper_diagnostics.m       # Generate diagnostic figures + compute metrics
├── LOOPER.m                       # Headless wrapper around LOOPER pipeline
│
├── atanas-data/                   # Atanas 2023 dataset (per-worm JSONs)
│   ├── baseline/                  #   21 freely moving worms
│   ├── heat/                      #   19 heat-stimulated worms (unused here)
│   ├── load_atanas_data.m         #   Data loader
│   └── DATASET.md                 #   Dataset documentation
│
├── kato_2015/                     # Kato 2015 dataset
│   ├── KATO_WT_NoStim.mat         #   5 immobilized worms (56 MB)
│   ├── load_kato_data.m           #   Data loader
│   └── DATASET.md                 #   Dataset documentation
│
├── LOOPER_github_2020/            # LOOPER source code (Brennan lab GitHub)
│   └── Functions/                 #   Core pipeline functions
│
├── LOOPER_OSF_hosted_2022/        # Original paper analysis scripts (OSF)
│   ├── wormLOOPER2.mat            #   Paper's saved LOOPER output (reference)
│   └── buildWormData.m            #   Paper's preprocessing script
│
└── results/                       # All outputs (see RUNBOOK.md for full layout)
    ├── kato_shared/               #   Shared-neuron concatenated run
    ├── kato_all/                  #   Per-worm Kato runs (fidelity + stationarity)
    ├── kato_single/               #   Single Kato worm runs
    ├── atanas_all/                #   Per-worm Atanas runs (fidelity + stationarity)
    ├── atanas_single/             #   Single Atanas worm runs
    └── osf_wormlooper2_diagnostics/  # Reference output from paper's analysis
```

---

## Datasets

**Kato 2015** (positive control): 5 immobilized worms, ~107-135 neurons each, ~18 min recordings at ~2.9 Hz. Data file: `kato_2015/KATO_WT_NoStim.mat`. Source: Kato et al. 2015, *Cell*.

**Atanas 2023** (target): 21 freely moving worms (baseline, no stimulus), ~109-153 NeuroPAL-identified neurons each, ~16 min recordings at ~1.7 Hz. Data files: `atanas-data/baseline/*.json`. Source: Atanas et al. 2023, *Cell*.

See [atanas-data/DATASET.md](https://github.com/vanshkumar/LOOPER-freely-moving-worms/blob/main/atanas-data/DATASET.md) and [kato_2015/DATASET.md](https://github.com/vanshkumar/LOOPER-freely-moving-worms/blob/main/kato_2015/DATASET.md) for full details on each dataset.

---

## Dependencies

- **MATLAB** (tested on R2023b+)
- **Statistics and Machine Learning Toolbox** (for clustering, linkage)
- **LOOPER source code** (included in `LOOPER_github_2020/`)
- No external toolboxes or Python dependencies

The LOOPER functions must be on the MATLAB path. The entry-point scripts handle `addpath` calls.

---

## Quick start

```matlab
% From the repo root in MATLAB:

% 1. Sanity check on one worm (~2-5 min)
kato_single_fidelity

% 2. All Kato worms, fidelity (~15-30 min)
kato_all_fidelity

% 3. All Atanas worms, fidelity (~30-60 min)
atanas_all_fidelity

% 4. Stationarity tests
kato_all_stationarity
atanas_all_stationarity
```

See [RUNBOOK.md](https://github.com/vanshkumar/LOOPER-freely-moving-worms/blob/main/RUNBOOK.md) for expected outputs, caching details, and troubleshooting.

---

## Notes

- **Scripts are the source of truth.** If documentation diverges from code, trust the scripts.
- **Results overwrite on re-run.** Scripts use fixed filenames with no timestamps.
- **Large files are not checked in.** `.mat` result files must be generated locally by running the scripts.
- Papers: `journal.pcbi.1010784.pdf` (LOOPER), `atanas-data/nihms-1958618.pdf` (Atanas), `kato_2015/PIIS0092867415011964.pdf` (Kato).
