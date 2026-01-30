# RUNBOOK.md

Minimal runbook for running LOOPER experiments in this repo.  
Pick a script below; each runs LOOPER and the evaluation.

---

**Quick definitions**
- **Fidelity** = train on the full trace, evaluate on the same trace.
- **Stationarity** = train on the first half, evaluate the full trace (stress test).

---

## Kato 2015 (positive control)

**Fidelity**
- Single worm: `kato_looper/kato_single_fidelity.m`
- All worms: `kato_looper/kato_all_fidelity.m`
- Shared‑neuron concatenation (fidelity‑only): `kato_looper/kato_shared_run.m` → `kato_looper/kato_shared_eval.m`

**Stationarity**
- Single worm: `kato_looper/kato_single_stationarity.m`
- All worms: `kato_looper/kato_all_stationarity.m`

---

## Atanas NeuroPAL (target test)

**Fidelity**
- Single worm: `atanas_single_fidelity.m`
- All worms: `atanas_all_fidelity.m`

**Stationarity**
- Single worm: `atanas_single_stationarity.m`
- All worms: `atanas_all_stationarity.m`

---

## Outputs

Run + eval scripts write:
- `results/<dataset>_<scope>/<mode>/summary.csv`
- `results/<dataset>_<scope>/<mode>/diagnostics/*.png`
- stationarity plots in `results/<dataset>_<scope>/<mode>/plots/`

---

## Notes

- Scripts use fixed filenames (no timestamps) and overwrite old results.
- Diffusion map caching is enabled for all‑worm runs. If you change preprocessing/params,
  clear the corresponding `results/<dataset>_<scope>/<mode>/cache/` folder to avoid stale reuse.
- If you later re‑enable heat‑pulse analyses, you will need event timing metadata.
