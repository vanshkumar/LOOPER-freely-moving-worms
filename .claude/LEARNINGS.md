# LEARNINGS (Claude)

- 2026-02-05: This repo applies LOOPER (Brennan et al. 2023) to freely moving C. elegans (Atanas NeuroPAL) with Kato 2015 immobilized worms as positive control. The core finding: fidelity passes for both datasets, but half-split stationarity fails for both. Papers are stored as PDFs in the repo root and subdirectories.
- 2026-02-05: LOOPER.m is a headless wrapper around the GUI app's pipeline (preprocess -> diffusion map -> reduce matrix -> find loops). It uses MATLAB scripts-as-functions (e.g., `buildDiffusionMap`, `reduceMatrix`, `buildMinimalModelFromMatrix`) that create variables in the caller workspace. This is fragile and requires careful attention to variable names.
- 2026-02-05: The .codex/LEARNINGS.md file captures MATLAB-specific gotchas (pcaBasis caching, Rcorr vs Rsquared, saveData field reattachment). Always check that file first when debugging pipeline issues.
- 2026-02-05: Atanas JSON files include behavioral covariates (velocity, angular_velocity, forwardness, etc.) that are not yet used in the analysis. These are the key to behavior-conditioned stationarity tests.
- 2026-02-05: The Kato shared-neuron run reports recon_corr=0.661, but the paper reports R^2=0.79. This gap is unexplained and may relate to neuron set (22 vs 15/18), detrending method, or MinimumReturnTime parameter.
