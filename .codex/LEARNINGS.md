# LEARNINGS

- 2026-01-22: pdftotext on atanas-data/nihms-1958618.pdf reports "Unknown filter Crypt" and only partially extracts text; use another PDF tool/viewer if full text is needed.
- 2026-01-22: Here-documents fail in this sandbox ("can't create temp file for here document"); use `python -c "code = '''...''' ; exec(code)"` instead.
- 2026-01-22: `numpy` is not installed in the environment; use pure-Python JSON inspection or add a dependency later if needed.
- 2026-01-22: `atanas-data/nihms-1958618.pdf` has a PDF filter that breaks `pdftotext`; use `gs` to rewrite a clean copy (e.g., to `/tmp`) before extracting text.
- 2026-01-22: Repo expects MATLAB-only workflows; avoid adding Python loaders unless explicitly approved.
- 2026-01-24: `mv` failed with "Operation not permitted" in this workspace; use `cp` with escalated permissions and then `rm` to rename files.
- 2026-01-26: MATLAB `addpath` does not recurse into subfolders; missing `tubeplot` in `plotLoops` required adding `LOOPER_github_2020/Functions/Library/Tubeplot` explicitly (or using `addpath(genpath(...))`).
- 2026-01-26: OSF `validateData.m` assumes equal-length trials; `wormlooper2.mat` has variable trial lengths, so `run_looper_diagnostics.m` now skips `validateData` for such cases instead of modifying OSF code.
- 2026-01-26: Diffusion-map caching in `LOOPER.m` must include `pcaBasis` (set in `buildDiffusionMap`), otherwise `reduceMatrix` errors; missing `pcaBasis` should invalidate the cache and force a full diffusion-map rerun.
- 2026-01-27: For Experiment 1 detrending with a train/test split, fit the trend on the training window only and apply it to the full trace to avoid leakage.
