# Atanas & Kim et al. (Cell 2023) — whole-brain imaging dataset (baseline + heat; NeuroPAL-labeled)

This file summarizes the experimental conditions and data structure described in **Atanas et al., Cell (author manuscript available in PMC 2024-09-14)**, focusing on the **freely moving whole-head calcium-imaging recordings** and the **subset with a single heat pulse**, including the **post-hoc NeuroPAL neuron-ID procedure**.

---

## 0) Repo JSON files (per-worm recordings from WormWideWeb)

Current repo contents include **per-worm JSON files** under:
- `atanas-data/baseline/` (21 worms)
- `atanas-data/heat/` (19 worms)

Total worms: **40**.

Each JSON file corresponds to **one worm recording**. Filenames follow:
`atanas_kim_2023-YYYY-MM-DD-XX.json`, and the internal `uid` matches that stamp.

### 0.1) File-level summary (observed across all files)

Baseline (`atanas-data/baseline/`):
- `num_neurons`: **109–153**
- `trace_array` shape: **num_neurons x T**, where **T = 1600–1615**
- `dataset_type`: `["baseline", "neuropal"]`

Heat (`atanas-data/heat/`):
- `num_neurons`: **120–147**
- `trace_array` shape: **num_neurons x T**, where **T = 1609–1615**
- `dataset_type`: `["heat", "neuropal"]`
- `events.heat`: single index per file (1-based), **802–804**

### 0.2) Common keys and shapes (per file)

Neural traces (neuron-major):
- `trace_array`: list length = `num_neurons`, each inner list length = `T` (timepoints).
  - Values are centered around 0 and include negatives; appears z-scored or normalized.
- `trace_original`: same shape as `trace_array`, but strictly positive and small-magnitude.

Time base:
- `timestamp_confocal`: list length `T`; time between frames is ~0.6 sec.
- `avg_timestep`: ~0.010; **minutes per frame** (0.010 min ~= 0.6 sec).
- `max_t`: integer used by the original analysis scripts (varies per file).

Behavior time series (length `T`):
- `velocity`, `angular_velocity`, `head_curvature`, `body_curvature`, `pumping`
- `reversal_events`: list of `[start, end]` frame indices (1-based)

Per-neuron metrics (length `num_neurons`):
- `dorsalness`, `feedingness`, `forwardness`
- `tau_vals`
- `rel_enc_str_v`, `rel_enc_str_θh`, `rel_enc_str_P`
- `encoding_changing_neurons`: list of neuron indices (1-based)

Neuron labels:
- `labeled`: dict keyed by neuron index as **string** (e.g., `"1"`, `"2"`).
  - Values contain fields like `label`, `neuron_class`, `LR`, `DV`, `confidence`, `roi_id`.

Neuron categorization:
- `neuron_categorization`: nested dict keyed by `"1"` and `"2"` (likely two data splits).
  - Within each split: groups for `v`, `θh`, and `P`, each containing lists of neuron indices.

Dataset tags:
- `dataset_type`: list of two strings:
  - baseline file: `["baseline", "neuropal"]`
  - heat file: `["heat", "neuropal"]`

Trace normalization (from paper + observed in JSONs):
- The paper’s Methods describe extracting ROI mean fluorescence, background subtraction,
  dividing activity-channel by marker-channel to reduce motion artifacts, and then
  bleach-correcting by dividing by a fitted exponential (these traces are referred to as **F**).
- The paper notes that some figures use **F / F_mean**, and that model fits z-score each neuron.
- In these JSONs:
  - `trace_original` is strictly positive with small magnitude and appears consistent with **F**
    (bleach-corrected activity/marker ratio).
  - `trace_array` is zero-mean, unit-variance per neuron (confirmed on sample files), i.e.
    **z-scored per neuron**.

### 0.3) Ranges and event windows (per worm)

Each JSON includes `ranges` with **two lists of 1-based frame indices**, used by
the original analysis as pre/post windows. In this repo, loaders expose this as
`ranges_raw`. The current stationarity pipeline uses a **half‑split** instead of
`ranges_raw`, but the ranges are still relevant for heat‑pulse analyses.

Observed in this repo:
- Baseline: `ranges` is consistently **[1..800]** and **[801..1600]**.
- Heat: `ranges` varies slightly across worms:
  - pre length: **799–800**
  - post length: **798–800**
  - post start: **810–812**

For heat recordings, `events.heat` is a single 1-based frame index (typically 802–804).
The post window often starts a few frames after the heat pulse; do not hardcode these
indices—use `ranges` and `events.heat` from each file.

### 0.4) How to use with LOOPER (this repo)

LOOPER in this repo expects **neurons × time** (channels × time). For these JSONs:
- `trace_array` is already neuron‑major (`num_neurons × T`)
- the loader returns `RawData = trace_array` for direct LOOPER use

Use `timestamp_confocal` (seconds) to derive `dt`, or compute:
- `dt_sec = mean(diff(timestamp_confocal))` (approximately 0.6 sec)
- or `dt_sec = avg_timestep * 60`

---

## 1) What was recorded

### A. Freely moving whole-head neural activity + behavior (main recordings)
- Animals were mounted on a thin agar pad under a coverslip with **microbeads (80 µm diameter)** placed at the pad corners to reduce coverslip pressure.
- During recording, the microscope stage tracked the worm so the head stayed in view while the animal moved.
- Neural activity was recorded via **volumetric spinning-disk confocal imaging** of the head.
- Behavior was recorded simultaneously in an **NIR (near-infrared) brightfield** channel, and behavioral features were later aligned to neural frames by timestamps.

### B. Post-hoc NeuroPAL imaging (subset of recordings)
For “NeuroPAL recordings,” the same animal was first recorded freely moving and then **immobilized by cooling** and imaged with multi-spectral channels to identify neurons via NeuroPAL color fingerprints.

### C. Heat-stimulation recordings (subset of recordings)
A subset of animals received a **single heat pulse to the head** during the freely moving recording.

---

## 2) Animals / strains (high level)

The Key Resources Table lists multiple strains used in the study. For the whole-brain freely moving recordings and the NeuroPAL-labeled subset, the relevant entries include:

- **SWF415**: pan-neuronal nuclear **GCaMP7F** plus several red landmark reporters and **NLS-mNeptune2.5** (used as a reference channel in the imaging system).
- **SWF702**: includes a **low-brightness NeuroPAL** transgene (otIs670) plus pan-neuronal imaging (flvIs17).
- Many strains are in a **lite-1(ce314); gur-3(ok2245)** background (as listed in the resources table).

(See the Key Resources Table for full genotypes and additional strains used for other perturbations/controls in the paper.)

---

## 3) Imaging hardware & acquisition parameters

### Confocal neural imaging
- Microscope: Nikon Eclipse Ti + Andor spinning disk confocal (Yokogawa CSU-X1) with dual-camera emission split.
- Objective: **40× water immersion, 1.15 NA**.
- Acquisition: **322×210 ROI**, **3×3 binning**.
- Volumetric sampling: **80 z-planes**, **0.54 µm** spacing.
- Volume rate: **~1.7 Hz** (some recordings at **~1.4 Hz** depending on objective piezo hardware).

### Behavior channel (NIR)
- Behavior video acquired in reflected **NIR brightfield**.
- Behavioral measurements were computed at **20 Hz** and then mapped to the confocal frame times using timestamps (averaging NIR samples within each confocal frame).

---

## 4) Heat stimulus protocol (head-localized pulse)

- Stimulus: **1436-nm infrared laser** (500 mW source) delivered **a single time** during a recording.
- Beam at sample: expanded to ~**600 µm** at the sample plane; fed into the NIR brightfield path via a dichroic.
- Calibration: temperature at the head region was calibrated using **rhodamine 110 (temperature-insensitive)** and **rhodamine B (temperature-sensitive)** ratio imaging on slides matching the worm-imaging configuration.
- Pulse shape (as characterized for analysis):
  - Mean temperature change during the first second of stimulation was set to **ΔT = 10.0°C**.
  - Exponential decay constant reported as **0.39 s**, returning to baseline within **~3 s**.

---

## 5) Recording structure useful for “pre vs post” analyses

The Methods describe how the authors split recordings into pre/post segments for some analyses (useful if you want standardized time windows across animals):

- **Baseline (no stimulation) datasets**: split in half for “first vs second half” comparisons.
- **NeuroPAL heat-stimulation datasets**:
  - Pre-stim: all timepoints up to just before stimulation (**799–800 timepoints**)
  - Post-stim: a fixed-length window of **800 timepoints** from **(stim+10) to (stim+809)**
- **SWF415 heat-stimulation datasets**:
  - Pre-stim: all timepoints up to just before stimulation
  - Post-stim: often **400 timepoints** from **(stim+10) to (stim+409)**, with alternative handling for recordings containing a gap.

The Methods note that much of the heat-stim analysis emphasizes the NeuroPAL heat dataset due to longer durations and matched pre/post window lengths.

---

## 6) NeuroPAL identification workflow (post-hoc)

After the freely moving recording, animals were cooled and multi-spectral images were collected to construct a NeuroPAL color atlas for that individual.

Key operational details:
- Cooling: agar cooled with a thermoelectric element; a closed-loop controller maintained a **1°C setpoint**, followed by a **~5 min** wait before imaging once immobilized.
- Channels: spectrally isolated image stacks were collected for **mTagBFP2**, **CyOFP1**, and **mNeptune2.5** (with stated excitation lasers and emission filters), plus a combined image used for segmentation/registration that included red markers and GCaMP.
- For each isolated image, **60 timepoints** were recorded and then registered/averaged to improve SNR.
- A 3D RGB composite was built by mapping **mTagBFP2→blue**, **CyOFP1→green**, **mNeptune2.5→red**.

---

## 7) Dataset sizes noted in the manuscript (selected)

The manuscript includes explicit “n” values for the heat stimulus analyses:
- Behavioral event-triggered averages: **32 animals**.
- Heat-stimulus neural responses pooled across **19 animals** (for neuron-class response summaries).

Additionally, the Methods mention **21 baseline NeuroPAL animals** used as training data for one decoding analysis (indicating the scale of the baseline NeuroPAL-labeled subset).

---

## 8) Data & code availability

The Key Resources Table points to:
- **Data**: Zenodo deposit **10.5281/zenodo.8150515** and the WormWideWeb portal.
- **Code**: Zenodo deposit **10.5281/zenodo.8151918** and the Flavell lab GitHub repository (**AtanasKim-Cell2023**).

---

## 9) Practical notes / caveats for downstream modeling

- Two temporal resolutions coexist: neural volumes at ~1.7 Hz vs behavior at 20 Hz; behavior is aligned to neural frames via timestamps and averaging.
- The heat stimulus is **brief and localized** (seconds), while many analyses consider **longer post-stim windows** (hundreds of neural timepoints), enabling studies of sustained state changes and/or recovery.
- NeuroPAL identification is **post-hoc** (after the freely moving recording), using multi-spectral imaging + registration.

---

## 10) Citation

Atanas A., Kim S.S., et al. **Cell (2023)**. Author manuscript available in PMC (2024-09-14). See STAR Methods / Key Resources Table for full experimental details.
