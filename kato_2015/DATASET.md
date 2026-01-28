# Kato et al. 2015 (Cell) — immobilized C. elegans whole-brain calcium dataset

This folder contains the **Kato 2015** immobilized-worm recordings used as the
canonical C. elegans dataset in the LOOPER paper.

Key reference:
- Kato S. et al. *Cell* (2015), “Global Brain Dynamics Embed the Motor Command
  Sequence of *C. elegans*.” PDF: `PIIS0092867415011964.pdf`.

---

## 1) What was recorded (paper summary)

From the Kato 2015 paper:
- **Whole‑brain, single‑cell calcium imaging** with a pan‑neuronal nuclear
  Ca²⁺ sensor.
- Animals were **immobilized in a microfluidic device** under **constant
  conditions** (no external stimulus).
- **18 minutes** per animal at **~2.85 volumes/sec**.
- **n = 5 worms** in the paper (and in `KATO_WT_NoStim.mat` in this repo).
- **107–131 neurons** detected per recording.
- The neural state trajectory is **low‑dimensional and cyclical**, forming a
  manifold in PCA space (the classic “loop” reference).
  - At 2.85 Hz, the frame interval is **~0.344 s**.

This dataset provides a canonical **immobilized, no‑stimulus** reference
recording with stable, cyclical population dynamics under constant conditions.

---

## 2) Files in this folder

- `KATO_WT_NoStim.mat` — main data file (wild‑type, no stimulus).
- `load_kato_data.m` — loader for this dataset:
  - `X` = time × neuron
  - `RawData` = neuron × time (LOOPER‑ready)
- `PIIS0092867415011964.pdf` — the paper.

---

## 3) MAT structure (as used by load_kato_data.m)

The loader loads:

```
S = load(matPath);
S = S.WT_NoStim;
```

Expected fields in `S`:
- `deltaFOverF_bc` — calcium traces (ΔF/F, bleach‑corrected).
  - Can be 2‑D (T × N) or 3‑D (T × N × W), or a cell array per worm.
- `tv` — time vector per worm (seconds).
- `NeuronNames` — neuron IDs (per worm).
- `dataset` — dataset name/ID (per worm).
- `fps` — sampling rate (per worm).

The helper `coerce_trace_to_t_by_n_` in `load_kato_data.m` converts traces to
**time × neuron** and determines `dt` from `tv` or `fps` (typically ~0.344 s).

---

## 4) Notes / caveats

- No behavior labels or event times are present; only continuous neural traces.
- The recordings are immobilized and stimulus‑free; dynamics reflect internal
  state progression rather than evoked responses.
