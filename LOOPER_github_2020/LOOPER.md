# LOOPER.md

This document describes **how LOOPER works in this repo**, grounded in the
PLOS Comp Bio paper (2023) and the actual MATLAB implementation under
`LOOPER_github_2020/Functions/`.

LOOPER is an unsupervised method for compressing noisy neural population
activity into a **set of interlocking 1‑D trajectories** ("loops") and a
**computational scaffold** that labels each time point by:
- a **trajectory/loop ID** (alpha), and
- a **phase position** along that trajectory (theta/phase bin).

The scaffold is a compact, interpretable representation of population dynamics
and can be used to compare computations across systems, predict behavior, and
simulate trajectories.

---

## Conceptual overview (paper)

- Start from high‑dimensional neural activity time series.
- Use **diffusion maps** (a Markov transition operator) to capture local
  transition structure under noisy dynamics.
- **Coarse‑grain** the transition matrix by clustering states with similar
  transition distributions.
- Assume the system stores information robustly, implying dynamics can be
  approximated by **stable 1‑D trajectories**.
- Identify and cluster these trajectories into a minimal model.
- Use the resulting model to assign a trajectory ID and phase to each time
  point (the **computational scaffold**).

---

## How the MATLAB code implements LOOPER

### 1) Preprocessing (`preprocessData.m`)
Inputs are raw neural data (optionally with input/output channels and trial IDs).
Key behavior:
- Accepts trialized data (cell arrays) and merges into a single stream.
- Optional smoothing.
- Optional z‑scoring (per neuron) using training data statistics.
- Optional delay embedding (`delayEmbed.m`).
- Scales input/output channels by a lambda weight relative to neural magnitude.
- Returns:
  - `tempData` (merged, processed stream),
  - `trialData` (cell array per trial),
  - `processedTrialSwitches` (trial end indices),
  - `dataMean`, `dataSTD` (for re‑use on validation data).

### 2) Local distances and diffusion map (`buildDiffusionMap.m`, `findBestAxes.m`)
LOOPER builds a Markov transition matrix from local similarity between states.

Core mechanics:
- `findBestAxes.m` computes a **local distance profile** for each time point by
  combining:
  - static similarity of activity vectors, and
  - similarity of local derivatives/velocity (cosine distance on a smoothed
    temporal derivative).
- `MinReturnTime` (called `MinimumReturnTime` in code) sets a minimum temporal
  separation when identifying candidate neighbors, which reduces trivial
  “nearest neighbor = immediate timepoint” matches.
- Distances are scaled by a locally estimated sigma and converted to a Gaussian
  kernel with thresholding (exp(-1) neighborhood).
- `buildDiffusionMap.m` assembles these local neighborhoods into a transition
  matrix and performs several cleanups:
  - removes trial terminal states from transitions,
  - optionally adds a terminal state for trial stitching,
  - raises the transition matrix until a target **repopulation density** is
    reached (ensures connectivity),
  - normalizes with a steady‑state density estimate.
- Output is an **asymmetric Markov transition matrix** representing dynamics.
  - Asymmetry comes from centering kernels on the *next* time step (xₜ₊₁),
    so the diffusion map is a transition operator rather than a purely
    geometric neighborhood graph.

### 3) State reduction / coarse graining (`reduceMatrix.m`)
- Treats each row of the transition matrix as a state’s transition distribution.
- Computes pairwise distances between rows (default: correlation distance).
- Performs hierarchical clustering for a range of cluster counts.
- Evaluates each clustering via **multi‑step KL divergence** between true and
  reduced transition dynamics and selects a best model using a BIC‑like score.
- Produces:
  - `clusterIDs` (state → cluster),
  - `reducedMatrix` (cluster‑level transitions).

### 4) Loop/trajectory discovery (`buildMinimalModelFromMatrix.m`)
- Builds transition counts between clusters and infers loop paths:
  - Uses shortest‑path search on a weighted graph derived from transition
    probabilities to find likely closed trajectories.
- Computes loop overlap and a loop‑similarity matrix, then clusters loops into
  a smaller set (either automatically or via GUI interaction).
- This yields a **set of candidate trajectories** that represent stable 1‑D
  dynamics in the data.

### 5) Loop alignment + phase binning (`buildMinimalModelFromLoops.m`)
- Aligns loops to a common phase reference.
- Splits each loop into a fixed number of **phase bins** (BINS_PER_LOOP).
- Assigns each time point to a (loop ID, phase bin) pair.
- Computes **emissions** (mean observed activity) for each bin.

Primary outputs written into the LOOPER data struct (names used in the scripts):
- `BestModel`: transition matrix between scaffold states (loop/phase bins).
- `BestEmission`: emission means per scaffold state.
- `BestStateMap`: per‑timepoint mapping to loop/phase.
- `BestLoopAssignments`: loop cluster membership mapping.

### 6) Validation (`validateModel.m`, `validateData.m`)
- Projects validation data onto the learned scaffold.
- Evaluates a log‑likelihood‑style score based on emissions and transition
  probabilities.
- Produces diagnostic plots in PCA space.

---

## What the "computational scaffold" is (in practice)

From the MATLAB outputs, the scaffold is effectively:
- A **discrete sequence of (loop ID, phase bin)** for every time point.
- A **transition graph** between these bins (BestModel).
- Emission means that map each bin back into neural space (BestEmission).

You can visualize it in two ways:
- Latent space trajectories (trajectories in diffusion‑map or PCA space).
- The **trajectory ID vs time** plot (the scaffold itself).

---

## Key parameters in this repo

These are the most important knobs you’ll see in scripts or the LOOPER app:
- `NearestNeighbors`: number of local neighbors in diffusion map construction.
- `RepopulateDensity`: how far the transition matrix is propagated to ensure
  connectivity (matrix powering step). In code this is implemented by
  repeatedly multiplying the transition matrix until the *max* reachable
  state count reaches a target fraction of all states.
- `MinReturnTime`: separation for repeated peaks when estimating local
  neighbors.
- `DelayTime`, `DelayCount`: delay embedding parameters (see below).
- `DistanceType` (in `reduceMatrix`): default `correlation`.
- `MaxCheckTime`: number of steps used for multi‑step KL evaluation.
- `TotalStates`: target total number of scaffold bins; sets
  `BINS_PER_LOOP = ceil(TotalStates / maxClusters)`, so it directly controls
  **phase resolution per loop**.
- `UseTerminalState`: whether to add a terminal state for trial stitching.
- `PutativeClusterCounts`: candidate cluster counts for coarse graining.
- `PutativeLoopCounts`: candidate loop counts for final model selection.

(See Table 2 of the paper for reference defaults; the OSF scripts often override
these via BayesOpt or grid search.)

### Reference settings from the paper (Table 2)

From Table 2 in the LOOPER paper, the **C. elegans** (continuous) row uses:
- Neighbor count = **8**
- Use local dims = **T**
- Repopulation density = **0.95**
- Min return time = **10**
- Distance measure = **correlation**
- Max check time = **10**
- Total state count = **25**
- Use terminal state = **F**

Note: `UseLocalDimensions` defaults to **false** in `LOOPER.m`, but the paper
used **true** for C. elegans and most other datasets. If you want to match the
paper closely, consider setting it explicitly.

### Delay embedding in this codebase (important detail)

`preprocessData.m` applies `delayEmbed` when `DelayTime > 0` and
`DelayCount > 0`. In this repo it calls:

```
delayEmbed(allTrials{i}, delayCount, delayTime, 0, true)
```

So **derivatives are NOT included** (useDerivatives = 0), and the first
`delayCount * delayTime` samples of each trial are dropped. This changes the
effective trial length and should be accounted for when aligning events.

### Parameter selection in practice (paper + OSF scripts)

Paper guidance:
- LOOPER has three core parameters: **NearestNeighbors**, **RepopulateDensity**,
  and **TotalStates** (bin count). Only **NearestNeighbors** has a strong effect
  on performance in their parameter sweeps; it should be roughly the number of
  *similar trials* in the dataset. Other parameters are comparatively robust.

OSF scripts in this repo show two practical strategies:
- **BayesOpt** (e.g., `ComapreToMarkovFitLOOPER.m`): optimizes
  `NearestNeighbors`, `RepopulateDensity`, and `TotalStates` by maximizing
  reconstruction correlation on held‑out trials.
- **Grid search** (e.g., `socialMoviesLOOPER.m`): sweeps over
  `DelayCount`, `Smoothing`, `NearestNeighbors`, and a single
  `PutativeClusterCounts`, then evaluates loop assignment consistency on
  validation subjects.

Takeaway: default parameters are acceptable for a first run, but for quantitative
comparisons the original pipelines **tune NN (and sometimes delay/smoothing)**
against a reconstruction/validation objective.

---

## Practical constraints / caveats

- **Heuristic clustering**: the loop reduction step is heuristic and may require
  manual intervention (GUI) for best results.
- **Parameter sensitivity**: nearest‑neighbor count, repopulation density, and
  cluster counts affect scaffold structure.
- **Missing neuron identity**: cross‑worm comparisons require careful handling
  of neuron identity and intersection sets.
- **Trial stitching**: if trials are present, transitions across trial boundaries
  are handled explicitly; be sure trial indices are correct.
- **Not a full dynamical system model**: LOOPER captures dominant 1‑D
  trajectories and their branching, but does not perform full stability
  analysis or parameterized dynamical modeling.

---

## Where to look in code

Core functions:
- `LOOPER_github_2020/Functions/preprocessData.m`
- `LOOPER_github_2020/Functions/buildDiffusionMap.m`
- `LOOPER_github_2020/Functions/findBestAxes.m`
- `LOOPER_github_2020/Functions/reduceMatrix.m`
- `LOOPER_github_2020/Functions/buildMinimalModelFromMatrix.m`
- `LOOPER_github_2020/Functions/buildMinimalModelFromLoops.m`
- `LOOPER_github_2020/Functions/validateModel.m`
- `LOOPER_github_2020/Functions/validateData.m`

Figure/analysis pipelines (OSF scripts) live in:
- `LOOPER_OSF_hosted_2022/`

The GUI app is:
- `LOOPER_github_2020/LOOPER.mlapp`
