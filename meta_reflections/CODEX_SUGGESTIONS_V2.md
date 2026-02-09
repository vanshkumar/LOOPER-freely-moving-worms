# CODEX_SUGGESTIONS (Top-Level, Timestamp-Grounded Revision)

Updated: 2026-02-09

## Scope Reviewed (with chat coverage details)

Priority order requested and followed:
1. `LOOPER-freely-moving-worms/CODEX_SUGGESTIONS.md`
2. `daily_logs/*.md`
3. `LOOPER-freely-moving-worms/` docs + paper (`journal.pcbi.1010784.pdf`)
4. `codex_chats/`
5. `chat_trajectory/`
6. other repos (primarily markdown)

For the chat corpora specifically:
- `codex_chats` raw rollouts scanned: 63 files (`rollout-*.jsonl`).
- `codex_chats/clean_llm` analyzed end-to-end for deduplicated behavioral signals: 21 files.
- `chat_trajectory/minimal_turns.jsonl` analyzed end-to-end: 35 conversations, 223 user turns.
- Combined non-boilerplate user turns used for timestamped pattern extraction: 1,095.

Why `clean_llm` is primary for behavior analysis:
- raw rollouts contain repeated snapshots/instruction boilerplate;
- `clean_llm` preserves chronology while deduping repeated policy wrappers, giving more faithful process signals.

## What the Timestamped Chats Show About Your Thinking Process

### 1) You repeatedly enforce an evidence gate before accepting interpretation

Representative turns:
- `2025-12-07T18:00:37Z`: “does this dataset contain enough data ... it feels underpowered?”
- `2025-12-13T16:38:26Z`: asks what each option means for what z-score “actually compute[s]”.
- `2026-01-06T01:34:29Z`: “explain it first before trying to fix”.

Interpretation:
- You do not accept narrative convenience without identifiability or mechanism checks.
- This is a core strength and should be formalized earlier in each branch.

### 2) You switch between abstraction and implementation rapidly (high leverage, high risk)

Representative turns:
- `2025-12-13T15:06:56Z`: asks for a simple conceptual statement of the analysis goal.
- `2026-01-23T14:50:40Z`: immediately surfaces freely-moving vs immobilized confound.
- `2026-01-27T16:19:32Z`: proposes a concrete missing metric (“phase continuity”).

Interpretation:
- You naturally bridge theory and concrete tests.
- Risk appears when this happens in long uninterrupted sessions: constraints shift mid-stream without explicit re-charter.

### 3) You actively orchestrate workflow and ordering, not just ask questions

Representative turns:
- `2026-01-09T14:11:10Z`: constrain eval to complex eigenvalues (exclude real-only).
- `2026-01-10T14:51:04Z`: remove spectral-gap item when it no longer matches the active framing.
- `2026-01-30T01:07:24Z`: explicit rerun instructions across Kato and Atanas before continuing.

Interpretation:
- You have strong steering behavior once a hypothesis/test frame is clear.
- Agent performance improves materially when this steering is externalized as fixed checklists.

### 4) You push continuity through uncertainty, then revisit branch validity

Representative turns:
- `2026-01-27T18:22:24Z`: considers pivot to negative-result framing.
- `2026-01-28T14:15:11Z`: “nvm keep going.”
- `2026-01-30T00:58:59Z`: “keep going” (again, after additional checks).

Interpretation:
- Strong persistence is present.
- The failure mode is local optimization in a stale objective if no checkpoint forces an explicit branch decision.

### 5) You treat documentation as part of scientific control, not post-hoc polish

Representative turns:
- `2025-12-13T14:42:29Z`: asks for a caveat-aware `METHODS.md`.
- `2025-12-13T19:11:19Z`: asks for exact coding-agent instructions.
- `2026-01-30` doc consistency passes across `LOOPER-freely-moving-worms`.

Interpretation:
- This is a major advantage: you already produce traceable reasoning artifacts.
- Missing piece: stronger coupling from docs to explicit claim status and decision checkpoints.

## Updated Diagnosis

Compared to the previous version, the larger timestamped corpus reinforces and sharpens the core diagnosis:

- Your bottleneck is not coding speed; it is objective stability across long multi-goal sessions.
- You think best in cycles of:
  - conceptual disambiguation,
  - concrete test design,
  - implementation,
  - audit/correction.
- Drift happens when those modes are mixed continuously without hard handoff markers.

## Updated Operating Recommendations (tailored to observed trajectory)

### 1) Use a mandatory `Question Gate` before any code branch

One-page pre-commit artifact:
- `Question`
- `Claim`
- `Data slice`
- `Metric`
- `Comparator/control`
- `Falsifier`
- `Stop condition`

No code before this exists.

### 2) Split sessions into explicit modes

Three-thread protocol:
- `SPEC`: question, identifiability, controls, confounds.
- `BUILD`: implementation only against frozen spec.
- `AUDIT`: challenge assumptions/interpretation and decide keep/pivot/kill.

Do not allow `BUILD` to redefine the question.

### 3) Enforce turn-budget handoffs

Hard rule:
- after 60-100 turns (or objective change), produce `SESSION_HANDOFF.md` and restart thread.

Required handoff fields:
- current claim
- strongest counter-explanation
- latest decisive result
- next discriminating test
- reason for continuing vs killing branch

### 4) Add a timestamped decision ledger

Create `docs/CLAIMS_LEDGER.md` with rows:
- `timestamp`
- `decision` (continue / pivot / kill)
- `trigger evidence`
- `confounds still open`
- `next test`

This converts your natural “keep going vs pivot” moments into auditable scientific state transitions.

### 5) Add a negative-result protocol

When uncertainty persists after planned checks:
- write explicit “negative-result candidate” note,
- list what was ruled out,
- list what remains unidentifiable,
- specify the new data required to resolve it.

This prevents losing signal from branches that do not confirm the initial story.

### 6) Keep metric governance strict

A metric enters `RESULTS.md` only if all pass:
- synthetic sanity,
- positive control,
- null/surrogate,
- plain-language definition in `METHODS.md`.

### 7) Keep dataset identifiability explicit

Per dataset-condition pair:
- regime (`immobilized` / `freely moving` / `stimulated`)
- what can be identified
- what cannot be identified
- likely confounds that mimic target behavior

This specifically addresses repeated freely-moving vs immobilized interpretation pressure.

## Suggested Docs To Add/Enforce

- `docs/QUESTION_GATE.md`
- `docs/DATASET_IDENTIFIABILITY.md`
- `docs/CLAIMS_LEDGER.md`
- `docs/METRIC_GOVERNANCE.md`
- `templates/SESSION_HANDOFF.md`

## Bottom Line

The expanded, timestamped chat record confirms a strong research pattern:
- high execution velocity,
- high skepticism,
- strong orchestration instincts.

The single biggest improvement lever is not “think harder” or “code faster”; it is **forcing explicit mode transitions and decision checkpoints so your strongest reasoning habits are captured before long-session drift appears**.

## Appendix: Weekly Timeline and Top Pattern Shifts

### Weekly timeline (timestamped chat patterns)

| Week | Build/Execute | Workflow | Operationalization | Evidence | Debug | Abstraction | Continuity |
|---|---:|---:|---:|---:|---:|---:|---:|
| 2025-W49 | 0 | 4 | 4 | 4 | 0 | 2 | 0 |
| 2025-W50 | 9 | 10 | 16 | 6 | 4 | 2 | 1 |
| 2025-W51 | 1 | 5 | 4 | 5 | 0 | 0 | 0 |
| 2025-W52 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 2026-W01 | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| 2026-W02 | 11 | 18 | 15 | 6 | 5 | 1 | 0 |
| 2026-W03 | 11 | 9 | 4 | 5 | 2 | 2 | 1 |
| 2026-W04 | 59 | 22 | 12 | 17 | 20 | 2 | 0 |
| 2026-W05 | 92 | 72 | 52 | 44 | 30 | 11 | 7 |

### Top 3 shifts for each week-to-week transition

- 2025-W49 -> 2025-W50: operationalization +12; build execute +9; workflow orchestration +6
- 2025-W50 -> 2025-W51: operationalization -12; build execute -8; workflow orchestration -5
- 2025-W51 -> 2025-W52: workflow orchestration -5; evidence seeking -5; operationalization -4
- 2025-W52 -> 2026-W01: evidence seeking +1; workflow orchestration +0; operationalization +0
- 2026-W01 -> 2026-W02: workflow orchestration +18; operationalization +15; build execute +11
- 2026-W02 -> 2026-W03: operationalization -11; workflow orchestration -9; debug fix -3
- 2026-W03 -> 2026-W04: build execute +48; debug fix +18; workflow orchestration +13
- 2026-W04 -> 2026-W05: workflow orchestration +50; operationalization +40; build execute +33
