# LET 'EM LOOSE: Where and How to Make LLMs More Autonomous

Prepared after reading: both CLAUDE_SUGGESTIONS.md files, CODEX_SUGGESTIONS.md, all 31 daily logs + Reflections, 21 Codex session logs, 35 ChatGPT conversation threads (~3,900 turns), AGENTS.md/METHODS.md/GOALS.md/DATASET.md across all repos, and both LEARNINGS.md files.

---

## The core tension

Your project revealed a specific version of a general problem: **LLMs are fast where you need them slow (idea generation, experiment design) and slow where you need them fast (code execution, data inspection, documentation)**. The Codex sessions show an agent that couldn't run MATLAB, creating a 20-turn debugging proxy loop for something that should have been 3 turns. The ChatGPT sessions show 300+ turn conversations that drifted because the model couldn't inspect the data it was reasoning about. Meanwhile, you spent your most valuable cognitive resource -- uninterrupted thinking time -- on tasks that should have been delegated.

The fix is not "more autonomy everywhere." It is **autonomy in the right places, with the right constraints, and mandatory human checkpoints where your project actually broke down**.

This document maps exactly where those boundaries should be, and gives you specific tools to enforce them.

---

## Part 1: The Autonomy Map

Based on what actually happened across 21 Codex sessions, 35 ChatGPT conversations, and 31 daily logs:

### Let them run (high autonomy, low risk)

| Task | Why it's safe | Evidence from your project |
|------|--------------|---------------------------|
| **Read & summarize codebases/papers** | Bounded output, easy to verify, 95% success rate | Every Codex session started with autonomous repo exploration. Consistently excellent. |
| **Write/update documentation** | Verifiable against source, low cost of errors | The Jan 30 "ya x15" session -- you rubber-stamped 15 consecutive doc edits because they were all correct. |
| **Code execution & debugging** | Errors are caught by the compiler/runtime, not by your judgment | The #1 bottleneck: you spent hours as a MATLAB execution proxy. An agent that can run code collapses 20-turn loops into 3. |
| **Data structure inspection** | Factual, verifiable | Jan 7: Codex loaded 1,425 neurons for a worm with ~107-131. You caught it by domain knowledge. But the *inspection* itself should be autonomous -- the agent should run `disp(fieldnames(S))` itself. |
| **Null model / control experiments** | Well-specified, bounded, the results are pass/fail | Your Priority 1 experiment (time-shuffled LOOPER) is a perfect autonomy target: clear specification, clear success criterion, no scientific judgment needed in execution. |
| **Code refactoring & cleanup** | Low scientific risk, high tedium for you | Dec 13: Codex extracted helpers into 8 .m files autonomously. Correct. |
| **Literature search & paper retrieval** | Breadth is the LLM's advantage; you verify relevance | ChatGPT's paper discovery was consistently useful. The failure was when you relied on its *interpretation* without reading the paper yourself. |
| **Cross-file consistency checks** | Mechanical, tedious for humans | Checking that METHODS.md matches RUNBOOK.md matches actual script arguments. |

### Keep the human in the loop (low autonomy, high risk)

| Task | Why it's dangerous to automate | Evidence from your project |
|------|-------------------------------|---------------------------|
| **Question formulation** | This is where your project actually broke down. Three times. | Heat pulse experiment: question didn't match data. Split-half stationarity: metric didn't match claim. QSimeon Stim vs NoStim: dataset didn't match question. All three were execution-before-understanding failures. |
| **Choosing what to build next** | "The danger of LLMs in research is that they make the wrong questions feel productive." | Day 24: you rejected ChatGPT's experiment suggestion because it "overindexed on speculation in the Moore paper." Your judgment was correct. The LLM's was plausible but wrong. |
| **Scientific interpretation** | The LLM cannot distinguish "surprising" from "expected" in your specific context | ChatGPT reframed your stationarity failure as a "stress test" -- a correct and important reframing. But *you* had to recognize it was correct. The LLM can generate interpretations; only you can evaluate them. |
| **Strategic pivots** | These require integrating your energy, timeline, skill gaps, and intellectual taste | Day 12: pivoting from mouse PFC to C. elegans. Day 21: pivoting from ADM to LOOPER. Day 24: rejecting the ChatGPT experiment. All were correct calls that required your judgment, not the LLM's. |
| **Email/communication to experts** | Tone, positioning, and intellectual honesty are high-stakes | The Brennan email was your sharpest output. ChatGPT helped calibrate tone, but you wrote the substance. |

### The gray zone (moderate autonomy, with guardrails)

| Task | The right level of autonomy | How to guardrail it |
|------|----------------------------|-------------------|
| **Experiment design (implementation)** | LLM writes the experiment contract; you approve before code | The experiment contract template (below) is the guardrail. |
| **Result summarization** | LLM writes first draft of RESULTS.md updates; you edit | Require the LLM to flag caveats and alternative interpretations alongside any positive claim. |
| **Dataset evaluation** | LLM runs basic stats and reports; you decide suitability | The dataset gate (below) structures what it reports. |
| **Code review** | LLM reviews for correctness; you review for scientific validity | Day 29: ChatGPT web found an indexing bug Codex missed. Use two models, triangulate. |

---

## Part 2: Specific Skills to Build

These are concrete Claude Code / Codex skills (or prompt templates) you can create. Each one operationalizes a lesson from the project.

### Skill 1: `/question-gate`

**What it does:** Before any new experiment or analysis direction, forces you through the three-gate system from CODEX_SUGGESTIONS.md. The LLM fills in what it can; you fill in what it can't.

**Template it generates:**

```markdown
## QUESTION GATE
- Question (one sentence, falsifiable): ___
- What would change my mind: ___
- What I actually believe right now (written by me, not the LLM): ___

## DATASET GATE
- Dataset: ___
- Recording regime (immobilized/freely moving/etc): ___
- Trial structure: ___
- Known nonstationarities: ___
- Missing metadata that matters: ___
- Can this dataset answer this question? (explicit yes/no with reasoning): ___

## METRIC GATE
- Metric: ___
- What it measures (in plain English): ___
- What a null model produces for this metric: ___
- What the positive control produces: ___
- What would constitute a meaningful effect size: ___

## KILL CONDITION
- Stop if: ___
```

**Why this skill specifically:** Your three most expensive mistakes (heat pulse, stationarity, QSimeon Stim/NoStim) would each have been caught by one of these gates. The skill makes the gates frictionless -- you type `/question-gate` and get the template. The 1-2 hours filling it in saves the multi-day pivots.

**Implementation:** A Claude Code skill that outputs this template pre-filled with whatever the LLM knows about the current repo context (dataset names, available metrics, etc.), with blanks for the parts that require your judgment.

### Skill 2: `/experiment-contract`

**What it does:** Generates a pre-registration document for a specific experiment. More detailed than the question gate -- this is the implementation spec.

**Template:**

```markdown
## EXPERIMENT CONTRACT: [name]
Date: ___
Status: proposed / approved / running / completed / killed

### Specification
- QUESTION: [one sentence, falsifiable]
- DATA: [specific files, fields, time windows]
- METHOD: [analysis, parameters, comparisons]
- EXPECTED IF TRUE: [specific numbers/patterns]
- EXPECTED IF FALSE: [specific numbers/patterns]
- ASSUMPTIONS: [list everything this assumes about the data]

### Pre-flight checks
- [ ] Paper methods section read for each dataset used
- [ ] Null model defined and ready to run
- [ ] Positive control identified
- [ ] Kill condition written

### Results (filled in after running)
- Outcome: ___
- Matches expected? ___
- Caveats: ___
- Next step: ___
```

**Why this skill:** Your EXPERIMENTS.md was strong on framing but weak on pre-registered expected results. This closes that gap. The LLM can pre-fill DATA, METHOD, and ASSUMPTIONS from the repo context; you fill in QUESTION, EXPECTED IF TRUE/FALSE, and KILL CONDITION.

### Skill 3: `/audit`

**What it does:** Opens a fresh context (new session/conversation) whose only job is to attack a result. The prompt is hard-coded:

> "Here are my results and my interpretation. Your job is to find every way this could be wrong, misleading, or uninterpretable. Do not suggest what to do next -- only identify problems."

**Why this skill:** You never had a consistent auditor. ChatGPT was your thinking partner (not adversarial enough). Codex was your builder (wrong role entirely). The audit skill enforces the Architect/Builder/Auditor separation from CODEX_SUGGESTIONS.md by making the auditor a separate session with a separate prompt.

**Implementation:** Reads the current RESULTS.md and the relevant experiment contract, then launches a new Claude context with the adversarial prompt. Returns a structured list of vulnerabilities.

### Skill 4: `/state-update`

**What it does:** Updates STATE.md (see Part 3 below) with the current state of your understanding. Prompts you for the parts the LLM can't fill in.

**Template it generates:**

```markdown
## STATE.md — Updated [date]

### Current hypothesis
[LLM pre-fills from latest RESULTS.md and experiment contracts]

### Evidence for
- [bullet list, auto-extracted from completed experiments]

### Evidence against
- [bullet list]

### Open questions
- [carried forward from last STATE.md + any new ones]

### What I actually believe (human-written)
___

### Next experiment
___

### Biggest risk to current direction
___
```

**Why this skill:** Your daily logs were chronological (good for reflection, bad for re-entry). STATE.md is a snapshot optimized for future-you and for LLM agents who need to understand where the project stands. The `/state-update` skill makes maintaining it frictionless.

### Skill 5: `/bridge`

**What it does:** Generates a handoff document when switching from thinking (ChatGPT) to building (Codex/Claude Code).

**Template:**

```markdown
## BRIDGE: [thinking tool] → [building tool]
Date: ___

### What was decided (from the thinking session)
- Question: ___
- Approach: ___
- Key constraints: ___

### What the builder needs to know
- Files to read first: ___
- Data to use: ___
- Expected output format: ___
- Known gotchas: ___

### What the builder should NOT do
- Do not choose a different question
- Do not add features beyond spec
- Do not use production-grade patterns (this is a research prototype)

### Success criterion
___
```

**Why this skill:** Your project showed a consistent pattern: insights developed in ChatGPT conversations didn't fully transfer to Codex sessions. The bridge document prevents the builder from re-deriving (or silently modifying) the architect's decisions.

### Skill 6: `/thinking-buffer`

**What it does:** After any major experiment result, blocks further coding for 20 minutes. Outputs a prompt:

> "Do not open any LLM. Blank page or blank screen. Write:
> 1. What does this result actually say?
> 2. What does it NOT say?
> 3. What is the cheapest next experiment that would disambiguate?
> 4. Am I still asking the right question?"

**Why this skill:** The 20 minutes of doing nothing that cracked your understanding of the Brennan paper was the highest-ROI activity of the entire project. This skill operationalizes it. It's a speed bump, not a roadblock -- but it prevents the execution-before-understanding failure mode that cost you the most time.

**Implementation:** Could be as simple as a skill that prints the prompt and sets a timer. The point is making it a habit, not a technology.

---

## Part 3: Documentation Architecture for Autonomous Agents

Your project already developed strong documentation patterns. Here's what to keep, what to add, and what to change.

### Keep (these worked)

1. **The "mission brief" AGENTS.md pattern** (from category-learning_mpfc). A single document stating: hypothesis, core analytical idea, what a positive result looks like. No process prescription. This is the most autonomy-enabling pattern because it constrains the *question* without constraining the *method*. The agent has freedom to choose implementation while staying aligned to the scientific goal.

2. **Tiered deliverables** (from AsymmetricDiffusionMapping GOALS.md). MVP / Good / Stretch creates natural checkpoints without requiring human intervention at each stage. The agent can self-assess whether it's reached each tier and decide whether to advance. This is an implicit autonomy dial.

3. **Anticipatory gotcha documentation** (from category-learning_mpfc METHODS.md). The 10 numbered caveats/gotchas are guard rails. The agent has freedom of movement within the lane but cannot silently drive off a cliff. Continue writing these for every new dataset.

4. **Disambiguating code examples** (from category-learning_mpfc METHODS.md). Short MATLAB snippets showing how to interpret ambiguous field names. These function as executable specifications. Continue embedding these in METHODS.md and DATASET.md.

5. **LEARNINGS.md** (from LOOPER). The 18 entries captured real friction. Continue, but restructure (see below).

### Add

1. **STATE.md** at the top level. Updated weekly via `/state-update`. This replaces daily logs as the primary LLM-readable context document. Daily logs are still valuable for your own reflection, but STATE.md is what the agent reads at session start.

2. **CLAIMS_LEDGER.md** (from CODEX_SUGGESTIONS.md). One row per claim:

   ```
   | claim_id | claim | status | evidence | controls_passed | confounds |
   ```

   This prevents the narrative drift that happens when RESULTS.md is updated incrementally. The ledger is the source of truth for what you've actually shown.

3. **Bridge documents** (via `/bridge` skill). One per thinking-to-building transition. Stored in a `bridges/` directory so agents can reference the chain of decisions.

4. **Experiment contracts** (via `/experiment-contract` skill). One per experiment. Stored in `experiments/` directory. The contract is written before code; the results section is filled in after.

### Change

1. **Restructure LEARNINGS.md by category, not chronology.** Your current LEARNINGS.md is a chronological list. An agent scanning for "caching issues" has to read all 18 entries. Restructure into sections:

   ```markdown
   ## Environment constraints
   - No pip/network access in sandbox
   - apply_patch must use dedicated tool, not shell

   ## MATLAB gotchas
   - addpath doesn't recurse; use genpath()
   - LOOPER overwrites saveData; reattach custom fields after

   ## Methodological guardrails
   - Detrend must fit on training window only
   - Split-half stationarity is a stress test, not a validation

   ## Rules (from human corrections)
   - ALWAYS explain what you think the issue is before implementing a fix
   - Keep code simple -- this is a research prototype, not production
   - Do not add optional name/value pair arguments
   ```

   The "Rules" section is critical. Your most important meta-cognitive instruction ("explain the issue before fixing") was in CLAUDE.md, not LEARNINGS.md. It should be in both.

2. **Add AGENTS.md to every repo.** Only category-learning_mpfc had one. AsymmetricDiffusionMapping had the documentation but not the explicit "you are the agent, here is your mission" framing. The mission brief pattern is what makes high-autonomy operation possible.

3. **Add explicit "Do NOT" sections to AGENTS.md.** Your Codex sessions showed recurring over-engineering (production-grade input validation, optional parameters, excessive error handling). A "Do NOT" section prevents this:

   ```markdown
   ## Do NOT
   - Add optional name/value pair arguments
   - Add input validation beyond basic type checks
   - Write production-grade error handling
   - Choose which experiment to run (ask the human)
   - Modify the scientific question without explicit approval
   ```

---

## Part 4: Process Architecture -- Keeping Directional Control

The failure modes from your project were not failures of execution. They were failures of direction. The LLM built what you asked for, quickly and competently. The problem was that what you asked for was sometimes the wrong thing. Here's how to prevent that without slowing down the parts that should be fast.

### The mandatory human checkpoints

These are non-negotiable speed bumps. They exist because your project showed that *removing* them leads to multi-day pivots.

1. **Before any new experiment: `/question-gate`** (15-30 min)
   - The LLM pre-fills what it can. You fill in: the question, the kill condition, what you actually believe.
   - If you can't fill in "what I actually believe" without consulting an LLM, you're not ready to run the experiment. Sit with it for 20 minutes first.

2. **After any major result: `/thinking-buffer`** (20 min, no LLM)
   - What does this result say? What doesn't it say?
   - This is where the "20 minutes of doing nothing" insight lives. Make it a habit, not an accident.

3. **Before switching from thinking to building: `/bridge`** (10 min)
   - Prevents the builder from re-deriving or silently changing the architect's decisions.

4. **Weekly: `/state-update`** (20 min)
   - Forces you to write what you actually believe, distinct from what the LLM has told you.
   - The "What I actually believe (human-written)" field is the critical one.

### The autonomous loops (let them run)

Between the checkpoints, the LLM should operate with high autonomy. Specifically:

**Loop 1: Build and debug**
After you approve an experiment contract, the agent:
1. Reads the contract, METHODS.md, DATASET.md, LEARNINGS.md
2. Writes the code
3. Runs it (this requires an execution environment -- see Part 5)
4. Debugs until it produces output
5. Updates RESULTS.md with raw output + caveats
6. Notifies you that results are ready for interpretation

You do NOT need to be in the loop for steps 1-5. You enter at step 6 for interpretation.

**Loop 2: Null model execution**
After you define a null model in the experiment contract, the agent:
1. Implements the null (time-shuffled, phase-randomized, etc.)
2. Runs it N times
3. Computes summary statistics
4. Compares to the real data
5. Reports: "Real data metric = X. Null distribution: mean = Y, 95% CI = [Z1, Z2]. Real data is [inside/outside] the null."

This entire loop is autonomous. You check the output.

**Loop 3: Documentation maintenance**
After any code change, the agent:
1. Checks if METHODS.md, RUNBOOK.md, or EXPERIMENTS.md need updating
2. Updates them
3. Runs consistency checks (does RUNBOOK reference scripts that exist? do METHODS parameters match script defaults?)

This is the "ya x15" pattern. Let it run.

**Loop 4: Literature retrieval**
When you identify a paper to read, the agent:
1. Finds and downloads the PDF
2. Extracts key methods, datasets, results, and limitations
3. Writes a structured summary
4. Flags connections to your current hypothesis (from STATE.md)
5. You read the paper yourself -- but the summary gives you scaffolding

The key difference from your project: in the Aliveline project, you sometimes read the LLM's summary *instead of* the paper. The summary should be a *companion* to reading, not a *replacement*.

### The conversation hygiene rules

From your ChatGPT usage patterns:

1. **New conversation every ~50 turns or when the topic shifts.** Your 300+ turn conversations accumulated stale context. The "take a step back, what are we doing here?" signal at turn 293 is diagnostic -- that reset should have been a new conversation.

2. **One question per conversation.** Not "explore this general area" but "is X true, and what evidence would distinguish it from Y?" Your shortest, most productive ChatGPT conversations had bounded questions. Your longest, least productive ones were open-ended explorations.

3. **Write the question before opening the conversation.** This is the "think before prompting" rule operationalized. If you can't write the question in one sentence, you need the thinking buffer, not ChatGPT.

---

## Part 5: Infrastructure for Actually Doing This

### The execution environment gap

The single biggest thing limiting your agents' autonomy was **inability to execute code**. The Codex sessions made this painfully clear: you were a human MATLAB proxy for 20-turn debugging loops.

**Options (in order of pragmatism):**

1. **Claude Code with MATLAB access.** If you can run Claude Code in an environment where MATLAB is available (or where MATLAB scripts can be executed via a shell command), the debugging proxy problem disappears entirely. Claude Code already has Bash access. The missing piece is MATLAB on the path.

2. **Codex with a MATLAB sandbox.** If using Codex, configure the sandbox to include MATLAB (or Octave as a free substitute for basic linear algebra / matrix operations). Many LOOPER operations would work in Octave.

3. **Python translation layer.** For some analyses (eigenvalue decomposition, covariance analysis, correlation metrics), a Python implementation might be faster to iterate on than MATLAB. The tradeoff is fidelity to the original LOOPER codebase.

4. **At minimum: a "run and report" skill.** A skill that takes a script path, runs it in MATLAB, captures stdout/stderr and any output files, and returns them to the agent. This collapses the proxy loop without requiring the agent to have interactive MATLAB access.

### The unified documentation pipeline

From your Reflections: "moving forward make sure all documentation (chats, codex sessions, markdown notes) is in a single easy to access place, with timestamps, so llms can read and directly reflect on it to improve iteratively."

**Concrete implementation:**

```
project/
  STATE.md                    # weekly snapshot (via /state-update)
  CLAIMS_LEDGER.md            # evidence tracking
  experiments/                # experiment contracts
    001_null_model.md
    002_kato_stationarity.md
  bridges/                    # thinking-to-building handoffs
    2026-02-10_chatgpt_to_claude.md
  daily_logs/                 # your reflective journal (unchanged)
  sessions/                   # auto-archived LLM session summaries
    2026-02-10_claude_session_1.md
    2026-02-10_codex_session_1.md
  LOOPER-freely-moving-worms/ # the actual code repo
    AGENTS.md
    METHODS.md
    DATASET.md
    RUNBOOK.md
    RESULTS.md
    LEARNINGS.md
```

The `sessions/` directory is new. After each LLM coding session, the agent writes a brief summary: what was attempted, what worked, what failed, what's next. This is the institutional memory that LEARNINGS.md doesn't capture (LEARNINGS captures gotchas; session summaries capture trajectory).

### The Claude Code harness pattern

From your Reflections: "claude code as harness for codex?"

This is the right instinct. The idea: Claude Code orchestrates, Codex (or another agent) executes. Specifically:

1. Claude Code reads STATE.md, the latest experiment contract, and LEARNINGS.md
2. Claude Code generates a specific, bounded task for the builder agent
3. The builder agent executes (writes code, runs it, reports results)
4. Claude Code reviews the output, updates RESULTS.md, flags anything that needs your attention
5. You review the flagged items

This gives you a two-layer system: Claude Code as the "senior researcher" (reads context, delegates, reviews) and the builder agent as the "research assistant" (writes and runs code). You are the PI who approves questions and interprets results.

**What this looks like in practice:**

```
You: /question-gate
  → You fill in the gates, approve the experiment

You: /bridge
  → Claude Code generates the handoff document

You: "run this experiment"
  → Claude Code reads the bridge + contract
  → Claude Code delegates code writing to a builder agent
  → Builder writes code, runs it, returns output
  → Claude Code reviews output, updates RESULTS.md
  → Claude Code: "Results ready. Recon_corr = 0.62 for real data, 0.31 ± 0.04 for null.
     Flagging: the null distribution has a long right tail. See RESULTS.md for details."

You: /thinking-buffer
  → 20 minutes, no LLM
  → You write what you think this means

You: /audit
  → Fresh Claude session attacks the result
  → Returns: "Three concerns: (1) ... (2) ... (3) ..."

You: /state-update
  → STATE.md updated with new evidence
```

Total human time in this loop: ~1 hour (question gate + thinking buffer + state update). Total agent time: however long the experiment takes to run. The ratio of your cognitive engagement to agent execution is high where it matters (question formulation, interpretation) and low where it doesn't (code, debugging, documentation).

---

## Part 6: Avoiding the Aliveline Failure Modes

Each specific failure mode from the project, mapped to the mechanism that prevents it:

| Failure mode | What happened | Prevention mechanism |
|---|---|---|
| **Question didn't match data** (heat pulse) | Designed experiment before reading Atanas methods section carefully | `/question-gate` dataset gate: "Can this dataset answer this question?" with explicit assumption listing |
| **Metric didn't match claim** (split-half stationarity) | Applied a strict test without understanding it was stricter than the original paper's validation | `/question-gate` metric gate: "What does the positive control produce for this metric?" |
| **Dataset mismatch** (QSimeon WT_Stim vs NoStim) | Assumed dataset compatibility without checking | `/question-gate` dataset gate: explicit recording regime and trial structure fields |
| **Execution before understanding** (general) | Code advanced faster than understanding | `/thinking-buffer` after every major result; `/bridge` before every build |
| **ChatGPT leading instead of stress-testing** | Open-ended brainstorming with ChatGPT produced plausible but untested ideas | Conversation hygiene: write the question before opening the conversation; use ChatGPT to *attack* your hypothesis, not generate one |
| **Builder choosing what to build** (Codex over-engineering) | Codex added production-grade features to a research prototype | AGENTS.md "Do NOT" section; `/bridge` with explicit constraints |
| **Context fragmentation across tools** | Insights from ChatGPT didn't transfer to Codex sessions | `/bridge` documents; STATE.md as single source of truth; `sessions/` directory for institutional memory |
| **No auditor** | Results were generated and interpreted in the same session/mindset | `/audit` skill: separate session, adversarial prompt, no "what to do next" |
| **Expert contact too late** (Brennan email on Day 31) | Activation energy too high, felt "not ready" | Week 1 rule: write the email as a thinking tool, then send it. The email IS the question gate for the human expert. |
| **Momentum loss after break** (holiday gap) | No re-entry document, had to re-bootstrap understanding | STATE.md + `/state-update` before break = instant re-entry. The agent reads STATE.md and knows where you are. |
| **Over-reliance on LLM for idea generation** | "LLM idea generation feels productive while being a substitute for sitting with confusion" | `/thinking-buffer` is the structural fix. The question gate's "What I actually believe (human-written)" field is the diagnostic: if you can't fill it in, the LLM has been doing your thinking. |

---

## Part 7: The Meta-System -- How to Iterate on This

This document is itself a hypothesis. It might be wrong in places. Here's how to test and iterate:

1. **Try it on one experiment first.** Run your Priority 1 (null model controls) through the full loop: question gate, experiment contract, bridge, autonomous execution, thinking buffer, audit, state update. Note where the process helps and where it creates friction.

2. **Log friction in LEARNINGS.md.** If a skill generates a template that doesn't match your thinking, update the template. If a checkpoint feels unnecessary for a particular task, note why -- it might reveal a category of task that doesn't need it.

3. **The 3-session rule.** After 3 coding sessions using this system, re-read this document and LET_EM_LOOSE.md. What's working? What's ceremony? Cut the ceremony. Keep what actually catches errors or saves time.

4. **The "scaffold removal" question.** From your Reflections: "how to remove the LLM scaffold as you get deeper in the field." The answer is: as your domain knowledge grows, the gates get faster (you can fill them in from memory instead of from papers), the thinking buffer gets shorter (you process results faster), and the bridge documents get simpler (you know what the builder needs). The system doesn't go away -- it gets lighter. The checkpoints that still feel necessary after 6 months are the ones that protect against fundamental cognitive biases, not just domain novice mistakes.

5. **Ask yourself monthly:** "Am I generating my own questions, or is the LLM generating them for me?" If the answer is the latter, increase thinking buffer time and decrease LLM brainstorming. If the answer is the former, you can safely increase agent autonomy because you're providing the direction.

---

## Summary: The One-Page Version

**More autonomous:** Code execution, debugging, null models, documentation, data inspection, literature retrieval, consistency checks. Let them run between checkpoints.

**Less autonomous:** Question formulation, experiment design decisions, scientific interpretation, strategic pivots, what to build next. These require your judgment and should be preceded by uninterrupted thinking time.

**The 6 skills to build:** `/question-gate`, `/experiment-contract`, `/audit`, `/state-update`, `/bridge`, `/thinking-buffer`

**The 3 documents to add:** STATE.md, CLAIMS_LEDGER.md, session summaries in `sessions/`

**The 3 changes to existing docs:** Restructure LEARNINGS.md by category, add AGENTS.md with "Do NOT" sections to every repo, add experiment contracts to an `experiments/` directory.

**The infrastructure fix:** Give the agent MATLAB execution access. This single change collapses the biggest bottleneck in the entire project.

**The deepest principle:** Autonomy is not a single dial. It is a map. The agent should be autonomous in *execution* and constrained in *direction*. You set the direction at checkpoints; the agent executes between them. The checkpoints are fast (15-30 min each) but mandatory. The execution loops can run for hours without you. The ratio of your cognitive time to agent compute time should be small -- but your cognitive time should be spent on the things only you can do: deciding what questions matter, interpreting what results mean, and sitting with confusion until understanding arrives.
