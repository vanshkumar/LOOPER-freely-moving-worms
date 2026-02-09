# Suggestions for the December 2025 Aliveline Project

Prepared after reading: the original LOOPER-specific CLAUDE_SUGGESTIONS.md and CODEX_SUGGESTIONS.md, all 31 daily log entries plus Reflections/Results/Email to Brennan, all LOOPER repo documentation (README, EXPERIMENTS, METHODS, RESULTS, RUNBOOK, LEARNINGS), 21 Codex CLI session logs, 35 ChatGPT conversation threads (~3,900 turns), the AsymmetricDiffusionMapping and category-learning_mpfc repos, the qsimeon_data documentation, and the project_chats_subset conversation index.

The original CLAUDE_SUGGESTIONS.md (in the LOOPER directory) offered strong process and science advice. This document builds on it with the additional context from the full multi-tool, multi-repo record of the project. Where the original is sufficient, I point back to it rather than repeating.

---

## Part 1: How you actually used LLMs -- a multi-tool picture

The original CLAUDE_SUGGESTIONS.md diagnosed a single pattern: "LLM used for idea generation where it shouldn't have been." The full record reveals something more nuanced. You used **three distinct LLM modalities** -- ChatGPT (web), Codex (CLI agent), and Claude (CLI agent) -- and the failure modes were different for each.

### 1.1 ChatGPT: thinking partner, paper companion, social coach

The 35 ChatGPT conversations (~3,900 turns, 11:1 assistant-to-user turn ratio) show ChatGPT used in several modes:

**Where it excelled:**
- **Paper interrogation.** Uploading PDFs and asking specific methodology questions ("did this paper collect their own data?", "explain their asymmetric diffusion map method"). The conversations on the Brennan-Proekt, Reinert, Kato, and DD-DC papers show genuine understanding being built through dialogue.
- **Critical feedback on your conclusions.** When you uploaded RESULTS.md, ChatGPT pushed back: "If your stationarity test fails on the positive control (Kato), then it's not a valid 'positive control stationarity test.' It becomes a *stress test*." That reframing was correct and important.
- **Social/professional coaching.** The email to Connor Brennan was drafted and refined with ChatGPT's help. ChatGPT's assessment -- "It's fine. Two things make it not read like 'please supervise me'" -- was calibrated and helpful.
- **Long, deep working sessions.** Your most productive ChatGPT conversations had 200+ turns ("Locomotion manifold comparison," "Dimensionality Reduction Techniques," "Worms Dataset Methodology Summary"). These were genuine extended thinking sessions where you built understanding incrementally.

**Where it failed or was misused:**
- **Idea generation from scratch.** This is the same finding as the original document, but the chat logs make it more concrete. Early conversations ("Neural trajectory curvature analysis") show ChatGPT producing plausible-sounding frameworks that you never returned to. The ideas were well-formed but not *yours*, so you didn't have the intuition to evaluate them.
- **Web search chains.** ChatGPT frequently performed extensive browsing to find paper details that could have been answered faster by reading the paper directly. The model's browsing was thorough but slow and sometimes circular.
- **Premature experiment design.** Day 24's rejection -- "I don't think the original ChatGPT experiment suggestion makes sense. It's overindexing on the speculation in the Moore paper" -- captures a recurring pattern: ChatGPT generates detailed experiment plans that sound rigorous but are built on assumptions you haven't verified.

**The most important pattern:** Your custom instructions to ChatGPT reveal your intent: *"Always pursue truth. No confirmation bias towards whatever I'm saying or suggesting, minimal editorializing unless asked, zero patronizing. Speak as simply as possible."* This is excellent. But even with these instructions, ChatGPT's confirmation bias operates at a subtler level: it *engages enthusiastically with whatever framing you bring*, which feels like validation even when it's just conversational responsiveness. The 20 minutes of doing nothing that cracked your understanding worked precisely because it offered zero social reinforcement.

### 1.2 Codex: code generation, repo comprehension, documentation

The 21 Codex sessions show a different tool used for different purposes:

**Where it excelled:**
- **Repository onboarding.** Codex consistently produced excellent structured summaries of repos it had never seen, identifying key files, data structures, and analysis pipelines. This saved hours of manual reading.
- **MATLAB pipeline construction.** The session auto-detection code, helper functions, eigenvalue analysis pipeline, and data loading code were all generated competently.
- **Documentation-as-specification.** Your pattern of having Codex create AGENTS.md, METHODS.md, GOALS.md, and DATASET.md files served a dual purpose: they guided Codex *and* documented the project for you. This is a strong workflow pattern worth keeping.
- **The LEARNINGS.md institutional memory system.** This carried context across stateless sessions. When Codex opened a new session on Jan 9, it found and read the LEARNINGS.md from Jan 8, picking up where it left off. The Codex LEARNINGS.md accumulated genuinely useful entries (pcaBasis caching, Rcorr vs Rsquared, saveData field reattachment).

**Where it failed:**
- **Repeated PDF extraction failures.** Codex tried PyPDF2/pypdf across multiple sessions (Jan 20, 22, 26, 27) and failed every time because the packages weren't installed. Despite the LEARNINGS.md system, it didn't learn to stop trying this approach. This suggests the LEARNINGS.md system works for *code-level* gotchas but not for *environment-level* constraints unless they're phrased as explicit prohibitions.
- **Context loss between sessions.** Each new Codex session re-scanned the repo from scratch. You had to tell Codex to "read the AGENTS.md" even after it had been created in the same long session. The LEARNINGS.md only partially addresses this.
- **Scientific reasoning.** Codex identified important caveats (boundary stimulus handling, Go/NoGo confounds) but needed to be redirected toward the actual scientific question. It is a builder, not an architect.

**Key difference from ChatGPT:** You used ChatGPT for *thinking* and Codex for *doing*. This is roughly the right division. The problem is that Codex sometimes tried to *think* (choosing what to build next) and ChatGPT sometimes tried to *do* (generating experiment code inline), and both were worse at the other's job.

### 1.3 What the multi-tool picture adds to the original recommendation

The original CLAUDE_SUGGESTIONS.md recommended: "You generate hypotheses, LLM stress-tests them." The chat logs refine this into a more specific protocol:

| Phase | Best tool | Your role | LLM role |
|-------|-----------|-----------|----------|
| Formulating the question | Blank paper + you | Think alone for 20-30 min | None yet |
| Stress-testing the question | ChatGPT (extended session) | Bring a written hypothesis | Attack it, find counterarguments, suggest papers |
| Designing the experiment | ChatGPT + you | Write the experiment contract first | Critique the contract, flag assumption gaps |
| Understanding the data | Codex / Claude (repo agent) | Point it at the right files | Read, summarize, flag gotchas |
| Building the pipeline | Codex / Claude (code agent) | Specify the interface and expected behavior | Write the code, run it, debug it |
| Interpreting results | You first, then ChatGPT | Write what you think the results mean | Push back, suggest alternative explanations |
| Writing / communicating | You first, then LLM for editing | Write the first draft (email, paper, notes) | Edit for clarity, check for overclaiming |

The key addition is that **the "think alone" phase is mandatory and comes before any LLM interaction**. The original document recommended this; the chat logs confirm it's the single highest-leverage change.

---

## Part 2: Process -- what the full record reveals beyond the original

### 2.1 Execution bandwidth > epistemic bandwidth (from CODEX_SUGGESTIONS)

The CODEX_SUGGESTIONS.md nailed the core diagnosis: your code advanced faster than your understanding. The chat logs provide the mechanism: **LLMs make execution feel frictionless, which removes the natural feedback signal that forces you to slow down and understand**.

When you write code by hand, getting stuck on a bug forces you to understand what the code is doing. When Codex writes the code, the bug still exists but you don't encounter it until the results are wrong -- and by then the distance between your mental model and the actual computation is large.

**Concrete manifestation from the logs:**
- The split-half stationarity test was designed, implemented, and run before you understood that it was a stricter criterion than the LOOPER paper's own validation. This was discovered only when *both datasets failed*, forcing a rethink.
- The heat-pulse experiment was designed before carefully reading the Atanas paper's methods section to verify that the heat pulse was a locomotion perturbation (it wasn't -- it was nociception).
- The QSimeon dataset was assumed to be directly comparable to the Kato paper's data before discovering the WT_Stim vs. NoStim mismatch.

All three could have been caught by the "experiment contract" recommended in the original document. The CODEX_SUGGESTIONS.md adds the "dataset suitability gate" and "metric validation ladder" which are complementary to the experiment contract. Together, they form a three-gate system:

1. **Question gate** (is the question well-formed and falsifiable?)
2. **Dataset gate** (can this dataset answer this question?)
3. **Metric gate** (does this metric measure what I think it measures?)

All three gates must pass before writing code. This adds maybe 1-2 hours per experiment but prevents the multi-day "the question didn't match the data" pivots that dominated the second half of the project.

### 2.2 The three-LLM-role separation

The CODEX_SUGGESTIONS.md proposed separating Architect / Builder / Auditor roles. The chat logs show this is approximately what you already did -- ChatGPT was the architect, Codex was the builder -- but you didn't enforce the boundaries:

- ChatGPT sometimes acted as builder (generating inline code for experiments).
- Codex sometimes acted as architect (choosing what to build next).
- Neither consistently acted as auditor.

**Concrete recommendation:** Make the auditor role explicit. After any major result, have a *separate* LLM session (or a separate ChatGPT conversation) whose only job is to attack the result. Prompt: "Here are my results and my interpretation. Your job is to find every way this could be wrong, misleading, or uninterpretable. Do not suggest what to do next -- only identify problems." This separation prevents the natural conversational tendency to be helpful (suggesting fixes) from short-circuiting the critical evaluation.

### 2.3 Context management for large projects

The project spanned multiple repos, multiple LLM tools, and multiple months. Context fragmentation was a real problem:

- ChatGPT conversations became unwieldy at 200+ turns (high branch depths, possible circular patterns).
- Codex sessions were stateless and required re-bootstrapping.
- Daily logs captured your thinking but were not systematically consulted during coding sessions.
- LEARNINGS.md captured technical gotchas but not scientific reasoning.

**Recommendations for next time:**

1. **Single source of truth for the scientific state.** Maintain a `STATE.md` (or similar) at the top level that is updated weekly with: current hypothesis, current evidence for/against, current blockers, next experiment. This is what you consult before starting any LLM session, and what you update after any major result. It is shorter and more actionable than daily logs.

2. **Conversation hygiene.** For ChatGPT, start a new conversation every ~50 turns or when the topic shifts significantly. Long conversations accumulate stale context that biases the model's responses. The 200+ turn conversations were productive but likely became less efficient over time.

3. **For very large chat corpora, use chunked map-reduce summarization with parallel workers, then synthesize from aggregated summaries** (per your CLAUDE.md update). This is the right approach and generalizes: any time you're asking an LLM to process more material than fits in a single context, chunk first.

4. **Bridge documents between tools.** When switching from ChatGPT (thinking) to Codex (building), write a brief handoff document: "Here is the question, the experiment contract, and the expected behavior. Build this." This prevents the builder from re-deriving (or worse, silently modifying) the architect's decisions.

### 2.4 The expert contact bottleneck (reinforced)

The original document recommended emailing paper authors in week 1. The full record makes this even more urgent:

- The email to Connor Brennan was the sharpest intellectual output of the entire project.
- It came on *Day 31* -- the last day.
- The follow-up question ("is distinguishing intrinsic attractors from closed-loop feedback a useful question?") could have redirected the entire second month if asked in week 2.
- ChatGPT's advice on the email tone was helpful but is no substitute for the actual expert response.

The activation energy for cold-emailing is real. A concrete trick: **write the email as if you're writing it, but don't send it yet. Use it as a thinking tool. Then send it.** Your daily logs show that the act of writing the email *was itself* a clarifying exercise: "an extremely useful exercise in thinking through wtf I'm actually doing."

### 2.5 The holiday break and momentum loss

The 2.5-week break (Dec 21 - Jan 9) caused significant momentum loss and is the single biggest structural problem in the project timeline. Your own diagnosis is correct: "splitting this over 2 months played a part" in losing focus.

The chat logs show the pre-break period (Days 1-17) was dominated by exploration and reading, while the post-break period (Days 18-31) was dominated by execution and coding. This is a reasonable arc, but the break meant that the transition from exploration to execution required *re-bootstrapping* your understanding, which cost about a week.

**For next time:** If a break is unavoidable, write a "re-entry document" before the break: what you currently believe, what you were about to do, what the key open questions are, and what files/papers you should read first upon return. This is different from daily logs (which are chronological) -- it's a snapshot of your mental state optimized for future-you.

### 2.6 Context switching (reinforced with evidence)

The daily logs mention context switching as a problem. The codex chat timestamps make it concrete: sessions were often short (30-60 minutes) with gaps of days between them. The ChatGPT conversations show a similar pattern -- many conversations were used in bursts rather than sustained sessions.

The highest-quality thinking in the logs correlates with the longest uninterrupted blocks. Day 21 (the detailed experiment taxonomy with falsification criteria) and the 20-minute "doing nothing" insight are both products of sustained focus.

**Concrete recommendation:** For the next project, designate 2-3 days per week as "deep work" days with no administrative tasks, no email, and no context-switching between projects. Use the remaining days for admin, correspondence, and exploratory reading. This is a well-known productivity pattern (Cal Newport's "Deep Work") but the logs provide direct evidence that it would have helped here.

---

## Part 3: The intellectual trajectory -- what the full record reveals

### 3.1 The actual arc of ideas

The daily logs show the *what* of your intellectual journey. The chat logs show the *how* -- the specific conversations, paper readings, and reasoning steps that drove each transition. Reconstructed from both sources:

**Phase 1: Casting the net (Days 1-6)**
- Started with Kolmogorov complexity as a lens on brain function
- Explored Tracy-Widom phase transitions, neural population geometry, the brain as a dynamical system
- Key ChatGPT conversations: "Neural trajectory curvature analysis," "Neural population transitions"
- Generated many interesting but unfocused ideas

**Phase 2: The mouse PFC detour (Days 7-12)**
- Settled on the Reinert mouse mPFC dataset as the most tractable test
- Built an eigenvalue analysis pipeline (via Codex) to look for BBP-style phase transitions during rule learning
- Discovered the confound problem is enormous ("it seems way harder in 'higher' animals like mice")
- Key realization: C. elegans is a better model system for these questions
- Key Codex sessions: Dec 12-13 (building the eigenvalue pipeline)

**Phase 3: Deepening into C. elegans (Days 13-20)**
- Read the aversive learning paper (Liang et al.), the Kato paper, the Brennan-Proekt paper
- Developed the "macroscopic universality, microscopic degeneracy" framework
- Connected locomotion manifolds, proprioception, and context-gated learning
- Key ChatGPT conversations: "C elegans locomotion learning," "Worm memory and behavior"
- Holiday break (Dec 21 - Jan 9)

**Phase 4: The DD-DC breakthrough (Days 21-23)**
- Read "The neuron as a direct data-driven controller" (Moore et al.)
- Integrated it into a multi-level hypothesis: genes -> feedback control laws -> computational scaffolds -> behavior
- Produced the detailed experiment taxonomy (Day 21) with falsification criteria
- Key ChatGPT conversations: "Neuron as DD-DC Controller," "DD-DC as Bridge Language," "Distributed Feedback Control Attractors"

**Phase 5: Execution and negative results (Days 24-31)**
- Applied LOOPER to Kato (positive control) and Atanas (target) datasets
- Discovered: fidelity passes for both, but stationarity fails for both
- Put heat-pulse experiments on hold
- Wrote the email to Connor Brennan
- Key Codex sessions: Jan 22-28 (running LOOPER pipeline, debugging MATLAB issues)
- Key ChatGPT conversations: "Code Review and Metrics," "LOOPER loop recovery comparison"

### 3.2 What the chat logs add to the original science suggestions

The original CLAUDE_SUGGESTIONS.md proposed 7 priority experiments. The chat logs and other repos provide additional context that refines these priorities:

**The mouse mPFC work is not a detour -- it's a parallel test of the same framework.** The category-learning_mpfc repo implements the same conceptual toolkit (eigenvalue analysis, covariance structure, phase-transition-like phenomena) applied to a different organism and cognitive domain. If the LOOPER/manifold framework works in C. elegans, testing the same spectral-emergence hypothesis in mouse mPFC during rule learning would strengthen the case for universality. This is a potential second paper or a section of a larger paper.

**The ADMM repo is a necessary comparison.** The AsymmetricDiffusionMapping repo provides the original Brennan-Proekt method as a benchmark. Running ADMM and LOOPER on the same QSimeon data and comparing their scaffold structures is a high-value experiment that the original CLAUDE_SUGGESTIONS doesn't mention. If both methods recover the same scaffold topology, the finding is robust to methodological choices. If they diverge, understanding why tells you something about what each method is actually measuring.

**The ChatGPT conversation "Kato stimulus and recovery" contains a specific, actionable insight:** The Kato paper's stimulus is O2 chemosensory stimulation, and the paper explicitly reports that the manifold topology doesn't change -- only occupancy and phase do. This means a "recovery" experiment on Kato WT_Stim should look for de-entrainment and phase distribution relaxation, not manifold disappearance/reappearance. This is a more precise experimental design than the original's Priority 6 (perturbation recovery).

**Revised priority ordering** (incorporating all sources):

| Priority | Action | Why it moved |
|----------|--------|-------------|
| 1 | Null model controls (time-shuffled) | Unchanged -- mandatory |
| 2 | ADMM vs. LOOPER comparison on same data | New -- validates methodology before interpreting results |
| 3 | Kato stationarity diagnosis (k-fold, windowed) | Unchanged |
| 4 | Phase velocity vs. behavior correlation | Unchanged |
| 5 | Behavior-conditioned LOOPER (using Atanas behavioral covariates) | Unchanged, but now explicitly use velocity, angular_velocity, forwardness fields |
| 6 | Cross-individual scaffold comparison | Unchanged |
| 7 | Kato WT_Stim de-entrainment analysis | Refined from original's "perturbation recovery" -- use O2 transitions, not heat pulse |
| 8 | Mouse mPFC eigenvalue emergence | New -- parallel test of spectral-emergence framework in different species |

### 3.3 The bigger question and where it stands

Your underlying question evolved across the project:

- Start: "What is the minimal description of brain computation?" (Kolmogorov complexity framing)
- Middle: "How does DNA bridge micro (neurons) to macro (behavior)?" (renormalization group framing)
- End: "Genes encode behavior as top-down constraints via modifying neurons' nonlinear feedback control algorithms" (DD-DC framing)

This is a genuine intellectual trajectory. The DD-DC framing is the most specific and testable version. It predicts:

1. Feedback control gains/parameters should be heritable (genetically specified).
2. Different individuals should share the same control law but with different internal coordinates (consistent with LOOPER's "conserved scaffold, variable coordinates" finding).
3. Perturbation recovery should follow control-theoretic dynamics (exponential return with characteristic time set by feedback gain), not passive attractor dynamics.
4. Eliminating sensory feedback (immobilization) should change the scaffold structure if feedback is essential, or leave it intact if the scaffold is largely intrinsic.

Predictions 2 and 4 are testable with existing data. Prediction 3 requires perturbation data. Prediction 1 requires genetic manipulation data (available in some C. elegans datasets but not in Atanas or Kato).

---

## Part 4: Process recommendations -- the consolidated list

These integrate the original CLAUDE_SUGGESTIONS, the CODEX_SUGGESTIONS, and the new insights from the full record.

### For the next research project

| # | Recommendation | Source |
|---|---------------|--------|
| 1 | **Think before prompting.** 20 min blank page before any LLM brainstorm. | Original, confirmed by chat logs |
| 2 | **Write experiment contracts.** Question + data + method + expected results + kill condition. | Original + CODEX_SUGGESTIONS |
| 3 | **Three gates before code.** Question gate, dataset gate, metric gate. | CODEX_SUGGESTIONS, new synthesis |
| 4 | **Email authors in week 1.** Write the email as a thinking tool, then send it. | Original, reinforced |
| 5 | **Protect deep work blocks.** 2-3 days/week with no admin or context-switching. | Original, reinforced by timestamp analysis |
| 6 | **Separate LLM roles.** Architect (ChatGPT), Builder (Codex/Claude), Auditor (separate session). | CODEX_SUGGESTIONS, refined |
| 7 | **Maintain a STATE.md.** Weekly snapshot: hypothesis, evidence, blockers, next experiment. | New, from context fragmentation analysis |
| 8 | **Write a re-entry document before any break.** Snapshot of mental state for future-you. | New, from holiday break analysis |
| 9 | **New ChatGPT conversation every ~50 turns.** Long conversations accumulate stale context. | New, from chat analysis |
| 10 | **Bridge documents between tools.** Handoff from thinking (ChatGPT) to building (Codex). | New, from multi-tool analysis |
| 11 | **Null models are mandatory.** Always know what chance looks like. | Original |
| 12 | **Study fundamentals between projects.** Strogatz (nonlinear dynamics), control theory. | Original |
| 13 | **LEARNINGS.md for environment constraints, not just code gotchas.** Include "do not try X in this sandbox." | New, from Codex PDF extraction failures |
| 14 | **For large corpora, use chunked map-reduce with parallel workers.** Don't attempt monolithic reads. | Your CLAUDE.md update |
| 15 | **Track claims, not narratives.** Use a CLAIMS_LEDGER with status, evidence, and known confounds. | CODEX_SUGGESTIONS |

### For LLM-assisted science specifically

The single most important insight, stated more precisely than the original:

**The danger of LLMs in research is not that they give wrong answers. It is that they make the wrong *questions* feel productive.** The feeling of progress from a detailed LLM-generated experiment plan is almost indistinguishable from the feeling of progress from actually understanding the problem. The only reliable way to tell the difference is to spend time thinking without any LLM at all, and to notice whether your understanding survives the silence.

Your Reflections entry captures this perfectly: "wow just sitting doing nothing for 20m and i thought through what the brennan paper's loops thing was actually doing, how they measured its quality, if my split half / heat pulse experiments actually make sense (they don't)."

That 20 minutes was the highest-ROI activity of the entire project. It should be the default starting point, not an exception.

---

## Part 5: What you should be proud of

The original document emphasized what went wrong. The full record also reveals what went right in ways that deserve acknowledgment:

1. **You entered a new field and developed a real hypothesis in 2 months.** The DD-DC + LOOPER + manifold synthesis is not trivial. Many people with formal training take longer to formulate a question this specific.

2. **You ran real experiments and got real results.** The LOOPER reproduction on Kato, the negative result on Atanas, the stationarity analysis -- these are genuine scientific outputs, not just reading and thinking.

3. **You built a documented, reproducible codebase.** The LOOPER-freely-moving-worms repo has EXPERIMENTS.md, METHODS.md, RESULTS.md, RUNBOOK.md, and LEARNINGS.md. This is better documentation than most academic code.

4. **You wrote a genuinely good email.** The Connor Brennan email compresses months of thinking into two precise questions. This is what scientific communication looks like.

5. **You maintained intellectual honesty throughout.** The daily logs show no instance of overclaiming or ignoring negative results. The Reflections document is unflinching. The negative stationarity result is presented as what it is -- a finding that the test is stricter than expected -- not spun into a false positive.

6. **You iterated on your own process.** The LEARNINGS.md system, the experiment contract idea, the "think before prompting" realization -- these are meta-improvements that compound across future projects.

7. **You used multiple LLM tools in complementary ways.** ChatGPT for thinking, Codex for building, custom GPTs for domain-specific conversations. This multi-tool approach is sophisticated and will become more valuable as the tools improve.

8. **You documented everything.** 31 daily logs, complete chat exports, Codex session recordings. This level of documentation makes retrospective analysis (like this document) possible and is itself a contribution to understanding how LLM-assisted research works.

The "ultimate result" is not a single email. It is: a refined hypothesis, a negative result that is itself informative, a documented codebase, a process for LLM-assisted research that you can iterate on, and a network connection with a domain expert. That is a solid foundation for whatever comes next.
