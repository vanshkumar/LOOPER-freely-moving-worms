# Suggestions for LOOPER-freely-moving-worms

Prepared after reading the full repository (code, docs, results, LEARNINGS.md), the LOOPER paper + supplement, and all 34 entries of the December 2025 Aliveline daily logs documenting the intellectual journey from start to finish.

---

## Part 1: Improving the process of co-working with LLMs to do science

The first version of this document offered generic process advice (pre-register decisions, add negative controls, version-pin code). That advice still stands, but after reading the Aliveline logs I think the important observations are at a different level. The logs reveal a specific person's actual struggle with a specific problem: **how do you use LLMs productively when entering a new field, under time pressure, with no mentor, trying to generate an original insight?** That is a much harder and more interesting question than "how do you organize your MATLAB pipeline."

### 1.1 What actually went well

**The question-first / research-first approach to learning a new field.** The Aliveline challenge forced you to learn computational neuroscience by trying to answer a question, not by working through a textbook. This is underrated. The sequence — Kolmogorov complexity as a lens on brain function -> C. elegans as model organism -> neural population geometry -> locomotion manifolds -> LOOPER -> DD-DC feedback control — is a genuine intellectual trajectory driven by curiosity and successive refinement. Many people with formal training never develop this kind of taste-driven exploration because the curriculum imposes the path for them.

**The evolving hypothesis.** Your thinking progressed from a vague "brain complexity" framing to the quite specific: "Genes encode behavior as top-down constraints on neuronal networks by modifying neurons' nonlinear feedback control algorithms, resulting in conserved computational scaffolds (loops) despite microscopic degeneracy across individual brains." That hypothesis connects LOOPER, DD-DC, the Atanas dataset, and the broader question of how a 302-neuron nervous system is specified by a genome. It is a real scientific hypothesis — testable, non-obvious, and grounded in existing work. Getting there in ~2 months as a newcomer to the field is not trivial.

**The email to Connor Brennan.** Your Reflections note that this was "extremely useful" and "the ultimate result." This is exactly right. The follow-up email — "is distinguishing intrinsic attractors from closed-loop feedback a useful question?" — is a precisely formulated question that compresses the entire project into something a domain expert can engage with. The fact that it took 2 months to get to that question is not a failure; it is the normal timescale of question formulation.

**LEARNINGS.md as institutional memory.** This pattern worked well for the mechanical parts of the project (MATLAB gotchas, caching bugs, struct reattachment). It is genuinely good practice.

**Incremental pipeline validation.** Single worm -> all worms -> shared neurons -> stationarity. Sound progression that caught bugs early.

### 1.2 What the logs reveal about where the process broke down

#### A. The LLM was used for idea generation where it shouldn't have been

This is the single most important observation from the logs.

You noticed this yourself in Reflections.md: "found myself relying on LLMs much more than is probably ideal for the idea generation parts of this." The logs show a pattern where ChatGPT was asked to generate hypotheses, suggest experiments, and evaluate ideas — and the results were often generic, surface-level, or actively misleading (e.g., the early suggestion to look at Kolmogorov complexity of neural recordings, which is not a productive direction; or the oscillation between mouse PFC data and C. elegans before settling on a question).

The problem is not that LLMs are bad at idea generation. It is that **LLM idea generation feels productive while often being a substitute for the hard, slow, uncomfortable work of sitting with confusion.** Your own log entry is the best evidence: "wow just sitting doing nothing for 20m and I thought through what the Brennan paper's loops thing was actually doing, how they measured its quality, if my split half / heat pulse experiments actually make sense (they don't)." That 20 minutes of staring at the wall produced more insight than many hours of LLM conversation.

**Concrete recommendation:** Institute a rule for yourself: before asking an LLM to brainstorm, spend 20 minutes with a blank page (paper, not screen) and write down what you actually think. Not what sounds smart, not what you think the LLM would say — what you actually believe and what confuses you. Then use the LLM to stress-test your thinking, not to generate it. The LLM is excellent at: "here is my hypothesis, what are the three strongest objections?" It is mediocre at: "what should I think about?"

A practical version of this for LLM-assisted science:

1. **Think first** (20-30 min, no screen): What do I actually believe? What confuses me? What would change my mind?
2. **Write a one-paragraph hypothesis** in your own words.
3. **Ask the LLM to attack it**: "Here is my hypothesis. What are the three strongest reasons it might be wrong? What experiment would most efficiently distinguish it from the two nearest alternative hypotheses?"
4. **Then** ask for literature, methods, code.

This ordering matters because it preserves your intellectual ownership of the question. The alternative — asking the LLM to generate the question — tends to produce questions that are well-formed but not deeply yours, which means you won't have the intuition to know when the results are surprising.

#### B. Insufficient time on question formulation, too fast to experiments

The logs show the following pattern repeating:

1. Interesting idea surfaces (e.g., "do locomotion loops survive in freely moving worms?")
2. Within hours/days, LLM is asked to write experiment code
3. Experiment runs, produces ambiguous results
4. Realization that the question wasn't well-formed (e.g., the heat pulse recovery idea, the half-split stationarity test)
5. Return to question formulation

You identified this yourself: "One major takeaway is to spend way more time formulating questions precisely. Solving problems is 90% about understanding them."

The specific instance that cost the most time was the heat pulse experiment design. The original plan was to test whether locomotion loops recover after perturbation (heat pulse), which would distinguish intrinsic attractors from feedback-stabilized dynamics. But the heat pulse data in the Atanas dataset turned out to be a thermal nociception response, not a perturbation to an ongoing locomotion loop. This mismatch between the data and the question could have been caught earlier by carefully reading the Atanas paper's methods section before writing any code.

**Concrete recommendation:** For each new experiment/analysis, write a one-page "experiment contract" BEFORE any code:

```
QUESTION: [one sentence, falsifiable]
DATA: [which specific data files, which fields, which time windows]
METHOD: [which analysis, which parameters, which comparisons]
EXPECTED RESULTS IF TRUE: [specific numbers/patterns]
EXPECTED RESULTS IF FALSE: [specific numbers/patterns]
ASSUMPTIONS: [list everything the question assumes about the data]
```

Then read the data paper's methods section to verify every assumption. Then have the LLM critique the contract. Only then write code. This adds maybe an hour of upfront work but saves days of "the question didn't match the data" pivots.

#### C. The expert-contact bottleneck

Your logs repeatedly reference wanting a mentor or advisor. The email to Connor Brennan happened in late January — roughly 6 weeks into the project. You note: "Step 1 after reading any paper should be emailing the author with (good) questions." This is exactly right.

LLMs are not a substitute for talking to the person who built the method. Connor Brennan can tell you in one email whether distinguishing intrinsic attractors from feedback control is a meaningful question, whether LOOPER is the right tool for it, and what the key pitfalls are. An LLM can generate plausible-sounding answers to these questions that may be wrong in subtle ways.

**Concrete recommendation:** Make expert outreach a first-week activity, not a last-week activity. The template:

1. Read the paper carefully (the first time, just read — don't code).
2. Write down your 2-3 best questions. Use the LLM to sharpen them if needed.
3. Email the corresponding author within the first week. Short email, specific questions, link to your preliminary thinking.
4. Continue working while waiting for a response, but hold major design decisions loosely.

The activation energy for this is high ("I don't know enough to ask good questions yet"), but the bar for "good question" from a genuinely curious newcomer is much lower than you think. Researchers are generally happy to hear from someone who read their paper carefully.

#### D. Context switching killed momentum

The logs show significant context switching — coworking admin, neuronal dynamics group logistics, email, and other projects. You note this yourself: "Frankly I did way too much context switching, especially in the back half of this work."

This matters more for question formulation than for code execution. Writing code with an LLM is interruptible; you can pick up where you left off because the state is in the codebase. Developing an understanding of a new field is not interruptible the same way. The 20-minute thinking session that cracked your understanding happened precisely because you had unbroken time.

**Concrete recommendation:** If doing another project like this, protect 2-3 hour blocks specifically for "thinking and reading, no coding." These blocks are the bottleneck — the LLM can write code 24/7, but you are the only one who can sit with the papers and develop intuition.

#### E. Where LLMs actually add the most value (and it's not idea generation)

Based on the logs, the highest-value LLM uses in this project were:

1. **Translating between your hypothesis and code** — once you knew what question to ask, the LLM could build the pipeline quickly. The LOOPER.m wrapper, the evaluation framework, the diagnostic plots — these would have taken weeks by hand and took days with an LLM.

2. **Literature search and paper interrogation** — asking "what does this equation mean" or "how does this method compare to X" while reading papers. The LLM is a tireless reading companion.

3. **Debugging MATLAB** — the LEARNINGS.md entries are evidence that the LLM caught and fixed non-trivial bugs efficiently.

4. **Stress-testing ideas** — when you brought a specific hypothesis and asked the LLM to critique it, the results were much better than when you asked it to generate ideas from scratch.

The lowest-value uses were: generating research directions from scratch, brainstorming without a specific hypothesis, and being asked "what should I do next" without constraints.

**Concrete recommendation for the next project:**

| Task | LLM role | Your role |
|------|----------|-----------|
| Idea generation | Stress-test your ideas | Generate the ideas yourself |
| Question formulation | Critique your questions | Write the first draft |
| Literature search | Find and summarize papers | Decide which are relevant |
| Experiment design | Write the code, critique the plan | Design the experiment |
| Debugging | Fix the bugs | Understand why they happened |
| Interpretation | Generate alternative explanations | Decide what you believe |
| Communication | Edit your writing | Write the first draft |

The common thread: **you provide direction and judgment, the LLM provides speed and breadth**. When the roles are reversed, the work feels productive but is often circular.

#### F. The explore-exploit tradeoff in entering a new field

The logs capture a real tension: "Viscerally felt the explore-exploit tradeoff / breadth vs depth thing when getting into a new field!" The path through Kolmogorov complexity -> population geometry -> manifold learning -> LOOPER -> DD-DC covered a lot of ground but left each topic shallow.

You suggest "alternating stupid objectives like this with free flowing exploration blocks would be very effective." I think this is right but incomplete. The missing piece is **deep study of fundamentals** — you note this too: "a couple months of deep study in a few relevant fields (linear algebra, nonlinear dynamics, neuronal dynamics, random matrix theory) would have accelerated this significantly."

The insight is that depth in fundamentals is what makes exploration efficient. Without a solid grasp of dynamical systems, it is hard to evaluate whether LOOPER's "loops" are the same kind of object as a limit cycle in a nonlinear ODE, or whether the DD-DC feedback control framing is mathematically compatible with the manifold learning framing. These connections are where original insights live, and they require depth.

**Concrete recommendation:** For the next phase, consider spending 4-6 weeks on Strogatz's "Nonlinear Dynamics and Chaos" (the textbook, not a summary) with an LLM as tutor. Not to read it cover-to-cover, but to work through the chapters on limit cycles, bifurcations, and stability — the math that LOOPER implicitly assumes. This would let you evaluate the feedback control hypothesis with real mathematical precision rather than analogy.

### 1.3 What the original version got right (retained)

The technical process recommendations from the first version remain valid:

- **Pre-register analysis decisions** in a DECISIONS.md before coding.
- **Add negative controls** (shuffled data) at every stage.
- **Parameter sensitivity sweeps** rather than single-point runs.
- **Version-pin external code** (record LOOPER commit hashes).
- **Composite summary figures** for cross-session continuity.

These are good hygiene. But they are secondary to the higher-level issues above.

---

## Part 2: Suggestions on the science

### 2.1 Stepping back: what question are you actually trying to answer?

The repo frames the question as: "Are locomotion loops in neural state space intrinsic attractors or feedback-stabilized?" This is the question in your second email to Connor Brennan, and it is well-posed. But reading the full arc of the Aliveline logs reveals a deeper question underneath:

**How does a genome (which is fixed) specify neural dynamics (which must be flexible) in a way that produces reliable behavior (which is conserved across individuals)?**

This is the question that connects your interest in Kolmogorov complexity ("what is the minimal description of brain computation?"), the DD-DC paper ("neurons implement feedback control without explicit system identification"), and LOOPER ("different brains produce the same computational scaffolds despite different neural coordinates").

The intrinsic-vs-feedback question is one specific test of this bigger question. If loops are intrinsic attractors, the genome might specify them directly (e.g., by setting synaptic weights that create a limit cycle). If loops are feedback-stabilized, the genome specifies a control law, and the loop is an emergent property of the controller interacting with the body/environment. These have different implications for how the genome encodes behavior and how robust behavior is to perturbation.

I raise this because the bigger question should guide which experiments you prioritize. Some of the technical suggestions below are more relevant to it than others.

### 2.2 What the current results actually tell you

**Fidelity (in-sample):** Both Kato (immobilized) and Atanas (freely moving) show comparable in-sample reconstruction (recon_corr ~0.57-0.66). This means LOOPER can fit 1-D loop-like structure to both datasets.

**Stationarity (half-split):** Both datasets fail. The scaffold learned from the first half does not generalize to the second half.

The fact that **both** datasets fail stationarity is the most important result and is somewhat under-emphasized. It means one of:

1. Half-split stationarity is too strict a test for any calcium imaging data of this length (i.e., it's a methodological artifact).
2. Neural dynamics are genuinely non-stationary on 5-10 minute timescales, even in immobilized worms (which would be a real finding — the Kato paper assumes stationarity implicitly).
3. LOOPER's scaffold is sensitive to initial conditions / noise in a way that makes cross-temporal generalization fragile even when the underlying dynamics are stationary.

The fact that Kato also fails is actually more informative than Atanas failing, because Kato is the "easy" case. If Kato passes but Atanas fails, you learn about freely moving worms. If both fail, you learn about the method or about neural dynamics in general.

### 2.3 The feedback control hypothesis and how to test it

Your core hypothesis (from the DD-DC connection) is roughly: "LOOPER strands are the invariant sets created by distributed negative-feedback control. Branch/merge reflect policy switching. Feedback control is how DNA bridges micro (neuronal rules) and macro (conserved behavior)."

This is a compelling framing, but it makes specific predictions that your current pipeline does not test:

**Prediction 1: Perturbation recovery.** If loops are feedback-stabilized, perturbing the system should cause transient deviation followed by return to the scaffold. If loops are intrinsic attractors, perturbation recovery should also occur but with different dynamics (exponential return to limit cycle vs. active correction via sensory feedback). The Atanas heat-pulse data cannot test this cleanly because the heat pulse is a nociception stimulus that elicits a specific behavioral response, not a perturbation to ongoing locomotion.

To actually test this, you would need either:
- Optogenetic perturbation data (brief activation/silencing of specific neurons during ongoing locomotion)
- Natural perturbations (unexpected mechanical stimuli during forward crawling), which some C. elegans datasets include

The Atanas data may have natural perturbation events buried in the behavioral time series (e.g., brief stalls or direction changes). Identifying these and checking whether the LOOPER phase trajectory recovers afterward would be a low-cost first test.

**Prediction 2: Feedback requires sensory input; intrinsic attractors don't.** The Kato data is from immobilized worms — proprioceptive feedback is eliminated or greatly reduced. If loops survive in immobilized worms (they do, per the original Kato paper and your reproduction), this suggests that the locomotion loop has a strong intrinsic component, at least for the basic oscillation. But "intrinsic oscillator + feedback modulation" is not the same as "pure attractor" — the freely moving worm may use feedback to shape, stabilize, or gate the intrinsic oscillation.

This framing suggests a specific experiment: **compare the LOOPER scaffold structure (not just fidelity) between immobilized and freely moving worms.** If feedback modulates the loop, you might see the same basic topology but with different branching patterns, phase velocity profiles, or strand structures in freely moving worms. LOOPER's "merging" operation specifically captures where trajectories converge to a common scaffold — more merging in freely moving worms could indicate feedback correcting drift.

**Prediction 3: Conserved scaffolds despite different coordinates.** This is LOOPER's headline finding from the paper (Fig. 6 — same scaffold across species). Your project could test this within C. elegans: do different individual worms in the Atanas dataset produce the same scaffold topology? The shared-neuron analysis already goes in this direction, but a more direct test would be to compare the scaffold structure (number of strands, branching pattern, loop topology) across individual worms.

### 2.4 Highest-priority experiments (revised)

The priority ordering changes substantially when informed by the bigger question (how does the genome specify dynamics?):

#### Priority 1: Null model controls

*(Unchanged from original version — this is mandatory before anything else.)*

Run LOOPER on time-shuffled data (each neuron shuffled independently). This tells you whether recon_corr ~0.6 is meaningful or whether LOOPER always finds something loop-like due to smoothing and delay embedding. Without this, no result is interpretable. This is a one-session task.

#### Priority 2: Understand why Kato fails stationarity

This is more important than behavior-conditioned LOOPER on Atanas, because if the positive control fails, the framework itself is in question.

- Run k-fold (5-fold) temporal cross-validation on Kato instead of half-split. If 5-fold passes, the half-split is too noisy (explanation 3 above). If 5-fold also fails, neural dynamics are genuinely non-stationary or LOOPER's scaffold is fragile.
- Compute reconstruction correlation as a function of time-lag from training set. Does it degrade gradually (drift) or abruptly (mode switch)?
- Check whether the Kato paper's trial-based validation (which does work) succeeds because trials are short enough to be within the stationarity horizon.

#### Priority 3: Phase velocity vs. behavior (connects loops to locomotion)

Before conditioning LOOPER on behavioral state, first check whether the LOOPER phase variable tracks behavior at all:

- Compute LOOPER phase (theta) from the Atanas in-sample scaffold.
- Compute d(theta)/dt (phase velocity).
- Correlate with velocity, angular_velocity, forwardness from the behavioral covariates.
- Check whether phase reverses during reversal events.

If phase velocity does not correlate with any behavioral variable, the "loop" may not be a locomotion loop. This would fundamentally change the interpretation. If it does correlate, you have a concrete link between the abstract scaffold and behavior, which strengthens everything else.

This is important to do early because it determines whether the rest of the project is about locomotion (as assumed) or about some other neural dynamic.

#### Priority 4: Behavior-conditioned LOOPER

If phase velocity does correlate with behavior (Priority 3), then segment by behavioral state and run LOOPER within each mode:

- Forward crawling epochs (forwardness > threshold)
- Non-forward epochs (reversals, turns, pauses)
- Run fidelity + stationarity on each subset separately.

If within-mode scaffolds are more stationary, you have evidence that mode-switching explains the stationarity failure. If not, the non-stationarity is intrinsic to each behavioral mode.

#### Priority 5: Cross-individual scaffold comparison

Run LOOPER on each Atanas worm individually and compare the scaffold structure:

- Number of strands, branching pattern, loop count.
- Phase velocity distribution.
- Reconstruction quality distribution across neurons.

If different worms produce similar scaffolds (as LOOPER predicts), this is evidence for the "genome specifies the computational strategy" hypothesis. If scaffolds are highly variable, the conserved-scaffold story does not hold for freely moving worms.

#### Priority 6: Perturbation recovery from natural events

Identify moments in the Atanas behavioral time series where the worm experiences apparent perturbations (sharp velocity changes, brief stalls, collision-like events). Track the LOOPER phase trajectory through these events. Does the phase return to the pre-perturbation scaffold trajectory? How quickly? This is a low-cost test of the feedback control prediction.

#### Priority 7: Immobilized vs. freely moving scaffold topology

Compare the Kato and Atanas scaffold structures (not just summary statistics):

- Do freely moving worms have more branching / merging?
- Are there additional strands in freely moving worms that might correspond to feedback-modulated dynamics?
- Is the basic oscillatory structure (the main loop) conserved, with differences only in secondary structure?

This directly tests whether feedback adds structure to the scaffold.

### 2.5 What probably isn't worth doing right now

- **Parameter sensitivity sweep** — useful for robustness but won't advance the scientific question. Do it later for a paper, not now for understanding.
- **Kato reproduction gap analysis** (recon_corr vs. paper R^2) — the gap is likely parameter/neuron-set differences and won't change the story. Document it but don't spend time closing it.
- **SNR-matched Kato degradation** — interesting but second-order. If Priorities 2-3 give clear results, this becomes unnecessary.
- **Neuron dropout analysis** — cool but requires deeper neuroscience knowledge to interpret. Better done with a collaborator (e.g., Brennan or someone in the Atanas group).

### 2.6 Framing for a potential paper (revised)

The original framing was:

> "LOOPER recovers loop-like structure in freely moving C. elegans but the scaffold is non-stationary."

That is a negative result dressed up. The logs reveal a much more interesting story if the experiments work out:

> "Neural population dynamics in freely moving C. elegans contain locomotion-correlated loop structure that is individually recoverable in each animal. Unlike immobilized preparations, the freely moving scaffold is modulated by behavioral state: forward crawling, reversals, and turns produce distinct but related scaffolds that share a common oscillatory core. Cross-individual comparison reveals conserved scaffold topology despite variable neural coordinates, consistent with a model in which the genome specifies a control law rather than a fixed dynamical landscape."

This is stronger because it:
1. Connects to behavior via phase velocity (testable, not assumed).
2. Explains non-stationarity as behavioral mode-switching (mechanistic).
3. Connects to the genome/behavior question (your actual intellectual interest, and what makes this more than a methods paper).
4. Leverages NeuroPAL identity and behavioral covariates (Atanas dataset's unique strengths).

But this framing requires that Priorities 1-5 produce positive results. If phase velocity doesn't correlate with behavior, or if null models show comparable fidelity, the story changes fundamentally. Which is why the priorities are ordered as they are.

### 2.7 The bigger picture: where this could go

If the feedback-control framing holds up, the natural next question is: what is the control variable? In the DD-DC framework, a neuron controls a specific variable without explicit system identification. For the locomotion loop:

- Is the controlled variable body curvature? Phase of the motor sequence? Forward velocity?
- Can you identify the "error signal" (deviation of controlled variable from reference) in the neural dynamics?
- Does the LOOPER scaffold correspond to the invariant set of the closed-loop system?

These questions connect LOOPER (a data analysis tool) to control theory (a mathematical framework), which is exactly where your hypothesis lives. They also connect to the broader question of how genetics shapes neural computation — if the "control law" is what's genetically specified, the parameters of that law (gains, time constants, reference values) should be heritable and species-specific.

This is a multi-year research direction, not a one-month project. But having the direction is the point.

---

## Summary: priority-ordered action items

| Priority | Action | Why | Estimated effort |
|----------|--------|-----|-----------------|
| 1 | Null model controls (time-shuffled) | Without this, nothing is interpretable | 1 session |
| 2 | Kato stationarity diagnosis (k-fold, time-lag curve) | Positive control also fails; need to understand why | 1-2 sessions |
| 3 | Phase velocity vs. behavior correlation | Tests whether loops = locomotion | 1 session |
| 4 | Behavior-conditioned LOOPER | Tests mode-switching explanation | 2-3 sessions |
| 5 | Cross-individual scaffold comparison | Tests conserved-scaffold hypothesis | 1-2 sessions |
| 6 | Natural perturbation recovery | Low-cost test of feedback control | 1 session |
| 7 | Immobilized vs. freely moving topology | Tests feedback-adds-structure prediction | 1-2 sessions |

---

## Process recommendations summary

| Recommendation | One-liner |
|----------------|-----------|
| Think before prompting | 20 min blank page before any LLM brainstorm |
| Write experiment contracts | Question + data + method + expected results before code |
| Email authors in week 1 | Expert feedback >> LLM speculation |
| Protect thinking blocks | 2-3 hour no-code blocks for reading and pondering |
| LLM for critique, not generation | You generate hypotheses, LLM stress-tests them |
| Study fundamentals between projects | Strogatz's nonlinear dynamics would pay dividends |
| Null models are mandatory | Always know what chance looks like |
