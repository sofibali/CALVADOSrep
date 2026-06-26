# What does this all mean? — TICA, MSMs, and the tests we ran

*A guide for someone comfortable with PCA who hasn't built a Markov State Model (MSM).*

This explains the machinery behind the PARP14 macrodomain clustering: **TICA**, **Markov
State Models**, and the three diagnostics we keep quoting — the **implied-timescale (ITS)
test**, the **Chapman–Kolmogorov (CK) test**, and **VAMP-2 cross-validation** — plus
**PCCA+**, which turns thousands of microstates into the handful of "representative
states" we report. Every claim is sourced; full citations are at the end.

---

## 0. The one-paragraph map

You ran a simulation and want to summarise it as **a few representative structures and
how the molecule moves between them**. PCA would give you the directions of largest
*variance*. But the biggest-variance motion is not necessarily the biggest-*timescale*
motion — a floppy loop can wag with huge amplitude and zero kinetic relevance.
**TICA** fixes this by finding the slowest-*decorrelating* directions instead. We then
chop the TICA space into many small **microstates**, count transitions between them at a
chosen time interval (the **lag time** τ) to build a **Markov State Model** — a transition
matrix whose **eigenvalues are relaxation timescales** and whose **eigenvectors are the
slow processes**. The ITS and CK tests check whether that matrix is actually Markovian
(memoryless); the VAMP-2 score tells us which input features capture the slow dynamics
best; and PCCA+ lumps the microstates into a few metastable **macrostates** — the states
we hand to BindCraft. The whole pipeline rests on a **variational principle** (Noé &
Nüske 2013) that makes "which model is better?" an objective, scoreable question rather
than a matter of expert taste — what Husic & Pande (2018) call the field's move *"from an
art to a science."*

---

## 1. From PCA to TICA: same shape, different objective

**PCA** diagonalises the covariance matrix `C = ⟨x xᵀ⟩`. Its eigenvectors are the
orthogonal directions of maximal *variance*. Nothing in PCA knows about *time* — shuffle
your frames and PCA is unchanged.

**TICA** (time-lagged independent component analysis) diagonalises the same kind of object
but with a **time lag** built in. It solves the generalised eigenvalue problem

```
C(τ) v = λ C(0) v
```

where `C(0) = ⟨x(t) x(t)ᵀ⟩` is the ordinary covariance and
`C(τ) = ⟨x(t) x(t+τ)ᵀ⟩` is the **time-lagged** covariance at lag τ. The eigenvectors are
the linear combinations of your input features whose **autocorrelation decays slowest** —
the slow collective coordinates. The eigenvalues `λᵢ ∈ (0,1]` are autocorrelations, which
convert directly to timescales (Section 3).

Why this matters: PCA's high-variance direction and the system's slowest motion are
**different things**. Pérez-Hernández et al. (2013), who introduced TICA for MSM building,
state the problem as identifying "*the slow subspace*" of configuration space and note
that "*a method to identify this slow subspace exists in statistics: the time-lagged
independent component analysis (TICA)*." Their paper shows time-lagged analyses are
better suited to detect slow modes than PCA. (The underlying transform dates to Molgedey
& Schuster 1994; its use for molecular kinetics is Pérez-Hernández 2013 and, independently,
Schwantes & Pande 2013.)

> **Takeaway:** TICA is "PCA that maximises autocorrelation-time instead of variance."
> Same eigenvalue machinery; the lag time τ is the new ingredient. Everything downstream
> is built on these slow coordinates because the slow coordinates *are* the interesting
> kinetics.

In our runs the input "features" are inter-domain CA–CA distances (`interface_ca`,
`ca_stride25`, …). TICA combines them into the few coordinates along which MD1/MD2/MD3
rearrange most slowly.

---

## 2. What a Markov State Model actually is

An MSM is almost embarrassingly simple once you have good coordinates:

1. **Discretise** the (TICA) space into many small **microstates** (we use k-means,
   200 microstates) — like a histogram of where the molecule goes.
2. **Count transitions** at a fixed lag τ: how often does a trajectory sitting in state *i*
   land in state *j* exactly τ later? Normalise the counts row-wise → a **transition
   probability matrix** `T(τ)`.
3. **Diagonalise `T(τ)`**. The eigenvalue `λ₁ = 1` is equilibrium (stationary
   distribution). Each remaining `λᵢ < 1` is a **slow relaxation process**; its
   eigenvector says *which states exchange* in that process.

The central modelling assumption is the **Markov property**: the probability of the next
state depends only on the current state, not on the history of how you got there. Bowman,
Pande & Noé's textbook (2014) is the standard reference for this construction; Husic &
Pande (2018) give a readable history.

The reason MSMs are powerful: they stitch together many **short** trajectories into one
**long**-timescale kinetic model. You never need a single trajectory that crosses every
barrier — you only need enough short trajectories that *each* barrier is crossed *somewhere*
in the dataset. (This is exactly the regime CALVADOS CG runs live in: many short replicates.)

---

## 3. Eigenvalues → timescales, and the implied-timescale (ITS) test

Each MSM eigenvalue converts to an **implied timescale**:

```
tᵢ(τ) = −τ / ln λᵢ(τ)
```

This is the relaxation time of the *i*-th slow process. Here is the crucial logic, due to
Swope et al. (2004) and formalised by Prinz et al. (2011):

> **If the dynamics are truly Markovian at lag τ, the implied timescale `tᵢ(τ)` must be
> independent of τ.** A real relaxation time is a property of the molecule, not of the
> bookkeeping interval you chose. So you compute `tᵢ(τ)` at several lags and plot it. When
> the curve **flattens (plateaus)**, you have found a lag long enough for the Markov
> assumption to hold — and the plateau value is the physical timescale.

Why a plateau, mechanistically? At very short lags the discretised model still "sees"
within-state recrossings (non-Markovian memory) and *underestimates* the timescale.
As τ grows past the memory time, that error vanishes and the estimate levels off. Prinz
et al. (2011) prove the approximation error is bounded and "*can be made arbitrarily small
with surprisingly little effort*" — provided the slow coordinates are resolved and τ is
chosen on the plateau. Picking the lag time this way is standard practice (Husic & Pande
2018; Chodera & Noé 2014).

**What a non-plateau means.** If `tᵢ(τ)` keeps **rising** as you increase τ and never
flattens, the model is **not Markovian at any accessible lag**. Two usual causes:
(a) the discretisation/features don't capture the true slow coordinate, or (b) the
trajectories are **shorter than the process you're trying to resolve**, so the slow
eigenvalue can't be estimated. You cannot fix this by re-picking τ; you fix it with better
coordinates or more/longer sampling.

> ⚠️ **This is exactly what we see for PARP14** (Section 7): even at 100 ns/replicate the
> ITS climb essentially linearly with τ and never plateau.

---

## 4. The Chapman–Kolmogorov (CK) test

The ITS test checks timescales; the CK test checks the **transition matrix itself**,
self-consistently. A genuine Markov process obeys the **Chapman–Kolmogorov equation**:

```
T(nτ) ≈ T(τ)ⁿ
```

In words: *propagating the model n steps of length τ should agree with the model estimated
directly at the longer lag nτ.* If the one-step model already captures the dynamics, you
can iterate it and it predicts longer-time behaviour correctly. Prinz et al. (2011)
introduced this as one of the two standard MSM validations; you overlay `T(τ)ⁿ` (prediction)
on `T(nτ)` (re-estimation) for each macrostate and check they fall within statistical error.

We report a single scalar summary, the mean `|T_obs − T_pred|`, with the convention
`<0.05` excellent, `<0.10` good, `>0.20` invalid. **Important caveat:** the CK test only
probes Markovianity *at the lag you built the model at*. Passing CK at τ = 0.8 ns does
**not** rescue a non-converging ITS — they test different things. (A model can be locally
self-consistent at short lag yet still fail to expose a slow process it never sampled.)

---

## 5. The variational principle and VAMP-2: choosing features objectively

We tried six featurizations (`interface_ca`, `ca_stride25`, `linker_ca`, `orient`, `com`,
`inter_exposed`). How do we say one is *better* without circular reasoning? This is where
the field's key theoretical advance comes in.

**The variational principle for conformational dynamics** (Noé & Nüske 2013; Nüske et al.
2014). The true dynamics have exact eigenfunctions/eigenvalues. Noé & Nüske prove that any
*approximate* set of slow coordinates yields eigenvalue estimates that are **lower bounds**
to the true ones — the better your coordinates, the **larger** (closer to the true) the
estimated slow eigenvalues. Their construction "*is based on the maximization of a Rayleigh
coefficient*" that "*can be estimated from statistical observables … obtained from short
distributed simulations*." This turns model selection into an optimisation: **maximise the
captured slow-eigenvalue content.**

**VAMP-2 score** (Wu & Noé 2020; the closely related GMRQ of McGibbon & Pande 2015). The
VAMP-2 score is essentially the sum of squared singular values of the lag-τ propagator
approximated by your features — a single number measuring "how much slow dynamics did these
coordinates capture." Higher = better. McGibbon & Pande define a "*generalized matrix
Rayleigh quotient (GMRQ), which measures the ability of a rank-m projection operator to
capture the slow subspace of the system.*"

**Why cross-validation is mandatory.** The variational bound holds for the *true* matrix
elements, but we estimate them from finite, noisy data — so a flexible feature set can
**overfit**, scoring high on the data it was fit to while generalising poorly. McGibbon &
Pande (2015) show the GMRQ bound "*can be violated when the requisite matrix elements are
estimated subject to statistical uncertainty*" and that "*this overfitting can be detected
and avoided through cross-validation*." That is why we run `--cv` (5-fold): we score each
featurization on **held-out** trajectories. (Caution: Husic & Pande's companion note, 2017,
warns the *lag time* itself must not be chosen by VAMP score — only features/coordinates at
fixed lag — which is exactly how we use it.)

**How to read our VAMP-2 table.** A good featurization has a **high** score that is
**stable as the lag grows** — meaning the variance it captures belongs to genuinely slow
processes, not fast noise that decorrelates quickly. A feature that scores high at short lag
but **collapses** at longer lag was capturing fast motion.

---

## 6. PCCA+: from 200 microstates to a few macrostates

The MSM has 200 microstates — too many to call "representative states." **PCCA+** (Robust
Perron Cluster Cluster Analysis; Deuflhard & Weber 2005; Röblitz & Weber 2013) lumps them
into a few **metastable macrostates** using the sign structure of the leading MSM
eigenvectors. The idea: the *m* slowest eigenvectors are nearly constant *within* a
metastable basin and change sign *between* basins, so their pattern reveals the natural
basins. Röblitz & Weber (2013) prove PCCA+ "*always delivers an optimal fuzzy clustering for
nearly uncoupled … Markov chains*" by transforming "*dominant right eigenvectors of the
transition matrix into membership functions.*"

The choice of *m* (our `--k`) should be informed by a **spectral gap** — a clear drop
between slow eigenvalues (states that interconvert slowly, worth separating) and fast ones
(within-basin wiggle). When there's no gap — as for a floppy, diffusive system — the choice
of *k* is arbitrary and the "macrostates" are slices of a continuum rather than true basins.
This is why, for PARP14, we chose to pick states by **structural diversity from the TICA
landscape** rather than by a (non-existent) kinetic gap.

---

## 7. What our PARP14 numbers actually mean

Putting the theory to work on the **complete** macrodomain construct (`md_full`, FL
790–1388 incl. the 1194–1206 linker), 25 replicates **extended to 100 ns** (2.5 µs total):

**(a) Feature ranking is clear and stable.** Cross-validated VAMP-2 (lags 0.8 / 2.0 / 4.0 ns):

| feature | 0.8 ns | 2.0 ns | 4.0 ns | reading |
|---|---|---|---|---|
| **interface_ca** | 2.18 | 1.74 | 1.56 | high **and stable** → captures real slow modes ✅ |
| ca_stride25 | 4.12 | 1.71 | 0.68 | high then **collapses** → much of it is fast |
| linker_ca | 0.41 | 0.007 | 0.006 | fast noise |
| com / orient | ≤0.24 | ~0 | ~0 | uninformative |

So **`interface_ca` (domain-edge contacts) is the right coordinate** — the variational/VAMP
logic of Section 5 picks it unambiguously, and the choice is robust to trajectory length
(it won the 5 ns sweep too).

**(b) But the kinetics do *not* converge — and that is itself the result.** The implied
timescales for `interface_ca` (Section 3):

| lag τ (ns) | 0.8 | 2.0 | 4.0 | 8.0 | 16 | 32 |
|---|---|---|---|---|---|---|
| IC1 (ns) | 6.7 | 16.7 | 33 | 66 | 132 | 265 |

`tᵢ(τ)` rises almost perfectly **linearly** with τ (IC1 ≈ 8.3 × τ; equivalently the slow
eigenvalue sits at λ₁ ≈ 0.89 at *every* lag). By the ITS logic this is the textbook
signature of **no Markovian regime**: the model never "forgets," even at 32 ns. Extending
from 5 ns to 100 ns did not produce a plateau — it just pushed the apparent timescale up
in lock-step with the lag.

**What this tells us physically.** Combined with the earlier finding that the macrodomains
have **no persistent inter-domain interface** (inter-domain contact frequency ≈ 0.005;
COM–COM distances fluctuate by ~0.4–0.8 nm), the non-plateau says the MD1–MD2–MD3 module
behaves less like a few **metastable states separated by barriers** and more like a
**flexible, diffusive multidomain chain** (macrodomains on flexible linkers, closer to a
semi-flexible polymer with a continuum of slow modes than to a two- or three-state switch).
The MSM framework *assumes* metastable basins; when the landscape is shallow/diffusive, the
diagnostics correctly refuse to certify states — no amount of re-tuning τ or *k* manufactures
basins that aren't there.

**Consequence for the design pipeline.** The clusters we extract are therefore best read as
**structural representatives sampled along a continuum**, not kinetically-distinct metastable
states. That is fine — and honest — for the BindCraft goal of *enforcing a new interface*:
we deliberately pick **diverse poses** to clamp, and structural diversity (not kinetic
metastability) is what that needs. We just don't claim these are long-lived states the
protein "sits in."

> **Bottom line.** The featurization question has a clean answer (`interface_ca`). The
> metastable-states question has an honest negative answer: with these data the module
> looks diffusive, the ITS don't plateau, and "representative states" means *geometric
> diversity*, not *kinetic basins*. The tests did their job — they told us which question
> the data can and cannot answer.

---

## Glossary (quick reference)

- **Lag time τ** — the time interval at which transitions are counted / autocorrelation is
  measured. The one free knob that ITS/CK exist to validate.
- **TICA** — finds slowest-decorrelating linear coordinates (PCA-like, but maximises
  autocorrelation time).
- **Microstate / macrostate** — fine k-means cells / a few PCCA+ lumps of them.
- **Transition matrix `T(τ)`** — row-stochastic matrix of microstate→microstate
  probabilities at lag τ. Its eigen-decomposition = the kinetics.
- **Implied timescale `tᵢ = −τ/ln λᵢ`** — relaxation time of slow process *i*. Should be
  τ-independent (plateau) if Markovian.
- **Chapman–Kolmogorov test** — checks `T(nτ) ≈ T(τ)ⁿ`.
- **Variational principle** — better slow coordinates give larger (true-er) slow
  eigenvalues; basis of objective model selection.
- **VAMP-2 / GMRQ** — score for how much slow dynamics a featurization captures; use
  **cross-validated** to avoid overfitting.
- **PCCA+** — spectral clustering of the MSM into metastable macrostates.
- **Spectral gap** — a clear separation between slow and fast eigenvalues; its *absence*
  means no well-defined number of states.

---

## Sources

1. **Prinz, J.-H.; Wu, H.; Sarich, M.; Keller, B.; Senne, M.; Held, M.; Chodera, J. D.;
   Schütte, C.; Noé, F.** "Markov models of molecular kinetics: Generation and validation."
   *J. Chem. Phys.* **134**, 174105 (2011). — *Implied-timescale and Chapman–Kolmogorov
   tests; bounded approximation error.*
   https://doi.org/10.1063/1.3565032
2. **Pérez-Hernández, G.; Paul, F.; Giorgino, T.; De Fabritiis, G.; Noé, F.**
   "Identification of slow molecular order parameters for Markov model construction."
   *J. Chem. Phys.* **139**, 015102 (2013). — *Introduces TICA for MSMs; slow subspace.*
   arXiv:1302.6614 · https://doi.org/10.1063/1.4811489
3. **Schwantes, C. R.; Pande, V. S.** "Improvements in Markov State Model Construction
   Reveal Many Non-Native Interactions in the Folding of NTL9." *J. Chem. Theory Comput.*
   **9**, 2000–2009 (2013). — *Independent introduction of TICA to MSMs.*
   https://doi.org/10.1021/ct300878a
4. **Molgedey, L.; Schuster, H. G.** "Separation of a mixture of independent signals using
   time delayed correlations." *Phys. Rev. Lett.* **72**, 3634 (1994). — *Original
   time-lagged decorrelation transform.* https://doi.org/10.1103/PhysRevLett.72.3634
5. **Noé, F.; Nüske, F.** "A Variational Approach to Modeling Slow Processes in Stochastic
   Dynamical Systems." *Multiscale Model. Simul.* **11**, 635–655 (2013). — *Variational
   principle; Rayleigh-coefficient maximisation.* arXiv:1211.7103 ·
   https://doi.org/10.1137/110858616
6. **Nüske, F.; Keller, B. G.; Pérez-Hernández, G.; Mey, A. S. J. S.; Noé, F.**
   "Variational Approach to Molecular Kinetics." *J. Chem. Theory Comput.* **10**,
   1739–1752 (2014). https://doi.org/10.1021/ct4009156
7. **McGibbon, R. T.; Pande, V. S.** "Variational cross-validation of slow dynamical modes
   in molecular kinetics." *J. Chem. Phys.* **142**, 124105 (2015). — *GMRQ score;
   overfitting detected via cross-validation.* arXiv:1407.8083 ·
   https://doi.org/10.1063/1.4916292
8. **Wu, H.; Noé, F.** "Variational Approach for Learning Markov Processes from Time Series
   Data." *J. Nonlinear Sci.* **30**, 23–66 (2020). — *VAMP / VAMP-2 score.*
   arXiv:1707.04659 · https://doi.org/10.1007/s00332-019-09567-y
9. **Husic, B. E.; Pande, V. S.** "Note: MSM lag time cannot be used for variational model
   selection." *J. Chem. Phys.* **147**, 176101 (2017). — *Use VAMP for features at fixed
   lag, not to choose the lag.* https://doi.org/10.1063/1.5002086
10. **Deuflhard, P.; Weber, M.** "Robust Perron cluster analysis in conformation dynamics."
    *Linear Algebra Appl.* **398**, 161–184 (2005). — *PCCA+.*
    https://doi.org/10.1016/j.laa.2004.10.026
11. **Röblitz, S.; Weber, M.** "Fuzzy spectral clustering by PCCA+: application to Markov
    state models and data classification." *Adv. Data Anal. Classif.* **7**, 147–179 (2013).
    https://doi.org/10.1007/s11634-013-0134-6
12. **Swope, W. C.; Pitera, J. W.; Suits, F.** "Describing Protein Folding Kinetics by
    Molecular Dynamics Simulations. 1. Theory." *J. Phys. Chem. B* **108**, 6571–6581
    (2004). — *Early implied-timescale convergence test.*
    https://doi.org/10.1021/jp037421y
13. **Bowman, G. R.; Pande, V. S.; Noé, F. (eds.)** *An Introduction to Markov State Models
    and Their Application to Long Timescale Molecular Simulation.* Adv. Exp. Med. Biol.
    **797**, Springer (2014). — *Standard textbook.*
    https://doi.org/10.1007/978-94-007-7606-7
14. **Husic, B. E.; Pande, V. S.** "Markov State Models: From an Art to a Science."
    *J. Am. Chem. Soc.* **140**, 2386–2396 (2018). — *Readable history; the variational
    principle as the turning point to objective model selection.*
    https://doi.org/10.1021/jacs.7b12191
15. **Chodera, J. D.; Noé, F.** "Markov state models of biomolecular conformational
    dynamics." *Curr. Opin. Struct. Biol.* **25**, 135–144 (2014). — *Concise review.*
    https://doi.org/10.1016/j.sbi.2014.04.002
16. **Scherer, M. K. et al.** "PyEMMA 2: A Software Package for Estimation, Validation, and
    Analysis of Markov Models." *J. Chem. Theory Comput.* **11**, 5525–5542 (2015). —
    *Reference implementation of these methods.* https://doi.org/10.1021/acs.jctc.5b00743

*Direct quotations above are taken from the abstracts / texts of refs 1, 2, 5, 7, 11, and
14. Where a statement paraphrases rather than quotes, no quotation marks are used.*
