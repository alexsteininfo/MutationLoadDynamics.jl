# Model Description

## Overview

We simulate somatic cell populations using a **non-Markovian birth–death process** implemented
with a global min-heap event queue. Each cell independently pre-schedules its next division or
death event by drawing waiting times from fitness-dependent Gamma distributions. The full binary
lineage tree is preserved throughout.

---

## Cell State

Each cell carries:

| Field | Symbol | Description |
|---|---|---|
| `birthtime` | $t_0$ | Absolute time of birth (or simulation start) |
| `mutations` | $j$ | Driver mutations acquired at this division (Poisson draw) |
| `fitness` | $f$ | Cumulative fitness (parent fitness updated per driver mutation) |

---

## Waiting-Time Distributions

We use **Gamma** distributions for both division and death times. The Gamma is the natural
generalization of the exponential that retains analytical tractability while capturing the
non-Markovian, clock-like nature of the mammalian cell cycle.

For a cell with fitness $f$, waiting times are drawn independently at the moment of birth:

$$T_\text{div} \sim \text{Gamma}(k,\, \theta_b / f), \qquad T_\text{die} \sim \text{Gamma}(k,\, \theta_d)$$

where:

- $k$ — shape parameter (shared for birth and death), controls timing variability
- $\theta_b = 1/(k b)$ — birth scale; mean division time $= k \cdot \theta_b / f = 1/(b f)$, inversely proportional to fitness
- $\theta_d = 1/(k d)$ — death scale; mean death time $= 1/d$, fitness-independent
- $b, d$ — baseline birth and death rates

The Gamma mean, variance, and coefficient of variation are:

$$\mathbb{E}[T] = k\theta, \qquad \mathrm{Var}[T] = k\theta^2, \qquad \mathrm{CV} = \frac{1}{\sqrt{k}}$$

The shape $k$ is the sole dimensionless parameter governing timing noise. Larger $k$ gives more
clock-like, less variable divisions.

### Why not Exponential?

The exponential distribution ($k = 1$) is the **Markovian limit** (Gillespie algorithm). It is
memoryless: the probability of dividing in the next instant is independent of how long the cell
has already been alive. This is biologically unrealistic — mammalian cells have a minimum cell
cycle duration of 8–12 h (fast cancer cells) to days (normal somatic cells), imposing a hard
refractory period incompatible with the exponential. Models with $k = 1$ give consistently poor
fits to experimental proliferation data (Zilman et al. 2010, PMID 20941358).

### Why not EMGD?

The most accurate empirical distribution for mammalian interdivision times is the
**Exponentially Modified Gamma Distribution** (EMGD = Gamma $*$ Exponential): the exponential
component captures the stochastic G1 restriction-point wait; the Gamma captures the
deterministic S/G2/M phases (Golubev 2016, PMID 26780652). EMGD outperforms plain Gamma on 77
published datasets across 16 cell types. We use plain Gamma as a tractable approximation that
captures the essential non-Markovian character without adding a third free parameter.

### Choice of $k$

| Cell type | Typical CV | $k = 1/\mathrm{CV}^2$ |
|---|---|---|
| Well-regulated somatic (fibroblasts, epithelial) | 0.10–0.20 | 25–100 |
| Fast cancer cells (leukemia, ovarian carcinoma) | 0.30–0.45 | 5–11 |
| Stimulated lymphocytes | 0.50–0.70 | 2–4 |

**Default: $k = 5$** (CV ≈ 0.45), appropriate for rapidly cycling cancer cells.
Empirical CV values from Yanagisawa et al. 1985 (PMID 4064838), Chiorino et al. 2001
(PMID 11162063), Zilman et al. 2010 (PMID 20941358).

### Note on "mean = variance"

Setting $\mathrm{Var}[T] = \mathbb{E}[T]$ for a Gamma requires $k\theta^2 = k\theta$, hence
$\theta = 1$. Combined with the constraint $\mathbb{E}[T_\text{div}] = 1/b$ we get $k = 1/b$.
For $b = 1$ this gives $k = 1$ — the exponential. Mean $=$ variance is therefore not a useful
simplification; it collapses back to Gillespie.

### Death times when $d = 0$

When $d = 0$ (pure birth, no apoptosis), we set $\theta_d = 10^8$ so death events essentially
never fire. For $d > 0$ we use $\theta_d = 1/(kd)$ with the same shape $k$ as birth.

---

## Competing Risks

At each birth event, both $T_\text{div}$ and $T_\text{die}$ are sampled. The cell's fate is
determined by the minimum:

$$t_\text{event} = \min\!\bigl(t_0 + T_\text{div},\; t_0 + T_\text{die}\bigr)$$

$$\text{event type} = \begin{cases} \text{division} & T_\text{div} \leq T_\text{die} \\ \text{death} & T_\text{die} < T_\text{div} \end{cases}$$

This is mathematically exact — no acceptance-rejection step is needed. For the exponential
special case ($k = 1$), the division probability simplifies to:

$$P(\text{divide}) = \frac{bf}{bf + d}$$

For general $k$, this probability is a function of the ratio $\theta_b f / \theta_d = d/(bf)$
and must be evaluated numerically.

---

## Algorithm

All living cells' next events are stored in a **global min-heap** ordered by absolute event time.
At each step:

1. Pop the earliest event — $O(\log N)$
2. Set simulation time $t \leftarrow t_\text{event}$
3. **Division**: remove parent from alive-cell registry; create two daughters with new fitness and
   mutation counts; draw their waiting times; push to heap
4. **Death**: prune cell from tree; remove from alive-cell registry

Total complexity to grow from 1 to $N$ cells (birth-dominated): $O(N \log N)$.

---

## Driver Mutation Model

At each division, each daughter independently:

1. Draws $j \sim \mathrm{Poisson}(\nu)$ — number of new driver mutations
2. For each mutation: draws $\delta \sim \pi_\delta$, then updates fitness:

$$f \;\leftarrow\; \text{fitness\_update}(f,\, \delta)$$

3. Stores $j$ (for SFS traversal) and the final $f$

Two natural update rules:

$$f \leftarrow f + \delta \qquad \text{(additive)}$$
$$f \leftarrow f \cdot (1 + \delta) \qquad \text{(multiplicative)}$$

For the **neutral case** ($\nu = 0$): no mutations are drawn, all cells maintain $f = f_0 = 1$,
and division times are identically distributed.

---

## Measurements

The simulation records a `TrajectoryPoint` every $\Delta t_\text{traj}$ time units:

| Field | Description |
|---|---|
| $t$ | Simulation time |
| $N$ | Population size |
| $\bar{f}$ | Mean fitness across alive cells |
| $\sigma^2_f$ | Variance of fitness |
| $\bar{k}$ | Mean total driver mutations per cell (summed to root) |
| $\sigma^2_k$ | Variance of driver count |

Full snapshots (tree, fitness distribution, driver SFS) can be recorded at triggers
`AtEnd`, `AtTime(t)`, `AtPopSize(N)`.

---

## Neutral Growth Parameters (`growth_neutral`)

| Parameter | Value | Notes |
|---|---|---|
| $b$ | 1.0 | Mean division time = 1.0 (time unit) |
| $d$ | 0.0, 0.5, 0.9 | Pure growth / moderate death / near-critical |
| $k$ | 5 | CV ≈ 0.45; cancer-cell regime |
| $\nu$ | 0 | Neutral — no driver mutations |
| $N_\text{target}$ | 100, 1 000 | Stop condition |
| $N_\text{sims}$ | 50 | Independent replicates per parameter set |

The net Malthusian growth rate under the Markovian approximation is $r \approx b - d$. In the
exact non-Markovian Gamma model, $r$ is the same for fixed means but the variance of offspring
number per unit time is reduced (compared to exponential), leading to lower effective drift and
slower coalescence in large populations.

---

## Selection Models

All four selection models share the same background parameters as the neutral model except that
driver mutations are now active (ν > 0) and fitness is no longer constant.

| Parameter | Value | Notes |
|---|---|---|
| $b$ | 1.0 | Baseline birth rate at fitness $f = 1$ |
| $d$ | 0.5 | Death rate (moderate; net growth rate ≈ 0.5 at baseline) |
| $k$ | 5 | Gamma shape for both birth and death waiting times |
| $\nu$ | 0.2 | Expected driver mutations per daughter cell (Poisson) |
| $N_\text{target}$ | 1 000 | Stop condition |
| $N_\text{sims}$ | 50 | Independent replicates per (model, $s$) pair |
| $s$ | 0.0 – 0.5 | Selection coefficient swept in {0.0, 0.1, 0.2, 0.3, 0.4, 0.5} |

Each driver mutation is drawn from a distribution specified per model. After drawing increment
$\delta$, the cell's fitness $f$ is updated according to a model-specific rule. The updated $f$
then sets the birth waiting-time distribution for all future events of that cell:

$$T_\text{div} \sim \text{Gamma}\!\left(k,\, \frac{1}{k \cdot f}\right), \quad \mathbb{E}[T_\text{div}] = \frac{1}{f}$$

---

### Model 1 — Additive Selection, Fixed Effect (`growth`)

**Fitness update rule:**
$$f \;\leftarrow\; f + s$$

Each driver mutation adds a fixed increment $s$ to the current fitness. After $n$ driver
mutations the cell's fitness is $f = 1 + n \cdot s$. There is no cap — fitness grows without
bound as mutations accumulate.

**Driver distribution:** $\delta = s$ (Dirac point mass; no randomness in effect size).

**Biological interpretation:** Models a scenario where every driver mutation confers the same
proliferative advantage, analogous to gain-of-function mutations in a single signalling pathway.
The linear accumulation rule corresponds to additive effects on the birth rate: each mutation
independently shortens mean division time by a fixed fraction.

**Key behaviour:**
- At $s = 0$: identical to neutral growth (no fitness change).
- At large $s$: early driver mutations sweep rapidly, compressing coalescence times.
- Fitness distribution at endpoint is approximately normal (sum of Poisson many $s$ increments).

**Output files:** `data/raw/growth/growth_s{s}_k5.0.jls`

---

### Model 2 — Max-Random Selection (`growth_rand`)

**Fitness update rule:**
$$f \;\leftarrow\; \max\!\left(f,\; 1 + \delta\right), \quad \delta \sim \text{Exp}(s)$$

Each driver mutation draws a candidate new fitness from $1 + \text{Exp}(s)$. The cell keeps the
_better_ of its current fitness and the new draw. Fitness can only increase; once a high-fitness
state is reached it is permanent.

**Driver distribution:** $\delta \sim \text{Exp}(s)$ (mean $s$, heavy right tail).  
For $s = 0$: $\text{Dirac}(0)$ is used so no fitness change occurs.

**Biological interpretation:** Models clonal competition where only the dominant pathway
determines division speed — acquiring a second oncogene in a pathway already activated does
nothing. Captures "winner-take-all" epistasis: subsequent mutations are only beneficial if they
exceed the current fitness level. The exponential distribution reflects uncertainty about the
magnitude of driver effects while preserving the observation that large effects are rare.

**Key behaviour:**
- Fitness converges to the order statistic $\max\{\delta_1, \ldots, \delta_n\}$ of the draws
  made across all mutations in the lineage; expected maximum of $n$ Exp($s$) draws scales as
  $s \ln n$.
- Additional mutations beyond the first beneficial one contribute little unless they exceed the
  current maximum — creating diminishing returns and a characteristic fitness plateau.
- The SFS is expected to show a more pronounced sweep signature than the additive model at the
  same mean effect $s$.

**Output files:** `data/raw/growth_rand/growth_rand_s{s}_k5.0.jls`

---

### Model 3 — Multiplicative Selection, Fixed Effect (`growth_multi`)

**Fitness update rule:**
$$f \;\leftarrow\; \min\!\left(f \cdot \left(1 + s \cdot \left(1 - \frac{f}{M}\right)\right),\; M\right), \quad M = 10$$

Each driver mutation multiplies the current fitness by a logistic gain factor
$(1 + s(1 - f/M))$. The factor $(1 - f/M)$ reduces the per-mutation gain as fitness approaches
the cap $M$, preventing unbounded acceleration. This is the discrete-time analogue of logistic
growth in fitness space.

**Driver distribution:** $X = 1$ (Dirac; fixed multiplicative magnitude, plugged in as $\delta = 1$).

**Biological interpretation:** Models tumour suppressor loss or oncogene dosage effects where
each additional mutation multiplies the division rate by a factor $(1 + s)$ reduced by the
fraction of remaining "room to grow" toward a physiological cap. The cap $M$ represents a
maximum achievable birth rate set by nutrient supply, checkpoint constraints, or physical
crowding. The logistic factor ensures diminishing returns and a stable attractor at $f = M$.

**Key behaviour:**
- At low $f$ (early tumour): mutations have near-full multiplicative effect.
- Near $M$: gain per mutation → 0, plateau emerges, further mutations are nearly neutral.
- Mean fitness at endpoint follows a logistic approach to $M$; variance narrows as $f \to M$.
- Unlike additive models, mutation load (number of drivers) matters non-linearly.

**Parameters:** $M = 10$, $s$ swept over {0.0, 0.1, 0.2, 0.3, 0.4, 0.5}.

**Output files:** `data/raw/growth_multi/growth_multi_s{s}_k5.0.jls`

---

### Model 4 — Multiplicative Selection, Random Effect (`growth_multi_rand`)

**Fitness update rule:**
$$f \;\leftarrow\; \min\!\left(f \cdot \left(1 + s \cdot X \cdot \left(1 - \frac{f}{M}\right)\right),\; M\right), \quad X \sim \text{Exp}(1),\; M = 10$$

Identical to Model 3 except the fixed effect magnitude $X = 1$ is replaced by a random draw
$X \sim \text{Exp}(1)$. The expected per-mutation gain is $\mathbb{E}[sX(1 - f/M)] = s(1 - f/M)$,
equal to the fixed-effect model, but with additional variance from $X$.

**Driver distribution:** $X \sim \text{Exp}(1)$ (mean 1, heavy right tail).

**Biological interpretation:** Extends Model 3 to the realistic setting where driver mutations
vary in their functional impact. Some mutations strongly activate a pathway ($X \gg 1$); others
are barely consequential ($X \ll 1$). The exponential distribution is a maximum-entropy prior
for positive random magnitudes given a fixed mean, and matches empirically measured distributions
of fitness effects of beneficial mutations in microbial evolution experiments.

**Key behaviour:**
- Overdispersion relative to Model 3: variance in $f$ at endpoint is higher.
- Right-tail draws ($X \gg 1$) can push fitness to $M$ in a single mutation when $s$ is large
  enough, creating a bimodal final fitness distribution (cells near $M$ vs. cells still near 1).
- At $s = 0$: $X$ is drawn but $s \cdot X = 0$, so no fitness change regardless of $X$.
- Comparison with Model 3 isolates the effect of effect-size heterogeneity at fixed mean $s$.

**Parameters:** $M = 10$, $s$ swept over {0.0, 0.1, 0.2, 0.3, 0.4, 0.5}.

**Output files:** `data/raw/growth_multi_rand/growth_multi_rand_s{s}_k5.0.jls`

---

## References

- **Golubev A (2016)** — Exponentially Modified Gamma distribution fits mammalian cell cycle
  data better than Gamma or lognormal alone; analysis of 77 datasets across 16 cell types.
  *J. Theor. Biol.* 393, 203–217. PMID 26780652.
  DOI: [10.1016/j.jtbi.2015.12.027](https://doi.org/10.1016/j.jtbi.2015.12.027)

- **Zilman A, Ganusov VV, Perelson AS (2010)** — Gamma distributions ($k=2$–3) needed to fit
  CD4⁺ T cell proliferation; $k=1$ (exponential) gives poor fits in all conditions tested.
  *PLoS ONE* 5(9): e12775. PMID 20941358.
  DOI: [10.1371/journal.pone.0012775](https://doi.org/10.1371/journal.pone.0012775)

- **Hahn GM (1966)** — CV of interdivision times governs desynchronization kinetics;
  exponential ($k=1$) predicts immediate loss of synchrony, contradicted by experiment.
  *Biophys. J.* 6(2), 197–207. PMID 5963460.
  DOI: [10.1016/S0006-3495(66)86656-0](https://doi.org/10.1016/S0006-3495(66)86656-0)

- **Sandler O et al. (2015)** — Cousin-cousin correlations in cell cycle duration dominate
  over mother-daughter; iid Gamma underestimates lineage-level clustering.
  *Nature* 519, 422–425. PMID 25762143.
  DOI: [10.1038/nature14318](https://doi.org/10.1038/nature14318)

- **Yanagisawa M et al. (1985)** — Direct measurement of CV for each cell cycle phase in CHO
  cells; G1 most variable, S/G2/M constrained.
  *Cytometry* 6(6), 550–558. PMID 4064838.
  DOI: [10.1002/cyto.990060609](https://doi.org/10.1002/cyto.990060609)

- **Chiorino G et al. (2001)** — Gamma-based desynchronization model for cancer cell lines;
  CV ≈ 0.2–0.4 for IGROV1 ovarian carcinoma and MOLT4 leukemia.
  *J. Theor. Biol.* 208(2), 185–199. PMID 11162063.
  DOI: [10.1006/jtbi.2000.2213](https://doi.org/10.1006/jtbi.2000.2213)

- **Fennell DA et al. (2005)** — Apoptosis kinetics modeled with exponential time-to-MOMP
  at population level; exponential is more defensible for death than for division.
  *Apoptosis* 10(3), 517–530. PMID 15843905.
  DOI: [10.1007/s10495-005-0818-2](https://doi.org/10.1007/s10495-005-0818-2)
