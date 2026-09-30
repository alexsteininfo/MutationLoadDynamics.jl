# Concepts

## A cell is a node in a binary tree

The unit of state is a [`NonMarkovCell`](@ref), carried as the `data` of a
[`BinaryNode`](@ref). It is immutable and holds five numbers:

| Field | Symbol | Meaning |
|:---|:---|:---|
| `id::Int64` | — | unique identifier, allocated in order of creation |
| `birthtime::Float64` | ``t_0`` | absolute simulation time at which this cell was born |
| `mutations::Int64` | ``j`` | driver mutations acquired **at this cell's own birth** |
| `total_mutations::Int64` | ``K`` | drivers on the whole path from the root, including ``j`` |
| `fitness::Float64` | ``f`` | the cell's fitness, inherited and then updated per mutation |

`mutations` is **local**: it is the Poisson draw made when this particular cell was
born, and nothing else. That is what makes a site-frequency spectrum a single
traversal: a mutation event sitting on one node is carried by exactly the leaves below
that node, so one post-order pass assigns every event to a frequency class.

`total_mutations` is the cell's whole burden, the parent's total plus its own
`mutations`. It is stored so that [`mutations_per_cell`](@ref), trajectory recording and
pairwise distances read a field instead of walking to the root. A hand-built tree must
keep it consistent; an `on_division` hook should change a cell only through
[`set_fitness!`](@ref).

`fitness` is **cumulative**: it is the parent's fitness with `fitness_update` applied
once per new mutation, stored outright. That makes the quantity the waiting-time
distributions actually need an ``O(1)`` field read at scheduling time rather than a walk
up the tree, which matters because scheduling happens twice per division for the life of
the run.

Ids are allocated strictly increasing, so an ancestor always has a smaller id than its
descendants. [`findMRCA`](@ref) and [`pairwisedistance`](@ref) both use that ordering to
climb toward a common ancestor without needing depths.

## The tree records the survivors

Cells that divide become internal nodes; cells alive now are leaves. A cell that dies is
removed from the tree, and so is any ancestor that is left childless by the removal —
that is [`prune_tree!`](@ref MutationLoadDynamics.prune_tree!)'s recursive step.

The tree you end up holding is therefore the **reduced tree of the survivors**, not a
complete record of everything that ever happened. Extinct side-lineages leave nothing
behind. This is a deliberate memory decision: under a supercritical process most of the
cells ever born are in lineages that die out, and keeping them would dominate storage
while contributing to no observable.

One structural consequence deserves attention, because it looks like corruption if you
meet it unprepared:

!!! note "Internal nodes can have one child"
    When a cell divides and one daughter's lineage later dies out entirely, the division
    itself still happened and stays in the tree — but the dead branch is gone, so that
    node now has a single child. The tree is a binary tree in which unary nodes are
    legitimate.

    This is the right behaviour, and statistics rely on it. A unary node is a real
    division, so counting it in [`leaf_depths`](@ref) correctly reports how many times a
    surviving cell's lineage actually divided. Collapsing unary nodes would turn depth
    into a count of *bifurcations that happen to have survived*, which is a property of
    the sampling, not of the cell. [`sample_leaves`](@ref) produces exactly the same
    shape, for exactly the same reason — see [Sampling](sampling.md).

## The population is the living cells only

[`Population`](@ref) holds alive cells in a `Dict{Int64, BinaryNode{NonMarkovCell}}`
keyed by id, so a birth or a death is ``O(1)``. It also carries the current simulation
time and the pending event queue.

```julia
pop = initialize_population(fitness_init = 1.0)   # one founder
pop = initialize_population(100; fitness_init = 1.0)  # 100 independent founders

popsize(pop)         # number of living cells
allcells(pop)        # Vector of their BinaryNodes
pop.t                # time of the most recently processed event
age(pop)             # === pop.t
```

The tree is reachable from the population and never stored separately: every living leaf
holds a `parent` chain back to its founder, so `getsingleroot(allcells(pop))` recovers
the root.

!!! warning "`initialize_population(N)` builds a forest, not a tree"
    The `N`-cell constructor seeds `N` **independent founders**, each the root of its own
    tree. That is the right model for a pre-existing pool of cells with no shared somatic
    ancestry, but several functions need a single root and say so:
    `getsingleroot` returns `nothing`, [`findMRCA`](@ref) returns `nothing`,
    [`clonal_mutations`](@ref) returns `0`, and [`sample_leaves`](@ref) /
    [`sample_trees`](@ref) throw with a message naming the number of roots.
    [`sitefrequencyspectrum`](@ref) handles the forest correctly, traversing each root.
    See [Populations with more than one root](statistics.md#Populations-with-more-than-one-root).

## Competing risks

At the moment a cell is born, **both** of its possible futures are drawn:

```math
T_\text{div} \sim \texttt{birth\_dist}(f), \qquad T_\text{die} \sim \texttt{death\_dist}(f)
```

and the earlier one is what will happen:

```math
t_\text{event} = t_0 + \min(T_\text{div},\, T_\text{die}), \qquad
\text{type} = \begin{cases} \text{division} & T_\text{div} \le T_\text{die} \\ \text{death} & \text{otherwise.}\end{cases}
```

This is [`schedule_cell!`](@ref MutationLoadDynamics.schedule_cell!), and it is exact —
no acceptance–rejection, no maximum-rate bound, no time discretisation. The losing draw
is simply discarded.

Two consequences worth stating explicitly:

**A cell's fate is fixed at birth and never revised.** Between a cell's birth and its
event, nothing can change what happens to it. Its own fitness cannot change (mutations
occur only at division, to daughters), which is why this is exact. But it also means the
simulator has no mechanism to invalidate a pending event when something *external*
changes — the population size, a neighbouring clone, a treatment. Density-dependent
rules can be written, but they are frozen as of each cell's birth; see
[Density-dependent and homeostatic growth](blocks.md#Density-dependent-and-homeostatic-growth).

**Both waiting times are always drawn.** `T_div` and `T_die` are both evaluated even
when one is certain to lose (a `Dirac` consumes no random numbers, every other family
does), so the random stream depends on the shape of the model, not
only on its outcomes. Two runs are draw-for-draw identical only when the distributions
are identical too — see [Reproducibility](#Reproducibility) below.

## The event queue

All pending events live in one global `BinaryMinHeap{CellEvent}` ordered by absolute
time. The loop is:

1. Pop the earliest event — ``O(\log N)``.
2. Advance `population.t` to that event's time.
3. **Division**: remove the parent from the living registry, create two daughters
   (drawing their mutations and fitness), fire `on_division` if given, then schedule
   both daughters.
4. **Death**: prune the cell from the tree and remove it from the registry.
5. Test `stopfunction`. If it is true, or the heap is empty, stop.

Growing from one cell to ``N`` costs ``O(N \log N)``.

The heap maintains a strict invariant: **exactly one pending event per living cell**. A
division pops one event and pushes two while the population grows by one; a death pops
one and pushes none while the population shrinks by one. `simulate!` checks this
invariant when it resumes on a population that already carries a queue, and throws
naming both counts if you have edited `population.cells` by hand —
[`reset_schedule!`](@ref) is the recovery.

### Why this is not Gillespie

The Gillespie algorithm draws a time to the *next event anywhere in the population* from
a single exponential whose rate is the sum over all cells, then picks which cell it was.
That factorisation is only valid because the exponential is memoryless: after an event
elsewhere, every other cell's remaining waiting time has the same distribution it had
before.

Drop the exponential and the factorisation fails. Cell ages then matter, and a correct
algorithm has to track the residual waiting time of every cell individually. Pre-drawing
each cell's own absolute event time and ordering them in a heap does exactly that — it
is the Next Reaction Method for a non-Markovian process — and it is exact for *any*
distributions, including the exponential.

So the exponential is not a separate code path here; it is the special case in which
this algorithm and Gillespie agree. See
[Exponential — the Markov special case](blocks.md#Exponential-—-the-Markov-special-case).

## Why not the exponential?

The exponential's defining property is that a cell's chance of dividing in the next
instant does not depend on how long it has already been alive. Mammalian cells do not
behave that way: a cell that has just divided cannot divide again for the length of a
cell cycle — 8–12 h for fast tumour cells, days for normal somatic tissue. The
exponential puts substantial probability mass in that refractory window, and gives
consistently poor fits to measured proliferation data (Zilman et al. 2010).

The natural replacement is the Gamma, which retains a closed-form mean and variance
while adding exactly one dimensionless parameter, the shape ``k``:

```math
\mathbb{E}[T] = k\theta, \qquad \mathrm{Var}[T] = k\theta^2, \qquad \mathrm{CV} = \frac{1}{\sqrt{k}}
```

``k = 1`` recovers the exponential; ``k \to \infty`` approaches a deterministic clock.
Reported coefficients of variation for interdivision time:

| Cell type | Typical CV | ``k = 1/\mathrm{CV}^2`` |
|:---|:---|:---|
| Well-regulated somatic (fibroblasts, epithelial) | 0.10–0.20 | 25–100 |
| Fast cancer cells (leukaemia, ovarian carcinoma) | 0.30–0.45 | 5–11 |
| Stimulated lymphocytes | 0.50–0.70 | 2–4 |

``k = 5`` (CV ≈ 0.45) is a reasonable default for rapidly cycling tumour cells.

The best-fitting empirical family is not the Gamma but the **exponentially modified
Gamma** — a Gamma convolved with an exponential, where the exponential captures the
stochastic G1 restriction-point wait and the Gamma the more deterministic S/G2/M phases.
It outperforms the plain Gamma across 77 published datasets and 16 cell types (Golubev
2016). Nothing in this package prevents you from using it: `birth_dist` returns any
`Distribution`, and the simulator only ever calls `rand(rng, d)` on it. Distributions.jl
has no built-in sum of a Gamma and an exponential, but a sampler is three lines:

```julia
struct ExpModGamma <: ContinuousUnivariateDistribution
    k::Float64; θ::Float64; λ::Float64       # Gamma shape and scale, exponential mean
end
Base.rand(rng::AbstractRNG, d::ExpModGamma) =
    rand(rng, Gamma(d.k, d.θ)) + rand(rng, Exponential(d.λ))

birth_dist = f -> ExpModGamma(k, θ / f, λ / f)
```

The Gamma is offered as the tractable default, not as a claim about biology.

!!! note "\"Mean equals variance\" is not a simplification here"
    Imposing ``\mathrm{Var}[T] = \mathbb{E}[T]`` on a Gamma forces ``k\theta^2 = k\theta``,
    hence ``\theta = 1``; combined with ``\mathbb{E}[T_\text{div}] = 1/b`` it gives
    ``k = 1/b``, which at ``b = 1`` is the exponential. The constraint collapses the model
    back to Gillespie rather than reducing it.

### What non-exponential timing changes

Cell-cycle variability is not a nuisance parameter — it changes the growth rate itself.
For a population where cells divide with density ``g`` and survive death with
probability ``S_d``, the Malthusian rate ``r`` solves the Euler–Lotka equation

```math
2\int_0^\infty e^{-rt}\, g(t)\, S_d(t)\, \mathrm{d}t = 1 .
```

For exponentials with rates ``b`` and ``d`` this gives ``2b/(r + b + d) = 1``, that is
``r = b - d``, the familiar answer. For pure birth with ``T_\text{div} \sim
\mathrm{Gamma}(k, 1/(kb))`` — mean division time ``1/b`` for every ``k`` — it gives

```math
r = k\,b\left(2^{1/k} - 1\right),
```

which decreases monotonically from ``r = b`` at ``k = 1`` to ``r = b\ln 2 \approx 0.693\,b``
as ``k \to \infty``. Holding the mean cell-cycle length fixed, a **more variable** cell
cycle grows the population **faster**, because early-dividing lineages compound. Any
calibration that matches a mean division time while changing ``k`` is therefore not
holding growth fixed, and any inference that assumes ``r = b - d`` while the data were
generated at ``k = 5`` is biased by about 25%.

The same decoupling is what makes the package's outputs informative: under exponential
timing, elapsed time and division count are hard to tell apart, and under Gamma timing
they are not. That is why [`leaf_depths`](@ref) (divisions) and coalescence times (real
time) are both retained as separate observables.

## Reproducibility

`simulate!` takes an `AbstractRNG` and consumes it in a fixed order. A run is
reproducible from `(initial population, block, seed)`, and a chained run is identical
draw-for-draw to the equivalent uninterrupted one — the test suite asserts this bit for
bit across seeds.

The order of consumption per division is: daughter 1's mutation count and its effect
sizes, then daughter 2's, then `on_division` if it draws, then daughter 1's two waiting
times, then daughter 2's two. Anything that changes how many values a step consumes
shifts everything after it, so the following are worth knowing:

- `Poisson(0.0)` **does** consume a draw. A neutral run (`ν = 0`) still burns one value
  per daughter, so it does not share a stream with a run that skips the mutation step.
- `Dirac(x)` consumes **nothing**. Swapping a `Dirac` driver effect for an `Exponential`
  one therefore shifts the stream, on top of changing the model.
- Both waiting times are always drawn, even when the death distribution is
  `Dirac(Inf)`.
- A fresh schedule (a new population, [`reset_schedule!`](@ref), an extinction restart)
  visits cells in increasing id order. A cell older than zero at that moment may need
  several pairs of draws, because its next event is conditioned on nothing having
  happened to it yet.

**Across Julia versions.** A seed reproduces a run only on the same Julia and
Distributions.jl versions: `MersenneTwister` streams and several samplers have changed
between releases. Passing a `StableRNG` from StableRNGs.jl removes the first source of
drift but not the second. Leaf sampling is the exception: it depends only on
`StableRNG` and is stable across versions (see [Sampling](sampling.md)).

The package holds no global mutable state. One population, one block and one closure per
simulation is safe to run under `Threads.@threads`, provided each thread has its own
rng. [`sample_leaves`](@ref) goes further and builds its own rng from a caller-supplied
seed, so draws are independent by construction.

## References

- **Golubev A (2016)** — the exponentially modified Gamma fits mammalian cell-cycle data
  better than Gamma or lognormal alone; 77 datasets, 16 cell types.
  *J. Theor. Biol.* 393, 203–217. PMID 26780652.
  DOI: [10.1016/j.jtbi.2015.12.027](https://doi.org/10.1016/j.jtbi.2015.12.027)
- **Zilman A, Ganusov VV, Perelson AS (2010)** — Gamma shapes ``k = 2``–3 are needed to fit
  CD4⁺ T-cell proliferation; ``k = 1`` fits poorly in every condition tested.
  *PLoS ONE* 5(9): e12775. PMID 20941358.
  DOI: [10.1371/journal.pone.0012775](https://doi.org/10.1371/journal.pone.0012775)
- **Hahn GM (1966)** — the CV of interdivision time governs desynchronisation kinetics;
  ``k = 1`` predicts immediate loss of synchrony, contradicted by experiment.
  *Biophys. J.* 6(2), 197–207. PMID 5963460.
  DOI: [10.1016/S0006-3495(66)86656-0](https://doi.org/10.1016/S0006-3495(66)86656-0)
- **Sandler O et al. (2015)** — cousin–cousin correlations in cell-cycle duration dominate
  over mother–daughter ones; i.i.d. Gamma underestimates lineage-level clustering.
  *Nature* 519, 422–425. PMID 25762143.
  DOI: [10.1038/nature14318](https://doi.org/10.1038/nature14318)
- **Yanagisawa M et al. (1985)** — direct per-phase CV measurement in CHO cells; G1 is the
  most variable phase, S/G2/M the most constrained.
  *Cytometry* 6(6), 550–558. PMID 4064838.
  DOI: [10.1002/cyto.990060609](https://doi.org/10.1002/cyto.990060609)
- **Chiorino G et al. (2001)** — Gamma-based desynchronisation for cancer cell lines;
  CV ≈ 0.2–0.4 for IGROV1 ovarian carcinoma and MOLT4 leukaemia.
  *J. Theor. Biol.* 208(2), 185–199. PMID 11162063.
  DOI: [10.1006/jtbi.2000.2213](https://doi.org/10.1006/jtbi.2000.2213)
- **Fennell DA et al. (2005)** — apoptosis kinetics modelled with an exponential
  time-to-MOMP at population level; the exponential is more defensible for death than for
  division.
  *Apoptosis* 10(3), 517–530. PMID 15843905.
  DOI: [10.1007/s10495-005-0818-2](https://doi.org/10.1007/s10495-005-0818-2)
