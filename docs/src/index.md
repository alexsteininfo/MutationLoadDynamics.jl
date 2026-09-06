# MutationLoadDynamics.jl

Forward simulation of a somatic cell population under **non-Markovian** birth–death
dynamics, with per-cell fitness and the complete lineage tree retained.

Division and death waiting times are drawn from arbitrary user-supplied distributions
that may depend on the cell's own fitness. Nothing in the package assumes the
exponential, so the cell cycle can be given a realistic refractory period — a minimum
duration below which a cell essentially never divides — which a memoryless model cannot
express. Each driver mutation draws its own fitness increment, so every cell carries an
individual fitness rather than inheriting a shared subclone label.

What comes back is the lineage tree itself: every division that led to a surviving cell,
with its time, its mutation count and its fitness. Site-frequency spectra, mutational
burdens, coalescence times and divisional depths are all computed from that tree, on the
full population or on a uniform sample of it.

This package is a **simulator**. It does not infer trees, fit parameters, or compute
distances between populations.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/alexsteininfo/MutationLoadDynamics.jl")
```

Dependencies are deliberately light — `Random`, `Distributions`, `StatsBase`,
`Statistics`, `AbstractTrees`, `DataStructures` — with no plotting or data-frame stack,
so depending on this package for its types and its trees stays cheap.

## Quickstart

```julia
using MutationLoadDynamics, Distributions, Random, Statistics

rng = MersenneTwister(42)

# One founding cell at fitness 1.
pop = initialize_population(fitness_init = 1.0)

# What the cells do.
block = NonMarkovBlock(
    birth_dist     = f -> Gamma(5.0, 1 / (5 * f)),  # mean division time 1/f, CV = 1/√5
    death_dist     = f -> Gamma(5.0, 1 / (5 * 0.3)),# mean death time 1/0.3, fitness-free
    stopfunction   = pop -> popsize(pop) >= 10_000, # grow to 10 000 cells
    driver_dist    = Exponential(0.05),             # each driver's effect size δ
    fitness_update = (f, δ) -> f + δ,               # additive selection
    ν              = 0.2,                           # mean drivers per daughter
    restart_on_extinction = true,                   # retry if the founder line dies
)

# What to record while it runs.
spec = MeasurementSpec(
    trajectory_dt     = 0.5,
    snapshot_triggers = [AtEnd()],
    snapshot_stats    = [SFS(), FitnessDistribution(), DriversPerCell()],
)
acc = MeasurementAccumulator(spec)

simulate!(pop, block, rng; accumulator = acc)
m = finalize_measurements(acc)

popsize(pop)                      # 10000
pop.t                             # simulation time reached
mean(fitness_per_cell(pop))       # mean fitness of the living cells
mean_k(pop)                       # mean driver burden per cell

m.trajectory                      # Vector{TrajectoryPoint}, one per 0.5 time units
m.snapshots[1].sfs                # driver site-frequency spectrum at the end

# The tree is the real output.
root = getsingleroot(allcells(pop))
mutations_per_cell(root)          # burden of every living cell
leaf_depths(root)                 # divisions from founder to each living cell
branch_spectrum(root)             # topology, separated from the mutation rate

# Sequence 1 000 of the 10 000 cells, the way an experiment would.
out = sample_trees(pop, SamplingSpec(1_000); seed = UInt64(0xBEEF))
sitefrequencyspectrum(out.samples[1].root, 1_000)
```

## The shape of a run

Three objects, each answering one question.

| Object | Question | Page |
|:---|:---|:---|
| [`Population`](@ref) | what exists now, and at what time | [Concepts](concepts.md) |
| [`NonMarkovBlock`](@ref) | what the cells do, and until when | [The simulation block](blocks.md) |
| [`MeasurementSpec`](@ref) | what to write down while it happens | [Output](output.md) |

[`simulate!`](@ref) advances the first under the rules of the second, recording through
the third. It mutates the population in place and returns it, so successive blocks chain
onto one lineage tree — see [Chaining blocks](blocks.md#Chaining-blocks).

## Where to go next

- [Concepts](concepts.md) — the cell, the tree, competing risks, and the event queue
  that makes non-exponential timing exact rather than approximate.
- [The simulation block](blocks.md) — every field of [`NonMarkovBlock`](@ref); how to
  write deterministic, exponential (Markov), Gamma or heavy-tailed waiting times; how to
  choose a stop condition; the division and restart hooks; chaining.
- [Mutations and selection](selection.md) — the driver channel, and the selection modes
  it expresses: neutral, additive, multiplicative, winner-takes-all, capped, deleterious.
- [Output](output.md) — trajectories, snapshot triggers, and reading the tree directly.
- [Tree statistics](statistics.md) — spectra, burdens, depths, distances, coalescence
  times, and which of them are co-indexed with which.
- [Sampling](sampling.md) — drawing `n` of `N` cells and the induced tree.
- [Limitations and open questions](limitations.md) — what is deliberately out of scope,
  and what is still undecided.

The scientific case for non-exponential cell-cycle timing, with references, is in
[Why not the exponential?](concepts.md#Why-not-the-exponential?).

## Related packages

For the **low** driver-mutation-rate regime, where clone-level bookkeeping is
meaningful and passenger mutations dominate the observable spectrum, see
[`BirthDeathMutation`](https://github.com/alexsteininfo/BirthDeathMutation).

To lay somatic **copy-number** alterations along a tree this package produces, see
[`CopyNumberEvolution.jl`](https://github.com/alexsteininfo/CopyNumberEvolution.jl),
which consumes a `MutationLoadDynamics.jl` tree directly.

The analyses, figures and theory that use this simulator live in
[`gITH-nonMarkovian`](https://github.com/alexsteininfo/gITH-nonMarkovian); this
repository holds only the simulation code.
