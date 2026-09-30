# MutationLoadDynamics.jl

Forward simulation of a somatic cell population under **non-Markovian** birth–death
dynamics, with per-cell fitness and the complete lineage tree retained.

Division and death waiting times are drawn from arbitrary user-supplied distributions
that may depend on the cell's own fitness. Nothing assumes the exponential, so the cell
cycle can be given a realistic refractory period — a minimum duration below which a cell
essentially never divides — which a memoryless model cannot express. Each driver mutation
draws its own fitness increment, so every cell carries an individual fitness rather than
inheriting a shared subclone label.

This is a **simulator**: it does not infer trees, fit parameters, or compare populations.

## Install

Requires Julia 1.10 or newer.

```julia
using Pkg
Pkg.add(url = "https://github.com/alexsteininfo/MutationLoadDynamics.jl")
```

## Quickstart

```julia
using MutationLoadDynamics, Distributions, Random, Statistics

rng = MersenneTwister(42)
pop = initialize_population(fitness_init = 1.0)

block = NonMarkovBlock(
    birth_dist     = f -> Gamma(5.0, 1 / (5 * f)),   # mean division time 1/f, CV = 1/√5
    death_dist     = f -> Gamma(5.0, 1 / (5 * 0.3)), # mean death time 1/0.3
    stopfunction   = pop -> popsize(pop) >= 10_000,
    driver_dist    = Exponential(0.05),              # effect size of one driver
    fitness_update = (f, δ) -> f + δ,                # additive selection
    ν              = 0.2,                            # mean drivers per daughter
    restart_on_extinction = true,
)

spec = MeasurementSpec(
    trajectory_dt     = 0.5,
    snapshot_triggers = [AtEnd()],
    snapshot_stats    = [SFS(), FitnessDistribution(), DriversPerCell()],
)
acc = MeasurementAccumulator(spec)

simulate!(pop, block, rng; accumulator = acc)
m = finalize_measurements(acc)

root = getsingleroot(allcells(pop))
mutations_per_cell(root)      # burden of every living cell
leaf_depths(root)             # divisions from founder to each living cell
branch_spectrum(root)         # topology, separated from the mutation rate

# Sequence 1 000 of the 10 000 cells, the way an experiment would.
out = sample_trees(pop, SamplingSpec(1_000); seed = UInt64(0xBEEF))
sitefrequencyspectrum(out.samples[1].root, 1_000)
```

Worked examples are in [`examples/`](examples): single-cell expansion, growth with
drivers, chained two-phase runs, and an arbitrary initial condition.

## What it models

**Competing risks.** At birth a cell draws both a division time and a death time; the
earlier one happens. This is exact for *any* pair of distributions — no
acceptance–rejection, no rate bound, no time discretisation — because a cell's fitness
cannot change between its birth and its own division.

**A global min-heap** orders every pending event by absolute time, so the next event
processed is always the globally earliest. This is the Next Reaction Method for a
non-Markovian process, and it costs `O(N log N)` to grow to `N` cells. The exponential is
not a separate code path; it is the special case in which this and Gillespie agree.

**Waiting-time modes** are whatever you write: deterministic (`Dirac`), exponential
(Markov), Gamma, Weibull, log-normal. Parameterised to hold the mean at `1/(bf)`, the
shape becomes a pure noise parameter — and it matters, because at fixed mean division
time the Malthusian rate falls from `b` at `k = 1` to `b·ln 2` as `k → ∞`.

**Selection modes** are two lines each — pick `driver_dist`, pick `fitness_update`:
neutral, additive (fixed or random effect), multiplicative, winner-takes-all
(`max(f, 1+δ)`), capped/logistic with diminishing returns, and deleterious load.

**The tree is the output.** Cells that divide become internal nodes; dead cells are
pruned along with any ancestor they leave childless, so what remains is the reduced tree
of the survivors. Site-frequency and branch spectra, per-cell burdens, divisional depths,
pairwise distances and coalescence times all come from it — on the full population or on
a uniform sample of it.

## Scope

| | |
|:---|:---|
| Here | the birth–death process, per-cell fitness, the lineage tree, tree statistics, leaf sampling |
| Not here | neutral passengers (derive them from `leaf_depths` and `branch_spectrum`), spatial structure, clone-level bookkeeping |
| [`CopyNumberEvolution.jl`](https://github.com/alexsteininfo/CopyNumberEvolution.jl) | copy-number alterations along a tree from this package |
| [`BirthDeathMutation`](https://github.com/alexsteininfo/BirthDeathMutation) | the low-`ν` regime, where clones are the natural unit |
| [`gITH-nonMarkovian`](https://github.com/alexsteininfo/gITH-nonMarkovian) | the analyses, figures and theory built on this simulator |

## Documentation

Build the manual locally:

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(path = "."); Pkg.instantiate()'
julia --project=docs docs/make.jl
# then open docs/build/index.html
```

It covers the concepts and the algorithm, every field of `NonMarkovBlock` with the
waiting-time and growth modes spelled out, the selection modes, measurement and output,
the tree statistics and which of them are co-indexed, leaf sampling and its
reproducibility guarantees, and the open questions.

## Open questions

Documented rather than silently settled — see the manual's Limitations page:

- **Rates are frozen at each cell's birth.** Density-dependent rules evaluate the density
  as of a cell's birth, so a homeostatic ceiling is overshot. Chaining blocks is the exact
  alternative.
- **Cell-cycle durations are independent between relatives**, which underestimates the
  lineage-level clustering seen in real data.
- **`leaf_depths` is not co-indexed** with the other per-leaf statistics, and its order is
  frozen because stored results depend on it.
- **Coalescence times across independent founders** return the time to the seeding moment,
  a convention rather than a settled answer.

## Changes

Breaking changes between releases are listed in [CHANGELOG.md](CHANGELOG.md).

## Tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## License

See [LICENSE](LICENSE).
