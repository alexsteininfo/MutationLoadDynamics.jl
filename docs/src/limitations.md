# Limitations and open questions

## Deliberately out of scope

**Neutral passenger mutations.** They are not simulated because they need no state. For a
neutral per-division rate ``m``, the per-cell burden is `Poisson(m * depth)` drawn against
[`leaf_depths`](@ref) and the expected spectrum is ``m`` times the
[`branch_spectrum`](@ref). Deriving them afterwards is exact in distribution, costs
nothing, and lets one simulated tree serve every value of ``m`` — whereas simulating them
would multiply memory by the mutation count for no gain. See
[Neutral evolution](selection.md#Neutral-evolution).

**Extinct lineages.** A cell that dies is pruned from the tree, together with any ancestor
it leaves childless. Under a supercritical process most cells ever born are in lineages
that die out; keeping them would dominate memory while contributing to no observable of
the surviving population. The consequence is that the package cannot answer questions
about lineages that failed — how many clones were lost to drift, say — unless you record
them yourself in an `on_division` hook as they happen.

**Spatial structure.** Cells have no position and no neighbours. Every cell competes with
every other only through the global stop condition. Spatially structured growth changes
tree shape substantially, and belongs in a different simulator rather than behind a flag
here.

**Clone-level bookkeeping.** There are no subclone labels, because in the regime this
package targets (``ν \gtrsim 0.1``) every cell has a distinct fitness history and a clone
label would be meaningless. For ``ν \ll 1``, where clones *are* the natural unit, see
[`BirthDeathMutation`](https://github.com/alexsteininfo/BirthDeathMutation).

**Copy-number alterations, sequencing noise, variant calling.** Everything between a
lineage tree and a data matrix lives downstream.
[`CopyNumberEvolution.jl`](https://github.com/alexsteininfo/CopyNumberEvolution.jl)
consumes a tree from this package and draws copy-number alterations along its edges.

**Inference of any kind.** No tree reconstruction, no parameter estimation, no distances
between populations. Keeping them apart is what lets estimators run on real patient data
without a simulator in their dependency chain.

## Known constraints of the model

### Rates are frozen at each cell's birth

A cell's waiting times are drawn once, when it is born, and the simulator has no mechanism
to invalidate a pending event when something external changes. Fitness cannot change
mid-cycle (mutations happen only at division, to daughters), which is exactly what makes
the competing-risks scheduling **exact** for the model as specified.

The cost is that anything genuinely time-varying is approximate. Density-dependent birth
rates evaluate the density as of each cell's birth, so a homeostatic ceiling is
overshot — see
[Density-dependent and homeostatic growth](blocks.md#Density-dependent-and-homeostatic-growth).
The same applies to a treatment that starts mid-run, or to any interaction between
cells. The sound alternative is to change the *regime* by chaining blocks, where each
block's rates are constant and the boundary is handled exactly.

Whether the package should grow an explicit event-invalidation mechanism — pop the
affected cells' events and redraw from the conditional law ``T \mid T > \text{age}`` —
is genuinely undecided. It would make true density dependence exact, at the cost of a
heap that supports deletion and of drawing from conditional distributions that most
`Distributions.jl` types do not expose directly.

### Chaining and `restart_on_extinction` do not combine

The extinction snapshot is taken at the start of a `simulate!` call, so restarting a
chained second block restores cells born under the first block and rewinds them to
birthtimes predating the boundary. Rescheduling them then has the very age-conditioning
defect that carrying the event queue exists to avoid. The snapshot is also a `deepcopy` of
a live tree, so on a large chained population it duplicates the entire history.

Use `restart_on_extinction` on the first block only, where the population is small and its
birthtimes are the initial ones; drive anything more selective with your own retry loop.

### Cell-cycle durations are independent between relatives

Every cell draws its waiting time independently. Real lineages do not behave that way:
cousin–cousin correlations in cell-cycle duration dominate over mother–daughter ones
(Sandler et al. 2015), so an i.i.d. model underestimates lineage-level clustering. Adding
inheritance of cycle length would mean carrying a per-cell latent variable through
`_make_daughter` and letting `birth_dist` see it. That is a small change to the cell type
and a large change to what the package claims, and it has not been made.

### The Gamma is a convenience, not a claim

The best-fitting empirical family for mammalian interdivision times is the exponentially
modified Gamma (Golubev 2016), not the plain Gamma. The plain Gamma is the recommended
default here because it captures the essential non-Markovian character with one extra
parameter. Nothing prevents you from supplying something better — `birth_dist` returns any
distribution — and if your conclusions turn on the shape of the cycle-time distribution,
you should.

## Open questions

These are genuinely undecided. The current behaviour is a documented placeholder, not a
settled position.

### Coalescence times across independent founders

For two cells descending from different founders of an
`initialize_population(N)` forest, [`coalescence_times`](@ref) returns the time back to
the seeding moment — treating independent founders as coalescing when the population was
initialised. That is a lower bound on the true divergence and it does not error, but it
puts an artefactual spike at ``t - t_0`` into any pooled histogram. Whether the right
answer is `Inf`, `missing`, or an error is unresolved; today it is a convention, and the
safe course is to sample each founder's tree separately.

### `leaf_depths` is not co-indexed, and cannot be fixed

It returns depths in its own stack order, which is neither `Leaves` order nor `Dict`
order, so it is valid only as a pooled distribution. The order is frozen because a large
volume of stored results was produced with it. The clean resolution is a second,
explicitly co-indexed function rather than a change to this one — it has not been added
because nothing has needed it yet.

### `pop.t` versus `age(root)`

`age(::Population)` is `pop.t`, the time of the last event of any kind. `age(::BinaryNode)`
is the birthtime of the last-born leaf, i.e. the last *division*. A run that ends on a
death has `age(root) < pop.t`, and the two `coalescence_times` methods inherit that
difference in their default reference time. Having two functions named `age` mean two
different things is a wart; changing either would break stored analyses, so both are
documented rather than reconciled.

### The realised cell-cycle distribution is doubly biased

[`celllifetimes`](@ref) is shifted earlier than the `birth_dist` you specified, for two
independent reasons. Competing risks mean a cell that dies never contributes a completed
lifetime, so what remains is the division time *conditioned on dividing first*. And
stopping at a fixed population size means the cells that have already divided are the ones
with short cycles, while the long ones are still in progress — a bias that survives even
with no death at all (a pure-birth neutral run with a true mean cycle of 1.0 returns a mean
completed lifetime near 0.87).

Both are correct behaviour, and together they are what an experiment measuring
interdivision times in a growing culture would actually see. But they mean
`celllifetimes` cannot serve as a check that `birth_dist` was implemented as intended, and
`excludeliving = false` does not fix it — those cells are right-censored at `age(root)`
rather than completed. Whether the package should offer an explicitly de-biased estimator,
or expose the censoring so a survival-analysis estimate is possible, is undecided.

### No cap on fitness

`fitness_update` is unconstrained, so an unbounded rule (additive or multiplicative with
no ceiling) lets fitness grow without limit and division times shrink toward zero. Nothing
in the package objects; the run simply becomes dominated by one lineage and, eventually,
numerically degenerate. Whether a soft guard belongs in the simulator — and what it should
do — is unsettled. The
[capped selection mode](selection.md#Capped-selection-(diminishing-returns)) is the
current answer, and it is opt-in.
