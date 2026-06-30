##
## Example 2: Two-phase simulation — growth then neutral drift
##
## Phase 1: fast growth with high driver rate (ν = 1.0) accumulates fitness diversity.
## Phase 2: continue the same tree with ν = 0 and a timed stop to observe how the
##           fitness landscape evolves by selection alone (no new mutations).
##
## This demonstrates chaining two successive simulate! calls on the same Population
## so the lineage tree spans both phases.
##

using Pkg
Pkg.activate(dirname(@__DIR__))

using MutationLoadDynamics
using Random
using Distributions
using Statistics: mean, std

rng = MersenneTwister(99)

# ── Phase 1: growth with driver accumulation ────────────────────────────────

block1 = NonMarkovBlock(
    birth_dist     = f -> Gamma(2.0, 1.0 / f),
    death_dist     = f -> Gamma(2.0, 10.0),
    stopfunction   = pop -> popsize(pop) >= 1_000,
    driver_dist    = Exponential(0.05),
    fitness_update = (f, δ) -> f + δ,
    ν              = 1.0,
    restart_on_extinction = true,
)

spec1 = MeasurementSpec(
    trajectory_dt     = 0.5,
    snapshot_triggers = [AtEnd()],
    snapshot_stats    = [SFS(), FitnessDistribution(), DriversPerCell()],
)
acc1 = MeasurementAccumulator(spec1)

pop = initialize_population(fitness_init = 1.0)
simulate!(pop, block1, rng; accumulator = acc1)
m1  = finalize_measurements(acc1)

println("=== Phase 1: growth with drivers (ν = 1.0) ===")
println("Population size : ", popsize(pop))
println("Simulation time : ", round(pop.t, digits = 2))
println("Mean fitness    : ", round(mean(m1.snapshots[1].fitness_distribution), digits = 4))
println("Std fitness     : ", round(std(m1.snapshots[1].fitness_distribution),  digits = 4))
println("Mean drivers    : ", round(mean(Float64.(m1.snapshots[1].drivers_per_cell)), digits = 2))
println("Trajectory pts  : ", length(m1.trajectory))

# ── Phase 2: second growth phase, no new mutations — selection sweeps ────────
#
# ν = 0: no new driver mutations. Cells with higher fitness (accumulated in
# phase 1) still divide faster, so the population continues to grow and the
# mean fitness increases by selection alone.

block2 = NonMarkovBlock(
    birth_dist     = f -> Gamma(2.0, 1.0 / f),
    death_dist     = f -> Gamma(2.0, 10.0),
    stopfunction   = pop -> popsize(pop) >= 5_000,
    driver_dist    = Exponential(0.05),            # unused (ν = 0)
    fitness_update = (f, δ) -> f + δ,
    ν              = 0.0,                          # no new driver mutations
)

spec2 = MeasurementSpec(
    trajectory_dt     = 0.5,
    snapshot_triggers = [AtEnd()],
    snapshot_stats    = [FitnessDistribution(), DriversPerCell()],
)
acc2 = MeasurementAccumulator(spec2)

simulate!(pop, block2, rng; accumulator = acc2)
m2  = finalize_measurements(acc2)

println("\n=== Phase 2: selection sweep (ν = 0, growth to 5 000) ===")
println("Population size : ", popsize(pop))
println("Simulation time : ", round(pop.t, digits = 2))
println("Mean fitness    : ", round(mean(m2.snapshots[1].fitness_distribution), digits = 4))
println("Std fitness     : ", round(std(m2.snapshots[1].fitness_distribution),  digits = 4))
println("Mean drivers    : ", round(mean(Float64.(m2.snapshots[1].drivers_per_cell)), digits = 2))
println("Trajectory pts  : ", length(m2.trajectory))
println()
println("(mean fitness rises in phase 2 without new mutations — fitter lineages expand)")
println("(mean drivers/cell also rises as fitter cells were those with more phase-1 drivers)")
