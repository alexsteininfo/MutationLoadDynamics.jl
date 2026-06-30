using MutationLoadDynamics
using Distributions
using Random
using Statistics

rng = MersenneTwister(42)

pop = initialize_population(fitness_init = 1.0)

block = NonMarkovBlock(
    # Division time ~ Gamma(2, 1/f): fitter cells divide faster on average
    birth_dist   = f -> Gamma(2.0, 1.0 / f),
    # Death time ~ Gamma(2, 5): mean death time = 10, fitness-independent
    death_dist   = f -> Gamma(2.0, 5.0),
    stopfunction = pop -> popsize(pop) >= 5_000,
    driver_dist  = Exponential(0.05),
    fitness_update = (f, δ) -> f + δ,   # additive fitness
    ν = 0.5,                             # mean 0.5 drivers per daughter per division
    restart_on_extinction = true,
)

spec = MeasurementSpec(
    trajectory_dt     = 1.0,
    snapshot_triggers = [AtEnd()],
    snapshot_stats    = [SFS(), FitnessDistribution(), DriversPerCell()],
)
acc = MeasurementAccumulator(spec)

simulate!(pop, block, rng; accumulator = acc)
m = finalize_measurements(acc)

println("Final population size : ", popsize(pop))
println("Simulation time       : ", round(pop.t, digits = 2))
println("Mean fitness          : ", round(mean(fitness_per_cell(pop)), digits = 4))
println("Std fitness           : ", round(std(fitness_per_cell(pop)), digits = 4))
println("Mean drivers/cell     : ", round(mean_k(pop), digits = 2))
println("Trajectory points     : ", length(m.trajectory))

if !isempty(m.snapshots)
    snap = m.snapshots[1]
    sfs  = snap.sfs
    println("SFS length (nonzero)  : ", count(>(0), sfs))
end
