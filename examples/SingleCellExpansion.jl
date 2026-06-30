##
## Example 1: Single-cell expansion with driver mutation accumulation
##
## A single founding cell grows to K = 5 000 cells. Division times follow
## Gamma(2, 1/f), so fitter cells divide faster. Each daughter acquires
## mutations at Poisson rate ν = 0.5, each incrementing fitness by
## δ ~ Exponential(0.05). restart_on_extinction retries automatically if
## the founding lineage goes extinct early.
##

using Pkg
Pkg.activate(dirname(@__DIR__))

using MutationLoadDynamics
using Random
using Distributions
using Statistics: mean, std

rng = MersenneTwister(12)

K = 5_000

block = NonMarkovBlock(
    birth_dist     = f -> Gamma(2.0, 1.0 / f),   # fitter cells divide faster
    death_dist     = f -> Gamma(2.0, 10.0),       # mean death time = 20 (rare relative to birth)
    stopfunction   = pop -> popsize(pop) >= K,
    driver_dist    = Exponential(0.05),
    fitness_update = (f, δ) -> f + δ,
    ν              = 0.5,
    restart_on_extinction = true,
)

pop = initialize_population(fitness_init = 1.0)
simulate!(pop, block, rng)

println("=== Single-Cell Expansion ===")
println("Population size : ", popsize(pop))
println("Simulation time : ", round(pop.t, digits = 2))

fits = fitness_per_cell(pop)
println("\nFitness across cells:")
println("  Mean   : ", round(mean(fits), digits = 4))
println("  Std    : ", round(std(fits),  digits = 4))
println("  Min    : ", round(minimum(fits), digits = 4))
println("  Max    : ", round(maximum(fits), digits = 4))

ks = mutations_per_cell(pop)
println("\nDriver mutations per cell:")
println("  Mean   : ", round(mean(Float64.(ks)), digits = 2))
println("  Min    : ", minimum(ks))
println("  Max    : ", maximum(ks))
println("  Clonal : ", clonal_mutations(pop), " (shared by all cells)")

sfs = sitefrequencyspectrum(pop)
sfs_nonzero = [(k, sfs[k]) for k in 1:length(sfs) if sfs[k] > 0]
println("\nDriver SFS (frequency, count) — top 10:")
for (freq, cnt) in first(sfs_nonzero, 10)
    println("  $freq cells: $cnt driver mutation(s)")
end
