##
## Example 3: Effect of initial condition on genetic diversity
##
## Compares two scenarios that both grow to the same final population size (N = 1 000)
## but differ in the number of founding cells:
##
##   Scenario A: 1 founding cell  — each mutation in the final population has a
##               well-defined lineage depth; many singletons, few clonal mutations.
##
##   Scenario B: 100 founding cells — independent lineages; driver mutations acquired
##               by one founding cell are never shared with cells from other lineages,
##               so the SFS has more singletons and the mean pairwise distance is larger.
##
## Both scenarios use identical block parameters (same rng seed) so differences are
## solely due to the initial condition.
##

using Pkg
Pkg.activate(dirname(@__DIR__))

using MutationLoadDynamics
using Random
using Distributions
using Statistics: mean, std

rng_a = MersenneTwister(77)
rng_b = MersenneTwister(77)   # same seed — only the initial condition differs

N_final = 1_000

shared_block = NonMarkovBlock(
    birth_dist     = f -> Gamma(2.0, 1.0 / f),
    death_dist     = f -> Gamma(2.0, 10.0),
    stopfunction   = pop -> popsize(pop) >= N_final,
    driver_dist    = Exponential(0.05),
    fitness_update = (f, δ) -> f + δ,
    ν              = 0.5,
    restart_on_extinction = true,
)

# ── Scenario A: single founding cell ────────────────────────────────────────

popA = initialize_population(fitness_init = 1.0)
simulate!(popA, shared_block, rng_a)

sfsA = sitefrequencyspectrum(popA)
mpcA = mutations_per_cell(popA)

println("=== Scenario A: 1 founding cell → $N_final cells ===")
println("Simulation time    : $(round(popA.t, digits = 2))")
println("Mean fitness       : $(round(mean(fitness_per_cell(popA)), digits = 4))")
println("Mean drivers/cell  : $(round(mean(Float64.(mpcA)), digits = 2))")
println("Clonal mutations   : $(clonal_mutations(popA))  (shared by all $N_final cells)")
singleton_frac_A = sfsA[1] / sum(sfsA)
println("Singleton fraction : $(round(100 * singleton_frac_A, digits = 1))%")
idx_a = 1:min(30, popsize(popA))
pd_a  = pairwisedistances(popA, collect(idx_a))
println("Mean pairwise dist : $(round(mean(pd_a), digits = 2))  (sample of $(length(idx_a)) cells)")

# ── Scenario B: 100 founding cells ──────────────────────────────────────────

popB = initialize_population(100; fitness_init = 1.0)
simulate!(popB, shared_block, rng_b)

sfsB = sitefrequencyspectrum(popB)
mpcB = mutations_per_cell(popB)

println("\n=== Scenario B: 100 founding cells → $N_final cells ===")
println("Simulation time    : $(round(popB.t, digits = 2))")
println("Mean fitness       : $(round(mean(fitness_per_cell(popB)), digits = 4))")
println("Mean drivers/cell  : $(round(mean(Float64.(mpcB)), digits = 2))")
println("Clonal mutations   : $(clonal_mutations(popB))  (shared by all $N_final cells)")
singleton_frac_B = sfsB[1] / sum(sfsB)
println("Singleton fraction : $(round(100 * singleton_frac_B, digits = 1))%")
idx_b = 1:min(30, popsize(popB))
pd_b  = pairwisedistances(popB, collect(idx_b))
println("Mean pairwise dist : $(round(mean(pd_b), digits = 2))  (sample of $(length(idx_b)) cells)")

println()
println("Expected: B has more singletons (100% private mutations) because each sub-lineage is shallow")
println("          (~10 cells per founder), so most mutations land on terminal branches.")
println("          B also has lower mean pairwise distance because cells carry fewer total mutations")
println("          (shallower trees) even though they come from independent lineages.")
