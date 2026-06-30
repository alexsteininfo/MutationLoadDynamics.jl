function simple_pop(; Nmax=50, ν=1.0, rng=MersenneTwister(1))
    block = NonMarkovBlock(
        birth_dist     = f -> Gamma(2.0, 1.0 / f),
        death_dist     = f -> Gamma(2.0, 20.0),
        stopfunction   = pop -> popsize(pop) >= Nmax,
        driver_dist    = Exponential(0.1),
        fitness_update = (f, δ) -> f + δ,
        ν              = ν,
    )
    pop = initialize_population(fitness_init = 1.0)
    simulate!(pop, block, rng)
    return pop
end

@testset "fitness_per_cell length and positivity" begin
    pop = simple_pop()
    fits = fitness_per_cell(pop)
    @test length(fits) == popsize(pop)
    @test all(f > 0 for f in fits)
end

@testset "mutations_per_cell is non-negative" begin
    pop = simple_pop(ν = 2.0)
    ks = mutations_per_cell(pop)
    @test length(ks) == popsize(pop)
    @test all(k >= 0 for k in ks)
    @test mean(Float64.(ks)) > 0
end

@testset "ν=0 gives zero mutations per cell" begin
    pop = simple_pop(ν = 0.0)
    @test all(k == 0 for k in mutations_per_cell(pop))
end

@testset "SFS sums correctly" begin
    pop = simple_pop(ν = 2.0)
    sfs = sitefrequencyspectrum(pop)
    @test length(sfs) == popsize(pop)
    @test all(s >= 0 for s in sfs)
    # each driver mutation appears in 1..N cells; total weighted count == total mutations
    N = popsize(pop)
    total_from_sfs = sum(sfs[k] * k for k in 1:N)
    total_from_cells = sum(mutations_per_cell(pop))
    @test total_from_sfs == total_from_cells
end

@testset "pairwisedistances are non-negative" begin
    pop = simple_pop(ν = 1.0, Nmax = 20)
    dists = pairwisedistances(pop)
    @test all(d >= 0 for d in dists)
    @test length(dists) == binomial(popsize(pop), 2)
end

@testset "pairwisedistance is symmetric" begin
    pop = simple_pop(ν = 1.0, Nmax = 15)
    cells = allcells(pop)
    i, j = cells[1], cells[2]
    @test pairwisedistance(i, j) == pairwisedistance(j, i)
end

@testset "findMRCA returns a node" begin
    pop = simple_pop(Nmax = 30)
    mrca = findMRCA(pop)
    @test !isnothing(mrca)
end

@testset "coalescence_times are positive" begin
    pop = simple_pop(Nmax = 20)
    ct = coalescence_times(pop)
    @test all(t >= 0 for t in ct)
    @test length(ct) == binomial(popsize(pop), 2)
end

@testset "mean_k and var_k" begin
    pop = simple_pop(ν = 2.0, Nmax = 100)
    mk = mean_k(pop)
    vk = var_k(pop)
    @test mk >= 0.0
    @test vk >= 0.0
end
