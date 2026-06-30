function make_block(; Nmax=50, ν=0.0, s_mean=0.1, restart=false)
    NonMarkovBlock(
        birth_dist     = f -> Gamma(2.0, 1.0 / f),
        death_dist     = f -> Gamma(2.0, 10.0),       # mean death time >> birth time
        stopfunction   = pop -> popsize(pop) >= Nmax,
        driver_dist    = Exponential(s_mean),
        fitness_update = (f, δ) -> f + δ,
        ν              = ν,
        restart_on_extinction = restart,
    )
end

@testset "population grows to target size" begin
    rng = MersenneTwister(1)
    pop = initialize_population(fitness_init = 1.0)
    simulate!(pop, make_block(Nmax = 20), rng)
    @test popsize(pop) >= 20
end

@testset "ν=0 keeps mean fitness constant" begin
    rng = MersenneTwister(7)
    pop = initialize_population(fitness_init = 1.0)
    simulate!(pop, make_block(Nmax = 100, ν = 0.0), rng)
    fits = fitness_per_cell(pop)
    @test all(f ≈ 1.0 for f in fits)
    @test mean(fits) ≈ 1.0
end

@testset "ν>0 increases mean fitness with positive selection" begin
    rng = MersenneTwister(3)
    pop = initialize_population(fitness_init = 1.0)
    simulate!(pop, make_block(Nmax = 200, ν = 1.0, s_mean = 0.1), rng)
    @test mean(fitness_per_cell(pop)) > 1.0
end

@testset "mutations_per_cell > 0 when ν>0" begin
    rng = MersenneTwister(5)
    pop = initialize_population(fitness_init = 1.0)
    simulate!(pop, make_block(Nmax = 50, ν = 2.0), rng)
    @test mean(mutations_per_cell(pop)) > 0
end

@testset "birth_dist is sampled correctly (single-cell draws)" begin
    # Verify schedule_cell! samples division times from the correct distribution
    # by observing the spacing between birth and first division in isolated cells.
    # Uses many independent single-cell simulations to avoid the population-level
    # length-biased sampling artefact that arises in a full growing tree.
    rng   = MersenneTwister(11)
    shape = 3.0
    scale = 0.5   # mean = shape * scale = 1.5
    n_obs = 400
    observed = Float64[]
    for _ in 1:n_obs
        pop = initialize_population(fitness_init = 1.0)
        block = NonMarkovBlock(
            birth_dist     = f -> Gamma(shape, scale),
            death_dist     = f -> Gamma(2.0, 1000.0),   # negligible death
            stopfunction   = pop -> popsize(pop) >= 3,   # stop after first division
            driver_dist    = Exponential(0.01),
            fitness_update = (f, δ) -> f + δ,
            ν              = 0.0,
        )
        simulate!(pop, block, rng)
        root = getsingleroot(allcells(pop))
        lt   = celllifetimes(root; excludeliving = true)
        isempty(lt) || push!(observed, lt[1])
    end
    @test length(observed) > 300
    expected_mean = shape * scale
    @test abs(mean(observed) - expected_mean) / expected_mean < 0.15
end

@testset "popsize matches allcells length" begin
    rng = MersenneTwister(99)
    pop = initialize_population(fitness_init = 1.0)
    simulate!(pop, make_block(Nmax = 80), rng)
    @test popsize(pop) == length(allcells(pop))
end

@testset "restart_on_extinction reaches target" begin
    rng = MersenneTwister(77)
    block = NonMarkovBlock(
        birth_dist     = f -> Exponential(1.0 / f),
        death_dist     = f -> Exponential(2.0 / f),   # death > birth, high extinction
        stopfunction   = pop -> popsize(pop) >= 10,
        driver_dist    = Exponential(0.1),
        fitness_update = (f, δ) -> f + δ,
        ν              = 0.0,
        restart_on_extinction = true,
    )
    pop = initialize_population(fitness_init = 1.0)
    simulate!(pop, block, rng)
    @test popsize(pop) >= 10
end

@testset "trajectory recording" begin
    rng = MersenneTwister(42)
    spec = MeasurementSpec(
        trajectory_dt     = 0.5,
        snapshot_triggers = [AtEnd()],
        snapshot_stats    = [FitnessDistribution()],
    )
    acc = MeasurementAccumulator(spec)
    pop = initialize_population(fitness_init = 1.0)
    simulate!(pop, make_block(Nmax = 100, ν = 1.0), rng; accumulator = acc)
    m = finalize_measurements(acc)
    @test length(m.trajectory) > 0
    @test length(m.snapshots) == 1
    @test m.snapshots[1].trigger isa AtEnd
    @test !isnothing(m.snapshots[1].fitness_distribution)
    for tp in m.trajectory
        @test tp.N_total > 0
        @test isfinite(tp.mean_fitness)
    end
end
