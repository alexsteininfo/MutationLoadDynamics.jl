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

# ── on_division / on_restart hooks ────────────────────────────────────────────

@testset "block without hooks still constructs and runs" begin
    block = make_block(Nmax = 20)
    @test isnothing(block.on_division)
    @test isnothing(block.on_restart)
    rng = MersenneTwister(1234)
    pop = initialize_population(fitness_init = 1.0)
    simulate!(pop, block, rng)
    @test popsize(pop) >= 20
end

@testset "on_division fires once per birth with the right arguments" begin
    rng     = MersenneTwister(2024)
    records = Tuple{Int, Int, Int, Int}[]
    parents_ok = Ref(true)
    block = NonMarkovBlock(
        birth_dist     = f -> Gamma(2.0, 1.0 / f),
        death_dist     = f -> Gamma(2.0, 1.0e6),   # effectively d = 0
        stopfunction   = pop -> popsize(pop) >= 40,
        driver_dist    = Exponential(0.05),
        fitness_update = (f, δ) -> f + δ,
        ν              = 0.0,
        on_division    = function (pop, parent, d1, d2)
            (d1.parent === parent && d2.parent === parent) || (parents_ok[] = false)
            push!(records, (popsize(pop), parent.data.id, d1.data.id, d2.data.id))
            return nothing
        end,
    )
    pop = initialize_population(fitness_init = 1.0)
    simulate!(pop, block, rng)

    @test parents_ok[]
    @test !isempty(records)
    # With no deaths every event is a birth, so the popsize seen inside the hook
    # walks 2, 3, 4, ... — i.e. pre-division size + 1, one per division.
    @test [r[1] for r in records] == collect(2:(1 + length(records)))
    @test popsize(pop) == last(records)[1]
    # Daughters are allocated consecutive fresh ids.
    @test all(r[4] == r[3] + 1 for r in records)
end

@testset "on_division fitness applies to the daughter's own first division" begin
    # Deterministic timing: a cell of fitness f divides exactly 1/f after its birth.
    # The root divides at t = 1. The hook doubles d1's fitness, so d1 must divide at
    # t = 1.5, not t = 2.0. A hook placed after schedule_cell! fails only this test.
    rng    = MersenneTwister(4242)
    boosted = Ref{Any}(nothing)
    block = NonMarkovBlock(
        birth_dist     = f -> Dirac(1.0 / f),
        death_dist     = f -> Dirac(1.0e6),
        stopfunction   = pop -> popsize(pop) >= 3,
        driver_dist    = Dirac(0.0),
        fitness_update = (f, δ) -> f,
        ν              = 0.0,
        on_division    = function (pop, parent, d1, d2)
            isnothing(boosted[]) || return nothing
            c = d1.data
            d1.data   = NonMarkovCell(c.id, c.birthtime, c.mutations, 2.0)
            boosted[] = d1
            return nothing
        end,
    )
    pop = initialize_population(fitness_init = 1.0)
    simulate!(pop, block, rng)

    d1 = boosted[]
    @test d1.data.fitness == 2.0
    @test d1.data.birthtime ≈ 1.0
    @test celllifetime(d1) ≈ 0.5
    @test pop.t ≈ 1.5
end

@testset "a single injected driver forms exactly one clade" begin
    rng      = MersenneTwister(31337)
    N_critic = 8
    injected = Ref(false)
    root_of_clade = Ref{Any}(nothing)
    block = NonMarkovBlock(
        birth_dist     = f -> Gamma(5.0, 1.0 / (5.0 * f)),
        death_dist     = f -> Gamma(5.0, 1.0e5),   # negligible death
        stopfunction   = pop -> popsize(pop) >= 120,
        driver_dist    = Dirac(0.0),
        fitness_update = (f, δ) -> f,
        ν              = 0.0,
        on_division    = function (pop, parent, d1, d2)
            injected[] && return nothing
            popsize(pop) == N_critic + 1 || return nothing
            injected[] = true
            c = d1.data
            d1.data = NonMarkovCell(c.id, c.birthtime, c.mutations, 1.5)
            root_of_clade[] = d1
            return nothing
        end,
    )
    pop = initialize_population(fitness_init = 1.0)
    simulate!(pop, block, rng)

    @test injected[]
    fit_ids   = Set(c.data.id for c in allcells(pop) if c.data.fitness > 1.0)
    clade_ids = Set(l.data.id for l in Leaves(root_of_clade[]) if haskey(pop.cells, l.data.id))
    @test !isempty(fit_ids)
    @test fit_ids == clade_ids
    @test all(c.data.fitness == 1.0 for c in allcells(pop) if !(c.data.id in fit_ids))
end

@testset "driver clone size increases with selection strength" begin
    N_critic = 8
    function _clone_fraction(s, seed)
        rng      = MersenneTwister(seed)
        injected = Ref(false)
        clade    = Ref{Any}(nothing)
        block = NonMarkovBlock(
            birth_dist     = f -> Gamma(5.0, 1.0 / (5.0 * f)),
            death_dist     = f -> Gamma(5.0, 1.0e5),
            stopfunction   = pop -> popsize(pop) >= 300,
            driver_dist    = Dirac(0.0),
            fitness_update = (f, δ) -> f,
            ν              = 0.0,
            on_division    = function (pop, parent, d1, d2)
                injected[] && return nothing
                popsize(pop) == N_critic + 1 || return nothing
                injected[] = true
                c = d1.data
                d1.data  = NonMarkovCell(c.id, c.birthtime, c.mutations, 1.0 + s)
                clade[]  = d1
                return nothing
            end,
        )
        pop = initialize_population(fitness_init = 1.0)
        simulate!(pop, block, rng)
        injected[] || return NaN
        n = count(l -> haskey(pop.cells, l.data.id), Leaves(clade[]))
        return n / popsize(pop)
    end
    seeds = 101:118
    weak   = filter(!isnan, [_clone_fraction(0.1, s) for s in seeds])
    strong = filter(!isnan, [_clone_fraction(1.5, s) for s in seeds])
    @test length(weak) > 10 && length(strong) > 10
    @test mean(strong) > mean(weak)
end

@testset "on_restart lets a stateful hook survive an extinction restart" begin
    # A closure holding an `injected` flag would otherwise burn its one injection on an
    # attempt that later goes extinct, and never inject again. `on_restart` resets it.
    # Seed 4 is chosen because it does exactly that: it injects, dies, and restarts.
    #
    # That narrative is a property of seed 4 under Julia's post-1.11 `rand`/randperm-style
    # range-sampling algorithm (the CI floor is 1.9, which predates that change and does not
    # reproduce the same draws from this seed), so the assertions below are skipped on
    # Julia < 1.11 rather than asserted against a stream that won't exhibit the same run.
    function _restart_run(reset_flag)
        rng        = MersenneTwister(4)
        restarts   = Ref(0)
        injected   = Ref(false)
        injections = Ref(0)
        block = NonMarkovBlock(
            birth_dist     = f -> Exponential(1.0 / f),
            death_dist     = f -> Exponential(2.0 / f),   # death > birth, high extinction
            stopfunction   = pop -> popsize(pop) >= 10,
            driver_dist    = Dirac(0.0),
            fitness_update = (f, δ) -> f,
            ν              = 0.0,
            restart_on_extinction = true,
            on_division    = function (pop, parent, d1, d2)
                injected[] && return nothing
                popsize(pop) == 3 || return nothing
                injected[]    = true
                injections[] += 1
                c = d1.data
                d1.data = NonMarkovCell(c.id, c.birthtime, c.mutations, 1.2)
                return nothing
            end,
            on_restart     = function (pop)
                restarts[] += 1
                reset_flag && (injected[] = false)
                return nothing
            end,
        )
        pop = initialize_population(fitness_init = 1.0)
        simulate!(pop, block, rng)
        return (popsize = popsize(pop), restarts = restarts[], injections = injections[],
                boosted = any(c.data.fitness > 1.0 for c in allcells(pop)))
    end

    if VERSION >= v"1.11"
        with_reset = _restart_run(true)
        @test with_reset.popsize >= 10
        @test with_reset.restarts >= 1        # this seed does go extinct
        @test with_reset.injections >= 2      # it re-injected after a restart
        @test with_reset.boosted              # the surviving population carries the driver

        # Same seed, same everything, except the hook does not reset its own flag: the one
        # injection is spent on a doomed attempt and the surviving population has no driver.
        # This is the failure mode `on_restart` exists to prevent.
        without_reset = _restart_run(false)
        @test without_reset.restarts == with_reset.restarts
        @test without_reset.injections == 1
        @test !without_reset.boosted
    else
        with_reset = _restart_run(true)
        @test_skip with_reset.restarts >= 1
        @test_skip with_reset.injections >= 2
        @test_skip with_reset.boosted
    end
end
