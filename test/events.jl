@testset "CellEvent ordering" begin
    pop  = initialize_population(fitness_init = 1.0)
    node = first(values(pop.cells))
    e1 = MutationLoadDynamics.CellEvent(0.5, node, :birth)
    e2 = MutationLoadDynamics.CellEvent(1.5, node, :death)
    @test e1 < e2
    @test !(e2 < e1)
end

@testset "schedule_cell! event times are after birthtime" begin
    rng  = MersenneTwister(7)
    pop  = initialize_population(fitness_init = 1.0)
    node = first(values(pop.cells))
    block = NonMarkovBlock(
        birth_dist     = f -> Gamma(2.0, 1.0 / f),
        death_dist     = f -> Gamma(2.0, 5.0),
        stopfunction   = pop -> false,
        driver_dist    = Exponential(0.1),
        fitness_update = (f, δ) -> f + δ,
        ν              = 0.0,
    )
    heap = DataStructures.BinaryMinHeap{MutationLoadDynamics.CellEvent}()
    for _ in 1:100
        MutationLoadDynamics.schedule_cell!(heap, node, block, rng)
    end
    while !isempty(heap)
        e = pop!(heap)
        @test e.time >= node.data.birthtime
        @test e.event_type in (:birth, :death)
    end
end

@testset "low death rate → births dominate events" begin
    rng  = MersenneTwister(3)
    pop  = initialize_population(fitness_init = 1.0)
    node = first(values(pop.cells))
    block = NonMarkovBlock(
        birth_dist     = f -> Exponential(1.0 / f),  # mean division time ≈ 1
        death_dist     = f -> Exponential(100.0),    # mean death time = 100 → rare
        stopfunction   = pop -> false,
        driver_dist    = Exponential(0.1),
        fitness_update = (f, δ) -> f + δ,
        ν              = 0.0,
    )
    n_birth = 0
    n_death = 0
    for _ in 1:500
        heap = DataStructures.BinaryMinHeap{MutationLoadDynamics.CellEvent}()
        MutationLoadDynamics.schedule_cell!(heap, node, block, rng)
        e = pop!(heap)
        e.event_type == :birth ? (n_birth += 1) : (n_death += 1)
    end
    @test n_birth > 450   # roughly 99% should be births
end
