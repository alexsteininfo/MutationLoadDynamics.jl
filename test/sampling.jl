# Both helpers must be defined at the TOP of this file, above every @testset.
# A `@testset` body executes eagerly while the file is being evaluated, so a
# testset that calls a function defined further down raises UndefVarError.
#
# Same fixture as test/statistics.jl. Duplicated deliberately: each test file must
# run standalone under `include`, and the tree is five lines.
function sampling_fixture()
    root = BinaryNode(NonMarkovCell(1, 0.0, 5, 1.0))
    L = leftchild!(root,  NonMarkovCell(2, 1.0, 1, 1.0))
    rightchild!(root,     NonMarkovCell(3, 1.2, 7, 1.0))
    leftchild!(L,         NonMarkovCell(4, 2.0, 2, 1.0))
    rightchild!(L,        NonMarkovCell(5, 2.1, 3, 1.0))
    return root
end

function simple_sampling_pop(; Nmax = 60, ν = 2.0, rng = MersenneTwister(11))
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

@testset "LeafSample records the draw" begin
    root = sampling_fixture()
    s = sample_leaves(root, 2; seed = UInt64(7))
    @test s isa LeafSample
    @test s.n == 2
    @test s.N_full == 3
    @test s.seed == UInt64(7)
    @test s.replicate == 1
    @test length(s.sampled_ids) == 2
    @test allunique(s.sampled_ids)
    @test all(id -> id in (3, 4, 5), s.sampled_ids)
end

@testset "GOLDEN: the draw recipe is frozen (HR-1)" begin
    # These ids are the output of MersenneTwister(seed) + randperm(rng, 3)[1:n]
    # over Leaves order [4, 5, 3]. 4.2 GB of sampled trees on disk were drawn with
    # this exact recipe. If this test fails, the recipe changed and the stored data
    # is no longer reproducible — do NOT update the expected values, fix the code.
    # Generated under Julia 1.12.1; treat a change of Julia minor version as a
    # reason to re-verify against real shards rather than to edit this test.
    root = sampling_fixture()
    leaves = [l.data.id for l in Leaves(root)]
    for (seed, n) in ((UInt64(1), 1), (UInt64(1), 2), (UInt64(42), 2), (UInt64(7), 3))
        rng      = MersenneTwister(seed)
        expected = leaves[randperm(rng, 3)[1:n]]
        @test sample_leaves(root, n; seed = seed).sampled_ids == expected
    end
end

@testset "n = N_full reproduces the source tree exactly" begin
    root = sampling_fixture()
    s = sample_leaves(root, 3; seed = UInt64(1))
    @test sort(s.sampled_ids) == [3, 4, 5]
    @test mutations_per_cell(s.root) == mutations_per_cell(root)
    @test leaf_depths(s.root)        == leaf_depths(root)
    @test sitefrequencyspectrum(s.root, 3) == sitefrequencyspectrum(root, 3)
    @test [l.data.id for l in Leaves(s.root)] == [l.data.id for l in Leaves(root)]
end

@testset "sampling is deterministic in the seed" begin
    root = sampling_fixture()
    a = sample_leaves(root, 2; seed = UInt64(42))
    b = sample_leaves(root, 2; seed = UInt64(42))
    @test a.sampled_ids == b.sampled_ids
end

@testset "sample_leaves errors on out-of-range n (HR-8)" begin
    root = sampling_fixture()
    @test_throws ArgumentError sample_leaves(root, 0; seed = UInt64(1))
    @test_throws ArgumentError sample_leaves(root, 4; seed = UInt64(1))
end

@testset "sample_leaves accepts a Population" begin
    pop = simple_sampling_pop()
    s = sample_leaves(pop, 10; seed = UInt64(3))
    @test s.N_full == popsize(pop)
    @test s.n == 10
    @test all(id -> haskey(pop.cells, id), s.sampled_ids)
end
