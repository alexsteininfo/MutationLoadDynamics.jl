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

@testset "GOLDEN LITERALS: frozen draw output (HR-1)" begin
    # Captured once from a passing run under Julia 1.12.1 and FROZEN. These are
    # literal expected values, not a re-derivation of the recipe — that is the
    # point. If this testset fails, the draw recipe changed and every serialized
    # sampled tree on disk became unreproducible. Fix the code; never edit these
    # numbers to match new behaviour.
    #
    # The companion "GOLDEN" testset above re-derives the recipe and catches a
    # different class of change (a wrong RNG or sampler). Keep both.
    #
    # Julia changed the algorithm behind `rand`/`randperm` range sampling in 1.11 (the CI
    # floor is 1.9), so a `MersenneTwister` seed does not reproduce the same draws across
    # that boundary. These literals are only meaningful on the Julia line they were
    # captured on, so they are skipped below 1.11 rather than pinned to a stream that
    # cannot match.
    root = sampling_fixture()
    if VERSION >= v"1.11"
        @test sample_leaves(root, 1; seed = UInt64(1)).sampled_ids  == Int64[3]
        @test sample_leaves(root, 2; seed = UInt64(1)).sampled_ids  == Int64[3, 5]
        @test sample_leaves(root, 2; seed = UInt64(42)).sampled_ids == Int64[5, 3]
        @test sample_leaves(root, 3; seed = UInt64(7)).sampled_ids  == Int64[5, 3, 4]
    else
        @test_skip sample_leaves(root, 1; seed = UInt64(1)).sampled_ids  == Int64[3]
        @test_skip sample_leaves(root, 2; seed = UInt64(1)).sampled_ids  == Int64[3, 5]
        @test_skip sample_leaves(root, 2; seed = UInt64(42)).sampled_ids == Int64[5, 3]
        @test_skip sample_leaves(root, 3; seed = UInt64(7)).sampled_ids  == Int64[5, 3, 4]
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

@testset "sample_leaves on a single-node tree (M-4)" begin
    # The bounds test the spec (§5) called for but the suite never got: a tree
    # with exactly one leaf, which is also its own root.
    root = BinaryNode(NonMarkovCell(1, 0.0, 0, 1.0))
    s = sample_leaves(root, 1; seed = UInt64(1))
    @test s.sampled_ids == Int64[1]
    @test s.N_full == 1
    @test s.root.data.id == 1
    @test_throws ArgumentError sample_leaves(root, 2; seed = UInt64(1))
end

@testset "sample_leaves accepts a Population" begin
    pop = simple_sampling_pop()
    s = sample_leaves(pop, 10; seed = UInt64(3))
    @test s.N_full == popsize(pop)
    @test s.n == 10
    @test all(id -> haskey(pop.cells, id), s.sampled_ids)
end

# ── Id-keyed reference maps ───────────────────────────────────────────────────
# leaf_depths is not co-indexed with anything, so comparing a sample cell-by-cell
# against the full tree needs id keys.
function id_depth_map(root::BinaryNode{NonMarkovCell})
    m = Dict{Int64, Int}()
    stack = Tuple{BinaryNode{NonMarkovCell}, Int}[(root, 0)]
    while !isempty(stack)
        node, d = pop!(stack)
        if isnothing(node.left) && isnothing(node.right)
            m[node.data.id] = d
        else
            isnothing(node.left)  || push!(stack, (node.left,  d + 1))
            isnothing(node.right) || push!(stack, (node.right, d + 1))
        end
    end
    return m
end

function id_burden_map(root::BinaryNode{NonMarkovCell})
    m = Dict{Int64, Int}()
    for leaf in Leaves(root)
        muts = leaf.data.mutations
        node = leaf
        while !isnothing(node.parent)
            node = node.parent
            muts += node.data.mutations
        end
        m[leaf.data.id] = muts
    end
    return m
end

@testset "INVARIANCE: burden and depth are the cell's full-tree values" begin
    # The test that catches an accidental unary collapse. If ancestors were merged
    # away, a sampled cell's depth would shrink and its burden would lose the
    # mutations on the collapsed edges.
    pop    = simple_sampling_pop(Nmax = 200)
    root   = getsingleroot(allcells(pop))
    burden = id_burden_map(root)
    depth  = id_depth_map(root)

    for n in (1, 5, 40, 200)
        s = sample_leaves(root, n; seed = UInt64(n))
        @test id_burden_map(s.root) == Dict(id => burden[id] for id in s.sampled_ids)
        @test id_depth_map(s.root)  == Dict(id => depth[id]  for id in s.sampled_ids)
        @test sort(leaf_depths(s.root)) == sort([depth[id] for id in s.sampled_ids])
    end
end

@testset "every retained node has a sampled leaf below it, and vice versa" begin
    root = sampling_fixture()
    s = sample_leaves(root, 2; seed = UInt64(5))
    drawn = Set(s.sampled_ids)

    # Direction 1: no node is retained without a sampled leaf beneath it.
    for node in PreOrderDFS(s.root)
        below = Set(l.data.id for l in Leaves(node))
        @test !isempty(intersect(below, drawn))
    end
    # Direction 2: the sampled leaves are exactly the leaves of the induced tree.
    @test Set(l.data.id for l in Leaves(s.root)) == drawn
end

@testset "n = 1 gives a path, and the founder stays the root" begin
    root  = sampling_fixture()
    depth = id_depth_map(root)
    seen  = Set{Int64}()
    for seed in UInt64(1):UInt64(60)
        s = sample_leaves(root, 1; seed = seed)
        push!(seen, s.sampled_ids[1])
        @test s.root.data.id == 1                     # founder retained as root
        @test isnothing(s.root.parent)
        @test length(collect(Leaves(s.root))) == 1
        @test leaf_depths(s.root) == [depth[s.sampled_ids[1]]]
        # A path: every node has exactly one child except the leaf.
        for node in PreOrderDFS(s.root)
            nkids = count(!isnothing, (node.left, node.right))
            @test nkids in (0, 1)
        end
    end
    @test seen == Set{Int64}([3, 4, 5])               # all leaves reachable
end

@testset "n = 2 keeps the branching node and drops the other side" begin
    root  = sampling_fixture()
    found = false
    # Draw until we get exactly {4, 5}: both under L, so R must be absent.
    for seed in UInt64(1):UInt64(200)
        s = sample_leaves(root, 2; seed = seed)
        if Set(s.sampled_ids) == Set{Int64}([4, 5])
            found = true
            @test isnothing(s.root.right)             # R dropped
            @test !isnothing(s.root.left)             # L retained, still binary
            @test !isnothing(s.root.left.left) && !isnothing(s.root.left.right)
            break
        end
    end
    # Without this the testset would silently assert nothing if no seed hit {4,5}.
    @test found
end

@testset "left/right slots are preserved, not re-balanced (HR-4)" begin
    root  = sampling_fixture()
    found = false
    for seed in UInt64(1):UInt64(200)
        s = sample_leaves(root, 2; seed = seed)
        if Set(s.sampled_ids) == Set{Int64}([3, 4])
            found = true
            # 4 is under L (left of root); 3 is R (right of root). L keeps its
            # left child and loses its right, and must NOT slide 4 into `right`.
            @test !isnothing(s.root.left) && !isnothing(s.root.right)
            @test s.root.left.left.data.id == 4
            @test isnothing(s.root.left.right)
            break
        end
    end
    @test found
end

@testset "parent links are consistent (HR-4)" begin
    pop  = simple_sampling_pop(Nmax = 80)
    root = getsingleroot(allcells(pop))
    s = sample_leaves(root, 20; seed = UInt64(9))
    @test isnothing(s.root.parent)
    for node in PreOrderDFS(s.root)
        isnothing(node.left)  || @test node.left.parent  === node
        isnothing(node.right) || @test node.right.parent === node
    end
end

@testset "the source tree is never mutated (HR-3)" begin
    root   = sampling_fixture()
    before = (mutations_per_cell(root), leaf_depths(root),
              sitefrequencyspectrum(root, 3), [l.data.id for l in Leaves(root)])
    sample_leaves(root, 1; seed = UInt64(3))
    sample_leaves(root, 2; seed = UInt64(4))
    sample_leaves(root, 3; seed = UInt64(5))
    @test (mutations_per_cell(root), leaf_depths(root),
           sitefrequencyspectrum(root, 3),
           [l.data.id for l in Leaves(root)]) == before
end

@testset "node data is shared, not copied (HR-5)" begin
    root = sampling_fixture()
    s = sample_leaves(root, 3; seed = UInt64(1))
    src = Dict(l.data.id => l.data for l in Leaves(root))
    for leaf in Leaves(s.root)
        @test leaf.data === src[leaf.data.id]
    end
end

@testset "draws at different n and different seeds are independent" begin
    root = sampling_fixture()
    a = sample_leaves(root, 2; seed = UInt64(42))
    c = sample_leaves(root, 2; seed = UInt64(2))
    @test a.sampled_ids != c.sampled_ids   # different seeds, different draw
    # A draw at n = 1 is not merely the first element of the n = 2 draw's ids
    # reversed or truncated in some coupled way — it is its own permutation.
    @test sample_leaves(root, 1; seed = UInt64(42)).sampled_ids == a.sampled_ids[1:1]
end

@testset "SFS of a sample sums to the sample's total burden" begin
    pop  = simple_sampling_pop(Nmax = 150)
    root = getsingleroot(allcells(pop))
    for n in (1, 10, 150)
        s   = sample_leaves(root, n; seed = UInt64(100 + n))
        sfs = sitefrequencyspectrum(s.root, n)
        @test sum(k * sfs[k] for k in 1:n) == sum(mutations_per_cell(s.root))
    end
end

@testset "SamplingSpec constructors express the three modes" begin
    @test SamplingSpec().sizes == Int[]
    @test SamplingSpec().retain_full
    @test SamplingSpec(100).sizes == [100]
    @test SamplingSpec([1000, 100]).sizes == [1000, 100]
    @test SamplingSpec(sizes = [10], retain_full = false).retain_full == false
    @test SamplingSpec().replicates == 1
end

@testset "SamplingSpec validates its inputs" begin
    @test_throws ArgumentError SamplingSpec(sizes = [10, 10])       # duplicate
    @test_throws ArgumentError SamplingSpec(sizes = [0])            # n < 1
    @test_throws ArgumentError SamplingSpec(sizes = [10], replicates = 0)
end

@testset "SamplingSpec cannot be constructed unvalidated" begin
    # The auto-generated positional constructor used to bypass every check,
    # which let two draws at the same (n, replicate) share a derived seed.
    @test_throws ArgumentError SamplingSpec([10, 10], 1, true)
    @test_throws ArgumentError SamplingSpec([0], 1, true)
    @test_throws ArgumentError SamplingSpec([10], 0, true)
end

@testset "sample_trees: full data only" begin
    root = sampling_fixture()
    out  = sample_trees(root, SamplingSpec(); seed = UInt64(1))
    @test out isa SampledTrees
    @test out.full === root
    @test isempty(out.samples)
end

@testset "sample_trees: one sample size, no full tree" begin
    root = sampling_fixture()
    out  = sample_trees(root, SamplingSpec(sizes = [2], retain_full = false);
                        seed = UInt64(1))
    @test isnothing(out.full)
    @test length(out.samples) == 1
    @test out.samples[1].n == 2
end

@testset "sample_trees: full tree plus several sizes" begin
    root = sampling_fixture()
    out  = sample_trees(root, SamplingSpec([3, 2, 1]); seed = UInt64(1))
    @test out.full === root
    @test [s.n for s in out.samples] == [3, 2, 1]
    @test allunique([s.seed for s in out.samples])
end

@testset "sample_trees: replicates give distinct draws and distinct seeds" begin
    pop = simple_sampling_pop(Nmax = 100)
    out = sample_trees(pop, SamplingSpec(sizes = [10], replicates = 3);
                       seed = UInt64(77))
    @test length(out.samples) == 3
    @test [s.replicate for s in out.samples] == [1, 2, 3]
    @test allunique([s.seed for s in out.samples])
    @test allunique([s.sampled_ids for s in out.samples])
    @test all(s -> s.n == 10, out.samples)
end

@testset "sample_trees derives seeds reproducibly from the base seed" begin
    root = sampling_fixture()
    a = sample_trees(root, SamplingSpec([2, 1]); seed = UInt64(5))
    b = sample_trees(root, SamplingSpec([2, 1]); seed = UInt64(5))
    @test [s.sampled_ids for s in a.samples] == [s.sampled_ids for s in b.samples]
    # And each stored seed replays its own draw in isolation.
    for s in a.samples
        @test sample_leaves(root, s.n; seed = s.seed).sampled_ids == s.sampled_ids
    end
end

@testset "sample_trees validates n against the tree before drawing" begin
    root = sampling_fixture()
    # All sizes are validated before any draw is made, regardless of order.
    @test_throws ArgumentError sample_trees(root, SamplingSpec([9, 1]);
                                           seed = UInt64(1))
end

@testset "END TO END: grow, sample, and check against the population" begin
    pop  = simple_sampling_pop(Nmax = 300, ν = 2.0, rng = MersenneTwister(2024))
    out  = sample_trees(pop, SamplingSpec([300, 30, 3]); seed = UInt64(0xC0FFEE))

    @test !isnothing(out.full)
    @test length(out.samples) == 3
    for s in out.samples
        @test s.N_full == popsize(pop)
        @test length(s.sampled_ids) == s.n
        @test allunique(s.sampled_ids)
        @test all(id -> haskey(pop.cells, id), s.sampled_ids)
        @test length(collect(Leaves(s.root))) == s.n
    end

    # n = N_full is the whole population: every statistic must match exactly.
    whole = out.samples[1]
    @test whole.n == popsize(pop)
    @test sort(mutations_per_cell(whole.root)) == sort(mutations_per_cell(pop))
    @test sitefrequencyspectrum(whole.root, popsize(pop)) ==
          sitefrequencyspectrum(pop)
    @test sort(leaf_depths(whole.root)) ==
          sort(leaf_depths(getsingleroot(allcells(pop))))
end

@testset "END TO END: sampling a tree that death has pruned" begin
    # d > 0 means prune_tree! has removed dead lineages, leaving unary nodes in the
    # stored tree already. A sampled tree must be indistinguishable in kind.
    block = NonMarkovBlock(
        birth_dist     = f -> Gamma(5.0, 1.0 / (5.0 * f)),
        death_dist     = f -> Gamma(5.0, 1.0 / (5.0 * 0.5)),
        stopfunction   = pop -> popsize(pop) >= 200,
        driver_dist    = Dirac(0.0),
        fitness_update = (f, δ) -> f,
        ν              = 2.0,
        restart_on_extinction = true,
    )
    pop = initialize_population(fitness_init = 1.0)
    simulate!(pop, block, MersenneTwister(7))

    root   = getsingleroot(allcells(pop))
    burden = id_burden_map(root)
    s = sample_leaves(root, 20; seed = UInt64(31))

    @test id_burden_map(s.root) == Dict(id => burden[id] for id in s.sampled_ids)
    @test length(collect(Leaves(s.root))) == 20
    @test isnothing(s.root.parent)
    for node in PreOrderDFS(s.root)
        isnothing(node.left)  || @test node.left.parent  === node
        isnothing(node.right) || @test node.right.parent === node
    end
end

@testset "sampling a forest is rejected by name (HR-8)" begin
    # Two independent founders: getsingleroot returns nothing, and both the
    # single-draw and the spec-driven entry points must say so rather than
    # silently sampling one tree of the forest.
    pop = initialize_population(fitness_init = 1.0)
    second = BinaryNode(NonMarkovCell(pop._next_id + 1, 0.0, 0, 1.0))
    pop.cells[second.data.id] = second
    pop._next_id += 1

    @test_throws ArgumentError sample_leaves(pop, 1; seed = UInt64(1))
    @test_throws ArgumentError sample_trees(pop, SamplingSpec(1); seed = UInt64(1))
    err = try sample_leaves(pop, 1; seed = UInt64(1)) catch e; e end
    @test occursin("2 independent roots", err.msg)
    err_trees = try sample_trees(pop, SamplingSpec(1); seed = UInt64(1)) catch e; e end
    @test occursin("2 independent roots", err_trees.msg)
end
