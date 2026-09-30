@isdefined(fixture_tree) || include("fixtures.jl")

simple_pop(; kwargs...) = grown_population(; kwargs...)

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

# ── Hand-built fixture (see fixtures.jl) ──────────────────────────────────────

@testset "fixture has the expected shape" begin
    root = fixture_tree()
    @test [l.data.id for l in Leaves(root)] == [4, 5, 3]
    @test mutations_per_cell(root; includeclonal = true) == [8, 9, 12]
end

@testset "includeclonal: false drops the root's own mutations, true keeps everything" begin
    root = fixture_tree()
    L    = root.left
    @test mutations_per_cell(root) == [3, 4, 7]        # root's 5 are clonal
    @test mutations_per_cell(L) == [2, 3]              # L's 1 and root's 5 are clonal
    @test mutations_per_cell(L; includeclonal = true) == [8, 9]
end

@testset "distances, MRCA and coalescence on the fixture" begin
    root = fixture_tree()
    LL, LR, R = root.left.left, root.left.right, root.right
    @test findMRCA(LL, LR) === root.left
    @test findMRCA(LL, R) === root
    @test findMRCA([LL, LR, R]) === root
    @test findMRCA(LL, LL) === LL
    @test pairwisedistance(LL, LR) == 5
    @test pairwisedistance(LL, R)  == 10
    @test pairwisedistance(LR, R)  == 11
    @test pairwisedistance(R, R)   == 0
    # Leaves order [LL, LR, R]; the MRCA of LL and LR divided at 2.0, the root at 1.0.
    @test coalescence_times(root; t = 3.0) == [1.0, 2.0, 2.0]
    @test coalescence_times(root, [1, 3]; t = 3.0) == [2.0]
    @test celllifetimes(root) == [1.0, 1.0]            # root 0→1, L 1→2
    @test age(root) == 2.1
end

@testset "tree statistics agree with independent walks on a simulated tree" begin
    pop  = simple_pop(ν = 2.0, Nmax = 150)
    root = getsingleroot(allcells(pop))
    @test MutationLoadDynamics._leaves(root) == collect(Leaves(root))
    @test getalivecells(root) == collect(Leaves(root))
    burden = id_burden_map(root)
    @test mutations_per_cell(root; includeclonal = true) ==
          [burden[l.data.id] for l in Leaves(root)]
    @test sort(mutations_per_cell(pop)) == sort(collect(values(burden)))
    @test sort(leaf_depths(root)) == sort(collect(values(id_depth_map(root))))
    @test clonal_mutations(pop) == findMRCA(pop).data.total_mutations
    cells = getalivecells(root)[1:12]
    for a in cells, b in cells
        m = findMRCA(a, b)
        @test pairwisedistance(a, b) == burden[a.data.id] + burden[b.data.id] -
                                         2 * m.data.total_mutations
    end
end

@testset "leaf_depths on the fixture" begin
    root = fixture_tree()
    # Stack order, not Leaves order: the function pushes left then right and pops
    # last-in-first-out, so R (depth 1) is emitted before LL and LR (depth 2).
    @test leaf_depths(root) == [1, 2, 2]
end

@testset "leaf_depths on a single-node tree" begin
    root = rootnode(1, 0.0, 3)
    @test leaf_depths(root) == [0]
end

@testset "leaf_depths length matches the leaf count" begin
    pop  = simple_pop(Nmax = 40)
    root = getsingleroot(allcells(pop))
    @test length(leaf_depths(root)) == popsize(pop)
    @test all(d >= 0 for d in leaf_depths(root))
end

@testset "root-based sitefrequencyspectrum on the fixture" begin
    root = fixture_tree()
    # sfs[1] = 2 + 3 + 7 (the three leaves' own mutations)
    # sfs[2] = 1         (L subtends LL and LR)
    # sfs[3] = 5         (the root subtends all three leaves)
    @test sitefrequencyspectrum(root, 3) == [12, 1, 5]
    @test sitefrequencyspectrum(root)    == [12, 1, 5]
    @test sitefrequencyspectrum(root.left) == [5, 1]   # a subtree is a tree too
end

@testset "SFS on a forest counts every root's tree" begin
    # `initialize_population(N)` seeds N independent founders, so the population is a
    # forest and `sitefrequencyspectrum` takes its multi-root branch. That branch has to
    # traverse the N *roots*; a version that traversed one representative alive *cell*
    # per tree instead would fill only `sfs[1]`, from those leaves' own mutations, and
    # silently lose almost every mutation in the population.
    block = NonMarkovBlock(
        birth_dist     = f -> Gamma(2.0, 1.0 / f),
        death_dist     = f -> Gamma(2.0, 20.0),
        stopfunction   = pop -> popsize(pop) >= 200,
        driver_dist    = Exponential(0.1),
        fitness_update = (f, δ) -> f + δ,
        ν              = 1.0,
    )
    pop = initialize_population(20; fitness_init = 1.0)
    simulate!(pop, block, MersenneTwister(11))

    @test isnothing(getsingleroot(allcells(pop)))    # it really is a forest
    sfs = sitefrequencyspectrum(pop)
    @test length(sfs) == popsize(pop)

    # Same identity the single-root test asserts: summing k-weighted spectrum entries
    # recovers the total burden carried by the living cells.
    N = popsize(pop)
    @test sum(sfs[k] * k for k in 1:N) == sum(mutations_per_cell(pop))
    @test sum(sfs) > popsize(pop) / 10   # not just the leaves' own mutations
end

@testset "root-based sitefrequencyspectrum agrees with the population method" begin
    pop  = simple_pop(ν = 2.0, Nmax = 40)
    root = getsingleroot(allcells(pop))
    @test sitefrequencyspectrum(root, popsize(pop)) == sitefrequencyspectrum(pop)
end

@testset "sitefrequencyspectrum accepts N larger than the leaf count" begin
    root = fixture_tree()
    # Padding with zeros is legitimate: a caller may want a fixed-length spectrum.
    @test sitefrequencyspectrum(root, 5) == [12, 1, 5, 0, 0]
end

@testset "sitefrequencyspectrum errors when N is too small" begin
    root = fixture_tree()
    @test_throws ArgumentError sitefrequencyspectrum(root, 2)
end

@testset "branch_spectrum on the fixture" begin
    root = fixture_tree()
    # L subtends 2 leaves, root subtends 3; leaves themselves are not counted.
    @test branch_spectrum(root, 3) == [0, 1, 1]
    @test branch_spectrum(root)    == [0, 1, 1]
end

@testset "branch_spectrum counts every internal node exactly once" begin
    pop  = simple_pop(Nmax = 40)
    root = getsingleroot(allcells(pop))
    N    = popsize(pop)
    bs   = branch_spectrum(root, N)
    n_internal = count(n -> !isnothing(n.left) || !isnothing(n.right), PreOrderDFS(root))
    @test sum(bs) == n_internal
    # >= rather than == : simple_pop has a non-zero death rate, so prune_tree!
    # can leave a unary node that subtends every leaf just as the root does.
    @test bs[N] >= 1          # the root subtends every leaf
end

@testset "branch_spectrum errors when N is too small" begin
    root = fixture_tree()
    @test_throws ArgumentError branch_spectrum(root, 2)
end

@testset "leaf_fitness on the fixture" begin
    root = fixture_tree()
    @test leaf_fitness(root) == [1.0, 1.0, 1.0]
end

@testset "leaf_fitness is co-indexed with mutations_per_cell" begin
    root = fixture_tree()
    # Both iterate getalivecells(root), so entry i is the same cell in both.
    ids = [l.data.id for l in getalivecells(root)]
    @test ids == [4, 5, 3]
    @test length(leaf_fitness(root)) == length(mutations_per_cell(root))
end

@testset "leaf_fitness picks up a driver" begin
    root = fixture_tree()
    set_fitness!(first(getalivecells(root)), 1.5)
    @test leaf_fitness(root) == [1.5, 1.0, 1.0]
end

@testset "filtered_mutations_per_cell with threshold 1.0 equals the full burden" begin
    root = fixture_tree()
    # floor(1.0 * 3) = 3, so no node is excluded.
    @test filtered_mutations_per_cell(root, 1.0) == mutations_per_cell(root; includeclonal = true)
end

@testset "filtered_mutations_per_cell excludes the root at a low threshold" begin
    root = fixture_tree()
    # floor(0.5 * 3) = 1, so only nodes with <= 1 live descendant contribute from
    # the ancestry: the root (3 descendants) and L (2 descendants) are excluded, so
    # each leaf keeps only its own mutations: [2, 3, 7].
    @test filtered_mutations_per_cell(root, 0.5) == [2, 3, 7]
end

@testset "filtered_mutations_per_cell is bounded by the full burden" begin
    pop  = simple_pop(ν = 2.0, Nmax = 40)
    root = getsingleroot(allcells(pop))
    full = mutations_per_cell(root; includeclonal = true)
    filt = filtered_mutations_per_cell(root, 0.3)
    @test length(filt) == length(full)
    @test all(filt .<= full)
end

@testset "filtered_mutations_per_cell on a subtree stops at the subtree root" begin
    root = fixture_tree()
    @test filtered_mutations_per_cell(root.left, 1.0) == [3, 4]   # L's 1 counted, root's not
    @test filtered_mutations_per_cell(root.left, 0.5) == [2, 3]   # L subtends 2 > floor(1)
end

@testset "forests: MRCA, distance and coalescence across trees" begin
    # Tree B's root has a larger id than a non-root node of tree A. The old climb
    # followed `parent` into `nothing` and crashed here.
    a  = rootnode(1, 0.0, 2)
    a1 = child!(:left, a, 2, 1.0, 1); child!(:right, a, 3, 1.0, 0)
    b  = rootnode(10, 0.5, 4)
    @test isnothing(findMRCA(a1, b))
    @test isnothing(findMRCA(b, a1))
    @test pairwisedistance(a1, b) == 3 + 4
    @test MutationLoadDynamics._coalescence_time(a1, b, 5.0) == 5.0   # back to t = 0
    pop = initialize_population(4)
    @test isnothing(findMRCA(pop))
    @test clonal_mutations(pop) == 0
    @test length(coalescence_times(pop)) == 6
    @test length(MutationLoadDynamics._roots(allcells(pop))) == 4
end

@testset "pairwise helpers handle fewer than two cells" begin
    pop = initialize_population()
    @test pairwisedistances(pop) == Int64[]
    @test coalescence_times(pop) == Float64[]
    @test isempty(pairwise_differences(pop))
end
