# ── Mutational burden ─────────────────────────────────────────────────────────

"""
    mutations_per_cell(population::Population) -> Vector{Int64}

Total accumulated driver mutations of every alive cell (its whole lineage back to the
root), in `allcells(population)` order.
"""
mutations_per_cell(population::Population) =
    Int64[node.data.total_mutations for node in values(population.cells)]

"""
    mutations_per_cell(root::BinaryNode; includeclonal = false) -> Vector{Int64}

Driver burden of every leaf under `root`, in `Leaves(root)` order.

- `includeclonal = false` (default): only mutations acquired strictly *below* `root`.
  Mutations on `root` itself and on its ancestors are carried by every leaf of the
  subtree, so they are clonal here and are left out.
- `includeclonal = true`: each leaf's full burden, back to the top of its tree.

For the root of a simulated tree the two differ only by the founder's own mutations,
which are `0` for a population from [`initialize_population`](@ref).
"""
function mutations_per_cell(root::BinaryNode{NonMarkovCell}; includeclonal::Bool = false)
    offset = includeclonal ? 0 : root.data.total_mutations
    return Int64[leaf.data.total_mutations - offset for leaf in _leaves(root)]
end

"""
    average_mutations(population::Population) -> Float64

Mean accumulated driver mutations across all alive cells.
"""
average_mutations(population::Population) = mean(mutations_per_cell(population))

"""
    clonal_mutations(population::Population) -> Int64

Number of driver mutations shared by every alive cell (the MRCA's total burden), or `0`
if the cells have no common ancestor.
"""
function clonal_mutations(population::Population)
    MRCA = findMRCA(population)
    isnothing(MRCA) && return 0
    return MRCA.data.total_mutations
end

"""
    fitness_per_cell(population::Population) -> Vector{Float64}

Fitness of each currently alive cell, in `allcells(population)` order.
"""
fitness_per_cell(population::Population) =
    [node.data.fitness for node in values(population.cells)]

"""
    fitness_distribution(population::Population) -> Vector{Float64}

Alias for `fitness_per_cell`; returns fitness values of all living cells.
"""
fitness_distribution(population::Population) = fitness_per_cell(population)

"""
    mean_k(population::Population) -> Float64

Mean total driver-mutation count per alive cell.
"""
mean_k(population::Population) = mean(Float64.(mutations_per_cell(population)))

"""
    var_k(population::Population) -> Float64

Variance in total driver-mutation count across alive cells.
"""
var_k(population::Population) = var(Float64.(mutations_per_cell(population)))

# ── Distances and coalescence ─────────────────────────────────────────────────

# `f(cells[i], cells[j])::T` for every pair i < j, in row-major order; `idx`, if given,
# restricts `cells` first.
function _pairwise(f, ::Type{T}, cells::AbstractVector, idx = nothing) where T
    isnothing(idx) || (cells = cells[idx])
    n   = length(cells)
    out = Vector{T}(undef, n * (n - 1) ÷ 2)
    k   = 0
    for i in 1:n, j in i+1:n
        out[k += 1] = f(cells[i], cells[j])
    end
    return out
end

"""
    pairwisedistance(node1, node2) -> Int64

Number of driver mutations that differ between two cells: those on the path from each
cell up to their MRCA, excluding the MRCA itself. Cells in different trees share
nothing, so their distance is the sum of both full burdens.
"""
function pairwisedistance(cellnode1::BinaryNode, cellnode2::BinaryNode)
    mrca   = findMRCA(cellnode1, cellnode2)
    shared = isnothing(mrca) ? 0 : mrca.data.total_mutations
    return cellnode1.data.total_mutations + cellnode2.data.total_mutations - 2shared
end

"""
    pairwisedistances(population::Population[, idx]) -> Vector{Int64}

All pairwise distances between alive cells as a flat vector. `idx` indexes into
`allcells(population)`.
"""
pairwisedistances(population::Population, idx = nothing) =
    _pairwise(pairwisedistance, Int64, allcells(population), idx)

"""
    pairwise_differences(population::Population[, idx]) -> Dict{Int64,Int64}

Histogram of pairwise mutational distances between alive cells.
"""
pairwise_differences(population::Population, idx = nothing) =
    countmap(pairwisedistances(population, idx))

function _treeroot(node::BinaryNode)
    while !isnothing(node.parent)
        node = node.parent
    end
    return node
end

# Time from `t` back to the division of the two cells' MRCA. For cells in different
# trees this is, by convention, the time back to the earlier of the two founders' births.
function _coalescence_time(node1::BinaryNode, node2::BinaryNode, t::Real)
    node1 === node2 && return 0.0
    mrca = findMRCA(node1, node2)
    if isnothing(mrca)
        return t - min(_treeroot(node1).data.birthtime, _treeroot(node2).data.birthtime)
    end
    return t - endtime(mrca)
end

"""
    coalescence_times(root[, idx]; t) -> Vector{Float64}
    coalescence_times(population[, idx]; t) -> Vector{Float64}

Time back to the MRCA's division for every pair of alive cells. `t` defaults to
`age(root)` (last division) for the root method and to `pop.t` for the population
method. `idx` indexes into `getalivecells(root)` or `allcells(population)`.
"""
function coalescence_times(root::BinaryNode, idx = nothing; t = nothing)
    tref = isnothing(t) ? age(root) : t
    return _pairwise((a, b) -> _coalescence_time(a, b, tref), Float64, getalivecells(root), idx)
end

function coalescence_times(population::Population, idx = nothing; t = nothing)
    tref = isnothing(t) ? age(population) : t
    return _pairwise((a, b) -> _coalescence_time(a, b, tref), Float64, allcells(population), idx)
end

# ── Spectra ───────────────────────────────────────────────────────────────────

"""
    sitefrequencyspectrum(population::Population) -> Vector{Int64}

Driver-mutation site-frequency spectrum: `sfs[k]` = number of driver mutation events
present in exactly `k` living cells. Length `popsize(population)`. On a forest every
tree contributes.
"""
function sitefrequencyspectrum(population::Population)
    sfs = zeros(Int64, popsize(population))
    for root in _population_roots(population)
        _sfs_fill!(sfs, _preorder(root)...)
    end
    return sfs
end

function _sfs_fill!(sfs::Vector{Int64}, nodes, parent)
    counts = _leafcounts(nodes, parent)
    for (node, c) in zip(nodes, counts)
        c > 0 && (sfs[c] += node.data.mutations)
    end
    return counts[1]
end

function _check_spectrum_length(nleaves::Int, N::Int, what::String)
    nleaves <= N || throw(ArgumentError(
        "$what: tree has $nleaves leaves but N = $N — the spectrum would " *
        "overflow. Pass N >= $nleaves."))
    return nothing
end

"""
    sitefrequencyspectrum(root::BinaryNode[, N::Int]) -> Vector{Int64}

Site-frequency spectrum of a lineage tree: `sfs[k]` is the number of mutation
events carried by exactly `k` of the tree's leaves. The returned vector has
length `N`, which defaults to the tree's leaf count `n`.

Pass `N` explicitly when the spectrum must have a particular length — notably for
a tree returned by [`sample_leaves`](@ref), where the meaningful length is the
sample size. `N > n` pads with zeros; `N < n` is an error.

Mutations on the root's own edge are accumulated into `sfs[n]`: they are clonal in
the given tree.
"""
function sitefrequencyspectrum(root::BinaryNode, N::Int)
    nodes, parent = _preorder(root)
    counts = _leafcounts(nodes, parent)
    _check_spectrum_length(counts[1], N, "sitefrequencyspectrum")
    sfs = zeros(Int64, N)
    for (node, c) in zip(nodes, counts)
        c > 0 && (sfs[c] += node.data.mutations)
    end
    return sfs
end

sitefrequencyspectrum(root::BinaryNode) = sitefrequencyspectrum(root, popsize(root))

"""
    branch_spectrum(root::BinaryNode[, N::Int]) -> Vector{Int}

Topological site-frequency spectrum: `bs[k]` is the number of internal nodes
subtending exactly `k` leaves. Leaves themselves are not counted.

For neutral mutations this encodes the full tree topology, separating it from the
mutation rate: under any neutral per-division rate `m`,

    E[sfs[k]] = m * bs[k]           for k >= 2
    E[sfs[1]] = m * (bs[1] + n)     for k == 1

where `n` is the tree's leaf count. The `k = 1` case needs the extra `n`
because a leaf's own edge also carries mutations into `sfs[1]`, while `bs`
counts only internal nodes — so `bs[1]` is the number of *unary* internal nodes
(a division whose other daughter's lineage died out), not the number of leaves.

Length and `N` semantics match [`sitefrequencyspectrum`](@ref).
"""
function branch_spectrum(root::BinaryNode, N::Int)
    nodes, parent = _preorder(root)
    counts = _leafcounts(nodes, parent)
    _check_spectrum_length(counts[1], N, "branch_spectrum")
    bs = zeros(Int, N)
    for (node, c) in zip(nodes, counts)
        haschildren(node) && c > 0 && (bs[c] += 1)
    end
    return bs
end

branch_spectrum(root::BinaryNode) = branch_spectrum(root, popsize(root))

# ── Leaf divisional depths ────────────────────────────────────────────────────

"""
    leaf_depths(root::BinaryNode) -> Vector{Int}

Number of division events on the path from `root` to each leaf.

For neutral simulations (`ν = 0`) this is the primary quantity: mutations per cell
follow by drawing `Poisson(m * depth)` per leaf afterwards.

!!! warning "Not co-indexed"
    The returned vector is in this function's own stack order (right subtrees first),
    which is *not* `Leaves(root)` order, so it is valid only as a pooled distribution.
    The order is kept fixed because stored results were produced with it.
"""
function leaf_depths(root::BinaryNode{T}) where {T <: AbstractTreeCell}
    depths = Int[]
    stack  = Tuple{BinaryNode{T}, Int}[(root, 0)]
    while !isempty(stack)
        node, d = pop!(stack)
        if isleaf(node)
            push!(depths, d)
        else
            isnothing(node.left)  || push!(stack, (node.left,  d + 1))
            isnothing(node.right) || push!(stack, (node.right, d + 1))
        end
    end
    return depths
end

# ── Filtered burdens and leaf fitness ────────────────────────────────────────

"""
    filtered_mutations_per_cell(root::BinaryNode, threshold::Float64) -> Vector{Int}

Per-leaf mutational burden, counted from the leaf up to and including `root`, excluding
mutations on any ancestral node whose leaf count exceeds `floor(threshold * N)`, where
`N` is the leaf count of the tree under `root`. A leaf's own mutations always count.

This is the tree-side analogue of dropping high-frequency variants before
estimating a mutation rate: a node subtending a large fraction of the population
contributes the same mutations to every cell below it, so it carries no
information about within-clone divergence. `threshold = 1.0` excludes nothing, and for
the top of a tree equals `mutations_per_cell(root; includeclonal = true)`.

Returned in `Leaves(root)` order.
"""
function filtered_mutations_per_cell(root::BinaryNode{T},
                                     threshold::Float64) where {T <: AbstractTreeCell}
    nodes, parent = _preorder(root)
    counts    = _leafcounts(nodes, parent)
    max_count = floor(Int, threshold * counts[1])

    # Forward sweep: parents precede children, so `above[p]` is final when read.
    above  = zeros(Int, length(nodes))   # filtered burden down to and including node i
    result = Int[]
    for i in eachindex(nodes)
        node    = nodes[i]
        inherit = parent[i] == 0 ? 0 : above[parent[i]]
        if isleaf(node)
            push!(result, inherit + node.data.mutations)
        else
            above[i] = inherit + (counts[i] <= max_count ? node.data.mutations : 0)
        end
    end
    return result
end

"""
    leaf_fitness(root::BinaryNode) -> Vector{Float64}

Fitness of every alive leaf, in `getalivecells(root)` order — the same order
[`mutations_per_cell`](@ref) uses, so the two are **co-indexed**: entry `i` is the
same cell in both. Note that [`leaf_depths`](@ref) is *not* co-indexed with either.
"""
leaf_fitness(root::BinaryNode) = [leaf.data.fitness for leaf in _leaves(root)]
