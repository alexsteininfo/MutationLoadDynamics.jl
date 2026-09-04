# Uniform leaf sampling of lineage trees.
#
# The induced tree keeps the sampled leaves plus *every ancestor* of a sampled
# leaf, and retains the resulting unary nodes rather than collapsing them. That is
# the load-bearing decision: because every division ancestral to a sampled cell is
# still a node, a sampled cell's root-to-leaf path is unchanged, so
# `mutations_per_cell` and `leaf_depths` return exactly that cell's *full-tree*
# burden and divisional depth. Collapsing would turn depth into a count of
# bifurcations that survived sampling — a property of the sample rather than of the
# cell. It is also the shape `prune_tree!` already leaves behind when a lineage
# dies out, so a sampled tree is indistinguishable in kind from a full one and
# every statistic in `statistics.jl` applies to it unchanged.

"""
    LeafSample

One uniform draw of `n` leaves from a lineage tree, with the induced tree.

# Fields
- `root::BinaryNode{NonMarkovCell}` — induced tree: the sampled leaves plus every
  ancestor of a sampled leaf, unary nodes retained. The original founder remains
  the root even when it has a single child, so `sitefrequencyspectrum` accumulates
  its mutations into `sfs[n]` exactly as it does into `sfs[N]` for a full tree.
  **The root is therefore the founder, not the MRCA of the sample** — call
  `findMRCA` on the sampled leaves if you need that.
- `n::Int` — cells drawn.
- `N_full::Int` — leaf count of the source tree.
- `seed::UInt64` — rng seed of this draw. The draw is a pure function of
  `(tree, n, seed)`, so any single draw replays in isolation from this record.
- `replicate::Int` — which independent draw this is at this `(tree, n)`. Defaults
  to 1; distinct replicates require distinct seeds, which
  [`sample_trees`](@ref) derives for you.
- `sampled_ids::Vector{Int64}` — `NonMarkovCell.id` of each drawn cell, in draw
  order.

The induced tree **shares** `node.data` with the source tree — `NonMarkovCell` is
immutable, so ids, birthtimes, mutation counts and fitness are identical by
construction rather than by copy. Replace a node's `data` wholesale if you need to
change it (that rebinds only the sampled tree); never mutate a cell in place.
"""
struct LeafSample
    root::BinaryNode{NonMarkovCell}
    n::Int
    N_full::Int
    seed::UInt64
    replicate::Int
    sampled_ids::Vector{Int64}
end

# Copy the marked part of the tree. Left/right slots are preserved, so a node whose
# left lineage was dropped keeps `left = nothing` — the tree stays a faithful
# sub-shape of the original rather than being silently re-balanced.
function _copy_marked(node::BinaryNode{NonMarkovCell},
                      marked::Set{BinaryNode{NonMarkovCell}})
    new = BinaryNode(node.data)
    if !isnothing(node.left) && node.left in marked
        new.left = _copy_marked(node.left, marked)
        new.left.parent = new
    end
    if !isnothing(node.right) && node.right in marked
        new.right = _copy_marked(node.right, marked)
        new.right.parent = new
    end
    return new
end

"""
    sample_leaves(root, n; seed, replicate = 1) -> LeafSample
    sample_leaves(population, n; seed, replicate = 1) -> LeafSample

Draw `n` of the tree's leaves uniformly without replacement and return the induced
lineage tree as a [`LeafSample`](@ref).

`seed` is required and supplied by the caller: the draw must be reproducible from
values the caller itself records, and callers derive seeds from their own
provenance (a filename stem, a simulation index). Two draws that should differ
must be given different seeds.

Non-destructive — `root` is left untouched, because the same tree is normally
drawn from again at other sample sizes and those draws must be independent.
Single-threaded by design: each call builds its own rng from `seed`, so
independent draws are safe to run concurrently from the caller's own threads.

Cost is proportional to the number of *retained* nodes, not to the size of the
tree.

!!! warning "The draw recipe is frozen"
    `MersenneTwister(seed)` and `randperm(rng, N_full)[1:n]` over `Leaves(root)`
    order are load-bearing, not incidental: serialized sampled trees produced by
    earlier versions of this algorithm must stay reproducible. Do not substitute a
    different rng, `StatsBase.sample`, or a direct `n`-index draw.
"""
function sample_leaves(root::BinaryNode{NonMarkovCell}, n::Int;
                       seed::UInt64, replicate::Int = 1)
    leaves = collect(Leaves(root))
    N_full = length(leaves)
    1 <= n <= N_full || throw(ArgumentError(
        "cannot draw n = $n cells from a tree with $N_full leaves"))
    replicate >= 1 || throw(ArgumentError("replicate must be >= 1, got $replicate"))

    # `Leaves` visits a fixed tree in a fixed order (`children` returns
    # `(left, right)`), so the draw is a pure function of (tree, seed, n).
    rng = MersenneTwister(seed)
    idx = randperm(rng, N_full)[1:n]

    # Mark each sampled leaf and its ancestors, stopping at the first node already
    # marked: total cost is the number of retained nodes, not the size of the tree.
    # `BinaryNode` is mutable, so `Set` compares by identity.
    marked = Set{BinaryNode{NonMarkovCell}}()
    for i in idx
        node = leaves[i]
        while !isnothing(node) && !(node in marked)
            push!(marked, node)
            node = node.parent
        end
    end

    new_root    = _copy_marked(root, marked)
    sampled_ids = Int64[leaves[i].data.id for i in idx]
    return LeafSample(new_root, n, N_full, seed, replicate, sampled_ids)
end

function sample_leaves(population::Population, n::Int;
                       seed::UInt64, replicate::Int = 1)
    root = getsingleroot(allcells(population))
    if isnothing(root)
        nroots = length(AbstractTrees.getroot(allcells(population)))
        throw(ArgumentError(
            "population has $nroots independent roots (it is a forest), but " *
            "sampling needs one — sample each root's tree separately"))
    end
    return sample_leaves(root, n; seed = seed, replicate = replicate)
end
