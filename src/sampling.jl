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

# ── Declaring what to produce ────────────────────────────────────────────────

"""
    SamplingSpec(; sizes = Int[], replicates = 1, retain_full = true)
    SamplingSpec(n::Int)
    SamplingSpec(sizes::Vector{Int})

What [`sample_trees`](@ref) should produce from one finished tree.

| intent | spec |
|---|---|
| full data only | `SamplingSpec()` |
| one sample of size `n` only | `SamplingSpec(sizes = [n], retain_full = false)` |
| full data plus several sizes | `SamplingSpec([1000, 100])` |

# Keyword arguments
- `sizes` — sample sizes to draw. Empty draws nothing. No duplicates: use
  `replicates` for repeated draws at the same size.
- `replicates` — independent draws per size, each with its own derived seed.
- `retain_full` — whether the returned [`SampledTrees`](@ref) carries the full
  tree.

!!! note "`retain_full = false` does not free memory by itself"
    It controls only what the returned bundle holds. It cannot release the
    caller's own reference to the tree, and it cannot release the tree reachable
    from a `Population` — every alive leaf holds a `parent` chain back to the
    founder, so while the population is alive the whole tree is alive. To
    actually bound memory across a sweep, drop both yourself:

    ```julia
    out = sample_trees(pop, SamplingSpec(sizes = [1000], retain_full = false);
                       seed = myseed)
    pop = nothing
    GC.gc()
    ```
"""
struct SamplingSpec
    sizes::Vector{Int}
    replicates::Int
    retain_full::Bool

    # Inner constructor: this is the ONLY way to build a `SamplingSpec`, so there is
    # no unvalidated path in — Julia would otherwise still expose the auto-generated
    # default positional constructor, which skips every check below.
    function SamplingSpec(sizes::AbstractVector{<:Integer},
                          replicates::Integer,
                          retain_full::Bool)
        sizes = Int[Int(n) for n in sizes]
        allunique(sizes) || throw(ArgumentError(
            "SamplingSpec: duplicate sample sizes in $sizes — use `replicates` for " *
            "repeated draws at the same size"))
        all(>=(1), sizes) || throw(ArgumentError(
            "SamplingSpec: every sample size must be >= 1, got $sizes"))
        replicates >= 1 || throw(ArgumentError(
            "SamplingSpec: replicates must be >= 1, got $replicates"))
        return new(sizes, Int(replicates), retain_full)
    end
end

function SamplingSpec(; sizes::AbstractVector{<:Integer} = Int[],
                        replicates::Integer = 1,
                        retain_full::Bool = true)
    return SamplingSpec(sizes, replicates, retain_full)
end

SamplingSpec(n::Integer) = SamplingSpec(sizes = [n])
SamplingSpec(sizes::AbstractVector{<:Integer}) = SamplingSpec(sizes = sizes)

"""
    SampledTrees

Result of [`sample_trees`](@ref): the full tree (or `nothing` when
`retain_full = false`) and every draw requested by the spec.

`samples` is ordered by the spec's `sizes`, and within a size by `replicate`.
"""
struct SampledTrees
    full::Union{BinaryNode{NonMarkovCell}, Nothing}
    samples::Vector{LeafSample}
end

# Per-draw seed. Depends on the base seed, the size and the replicate index, so no
# two draws in one spec share a seed, and each stored seed replays its own draw in
# isolation.
_draw_seed(base::UInt64, n::Int, replicate::Int) = hash((base, n, replicate))

"""
    sample_trees(root, spec; seed) -> SampledTrees
    sample_trees(population, spec; seed) -> SampledTrees

Apply a [`SamplingSpec`](@ref) to one finished tree.

Sizes are validated against the tree, largest first, *before* any draw is made, so
a size larger than the tree fails immediately rather than after minutes of work.
Per-draw seeds are derived from `seed`, the size and the replicate index and are
recorded on each [`LeafSample`](@ref).

Sampling is post-hoc: it applies to a tree that has finished growing. It is
deliberately not part of `NonMarkovBlock` or `MeasurementSpec`, both of which
describe things that happen *during* `simulate!`.
"""
function sample_trees(root::BinaryNode{NonMarkovCell}, spec::SamplingSpec;
                      seed::UInt64)
    N_full = length(collect(Leaves(root)))
    for n in sort(spec.sizes; rev = true)
        n <= N_full || throw(ArgumentError(
            "SamplingSpec asks for n = $n cells but the tree has $N_full leaves"))
    end

    samples = LeafSample[]
    for n in spec.sizes, r in 1:spec.replicates
        push!(samples, sample_leaves(root, n;
                                     seed = _draw_seed(seed, n, r), replicate = r))
    end
    return SampledTrees(spec.retain_full ? root : nothing, samples)
end

function sample_trees(population::Population, spec::SamplingSpec; seed::UInt64)
    root = getsingleroot(allcells(population))
    if isnothing(root)
        nroots = length(AbstractTrees.getroot(allcells(population)))
        throw(ArgumentError(
            "population has $nroots independent roots (it is a forest), but " *
            "sampling needs one — sample each root's tree separately"))
    end
    return sample_trees(root, spec; seed = seed)
end
