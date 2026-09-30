# ── Abstract cell hierarchy ────────────────────────────────────────────────────

abstract type AbstractCell end
abstract type AbstractTreeCell <: AbstractCell end

# ── BinaryNode ─────────────────────────────────────────────────────────────────

"""
    BinaryNode{T}

Basic unit of a binary tree, used to represent cell lineages.

# Fields
- `data::T`
- `parent::Union{Nothing, BinaryNode{T}}`
- `left::Union{Nothing, BinaryNode{T}}`
- `right::Union{Nothing, BinaryNode{T}}`
"""
mutable struct BinaryNode{T}
    data::T
    parent::Union{Nothing, BinaryNode{T}}
    left::Union{Nothing, BinaryNode{T}}
    right::Union{Nothing, BinaryNode{T}}

    function BinaryNode{T}(data, parent=nothing, l=nothing, r=nothing) where T
        new{T}(data, parent, l, r)
    end
end
BinaryNode(data) = BinaryNode{typeof(data)}(data)

"""
    leftchild!(parent::BinaryNode, data)

Create a new `BinaryNode` from `data` and assign it to `parent.left`.
"""
function leftchild!(parent::BinaryNode, data)
    isnothing(parent.left) || error("left child is already assigned")
    node = typeof(parent)(data, parent)
    parent.left = node
end

"""
    rightchild!(parent::BinaryNode, data)

Create a new `BinaryNode` from `data` and assign it to `parent.right`.
"""
function rightchild!(parent::BinaryNode, data)
    isnothing(parent.right) || error("right child is already assigned")
    node = typeof(parent)(data, parent)
    parent.right = node
end

function AbstractTrees.children(node::BinaryNode)
    if isnothing(node.left) && isnothing(node.right)
        ()
    elseif isnothing(node.left) && !isnothing(node.right)
        (node.right,)
    elseif !isnothing(node.left) && isnothing(node.right)
        (node.left,)
    else
        (node.left, node.right)
    end
end

function AbstractTrees.nextsibling(child::BinaryNode)
    isnothing(child.parent) && return nothing
    p = child.parent
    if !isnothing(p.right)
        child === p.right && return nothing
        return p.right
    end
    return nothing
end

function AbstractTrees.prevsibling(child::BinaryNode)
    isnothing(child.parent) && return nothing
    p = child.parent
    if !isnothing(p.left)
        child === p.left && return nothing
        return p.left
    end
    return nothing
end

AbstractTrees.nodevalue(n::BinaryNode) = n.data
AbstractTrees.ParentLinks(::Type{<:BinaryNode}) = StoredParents()
AbstractTrees.parent(n::BinaryNode) = n.parent
AbstractTrees.NodeType(::Type{<:BinaryNode{T}}) where {T} = HasNodeType()
AbstractTrees.nodetype(::Type{<:BinaryNode{T}}) where {T} = BinaryNode{T}

Base.eltype(::Type{<:TreeIterator{BinaryNode{T}}}) where T = BinaryNode{T}
Base.IteratorEltype(::Type{<:TreeIterator{BinaryNode{T}}}) where T = Base.HasEltype()

AbstractTrees.printnode(io::IO, node::BinaryNode) = print(io, node.data)

haschildren(node::BinaryNode) = !isnothing(node.left) || !isnothing(node.right)
isleaf(node::BinaryNode) = isnothing(node.left) && isnothing(node.right)

function Base.show(io::IO, node::BinaryNode)
    show(io, node.data)
end

# ── Fast traversals ────────────────────────────────────────────────────────────
#
# AbstractTrees' generic iterators cost about 2 µs per node on these trees, which made
# every post-hoc statistic slower than the simulation that produced the tree. All
# package code traverses with the explicit stacks below instead. `Leaves` and
# `PreOrderDFS` still work on a `BinaryNode` for users.

"""
    _leaves(root) -> Vector{BinaryNode}

The leaves under `root`, in exactly `AbstractTrees.Leaves(root)` order (left before
right at every node). Sampling draws index into this order, so it must never change.
"""
function _leaves(root::BinaryNode{T}) where T
    leaves = BinaryNode{T}[]
    stack  = BinaryNode{T}[root]
    while !isempty(stack)
        node = pop!(stack)
        if isleaf(node)
            push!(leaves, node)
        else
            isnothing(node.right) || push!(stack, node.right)
            isnothing(node.left)  || push!(stack, node.left)
        end
    end
    return leaves
end

"""
    _preorder(root) -> (nodes, parent_index)

Every node under `root` in pre-order (left before right), with `parent_index[i]` the
position of `nodes[i]`'s parent in `nodes` (`0` for `root`). Children always come after
their parent, so a reverse sweep is a post-order pass and a forward sweep a pre-order
pass — which is how the spectra and the filtered burden avoid recursion.
"""
function _preorder(root::BinaryNode{T}) where T
    nodes  = BinaryNode{T}[]
    parent = Int[]
    stack  = Tuple{BinaryNode{T}, Int}[(root, 0)]
    while !isempty(stack)
        node, p = pop!(stack)
        push!(nodes, node)
        push!(parent, p)
        i = length(nodes)
        isnothing(node.right) || push!(stack, (node.right, i))
        isnothing(node.left)  || push!(stack, (node.left,  i))
    end
    return nodes, parent
end

# Number of leaves below each node of a `_preorder` listing.
function _leafcounts(nodes::Vector{<:BinaryNode}, parent::Vector{Int})
    counts = zeros(Int, length(nodes))
    for i in length(nodes):-1:1
        isleaf(nodes[i]) && (counts[i] = 1)
        parent[i] == 0 || (counts[parent[i]] += counts[i])
    end
    return counts
end

"""
    _roots(nodes) -> Vector{BinaryNode}

The distinct roots of the trees that `nodes` belong to, in order of first appearance.
Climbs from each node but stops at the first node already visited, so the total cost is
the number of distinct ancestors rather than `length(nodes) × depth`.
"""
function _roots(nodes::AbstractVector{BinaryNode{T}}) where T
    visited = Base.IdSet{BinaryNode{T}}()
    roots   = BinaryNode{T}[]
    for start in nodes
        node = start
        while !(node in visited)
            push!(visited, node)
            if isnothing(node.parent)
                push!(roots, node)
                break
            end
            node = node.parent
        end
    end
    return roots
end

"""
    popsize(root::BinaryNode)

Number of living cells (leaves) under `root`.
"""
popsize(root::BinaryNode) = length(_leaves(root))

# ── MRCA helpers ───────────────────────────────────────────────────────────────

AbstractTrees.getroot(nodevec::Vector{BinaryNode{T}}) where T = _roots(nodevec)

"""
    getsingleroot(nodevec)

Return the unique root of `nodevec`, or `nothing` if there is more than one.
"""
function getsingleroot(nodevec::Vector{BinaryNode{T}}) where T
    roots = _roots(nodevec)
    return length(roots) == 1 ? roots[1] : nothing
end

"""
    findMRCA(population)
    findMRCA(nodes)
    findMRCA(node1, node2)

Find the most recent common ancestor, or `nothing` if the nodes lie in different trees.

Relies on ids increasing along every lineage (a parent's id is smaller than its
children's), which holds for every tree `simulate!` builds. A hand-built tree must keep
that ordering too.
"""
function findMRCA end

function findMRCA(node1, node2)
    (isnothing(node1) || isnothing(node2)) && return nothing
    while node1 !== node2
        node1.data.id > node2.data.id && ((node1, node2) = (node2, node1))
        # node2 has the larger id, so it cannot be an ancestor of node1: climb it. A root
        # holds the smallest id of its tree, so a root with the larger id means the two
        # nodes are in different trees.
        isnothing(node2.parent) && return nothing
        node2 = node2.parent
    end
    return node1
end

function findMRCA(nodes::Vector)
    isempty(nodes) && return nothing
    nodes = copy(nodes)
    node1 = pop!(nodes)
    while length(nodes) > 0
        node2 = pop!(nodes)
        node1 = findMRCA(node1, node2)
        isnothing(node1) && return nothing   # forest: no shared ancestor
    end
    return node1
end


# ── NonMarkovCell ──────────────────────────────────────────────────────────────

"""
    NonMarkovCell <: AbstractTreeCell

Represents a single cell in the non-Markovian birth-death simulation.

# Fields
- `id::Int64` — unique cell identifier; always larger than the parent's id
- `birthtime::Float64` — simulation time at which the cell was born
- `mutations::Int64` — driver mutations acquired at this cell's own birth
- `total_mutations::Int64` — drivers on the whole path from the root to this cell,
  including `mutations`: the parent's `total_mutations + mutations`
- `fitness::Float64` — cumulative fitness (parent fitness updated once per driver mutation)

`total_mutations` is stored so that a cell's burden is a field read, not a walk to the
root. A hand-built tree must keep it consistent: a root has
`total_mutations == mutations`, and every child adds its own `mutations` to its
parent's total. To change a cell's fitness from a hook, use [`set_fitness!`](@ref).
"""
struct NonMarkovCell <: AbstractTreeCell
    id::Int64
    birthtime::Float64
    mutations::Int64
    total_mutations::Int64
    fitness::Float64
end

"""
    set_fitness!(node::BinaryNode{NonMarkovCell}, fitness) -> node

Replace `node`'s cell with an identical one of the given fitness. This is the supported
way for an `on_division` hook to change a daughter: it keeps `id`, `birthtime` and both
mutation counts intact. The new fitness is inherited by all later descendants.
"""
function set_fitness!(node::BinaryNode{NonMarkovCell}, fitness::Real)
    c = node.data
    node.data = NonMarkovCell(c.id, c.birthtime, c.mutations, c.total_mutations, fitness)
    return node
end

id(cellnode::BinaryNode{<:AbstractTreeCell}) = cellnode.data.id

# ── CellEvent ──────────────────────────────────────────────────────────────────

"""
    CellEvent

An entry in the global min-heap event queue.

Defined here rather than in `events.jl` because `Population` carries a
`BinaryMinHeap{CellEvent}` of pending events and therefore needs the type to
already exist. The scheduling logic itself lives in `events.jl`.

# Fields
- `time::Float64` — absolute simulation time at which the event fires
- `node::BinaryNode{NonMarkovCell}` — the cell whose event this is
- `event_type::Symbol` — `:birth` (cell divides) or `:death` (cell dies)
"""
struct CellEvent
    time::Float64
    node::BinaryNode{NonMarkovCell}
    event_type::Symbol
end

Base.isless(a::CellEvent, b::CellEvent) = a.time < b.time

# ── Population ─────────────────────────────────────────────────────────────────

"""
    Population

Holds the set of currently alive cells in a `Dict` keyed by cell id, enabling O(1)
insertion and removal. The simulation time `t` is updated after each event.

# Fields
- `cells` — alive cells, keyed by cell id.
- `t::Float64` — current simulation time.
- `_next_id::Int64` — id allocated to the most recently created cell.
- `_pending` — the queue of already-drawn, not-yet-fired events, one per alive cell,
  or `nothing` if no events have been scheduled yet. `simulate!` carries this across
  calls so that chained blocks resume from the event times already drawn. Use
  [`reset_schedule!`](@ref) to discard it.
"""
mutable struct Population
    cells::Dict{Int64, BinaryNode{NonMarkovCell}}
    t::Float64
    _next_id::Int64
    _pending::Union{Nothing, BinaryMinHeap{CellEvent}}
end

"""
    Population(cells, t, _next_id)

Construct a `Population` with no carried event queue. Preserves the three-argument
positional form used before `_pending` was added.
"""
Population(cells, t, _next_id) = Population(cells, t, _next_id, nothing)

"""
    allcells(population) -> Vector{BinaryNode{NonMarkovCell}}

Return all currently alive cells as a vector.
"""
allcells(population::Population) = collect(values(population.cells))

"""
    getsingleroot(population) -> Union{BinaryNode, Nothing}

The unique root of the population's lineage tree, or `nothing` if it is a forest or
empty. Faster than `getsingleroot(allcells(population))`, with the same result.
"""
getsingleroot(population::Population) = _single_root_or_nothing(population)

# Every leaf of a pruned tree is a living cell. So if the root reached from any one cell
# has as many leaves as the population has cells, those leaves *are* the population and
# that root is the only one. One tree traversal instead of a climb from every cell.
function _single_root_or_nothing(population::Population)
    isempty(population.cells) && return nothing
    root = first(values(population.cells))
    while !isnothing(root.parent)
        root = root.parent
    end
    return popsize(root) == popsize(population) ? root : nothing
end

# The distinct roots of a population's trees, with the single-tree case fast.
function _population_roots(population::Population)
    root = _single_root_or_nothing(population)
    isnothing(root) || return [root]
    return _roots(allcells(population))
end

"""
    popsize(population) -> Int

Return the number of currently alive cells.
"""
popsize(population::Population) = length(population.cells)

function Base.show(io::IO, pop::Population)
    print(io, "Population: $(popsize(pop)) cells (t = $(round(pop.t, digits=3)))")
end

function findMRCA(population::Population)
    # On a single tree the MRCA of all leaves is the first node below the root with
    # two children (or the only leaf); no pairwise climbing needed.
    root = _single_root_or_nothing(population)
    isnothing(root) && return nothing
    node = root
    while true
        isnothing(node.left)  && !isnothing(node.right) && (node = node.right; continue)
        isnothing(node.right) && !isnothing(node.left)  && (node = node.left;  continue)
        return node
    end
end
