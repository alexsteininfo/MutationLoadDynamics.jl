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

haschildren(node::BinaryNode) = length(children(node)) != 0

function Base.show(io::IO, node::BinaryNode)
    show(io, node.data)
end

"""
    popsize(root::BinaryNode)

Count alive descendant leaf cells under `root`.
"""
function popsize(root::BinaryNode)
    N = 0
    for l in Leaves(root)
        if isalive(nodevalue(l))
            N += 1
        end
    end
    return N
end

# ── MRCA helpers ───────────────────────────────────────────────────────────────

AbstractTrees.getroot(nodevec::Vector{BinaryNode{T}}) where T =
    collect(Set(getroot.(nodevec)))

"""
    getsingleroot(nodevec)

Return the unique root of `nodevec`, or `nothing` if there is more than one.
"""
function getsingleroot(nodevec::Vector{BinaryNode{T}}) where T
    roots = AbstractTrees.getroot(nodevec)
    return length(roots) == 1 ? roots[1] : nothing
end

"""
    findMRCA(population)
    findMRCA(nodes)
    findMRCA(node1, node2)

Find the most recent common ancestor.
"""
function findMRCA end

function findMRCA(node1, node2)
    (isnothing(node1) || isnothing(node2)) && return nothing
    node1 == node2 && return node1
    while true
        if node1.data.id > node2.data.id
            node1, node2 = node2, node1
        end
        if node1.parent == node2.parent
            return node1.parent
        elseif isnothing(node1.parent) && isnothing(node2.parent)
            return nothing
        elseif node1 == node2.parent
            return node1
        else
            node2 = node2.parent
        end
    end
end

function findMRCA(nodes::Vector)
    nodes = copy(nodes)
    node1 = pop!(nodes)
    while length(nodes) > 0
        node2 = pop!(nodes)
        node1 = findMRCA(node1, node2)
        isnothing(node1) && return nothing   # forest: no shared ancestor
    end
    return node1
end

function findMRCA(population)
    return findMRCA(allcells(population))
end

# ── NonMarkovCell ──────────────────────────────────────────────────────────────

"""
    NonMarkovCell <: AbstractTreeCell

Represents a single cell in the non-Markovian birth-death simulation.

# Fields
- `id::Int64` — unique cell identifier
- `birthtime::Float64` — simulation time at which the cell was born
- `mutations::Int64` — number of new driver mutations acquired at this cell's birth
- `fitness::Float64` — cumulative fitness (parent fitness updated once per driver mutation)
"""
struct NonMarkovCell <: AbstractTreeCell
    id::Int64
    birthtime::Float64
    mutations::Int64
    fitness::Float64
end

id(cellnode::BinaryNode{<:AbstractTreeCell}) = cellnode.data.id

# ── Population ─────────────────────────────────────────────────────────────────

"""
    Population

Holds the set of currently alive cells in a `Dict` keyed by cell id, enabling O(1)
insertion and removal. The simulation time `t` is updated after each event.
"""
mutable struct Population
    cells::Dict{Int64, BinaryNode{NonMarkovCell}}
    t::Float64
    _next_id::Int64
end

"""
    allcells(population) -> Vector{BinaryNode{NonMarkovCell}}

Return all currently alive cells as a vector.
"""
allcells(population::Population) = collect(values(population.cells))

"""
    popsize(population) -> Int

Return the number of currently alive cells.
"""
popsize(population::Population) = length(population.cells)

function Base.show(io::IO, pop::Population)
    print(io, "Population: $(popsize(pop)) cells (t = $(round(pop.t, digits=3)))")
end
