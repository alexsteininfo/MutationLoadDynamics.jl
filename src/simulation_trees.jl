# Tree utilities: pruning, roots, MRCAs, division times and lifetimes.

"""
    prune_tree!(cellnode)

Remove `cellnode` from the tree and then every ancestor that becomes childless. Called
when a leaf cell dies without dividing, so the tree keeps only the survivors' ancestry.
"""
function prune_tree!(cellnode::BinaryNode)
    while true
        parent = cellnode.parent
        isnothing(parent) && return
        cellnode.parent = nothing
        if parent.left === cellnode
            parent.left = nothing
        elseif parent.right === cellnode
            parent.right = nothing
        else
            error("dead cell is neither left nor right child of its parent")
        end
        haschildren(parent) && return
        cellnode = parent
    end
end

# ── Cells, roots and MRCAs ────────────────────────────────────────────────────

"""
    alive_cells(population) -> Vector{BinaryNode{NonMarkovCell}}
    alive_cells(root::BinaryNode) -> Vector{BinaryNode{NonMarkovCell}}

The living cells: of the population in increasing id order, or the leaves under `root`
in `Leaves(root)` order. Dead cells are pruned from the tree, so every leaf is alive.
"""
alive_cells(population::Population) = _sorted_cells(population)
alive_cells(root::BinaryNode) = _leaves(root)

"""
    roots(population) -> Vector{BinaryNode}
    roots(nodes::Vector{<:BinaryNode}) -> Vector{BinaryNode}

The distinct roots of the trees the population's cells (or the given nodes) belong to.
A population from `initialize_population(N)` with `N > 1` is a forest of `N` trees.
"""
roots(population::Population) = _population_roots(population)
roots(nodes::AbstractVector{<:BinaryNode}) = _roots(nodes)

"""
    single_root(population) -> Union{BinaryNode, Nothing}
    single_root(nodes::Vector{<:BinaryNode}) -> Union{BinaryNode, Nothing}

The unique root of the population's (or the nodes') lineage tree, or `nothing` if they
form a forest or are empty.
"""
single_root(population::Population) = _single_root_or_nothing(population)
function single_root(nodes::AbstractVector{<:BinaryNode})
    found = _roots(nodes)
    return length(found) == 1 ? only(found) : nothing
end

# Every leaf of a pruned tree is a living cell. So if the root reached from any one cell
# has as many leaves as the population has cells, those leaves *are* the population and
# that root is the only one: one tree traversal instead of a climb from every cell.
function _single_root_or_nothing(population::Population)
    isempty(population.cells) && return nothing
    root = _treeroot(first(values(population.cells)))
    return popsize(root) == popsize(population) ? root : nothing
end

function _population_roots(population::Population)
    root = _single_root_or_nothing(population)
    isnothing(root) || return [root]
    return _roots(collect(values(population.cells)))
end

"""
    find_mrca(node1, node2)
    find_mrca(nodes::Vector)
    find_mrca(population)

The most recent common ancestor, or `nothing` if the nodes lie in different trees.

Relies on ids increasing along every lineage (a parent's id is smaller than its
children's), which holds for every tree `simulate!` builds. A hand-built tree must keep
that ordering too.
"""
function find_mrca(node1::BinaryNode, node2::BinaryNode)
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

function find_mrca(nodes::AbstractVector{<:BinaryNode})
    isempty(nodes) && return nothing
    mrca = first(nodes)
    for node in Iterators.drop(nodes, 1)
        mrca = find_mrca(mrca, node)
        isnothing(mrca) && return nothing   # forest: no shared ancestor
    end
    return mrca
end

function find_mrca(population::Population)
    # On a single tree the MRCA of all leaves is the first node below the root with two
    # children (or the only leaf): everything above it is a chain of unary nodes.
    node = _single_root_or_nothing(population)
    isnothing(node) && return nothing
    while true
        if isnothing(node.left) && !isnothing(node.right)
            node = node.right
        elseif isnothing(node.right) && !isnothing(node.left)
            node = node.left
        else
            return node
        end
    end
end

# ── Times ─────────────────────────────────────────────────────────────────────

"""
    division_time(node::BinaryNode) -> Union{Float64, Nothing}

The time at which the cell divided (its daughters' birthtime), or `nothing` if it is
still alive (a leaf).
"""
function division_time(node::BinaryNode)
    isnothing(node.left) || return node.left.data.birthtime
    isnothing(node.right) || return node.right.data.birthtime
    return nothing
end

"""
    last_division_time(root::BinaryNode) -> Float64

The birthtime of the most recently born leaf under `root`, i.e. the time of the last
division in that tree. This is earlier than the population's clock `pop.t` if the run
ended on a death.
"""
last_division_time(root::BinaryNode) = maximum(leaf.data.birthtime for leaf in _leaves(root))

"""
    cell_lifetime(node::BinaryNode, tnow) -> Float64

Time from the cell's birth to its division, or to `tnow` if it is still alive. Pass the
population's clock `pop.t` as `tnow` for "how old is this cell right now".
"""
function cell_lifetime(node::BinaryNode, tnow::Real)
    t_end = division_time(node)
    return (isnothing(t_end) ? tnow : t_end) - node.data.birthtime
end

"""
    cell_lifetimes(root; include_alive = false, tnow = last_division_time(root))

The lifetime of every cell in the tree under `root`, in pre-order. By default only cells
that have divided contribute. With `include_alive = true`, living cells contribute their
age at `tnow`, which is a right-censored lifetime rather than a completed one.
"""
function cell_lifetimes(root::BinaryNode; include_alive::Bool = false,
                        tnow::Real = last_division_time(root))
    nodes, _ = _preorder(root)
    include_alive && return Float64[cell_lifetime(n, tnow) for n in nodes]
    return Float64[division_time(n) - n.data.birthtime for n in nodes if haschildren(n)]
end
